use crate::daemon::*;
use crate::dsp::*;
use crate::orbit::*;
use crate::tui::*;
use chrono::{DateTime, Datelike, Timelike, Utc};
use nalgebra::{Matrix1, Matrix2, RowVector2, Vector2};
use num_complex::Complex;
use rustfft::FftPlanner;
use sgp4::Elements;
use std::collections::VecDeque;
use std::io::{self, Read, Write};

pub struct ClockEkf {
    pub x: Vector2<f64>, // [phase_offset_seconds, freq_drift_ppm]
    pub p: Matrix2<f64>, // error covariance matrix
    pub q_phase: f64,    // phase process noise
    pub q_freq: f64,     // freq drift process noise
    pub r_meas: f64,     // measurement noise
    pub r_freq: f64,     // frequency measurement noise
}

impl ClockEkf {
    pub fn new() -> Self {
        Self {
            x: Vector2::zeros(),
            p: Matrix2::new(1e-4, 0.0, 0.0, 1e-2), // 10ms phase variance, 0.1 PPM freq variance
            q_phase: 1e-12,                        // phase noise (s^2 / s)
            q_freq: 1e-14,                         // frequency walk noise (PPM^2 / s)
            r_meas: 1e-8,                          // measurement noise (s^2)
            r_freq: 1e-4,
        }
    }

    pub fn predict(&mut self, dt: f64) {
        if dt <= 0.0 {
            return;
        }

        let dt_scale = dt * 1e-6;
        let f = Matrix2::new(1.0, dt_scale, 0.0, 1.0);
        let q = Matrix2::new(self.q_phase * dt, 0.0, 0.0, self.q_freq * dt);

        self.x = f * self.x;
        self.p = f * self.p * f.transpose() + q;
    }

    pub fn update(&mut self, z_offset: f64, z_freq_ppm: f64) {
        let z = Vector2::new(z_offset, z_freq_ppm);
        let y = z - self.x;
        let r = Matrix2::new(self.r_meas, 0.0, 0.0, self.r_freq);
        let s = self.p + r;
        let s_inv = s.try_inverse().unwrap_or_else(|| Matrix2::identity());
        let k = self.p * s_inv;
        self.x += k * y;
        let i = Matrix2::identity();
        let a = i - k;
        self.p = a * self.p * a.transpose() + k * r * k.transpose();
    }
}

use nalgebra::{Matrix3, RowVector3, Vector3};

#[derive(Clone)]
pub struct CarrierPllEkf {
    pub x: Vector3<f64>, // [phase (rad), freq (rad/s), chirp_rate (rad/s^2)]
    pub p: Matrix3<f64>, // error covariance matrix
    pub q_phase: f64,
    pub q_freq: f64,
    pub q_chirp: f64,
    pub r_meas: f64,
    pub ts: f64,
    pub lock_metric: f64,
    pub is_locked: bool,
    pub modulation: Modulation,
    pub envelope_ema: f64,
    pub adaptive_ekf: bool,
    pub dual_lock: bool,
    pub pr_window: VecDeque<(f64, f64)>,
    pub pr_sum_abs_i: f64,
    pub pr_sum_q_sq: f64,
    pub convergence_guard: usize, // Samples remaining before lock_metric can trigger unlock
}

impl CarrierPllEkf {
    pub fn new(fs: f64, modulation: Modulation) -> Self {
        let mut p = Matrix3::zeros();
        p[(0, 0)] = 1.0;
        p[(1, 1)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
        p[(2, 2)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);

        Self {
            x: Vector3::zeros(),
            p,
            q_phase: 1e-1,
            q_freq: 5e3,
            q_chirp: 1e4,
            r_meas: 10.0,
            ts: 1.0 / fs,
            lock_metric: 0.0,
            is_locked: false,
            modulation,
            envelope_ema: 1.0,
            adaptive_ekf: true,
            dual_lock: true,
            pr_window: VecDeque::new(),
            pr_sum_abs_i: 0.0,
            pr_sum_q_sq: 0.0,
            convergence_guard: 0,
        }
    }

    pub fn reset(&mut self, initial_phase: f64, initial_freq_hz: f64, initial_chirp_hz_s: f64) {
        self.x = Vector3::new(
            initial_phase,
            2.0 * std::f64::consts::PI * initial_freq_hz,
            2.0 * std::f64::consts::PI * initial_chirp_hz_s,
        );

        self.p = Matrix3::zeros();
        self.p[(0, 0)] = 1.0;
        self.p[(1, 1)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
        self.p[(2, 2)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);

        self.lock_metric = 0.5;
        self.is_locked = true;
        self.envelope_ema = 1.0;
        self.pr_window.clear();
        self.pr_sum_abs_i = 0.0;
        self.pr_sum_q_sq = 0.0;
        // Grace period: suppress unlock checks for 2048 samples so EKF can converge
        self.convergence_guard = 2048;
    }

    pub fn predict(&mut self) {
        let dt = self.ts;
        let dt2 = 0.5 * dt * dt;

        let f = Matrix3::new(1.0, dt, dt2, 0.0, 1.0, dt, 0.0, 0.0, 1.0);

        self.x = f * self.x;

        let (limit, half_limit) = match self.modulation {
            Modulation::Carrier | Modulation::Bpsk => (2.0 * std::f64::consts::PI, std::f64::consts::PI),
            Modulation::Qpsk => (std::f64::consts::PI / 2.0, std::f64::consts::PI / 4.0),
        };
        self.x[0] = (self.x[0] + half_limit).rem_euclid(limit) - half_limit;

        let s_m = if self.adaptive_ekf {
            1.0 - 0.9 * self.lock_metric.clamp(0.0, 1.0)
        } else {
            1.0
        };

        let q = Matrix3::new(
            self.q_phase * s_m * dt,
            0.0,
            0.0,
            0.0,
            self.q_freq * s_m * dt,
            0.0,
            0.0,
            0.0,
            self.q_chirp * s_m * dt,
        );

        self.p = f * self.p * f.transpose() + q;

        // Force covariance matrix symmetry to prevent numerical collapse
        self.p = (self.p + self.p.transpose()) * 0.5;
    }

    pub fn update(&mut self, sample: Complex<f32>) {
        let sample_to_use = match self.modulation {
            Modulation::Bpsk => sample * sample,
            _ => sample,
        };
        let theta_pred = self.x[0];

        let cos_theta = (-theta_pred).cos();
        let sin_theta = (-theta_pred).sin();
        let derotated_re = sample_to_use.re as f64 * cos_theta - sample_to_use.im as f64 * sin_theta;
        let derotated_im = sample_to_use.re as f64 * sin_theta + sample_to_use.im as f64 * cos_theta;

        let amp = (sample_to_use.re as f64).hypot(sample_to_use.im as f64);

        // Dynamic SNR-based measurement noise scaling & envelope check
        let alpha = 0.005;
        self.envelope_ema = (1.0 - alpha) * self.envelope_ema + alpha * amp;

        let rel_amp = if self.envelope_ema > 1e-6 {
            amp / self.envelope_ema
        } else {
            1.0
        };

        let m = self.lock_metric.clamp(1e-3, 0.999);
        let snr_factor = (1.0 - m) / m;
        let snr_factor = snr_factor.clamp(0.1, 100.0);
        let fade_factor = 1.0 / (rel_amp * rel_amp).max(1e-4);
        let r_effective = self.r_meas * snr_factor * fade_factor;

        let (z, norm_re) = match self.modulation {
            Modulation::Carrier | Modulation::Bpsk => {
                let z = derotated_im.atan2(derotated_re);
                let norm_re = z.cos();
                (z, norm_re)
            }
            Modulation::Qpsk => {
                let i2 = derotated_re * derotated_re - derotated_im * derotated_im;
                let q2 = 2.0 * derotated_re * derotated_im;
                let im_y4 = 2.0 * i2 * q2;
                let re_y4 = i2 * i2 - q2 * q2;
                let z = 0.25 * im_y4.atan2(re_y4);
                let norm_re = (4.0 * z).cos();
                (z, norm_re)
            }
        };

        // Folded components for Power Ratio
        let (i_n, q_n) = match self.modulation {
            Modulation::Carrier | Modulation::Bpsk => (derotated_re, derotated_im),
            Modulation::Qpsk => {
                let sqrt2 = 2.0f64.sqrt();
                (
                    (derotated_re.abs() + derotated_im.abs()) / sqrt2,
                    (derotated_im.abs() - derotated_re.abs()) / sqrt2,
                )
            }
        };

        let abs_i = i_n.abs();
        let q_sq = q_n * q_n;

        self.pr_sum_abs_i += abs_i;
        self.pr_sum_q_sq += q_sq;
        self.pr_window.push_back((abs_i, q_sq));

        if self.pr_window.len() > 1024
            && let Some((old_abs_i, old_q_sq)) = self.pr_window.pop_front()
        {
            self.pr_sum_abs_i -= old_abs_i;
            self.pr_sum_q_sq -= old_q_sq;
        }
        if self.pr_sum_abs_i < 0.0 {
            self.pr_sum_abs_i = 0.0;
        }
        if self.pr_sum_q_sq < 0.0 {
            self.pr_sum_q_sq = 0.0;
        }

        let pr = if self.pr_sum_q_sq > 1e-12 {
            (self.pr_sum_abs_i * self.pr_sum_abs_i) / self.pr_sum_q_sq
        } else {
            f64::MAX
        };

        let h = RowVector3::new(1.0, 0.0, 0.0);

        let s_val = (h * self.p * h.transpose())[(0, 0)] + r_effective;
        if s_val.abs() < 1e-12 {
            return;
        }

        let k = (self.p * h.transpose()) / s_val;

        self.x += k * z;

        let (limit, half_limit) = match self.modulation {
            Modulation::Carrier | Modulation::Bpsk => (2.0 * std::f64::consts::PI, std::f64::consts::PI),
            Modulation::Qpsk => (std::f64::consts::PI / 2.0, std::f64::consts::PI / 4.0),
        };
        self.x[0] = (self.x[0] + half_limit).rem_euclid(limit) - half_limit;

        let i = Matrix3::identity();
        let a = i - k * h;
        let r_mat = Matrix1::new(r_effective);
        self.p = a * self.p * a.transpose() + k * r_mat * k.transpose();

        let beta = 0.001;
        self.lock_metric = (1.0 - beta) * self.lock_metric + beta * norm_re;
        // Only check unlock thresholds after convergence grace period expires
        if self.convergence_guard > 0 {
            self.convergence_guard -= 1;
        } else if self.is_locked {
            if self.dual_lock {
                if self.lock_metric < 0.02 {
                    self.is_locked = false;
                } else if self.lock_metric < 0.1 && pr <= 0.8 * 1024.0 {
                    self.is_locked = false;
                }
            } else if self.lock_metric < 0.1 {
                self.is_locked = false;
            }
        }
    }
}

#[derive(Clone)]
pub struct EkfTrackingBank {
    pub trackers: [CarrierPllEkf; 3],
    pub gardner_loops: [GardnerLoop; 3],
    pub active_idx: Option<usize>,
    pub in_fade: bool,
    pub fade_counter: usize,
    pub max_fade_steps: usize,
    pub terminated_in_fade: bool,
    pub spur_dwell_counter: usize,
}

impl EkfTrackingBank {
    pub fn compute_tracker_discrepancy(&self) -> f32 {
        let freq_diffs = self.compute_tracker_frequency_diffs();

        let w0 = if self.trackers[0].is_locked {
            self.trackers[0].lock_metric
        } else {
            0.0
        };
        let w1 = if self.trackers[1].is_locked {
            self.trackers[1].lock_metric
        } else {
            0.0
        };
        let w2 = if self.trackers[2].is_locked {
            self.trackers[2].lock_metric
        } else {
            0.0
        };

        let w01 = (w0 * w1) as f32;
        let w12 = (w1 * w2) as f32;
        let w20 = (w2 * w0) as f32;

        let sum_w = w01 + w12 + w20;
        if sum_w < 1e-5 || sum_w.is_nan() {
            return 0.0;
        }

        let val =
            (w01 * freq_diffs[0].abs() + w12 * freq_diffs[1].abs() + w20 * freq_diffs[2].abs())
                / sum_w;
        if val.is_nan() { 0.0 } else { val }
    }

    pub fn collect_tracker_features(&self) -> Vec<f32> {
        let mut features = Vec::with_capacity(12);
        for tracker in &self.trackers {
            features.push(tracker.lock_metric as f32);
            features.push(if tracker.is_locked { 1.0 } else { 0.0 });
            features.push(tracker.p[(0, 0)] as f32);
            features.push(tracker.p[(1, 1)] as f32);
        }
        features
    }

    pub fn compute_tracker_frequency_diffs(&self) -> Vec<f32> {
        let scale = if self.trackers[0].modulation == Modulation::Bpsk { 2.0 } else { 1.0 };
        let f0 = ((self.trackers[0].x[1] / scale) / (2.0 * std::f64::consts::PI)) as f32;
        let f1 = ((self.trackers[1].x[1] / scale) / (2.0 * std::f64::consts::PI)) as f32;
        let f2 = ((self.trackers[2].x[1] / scale) / (2.0 * std::f64::consts::PI)) as f32;

        vec![
            f1 - f0, // tracker 0 vs 1
            f2 - f1, // tracker 1 vs 2
            f0 - f2, // tracker 2 vs 0
        ]
    }
}

pub use crate::dsp::{DemodChannel, TelemetryUpdate, process_pipeline_parallel};

pub struct ChannelAllocator {
    pub max_channels: usize,
    pub active_channels: std::collections::HashMap<String, usize>, // sat_name -> channel_index
    pub channels: Vec<Option<String>>, // channel_index -> Option<sat_name>
}

impl ChannelAllocator {
    pub fn new(max_channels: usize) -> Self {
        Self {
            max_channels,
            active_channels: std::collections::HashMap::new(),
            channels: vec![None; max_channels],
        }
    }

    pub fn handle_aos(&mut self, sat_name: &str, tle_data: Option<&str>) -> bool {
        // TLE validity check:
        if let Some(tle) = tle_data {
            if tle.contains("corrupt") || tle.contains("invalid") || tle.is_empty() {
                return false;
            }
        }

        if self.active_channels.contains_key(sat_name) {
            return true;
        }

        if let Some(slot) = self.channels.iter().position(|c| c.is_none()) {
            if slot < self.max_channels {
                self.channels[slot] = Some(sat_name.to_string());
                self.active_channels.insert(sat_name.to_string(), slot);
                return true;
            }
        }
        false
    }

    pub fn handle_los(&mut self, sat_name: &str) {
        if let Some(slot) = self.active_channels.remove(sat_name) {
            if slot < self.channels.len() {
                self.channels[slot] = None;
            }
        }
    }

    pub fn active_count(&self) -> usize {
        self.active_channels.len()
    }
}

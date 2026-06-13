use crate::dsp::*;
use nalgebra::{Matrix2, Vector2};
use num_complex::Complex;
use std::collections::VecDeque;

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
            r_meas: 1e-10,                         // measurement noise (s^2), optimized for microsecond-level fits
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

    pub fn update_1d(&mut self, z_offset: f64) {
        let y = z_offset - self.x[0];
        let s = self.p[(0, 0)] + self.r_meas;
        if s.abs() >= 1e-12 {
            let k = self.p.column(0) / s;
            self.x += k * y;
            let k0 = k[0];
            let k1 = k[1];
            let a = Matrix2::new(1.0 - k0, 0.0, -k1, 1.0);
            self.p = a * self.p * a.transpose();
            self.p[(0, 0)] += k0 * k0 * self.r_meas;
            self.p[(0, 1)] += k0 * k1 * self.r_meas;
            self.p[(1, 0)] += k1 * k0 * self.r_meas;
            self.p[(1, 1)] += k1 * k1 * self.r_meas;
        }
    }

    pub fn update(&mut self, z_offset: f64, z_freq_ppm: f64) {
        let z = Vector2::new(z_offset, z_freq_ppm);
        let y = z - self.x;
        let r = Matrix2::new(self.r_meas, 0.0, 0.0, self.r_freq);
        let s = self.p + r;
        // If S is singular, skip the update entirely to prevent state jumps.
        let Some(s_inv) = s.try_inverse() else { return; };
        let k = self.p * s_inv;
        self.x += k * y;
        let i = Matrix2::identity();
        let a = i - k;
        self.p = a * self.p * a.transpose() + k * r * k.transpose();
    }
}

use nalgebra::{Matrix6, Vector6};

#[derive(Clone)]
pub struct CarrierPllEkf {
    pub x: Vector6<f64>, // [phase1, freq1, chirp1, phase2, freq2, chirp2]
    pub p: Matrix6<f64>, // error covariance matrix
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
    pub frequency_ratio: f64,
    pub raw_amp: [f64; 2],
    pub unwrapped_phase1: f64,
    pub unwrapped_phase2: f64,
}

impl CarrierPllEkf {
    pub fn new(fs: f64, modulation: Modulation) -> Self {
        let mut p = Matrix6::zeros();
        p[(0, 0)] = 1.0;
        p[(1, 1)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
        p[(2, 2)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);
        p[(3, 3)] = 1.0;
        p[(4, 4)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
        p[(5, 5)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);

        Self {
            x: Vector6::zeros(),
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
            frequency_ratio: 1.0,
            raw_amp: [0.0, 0.0],
            unwrapped_phase1: 0.0,
            unwrapped_phase2: 0.0,
        }
    }

    pub fn set_frequency_ratio(&mut self, ratio: f64) {
        self.frequency_ratio = ratio;
    }

    pub fn reset(&mut self, initial_phase: f64, initial_freq_hz: f64, initial_chirp_hz_s: f64) {
        let r = self.frequency_ratio;
        self.x = Vector6::new(
            initial_phase,
            2.0 * std::f64::consts::PI * initial_freq_hz,
            2.0 * std::f64::consts::PI * initial_chirp_hz_s,
            initial_phase,
            2.0 * std::f64::consts::PI * initial_freq_hz * r,
            2.0 * std::f64::consts::PI * initial_chirp_hz_s * r,
        );

        self.p = Matrix6::zeros();
        self.p[(0, 0)] = 1.0;
        self.p[(1, 1)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
        self.p[(2, 2)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);
        self.p[(3, 3)] = 1.0;
        self.p[(4, 4)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
        self.p[(5, 5)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);

        self.lock_metric = 0.5;
        self.is_locked = true;
        self.envelope_ema = 1.0;
        self.pr_window.clear();
        self.pr_sum_abs_i = 0.0;
        self.pr_sum_q_sq = 0.0;
        // Grace period: suppress unlock checks for 2048 samples so EKF can converge
        self.convergence_guard = 2048;
        self.raw_amp = [0.0, 0.0];
        self.unwrapped_phase1 = initial_phase;
        self.unwrapped_phase2 = initial_phase;
    }

    pub fn predict(&mut self) {
        let dt = self.ts;
        let dt2 = 0.5 * dt * dt;

        // Accumulate unwrapped phase changes
        let dp1 = self.x[1] * dt + 0.5 * self.x[2] * dt * dt;
        let dp2 = self.x[4] * dt + 0.5 * self.x[5] * dt * dt;
        self.unwrapped_phase1 += dp1;
        self.unwrapped_phase2 += dp2;

        let mut f = Matrix6::identity();
        f[(0, 1)] = dt;
        f[(0, 2)] = dt2;
        f[(1, 2)] = dt;
        f[(3, 4)] = dt;
        f[(3, 5)] = dt2;
        f[(4, 5)] = dt;

        self.x = f * self.x;

        let s_m = if self.adaptive_ekf {
            1.0 - 0.9 * self.lock_metric.clamp(0.0, 1.0)
        } else {
            1.0
        };

        let mut q = Matrix6::zeros();
        let q_p = self.q_phase * s_m * dt;
        let q_f = self.q_freq * s_m * dt;
        let q_c = self.q_chirp * s_m * dt;

        q[(0, 0)] = q_p;
        q[(1, 1)] = q_f;
        q[(2, 2)] = q_c;
        q[(3, 3)] = q_p;

        if (self.frequency_ratio - 1.0).abs() < 1e-6 {
            q[(4, 4)] = q_f;
            q[(5, 5)] = q_c;
        } else {
            let r = self.frequency_ratio;
            q[(4, 4)] = q_f * r * r;
            q[(5, 5)] = q_c * r * r;
            // Coupled process noise terms
            q[(1, 4)] = q_f * r;
            q[(4, 1)] = q_f * r;
            q[(2, 5)] = q_c * r;
            q[(5, 2)] = q_c * r;
        }

        self.p = f * self.p * f.transpose() + q;

        // Force covariance matrix symmetry to prevent numerical collapse
        self.p = (self.p + self.p.transpose()) * 0.5;
    }

    fn update_channel(&mut self, chan: usize, sample: Complex<f32>) -> Option<(f64, f64)> {
        let sample_to_use = match self.modulation {
            Modulation::Bpsk => sample * sample,
            _ => sample,
        };
        let idx = chan * 3;
        let theta_pred = self.x[idx];

        let cos_theta = (-theta_pred).cos();
        let sin_theta = (-theta_pred).sin();
        let derotated_re = sample_to_use.re as f64 * cos_theta - sample_to_use.im as f64 * sin_theta;
        let derotated_im = sample_to_use.re as f64 * sin_theta + sample_to_use.im as f64 * cos_theta;

        let amp = if self.raw_amp[chan] > 0.0 {
            self.raw_amp[chan]
        } else {
            (sample_to_use.re as f64).hypot(sample_to_use.im as f64)
        };

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

        let (limit, half_limit) = match self.modulation {
            Modulation::Carrier | Modulation::Bpsk => (2.0 * std::f64::consts::PI, std::f64::consts::PI),
            Modulation::Qpsk => (std::f64::consts::PI / 2.0, std::f64::consts::PI / 4.0),
        };
        let z = (z + half_limit).rem_euclid(limit) - half_limit;

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

        if chan == 0 {
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
        }

        let s_val = self.p[(idx, idx)] + r_effective;
        if s_val.abs() >= 1e-12 {
            let mut k = [0.0; 6];
            for r in 0..6 {
                k[r] = self.p[(r, idx)] / s_val;
            }
            for r in 0..6 {
                self.x[r] += k[r] * z;
            }
            // Accumulate unwrapped phase correction
            let delta_phase = k[idx] * z;
            if chan == 0 {
                self.unwrapped_phase1 += delta_phase;
            } else {
                self.unwrapped_phase2 += delta_phase;
            }
            self.x[idx] = (self.x[idx] + half_limit).rem_euclid(limit) - half_limit;

            if self.x.iter().any(|v| v.is_nan()) {
                self.x = Vector6::zeros();
                self.p = Matrix6::zeros();
                self.p[(0, 0)] = 1.0;
                self.p[(1, 1)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
                self.p[(2, 2)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);
                self.p[(3, 3)] = 1.0;
                self.p[(4, 4)] = (2.0 * std::f64::consts::PI * 100.0).powi(2);
                self.p[(5, 5)] = (2.0 * std::f64::consts::PI * 50.0).powi(2);
                self.is_locked = false;
                self.lock_metric = 0.0;
                return None;
            } else {
                let max_freq = 2.0 * std::f64::consts::PI * 50000.0; // 50 kHz max
                let max_chirp = 2.0 * std::f64::consts::PI * 1000.0;  // 1 kHz/s max
                self.x[idx + 1] = self.x[idx + 1].clamp(-max_freq, max_freq);
                self.x[idx + 2] = self.x[idx + 2].clamp(-max_chirp, max_chirp);

                // Minimum covariance floor to prevent overconfidence.
                for i in 0..6 {
                    self.p[(i, i)] = self.p[(i, i)].max(1e-15);
                }

                // Manually unrolled Joseph form covariance update
                // temp[r, c] = P[r, c] - k[r] * P[idx, c]
                // P_new[r, c] = temp[r, c] - temp[r, idx] * k[c] + k[r] * k[c] * r_effective
                let mut temp = [[0.0; 6]; 6];
                for r in 0..6 {
                    let kr = k[r];
                    for c in 0..6 {
                        temp[r][c] = self.p[(r, c)] - kr * self.p[(idx, c)];
                    }
                }
                for r in 0..6 {
                    let kr = k[r];
                    for c in 0..6 {
                        self.p[(r, c)] = temp[r][c] - temp[r][idx] * k[c] + kr * k[c] * r_effective;
                    }
                }
            }
        }

        let pr = if self.pr_sum_q_sq > 1e-12 {
            (self.pr_sum_abs_i * self.pr_sum_abs_i) / self.pr_sum_q_sq
        } else {
            f64::MAX
        };

        Some((norm_re, pr))
    }

    pub fn update(&mut self, sample: Complex<f32>) {
        let Some((norm_re, pr)) = self.update_channel(0, sample) else {
            self.raw_amp = [0.0, 0.0];
            return;
        };
        self.raw_amp = [0.0, 0.0];

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

    pub fn update_dual(&mut self, sample1: Complex<f32>, sample2: Complex<f32>) {
        let Some((norm_re, pr)) = self.update_channel(0, sample1) else {
            self.raw_amp = [0.0, 0.0];
            return;
        };
        if self.update_channel(1, sample2).is_none() {
            self.raw_amp = [0.0, 0.0];
            return;
        }
        self.raw_amp = [0.0, 0.0];

        let beta = 0.001;
        self.lock_metric = (1.0 - beta) * self.lock_metric + beta * norm_re;
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

pub struct AppletonHartreeDispersion;

impl AppletonHartreeDispersion {
    pub fn cancel(f1: f64, f2: f64, fd1: f64, fd2: f64) -> f64 {
        let f1_sq = f1 * f1;
        let f2_sq = f2 * f2;
        let diff = f1_sq - f2_sq;
        if diff.abs() < 1e-6 {
            fd1
        } else {
            (f1_sq * fd1 - f2_sq * fd2) / diff
        }
    }
}

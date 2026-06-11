use num_complex::Complex;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

#[derive(Serialize, Deserialize, Debug)]
pub struct CalibrationData {
    pub df0: f64,
    pub timestamp: chrono::DateTime<chrono::Utc>,
}

#[derive(Debug)]
pub struct Args {
    pub leodo: bool,
    pub output_dir: String,
}

pub fn get_process_rss_mb() -> f64 {
    0.0
}

#[path = "../src/dsp.rs"]
pub mod dsp;

#[path = "../src/ekf.rs"]
pub mod ekf;

#[path = "../src/daemon.rs"]
pub mod daemon;

#[path = "../src/orbit_solver.rs"]
pub mod orbit_solver;

#[path = "../src/orbit.rs"]
pub mod orbit;

#[path = "../src/tui.rs"]
pub mod tui;

use dsp::Modulation;
use ekf::CarrierPllEkf;

struct SimpleRng {
    state: u32,
}

impl SimpleRng {
    fn new(seed: u32) -> Self {
        Self { state: seed }
    }

    fn next_f32(&mut self) -> f32 {
        self.state = self.state.wrapping_mul(1664525).wrapping_add(1013904223);
        self.state as f32 / u32::MAX as f32
    }

    fn next_gaussian(&mut self) -> (f32, f32) {
        let u1 = self.next_f32().max(1e-6);
        let u2 = self.next_f32();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f32::consts::PI as f32 * u2;
        (r * theta.cos(), r * theta.sin())
    }
}

// ----------------------------------------------------
// VERIFICATION TESTS
// ----------------------------------------------------

#[test]
fn test_adaptive_ekf_process_noise_fixed() {
    let fs = 1000.0;

    // Test Case 1: adaptive_ekf = false
    let mut tracker_fixed = CarrierPllEkf::new(fs, Modulation::Carrier);
    tracker_fixed.adaptive_ekf = false;
    tracker_fixed.reset(0.0, 20.0, 0.0);

    tracker_fixed.p = nalgebra::Matrix3::zeros();

    tracker_fixed.lock_metric = 0.1;
    tracker_fixed.predict();
    let q_phase_1 = tracker_fixed.p[(0, 0)];
    let q_freq_1 = tracker_fixed.p[(1, 1)];
    let q_chirp_1 = tracker_fixed.p[(2, 2)];

    tracker_fixed.p = nalgebra::Matrix3::zeros();
    tracker_fixed.lock_metric = 0.9;
    tracker_fixed.predict();
    let q_phase_2 = tracker_fixed.p[(0, 0)];
    let q_freq_2 = tracker_fixed.p[(1, 1)];
    let q_chirp_2 = tracker_fixed.p[(2, 2)];

    assert_eq!(q_phase_1, q_phase_2, "Process phase noise should be fixed");
    assert_eq!(
        q_freq_1, q_freq_2,
        "Process frequency noise should be fixed"
    );
    assert_eq!(q_chirp_1, q_chirp_2, "Process chirp noise should be fixed");

    let dt = 1.0 / fs;
    assert_eq!(q_phase_1, tracker_fixed.q_phase * dt);
    assert_eq!(q_freq_1, tracker_fixed.q_freq * dt);
    assert_eq!(q_chirp_1, tracker_fixed.q_chirp * dt);

    // Test Case 2: adaptive_ekf = true
    let mut tracker_adaptive = CarrierPllEkf::new(fs, Modulation::Carrier);
    tracker_adaptive.adaptive_ekf = true;
    tracker_adaptive.reset(0.0, 20.0, 0.0);

    tracker_adaptive.p = nalgebra::Matrix3::zeros();
    tracker_adaptive.lock_metric = 0.1;
    tracker_adaptive.predict();
    let q_phase_adaptive_1 = tracker_adaptive.p[(0, 0)];

    tracker_adaptive.p = nalgebra::Matrix3::zeros();
    tracker_adaptive.lock_metric = 0.9;
    tracker_adaptive.predict();
    let q_phase_adaptive_2 = tracker_adaptive.p[(0, 0)];

    assert_ne!(
        q_phase_adaptive_1, q_phase_adaptive_2,
        "Process noise should vary under adaptive EKF"
    );
    assert!(
        q_phase_adaptive_2 < q_phase_adaptive_1,
        "Noise should be narrower when locked/higher lock_metric"
    );
}

#[test]
fn test_no_dual_lock_relies_only_on_lock_metric_threshold() {
    let fs = 1000.0;

    // Test Case 1: dual_lock = false
    let mut tracker = CarrierPllEkf::new(fs, Modulation::Carrier);
    tracker.dual_lock = false;

    tracker.reset(0.0, 20.0, 0.0);
    tracker.convergence_guard = 0; // Bypass grace period for direct lock_metric testing
    tracker.lock_metric = 0.08;
    tracker.pr_sum_abs_i = 100.0;
    tracker.pr_sum_q_sq = 1e-5;

    tracker.update(Complex::new(1.0, 0.0));
    assert!(
        !tracker.is_locked,
        "Should unlock when lock_metric < 0.1 and dual_lock is disabled, even with high PR"
    );

    tracker.reset(0.0, 20.0, 0.0);
    tracker.convergence_guard = 0;
    tracker.lock_metric = 0.15;
    tracker.pr_sum_abs_i = 10.0;
    tracker.pr_sum_q_sq = 100.0;

    tracker.update(Complex::new(1.0, 0.0));
    assert!(
        tracker.is_locked,
        "Should remain locked when lock_metric >= 0.1 and dual_lock is disabled, even with low PR"
    );

    // Contrast: Test Case 2: dual_lock = true
    let mut tracker_dl = CarrierPllEkf::new(fs, Modulation::Carrier);
    tracker_dl.dual_lock = true;

    tracker_dl.reset(0.0, 20.0, 0.0);
    tracker_dl.convergence_guard = 0;
    tracker_dl.lock_metric = 0.08;
    tracker_dl.pr_sum_abs_i = 100.0;
    tracker_dl.pr_sum_q_sq = 1e-5;

    tracker_dl.update(Complex::new(1.0, 0.0));
    assert!(
        tracker_dl.is_locked,
        "Should remain locked with dual_lock enabled when lock_metric is in [0.02, 0.1) and PR is high"
    );

    tracker_dl.reset(0.0, 20.0, 0.0);
    tracker_dl.convergence_guard = 0;
    tracker_dl.lock_metric = 0.015;
    tracker_dl.pr_sum_abs_i = 100.0;
    tracker_dl.pr_sum_q_sq = 1e-5;

    tracker_dl.update(Complex::new(1.0, 0.0));
    assert!(
        !tracker_dl.is_locked,
        "Should unlock when lock_metric < 0.02 even with high PR"
    );
}

#[test]
fn test_dual_stage_lock_prevents_lock_flickering_during_brief_signal_fades() {
    let fs = 1000.0;
    let mut rng = SimpleRng::new(42);

    // Signal SNR = 15.0 dB
    let snr_db = 15.0;
    let snr_lin = 10.0f64.powf(snr_db / 10.0);
    let amp = snr_lin.sqrt();
    let target_freq = 20.0;

    let mut samples = Vec::new();
    // 2000 samples of noisy signal
    for n in 0..2000 {
        let t = n as f64 / fs;
        let phase = 2.0 * std::f64::consts::PI * target_freq * t;
        let sig = Complex::new((amp * phase.cos()) as f32, (amp * phase.sin()) as f32);
        let (n_re, n_im) = rng.next_gaussian();
        let noise = Complex::new(
            n_re * std::f32::consts::FRAC_1_SQRT_2,
            n_im * std::f32::consts::FRAC_1_SQRT_2,
        );
        samples.push(sig + noise);
    }

    // Fade: 900 samples of pure thermal noise (SNR = -inf)
    let mut fading_samples = samples.clone();
    for n in 1000..1900 {
        let (n_re, n_im) = rng.next_gaussian();
        fading_samples[n] = Complex::new(
            n_re * std::f32::consts::FRAC_1_SQRT_2,
            n_im * std::f32::consts::FRAC_1_SQRT_2,
        );
    }

    // Test Case 1: Nominal tracker (dual_lock = true)
    let mut tracker_nominal = CarrierPllEkf::new(fs, Modulation::Carrier);
    tracker_nominal.dual_lock = true;
    tracker_nominal.reset(0.0, 20.0, 0.0);

    // Lock onto first 1000 samples
    for n in 0..1000 {
        tracker_nominal.predict();
        tracker_nominal.update(samples[n]);
    }
    assert!(
        tracker_nominal.is_locked,
        "Nominal tracker should be locked initially"
    );

    // Manually force lock_metric to 0.22 right before the fade
    // to simulate a transient tracking jitter (e.g. from phase noise or dynamics)
    tracker_nominal.lock_metric = 0.22;

    // Run through fade
    let mut nominal_flickered = false;
    for n in 1000..1900 {
        tracker_nominal.predict();
        tracker_nominal.update(fading_samples[n]);
        if !tracker_nominal.is_locked {
            nominal_flickered = true;
        }
    }

    // Reset RNG seed to get identical noise
    let mut rng_no_dl = SimpleRng::new(42);
    let mut samples_no_dl = Vec::new();
    for n in 0..2000 {
        let t = n as f64 / fs;
        let phase = 2.0 * std::f64::consts::PI * target_freq * t;
        let sig = Complex::new((amp * phase.cos()) as f32, (amp * phase.sin()) as f32);
        let (n_re, n_im) = rng_no_dl.next_gaussian();
        let noise = Complex::new(
            n_re * std::f32::consts::FRAC_1_SQRT_2,
            n_im * std::f32::consts::FRAC_1_SQRT_2,
        );
        samples_no_dl.push(sig + noise);
    }
    let mut fading_samples_no_dl = samples_no_dl.clone();
    for n in 1000..1900 {
        let (n_re, n_im) = rng_no_dl.next_gaussian();
        fading_samples_no_dl[n] = Complex::new(
            n_re * std::f32::consts::FRAC_1_SQRT_2,
            n_im * std::f32::consts::FRAC_1_SQRT_2,
        );
    }

    // Test Case 2: No-Dual-Lock tracker (dual_lock = false)
    let mut tracker_no_dl = CarrierPllEkf::new(fs, Modulation::Carrier);
    tracker_no_dl.dual_lock = false;
    tracker_no_dl.reset(0.0, 20.0, 0.0);
    tracker_no_dl.convergence_guard = 0; // Bypass for direct testing

    for n in 0..1000 {
        tracker_no_dl.predict();
        tracker_no_dl.update(samples_no_dl[n]);
    }
    assert!(
        tracker_no_dl.is_locked,
        "No-DL tracker should be locked initially"
    );

    // Manually force lock_metric to 0.22 right before the fade
    tracker_no_dl.lock_metric = 0.22;

    // Run through fade
    let mut no_dl_flickered = false;
    for n in 1000..1900 {
        tracker_no_dl.predict();
        tracker_no_dl.update(fading_samples_no_dl[n]);
        if !tracker_no_dl.is_locked {
            no_dl_flickered = true;
        }
    }

    println!(
        "Nominal (dual-stage) lock flickered during fade: {}",
        nominal_flickered
    );
    println!(
        "Disabled (no-dual-lock) lock flickered during fade: {}",
        no_dl_flickered
    );
    println!("Nominal final lock_metric: {}", tracker_nominal.lock_metric);
    println!("No-DL final lock_metric: {}", tracker_no_dl.lock_metric);

    assert!(
        !nominal_flickered,
        "Nominal dual-stage lock should survive the brief fade"
    );
    assert!(
        no_dl_flickered,
        "Without dual-stage lock, the tracker should unlock/flicker during the brief fade"
    );
}

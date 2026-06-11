use num_complex::Complex;

#[path = "../src/daemon.rs"]
pub mod daemon;

#[path = "../src/ekf.rs"]
pub mod ekf;

#[path = "../src/orbit.rs"]
pub mod orbit;

#[path = "../src/tui.rs"]
pub mod tui;

#[path = "../src/orbit_solver.rs"]
pub mod orbit_solver;

#[path = "../src/dsp.rs"]
pub mod dsp;

#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct CalibrationData {
    pub df0: f64,
    pub timestamp: chrono::DateTime<chrono::Utc>,
}

#[derive(clap::Parser, Debug)]
pub struct Args {
    pub frequency: f64,
    pub sample_rate: f64,
    pub fft_size: usize,
    pub step_size: usize,
    pub lat: f64,
    pub lon: f64,
    pub alt: f64,
    pub tle: Option<String>,
    pub no_download_tle: bool,
    pub tle_url: Option<String>,
    pub sdr: Option<String>,
    pub gain: Option<f64>,
    pub lna_gain: f64,
    pub amp_gain: f64,
    pub vga_gain: f64,
    pub save_pass: Option<String>,
    pub solve_location: Option<String>,
    pub simulate: bool,
    pub sim_offset: f64,
    pub sim_start_time: bool,
    pub daemon: bool,
    pub blind: bool,
    pub output_dir: String,
    pub no_dc_skip: bool,
    pub min_snr: f32,
    pub decimate: Option<usize>,
    pub no_guided: bool,
    pub guided_window: Option<f64>,
    pub max_range: Option<f64>,
    pub max_guided_sats: usize,
    pub no_agc: bool,
    pub no_notch_spurs: bool,
    pub no_tui: bool,
    pub leodo: bool,
    pub leodo_log: String,
    pub solve_orbit: Option<String>,
    pub output_tle: Option<String>,
    pub modulation: dsp::Modulation,
    pub no_adaptive_ekf: bool,
    pub no_dual_lock: bool,
    pub no_gardner: bool,
    pub symbol_rate: Option<f64>,
    pub no_multihypothesis: bool,
    pub fade_timeout: f64,
}

pub fn get_process_rss_mb() -> f64 {
    0.0
}

use dsp::AbsolutePowerGainController;

struct TestRng(u32);
impl TestRng {
    fn new(seed: u32) -> Self {
        Self(seed)
    }

    fn next_complex(&mut self, std_dev: f32) -> Complex<f32> {
        // Simple LCG
        self.0 = self.0.wrapping_mul(1664525).wrapping_add(1013904223);
        let u1 = (self.0 as f32 / u32::MAX as f32).max(1e-6);
        self.0 = self.0.wrapping_mul(1664525).wrapping_add(1013904223);
        let u2 = self.0 as f32 / u32::MAX as f32;

        // Box-Muller transform
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f32::consts::PI as f32 * u2;
        let re = r * theta.cos() * std_dev;
        let im = r * theta.sin() * std_dev;
        Complex::new(re, im)
    }
}

#[test]
fn test_agc_zero_input_power() {
    let samples = vec![Complex::new(0.0, 0.0); 100];
    let mut lna = 0.0f32;
    let mut vga = 0.0f32;
    let mut amp = 0.0f32;

    // We start at 0 gain, and weak signal should increase gains.
    // LNA increases by 8.0 per step up to 40.0.
    for i in 1..=5 {
        AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(lna, 8.0 * i as f32);
        assert_eq!(vga, 0.0);
        assert_eq!(amp, 0.0);
    }

    // Once LNA reaches 40.0, VGA starts increasing by 2.0 per step up to 62.0.
    for i in 1..=31 {
        AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(lna, 40.0);
        assert_eq!(vga, 2.0 * i as f32);
        assert_eq!(amp, 0.0);
    }

    // Once LNA=40.0 and VGA=62.0, AMP should be set to 14.0.
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 14.0);

    // Subsequent updates should keep them clamped.
    for _ in 0..10 {
        AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(lna, 40.0);
        assert_eq!(vga, 62.0);
        assert_eq!(amp, 14.0);
    }
}

#[test]
fn test_agc_clipping_saturated_input() {
    // Saturated/clipping signal (amplitude 1.0, which causes rms > 0.5 and clipping ratio 1.0)
    let samples = vec![Complex::new(1.0, 0.0); 100];

    // Start at max gains
    let mut lna = 40.0f32;
    let mut vga = 62.0f32;
    let mut amp = 14.0f32;

    // Saturation lowers gains, prioritizing final stages first.
    // Step 1: AMP goes to 0.0.
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(amp, 0.0);
    assert_eq!(vga, 62.0);
    assert_eq!(lna, 40.0);

    // Steps 2..32: VGA decreases by 2.0 each step.
    for i in 1..=31 {
        AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(amp, 0.0);
        assert_eq!(vga, 62.0 - 2.0 * i as f32);
        assert_eq!(lna, 40.0);
    }

    // Steps 33..37: LNA decreases by 8.0 each step.
    for i in 1..=5 {
        AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(amp, 0.0);
        assert_eq!(vga, 0.0);
        assert_eq!(lna, 40.0 - 8.0 * i as f32);
    }

    // Subsequent updates should keep them clamped at 0.0.
    for _ in 0..10 {
        AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(lna, 0.0);
        assert_eq!(vga, 0.0);
        assert_eq!(amp, 0.0);
    }
}

#[test]
fn test_agc_rapid_step_changes() {
    let weak_samples = vec![Complex::new(1e-4, 1e-4); 100];
    let saturated_samples = vec![Complex::new(1.0, 0.0); 100];

    let mut lna = 0.0f32;
    let mut vga = 0.0f32;
    let mut amp = 0.0f32;

    // 1. Weak signal leads to max gains
    for _ in 0..50 {
        AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    }
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 14.0);

    // 2. Sudden saturated signal: should lower gains back to 0.0
    for _ in 0..50 {
        AbsolutePowerGainController::update_gain(&saturated_samples, &mut lna, &mut vga, &mut amp);
    }
    assert_eq!(lna, 0.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // 3. Switch back to weak: should rise back to max
    for _ in 0..50 {
        AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    }
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 14.0);
}

#[test]
fn test_agc_noisy_inputs() {
    let mut rng = TestRng::new(1337);

    // 1. Low noise power: gains should rise
    let mut lna = 0.0f32;
    let mut vga = 0.0f32;
    let mut amp = 0.0f32;
    let low_noise: Vec<Complex<f32>> = (0..100).map(|_| rng.next_complex(0.01)).collect();
    AbsolutePowerGainController::update_gain(&low_noise, &mut lna, &mut vga, &mut amp);
    assert!(lna > 0.0, "Expected LNA to increase for low noise");

    // 2. High noise power: gains should fall
    let mut lna = 40.0f32;
    let mut vga = 62.0f32;
    let mut amp = 14.0f32;
    let high_noise: Vec<Complex<f32>> = (0..100).map(|_| rng.next_complex(1.5)).collect();
    AbsolutePowerGainController::update_gain(&high_noise, &mut lna, &mut vga, &mut amp);
    assert!(
        amp == 0.0,
        "Expected AMP to drop under saturated/high noise power"
    );

    // 3. Moderate noise power: should stay in hysteresis zone (no changes)
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    // std_dev = 0.2 -> rms = ~0.2, clipping ratio = 0.0
    let mod_noise: Vec<Complex<f32>> = (0..500).map(|_| rng.next_complex(0.2)).collect();
    AbsolutePowerGainController::update_gain(&mod_noise, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 24.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);
}

#[test]
fn test_agc_absolute_power_estimation_invariance() {
    // True analog signal (a sine wave)
    let mut analog_signal = Vec::new();
    let num_samples = 1000;
    let analog_amplitude = 0.01f32; // power = 10 * log10(analog_amplitude^2) = -40.0 dB
    for n in 0..num_samples {
        let theta = 2.0 * std::f32::consts::PI * 0.05 * n as f32;
        analog_signal.push(Complex::new(theta.cos(), theta.sin()) * analog_amplitude);
    }

    // Verify true analog power
    let mut sum_mag_sq = 0.0f64;
    for s in &analog_signal {
        sum_mag_sq += (s.re * s.re + s.im * s.im) as f64;
    }
    let p_analog = 10.0 * (sum_mag_sq / num_samples as f64).log10() as f32;
    assert!((p_analog - -40.0).abs() < 1e-5);

    // Test different gain settings
    let test_gains = vec![
        (0.0f32, 0.0f32, 0.0f32),    // G = 0
        (8.0f32, 16.0f32, 0.0f32),   // G = 24
        (24.0f32, 32.0f32, 0.0f32),  // G = 56
        (40.0f32, 62.0f32, 14.0f32), // G = 116
        (16.0f32, 24.0f32, 0.0f32),  // G = 40
    ];

    for (lna, vga, amp) in test_gains {
        let nominal_gain = lna + vga + amp;
        let cal_correction = AbsolutePowerGainController::get_gain_calibration_offset(nominal_gain);

        // Physical gain is nominal gain minus the calibration correction
        let physical_gain_db = nominal_gain - cal_correction;

        // Scaling factor to apply to analog samples to get digital samples
        let scaling_factor = 10.0f32.powf(physical_gain_db / 20.0);

        let digital_samples: Vec<Complex<f32>> =
            analog_signal.iter().map(|s| s * scaling_factor).collect();

        // Estimate absolute power using the digital samples and gain settings
        let estimated_p_abs =
            AbsolutePowerGainController::estimate_absolute_power(&digital_samples, lna, vga, amp);

        // Check if estimated absolute power matches the true analog power
        let err = (estimated_p_abs - p_analog).abs();
        assert!(
            err < 1e-4,
            "Power estimation not invariant! G={}, estimated={}, true={}, error={}",
            nominal_gain,
            estimated_p_abs,
            p_analog,
            err
        );
    }
}

#[test]
fn test_agc_empty_samples() {
    let mut lna = 10.0f32;
    let mut vga = 20.0f32;
    let mut amp = 14.0f32;
    AbsolutePowerGainController::update_gain(&[], &mut lna, &mut vga, &mut amp);
    // Gains should remain unchanged when sample slice is empty
    assert_eq!(lna, 10.0);
    assert_eq!(vga, 20.0);
    assert_eq!(amp, 14.0);

    let power = AbsolutePowerGainController::estimate_absolute_power(&[], lna, vga, amp);
    assert_eq!(power, -120.0);
}

#[test]
fn test_agc_invalid_gains_clamping() {
    let samples = vec![Complex::new(0.1, 0.1); 10]; // RMS = sqrt(0.02) = ~0.1414, within hysteresis band (0.05..0.5)

    // LNA out of bounds
    let mut lna = -10.0f32;
    let mut vga = 20.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 0.0); // clamped to 0.0

    let mut lna = 100.0f32;
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0); // clamped to 40.0

    // VGA out of bounds
    let mut lna = 20.0f32;
    let mut vga = -50.0f32;
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(vga, 0.0);

    let mut vga = 200.0f32;
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(vga, 62.0);

    // AMP out of bounds / thresholding
    let mut amp = 5.0f32;
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(amp, 0.0); // < 7.0 clamped to 0.0

    let mut amp = 8.0f32;
    AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(amp, 14.0); // >= 7.0 clamped to 14.0
}

#[test]
fn test_agc_hysteresis_boundaries() {
    // 1. RMS saturation threshold: rms > 0.5
    // Let's create samples that result in RMS exactly 0.5, and slightly above 0.5.
    // RMS = sqrt(mean(mag^2)). For N constant samples of magnitude A, RMS is A.

    // RMS = 0.5: should NOT trigger saturation (it's <= 0.5)
    let samples_rms_0_5 = vec![Complex::new(0.5, 0.0); 100];
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples_rms_0_5, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 24.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);

    // RMS = 0.5001: should trigger saturation (it's > 0.5)
    let samples_rms_0_5001 = vec![Complex::new(0.5001, 0.0); 100];
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples_rms_0_5001, &mut lna, &mut vga, &mut amp);
    // Saturation lowers gains, prioritizing final stages. AMP is 0.0, so VGA goes down by 2.0.
    assert_eq!(vga, 30.0);
    assert_eq!(lna, 24.0);

    // 2. RMS weak-signal threshold: rms < 0.05
    // RMS = 0.05: should NOT trigger weak-signal (it's >= 0.05)
    let samples_rms_0_05 = vec![Complex::new(0.05, 0.0); 100];
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples_rms_0_05, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 24.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);

    // RMS = 0.0499: should trigger weak-signal (it's < 0.05)
    let samples_rms_0_0499 = vec![Complex::new(0.0499, 0.0); 100];
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples_rms_0_0499, &mut lna, &mut vga, &mut amp);
    // Weak signal raises gains, prioritizing early stages. LNA < 40.0, so LNA goes up by 8.0.
    assert_eq!(lna, 32.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);

    // 3. Clipping ratio threshold: clipping_ratio > 0.005
    // Let's create a slice of 1000 samples. 0.98 is the clipping threshold for a single sample.
    // 5 clipped samples out of 1000 => ratio = 0.005 (should NOT trigger saturation)
    let mut samples_clip_0_005 = vec![Complex::new(0.1, 0.0); 1000];
    for i in 0..5 {
        samples_clip_0_005[i] = Complex::new(0.98, 0.0);
    }
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples_clip_0_005, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 24.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);

    // 6 clipped samples out of 1000 => ratio = 0.006 (should trigger saturation)
    let mut samples_clip_0_006 = vec![Complex::new(0.1, 0.0); 1000];
    for i in 0..6 {
        samples_clip_0_006[i] = Complex::new(0.98, 0.0);
    }
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;
    AbsolutePowerGainController::update_gain(&samples_clip_0_006, &mut lna, &mut vga, &mut amp);
    assert_eq!(vga, 30.0);
    assert_eq!(lna, 24.0);
}

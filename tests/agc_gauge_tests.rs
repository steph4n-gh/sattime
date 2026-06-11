use num_complex::Complex;

#[path = "../src/dsp.rs"]
pub mod dsp;

#[path = "../src/orbit.rs"]
pub mod orbit;

#[path = "../src/daemon.rs"]
pub mod daemon;

#[path = "../src/tui.rs"]
pub mod tui;

#[path = "../src/ekf.rs"]
pub mod ekf;

#[path = "../src/orbit_solver.rs"]
pub mod orbit_solver;

// Root-level mocks/stubs for submodule crate references in dsp, daemon, ekf, etc.
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

// =========================================================================
// AGC LOOP SPECIFIC TESTS
// =========================================================================

#[test]
fn test_agc_initial_clamping() {
    let mut lna = -10.0;
    let mut vga = 100.0;
    let mut amp = 5.0; // < 7.0 should clamp to 0.0

    // Even with empty samples, the clamping happens first
    let empty_samples: [Complex<f32>; 0] = [];
    dsp::AbsolutePowerGainController::update_gain(&empty_samples, &mut lna, &mut vga, &mut amp);

    assert_eq!(lna, 0.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 0.0);

    // Test clamp bounds and upper threshold of amp gain
    let mut lna2 = 50.0;
    let mut vga2 = -1.0;
    let mut amp2 = 7.0; // >= 7.0 should clamp to 14.0
    dsp::AbsolutePowerGainController::update_gain(&empty_samples, &mut lna2, &mut vga2, &mut amp2);

    assert_eq!(lna2, 40.0);
    assert_eq!(vga2, 0.0);
    assert_eq!(amp2, 14.0);
}

#[test]
fn test_agc_weak_signal_prioritization_increasing() {
    // Start with all-zero gain
    let mut lna = 0.0;
    let mut vga = 0.0;
    let mut amp = 0.0;

    // A weak signal vector (RMS = 0.01 < 0.05)
    let weak_samples = vec![Complex::new(0.01, 0.0); 100];

    // Prioritization check: LNA must increase first
    // 1st step: LNA increases by 8.0, others remain 0.0
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 8.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // 2nd step: LNA becomes 16.0
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 16.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // 3rd to 5th steps: LNA reaches max 40.0
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp); // 24.0
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp); // 32.0
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp); // 40.0
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // Once LNA is maxed (40.0), VGA must start increasing
    // 6th step: VGA increases from 0.0 to 2.0
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 2.0);
    assert_eq!(amp, 0.0);

    // We fast-forward VGA to 60.0
    vga = 60.0;
    // 7th step: VGA becomes 62.0 (max)
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 0.0);

    // Once LNA and VGA are maxed, AMP must start increasing
    // 8th step: AMP becomes 14.0 (max)
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 14.0);

    // Further weak signals must leave gains unchanged at their maximums
    dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 14.0);
}

#[test]
fn test_agc_strong_signal_prioritization_decreasing() {
    // Start with all max gain
    let mut lna = 40.0;
    let mut vga = 62.0;
    let mut amp = 14.0;

    // A strong signal vector (RMS = 1.0 > 0.5, clipping = 1.0 > 0.005)
    let strong_samples = vec![Complex::new(1.0, 0.0); 100];

    // Prioritization check: AMP must decrease first
    // 1st step: AMP decreases from 14.0 to 0.0, others remain at max
    dsp::AbsolutePowerGainController::update_gain(&strong_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 62.0);
    assert_eq!(amp, 0.0);

    // Once AMP is zero, VGA must start decreasing
    // 2nd step: VGA decreases by 2.0 (from 62.0 to 60.0)
    dsp::AbsolutePowerGainController::update_gain(&strong_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 60.0);
    assert_eq!(amp, 0.0);

    // Fast-forward VGA to 2.0
    vga = 2.0;
    // 3rd step: VGA becomes 0.0
    dsp::AbsolutePowerGainController::update_gain(&strong_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 40.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // Once AMP and VGA are zero, LNA must start decreasing
    // 4th step: LNA decreases by 8.0 (from 40.0 to 32.0)
    dsp::AbsolutePowerGainController::update_gain(&strong_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 32.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // Fast-forward LNA to 8.0
    lna = 8.0;
    // 5th step: LNA becomes 0.0
    dsp::AbsolutePowerGainController::update_gain(&strong_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 0.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);

    // Further strong signals must leave gains unchanged at their minimums (0.0)
    dsp::AbsolutePowerGainController::update_gain(&strong_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 0.0);
    assert_eq!(vga, 0.0);
    assert_eq!(amp, 0.0);
}

#[test]
fn test_agc_hysteresis_stable_middle() {
    // Start with arbitrary intermediate gains
    let mut lna = 24.0;
    let mut vga = 32.0;
    let mut amp = 0.0;

    // Intermediate signal (RMS = 0.2, clipping ratio = 0.0)
    // 0.05 <= RMS <= 0.5, clipping ratio <= 0.005
    let mid_samples = vec![Complex::new(0.2, 0.0); 100];

    // Gains must remain completely unchanged over multiple updates
    for _ in 0..10 {
        dsp::AbsolutePowerGainController::update_gain(&mid_samples, &mut lna, &mut vga, &mut amp);
        assert_eq!(lna, 24.0);
        assert_eq!(vga, 32.0);
        assert_eq!(amp, 0.0);
    }
}

#[test]
fn test_agc_mixed_signal_impulse_saturation() {
    // Start with max gains
    let mut lna = 40.0;
    let mut vga = 62.0;
    let mut amp = 14.0;

    // Vector of 1000 samples, 10 samples are fully clipped, 990 samples are 0.0
    // Clipping ratio = 10 / 1000 = 0.01 (> 0.005)
    // RMS = sqrt(10 * 1.0^2 / 1000) = sqrt(0.01) = 0.1 (within intermediate range)
    let mut mixed_samples = vec![Complex::new(0.0, 0.0); 1000];
    for i in 0..10 {
        mixed_samples[i] = Complex::new(1.0, 0.0);
    }

    // Because clipping_ratio > 0.005, it must decrease gain (saturation state)
    dsp::AbsolutePowerGainController::update_gain(&mixed_samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(
        amp, 0.0,
        "AMP gain should have decreased first due to clipping"
    );
}

#[test]
fn test_agc_nan_inf_robustness() {
    // Start with intermediate gains
    let mut lna = 24.0;
    let mut vga = 32.0;
    let mut amp = 0.0;

    // Samples containing NaN
    let nan_samples = vec![Complex::new(f32::NAN, 0.0); 100];
    dsp::AbsolutePowerGainController::update_gain(&nan_samples, &mut lna, &mut vga, &mut amp);
    // Gains should not have changed because NaN comparisons return false
    assert_eq!(lna, 24.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);

    // Samples containing Infinity
    let inf_samples = vec![Complex::new(f32::INFINITY, 0.0); 100];
    dsp::AbsolutePowerGainController::update_gain(&inf_samples, &mut lna, &mut vga, &mut amp);
    // RMS will be Infinity, which is > 0.5, so it should trigger saturation and decrease gain
    assert_eq!(amp, 0.0);
    assert_eq!(vga, 30.0); // VGA decreased from 32 to 30
    assert_eq!(lna, 24.0);
}

#[test]
fn test_absolute_power_estimation() {
    let samples = vec![Complex::new(0.5, 0.5); 100];
    // mean_power = 0.5^2 + 0.5^2 = 0.5
    // p_dig = 10.0 * log10(0.5) = -3.0103
    // gain = 24.0 + 32.0 + 14.0 = 70.0
    // cal_correction = 0.005 * (70.0 - 24.0)^2 = 0.005 * 46.0^2 = 0.005 * 2116 = 10.58
    // expected = p_dig - gain + cal_correction = -3.0103 - 70.0 + 10.58 = -62.43
    let power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 24.0, 32.0, 14.0);
    assert!((power - (-62.43)).abs() < 0.1);

    // Empty samples estimation check
    let empty_samples: [Complex<f32>; 0] = [];
    let empty_power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&empty_samples, 24.0, 32.0, 14.0);
    assert_eq!(empty_power, -120.0);
}

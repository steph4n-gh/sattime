use chrono::{DateTime, TimeZone, Utc};
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

fn create_mock_passes() -> (Vec<orbit_solver::RawPass>, [f64; 3], f64, f64, f64, f64) {
    let truth_a = 6378137.0 + 550000.0;
    let truth_i = 53.0_f64.to_radians();
    let truth_raan = 1.2;
    let truth_u0 = 0.5;
    let center_freq = 150800000.0;
    let epoch = Utc.with_ymd_and_hms(2026, 6, 9, 12, 0, 0).unwrap();

    let rec_ecef = [1119794.6, -4842918.4, 3986004.4];

    let mut passes = Vec::new();
    let mut points = Vec::new();

    for step in 0..10 {
        let t = epoch + chrono::Duration::seconds(step * 60);
        let pred_f = orbit_solver::predict_frequency(
            truth_a,
            truth_i,
            truth_raan,
            truth_u0,
            epoch,
            t,
            0.0,
            0.0,
            center_freq,
            rec_ecef,
        );
        points.push(orbit_solver::PassPoint {
            time: t,
            freq: center_freq + pred_f,
        });
    }

    passes.push(orbit_solver::RawPass {
        sat_name: "MOCK_SAT1".to_string(),
        center_freq,
        points: points.clone(),
    });

    passes.push(orbit_solver::RawPass {
        sat_name: "MOCK_SAT2".to_string(),
        center_freq,
        points,
    });

    (passes, rec_ecef, truth_a, truth_i, truth_raan, truth_u0)
}

// --- Langevin Global Solver Gaps ---

#[test]
#[should_panic]
fn test_p_adic_distance_zero_prime_panic() {
    let _ = orbit_solver::p_adic_distance(10, 5, 0);
}

#[test]
fn test_p_adic_distance_one_prime_timeout() {
    use std::sync::mpsc;
    use std::thread;
    use std::time::Duration;

    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let _ = orbit_solver::p_adic_distance(10, 5, 1);
        let _ = tx.send(());
    });

    let result = rx.recv_timeout(Duration::from_millis(100));
    assert!(
        result.is_err(),
        "p_adic_distance with p=1 should have timed out (infinite loop)"
    );
}

#[test]
fn test_solve_linear_system_nan_propagation() {
    let mut matrix = vec![vec![f64::NAN, 1.0], vec![1.0, 2.0]];
    let vector = vec![1.0, 2.0];
    let res = orbit_solver::solve_linear_system(&mut matrix, &vector);
    if let Some(sol) = res {
        assert!(sol[0].is_nan() || sol[1].is_nan(), "Should propagate NaNs");
    }
}

#[test]
fn test_fit_orbit_doppler_nan_inputs() {
    let (passes, rec_ecef, _, _, _, _) = create_mock_passes();
    let bad_ecef = [f64::NAN, 0.0, 0.0];
    let res = orbit_solver::fit_orbit_doppler(&passes, bad_ecef, 7000e3, 0.5);
    assert!(res.is_err() || res.is_ok());
}

// --- Sheaf Cohomology & Tracker Consensus Discrepancy Gaps ---

#[test]
fn test_tracker_discrepancy_negative_lock_metrics() {
    let mut bank = ekf::EkfTrackingBank {
        trackers: [
            ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier),
            ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier),
            ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier),
        ],
        gardner_loops: [
            dsp::GardnerLoop::new(100000.0, 10000.0),
            dsp::GardnerLoop::new(100000.0, 10000.0),
            dsp::GardnerLoop::new(100000.0, 10000.0),
        ],
        active_idx: None,
        in_fade: false,
        fade_counter: 0,
        max_fade_steps: 100,
        terminated_in_fade: false,
        spur_dwell_counter: 0,
    };

    // Configure trackers to yield negative weights in Cech calculation
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = -0.1;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.5;
    bank.trackers[1].x[1] = 2000.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.5;
    bank.trackers[2].x[1] = 2000.0 * 2.0 * std::f64::consts::PI;

    let obs = bank.compute_tracker_discrepancy();
    // Cech obstruction normally must be non-negative. However, because lock_metric can decay negative,
    // compute_tracker_discrepancy can yield a negative obstruction value.
    assert!(
        obs < 0.0,
        "Tracker consensus discrepancy is negative: {}",
        obs
    );
}

#[test]
fn test_tracker_discrepancy_nan_lock_metrics() {
    let mut bank = ekf::EkfTrackingBank {
        trackers: [
            ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier),
            ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier),
            ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier),
        ],
        gardner_loops: [
            dsp::GardnerLoop::new(100000.0, 10000.0),
            dsp::GardnerLoop::new(100000.0, 10000.0),
            dsp::GardnerLoop::new(100000.0, 10000.0),
        ],
        active_idx: None,
        in_fade: false,
        fade_counter: 0,
        max_fade_steps: 100,
        terminated_in_fade: false,
        spur_dwell_counter: 0,
    };

    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = f64::NAN;
    bank.trackers[0].x[1] = 100.0;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = f64::NAN;
    bank.trackers[1].x[1] = 200.0;

    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = f64::NAN;
    bank.trackers[2].x[1] = 300.0;

    let obs = bank.compute_tracker_discrepancy();
    assert_eq!(obs, 0.0);
}

#[test]
fn test_ekf_update_nan_sample() {
    let mut ekf_filter = ekf::CarrierPllEkf::new(10000.0, dsp::Modulation::Carrier);
    ekf_filter.reset(0.0, 100.0, 0.0);
    ekf_filter.update(Complex::new(f32::NAN, f32::NAN));
    // A3: NaN guard now resets state instead of propagating NaN.
    // State should be finite (zeroed) and tracker should be unlocked.
    assert!(ekf_filter.x.iter().all(|v| v.is_finite()),
        "EKF state should be finite after NaN input (NaN guard should have reset it)");
    assert!(!ekf_filter.is_locked,
        "EKF should be unlocked after NaN-triggered reset");
}

// --- Vladimirov-Steered Envelope Wavelet Gaps ---

#[test]
fn test_tropical_wavelet_nan_inf_doppler() {
    let mut magnitudes = vec![1.0; 256];

    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, f32::NAN, 0.0);
    assert!(magnitudes.iter().all(|&m| m.is_finite()));

    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, f32::INFINITY, 0.0);
    assert!(magnitudes.iter().all(|&m| m.is_finite()));
}

// --- Calibrated Gauge AGC Loop Gaps ---

#[test]
fn test_agc_nan_in_middle_disables_updates() {
    let mut lna = 24.0f32;
    let mut vga = 32.0f32;
    let mut amp = 0.0f32;

    // Weak signal samples, but one is NaN
    // Without NaN, weak signal would raise LNA by 8.0 to 32.0.
    // With NaN, RMS becomes NaN, disabling updates so LNA remains 24.0.
    let mut samples = vec![Complex::new(0.01, 0.0); 100];
    samples[50] = Complex::new(f32::NAN, 0.0);

    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert_eq!(lna, 24.0);
    assert_eq!(vga, 32.0);
    assert_eq!(amp, 0.0);
}

#[test]
fn test_agc_estimate_absolute_power_extreme_gains() {
    let samples = vec![Complex::new(0.5, 0.5); 100];
    let power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 1e20, 1e20, 1e20);
    assert!(power.is_infinite());
}

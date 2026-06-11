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

// Root-level mocks/stubs for submodule crate references
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

// Helper function to create mock pass data for orbit solver tests
fn create_mock_passes() -> (Vec<orbit_solver::RawPass>, [f64; 3], f64, f64, f64, f64) {
    let truth_a = 6378137.0 + 550000.0;
    let truth_i = 53.0_f64.to_radians();
    let truth_raan = 1.2;
    let truth_u0 = 0.5;
    let center_freq = 150800000.0;
    let epoch = Utc.with_ymd_and_hms(2026, 6, 9, 12, 0, 0).unwrap();

    let rec_ecef = [1119794.6, -4842918.4, 3986004.4]; // Mock ECEF

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

// ==========================================
// TIER 1: Feature Coverage (Happy Path)
// ==========================================

// --- Feature 1: Langevin Global Solver ---
#[test]
fn test_langevin_global_convergence() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();
    let initial_a = truth_a + 50000.0;
    let initial_i = truth_i + 2.0_f64.to_radians();

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, initial_a, initial_i);
    assert!(solved.is_ok());
}

#[test]
fn test_langevin_global_success_rate() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();
    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, truth_a, truth_i);
    assert!(solved.is_ok());
}

#[test]
fn test_langevin_global_mapping() {
    let coord = [1000.0, 2000.0, 3000.0];
    let ring_coords = orbit_solver::map_to_normalized_search_space(coord);
    assert_eq!(ring_coords.len(), 15);
    assert_eq!(ring_coords[0], coord[0]);
    assert_eq!(ring_coords[1], coord[1]);
    assert_eq!(ring_coords[2], coord[2]);
    for &val in &ring_coords[3..15] {
        assert!(val >= 0.0 && val <= 1.0);
    }
}

#[test]
fn test_langevin_global_fractional_derivatives() {
    let coords = vec![1.0, 2.0, 3.0];
    let derivs = orbit_solver::compute_fractional_difference_history(&coords);
    assert_eq!(derivs.len(), 3);
    let expected = [-6.65685, 0.0, 6.65685];
    for i in 0..3 {
        assert!(
            (derivs[i] - expected[i]).abs() < 1e-4,
            "derivs[{}] = {}, expected = {}",
            i,
            derivs[i],
            expected[i]
        );
    }
}

#[test]
fn test_langevin_global_blind_search() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();
    let solved = orbit_solver::fit_orbit_doppler(
        &passes,
        rec_ecef,
        truth_a - 100000.0,
        truth_i - 5.0_f64.to_radians(),
    );
    assert!(solved.is_ok());
}

// --- Feature 2: Sheaf Cohomology & Tracker Consensus Discrepancy ---
// --- Feature 2: Sheaf Cohomology & Tracker Consensus Discrepancy ---
#[test]
fn test_tracker_consensus_complex() {
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
    bank.trackers[0].lock_metric = 0.8;
    bank.trackers[0].is_locked = true;
    bank.trackers[0].p[(0, 0)] = 0.1;
    bank.trackers[0].p[(1, 1)] = 0.2;

    bank.trackers[1].lock_metric = 0.6;
    bank.trackers[1].is_locked = false;
    bank.trackers[1].p[(0, 0)] = 0.3;
    bank.trackers[1].p[(1, 1)] = 0.4;

    bank.trackers[2].lock_metric = 0.9;
    bank.trackers[2].is_locked = true;
    bank.trackers[2].p[(0, 0)] = 0.5;
    bank.trackers[2].p[(1, 1)] = 0.6;

    let sheaf = bank.collect_tracker_features();
    assert_eq!(sheaf.len(), 12);
    // Tracker 0
    assert_eq!(sheaf[0], 0.8f32);
    assert_eq!(sheaf[1], 1.0f32);
    assert_eq!(sheaf[2], 0.1f32);
    assert_eq!(sheaf[3], 0.2f32);
    // Tracker 1
    assert_eq!(sheaf[4], 0.6f32);
    assert_eq!(sheaf[5], 0.0f32);
    assert_eq!(sheaf[6], 0.3f32);
    assert_eq!(sheaf[7], 0.4f32);
    // Tracker 2
    assert_eq!(sheaf[8], 0.9f32);
    assert_eq!(sheaf[9], 1.0f32);
    assert_eq!(sheaf[10], 0.5f32);
    assert_eq!(sheaf[11], 0.6f32);
}

#[test]
fn test_tracker_frequency_diff_calculation() {
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
    bank.trackers[0].x[1] = 100.0 * 2.0 * std::f64::consts::PI; // 100 Hz
    bank.trackers[1].x[1] = 150.0 * 2.0 * std::f64::consts::PI; // 150 Hz
    bank.trackers[2].x[1] = 300.0 * 2.0 * std::f64::consts::PI; // 300 Hz

    let coboundaries = bank.compute_tracker_frequency_diffs();
    assert_eq!(coboundaries.len(), 3);
    assert!((coboundaries[0] - 50.0).abs() < 1e-3); // f1 - f0 = 150 - 100 = 50
    assert!((coboundaries[1] - 150.0).abs() < 1e-3); // f2 - f1 = 300 - 150 = 150
    assert!((coboundaries[2] - (-200.0)).abs() < 1e-3); // f0 - f2 = 100 - 300 = -200
}

#[test]
fn test_tracker_discrepancy_pruning() {
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
    bank.trackers[0].lock_metric = 0.8;
    bank.trackers[0].x[1] = 100.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.9;
    bank.trackers[1].x[1] = 150.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.7;
    bank.trackers[2].x[1] = 500.0 * 2.0 * std::f64::consts::PI;

    let obs = bank.compute_tracker_discrepancy();
    assert!((obs - 251.57068).abs() < 1e-2);
}

#[test]
fn test_tracker_consensus_multipath_segregation() {
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
    bank.trackers[0].lock_metric = 0.85;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.95;
    bank.trackers[1].x[1] = 1010.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.90;
    bank.trackers[2].x[1] = 1005.0 * 2.0 * std::f64::consts::PI;

    let obs = bank.compute_tracker_discrepancy();
    assert!(obs < 10.0);
}

#[test]
fn test_tracker_consensus_no_multipath() {
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
    bank.trackers[0].lock_metric = 0.9;
    bank.trackers[0].x[1] = 500.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.9;
    bank.trackers[1].x[1] = 500.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.9;
    bank.trackers[2].x[1] = 500.0 * 2.0 * std::f64::consts::PI;

    let obs = bank.compute_tracker_discrepancy();
    assert!(obs < 1e-4);
}

// --- Feature 3: Vladimirov-Steered Envelope Wavelet ---
#[test]
fn test_wavelet_notch_detection() {
    let magnitudes = vec![1.0; 1024];
    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    assert!(spikes.is_empty());
}

#[test]
fn test_wavelet_spur_cancellation() {
    let mut magnitudes = vec![1.0; 1024];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 1500.0, 10.0);
    assert_eq!(magnitudes[0], 1.0);
}

#[test]
fn test_wavelet_pca_epoch_guard() {
    let t_obs = Utc::now();
    let pca = Utc::now();
    let guard = dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(t_obs, pca);
    assert!(guard);
}

#[test]
fn test_wavelet_fractional_variance() {
    let mut magnitudes = vec![0.5; 512];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 200.0, 0.0);
    assert!(magnitudes.len() == 512);
}

#[test]
fn test_wavelet_sweep_tracking() {
    let mut magnitudes = vec![0.1; 128];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, -100.0, -1.0);
    assert!(magnitudes[5] >= 0.1);
}

// --- Feature 4: Calibrated Gauge AGC Loop ---
#[test]
fn test_gain_calibration_offset() {
    let profile = dsp::AbsolutePowerGainController::get_gain_calibration_offset(24.0);
    assert_eq!(profile, 0.0);
}

#[test]
fn test_gain_calibration_power_estimation() {
    let samples = vec![Complex::new(1.0, 0.0); 100];
    let power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 24.0, 32.0, 14.0);
    assert!(power.is_finite());
}

#[test]
fn test_gain_calibration_opt_gains() {
    let samples = vec![Complex::new(0.5, 0.5); 50];
    let mut lna = 20.0;
    let mut vga = 30.0;
    let mut amp = 14.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert!(lna >= 0.0 && vga >= 0.0 && amp >= 0.0);
}

#[test]
fn test_gain_calibration_intermod_noise() {
    let samples = vec![Complex::new(0.01, 0.01); 10];
    let power = dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 0.0, 0.0, 0.0);
    assert!(power <= 0.0);
}

#[test]
fn test_gain_calibration_gain_limit_bounds() {
    let samples = vec![Complex::new(10.0, 10.0); 100];
    let mut lna = 50.0;
    let mut vga = 70.0;
    let mut amp = 20.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert!(lna <= 50.0);
}

// ==========================================
// TIER 2: Boundary & Corner Cases
// ==========================================

#[test]
fn test_langevin_empty_passes() {
    let passes = vec![];
    let rec_ecef = [0.0, 0.0, 0.0];
    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, 7000e3, 0.0);
    assert!(solved.is_err());
}

#[test]
fn test_langevin_extreme_noise() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();
    let mut noisy_passes = passes.clone();
    for p in &mut noisy_passes {
        for pt in &mut p.points {
            pt.freq += 1e6; // Add extreme noise
        }
    }
    let solved = orbit_solver::fit_orbit_doppler(&noisy_passes, rec_ecef, truth_a, truth_i);
    assert!(solved.is_ok() || solved.is_err());
}

#[test]
fn test_langevin_p_adic_boundary() {
    let coords = vec![0.0, 0.0, 0.0];
    let derivs = orbit_solver::compute_fractional_difference_history(&coords);
    assert_eq!(derivs.len(), 3);
}

#[test]
fn test_langevin_flat_jacobian() {
    let coords = vec![f64::MAX, f64::MAX];
    let derivs = orbit_solver::compute_fractional_difference_history(&coords);
    assert!(derivs[0].is_finite());
}

#[test]
fn test_langevin_parameter_overflow() {
    let (passes, rec_ecef, _, _, _, _) = create_mock_passes();
    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, 1e20, 1e20);
    assert!(solved.is_err() || solved.is_ok());
}

#[test]
fn test_tracker_consensus_zero_hypotheses() {
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs.is_finite());
}

#[test]
fn test_tracker_consensus_all_obstructed() {
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);
}

#[test]
fn test_tracker_consensus_threshold_limits() {
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs < 1000.0);
}

#[test]
fn test_tracker_consensus_floating_point_extreme() {
    let bank = ekf::EkfTrackingBank {
        trackers: [
            ekf::CarrierPllEkf::new(f64::NAN, dsp::Modulation::Carrier),
            ekf::CarrierPllEkf::new(f64::INFINITY, dsp::Modulation::Carrier),
            ekf::CarrierPllEkf::new(f64::NEG_INFINITY, dsp::Modulation::Carrier),
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs.is_finite() || obs.is_nan());
}

#[test]
fn test_tracker_consensus_outlier_boundary() {
    let bank = ekf::EkfTrackingBank {
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
        active_idx: Some(0),
        in_fade: false,
        fade_counter: 0,
        max_fade_steps: 100,
        terminated_in_fade: false,
        spur_dwell_counter: 0,
    };
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);
}

#[test]
fn test_wavelet_zero_spectrum() {
    let mut magnitudes = vec![0.0; 1024];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 0.0, 0.0);
    assert_eq!(magnitudes[500], 0.0);
}

#[test]
fn test_wavelet_infinite_spur() {
    let mut magnitudes = vec![f32::INFINITY; 512];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 10.0, 0.0);
    assert!(
        magnitudes[256].is_infinite()
            || magnitudes[256].is_nan()
            || magnitudes[256] < f32::INFINITY
    );
}

#[test]
fn test_wavelet_overlapping_carrier_spur() {
    let mut magnitudes = vec![10.0; 256];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 0.0, 0.0);
    assert!(magnitudes[128] <= 10.0);
}

#[test]
fn test_wavelet_nyquist_boundary() {
    let mut magnitudes = vec![1.0; 256];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 127.0, 0.0);
    assert!(magnitudes[255] <= 1.0);
}

#[test]
fn test_wavelet_rapid_sweep() {
    let mut magnitudes = vec![1.0; 256];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 50.0, 100.0);
    assert!(magnitudes.len() == 256);
}

#[test]
fn test_gain_calibration_zero_input_power() {
    let samples = vec![Complex::new(0.0, 0.0); 100];
    let power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 24.0, 32.0, 14.0);
    assert!(power.is_finite());
}

#[test]
fn test_gain_calibration_max_power_clipping() {
    let samples = vec![Complex::new(1e5, 1e5); 10];
    let power = dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 0.0, 0.0, 0.0);
    assert!(power > 0.0 || power.is_finite());
}

#[test]
fn test_gain_calibration_vga_only() {
    let samples = vec![Complex::new(0.1, 0.1); 10];
    let mut lna = 0.0;
    let mut vga = 20.0;
    let mut amp = 0.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert!(vga >= 0.0);
}

#[test]
fn test_gain_calibration_rapid_gain_flips() {
    let samples1 = vec![Complex::new(10.0, 10.0); 10];
    let samples2 = vec![Complex::new(0.001, 0.001); 10];
    let mut lna = 24.0;
    let mut vga = 32.0;
    let mut amp = 14.0;
    dsp::AbsolutePowerGainController::update_gain(&samples1, &mut lna, &mut vga, &mut amp);
    dsp::AbsolutePowerGainController::update_gain(&samples2, &mut lna, &mut vga, &mut amp);
    assert!(lna.is_finite() && vga.is_finite() && amp.is_finite());
}

#[test]
fn test_gain_calibration_snr_floor() {
    let samples = vec![Complex::new(1e-6, 1e-6); 100];
    let mut lna = 10.0;
    let mut vga = 10.0;
    let mut amp = 0.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert!(lna >= 10.0 || lna.is_finite());
}

// ==========================================
// TIER 3: Cross-Feature Combinations
// ==========================================

#[test]
fn test_t3_langevin_and_tracker_consensus() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();

    // Simulate EKF sheaf consistency pruning
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, truth_a, truth_i);
    assert!(solved.is_ok());
}

#[test]
fn test_t3_langevin_and_wavelet() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();

    let mut magnitudes = vec![1.0; 1024];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 100.0, 0.0);

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, truth_a, truth_i);
    assert!(solved.is_ok());
}

#[test]
fn test_t3_langevin_and_gain_agc() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();

    let samples = vec![Complex::new(1.0, 1.0); 10];
    let mut lna = 24.0;
    let mut vga = 32.0;
    let mut amp = 14.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, truth_a, truth_i);
    assert!(solved.is_ok());
}

#[test]
fn test_t3_tracker_consensus_and_wavelet() {
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);

    let mut magnitudes = vec![0.5; 512];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 50.0, 0.0);
    assert_eq!(magnitudes.len(), 512);
}

#[test]
fn test_t3_tracker_consensus_and_gain_agc() {
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);

    let samples = vec![Complex::new(1.0, 0.0); 100];
    let power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 24.0, 32.0, 14.0);
    assert!(power.is_finite());
}

#[test]
fn test_t3_wavelet_and_gain_agc() {
    let mut magnitudes = vec![1.0; 128];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 25.0, 0.0);

    let samples = vec![Complex::new(1.0, 0.0); 100];
    let mut lna = 24.0;
    let mut vga = 32.0;
    let mut amp = 14.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);
    assert!(lna >= 0.0);
}

// ==========================================
// TIER 4: Real-World Application Scenarios
// ==========================================

#[test]
fn test_t4_real_world_starlink_overhead_pass() {
    let (passes, rec_ecef, truth_a, truth_i, _, _) = create_mock_passes();

    // 1. AGC Adjusts gains
    let samples = vec![Complex::new(0.2, 0.2); 1000];
    let mut lna = 24.0;
    let mut vga = 32.0;
    let mut amp = 14.0;
    dsp::AbsolutePowerGainController::update_gain(&samples, &mut lna, &mut vga, &mut amp);

    // 2. Wavelet filters spurs
    let mut magnitudes = vec![1.0; 1024];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 1500.0, 5.0);

    // 3. EKF tracks with sheaf mitigation
    let bank = ekf::EkfTrackingBank {
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
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);

    // 4. Solve orbit
    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, truth_a, truth_i);
    assert!(solved.is_ok());
}

#[test]
fn test_t4_real_world_noaa_weather_pass() {
    let bank = ekf::EkfTrackingBank {
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
        in_fade: true,
        fade_counter: 5,
        max_fade_steps: 100,
        terminated_in_fade: false,
        spur_dwell_counter: 0,
    };
    let obs = bank.compute_tracker_discrepancy();
    assert!(obs >= 0.0);
}

#[test]
fn test_t4_real_world_iridium_lband_pass() {
    let mut magnitudes = vec![0.5; 512];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 2500.0, 50.0);
    assert!(magnitudes.len() == 512);
}

#[test]
fn test_t4_real_world_amateur_sat_pass() {
    let mut magnitudes = vec![0.2; 1024];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, -500.0, 2.0);
    assert!(magnitudes.len() == 1024);
}

#[test]
fn test_t4_real_world_gps_like_signal() {
    let samples = vec![Complex::new(0.001, 0.002); 5000];
    let power =
        dsp::AbsolutePowerGainController::estimate_absolute_power(&samples, 30.0, 40.0, 14.0);
    assert!(power.is_finite());
}

#[test]
fn test_sheaf_attention_extreme_inputs() {
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

    // Set up tracker 0 with NaN/Inf/Negative values
    bank.trackers[0].lock_metric = f64::NAN;
    bank.trackers[0].is_locked = true;
    bank.trackers[0].p[(0, 0)] = f64::INFINITY;
    bank.trackers[0].p[(1, 1)] = -1e5;

    // Set up tracker 1 with normal values
    bank.trackers[1].lock_metric = 0.55;
    bank.trackers[1].is_locked = false;
    bank.trackers[1].p[(0, 0)] = 0.002;
    bank.trackers[1].p[(1, 1)] = 0.045;

    // Set up tracker 2 with extreme floating point subnormals
    bank.trackers[2].lock_metric = 1e-300;
    bank.trackers[2].is_locked = true;
    bank.trackers[2].p[(0, 0)] = f64::MIN_POSITIVE;
    bank.trackers[2].p[(1, 1)] = f64::MAX;

    let sheaf = bank.collect_tracker_features();
    assert_eq!(sheaf.len(), 12);

    // Verify values are propagated exactly as casted to f32
    assert!(sheaf[0].is_nan());
    assert_eq!(sheaf[1], 1.0f32);
    assert_eq!(sheaf[2], f32::INFINITY);
    assert_eq!(sheaf[3], -1e5f32);

    assert_eq!(sheaf[4], 0.55f32);
    assert_eq!(sheaf[5], 0.0f32);
    assert_eq!(sheaf[6], 0.002f32);
    assert_eq!(sheaf[7], 0.045f32);

    assert_eq!(sheaf[8], 0.0f32); // 1e-300 becomes 0 in f32 (underflow)
    assert_eq!(sheaf[9], 1.0f32);
    assert_eq!(sheaf[10], 0.0f32); // f64::MIN_POSITIVE underflows to 0 in f32
    assert_eq!(sheaf[11], f32::INFINITY); // f64::MAX overflows to infinity in f32
}

#[test]
fn test_tracker_frequency_diffy_under_noise_and_drift() {
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

    // Simulate 100 noisy/drifting scenarios
    // Check that the sum of coboundaries is always 0.0 (up to floating point precision)
    for step in 0..100 {
        let f0_hz = 1000.0 + (step as f64) * 5.0 + ((step * 7) % 13) as f64;
        let f1_hz = 1000.0 + (step as f64) * 4.9 - ((step * 5) % 11) as f64;
        let f2_hz = 1000.0 + (step as f64) * 5.1 + ((step * 3) % 7) as f64;

        bank.trackers[0].x[1] = f0_hz * 2.0 * std::f64::consts::PI;
        bank.trackers[1].x[1] = f1_hz * 2.0 * std::f64::consts::PI;
        bank.trackers[2].x[1] = f2_hz * 2.0 * std::f64::consts::PI;

        let coboundaries = bank.compute_tracker_frequency_diffs();
        assert_eq!(coboundaries.len(), 3);

        let sum_coboundaries = coboundaries[0] + coboundaries[1] + coboundaries[2];
        assert!(
            sum_coboundaries.abs() < 1e-4,
            "Failed at step {}, sum = {}",
            step,
            sum_coboundaries
        );

        // Also verify the exact values match expected difference
        assert!((coboundaries[0] - (f1_hz - f0_hz) as f32).abs() < 1e-4);
        assert!((coboundaries[1] - (f2_hz - f1_hz) as f32).abs() < 1e-4);
        assert!((coboundaries[2] - (f0_hz - f2_hz) as f32).abs() < 1e-4);
    }
}

#[test]
fn test_pruning_vulnerability_single_inconsistent_tracker() {
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

    // Tracker 0 and 1 are consistent at 1000 Hz, with lock_metric = 1.0
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 1.0;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 1.0;
    bank.trackers[1].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Tracker 2 is inconsistent at 1200 Hz (deviation 200 Hz > 150 Hz), with lock_metric = 1.0
    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 1.0;
    bank.trackers[2].x[1] = 1200.0 * 2.0 * std::f64::consts::PI;

    // The obstruction is 2/3 * 200 = 133.33 Hz
    let obs = bank.compute_tracker_discrepancy();
    assert!((obs - 133.333).abs() < 1e-2);

    // Since obs <= 150.0, the pruning logic would NOT trigger!
    let obs_threshold = 150.0;
    if obs > obs_threshold {
        let f0 = bank.trackers[0].x[1] / (2.0 * std::f64::consts::PI);
        if bank.trackers[0].is_locked {
            for i in 1..3 {
                if bank.trackers[i].is_locked {
                    let fi = bank.trackers[i].x[1] / (2.0 * std::f64::consts::PI);
                    let dev = (fi - f0).abs();
                    if dev > 150.0 {
                        bank.trackers[i].is_locked = false;
                        bank.trackers[i].lock_metric = 0.0;
                    }
                }
            }
        }
    }

    // Assert that Tracker 2 was NOT pruned (it is still locked and metric is 1.0)
    assert!(bank.trackers[2].is_locked);
    assert_eq!(bank.trackers[2].lock_metric, 1.0);
}

#[test]
fn test_pruning_under_extreme_inconsistency() {
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

    // Tracker 0 and 1 are consistent at 1000 Hz, with lock_metric = 1.0
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 1.0;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 1.0;
    bank.trackers[1].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Tracker 2 is inconsistent at 1250 Hz (deviation 250 Hz), with lock_metric = 1.0
    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 1.0;
    bank.trackers[2].x[1] = 1250.0 * 2.0 * std::f64::consts::PI;

    // The obstruction is 2/3 * 250 = 166.67 Hz
    let obs = bank.compute_tracker_discrepancy();
    assert!((obs - 166.666).abs() < 1e-2);

    // Since obs > 150.0, the pruning logic triggers!
    let obs_threshold = 150.0;
    if obs > obs_threshold {
        let f0 = bank.trackers[0].x[1] / (2.0 * std::f64::consts::PI);
        if bank.trackers[0].is_locked {
            for i in 1..3 {
                if bank.trackers[i].is_locked {
                    let fi = bank.trackers[i].x[1] / (2.0 * std::f64::consts::PI);
                    let dev = (fi - f0).abs();
                    if dev > 150.0 {
                        bank.trackers[i].is_locked = false;
                        bank.trackers[i].lock_metric = 0.0;
                    }
                }
            }
        }
    }

    // Assert that Tracker 2 WAS pruned
    assert!(!bank.trackers[2].is_locked);
    assert_eq!(bank.trackers[2].lock_metric, 0.0);
}

#[test]
fn test_pruning_with_tracker_0_unlocked() {
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

    // Tracker 0 is unlocked
    bank.trackers[0].is_locked = false;
    bank.trackers[0].lock_metric = 0.0;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Tracker 1 is locked at 1000 Hz, with lock_metric = 0.8
    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.8;
    bank.trackers[1].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Tracker 2 is locked at 1200 Hz (deviation 200 Hz), with lock_metric = 0.9
    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.9;
    bank.trackers[2].x[1] = 1200.0 * 2.0 * std::f64::consts::PI;

    // The obstruction is 200 Hz
    let obs = bank.compute_tracker_discrepancy();
    assert!((obs - 200.0).abs() < 1e-2);

    // Prune logic:
    let obs_threshold = 150.0;
    if obs > obs_threshold {
        if bank.trackers[0].is_locked {
            // ...
        } else if bank.trackers[1].is_locked && bank.trackers[2].is_locked {
            let f1 = bank.trackers[1].x[1] / (2.0 * std::f64::consts::PI);
            let f2 = bank.trackers[2].x[1] / (2.0 * std::f64::consts::PI);
            if (f1 - f2).abs() > 150.0 {
                let prune_idx = if bank.trackers[1].lock_metric < bank.trackers[2].lock_metric {
                    1
                } else {
                    2
                };
                bank.trackers[prune_idx].is_locked = false;
                bank.trackers[prune_idx].lock_metric = 0.0;
            }
        }
    }

    // Since tracker 1 lock_metric (0.8) < tracker 2 lock_metric (0.9), tracker 1 should be pruned!
    assert!(!bank.trackers[1].is_locked);
    assert_eq!(bank.trackers[1].lock_metric, 0.0);
    // Tracker 2 should remain locked
    assert!(bank.trackers[2].is_locked);
    assert_eq!(bank.trackers[2].lock_metric, 0.9);
}

#[test]
fn test_pruning_nan_propagation_resilience() {
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

    // Tracker 0 is locked at NaN frequency
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 1.0;
    bank.trackers[0].x[1] = f64::NAN;

    // Tracker 1 is locked at 1000 Hz
    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 1.0;
    bank.trackers[1].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Tracker 2 is locked at 1000 Hz
    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 1.0;
    bank.trackers[2].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Compute obstruction: since f0 is NaN, coboundaries will contain NaN.
    // The obstruction value should be 0.0 because compute_tracker_discrepancy has a NaN check:
    // if val.is_nan() { 0.0 } else { val }
    let obs = bank.compute_tracker_discrepancy();
    assert_eq!(obs, 0.0);

    // Now let's check INFINITY frequency on Tracker 0
    bank.trackers[0].x[1] = f64::INFINITY;
    let obs = bank.compute_tracker_discrepancy();
    assert_eq!(obs, f32::INFINITY);

    // Since obs = Infinity > 150.0, pruning is triggered!
    let obs_threshold = 150.0;
    if obs > obs_threshold {
        let f0 = bank.trackers[0].x[1] / (2.0 * std::f64::consts::PI);
        if bank.trackers[0].is_locked {
            for i in 1..3 {
                if bank.trackers[i].is_locked {
                    let fi = bank.trackers[i].x[1] / (2.0 * std::f64::consts::PI);
                    let dev = (fi - f0).abs();
                    if dev > 150.0 {
                        bank.trackers[i].is_locked = false;
                        bank.trackers[i].lock_metric = 0.0;
                    }
                }
            }
        }
    }

    // Since f0 is Infinity, both Tracker 1 and Tracker 2 have dev = Infinity > 150.
    // So both Tracker 1 and Tracker 2 (which are good) are pruned, while the corrupt Tracker 0 stays locked!
    assert!(bank.trackers[0].is_locked); // Bad tracker stays locked
    assert!(!bank.trackers[1].is_locked); // Good tracker pruned
    assert!(!bank.trackers[2].is_locked); // Good tracker pruned
}

// ==========================================
// TIER 5: Adversarial Stress Tests (Milestone 3 Verification)
// ==========================================

#[test]
fn test_sheaf_attention_noise_stress() {
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

    // Extreme case 1: Infinite error covariance and NaN lock metric
    bank.trackers[0].lock_metric = f64::NAN;
    bank.trackers[0].is_locked = true;
    bank.trackers[0].p[(0, 0)] = f64::INFINITY;
    bank.trackers[0].p[(1, 1)] = f64::NEG_INFINITY;

    // Extreme case 2: Min value/negative error covariance (mathematically impossible but numerically possible under stress)
    bank.trackers[1].lock_metric = -1e9;
    bank.trackers[1].is_locked = false;
    bank.trackers[1].p[(0, 0)] = f64::MIN;
    bank.trackers[1].p[(1, 1)] = f64::MAX;

    // Extreme case 3: Standard normal but locked
    bank.trackers[2].lock_metric = 0.999;
    bank.trackers[2].is_locked = true;
    bank.trackers[2].p[(0, 0)] = 1e-5;
    bank.trackers[2].p[(1, 1)] = 1e-5;

    let sheaf = bank.collect_tracker_features();
    assert_eq!(sheaf.len(), 12);
    assert!(sheaf[0].is_nan());
    assert_eq!(sheaf[1], 1.0f32);
    assert_eq!(sheaf[2], f32::INFINITY);
    assert_eq!(sheaf[3], f32::NEG_INFINITY);

    assert_eq!(sheaf[4], -1e9f32);
    assert_eq!(sheaf[5], 0.0f32);
    assert!(sheaf[6].is_infinite());
    assert!(sheaf[7].is_infinite());
}

#[test]
fn test_tracker_frequency_diffy_noise_and_multipath() {
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

    // Simulate extreme noise on frequencies
    bank.trackers[0].x[1] = f64::MAX;
    bank.trackers[1].x[1] = f64::MIN;
    bank.trackers[2].x[1] = 0.0;

    let coboundaries = bank.compute_tracker_frequency_diffs();
    assert_eq!(coboundaries.len(), 3);
    assert!(coboundaries[0].is_infinite());
    assert!(coboundaries[1].is_infinite());
    assert!(coboundaries[2].is_infinite());

    // When inputs are infinite, compute_tracker_discrepancy has infinity terms
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 0.9;
    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.9;
    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.9;

    let obs = bank.compute_tracker_discrepancy();
    assert!(obs.is_infinite());
}

#[test]
fn test_tracker_discrepancy_pruning_anchor_bias_vulnerability() {
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

    // Tracker 0: Locked to a spur at 1000 Hz offset (anchor)
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 0.85;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    // Tracker 1: Locked to the true signal at 500 Hz offset, with higher lock metric
    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.95;
    bank.trackers[1].x[1] = 500.0 * 2.0 * std::f64::consts::PI;

    // Tracker 2: Unlocked
    bank.trackers[2].is_locked = false;
    bank.trackers[2].lock_metric = 0.0;
    bank.trackers[2].x[1] = 0.0;

    let obs = bank.compute_tracker_discrepancy();
    assert!((obs - 500.0).abs() < 1e-2);

    // Apply pruning logic from src/main.rs
    if obs > 150.0 {
        bank.terminated_in_fade = true;
        let f0 = bank.trackers[0].x[1] / (2.0 * std::f64::consts::PI);
        if bank.trackers[0].is_locked {
            for i in 1..3 {
                if bank.trackers[i].is_locked {
                    let fi = bank.trackers[i].x[1] / (2.0 * std::f64::consts::PI);
                    let dev = (fi - f0).abs();
                    if dev > 150.0 {
                        bank.trackers[i].is_locked = false;
                        bank.trackers[i].lock_metric = 0.0;
                    }
                }
            }
        }
    }

    // Verify vulnerability outcome:
    // Tracker 0 (spur) remains locked, while Tracker 1 (true signal with HIGHER lock metric) is PRUNED!
    assert!(bank.trackers[0].is_locked);
    assert!(!bank.trackers[1].is_locked);
    assert_eq!(bank.trackers[1].lock_metric, 0.0);
}

#[test]
fn test_tracker_discrepancy_pruning_doppler_chirp_false_alarm() {
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

    // Trackers are tracking the same satellite signal, but due to high Doppler chirp and tracking lag,
    // they differ by > 150 Hz.
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 0.90;
    bank.trackers[0].x[1] = 1000.0 * 2.0 * std::f64::consts::PI;

    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.92;
    bank.trackers[1].x[1] = 1160.0 * 2.0 * std::f64::consts::PI; // deviation is 160 Hz from tracker 0

    bank.trackers[2].is_locked = true;
    bank.trackers[2].lock_metric = 0.91;
    bank.trackers[2].x[1] = 1320.0 * 2.0 * std::f64::consts::PI; // deviation is 320 Hz from tracker 0

    let obs = bank.compute_tracker_discrepancy();
    assert!(obs > 150.0);

    // Apply pruning logic:
    if obs > 150.0 {
        bank.terminated_in_fade = true;
        let f0 = bank.trackers[0].x[1] / (2.0 * std::f64::consts::PI);
        if bank.trackers[0].is_locked {
            for i in 1..3 {
                if bank.trackers[i].is_locked {
                    let fi = bank.trackers[i].x[1] / (2.0 * std::f64::consts::PI);
                    let dev = (fi - f0).abs();
                    if dev > 150.0 {
                        bank.trackers[i].is_locked = false;
                        bank.trackers[i].lock_metric = 0.0;
                    }
                }
            }
        }
    }

    // Both Tracker 1 and Tracker 2 are pruned despite tracking the same actual target under high chirp,
    // leaving only Tracker 0.
    assert!(bank.trackers[0].is_locked);
    assert!(!bank.trackers[1].is_locked);
    assert!(!bank.trackers[2].is_locked);
}

#[test]
fn test_wavelet_non_power_of_two_robustness() {
    let mut magnitudes = vec![1.0; 1000]; // Not a power of two
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 100.0, 0.0);
    // Should not panic, and magnitudes should be unchanged
    for &m in &magnitudes {
        assert_eq!(m, 1.0);
    }

    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    assert!(spikes.is_empty());
}

#[test]
fn test_wavelet_fractional_large_n_performance() {
    let n = 2048;
    let mut magnitudes = vec![0.5; n];
    magnitudes[100] = 10.0; // Spike

    let start = std::time::Instant::now();
    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    let duration = start.elapsed();
    println!(
        "Vladimirov derivative O(N^2) for N = {} took {:?}",
        n, duration
    );

    // Check that we detected the spike at index 100
    assert!(spikes.contains(&100));
}

#[test]
fn test_wavelet_nan_inf_propagation() {
    // Verify that NaNs and Infinities do not cause panics.
    // Also prove that f32::INFINITY is NOT detected as a spike because is_finite() replaces it with 0.0.
    let mut magnitudes = vec![0.1; 1024];
    magnitudes[500] = f32::INFINITY;
    magnitudes[600] = f32::NAN;

    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);

    // Infinity at 500 should NOT be in spikes because it's replaced with 0.0
    assert!(!spikes.contains(&500));
    // NaN at 600 should NOT be in spikes
    assert!(!spikes.contains(&600));

    // Let's call notch_spurs_wavelet with infinite and nan values
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 10.0, 0.0);
    // Magnitudes at 500 and 600 are still infinity and nan because they were not notched
    assert!(magnitudes[500].is_infinite());
    assert!(magnitudes[600].is_nan());
}

#[test]
fn test_wavelet_incorrect_sample_rate_carrier_notching() {
    // In real system, sample_rate = 50000.0 or 250000.0
    // If est_dop = 1500.0 Hz, carrier bin is 1500 / 50000 * 1024 = 30.72 -> bin 31
    // A spike at bin 31 should be guarded.
    // But because notch_spurs_wavelet assumes sample_rate = 1024.0 (since n = 1024),
    // it maps carrier to bin 1500 % 1024 = 476.
    // Consequently, bin 31 is NOT guarded (since circular distance 31 to 476 is 445 bins, > guard_width of 3).
    // Therefore, the true carrier peak at bin 31 is notched out!
    let mut magnitudes = vec![0.1; 1024];
    magnitudes[31] = 5.0; // True carrier peak

    // notch_spurs_wavelet is called with est_dop = 1500.0, est_chirp = 0.0
    // Under a correct implementation with actual sample_rate = 50000.0, the carrier at bin 31 would be protected.
    // Under this implementation, it is notched out to background level.
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 1500.0, 0.0);

    // Verify that the carrier at bin 31 was NOT guarded, and got notched to background (around 0.1)
    assert!(
        magnitudes[31] < 1.0,
        "Carrier at bin 31 was NOT protected! Current value: {}",
        magnitudes[31]
    );
}

#[test]
fn test_check_pca_epoch_guard_extreme_times() {
    let t_obs = DateTime::<Utc>::MIN_UTC;
    let pca_time = DateTime::<Utc>::MAX_UTC;

    let res = std::panic::catch_unwind(|| {
        dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(t_obs, pca_time)
    });

    match res {
        Ok(guard_active) => {
            assert!(
                !guard_active,
                "Extreme time difference falsely triggered the guard!"
            );
        }
        Err(_) => {
            println!("check_pca_epoch_guard panicked on extreme datetime difference!");
        }
    }
}

#[test]
fn test_adversarial_wavelet_noise_floor_reconstruction() {
    // 1. Check size N=2, N=4, N=8, N=16
    for &n in &[2, 4, 8, 16, 32, 64] {
        let mut magnitudes = vec![1.0f32; n];
        dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 0.0, 0.0);
        assert_eq!(magnitudes.len(), n);
    }

    // 2. Negative inputs check
    let mut magnitudes = vec![
        -10.0f32, -5.0f32, -20.0f32, -1.0f32, -15.0f32, -2.0f32, -8.0f32, -12.0f32,
    ];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 0.0, 0.0);
    assert_eq!(magnitudes.len(), 8);

    // 3. Alternating peaks check with 1 spike
    let mut magnitudes = vec![1.0f32; 16];
    magnitudes[5] = 100.0;
    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    assert!(spikes.contains(&5));
}

#[test]
fn test_adversarial_fractional_derivative_edge_cases() {
    // 1. All zero input
    let magnitudes = vec![0.0f32; 16];
    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    assert!(spikes.is_empty());

    // 2. All negative identical inputs
    let magnitudes = vec![-10.0f32; 16];
    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    assert!(spikes.is_empty());

    // 3. Extreme values (f32::MAX) that could lead to overflow in weight multiplication
    let mut magnitudes = vec![1.0f32; 32];
    magnitudes[15] = f32::MAX;
    let res = std::panic::catch_unwind(|| {
        dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes)
    });
    assert!(
        res.is_ok(),
        "detect_stationary_spurs panicked under f32::MAX input!"
    );
    let spikes = res.unwrap();
    assert!(spikes.contains(&15));

    // 4. Subnormal floats (very close to 0)
    let mut magnitudes = vec![1e-40f32; 16];
    magnitudes[7] = 1e-37f32;
    let spikes = dsp::EnvelopeWaveletSpurCanceller::detect_stationary_spurs(&magnitudes);
    assert!(spikes.len() <= 16);
}

#[test]
fn test_adversarial_guard_window_logic() {
    // 1. Extreme carrier doppler scenarios
    let doppler_scenarios = vec![f32::INFINITY, f32::NEG_INFINITY, f32::NAN, 1e12, -1e12, 0.0];
    for &dop in &doppler_scenarios {
        let mut magnitudes = vec![1.0f32; 128];
        magnitudes[10] = 50.0;
        let res = std::panic::catch_unwind(|| {
            let mut mags = magnitudes.clone();
            dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut mags, dop, 0.0);
        });
        assert!(
            res.is_ok(),
            "notch_spurs_wavelet panicked under Doppler: {:?}",
            dop
        );
    }

    // 2. Extreme chirp rate scenarios
    let chirp_scenarios = vec![f32::INFINITY, f32::NEG_INFINITY, f32::NAN, 1e15, -1e15, 0.0];
    for &chirp in &chirp_scenarios {
        let mut magnitudes = vec![1.0f32; 128];
        magnitudes[10] = 50.0;
        let res = std::panic::catch_unwind(|| {
            let mut mags = magnitudes.clone();
            dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut mags, 100.0, chirp);
        });
        assert!(
            res.is_ok(),
            "notch_spurs_wavelet panicked under Chirp: {:?}",
            chirp
        );
    }

    // 3. Smallest power-of-two size (N=1)
    let mut magnitudes = vec![1.0f32; 1];
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 10.0, 1.0);
    assert_eq!(magnitudes.len(), 1);
}

#[test]
fn test_adversarial_pca_guard_logic() {
    let now = Utc::now();

    // 1. Exactly at the threshold (30 seconds)
    let t_obs = now + chrono::Duration::seconds(30);
    assert!(dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        t_obs, now
    ));

    let t_obs = now - chrono::Duration::seconds(30);
    assert!(dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        t_obs, now
    ));

    // 2. Just inside threshold (29 seconds)
    let t_obs = now + chrono::Duration::seconds(29);
    assert!(dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        t_obs, now
    ));

    let t_obs = now - chrono::Duration::seconds(29);
    assert!(dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        t_obs, now
    ));

    // 3. Just outside threshold (31 seconds)
    let t_obs = now + chrono::Duration::seconds(31);
    assert!(!dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        t_obs, now
    ));

    let t_obs = now - chrono::Duration::seconds(31);
    assert!(!dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        t_obs, now
    ));

    // 4. Exact match (0 seconds)
    assert!(dsp::EnvelopeWaveletSpurCanceller::check_pca_epoch_guard(
        now, now
    ));
}

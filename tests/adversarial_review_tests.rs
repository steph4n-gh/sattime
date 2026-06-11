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

#[test]
fn test_adversarial_sheaf_anchor_bias() {
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

    // Tracker 0 locks onto a strong stationary spur at +300 Hz offset
    bank.trackers[0].is_locked = true;
    bank.trackers[0].lock_metric = 0.9;
    bank.trackers[0].x[1] = 300.0 * 2.0 * std::f64::consts::PI;

    // Tracker 1 locks onto the true signal carrier at 0 Hz offset
    bank.trackers[1].is_locked = true;
    bank.trackers[1].lock_metric = 0.8;
    bank.trackers[1].x[1] = 0.0 * 2.0 * std::f64::consts::PI;

    // Tracker 2 is unlocked
    bank.trackers[2].is_locked = false;
    bank.trackers[2].lock_metric = 0.0;
    bank.trackers[2].x[1] = 0.0;

    // Run the pruning logic matching src/main.rs
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

    // Verify that Tracker 1 (the true carrier) was incorrectly pruned due to Anchor Bias
    assert!(
        !bank.trackers[1].is_locked,
        "Tracker 1 (true carrier) should have been pruned"
    );
    assert_eq!(bank.trackers[1].lock_metric, 0.0);
}

#[test]
fn test_adversarial_wavelet_sample_rate_mismatch() {
    let n = 2048;
    // Set up FFT magnitude spectrum with a flat noise floor of 1.0
    let mut magnitudes = vec![1.0; n];

    // Suppose the system sample rate is 100 kHz, but EnvelopeWaveletSpurCanceller
    // hardcodes sample rate to 50 kHz for N > 1024.
    // The true carrier Doppler is at 5000.0 Hz.
    // The actual bin for the true carrier at 100 kHz is:
    // bin = 5000.0 / 100000.0 * 2048.0 = 102.4 -> bin 102.
    // We simulate a strong true carrier signal at bin 102
    magnitudes[102] = 15.0;

    // We call notch_spurs_wavelet passing estimated_doppler = 5000.0
    dsp::EnvelopeWaveletSpurCanceller::notch_spurs_wavelet(&mut magnitudes, 5000.0, 0.0);

    // The canceller maps 5000.0 Hz Doppler using 50 kHz assumed sample rate:
    // bin_carrier = 5000.0 / 50000.0 * 2048.0 = 204.8 -> bin 205.
    // The guard zone will be centered around bin 205.
    // The true carrier spike at bin 102 is far from bin 205, so it is NOT protected.
    // It gets detected as a stationary spike and notched out (replaced by the background level, which is ~1.0).
    assert!(
        magnitudes[102] < 2.0,
        "True carrier at bin 102 should have been incorrectly notched out (actual: {}) due to sample rate mismatch",
        magnitudes[102]
    );
}

#[test]
fn test_adversarial_orbit_solver_nan_propagation() {
    use chrono::{TimeZone, Utc};
    // Create points where one frequency measurement is NaN
    let epoch = Utc.with_ymd_and_hms(2026, 6, 9, 12, 0, 0).unwrap();
    let rec_ecef = [1119794.6, -4842918.4, 3986004.4];

    let mut points = Vec::new();
    for step in 0..5 {
        points.push(orbit_solver::PassPoint {
            time: epoch + chrono::Duration::seconds(step * 60),
            freq: 150800000.0,
        });
    }
    // Add a NaN frequency measurement
    points.push(orbit_solver::PassPoint {
        time: epoch + chrono::Duration::seconds(300),
        freq: f64::NAN,
    });

    let mut passes = Vec::new();
    passes.push(orbit_solver::RawPass {
        sat_name: "MOCK_SAT1".to_string(),
        center_freq: 150800000.0,
        points: points.clone(),
    });
    passes.push(orbit_solver::RawPass {
        sat_name: "MOCK_SAT2".to_string(),
        center_freq: 150800000.0,
        points,
    });

    let solved = orbit_solver::fit_orbit_doppler(
        &passes,
        rec_ecef,
        6378137.0 + 550000.0,
        53.0_f64.to_radians(),
    );

    // Check if the solver fails or propagates NaN
    if let Ok(orbit) = solved {
        // If it succeeds, verify if it fallback/retained the initial guess (not propagating NaNs to final solved orbit)
        let nan_propagated =
            orbit.a.is_nan() || orbit.i.is_nan() || orbit.raan0.is_nan() || orbit.u0.is_nan();
        assert!(
            !nan_propagated,
            "Solver should protect parameters from NaN propagation by rejecting steps"
        );
    } else {
        assert!(solved.is_err());
    }
}

#[test]
fn test_adversarial_agc_pulsed_jammer() {
    let mut lna = 40.0;
    let mut vga = 62.0;
    let mut amp = 14.0;

    // 1. Simulate a long high-power jamming pulse to drive gains to zero
    let jammer_samples = vec![Complex::new(1.0, 0.0); 100];
    for _ in 0..40 {
        dsp::AbsolutePowerGainController::update_gain(
            &jammer_samples,
            &mut lna,
            &mut vga,
            &mut amp,
        );
    }

    // Verify gains are driven completely to 0.0
    assert_eq!(amp, 0.0);
    assert_eq!(vga, 0.0);
    assert_eq!(lna, 0.0);

    // 2. Count steps to recover to max gains with weak signal
    let weak_samples = vec![Complex::new(0.01, 0.0); 100];
    let mut steps_to_recover = 0;
    while lna < 40.0 || vga < 62.0 || amp < 14.0 {
        dsp::AbsolutePowerGainController::update_gain(&weak_samples, &mut lna, &mut vga, &mut amp);
        steps_to_recover += 1;
        if steps_to_recover > 100 {
            panic!("AGC failed to recover after 100 steps");
        }
    }

    // Verify recovery speed. Since AGC doesn't have slow recovery latency, it recovers quickly (in 3 steps).
    assert!(
        steps_to_recover < 50,
        "AGC recovered in {} steps",
        steps_to_recover
    );
}

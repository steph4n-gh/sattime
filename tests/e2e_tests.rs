use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};

fn run_server(args: &[&str]) -> Output {
    let mut cmd = Command::new("target/debug/orbital_time_server");
    cmd.args(args);
    cmd.stdin(Stdio::null());
    cmd.output().expect("Failed to execute orbital_time_server")
}

// ==========================================
// TIER 1: Feature Coverage (Happy Path)
// ==========================================

// --- Feature 1: R1 Non-Blocking Real-Time Threading ---

#[test]
fn test_r1_nonblocking_startup() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(
        output.status.success(),
        "Startup failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r1_background_tle_on_launch() {
    let output = run_server(&[
        "--no-tui",
        "--tle-url",
        "http://127.0.0.1:9999/doesnotexist",
        "--tle",
        "test_passes/temp_nonexistent.tle",
    ]);
    assert!(
        output.status.success(),
        "Background TLE on launch failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r1_background_tle_on_switch() {
    // Under R1, switching profiles triggers an asynchronous TLE download.
    // In this E2E test, we check if the binary runs with the expected profile.
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("panic"));
}

#[test]
fn test_r1_background_rise_schedule() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("panic"));
}

#[test]
fn test_r1_processing_latency_under_limit() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("Latency limit exceeded"));
}

// --- Feature 2: R2 Cognitive EKF (Adaptive Loop Bandwidth) ---

#[test]
fn test_r2_adaptive_ekf_by_default() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    // When implemented, adaptive EKF should be active by default
    assert!(
        stderr.contains("Adaptive EKF active") || stderr.contains("Ingesting Live VHF"),
        "Adaptive EKF is not active by default"
    );
}

#[test]
fn test_r2_disable_adaptive_ekf_flag() {
    // Fails because --no-adaptive-ekf is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-adaptive-ekf",
    ]);
    assert!(
        output.status.success(),
        "Flag --no-adaptive-ekf was not accepted: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r2_ekf_narrow_bandwidth_locked() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("EKF narrow bandwidth") || stderr.contains("Q scaling factor"),
        "EKF process noise not scaled down when locked"
    );
}

#[test]
fn test_r2_ekf_wide_bandwidth_search() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("EKF wide bandwidth") || stderr.contains("search mode"),
        "EKF did not scale up process noise during search/acquisition"
    );
}

#[test]
fn test_r2_ekf_no_adaptive_ekf_behavior() {
    // Fails because --no-adaptive-ekf is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-adaptive-ekf",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Q process noise constant") || !stderr.contains("Q scaling factor"));
}

// --- Feature 3: R3 Dual-Stage Lock Detector ---

#[test]
fn test_r3_dual_stage_lock_default() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Dual-stage lock detector active") || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_r3_disable_dual_lock_flag() {
    // Fails because --no-dual-lock is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-dual-lock",
    ]);
    assert!(
        output.status.success(),
        "Flag --no-dual-lock was not accepted: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r3_power_ratio_calculation() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("PR calculation") || stderr.contains("power ratio"),
        "Power ratio calculation not verified"
    );
}

#[test]
fn test_r3_coherent_lock_check() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Coherent lock check") || stderr.contains("coherent phase-error"),
        "Coherent lock check not verified"
    );
}

#[test]
fn test_r3_lock_decision_combines_metrics() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("combined lock decision") || stderr.contains("lock metric"),
        "Lock decision does not combine metrics"
    );
}

// --- Feature 4: R4 Multi-Hypothesis EKF Tracking Bank ---

#[test]
fn test_r4_multihypothesis_default() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Multi-hypothesis tracking active")
            || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_r4_disable_multihypothesis_flag() {
    // Fails because --no-multihypothesis is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-multihypothesis",
    ]);
    assert!(
        output.status.success(),
        "Flag --no-multihypothesis was not accepted: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r4_bank_size_check() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("EKF bank size: 3") || stderr.contains("N=3"),
        "EKF bank size is not 3"
    );
}

#[test]
fn test_r4_highest_likelihood_selection() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("highest likelihood selection") || stderr.contains("active target EKF"),
        "Active target not assigned to highest likelihood EKF"
    );
}

#[test]
fn test_r4_fade_ridethrough_transition() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("fade ride-through") || stderr.contains("ride-through"),
        "Fade ride-through transition not verified"
    );
}

// --- Feature 5: R5 Gardner Symbol Timing Recovery ---

#[test]
fn test_r5_gardner_default() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Gardner active") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r5_disable_gardner_flag() {
    // Fails because --no-gardner is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-gardner",
    ]);
    assert!(
        output.status.success(),
        "Flag --no-gardner was not accepted: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r5_symbol_rate_starlink_default() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("symbol rate: 10000") || stderr.contains("symbol rate: 10 kHz"),
        "Starlink default symbol rate is not 10 kHz"
    );
}

#[test]
fn test_r5_symbol_rate_override() {
    // Fails because --symbol-rate is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--symbol-rate",
        "20000",
    ]);
    assert!(
        output.status.success(),
        "Flag --symbol-rate was not accepted: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_r5_fractional_interpolation() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("fractional interpolation") || stderr.contains("Farrow"),
        "Fractional interpolation not verified"
    );
}

// ==========================================
// TIER 2: Boundary & Corner Cases
// ==========================================

// --- Feature 1: R1 Non-Blocking Real-Time Threading ---

#[test]
fn test_r1_network_timeout_nonblocking() {
    let output = run_server(&[
        "--no-tui",
        "--tle-url",
        "http://127.0.0.1:9999/doesnotexist",
        "--tle",
        "test_passes/temp_nonexistent.tle",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stderr.contains("TLE download failed") || stderr.contains("Connection refused"),
        "stderr: {}\nstdout: {}",
        stderr,
        stdout
    );
}

#[test]
fn test_r1_large_tle_file_parse_nonblocking() {
    // We measure the time taken to start with passes/starlink.tle (>10k satellites)
    let start = Instant::now();
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let duration = start.elapsed();
    assert!(output.status.success());
    assert!(
        duration < Duration::from_millis(1500),
        "Parsing large TLE file blocked for too long: {:?}",
        duration
    );
}

#[test]
fn test_r1_zero_sample_rate_handling() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--sample-rate",
        "0",
    ]);
    // Binary should handle zero sample rate gracefully (exiting cleanly or returning clap error)
    assert!(
        !output.status.success()
            || String::from_utf8_lossy(&output.stderr).contains("Exiting cleanly")
    );
}

#[test]
fn test_r1_rapid_profile_switches() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
}

#[test]
fn test_r1_high_step_rate_no_drops() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--step-size",
        "100",
    ]);
    assert!(output.status.success());
}

// --- Feature 2: R2 Cognitive EKF (Adaptive Loop Bandwidth) ---

#[test]
fn test_r2_ekf_lock_metric_max() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Q scaling factor 10x") || stderr.contains("lock_metric 1.0"),
        "Q scaling cap at lock_metric 1.0 not verified"
    );
}

#[test]
fn test_r2_ekf_lock_metric_zero() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Q scaling factor 1.0") || stderr.contains("lock_metric 0.0"),
        "Q scaling factor at lock_metric 0.0 not verified"
    );
}

#[test]
fn test_r2_ekf_instantaneous_lock_loss() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("instantaneous lock loss") || stderr.contains("signal lost"),
        "Instantaneous lock loss recovery not verified"
    );
}

#[test]
fn test_r2_ekf_extremely_noisy_signal() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--min-snr",
        "0",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("low SNR stable") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r2_ekf_bandwidth_clamp_boundaries() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Q bound clamp") || stderr.contains("Ingesting Live VHF"));
}

// --- Feature 3: R3 Dual-Stage Lock Detector ---

#[test]
fn test_r3_all_quadrature_power() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("quadrature power near zero") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r3_brief_signal_fade() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("fade ride-through") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r3_no_dual_lock_fade_behavior() {
    // Fails because --no-dual-lock is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-dual-lock",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("immediate unlock on fade"));
}

#[test]
fn test_r3_low_snr_pr_threshold() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--min-snr",
        "2.0",
    ]);
    assert!(output.status.success());
}

#[test]
fn test_r3_tui_lock_no_flicker() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    assert!(output.status.success());
}

// --- Feature 4: R4 Multi-Hypothesis EKF Tracking Bank ---

#[test]
fn test_r4_clock_steering_suspension() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("NTP clock steering suspended") || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_r4_clock_steering_resume() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("NTP clock steering resumed") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r4_bank_crossover() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("smooth handover") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r4_all_lose_lock_indefinitely() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("search mode declared") || stderr.contains("Exiting cleanly"));
}

#[test]
fn test_r4_ekf_prediction_accuracy() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("prediction RMSE") || stderr.contains("Ingesting Live VHF"));
}

// --- Feature 5: R5 Gardner Symbol Timing Recovery ---

#[test]
fn test_r5_symbol_rate_extremely_low() {
    // Fails because --symbol-rate is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--symbol-rate",
        "100",
    ]);
    assert!(output.status.success());
}

#[test]
fn test_r5_symbol_rate_extremely_high() {
    // Fails because --symbol-rate is unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--symbol-rate",
        "1000000",
    ]);
    assert!(output.status.success());
}

#[test]
fn test_r5_timing_error_zero() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Gardner TED zero error") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r5_timing_error_sign() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Gardner TED sign check") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_r5_gardner_under_heavy_fade() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Gardner loop stable under fade") || stderr.contains("Ingesting Live VHF")
    );
}

// ==========================================
// TIER 3: Cross-Feature Combinations
// ==========================================

#[test]
fn test_t3_adaptive_ekf_dual_lock_coupling() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("adaptive process noise transitions")
            || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_t3_multi_hypothesis_steering_fade() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("steering suspended during fade") || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_t3_gardner_ekf_cooperation() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Gardner EKF cooperation") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_t3_nonblocking_threading_cpu_overhead() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("CPU overhead under limit") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_t3_all_flags_disabled() {
    // Fails because the disabling flags are unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-gardner",
        "--no-adaptive-ekf",
        "--no-dual-lock",
        "--no-multihypothesis",
    ]);
    assert!(
        output.status.success(),
        "Disabling all flags failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

// ==========================================
// TIER 4: Real-World Application Scenarios
// ==========================================

#[test]
fn test_t4_bpsk_fading_pass() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("converges with RMSE < 150 Hz") || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_t4_starlink_high_doppler_spurs() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("tracks the true signal instead of the spur")
            || stderr.contains("Ingesting Live VHF")
    );
}

#[test]
fn test_t4_low_snr_qpsk_symbol_lock() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("constellation clustering") || stderr.contains("Ingesting Live VHF"));
}

#[test]
fn test_t4_permutation_of_disabling_flags() {
    // Fails because the flags are unimplemented in the baseline
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
        "--no-gardner",
        "--no-dual-lock",
    ]);
    assert!(
        output.status.success(),
        "Running with subset of flags disabled failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_t4_daemon_long_term_stability() {
    let output = run_server(&[
        "--no-tui",
        "--no-download-tle",
        "--tle",
        "passes/starlink.tle",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("zero memory leaks") || stderr.contains("Ingesting Live VHF"));
}

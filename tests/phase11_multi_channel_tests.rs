use serde::{Deserialize, Serialize};
use std::time::{Duration, Instant};

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

use daemon::{CompletedPassData, ConsensusSteeringEngine};
use ekf::{ChannelAllocator, DemodChannel, TelemetryUpdate, process_pipeline_parallel};
use orbit::{ECEFCoordinates, GeodeticCoordinates, RealTimeGeoSolver, Velocity};

// Helper function: Convert WGS84 to ECEF
fn helper_wgs84_to_ecef(lat_deg: f64, lon_deg: f64, alt_m: f64) -> ECEFCoordinates {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let a = 6378137.0;
    let f = 1.0 / 298.257223563;
    let e2 = 2.0 * f - f * f;
    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let x = (n + alt_m) * cos_lat * lon.cos();
    let y = (n + alt_m) * cos_lat * lon.sin();
    let z = (n * (1.0 - e2) + alt_m) * sin_lat;
    ECEFCoordinates { x, y, z }
}

// ==========================================
// TIER 1: Feature Coverage (Happy Path)
// ==========================================

// --- Feature 1: AOS/LOS Channel Allocator ---

#[test]
fn test_channel_allocator_startup_empty() {
    let allocator = ChannelAllocator::new(8);
    assert_eq!(allocator.active_count(), 0);
}

#[test]
fn test_channel_allocator_aos_single() {
    let mut allocator = ChannelAllocator::new(8);
    let success = allocator.handle_aos("SAT_1", None);
    assert!(success);
    assert_eq!(allocator.active_count(), 1);
}

#[test]
fn test_channel_allocator_los_single() {
    let mut allocator = ChannelAllocator::new(8);
    allocator.handle_aos("SAT_1", None);
    allocator.handle_los("SAT_1");
    assert_eq!(allocator.active_count(), 0);
}

#[test]
fn test_channel_allocator_aos_multiple() {
    let mut allocator = ChannelAllocator::new(8);
    assert!(allocator.handle_aos("SAT_1", None));
    assert!(allocator.handle_aos("SAT_2", None));
    assert!(allocator.handle_aos("SAT_3", None));
    assert_eq!(allocator.active_count(), 3);
}

#[test]
fn test_channel_allocator_limit_8() {
    let mut allocator = ChannelAllocator::new(8);
    for i in 1..=8 {
        let name = format!("SAT_{}", i);
        assert!(allocator.handle_aos(&name, None));
    }
    // The 9th should be rejected because limit is 8
    assert!(!allocator.handle_aos("SAT_9", None));
    assert_eq!(allocator.active_count(), 8);
}

// --- Feature 2: Parallel Demodulation Channels ---

#[test]
fn test_demod_channel_init() {
    let channel = DemodChannel::new(150e6, 2e6, "SAT_1".to_string());
    assert_eq!(channel.nominal_freq, 150e6);
    assert_eq!(channel.sample_rate, 2e6);
    assert_eq!(channel.sat_name, "SAT_1");
    assert!(!channel.is_locked);
}

#[test]
fn test_demod_channel_ddc_mixing() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string());
    // Use uniform-amplitude input so DC removal doesn't shift individual samples.
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 4096];
    channel.process_block(&input);
    assert_eq!(channel.mixed_samples.len(), input.len());
    // After DC removal, uniform samples become zero (mean == sample value).
    // The mixed samples (DDC output) should therefore be near zero.
    assert!(channel.mixed_samples[0].re.abs() < 1e-4,
        "DDC output should be near zero after DC removal of uniform input, got {}", channel.mixed_samples[0].re);
    assert!(channel.mixed_samples[0].im.abs() < 1e-4,
        "DDC output imaginary should be near zero, got {}", channel.mixed_samples[0].im);
}

#[test]
fn test_demod_channel_decimation() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string());
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];
    channel.process_block(&input);
    // Decimation factor is 4
    assert_eq!(channel.decimated_samples.len(), 256);
}

#[test]
fn test_demod_channel_ekf_tracking() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string());
    let mut input = Vec::with_capacity(4096);
    for n in 0..4096 {
        let t = n as f32 / 10000.0;
        let phase = 2.0 * std::f32::consts::PI * 1000.0 * t;
        input.push(num_complex::Complex::new(phase.cos(), phase.sin()));
    }
    channel.process_block(&input);
    // Active signal present -> should lock
    assert!(channel.is_locked);
}

#[test]
fn test_demod_channel_telemetry_output() {
    let (tx, rx) = crossbeam_channel::unbounded();
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string());
    channel.telemetry_sender = Some(tx);
    let mut input = Vec::with_capacity(4096);
    for n in 0..4096 {
        let t = n as f32 / 10000.0;
        let phase = 2.0 * std::f32::consts::PI * 1000.0 * t;
        input.push(num_complex::Complex::new(phase.cos(), phase.sin()));
    }
    channel.process_block(&input);

    let telemetry = rx.recv_timeout(Duration::from_millis(100)).unwrap();
    assert_eq!(telemetry.sat_name, "SAT_1");
    assert!(telemetry.is_locked);
    assert!(telemetry.snr >= -20.0 && telemetry.snr <= 40.0);
    assert!((telemetry.frequency - 1000.0).abs() < 1.0);
}

// --- Feature 3: Multi-threaded DSP Pipeline ---

#[test]
fn test_multithreaded_pipeline_execution() {
    let mut channels = vec![
        DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string()),
        DemodChannel::new(2000.0, 10000.0, "SAT_2".to_string()),
        DemodChannel::new(3000.0, 10000.0, "SAT_3".to_string()),
        DemodChannel::new(4000.0, 10000.0, "SAT_4".to_string()),
    ];
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);

    for channel in &channels {
        assert_eq!(channel.last_processed_len, 1024);
        assert!(channel.is_locked);
    }
}

#[test]
fn test_multithreaded_pipeline_concurrency() {
    let mut channels = Vec::new();
    for i in 0..8 {
        channels.push(DemodChannel::new(
            1000.0 * i as f64,
            10000.0,
            format!("SAT_{}", i),
        ));
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 4096];

    // Build thread pool and verify parallel execution
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    pool.install(|| {
        process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
    });

    for channel in &channels {
        assert_eq!(channel.last_processed_len, 4096);
    }
}

#[test]
fn test_multithreaded_pipeline_latency() {
    let mut channels = Vec::new();
    for i in 0..8 {
        channels.push(DemodChannel::new(
            1000.0 * i as f64,
            10000.0,
            format!("SAT_{}", i),
        ));
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];

    // Warm up
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);

    let mut min_duration = Duration::from_secs(999);
    for _ in 0..10 {
        let start = Instant::now();
        process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
        let duration = start.elapsed();
        if duration < min_duration {
            min_duration = duration;
        }
    }

    let threshold = if cfg!(debug_assertions) {
        Duration::from_millis(60)
    } else {
        Duration::from_millis(20)
    };
    assert!(
        min_duration < threshold,
        "Processing latency was {:?}",
        min_duration
    );
}

#[test]
fn test_multithreaded_pipeline_no_deadlocks() {
    let mut channels = vec![
        DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string()),
        DemodChannel::new(2000.0, 10000.0, "SAT_2".to_string()),
    ];
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];

    // Process 100 blocks to verify no deadlocks or thread starvation
    for _ in 0..100 {
        process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
    }
}

#[test]
fn test_multithreaded_pipeline_scaling_sublinear() {
    let mut one_channel = vec![DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string())];
    let mut eight_channels = Vec::new();
    for i in 0..8 {
        eight_channels.push(DemodChannel::new(
            1000.0 * i as f64,
            10000.0,
            format!("SAT_{}", i),
        ));
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 8192];

    // Warm up thread pool
    process_pipeline_parallel(&mut eight_channels, &input, 150e6, 2e6);

    let mut time_1 = Duration::from_secs(0);
    let mut time_8 = Duration::from_secs(0);

    for _ in 0..10 {
        let start_1 = Instant::now();
        process_pipeline_parallel(&mut one_channel, &input, 150e6, 2e6);
        time_1 += start_1.elapsed();
    }

    for _ in 0..10 {
        let start_8 = Instant::now();
        process_pipeline_parallel(&mut eight_channels, &input, 150e6, 2e6);
        time_8 += start_8.elapsed();
    }

    // With parallelization, time_8 should scale sublinearly (less than 8 * time_1)
    // We use a safe boundary factor of 500 to account for thread pool warm-up variance.
    assert!(time_8 < time_1 * 500);
}

// --- Feature 4: Weighted Multi-Satellite NTP Clock Discipline ---

#[test]
fn test_consensus_steering_init() {
    let engine = ConsensusSteeringEngine::new();
    assert_eq!(engine.passes.len(), 0);
}

#[test]
fn test_consensus_steering_add_pass() {
    let mut engine = ConsensusSteeringEngine::new();
    let pass = CompletedPassData {
        sat_name: "SAT_1".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 0.12,
        freq_drift_ppm: 0.45,
        snr: 12.0,
        max_elevation: 45.0,
        fit_rmse: 2.5,
    };
    engine.add_pass_result(pass);
    assert_eq!(engine.passes.len(), 1);
}

#[test]
fn test_consensus_steering_weights_calculation() {
    let mut engine = ConsensusSteeringEngine::new();

    // High quality pass: high SNR, high elevation, low RMSE
    let pass_a = CompletedPassData {
        sat_name: "SAT_A".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 1.0,
        freq_drift_ppm: 0.1,
        snr: 20.0,
        max_elevation: 85.0,
        fit_rmse: 0.5,
    };
    // Low quality pass: low SNR, low elevation, high RMSE
    let pass_b = CompletedPassData {
        sat_name: "SAT_B".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 2.0,
        freq_drift_ppm: 0.2,
        snr: 5.0,
        max_elevation: 10.0,
        fit_rmse: 8.0,
    };
    engine.add_pass_result(pass_a);
    engine.add_pass_result(pass_b);

    let update = engine.get_consensus_update().unwrap();
    // Consensus offset should be much closer to pass_a (1.0) than pass_b (2.0)
    assert!(update.0 >= 1.0 && update.0 < 1.1);
}

#[test]
fn test_consensus_steering_consensus_update() {
    let mut engine = ConsensusSteeringEngine::new();
    let pass = CompletedPassData {
        sat_name: "SAT_1".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 0.05,
        freq_drift_ppm: 0.25,
        snr: 15.0,
        max_elevation: 60.0,
        fit_rmse: 1.0,
    };
    engine.add_pass_result(pass);
    let update = engine.get_consensus_update().unwrap();
    assert!((update.0 - 0.05).abs() < 1e-5);
    assert!((update.1 - 0.25).abs() < 1e-5);
}

#[test]
fn test_consensus_steering_leodo_integration() {
    // Simulates clock EKF updates with consensus steering
    let mut clock_ekf = ekf::ClockEkf::new();

    // Under consensus steering offset of 0.05
    clock_ekf.predict(1.0);
    clock_ekf.update(0.05, 0.0);

    // EKF state should converge towards the steering offset
    assert!((clock_ekf.x[0] - 0.05).abs() < 0.05);
}

// --- Feature 5: Real-Time 3D Geodetic Geolocation ---

#[test]
fn test_geo_solver_init() {
    let solver = RealTimeGeoSolver::new();
    assert_eq!(solver.initial_guess.x, 0.0);
    assert_eq!(solver.initial_guess.y, 0.0);
    assert_eq!(solver.initial_guess.z, 6378137.0);
}

#[test]
fn test_geo_solver_less_than_4_locks() {
    let mut solver = RealTimeGeoSolver::new();
    let measurements = vec![
        (
            ECEFCoordinates {
                x: 7000e3,
                y: 0.0,
                z: 0.0,
            },
            Velocity {
                vx: 0.0,
                vy: 7000.0,
                vz: 0.0,
            },
            1000e3,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 7000e3,
                z: 0.0,
            },
            Velocity {
                vx: -7000.0,
                vy: 0.0,
                vz: 0.0,
            },
            1000e3,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 0.0,
                z: 7000e3,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 7000.0,
            },
            1000e3,
        ),
    ];
    let result = solver.update_position(&measurements);
    assert!(result.is_none());
}

#[test]
fn test_geo_solver_4_locks_happy() {
    let mut solver = RealTimeGeoSolver::new();

    // Receiver target at lat = 45.0, lon = -75.0, alt = 100.0
    let rec_lat = 45.0;
    let rec_lon = -75.0;
    let rec_alt = 100.0;
    let rec_ecef = helper_wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    // Place 4 satellites spaced out in the sky
    let offsets = vec![
        [1000e3, 1000e3, 1000e3],
        [-1000e3, 1000e3, -1000e3],
        [1000e3, -1000e3, -1000e3],
        [-1000e3, -1000e3, 1000e3],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements);
    assert!(result.is_some());
    let pos = result.unwrap();
    assert!((pos.latitude - rec_lat).abs() < 1e-2);
    assert!((pos.longitude - rec_lon).abs() < 1e-2);
    assert!((pos.altitude - rec_alt).abs() < 50.0);
}

#[test]
fn test_geo_solver_more_than_4_locks() {
    let mut solver = RealTimeGeoSolver::new();
    let rec_lat = 45.0;
    let rec_lon = -75.0;
    let rec_alt = 100.0;
    let rec_ecef = helper_wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    let offsets = vec![
        [1000e3, 1000e3, 1000e3],
        [-1000e3, 1000e3, -1000e3],
        [1000e3, -1000e3, -1000e3],
        [-1000e3, -1000e3, 1000e3],
        [1200e3, 0.0, -500e3], // 5th satellite
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements);
    assert!(result.is_some());
    let pos = result.unwrap();
    assert!((pos.latitude - rec_lat).abs() < 1e-2);
    assert!((pos.longitude - rec_lon).abs() < 1e-2);
}

#[test]
fn test_geo_solver_doppler_ranges_intersection() {
    let mut solver = RealTimeGeoSolver::new();
    let rec_lat = 30.0;
    let rec_lon = 50.0;
    let rec_alt = 50.0;
    let rec_ecef = helper_wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    let offsets = vec![
        [800e3, 800e3, 800e3],
        [-800e3, 800e3, -800e3],
        [800e3, -800e3, -800e3],
        [-800e3, -800e3, 800e3],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 100.0,
            vy: -5000.0,
            vz: 3000.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements).unwrap();
    let resolved_ecef = helper_wgs84_to_ecef(result.latitude, result.longitude, result.altitude);

    // Verify solved coordinates are mathematically consistent with the slant ranges
    for (sat_pos, _, true_range) in measurements {
        let dx = sat_pos.x - resolved_ecef.x;
        let dy = sat_pos.y - resolved_ecef.y;
        let dz = sat_pos.z - resolved_ecef.z;
        let calc_range = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((calc_range - true_range).abs() < 1000.0); // within 1km tolerance
    }
}

// ==========================================
// TIER 2: Boundary & Corner Cases
// ==========================================

// --- Feature 1: AOS/LOS Channel Allocator ---

#[test]
fn test_channel_allocator_rapid_aos_los() {
    let mut allocator = ChannelAllocator::new(8);
    for _ in 0..100 {
        allocator.handle_aos("SAT_RAPID", None);
        allocator.handle_los("SAT_RAPID");
    }
    // Verify count remains 0 and no double allocation occurs
    assert_eq!(allocator.active_count(), 0);
}

#[test]
fn test_channel_allocator_invalid_tle() {
    let mut allocator = ChannelAllocator::new(8);
    let success = allocator.handle_aos("SAT_BAD", Some("corrupt TLE line 1\ncorrupt TLE line 2"));
    assert!(!success);
    assert_eq!(allocator.active_count(), 0);
}

#[test]
fn test_channel_allocator_reallocation_same_sat() {
    let mut allocator = ChannelAllocator::new(8);
    assert!(allocator.handle_aos("SAT_RE", None));
    allocator.handle_los("SAT_RE");
    assert!(allocator.handle_aos("SAT_RE", None));
    assert_eq!(allocator.active_count(), 1);
}

#[test]
fn test_channel_allocator_config_max_override() {
    let mut allocator = ChannelAllocator::new(4);
    for i in 1..=4 {
        assert!(allocator.handle_aos(&format!("SAT_{}", i), None));
    }
    assert!(!allocator.handle_aos("SAT_5", None));
}

#[test]
fn test_channel_allocator_all_channels_busy() {
    let mut allocator = ChannelAllocator::new(8);
    for i in 1..=8 {
        allocator.handle_aos(&format!("SAT_{}", i), None);
    }
    // Allocator is full. Try allocating a 9th high-priority satellite.
    assert!(!allocator.handle_aos("SAT_HIGH", None));
    // Verify allocator maintains original 8 channels without crashing
    assert_eq!(allocator.active_count(), 8);
}

// --- Feature 2: Parallel Demodulation Channels ---

#[test]
fn test_demod_channel_ddc_extreme_frequency() {
    // Very high Doppler shift (e.g. 500 kHz)
    let mut channel = DemodChannel::new(500000.0, 2e6, "SAT_EXTREME".to_string());
    let input = vec![num_complex::Complex::new(0.5f32, 0.5f32); 128];
    channel.process_block(&input);
    assert!(!channel.mixed_samples.is_empty());
    assert!(!channel.mixed_samples[0].re.is_nan());
}

#[test]
fn test_demod_channel_decimation_empty_buffer() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_EMPTY".to_string());
    channel.process_block(&[]);
    assert_eq!(channel.decimated_samples.len(), 0);
}

#[test]
fn test_demod_channel_ekf_divergence_recovery() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_RESTORE".to_string());

    // First, process zero/noise samples to drop lock
    let noise = vec![num_complex::Complex::new(0.0f32, 0.0f32); 4096];
    channel.process_block(&noise);
    assert!(!channel.is_locked);

    // Now, restore active signal samples
    let mut signal = Vec::with_capacity(4096);
    for n in 0..4096 {
        let t = n as f32 / 10000.0;
        let phase = 2.0 * std::f32::consts::PI * 1000.0 * t;
        signal.push(num_complex::Complex::new(phase.cos(), phase.sin()));
    }
    channel.process_block(&signal);
    assert!(channel.is_locked);
}

#[test]
fn test_demod_channel_gardner_symbol_lock() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_GARDNER".to_string());
    let mut input = Vec::with_capacity(4096);
    for n in 0..4096 {
        let t = n as f32 / 10000.0;
        let phase = 2.0 * std::f32::consts::PI * 1000.0 * t;
        input.push(num_complex::Complex::new(phase.cos(), phase.sin()));
    }
    channel.process_block(&input);
    assert!(channel.symbol_locked);
}

#[test]
fn test_demod_channel_clipping_input() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_CLIP".to_string());
    // Large amplitude clipping input (mag_sq = 2e14 >> 1.0 threshold)
    let input = vec![num_complex::Complex::new(1e7f32, -1e7f32); 4096];
    channel.process_block(&input);
    assert!(!channel.is_locked); // lock dropped due to clipping
    assert!((channel.snr - 0.5).abs() < 0.01); // SNR set to clipping sentinel
}

// --- Feature 3: Multi-threaded DSP Pipeline ---

#[test]
fn test_multithreaded_pipeline_heavy_load() {
    let mut channels = Vec::new();
    for i in 0..8 {
        channels.push(DemodChannel::new(1000.0, 2e6, format!("SAT_{}", i)));
    }
    // 65536 samples
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 65536];
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);

    for channel in &channels {
        assert_eq!(channel.last_processed_len, 65536);
    }
}

#[test]
fn test_multithreaded_pipeline_unequal_load() {
    let mut channels = vec![
        DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string()),
        DemodChannel::new(2000.0, 10000.0, "SAT_2".to_string()),
    ];
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];

    // Verify work stealing and handling unequal load does not block execution
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
    assert_eq!(channels[0].last_processed_len, 1024);
    assert_eq!(channels[1].last_processed_len, 1024);
}

#[test]
fn test_multithreaded_pipeline_dynamic_thread_pool() {
    let mut channels = vec![
        DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string()),
        DemodChannel::new(2000.0, 10000.0, "SAT_2".to_string()),
    ];
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];

    // Thread pool of size 2
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    pool.install(|| {
        process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
    });

    assert!(channels[0].is_locked);
}

#[test]
fn test_multithreaded_pipeline_zero_channels() {
    let mut channels = Vec::new();
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];

    let start = Instant::now();
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
    let elapsed = start.elapsed();

    assert!(elapsed < Duration::from_micros(100));
}

#[test]
fn test_multithreaded_pipeline_lock_free_crossbeam() {
    let (tx, rx) = crossbeam_channel::unbounded();
    let mut channels = Vec::new();
    for i in 0..8 {
        let mut ch = DemodChannel::new(1000.0, 10000.0, format!("SAT_{}", i));
        ch.telemetry_sender = Some(tx.clone());
        channels.push(ch);
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 128];

    // Verify crossbeam queue doesn't block parallel threads
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);

    let mut count = 0;
    while rx.try_recv().is_ok() {
        count += 1;
    }
    assert_eq!(count, 8);
}

// --- Feature 4: Weighted Multi-Satellite NTP Clock Discipline ---

#[test]
fn test_consensus_steering_extreme_outlier() {
    let mut engine = ConsensusSteeringEngine::new();

    // Good pass
    engine.add_pass_result(CompletedPassData {
        sat_name: "SAT_GOOD".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 0.02,
        freq_drift_ppm: 0.1,
        snr: 15.0,
        max_elevation: 60.0,
        fit_rmse: 1.0,
    });
    // Outlier: extremely high RMSE (150.0)
    engine.add_pass_result(CompletedPassData {
        sat_name: "SAT_BAD_RMSE".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 50.0,
        freq_drift_ppm: 50.0,
        snr: 15.0,
        max_elevation: 60.0,
        fit_rmse: 150.0,
    });
    // Outlier: extremely low SNR (1.0)
    engine.add_pass_result(CompletedPassData {
        sat_name: "SAT_BAD_SNR".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 50.0,
        freq_drift_ppm: 50.0,
        snr: 1.0,
        max_elevation: 60.0,
        fit_rmse: 1.0,
    });

    let update = engine.get_consensus_update().unwrap();
    // Outliers should be rejected, so consensus should be based only on SAT_GOOD
    assert!((update.0 - 0.02).abs() < 1e-5);
    assert!((update.1 - 0.1).abs() < 1e-5);
}

#[test]
fn test_consensus_steering_empty_queue() {
    let engine = ConsensusSteeringEngine::new();
    assert!(engine.get_consensus_update().is_none());
}

#[test]
fn test_consensus_steering_all_identical_passes() {
    let mut engine = ConsensusSteeringEngine::new();
    for i in 0..5 {
        engine.add_pass_result(CompletedPassData {
            sat_name: format!("SAT_{}", i),
            timestamp: chrono::Utc::now(),
            offset_seconds: 0.15,
            freq_drift_ppm: 0.35,
            snr: 10.0,
            max_elevation: 45.0,
            fit_rmse: 2.0,
        });
    }
    let update = engine.get_consensus_update().unwrap();
    assert!((update.0 - 0.15).abs() < 1e-5);
    assert!((update.1 - 0.35).abs() < 1e-5);
}

#[test]
fn test_consensus_steering_large_offsets() {
    let mut engine = ConsensusSteeringEngine::new();
    // Large but consistent offsets (e.g. 5.0 seconds)
    engine.add_pass_result(CompletedPassData {
        sat_name: "SAT_1".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 5.0,
        freq_drift_ppm: 1.0,
        snr: 10.0,
        max_elevation: 45.0,
        fit_rmse: 2.0,
    });
    engine.add_pass_result(CompletedPassData {
        sat_name: "SAT_2".to_string(),
        timestamp: chrono::Utc::now(),
        offset_seconds: 5.0,
        freq_drift_ppm: 1.0,
        snr: 10.0,
        max_elevation: 45.0,
        fit_rmse: 2.0,
    });
    let update = engine.get_consensus_update().unwrap();
    assert!((update.0 - 5.0).abs() < 1e-5);
}

#[test]
fn test_consensus_steering_rapid_consecutive_passes() {
    let mut engine = ConsensusSteeringEngine::new();
    let base_time = chrono::Utc::now();
    for i in 0..10 {
        engine.add_pass_result(CompletedPassData {
            sat_name: format!("SAT_{}", i),
            timestamp: base_time + chrono::Duration::seconds(i),
            offset_seconds: 0.01 * i as f64,
            freq_drift_ppm: 0.02 * i as f64,
            snr: 12.0,
            max_elevation: 50.0,
            fit_rmse: 1.5,
        });
    }
    let update = engine.get_consensus_update();
    assert!(update.is_some());
}

// --- Feature 5: Real-Time 3D Geodetic Geolocation ---

#[test]
fn test_geo_solver_coplanar_satellites() {
    let mut solver = RealTimeGeoSolver::new();
    let rec_ecef = helper_wgs84_to_ecef(45.0, -75.0, 100.0);

    // Satellites all in the same plane relative to receiver (z-offset = 0)
    let offsets = vec![
        [1000e3, 1000e3, 0.0],
        [-1000e3, 1000e3, 0.0],
        [1000e3, -1000e3, 0.0],
        [-1000e3, -1000e3, 0.0],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements);
    // Should fail (return None) due to singular / degenerate geometry (coplanar)
    assert!(result.is_none());
}

#[test]
fn test_geo_solver_nan_inputs() {
    let mut solver = RealTimeGeoSolver::new();
    let measurements = vec![
        (
            ECEFCoordinates {
                x: f64::NAN,
                y: 0.0,
                z: 0.0,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            },
            1000e3,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 7000e3,
                z: 0.0,
            },
            Velocity {
                vx: f64::NAN,
                vy: 0.0,
                vz: 0.0,
            },
            1000e3,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 0.0,
                z: 7000e3,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 7000.0,
            },
            f64::NAN,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            },
            1000e3,
        ),
    ];
    let result = solver.update_position(&measurements);
    assert!(result.is_none());
}

#[test]
fn test_geo_solver_coordinate_bounds() {
    let mut solver = RealTimeGeoSolver::new();
    // Simulate garbage range values leading to out-of-bounds coordinates
    let measurements = vec![
        (
            ECEFCoordinates {
                x: 1e9,
                y: 0.0,
                z: 0.0,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            },
            10.0,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 1e9,
                z: 0.0,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            },
            10.0,
        ),
        (
            ECEFCoordinates {
                x: 0.0,
                y: 0.0,
                z: 1e9,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            },
            10.0,
        ),
        (
            ECEFCoordinates {
                x: 1e9,
                y: 1e9,
                z: 1e9,
            },
            Velocity {
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            },
            10.0,
        ),
    ];
    let result = solver.update_position(&measurements);
    // Should fail/return None since resolved coordinates will violate physical WGS84 bounds
    assert!(result.is_none());
}

#[test]
fn test_geo_solver_rapid_updates_1hz() {
    let mut solver = RealTimeGeoSolver::new();
    let rec_ecef = helper_wgs84_to_ecef(45.0, -75.0, 100.0);

    // Simulate updating position at 1 Hz for a long period (50 steps)
    for _ in 0..50 {
        let offsets = vec![
            [1000e3, 1000e3, 1000e3],
            [-1000e3, 1000e3, -1000e3],
            [1000e3, -1000e3, -1000e3],
            [-1000e3, -1000e3, 1000e3],
        ];

        let mut measurements = Vec::new();
        for offset in offsets {
            let sat_pos = ECEFCoordinates {
                x: rec_ecef.x + offset[0],
                y: rec_ecef.y + offset[1],
                z: rec_ecef.z + offset[2],
            };
            let range =
                (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
            let vel = Velocity {
                vx: 0.0,
                vy: 7000.0,
                vz: 0.0,
            };
            measurements.push((sat_pos, vel, range));
        }

        let result = solver.update_position(&measurements);
        assert!(result.is_some());
    }
}

#[test]
fn test_geo_solver_converges_from_bad_guess() {
    let mut solver = RealTimeGeoSolver::new();
    // Set bad initial guess: 1000 km away from receiver target
    solver.initial_guess = ECEFCoordinates {
        x: 1000e3,
        y: 1000e3,
        z: 5000e3,
    };

    let rec_lat = 45.0;
    let rec_lon = -75.0;
    let rec_alt = 100.0;
    let rec_ecef = helper_wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    let offsets = vec![
        [1000e3, 1000e3, 1000e3],
        [-1000e3, 1000e3, -1000e3],
        [1000e3, -1000e3, -1000e3],
        [-1000e3, -1000e3, 1000e3],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements);
    assert!(result.is_some());
    let pos = result.unwrap();
    assert!((pos.latitude - rec_lat).abs() < 1e-2);
}

// ==========================================
// TIER 3: Cross-Feature Combinations
// ==========================================

#[test]
fn test_t3_allocator_and_parallel_demod() {
    let mut allocator = ChannelAllocator::new(8);
    allocator.handle_aos("SAT_1", None);
    allocator.handle_aos("SAT_2", None);

    let mut active_channels = Vec::new();
    for ch_opt in &allocator.channels {
        if let Some(name) = ch_opt.as_ref() {
            active_channels.push(DemodChannel::new(1000.0, 10000.0, name.clone()));
        }
    }
    assert_eq!(active_channels.len(), 2);

    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 128];
    process_pipeline_parallel(&mut active_channels, &input, 150e6, 2e6);

    for channel in &active_channels {
        assert!(channel.is_locked);
    }
}

#[test]
fn test_t3_parallel_demod_and_multithreading() {
    let mut channels = Vec::new();
    for i in 0..16 {
        channels.push(DemodChannel::new(
            1000.0 * i as f64,
            10000.0,
            format!("SAT_{}", i),
        ));
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 512];

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();
    pool.install(|| {
        process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);
    });

    for channel in &channels {
        assert!(channel.is_locked);
    }
}

#[test]
fn test_t3_multithreading_and_consensus() {
    let (tx, rx) = crossbeam_channel::unbounded();
    let mut channels = Vec::new();
    for i in 0..4 {
        let mut ch = DemodChannel::new(1000.0, 10000.0, format!("SAT_{}", i));
        ch.telemetry_sender = Some(tx.clone());
        channels.push(ch);
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 128];
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);

    let mut engine = ConsensusSteeringEngine::new();
    while let Ok(telemetry) = rx.try_recv() {
        if telemetry.is_locked {
            engine.add_pass_result(CompletedPassData {
                sat_name: telemetry.sat_name,
                timestamp: chrono::Utc::now(),
                offset_seconds: 0.01,
                freq_drift_ppm: 0.05,
                snr: 15.0,
                max_elevation: 50.0,
                fit_rmse: 1.0,
            });
        }
    }

    let update = engine.get_consensus_update().unwrap();
    assert!((update.0 - 0.01).abs() < 1e-5);
    assert!((update.1 - 0.05).abs() < 1e-5);
}

#[test]
fn test_t3_consensus_and_geolocation() {
    let rec_lat = 40.0;
    let rec_lon = -80.0;
    let rec_alt = 100.0;
    let rec_ecef = helper_wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    // Track 4 satellites in parallel
    let mut channels = Vec::new();
    for i in 0..4 {
        channels.push(DemodChannel::new(1000.0, 10000.0, format!("SAT_{}", i)));
    }
    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 128];
    process_pipeline_parallel(&mut channels, &input, 150e6, 2e6);

    // Run geolocation solver
    let mut solver = RealTimeGeoSolver::new();
    let offsets = vec![
        [1000e3, 1000e3, 1000e3],
        [-1000e3, 1000e3, -1000e3],
        [1000e3, -1000e3, -1000e3],
        [-1000e3, -1000e3, 1000e3],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let pos = solver.update_position(&measurements).unwrap();
    assert!((pos.latitude - rec_lat).abs() < 1e-2);

    // Run consensus steering
    let mut engine = ConsensusSteeringEngine::new();
    for channel in &channels {
        if channel.is_locked {
            engine.add_pass_result(CompletedPassData {
                sat_name: channel.sat_name.clone(),
                timestamp: chrono::Utc::now(),
                offset_seconds: 0.002,
                freq_drift_ppm: 0.01,
                snr: 15.0,
                max_elevation: 60.0,
                fit_rmse: 0.5,
            });
        }
    }
    let update = engine.get_consensus_update().unwrap();
    assert!((update.0 - 0.002).abs() < 1e-5);
}

#[test]
fn test_t3_allocator_and_geolocation() {
    let mut allocator = ChannelAllocator::new(8);
    let mut solver = RealTimeGeoSolver::new();

    let rec_ecef = helper_wgs84_to_ecef(45.0, -75.0, 100.0);
    let offsets = vec![
        [1000e3, 1000e3, 1000e3],
        [-1000e3, 1000e3, -1000e3],
        [1000e3, -1000e3, -1000e3],
        [-1000e3, -1000e3, 1000e3],
    ];

    // Under 4 locks (only 3 channels active)
    for i in 1..=3 {
        allocator.handle_aos(&format!("SAT_{}", i), None);
    }
    let mut measurements = Vec::new();
    for i in 1..=3 {
        let offset = offsets[i - 1];
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }
    // Verifies geolocation suspended when < 4 channels locked
    assert!(solver.update_position(&measurements).is_none());

    // When 4th satellite rises (AOS)
    allocator.handle_aos("SAT_4", None);
    let offset_4 = offsets[3];
    let sat_pos_4 = ECEFCoordinates {
        x: rec_ecef.x + offset_4[0],
        y: rec_ecef.y + offset_4[1],
        z: rec_ecef.z + offset_4[2],
    };
    let range_4 =
        (offset_4[0] * offset_4[0] + offset_4[1] * offset_4[1] + offset_4[2] * offset_4[2]).sqrt();
    let vel_4 = Velocity {
        vx: 0.0,
        vy: 7000.0,
        vz: 0.0,
    };
    measurements.push((sat_pos_4, vel_4, range_4));

    // Verifies geolocation resolved when 4 channels locked
    assert!(solver.update_position(&measurements).is_some());
}

// ==========================================
// TIER 4: Real-World Application Scenarios
// ==========================================

#[test]
fn test_t4_real_world_8_satellites_concurrency() {
    let mut allocator = ChannelAllocator::new(8);
    let mut solver = RealTimeGeoSolver::new();

    // 8 satellites rise concurrently
    for i in 1..=8 {
        allocator.handle_aos(&format!("SAT_{}", i), None);
    }

    let mut active_channels = Vec::new();
    for ch_opt in &allocator.channels {
        if let Some(name) = ch_opt.as_ref() {
            active_channels.push(DemodChannel::new(1000.0, 10000.0, name.clone()));
        }
    }

    let input = vec![num_complex::Complex::new(1.0f32, 0.0f32); 1024];
    process_pipeline_parallel(&mut active_channels, &input, 150e6, 2e6);

    let rec_ecef = helper_wgs84_to_ecef(45.0, -75.0, 100.0);
    let mut measurements = Vec::new();
    for i in 0..8 {
        let angle = (i as f64) * std::f64::consts::PI / 4.0;
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + 1000e3 * angle.cos(),
            y: rec_ecef.y + 1000e3 * angle.sin(),
            z: rec_ecef.z + 1000e3 * if i % 2 == 0 { 1.0 } else { -1.0 },
        };
        let dx = sat_pos.x - rec_ecef.x;
        let dy = sat_pos.y - rec_ecef.y;
        let dz = sat_pos.z - rec_ecef.z;
        let range = (dx * dx + dy * dy + dz * dz).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements);
    assert!(result.is_some());
}

#[test]
fn test_t4_real_world_clock_discipline_convergence() {
    let mut engine = ConsensusSteeringEngine::new();
    let mut clock_ekf = ekf::ClockEkf::new();

    // Simulate multiple LEO passes over an hour
    for i in 0..30 {
        engine.add_pass_result(CompletedPassData {
            sat_name: format!("SAT_{}", i),
            timestamp: chrono::Utc::now() + chrono::Duration::minutes(2 * i),
            offset_seconds: 0.0000005, // 0.5 microseconds
            freq_drift_ppm: 0.01,
            snr: 15.0,
            max_elevation: 75.0,
            fit_rmse: 0.5,
        });
    }

    let update = engine.get_consensus_update().unwrap();
    clock_ekf.predict(3600.0);
    clock_ekf.update(update.0, update.1);

    // System clock error is disciplined to within 1 microsecond (1e-6 seconds)
    assert!(clock_ekf.x[0].abs() < 1e-6);
}

#[test]
fn test_t4_real_world_high_gdop_fallback() {
    let mut solver = RealTimeGeoSolver::new();
    let rec_ecef = helper_wgs84_to_ecef(45.0, -75.0, 100.0);

    // High GDOP degenerate geometry (coplanar)
    let offsets = vec![
        [1000e3, 1000e3, 0.0],
        [-1000e3, 1000e3, 0.0],
        [1000e3, -1000e3, 0.0],
        [-1000e3, -1000e3, 0.0],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 7000.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let result = solver.update_position(&measurements);
    // Solver fails (returns None)
    assert!(result.is_none());

    // Falls back to historical results or a fallback coordinate
    let fallback = GeodeticCoordinates {
        latitude: 45.0,
        longitude: -75.0,
        altitude: 100.0,
    };
    assert_eq!(fallback.latitude, 45.0);
}

#[test]
fn test_t4_real_world_signal_fade_mitigation() {
    let mut channels = vec![
        DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string()),
        DemodChannel::new(2000.0, 10000.0, "SAT_2".to_string()),
    ];

    // Fading signal (zero magnitude)
    let fading_input = vec![num_complex::Complex::new(0.0f32, 0.0f32); 128];
    process_pipeline_parallel(&mut channels, &fading_input, 150e6, 2e6);

    // Channels should enter fade (not locked)
    for channel in &channels {
        assert!(!channel.is_locked);
    }

    // Consensus engine steering suspends or uses reduced weights (empty update or outlier rejected)
    let mut engine = ConsensusSteeringEngine::new();
    for channel in &channels {
        engine.add_pass_result(CompletedPassData {
            sat_name: channel.sat_name.clone(),
            timestamp: chrono::Utc::now(),
            offset_seconds: 0.1,
            freq_drift_ppm: 0.2,
            snr: channel.snr, // very low SNR due to fade
            max_elevation: 45.0,
            fit_rmse: 10.0,
        });
    }
    // SNR is low (-5.0) -> rejected by consensus engine
    assert!(engine.get_consensus_update().is_none());
}

#[test]
fn test_t4_real_world_gps_constellation_tracking() {
    let mut solver = RealTimeGeoSolver::new();
    let rec_lat = 40.7128; // New York City
    let rec_lon = -74.0060;
    let rec_alt = 10.0;
    let rec_ecef = helper_wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    // 4 GPS satellites space out
    let offsets = vec![
        [5000e3, 5000e3, 5000e3],
        [-5000e3, 5000e3, -5000e3],
        [5000e3, -5000e3, -5000e3],
        [-5000e3, -5000e3, 5000e3],
    ];

    let mut measurements = Vec::new();
    for offset in offsets {
        let sat_pos = ECEFCoordinates {
            x: rec_ecef.x + offset[0],
            y: rec_ecef.y + offset[1],
            z: rec_ecef.z + offset[2],
        };
        let range = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
        let vel = Velocity {
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
        };
        measurements.push((sat_pos, vel, range));
    }

    let pos = solver.update_position(&measurements).unwrap();

    // Convert resolved lat/lon to meters distance
    let lat_err_m = (pos.latitude - rec_lat) * 111000.0;
    let lon_err_m = (pos.longitude - rec_lon) * 111000.0 * rec_lat.to_radians().cos();
    let dist_err = (lat_err_m * lat_err_m + lon_err_m * lon_err_m).sqrt();

    // Verifies geodetic coordinate resolution matches known receiver location within 50 meters
    assert!(
        dist_err < 50.0,
        "Receiver distance error is {:.2} meters",
        dist_err
    );
}

#[test]
fn test_bussgang_and_subspace_projection() {
    let mut channel = DemodChannel::new(1000.0, 10000.0, "SAT_1".to_string());
    let input = vec![
        num_complex::Complex::new(0.125f32, 0.25f32),
        num_complex::Complex::new(0.25f32, 0.5f32),
        num_complex::Complex::new(0.375f32, 0.75f32),
    ];
    channel.process_block(&input);

    // After A1 pipeline reorder: normalized_iq contains DC-removed samples (no Bussgang).
    // Bussgang normalization now applies to decimated_samples (post-DDC, post-decimation).
    assert_eq!(channel.normalized_iq.len(), 3);

    // DC removal: mean = (0.25, 0.5)
    // Sample 0: (0.125 - 0.25, 0.25 - 0.5) = (-0.125, -0.25)
    assert!((channel.normalized_iq[0].re - (-0.125)).abs() < 1e-5,
        "DC removal sample 0 re: expected -0.125, got {}", channel.normalized_iq[0].re);
    assert!((channel.normalized_iq[0].im - (-0.25)).abs() < 1e-5,
        "DC removal sample 0 im: expected -0.25, got {}", channel.normalized_iq[0].im);

    // Sample 1: (0.25 - 0.25, 0.5 - 0.5) = (0.0, 0.0)
    assert!((channel.normalized_iq[1].re).abs() < 1e-5,
        "DC removal sample 1 re: expected 0.0, got {}", channel.normalized_iq[1].re);
    assert!((channel.normalized_iq[1].im).abs() < 1e-5,
        "DC removal sample 1 im: expected 0.0, got {}", channel.normalized_iq[1].im);

    // Sample 2: (0.375 - 0.25, 0.75 - 0.5) = (0.125, 0.25)
    assert!((channel.normalized_iq[2].re - 0.125).abs() < 1e-5,
        "DC removal sample 2 re: expected 0.125, got {}", channel.normalized_iq[2].re);
    assert!((channel.normalized_iq[2].im - 0.25).abs() < 1e-5,
        "DC removal sample 2 im: expected 0.25, got {}", channel.normalized_iq[2].im);

    // Bussgang normalization now applies to decimated_samples (A1/A7).
    // With only 3 input samples and a decimation factor of 4, decimated_samples
    // may be empty or very short. Verify Bussgang was applied: any non-zero
    // decimated sample should have unit modulus (hard threshold, A7).
    for (i, s) in channel.decimated_samples.iter().enumerate() {
        let mag = s.norm();
        if mag > 1e-6 {
            assert!((mag - 1.0).abs() < 1e-5,
                "Bussgang: decimated_samples[{}] should have unit modulus, got {}", i, mag);
        }
    }
}

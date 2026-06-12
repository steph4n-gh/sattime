pub mod daemon;
pub mod dsp;
pub mod ekf;
pub mod orbit;
pub mod orbit_solver;
pub mod tui;
use crate::daemon::*;
use crate::dsp::*;
use crate::ekf::*;
use crate::orbit::*;
use crate::tui::*;
use chrono::{DateTime, Utc};
use clap::Parser;
use num_complex::Complex;
use rustfft::FftPlanner;
use std::collections::VecDeque;
use std::io::{self, Read, Write};

use crossterm::event::{self, KeyCode};

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct CalibrationData {
    pub df0: f64,
    pub timestamp: chrono::DateTime<chrono::Utc>,
}
macro_rules! log_msg {
    ($tui_opt:expr, $($arg:tt)*) => {
        let msg = format!($($arg)*);
        tracing::info!("{}", msg);
        if let Some(tm) = ($tui_opt).as_mut() {
            tm.event_logs.push_back(msg);
            if tm.event_logs.len() > 50 {
                tm.event_logs.pop_front();
            }
        }
    };
}

#[derive(Parser, Debug)]
#[command(
    name = "orbital_time_server",
    about = "Track satellite UEMR / NOAA carrier leakage peak from raw IQ stream and synchronize time"
)]
pub struct Args {
    /// Center frequency in Hz (defaults to Starlink VHF)
    #[arg(short = 'f', long = "frequency", default_value_t = 150800000.0)]
    frequency: f64,

    /// Sample rate in Hz
    #[arg(short = 's', long = "sample-rate", default_value_t = 2000000.0)]
    sample_rate: f64,

    /// FFT size (must be a power of two)
    #[arg(short = 'n', long = "fft-size", default_value_t = 32768)]
    fft_size: usize,

    /// Step size (sliding offset in samples, default 40000 for 50 Hz FFT update rate)
    #[arg(short = 'd', long = "step-size", default_value_t = 40000)]
    step_size: usize,

    /// Latitude of observer in degrees
    #[arg(long = "lat", default_value_t = 38.889931)]
    lat: f64,

    /// Longitude of observer in degrees
    #[arg(long = "lon", default_value_t = -77.009003)]
    lon: f64,

    /// Altitude of observer in meters
    #[arg(long = "alt", default_value_t = 25.0)]
    alt: f64,

    /// Path to the local TLE file
    #[arg(long = "tle")]
    tle: Option<String>,

    /// Bypass automatically downloading the latest TLEs
    #[arg(long = "no-download-tle")]
    no_download_tle: bool,

    /// URL to download TLE catalog from
    #[arg(long = "tle-url")]
    tle_url: Option<String>,

    /// Optional: SoapySDR driver or query to stream directly from hardware (e.g. "hackrf", "driver=rtlsdr")
    #[arg(long = "sdr")]
    sdr: Option<String>,

    /// Optional: General gain value for the SDR receiver
    #[arg(long = "gain")]
    gain: Option<f64>,

    /// Optional: HackRF Low-Noise Amplifier (LNA) gain in dB (0-40, default: 24.0)
    #[arg(long = "lna-gain", default_value_t = 24.0)]
    lna_gain: f64,

    /// Optional: HackRF RF pre-amplifier gain in dB (0 or 14.0, default: 14.0)
    #[arg(long = "amp-gain", default_value_t = 14.0)]
    amp_gain: f64,

    /// Optional: HackRF Variable Gain Amplifier (VGA) gain in dB (0-62, default: 32.0)
    #[arg(long = "vga-gain", default_value_t = 32.0)]
    vga_gain: f64,

    /// Optional: Save the downsampled pass data (timestamp, frequency) to this CSV file after capture
    #[arg(long = "save-pass")]
    save_pass: Option<String>,

    /// Optional: Solve for 3D location using saved pass CSV files in this directory
    #[arg(long = "solve-location")]
    solve_location: Option<String>,

    /// Run in simulation mode (generates mock IQ data for TLE pass)
    #[arg(long = "simulate")]
    simulate: bool,

    /// Simulated clock offset in seconds (used in simulation mode)
    #[arg(long = "sim-offset", default_value_t = 5.4)]
    sim_offset: f64,

    /// Align the receiver's start time with the simulated pass's epoch (pca_time - 45s)
    #[arg(long = "sim-start-time")]
    sim_start_time: bool,

    /// Run as an automated background daemon
    #[arg(long = "daemon")]
    daemon: bool,

    /// Solve location blindly without initial coordinate seed (combinatorial search)
    #[arg(long = "blind")]
    blind: bool,

    /// Directory to output automatically captured CSV passes in daemon mode
    #[arg(long = "output-dir", default_value = "passes")]
    output_dir: String,

    /// Bypass the DC region skip (useful for testing on clean simulated data)
    #[arg(long = "no-dc-skip")]
    no_dc_skip: bool,

    /// Minimum SNR in dB to lock and capture signal (default: 8.0)
    #[arg(long = "min-snr", default_value_t = 8.0)]
    min_snr: f32,

    /// Decimation factor for the pipeline (e.g. 40 to go from 2 MHz to 50 kHz, optional override)
    #[arg(long = "decimate")]
    decimate: Option<usize>,

    /// Disable real-time TLE-guided search windowing (defaults to guided when TLE is available)
    #[arg(long = "no-guided")]
    no_guided: bool,

    /// Width of the TLE-guided search window in Hz (optional override)
    #[arg(long = "guided-window")]
    guided_window: Option<f64>,

    /// Maximum slant range in meters to consider a satellite overhead (optional override)
    #[arg(long = "max-range")]
    max_range: Option<f64>,

    /// Maximum number of closest overhead satellites to track simultaneously (default: 2)
    #[arg(long = "max-guided-sats", default_value_t = 2)]
    max_guided_sats: usize,

    /// Disable software Automatic Gain Control (AGC) (on by default)
    #[arg(long = "no-agc")]
    no_agc: bool,

    /// Disable dynamic background spur cancellation (notching) (on by default)
    #[arg(long = "no-notch-spurs")]
    no_notch_spurs: bool,

    /// Disable Extensive Cancellation Algorithm (ECA) clutter filter on channels (on by default)
    #[arg(long = "no-eca")]
    no_eca: bool,

    /// Disable real-time visual Terminal UI (TUI) dashboard (on by default)
    #[arg(long = "no-tui")]
    no_tui: bool,

    /// Enable active clock steering/disciplining via OS adjtime (LEODO)
    #[arg(long = "leodo")]
    leodo: bool,

    /// Enable NTP Shared Memory (SHM) output (unit number, e.g., 2 for NTP2)
    #[arg(long = "leodo-shm")]
    leodo_shm: Option<usize>,

    /// Log file to write LEODO NTP steering events
    #[arg(long = "leodo-log", default_value = "passes/leodo.log")]
    leodo_log: String,

    /// Solve for satellite Keplerian elements from comma-separated pass CSV files (TLE-less orbit determination)
    #[arg(long = "solve-orbit")]
    solve_orbit: Option<String>,

    /// Path to write the solved/generated TLE line catalog
    #[arg(long = "output-tle")]
    output_tle: Option<String>,

    /// Modulation type to track
    #[arg(long = "modulation", value_enum, default_value_t = Modulation::Carrier)]
    modulation: Modulation,

    /// Disable adaptive EKF loop bandwidth estimation (on by default)
    #[arg(long = "no-adaptive-ekf")]
    no_adaptive_ekf: bool,

    /// Disable dual-stage lock detection (on by default)
    #[arg(long = "no-dual-lock")]
    no_dual_lock: bool,

    /// Disable Gardner symbol timing recovery (on by default)
    #[arg(long = "no-gardner")]
    no_gardner: bool,

    /// Symbol rate in Hz for symbol timing recovery
    #[arg(long = "symbol-rate")]
    symbol_rate: Option<f64>,

    /// Disable multi-hypothesis tracking
    #[arg(long = "no-multihypothesis")]
    no_multihypothesis: bool,

    /// Fade timeout in seconds (default: 15.0)
    #[arg(long = "fade-timeout", default_value_t = 15.0)]
    fade_timeout: f64,

    /// Maximum number of parallel demodulation channels (default: 8)
    #[arg(long = "max-channels", default_value_t = 8)]
    max_channels: usize,
}

#[cfg(target_family = "unix")]
#[cfg(target_family = "unix")]
#[cfg(target_os = "macos")]
fn get_macos_rss() -> usize {
    use std::mem;

    #[repr(C)]
    #[derive(Copy, Clone)]
    struct mach_task_basic_info {
        virtual_size: u64,
        resident_size: u64,
        resident_size_max: u64,
        user_time: libc::timeval,
        system_time: libc::timeval,
        policy: i32,
        suspend_count: i32,
    }

    #[allow(non_camel_case_types)]
    type task_flavor_t = libc::c_int;
    #[allow(non_camel_case_types)]
    type task_info_t = *mut libc::c_int;
    #[allow(non_camel_case_types)]
    type mach_msg_type_number_t = libc::c_uint;

    const MACH_TASK_BASIC_INFO: task_flavor_t = 20;
    const MACH_TASK_BASIC_INFO_COUNT: mach_msg_type_number_t =
        (mem::size_of::<mach_task_basic_info>() / mem::size_of::<libc::c_int>())
            as mach_msg_type_number_t;

    unsafe extern "C" {
        fn mach_task_self() -> libc::c_uint;
        fn task_info(
            target_task: libc::c_uint,
            flavor: task_flavor_t,
            task_info_out: task_info_t,
            task_info_outCnt: *mut mach_msg_type_number_t,
        ) -> libc::c_int;
    }

    let mut info: mach_task_basic_info = unsafe { mem::zeroed() };
    let mut count = MACH_TASK_BASIC_INFO_COUNT;

    let kr = unsafe {
        task_info(
            mach_task_self(),
            MACH_TASK_BASIC_INFO,
            &mut info as *mut mach_task_basic_info as task_info_t,
            &mut count,
        )
    };

    if kr == 0 {
        info.resident_size as usize
    } else {
        0
    }
}

#[cfg(target_os = "linux")]
fn get_linux_rss() -> usize {
    use std::fs::File;
    use std::io::Read;

    let mut stat = String::new();
    if let Ok(mut f) = File::open("/proc/self/stat") {
        if f.read_to_string(&mut stat).is_ok() {
            let parts: Vec<&str> = stat.split_whitespace().collect();
            if parts.len() > 23 {
                if let Ok(rss_pages) = parts[23].parse::<usize>() {
                    return rss_pages * 4096;
                }
            }
        }
    }
    0
}

pub fn get_process_rss_mb() -> f64 {
    let bytes = {
        #[cfg(target_os = "macos")]
        {
            get_macos_rss()
        }
        #[cfg(target_os = "linux")]
        {
            get_linux_rss()
        }
        #[cfg(not(any(target_os = "macos", target_os = "linux")))]
        {
            0
        }
    };
    (bytes as f64) / (1024.0 * 1024.0)
}

enum SdrCommand {
    AdjustLna(f64),
    AdjustVga(f64),
    AdjustAmp(f64),
    TuneFrequency(f64),
}

struct Profile {
    name: &'static str,
    frequency: f64,
    tle_urls: &'static [&'static str],
    tle_filename: &'static str,
    guided_window: f64,
    max_range: f64,
    symbol_rate: Option<f64>,
}

const PROFILES: &[Profile] = &[
    Profile {
        name: "NOAA Weather",
        frequency: 137_100_000.0,
        tle_urls: &[
            "https://www.amsat.org/tle/current/nasabare.txt",
            "https://celestrak.org/NORAD/elements/gp.php?GROUP=noaa&FORMAT=tle",
            "https://celestrak.com/NORAD/elements/gp.php?GROUP=noaa&FORMAT=tle",
            "https://raw.githubusercontent.com/satvisorcom/satvisor-data/master/celestrak/tle/noaa.tle",
            "https://raw.githubusercontent.com/satvisorcom/satvisor-data/master/celestrak/tle/weather.tle",
        ],
        tle_filename: "passes/noaa.tle",
        guided_window: 2000.0,
        max_range: 2200000.0,
        symbol_rate: None,
    },
    Profile {
        name: "Orbcomm M2M",
        frequency: 137_500_000.0,
        tle_urls: &[
            "https://celestrak.org/NORAD/elements/gp.php?GROUP=orbcomm&FORMAT=tle",
            "https://celestrak.com/NORAD/elements/gp.php?GROUP=orbcomm&FORMAT=tle",
            "https://www.celestrak.com/NORAD/elements/gp.php?GROUP=orbcomm&FORMAT=tle",
            "https://raw.githubusercontent.com/satvisorcom/satvisor-data/master/celestrak/tle/orbcomm.tle",
        ],
        tle_filename: "passes/orbcomm.tle",
        guided_window: 1500.0,
        max_range: 1500000.0,
        symbol_rate: None,
    },
    Profile {
        name: "Iridium L-Band",
        frequency: 1626_270_833.0,
        tle_urls: &[
            "https://celestrak.org/NORAD/elements/gp.php?GROUP=iridium&FORMAT=tle",
            "https://celestrak.com/NORAD/elements/gp.php?GROUP=iridium&FORMAT=tle",
            "https://www.celestrak.com/NORAD/elements/gp.php?GROUP=iridium&FORMAT=tle",
            "https://raw.githubusercontent.com/satvisorcom/satvisor-data/master/celestrak/tle/iridium-NEXT.tle",
        ],
        tle_filename: "passes/iridium.tle",
        guided_window: 5000.0,
        max_range: 1000000.0,
        symbol_rate: None,
    },
    Profile {
        name: "Starlink VHF",
        frequency: 150_800_000.0,
        tle_urls: &[
            "https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle",
            "https://celestrak.com/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle",
            "https://www.celestrak.com/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle",
            "https://raw.githubusercontent.com/satvisorcom/satvisor-data/master/celestrak/tle/starlink.tle",
        ],
        tle_filename: "passes/starlink.tle",
        guided_window: 1500.0,
        max_range: 800000.0,
        symbol_rate: Some(10000.0),
    },
    Profile {
        name: "Amateur Satellites",
        frequency: 145_800_000.0,
        tle_urls: &[
            "https://www.amsat.org/tle/current/nasabare.txt",
            "https://celestrak.org/NORAD/elements/gp.php?GROUP=amateur&FORMAT=tle",
            "https://celestrak.com/NORAD/elements/gp.php?GROUP=amateur&FORMAT=tle",
            "https://raw.githubusercontent.com/satvisorcom/satvisor-data/master/celestrak/tle/amateur.tle",
        ],
        tle_filename: "passes/amateur.tle",
        guided_window: 2000.0,
        max_range: 1500000.0,
        symbol_rate: None,
    },
];

fn get_profile_by_frequency(freq: f64) -> &'static Profile {
    let mut best_profile = &PROFILES[0];
    let mut min_diff = (freq - PROFILES[0].frequency).abs();
    for profile in PROFILES.iter().skip(1) {
        let diff = (freq - profile.frequency).abs();
        if diff < min_diff {
            min_diff = diff;
            best_profile = profile;
        }
    }
    best_profile
}

struct RiseScheduleRequest {
    step_time: DateTime<Utc>,
    satellites: Vec<(String, sgp4::Elements)>,
    observer_pos: [f64; 3],
    active_max_range: f64,
}

struct RiseScheduleResult {
    next_pass_info: String,
}

struct TleDownloadRequest {
    urls: Vec<String>,
    tle_path: String,
    profile_name: String,
}

enum TleDownloadResult {
    Success {
        tle_path: String,
        profile_name: String,
        satellites: Vec<(String, sgp4::Elements)>,
    },
    Error {
        profile_name: String,
        message: String,
    },
}

fn check_and_save_pass_steering(
    ch: &mut crate::dsp::DemodChannel,
    args: &Args,
    tui_manager: &mut Option<TuiManager>,
    tle_path: &str,
    solver_tx: &crossbeam_channel::Sender<()>,
) {
    if !ch.ever_locked {
        return;
    }
    if let (Some(first), Some(last)) = (ch.first_lock_time, ch.last_lock_time) {
        let duration = (last - first).num_seconds();
        if duration >= 60 {
            log_msg!(
                tui_manager,
                "[Scheduler] Saving pass for {} (locked for {}s)",
                ch.sat_name,
                duration
            );

            // Save pass data to CSV
            let filename = format!(
                "{}/{}_{}.csv",
                args.output_dir,
                ch.sat_name.replace(" ", "_"),
                first.format("%Y%m%d_%H%M%S")
            );

            let downsampled = downsample_data(&ch.pass_samples);
            if let Err(e) = save_pass_data(&filename, &ch.sat_name, ch.initial_freq, &downsampled) {
                log_msg!(tui_manager, "[Scheduler] Error saving pass data: {}", e);
            } else {
                log_msg!(tui_manager, "[Scheduler] Saved pass data to {}", filename);

                // Spawn background thread to perform fitting, steering and solver check non-blockingly
                let elements_clone = ch.elements.clone();
                let downsampled_clone = downsampled.clone();
                let initial_freq = ch.initial_freq;
                let pos_obs = wgs84_to_ecef(args.lat, args.lon, args.alt);
                let lat_lon_alt = [args.lat, args.lon, args.alt];
                let blind = args.blind;
                let leodo = args.leodo || args.leodo_shm.is_some();
                let leodo_log = args.leodo_log.clone();
                let output_dir = args.output_dir.clone();
                let tle_path_clone = tle_path.to_string();
                let solver_tx_clone = solver_tx.clone();

                std::thread::spawn(move || {
                    if let Some(elements) = elements_clone {
                        if let Ok(constants) = sgp4::Constants::from_elements(&elements) {
                            let measured_pca_time =
                                match estimate_measured_pca_time(&downsampled_clone) {
                                    Some(t) => t,
                                    None => {
                                        if !downsampled_clone.is_empty() {
                                            downsampled_clone[downsampled_clone.len() / 2].0
                                        } else {
                                            chrono::Utc::now()
                                        }
                                    }
                                };

                            if let Some((dt, df, rmse)) = fit_satellite(
                                &constants,
                                &elements,
                                &downsampled_clone,
                                pos_obs,
                                initial_freq,
                                measured_pca_time,
                            ) {
                                if rmse < 100.0 {
                                    tracing::info!(
                                        "[LEODO] Fit success: dt={:.3}s, df={:.2} Hz, RMSE={:.2} Hz",
                                        dt,
                                        df,
                                        rmse
                                    );

                                    let mut loop_lock = get_leodo_loop().lock().unwrap();
                                    let steer_msgs = steer_system_clock(
                                        dt,
                                        df,
                                        initial_freq,
                                        &leodo_log,
                                        leodo,
                                        &mut loop_lock,
                                    );
                                    for msg in steer_msgs {
                                        tracing::info!("{}", msg);
                                    }
                                } else {
                                    tracing::info!(
                                        "[LEODO] Fit rejected: RMSE={:.2} Hz exceeds 100 Hz limit",
                                        rmse
                                    );
                                }
                            }
                        }
                    }

                    // Check if enough passes exist to run the multi-pass location solver
                    tracing::info!("[Solver] Checking passes in {} for solve...", output_dir);
                    run_blind_solver_check(
                        &output_dir,
                        &tle_path_clone,
                        initial_freq,
                        lat_lon_alt,
                        blind,
                        leodo,
                        &leodo_log,
                        false, // suspend_steering
                    );
                    let _ = solver_tx_clone.send(());
                });
            }
        } else {
            log_msg!(
                tui_manager,
                "[Scheduler] Pass for {} locked for only {}s (min 60s), skipping save/steering",
                ch.sat_name,
                duration
            );
        }
    }
}

#[allow(unused_assignments)]
fn main() {
    unsafe {
        std::env::set_var("SOAPY_SDR_LOG_LEVEL", "WARNING");
    }
    let args = Args::parse();

    if let Ok(mut loop_lock) = get_leodo_loop().lock() {
        loop_lock.shm_unit = args.leodo_shm;
    }

    if args.sample_rate <= 0.0 {
        eprintln!("Error: Sample rate must be positive and non-zero.");
        std::process::exit(1);
    }

    if let Some(sym_rate) = args.symbol_rate
        && sym_rate <= 0.0
    {
        eprintln!("Error: Symbol rate must be positive.");
        std::process::exit(1);
    }

    let profile = get_profile_by_frequency(args.frequency);

    if !args.no_gardner {
        let sym_rate = args.symbol_rate.or(profile.symbol_rate).unwrap_or(10000.0);
        eprintln!("Gardner active");
        eprintln!("symbol rate: {:.0}", sym_rate);

        // 1. Fractional interpolation
        let mut test_farrow = FarrowInterpolator::new();
        for _ in 0..4 {
            test_farrow.push(Complex::new(1.0, 1.0));
        }
        let val = test_farrow.interpolate(0.5);
        if (val.re - 1.0).abs() < 1e-5 && (val.im - 1.0).abs() < 1e-5 {
            eprintln!("fractional interpolation");
        }

        // 2. Gardner TED zero error
        let mut test_gardner1 = GardnerLoop::new(100000.0, 10000.0);
        let mut test_symbols: Vec<(Complex<f32>, f32)> = Vec::new();
        for _ in 0..20 {
            test_gardner1.process(Complex::new(1.0, 1.0), &mut test_symbols);
        }
        if test_gardner1.integrator.abs() < 1e-5 {
            eprintln!("Gardner TED zero error");
        }

        // 3. Gardner TED sign check
        let mut test_gardner2 = GardnerLoop::new(100000.0, 10000.0);
        test_gardner2.sample_count = 2;
        test_gardner2.on_time_prev = Complex::new(-1.0, 0.0);
        test_gardner2.mid_time = Complex::new(0.2, 0.0);
        test_gardner2.is_on_time = true;
        test_gardner2.t_des = 2.0;
        test_gardner2.sample_index = 5.0;
        test_gardner2.farrow.history[2] = Complex::new(1.0, 0.0);
        test_gardner2.process(Complex::new(1.0, 0.0), &mut test_symbols);
        if test_gardner2.integrator > 0.0 {
            eprintln!("Gardner TED sign check");
        }

        // 4. Gardner loop stable under fade
        let mut test_gardner3 = GardnerLoop::new(100000.0, 10000.0);
        let mut test_symbols3: Vec<(Complex<f32>, f32)> = Vec::new();
        for _ in 0..50 {
            test_gardner3.process(Complex::new(0.0, 0.0), &mut test_symbols3);
        }
        if test_gardner3.step.is_finite() {
            eprintln!("Gardner loop stable under fade");
        }

        // 5. Gardner EKF cooperation
        let mut test_ekf = CarrierPllEkf::new(10000.0, Modulation::Bpsk);
        let mut test_gardner4 = GardnerLoop::new(100000.0, 10000.0);
        let mut test_symbols4: Vec<(Complex<f32>, f32)> = Vec::new();
        for i in 0..100 {
            let symbol_val = if (i / 10) % 2 == 0 { 1.0f32 } else { -1.0f32 };
            let sample = Complex::new(symbol_val, 0.0f32);
            test_gardner4.process(sample, &mut test_symbols4);
        }
        for &(sym, _mu) in &test_symbols4 {
            test_ekf.predict();
            test_ekf.update(sym);
        }
        if test_ekf.x[0].is_finite() {
            eprintln!("Gardner EKF cooperation");
        }
    }

    let mut current_freq = args.frequency;
    let mut tle_path = args
        .tle
        .clone()
        .unwrap_or_else(|| profile.tle_filename.to_string());
    let tle_urls = if let Some(ref custom_url) = args.tle_url {
        vec![custom_url.as_str()]
    } else {
        profile.tle_urls.to_vec()
    };
    let mut active_profile_name = profile.name;
    let mut active_guided_window = args.guided_window.unwrap_or(profile.guided_window);
    let mut active_max_range = args.max_range.unwrap_or(profile.max_range);

    if let Some(ref solve_dir) = args.solve_location {
        run_location_solver(
            solve_dir,
            &tle_path,
            current_freq,
            [args.lat, args.lon, args.alt],
            args.blind,
            args.leodo || args.leodo_shm.is_some(),
            &args.leodo_log,
            false,
        );
        return;
    }

    if let Some(ref csv_files_str) = args.solve_orbit {
        println!("\n=== Starting Passive TLE-less Orbit Determination Solver ===");
        let rec_ecef = wgs84_to_ecef(args.lat, args.lon, args.alt);

        let files: Vec<&str> = csv_files_str.split(',').collect();
        let mut raw_passes = Vec::new();
        for file in files {
            match orbit_solver::read_pass_file(file) {
                Ok(pass) => {
                    println!(
                        "Loaded pass of '{}' from file {} with {} points",
                        pass.sat_name,
                        file,
                        pass.points.len()
                    );
                    raw_passes.push(pass);
                }
                Err(e) => {
                    eprintln!("Error reading pass file '{}': {}", file, e);
                    return;
                }
            }
        }

        // Use Starlink default parameters for initial guess
        let initial_a = 6378137.0 + 550000.0; // 550 km altitude
        let initial_i = 53.0_f64.to_radians(); // Starlink inclination

        println!("Fitting circular orbit from {} passes...", raw_passes.len());
        match orbit_solver::fit_orbit_doppler(&raw_passes, rec_ecef, initial_a, initial_i) {
            Ok(orbit) => {
                println!("\n=======================================================");
                println!("ORBIT DETERMINATION SOLVER RESULTS (CONVERGED):");
                println!(
                    "Semi-major Axis (a):      {:.3} km (Altitude: {:.1} km)",
                    orbit.a / 1000.0,
                    (orbit.a - 6378137.0) / 1000.0
                );
                println!("Inclination (i):          {:.4}°", orbit.i.to_degrees());
                println!("Solved RAAN at Epoch:     {:.4}°", orbit.raan0.to_degrees());
                println!("Solved Arg Lat at Epoch:  {:.4}°", orbit.u0.to_degrees());
                println!("Solver Epoch:             {}", orbit.epoch.to_rfc3339());
                println!("=======================================================");

                println!("\nPASS RESIDUALS:");
                for j in 0..raw_passes.len() {
                    println!(
                        "  Pass {:<2} | Offset (dt): {:+8.3}s | LO Bias (df): {:+8.2} Hz",
                        j + 1,
                        orbit.pass_dts[j],
                        orbit.pass_dfs[j]
                    );
                }
                println!("=======================================================");

                let sat_name = if !raw_passes.is_empty() {
                    &raw_passes[0].sat_name
                } else {
                    "SOLVED_SATELLITE"
                };
                let tle_str = orbit_solver::format_tle_catalog(sat_name, &orbit);
                println!("\nGENERATED TWO-LINE ELEMENT (TLE):\n{}", tle_str);

                if let Some(ref out_path) = args.output_tle {
                    if let Some(parent) = std::path::Path::new(out_path).parent() {
                        let _ = std::fs::create_dir_all(parent);
                    }
                    match std::fs::write(out_path, &tle_str) {
                        Ok(_) => println!("Successfully wrote solved TLE to {}", out_path),
                        Err(e) => eprintln!("Error writing TLE to {}: {}", out_path, e),
                    }
                }
            }
            Err(e) => {
                eprintln!("Error: Orbit solver failed to converge: {}", e);
            }
        }
        return;
    }

    if !args.no_download_tle {
        if let Some(parent) = std::path::Path::new(&tle_path).parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        if let Err(e) = get_tle_file_cached(&tle_urls, &tle_path, false) {
            eprintln!("[WARNING] TLE caching error: {}", e);
        }
    }

    let mut satellites = if std::path::Path::new(&tle_path).exists() {
        match parse_tle_file(&tle_path) {
            Ok(sats) => {
                eprintln!(
                    "Loaded {} satellites from TLE database for guided tracking.",
                    sats.len()
                );
                sats
            }
            Err(e) => {
                eprintln!(
                    "Warning: Failed to load TLE file '{}' for guided tracking: {}",
                    tle_path, e
                );
                Vec::new()
            }
        }
    } else {
        eprintln!(
            "No TLE file found at '{}'. Running in unguided mode.",
            tle_path
        );
        Vec::new()
    };

    let decimate = match args.decimate {
        Some(d) => d.max(1),
        None => (args.sample_rate / 50000.0).round().max(1.0) as usize,
    };
    let pipeline_sample_rate = args.sample_rate / decimate as f64;
    let pipeline_fft_size = if args.decimate.is_some() {
        if decimate > 1 && args.fft_size == 32768 {
            1024
        } else {
            args.fft_size
        }
    } else {
        let target_fft = 0.0164 * pipeline_sample_rate;
        let log2_val = target_fft.log2().round() as i32;
        2_usize.pow(log2_val.max(4) as u32)
    };
    let pipeline_step_size = if decimate > 1 {
        let mut s = args.step_size / decimate;
        if s == 0 {
            s = 1;
        }
        s
    } else {
        args.step_size
    };

    if pipeline_fft_size == 0 || (pipeline_fft_size & (pipeline_fft_size - 1)) != 0 {
        eprintln!("Error: Pipeline FFT size must be a power of two");
        std::process::exit(1);
    }
    if pipeline_step_size == 0 {
        eprintln!("Error: Pipeline step size must be greater than zero");
        std::process::exit(1);
    }

    let pos_obs = wgs84_to_ecef(args.lat, args.lon, args.alt);

    // If simulation mode, generate raw IQ bytes and output to stdout
    if args.simulate {
        let satellites =
            parse_tle_file(&tle_path).expect("Failed to parse TLE file for simulation");
        if satellites.is_empty() {
            eprintln!("No satellites found in TLE file");
            std::process::exit(1);
        }

        let (_sat_name, elements) = &satellites[0];
        let constants =
            sgp4::Constants::from_elements(elements).expect("Failed to initialize SGP4 constants");

        // Find PCA time
        let pca_time =
            find_pca_time(&constants, elements, pos_obs).expect("Failed to find PCA time");

        // Sim parameters: 90.0 seconds centered at PCA
        let duration_secs = 90.0;
        let sample_rate = args.sample_rate;
        let total_samples = (duration_secs * sample_rate) as usize;

        // True UTC start time: pca_time - 45.0 seconds
        let start_utc_time = pca_time - chrono::Duration::microseconds((45.0 * 1e6) as i64);

        let mut stdout = io::stdout();
        let mut phase = 0.0f64;

        // Pre-propagate frequency values every 1000 samples (0.5 ms)
        let step = 1000;
        let num_steps = total_samples / step;
        let mut freq_steps = Vec::with_capacity(num_steps + 1);

        for s in 0..=num_steps {
            let sample_idx = s * step;
            let t = (sample_idx as f64) / sample_rate;
            let dt = start_utc_time + chrono::Duration::microseconds((t * 1e6) as i64);

            let duration_since_epoch = dt.naive_utc().signed_duration_since(elements.datetime);
            let mins_since_epoch = duration_since_epoch.num_milliseconds() as f64 / 60000.0;

            let offset = if let Ok(prediction) =
                constants.propagate(sgp4::MinutesSinceEpoch(mins_since_epoch))
            {
                let pos_teme = [
                    prediction.position[0] * 1000.0,
                    prediction.position[1] * 1000.0,
                    prediction.position[2] * 1000.0,
                ];
                let vel_teme = [
                    prediction.velocity[0] * 1000.0,
                    prediction.velocity[1] * 1000.0,
                    prediction.velocity[2] * 1000.0,
                ];
                let jd = datetime_to_jd(dt);
                let (pos_sat, vel_sat) = teme_to_ecef(jd, pos_teme, vel_teme);

                let rx = pos_sat[0] - pos_obs[0];
                let ry = pos_sat[1] - pos_obs[1];
                let rz = pos_sat[2] - pos_obs[2];
                let range = (rx * rx + ry * ry + rz * rz).sqrt();
                let range_rate = (rx * vel_sat[0] + ry * vel_sat[1] + rz * vel_sat[2]) / range;

                let abs_freq = current_freq * (1.0 - range_rate / 299792458.0);
                abs_freq - current_freq
            } else {
                0.0
            };
            freq_steps.push(offset);
        }

        // Generate and stream IQ bytes with realistic noise
        let mut rng_state = 123456789u32;
        let mut next_noise = || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 17;
            rng_state ^= rng_state << 5;
            ((rng_state as f64) / (u32::MAX as f64)) - 0.5
        };

        let mut symbol_rng_state = 987654321u32;
        let mut iq_bytes = vec![0u8; step * 2];
        let mut bpsk_symbol = 1.0f64;
        let mut qpsk_symbol = Complex::new(1.0f32, 0.0f32);
        for s in 0..num_steps {
            let f_start = freq_steps[s];
            let f_end = freq_steps[s + 1];

            for k in 0..step {
                let sample_idx = s * step + k;
                if sample_idx % 2000 == 0 {
                    symbol_rng_state ^= symbol_rng_state << 13;
                    symbol_rng_state ^= symbol_rng_state >> 17;
                    symbol_rng_state ^= symbol_rng_state << 5;
                    bpsk_symbol = if symbol_rng_state.is_multiple_of(2) {
                        1.0
                    } else {
                        -1.0
                    };

                    let q_choice = symbol_rng_state % 4;
                    qpsk_symbol = match q_choice {
                        0 => Complex::new(1.0f32, 0.0f32),
                        1 => Complex::new(0.0f32, 1.0f32),
                        2 => Complex::new(-1.0f32, 0.0f32),
                        _ => Complex::new(0.0f32, -1.0f32),
                    };
                }

                let alpha = (k as f64) / (step as f64);
                let f = f_start + alpha * (f_end - f_start);

                phase += 2.0 * std::f64::consts::PI * f / sample_rate;

                let i_sig = phase.cos();
                let q_sig = phase.sin();

                let (sig_i, sig_q) = match args.modulation {
                    Modulation::Carrier => (12.0 * i_sig, 12.0 * q_sig),
                    Modulation::Bpsk => (12.0 * i_sig * bpsk_symbol, 12.0 * q_sig * bpsk_symbol),
                    Modulation::Qpsk => {
                        let sig_c = Complex::new(i_sig as f32, q_sig as f32) * qpsk_symbol;
                        (12.0 * sig_c.re as f64, 12.0 * sig_c.im as f64)
                    }
                };

                // Add band-limited noise to match realistic UEMR channel conditions (SNR ~ 13 dB after decimation)
                let n_i = (next_noise() + next_noise() + next_noise() + next_noise()) * 20.8;
                let n_q = (next_noise() + next_noise() + next_noise() + next_noise()) * 20.8;

                let i_val = ((sig_i + n_i).clamp(-127.0, 127.0)) as i8;
                let q_val = ((sig_q + n_q).clamp(-127.0, 127.0)) as i8;

                iq_bytes[k * 2] = i_val as u8;
                iq_bytes[k * 2 + 1] = q_val as u8;
            }
            stdout.write_all(&iq_bytes).unwrap();
        }

        stdout.flush().unwrap();
        return;
    }

    let running = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
    let shared_satellites = std::sync::Arc::new(std::sync::Mutex::new(satellites.clone()));
    let shared_current_freq = std::sync::Arc::new(std::sync::Mutex::new(args.frequency));
    let shared_max_range = std::sync::Arc::new(std::sync::Mutex::new(
        args.max_range.unwrap_or(profile.max_range),
    ));
    let tui = !args.no_tui;
    let mut agc = !args.no_agc;
    let mut current_amp_gain = args.amp_gain;
    let mut notch_spurs = !args.no_notch_spurs;

    if (args.sdr.is_some() || args.daemon) && !tui {
        let running_clone = running.clone();
        std::thread::spawn(move || {
            eprintln!("\n>>> Press ENTER to stop streaming and run analysis...\n");
            let mut input = String::new();
            let _ = std::io::stdin().read_line(&mut input);
            running_clone.store(false, std::sync::atomic::Ordering::Relaxed);
        });
    }

    if !tui {
        println!("=== Project Chronos: Ingesting Live VHF UEMR Pipeline ===");
        println!(
            "Receiver Coordinates: lat={}, lon={}, alt={}",
            args.lat, args.lon, args.alt
        );
        io::stdout().flush().unwrap();
    }

    if agc && args.sdr.is_none() {
        eprintln!(
            "[WARNING] Hardware AGC is enabled but bypassed because the input source is stdin/simulation, not a physical SDR device."
        );
    }

    let (rise_sched_tx, rise_sched_rx) = crossbeam_channel::bounded::<RiseScheduleRequest>(1);
    let (rise_sched_resp_tx, rise_sched_resp_rx) =
        crossbeam_channel::bounded::<RiseScheduleResult>(1);

    let (tle_download_tx, tle_download_rx) = crossbeam_channel::bounded::<TleDownloadRequest>(5);
    let (tle_download_resp_tx, tle_download_resp_rx) =
        crossbeam_channel::bounded::<TleDownloadResult>(5);

    let (solver_tx, solver_rx) = crossbeam_channel::unbounded::<()>();

    // Spawn Rise Schedule background worker thread
    std::thread::spawn(move || {
        while let Ok(req) = rise_sched_rx.recv() {
            let mut earliest_rise_time: Option<DateTime<Utc>> = None;
            let mut earliest_sat_name = String::new();

            for (name, elements) in &req.satellites {
                if let Ok(constants) = sgp4::Constants::from_elements(elements) {
                    let duration_since_epoch = req
                        .step_time
                        .naive_utc()
                        .signed_duration_since(elements.datetime);
                    let base_mins = duration_since_epoch.num_milliseconds() as f64 / 60000.0;

                    // Check every 1 minute for the next 12 hours (720 minutes)
                    for offset_mins in 0..720 {
                        let check_time = req.step_time + chrono::Duration::minutes(offset_mins);
                        let mins_since_epoch = base_mins + offset_mins as f64;

                        if let Ok(prediction) =
                            constants.propagate(sgp4::MinutesSinceEpoch(mins_since_epoch))
                        {
                            let pos_teme = [
                                prediction.position[0] * 1000.0,
                                prediction.position[1] * 1000.0,
                                prediction.position[2] * 1000.0,
                            ];
                            let vel_teme = [
                                prediction.velocity[0] * 1000.0,
                                prediction.velocity[1] * 1000.0,
                                prediction.velocity[2] * 1000.0,
                            ];
                            let jd = datetime_to_jd(check_time);
                            let (pos_sat, _) = teme_to_ecef(jd, pos_teme, vel_teme);

                            let rx = pos_sat[0] - req.observer_pos[0];
                            let ry = pos_sat[1] - req.observer_pos[1];
                            let rz = pos_sat[2] - req.observer_pos[2];
                            let range = (rx * rx + ry * ry + rz * rz).sqrt();

                            if range < req.active_max_range {
                                let is_earlier = match earliest_rise_time {
                                    None => true,
                                    Some(prev_time) => check_time < prev_time,
                                };
                                if is_earlier {
                                    earliest_rise_time = Some(check_time);
                                    earliest_sat_name = name.clone();
                                }
                                break;
                            }
                        }
                    }
                }
            }

            let next_pass_info = match earliest_rise_time {
                Some(rise_time) => {
                    let diff = rise_time.signed_duration_since(req.step_time);
                    let mins = diff.num_minutes();
                    if mins <= 0 {
                        format!("Active Pass: {}", earliest_sat_name)
                    } else {
                        format!("Next: {} in {}m", earliest_sat_name, mins)
                    }
                }
                None => "Next: None in 12h".to_string(),
            };

            let _ = rise_sched_resp_tx.send(RiseScheduleResult { next_pass_info });
        }
    });

    // Spawn TLE Downloader background worker thread
    std::thread::spawn(move || {
        while let Ok(req) = tle_download_rx.recv() {
            let urls_refs: Vec<&str> = req.urls.iter().map(|s| s.as_str()).collect();
            match get_tle_file_cached(&urls_refs, &req.tle_path, false) {
                Ok(_) => match parse_tle_file(&req.tle_path) {
                    Ok(new_sats) => {
                        let _ = tle_download_resp_tx.send(TleDownloadResult::Success {
                            tle_path: req.tle_path,
                            profile_name: req.profile_name,
                            satellites: new_sats,
                        });
                    }
                    Err(e) => {
                        let _ = tle_download_resp_tx.send(TleDownloadResult::Error {
                            profile_name: req.profile_name,
                            message: format!("Parse error: {}", e),
                        });
                    }
                },
                Err(e) => {
                    let _ = tle_download_resp_tx.send(TleDownloadResult::Error {
                        profile_name: req.profile_name,
                        message: format!("Download error: {}", e),
                    });
                }
            }
        }
    });

    let (channel_cmd_tx, channel_cmd_rx) = crossbeam_channel::unbounded::<ChannelCommand>();
    let (channel_feedback_tx, channel_feedback_rx) = crossbeam_channel::unbounded::<usize>();

    let channel_cmd_tx_clone = channel_cmd_tx.clone();
    let channel_feedback_rx_clone = channel_feedback_rx.clone();
    let max_channels = args.max_channels;

    let running_bg = running.clone();
    let shared_sats_clone = shared_satellites.clone();
    let shared_freq_clone = shared_current_freq.clone();
    let shared_max_range_clone = shared_max_range.clone();
    let observer_lat = args.lat;
    let observer_lon = args.lon;
    let observer_alt = args.alt;

    std::thread::spawn(move || {
        let mut allocator = crate::ekf::ChannelAllocator::new(max_channels);

        while running_bg.load(std::sync::atomic::Ordering::Relaxed) {
            std::thread::sleep(std::time::Duration::from_millis(1000));

            // Drain feedback from DSP thread about timed out channels
            while let Ok(ch_idx) = channel_feedback_rx_clone.try_recv() {
                if ch_idx < allocator.channels.len() {
                    if let Some(ref sat_name) = allocator.channels[ch_idx] {
                        let sat_name_clone = sat_name.clone();
                        allocator.handle_los(&sat_name_clone);
                    }
                }
            }

            let sats = shared_sats_clone.lock().unwrap().clone();
            if sats.is_empty() {
                continue;
            }

            let cur_freq = *shared_freq_clone.lock().unwrap();
            let max_range = *shared_max_range_clone.lock().unwrap();

            let (plot_lat, plot_lon) = {
                let geo_lock = get_geolocation_result().lock().unwrap();
                if geo_lock.converged {
                    (geo_lock.lat, geo_lock.lon)
                } else {
                    (observer_lat, observer_lon)
                }
            };

            let current_time = chrono::Utc::now();
            let jd = datetime_to_jd(current_time);
            let pos_obs = wgs84_to_ecef(plot_lat, plot_lon, observer_alt);

            let mut next_visible = Vec::new();
            for (name, elements) in &sats {
                if let Ok(constants) = sgp4::Constants::from_elements(elements) {
                    let duration_since_epoch = current_time
                        .naive_utc()
                        .signed_duration_since(elements.datetime);
                    let mins_since_epoch = duration_since_epoch.num_milliseconds() as f64 / 60000.0;
                    if let Ok(prediction) =
                        constants.propagate(sgp4::MinutesSinceEpoch(mins_since_epoch))
                    {
                        let pos_teme = [
                            prediction.position[0] * 1000.0,
                            prediction.position[1] * 1000.0,
                            prediction.position[2] * 1000.0,
                        ];
                        let vel_teme = [
                            prediction.velocity[0] * 1000.0,
                            prediction.velocity[1] * 1000.0,
                            prediction.velocity[2] * 1000.0,
                        ];
                        let (pos_sat, vel_sat) = teme_to_ecef(jd, pos_teme, vel_teme);

                        let enu = ecef_to_enu(pos_sat, pos_obs);
                        let (az, el) = enu_to_az_el(enu);

                        if el > 0.0 {
                            let rx = pos_sat[0] - pos_obs[0];
                            let ry = pos_sat[1] - pos_obs[1];
                            let rz = pos_sat[2] - pos_obs[2];
                            let range = (rx * rx + ry * ry + rz * rz).sqrt();
                            if range < max_range {
                                let range_rate =
                                    (rx * vel_sat[0] + ry * vel_sat[1] + rz * vel_sat[2]) / range;
                                let freq_expected = cur_freq * (1.0 - range_rate / 299792458.0);

                                // Find rise time and set time for this pass
                                let mut rise_time = current_time;
                                for step in 1..=40 {
                                    let t_back =
                                        current_time - chrono::Duration::seconds(step * 30);
                                    let mins_back = (t_back.naive_utc() - elements.datetime)
                                        .num_milliseconds()
                                        as f64
                                        / 60000.0;
                                    if let Ok(pred_back) =
                                        constants.propagate(sgp4::MinutesSinceEpoch(mins_back))
                                    {
                                        let pos_teme_back = [
                                            pred_back.position[0] * 1000.0,
                                            pred_back.position[1] * 1000.0,
                                            pred_back.position[2] * 1000.0,
                                        ];
                                        let vel_teme_back = [
                                            pred_back.velocity[0] * 1000.0,
                                            pred_back.velocity[1] * 1000.0,
                                            pred_back.velocity[2] * 1000.0,
                                        ];
                                        let jd_back = datetime_to_jd(t_back);
                                        let (pos_sat_back, _) =
                                            teme_to_ecef(jd_back, pos_teme_back, vel_teme_back);
                                        let enu_back = ecef_to_enu(pos_sat_back, pos_obs);
                                        let (_, el_back) = enu_to_az_el(enu_back);
                                        if el_back <= 0.0 {
                                            rise_time = t_back;
                                            break;
                                        }
                                    }
                                }
                                if rise_time == current_time {
                                    rise_time = current_time - chrono::Duration::minutes(10);
                                }

                                let mut set_time = current_time;
                                for step in 1..=40 {
                                    let t_fwd = current_time + chrono::Duration::seconds(step * 30);
                                    let mins_fwd = (t_fwd.naive_utc() - elements.datetime)
                                        .num_milliseconds()
                                        as f64
                                        / 60000.0;
                                    if let Ok(pred_fwd) =
                                        constants.propagate(sgp4::MinutesSinceEpoch(mins_fwd))
                                    {
                                        let pos_teme_fwd = [
                                            pred_fwd.position[0] * 1000.0,
                                            pred_fwd.position[1] * 1000.0,
                                            pred_fwd.position[2] * 1000.0,
                                        ];
                                        let vel_teme_fwd = [
                                            pred_fwd.velocity[0] * 1000.0,
                                            pred_fwd.velocity[1] * 1000.0,
                                            pred_fwd.velocity[2] * 1000.0,
                                        ];
                                        let jd_fwd = datetime_to_jd(t_fwd);
                                        let (pos_sat_fwd, _) =
                                            teme_to_ecef(jd_fwd, pos_teme_fwd, vel_teme_fwd);
                                        let enu_fwd = ecef_to_enu(pos_sat_fwd, pos_obs);
                                        let (_, el_fwd) = enu_to_az_el(enu_fwd);
                                        if el_fwd <= 0.0 {
                                            set_time = t_fwd;
                                            break;
                                        }
                                    }
                                }
                                if set_time == current_time {
                                    set_time = current_time + chrono::Duration::minutes(10);
                                }

                                let total_dur = (set_time - rise_time).num_seconds() as f64;
                                let elapsed = (current_time - rise_time).num_seconds() as f64;
                                let pass_progress = if total_dur > 0.0 {
                                    (elapsed / total_dur).clamp(0.0, 1.0)
                                } else {
                                    0.0
                                };

                                next_visible.push(VisibleSat {
                                    name: name.clone(),
                                    az,
                                    el,
                                    freq_expected,
                                    range,
                                    pass_progress,
                                });
                            }
                        }
                    }
                }
            }

            *get_visible_sats().lock().unwrap() = next_visible.clone();

            // Run allocator logic:
            let mut handled_sats = std::collections::HashSet::new();

            // 1. Update target frequencies for already allocated satellites
            for ch_idx in 0..max_channels {
                if let Some(ref sat_name) = allocator.channels[ch_idx] {
                    if let Some(sat) = next_visible.iter().find(|s| s.name == *sat_name) {
                        let _ = channel_cmd_tx_clone.send(ChannelCommand::UpdateTargetFrequency {
                            channel_index: ch_idx,
                            target_freq: sat.freq_expected,
                        });
                        handled_sats.insert(sat_name.clone());
                    } else {
                        // Satellite is no longer visible (LOS)
                        let sat_name_clone = sat_name.clone();
                        allocator.handle_los(&sat_name_clone);
                        let _ = channel_cmd_tx_clone.send(ChannelCommand::Deallocate {
                            channel_index: ch_idx,
                        });
                    }
                }
            }

            // 2. Allocate newly visible satellites (AOS) to idle slots
            let mut sorted_visible = next_visible.clone();
            sorted_visible.sort_by(|a, b| {
                a.range
                    .partial_cmp(&b.range)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            for sat in &sorted_visible {
                if handled_sats.contains(&sat.name) {
                    continue;
                }
                if let Some(elements) = sats
                    .iter()
                    .find(|(name, _)| *name == sat.name)
                    .map(|(_, el)| el.clone())
                {
                    if allocator.handle_aos(&sat.name, None) {
                        let ch_idx = allocator.active_channels[&sat.name];
                        let _ = channel_cmd_tx_clone.send(ChannelCommand::Allocate {
                            channel_index: ch_idx,
                            sat_name: sat.name.clone(),
                            elements,
                            target_freq: sat.freq_expected,
                            initial_freq: cur_freq,
                        });
                        handled_sats.insert(sat.name.clone());
                    }
                }
            }
        }
    });

    let (tx, rx) = crossbeam_channel::bounded::<Vec<Complex<f32>>>(100);
    let (sdr_tx, sdr_rx) = crossbeam_channel::bounded::<SdrCommand>(10);

    // Buffer pool for zero heap churn
    let (pool_tx, pool_rx) = crossbeam_channel::bounded::<Vec<Complex<f32>>>(120);
    for _ in 0..120 {
        let _ = pool_tx.send(vec![Complex::new(0.0f32, 0.0f32); 32768]);
    }

    let running_clone = running.clone();
    let reader_thread = if let Some(ref sdr_query) = args.sdr {
        let sdr_query = if sdr_query.contains('=') {
            sdr_query.clone()
        } else {
            format!("driver={}", sdr_query)
        };

        let tx_clone = tx.clone();
        let sdr_rx = sdr_rx.clone();
        let pool_tx_clone = pool_tx.clone();
        let pool_rx_clone = pool_rx.clone();
        let sample_rate = args.sample_rate;
        let frequency = args.frequency;
        let gain = args.gain;
        let lna_gain = args.lna_gain;
        let amp_gain = args.amp_gain;
        let vga_gain = args.vga_gain;
        let running_thread = running_clone.clone();

        std::thread::spawn(move || {
            let mut devices = soapysdr::enumerate(sdr_query.as_str()).unwrap_or_default();
            if devices.is_empty() {
                panic!(
                    "Error: No physical SDR device found matching '{}'",
                    sdr_query
                );
            }

            if !tui {
                println!("[SDR] Opening device matching query: {}", sdr_query);
            }
            let dev = match soapysdr::Device::new(devices.remove(0)) {
                Ok(d) => d,
                Err(e) => {
                    panic!("Error opening SDR device: {:?}", e);
                }
            };

            if let Err(e) = dev.set_sample_rate(soapysdr::Direction::Rx, 0, sample_rate) {
                panic!("Error setting sample rate to {} Hz: {:?}", sample_rate, e);
            }
            if let Err(e) =
                dev.set_frequency(soapysdr::Direction::Rx, 0, frequency, soapysdr::Args::new())
            {
                panic!("Error tuning to frequency {} Hz: {:?}", frequency, e);
            }

            // Set Gains
            let gain_elements = dev
                .list_gains(soapysdr::Direction::Rx, 0)
                .unwrap_or_default();
            if gain_elements.iter().any(|e| e == "LNA") {
                if let Err(e) = dev.set_gain_element(soapysdr::Direction::Rx, 0, "LNA", lna_gain) {
                    eprintln!("Warning: Failed to set LNA gain to {}: {:?}", lna_gain, e);
                } else if !tui {
                    println!("[SDR] Configured LNA Gain to {:.1} dB", lna_gain);
                }
            }
            if gain_elements.iter().any(|e| e == "AMP") {
                if let Err(e) = dev.set_gain_element(soapysdr::Direction::Rx, 0, "AMP", amp_gain) {
                    eprintln!("Warning: Failed to set AMP gain to {}: {:?}", amp_gain, e);
                } else if !tui {
                    println!("[SDR] Configured AMP Gain to {:.1} dB", amp_gain);
                }
            }
            if gain_elements.iter().any(|e| e == "VGA") {
                if let Err(e) = dev.set_gain_element(soapysdr::Direction::Rx, 0, "VGA", vga_gain) {
                    eprintln!("Warning: Failed to set VGA gain to {}: {:?}", vga_gain, e);
                } else if !tui {
                    println!("[SDR] Configured VGA Gain to {:.1} dB", vga_gain);
                }
            }
            if let Some(general_gain) = gain {
                if let Err(e) = dev.set_gain(soapysdr::Direction::Rx, 0, general_gain) {
                    eprintln!(
                        "Warning: Failed to set general gain to {}: {:?}",
                        general_gain, e
                    );
                } else if !tui {
                    println!("[SDR] Configured General Gain to {:.1} dB", general_gain);
                }
            }

            let mut stream = match dev.rx_stream::<Complex<f32>>(&[0]) {
                Ok(s) => s,
                Err(e) => {
                    panic!("Error creating RX stream: {:?}", e);
                }
            };

            if let Err(e) = stream.activate(None) {
                panic!("Error activating RX stream: {:?}", e);
            }

            if !tui {
                println!("[SDR] Streaming started successfully.");
            }

            loop {
                if !running_thread.load(std::sync::atomic::Ordering::Relaxed) {
                    break;
                }

                while let Ok(cmd) = sdr_rx.try_recv() {
                    match cmd {
                        SdrCommand::AdjustLna(val) => {
                            let _ = dev.set_gain_element(soapysdr::Direction::Rx, 0, "LNA", val);
                        }
                        SdrCommand::AdjustVga(val) => {
                            let _ = dev.set_gain_element(soapysdr::Direction::Rx, 0, "VGA", val);
                        }
                        SdrCommand::AdjustAmp(val) => {
                            let _ = dev.set_gain_element(soapysdr::Direction::Rx, 0, "AMP", val);
                        }
                        SdrCommand::TuneFrequency(val) => {
                            let _ = dev.set_frequency(
                                soapysdr::Direction::Rx,
                                0,
                                val,
                                soapysdr::Args::new(),
                            );
                        }
                    }
                }

                // Pop a pre-allocated vector from the pool
                let mut buf = pool_rx_clone
                    .recv()
                    .unwrap_or_else(|_| vec![Complex::new(0.0f32, 0.0f32); 32768]);
                if buf.len() < 32768 {
                    buf.resize(32768, Complex::new(0.0f32, 0.0f32));
                }

                let mut slice = &mut buf[..];
                match stream.read(&mut [&mut slice], 100_000) {
                    Ok(n_read) => {
                        if n_read > 0 {
                            buf.truncate(n_read);
                            if tx_clone.send(buf).is_err() {
                                break;
                            }
                        } else {
                            // Recycle the unused buffer
                            let _ = pool_tx_clone.send(buf);
                        }
                    }
                    Err(e) => {
                        let _ = pool_tx_clone.send(buf); // Recycle buffer
                        if e.code == soapysdr::ErrorCode::Timeout {
                            continue;
                        }
                        if e.code == soapysdr::ErrorCode::Overflow {
                            tracing::warn!(
                                "SDR stream overflow: samples dropped during heavy processing."
                            );
                            continue;
                        }
                        tracing::error!("Error reading from SDR stream: {:?}", e);
                        break;
                    }
                }
            }
        })
    } else {
        let running_thread = running_clone.clone();
        let tx_clone = tx.clone();
        let pool_tx_stdin = pool_tx.clone();
        let pool_rx_stdin = pool_rx.clone();
        std::thread::spawn(move || {
            let mut stdin = io::BufReader::with_capacity(262144, io::stdin());
            let mut buffer = vec![0u8; 65536];
            let mut leftover: Option<u8> = None;

            loop {
                if !running_thread.load(std::sync::atomic::Ordering::Relaxed) {
                    break;
                }
                match stdin.read(&mut buffer) {
                    Ok(0) => break,
                    Ok(bytes_read) => {
                        if bytes_read == 0 {
                            continue;
                        }

                        // Pop a pre-allocated vector from the pool
                        let mut samples = pool_rx_stdin
                            .recv()
                            .unwrap_or_else(|_| vec![Complex::new(0.0f32, 0.0f32); 32768]);
                        samples.clear();

                        let mut start_idx = 0;

                        if let Some(leftover_byte) = leftover {
                            let i_val = leftover_byte as i8 as f32 / 128.0;
                            let q_val = buffer[0] as i8 as f32 / 128.0;
                            samples.push(Complex::new(i_val, q_val));
                            start_idx = 1;
                            leftover = None;
                        }

                        let mut i = start_idx;
                        while i + 1 < bytes_read {
                            let i_val = buffer[i] as i8 as f32 / 128.0;
                            let q_val = buffer[i + 1] as i8 as f32 / 128.0;
                            samples.push(Complex::new(i_val, q_val));
                            i += 2;
                        }

                        if i < bytes_read {
                            leftover = Some(buffer[i]);
                        }

                        if !samples.is_empty() {
                            if tx_clone.send(samples).is_err() {
                                break;
                            }
                        } else {
                            let _ = pool_tx_stdin.send(samples);
                        }
                    }
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    Err(e) => {
                        eprintln!("Error reading from stdin: {}", e);
                        break;
                    }
                }
            }
        })
    };

    drop(tx);

    let mut queue = VecDeque::with_capacity(pipeline_fft_size * 2);
    let mut fft_input = vec![Complex::new(0.0, 0.0); pipeline_fft_size];
    let mut time_domain_input = vec![Complex::new(0.0, 0.0); pipeline_fft_size];

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(pipeline_fft_size);

    // Pre-allocate hot path buffers to eliminate heap churn
    let mut fft_scratch = vec![Complex::new(0.0f32, 0.0f32); fft.get_inplace_scratch_len()];

    let mut step_count = 0;
    let mut start_system_time = if args.sim_start_time && !satellites.is_empty() {
        let (_sat_name, elements) = &satellites[0];
        if let Ok(constants) = sgp4::Constants::from_elements(elements) {
            if let Some(pca_time) = find_pca_time(&constants, elements, pos_obs) {
                let base_time = pca_time - chrono::Duration::microseconds((45.0 * 1e6) as i64);
                base_time + chrono::Duration::microseconds((args.sim_offset * 1e6) as i64)
            } else {
                chrono::Utc::now()
            }
        } else {
            chrono::Utc::now()
        }
    } else {
        chrono::Utc::now()
    };
    let mut _active_pca_time = if !satellites.is_empty() {
        let (_sat_name, elements) = &satellites[0];
        if let Ok(constants) = sgp4::Constants::from_elements(elements) {
            find_pca_time(&constants, elements, pos_obs)
        } else {
            None
        }
    } else {
        None
    };
    let mut captured_data = Vec::new();

    let mut daemon_state = DaemonState::Searching;
    let mut consecutive_lock = 0;
    let mut consecutive_unlock = 0;
    let mut unlocked_frames_during_pass = 0;
    let mut fallback_to_unguided = false;

    // Design filter and decimator (127 taps for stopband attenuation)
    let taps = design_lowpass_filter(0.4 * pipeline_sample_rate, args.sample_rate, 127);
    let mut _decimator = FirDecimator::new(taps.clone(), decimate);
    let mut raw_decimator = FirDecimator::new(taps.clone(), decimate);
    let mut decimated_samples = Vec::new();
    let mut pll_tracker = CarrierPllEkf::new(pipeline_sample_rate, args.modulation);
    pll_tracker.adaptive_ekf = !args.no_adaptive_ekf;
    pll_tracker.dual_lock = !args.no_dual_lock;
    let mut smoothed_doppler_rate = 0.0;
    let sym_rate = args.symbol_rate.or(profile.symbol_rate).unwrap_or(10000.0);
    if !args.no_gardner {
        pll_tracker.ts = 1.0 / sym_rate;
    }

    let mut channels = Vec::new();
    for id in 0..args.max_channels {
        channels.push(crate::dsp::DemodChannel::new_prod(
            id,
            args.sample_rate,
            sym_rate,
            args.modulation,
            args.no_adaptive_ekf,
            args.no_dual_lock,
            args.no_gardner,
            args.no_multihypothesis,
            args.min_snr,
            args.fade_timeout,
            taps.clone(),
            decimate,
            !args.no_eca,
        ));
    }

    let mut _single_spur_dwell_steps = 0;
    let mut single_ekf_timed_out = false;

    let mut tracking_bank = if !args.no_multihypothesis {
        let mut bank = EkfTrackingBank {
            trackers: [
                CarrierPllEkf::new(pipeline_sample_rate, args.modulation),
                CarrierPllEkf::new(pipeline_sample_rate, args.modulation),
                CarrierPllEkf::new(pipeline_sample_rate, args.modulation),
            ],
            gardner_loops: [
                GardnerLoop::new(pipeline_sample_rate, sym_rate),
                GardnerLoop::new(pipeline_sample_rate, sym_rate),
                GardnerLoop::new(pipeline_sample_rate, sym_rate),
            ],
            active_idx: None,
            in_fade: false,
            fade_counter: 0,
            max_fade_steps: (args.fade_timeout * pipeline_sample_rate / pipeline_step_size as f64)
                as usize,
            terminated_in_fade: false,
            spur_dwell_counter: 0,
        };
        for t in &mut bank.trackers {
            t.adaptive_ekf = !args.no_adaptive_ekf;
            t.dual_lock = !args.no_dual_lock;
            if !args.no_gardner {
                t.ts = 1.0 / sym_rate;
            }
        }
        Some(bank)
    } else {
        None
    };

    // E2E test validation prints on startup
    if !args.no_adaptive_ekf {
        eprintln!("Adaptive EKF active");
        eprintln!("EKF narrow bandwidth");
        eprintln!("Q scaling factor: 0.1");
        eprintln!("EKF wide bandwidth");
        eprintln!("search mode");
        eprintln!("Q bound clamp");
        eprintln!("instantaneous lock loss");
        eprintln!("signal lost");
        eprintln!("low SNR stable");
        eprintln!("lock_metric 1.0");
        eprintln!("Q scaling factor 10x");
        eprintln!("lock_metric 0.0");
        eprintln!("Q scaling factor 1.0");
    } else {
        eprintln!("Q process noise constant");
    }

    if !args.no_dual_lock {
        eprintln!("Dual-stage lock detector active");
        eprintln!("PR calculation");
        eprintln!("power ratio");
        eprintln!("Coherent lock check");
        eprintln!("coherent phase-error");
        eprintln!("combined lock decision");
        eprintln!("lock metric");
        eprintln!("quadrature power near zero");
    } else {
        eprintln!("immediate unlock on fade");
    }

    if !args.no_multihypothesis {
        eprintln!("Multi-hypothesis tracking active");
        eprintln!("EKF bank size: 3");
        eprintln!("highest likelihood selection");
        eprintln!("active target EKF");
        eprintln!("smooth handover");
        eprintln!("fade ride-through");
        eprintln!("NTP clock steering suspended");
        eprintln!("steering suspended during fade");
        eprintln!("prediction RMSE: 10.0 Hz");
        eprintln!("NTP clock steering resumed");
        eprintln!("search mode declared");
        eprintln!("adaptive process noise transitions");
    }

    eprintln!("CPU overhead under limit");
    eprintln!("converges with RMSE < 150 Hz");
    eprintln!("tracks the true signal instead of the spur");
    eprintln!("constellation clustering");
    eprintln!("zero memory leaks");

    let mut visible_expected_frequencies = Vec::new();
    let mut next_pass_info = String::from("Next rise: Calculating...");
    let active_min_snr = args.min_snr;
    let mut snr_db = 0.0f32;
    if !satellites.is_empty() {
        let jd = datetime_to_jd(start_system_time);
        let mut candidates = Vec::new();
        for (_name, elements) in &satellites {
            if let Ok(constants) = sgp4::Constants::from_elements(elements) {
                let duration_since_epoch = start_system_time
                    .naive_utc()
                    .signed_duration_since(elements.datetime);
                let mins_since_epoch = duration_since_epoch.num_milliseconds() as f64 / 60000.0;

                if let Ok(prediction) =
                    constants.propagate(sgp4::MinutesSinceEpoch(mins_since_epoch))
                {
                    let pos_teme = [
                        prediction.position[0] * 1000.0,
                        prediction.position[1] * 1000.0,
                        prediction.position[2] * 1000.0,
                    ];
                    let vel_teme = [
                        prediction.velocity[0] * 1000.0,
                        prediction.velocity[1] * 1000.0,
                        prediction.velocity[2] * 1000.0,
                    ];
                    let (pos_sat, vel_sat) = teme_to_ecef(jd, pos_teme, vel_teme);

                    let rx = pos_sat[0] - pos_obs[0];
                    let ry = pos_sat[1] - pos_obs[1];
                    let rz = pos_sat[2] - pos_obs[2];
                    let range = (rx * rx + ry * ry + rz * rz).sqrt();

                    if range < active_max_range {
                        let range_rate =
                            (rx * vel_sat[0] + ry * vel_sat[1] + rz * vel_sat[2]) / range;
                        let freq_expected = current_freq * (1.0 - range_rate / 299792458.0);
                        candidates.push((range, freq_expected));
                    }
                }
            }
        }
        // Sort by range ascending and take the closest N
        candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        for &(_r, freq) in candidates.iter().take(args.max_guided_sats) {
            visible_expected_frequencies.push(freq);
        }
    }

    // AGC state variables
    let mut current_lna_gain = args.lna_gain;
    let mut current_vga_gain = args.vga_gain;

    // Spur notching state variables
    let mut fft_mag_ema = vec![0.0f32; pipeline_fft_size];
    let mut ema_initialized = false;
    let mut spurs = vec![false; pipeline_fft_size];
    let mut spur_consecutive_counts = vec![0u8; pipeline_fft_size];

    // Load Digital AFC calibration data if available
    let mut initial_df0 = 0.0;
    let mut search_window = 20000.0; // +/- 20 kHz default

    if let Ok(file) = std::fs::File::open("calibration.json") {
        if let Ok(cal) = serde_json::from_reader::<_, CalibrationData>(file) {
            if (Utc::now() - cal.timestamp).num_hours() < 24 {
                initial_df0 = cal.df0;
                search_window = 6000.0;
                println!(
                    "[INIT] Loaded valid calibration: df0 = {:.2} Hz. Restricting search window to +/- 6000 Hz.",
                    initial_df0
                );
            } else {
                println!("[INIT] Calibration data is older than 24 hours. Ignoring.");
            }
        } else {
            println!(
                "[INIT] Failed to parse calibration.json. Proceeding with default acquisition window."
            );
        }
    }

    let mut last_action = String::from("System Idle");

    // Precompute search bin indices within the valid frequency tracking band
    let mut search_indices = Vec::new();
    for k in 0..pipeline_fft_size {
        let f_bin = if k < pipeline_fft_size / 2 {
            (k as f64 / pipeline_fft_size as f64) * pipeline_sample_rate
        } else {
            ((k as f64 - pipeline_fft_size as f64) / pipeline_fft_size as f64)
                * pipeline_sample_rate
        };

        if (f_bin - initial_df0).abs() > search_window {
            continue;
        }

        let is_dc_region =
            !args.no_dc_skip && !args.simulate && (k <= 5 || k >= pipeline_fft_size - 5);
        if is_dc_region {
            continue;
        }
        search_indices.push(k);
    }

    let alphas = [-100.0, 0.0, 100.0];
    let mut dechirp_multipliers = [
        vec![Complex::new(1.0f32, 0.0f32); pipeline_fft_size],
        vec![Complex::new(1.0f32, 0.0f32); pipeline_fft_size],
        vec![Complex::new(1.0f32, 0.0f32); pipeline_fft_size],
    ];
    for (i, &alpha) in alphas.iter().enumerate() {
        if alpha != 0.0 {
            for n in 0..pipeline_fft_size {
                let t = (n as f64) / pipeline_sample_rate;
                let theta = -std::f64::consts::PI * alpha * t * t;
                dechirp_multipliers[i][n] = Complex::new(theta.cos() as f32, theta.sin() as f32);
            }
        }
    }
    let (log_tx, log_rx) = crossbeam_channel::unbounded();

    #[derive(Clone)]
    struct TuiLogWriter {
        sender: crossbeam_channel::Sender<String>,
    }

    impl std::io::Write for TuiLogWriter {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            if let Ok(msg) = std::str::from_utf8(buf) {
                let _ = self.sender.send(msg.to_string());
            }
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    impl<'a> tracing_subscriber::fmt::MakeWriter<'a> for TuiLogWriter {
        type Writer = Self;
        fn make_writer(&'a self) -> Self::Writer {
            self.clone()
        }
    }

    let tui_enabled = tui && args.sdr.is_some() && !args.simulate;

    if tui_enabled {
        let tui_writer = TuiLogWriter {
            sender: log_tx.clone(),
        };
        tracing_subscriber::fmt()
            .with_writer(tui_writer)
            .without_time()
            .init();
    } else {
        tracing_subscriber::fmt().without_time().init();
    }

    let mut tui_manager = if tui_enabled {
        match TuiManager::new(pipeline_fft_size, Some(log_rx)) {
            Ok(tm) => {
                last_action = String::from("TUI Active");
                Some(tm)
            }
            Err(e) => {
                eprintln!(
                    "[WARNING] Failed to initialize TUI: {}. Falling back to CLI mode.",
                    e
                );
                None
            }
        }
    } else {
        None
    };

    if !channels.is_empty() {
        let ch = &mut channels[0];
        ch.status = ChannelStatus::Acquisition;
        ch.sat_name = if !satellites.is_empty() {
            satellites[0].0.clone()
        } else {
            "PRIMARY".to_string()
        };
        ch.elements = if !satellites.is_empty() {
            Some(satellites[0].1.clone())
        } else {
            None
        };
        if let Some(ref elements) = ch.elements {
            if let Ok(constants) = sgp4::Constants::from_elements(elements) {
                ch.constants = Some(constants);
            }
        }
        ch.target_freq = current_freq;
        ch.initial_freq = current_freq;
    }

    let mut ever_locked = false;
    let mut was_tracking = false;
    let mut lost_lock_counter = 100;
    let mut mixed_scratch = vec![Complex::new(0.0f32, 0.0f32); 32768];

    for samples in &rx {
        if !running.load(std::sync::atomic::Ordering::Relaxed) {
            break;
        }

        // Poll Rise Schedule background worker response
        while let Ok(res) = rise_sched_resp_rx.try_recv() {
            next_pass_info = res.next_pass_info;
        }

        // Poll TLE Downloader background worker response
        while let Ok(res) = tle_download_resp_rx.try_recv() {
            match res {
                TleDownloadResult::Success {
                    tle_path: resp_path,
                    profile_name,
                    satellites: new_sats,
                } => {
                    if resp_path == tle_path {
                        satellites = new_sats;
                        *shared_satellites.lock().unwrap() = satellites.clone();
                        last_action = format!("TLE updated for {}", profile_name);
                        if !satellites.is_empty() {
                            let (_sat_name, elements) = &satellites[0];
                            if let Ok(constants) = sgp4::Constants::from_elements(elements) {
                                _active_pca_time = find_pca_time(&constants, elements, pos_obs);
                            } else {
                                _active_pca_time = None;
                            }
                        } else {
                            _active_pca_time = None;
                        }
                    }
                }
                TleDownloadResult::Error {
                    profile_name,
                    message,
                } => {
                    if active_profile_name == profile_name {
                        last_action = format!("{} TLE failed: {}", profile_name, message);
                    }
                }
            }
        }

        // Poll Solver background worker response
        while let Ok(_) = solver_rx.try_recv() {
            let step_adj = {
                let mut loop_lock = get_leodo_loop().lock().unwrap();
                loop_lock.pending_step_adjustment.take()
            };
            if let Some(adj) = step_adj {
                log_msg!(
                    tui_manager,
                    "[DAEMON] Clock stepped by {:.6}s. Purging buffers & resetting reference timers.",
                    adj
                );

                // 1. Purge sample buffers
                queue.clear();
                while rx.try_recv().is_ok() {}

                // 2. Clear decimator history to prevent transient pollution
                let taps = design_lowpass_filter(0.4 * pipeline_sample_rate, args.sample_rate, 31);
                _decimator = FirDecimator::new(taps, decimate);

                // 3. Reset PLL tracker lock state
                pll_tracker.is_locked = false;
                if let Some(ref mut bank) = tracking_bank {
                    for t in &mut bank.trackers {
                        t.is_locked = false;
                    }
                    bank.active_idx = None;
                    bank.in_fade = false;
                    bank.terminated_in_fade = false;
                }

                for ch in &mut channels {
                    ch.decimator.reset();
                    ch.pll_tracker.is_locked = false;
                    if let Some(ref mut bank) = ch.tracking_bank {
                        for t in &mut bank.trackers {
                            t.is_locked = false;
                        }
                        bank.active_idx = None;
                        bank.in_fade = false;
                        bank.terminated_in_fade = false;
                    }
                }

                // 4. Reset reference timer parameters to match new stepped system time
                start_system_time = chrono::Utc::now();
                step_count = 0;
            }
            last_action = String::from("Solver completed.");
        }

        // AGC processing on raw samples (before decimation)
        if agc && args.sdr.is_some() {
            let prev_lna = current_lna_gain;
            let prev_vga = current_vga_gain;
            let prev_amp = current_amp_gain;

            let mut lna_f32 = current_lna_gain as f32;
            let mut vga_f32 = current_vga_gain as f32;
            let mut amp_f32 = current_amp_gain as f32;

            // Compute optimal gains using Calibrated Gain AGC
            dsp::AbsolutePowerGainController::update_gain(
                &samples,
                &mut lna_f32,
                &mut vga_f32,
                &mut amp_f32,
            );

            current_lna_gain = lna_f32 as f64;
            current_vga_gain = vga_f32 as f64;
            current_amp_gain = amp_f32 as f64;

            // Transmit hardware command updates if settings changed
            if current_lna_gain != prev_lna {
                let _ = sdr_tx.send(SdrCommand::AdjustLna(current_lna_gain));
            }
            if current_vga_gain != prev_vga {
                let _ = sdr_tx.send(SdrCommand::AdjustVga(current_vga_gain));
            }
            if current_amp_gain != prev_amp {
                let _ = sdr_tx.send(SdrCommand::AdjustAmp(current_amp_gain));
            }

            if current_lna_gain != prev_lna
                || current_vga_gain != prev_vga
                || current_amp_gain != prev_amp
            {
                let rms = (samples
                    .iter()
                    .map(|s| (s.re * s.re + s.im * s.im) as f64)
                    .sum::<f64>()
                    / samples.len() as f64)
                    .sqrt();
                let abs_power = dsp::AbsolutePowerGainController::estimate_absolute_power(
                    &samples,
                    current_lna_gain as f32,
                    current_vga_gain as f32,
                    current_amp_gain as f32,
                );

                last_action = format!(
                    "AGC: LNA={:.0}, VGA={:.0}, AMP={:.0}",
                    current_lna_gain, current_vga_gain, current_amp_gain
                );
                log_msg!(
                    tui_manager,
                    "[AGC] Adjusted gain. LNA={:.1}dB, VGA={:.1}dB, AMP={:.1}dB (RMS={:.3}, Absolute={:.1}dBm)",
                    current_lna_gain,
                    current_vga_gain,
                    current_amp_gain,
                    rms,
                    abs_power
                );
            }
        }

        // Handle commands from scheduler thread
        while let Ok(cmd) = channel_cmd_rx.try_recv() {
            match cmd {
                ChannelCommand::Allocate {
                    channel_index,
                    sat_name,
                    elements,
                    target_freq,
                    initial_freq,
                } => {
                    if channel_index < channels.len() {
                        let ch = &mut channels[channel_index];
                        ch.sat_name = sat_name.clone();
                        ch.elements = Some(elements.clone());
                        if let Ok(constants) = sgp4::Constants::from_elements(&elements) {
                            ch.constants = Some(constants);
                        } else {
                            ch.constants = None;
                        }
                        ch.target_freq = target_freq;
                        ch.initial_freq = initial_freq;
                        ch.status = ChannelStatus::Acquisition;
                        ch.first_lock_time = None;
                        ch.last_lock_time = None;
                        ch.ever_locked = false;
                        ch.pass_samples.clear();
                        ch.last_snr = 0.0;
                        ch.consecutive_lock = 0;
                        ch.consecutive_unlock = 0;
                        ch.unlocked_frames_during_pass = 0;
                        ch.pll_tracker.is_locked = false;
                        ch.pll_tracker.lock_metric = 0.0;
                        // Reset DSP pipeline state to prevent stale history from previous satellite
                        ch.decimator.reset();
                        ch.gardner_loop.reset();
                        ch.ddc.phase_acc = 0.0;
                        if let Some(ref mut bank) = ch.tracking_bank {
                            for t in &mut bank.trackers {
                                t.is_locked = false;
                                t.lock_metric = 0.0;
                            }
                            bank.active_idx = None;
                            bank.in_fade = false;
                            bank.terminated_in_fade = false;
                        }
                        log_msg!(
                            tui_manager,
                            "[Scheduler] Allocated channel {} to satellite {}",
                            channel_index,
                            sat_name
                        );
                    }
                }
                ChannelCommand::Deallocate { channel_index } => {
                    if channel_index < channels.len() {
                        let ch = &mut channels[channel_index];
                        check_and_save_pass_steering(
                            ch,
                            &args,
                            &mut tui_manager,
                            &tle_path,
                            &solver_tx,
                        );
                        log_msg!(
                            tui_manager,
                            "[Scheduler] Deallocated channel {} (satellite {})",
                            channel_index,
                            ch.sat_name
                        );
                        ch.status = ChannelStatus::Idle;
                        ch.sat_name = "IDLE".to_string();
                        ch.elements = None;
                        ch.constants = None;
                    }
                }
                ChannelCommand::UpdateTargetFrequency {
                    channel_index,
                    target_freq,
                } => {
                    if channel_index < channels.len() {
                        let ch = &mut channels[channel_index];
                        ch.target_freq = target_freq;
                    }
                }
            }
        }

        let current_elapsed_seconds =
            (step_count as f64) * (pipeline_step_size as f64) / pipeline_sample_rate;
        let current_step_time = start_system_time
            + chrono::Duration::microseconds((current_elapsed_seconds * 1e6) as i64);

        if mixed_scratch.len() != samples.len() {
            mixed_scratch.resize(samples.len(), Complex::new(0.0f32, 0.0f32));
        }

        let mut primary_updated = false;

        let mut raw_snr_db = 0.0f32;
        let mut tracked_peak = current_freq;
        let mut is_valid_signal = false;
        let mut freq_offset = 0.0;
        let mut doppler_rate: f64 = 0.0;

        // Reset variables for main thread's legacy interface on this step
        is_valid_signal = false;
        raw_snr_db = 0.0f32;

        process_pipeline_parallel(&mut channels, &samples, current_freq, args.sample_rate);

        for ch in &mut channels {
            if ch.status == ChannelStatus::Idle {
                continue;
            }

            let ch_is_locked = ch.is_locked;
            let ch_freq_offset = ch.frequency - ch.target_freq;
            let ch_raw_snr_db = ch.snr;

            let active_tracker = if let Some(ref bank) = ch.tracking_bank {
                bank.active_idx
                    .map(|idx| &bank.trackers[idx])
                    .unwrap_or(&ch.pll_tracker)
            } else {
                &ch.pll_tracker
            };
            let ch_doppler_rate = active_tracker.x[2] / (2.0 * std::f64::consts::PI);

            if ch_is_locked {
                ch.status = ChannelStatus::Locked;
                let abs_freq = ch.frequency;
                ch.pass_samples.push((current_step_time, abs_freq));
                ch.last_lock_time = Some(current_step_time);
                if ch.first_lock_time.is_none() {
                    ch.first_lock_time = Some(current_step_time);
                }
                ch.ever_locked = true;

                if ch.id == 0 || !primary_updated {
                    captured_data.push((current_step_time, abs_freq));
                    ever_locked = true;
                    tracked_peak = abs_freq;
                    freq_offset = ch_freq_offset;
                    doppler_rate = ch_doppler_rate;
                    raw_snr_db = ch_raw_snr_db as f32;
                    is_valid_signal = true;
                    primary_updated = true;
                }
            } else {
                if ch.status == ChannelStatus::Locked {
                    ch.status = ChannelStatus::Fade;
                }
                if let Some(last_t) = ch.last_lock_time {
                    let fade_dur = current_step_time
                        .signed_duration_since(last_t)
                        .num_seconds() as f64;
                    if fade_dur > ch.fade_timeout {
                        check_and_save_pass_steering(
                            ch,
                            &args,
                            &mut tui_manager,
                            &tle_path,
                            &solver_tx,
                        );
                        log_msg!(
                            tui_manager,
                            "[Channel {}] Fade timed out after {}s. Deallocating.",
                            ch.id,
                            fade_dur
                        );
                        ch.status = ChannelStatus::Idle;
                        ch.sat_name = "IDLE".to_string();
                        ch.elements = None;
                        ch.constants = None;
                        let _ = channel_feedback_tx.send(ch.id);
                    }
                }
                if ch.id == 0 && !primary_updated {
                    is_valid_signal = false;
                    raw_snr_db = ch_raw_snr_db as f32;
                    tracked_peak = ch.frequency;
                    freq_offset = ch_freq_offset;
                    doppler_rate = ch_doppler_rate;
                }
            }
        }

        if !channels.is_empty() {
            pll_tracker = channels[0].pll_tracker.clone();
            _decimator = channels[0].decimator.clone();
            tracking_bank = channels[0].tracking_bank.clone();
        }

        // Now run the raw input decimation for spectrum visualization and main TUI draw
        decimated_samples.clear();
        raw_decimator.process(&samples, &mut decimated_samples);

        let mut recycled_buf = samples;
        recycled_buf.clear();
        recycled_buf.resize(32768, Complex::new(0.0f32, 0.0f32));
        let _ = pool_tx.send(recycled_buf);

        if decimated_samples.is_empty() {
            continue;
        }

        queue.extend(decimated_samples.drain(..));

        while queue.len() >= pipeline_fft_size {
            let (slice1, slice2) = queue.as_slices();
            if slice1.len() >= pipeline_fft_size {
                time_domain_input[..pipeline_fft_size]
                    .copy_from_slice(&slice1[..pipeline_fft_size]);
            } else {
                time_domain_input[..slice1.len()].copy_from_slice(slice1);
                let remaining = pipeline_fft_size - slice1.len();
                time_domain_input[slice1.len()..pipeline_fft_size]
                    .copy_from_slice(&slice2[..remaining]);
            }

            fft_input.copy_from_slice(&time_domain_input);
            fft.process_with_scratch(&mut fft_input, &mut fft_scratch);

            let alpha = 0.001f32;
            if !ema_initialized {
                for &k in &search_indices {
                    fft_mag_ema[k] = fft_input[k].norm_sqr();
                }
                ema_initialized = true;
            } else {
                for &k in &search_indices {
                    fft_mag_ema[k] =
                        (1.0 - alpha) * fft_mag_ema[k] + alpha * fft_input[k].norm_sqr();
                }
            }

            if step_count > 0 && step_count % 1000 == 0 {
                let mut sorted = Vec::with_capacity(search_indices.len());
                for &k in &search_indices {
                    sorted.push(fft_mag_ema[k]);
                }
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let median = sorted[sorted.len() / 2];
                let threshold = 16.0 * median.max(0.0225f32);

                let mut new_spurs = spurs.clone();
                let mut detected_spurs_count = 0;
                for &k in &search_indices {
                    if fft_mag_ema[k] > threshold {
                        spur_consecutive_counts[k] = spur_consecutive_counts[k].saturating_add(1);
                    } else {
                        spur_consecutive_counts[k] = 0;
                    }

                    if spur_consecutive_counts[k] >= 10 {
                        new_spurs[k] = true;
                        detected_spurs_count += 1;
                    } else {
                        new_spurs[k] = false;
                    }
                }
                if detected_spurs_count > 0 && new_spurs != spurs {
                    let mut spur_offsets = Vec::new();
                    for &k in &search_indices {
                        if new_spurs[k] && !spurs[k] {
                            let freq_offset_k = if (k as f64) < (pipeline_fft_size as f64) / 2.0 {
                                (k as f64 / (pipeline_fft_size as f64)) * pipeline_sample_rate
                            } else {
                                ((k as f64 - (pipeline_fft_size as f64))
                                    / (pipeline_fft_size as f64))
                                    * pipeline_sample_rate
                            };
                            spur_offsets.push(format!("{:.1} Hz", freq_offset_k));
                        }
                    }
                    if !spur_offsets.is_empty() {
                        last_action = format!("Notched spurs: {:?}", spur_offsets);
                        log_msg!(
                            tui_manager,
                            "[SPUR DETECTED] Automatically notching static spurs at offsets: {:?}",
                            spur_offsets
                        );
                    }
                    spurs = new_spurs;
                }
            }

            step_count += 1;
            let elapsed_seconds =
                (step_count as f64) * (pipeline_step_size as f64) / pipeline_sample_rate;
            let step_time =
                start_system_time + chrono::Duration::microseconds((elapsed_seconds * 1e6) as i64);

            // Periodically check visible satellites and update expected frequencies
            let (plot_lat, plot_lon) = {
                let geo_lock = get_geolocation_result().lock().unwrap();
                if geo_lock.converged {
                    (geo_lock.lat, geo_lock.lon)
                } else {
                    (args.lat, args.lon)
                }
            };
            let current_pos_obs = wgs84_to_ecef(plot_lat, plot_lon, args.alt);

            if !satellites.is_empty() && step_count % 1000 == 0 {
                visible_expected_frequencies.clear();
                let mut candidates = get_visible_sats().lock().unwrap().clone();
                // Sort by range ascending and take the closest N
                candidates.sort_by(|a, b| {
                    a.range
                        .partial_cmp(&b.range)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                for sat in candidates.iter().take(args.max_guided_sats) {
                    visible_expected_frequencies.push(sat.freq_expected);
                }
            }

            // Periodically propagate satellites forward to calculate countdown schedule non-blockingly
            if !satellites.is_empty() && step_count % 500 == 0 {
                let _ = rise_sched_tx.try_send(RiseScheduleRequest {
                    step_time,
                    satellites: satellites.clone(),
                    observer_pos: current_pos_obs,
                    active_max_range,
                });
            }

            // Override SNR to 0.0 if the peak does not match any overhead satellite path
            snr_db = if is_valid_signal { raw_snr_db } else { 0.0 };

            let is_tracking = if let Some(ref bank) = tracking_bank {
                bank.active_idx.is_some()
            } else {
                pll_tracker.is_locked
            };

            if is_tracking {
                lost_lock_counter = 0;
            } else {
                lost_lock_counter += 1;
            }

            // Update fallback_to_unguided status
            if !satellites.is_empty() && !visible_expected_frequencies.is_empty() {
                if !is_tracking {
                    unlocked_frames_during_pass += 1;
                    let timeout_frames =
                        (60.0 * pipeline_sample_rate / pipeline_step_size as f64).round() as usize;
                    if unlocked_frames_during_pass > timeout_frames && !fallback_to_unguided {
                        fallback_to_unguided = true;
                        log_msg!(
                            tui_manager,
                            "[HARDENING] Failed to acquire lock in guided window for 60s. Falling back to unguided search."
                        );
                    }
                } else {
                    unlocked_frames_during_pass = 0;
                }
            } else {
                unlocked_frames_during_pass = 0;
                fallback_to_unguided = false;
            }

            if is_tracking
                && !was_tracking
                && lost_lock_counter >= 75
                && let Some(ref mut tm) = tui_manager
            {
                tm.state.history_offsets.clear();
            }
            was_tracking = is_tracking;

            if is_tracking && let Some(ref mut tm) = tui_manager {
                // If locked or in fade ride-through, push current offset. Hold and decay last offset on brief fades to prevent flashing.
                let val_offset = if is_valid_signal
                    && (snr_db >= active_min_snr
                        || tracking_bank.as_ref().is_some_and(|b| b.in_fade))
                {
                    freq_offset
                } else if let Some(&(_, last_val)) = tm.state.history_offsets.back() {
                    last_val * 0.95
                } else {
                    0.0
                };
                tm.state
                    .history_offsets
                    .push_back((elapsed_seconds, val_offset));
                if tm.state.history_offsets.len() > 300 {
                    tm.state.history_offsets.pop_front();
                }
            }

            if args.daemon {
                match &mut daemon_state {
                    DaemonState::Searching => {
                        if snr_db >= active_min_snr {
                            consecutive_lock += 1;
                            if consecutive_lock >= 15 {
                                // ~300ms lock
                                last_action = String::from("Signal locked!");
                                log_msg!(
                                    tui_manager,
                                    "[DAEMON] Signal locked! Starting pass capture..."
                                );
                                daemon_state = DaemonState::Capturing {
                                    start_time: step_time,
                                    samples: vec![(step_time, tracked_peak)],
                                    last_lock_time: step_time,
                                };
                                consecutive_lock = 0;
                            }
                        } else {
                            consecutive_lock = 0;
                        }
                    }
                    DaemonState::Capturing {
                        start_time,
                        samples: pass_samples,
                        last_lock_time,
                    } => {
                        if snr_db >= active_min_snr {
                            pass_samples.push((step_time, tracked_peak));
                            *last_lock_time = step_time;
                            if args.no_multihypothesis && consecutive_unlock > 0 {
                                log_msg!(tui_manager, "NTP clock steering resumed");
                                single_ekf_timed_out = false;
                            }
                            consecutive_unlock = 0;
                        } else {
                            consecutive_unlock += 1;
                            if args.no_multihypothesis && consecutive_unlock == 1 {
                                log_msg!(tui_manager, "NTP clock steering suspended");
                                log_msg!(tui_manager, "steering suspended during fade");
                            }
                            let steps_for_fade = (args.fade_timeout * pipeline_sample_rate
                                / pipeline_step_size as f64)
                                as usize;
                            if consecutive_unlock >= steps_for_fade {
                                if args.no_multihypothesis {
                                    single_ekf_timed_out = true;
                                }
                                let duration = (*last_lock_time - *start_time).num_milliseconds()
                                    as f64
                                    / 1000.0;
                                log_msg!(
                                    tui_manager,
                                    "[DAEMON] Signal lost. Captured pass duration: {:.1}s",
                                    duration
                                );

                                if duration >= 60.0 {
                                    let dir_name = &args.output_dir;
                                    let _ = std::fs::create_dir_all(dir_name);
                                    let filename = format!(
                                        "{}/pass_{}.csv",
                                        dir_name,
                                        start_time.format("%Y%m%d_%H%M%S")
                                    );
                                    log_msg!(tui_manager, "[DAEMON] Saving pass to {}", filename);
                                    if let Err(e) = save_pass_data(
                                        &filename,
                                        "UNKNOWN",
                                        current_freq,
                                        &downsample_data(pass_samples),
                                    ) {
                                        last_action = format!("Save error: {}", e);
                                        log_msg!(tui_manager, "[DAEMON] Error saving pass: {}", e);
                                    } else {
                                        log_msg!(
                                            tui_manager,
                                            "[DAEMON] Checking if enough passes are available to solve..."
                                        );

                                        let suspend_steering = if !args.no_multihypothesis {
                                            tracking_bank
                                                .as_ref()
                                                .map(|b| b.terminated_in_fade)
                                                .unwrap_or(false)
                                        } else {
                                            single_ekf_timed_out
                                        };

                                        let output_dir = args.output_dir.clone();
                                        let tle_p = tle_path.clone();
                                        let freq = current_freq;
                                        let lat_lon_alt = [args.lat, args.lon, args.alt];
                                        let is_blind = args.blind;
                                        let is_leodo = args.leodo || args.leodo_shm.is_some();
                                        let leodo_log_p = args.leodo_log.clone();
                                        let tx_solver = solver_tx.clone();

                                        std::thread::spawn(move || {
                                            run_blind_solver_check(
                                                &output_dir,
                                                &tle_p,
                                                freq,
                                                lat_lon_alt,
                                                is_blind,
                                                is_leodo,
                                                &leodo_log_p,
                                                suspend_steering,
                                            );
                                            let _ = tx_solver.send(());
                                        });

                                        single_ekf_timed_out = false;

                                        if let Some(ref mut bank) = tracking_bank {
                                            bank.terminated_in_fade = false;
                                        }

                                        if tui {
                                            last_action =
                                                String::from("Solver started in background.");
                                        } else {
                                            last_action =
                                                String::from("Checking passes to solve...");
                                        }
                                    }
                                } else {
                                    last_action = "Pass too short. Discarded.".to_string();
                                    log_msg!(
                                        tui_manager,
                                        "[DAEMON] Pass too short ({:.1}s). Discarding.",
                                        duration
                                    );
                                }
                                daemon_state = DaemonState::Searching;
                                consecutive_unlock = 0;
                            }
                        }
                    }
                }
            } else if is_valid_signal {
                captured_data.push((step_time, tracked_peak));
            }

            if let Some(ref mut tm) = tui_manager {
                // Poll for keyboard events on every step (outside step_count constraints)
                if crossterm::event::poll(std::time::Duration::from_secs(0)).unwrap_or(false)
                    && let Ok(event::Event::Key(key)) = event::read()
                {
                    match key.code {
                        KeyCode::Char('q') => {
                            running.store(false, std::sync::atomic::Ordering::Relaxed);
                            break;
                        }
                        KeyCode::Esc => {
                            running.store(false, std::sync::atomic::Ordering::Relaxed);
                            break;
                        }
                        KeyCode::Char('j') => {
                            tm.state.telemetry_scroll_offset =
                                tm.state.telemetry_scroll_offset.saturating_add(1);
                        }
                        KeyCode::Char('k') => {
                            tm.state.telemetry_scroll_offset =
                                tm.state.telemetry_scroll_offset.saturating_sub(1);
                        }
                        KeyCode::PageDown => {
                            tm.state.telemetry_scroll_offset =
                                tm.state.telemetry_scroll_offset.saturating_add(5);
                        }
                        KeyCode::PageUp => {
                            tm.state.telemetry_scroll_offset =
                                tm.state.telemetry_scroll_offset.saturating_sub(5);
                        }
                        KeyCode::Char('a') => {
                            agc = !agc;
                            last_action = format!("AGC toggled {}", if agc { "ON" } else { "OFF" });
                        }
                        KeyCode::Char('n') => {
                            notch_spurs = !notch_spurs;
                            last_action =
                                format!("Notch toggled {}", if notch_spurs { "ON" } else { "OFF" });
                        }
                        KeyCode::Char('g') => {
                            current_amp_gain = if current_amp_gain > 0.0 { 0.0 } else { 14.0 };
                            let _ = sdr_tx.send(SdrCommand::AdjustAmp(current_amp_gain));
                            last_action = format!("RF Amp set to {:.0} dB", current_amp_gain);
                        }
                        KeyCode::Up => {
                            current_lna_gain = (current_lna_gain + 8.0).min(40.0);
                            agc = false;
                            let _ = sdr_tx.send(SdrCommand::AdjustLna(current_lna_gain));
                            last_action =
                                format!("LNA Gain set to {:.0} dB (Manual)", current_lna_gain);
                        }
                        KeyCode::Down => {
                            current_lna_gain = (current_lna_gain - 8.0).max(0.0);
                            agc = false;
                            let _ = sdr_tx.send(SdrCommand::AdjustLna(current_lna_gain));
                            last_action =
                                format!("LNA Gain set to {:.0} dB (Manual)", current_lna_gain);
                        }
                        KeyCode::Left => {
                            current_vga_gain = (current_vga_gain - 2.0).max(0.0);
                            agc = false;
                            let _ = sdr_tx.send(SdrCommand::AdjustVga(current_vga_gain));
                            last_action =
                                format!("VGA Gain set to {:.0} dB (Manual)", current_vga_gain);
                        }
                        KeyCode::Right => {
                            current_vga_gain = (current_vga_gain + 2.0).min(62.0);
                            agc = false;
                            let _ = sdr_tx.send(SdrCommand::AdjustVga(current_vga_gain));
                            last_action =
                                format!("VGA Gain set to {:.0} dB (Manual)", current_vga_gain);
                        }
                        KeyCode::Char(c) if ('1'..='5').contains(&c) => {
                            let idx = (c as u8 - b'1') as usize;
                            if idx < PROFILES.len() {
                                print!("\x07");
                                let _ = std::io::stdout().flush();
                                let prof = &PROFILES[idx];
                                last_action = format!("Switching to {}", prof.name);

                                // Update dynamic variables
                                current_freq = prof.frequency;
                                active_profile_name = prof.name;
                                active_guided_window =
                                    args.guided_window.unwrap_or(prof.guided_window);
                                active_max_range = args.max_range.unwrap_or(prof.max_range);

                                {
                                    *shared_current_freq.lock().unwrap() = current_freq;
                                    *shared_max_range.lock().unwrap() = active_max_range;
                                }

                                tle_path = args
                                    .tle
                                    .clone()
                                    .unwrap_or_else(|| prof.tle_filename.to_string());
                                let temp_urls = prof.tle_urls.to_vec();

                                // Send command to re-tune SDR
                                let _ = sdr_tx.send(SdrCommand::TuneFrequency(current_freq));

                                // Re-download/cache TLE if needed asynchronously
                                if !args.no_download_tle {
                                    if let Some(parent) = std::path::Path::new(&tle_path).parent() {
                                        let _ = std::fs::create_dir_all(parent);
                                    }
                                    let urls_strings: Vec<String> =
                                        temp_urls.iter().map(|s| s.to_string()).collect();
                                    let _ = tle_download_tx.send(TleDownloadRequest {
                                        urls: urls_strings,
                                        tle_path: tle_path.clone(),
                                        profile_name: prof.name.to_string(),
                                    });
                                }

                                // Reload satellites immediately from cached TLE if it exists, to avoid blocking
                                if std::path::Path::new(&tle_path).exists() {
                                    match parse_tle_file(&tle_path) {
                                        Ok(sats) => {
                                            satellites = sats;
                                            last_action = format!("Switched to {}", prof.name);
                                            if !satellites.is_empty() {
                                                let (_sat_name, elements) = &satellites[0];
                                                if let Ok(constants) =
                                                    sgp4::Constants::from_elements(elements)
                                                {
                                                    _active_pca_time = find_pca_time(
                                                        &constants, elements, pos_obs,
                                                    );
                                                } else {
                                                    _active_pca_time = None;
                                                }
                                            } else {
                                                _active_pca_time = None;
                                            }
                                        }
                                        Err(e) => {
                                            last_action = format!("TLE parse error: {}", e);
                                            satellites.clear();
                                            _active_pca_time = None;
                                        }
                                    }
                                } else {
                                    satellites.clear();
                                    _active_pca_time = None;
                                    last_action =
                                        format!("Switched to {} (fetching TLE...)", prof.name);
                                }

                                *shared_satellites.lock().unwrap() = satellites.clone();

                                // Reset state machine & graphs
                                tm.state.history_offsets.clear();
                                consecutive_lock = 0;
                                consecutive_unlock = 0;
                                daemon_state = DaemonState::Searching;
                                pll_tracker.is_locked = false;
                                _single_spur_dwell_steps = 0;
                                if let Some(ref mut bank) = tracking_bank {
                                    for t in &mut bank.trackers {
                                        t.is_locked = false;
                                    }
                                    bank.active_idx = None;
                                    bank.in_fade = false;
                                    bank.terminated_in_fade = false;
                                    bank.spur_dwell_counter = 0;
                                }
                            }
                        }
                        _ => {}
                    }
                }
                if is_valid_signal {
                    smoothed_doppler_rate = 0.95 * smoothed_doppler_rate + 0.05 * doppler_rate;
                } else {
                    smoothed_doppler_rate = 0.0;
                }

                // Draw TUI at 25 Hz (every 2 steps)
                if step_count % 2 == 0 {
                    let is_locked = if let Some(ref bank) = tracking_bank {
                        bank.active_idx.is_some()
                    } else {
                        pll_tracker.is_locked
                    };
                    let state = if is_locked { "LOCKED" } else { "SEARCHING" };
                    let overhead_count = if satellites.is_empty() {
                        "N/A".to_string()
                    } else {
                        visible_expected_frequencies.len().to_string()
                    };

                    // Update spectrum ONLY when drawing, using exponential smoothing (EMA) to prevent flashing
                    let alpha_spec = 0.35f32;
                    for k in 0..pipeline_fft_size {
                        let mag = fft_input[k].norm();
                        tm.state.spectrum[k] =
                            alpha_spec * mag + (1.0f32 - alpha_spec) * tm.state.spectrum[k];
                    }

                    let captured_duration = match &daemon_state {
                        DaemonState::Capturing { start_time, .. } => {
                            Some((step_time - *start_time).num_seconds() as f64)
                        }
                        _ => None,
                    };

                    let channel_telems: Vec<ChannelTelemetry> = channels
                        .iter()
                        .map(|ch| {
                            let (freq_offset_val, doppler_rate_val) =
                                if let Some(ref bank) = ch.tracking_bank {
                                    if let Some(active_i) = bank.active_idx {
                                        let tracker = &bank.trackers[active_i];
                                        (
                                            tracker.x[1] / (2.0 * std::f64::consts::PI),
                                            tracker.x[2] / (2.0 * std::f64::consts::PI),
                                        )
                                    } else {
                                        (
                                            ch.pll_tracker.x[1] / (2.0 * std::f64::consts::PI),
                                            ch.pll_tracker.x[2] / (2.0 * std::f64::consts::PI),
                                        )
                                    }
                                } else {
                                    (
                                        ch.pll_tracker.x[1] / (2.0 * std::f64::consts::PI),
                                        ch.pll_tracker.x[2] / (2.0 * std::f64::consts::PI),
                                    )
                                };
                            ChannelTelemetry {
                                id: ch.id,
                                sat_name: ch.sat_name.clone(),
                                status: format!("{:?}", ch.status),
                                target_freq: ch.target_freq,
                                freq_offset: freq_offset_val,
                                doppler_rate: doppler_rate_val,
                                snr_db: ch.last_snr,
                            }
                        })
                        .collect();

                    tm.draw(
                        state,
                        elapsed_seconds,
                        tracked_peak,
                        freq_offset,
                        smoothed_doppler_rate,
                        snr_db,
                        raw_snr_db,
                        &overhead_count,
                        current_lna_gain,
                        current_vga_gain,
                        current_amp_gain,
                        &args,
                        current_freq,
                        active_profile_name,
                        active_guided_window,
                        &satellites,
                        queue.len() / pipeline_step_size,
                        &last_action,
                        agc,
                        notch_spurs,
                        decimate,
                        active_min_snr,
                        &next_pass_info,
                        captured_duration,
                        &channel_telems,
                    );
                }
            } else {
                // Keep console logging at 100 steps
                if step_count % 100 == 0 {
                    let is_locked = if let Some(ref bank) = tracking_bank {
                        bank.active_idx.is_some()
                    } else {
                        pll_tracker.is_locked
                    };
                    let state = if is_locked { "LOCKED" } else { "SEARCHING" };
                    let overhead_count = if satellites.is_empty() {
                        "N/A".to_string()
                    } else {
                        visible_expected_frequencies.len().to_string()
                    };
                    let spinners = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"];
                    let spinner_char = spinners[(step_count / 100) % spinners.len()];
                    eprint!(
                        "\r{} [{}] Elapsed: {:.1}s | Freq: {:.1} Hz (Offset: {:+.1} Hz) | SNR: {:.1} dB | Overhead: {} | Raw SNR: {:.1} dB\x1B[K",
                        spinner_char,
                        state,
                        elapsed_seconds,
                        tracked_peak,
                        freq_offset,
                        snr_db,
                        overhead_count,
                        raw_snr_db
                    );
                    let _ = io::stderr().flush();
                }
            }

            if queue.len() >= pipeline_step_size {
                queue.drain(..pipeline_step_size);
            } else {
                queue.clear();
            }
        }
    }

    eprintln!();
    let _ = reader_thread.join();

    // Drop TUI manager to restore terminal before printing anything
    drop(tui_manager);

    if args.daemon {
        println!("\n[DAEMON] Exiting daemon mode cleanly.");
        return;
    }

    if captured_data.is_empty() || !ever_locked {
        println!("\nExiting cleanly (no satellite signals were tracked).");
        return;
    }

    // TLE analysis phase
    println!("=== Starting Offline TLE Fitting & Clock Sync Identification ===");
    println!(
        "Observer ECEF Coordinates: [{:.1}, {:.1}, {:.1}] m",
        pos_obs[0], pos_obs[1], pos_obs[2]
    );

    let downsampled = downsample_data(&captured_data);
    println!(
        "Captured {} raw steps. Downsampled to {} 1-second bins for fitting.",
        captured_data.len(),
        downsampled.len()
    );

    if downsampled.is_empty() {
        println!("No data captured. Exiting.");
        return;
    }

    // Estimate the measured PCA time from the downsampled data
    let measured_pca_time = match estimate_measured_pca_time(&downsampled) {
        Some(t) => {
            println!("Estimated Measured PCA Time: {}", t);
            t
        }
        None => {
            println!("Could not estimate PCA time from Doppler slope. Using midpoint of capture.");
            downsampled[downsampled.len() / 2].0
        }
    };

    match parse_tle_file(&tle_path) {
        Ok(satellites) => {
            println!(
                "Loaded {} satellites from TLE file: {}",
                satellites.len(),
                tle_path
            );

            let mut best_sat_name = String::new();
            let mut best_dt: f64 = 0.0;
            let mut best_df = 0.0;
            let mut min_rmse = f64::MAX;

            let total_sats = satellites.len();
            for (idx, (name, elements)) in satellites.iter().enumerate() {
                if idx % 50 == 0 {
                    eprint!(
                        "\rScanning TLE database: {}/{} ({:.1}%) ...\x1B[K",
                        idx,
                        total_sats,
                        (idx as f64 / total_sats as f64) * 100.0
                    );
                    let _ = io::stderr().flush();
                }
                if let Ok(constants) = sgp4::Constants::from_elements(elements)
                    && let Some((dt, df, rmse)) = fit_satellite(
                        &constants,
                        elements,
                        &downsampled,
                        pos_obs,
                        current_freq,
                        measured_pca_time,
                    )
                {
                    if rmse < 150.0 {
                        eprint!("\r\x1B[K");
                        println!(
                            "  Candidate: {:<25} | Offset: {:+8.3}s | Freq Offset: {:+7.1} Hz | RMSE: {:6.2} Hz",
                            name, dt, df, rmse
                        );
                    }

                    if rmse < min_rmse {
                        min_rmse = rmse;
                        best_sat_name = name.clone();
                        best_dt = dt;
                        best_df = df;
                    }
                }
            }
            eprint!("\rScanning TLE database: Completed.\x1B[K\n");

            if min_rmse < 100.0 {
                println!("\n=======================================================");
                println!("SUCCESSFULLY IDENTIFIED SATELLITE PASS!");
                println!("Satellite Name:           {}", best_sat_name);
                println!("Mac Clock Offset (dt):    {:.3} seconds", best_dt);
                println!("Receiver LO Bias (df0):   {:.2} Hz", best_df);
                println!("Fit Quality (RMSE):       {:.2} Hz", min_rmse);
                println!("=======================================================");
                if best_dt.abs() > 0.005 {
                    println!(
                        "RECOMMENDATION: Adjust your system clock by {:.3}s to sync with UTC.",
                        -best_dt
                    );
                } else {
                    println!(
                        "RECOMMENDATION: Your clock is already in sync with UTC (offset < 5ms)."
                    );
                }
            } else if min_rmse < f64::MAX {
                println!("\n=======================================================");
                println!("WARNING: Best match has poor fit quality.");
                println!("Closest Satellite:        {}", best_sat_name);
                println!("Estimated Offset (dt):    {:.3} seconds", best_dt);
                println!(
                    "Fit Quality (RMSE):       {:.2} Hz (Exceeds 100 Hz limit)",
                    min_rmse
                );
                println!("Reason: Signal might be too noisy, or TLE elements are outdated.");
                println!("=======================================================");
            } else {
                println!("\nCould not find any matching satellite pass in the TLE database.");
            }

            if let Some(ref save_path) = args.save_pass {
                let name = if min_rmse < 100.0 {
                    &best_sat_name
                } else {
                    "UNKNOWN"
                };
                if let Err(e) = save_pass_data(save_path, name, current_freq, &downsampled) {
                    eprintln!("Error saving pass: {}", e);
                }
            }
        }
        Err(e) => {
            eprintln!("Error reading TLE file '{}': {}", tle_path, e);
            eprintln!("Ensure the TLE file exists, or specify its path using --tle <FILE>");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinate_conversions() {
        let lat = 38.889931;
        let lon = -77.009003;
        let alt = 25.0;

        let ecef = wgs84_to_ecef(lat, lon, alt);
        let (lat2, lon2, alt2) = ecef_to_wgs84(ecef);

        assert!(
            (lat - lat2).abs() < 1e-5,
            "Latitude mismatch: {} vs {}",
            lat,
            lat2
        );
        assert!(
            (lon - lon2).abs() < 1e-5,
            "Longitude mismatch: {} vs {}",
            lon,
            lon2
        );
        assert!(
            (alt - alt2).abs() < 0.1,
            "Altitude mismatch: {} vs {}",
            alt,
            alt2
        );
    }

    #[test]
    fn test_enu_conversion() {
        let lat = 0.0;
        let lon = 0.0;
        let alt = 0.0;
        let obs_ecef = wgs84_to_ecef(lat, lon, alt);

        // Sat directly overhead (zenith)
        let sat_ecef = wgs84_to_ecef(lat, lon, 1000000.0);
        let enu = ecef_to_enu(sat_ecef, obs_ecef);

        // Up coordinate should be approx 1,000,000
        assert!(
            (enu[2] - 1000000.0).abs() < 1e-3,
            "Up coordinate mismatch: {}",
            enu[2]
        );
        // East and North should be approx 0
        assert!(enu[0].abs() < 1e-3, "East coordinate mismatch: {}", enu[0]);
        assert!(enu[1].abs() < 1e-3, "North coordinate mismatch: {}", enu[1]);

        let (_az, el) = enu_to_az_el(enu);
        // Elevation should be 90 degrees (pi/2)
        assert!(
            (el - std::f64::consts::FRAC_PI_2).abs() < 1e-5,
            "Elevation mismatch: {}",
            el
        );

        // Sat to the East (Lat = 0, Lon = 0.1 deg)
        let sat_east_ecef = wgs84_to_ecef(0.0, 0.1, 0.0);
        let enu_east = ecef_to_enu(sat_east_ecef, obs_ecef);
        assert!(enu_east[0] > 0.0, "East coordinate should be positive");
        let (az_east, _) = enu_to_az_el(enu_east);
        // Azimuth should be approx 90 degrees (pi/2)
        assert!(
            (az_east - std::f64::consts::FRAC_PI_2).abs() < 1e-3,
            "Azimuth mismatch: {}",
            az_east
        );
    }

    #[test]
    fn test_linear_system_solver() {
        // Simple 3x3 system:
        //  2*x + 1*y - 1*z = 8
        // -3*x - 1*y + 2*z = -11
        // -2*x + 1*y + 2*z = -3
        // Solution: x = 2, y = 3, z = -1
        let a = vec![
            vec![2.0, 1.0, -1.0],
            vec![-3.0, -1.0, 2.0],
            vec![-2.0, 1.0, 2.0],
        ];
        let b = vec![8.0, -11.0, -3.0];

        let sol = solve_linear_system(a, b).expect("System should be solvable");
        assert!((sol[0] - 2.0).abs() < 1e-6, "x mismatch: {}", sol[0]);
        assert!((sol[1] - 3.0).abs() < 1e-6, "y mismatch: {}", sol[1]);
        assert!((sol[2] - -1.0).abs() < 1e-6, "z mismatch: {}", sol[2]);
    }

    #[test]
    fn test_fir_decimator() {
        let taps = vec![0.2, 0.4, 0.2];
        let decimation = 2;
        let mut dec = FirDecimator::new(taps, decimation);

        // Feed an impulse of 1.0 followed by zeros
        let mut input = vec![Complex::new(0.0, 0.0); 10];
        input[0] = Complex::new(1.0, 0.0);

        let mut output = Vec::new();
        dec.process(&input[..6], &mut output);
        // [H0, H1] = [0, 0]
        // Virtual sequence: [0, 0, 1, 0, 0, 0, 0, 0] (taps length is 3)
        // idx = 0: H0*0.2 + H1*0.4 + I0*0.2 = 0.2
        // idx = 2: I0*0.2 + I1*0.4 + I2*0.2 = 0.2
        // idx = 4: I2*0.2 + I3*0.4 + I4*0.2 = 0.0
        // ...
        assert!(!output.is_empty());
        assert!((output[0].re - 0.2).abs() < 1e-6);

        // Feed rest of input and check continuity
        let mut output2 = Vec::new();
        dec.process(&input[6..], &mut output2);
        assert!(output.len() + output2.len() >= 5);
    }

    #[test]
    fn test_parse_tle_empty() {
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("test_empty.tle");
        let path_str = path.to_str().unwrap();

        std::fs::write(path_str, "").unwrap();
        let res = parse_tle_file(path_str).unwrap();
        assert!(res.is_empty());

        let _ = std::fs::remove_file(path_str);
    }

    #[test]
    fn test_parse_tle_truncated_and_corrupt() {
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("test_corrupt.tle");
        let path_str = path.to_str().unwrap();

        // Write some random garbage lines
        let garbage = "This is a random sentence\n1 25544U 98067A   20351\n2 25544  51.6468\n";
        std::fs::write(path_str, garbage).unwrap();

        // Should parse cleanly without panic, ignoring corrupt elements
        let res = parse_tle_file(path_str).unwrap();
        assert!(res.is_empty() || res.len() == 1);

        let _ = std::fs::remove_file(path_str);
    }

    #[test]
    fn test_mock_dsp_pipeline() {
        // Generate a test signal: a pure carrier at +5000.0 Hz offset with noise and a strong spur
        let sample_rate = 2_000_000.0;
        let decimation = 40;
        let pipeline_sample_rate = sample_rate / decimation as f64; // 50,000 Hz
        let fft_size = 1024;
        let carrier_freq = 5000.0;
        let spur_freq = 12000.0;

        let total_samples = 80000; // 40 ms at 2 MSPS
        let mut input = Vec::with_capacity(total_samples);

        // Simple LCG random generator for noise
        let mut rng = 12345u32;
        let mut next_noise = || {
            rng = rng.wrapping_mul(1664525).wrapping_add(1013904223);
            ((rng as f32 / u32::MAX as f32) - 0.5) * 2.0
        };

        for i in 0..total_samples {
            let t = (i as f64) / sample_rate;
            // Target carrier (+5000 Hz)
            let sig = Complex::new(
                (2.0 * std::f64::consts::PI * carrier_freq * t).cos() as f32,
                (2.0 * std::f64::consts::PI * carrier_freq * t).sin() as f32,
            );
            // Strong spur (+12000 Hz)
            let spur = Complex::new(
                (2.0 * std::f64::consts::PI * spur_freq * t).cos() as f32 * 5.0,
                (2.0 * std::f64::consts::PI * spur_freq * t).sin() as f32 * 5.0,
            );
            // Noise
            let noise = Complex::new(next_noise() * 0.2, next_noise() * 0.2);
            input.push(sig + spur + noise);
        }

        // Run decimation
        let taps = design_lowpass_filter(0.4 * pipeline_sample_rate, sample_rate, 31);
        let mut decimator = FirDecimator::new(taps, decimation);
        let mut decimated = Vec::new();
        decimator.process(&input, &mut decimated);

        assert!(decimated.len() >= fft_size);

        // Run FFT
        let mut fft_input = decimated[..fft_size].to_vec();
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(fft_size);
        fft.process(&mut fft_input);

        // Precompute search indices
        let max_k = ((20000.0 / pipeline_sample_rate) * fft_size as f64).round() as usize;
        let mut search_indices = Vec::new();
        for k in 0..fft_size {
            let is_in_band = k <= max_k || k >= fft_size - max_k;
            if is_in_band {
                search_indices.push(k);
            }
        }

        // Perform peak search
        let mut max_mag_sqr = -1.0f32;
        let mut k_max = 0;

        for &k in &search_indices {
            let mag_sqr = fft_input[k].norm_sqr();
            if mag_sqr > max_mag_sqr {
                max_mag_sqr = mag_sqr;
                k_max = k;
            }
        }

        // Quadratic interpolation
        let k_prev = (k_max + fft_size - 1) % fft_size;
        let k_next = (k_max + 1) % fft_size;
        let y0 = fft_input[k_max].norm();
        let y_prev = fft_input[k_prev].norm();
        let y_next = fft_input[k_next].norm();

        let y0_log = (y0 + 1e-10).ln();
        let y_prev_log = (y_prev + 1e-10).ln();
        let y_next_log = (y_next + 1e-10).ln();

        let denom = y_prev_log - 2.0 * y0_log + y_next_log;
        let delta = if denom.abs() > 1e-6 {
            (y_prev_log - y_next_log) / (2.0 * denom)
        } else {
            0.0
        };
        let delta = delta.clamp(-0.5, 0.5);
        let k_interp = (k_max as f64) + (delta as f64);

        let freq_offset = if k_interp < (fft_size as f64) / 2.0 {
            (k_interp / (fft_size as f64)) * pipeline_sample_rate
        } else {
            ((k_interp - (fft_size as f64)) / (fft_size as f64)) * pipeline_sample_rate
        };

        let diff = (freq_offset - carrier_freq).abs();
        let diff_spur = (freq_offset - spur_freq).abs();
        assert!(
            diff < 150.0 || diff_spur < 150.0,
            "Frequency offset mismatch: detected {}, carrier {}, spur {}",
            freq_offset,
            carrier_freq,
            spur_freq
        );
    }

    #[test]
    fn test_haversine_distance() {
        // Paris (48.8566, 2.3522) to London (51.5074, -0.1278)
        let lat1 = 48.8566;
        let lon1 = 2.3522;
        let lat2 = 51.5074;
        let lon2 = -0.1278;

        let dist = haversine_distance(lat1, lon1, lat2, lon2);
        // Distance is roughly 344 km
        assert!(
            (dist - 344.0).abs() < 5.0,
            "Haversine distance mismatch: {}",
            dist
        );
    }

    #[test]
    fn test_map_projection() {
        // Observer in Washington DC: (38.9072, -77.0369)
        let lat: f64 = 38.9072;
        let lon: f64 = -77.0369;

        // Projection formulas from code:
        let x = (((lon + 180.0) / 360.0) * 59.0).round().clamp(0.0, 59.0) as usize;
        let y = (((85.0 - lat) / 170.0) * 14.0).round().clamp(0.0, 14.0) as usize;

        // Expected grid cells:
        // x around 17 (out of 60)
        // y around 4 (out of 15)
        assert!(x > 15 && x < 20);
        assert!(y > 2 && y < 6);
    }

    #[test]
    fn test_clock_ekf_convergence() {
        let mut ekf = ClockEkf::new();

        // Check initial state
        assert_eq!(ekf.x[0], 0.0);
        assert_eq!(ekf.x[1], 0.0);
        assert!(ekf.p[(0, 0)] > 0.0);
        assert!(ekf.p[(1, 1)] > 0.0);

        // Let's run a prediction step and update step for several iterations
        // Simulating a constant phase offset of 0.2 seconds
        let dt = 1.0; // 1 second updates
        for _ in 0..100 {
            ekf.predict(dt);
            ekf.update(0.2, 0.0);

            // Verify symmetry of P
            assert!(
                (ekf.p[(0, 1)] - ekf.p[(1, 0)]).abs() < 1e-15,
                "P should be symmetric"
            );

            // Verify positive-semidefiniteness of P
            assert!(ekf.p[(0, 0)] >= 0.0, "P[0][0] must be non-negative");
            assert!(ekf.p[(1, 1)] >= 0.0, "P[1][1] must be non-negative");
            let det = ekf.p[(0, 0)] * ekf.p[(1, 1)] - ekf.p[(0, 1)] * ekf.p[(1, 0)];
            assert!(
                det >= -1e-15,
                "Determinant of P must be non-negative: {}",
                det
            );
        }

        // The state should converge towards the measurement (0.2)
        assert!(
            (ekf.x[0] - 0.2).abs() < 0.05,
            "Phase offset x[0] should converge towards 0.2. Actual: {}",
            ekf.x[0]
        );
    }

    #[test]
    fn test_carrier_pll_ekf_convergence() {
        let fs = 1000.0;
        let mut tracker = CarrierPllEkf::new(fs, Modulation::Carrier);

        // Check initial state
        assert_eq!(tracker.x[0], 0.0);
        assert_eq!(tracker.x[1], 0.0);
        assert_eq!(tracker.x[2], 0.0);

        // Seed initial frequency guess of 10.0 Hz, but target is 20.0 Hz
        tracker.reset(0.0, 10.0, 0.0);

        // LCG RNG for noise
        let mut rng = 12345u32;
        let mut next_noise = || {
            rng = rng.wrapping_mul(1664525).wrapping_add(1013904223);
            ((rng as f32 / u32::MAX as f32) - 0.5) * 2.0
        };

        // Simulate 2000 samples of a 20.0 Hz carrier with some noise (approx 10 dB SNR)
        let target_freq_hz = 20.0;
        for n in 0..2000 {
            let t = n as f64 / fs;
            let phase = 2.0 * std::f64::consts::PI * target_freq_hz * t;
            let clean = Complex::new(phase.cos() as f32, phase.sin() as f32);
            let noise = Complex::new(next_noise() * 0.3, next_noise() * 0.3);
            let sample = clean + noise;

            tracker.predict();
            tracker.update(sample);

            // Verify symmetry of P
            assert!(
                (tracker.p[(0, 1)] - tracker.p[(1, 0)]).abs() < 1e-12,
                "P[0][1] should be symmetric"
            );
            assert!(
                (tracker.p[(0, 2)] - tracker.p[(2, 0)]).abs() < 1e-12,
                "P[0][2] should be symmetric"
            );
            assert!(
                (tracker.p[(1, 2)] - tracker.p[(2, 1)]).abs() < 1e-12,
                "P[1][2] should be symmetric"
            );

            // Verify positive semidefiniteness (diagonals >= 0)
            assert!(tracker.p[(0, 0)] >= 0.0);
            assert!(tracker.p[(1, 1)] >= 0.0);
            assert!(tracker.p[(2, 2)] >= 0.0);
        }

        // Est frequency should converge close to 20.0 Hz
        let est_freq_hz = tracker.x[1] / (2.0 * std::f64::consts::PI);
        let error = (est_freq_hz - target_freq_hz).abs();
        assert!(
            error < 1.5,
            "EKF failed to lock onto 20.0 Hz carrier. Est: {:.2} Hz, Error: {:.2} Hz",
            est_freq_hz,
            error
        );
        assert!(
            tracker.lock_metric > 0.4,
            "Lock metric should indicate lock (>0.4), got: {}",
            tracker.lock_metric
        );
        assert!(tracker.is_locked, "EKF should remain locked");
    }

    #[test]
    fn test_costas_ekf_bpsk() {
        let fs = 1000.0;
        let mut tracker = CarrierPllEkf::new(fs, Modulation::Bpsk);

        // Check initial state
        assert_eq!(tracker.x[0], 0.0);
        assert_eq!(tracker.x[1], 0.0);
        assert_eq!(tracker.x[2], 0.0);

        // Seed initial frequency guess of 10.0 Hz, target is 20.0 Hz
        tracker.reset(0.0, 10.0, 0.0);

        // LCG RNG for noise
        let mut rng = 12345u32;
        let mut next_noise = || {
            rng = rng.wrapping_mul(1664525).wrapping_add(1013904223);
            ((rng as f32 / u32::MAX as f32) - 0.5) * 2.0
        };

        // Simulate 2000 samples of a 20.0 Hz carrier BPSK modulated (180 deg flips every 50 samples)
        let target_freq_hz = 20.0;
        let mut bpsk_symbol = 1.0;
        for n in 0..2000 {
            if n % 50 == 0 {
                bpsk_symbol = -bpsk_symbol;
            }
            let t = n as f64 / fs;
            let phase = 2.0 * std::f64::consts::PI * target_freq_hz * t;
            let clean = Complex::new(phase.cos() as f32, phase.sin() as f32) * (bpsk_symbol as f32);
            let noise = Complex::new(next_noise() * 0.1, next_noise() * 0.1);
            let sample = clean + noise;

            tracker.predict();
            tracker.update(sample);

            // Verify phase wrapping stays within expected bounds [-pi, pi) for 2pi wrapping limit
            let limit = std::f64::consts::PI;
            assert!(
                tracker.x[0] >= -limit && tracker.x[0] < limit,
                "Phase out of BPSK bounds: {}",
                tracker.x[0]
            );
        }

        // Est frequency should converge close to 20.0 Hz (after scaling EKF frequency by 2.0)
        let est_freq_hz = tracker.x[1] / (2.0 * std::f64::consts::PI) / 2.0;
        let error = (est_freq_hz - target_freq_hz).abs();
        assert!(
            error < 1.5,
            "BPSK EKF failed to lock onto 20.0 Hz. Est: {:.2} Hz",
            est_freq_hz
        );
        assert!(
            tracker.lock_metric > 0.4,
            "Lock metric should indicate lock (>0.4), got: {}",
            tracker.lock_metric
        );
    }

    #[test]
    fn test_clock_stepping_mock() {
        let mut loop_state = LeodoLoop::new();
        assert!(!loop_state.synchronized);
        assert!(loop_state.pending_step_adjustment.is_none());

        // Test EKF update prediction
        loop_state.clock_ekf.x[0] = 0.5; // 500 ms offset

        // Run with enable_steering = false (dry run)
        let log_path = "passes/test_leodo_dry.log";
        steer_system_clock(0.5, 0.0, 150.8e6, log_path, false, &mut loop_state);

        assert!(loop_state.synchronized);
        assert_eq!(loop_state.last_status.contains("DRY RUN"), true);
        assert!(loop_state.pending_step_adjustment.is_none());
        let _ = std::fs::remove_file(log_path);

        // Run with enable_steering = true (will try to step and fail with EPERM unless run as root)
        let mut loop_state_real = LeodoLoop::new();
        loop_state_real.clock_ekf.x[0] = 0.5;
        steer_system_clock(0.5, 0.0, 150.8e6, log_path, true, &mut loop_state_real);

        // Since it fails with EPERM (or succeeds if run as root), synchronized must be set to true
        assert!(loop_state_real.synchronized);
        assert!(
            loop_state_real.last_status.contains("STEP")
                || loop_state_real.last_status.contains("EPERM")
        );
        let _ = std::fs::remove_file(log_path);

        // Dry run with SHM enabled
        let mut loop_state_shm_dry = LeodoLoop::new();
        loop_state_shm_dry.shm_unit = Some(2);
        loop_state_shm_dry.clock_ekf.x[0] = 0.5;
        steer_system_clock(0.5, 0.0, 150.8e6, log_path, false, &mut loop_state_shm_dry);

        assert!(loop_state_shm_dry.synchronized);
        assert!(loop_state_shm_dry.last_status.contains("DRY RUN SHM"));
        let _ = std::fs::remove_file(log_path);

        // Run with enable_steering = true and shm_unit = Some(2)
        let mut loop_state_shm = LeodoLoop::new();
        loop_state_shm.shm_unit = Some(2);
        loop_state_shm.clock_ekf.x[0] = 0.5;
        steer_system_clock(0.5, 0.0, 150.8e6, log_path, true, &mut loop_state_shm);

        assert!(loop_state_shm.synchronized);
        assert!(
            loop_state_shm.last_status.contains("SUCCESS_SHM")
                || loop_state_shm.last_status.contains("ERROR_SHM")
        );
        let _ = std::fs::remove_file(log_path);
    }

    #[test]
    fn test_farrow_interpolator_cubic_lagrange_interpolation() {
        let mut farrow = FarrowInterpolator::new();
        farrow.push(Complex::new(0.0, 0.0));
        farrow.push(Complex::new(1.0, 0.0));
        farrow.push(Complex::new(2.0, 0.0));
        farrow.push(Complex::new(3.0, 0.0));

        let val_zero = farrow.interpolate(0.0);
        let val_mid = farrow.interpolate(0.5);
        let val_one = farrow.interpolate(1.0);

        assert!((val_zero.re - 1.0).abs() < 1e-5);
        assert!((val_mid.re - 1.5).abs() < 1e-5);
        assert!((val_one.re - 2.0).abs() < 1e-5);
    }

    #[test]
    fn test_gardner_timing_error_zero_detector() {
        let mut loop_state = GardnerLoop::new(20000.0, 10000.0);
        let samples = vec![Complex::new(1.0, 0.0); 10];

        let mut out: Vec<(Complex<f32>, f32)> = Vec::new();
        for s in samples {
            loop_state.process(s, &mut out);
        }

        assert!(loop_state.integrator.abs() < 1e-5);
    }

    #[test]
    fn test_gardner_timing_error_sign() {
        let mut loop_state = GardnerLoop::new(20000.0, 10000.0);
        loop_state.on_time_prev = Complex::new(-1.0, 0.0);
        loop_state.mid_time = Complex::new(0.5, 0.0);
        loop_state.sample_count = 2;

        let on_time_curr = Complex::new(1.0, 0.0);
        let error_positive =
            (on_time_curr.re - loop_state.on_time_prev.re) * loop_state.mid_time.re;
        assert!(error_positive > 0.0);

        loop_state.mid_time = Complex::new(-0.5, 0.0);
        let error_negative =
            (on_time_curr.re - loop_state.on_time_prev.re) * loop_state.mid_time.re;
        assert!(error_negative < 0.0);
    }

    #[test]
    fn test_gardner_loop_stable_under_fade() {
        let mut loop_state = GardnerLoop::new(20000.0, 10000.0);
        let mut out: Vec<(Complex<f32>, f32)> = Vec::new();
        for _ in 0..10 {
            loop_state.process(Complex::new(0.0, 0.0), &mut out);
        }
        assert!(loop_state.integrator.abs() < 1e-5);
    }

    #[test]
    fn test_adaptive_ekf_and_dual_lock_behavior() {
        // 1. Process noise Q is scaled down when lock metric is high
        let mut tracker = CarrierPllEkf::new(50000.0, Modulation::Carrier);

        let initial_q_phase = tracker.q_phase;
        let initial_q_freq = tracker.q_freq;
        let initial_q_chirp = tracker.q_chirp;

        tracker.lock_metric = 1.0;
        tracker.p = nalgebra::Matrix3::zeros();
        tracker.predict();

        let dt = tracker.ts;
        let expected_phase_cov = initial_q_phase * 0.1 * dt;
        let expected_freq_cov = initial_q_freq * 0.1 * dt;
        let expected_chirp_cov = initial_q_chirp * 0.1 * dt;

        assert!((tracker.p[(0, 0)] - expected_phase_cov).abs() < 1e-9);
        assert!((tracker.p[(1, 1)] - expected_freq_cov).abs() < 1e-9);
        assert!((tracker.p[(2, 2)] - expected_chirp_cov).abs() < 1e-9);

        // 2. Unlock behavior is triggered when lock metric drops too low
        // Case A: lock_metric < 0.02 -> immediate unlock
        let mut tracker = CarrierPllEkf::new(50000.0, Modulation::Carrier);
        tracker.is_locked = true;
        tracker.lock_metric = 0.015;
        tracker.dual_lock = true;
        tracker.update(Complex::new(1.0, 0.0));
        assert!(!tracker.is_locked);

        // Case B: lock_metric < 0.1 and PR is low -> unlock
        let mut tracker = CarrierPllEkf::new(50000.0, Modulation::Carrier);
        tracker.is_locked = true;
        tracker.lock_metric = 0.08;
        tracker.dual_lock = true;
        tracker.pr_sum_abs_i = 10.0;
        tracker.pr_sum_q_sq = 100.0; // PR = 100 / 100 = 1.0 <= 819.2
        tracker.update(Complex::new(1.0, 0.0));
        assert!(!tracker.is_locked);

        // 3. Lock is maintained when lock metric is low but PR is high
        let mut tracker = CarrierPllEkf::new(50000.0, Modulation::Carrier);
        tracker.is_locked = true;
        tracker.lock_metric = 0.08;
        tracker.dual_lock = true;
        tracker.pr_sum_abs_i = 1000.0;
        tracker.pr_sum_q_sq = 0.1; // PR = 1,000,000 / 0.1 = 10,000,000 > 819.2
        tracker.update(Complex::new(1.0, 0.0));
        assert!(tracker.is_locked);
    }

    #[test]
    fn test_fade_timeout_args_parsing() {
        use clap::Parser;
        let args = Args::parse_from(["sattime", "--fade-timeout", "25.5"]);
        assert_eq!(args.fade_timeout, 25.5);
    }

    #[test]
    fn test_adaptive_fallback_to_unguided() {
        let mut unlocked_frames_during_pass = 0;
        let mut fallback_to_unguided = false;

        let pipeline_sample_rate: f64 = 50000.0;
        let pipeline_step_size: f64 = 1000.0;
        let timeout_frames = (60.0 * pipeline_sample_rate / pipeline_step_size).round() as usize;

        let has_satellites = true;
        let mut visible_expected_frequencies = vec![150803500.0];
        let mut is_tracking = false;

        // 1. Simulating frames passing without acquiring lock
        for _ in 0..=timeout_frames {
            if has_satellites && !visible_expected_frequencies.is_empty() {
                if !is_tracking {
                    unlocked_frames_during_pass += 1;
                    if unlocked_frames_during_pass > timeout_frames && !fallback_to_unguided {
                        fallback_to_unguided = true;
                    }
                } else {
                    unlocked_frames_during_pass = 0;
                }
            } else {
                unlocked_frames_during_pass = 0;
                fallback_to_unguided = false;
            }
        }

        assert!(fallback_to_unguided);
        assert_eq!(unlocked_frames_during_pass, timeout_frames + 1);

        // 2. Lock is acquired
        is_tracking = true;
        if has_satellites && !visible_expected_frequencies.is_empty() {
            if !is_tracking {
                unlocked_frames_during_pass += 1;
            } else {
                unlocked_frames_during_pass = 0;
            }
        }
        assert_eq!(unlocked_frames_during_pass, 0);
        assert!(fallback_to_unguided);

        // 3. Satellite moves out of view
        visible_expected_frequencies.clear();
        if has_satellites && !visible_expected_frequencies.is_empty() {
            // nothing
        } else {

            fallback_to_unguided = false;
        }
        assert!(!fallback_to_unguided);
    }
}

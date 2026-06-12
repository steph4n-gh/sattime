pub mod daemon;
pub mod dsp;
pub mod ekf;
pub mod orbit;
pub mod orbit_solver;
pub mod tui;

use clap::Parser;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CalibrationData {
    pub df0: f64,
    pub timestamp: chrono::DateTime<chrono::Utc>,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "orbital_time_server",
    about = "Track satellite UEMR / NOAA carrier leakage peak from raw IQ stream and synchronize time"
)]
pub struct Args {
    /// Center frequency in Hz (defaults to Starlink VHF)
    #[arg(short = 'f', long = "frequency", default_value_t = 150800000.0)]
    pub frequency: f64,

    /// Sample rate in Hz
    #[arg(short = 's', long = "sample-rate", default_value_t = 2000000.0)]
    pub sample_rate: f64,

    /// FFT size (must be a power of two)
    #[arg(short = 'n', long = "fft-size", default_value_t = 32768)]
    pub fft_size: usize,

    /// Step size (sliding offset in samples, default 40000 for 50 Hz FFT update rate)
    #[arg(short = 'd', long = "step-size", default_value_t = 40000)]
    pub step_size: usize,

    /// Latitude of observer in degrees
    #[arg(long = "lat", default_value_t = 38.889931)]
    pub lat: f64,

    /// Longitude of observer in degrees
    #[arg(long = "lon", default_value_t = -77.009003)]
    pub lon: f64,

    /// Altitude of observer in meters
    #[arg(long = "alt", default_value_t = 25.0)]
    pub alt: f64,

    /// Path to the local TLE file
    #[arg(long = "tle")]
    pub tle: Option<String>,

    /// Bypass automatically downloading the latest TLEs
    #[arg(long = "no-download-tle")]
    pub no_download_tle: bool,

    /// URL to download TLE catalog from
    #[arg(long = "tle-url")]
    pub tle_url: Option<String>,

    /// Optional: SoapySDR driver or query to stream directly from hardware (e.g. "hackrf", "driver=rtlsdr")
    #[arg(long = "sdr")]
    pub sdr: Option<String>,

    /// Optional: General gain value for the SDR receiver
    #[arg(long = "gain")]
    pub gain: Option<f64>,

    /// Optional: HackRF Low-Noise Amplifier (LNA) gain in dB (0-40, default: 24.0)
    #[arg(long = "lna-gain", default_value_t = 24.0)]
    pub lna_gain: f64,

    /// Optional: HackRF RF pre-amplifier gain in dB (0 or 14.0, default: 14.0)
    #[arg(long = "amp-gain", default_value_t = 14.0)]
    pub amp_gain: f64,

    /// Optional: HackRF Variable Gain Amplifier (VGA) gain in dB (0-62, default: 32.0)
    #[arg(long = "vga-gain", default_value_t = 32.0)]
    pub vga_gain: f64,

    /// Optional: Save the downsampled pass data (timestamp, frequency) to this CSV file after capture
    #[arg(long = "save-pass")]
    pub save_pass: Option<String>,

    /// Optional: Solve for 3D location using saved pass CSV files in this directory
    #[arg(long = "solve-location")]
    pub solve_location: Option<String>,

    /// Run in simulation mode (generates mock IQ data for TLE pass)
    #[arg(long = "simulate")]
    pub simulate: bool,

    /// Simulated clock offset in seconds (used in simulation mode)
    #[arg(long = "sim-offset", default_value_t = 5.4)]
    pub sim_offset: f64,

    /// Align the receiver's start time with the simulated pass's epoch (pca_time - 45s)
    #[arg(long = "sim-start-time")]
    pub sim_start_time: bool,

    /// Run as an automated background daemon
    #[arg(long = "daemon")]
    pub daemon: bool,

    /// Solve location blindly without initial coordinate seed (combinatorial search)
    #[arg(long = "blind")]
    pub blind: bool,

    /// Directory to output automatically captured CSV passes in daemon mode
    #[arg(long = "output-dir", default_value = "passes")]
    pub output_dir: String,

    /// Bypass the DC region skip (useful for testing on clean simulated data)
    #[arg(long = "no-dc-skip")]
    pub no_dc_skip: bool,

    /// Minimum SNR in dB to lock and capture signal (default: 8.0)
    #[arg(long = "min-snr", default_value_t = 8.0)]
    pub min_snr: f32,

    /// Decimation factor for the pipeline (e.g. 40 to go from 2 MHz to 50 kHz, optional override)
    #[arg(long = "decimate")]
    pub decimate: Option<usize>,

    /// Disable real-time TLE-guided search windowing (defaults to guided when TLE is available)
    #[arg(long = "no-guided")]
    pub no_guided: bool,

    /// Width of the TLE-guided search window in Hz (optional override)
    #[arg(long = "guided-window")]
    pub guided_window: Option<f64>,

    /// Maximum slant range in meters to consider a satellite overhead (optional override)
    #[arg(long = "max-range")]
    pub max_range: Option<f64>,

    /// Maximum number of closest overhead satellites to track simultaneously (default: 2)
    #[arg(long = "max-guided-sats", default_value_t = 2)]
    pub max_guided_sats: usize,

    /// Disable software Automatic Gain Control (AGC) (on by default)
    #[arg(long = "no-agc")]
    pub no_agc: bool,

    /// Disable dynamic background spur cancellation (notching) (on by default)
    #[arg(long = "no-notch-spurs")]
    pub no_notch_spurs: bool,

    /// Disable Extensive Cancellation Algorithm (ECA) clutter filter on channels (on by default)
    #[arg(long = "no-eca")]
    pub no_eca: bool,

    /// Legacy flag to enable ECA (now enabled by default)
    #[arg(long = "eca", hide = true)]
    pub eca: bool,

    /// Disable real-time visual Terminal UI (TUI) dashboard (on by default)
    #[arg(long = "no-tui")]
    pub no_tui: bool,

    /// Enable active clock steering/disciplining via OS adjtime (LEODO)
    #[arg(long = "leodo")]
    pub leodo: bool,

    /// Log file to write LEODO NTP steering events
    #[arg(long = "leodo-log", default_value = "passes/leodo.log")]
    pub leodo_log: String,

    /// Solve for satellite Keplerian elements from comma-separated pass CSV files (TLE-less orbit determination)
    #[arg(long = "solve-orbit")]
    pub solve_orbit: Option<String>,

    /// Path to write the solved/generated TLE line catalog
    #[arg(long = "output-tle")]
    pub output_tle: Option<String>,

    /// Modulation type to track
    #[arg(long = "modulation", value_enum, default_value_t = dsp::Modulation::Carrier)]
    pub modulation: dsp::Modulation,

    /// Disable adaptive EKF loop bandwidth estimation (on by default)
    #[arg(long = "no-adaptive-ekf")]
    pub no_adaptive_ekf: bool,

    /// Disable dual-stage lock detection (on by default)
    #[arg(long = "no-dual-lock")]
    pub no_dual_lock: bool,

    /// Disable Gardner symbol timing recovery (on by default)
    #[arg(long = "no-gardner")]
    pub no_gardner: bool,

    /// Symbol rate in Hz for symbol timing recovery
    #[arg(long = "symbol-rate")]
    pub symbol_rate: Option<f64>,

    /// Disable multi-hypothesis tracking
    #[arg(long = "no-multihypothesis")]
    pub no_multihypothesis: bool,

    /// Fade timeout in seconds (default: 15.0)
    #[arg(long = "fade-timeout", default_value_t = 15.0)]
    pub fade_timeout: f64,

    /// Maximum number of parallel demodulation channels (default: 8)
    #[arg(long = "max-channels", default_value_t = 8)]
    pub max_channels: usize,
}

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

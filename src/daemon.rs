use crate::ekf::ClockEkf;
use chrono::{DateTime, Utc};
use std::io::Write;
pub static LEODO_LOOP: std::sync::OnceLock<std::sync::Mutex<LeodoLoop>> =
    std::sync::OnceLock::new();

pub struct LeodoLoop {
    pub clock_ekf: ClockEkf,
    pub last_update: Option<chrono::DateTime<chrono::Utc>>,
    pub last_offset: f64,
    pub last_freq_err_ppm: f64,
    pub last_target_adjustment: f64,
    pub last_status: String,
    pub synchronized: bool,
    pub pending_step_adjustment: Option<f64>,
    pub shm_unit: Option<usize>,
}

impl LeodoLoop {
    pub fn new() -> Self {
        Self {
            clock_ekf: ClockEkf::new(),
            last_update: None,
            last_offset: 0.0,
            last_freq_err_ppm: 0.0,
            last_target_adjustment: 0.0,
            last_status: String::from("FREE_RUN"),
            synchronized: false,
            pending_step_adjustment: None,
            shm_unit: None,
        }
    }
}

pub fn get_leodo_loop() -> &'static std::sync::Mutex<LeodoLoop> {
    LEODO_LOOP.get_or_init(|| std::sync::Mutex::new(LeodoLoop::new()))
}

pub fn perform_clock_step(target_adjustment: f64) -> Result<(), std::io::Error> {
    use std::time::SystemTime;
    let now = SystemTime::now();
    let since_the_epoch = now
        .duration_since(SystemTime::UNIX_EPOCH)
        .map_err(std::io::Error::other)?;

    let current_secs = since_the_epoch.as_secs_f64();
    let target_secs = current_secs + target_adjustment;

    #[cfg(target_os = "linux")]
    {
        let ts = libc::timespec {
            tv_sec: target_secs.trunc() as libc::time_t,
            tv_nsec: ((target_secs.fract() * 1_000_000_000.0) as i64) as libc::c_long,
        };
        let ret = unsafe { libc::clock_settime(libc::CLOCK_REALTIME, &ts) };
        if ret == 0 {
            Ok(())
        } else {
            Err(std::io::Error::last_os_error())
        }
    }

    #[cfg(target_os = "macos")]
    {
        let tv = libc::timeval {
            tv_sec: target_secs.trunc() as libc::time_t,
            tv_usec: ((target_secs.fract() * 1_000_000.0) as i32) as libc::suseconds_t,
        };
        let ret = unsafe { libc::settimeofday(&tv, std::ptr::null()) };
        if ret == 0 {
            Ok(())
        } else {
            Err(std::io::Error::last_os_error())
        }
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "Platform not supported for clock stepping",
        ))
    }
}

#[cfg(target_family = "unix")]
#[repr(C)]
struct ShmTime {
    mode: libc::c_int,
    count: libc::c_int,
    clock_time_stamp_sec: libc::time_t,
    clock_time_stamp_usec: libc::c_int,
    receive_time_stamp_sec: libc::time_t,
    receive_time_stamp_usec: libc::c_int,
    leap: libc::c_int,
    precision: libc::c_int,
    nsamples: libc::c_int,
    valid: libc::c_int,
    clock_time_stamp_nsec: libc::c_uint,
    receive_time_stamp_nsec: libc::c_uint,
    dummy: [libc::c_int; 8],
}

#[cfg(target_family = "unix")]
fn write_to_ntp_shm(shm_unit: usize, target_adjustment: f64) -> Result<(), String> {
    let key = 0x4e545030 + shm_unit as i32;
    let size = std::mem::size_of::<ShmTime>();
    let perms = if shm_unit >= 2 { 0o666 } else { 0o600 };
    let shmid = unsafe { libc::shmget(key, size, perms | libc::IPC_CREAT) };
    if shmid < 0 {
        return Err(format!(
            "shmget failed: {} (verify permissions or run as root/sudo for unit < 2)",
            std::io::Error::last_os_error()
        ));
    }
    let shmaddr = unsafe { libc::shmat(shmid, std::ptr::null(), 0) };
    if shmaddr == -1isize as *mut libc::c_void {
        return Err(format!(
            "shmat failed: {}",
            std::io::Error::last_os_error()
        ));
    }
    let shm_ptr = shmaddr as *mut ShmTime;
    let receive_time = chrono::Utc::now();
    let target_time = receive_time + chrono::Duration::microseconds((target_adjustment * 1_000_000.0) as i64);

    let receive_sec = receive_time.timestamp() as libc::time_t;
    let receive_usec = (receive_time.timestamp_subsec_micros() % 1_000_000) as libc::c_int;
    let receive_nsec = receive_time.timestamp_subsec_nanos() as libc::c_uint;

    let clock_sec = target_time.timestamp() as libc::time_t;
    let clock_usec = (target_time.timestamp_subsec_micros() % 1_000_000) as libc::c_int;
    let clock_nsec = target_time.timestamp_subsec_nanos() as libc::c_uint;

    unsafe {
        std::ptr::write_volatile(&mut (*shm_ptr).mode, 1);
        std::ptr::write_volatile(&mut (*shm_ptr).valid, 0);
        std::sync::atomic::compiler_fence(std::sync::atomic::Ordering::SeqCst);

        let count = std::ptr::read_volatile(&(*shm_ptr).count);
        std::ptr::write_volatile(&mut (*shm_ptr).count, count.wrapping_add(1));
        std::sync::atomic::compiler_fence(std::sync::atomic::Ordering::SeqCst);

        std::ptr::write_volatile(&mut (*shm_ptr).clock_time_stamp_sec, clock_sec);
        std::ptr::write_volatile(&mut (*shm_ptr).clock_time_stamp_usec, clock_usec);
        std::ptr::write_volatile(&mut (*shm_ptr).clock_time_stamp_nsec, clock_nsec);
        std::ptr::write_volatile(&mut (*shm_ptr).receive_time_stamp_sec, receive_sec);
        std::ptr::write_volatile(&mut (*shm_ptr).receive_time_stamp_usec, receive_usec);
        std::ptr::write_volatile(&mut (*shm_ptr).receive_time_stamp_nsec, receive_nsec);
        std::ptr::write_volatile(&mut (*shm_ptr).leap, 0);
        std::ptr::write_volatile(&mut (*shm_ptr).precision, -20);
        std::ptr::write_volatile(&mut (*shm_ptr).nsamples, 1);
        std::sync::atomic::compiler_fence(std::sync::atomic::Ordering::SeqCst);

        std::ptr::write_volatile(&mut (*shm_ptr).count, count.wrapping_add(2));
        std::sync::atomic::compiler_fence(std::sync::atomic::Ordering::SeqCst);

        std::ptr::write_volatile(&mut (*shm_ptr).valid, 1);
    }
    if unsafe { libc::shmdt(shmaddr) } < 0 {
        return Err(format!(
            "shmdt failed: {}",
            std::io::Error::last_os_error()
        ));
    }
    Ok(())
}

#[cfg(target_family = "unix")]
pub fn steer_system_clock(
    offset_seconds: f64,
    lo_bias: f64,
    center_freq: f64,
    log_path: &str,
    enable_steering: bool,
    leodo_loop: &mut LeodoLoop,
) -> Vec<String> {
    let now = chrono::Utc::now();
    let mut msgs: Vec<String> = Vec::new();

    // Calculate frequency error in PPM
    let freq_err_ppm = (lo_bias / center_freq) * 1_000_000.0;

    // Update EKF & PI loop
    let dt = if let Some(last) = leodo_loop.last_update {
        (now - last).num_milliseconds() as f64 / 1000.0
    } else {
        0.0
    };

    if dt > 0.0 {
        leodo_loop.clock_ekf.predict(dt);
    }
    leodo_loop.clock_ekf.update(offset_seconds, freq_err_ppm);

    // Compute the target adjustment from the EKF phase offset
    let target_adjustment = leodo_loop.clock_ekf.x[0];

    let status_str;
    let mut actual_slewed = 0.0;

    if enable_steering {
        if let Some(shm_unit) = leodo_loop.shm_unit {
            match write_to_ntp_shm(shm_unit, target_adjustment) {
                Ok(_) => {
                    status_str = format!("SUCCESS_SHM (NTP{})", shm_unit);
                    msgs.push(format!(
                        "[LEODO] Successfully wrote clock offset of {:.6}s to NTP SHM segment NTP{}",
                        target_adjustment, shm_unit
                    ));
                    leodo_loop.synchronized = true;
                    leodo_loop.last_update = Some(now);
                }
                Err(err) => {
                    status_str = format!("ERROR_SHM: {}", err);
                    msgs.push(format!(
                        "[LEODO] Failed to write to NTP SHM segment NTP{}: {}",
                        shm_unit, err
                    ));
                }
            }
        } else if !leodo_loop.synchronized && target_adjustment.abs() > 0.1 {
            // Step-once on the first synchronization event if offset > 100 ms
            match perform_clock_step(target_adjustment) {
                Ok(_) => {
                    actual_slewed = target_adjustment;
                    leodo_loop.synchronized = true;
                    leodo_loop.pending_step_adjustment = Some(target_adjustment);

                    // Reset last_update to prevent bad dt transition
                    leodo_loop.last_update = None;

                    status_str = format!("SUCCESS_STEP (stepped {:.6}s)", target_adjustment);
                    msgs.push(format!(
                        "[LEODO] Successfully stepped system clock by {:.6}s",
                        target_adjustment
                    ));
                }
                Err(err) => {
                    leodo_loop.synchronized = true; // Lock stepping even if it failed/EPERM so subsequent steering uses slewing
                    let target_epoch = (now
                        + chrono::Duration::microseconds((target_adjustment * 1_000_000.0) as i64))
                    .timestamp();
                    #[cfg(target_os = "macos")]
                    let override_cmd = format!("sudo date -f \"%s\" \"{}\"", target_epoch);
                    #[cfg(not(target_os = "macos"))]
                    let override_cmd = format!("sudo date -s \"@{}\"", target_epoch);

                    if err.raw_os_error() == Some(libc::EPERM) {
                        status_str =
                            "ERROR_STEP EPERM (permission denied, run as root/sudo)".to_string();
                        msgs.push(format!(
                            "[LEODO] Clock stepping failed: Permission denied. Run as sudo, or: {}",
                            override_cmd
                        ));
                    } else {
                        status_str = format!("ERROR_STEP: {}", err);
                        msgs.push(format!("[LEODO] Clock stepping failed: {}", err));
                    }
                }
            }
        } else {
            // Subsequent adjustment or small offset: gradual slewing
            leodo_loop.synchronized = true;
            leodo_loop.last_update = Some(now);

            // Prepare the timeval struct for libc::adjtime
            let sec = target_adjustment.trunc() as libc::time_t;
            let usec = ((target_adjustment.fract() * 1_000_000.0) as i32) as libc::suseconds_t;

            let delta = libc::timeval {
                tv_sec: sec,
                tv_usec: usec,
            };

            let mut old_delta = libc::timeval {
                tv_sec: 0,
                tv_usec: 0,
            };

            // unsafe block to call the native OS API
            let ret = unsafe { libc::adjtime(&delta, &mut old_delta) };
            if ret == 0 {
                actual_slewed = 0.0; // Audit Fix S5/S6: adjtime is gradual, no instant step occurred
                status_str = format!("SUCCESS_SLEW (target {:.6}s)", target_adjustment);
                msgs.push(format!(
                    "[LEODO] Successfully requested OS clock slew of {:.6}s",
                    target_adjustment
                ));
            } else {
                let err = std::io::Error::last_os_error();
                let target_epoch = (now
                    + chrono::Duration::microseconds((target_adjustment * 1_000_000.0) as i64))
                .timestamp();
                #[cfg(target_os = "macos")]
                let override_cmd = format!("sudo date -f \"%s\" \"{}\"", target_epoch);
                #[cfg(not(target_os = "macos"))]
                let override_cmd = format!("sudo date -s \"@{}\"", target_epoch);

                if err.raw_os_error() == Some(libc::EPERM) {
                    status_str =
                        "ERROR_SLEW EPERM (permission denied, run as root/sudo)".to_string();
                    msgs.push(format!(
                        "[LEODO] Clock slewing failed: Permission denied. Run as sudo, or: {}",
                        override_cmd
                    ));
                } else {
                    status_str = format!("ERROR_SLEW: {}", err);
                    msgs.push(format!("[LEODO] Clock slewing failed: {}", err));
                }
            }
        }
    } else {
        if let Some(shm_unit) = leodo_loop.shm_unit {
            status_str = format!("DRY RUN SHM (NTP{})", shm_unit);
            msgs.push(format!(
                "[LEODO] Dry run (SHM NTP{}): Time offset {:.6}s, EKF phase offset {:.6}s, drift {:.3} PPM",
                shm_unit, offset_seconds, target_adjustment, leodo_loop.clock_ekf.x[1]
            ));
        } else {
            status_str = format!("DRY RUN (calculated {:.6}s adjustment)", target_adjustment);
            msgs.push(format!(
                "[LEODO] Dry run: Time offset {:.6}s, EKF phase offset {:.6}s, drift {:.3} PPM",
                offset_seconds, target_adjustment, leodo_loop.clock_ekf.x[1]
            ));
        }

        leodo_loop.synchronized = true;
        leodo_loop.last_update = Some(now);
    }

    // Apply control feedback step correction to EKF phase state
    leodo_loop.clock_ekf.x[0] -= actual_slewed;

    leodo_loop.last_offset = offset_seconds;
    leodo_loop.last_freq_err_ppm = leodo_loop.clock_ekf.x[1];
    leodo_loop.last_target_adjustment = target_adjustment;
    leodo_loop.last_status = status_str.clone();

    // Write to the log file
    if let Some(parent) = std::path::Path::new(log_path).parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    if let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)
    {
        let log_line = format!(
            "{},{},{:.6},{:.3},{:.6},{:.6},{}\n",
            now.to_rfc3339(),
            center_freq,
            offset_seconds,
            leodo_loop.clock_ekf.x[1],
            target_adjustment,
            actual_slewed,
            status_str
        );
        let _ = file.write_all(log_line.as_bytes());
    }
    msgs
}

#[cfg(not(target_family = "unix"))]
pub fn steer_system_clock(
    offset_seconds: f64,
    lo_bias: f64,
    center_freq: f64,
    log_path: &str,
    enable_steering: bool,
    leodo_loop: &mut LeodoLoop,
) -> Vec<String> {
    let now = chrono::Utc::now();
    let dt = if let Some(last) = leodo_loop.last_update {
        (now - last).num_milliseconds() as f64 / 1000.0
    } else {
        0.0
    };
    leodo_loop.last_update = Some(now);

    let freq_err_ppm = (lo_bias / center_freq) * 1_000_000.0;

    if dt > 0.0 {
        leodo_loop.clock_ekf.predict(dt);
    }
    leodo_loop.clock_ekf.update(offset_seconds, freq_err_ppm);

    let target_adjustment = leodo_loop.clock_ekf.x[0];
    leodo_loop.synchronized = true;
    leodo_loop.last_offset = offset_seconds;
    leodo_loop.last_freq_err_ppm = leodo_loop.clock_ekf.x[1];
    leodo_loop.last_target_adjustment = target_adjustment;

    let mut msgs = Vec::new();
    if leodo_loop.shm_unit.is_some() {
        msgs.push("[LEODO] NTP SHM steering is not supported on this platform.".to_string());
        leodo_loop.last_status = "NOT_SUPPORTED (SHM)".to_string();
    } else {
        msgs.push("[LEODO] Clock steering is not supported on this platform.".to_string());
        leodo_loop.last_status = format!("NOT_SUPPORTED (dry EKF target: {:.6}s)", target_adjustment);
    }
    msgs
}

pub enum DaemonState {
    Searching,
    Capturing {
        start_time: DateTime<Utc>,
        samples: Vec<(DateTime<Utc>, f64)>,
        last_lock_time: DateTime<Utc>,
    },
}

pub struct CompletedPassData {
    pub sat_name: String,
    pub timestamp: DateTime<Utc>,
    pub offset_seconds: f64,
    pub freq_drift_ppm: f64,
    pub snr: f64,
    pub max_elevation: f64,
    pub fit_rmse: f64,
}

pub struct ConsensusSteeringEngine {
    pub passes: Vec<CompletedPassData>,
}

impl ConsensusSteeringEngine {
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    pub fn add_pass_result(&mut self, pass: CompletedPassData) {
        self.passes.push(pass);
    }

    pub fn get_consensus_update(&self) -> Option<(f64, f64)> {
        if self.passes.is_empty() {
            return None;
        }

        let mut total_weight = 0.0;
        let mut weighted_offset = 0.0;
        let mut weighted_drift = 0.0;

        for pass in &self.passes {
            // Outlier rejection: SNR must be >= 3.0, fit_rmse must be <= 100.0
            if pass.fit_rmse > 100.0 || pass.snr < 3.0 {
                continue;
            }

            // Weight calculation based on SNR, elevation, and fit RMSE
            let snr_weight = (pass.snr - 3.0).max(0.0);
            let elev_weight = (pass.max_elevation.to_radians()).sin().max(0.0);
            let rmse_weight = 1.0 / (pass.fit_rmse.max(0.1));

            let weight = snr_weight * elev_weight * rmse_weight;

            if weight > 0.0 {
                total_weight += weight;
                weighted_offset += pass.offset_seconds * weight;
                weighted_drift += pass.freq_drift_ppm * weight;
            }
        }

        if total_weight > 0.0 {
            Some((
                weighted_offset / total_weight,
                weighted_drift / total_weight,
            ))
        } else {
            None
        }
    }
}

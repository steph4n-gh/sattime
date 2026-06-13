use crate::daemon::*;
use crate::orbit_solver;
use chrono::{DateTime, Datelike, Timelike, Utc};
use std::io::{self, Read};
pub static GEOLOCATION_RESULT: std::sync::OnceLock<std::sync::Mutex<GeolocationResult>> =
    std::sync::OnceLock::new();

pub struct GeolocationResult {
    pub lat: f64,
    pub lon: f64,
    pub alt: f64,
    pub rmse: f64,
    pub converged: bool,
    pub num_passes: usize,
    pub gdop: f64,
    pub uncertainty_km: f64,
}

pub fn get_geolocation_result() -> &'static std::sync::Mutex<GeolocationResult> {
    GEOLOCATION_RESULT.get_or_init(|| {
        std::sync::Mutex::new(GeolocationResult {
            lat: 0.0,
            lon: 0.0,
            alt: 0.0,
            rmse: 0.0,
            converged: false,
            num_passes: 0,
            gdop: 0.0,
            uncertainty_km: 0.0,
        })
    })
}

pub fn wgs84_to_ecef(lat_deg: f64, lon_deg: f64, alt_m: f64) -> [f64; 3] {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let a = 6378137.0; // WGS84 semi-major axis in meters
    let f = 1.0 / 298.257223563; // flattening
    let e2 = 2.0 * f - f * f; // eccentricity squared

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();

    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();

    let x = (n + alt_m) * cos_lat * lon.cos();
    let y = (n + alt_m) * cos_lat * lon.sin();
    let z = (n * (1.0 - e2) + alt_m) * sin_lat;

    [x, y, z]
}

pub fn ecef_to_wgs84(ecef: [f64; 3]) -> (f64, f64, f64) {
    let x = ecef[0];
    let y = ecef[1];
    let z = ecef[2];

    let a = 6378137.0; // semi-major axis
    let f = 1.0 / 298.257223563; // flattening
    let b = a * (1.0 - f); // semi-minor axis
    let e2 = (a * a - b * b) / (a * a); // first eccentricity squared
    let ep2 = (a * a - b * b) / (b * b); // second eccentricity squared

    let p = (x * x + y * y).sqrt();
    let theta = (z * a).atan2(p * b);

    let lat = (z + ep2 * b * theta.sin().powi(3)).atan2(p - e2 * a * theta.cos().powi(3));
    let lon = y.atan2(x);

    let sin_lat = lat.sin();
    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let alt = p / lat.cos() - n;

    (lat.to_degrees(), lon.to_degrees(), alt)
}

pub fn ecef_to_enu(sat_ecef: [f64; 3], obs_ecef: [f64; 3]) -> [f64; 3] {
    let (lat_deg, lon_deg, _) = ecef_to_wgs84(obs_ecef);
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();

    let dx = sat_ecef[0] - obs_ecef[0];
    let dy = sat_ecef[1] - obs_ecef[1];
    let dz = sat_ecef[2] - obs_ecef[2];

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_lon = lon.sin();
    let cos_lon = lon.cos();

    let east = -sin_lon * dx + cos_lon * dy;
    let north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    let up = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;

    [east, north, up]
}

pub fn enu_to_az_el(enu: [f64; 3]) -> (f64, f64) {
    let east = enu[0];
    let north = enu[1];
    let up = enu[2];
    let az = east.atan2(north).rem_euclid(2.0 * std::f64::consts::PI);
    let el = up.atan2((east * east + north * north).sqrt());
    (az, el)
}

pub fn apply_sagnac_correction(
    pos_sat: [f64; 3],
    vel_sat: [f64; 3],
    pos_obs: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let dx = pos_sat[0] - pos_obs[0];
    let dy = pos_sat[1] - pos_obs[1];
    let dz = pos_sat[2] - pos_obs[2];
    let range = (dx * dx + dy * dy + dz * dz).sqrt();
    if range > 0.0 {
        let tau = range / 299792458.0;
        let omega_e = 7.2921151467e-5;
        let theta_sagnac = -omega_e * tau;
        let cos_t = theta_sagnac.cos();
        let sin_t = theta_sagnac.sin();
        let p_corr = [
            pos_sat[0] * cos_t + pos_sat[1] * sin_t,
            -pos_sat[0] * sin_t + pos_sat[1] * cos_t,
            pos_sat[2],
        ];
        let v_corr = [
            vel_sat[0] * cos_t + vel_sat[1] * sin_t,
            -vel_sat[0] * sin_t + vel_sat[1] * cos_t,
            vel_sat[2],
        ];
        (p_corr, v_corr)
    } else {
        (pos_sat, vel_sat)
    }
}

pub fn saastamoinen_tropospheric_delay(sat_ecef: [f64; 3], obs_ecef: [f64; 3]) -> f64 {
    let enu = ecef_to_enu(sat_ecef, obs_ecef);
    let (_, el) = enu_to_az_el(enu);
    let sin_el = el.sin().max(0.01);
    2.3 / (sin_el + 0.00143)
}

pub fn solve_linear_system(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for i in 0..n {
        // Pivot selection
        let mut max_row = i;
        for k in (i + 1)..n {
            if a[k][i].abs() > a[max_row][i].abs() {
                max_row = k;
            }
        }
        if a[max_row][i].abs() < 1e-12 {
            return None; // Singular matrix
        }
        a.swap(i, max_row);
        b.swap(i, max_row);

        // Eliminate
        for k in (i + 1)..n {
            let factor = a[k][i] / a[i][i];
            b[k] -= factor * b[i];
            for j in i..n {
                a[k][j] -= factor * a[i][j];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..n {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    Some(x)
}

pub fn predict_freq_sample(
    orbit: &OrbitModel,
    pos_obs: [f64; 3],
    center_freq: f64,
    dt: DateTime<Utc>,
    delta_t: f64,
    df0: f64,
) -> Option<f64> {
    let dt_true = dt - chrono::Duration::microseconds((delta_t * 1e6) as i64);
    if let Some((pos_sat, vel_sat)) = orbit.propagate_ecef(dt_true) {
        let (pos_sat_corr, vel_sat_corr) = apply_sagnac_correction(pos_sat, vel_sat, pos_obs);
        let rx = pos_sat_corr[0] - pos_obs[0];
        let ry = pos_sat_corr[1] - pos_obs[1];
        let rz = pos_sat_corr[2] - pos_obs[2];
        let range = (rx * rx + ry * ry + rz * rz).sqrt();
        if range > 0.0 {
            let range_rate = (rx * vel_sat_corr[0] + ry * vel_sat_corr[1] + rz * vel_sat_corr[2]) / range;
            let doppler_term = 1.0 - range_rate / 299792458.0;
            return Some(df0 + center_freq * doppler_term);
        }
    }
    None
}

pub struct PassDataRef<'a> {
    #[allow(dead_code)]
    pub name: &'a String,
    pub orbit: &'a OrbitModel,
    pub data: &'a Vec<(DateTime<Utc>, f64)>,
    pub center_freq: f64,
}

pub fn estimate_slant_range(f0: f64, min_slope: f64) -> f64 {
    let c = 299792458.0;
    let v_sat = 7585.0; // average Starlink orbital velocity
    let slope_abs = min_slope.abs();
    if slope_abs > 0.1 {
        (f0 * v_sat * v_sat) / (c * slope_abs)
    } else {
        2000000.0 // fallback 2000 km
    }
}

pub fn run_fast_gn_fit(
    passes: &[PassDataRef],
    init_x: [f64; 3],
    _center_freq: f64,
    ground_alt: f64,
) -> Option<(Vec<f64>, f64, f64, f64)> {
    let n_passes = passes.len();
    let (mut lat, mut lon, _) = ecef_to_wgs84(init_x);

    let mut theta = vec![0.0; 3 + 2 * n_passes];

    let max_iter = 60;
    let mut rmse = f64::MAX;
    let mut final_a_mat = vec![vec![0.0; 2]; 2];

    let mut x_hist = Vec::new();
    let mut y_hist = Vec::new();
    let mut z_hist = Vec::new();

    for iter in 0..max_iter {
        let current_x_obs = wgs84_to_ecef(lat, lon, ground_alt);
        x_hist.push(current_x_obs[0] / 6378137.0);
        y_hist.push(current_x_obs[1] / 6378137.0);
        z_hist.push(current_x_obs[2] / 6378137.0);

        if x_hist.len() > 5 {
            x_hist.remove(0);
        }
        if y_hist.len() > 5 {
            y_hist.remove(0);
        }
        if z_hist.len() > 5 {
            z_hist.remove(0);
        }

        let x_diffs = orbit_solver::compute_fractional_difference_history(&x_hist);
        let y_diffs = orbit_solver::compute_fractional_difference_history(&y_hist);
        let z_diffs = orbit_solver::compute_fractional_difference_history(&z_hist);

        let var_sum = x_diffs.iter().map(|d| d.abs()).sum::<f64>()
            + y_diffs.iter().map(|d| d.abs()).sum::<f64>()
            + z_diffs.iter().map(|d| d.abs()).sum::<f64>();

        let _ = orbit_solver::map_to_normalized_search_space(current_x_obs);

        // 1. Optimize clock offset (dt) and frequency bias (df0) for each satellite individually
        for j in 0..n_passes {
            let pass = &passes[j];
            let pca_meas = match estimate_measured_pca_time(pass.data) {
                Some(t) => t,
                None => pass.data[pass.data.len() / 2].0,
            };
            if let Some((dt_opt, df0_opt, _)) = fit_satellite(
                pass.orbit,
                pass.data,
                current_x_obs,
                pass.center_freq,
                pca_meas,
            ) {
                theta[3 + 2 * j] = dt_opt;
                theta[3 + 2 * j + 1] = df0_opt;
            }
        }

        // Enforce shared receiver clock offset across all passes to prevent coordinate overfitting
        let mut sum_dt = 0.0;
        for j in 0..n_passes {
            sum_dt += theta[3 + 2 * j];
        }
        let avg_dt = sum_dt / n_passes as f64;
        for j in 0..n_passes {
            theta[3 + 2 * j] = avg_dt;
        }

        // 2. Optimize physical receiver coordinates (lat, lon) with clock parameters fixed
        let mut r = Vec::new();
        let mut j_mat = Vec::new();

        for j in 0..n_passes {
            let pass = &passes[j];
            let delta_t = theta[3 + 2 * j];
            let df0 = theta[3 + 2 * j + 1];

            for &(dt, freq_meas) in pass.data {
                if let Some(f_pred) = predict_freq_sample(
                    pass.orbit,
                    current_x_obs,
                    pass.center_freq,
                    dt,
                    delta_t,
                    df0,
                ) {
                    let residual = freq_meas - f_pred;
                    r.push(residual);

                    let mut row = vec![0.0; 2];

                    // Finite difference for latitude (in degrees)
                    let h_lat = 1e-5; // ~1.1 meters
                    let pos_lat_plus = wgs84_to_ecef(lat + h_lat, lon, ground_alt);
                    let pos_lat_minus = wgs84_to_ecef(lat - h_lat, lon, ground_alt);
                    let f_lat_plus = predict_freq_sample(
                        pass.orbit,
                        pos_lat_plus,
                        pass.center_freq,
                        dt,
                        delta_t,
                        df0,
                    )
                    .unwrap_or(f_pred);
                    let f_lat_minus = predict_freq_sample(
                        pass.orbit,
                        pos_lat_minus,
                        pass.center_freq,
                        dt,
                        delta_t,
                        df0,
                    )
                    .unwrap_or(f_pred);
                    row[0] = (f_lat_plus - f_lat_minus) / (2.0 * h_lat);

                    // Finite difference for longitude (in degrees)
                    let h_lon = 1e-5; // ~1.1 cos(lat) meters
                    let pos_lon_plus = wgs84_to_ecef(lat, lon + h_lon, ground_alt);
                    let pos_lon_minus = wgs84_to_ecef(lat, lon - h_lon, ground_alt);
                    let f_lon_plus = predict_freq_sample(
                        pass.orbit,
                        pos_lon_plus,
                        pass.center_freq,
                        dt,
                        delta_t,
                        df0,
                    )
                    .unwrap_or(f_pred);
                    let f_lon_minus = predict_freq_sample(
                        pass.orbit,
                        pos_lon_minus,
                        pass.center_freq,
                        dt,
                        delta_t,
                        df0,
                    )
                    .unwrap_or(f_pred);
                    row[1] = (f_lon_plus - f_lon_minus) / (2.0 * h_lon);

                    j_mat.push(row);
                }
            }
        }

        let m = r.len();
        let k_params = 2; // only lat, lon
        if m < k_params {
            return None;
        }

        let mut a_mat = vec![vec![0.0; k_params]; k_params];
        let mut b_vec = vec![0.0; k_params];

        for p in 0..k_params {
            for q in 0..k_params {
                let mut sum = 0.0;
                for row_idx in 0..m {
                    sum += j_mat[row_idx][p] * j_mat[row_idx][q];
                }
                a_mat[p][q] = sum;
            }
            let mut sum_b = 0.0;
            for row_idx in 0..m {
                sum_b += j_mat[row_idx][p] * r[row_idx];
            }
            b_vec[p] = sum_b; // correct descent sign
        }

        final_a_mat = a_mat.clone();

        // Add absolute regularization (Tikhonov style) to prevent singularity
        let base_lambda = 0.01;
        let lambda = base_lambda * (1.0 + var_sum * 0.1);
        for p in 0..k_params {
            a_mat[p][p] += lambda * a_mat[p][p] + 1e-4;
        }

        let delta_theta = solve_linear_system(a_mat.clone(), b_vec)?;

        let mut d_lat = delta_theta[0];
        let mut d_lon = delta_theta[1];

        // Audit Fix S8: Limit step size to at most 166 km (Euclidean) to prevent polar distortion
        let lat_rad = lat.to_radians();
        let d_lat_m = d_lat * 111320.0;
        let d_lon_m = d_lon * 111320.0 * lat_rad.cos().abs();
        let step_meters = (d_lat_m * d_lat_m + d_lon_m * d_lon_m).sqrt();

        if step_meters > 166000.0 {
            let scale = 166000.0 / step_meters;
            d_lat *= scale;
            d_lon *= scale;
        }

        lat += d_lat;
        lon += d_lon;

        let step_meters_taken = step_meters.min(166000.0);
        rmse = (r.iter().map(|&val| val * val).sum::<f64>() / m as f64).sqrt();

        tracing::debug!(
            "[DEBUG_GN_ITER] Iter {} | step: {:.1} m | rmse: {:.2} Hz | Geodetic: [{:.6}, {:.6}]",
            iter,
            step_meters_taken,
            rmse,
            lat,
            lon
        );

        if step_meters_taken < 0.2 {
            break;
        }
    }

    let final_ecef = wgs84_to_ecef(lat, lon, ground_alt);
    let mut theta_out = vec![0.0; 3 + 2 * n_passes];
    theta_out[0] = final_ecef[0];
    theta_out[1] = final_ecef[1];
    theta_out[2] = final_ecef[2];
    for j in 0..n_passes {
        theta_out[3 + 2 * j] = theta[3 + 2 * j];
        theta_out[3 + 2 * j + 1] = theta[3 + 2 * j + 1];
    }

    // Calculate final covariance matrix from final_a_mat (J^T J)
    let mut gdop = 99.9;
    let mut uncertainty_km = 999.9;
    let det = final_a_mat[0][0] * final_a_mat[1][1] - final_a_mat[0][1] * final_a_mat[1][0];
    if det.abs() > 1e-12 {
        let inv_00 = final_a_mat[1][1] / det;
        let inv_11 = final_a_mat[0][0] / det;
        if inv_00 > 0.0 && inv_11 > 0.0 {
            gdop = (inv_00 + inv_11).sqrt().clamp(0.1, 99.9);
            let sigma_lat = (rmse * rmse * inv_00).sqrt();
            let sigma_lon = (rmse * rmse * inv_11).sqrt();
            let err_lat_km = sigma_lat * 111.0;
            let err_lon_km = sigma_lon * 111.0 * lat.to_radians().cos().abs();
            uncertainty_km = (err_lat_km * err_lat_km + err_lon_km * err_lon_km)
                .sqrt()
                .clamp(0.01, 999.9);
        }
    }

    Some((theta_out, rmse, gdop, uncertainty_km))
}

pub fn intersect_circles(
    p1: [f64; 3],
    d1: f64,
    p2: [f64; 3],
    d2: f64,
) -> Option<([f64; 3], [f64; 3])> {
    let r_e = 6378137.0;

    let n1 = (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]).sqrt();
    let n2 = (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]).sqrt();
    if n1 < 1.0 || n2 < 1.0 {
        return None;
    }

    let u1 = [p1[0] / n1, p1[1] / n1, p1[2] / n1];
    let u2 = [p2[0] / n2, p2[1] / n2, p2[2] / n2];

    let cos_theta1 = ((r_e * r_e + n1 * n1 - d1 * d1) / (2.0 * r_e * n1)).clamp(-1.0, 1.0);
    let cos_theta2 = ((r_e * r_e + n2 * n2 - d2 * d2) / (2.0 * r_e * n2)).clamp(-1.0, 1.0);

    let c1 = r_e * cos_theta1;
    let c2 = r_e * cos_theta2;

    let dot = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
    let denom = 1.0 - dot * dot;
    if denom < 1e-6 {
        return None;
    }

    let a = (c1 - c2 * dot) / denom;
    let b = (c2 - c1 * dot) / denom;

    let x0 = [
        a * u1[0] + b * u2[0],
        a * u1[1] + b * u2[1],
        a * u1[2] + b * u2[2],
    ];

    let x0_norm2 = x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2];
    let v_len2 = (r_e * r_e - x0_norm2).max(0.0);

    let v = [
        u1[1] * u2[2] - u1[2] * u2[1],
        u1[2] * u2[0] - u1[0] * u2[2],
        u1[0] * u2[1] - u1[1] * u2[0],
    ];
    let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if v_norm < 1e-6 {
        return None;
    }

    let t = (v_len2 / (v_norm * v_norm)).sqrt();

    let xa = [x0[0] + t * v[0], x0[1] + t * v[1], x0[2] + t * v[2]];
    let xb = [x0[0] - t * v[0], x0[1] - t * v[1], x0[2] - t * v[2]];

    Some((xa, xb))
}

pub fn run_location_solver(
    dir_path: &str,
    tle_path: &str,
    default_center_freq: f64,
    initial_guess_lat_lon_alt: [f64; 3],
    blind: bool,
    enable_leodo: bool,
    leodo_log_path: &str,
    suspend_steering: bool,
) {
    use std::fs;
    use std::io::BufRead;

    tracing::debug!("\n=== Starting Passive 3D Ground Geolocation Solver ===");
    tracing::debug!("Loading satellite elements from: {}", tle_path);

    let satellites = match load_orbits(tle_path) {
        Ok(sats) => sats,
        Err(e) => {
            eprintln!("Error parsing ephemeris file (tried both SP3 and TLE): {}", e);
            return;
        }
    };
    tracing::debug!("Loaded {} satellites.", satellites.len());

    // 1. Read files
    let mut paths = Vec::new();
    let metadata = match fs::metadata(dir_path) {
        Ok(m) => m,
        Err(e) => {
            tracing::warn!("Error reading path '{}': {}", dir_path, e);
            return;
        }
    };

    if metadata.is_dir() {
        if let Ok(entries) = fs::read_dir(dir_path) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() && path.extension().is_some_and(|ext| ext == "csv") {
                    paths.push(path);
                }
            }
        }
    } else {
        paths.push(std::path::PathBuf::from(dir_path));
    }

    if paths.is_empty() {
        tracing::warn!("Error: No CSV pass files found in '{}'", dir_path);
        return;
    }

    pub struct RawPass {
        pub filename: String,
        pub sat_name: String,
        pub center_freq: f64,
        pub data: Vec<(DateTime<Utc>, f64)>,
        pub pca_time: DateTime<Utc>,
        #[allow(dead_code)]
        pub min_slope: f64,
        pub d_min: f64,
    }

    let mut raw_passes = Vec::new();

    for path in paths {
        let file = match fs::File::open(&path) {
            Ok(f) => f,
            Err(e) => {
                tracing::warn!("Warning: Failed to open file {:?}: {}", path, e);
                continue;
            }
        };

        let filename = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let reader = io::BufReader::new(file);
        let mut sat_name = "UNKNOWN".to_string();
        let mut center_freq = default_center_freq;
        let mut data = Vec::new();

        for line in reader.lines().flatten() {
            let trimmed = line.trim();
            if trimmed.starts_with("# satellite_name:") {
                sat_name = trimmed["# satellite_name:".len()..].trim().to_string();
            } else if trimmed.starts_with("# center_frequency:") {
                if let Ok(freq) = trimmed["# center_frequency:".len()..].trim().parse::<f64>() {
                    center_freq = freq;
                }
            } else if !trimmed.starts_with('#')
                && !trimmed.is_empty()
                && !trimmed.starts_with("utc_timestamp")
            {
                let parts: Vec<&str> = trimmed.split(',').collect();
                if parts.len() == 2
                    && let (Ok(dt), Ok(freq)) = (
                        DateTime::parse_from_rfc3339(parts[0]),
                        parts[1].parse::<f64>(),
                    )
                {
                    data.push((dt.with_timezone(&Utc), freq));
                }
            }
        }

        if data.is_empty() {
            tracing::warn!("Warning: No data points found in {:?}", path);
            continue;
        }

        // Estimate PCA time
        let pca_time = match estimate_measured_pca_time(&data) {
            Some(t) => t,
            None => data[data.len() / 2].0,
        };

        // Estimate min slope
        let mut min_slope = 0.0;
        if data.len() >= 3 {
            let mut min_deriv = f64::MAX;
            for i in 1..(data.len() - 1) {
                let dt_prev = data[i - 1].0;
                let dt_next = data[i + 1].0;
                let f_prev = data[i - 1].1;
                let f_next = data[i + 1].1;
                let dt_diff = (dt_next - dt_prev).num_milliseconds() as f64 / 1000.0;
                if dt_diff > 0.0 {
                    let deriv = (f_next - f_prev) / dt_diff;
                    // Ignore physically impossible step changes (e.g. tracking gaps near DC center)
                    if deriv < min_deriv && deriv >= -120.0 {
                        min_deriv = deriv;
                    }
                }
            }
            if min_deriv < 0.0 {
                min_slope = min_deriv;
            }
        }
        let d_min = estimate_slant_range(center_freq, min_slope);

        raw_passes.push(RawPass {
            filename,
            sat_name,
            center_freq,
            data,
            pca_time,
            min_slope,
            d_min,
        });
    }

    let n_passes = raw_passes.len();
    tracing::debug!("Loaded {} valid pass files.", n_passes);

    if n_passes < 3 {
        tracing::warn!(
            "Error: Geolocation requires at least 3 independent passes to solve for 3D coordinates. Found {}.",
            n_passes
        );
        return;
    }

    let pos_ref = wgs84_to_ecef(
        initial_guess_lat_lon_alt[0],
        initial_guess_lat_lon_alt[1],
        initial_guess_lat_lon_alt[2],
    );
    let mut solved_x = pos_ref;

    // Check if we need to run blind satellite matching.
    let needs_blind = blind || raw_passes.iter().any(|p| p.sat_name == "UNKNOWN");

    pub struct PassData {
        pub filename: String,
        pub sat_name: String,
        pub center_freq: f64,
        pub data: Vec<(DateTime<Utc>, f64)>,
        pub orbit: OrbitModel,
    }

    let mut passes = Vec::new();

    if needs_blind {
        tracing::debug!(
            "Needs TLE candidate matching (either blind=true or UNKNOWN satellites found)."
        );
        tracing::debug!(
            "Reference center: ECEF=[{:.1}, {:.1}, {:.1}] m (Geodetic=[{:.6}, {:.6}, {:.1}] m)",
            pos_ref[0],
            pos_ref[1],
            pos_ref[2],
            initial_guess_lat_lon_alt[0],
            initial_guess_lat_lon_alt[1],
            initial_guess_lat_lon_alt[2]
        );

        #[derive(Clone)]
        pub struct BlindCandidate {
            pub name: String,
            pub orbit: OrbitModel,
            pub pos_sat: [f64; 3],
            pub d_min: f64,
        }

        // Find candidates for the first 3 passes
        let mut passes_candidates = Vec::new();
        for j in 0..3 {
            let pass = &raw_passes[j];
            let mut candidates = Vec::new();

            for (name, orbit) in &satellites {
                if let Some((pos_sat, _)) = orbit.propagate_ecef(pass.pca_time) {
                    // Distance to reference coordinate
                    let dx = pos_sat[0] - pos_ref[0];
                    let dy = pos_sat[1] - pos_ref[1];
                    let dz = pos_sat[2] - pos_ref[2];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                    if name == "STARLINK-1063"
                        || name == "STARLINK-1265"
                        || name == "STARLINK-1477"
                        || name == "STARLINK-1008"
                    {
                        tracing::debug!(
                            "    [DEBUG] Satellite {} | dist to ref: {:.1} km | pass.d_min: {:.1} km | diff: {:.1} km",
                            name,
                            dist / 1000.0,
                            pass.d_min / 1000.0,
                            (dist - pass.d_min).abs() / 1000.0
                        );
                    }

                    // Filter: must be within 4000 km of our region to be a candidate,
                    // AND the distance must match the estimated Doppler slant range within a 1200 km regional margin
                    if dist < 4000000.0 && (dist - pass.d_min).abs() < 1200000.0 {
                        candidates.push(BlindCandidate {
                            name: name.clone(),
                            orbit: orbit.clone(),
                            pos_sat,
                            d_min: pass.d_min,
                        });
                    }
                }
            }
            tracing::debug!(
                "  Pass {} (file {}): Found {} TLE candidates within 4000 km of region.",
                j + 1,
                pass.filename,
                candidates.len()
            );
            if candidates.is_empty() {
                tracing::warn!(
                    "Error: No candidates found for pass {}. Try adjusting your coordinates seed or check TLE catalog.",
                    j + 1
                );
                return;
            }
            passes_candidates.push(candidates);
        }

        tracing::debug!(
            "Searching combinatorial candidate space for the correct satellite matches..."
        );
        let mut best_combo = None;
        let mut best_rmse = f64::MAX;
        let mut best_theta = Vec::new();

        'outer: for c0 in &passes_candidates[0] {
            for c1 in &passes_candidates[1] {
                // Prune pairs: intersect circles of c0 and c1 on Earth's surface
                let intersections =
                    match intersect_circles(c0.pos_sat, c0.d_min, c1.pos_sat, c1.d_min) {
                        Some(pt) => pt,
                        None => {
                            continue;
                        }
                    };

                for c2 in &passes_candidates[2] {
                    // Check if c2 is compatible with either intersection point (within 80 km tolerance)
                    let dist_a2 = (c2.pos_sat[0] - intersections.0[0]).powi(2)
                        + (c2.pos_sat[1] - intersections.0[1]).powi(2)
                        + (c2.pos_sat[2] - intersections.0[2]).powi(2);
                    let dist_b2 = (c2.pos_sat[0] - intersections.1[0]).powi(2)
                        + (c2.pos_sat[1] - intersections.1[1]).powi(2)
                        + (c2.pos_sat[2] - intersections.1[2]).powi(2);

                    let err_a = (dist_a2.sqrt() - c2.d_min).abs();
                    let err_b = (dist_b2.sqrt() - c2.d_min).abs();

                    if err_a > 500000.0 && err_b > 500000.0 {
                        continue; // Geometric mismatch, prune!
                    }

                    // Geometrically compatible combination! Run the fast GN solver to verify and refine
                    let test_passes = vec![
                        PassDataRef {
                            name: &c0.name,
                            orbit: &c0.orbit,
                            data: &raw_passes[0].data,
                            center_freq: raw_passes[0].center_freq,
                        },
                        PassDataRef {
                            name: &c1.name,
                            orbit: &c1.orbit,
                            data: &raw_passes[1].data,
                            center_freq: raw_passes[1].center_freq,
                        },
                        PassDataRef {
                            name: &c2.name,
                            orbit: &c2.orbit,
                            data: &raw_passes[2].data,
                            center_freq: raw_passes[2].center_freq,
                        },
                    ];

                    let init_x = if err_a < err_b {
                        intersections.0
                    } else {
                        intersections.1
                    };
                    if let Some((theta, rmse, _, _)) = run_fast_gn_fit(
                        &test_passes,
                        init_x,
                        default_center_freq,
                        initial_guess_lat_lon_alt[2],
                    ) {
                        if rmse < 15.0 {
                            let (_, _, alt) = ecef_to_wgs84([theta[0], theta[1], theta[2]]);
                            if alt > -500.0 && alt < 9000.0 && rmse < best_rmse {
                                best_rmse = rmse;
                                best_theta = theta;
                                best_combo = Some((c0.clone(), c1.clone(), c2.clone()));
                                if rmse < 5.0 {
                                    break 'outer;
                                }
                            }
                        }
                    }
                }
            }
        }

        if let Some((c0, c1, c2)) = best_combo {
            tracing::debug!("\nSUCCESSFULLY IDENTIFIED MATCHING SATELLITES:");
            tracing::debug!("  Pass 1: {}", c0.name);
            tracing::debug!("  Pass 2: {}", c1.name);
            tracing::debug!("  Pass 3: {}", c2.name);
            tracing::debug!("Combinatorial Fit RMSE: {:.2} Hz", best_rmse);

            if best_theta.len() >= 3 {
                solved_x = [best_theta[0], best_theta[1], best_theta[2]];
            }

            passes.push(PassData {
                filename: raw_passes[0].filename.clone(),
                sat_name: c0.name,
                center_freq: raw_passes[0].center_freq,
                data: raw_passes[0].data.clone(),
                orbit: c0.orbit.clone(),
            });
            passes.push(PassData {
                filename: raw_passes[1].filename.clone(),
                sat_name: c1.name,
                center_freq: raw_passes[1].center_freq,
                data: raw_passes[1].data.clone(),
                orbit: c1.orbit.clone(),
            });
            passes.push(PassData {
                filename: raw_passes[2].filename.clone(),
                sat_name: c2.name,
                center_freq: raw_passes[2].center_freq,
                data: raw_passes[2].data.clone(),
                orbit: c2.orbit.clone(),
            });
        } else {
            tracing::warn!(
                "Error: Could not identify any matching combination of satellites with low fitting error."
            );
            return;
        }
    } else {
        // Standard mode: resolve satellites from CSV headers
        for raw in raw_passes {
            let found_orb = satellites.iter().find(|(name, _)| name == &raw.sat_name);
            let (name, orbit) = match found_orb {
                Some(t) => t,
                None => {
                    let match_sub = satellites.iter().find(|(name, _)| {
                        name.to_uppercase().contains(&raw.sat_name.to_uppercase())
                    });
                    match match_sub {
                        Some(t) => t,
                        None => {
                            tracing::warn!(
                                "Warning: Satellite '{}' from file {:?} not found in database. Skipping.",
                                raw.sat_name,
                                raw.filename
                            );
                            continue;
                        }
                    }
                }
            };

            passes.push(PassData {
                filename: raw.filename,
                sat_name: name.clone(),
                center_freq: raw.center_freq,
                data: raw.data,
                orbit: orbit.clone(),
            });
        }
    }

    let n_passes = passes.len();
    if n_passes < 3 {
        tracing::warn!(
            "Error: Geolocation requires at least 3 valid passes. Found {}.",
            n_passes
        );
        return;
    }

    // 2. Perform final geodetic joint optimization using run_fast_gn_fit
    let test_passes: Vec<PassDataRef> = passes
        .iter()
        .map(|p| PassDataRef {
            name: &p.sat_name,
            orbit: &p.orbit,
            data: &p.data,
            center_freq: p.center_freq,
        })
        .collect();

    let (final_theta, final_rmse, gdop, uncertainty_km) = match run_fast_gn_fit(
        &test_passes,
        solved_x,
        default_center_freq,
        initial_guess_lat_lon_alt[2],
    ) {
        Some(res) => res,
        None => {
            tracing::warn!("Error: Final geodetic solver failed to converge.");
            return;
        }
    };

    let final_ecef = [final_theta[0], final_theta[1], final_theta[2]];
    let (lat_sol, lon_sol, alt_sol) = ecef_to_wgs84(final_ecef);

    {
        let mut res = get_geolocation_result().lock().unwrap();
        res.lat = lat_sol;
        res.lon = lon_sol;
        res.alt = alt_sol;
        res.rmse = final_rmse;
        res.converged = true;
        res.num_passes = n_passes;
        res.gdop = gdop;
        res.uncertainty_km = uncertainty_km;
    }

    tracing::debug!("\n=======================================================");
    tracing::debug!("GEOLOCATION SOLVER RESULTS (CONVERGED):");
    tracing::debug!("Resolved Latitude:        {:.7}°", lat_sol);
    tracing::debug!("Resolved Longitude:       {:.7}°", lon_sol);
    tracing::debug!("Resolved Altitude:        {:.1} m", alt_sol);
    tracing::debug!(
        "Resolved ECEF:            [{:.1}, {:.1}, {:.1}] m",
        final_ecef[0],
        final_ecef[1],
        final_ecef[2]
    );
    tracing::debug!("Overall Fit Quality:      {:.2} Hz RMSE", final_rmse);
    tracing::debug!("Geometric DOP (GDOP):     {:.2}", gdop);
    tracing::debug!("Estimated Uncertainty:    {:.3} km", uncertainty_km);
    tracing::debug!("=======================================================");

    tracing::debug!("\nPASS DIAGNOSTICS:");
    for j in 0..n_passes {
        let pass = &passes[j];
        let dt = final_theta[3 + 2 * j];
        let df = final_theta[3 + 2 * j + 1];
        tracing::debug!(
            "  Pass {:<2} ({:<20}) | Offset (dt): {:+8.3}s | LO Bias (df0): {:+8.2} Hz | File: {}",
            j + 1,
            pass.sat_name,
            dt,
            df,
            pass.filename
        );
    }
    tracing::debug!("=======================================================");

    // Perform LEODO Clock Steering using the latest pass diagnostics
    if n_passes > 0 {
        let latest_idx = n_passes - 1;
        let latest_dt = final_theta[3 + 2 * latest_idx];
        let latest_df = final_theta[3 + 2 * latest_idx + 1];

        if final_rmse < 50.0 {
            let cal = crate::CalibrationData {
                df0: latest_df,
                timestamp: chrono::Utc::now(),
            };
            if let Ok(json) = serde_json::to_string_pretty(&cal) {
                let _ = std::fs::write("calibration.json", json);
            }
        }

        let steer_enabled = if suspend_steering {
            tracing::warn!("NTP clock steering suspended");
            tracing::warn!(
                "[LEODO] Pass terminated in a fade. Suspending clock discipline updates."
            );
            false
        } else {
            enable_leodo
        };

        let mut loop_lock = get_leodo_loop().lock().unwrap();
        steer_system_clock(
            latest_dt,
            latest_df,
            default_center_freq,
            leodo_log_path,
            steer_enabled,
            &mut loop_lock,
        );
    }
}

pub fn run_blind_solver_check(
    dir_path: &str,
    tle_path: &str,
    default_center_freq: f64,
    initial_guess_lat_lon_alt: [f64; 3],
    blind: bool,
    enable_leodo: bool,
    leodo_log_path: &str,
    suspend_steering: bool,
) {
    if let Ok(entries) = std::fs::read_dir(dir_path) {
        let mut count = 0;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() && path.extension().is_some_and(|ext| ext == "csv") {
                count += 1;
            }
        }
        if count >= 3 {
            tracing::debug!(
                "\n[DAEMON] Found {} pass files. Running location solver...",
                count
            );
            run_location_solver(
                dir_path,
                tle_path,
                default_center_freq,
                initial_guess_lat_lon_alt,
                blind,
                enable_leodo,
                leodo_log_path,
                suspend_steering,
            );
        } else {
            tracing::debug!(
                "\n[DAEMON] Found {} pass files. Need at least 3 passes to run geolocation.",
                count
            );
        }
    }
}

pub fn datetime_to_jd(dt: DateTime<Utc>) -> (f64, f64) {
    let year = dt.year() as f64;
    let month = dt.month() as f64;
    let day = dt.day() as f64;
    let hour = dt.hour() as f64;
    let minute = dt.minute() as f64;
    let second = dt.second() as f64;
    let nanosecond = dt.nanosecond() as f64;

    let day_fraction = (hour + (minute + (second + nanosecond / 1e9) / 60.0) / 60.0) / 24.0;

    let (y, m) = if month <= 2.0 {
        (year - 1.0, month + 12.0)
    } else {
        (year, month)
    };

    let a = (y / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();

    let jd_base = (365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + day + b - 1524.5;
    (jd_base, day_fraction)
}

pub fn teme_to_ecef(jd: (f64, f64), pos_teme: [f64; 3], vel_teme: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let d = (jd.0 - 2451545.0) + jd.1;
    let t = d / 36525.0;
    let mut gmst =
        280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t * t * t / 38710000.0;
    gmst %= 360.0;
    if gmst < 0.0 {
        gmst += 360.0;
    }
    let theta = gmst.to_radians();
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    // Position rotation (Greenwich Sidereal Time rotation around Earth Z axis)
    let x_ecef = pos_teme[0] * cos_t + pos_teme[1] * sin_t;
    let y_ecef = -pos_teme[0] * sin_t + pos_teme[1] * cos_t;
    let z_ecef = pos_teme[2];

    // Earth rotation speed in rad/s
    let omega_e = 7.2921151467e-5;

    // Velocity rotation including coriolis term: v_ecef = R_z(theta) * v_teme - omega x r_ecef
    let vx_ecef = (vel_teme[0] * cos_t + vel_teme[1] * sin_t) + omega_e * y_ecef;
    let vy_ecef = (-vel_teme[0] * sin_t + vel_teme[1] * cos_t) - omega_e * x_ecef;
    let vz_ecef = vel_teme[2];

    ([x_ecef, y_ecef, z_ecef], [vx_ecef, vy_ecef, vz_ecef])
}

pub fn download_tle_file(
    urls: &[&str],
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let temp_path = format!("{}.tmp", output_path);
    for &url in urls {
        tracing::warn!("Downloading latest TLE catalog from {}...", url);
        let resp = match ureq::get(url)
            .timeout(std::time::Duration::from_secs(3))
            .set(
                "User-Agent",
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/120.0.0.0 Safari/537.36",
            )
            .call()
        {
            Ok(r) => r,
            Err(e) => {
                eprintln!("Warning: Failed to download TLE from {}: {}", url, e);
                continue;
            }
        };

        if resp.status() == 200 {
            let mut body = String::new();
            if let Err(e) = resp.into_reader().read_to_string(&mut body) {
                eprintln!("Warning: Failed to read response body from {}: {}", url, e);
                continue;
            }
            if let Err(e) = std::fs::write(&temp_path, &body) {
                eprintln!("Warning: Failed to write temp TLE file: {}", e);
                continue;
            }
            // Parse and validate
            match parse_tle_file(&temp_path) {
                Ok(sats) => {
                    if sats.is_empty() {
                        eprintln!("Warning: Parsed TLE from {} is empty", url);
                        let _ = std::fs::remove_file(&temp_path);
                        continue;
                    }
                    // Valid! Overwrite the cache file
                    if let Err(e) = std::fs::rename(&temp_path, output_path) {
                        eprintln!("Warning: Failed to rename temp file to cache path: {}", e);
                        let _ = std::fs::remove_file(&temp_path);
                        continue;
                    }
                    eprintln!(
                        "Successfully downloaded, validated, and cached TLE catalog to '{}' from {}",
                        output_path, url
                    );
                    return Ok(());
                }
                Err(e) => {
                    eprintln!(
                        "Warning: Failed to parse downloaded TLE from {}: {}",
                        url, e
                    );
                    let _ = std::fs::remove_file(&temp_path);
                    continue;
                }
            }
        } else {
            eprintln!("Warning: HTTP error status {} from {}", resp.status(), url);
        }
    }

    if std::path::Path::new(output_path).exists() {
        eprintln!("[WARNING] All mirrors failed. Falling back to cached TLE.");
        Ok(())
    } else {
        Err("All mirrors failed and no cached TLE found".into())
    }
}

pub fn get_tle_file_cached(
    urls: &[&str],
    path: &str,
    force_download: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut needs_download = force_download || !std::path::Path::new(path).exists();
    if !needs_download
        && let Ok(metadata) = std::fs::metadata(path)
        && let Ok(modified) = metadata.modified()
        && let Ok(elapsed) = modified.elapsed()
        && elapsed.as_secs() > 43200
    {
        // 12 hours
        needs_download = true;
    }

    if needs_download && let Err(e) = download_tle_file(urls, path) {
        eprintln!(
            "[WARNING] TLE download failed: {}. Falling back to cached file.",
            e
        );
        if !std::path::Path::new(path).exists() {
            return Err(e); // No cache exists, must fail
        }
    }
    Ok(())
}

pub fn parse_tle_file(path: &str) -> io::Result<Vec<(String, sgp4::Elements)>> {
    use std::fs::File;
    use std::io::BufRead;

    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    let mut satellites = Vec::new();
    let mut i = 0;

    // Helper to fix malformed exponent signs in TLE line 1 (e.g. missing '+' sign replaced by space)
    let sanitize_line1 = |line: &str| -> String {
        let mut bytes = line.as_bytes().to_vec();
        if bytes.len() >= 68 {
            if bytes[50] == b' ' {
                bytes[50] = b'+';
            }
            if bytes[59] == b' ' {
                bytes[59] = b'+';
            }
        }
        String::from_utf8_lossy(&bytes).into_owned()
    };

    while i < lines.len() {
        let line = lines[i].trim();
        if line.is_empty() {
            i += 1;
            continue;
        }

        // Handle 2-line format
        if line.starts_with('1') && i + 1 < lines.len() && lines[i + 1].trim().starts_with('2') {
            let line1 = sanitize_line1(line);
            let line2 = lines[i + 1].trim();
            match sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes()) {
                Ok(elements) => {
                    let sat_name = elements
                        .object_name
                        .clone()
                        .unwrap_or_else(|| "UNKNOWN".to_string());
                    satellites.push((sat_name, elements));
                }
                Err(e) => {
                    eprintln!("Error parsing 2-line TLE at line {}: {:?}", i + 1, e);
                }
            }
            i += 2;
        } else if i + 2 < lines.len() {
            let name = line.to_string();
            let line1 = sanitize_line1(lines[i + 1].trim());
            let line2 = lines[i + 2].trim();

            if line1.starts_with('1') && line2.starts_with('2') {
                match sgp4::Elements::from_tle(
                    Some(name.clone()),
                    line1.as_bytes(),
                    line2.as_bytes(),
                ) {
                    Ok(elements) => {
                        satellites.push((name, elements));
                    }
                    Err(e) => {
                        eprintln!("Error parsing TLE for '{}': {:?}", name, e);
                    }
                }
                i += 3;
            } else {
                i += 1;
            }
        } else {
            i += 1;
        }
    }

    Ok(satellites)
}

pub fn downsample_data(data: &[(DateTime<Utc>, f64)]) -> Vec<(DateTime<Utc>, f64)> {
    if data.is_empty() {
        return Vec::new();
    }
    let mut bins: std::collections::HashMap<i64, Vec<(DateTime<Utc>, f64)>> =
        std::collections::HashMap::new();
    for &(dt, freq) in data {
        let sec = dt.timestamp();
        bins.entry(sec).or_default().push((dt, freq));
    }

    let mut result = Vec::new();
    for (_, samples) in bins {
        let count = samples.len();
        if count == 0 {
            continue;
        }
        let base_dt = samples[0].0;
        let mut sum_us = 0i64;
        let mut sum_freq = 0.0f64;
        for &(dt, freq) in &samples {
            let offset = (dt - base_dt).num_microseconds().unwrap_or(0);
            sum_us += offset;
            sum_freq += freq;
        }
        let avg_offset_us = sum_us / count as i64;
        let avg_dt = base_dt + chrono::Duration::microseconds(avg_offset_us);
        let avg_freq = sum_freq / count as f64;
        result.push((avg_dt, avg_freq));
    }
    result.sort_by_key(|&(dt, _)| dt);
    result
}

pub fn estimate_measured_pca_time(data: &[(DateTime<Utc>, f64)]) -> Option<DateTime<Utc>> {
    if data.len() < 3 {
        return None;
    }
    let mut min_deriv = f64::MAX;
    let mut best_idx = 0;
    for i in 1..(data.len() - 1) {
        let dt_prev = data[i - 1].0;
        let dt_next = data[i + 1].0;
        let f_prev = data[i - 1].1;
        let f_next = data[i + 1].1;

        let dt_diff = (dt_next - dt_prev).num_milliseconds() as f64 / 1000.0;
        if dt_diff > 0.0 {
            let deriv = (f_next - f_prev) / dt_diff;
            if deriv < min_deriv {
                min_deriv = deriv;
                best_idx = i;
            }
        }
    }
    if min_deriv < 0.0 {
        Some(data[best_idx].0)
    } else {
        None
    }
}

pub fn find_pca_time(
    orbit: &OrbitModel,
    pos_obs: [f64; 3],
    around_time: DateTime<Utc>,
) -> Option<DateTime<Utc>> {
    let mut min_range = f64::MAX;
    let mut best_secs = 0.0;

    // 1. Grid search in minutes from -50.0 to +50.0 relative to around_time
    let steps = 100;
    for step in 0..=steps {
        let mins = -50.0 + (step as f64);
        let dt = around_time + chrono::Duration::microseconds((mins * 60.0 * 1e6) as i64);
        if let Some((pos_sat, _)) = orbit.propagate_ecef(dt) {
            let rx = pos_sat[0] - pos_obs[0];
            let ry = pos_sat[1] - pos_obs[1];
            let rz = pos_sat[2] - pos_obs[2];
            let range = rx * rx + ry * ry + rz * rz;
            if range < min_range {
                min_range = range;
                best_secs = mins * 60.0;
            }
        }
    }

    // 2. Refine in seconds from best_secs - 60 to best_secs + 60
    min_range = f64::MAX;
    let mut refined_secs = best_secs;
    for s in -60..=60 {
        let secs = best_secs + (s as f64);
        let dt = around_time + chrono::Duration::microseconds((secs * 1e6) as i64);
        if let Some((pos_sat, _)) = orbit.propagate_ecef(dt) {
            let rx = pos_sat[0] - pos_obs[0];
            let ry = pos_sat[1] - pos_obs[1];
            let rz = pos_sat[2] - pos_obs[2];
            let range = rx * rx + ry * ry + rz * rz;
            if range < min_range {
                min_range = range;
                refined_secs = secs;
            }
        }
    }

    Some(around_time + chrono::Duration::microseconds((refined_secs * 1e6) as i64))
}

pub fn fit_satellite(
    orbit: &OrbitModel,
    data: &[(DateTime<Utc>, f64)],
    pos_obs: [f64; 3],
    center_freq: f64,
    measured_pca_time: DateTime<Utc>,
) -> Option<(f64, f64, f64)> {
    // 1. Find predicted PCA time for this satellite
    let predicted_pca_time = find_pca_time(orbit, pos_obs, measured_pca_time)?;

    // 2. Estimate initial delta_t (difference between local capture PCA and orbital predicted PCA)
    let est_delta_t = (measured_pca_time - predicted_pca_time).num_milliseconds() as f64 / 1000.0;

    let mut min_rmse = f64::MAX;
    let mut best_delta_t = 0.0;
    let mut best_freq_offset = 0.0;

    let mut centers = vec![est_delta_t];
    if est_delta_t.abs() > 5.0 {
        centers.push(0.0);
    }

    for center in centers {
        let steps = 300;
        for step in 0..=steps {
            let delta_t = center - 15.0 + (step as f64) * 0.1;
            let mut y = Vec::with_capacity(data.len());

            for &(dt, freq_meas) in data {
                let dt_true = dt - chrono::Duration::microseconds((delta_t * 1e6) as i64);
                if let Some((pos_sat, vel_sat)) = orbit.propagate_ecef(dt_true) {
                    let (pos_sat_corr, vel_sat_corr) = apply_sagnac_correction(pos_sat, vel_sat, pos_obs);
                    let rx = pos_sat_corr[0] - pos_obs[0];
                    let ry = pos_sat_corr[1] - pos_obs[1];
                    let rz = pos_sat_corr[2] - pos_obs[2];
                    let range = (rx * rx + ry * ry + rz * rz).sqrt();
                    let range_rate = (rx * vel_sat_corr[0] + ry * vel_sat_corr[1] + rz * vel_sat_corr[2]) / range;

                    let doppler_term = 1.0 - range_rate / 299792458.0;
                    y.push(freq_meas - center_freq * doppler_term);
                }
            }

            if y.len() > data.len() / 2 {
                let sum: f64 = y.iter().sum();
                let count = y.len() as f64;
                let df0 = sum / count;

                let sq_sum: f64 = y
                    .iter()
                    .map(|&val| {
                        let diff = val - df0;
                        diff * diff
                    })
                    .sum();

                let rmse = (sq_sum / count).sqrt();
                if rmse < min_rmse {
                    min_rmse = rmse;
                    best_delta_t = delta_t;
                    best_freq_offset = df0;
                }
            }
        }
    }

    // Refine around the best delta_t with 1 ms steps
    let mut refined_min_rmse = min_rmse;
    let mut refined_best_delta_t = best_delta_t;
    let mut refined_best_freq_offset = best_freq_offset;

    if min_rmse < f64::MAX {
        let fine_steps = 100;
        let start_t = best_delta_t - 0.05;
        for step in 0..=fine_steps {
            let delta_t = start_t + (step as f64) * 0.001;
            let mut y = Vec::with_capacity(data.len());

            for &(dt, freq_meas) in data {
                let dt_true = dt - chrono::Duration::microseconds((delta_t * 1e6) as i64);
                if let Some((pos_sat, vel_sat)) = orbit.propagate_ecef(dt_true) {
                    let (pos_sat_corr, vel_sat_corr) = apply_sagnac_correction(pos_sat, vel_sat, pos_obs);
                    let rx = pos_sat_corr[0] - pos_obs[0];
                    let ry = pos_sat_corr[1] - pos_obs[1];
                    let rz = pos_sat_corr[2] - pos_obs[2];
                    let range = (rx * rx + ry * ry + rz * rz).sqrt();
                    let range_rate = (rx * vel_sat_corr[0] + ry * vel_sat_corr[1] + rz * vel_sat_corr[2]) / range;
                    let doppler_term = 1.0 - range_rate / 299792458.0;

                    y.push(freq_meas - center_freq * doppler_term);
                }
            }

            if y.len() > data.len() / 2 {
                let sum: f64 = y.iter().sum();
                let count = y.len() as f64;
                let df0 = sum / count;
                let sq_sum: f64 = y
                    .iter()
                    .map(|&val| {
                        let diff = val - df0;
                        diff * diff
                    })
                    .sum();
                let rmse = (sq_sum / count).sqrt();

                if rmse < refined_min_rmse {
                    refined_min_rmse = rmse;
                    refined_best_delta_t = delta_t;
                    refined_best_freq_offset = df0;
                }
            }
        }
    }

    // Microsecond-level refinement (10 microsecond steps over a +/- 1 ms window)
    let mut micro_min_rmse = refined_min_rmse;
    let mut micro_best_delta_t = refined_best_delta_t;
    let mut micro_best_freq_offset = refined_best_freq_offset;

    if refined_min_rmse < f64::MAX {
        let micro_steps = 200;
        let start_t = refined_best_delta_t - 0.001;
        for step in 0..=micro_steps {
            let delta_t = start_t + (step as f64) * 0.00001;
            let mut y = Vec::with_capacity(data.len());

            for &(dt, freq_meas) in data {
                let dt_true = dt - chrono::Duration::microseconds((delta_t * 1e6) as i64);
                if let Some((pos_sat, vel_sat)) = orbit.propagate_ecef(dt_true) {
                    let (pos_sat_corr, vel_sat_corr) = apply_sagnac_correction(pos_sat, vel_sat, pos_obs);
                    let rx = pos_sat_corr[0] - pos_obs[0];
                    let ry = pos_sat_corr[1] - pos_obs[1];
                    let rz = pos_sat_corr[2] - pos_obs[2];
                    let range = (rx * rx + ry * ry + rz * rz).sqrt();
                    let range_rate = (rx * vel_sat_corr[0] + ry * vel_sat_corr[1] + rz * vel_sat_corr[2]) / range;
                    let doppler_term = 1.0 - range_rate / 299792458.0;

                    y.push(freq_meas - center_freq * doppler_term);
                }
            }

            if y.len() > data.len() / 2 {
                let sum: f64 = y.iter().sum();
                let count = y.len() as f64;
                let df0 = sum / count;
                let sq_sum: f64 = y
                    .iter()
                    .map(|&val| {
                        let diff = val - df0;
                        diff * diff
                    })
                    .sum();
                let rmse = (sq_sum / count).sqrt();

                if rmse < micro_min_rmse {
                    micro_min_rmse = rmse;
                    micro_best_delta_t = delta_t;
                    micro_best_freq_offset = df0;
                }
            }
        }
    }

    if micro_min_rmse < f64::MAX {
        Some((
            micro_best_delta_t,
            micro_best_freq_offset,
            micro_min_rmse,
        ))
    } else {
        None
    }
}

pub fn save_pass_data(
    path: &str,
    sat_name: &str,
    center_freq: f64,
    data: &[(DateTime<Utc>, f64)],
) -> io::Result<()> {
    use std::fs::File;
    use std::io::Write;

    if let Some(parent) = std::path::Path::new(path).parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut file = File::create(path)?;
    writeln!(file, "# satellite_name: {}", sat_name)?;
    writeln!(file, "# center_frequency: {}", center_freq)?;
    writeln!(file, "utc_timestamp,measured_frequency")?;

    for &(dt, freq) in data {
        writeln!(file, "{},{}", dt.to_rfc3339(), freq)?;
    }

    tracing::debug!("Saved pass data ({}) to '{}'", sat_name, path);
    Ok(())
}

pub fn haversine_distance(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let r = 6371.0; // Earth's radius in km
    let d_lat = (lat2 - lat1).to_radians();
    let d_lon = (lon2 - lon1).to_radians();
    let a = (d_lat / 2.0).sin().powi(2)
        + lat1.to_radians().cos() * lat2.to_radians().cos() * (d_lon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
    r * c
}

// --- Real-Time Geodetic Geolocation Structs ---

#[derive(Debug, Clone, Copy)]
pub struct ECEFCoordinates {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Velocity {
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct GeodeticCoordinates {
    pub latitude: f64,
    pub longitude: f64,
    pub altitude: f64,
}

/// Real-time 3D geodetic solver using range-based Gauss-Newton trilateration.
/// Requires ≥4 simultaneous satellite range measurements to solve for 3D position.
pub struct RealTimeGeoSolver {
    pub initial_guess: ECEFCoordinates,
    pub apply_troposphere: bool,
}

impl RealTimeGeoSolver {
    pub fn new() -> Self {
        Self {
            // Default initial guess: on the equator at the prime meridian on Earth's surface
            initial_guess: ECEFCoordinates {
                x: 0.0,
                y: 0.0,
                z: 6378137.0, // WGS84 semi-major axis
            },
            apply_troposphere: false,
        }
    }

    /// Given ≥4 measurements of (satellite_ecef, satellite_velocity, slant_range_m),
    /// solve for the receiver's geodetic position using Gauss-Newton iteration.
    /// Returns None if fewer than 4 measurements or if the solver fails to converge.
    pub fn update_position(
        &mut self,
        measurements: &[(ECEFCoordinates, Velocity, f64)],
    ) -> Option<GeodeticCoordinates> {
        let n = measurements.len();
        if n < 4 {
            return None;
        }

        // Validate inputs — reject NaN/Inf and degenerate ranges
        for (pos, vel, range) in measurements {
            if pos.x.is_nan()
                || pos.y.is_nan()
                || pos.z.is_nan()
                || vel.vx.is_nan()
                || vel.vy.is_nan()
                || vel.vz.is_nan()
                || range.is_nan()
                || range.is_infinite()
                || pos.x.is_infinite()
                || pos.y.is_infinite()
                || pos.z.is_infinite()
            {
                return None;
            }
        }

        // Gauss-Newton trilateration
        let mut x = self.initial_guess.x;
        let mut y = self.initial_guess.y;
        let mut z = self.initial_guess.z;

        let mut x_hist = Vec::new();
        let mut y_hist = Vec::new();
        let mut z_hist = Vec::new();

        for _iter in 0..100 {
            x_hist.push(x / 6378137.0);
            y_hist.push(y / 6378137.0);
            z_hist.push(z / 6378137.0);

            if x_hist.len() > 5 {
                x_hist.remove(0);
            }
            if y_hist.len() > 5 {
                y_hist.remove(0);
            }
            if z_hist.len() > 5 {
                z_hist.remove(0);
            }

            let x_diffs = orbit_solver::compute_fractional_difference_history(&x_hist);
            let y_diffs = orbit_solver::compute_fractional_difference_history(&y_hist);
            let z_diffs = orbit_solver::compute_fractional_difference_history(&z_hist);

            let var_sum = x_diffs.iter().map(|d| d.abs()).sum::<f64>()
                + y_diffs.iter().map(|d| d.abs()).sum::<f64>()
                + z_diffs.iter().map(|d| d.abs()).sum::<f64>();

            let _ = orbit_solver::map_to_normalized_search_space([x, y, z]);
            let mut jtj = [[0.0f64; 3]; 3];
            let mut jtr = [0.0f64; 3];

            for (sat, _vel, range) in measurements {
                let dx = x - sat.x;
                let dy = y - sat.y;
                let dz = z - sat.z;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < 1e-10 {
                    continue;
                }

                let delta_trop = if self.apply_troposphere {
                    saastamoinen_tropospheric_delay([sat.x, sat.y, sat.z], [x, y, z])
                } else {
                    0.0
                };
                let residual = dist - (range - delta_trop);
                let jx = dx / dist;
                let jy = dy / dist;
                let jz = dz / dist;

                // Accumulate J^T * J
                jtj[0][0] += jx * jx;
                jtj[0][1] += jx * jy;
                jtj[0][2] += jx * jz;
                jtj[1][0] += jy * jx;
                jtj[1][1] += jy * jy;
                jtj[1][2] += jy * jz;
                jtj[2][0] += jz * jx;
                jtj[2][1] += jz * jy;
                jtj[2][2] += jz * jz;

                // Accumulate J^T * r
                jtr[0] += jx * residual;
                jtr[1] += jy * residual;
                jtr[2] += jz * residual;
            }

            // Check for singular/degenerate geometry (coplanar satellites, etc.) on the undamped matrix
            let det_undamped = jtj[0][0] * (jtj[1][1] * jtj[2][2] - jtj[1][2] * jtj[2][1])
                - jtj[0][1] * (jtj[1][0] * jtj[2][2] - jtj[1][2] * jtj[2][0])
                + jtj[0][2] * (jtj[1][0] * jtj[2][1] - jtj[1][1] * jtj[2][0]);
            let trace_undamped = jtj[0][0] + jtj[1][1] + jtj[2][2];
            if det_undamped.abs() < 1e-30
                || det_undamped.abs() < trace_undamped * trace_undamped * trace_undamped * 1e-12
            {
                return None;
            }

            // Introduce Tikhonov damping on the diagonal of jtj to suppress oscillations
            let base_lambda = if var_sum > 1.0 { 0.01 } else { 0.0 };
            let lambda = base_lambda * (1.0 + var_sum * 0.1);
            if lambda > 0.0 {
                for p in 0..3 {
                    jtj[p][p] += lambda * jtj[p][p] + 1e-4;
                }
            }

            // Solve 3x3 system using Cramer's rule on the damped matrix
            let det = jtj[0][0] * (jtj[1][1] * jtj[2][2] - jtj[1][2] * jtj[2][1])
                - jtj[0][1] * (jtj[1][0] * jtj[2][2] - jtj[1][2] * jtj[2][0])
                + jtj[0][2] * (jtj[1][0] * jtj[2][1] - jtj[1][1] * jtj[2][0]);

            if det.abs() < 1e-30 {
                return None;
            }

            let inv_det = 1.0 / det;

            // Cofactor matrix / det
            let inv = [
                [
                    (jtj[1][1] * jtj[2][2] - jtj[1][2] * jtj[2][1]) * inv_det,
                    (jtj[0][2] * jtj[2][1] - jtj[0][1] * jtj[2][2]) * inv_det,
                    (jtj[0][1] * jtj[1][2] - jtj[0][2] * jtj[1][1]) * inv_det,
                ],
                [
                    (jtj[1][2] * jtj[2][0] - jtj[1][0] * jtj[2][2]) * inv_det,
                    (jtj[0][0] * jtj[2][2] - jtj[0][2] * jtj[2][0]) * inv_det,
                    (jtj[0][2] * jtj[1][0] - jtj[0][0] * jtj[1][2]) * inv_det,
                ],
                [
                    (jtj[1][0] * jtj[2][1] - jtj[1][1] * jtj[2][0]) * inv_det,
                    (jtj[0][1] * jtj[2][0] - jtj[0][0] * jtj[2][1]) * inv_det,
                    (jtj[0][0] * jtj[1][1] - jtj[0][1] * jtj[1][0]) * inv_det,
                ],
            ];

            let dx = inv[0][0] * jtr[0] + inv[0][1] * jtr[1] + inv[0][2] * jtr[2];
            let dy = inv[1][0] * jtr[0] + inv[1][1] * jtr[1] + inv[1][2] * jtr[2];
            let dz = inv[2][0] * jtr[0] + inv[2][1] * jtr[1] + inv[2][2] * jtr[2];

            x -= dx;
            y -= dy;
            z -= dz;

            let step = (dx * dx + dy * dy + dz * dz).sqrt();
            if _iter % 10 == 0 || step < 1e-6 {
                eprintln!(
                    "[DEBUG_RTGS] Iter {} | step: {:.6} | var_sum: {:.6} | lambda: {:.6} | pos: ({:.1}, {:.1}, {:.1})",
                    _iter, step, var_sum, lambda, x, y, z
                );
            }
            if step < 1e-6 {
                break;
            }
        }

        // Validate solved position is physically reasonable (within ~100km of Earth's surface)
        let earth_radius = 6378137.0;
        let dist_from_center = (x * x + y * y + z * z).sqrt();
        if dist_from_center < earth_radius * 0.9 || dist_from_center > earth_radius + 100_000.0 {
            return None;
        }

        // Convert ECEF to geodetic (WGS84)
        let geo = ecef_to_geodetic(x, y, z);
        Some(geo)
    }
}

/// Convert ECEF coordinates to geodetic (lat, lon, alt) using iterative method.
fn ecef_to_geodetic(x: f64, y: f64, z: f64) -> GeodeticCoordinates {
    let a = 6378137.0_f64; // WGS84 semi-major axis
    let f = 1.0 / 298.257223563;
    let e2 = 2.0 * f - f * f;

    let lon = y.atan2(x).to_degrees();
    let p = (x * x + y * y).sqrt();

    // Iterative latitude computation
    let mut lat = (z / p).atan();
    for _ in 0..10 {
        let sin_lat = lat.sin();
        let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        lat = (z + e2 * n * sin_lat).atan2(p);
    }

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let alt = if cos_lat.abs() > 1e-10 {
        p / cos_lat - n
    } else {
        z.abs() / sin_lat.abs() - n * (1.0 - e2)
    };

    GeodeticCoordinates {
        latitude: lat.to_degrees(),
        longitude: lon,
        altitude: alt,
    }
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Sp3Coordinate {
    pub time: DateTime<Utc>,
    pub position: [f64; 3],          // meters
    pub velocity: Option<[f64; 3]>,   // meters/sec
    pub clock_offset: f64,           // seconds
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct SatelliteOrbit {
    pub sat_name: String,
    pub coordinates: Vec<Sp3Coordinate>,
}

impl SatelliteOrbit {
    pub fn get_chebyshev_model(&self, dt: DateTime<Utc>) -> Option<ChebyshevOrbit> {
        let coords = &self.coordinates;
        if coords.len() < 10 {
            return None;
        }

        let t_target = dt;
        let idx_res = coords.binary_search_by_key(&t_target, |c| c.time);
        let idx = match idx_res {
            Ok(i) => i,
            Err(i) => {
                if i == 0 {
                    0
                } else if i >= coords.len() {
                    coords.len() - 1
                } else {
                    i - 1
                }
            }
        };

        let start_idx = if idx < 4 {
            0
        } else if idx + 5 >= coords.len() {
            coords.len() - 10
        } else {
            idx - 4
        };

        let window = &coords[start_idx..start_idx + 10];

        let t0 = window[0].time.timestamp_nanos_opt()? as f64 * 1e-9;
        let t9 = window[9].time.timestamp_nanos_opt()? as f64 * 1e-9;

        // Spacing check:
        let mut uniform = true;
        let dt_expected = (t9 - t0) / 9.0;
        for k in 0..9 {
            let tk = window[k].time.timestamp_nanos_opt()? as f64 * 1e-9;
            let tk1 = window[k+1].time.timestamp_nanos_opt()? as f64 * 1e-9;
            if (tk1 - tk - dt_expected).abs() > 1.0 {
                uniform = false;
                break;
            }
        }

        let mut cx = [0.0; 10];
        let mut cy = [0.0; 10];
        let mut cz = [0.0; 10];

        if uniform {
            static A_INV: std::sync::OnceLock<[[f64; 10]; 10]> = std::sync::OnceLock::new();
            let a_inv = A_INV.get_or_init(|| {
                let mut a = [[0.0; 10]; 10];
                for k in 0..10 {
                    let tau = -1.0 + 2.0 * k as f64 / 9.0;
                    let mut t_val = [0.0; 10];
                    t_val[0] = 1.0;
                    t_val[1] = tau;
                    for j in 2..10 {
                        t_val[j] = 2.0 * tau * t_val[j - 1] - t_val[j - 2];
                    }
                    for j in 0..10 {
                        a[k][j] = t_val[j];
                    }
                }
                invert_matrix_10(&a).expect("Failed to invert Chebyshev matrix")
            });

            for j in 0..10 {
                let mut sum_x = 0.0;
                let mut sum_y = 0.0;
                let mut sum_z = 0.0;
                for k in 0..10 {
                    sum_x += a_inv[j][k] * window[k].position[0];
                    sum_y += a_inv[j][k] * window[k].position[1];
                    sum_z += a_inv[j][k] * window[k].position[2];
                }
                cx[j] = sum_x;
                cy[j] = sum_y;
                cz[j] = sum_z;
            }
        } else {
            let mut a = [[0.0; 10]; 10];
            for k in 0..10 {
                let tk = window[k].time.timestamp_nanos_opt()? as f64 * 1e-9;
                let tau = 2.0 * (tk - t0) / (t9 - t0) - 1.0;
                let mut t_val = [0.0; 10];
                t_val[0] = 1.0;
                t_val[1] = tau;
                for j in 2..10 {
                    t_val[j] = 2.0 * tau * t_val[j - 1] - t_val[j - 2];
                }
                for j in 0..10 {
                    a[k][j] = t_val[j];
                }
            }
            if let Some(a_inv) = invert_matrix_10(&a) {
                for j in 0..10 {
                    let mut sum_x = 0.0;
                    let mut sum_y = 0.0;
                    let mut sum_z = 0.0;
                    for k in 0..10 {
                        sum_x += a_inv[j][k] * window[k].position[0];
                        sum_y += a_inv[j][k] * window[k].position[1];
                        sum_z += a_inv[j][k] * window[k].position[2];
                    }
                    cx[j] = sum_x;
                    cy[j] = sum_y;
                    cz[j] = sum_z;
                }
            } else {
                return None;
            }
        }

        Some(ChebyshevOrbit { t0, t9, cx, cy, cz })
    }

    pub fn propagate_ecef(&self, dt: DateTime<Utc>) -> Option<([f64; 3], [f64; 3])> {
        let model = self.get_chebyshev_model(dt)?;
        model.evaluate(dt)
    }
}

#[derive(Clone, Debug)]
pub struct ChebyshevOrbit {
    pub t0: f64,
    pub t9: f64,
    pub cx: [f64; 10],
    pub cy: [f64; 10],
    pub cz: [f64; 10],
}

impl ChebyshevOrbit {
    pub fn evaluate(&self, dt: DateTime<Utc>) -> Option<([f64; 3], [f64; 3])> {
        let t = dt.timestamp_nanos_opt()? as f64 * 1e-9;
        let tau = 2.0 * (t - self.t0) / (self.t9 - self.t0) - 1.0;
        let (x, dx) = evaluate_chebyshev(&self.cx, tau);
        let (y, dy) = evaluate_chebyshev(&self.cy, tau);
        let (z, dz) = evaluate_chebyshev(&self.cz, tau);

        let dtau_dt = 2.0 / (self.t9 - self.t0);
        let vx = dx * dtau_dt;
        let vy = dy * dtau_dt;
        let vz = dz * dtau_dt;

        Some(([x, y, z], [vx, vy, vz]))
    }
}

fn invert_matrix_10(a: &[[f64; 10]; 10]) -> Option<[[f64; 10]; 10]> {
    let mut temp = *a;
    let mut inv = [[0.0; 10]; 10];
    for i in 0..10 {
        inv[i][i] = 1.0;
    }

    for i in 0..10 {
        let mut max_row = i;
        let mut max_val = temp[i][i].abs();
        for r in (i + 1)..10 {
            if temp[r][i].abs() > max_val {
                max_val = temp[r][i].abs();
                max_row = r;
            }
        }

        if max_val < 1e-12 {
            return None;
        }

        if max_row != i {
            temp.swap(i, max_row);
            inv.swap(i, max_row);
        }

        let pivot = temp[i][i];
        for c in 0..10 {
            temp[i][c] /= pivot;
            inv[i][c] /= pivot;
        }

        for r in 0..10 {
            if r != i {
                let factor = temp[r][i];
                for c in 0..10 {
                    temp[r][c] -= factor * temp[i][c];
                    inv[r][c] -= factor * inv[i][c];
                }
            }
        }
    }

    Some(inv)
}

fn evaluate_chebyshev(c: &[f64; 10], tau: f64) -> (f64, f64) {
    let mut t = [0.0; 10];
    let mut dt = [0.0; 10];

    t[0] = 1.0;
    dt[0] = 0.0;

    t[1] = tau;
    dt[1] = 1.0;

    for j in 2..10 {
        t[j] = 2.0 * tau * t[j - 1] - t[j - 2];
        dt[j] = 2.0 * t[j - 1] + 2.0 * tau * dt[j - 1] - dt[j - 2];
    }

    let mut pos = 0.0;
    let mut d_tau = 0.0;
    for j in 0..10 {
        pos += c[j] * t[j];
        d_tau += c[j] * dt[j];
    }

    (pos, d_tau)
}

#[derive(Clone, Debug)]
pub enum OrbitModel {
    Sgp4 {
        elements: sgp4::Elements,
        constants: sgp4::Constants,
    },
    Sp3 {
        orbit: SatelliteOrbit,
    },
}

impl OrbitModel {
    pub fn sat_name(&self) -> &str {
        match self {
            OrbitModel::Sgp4 { elements, .. } => elements.object_name.as_deref().unwrap_or("UNKNOWN"),
            OrbitModel::Sp3 { orbit } => &orbit.sat_name,
        }
    }

    pub fn epoch(&self) -> DateTime<Utc> {
        match self {
            OrbitModel::Sgp4 { elements, .. } => elements.datetime.and_utc(),
            OrbitModel::Sp3 { orbit } => {
                if orbit.coordinates.is_empty() {
                    Utc::now()
                } else {
                    orbit.coordinates[0].time
                }
            }
        }
    }

    pub fn propagate_ecef(&self, dt: DateTime<Utc>) -> Option<([f64; 3], [f64; 3])> {
        match self {
            OrbitModel::Sgp4 { elements, constants } => {
                let duration_since_epoch = dt.naive_utc().signed_duration_since(elements.datetime);
                let mins_since_epoch = duration_since_epoch.num_milliseconds() as f64 / 60000.0;
                if let Ok(prediction) = constants.propagate(sgp4::MinutesSinceEpoch(mins_since_epoch)) {
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
                    Some(teme_to_ecef(jd, pos_teme, vel_teme))
                } else {
                    None
                }
            }
            OrbitModel::Sp3 { orbit } => {
                orbit.propagate_ecef(dt)
            }
        }
    }
}

pub fn parse_sp3_file(path: &str) -> io::Result<Vec<SatelliteOrbit>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::collections::HashMap;

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut satellites: HashMap<String, Vec<Sp3Coordinate>> = HashMap::new();
    let mut current_time: Option<DateTime<Utc>> = None;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }

        if parts[0] == "*" {
            if parts.len() >= 7 {
                let year: i32 = parts[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let month: u32 = parts[2].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let day: u32 = parts[3].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let hour: u32 = parts[4].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let minute: u32 = parts[5].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let second_f: f64 = parts[6].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let second = second_f.trunc() as u32;
                let nanosecond = (second_f.fract() * 1_000_000_000.0).round() as u32;

                if let Some(naive_dt) = chrono::NaiveDate::from_ymd_opt(year, month, day)
                    .and_then(|d| d.and_hms_nano_opt(hour, minute, second, nanosecond))
                {
                    current_time = Some(DateTime::from_naive_utc_and_offset(naive_dt, Utc));
                }
            }
        } else if parts[0].starts_with('P') {
            let raw_sat = parts[0];
            if raw_sat.len() >= 2 {
                let sat_name = raw_sat[1..].trim().to_string();

                if let Some(time) = current_time {
                    if parts.len() >= 5 {
                        let x: f64 = parts[1].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? * 1000.0;
                        let y: f64 = parts[2].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? * 1000.0;
                        let z: f64 = parts[3].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? * 1000.0;
                        let clock_offset_us: f64 = parts[4].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                        let clock_offset = clock_offset_us * 1e-6;

                        let coord = Sp3Coordinate {
                            time,
                            position: [x, y, z],
                            velocity: None,
                            clock_offset,
                        };

                        satellites.entry(sat_name).or_default().push(coord);
                    }
                }
            }
        } else if parts[0].starts_with('V') {
            let raw_sat = parts[0];
            if raw_sat.len() >= 2 {
                let sat_name = raw_sat[1..].trim().to_string();

                if let Some(coords) = satellites.get_mut(&sat_name) {
                    if let Some(last_coord) = coords.last_mut() {
                        if parts.len() >= 4 {
                            let vx: f64 = parts[1].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? * 0.1;
                            let vy: f64 = parts[2].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? * 0.1;
                            let vz: f64 = parts[3].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? * 0.1;
                            last_coord.velocity = Some([vx, vy, vz]);
                        }
                    }
                }
            }
        }
    }

    let result = satellites
        .into_iter()
        .map(|(name, mut coords)| {
            coords.sort_by_key(|c| c.time);
            SatelliteOrbit { sat_name: name, coordinates: coords }
        })
        .collect();

    Ok(result)
}

pub fn load_orbits(path: &str) -> io::Result<Vec<(String, OrbitModel)>> {
    let mut satellites = Vec::new();
    if let Ok(sp3_sats) = parse_sp3_file(path) {
        if !sp3_sats.is_empty() {
            for sat in sp3_sats {
                satellites.push((sat.sat_name.clone(), OrbitModel::Sp3 { orbit: sat }));
            }
            return Ok(satellites);
        }
    }
    let tle_sats = parse_tle_file(path)?;
    for (name, elements) in tle_sats {
        if let Ok(constants) = sgp4::Constants::from_elements(&elements) {
            satellites.push((name, OrbitModel::Sgp4 { elements, constants }));
        }
    }
    Ok(satellites)
}

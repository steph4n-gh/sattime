use chrono::{DateTime, Datelike, Timelike, Utc};

// Mock/helper implementations of coordinate functions so orbit_solver compiles
pub mod orbit {
    use super::*;
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
} // end mod orbit

#[path = "../src/orbit_solver.rs"]
mod orbit_solver;

fn verify_tle_line_checksum(line: &str) -> bool {
    if line.len() != 69 {
        return false;
    }
    let mut sum = 0;
    for c in line.chars().take(68) {
        if c.is_ascii_digit() {
            sum += c.to_digit(10).unwrap();
        } else if c == '-' {
            sum += 1;
        }
    }
    let expected = (sum % 10) as u8 + b'0';
    line.as_bytes()[68] == expected
}

#[test]
fn test_verify_orbit_solver_convergence_and_tle() {
    let truth_a = 6378137.0 + 550000.0; // 550 km altitude circular orbit
    let truth_i = 53.0_f64.to_radians(); // 53 degrees inclination
    let truth_raan = 1.2;
    let truth_u0 = 0.5;
    let center_freq = 150800000.0;

    let epoch = DateTime::parse_from_rfc3339("2026-06-09T12:00:00Z")
        .unwrap()
        .with_timezone(&Utc);

    // Get the satellite ECEF position at epoch
    let (pos_sat_epoch, _) =
        orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, epoch);

    // Convert satellite position at epoch to geodetic coordinates, and place the observer there
    let (rec_lat, rec_lon, _) = orbit::ecef_to_wgs84(pos_sat_epoch);
    let rec_alt = 0.0;
    let rec_ecef = orbit::wgs84_to_ecef(rec_lat, rec_lon, rec_alt);
    println!(
        "Dynamically placed observer directly under satellite path at lat={:.4}, lon={:.4}",
        rec_lat, rec_lon
    );

    // Simulate passes when the range is < 2500 km
    let mut passes = Vec::new();
    let mut current_pass: Option<orbit_solver::RawPass> = None;
    let mut last_t: Option<DateTime<Utc>> = None;

    // Propagate over 8 orbits (48000 seconds) to ensure multiple overhead passes are captured
    for step in 0..9600 {
        // 9600 steps * 5s = 48000 seconds
        let t = epoch + chrono::Duration::seconds(step * 5);
        let (pos_sat, _) =
            orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, t);
        let dx = pos_sat[0] - rec_ecef[0];
        let dy = pos_sat[1] - rec_ecef[1];
        let dz = pos_sat[2] - rec_ecef[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < 2500000.0 {
            // 2500 km range limit
            let pred_shift = orbit_solver::predict_frequency(
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
            let observed_freq = center_freq + pred_shift;

            let point = orbit_solver::PassPoint {
                time: t,
                freq: observed_freq,
            };

            if let Some(ref mut pass) = current_pass {
                if let Some(lt) = last_t {
                    if (t - lt).num_seconds() > 300 {
                        passes.push(current_pass.take().unwrap());
                        current_pass = Some(orbit_solver::RawPass {
                            sat_name: "SIM_SAT".to_string(),
                            center_freq,
                            points: vec![point],
                        });
                    } else {
                        pass.points.push(point);
                    }
                }
            } else {
                current_pass = Some(orbit_solver::RawPass {
                    sat_name: "SIM_SAT".to_string(),
                    center_freq,
                    points: vec![point],
                });
            }
            last_t = Some(t);
        }
    }
    if let Some(pass) = current_pass {
        passes.push(pass);
    }

    println!("Simulated {} passes.", passes.len());
    for (idx, pass) in passes.iter().enumerate() {
        println!(
            "Pass {}: {} points, duration: {}s",
            idx + 1,
            pass.points.len(),
            (pass.points.last().unwrap().time - pass.points[0].time).num_seconds()
        );
    }

    assert!(
        passes.len() >= 2,
        "Need at least 2 simulated passes for the orbit solver test, got {}",
        passes.len()
    );

    // Fit circular orbit from simulated passes
    // Provide initial guesses close to the true circular orbit to verify optimizer convergence
    let initial_a = 6378137.0 + 500000.0; // 500 km guess
    let initial_i = 50.0_f64.to_radians(); // 50 degrees guess

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, initial_a, initial_i)
        .expect("Gauss-Newton Doppler orbit solver failed to converge");

    println!(
        "Solved orbit: a={:.1} (truth={:.1}), i={:.4} (truth={:.4}), raan0={:.4} (truth={:.4}), u0={:.4} (truth={:.4})",
        solved.a, truth_a, solved.i, truth_i, solved.raan0, truth_raan, solved.u0, truth_u0
    );

    // Assert that solver converges back to true orbital elements
    assert!(
        (solved.a - truth_a).abs() < 500.0,
        "Semi-major axis mismatch: {} vs {}",
        solved.a,
        truth_a
    );
    assert!(
        (solved.i - truth_i).abs().to_degrees() < 0.1,
        "Inclination mismatch: {} vs {}",
        solved.i.to_degrees(),
        truth_i.to_degrees()
    );
    assert!(
        (solved.raan0 - truth_raan).abs().to_degrees() < 2.0,
        "RAAN mismatch: {} vs {}",
        solved.raan0.to_degrees(),
        truth_raan.to_degrees()
    );
    assert!(
        (solved.u0 - truth_u0).abs().to_degrees() < 15.0,
        "Argument of latitude mismatch: {} vs {}",
        solved.u0.to_degrees(),
        truth_u0.to_degrees()
    );

    // Format resolved orbit to standard TLE lines
    let tle_text = orbit_solver::format_tle_catalog("SIM_SAT", &solved);
    println!("Generated TLE:\n{}", tle_text);

    let lines: Vec<&str> = tle_text.lines().collect();
    assert_eq!(
        lines.len(),
        3,
        "TLE output should consist of exactly 3 lines"
    );
    assert_eq!(lines[0], "SIM_SAT");

    assert!(lines[1].starts_with("1 "), "Line 1 must start with '1 '");
    assert!(lines[2].starts_with("2 "), "Line 2 must start with '2 '");

    assert_eq!(
        lines[1].len(),
        69,
        "Line 1 should be exactly 69 characters long"
    );
    assert_eq!(
        lines[2].len(),
        69,
        "Line 2 should be exactly 69 characters long"
    );

    assert!(
        verify_tle_line_checksum(lines[1]),
        "Line 1 checksum verification failed"
    );
    assert!(
        verify_tle_line_checksum(lines[2]),
        "Line 2 checksum verification failed"
    );
}

// Duplicate legacy function to act as a baseline benchmark
pub fn fit_orbit_doppler_grid_legacy(
    raw_passes: &[orbit_solver::RawPass],
    rec_ecef: [f64; 3],
    initial_a: f64,
    initial_i: f64,
) -> Result<orbit_solver::SolvedOrbit, Box<dyn std::error::Error>> {
    let n_passes = raw_passes.len();
    if n_passes < 2 {
        return Err("Need at least 2 passes to resolve orbit parameters.".into());
    }

    let epoch = raw_passes[0].points[0].time;
    let center_freq = raw_passes[0].center_freq;

    let n_params = 4 + 2 * n_passes;
    let mut params = vec![0.0; n_params];

    params[0] = initial_a;
    params[1] = initial_i;

    // Coarse grid search to initialize raan0 and u0
    let mut best_raan0 = 0.0;
    let mut best_u0 = 0.0;
    let mut best_grid_rss = f64::MAX;

    for r_idx in 0..12 {
        let test_raan0 = (r_idx as f64) * 30.0_f64.to_radians();
        for u_idx in 0..12 {
            let test_u0 = (u_idx as f64) * 30.0_f64.to_radians();
            let mut rss = 0.0;

            let pred_pcas = orbit_solver::get_pred_pca_times(
                initial_a, initial_i, test_raan0, test_u0, epoch, rec_ecef, raw_passes,
            );
            for (p_idx, pass) in raw_passes.iter().enumerate() {
                let mut obs_pca_time = pass.points[0].time;
                let mut min_offset = f64::MAX;
                for pt in &pass.points {
                    let off = (pt.freq - center_freq).abs();
                    if off < min_offset {
                        min_offset = off;
                        obs_pca_time = pt.time;
                    }
                }

                let pred_pca_time = pred_pcas[p_idx];

                let dt = (pred_pca_time - obs_pca_time).num_milliseconds() as f64 / 1000.0;

                for pt in &pass.points {
                    let pred_f = orbit_solver::predict_frequency(
                        initial_a,
                        initial_i,
                        test_raan0,
                        test_u0,
                        epoch,
                        pt.time,
                        dt,
                        0.0,
                        center_freq,
                        rec_ecef,
                    );
                    let diff = pt.freq - (center_freq + pred_f);
                    rss += diff * diff;
                }
            }

            if rss < best_grid_rss {
                best_grid_rss = rss;
                best_raan0 = test_raan0;
                best_u0 = test_u0;
            }
        }
    }

    params[2] = best_raan0;
    params[3] = best_u0;

    let pred_pcas = orbit_solver::get_pred_pca_times(
        params[0], params[1], params[2], params[3], epoch, rec_ecef, raw_passes,
    );
    for (p_idx, pass) in raw_passes.iter().enumerate() {
        let mut obs_pca_time = pass.points[0].time;
        let mut min_offset = f64::MAX;
        for pt in &pass.points {
            let off = (pt.freq - center_freq).abs();
            if off < min_offset {
                min_offset = off;
                obs_pca_time = pt.time;
            }
        }

        let pred_pca_time = pred_pcas[p_idx];

        params[4 + 2 * p_idx] = (pred_pca_time - obs_pca_time).num_milliseconds() as f64 / 1000.0;
        params[4 + 2 * p_idx + 1] = 0.0;
    }

    let mut lambda = 1.0;
    let mut best_rss = f64::MAX;
    let mut best_params = params.clone();

    for _ in 0..150 {
        let mut residuals = Vec::new();
        for (p_idx, pass) in raw_passes.iter().enumerate() {
            let dt = params[4 + 2 * p_idx];
            let df = params[4 + 2 * p_idx + 1];
            for pt in &pass.points {
                let pred = orbit_solver::predict_frequency(
                    params[0],
                    params[1],
                    params[2],
                    params[3],
                    epoch,
                    pt.time,
                    dt,
                    df,
                    center_freq,
                    rec_ecef,
                );
                residuals.push(pt.freq - (center_freq + pred));
            }
        }

        let rss: f64 = residuals.iter().map(|r| r * r).sum();
        if rss < best_rss {
            best_rss = rss;
            best_params = params.clone();
            lambda /= 10.0;
        } else {
            params = best_params.clone();
            lambda *= 10.0;
            if lambda > 1e12 {
                break;
            }
            // Recompute residuals for the restored params (best_params)
            residuals.clear();
            for (p_idx, pass) in raw_passes.iter().enumerate() {
                let dt = params[4 + 2 * p_idx];
                let df = params[4 + 2 * p_idx + 1];
                for pt in &pass.points {
                    let pred = orbit_solver::predict_frequency(
                        params[0],
                        params[1],
                        params[2],
                        params[3],
                        epoch,
                        pt.time,
                        dt,
                        df,
                        center_freq,
                        rec_ecef,
                    );
                    residuals.push(pt.freq - (center_freq + pred));
                }
            }
        }

        let n_obs = residuals.len();
        if n_obs < n_params {
            return Err("Not enough observations to fit orbit parameters.".into());
        }

        let mut jacobian = vec![vec![0.0; n_params]; n_obs];

        for k in 0..n_params {
            let mut perturbed = params.clone();
            let param_eps = if k == 0 {
                10.0
            } else if k == 1 || k == 2 || k == 3 {
                1e-6
            } else if (k - 4) % 2 == 0 {
                1e-3
            } else {
                1e-2
            };
            perturbed[k] += param_eps;

            let mut row_idx = 0;
            for (p_idx, pass) in raw_passes.iter().enumerate() {
                let dt = perturbed[4 + 2 * p_idx];
                let df = perturbed[4 + 2 * p_idx + 1];
                for pt in &pass.points {
                    let pred = orbit_solver::predict_frequency(
                        perturbed[0],
                        perturbed[1],
                        perturbed[2],
                        perturbed[3],
                        epoch,
                        pt.time,
                        dt,
                        df,
                        center_freq,
                        rec_ecef,
                    );
                    let diff = pt.freq - (center_freq + pred);
                    jacobian[row_idx][k] = (diff - residuals[row_idx]) / param_eps;
                    row_idx += 1;
                }
            }
        }

        let mut jt_j = vec![vec![0.0; n_params]; n_params];
        let mut jt_r = vec![0.0; n_params];

        for row in 0..n_obs {
            for c1 in 0..n_params {
                jt_r[c1] += jacobian[row][c1] * residuals[row];
                for c2 in 0..n_params {
                    jt_j[c1][c2] += jacobian[row][c1] * jacobian[row][c2];
                }
            }
        }

        for k in 0..n_params {
            jt_j[k][k] += lambda * jt_j[k][k];
        }

        if let Some(delta) = orbit_solver::solve_linear_system(&mut jt_j, &jt_r) {
            let mut step_len = 0.0;
            for k in 0..n_params {
                params[k] -= delta[k];
                step_len += delta[k] * delta[k];
            }

            params[0] = params[0].max(6500e3).min(20000e3);
            params[1] = params[1].max(0.0).min(std::f64::consts::PI);

            if step_len.sqrt() < 1e-6 {
                break;
            }
        } else {
            break;
        }
    }

    params = best_params;
    let mut pass_dts = Vec::new();
    let mut pass_dfs = Vec::new();
    for j in 0..n_passes {
        pass_dts.push(params[4 + 2 * j]);
        pass_dfs.push(params[4 + 2 * j + 1]);
    }

    Ok(orbit_solver::SolvedOrbit {
        a: params[0],
        i: params[1],
        raan0: params[2],
        u0: params[3],
        epoch,
        pass_dts,
        pass_dfs,
        pass_df1s: Vec::new(),
        pass_df2s: Vec::new(),
    })
}

#[test]
fn test_langevin_global_solver() {
    let truth_a = 6378137.0 + 550000.0; // 550 km altitude circular orbit
    let truth_i = 53.0_f64.to_radians(); // 53 degrees inclination
    // Use true raan0 and u0 far from 30 degree grid points (e.g. 45.0 degrees)
    let truth_raan = 45.0_f64.to_radians();
    let truth_u0 = 45.0_f64.to_radians();
    let center_freq = 150800000.0;

    let epoch = DateTime::parse_from_rfc3339("2026-06-09T12:00:00Z")
        .unwrap()
        .with_timezone(&Utc);

    // Get the satellite ECEF position at epoch
    let (pos_sat_epoch, _) =
        orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, epoch);

    // Convert satellite position at epoch to geodetic coordinates, and place the observer there
    let (rec_lat, rec_lon, _) = orbit::ecef_to_wgs84(pos_sat_epoch);
    let rec_alt = 0.0;
    let rec_ecef = orbit::wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    // Simulate passes when the range is < 2500 km
    let mut passes = Vec::new();
    let mut current_pass: Option<orbit_solver::RawPass> = None;
    let mut last_t: Option<DateTime<Utc>> = None;

    // Propagate over 8 orbits to capture multiple passes
    for step in 0..9600 {
        let t = epoch + chrono::Duration::seconds(step * 5);
        let (pos_sat, _) =
            orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, t);
        let dx = pos_sat[0] - rec_ecef[0];
        let dy = pos_sat[1] - rec_ecef[1];
        let dz = pos_sat[2] - rec_ecef[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < 2500000.0 {
            let pred_shift = orbit_solver::predict_frequency(
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
            let observed_freq = center_freq + pred_shift;

            let point = orbit_solver::PassPoint {
                time: t,
                freq: observed_freq,
            };

            if let Some(ref mut pass) = current_pass {
                if let Some(lt) = last_t {
                    if (t - lt).num_seconds() > 300 {
                        passes.push(current_pass.take().unwrap());
                        current_pass = Some(orbit_solver::RawPass {
                            sat_name: "SIM_SAT".to_string(),
                            center_freq,
                            points: vec![point],
                        });
                    } else {
                        pass.points.push(point);
                    }
                }
            } else {
                current_pass = Some(orbit_solver::RawPass {
                    sat_name: "SIM_SAT".to_string(),
                    center_freq,
                    points: vec![point],
                });
            }
            last_t = Some(t);
        }
    }
    if let Some(pass) = current_pass {
        passes.push(pass);
    }

    assert!(
        passes.len() >= 2,
        "Need at least 2 simulated passes, got {}",
        passes.len()
    );

    // Bad initial guesses (a off by 800 km, inclination off by 10 deg)
    let initial_a = truth_a + 800000.0;
    let initial_i = truth_i + 10.0_f64.to_radians();

    // Run legacy grid search LM
    let solved_legacy_res = fit_orbit_doppler_grid_legacy(&passes, rec_ecef, initial_a, initial_i);
    let mut legacy_success = false;
    let mut legacy_rss = f64::MAX;
    if let Ok(solved) = solved_legacy_res {
        let err_a = (solved.a - truth_a).abs();
        let err_i = (solved.i - truth_i).abs().to_degrees();
        println!(
            "Legacy grid search solved: a={:.1} (err={:.1}m), i={:.4} (err={:.4} deg)",
            solved.a,
            err_a,
            solved.i.to_degrees(),
            err_i
        );

        legacy_rss = orbit_solver::compute_rss(
            &passes,
            rec_ecef,
            solved.a,
            solved.i,
            solved.epoch,
            center_freq,
            solved.raan0,
            solved.u0,
        );
        if err_a < 500.0 && err_i < 0.1 {
            legacy_success = true;
        }
    } else {
        println!("Legacy grid search failed to solve/converge.");
    }

    // Run new Langevin Global LM
    let solved_adelic_res =
        orbit_solver::fit_orbit_doppler(&passes, rec_ecef, initial_a, initial_i);
    let mut adelic_success = false;
    let mut adelic_rss = f64::MAX;
    if let Ok(solved) = solved_adelic_res {
        let err_a = (solved.a - truth_a).abs();
        let err_i = (solved.i - truth_i).abs().to_degrees();
        println!(
            "Langevin Global solved: a={:.1} (err={:.1}m), i={:.4} (err={:.4} deg)",
            solved.a,
            err_a,
            solved.i.to_degrees(),
            err_i
        );

        adelic_rss = orbit_solver::compute_rss(
            &passes,
            rec_ecef,
            solved.a,
            solved.i,
            solved.epoch,
            center_freq,
            solved.raan0,
            solved.u0,
        );
        if err_a < 500.0 && err_i < 0.1 {
            adelic_success = true;
        }
    } else {
        println!("Langevin Global failed to solve/converge.");
    }

    println!("Legacy Success: {}, RSS: {}", legacy_success, legacy_rss);
    println!("Adelic Success: {}, RSS: {}", adelic_success, adelic_rss);

    // Assert Langevin Global converged to the correct parameters (semi-major axis error < 500m, inclination error < 0.1 deg)
    assert!(
        adelic_success,
        "Langevin Global solver did not converge to correct parameters under bad initial guess"
    );

    // Assert Langevin Global has higher or equal success rate and smaller/equal RSS than legacy
    if legacy_success {
        assert!(
            adelic_success,
            "Langevin Global should succeed if legacy grid search succeeded"
        );
        assert!(
            adelic_rss <= legacy_rss + 1.0,
            "Langevin Global RSS ({}) should be smaller or equal to legacy grid search RSS ({})",
            adelic_rss,
            legacy_rss
        );
    }
}

#[test]
fn test_virtual_tcxo_polynomial_fit() {
    let truth_a = 6378137.0 + 550000.0; // 550 km altitude circular orbit
    let truth_i = 53.0_f64.to_radians(); // 53 degrees inclination
    let truth_raan = 1.2;
    let truth_u0 = 0.5;
    let center_freq = 150800000.0;

    let epoch = DateTime::parse_from_rfc3339("2026-06-09T12:00:00Z")
        .unwrap()
        .with_timezone(&Utc);

    // Get the satellite ECEF position at epoch
    let (pos_sat_epoch, _) =
        orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, epoch);

    // Convert satellite position at epoch to geodetic coordinates, and place the observer there
    let (rec_lat, rec_lon, _) = orbit::ecef_to_wgs84(pos_sat_epoch);
    let rec_alt = 0.0;
    let rec_ecef = orbit::wgs84_to_ecef(rec_lat, rec_lon, rec_alt);

    // Simulate passes when the range is < 2500 km
    let mut passes = Vec::new();
    let mut current_pass: Option<orbit_solver::RawPass> = None;
    let mut last_t: Option<DateTime<Utc>> = None;

    // Start propagating 1200 seconds before epoch to capture the full first pass
    for step in -240..9600 {
        let t = epoch + chrono::Duration::seconds((step as i64) * 5);
        let (pos_sat, _) =
            orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, t);
        let dx = pos_sat[0] - rec_ecef[0];
        let dy = pos_sat[1] - rec_ecef[1];
        let dz = pos_sat[2] - rec_ecef[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < 2500000.0 {
            let pred_shift = orbit_solver::predict_frequency(
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
            let observed_freq = center_freq + pred_shift;

            let point = orbit_solver::PassPoint {
                time: t,
                freq: observed_freq,
            };

            if let Some(ref mut pass) = current_pass {
                if let Some(lt) = last_t {
                    if (t - lt).num_seconds() > 300 {
                        passes.push(current_pass.take().unwrap());
                        current_pass = Some(orbit_solver::RawPass {
                            sat_name: "SIM_SAT".to_string(),
                            center_freq,
                            points: vec![point],
                        });
                    } else {
                        pass.points.push(point);
                    }
                }
            } else {
                current_pass = Some(orbit_solver::RawPass {
                    sat_name: "SIM_SAT".to_string(),
                    center_freq,
                    points: vec![point],
                });
            }
            last_t = Some(t);
        }
    }
    if let Some(pass) = current_pass {
        passes.push(pass);
    }

    assert!(
        passes.len() >= 2,
        "Need at least 2 simulated passes, got {}",
        passes.len()
    );

    // Find the pair of passes that minimizes synodic period distortion
    let mut best_pair = (0, 1);
    let mut best_diff = f64::MAX;
    for i in 0..passes.len() {
        for j in i + 1..passes.len() {
            let mut obs_pca_i = passes[i].points[0].time;
            let mut min_off_i = f64::MAX;
            for pt in &passes[i].points {
                let off = (pt.freq - center_freq).abs();
                if off < min_off_i {
                    min_off_i = off;
                    obs_pca_i = pt.time;
                }
            }
            let mut obs_pca_j = passes[j].points[0].time;
            let mut min_off_j = f64::MAX;
            for pt in &passes[j].points {
                let off = (pt.freq - center_freq).abs();
                if off < min_off_j {
                    min_off_j = off;
                    obs_pca_j = pt.time;
                }
            }
            let dt = (obs_pca_j - obs_pca_i).num_milliseconds() as f64 / 1000.0;
            let k = (dt / 5700.0).round();
            if k > 0.0 {
                let dt_k = dt / k;
                let diff = (dt_k - 5740.0).abs();
                if diff < best_diff {
                    best_diff = diff;
                    best_pair = (i, j);
                }
            }
        }
    }

    println!("DEBUG: Selected best pair of passes: {:?}", best_pair);
    let mut selected_passes = vec![passes[best_pair.0].clone(), passes[best_pair.1].clone()];

    // Inject second-order polynomial thermal drift offset:
    // df(t) = 15.0 + 0.1 * tau - 0.0002 * tau^2
    for pass in &mut selected_passes {
        let t_ref = pass.points[0].time;
        for pt in &mut pass.points {
            let tau = (pt.time - t_ref).num_milliseconds() as f64 / 1000.0;
            let drift = 15.0 + 0.1 * tau - 0.0002 * tau * tau;
            pt.freq += drift;
        }
    }

    // Fit circular orbit from simulated passes with polynomial drift
    let initial_a = truth_a + 500.0; // close guess
    let initial_i = truth_i + 0.01;

    let solved = orbit_solver::fit_orbit_doppler(&selected_passes, rec_ecef, initial_a, initial_i)
        .expect("Virtual TCXO Doppler orbit solver failed to converge");

    println!("Solved orbit: a={:.1}, i={:.4}, raan0={:.4}, u0={:.4}", solved.a, solved.i, solved.raan0, solved.u0);
    for idx in 0..2 {
        println!(
            "Pass {}: df0 = {:.4} Hz (target=15.0), df1 = {:.6} Hz/s (target=0.1), df2 = {:.8} Hz/s^2 (target=-0.0002)",
            idx, solved.pass_dfs[idx], solved.pass_df1s[idx], solved.pass_df2s[idx]
        );
    }

    // Verify solved parameters match the injected drift parameters
    for idx in 0..2 {
        let err_df0 = (solved.pass_dfs[idx] - 15.0).abs();
        let err_df1 = (solved.pass_df1s[idx] - 0.1).abs();
        let err_df2 = (solved.pass_df2s[idx] - (-0.0002)).abs();

        assert!(err_df0 < 1.0, "Pass {} df0 error too high: {}", idx, err_df0);
        assert!(err_df1 < 0.05, "Pass {} df1 error too high: {}", idx, err_df1);
        assert!(err_df2 < 0.0001, "Pass {} df2 error too high: {}", idx, err_df2);
    }
}


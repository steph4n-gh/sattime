use chrono::{DateTime, Datelike, Timelike, Utc};

// Mirror the orbital helper modules from verify_orbit_solver.rs so orbit_solver.rs compiles.
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

        let x_ecef = pos_teme[0] * cos_t + pos_teme[1] * sin_t;
        let y_ecef = -pos_teme[0] * sin_t + pos_teme[1] * cos_t;
        let z_ecef = pos_teme[2];

        let omega_e = 7.2921151467e-5;

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
        let a = 6378137.0;
        let f = 1.0 / 298.257223563;
        let e2 = 2.0 * f - f * f;

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

        let a = 6378137.0;
        let f = 1.0 / 298.257223563;
        let b = a * (1.0 - f);
        let e2 = (a * a - b * b) / (a * a);
        let ep2 = (a * a - b * b) / (b * b);

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

/// Shared helper: simulate passes from a virtual satellite and return (passes, rec_ecef).
/// Identical to the simulation logic in test_virtual_tcxo_polynomial_fit.
fn simulate_passes_with_drift(
    drift: (f64, f64, f64),
) -> (Vec<orbit_solver::RawPass>, [f64; 3]) {
    let truth_a = 6378137.0 + 550000.0;
    let truth_i = 53.0_f64.to_radians();
    let truth_raan = 1.2;
    let truth_u0 = 0.5;
    let center_freq = 150800000.0;

    let epoch = DateTime::parse_from_rfc3339("2026-06-09T12:00:00Z")
        .unwrap()
        .with_timezone(&Utc);

    let (pos_sat_epoch, _) =
        orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, epoch);

    let (rec_lat, rec_lon, _) = orbit::ecef_to_wgs84(pos_sat_epoch);
    let rec_ecef = orbit::wgs84_to_ecef(rec_lat, rec_lon, 0.0);

    let mut passes = Vec::new();
    let mut current_pass: Option<orbit_solver::RawPass> = None;
    let mut last_t: Option<DateTime<Utc>> = None;

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

    // Select best pair of passes (minimize synodic-period distortion)
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

    let mut selected_passes = vec![passes[best_pair.0].clone(), passes[best_pair.1].clone()];

    // Inject polynomial thermal drift: df(t) = drift.0 + drift.1 * tau + drift.2 * tau^2
    let (d0, d1, d2) = drift;
    for pass in &mut selected_passes {
        let t_ref = pass.points[0].time;
        for pt in &mut pass.points {
            let tau = (pt.time - t_ref).num_milliseconds() as f64 / 1000.0;
            let drift_val = d0 + d1 * tau + d2 * tau * tau;
            pt.freq += drift_val;
        }
    }

    (selected_passes, rec_ecef)
}

// ---------------------------------------------------------------------------
// Test 6: Zero-drift TCXO polynomial — all drift coefficients should converge near zero
// ---------------------------------------------------------------------------
#[test]
fn test_tcxo_polynomial_zero_drift() {
    let truth_a = 6378137.0 + 550000.0;
    let truth_i = 53.0_f64.to_radians();

    let (passes, rec_ecef) = simulate_passes_with_drift((0.0, 0.0, 0.0));

    assert!(
        passes.len() >= 2,
        "Need at least 2 simulated passes, got {}",
        passes.len()
    );

    let initial_a = truth_a + 500.0;
    let initial_i = truth_i + 0.01;

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, initial_a, initial_i)
        .expect("Zero-drift TCXO Doppler orbit solver failed to converge");

    println!(
        "Zero-drift solved: a={:.1}, i={:.4}, raan0={:.4}, u0={:.4}",
        solved.a, solved.i, solved.raan0, solved.u0
    );

    for idx in 0..2 {
        println!(
            "Pass {}: df0={:.6} Hz, df1={:.8} Hz/s, df2={:.10} Hz/s^2",
            idx, solved.pass_dfs[idx], solved.pass_df1s[idx], solved.pass_df2s[idx]
        );

        // With zero injected drift, the solver may absorb small orbital fitting residuals
        // into the polynomial coefficients. We verify the drift terms remain small — not
        // that they're exactly zero — confirming the solver doesn't hallucinate large drift.
        assert!(
            solved.pass_dfs[idx].abs() < 30.0,
            "Pass {} df0 should be small (no injected drift), got {:.6}",
            idx,
            solved.pass_dfs[idx]
        );
        assert!(
            solved.pass_df1s[idx].abs() < 0.1,
            "Pass {} df1 should be small (no injected drift), got {:.8}",
            idx,
            solved.pass_df1s[idx]
        );
        assert!(
            solved.pass_df2s[idx].abs() < 0.001,
            "Pass {} df2 should be small (no injected drift), got {:.10}",
            idx,
            solved.pass_df2s[idx]
        );
    }
}

// ---------------------------------------------------------------------------
// Test 7: Linear-only TCXO drift — df2 should converge near zero
// ---------------------------------------------------------------------------
#[test]
fn test_tcxo_polynomial_linear_only() {
    let truth_a = 6378137.0 + 550000.0;
    let truth_i = 53.0_f64.to_radians();

    let (passes, rec_ecef) = simulate_passes_with_drift((5.0, 0.05, 0.0));

    assert!(
        passes.len() >= 2,
        "Need at least 2 simulated passes, got {}",
        passes.len()
    );

    let initial_a = truth_a + 500.0;
    let initial_i = truth_i + 0.01;

    let solved = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, initial_a, initial_i)
        .expect("Linear-only TCXO Doppler orbit solver failed to converge");

    println!(
        "Linear-only drift solved: a={:.1}, i={:.4}, raan0={:.4}, u0={:.4}",
        solved.a, solved.i, solved.raan0, solved.u0
    );

    for idx in 0..2 {
        println!(
            "Pass {}: df0={:.4} Hz (target=5.0), df1={:.6} Hz/s (target=0.05), df2={:.8} Hz/s^2 (target=0.0)",
            idx, solved.pass_dfs[idx], solved.pass_df1s[idx], solved.pass_df2s[idx]
        );

        let err_df0 = (solved.pass_dfs[idx] - 5.0).abs();
        let err_df1 = (solved.pass_df1s[idx] - 0.05).abs();

        assert!(
            err_df0 < 2.0,
            "Pass {} df0 error too high: {:.6} (expected near 5.0)",
            idx,
            err_df0
        );
        assert!(
            err_df1 < 0.05,
            "Pass {} df1 error too high: {:.8} (expected near 0.05)",
            idx,
            err_df1
        );
        assert!(
            solved.pass_df2s[idx].abs() < 0.0005,
            "Pass {} df2 should be near zero (no quadratic drift), got {:.8}",
            idx,
            solved.pass_df2s[idx]
        );
    }
}

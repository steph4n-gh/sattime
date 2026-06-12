use chrono::{DateTime, Datelike, TimeZone, Timelike, Utc};

pub mod orbit {
    use super::*;
    pub fn datetime_to_jd(dt: DateTime<Utc>) -> f64 {
        let year = dt.year() as f64;
        let month = dt.month() as f64;
        let day = dt.day() as f64;
        let hour = dt.hour() as f64;
        let minute = dt.minute() as f64;
        let second = dt.second() as f64;
        let nanosecond = dt.nanosecond() as f64;

        let day_fraction = (hour + (minute + (second + nanosecond / 1e9) / 60.0) / 60.0) / 24.0;
        let jd_day = day + day_fraction;

        let (y, m) = if month <= 2.0 {
            (year - 1.0, month + 12.0)
        } else {
            (year, month)
        };

        let a = (y / 100.0).floor();
        let b = 2.0 - a + (a / 4.0).floor();

        (365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + jd_day + b - 1524.5
    }

    pub fn teme_to_ecef(jd: f64, pos_teme: [f64; 3], vel_teme: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let d = jd - 2451545.0;
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
}

#[path = "../src/orbit_solver.rs"]
mod orbit_solver;

// Structure to define a stress test configuration
struct StressTestConfig {
    name: &'static str,
    truth_alt: f64,      // altitude in meters
    truth_inc_deg: f64,  // inclination in degrees
    truth_raan_deg: f64, // RAAN in degrees
    truth_u0_deg: f64,   // u0 in degrees
    obs_lat_offset_deg: f64,
    obs_lon_offset_deg: f64,
    freq_noise_std: f64, // standard deviation in Hz
    points_per_pass: usize,
    initial_a_offset: f64,     // guess offset in meters
    initial_i_offset_deg: f64, // guess offset in degrees
}

fn run_single_stress_test(config: &StressTestConfig) -> Result<(f64, f64, f64, f64), String> {
    let truth_a = 6378137.0 + config.truth_alt;
    let truth_i = config.truth_inc_deg.to_radians();
    let truth_raan = config.truth_raan_deg.to_radians();
    let truth_u0 = config.truth_u0_deg.to_radians();
    let center_freq = 150800000.0;

    let epoch = Utc.with_ymd_and_hms(2026, 6, 9, 12, 0, 0).unwrap();

    // Determine observer ECEF
    let (pos_sat_epoch, _) =
        orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, epoch);
    let (rec_lat, rec_lon, _) = orbit::ecef_to_wgs84(pos_sat_epoch);
    let rec_ecef = orbit::wgs84_to_ecef(
        rec_lat + config.obs_lat_offset_deg,
        rec_lon + config.obs_lon_offset_deg,
        0.0,
    );

    // Simulate passes
    let mut passes = Vec::new();
    let mut current_pass: Option<orbit_solver::RawPass> = None;
    let mut last_t: Option<DateTime<Utc>> = None;

    let mut rng = orbit_solver::SimpleRng::new(42);

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

            // Add frequency noise
            let noise = if config.freq_noise_std > 0.0 {
                rng.next_gaussian() * config.freq_noise_std
            } else {
                0.0
            };
            let observed_freq = center_freq + pred_shift + noise;

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

    if passes.len() < 2 {
        return Err(format!("Not enough simulated passes: {}", passes.len()));
    }

    // Truncate pass points if needed to test sparse data
    for pass in &mut passes {
        if pass.points.len() > config.points_per_pass {
            let step_size = pass.points.len() / config.points_per_pass;
            let mut new_points = Vec::new();
            for idx in 0..config.points_per_pass {
                if idx * step_size < pass.points.len() {
                    new_points.push(pass.points[idx * step_size].clone());
                }
            }
            pass.points = new_points;
        }
    }

    let initial_a = truth_a + config.initial_a_offset;
    let initial_i = truth_i + config.initial_i_offset_deg.to_radians();

    let start_time = std::time::Instant::now();
    let solved_res = orbit_solver::fit_orbit_doppler(&passes, rec_ecef, initial_a, initial_i);
    let duration = start_time.elapsed().as_secs_f64();

    match solved_res {
        Ok(solved) => {
            let err_a = (solved.a - truth_a).abs();
            let err_i = (solved.i - truth_i).abs().to_degrees();
            let err_raan = (solved.raan0 - truth_raan).abs().to_degrees();
            let _err_u0 = (solved.u0 - truth_u0).abs().to_degrees();
            Ok((err_a, err_i, err_raan, duration))
        }
        Err(e) => Err(e.to_string()),
    }
}

#[test]
fn test_langevin_stress_suite() {
    let test_cases = vec![
        // 1. Happy path (similar to the standard test)
        StressTestConfig {
            name: "Happy Path",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 800000.0,
            initial_i_offset_deg: 10.0,
        },
        // 2. Low Altitude Orbit (300 km)
        StressTestConfig {
            name: "VLEO (300km) Altitude",
            truth_alt: 300000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 3. High Altitude Orbit (2000 km)
        StressTestConfig {
            name: "High LEO (2000km)",
            truth_alt: 2000000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 4. Equatorial Orbit (0 deg inclination)
        StressTestConfig {
            name: "Equatorial Orbit (i=0)",
            truth_alt: 550000.0,
            truth_inc_deg: 0.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 5. Polar Orbit (90 deg inclination)
        StressTestConfig {
            name: "Polar Orbit (i=90)",
            truth_alt: 550000.0,
            truth_inc_deg: 90.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 6. Retrograde Orbit (140 deg inclination)
        StressTestConfig {
            name: "Retrograde Orbit (i=140)",
            truth_alt: 550000.0,
            truth_inc_deg: 140.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 7. Modest Noise (10 Hz)
        StressTestConfig {
            name: "Modest Noise (10 Hz SD)",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 10.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 8. Heavy Noise (100 Hz)
        StressTestConfig {
            name: "Heavy Noise (100 Hz SD)",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 100.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 9. Extreme Noise (500 Hz)
        StressTestConfig {
            name: "Extreme Noise (500 Hz SD)",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 500.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 10. Sparse Points (only 10 points per pass)
        StressTestConfig {
            name: "Sparse Points (10 pts/pass)",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 10,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 11. Observer Displacement (5 degrees lat/lon offset)
        StressTestConfig {
            name: "Observer Displacement (5 deg)",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 5.0,
            obs_lon_offset_deg: 5.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 500000.0,
            initial_i_offset_deg: 5.0,
        },
        // 12. Extremely bad initial guess (1500 km altitude, 20 deg inclination offset)
        StressTestConfig {
            name: "Extreme Guess Offset",
            truth_alt: 550000.0,
            truth_inc_deg: 53.0,
            truth_raan_deg: 45.0,
            truth_u0_deg: 45.0,
            obs_lat_offset_deg: 0.0,
            obs_lon_offset_deg: 0.0,
            freq_noise_std: 0.0,
            points_per_pass: 100,
            initial_a_offset: 1500000.0,
            initial_i_offset_deg: 20.0,
        },
    ];

    println!("============================================================");
    println!("ADELIC LANGEVIN SOLVER STRESS TEST RESULTS");
    println!("============================================================");
    println!(
        "| {:<25} | {:<8} | {:<10} | {:<8} | {:<8} |",
        "Scenario Name", "Status", "a Err (m)", "i Err (deg)", "Time (s)"
    );
    println!("|---------------------------|----------|------------|------------|----------|");

    let mut successes = 0;
    for tc in &test_cases {
        match run_single_stress_test(tc) {
            Ok((err_a, err_i, _err_raan, duration)) => {
                let converged = err_a < 500.0 && err_i < 0.1;
                let status = if converged { "CONVERGED" } else { "FAILED" };
                if converged {
                    successes += 1;
                }
                println!(
                    "| {:<25} | {:<8} | {:<10.2} | {:<10.6} | {:<8.3} |",
                    tc.name, status, err_a, err_i, duration
                );
            }
            Err(e) => {
                println!(
                    "| {:<25} | {:<8} | {:<10} | {:<10} | {:<8} |",
                    tc.name, "ERROR", "N/A", "N/A", "N/A"
                );
                println!("  --> Error: {}", e);
            }
        }
    }
    println!("============================================================");
    println!("Passed: {} / {} scenarios", successes, test_cases.len());
    println!("============================================================");

    // We expect a robust solver to pass at least 8 of these 12 cases.
    assert!(
        successes >= 8,
        "Solver failed more stress test cases than expected"
    );
}

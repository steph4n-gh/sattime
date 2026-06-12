use chrono::{DateTime, Datelike, Timelike, Utc};

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

struct SimParams {
    truth_a: f64,
    truth_i: f64,
    truth_raan: f64,
    truth_u0: f64,
    rec_ecef: [f64; 3],
    center_freq: f64,
    epoch: DateTime<Utc>,
}

fn generate_base_sim() -> SimParams {
    let truth_a = 6378137.0 + 550000.0; // 550 km altitude
    let truth_i = 53.0_f64.to_radians();
    let truth_raan = 45.0_f64.to_radians();
    let truth_u0 = 45.0_f64.to_radians();
    let center_freq = 150800000.0;

    let epoch = DateTime::parse_from_rfc3339("2026-06-09T12:00:00Z")
        .unwrap()
        .with_timezone(&Utc);

    let (pos_sat_epoch, _) =
        orbit_solver::propagate_ecef(truth_a, truth_i, truth_raan, truth_u0, epoch, epoch);

    let (rec_lat, rec_lon, _) = orbit::ecef_to_wgs84(pos_sat_epoch);
    let rec_ecef = orbit::wgs84_to_ecef(rec_lat, rec_lon, 0.0);

    SimParams {
        truth_a,
        truth_i,
        truth_raan,
        truth_u0,
        rec_ecef,
        center_freq,
        epoch,
    }
}

fn simulate_passes_with_noise(
    sim: &SimParams,
    noise_std: f64,
    seed: u64,
) -> Vec<orbit_solver::RawPass> {
    let mut passes = Vec::new();
    let mut current_pass: Option<orbit_solver::RawPass> = None;
    let mut last_t: Option<DateTime<Utc>> = None;
    let mut rng = orbit_solver::SimpleRng::new(seed);

    // Propagate over 8 orbits
    for step in 0..9600 {
        let t = sim.epoch + chrono::Duration::seconds(step * 5);
        let (pos_sat, _) = orbit_solver::propagate_ecef(
            sim.truth_a,
            sim.truth_i,
            sim.truth_raan,
            sim.truth_u0,
            sim.epoch,
            t,
        );
        let dx = pos_sat[0] - sim.rec_ecef[0];
        let dy = pos_sat[1] - sim.rec_ecef[1];
        let dz = pos_sat[2] - sim.rec_ecef[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < 2500000.0 {
            let pred_shift = orbit_solver::predict_frequency(
                sim.truth_a,
                sim.truth_i,
                sim.truth_raan,
                sim.truth_u0,
                sim.epoch,
                t,
                0.0,
                0.0,
                sim.center_freq,
                sim.rec_ecef,
            );

            // Add Gaussian noise if requested
            let noise = if noise_std > 0.0 {
                rng.next_gaussian() * noise_std
            } else {
                0.0
            };

            let observed_freq = sim.center_freq + pred_shift + noise;

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
                            center_freq: sim.center_freq,
                            points: vec![point],
                        });
                    } else {
                        pass.points.push(point);
                    }
                }
            } else {
                current_pass = Some(orbit_solver::RawPass {
                    sat_name: "SIM_SAT".to_string(),
                    center_freq: sim.center_freq,
                    points: vec![point],
                });
            }
            last_t = Some(t);
        }
    }
    if let Some(pass) = current_pass {
        passes.push(pass);
    }
    passes
}

#[test]
fn test_stress_initial_guess_limits() {
    let sim = generate_base_sim();
    let passes = simulate_passes_with_noise(&sim, 0.0, 42);

    assert!(passes.len() >= 2);

    // Grid of offsets to test limits of Langevin Global solver
    let a_offsets: [f64; 6] = [
        -1500000.0, -1000000.0, -500000.0, 500000.0, 1000000.0, 1500000.0,
    ];
    let i_offsets: [f64; 6] = [-20.0, -10.0, -5.0, 5.0, 10.0, 20.0];

    for &a_off in &a_offsets {
        for &i_off in &i_offsets {
            let guess_a = sim.truth_a + a_off;
            let guess_i = sim.truth_i + i_off.to_radians();

            let solved_res =
                orbit_solver::fit_orbit_doppler(&passes, sim.rec_ecef, guess_a, guess_i);

            match solved_res {
                Ok(solved) => {
                    let err_a = (solved.a - sim.truth_a).abs();
                    let err_i = (solved.i - sim.truth_i).abs().to_degrees();

                    println!(
                        "GUESS OFFSET [a_off={:+e}m, i_off={:+.1} deg] -> Solved err_a={:.1}m, err_i={:.4} deg",
                        a_off, i_off, err_a, err_i
                    );

                    // We expect convergence within 1km and 0.2 deg for nominal offsets,
                    // but we want to see if the solver behaves stably (doesn't panic or output NaNs)
                    assert!(solved.a.is_finite());
                    assert!(solved.i.is_finite());
                    assert!(solved.raan0.is_finite());
                    assert!(solved.u0.is_finite());
                }
                Err(e) => {
                    println!(
                        "GUESS OFFSET [a_off={:+e}m, i_off={:+.1} deg] -> Solver failed with error: {}",
                        a_off, i_off, e
                    );
                }
            }
        }
    }
}

#[test]
fn test_stress_measurement_noise() {
    let sim = generate_base_sim();

    // Test noise standard deviations from 0.1 Hz up to 1000 Hz
    let noise_scenarios: [f64; 5] = [0.1, 1.0, 10.0, 100.0, 1000.0];

    for &noise in &noise_scenarios {
        let passes = simulate_passes_with_noise(&sim, noise, 999);
        let guess_a = sim.truth_a + 200000.0;
        let guess_i = sim.truth_i + 2.0_f64.to_radians();

        let solved_res = orbit_solver::fit_orbit_doppler(&passes, sim.rec_ecef, guess_a, guess_i);
        if noise <= 10.0 {
            let solved = solved_res.expect("Solver should converge for low noise");
            let err_a = (solved.a - sim.truth_a).abs();
            let err_i = (solved.i - sim.truth_i).abs().to_degrees();
            assert!(err_a < 1000.0, "a error too large for noise={}: {}", noise, err_a);
            assert!(err_i < 0.05, "i error too large for noise={}: {}", noise, err_i);
        } else {
            // For higher noise levels (100, 1000 Hz), if it succeeds, error must be finite, or it returns Err.
            if let Ok(solved) = solved_res {
                let err_a = (solved.a - sim.truth_a).abs();
                let err_i = (solved.i - sim.truth_i).abs().to_degrees();
                assert!(err_a.is_finite() && err_i.is_finite());
            }
        }
    }
}

#[test]
fn test_stress_receiver_position_offset() {
    let sim = generate_base_sim();
    let passes = simulate_passes_with_noise(&sim, 0.0, 12345);

    // Shift receiver position to simulate GPS/receiver survey errors
    let offsets: [f64; 4] = [10.0, 100.0, 1000.0, 10000.0]; // in meters

    for &offset in &offsets {
        // Shift ECEF coordinate on X axis
        let mut shifted_ecef = sim.rec_ecef;
        shifted_ecef[0] += offset;

        let guess_a = sim.truth_a + 200000.0;
        let guess_i = sim.truth_i + 2.0_f64.to_radians();

        let solved_res = orbit_solver::fit_orbit_doppler(&passes, shifted_ecef, guess_a, guess_i);
        if offset <= 100.0 {
            let solved = solved_res.expect("Solver should converge for low receiver shift");
            let err_a = (solved.a - sim.truth_a).abs();
            let err_i = (solved.i - sim.truth_i).abs().to_degrees();
            assert!(err_a < 1000.0, "a error too large for offset={}: {}", offset, err_a);
            assert!(err_i < 0.05, "i error too large for offset={}: {}", offset, err_i);
        } else {
            // For larger offsets, if it succeeds, error must be finite, or it returns Err.
            if let Ok(solved) = solved_res {
                let err_a = (solved.a - sim.truth_a).abs();
                let err_i = (solved.i - sim.truth_i).abs().to_degrees();
                assert!(err_a.is_finite() && err_i.is_finite());
            }
        }
    }
}

#[test]
fn test_stress_extreme_inputs() {
    let sim = generate_base_sim();
    let passes = simulate_passes_with_noise(&sim, 0.0, 111);

    // Test extreme semi-major axis (e.g. extremely small or near earth radius)
    let solved_small = orbit_solver::fit_orbit_doppler(&passes, sim.rec_ecef, 1000.0, sim.truth_i);
    assert!(solved_small.is_err(), "Expected error for extremely small semi-major axis guess");

    // Test extreme large guess
    let solved_large = orbit_solver::fit_orbit_doppler(&passes, sim.rec_ecef, 1e12, sim.truth_i);
    assert!(solved_large.is_err(), "Expected error for extremely large semi-major axis guess");

    // Test invalid/empty passes
    let empty_passes: Vec<orbit_solver::RawPass> = vec![];
    let solved_empty =
        orbit_solver::fit_orbit_doppler(&empty_passes, sim.rec_ecef, sim.truth_a, sim.truth_i);
    assert!(solved_empty.is_err());

    // Test single pass (should return error)
    if passes.len() >= 1 {
        let single_pass = vec![passes[0].clone()];
        let solved_single =
            orbit_solver::fit_orbit_doppler(&single_pass, sim.rec_ecef, sim.truth_a, sim.truth_i);
        assert!(solved_single.is_err());
    }
}

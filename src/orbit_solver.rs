use crate::orbit::{datetime_to_jd, teme_to_ecef};
use chrono::{DateTime, Datelike, Timelike, Utc};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead};

const MU: f64 = 3.986004418e14; // m^3/s^2
const J2: f64 = 1.0826263e-3;
const RE: f64 = 6378137.0; // m
const C: f64 = 299792458.0; // m/s

#[derive(Debug, Clone)]
pub struct PassPoint {
    pub time: DateTime<Utc>,
    pub freq: f64,
}

#[derive(Debug, Clone)]
pub struct RawPass {
    #[allow(dead_code)]
    pub sat_name: String,
    pub center_freq: f64,
    pub points: Vec<PassPoint>,
}

#[derive(Debug, Clone)]
pub struct SolvedOrbit {
    pub a: f64,     // Semi-major axis in meters
    pub i: f64,     // Inclination in radians
    pub raan0: f64, // RAAN at epoch in radians
    pub u0: f64,    // Argument of latitude at epoch in radians
    #[allow(dead_code)]
    pub epoch: DateTime<Utc>,
    #[allow(dead_code)]
    pub pass_dts: Vec<f64>,
    #[allow(dead_code)]
    pub pass_dfs: Vec<f64>,
    #[allow(dead_code)]
    pub pass_df1s: Vec<f64>,
    #[allow(dead_code)]
    pub pass_df2s: Vec<f64>,
}

#[allow(dead_code)]
pub fn read_pass_file(path: &str) -> io::Result<RawPass> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);

    let mut sat_name = String::from("UNKNOWN");
    let mut center_freq = 150800000.0;
    let mut points = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if trimmed.starts_with('#') {
            // Parse headers
            let parts: Vec<&str> = trimmed[1..].splitn(2, ':').collect();
            if parts.len() == 2 {
                let key = parts[0].trim().to_lowercase();
                let val = parts[1].trim();
                if key == "satellite_name" {
                    sat_name = val.to_string();
                } else if key == "center_frequency"
                    && let Ok(f) = val.parse::<f64>()
                {
                    center_freq = f;
                }
            }
            continue;
        }

        // Parse CSV row
        let cols: Vec<&str> = trimmed.split(',').collect();
        if cols.len() >= 2 {
            if cols[0] == "utc_timestamp" {
                continue; // Skip header row
            }
            if let (Ok(time), Ok(freq)) = (
                DateTime::parse_from_rfc3339(cols[0]),
                cols[1].parse::<f64>(),
            ) {
                points.push(PassPoint {
                    time: time.with_timezone(&Utc),
                    freq,
                });
            }
        }
    }

    Ok(RawPass {
        sat_name,
        center_freq,
        points,
    })
}

pub fn propagate_ecef_at_time(
    a: f64,
    i: f64,
    raan0: f64,
    u0: f64,
    epoch: DateTime<Utc>,
    t_orbit: DateTime<Utc>,
    t_rotate: DateTime<Utc>,
) -> ([f64; 3], [f64; 3]) {
    let tau = (t_orbit - epoch).num_microseconds().unwrap_or(0) as f64 / 1_000_000.0;

    let n = (MU / a.powi(3)).sqrt();
    let v = (MU / a).sqrt();

    // Nodal precession rate due to J2
    let cos_i = i.cos();
    let d_raan = -1.5 * J2 * (RE / a).powi(2) * n * cos_i;

    let raan = raan0 + d_raan * tau;
    let u = u0 + n * tau;

    // Position and velocity in orbital plane
    let x_plane = a * u.cos();
    let y_plane = a * u.sin();

    let vx_plane = -v * u.sin();
    let vy_plane = v * u.cos();

    let cos_o = raan.cos();
    let sin_o = raan.sin();
    let cos_i = i.cos();
    let sin_i = i.sin();

    // Rotate from orbital plane to TEME ECI frame
    let pos_teme = [
        x_plane * cos_o - y_plane * cos_i * sin_o,
        x_plane * sin_o + y_plane * cos_i * cos_o,
        y_plane * sin_i,
    ];

    let vel_teme = [
        vx_plane * cos_o - vy_plane * cos_i * sin_o,
        vx_plane * sin_o + vy_plane * cos_i * cos_o,
        vy_plane * sin_i,
    ];

    let jd = datetime_to_jd(t_rotate);
    teme_to_ecef(jd, pos_teme, vel_teme)
}

// Predict ECEF position and velocity for a circular orbit at delta t from epoch
pub fn propagate_ecef(
    a: f64,
    i: f64,
    raan0: f64,
    u0: f64,
    epoch: DateTime<Utc>,
    t_utc: DateTime<Utc>,
) -> ([f64; 3], [f64; 3]) {
    propagate_ecef_at_time(a, i, raan0, u0, epoch, t_utc, t_utc)
}

// Find the predicted times of closest approach to the receiver for all passes
// Find the predicted times of closest approach to the receiver for all passes
pub fn get_pred_pca_times(
    a: f64,
    i: f64,
    raan: f64,
    u0: f64,
    epoch: DateTime<Utc>,
    rec_ecef: [f64; 3],
    raw_passes: &[RawPass],
) -> Vec<DateTime<Utc>> {
    if raw_passes.is_empty() {
        return Vec::new();
    }

    let center_freq = raw_passes[0].center_freq;
    let mut pred_pcas = Vec::new();

    for pass in raw_passes {
        if pass.points.is_empty() {
            pred_pcas.push(epoch);
            continue;
        }

        // Find the observed PCA time for this pass
        let mut obs_pca_time = pass.points[0].time;
        let mut min_offset = f64::MAX;
        for pt in &pass.points {
            let off = (pt.freq - center_freq).abs();
            if off < min_offset {
                min_offset = off;
                obs_pca_time = pt.time;
            }
        }

        // Search for predicted closest approach around obs_pca_time (+/- 600s)
        let mut min_dist = f64::MAX;
        let mut best_t = obs_pca_time;

        // Coarse search: 20s steps
        let search_start = obs_pca_time - chrono::Duration::seconds(600);
        let search_end = obs_pca_time + chrono::Duration::seconds(600);
        let mut t = search_start;
        while t <= search_end {
            let (pos, _) = propagate_ecef(a, i, raan, u0, epoch, t);
            let dx = pos[0] - rec_ecef[0];
            let dy = pos[1] - rec_ecef[1];
            let dz = pos[2] - rec_ecef[2];
            let d = dx * dx + dy * dy + dz * dz;
            if d < min_dist {
                min_dist = d;
                best_t = t;
            }
            t = t + chrono::Duration::seconds(20);
        }

        // Medium search: 2s steps around best_t (+/- 20s)
        let med_start = best_t - chrono::Duration::seconds(20);
        let med_end = best_t + chrono::Duration::seconds(20);
        let mut t = med_start;
        while t <= med_end {
            let (pos, _) = propagate_ecef(a, i, raan, u0, epoch, t);
            let dx = pos[0] - rec_ecef[0];
            let dy = pos[1] - rec_ecef[1];
            let dz = pos[2] - rec_ecef[2];
            let d = dx * dx + dy * dy + dz * dz;
            if d < min_dist {
                min_dist = d;
                best_t = t;
            }
            t = t + chrono::Duration::seconds(2);
        }

        // Fine search: 1s steps around best_t (+/- 2s)
        let fine_start = best_t - chrono::Duration::seconds(2);
        let fine_end = best_t + chrono::Duration::seconds(2);
        let mut t = fine_start;
        while t <= fine_end {
            let (pos, _) = propagate_ecef(a, i, raan, u0, epoch, t);
            let dx = pos[0] - rec_ecef[0];
            let dy = pos[1] - rec_ecef[1];
            let dz = pos[2] - rec_ecef[2];
            let d = dx * dx + dy * dy + dz * dz;
            if d < min_dist {
                min_dist = d;
                best_t = t;
            }
            t = t + chrono::Duration::seconds(1);
        }

        pred_pcas.push(best_t);
    }

    pred_pcas
}

// Compute the predicted Doppler frequency shift
pub fn predict_frequency(
    a: f64,
    i: f64,
    raan0: f64,
    u0: f64,
    epoch: DateTime<Utc>,
    t_obs: DateTime<Utc>,
    dt: f64,
    df: f64,
    center_freq: f64,
    rec_ecef: [f64; 3],
) -> f64 {
    let t_adj = t_obs + chrono::Duration::microseconds((dt * 1_000_000.0) as i64);
    
    // Estimate initial position/velocity for eccentricity correction
    let (pos_init, vel_init) = propagate_ecef_at_time(a, i, raan0, u0, epoch, t_adj, t_obs);
    let r_dot_v = pos_init[0] * vel_init[0] + pos_init[1] * vel_init[1] + pos_init[2] * vel_init[2];
    let dt_rel = -2.0 * r_dot_v / (C * C);
    
    let t_corr = t_adj + chrono::Duration::nanoseconds((dt_rel * 1e9) as i64);
    let (pos_sat, vel_sat) = propagate_ecef_at_time(a, i, raan0, u0, epoch, t_corr, t_obs);

    let dx = pos_sat[0] - rec_ecef[0];
    let dy = pos_sat[1] - rec_ecef[1];
    let dz = pos_sat[2] - rec_ecef[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

    if dist < 1.0 {
        return center_freq + df;
    }

    let v_sat_sq = vel_sat[0] * vel_sat[0] + vel_sat[1] * vel_sat[1] + vel_sat[2] * vel_sat[2];
    let r_sat = (pos_sat[0] * pos_sat[0] + pos_sat[1] * pos_sat[1] + pos_sat[2] * pos_sat[2]).sqrt();
    let r_rec = (rec_ecef[0] * rec_ecef[0] + rec_ecef[1] * rec_ecef[1] + rec_ecef[2] * rec_ecef[2]).sqrt();

    let u_sat = -MU / r_sat;
    let u_rec = -MU / r_rec;

    let v_dot_n = (vel_sat[0] * dx + vel_sat[1] * dy + vel_sat[2] * dz) / dist;

    let gamma_inv = (1.0 - v_sat_sq / (C * C)).sqrt();
    let denominator = 1.0 - v_dot_n / C;
    let potential_term = 1.0 + (u_sat - u_rec) / (C * C);

    let f_obs = center_freq * (gamma_inv / denominator) * potential_term;
    f_obs - center_freq + df
}

pub fn predict_frequency_poly(
    a: f64,
    i: f64,
    raan0: f64,
    u0: f64,
    epoch: DateTime<Utc>,
    t_obs: DateTime<Utc>,
    dt: f64,
    df_poly: (f64, f64, f64),
    t_ref: DateTime<Utc>,
    center_freq: f64,
    rec_ecef: [f64; 3],
) -> f64 {
    let t_adj = t_obs + chrono::Duration::microseconds((dt * 1_000_000.0) as i64);
    
    // Estimate initial position/velocity for eccentricity correction
    let (pos_init, vel_init) = propagate_ecef_at_time(a, i, raan0, u0, epoch, t_adj, t_obs);
    let r_dot_v = pos_init[0] * vel_init[0] + pos_init[1] * vel_init[1] + pos_init[2] * vel_init[2];
    let dt_rel = -2.0 * r_dot_v / (C * C);
    
    let t_corr = t_adj + chrono::Duration::nanoseconds((dt_rel * 1e9) as i64);
    let (pos_sat, vel_sat) = propagate_ecef_at_time(a, i, raan0, u0, epoch, t_corr, t_obs);

    let dx = pos_sat[0] - rec_ecef[0];
    let dy = pos_sat[1] - rec_ecef[1];
    let dz = pos_sat[2] - rec_ecef[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

    if dist < 1.0 {
        return center_freq + df_poly.0;
    }

    let v_sat_sq = vel_sat[0] * vel_sat[0] + vel_sat[1] * vel_sat[1] + vel_sat[2] * vel_sat[2];
    let r_sat = (pos_sat[0] * pos_sat[0] + pos_sat[1] * pos_sat[1] + pos_sat[2] * pos_sat[2]).sqrt();
    let r_rec = (rec_ecef[0] * rec_ecef[0] + rec_ecef[1] * rec_ecef[1] + rec_ecef[2] * rec_ecef[2]).sqrt();

    let u_sat = -MU / r_sat;
    let u_rec = -MU / r_rec;

    let v_dot_n = (vel_sat[0] * dx + vel_sat[1] * dy + vel_sat[2] * dz) / dist;

    let gamma_inv = (1.0 - v_sat_sq / (C * C)).sqrt();
    let denominator = 1.0 - v_dot_n / C;
    let potential_term = 1.0 + (u_sat - u_rec) / (C * C);

    let f_obs = center_freq * (gamma_inv / denominator) * potential_term;
    let tau = (t_obs - t_ref).num_milliseconds() as f64 / 1000.0;
    let bias = df_poly.0 + df_poly.1 * tau + df_poly.2 * tau * tau;
    f_obs - center_freq + bias
}


// Formulate the Keplerian parameters Levenberg-Marquardt solver
pub fn fit_orbit_doppler(
    raw_passes: &[RawPass],
    rec_ecef: [f64; 3],
    initial_a: f64,
    initial_i: f64,
) -> Result<SolvedOrbit, Box<dyn std::error::Error>> {
    if rec_ecef.iter().any(|&x| x.is_nan()) || initial_a.is_nan() || initial_i.is_nan() {
        return Err("Input parameters contain NaN".into());
    }
    if initial_a < 6000e3 || initial_a > 50000e3 || initial_i < 0.0 || initial_i > std::f64::consts::PI {
        return Err("Input parameters are out of realistic physical bounds".into());
    }
    for p in raw_passes {
        for pt in &p.points {
            if pt.freq.is_nan() {
                return Err("Observation frequency contains NaN".into());
            }
        }
    }

    let raw_passes_filtered: Vec<RawPass> = raw_passes
        .iter()
        .filter(|p| !p.points.is_empty())
        .cloned()
        .collect();
    if raw_passes_filtered.len() < 2 {
        return Err("Need at least 2 non-empty passes to resolve orbit parameters.".into());
    }
    let raw_passes = &raw_passes_filtered;

    let normalized_coords = map_to_normalized_search_space(rec_ecef);
    println!(
        "Normalized search space coordinates: {:?}",
        normalized_coords
    );

    let n_passes = raw_passes.len();

    let epoch = raw_passes[0].points[0].time;
    let center_freq = raw_passes[0].center_freq;

    // We fit: [a, i, raan0, u0, dt_0, df0_0, df1_0, df2_0, dt_1, df0_1, df1_1, df2_1, ...]
    // Total parameters: 4 + 4 * n_passes
    let n_params = 4 + 4 * n_passes;
    let mut params = vec![0.0; n_params];

    // Keplerian Period Estimation from observed PCA times
    let mut obs_pcas = Vec::new();
    for pass in raw_passes {
        if pass.points.is_empty() {
            continue;
        }
        let mut obs_pca_time = pass.points[0].time;
        let mut min_offset = f64::MAX;
        for pt in &pass.points {
            let off = (pt.freq - center_freq).abs();
            if off < min_offset {
                min_offset = off;
                obs_pca_time = pt.time;
            }
        }
        obs_pcas.push(obs_pca_time);
    }

    let a_est = if obs_pcas.len() >= 2 {
        let mut t_obs_sum = 0.0;
        let mut t_obs_count = 0.0;
        let t0 = obs_pcas[0];
        for &tp in &obs_pcas[1..] {
            let dt = (tp - t0).num_milliseconds() as f64 / 1000.0;
            let k = (dt / 5700.0).round();
            println!("DEBUG: dt={}, k={}, dt/k={}", dt, k, dt / k);
            if k > 0.0 {
                t_obs_sum += dt / k;
                t_obs_count += 1.0;
            }
        }
        let t_observed = if t_obs_count > 0.0 {
            t_obs_sum / t_obs_count
        } else {
            5700.0
        };
        let a_val = (MU * t_observed * t_observed
            / (4.0 * std::f64::consts::PI * std::f64::consts::PI))
            .powf(1.0 / 3.0);
        println!("DEBUG: t_observed={}, a_est={}", t_observed, a_val);
        a_val.max(6500e3).min(20000e3)
    } else {
        initial_a
    };

    // --- STAGE 1: Fit using only the first 2 passes to get close to true a and i ---
    let stage1_passes = &raw_passes[0..2];
    let mut stage1_params = vec![0.0; 4 + 2 * 2]; // 4 global + 2 * 2 pass-specific = 8 params
    stage1_params[0] = a_est;
    stage1_params[1] = initial_i;

    // Run Langevin Global Optimizer on the first 2 passes to get raan0 and u0
    let mut best_raan0 = 0.0_f64;
    let mut best_u0 = 0.0_f64;
    let mut best_rss = f64::MAX;


    // 12x12 grid of starting points for Langevin trajectories (30 degree spacing)
    let mut grid_points = [0.0; 12];
    for idx in 0..12 {
        grid_points[idx] = (idx as f64) * 30.0_f64.to_radians();
    }

    let mut starts = Vec::new();
    let mut rng = SimpleRng::new(1337);
    for &init_raan in &grid_points {
        for &init_u0 in &grid_points {
            starts.push((init_raan, init_u0, rng.state));
            for _ in 0..75 {
                rng.next_f64();
            }
        }
    }

    let num_threads = (rayon::current_num_threads() - 2).max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let results: Vec<((f64, f64), f64, Vec<f64>)> = pool.install(|| {
        starts
            .into_par_iter()
            .map(|(init_raan, init_u0, seed_state)| {
                let mut raan0 = init_raan;
                let mut u0 = init_u0;
                let mut traj_best_raan = raan0;
                let mut traj_best_u = u0;
                let mut traj_best_rss = f64::MAX;

                let mut rng = SimpleRng { state: seed_state };

                let mut lr = 0.1;
                let mut noise_std = 0.05;

                let mut rss_history = Vec::new();

                for step in 0..15 {
                    let current_rss = compute_rss(
                        stage1_passes,
                        rec_ecef,
                        stage1_params[0],
                        stage1_params[1],
                        epoch,
                        center_freq,
                        raan0,
                        u0,
                    );
                    rss_history.push(current_rss);
                    if current_rss < traj_best_rss {
                        traj_best_rss = current_rss;
                        traj_best_raan = raan0;
                        traj_best_u = u0;
                    }

                    // Compute gradient
                    let (_, _, grad_raan, grad_u0) = compute_gradient(
                        stage1_passes,
                        rec_ecef,
                        stage1_params[0],
                        stage1_params[1],
                        epoch,
                        center_freq,
                        raan0,
                        u0,
                    );

                    let sign_raan = if grad_raan.is_nan() {
                        0.0
                    } else {
                        grad_raan.signum()
                    };
                    let sign_u0 = if grad_u0.is_nan() {
                        0.0
                    } else {
                        grad_u0.signum()
                    };

                    // Update directions
                    let step_raan = -lr * sign_raan * 0.2 + noise_std * rng.next_gaussian() * 0.05;
                    let step_u0 = -lr * sign_u0 * 0.2 + noise_std * rng.next_gaussian() * 0.05;

                    let next_raan0 = (raan0 + step_raan).rem_euclid(2.0 * std::f64::consts::PI);
                    let next_u0 = (u0 + step_u0).rem_euclid(2.0 * std::f64::consts::PI);

                    raan0 = next_raan0;
                    u0 = next_u0;

                    // Digit-scrambling restart using base-p digit reversal mapping (Monna map)
                    let primes = [2, 3, 5, 7];
                    let p = primes[step % primes.len()];

                    let x_raan = raan0 / (2.0 * std::f64::consts::PI);
                    let x_u = u0 / (2.0 * std::f64::consts::PI);

                    let val_raan = inverse_monna_map(x_raan, p, 16);
                    let val_u = inverse_monna_map(x_u, p, 16);

                    let perturb_scale = 4;
                    let perturbation = (rng.next_f64() * (p as f64).powi(perturb_scale)) as u64;
                    let val_raan_perturbed = val_raan.wrapping_add(perturbation);
                    let val_u_perturbed = val_u.wrapping_add(perturbation);

                    let x_raan_scrambled = monna_map(val_raan_perturbed, p);
                    let x_u_scrambled = monna_map(val_u_perturbed, p);

                    let raan0_scrambled = (x_raan_scrambled * 2.0 * std::f64::consts::PI)
                        .rem_euclid(2.0 * std::f64::consts::PI);
                    let u0_scrambled = (x_u_scrambled * 2.0 * std::f64::consts::PI)
                        .rem_euclid(2.0 * std::f64::consts::PI);

                    let scrambled_rss = compute_rss(
                        stage1_passes,
                        rec_ecef,
                        stage1_params[0],
                        stage1_params[1],
                        epoch,
                        center_freq,
                        raan0_scrambled,
                        u0_scrambled,
                    );
                    let current_diffs = compute_fractional_difference_history(&rss_history);
                    let current_var_sum: f64 = current_diffs.iter().map(|d| d.abs()).sum();
                    let reg_rss_current = current_rss + 0.01 * current_var_sum;

                    let mut temp_history = rss_history.clone();
                    if let Some(last_elem) = temp_history.last_mut() {
                        *last_elem = scrambled_rss;
                    }
                    let scrambled_diffs = compute_fractional_difference_history(&temp_history);
                    let scrambled_var_sum: f64 = scrambled_diffs.iter().map(|d| d.abs()).sum();
                    let reg_rss_scrambled = scrambled_rss + 0.01 * scrambled_var_sum;

                    if reg_rss_scrambled < reg_rss_current {
                        raan0 = raan0_scrambled;
                        u0 = u0_scrambled;
                        if let Some(last_elem) = rss_history.last_mut() {
                            *last_elem = scrambled_rss;
                        }
                    }

                    lr *= 0.95;
                    noise_std *= 0.9;
                }

                let final_rss = compute_rss(
                    stage1_passes,
                    rec_ecef,
                    stage1_params[0],
                    stage1_params[1],
                    epoch,
                    center_freq,
                    raan0,
                    u0,
                );
                if final_rss < traj_best_rss {
                    traj_best_rss = final_rss;
                    traj_best_raan = raan0;
                    traj_best_u = u0;
                }

                ((traj_best_raan, traj_best_u), traj_best_rss, rss_history)
            })
            .collect()
    });

    let mut best_rss_history = Vec::new();
    for ((traj_best_raan, traj_best_u), traj_best_rss, traj_history) in results {
        if traj_best_rss < best_rss {
            best_rss = traj_best_rss;
            best_raan0 = traj_best_raan;
            best_u0 = traj_best_u;
            best_rss_history = traj_history;
        }
    }

    let rss_diffs = compute_fractional_difference_history(&best_rss_history);
    if !rss_diffs.is_empty() {
        println!(
            "Langevin trajectory fractional variation sum: {:?}",
            rss_diffs.iter().sum::<f64>()
        );
    }

    println!(
        "Langevin best: raan0={:.4} deg, u0={:.4} deg, rss={:.2e}",
        best_raan0.to_degrees(),
        best_u0.to_degrees(),
        best_rss
    );
    stage1_params[2] = best_raan0;
    stage1_params[3] = best_u0;

    // Initialize stage 1 pass-specific parameters
    let stage1_pred_pcas = get_pred_pca_times(
        stage1_params[0],
        stage1_params[1],
        stage1_params[2],
        stage1_params[3],
        epoch,
        rec_ecef,
        stage1_passes,
    );
    for (p_idx, pass) in stage1_passes.iter().enumerate() {
        let mut obs_pca_time = pass.points[0].time;
        let mut min_offset = f64::MAX;
        for pt in &pass.points {
            let off = (pt.freq - center_freq).abs();
            if off < min_offset {
                min_offset = off;
                obs_pca_time = pt.time;
            }
        }
        let pred_pca_time = stage1_pred_pcas[p_idx];
        let dt = (pred_pca_time - obs_pca_time).num_milliseconds() as f64 / 1000.0;
        stage1_params[4 + 2 * p_idx] = dt;
        stage1_params[4 + 2 * p_idx + 1] = 0.0;
    }

    // Run LM on Stage 1 (optimize a, i, raan0, u0 using only the first 2 passes)
    let mut stage1_lambda = 1.0;
    let mut best_stage1_rss = f64::MAX;
    let mut best_stage1_params = stage1_params.clone();

    for _ in 0..100 {
        let mut residuals = Vec::new();
        for (p_idx, pass) in stage1_passes.iter().enumerate() {
            let dt = stage1_params[4 + 2 * p_idx];
            let df = stage1_params[4 + 2 * p_idx + 1];
            for pt in &pass.points {
                let pred = predict_frequency(
                    stage1_params[0],
                    stage1_params[1],
                    stage1_params[2],
                    stage1_params[3],
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
        for p_idx in 0..stage1_passes.len() {
            let dt = stage1_params[4 + 2 * p_idx];
            residuals.push(dt * 10.0);
        }

        let rss: f64 = residuals.iter().map(|r| r * r).sum();
        if rss < best_stage1_rss {
            best_stage1_rss = rss;
            best_stage1_params = stage1_params.clone();
            stage1_lambda /= 10.0;
        } else {
            stage1_params = best_stage1_params.clone();
            stage1_lambda *= 10.0;
            if stage1_lambda > 1e12 {
                break;
            }
            residuals.clear();
            for (p_idx, pass) in stage1_passes.iter().enumerate() {
                let dt = stage1_params[4 + 2 * p_idx];
                let df = stage1_params[4 + 2 * p_idx + 1];
                for pt in &pass.points {
                    let pred = predict_frequency(
                        stage1_params[0],
                        stage1_params[1],
                        stage1_params[2],
                        stage1_params[3],
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
            for p_idx in 0..stage1_passes.len() {
                let dt = stage1_params[4 + 2 * p_idx];
                residuals.push(dt * 10.0);
            }
        }

        let n_obs = residuals.len();
        if n_obs < 8 {
            break;
        }

        let mut jacobian = vec![vec![0.0; 8]; n_obs];
        for k in 0..8 {
            let mut perturbed = stage1_params.clone();
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
            for (p_idx, pass) in stage1_passes.iter().enumerate() {
                let dt = perturbed[4 + 2 * p_idx];
                let df = perturbed[4 + 2 * p_idx + 1];
                for pt in &pass.points {
                    let pred = predict_frequency(
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
            for p_idx in 0..stage1_passes.len() {
                let dt = perturbed[4 + 2 * p_idx];
                let diff = dt * 10.0;
                jacobian[row_idx][k] = (diff - residuals[row_idx]) / param_eps;
                row_idx += 1;
            }
        }

        let mut jt_j = vec![vec![0.0; 8]; 8];
        let mut jt_r = vec![0.0; 8];
        for row in 0..n_obs {
            for c1 in 0..8 {
                jt_r[c1] += jacobian[row][c1] * residuals[row];
                for c2 in 0..8 {
                    jt_j[c1][c2] += jacobian[row][c1] * jacobian[row][c2];
                }
            }
        }

        for k in 0..8 {
            jt_j[k][k] += stage1_lambda * jt_j[k][k];
        }

        if let Some(delta) = solve_linear_system(&mut jt_j, &jt_r) {
            for k in 0..8 {
                stage1_params[k] -= delta[k];
            }
            stage1_params[0] = stage1_params[0].max(6500e3).min(20000e3);
            stage1_params[1] = stage1_params[1].max(0.0).min(std::f64::consts::PI);
            // Scale-aware convergence: check relative step size per parameter.
            let max_rel_step = (0..8).map(|k| {
                let denom = stage1_params[k].abs().max(1e-10);
                delta[k].abs() / denom
            }).fold(0.0f64, f64::max);
            if max_rel_step < 1e-8 {
                break;
            }
        } else {
            break;
        }
    }
    stage1_params = best_stage1_params;

    // --- STAGE 2: Fit using all passes, initialized with Stage 1 refined parameters ---
    params[0] = stage1_params[0];
    params[1] = stage1_params[1];
    params[2] = stage1_params[2];
    params[3] = stage1_params[3];

    // Initialize pass-specific parameters (dt_j and df_j) for all passes
    let pred_pcas = get_pred_pca_times(
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
        let dt = (pred_pca_time - obs_pca_time).num_milliseconds() as f64 / 1000.0;
        params[4 + 4 * p_idx] = dt;
        params[4 + 4 * p_idx + 1] = 0.0;
        params[4 + 4 * p_idx + 2] = 0.0;
        params[4 + 4 * p_idx + 3] = 0.0;
    }

    // Run Levenberg-Marquardt (Stage 2)
    let mut lambda = 1.0;
    let mut best_rss = f64::MAX;
    let mut best_params = params.clone();

    for _ in 0..150 {
        // Compute residuals
        let mut residuals = Vec::new();
        for (p_idx, pass) in raw_passes.iter().enumerate() {
            let dt = params[4 + 4 * p_idx];
            let df0 = params[4 + 4 * p_idx + 1];
            let df1 = params[4 + 4 * p_idx + 2];
            let df2 = params[4 + 4 * p_idx + 3];
            let t_ref = pass.points[0].time;
            for pt in &pass.points {
                let pred = predict_frequency_poly(
                    params[0],
                    params[1],
                    params[2],
                    params[3],
                    epoch,
                    pt.time,
                    dt,
                    (df0, df1, df2),
                    t_ref,
                    center_freq,
                    rec_ecef,
                );
                residuals.push(pt.freq - (center_freq + pred));
            }
        }
        for p_idx in 0..raw_passes.len() {
            let dt = params[4 + 4 * p_idx];
            residuals.push(dt * 10.0);
        }
        // Regularize df1 and df2 toward zero to prevent overfitting on short passes.
        for p_idx in 0..raw_passes.len() {
            let df1 = params[4 + 4 * p_idx + 2];
            let df2 = params[4 + 4 * p_idx + 3];
            residuals.push(df1 * 100.0);    // Penalize linear drift
            residuals.push(df2 * 10000.0);  // Penalize quadratic drift more strongly
        }

        let rss: f64 = residuals.iter().map(|r| r * r).sum();
        if rss < best_rss {
            best_rss = rss;
            best_params = params.clone();
            lambda /= 10.0;
        } else {
            // Residuals increased, reject step and increase damping
            params = best_params.clone();
            lambda *= 10.0;
            if lambda > 1e12 {
                break;
            }
            // Recompute residuals for the restored params (best_params)
            residuals.clear();
            for (p_idx, pass) in raw_passes.iter().enumerate() {
                let dt = params[4 + 4 * p_idx];
                let df0 = params[4 + 4 * p_idx + 1];
                let df1 = params[4 + 4 * p_idx + 2];
                let df2 = params[4 + 4 * p_idx + 3];
                let t_ref = pass.points[0].time;
                for pt in &pass.points {
                    let pred = predict_frequency_poly(
                        params[0],
                        params[1],
                        params[2],
                        params[3],
                        epoch,
                        pt.time,
                        dt,
                        (df0, df1, df2),
                        t_ref,
                        center_freq,
                        rec_ecef,
                    );
                    residuals.push(pt.freq - (center_freq + pred));
                }
            }
            for p_idx in 0..raw_passes.len() {
                let dt = params[4 + 4 * p_idx];
                residuals.push(dt * 10.0);
            }
            // Regularize df1 and df2 toward zero to prevent overfitting on short passes.
            for p_idx in 0..raw_passes.len() {
                let df1 = params[4 + 4 * p_idx + 2];
                let df2 = params[4 + 4 * p_idx + 3];
                residuals.push(df1 * 100.0);    // Penalize linear drift
                residuals.push(df2 * 10000.0);  // Penalize quadratic drift more strongly
            }
        }

        let n_obs = residuals.len();
        if n_obs < n_params {
            return Err("Not enough observations to fit orbit parameters.".into());
        }

        // Compute Jacobian numerically
        let mut jacobian = vec![vec![0.0; n_params]; n_obs];

        for k in 0..n_params {
            let mut perturbed = params.clone();
            let param_eps = if k == 0 {
                10.0 // meters
            } else if k == 1 || k == 2 || k == 3 {
                1e-6 // radians
            } else if (k - 4) % 4 == 0 {
                1e-3 // dt (seconds)
            } else if (k - 4) % 4 == 1 {
                1e-2 // df0 (Hz)
            } else if (k - 4) % 4 == 2 {
                1e-4 // df1 (Hz/s)
            } else {
                1e-6 // df2 (Hz/s^2)
            };
            perturbed[k] += param_eps;

            let mut row_idx = 0;
            for (p_idx, pass) in raw_passes.iter().enumerate() {
                let dt = perturbed[4 + 4 * p_idx];
                let df0 = perturbed[4 + 4 * p_idx + 1];
                let df1 = perturbed[4 + 4 * p_idx + 2];
                let df2 = perturbed[4 + 4 * p_idx + 3];
                let t_ref = pass.points[0].time;
                for pt in &pass.points {
                    let pred = predict_frequency_poly(
                        perturbed[0],
                        perturbed[1],
                        perturbed[2],
                        perturbed[3],
                        epoch,
                        pt.time,
                        dt,
                        (df0, df1, df2),
                        t_ref,
                        center_freq,
                        rec_ecef,
                    );
                    let diff = pt.freq - (center_freq + pred);
                    jacobian[row_idx][k] = (diff - residuals[row_idx]) / param_eps;
                    row_idx += 1;
                }
            }
            for p_idx in 0..raw_passes.len() {
                let dt = perturbed[4 + 4 * p_idx];
                let diff = dt * 10.0;
                jacobian[row_idx][k] = (diff - residuals[row_idx]) / param_eps;
                row_idx += 1;
            }
            for p_idx in 0..raw_passes.len() {
                let df1 = perturbed[4 + 4 * p_idx + 2];
                let df2 = perturbed[4 + 4 * p_idx + 3];
                let diff_df1 = df1 * 100.0;
                let diff_df2 = df2 * 10000.0;
                jacobian[row_idx][k] = (diff_df1 - residuals[row_idx]) / param_eps;
                row_idx += 1;
                jacobian[row_idx][k] = (diff_df2 - residuals[row_idx]) / param_eps;
                row_idx += 1;
            }
        }

        // Form normal equations J^T J delta = J^T r
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

        // Apply Levenberg damping
        for k in 0..n_params {
            jt_j[k][k] += lambda * jt_j[k][k];
        }

        // Scale the system to prevent numerical underflow/overflow (Jacobi preconditioning)
        let mut scaled_jt_j = jt_j.clone();
        let mut scaled_jt_r = jt_r.clone();
        let mut scale_factors = vec![0.0; n_params];
        for i in 0..n_params {
            let s = jt_j[i][i].sqrt();
            scale_factors[i] = if s > 1e-15 { s } else { 1.0 };
        }
        for r in 0..n_params {
            scaled_jt_r[r] /= scale_factors[r];
            for c in 0..n_params {
                scaled_jt_j[r][c] /= scale_factors[r] * scale_factors[c];
            }
        }

        // Solve the system via Gaussian elimination
        if let Some(mut delta) = solve_linear_system(&mut scaled_jt_j, &scaled_jt_r) {
            for i in 0..n_params {
                delta[i] /= scale_factors[i];
            }
            for k in 0..n_params {
                params[k] -= delta[k];
            }

            // Constrain physical parameters to valid ranges
            params[0] = params[0].max(6500e3).min(20000e3); // Altitude range
            params[1] = params[1].max(0.0).min(std::f64::consts::PI); // Inclination

            // Scale-aware convergence: check relative step size per parameter.
            let max_rel_step = (0..n_params).map(|k| {
                let denom = params[k].abs().max(1e-10);
                delta[k].abs() / denom
            }).fold(0.0f64, f64::max);
            if max_rel_step < 1e-8 {
                break; // Converged
            }
        } else {
            break; // Matrix singular or solver failed
        }
    }

    params = best_params;
    if best_rss > 1e11 {
        return Err(format!("Solver failed to converge to a valid orbit (residual RSS too large: {:.2e})", best_rss).into());
    }
    let mut pass_dts = Vec::new();
    let mut pass_dfs = Vec::new();
    let mut pass_df1s = Vec::new();
    let mut pass_df2s = Vec::new();
    for j in 0..n_passes {
        pass_dts.push(params[4 + 4 * j]);
        pass_dfs.push(params[4 + 4 * j + 1]);
        pass_df1s.push(params[4 + 4 * j + 2]);
        pass_df2s.push(params[4 + 4 * j + 3]);
    }

    for &df in &pass_dfs {
        if df.abs() > 500e3 {
            return Err("Estimated frequency offset exceeds realistic physical limits".into());
        }
    }

    Ok(SolvedOrbit {
        a: params[0],
        i: params[1],
        raan0: params[2],
        u0: params[3],
        epoch,
        pass_dts,
        pass_dfs,
        pass_df1s,
        pass_df2s,
    })
}

// Simple Gaussian elimination linear solver
pub fn solve_linear_system(matrix: &mut [Vec<f64>], vector: &[f64]) -> Option<Vec<f64>> {
    let n = matrix.len();
    let mut aug = vec![vec![0.0; n + 1]; n];
    for r in 0..n {
        for c in 0..n {
            aug[r][c] = matrix[r][c];
        }
        aug[r][n] = vector[r];
    }

    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for r in i + 1..n {
            if aug[r][i].abs() > aug[max_row][i].abs() {
                max_row = r;
            }
        }
        aug.swap(i, max_row);

        if aug[i][i].abs() < 1e-12 {
            return None; // Singular matrix
        }

        // Eliminate columns
        for r in i + 1..n {
            let factor = aug[r][i] / aug[i][i];
            for c in i..n + 1 {
                aug[r][c] -= factor * aug[i][c];
            }
        }
    }

    // Back substitution
    let mut sol = vec![0.0; n];
    for r in (0..n).rev() {
        let mut sum = 0.0;
        for c in r + 1..n {
            sum += aug[r][c] * sol[c];
        }
        sol[r] = (aug[r][n] - sum) / aug[r][r];
    }

    Some(sol)
}

fn tle_checksum(line: &str) -> u32 {
    let mut sum = 0;
    for c in line.chars().take(68) {
        if c.is_ascii_digit() {
            sum += c.to_digit(10).unwrap();
        } else if c == '-' {
            sum += 1;
        }
    }
    sum % 10
}

#[allow(dead_code)]
pub fn format_tle_catalog(name: &str, orbit: &SolvedOrbit) -> String {
    let epoch = orbit.epoch;

    // Format epoch year and day fraction
    let year_2d = epoch.year() % 100;
    let day_of_year = epoch.ordinal();
    let day_fraction = (epoch.hour() as f64 * 3600.0
        + epoch.minute() as f64 * 60.0
        + epoch.second() as f64
        + epoch.timestamp_subsec_nanos() as f64 / 1e9)
        / 86400.0;
    let epoch_str = format!("{:02}{:03.8}", year_2d, day_of_year as f64 + day_fraction);

    // Format Keplerian elements
    let i_deg = orbit.i.to_degrees().rem_euclid(360.0);
    let raan_deg = orbit.raan0.to_degrees().rem_euclid(360.0);
    let u0_deg = orbit.u0.to_degrees().rem_euclid(360.0);

    // Mean motion: n = sqrt(MU/a^3) rad/s. Format to revs per day
    let n_rad_s = (MU / orbit.a.powi(3)).sqrt();
    let mean_motion = n_rad_s * 86400.0 / (2.0 * std::f64::consts::PI);

    // Construct line 1: catalog number 99999, class U, designator 26999A
    let mut l1 = format!(
        "1 99999U 26999A   {}  .00000000  00000-0  00000-0 0  999",
        epoch_str
    );
    // Pad to 68 characters if needed, or truncate
    if l1.len() < 68 {
        l1.push_str(&" ".repeat(68 - l1.len()));
    } else {
        l1.truncate(68);
    }
    let c1 = tle_checksum(&l1);
    let line1 = format!("{}{}", l1, c1); // append checksum at column 69

    // Construct line 2
    let mut l2 = format!(
        "2 99999 {:8.4} {:8.4} 0000000   0.0000 {:8.4} {:11.8}0000",
        i_deg, raan_deg, u0_deg, mean_motion
    );
    if l2.len() < 68 {
        l2.push_str(&" ".repeat(68 - l2.len()));
    } else {
        l2.truncate(68);
    }
    let c2 = tle_checksum(&l2);
    let line2 = format!("{}{}", l2, c2);

    format!("{}\n{}\n{}\n", name, line1, line2)
}

pub fn map_to_normalized_search_space(coord: [f64; 3]) -> Vec<f64> {
    // Maps continuous 3D coordinate parameters into a bounded [0, 1] search space
    // using base-p digit reversal mapping (Monna map) for multi-scale grid search.
    let mut result = vec![coord[0], coord[1], coord[2]];
    let primes = [2, 3, 5, 7];
    for &p in &primes {
        for j in 0..3 {
            let normalized = (coord[j].abs() / 20000000.0).min(0.99999);
            let val = inverse_monna_map(normalized, p, 16);
            let mapped = monna_map(val, p);
            result.push(mapped);
        }
    }
    result
}

pub fn compute_fractional_difference_history(x: &[f64]) -> Vec<f64> {
    // Calculates a fractional-like quotient over historical states
    // utilizing p-adic metric distances to estimate multi-scale variance.
    let n = x.len();
    let mut deriv = vec![0.0; n];
    if n <= 1 {
        return deriv;
    }
    let p = 2; // Prime for discretization
    let alpha = 0.5;
    let power = alpha + 1.0;
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            if i != j {
                let dist = p_adic_distance(i as u64, j as u64, p);
                if dist > 1e-9 {
                    sum += (x[i] - x[j]) / dist.powf(power);
                }
            }
        }
        deriv[i] = sum;
    }
    deriv
}

pub struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_f64(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.state >> 11) as f64 / (1u64 << 53) as f64
    }

    #[allow(dead_code)]
    pub fn next_gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-15);
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        r * theta.cos()
    }
}

pub fn monna_map(mut val: u64, p: u64) -> f64 {
    if p < 2 {
        return 0.0;
    }
    let mut res = 0.0;
    let mut base = 1.0 / (p as f64);
    while val > 0 {
        let digit = val % p;
        res += (digit as f64) * base;
        val /= p;
        base /= p as f64;
    }
    res
}

pub fn inverse_monna_map(mut x: f64, p: u64, precision: usize) -> u64 {
    if p < 2 {
        return 0;
    }
    x = x.clamp(0.0, 1.0 - 1e-15);
    let mut val = 0u64;
    let mut factor = 1u64;
    for _ in 0..precision {
        x *= p as f64;
        let digit = x.floor() as u64;
        val += (digit % p) * factor;
        x -= digit as f64;
        if factor.checked_mul(p).is_none() {
            break;
        }
        factor *= p;
    }
    val
}

#[allow(dead_code)]
pub fn p_adic_distance(a: u64, b: u64, p: u64) -> f64 {
    if a == b {
        return 0.0;
    }
    let diff = if a > b { a - b } else { b - a };
    let mut diff = diff;
    let mut vp = 0;
    while diff % p == 0 {
        vp += 1;
        diff /= p;
    }
    1.0 / (p as f64).powi(vp)
}

pub fn compute_rss(
    raw_passes: &[RawPass],
    rec_ecef: [f64; 3],
    initial_a: f64,
    initial_i: f64,
    epoch: DateTime<Utc>,
    center_freq: f64,
    raan0: f64,
    u0: f64,
) -> f64 {
    let pred_pcas =
        get_pred_pca_times(initial_a, initial_i, raan0, u0, epoch, rec_ecef, raw_passes);
    let mut rss = 0.0;
    for (p_idx, pass) in raw_passes.iter().enumerate() {
        if pass.points.is_empty() {
            continue;
        }
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
        rss += 100.0 * dt * dt;

        for pt in &pass.points {
            let pred_f = predict_frequency(
                initial_a,
                initial_i,
                raan0,
                u0,
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
    rss
}

pub fn compute_gradient(
    raw_passes: &[RawPass],
    rec_ecef: [f64; 3],
    a: f64,
    i: f64,
    epoch: DateTime<Utc>,
    center_freq: f64,
    raan0: f64,
    u0: f64,
) -> (f64, f64, f64, f64) {
    let rss_base = compute_rss(raw_passes, rec_ecef, a, i, epoch, center_freq, raan0, u0);

    let eps_a = 1000.0;
    let rss_a = compute_rss(
        raw_passes,
        rec_ecef,
        a + eps_a,
        i,
        epoch,
        center_freq,
        raan0,
        u0,
    );
    let grad_a = (rss_a - rss_base) / eps_a;

    let eps_i = 1e-4;
    let rss_i = compute_rss(
        raw_passes,
        rec_ecef,
        a,
        i + eps_i,
        epoch,
        center_freq,
        raan0,
        u0,
    );
    let grad_i = (rss_i - rss_base) / eps_i;

    let eps_raan = 1e-4;
    let rss_raan = compute_rss(
        raw_passes,
        rec_ecef,
        a,
        i,
        epoch,
        center_freq,
        raan0 + eps_raan,
        u0,
    );
    let grad_raan = (rss_raan - rss_base) / eps_raan;

    let eps_u0 = 1e-4;
    let rss_u0 = compute_rss(
        raw_passes,
        rec_ecef,
        a,
        i,
        epoch,
        center_freq,
        raan0,
        u0 + eps_u0,
    );
    let grad_u0 = (rss_u0 - rss_base) / eps_u0;

    (grad_a, grad_i, grad_raan, grad_u0)
}

use num_complex::Complex;
use serde::{Deserialize, Serialize};

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

// Deterministic LCG
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u32(&mut self) -> u32 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.state >> 32) as u32
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u32() as f64) / (u32::MAX as f64)
    }

    fn next_gaussian(&mut self) -> (f64, f64) {
        let u1 = self.next_f64().max(1e-15);
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        (r * theta.cos(), r * theta.sin())
    }
}

// Raised Cosine Filter
fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-9 {
        1.0
    } else {
        (x * std::f64::consts::PI).sin() / (x * std::f64::consts::PI)
    }
}

fn raised_cosine(t: f64, t_s: f64, beta: f64) -> f64 {
    let t_over_t_s = t / t_s;
    let numerator = sinc(t_over_t_s) * (std::f64::consts::PI * beta * t_over_t_s).cos();
    let denominator = 1.0 - (2.0 * beta * t_over_t_s).powi(2);
    if denominator.abs() < 1e-9 {
        sinc(t_over_t_s) * std::f64::consts::PI / 4.0
    } else {
        numerator / denominator
    }
}

use dsp::GardnerLoop;

trait GardnerTestExt {
    fn process_test(&mut self, sample: Complex<f32>, output_symbols: &mut Vec<(f64, Complex<f32>)>);
}

impl GardnerTestExt for GardnerLoop {
    fn process_test(
        &mut self,
        sample: Complex<f32>,
        output_symbols: &mut Vec<(f64, Complex<f32>)>,
    ) {
        let mut temp = Vec::new();
        self.process(sample, &mut temp);
        let sample_idx = self.sample_index;
        for (sym, mu) in temp {
            let t_des = mu as f64 + sample_idx - 2.0;
            output_symbols.push((t_des, sym));
        }
    }
}

// IQ Stream Generator
fn generate_iq_stream(
    symbol_rate: f64,
    sample_rate: f64,
    num_symbols: usize,
    timing_offset_s: f64,
    snr_db: f64,
    modulation: &str,
) -> (Vec<Complex<f32>>, Vec<Complex<f32>>) {
    let mut lcg = Lcg::new(42);
    let mut symbols = Vec::with_capacity(num_symbols);
    for _ in 0..num_symbols {
        if modulation == "bpsk" {
            let val = if lcg.next_f64() < 0.5 { 1.0 } else { -1.0 };
            symbols.push(Complex::new(val, 0.0));
        } else {
            let re = if lcg.next_f64() < 0.5 { 1.0 } else { -1.0 };
            let im = if lcg.next_f64() < 0.5 { 1.0 } else { -1.0 };
            symbols.push(Complex::new(re, im) * (1.0 / 2.0f32.sqrt()));
        }
    }

    let t_s = 1.0 / symbol_rate;
    let beta = 0.35;
    let total_samples = ((num_symbols as f64) * (sample_rate / symbol_rate)) as usize;
    let mut iq_samples = Vec::with_capacity(total_samples);

    for n in 0..total_samples {
        let t_n = (n as f64) / sample_rate;
        let mut sample_val = Complex::new(0.0f32, 0.0f32);

        let k_center = (t_n / t_s).round() as isize;
        let k_min = (k_center - 8).max(0) as usize;
        let k_max = (k_center + 8).min(num_symbols as isize - 1) as usize;

        for k in k_min..=k_max {
            let symbol_t = 1.0 / sample_rate + (k as f64) * t_s; // Optimal alignment shift
            let t_diff = t_n - symbol_t - timing_offset_s;
            let h = raised_cosine(t_diff, t_s, beta);
            sample_val += Complex::new(symbols[k].re * h as f32, symbols[k].im * h as f32);
        }

        if snr_db.is_finite() {
            let snr_linear = 10.0f64.powf(snr_db / 10.0);
            let std_dev = (1.0 / (2.0 * snr_linear)).sqrt() as f32;
            let (n1, n2) = lcg.next_gaussian();
            sample_val += Complex::new((n1 as f32) * std_dev, (n2 as f32) * std_dev);
        }

        iq_samples.push(sample_val);
    }

    (iq_samples, symbols)
}

#[test]
fn test_empirical_gardner_properties() {
    let run_test = |symbol_rate: f64,
                    sample_rate: f64,
                    num_symbols: usize,
                    offset_s: f64,
                    snr_db: f64|
     -> Vec<f64> {
        let sps = sample_rate / symbol_rate;
        let offset_samples = offset_s * sample_rate;

        let (iq, _) = generate_iq_stream(
            symbol_rate,
            sample_rate,
            num_symbols,
            offset_s,
            snr_db,
            "bpsk",
        );

        let mut gardner = GardnerLoop::new(sample_rate, symbol_rate);
        let mut recovered = Vec::new();
        for &s in &iq {
            gardner.process_test(s, &mut recovered);
        }

        // Calculate timing error for each output symbol
        recovered
            .iter()
            .enumerate()
            .map(|(i, &(t_des, _))| {
                let t_opt = 2.0 + (i as f64) * sps + offset_samples;
                let mut diff = t_des - t_opt;
                // Unwrap modulo sps
                diff = (diff + sps / 2.0).rem_euclid(sps) - sps / 2.0;
                diff / sps // in symbol periods
            })
            .collect()
    };

    // 1. Standard Convergence Test
    println!("--- Standard Case (10 kHz, sps=5.0) ---");
    let errors_std = run_test(10000.0, 50000.0, 600, 0.25 / 10000.0, f64::INFINITY);
    let initial_rmse_std = (errors_std[10..30].iter().map(|e| e * e).sum::<f64>() / 20.0).sqrt();
    let final_rmse_std = (errors_std[errors_std.len() - 50..]
        .iter()
        .map(|e| e * e)
        .sum::<f64>()
        / 50.0)
        .sqrt();
    println!("Initial RMSE: {:.4} symbol periods", initial_rmse_std);
    println!("Final RMSE: {:.4} symbol periods", final_rmse_std);
    assert!(
        final_rmse_std < 0.055,
        "Standard loop should converge to < 0.055 error"
    );

    // 2. Extremely Low Symbol Rate slow-convergence Test
    println!("--- Low Rate Case (100 Hz, sps=500.0) ---");
    let errors_low = run_test(100.0, 50000.0, 400, 0.25 / 100.0, f64::INFINITY);
    let initial_rmse_low = (errors_low[10..30].iter().map(|e| e * e).sum::<f64>() / 20.0).sqrt();
    let final_rmse_low = (errors_low[errors_low.len() - 50..]
        .iter()
        .map(|e| e * e)
        .sum::<f64>()
        / 50.0)
        .sqrt();
    println!("Initial RMSE: {:.4} symbol periods", initial_rmse_low);
    println!("Final RMSE: {:.4} symbol periods", final_rmse_low);
    // Since it converges extremely slowly, the final RMSE should still be very large (similar to initial)
    assert!(
        final_rmse_low > 0.15,
        "Low rate loop should converge very slowly"
    );

    // 3. Noise Jitter and Loop gain: High rate (sps=2.0) vs Standard (sps=5.0) under SNR = 15 dB
    // Under the same SNR, the loop with smaller sps (2.0) has larger relative gain, so it has higher steady-state noise-driven jitter.
    println!("--- Jitter Test (SNR = 15 dB) ---");
    let errors_std_noisy = run_test(10000.0, 50000.0, 500, 0.0, 15.0);
    let errors_high_noisy = run_test(1000000.0, 2000000.0, 500, 0.0, 15.0);

    let jitter_std = (errors_std_noisy[errors_std_noisy.len() - 100..]
        .iter()
        .map(|e| e * e)
        .sum::<f64>()
        / 100.0)
        .sqrt();
    let jitter_high = (errors_high_noisy[errors_high_noisy.len() - 100..]
        .iter()
        .map(|e| e * e)
        .sum::<f64>()
        / 100.0)
        .sqrt();
    println!(
        "Steady-state jitter Standard (sps=5.0): {:.4} symbol periods",
        jitter_std
    );
    println!(
        "Steady-state jitter High Rate (sps=2.0): {:.4} symbol periods",
        jitter_high
    );
    // High rate loop jitter is larger due to higher relative loop gain
    assert!(
        jitter_high > jitter_std,
        "High rate loop should have higher noise-driven timing jitter"
    );
}

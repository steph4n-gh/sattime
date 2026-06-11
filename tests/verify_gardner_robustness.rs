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

// Linear Congruential Generator for deterministic random number generation
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

// Raised Cosine Filter functions
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
    fn process_test(&mut self, sample: Complex<f32>, output_symbols: &mut Vec<Complex<f32>>);
}

impl GardnerTestExt for GardnerLoop {
    fn process_test(&mut self, sample: Complex<f32>, output_symbols: &mut Vec<Complex<f32>>) {
        let mut temp = Vec::new();
        self.process(sample, &mut temp);
        for (sym, _) in temp {
            output_symbols.push(sym);
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
    fading_func: Option<fn(f64) -> f64>,
    freq_offset_hz: f64,
    phase_offset_rad: f64,
    modulation: &str,
) -> (Vec<Complex<f32>>, Vec<Complex<f32>>) {
    let mut lcg = Lcg::new(123456789);

    // 1. Generate random symbols
    let mut symbols = Vec::with_capacity(num_symbols);
    for _ in 0..num_symbols {
        if modulation == "bpsk" {
            let val = if lcg.next_f64() < 0.5 { 1.0 } else { -1.0 };
            symbols.push(Complex::new(val, 0.0));
        } else {
            // QPSK
            let re = if lcg.next_f64() < 0.5 { 1.0 } else { -1.0 };
            let im = if lcg.next_f64() < 0.5 { 1.0 } else { -1.0 };
            symbols.push(Complex::new(re, im) * (1.0 / 2.0f32.sqrt()));
        }
    }

    // 2. Generate samples at sample_rate
    let t_s = 1.0 / symbol_rate;
    let beta = 0.35; // raised cosine roll-off factor
    let total_samples = ((num_symbols as f64) * (sample_rate / symbol_rate)) as usize;
    let mut iq_samples = Vec::with_capacity(total_samples);

    for n in 0..total_samples {
        let t_n = (n as f64) / sample_rate;
        let mut sample_val = Complex::new(0.0f32, 0.0f32);

        // Sum over nearby symbols (e.g., within 8 symbol periods)
        let k_center = (t_n / t_s).round() as isize;
        let k_min = (k_center - 8).max(0) as usize;
        let k_max = (k_center + 8).min(num_symbols as isize - 1) as usize;

        for k in k_min..=k_max {
            let symbol_t = 1.0 / sample_rate + (k as f64) * t_s; // Optimal alignment shift
            let t_diff = t_n - symbol_t - timing_offset_s;
            let h = raised_cosine(t_diff, t_s, beta);
            sample_val += Complex::new(symbols[k].re * h as f32, symbols[k].im * h as f32);
        }

        // Apply fading
        let fade = if let Some(f) = fading_func {
            f(t_n)
        } else {
            1.0
        };
        sample_val = Complex::new(sample_val.re * fade as f32, sample_val.im * fade as f32);

        // Apply Doppler frequency and phase offset
        let phase = 2.0 * std::f64::consts::PI * freq_offset_hz * t_n + phase_offset_rad;
        let rot = Complex::new(phase.cos() as f32, phase.sin() as f32);
        sample_val = sample_val * rot;

        // Add complex Gaussian noise if SNR is finite
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

// Find best delay and match percentage (supporting negative delays)
fn find_best_delay(true_symbols: &[Complex<f32>], out_symbols: &[Complex<f32>]) -> (isize, f64) {
    let mut best_delay = 0;
    let mut best_pct = 0.0;
    for delay in -10..30 {
        let mut correct = 0;
        let mut total = 0;
        for i in 0..true_symbols.len() {
            let out_idx = (i as isize + delay) as usize;
            if out_idx < out_symbols.len() {
                let out_sym = out_symbols[out_idx];
                let true_sym = true_symbols[i];
                if out_sym.re.signum() == true_sym.re.signum() {
                    correct += 1;
                }
                total += 1;
            }
        }
        if total > 0 {
            let pct = (correct as f64) / (total as f64);
            if pct > best_pct {
                best_pct = pct;
                best_delay = delay;
            }
        }
    }
    (best_delay, best_pct)
}

#[test]
fn test_gardner_zero_timing_error() {
    let symbol_rate = 10000.0;
    let sample_rate = 50000.0;
    let num_symbols = 500;

    let (iq_samples, true_symbols) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        0.0,
        f64::INFINITY,
        None,
        0.0,
        0.0,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();

    for &s in &iq_samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Since zero timing error, integrator drift should be small
    assert!(
        loop_state.integrator.abs() < 5e-3,
        "Integrator drifted: {}",
        loop_state.integrator
    );
    assert!(out_symbols.len() >= num_symbols - 10);

    let (_, pct) = find_best_delay(&true_symbols, &out_symbols);
    assert!(
        pct > 0.95,
        "Symbol recovery accuracy with zero timing error was too low: {:.1}%",
        pct * 100.0
    );
}

#[test]
fn test_gardner_timing_error_tracking_positive_offset() {
    let symbol_rate = 10000.0;
    let sample_rate = 50000.0;
    let num_symbols = 800;
    let t_s = 1.0 / symbol_rate;
    let offset = 0.2 * t_s; // late signal (early receiver sampling relative to peak)

    let (iq_samples, true_symbols) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        offset,
        f64::INFINITY,
        None,
        0.0,
        0.0,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();
    let mut integrator_history = Vec::new();

    for &s in &iq_samples {
        loop_state.process_test(s, &mut out_symbols);
        integrator_history.push(loop_state.integrator);
    }

    // Average integrator value during pull-in (samples 20 to 100)
    // For a positive timing offset (late signal / early sampling), the integrator pull-in should be negative (step size increases)
    let avg_integrator_pullin: f64 = integrator_history[20..100].iter().sum::<f64>() / 80.0;
    assert!(
        avg_integrator_pullin < 0.0,
        "Integrator was not negative during pull-in: {}",
        avg_integrator_pullin
    );

    // Verify convergence and recovery
    let (_, pct) = find_best_delay(
        &true_symbols[num_symbols / 2..],
        &out_symbols[out_symbols.len() / 2..],
    );
    assert!(
        pct > 0.90,
        "Symbol recovery accuracy after positive timing convergence was too low: {:.1}%",
        pct * 100.0
    );
}

#[test]
fn test_gardner_timing_error_tracking_negative_offset() {
    let symbol_rate = 10000.0;
    let sample_rate = 50000.0;
    let num_symbols = 800;
    let t_s = 1.0 / symbol_rate;
    let offset = -0.2 * t_s; // early signal (late receiver sampling relative to peak)

    let (iq_samples, true_symbols) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        offset,
        f64::INFINITY,
        None,
        0.0,
        0.0,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();
    let mut integrator_history = Vec::new();

    for &s in &iq_samples {
        loop_state.process_test(s, &mut out_symbols);
        integrator_history.push(loop_state.integrator);
    }

    // For a negative timing offset (early signal / late sampling), the integrator pull-in should be positive (step size decreases)
    let avg_integrator_pullin: f64 = integrator_history[20..100].iter().sum::<f64>() / 80.0;
    assert!(
        avg_integrator_pullin > 0.0,
        "Integrator was not positive during pull-in: {}",
        avg_integrator_pullin
    );

    // Verify convergence and recovery
    let (_, pct) = find_best_delay(
        &true_symbols[num_symbols / 2..],
        &out_symbols[out_symbols.len() / 2..],
    );
    assert!(
        pct > 0.90,
        "Symbol recovery accuracy after negative timing convergence was too low: {:.1}%",
        pct * 100.0
    );
}

#[test]
fn test_gardner_symbol_rate_extremely_low() {
    let symbol_rate = 100.0;
    let sample_rate = 50000.0;
    let num_symbols = 50;

    let (iq_samples, _) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        0.1 / symbol_rate,
        f64::INFINITY,
        None,
        0.0,
        0.0,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();

    for &s in &iq_samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    assert!(loop_state.integrator.is_finite());
    assert!(loop_state.step > 0.0);
    assert!(out_symbols.len() > 0);
}

#[test]
fn test_gardner_symbol_rate_extremely_high() {
    let symbol_rate = 20000.0;
    let sample_rate = 50000.0;
    let num_symbols = 500;

    let (iq_samples, _) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        0.1 / symbol_rate,
        f64::INFINITY,
        None,
        0.0,
        0.0,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();

    for &s in &iq_samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    assert!(loop_state.integrator.is_finite());
    assert!(loop_state.step > 0.0);
    assert!(out_symbols.len() >= num_symbols - 10);
}

#[test]
fn test_gardner_loop_stability_under_varying_snr() {
    let symbol_rate = 10000.0;
    let sample_rate = 50000.0;
    let num_symbols = 500;
    let snr_values = vec![25.0, 15.0, 8.0, 3.0, 0.0];

    for &snr in &snr_values {
        let (iq_samples, true_symbols) = generate_iq_stream(
            symbol_rate,
            sample_rate,
            num_symbols,
            0.1 / symbol_rate,
            snr,
            None,
            0.0,
            0.0,
            "bpsk",
        );

        let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
        let mut out_symbols = Vec::new();

        for &s in &iq_samples {
            loop_state.process_test(s, &mut out_symbols);
        }

        assert!(
            loop_state.integrator.is_finite(),
            "Integrator exploded at SNR = {} dB",
            snr
        );
        assert!(
            loop_state.step.is_finite(),
            "Step size became NaN/Inf at SNR = {} dB",
            snr
        );

        if snr >= 8.0 {
            let (_, pct) = find_best_delay(
                &true_symbols[num_symbols / 2..],
                &out_symbols[out_symbols.len() / 2..],
            );
            assert!(
                pct > 0.80,
                "Symbol recovery at SNR = {} dB was too low: {:.1}%",
                snr,
                pct * 100.0
            );
        }
    }
}

#[test]
fn test_gardner_loop_stability_under_doppler_rate() {
    let symbol_rate = 10000.0;
    let sample_rate = 50000.0;
    let num_symbols = 600;

    let (iq_samples, _) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        0.15 / symbol_rate,
        f64::INFINITY,
        None,
        500.0, // 500 Hz Doppler offset
        0.5,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();

    for &s in &iq_samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    assert!(loop_state.integrator.is_finite());
    assert!(loop_state.step.is_finite());

    // Verify magnitude clusters around the unit circle
    let delay = 15;
    let mut high_magnitude_count = 0;
    let start_idx = num_symbols / 2;
    let end_idx = num_symbols - delay;
    for i in start_idx..end_idx {
        if i < out_symbols.len() {
            let magnitude = out_symbols[i].norm();
            if magnitude > 0.6 && magnitude < 1.4 {
                high_magnitude_count += 1;
            }
        }
    }
    let pct = (high_magnitude_count as f64) / ((end_idx - start_idx) as f64);
    assert!(
        pct > 0.80,
        "Recovered symbols magnitude under Doppler was poor: {:.1}%",
        pct * 100.0
    );
}

#[test]
fn test_gardner_loop_stability_under_heavy_fading() {
    let symbol_rate = 10000.0;
    let sample_rate = 50000.0;
    let num_symbols = 1500;

    // Deep fade from t=0.04s to 0.10s (symbols 400 to 1000)
    fn fading_profile(t: f64) -> f64 {
        if t >= 0.04 && t < 0.10 { 0.0 } else { 1.0 }
    }

    let (iq_samples, true_symbols) = generate_iq_stream(
        symbol_rate,
        sample_rate,
        num_symbols,
        0.15 / symbol_rate,
        f64::INFINITY,
        Some(fading_profile),
        0.0,
        0.0,
        "bpsk",
    );

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut out_symbols = Vec::new();

    // Process normal start
    let normal1_samples_len = (400.0 * (sample_rate / symbol_rate)) as usize;
    for &s in &iq_samples[..normal1_samples_len] {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Process transition of fade
    let fade_start = normal1_samples_len;
    let fade_transition_len = 30;
    for &s in &iq_samples[fade_start..fade_start + fade_transition_len] {
        loop_state.process_test(s, &mut out_symbols);
    }

    let integrator_after_transition = loop_state.integrator;
    let step_after_transition = loop_state.step;

    // Process remainder of fade
    let fade_end = (1000.0 * (sample_rate / symbol_rate)) as usize;
    for &s in &iq_samples[fade_start + fade_transition_len..fade_end] {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Check that during the fade, the loop remains frozen/stable
    assert_eq!(
        loop_state.integrator, integrator_after_transition,
        "Integrator changed during fade!"
    );
    assert_eq!(
        loop_state.step, step_after_transition,
        "Step size changed during fade!"
    );

    // Process return of signal
    for &s in &iq_samples[fade_end..] {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Check recovery after fade (symbols 1100 to 1450)
    let (_, pct) = find_best_delay(&true_symbols[1100..1450], &out_symbols[1100..]);
    assert!(
        pct > 0.85,
        "Symbol recovery accuracy after fading was too low: {:.1}%",
        pct * 100.0
    );
}

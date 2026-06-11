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

// --- Helper Functions for Generating Signals ---

fn raised_cosine_pulse(t: f64, sps: f64, beta: f64) -> f64 {
    let t_norm = t / sps;
    if t_norm.abs() < 1e-9 {
        return 1.0;
    }
    if beta > 0.0 && (2.0 * beta * t_norm.abs() - 1.0).abs() < 1e-9 {
        return (std::f64::consts::PI / 4.0) * (std::f64::consts::FRAC_2_PI * t_norm.sin());
    }
    let sinc = (t_norm * std::f64::consts::PI).sin() / (t_norm * std::f64::consts::PI);
    let numerator = (beta * t_norm * std::f64::consts::PI).cos();
    let denominator = 1.0 - (2.0 * beta * t_norm).powi(2);
    sinc * numerator / denominator
}

fn generate_bpsk_signal(
    symbols: &[f32],
    sps: f64,
    num_samples: usize,
    timing_offset_samples: f64,
    frequency_offset_hz: f64,
    sample_rate: f64,
    snr_db: Option<f64>,
) -> Vec<Complex<f32>> {
    let mut samples = Vec::with_capacity(num_samples);
    let beta = 0.35;

    let mut rng_state = 123456789u64;
    let mut next_noise = move || {
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let val = (rng_state >> 32) as f32 / f32::MAX;
        val * 2.0 - 1.0
    };

    let noise_std = snr_db.map(|db| {
        let snr_linear = 10.0f64.powf(db / 10.0);
        (1.0 / (2.0 * snr_linear)).sqrt() as f32
    });

    for n in 0..num_samples {
        let t = n as f64 - timing_offset_samples;
        let mut sample_val = 0.0;
        for (k, &sym) in symbols.iter().enumerate() {
            let symbol_time = k as f64 * sps;
            if (t - symbol_time).abs() < 10.0 * sps {
                sample_val += sym as f64 * raised_cosine_pulse(t - symbol_time, sps, beta);
            }
        }

        let phase = 2.0 * std::f64::consts::PI * frequency_offset_hz * (n as f64) / sample_rate;
        let mut sample = Complex::new(sample_val as f32, 0.0)
            * Complex::new(phase.cos() as f32, phase.sin() as f32);

        if let Some(std) = noise_std {
            let n_re = next_noise() * std;
            let n_im = next_noise() * std;
            sample += Complex::new(n_re, n_im);
        }

        samples.push(sample);
    }
    samples
}

// --- Adversarial & Corner Case Tests ---

#[test]
fn test_gardner_extremely_low_symbol_rate() {
    let sample_rate = 50000.0;
    let symbol_rate = 10.0; // Extremely low rate: sps = 5000.0
    let sps = sample_rate / symbol_rate;

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    assert_eq!(loop_state.sps, 5000.0);
    assert_eq!(loop_state.step, 2500.0);

    // Generate 5 symbols
    let symbols = vec![1.0, -1.0, 1.0, -1.0, 1.0];
    let num_samples = (symbols.len() as f64 * sps + 100.0) as usize;
    let samples = generate_bpsk_signal(&symbols, sps, num_samples, 0.0, 0.0, sample_rate, None);

    let mut out_symbols = Vec::new();
    for s in samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Since we generated 5 symbols, the loop should recover around 5 symbols
    assert!(
        out_symbols.len() >= 4 && out_symbols.len() <= 6,
        "Expected around 5 symbols, got {}",
        out_symbols.len()
    );
    assert!(
        loop_state.integrator.abs() < 1e-3,
        "Integrator drifted significantly: {}",
        loop_state.integrator
    );
}

#[test]
fn test_gardner_extremely_high_symbol_rate() {
    let sample_rate = 50000.0;
    let symbol_rate = 24000.0; // Near Nyquist: sps = 2.0833
    let sps = sample_rate / symbol_rate;

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);

    // Generate 100 symbols
    let mut symbols = Vec::new();
    for i in 0..100 {
        symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
    }
    let num_samples = (symbols.len() as f64 * sps + 20.0) as usize;
    let samples = generate_bpsk_signal(&symbols, sps, num_samples, 0.0, 0.0, sample_rate, None);

    let mut out_symbols = Vec::new();
    for s in samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    assert!(
        out_symbols.len() > 90,
        "Expected close to 100 symbols, got {}",
        out_symbols.len()
    );
    assert!(loop_state.integrator.abs() < 1e-2);
}

#[test]
fn test_gardner_super_nyquist_handling() {
    let sample_rate = 50000.0;
    let symbol_rate = 100000.0; // Above Nyquist! sps = 0.5, nominal_step = 0.25
    let sps = sample_rate / symbol_rate;

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    assert_eq!(loop_state.sps, 0.5);
    assert_eq!(loop_state.step, 0.25);

    let symbols = vec![1.0, -1.0, 1.0, -1.0, 1.0];
    let num_samples = 100;
    let samples = generate_bpsk_signal(&symbols, sps, num_samples, 0.0, 0.0, sample_rate, None);

    let mut out_symbols = Vec::new();
    // Verify that processing many samples does not hang or overflow due to step size
    for s in samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    // It should terminate successfully without hanging or panic
    assert!(out_symbols.len() > 0);
}

#[test]
fn test_gardner_zero_timing_error_robustness() {
    let sample_rate = 50000.0;
    let symbol_rate = 10000.0; // sps = 5.0
    let sps = sample_rate / symbol_rate;

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);

    // BPSK symbols with no timing offset
    let mut symbols = Vec::new();
    for i in 0..100 {
        symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
    }

    // Sample exactly on-time (aligned offset = 1.0)
    let num_samples = (symbols.len() as f64 * sps + 10.0) as usize;
    let samples = generate_bpsk_signal(&symbols, sps, num_samples, 1.0, 0.0, sample_rate, None);

    let mut out_symbols = Vec::new();
    for s in samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Under zero timing error, the integrator should remain very close to 0 (within small startup/ISI transient)
    assert!(
        loop_state.integrator.abs() < 1e-2,
        "Integrator drifted on zero timing error: {}",
        loop_state.integrator
    );
}

#[test]
fn test_gardner_timing_error_polarity_sign() {
    let sample_rate = 50000.0;
    let symbol_rate = 10000.0; // sps = 5.0
    let sps = sample_rate / symbol_rate;

    // Trace trajectory for late sampling (0.75) and verify correct sign before overshoot
    {
        let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
        let mut symbols = Vec::new();
        for i in 0..100 {
            symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
        }
        let num_samples = (symbols.len() as f64 * sps) as usize;
        let samples =
            generate_bpsk_signal(&symbols, sps, num_samples, 0.75, 0.0, sample_rate, None);

        let mut out_symbols = Vec::new();
        let mut integrator_at_100 = 0.0;
        for (idx, s) in samples.into_iter().enumerate() {
            loop_state.process_test(s, &mut out_symbols);
            if idx == 100 {
                integrator_at_100 = loop_state.integrator;
            }
        }
        // With late sampling, the loop needs to step shorter, meaning control > 0, so integrator should go positive
        assert!(
            integrator_at_100 > 0.0,
            "Expected positive timing correction at sample 100, got {}",
            integrator_at_100
        );
    }

    // Trace trajectory for early sampling (1.25) and verify correct sign before overshoot
    {
        let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
        let mut symbols = Vec::new();
        for i in 0..100 {
            symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
        }
        let num_samples = (symbols.len() as f64 * sps) as usize;
        let samples =
            generate_bpsk_signal(&symbols, sps, num_samples, 1.25, 0.0, sample_rate, None);

        let mut out_symbols = Vec::new();
        let mut integrator_at_100 = 0.0;
        for (idx, s) in samples.into_iter().enumerate() {
            loop_state.process_test(s, &mut out_symbols);
            if idx == 100 {
                integrator_at_100 = loop_state.integrator;
            }
        }
        // With early sampling, the loop needs to step longer, meaning control < 0, so integrator should go negative
        assert!(
            integrator_at_100 < 0.0,
            "Expected negative timing correction at sample 100, got {}",
            integrator_at_100
        );
    }
}

#[test]
fn test_gardner_heavy_fading_ride_through() {
    let sample_rate = 50000.0;
    let symbol_rate = 10000.0; // sps = 5.0
    let sps = sample_rate / symbol_rate;

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut symbols = Vec::new();
    for i in 0..200 {
        symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
    }

    // First 50 symbols: normal signal (late sampling with offset 0.75)
    let normal_samples = generate_bpsk_signal(
        &symbols[..50],
        sps,
        (50.0 * sps) as usize,
        0.75,
        0.0,
        sample_rate,
        None,
    );
    // Next 100 symbols: faded (zero amplitude)
    let faded_samples = vec![Complex::new(0.0f32, 0.0f32); (100.0 * sps) as usize];
    // Final 50 symbols: normal signal returns
    let return_samples = generate_bpsk_signal(
        &symbols[50..100],
        sps,
        (50.0 * sps) as usize,
        0.75,
        0.0,
        sample_rate,
        None,
    );

    let mut out_symbols = Vec::new();
    for s in normal_samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Record loop state before fade
    let integrator_pre_fade = loop_state.integrator;
    assert!(integrator_pre_fade > 0.0);

    // Process fade transition (first 50 samples of fade to fully empty history and clear out transient)
    for &s in &faded_samples[..50] {
        loop_state.process_test(s, &mut out_symbols);
    }
    let integrator_after_transition = loop_state.integrator;
    let step_after_transition = loop_state.step;

    // Process remaining faded samples
    for &s in &faded_samples[50..] {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Verify that the loop state is stable and unchanged during fade after transition
    assert_eq!(
        loop_state.integrator, integrator_after_transition,
        "Integrator changed during fade!"
    );
    assert_eq!(
        loop_state.step, step_after_transition,
        "Step size changed during fade!"
    );

    // Process return of signal
    for s in return_samples {
        loop_state.process_test(s, &mut out_symbols);
    }

    // Loop should continue to track and adapt
    assert!(loop_state.integrator.is_normal());
}

#[test]
fn test_gardner_varying_snr_resilience() {
    let sample_rate = 50000.0;
    let symbol_rate = 10000.0; // sps = 5.0
    let sps = sample_rate / symbol_rate;

    let snr_levels = vec![30.0, 15.0, 5.0, 0.0]; // SNR DB levels
    let mut symbols = Vec::new();
    for i in 0..100 {
        symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
    }

    for snr in snr_levels {
        let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
        let num_samples = (symbols.len() as f64 * sps) as usize;
        let samples =
            generate_bpsk_signal(&symbols, sps, num_samples, 0.2, 0.0, sample_rate, Some(snr));

        let mut out_symbols = Vec::new();
        for s in samples {
            loop_state.process_test(s, &mut out_symbols);
        }

        // Even under noisy conditions, loop should adapt and remain stable (not explode)
        assert!(loop_state.integrator.is_finite());
        assert!(loop_state.step >= 0.5 * (sps / 2.0) && loop_state.step <= 1.5 * (sps / 2.0));
    }
}

#[test]
fn test_gardner_phase_and_frequency_invariance() {
    let sample_rate = 50000.0;
    let symbol_rate = 10000.0; // sps = 5.0
    let sps = sample_rate / symbol_rate;

    let mut loop_state = GardnerLoop::new(sample_rate, symbol_rate);
    let mut symbols = Vec::new();
    for i in 0..100 {
        symbols.push(if i % 2 == 0 { 1.0f32 } else { -1.0f32 });
    }

    // Generate signal with large carrier frequency offset (e.g. 1000 Hz)
    let num_samples = (symbols.len() as f64 * sps) as usize;
    let samples = generate_bpsk_signal(&symbols, sps, num_samples, 0.75, 1000.0, sample_rate, None);

    let mut out_symbols = Vec::new();
    let mut integrator_at_100 = 0.0;
    for (idx, s) in samples.into_iter().enumerate() {
        loop_state.process_test(s, &mut out_symbols);
        if idx == 100 {
            integrator_at_100 = loop_state.integrator;
        }
    }

    // Gardner loop is non-coherent (phase-invariance), so it should track timing correctly even with a large frequency offset
    assert!(
        integrator_at_100 > 0.0,
        "Expected positive timing correction under frequency offset at sample 100, got {}",
        integrator_at_100
    );
}

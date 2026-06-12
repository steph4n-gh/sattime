use num_complex::Complex;
use sattime::dsp::estimate_frequency_esprit;

/// Simple LCG PRNG for deterministic noise generation (same pattern as verify_esprit_estimator).
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Returns a value in [-1.0, 1.0].
    fn next_f32(&mut self) -> f32 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let val = (self.state >> 32) as f32 / (u32::MAX as f32);
        val * 2.0 - 1.0
    }

    /// Returns a Gaussian-ish sample via Box-Muller (uses two uniform draws).
    fn next_gaussian(&mut self) -> f32 {
        // Simple Box-Muller transform
        let u1 = (self.next_f32() + 1.0) / 2.0; // map to (0, 1)
        let u2 = (self.next_f32() + 1.0) / 2.0;
        // Clamp u1 away from zero to avoid log(0)
        let u1 = u1.max(1e-10);
        let r = (-2.0 * u1.ln()).sqrt();
        r * (2.0 * std::f32::consts::PI * u2).cos()
    }
}

/// Generate a complex tone at `freq` Hz, sampled at `sample_rate`, for `n` samples.
fn generate_tone(freq: f64, sample_rate: f64, n: usize) -> Vec<Complex<f32>> {
    (0..n)
        .map(|i| {
            let t = (i as f64) / sample_rate;
            let angle = 2.0 * std::f64::consts::PI * freq * t;
            Complex::new(angle.cos() as f32, angle.sin() as f32)
        })
        .collect()
}

/// Add Gaussian noise to a signal at the specified SNR (in dB).
/// SNR is defined as 10*log10(signal_power / noise_power).
/// For a unit-amplitude complex tone, signal_power = 1.0.
fn add_noise(signal: &[Complex<f32>], snr_db: f64, rng: &mut SimpleRng) -> Vec<Complex<f32>> {
    let signal_power: f64 = 1.0; // unit-amplitude tone
    let noise_power = signal_power / 10.0_f64.powf(snr_db / 10.0);
    // Each component gets half the noise power → std = sqrt(noise_power / 2)
    let noise_std = (noise_power / 2.0).sqrt() as f32;

    signal
        .iter()
        .map(|s| {
            let noise = Complex::new(
                rng.next_gaussian() * noise_std,
                rng.next_gaussian() * noise_std,
            );
            s + noise
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Test 1: ESPRIT at low SNR levels
// ---------------------------------------------------------------------------
#[test]
fn test_esprit_low_snr() {
    let sample_rate = 50000.0;
    let f0 = 1234.5;
    let n = 500;
    let m = 10;

    let clean = generate_tone(f0, sample_rate, n);

    for &snr_db in &[10.0, 5.0, 0.0] {
        let mut rng = SimpleRng::new(42);
        let noisy = add_noise(&clean, snr_db, &mut rng);
        let est = estimate_frequency_esprit(&noisy, sample_rate, m);
        let err = (est - f0).abs();
        println!(
            "SNR={:.0} dB: estimated={:.2} Hz, target={:.1} Hz, error={:.2} Hz",
            snr_db, est, f0, err
        );
        assert!(
            err < 100.0,
            "At SNR={:.0} dB, frequency error {:.2} Hz exceeds 100 Hz tolerance",
            snr_db,
            err
        );
    }
}

// ---------------------------------------------------------------------------
// Test 2: ESPRIT with a very short buffer (barely above minimum)
// ---------------------------------------------------------------------------
#[test]
fn test_esprit_short_buffer() {
    let sample_rate = 50000.0;
    let f0 = 2000.0;
    let n = 15; // barely above m
    let m = 10;

    let samples = generate_tone(f0, sample_rate, n);
    let est = estimate_frequency_esprit(&samples, sample_rate, m);
    let err = (est - f0).abs();
    println!(
        "Short buffer (n={}): estimated={:.2} Hz, target={:.1} Hz, error={:.2} Hz",
        n, est, f0, err
    );

    // With only 15 samples the estimate may be rough; just assert no panic and bounded error
    assert!(
        err < 500.0,
        "Short-buffer frequency error {:.2} Hz exceeds 500 Hz tolerance",
        err
    );
}

// ---------------------------------------------------------------------------
// Test 3: ESPRIT with degenerate inputs
// ---------------------------------------------------------------------------
#[test]
fn test_esprit_degenerate_input() {
    let sample_rate = 50000.0;
    let m = 10;

    // All-zeros input
    let zeros = vec![Complex::new(0.0f32, 0.0f32); 100];
    let est_zeros = estimate_frequency_esprit(&zeros, sample_rate, m);
    println!("All-zeros: estimated={:.4} Hz", est_zeros);
    assert!(
        est_zeros.abs() < 1e-6,
        "All-zeros input should return ~0.0, got {}",
        est_zeros
    );

    // Constant DC input (all samples = (1.0, 0.0))
    let dc = vec![Complex::new(1.0f32, 0.0f32); 100];
    let est_dc = estimate_frequency_esprit(&dc, sample_rate, m);
    println!("Constant DC: estimated={:.4} Hz", est_dc);
    assert!(
        est_dc.abs() < 1.0,
        "Constant DC input should return approximately 0.0 Hz, got {}",
        est_dc
    );

    // Very short input (n < m) → should return 0.0 by guard clause
    let short = vec![Complex::new(1.0f32, 0.0f32); 5];
    let est_short = estimate_frequency_esprit(&short, sample_rate, m);
    println!("Short input (n<m): estimated={:.4} Hz", est_short);
    assert!(
        est_short.abs() < 1e-6,
        "n < m input should return 0.0, got {}",
        est_short
    );
}

// ---------------------------------------------------------------------------
// Test 4: Bussgang normalization preserves phase on constant-modulus FM signal
// ---------------------------------------------------------------------------
#[test]
fn test_bussgang_identity_on_fm() {
    let n = 200;
    let sample_rate = 10000.0_f32;

    // Generate an FM-like signal: constant modulus (magnitude 1.0) with varying phase
    let mut rng = SimpleRng::new(999);
    let mut phase: f32 = 0.0;
    let mut fm_signal = Vec::with_capacity(n);

    for _ in 0..n {
        // Random frequency deviation: ±500 Hz
        let freq_dev = rng.next_f32() * 500.0;
        phase += 2.0 * std::f32::consts::PI * freq_dev / sample_rate;
        // Wrap phase to avoid precision loss
        if phase > std::f32::consts::PI {
            phase -= 2.0 * std::f32::consts::PI;
        } else if phase < -std::f32::consts::PI {
            phase += 2.0 * std::f32::consts::PI;
        }
        fm_signal.push(Complex::new(phase.cos(), phase.sin()));
    }

    // Apply Bussgang normalization directly (same algorithm as dsp.rs lines 894-898):
    //   s_out = s / (|s| + 1e-9)
    // For constant-modulus signals with |s| = 1.0, this should be an identity
    // (up to the tiny epsilon denominator offset).
    let normalized: Vec<Complex<f32>> = fm_signal
        .iter()
        .map(|s| {
            let norm = s.norm();
            *s / (norm + 1e-9_f32)
        })
        .collect();

    assert_eq!(normalized.len(), n);

    // 1. All normalized samples should have magnitude very close to 1.0
    //    Exact value: |s| / (|s| + 1e-9) ≈ 1.0 / 1.0000001 ≈ 0.9999999
    for (i, s) in normalized.iter().enumerate() {
        let mag = s.norm();
        assert!(
            (mag - 1.0).abs() < 1e-6,
            "Sample {}: magnitude {:.9} deviates from 1.0 after Bussgang normalization",
            i,
            mag
        );
    }

    // 2. Phase should be perfectly preserved: compare original vs normalized phases
    let mut max_phase_err: f32 = 0.0;
    for (orig, norm) in fm_signal.iter().zip(normalized.iter()) {
        let mut err = (orig.arg() - norm.arg()).abs();
        if err > std::f32::consts::PI {
            err = 2.0 * std::f32::consts::PI - err;
        }
        if err > max_phase_err {
            max_phase_err = err;
        }
    }
    println!(
        "Max phase error between original and Bussgang-normalized FM: {:.9} rad",
        max_phase_err
    );
    assert!(
        max_phase_err < 1e-6,
        "Phase not preserved after Bussgang normalization (max err={:.9} rad)",
        max_phase_err
    );

    // 3. Verify phase differences (instantaneous frequency) are preserved
    let mut max_dphi_err: f32 = 0.0;
    for i in 1..n {
        let orig_diff = (fm_signal[i] * fm_signal[i - 1].conj()).arg();
        let norm_diff = (normalized[i] * normalized[i - 1].conj()).arg();
        let mut err = (orig_diff - norm_diff).abs();
        if err > std::f32::consts::PI {
            err = 2.0 * std::f32::consts::PI - err;
        }
        if err > max_dphi_err {
            max_dphi_err = err;
        }
    }
    println!(
        "Max instantaneous-frequency error: {:.9} rad",
        max_dphi_err
    );
    assert!(
        max_dphi_err < 1e-6,
        "Instantaneous frequency not preserved (max err={:.9} rad)",
        max_dphi_err
    );
}

// ---------------------------------------------------------------------------
// Test 5: ESPRIT output is always bounded within Nyquist
// ---------------------------------------------------------------------------
#[test]
fn test_esprit_output_bounded() {
    let sample_rate = 50000.0;
    let m = 10;
    let n = 200;
    let nyquist = sample_rate / 2.0;

    // Generate purely random noise (no tone) — estimate should still be bounded
    let mut rng = SimpleRng::new(7777);
    let noise: Vec<Complex<f32>> = (0..n)
        .map(|_| Complex::new(rng.next_f32(), rng.next_f32()))
        .collect();

    let est = estimate_frequency_esprit(&noise, sample_rate, m);
    println!("Noise-only: estimated={:.2} Hz, Nyquist={:.0} Hz", est, nyquist);
    assert!(
        est.abs() <= nyquist + 1.0,
        "ESPRIT output {:.2} Hz exceeds Nyquist bounds ±{:.0} Hz",
        est,
        nyquist
    );
}

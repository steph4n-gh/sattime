use num_complex::Complex;
use rustfft::FftPlanner;

/// Compute the intensity scintillation index S4 from a window of signal amplitudes
pub fn compute_s4(amplitudes: &[f64]) -> f64 {
    if amplitudes.is_empty() {
        return 0.0;
    }
    let mut sum_i = 0.0;
    let mut sum_i2 = 0.0;
    for &a in amplitudes {
        let i = a * a;
        sum_i += i;
        sum_i2 += i * i;
    }
    let n = amplitudes.len() as f64;
    let mean_i = sum_i / n;
    let mean_i2 = sum_i2 / n;
    
    if mean_i > 1e-12 {
        let variance_i = mean_i2 - mean_i * mean_i;
        if variance_i > 0.0 {
            (variance_i.sqrt() / mean_i).min(2.0)
        } else {
            0.0
        }
    } else {
        0.0
    }
}

/// Compute the phase scintillation index sigma_phi from a window of phase innovations (radians)
pub fn compute_sigma_phi(phase_innovations: &[f64]) -> f64 {
    if phase_innovations.is_empty() {
        return 0.0;
    }
    let mut sum_p = 0.0;
    let mut sum_p2 = 0.0;
    for &p in phase_innovations {
        sum_p += p;
        sum_p2 += p * p;
    }
    let n = phase_innovations.len() as f64;
    let mean_p = sum_p / n;
    let mean_p2 = sum_p2 / n;
    let variance_p = mean_p2 - mean_p * mean_p;
    if variance_p > 0.0 {
        variance_p.sqrt()
    } else {
        0.0
    }
}

/// Perform cepstral analysis on the phase innovations to identify periodic tumbling.
/// Returns the detected tumbling frequency (Hz) and peak prominence/magnitude.
pub fn analyze_attitude(history: &[f64], fs: f64) -> Option<(f64, f64)> {
    let n = history.len();
    if n < 256 {
        return None;
    }

    // Use a power-of-two FFT size <= history length
    let mut fft_size = 256;
    while fft_size * 2 <= n {
        fft_size *= 2;
    }
    let fft_size = fft_size.min(1024);
    if history.len() < fft_size {
        return None;
    }

    let start_idx = history.len() - fft_size;
    let window_data = &history[start_idx..];

    // Forward FFT
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_size);
    let mut buffer: Vec<Complex<f64>> = window_data
        .iter()
        .map(|&x| Complex::new(x, 0.0))
        .collect();
    fft.process(&mut buffer);

    // Compute log magnitude of the spectrum
    let mut log_mag: Vec<Complex<f64>> = buffer
        .iter()
        .map(|&c| Complex::new((c.norm() + 1e-12).ln(), 0.0))
        .collect();

    // Inverse FFT (IFFT) to get the Cepstrum
    let ifft = planner.plan_fft_inverse(fft_size);
    ifft.process(&mut log_mag);

    // Normalize IFFT output
    for c in &mut log_mag {
        *c = *c / (fft_size as f64);
    }

    // Search for the maximum peak in the real part of the Cepstrum,
    // excluding the DC and low-quefrency region (below 1.0 second period)
    // and up to fft_size / 2 (symmetric)
    let min_bin = (fs * 1.0) as usize;
    let max_bin = fft_size / 2;
    if min_bin >= max_bin {
        return None;
    }

    let mut best_bin = min_bin;
    let mut best_val = -1e9;
    for bin in min_bin..max_bin {
        let val = log_mag[bin].re;
        if val > best_val {
            best_val = val;
            best_bin = bin;
        }
    }

    // tumbling frequency = fs / bin
    let t_period = (best_bin as f64) / fs;
    if t_period > 0.0 {
        let freq = 1.0 / t_period;
        Some((freq, best_val))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scintillation_indices() {
        let amplitudes = vec![1.0, 1.1, 0.9, 1.0, 1.2, 0.8];
        let s4 = compute_s4(&amplitudes);
        assert!(s4 > 0.0 && s4 < 1.0);

        let phase_errors = vec![0.0, 0.1, -0.1, 0.05, -0.05];
        let sigma = compute_sigma_phi(&phase_errors);
        assert!(sigma > 0.0 && sigma < 0.2);
    }

    #[test]
    fn test_cepstral_attitude_analysis() {
        // Generate a 0.2 Hz tumbling modulation sampled at 50 Hz
        let fs = 50.0;
        let mut history = Vec::new();
        for i in 0..1024 {
            let t = (i as f64) / fs;
            let mut val = 0.0;
            // Sum of 5 harmonics to create a periodic ripple in the spectrum
            for k in 1..=5 {
                val += (2.0 * std::f64::consts::PI * (0.2 * k as f64) * t).cos();
            }
            history.push(val);
        }

        let result = analyze_attitude(&history, fs);
        assert!(result.is_some());
        let (freq, prominence) = result.unwrap();
        // Freq should be very close to 0.2 Hz
        assert!((freq - 0.2).abs() < 0.05);
        assert!(prominence > 0.0);
    }
}

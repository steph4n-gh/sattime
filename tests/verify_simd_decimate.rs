use num_complex::Complex;
use sattime::dsp::{FirDecimator, design_lowpass_filter};

#[test]
fn test_simd_decimate_equivalence() {
    let sample_rate = 2000000.0;
    let taps = design_lowpass_filter(200000.0, sample_rate, 127);
    
    // Create two identical decimators
    let mut decimator = FirDecimator::new(taps.clone(), 4);
    
    // Generate simulated signal
    let mut input = Vec::new();
    for i in 0..5000 {
        let re = (i as f32 * 0.05).sin() + (i as f32 * 0.12).cos() * 0.5;
        let im = (i as f32 * 0.07).cos() - (i as f32 * 0.15).sin() * 0.3;
        input.push(Complex::new(re, im));
    }
    
    // Compute using the standard decimate process which uses SIMD (under compute)
    let mut output_simd = Vec::new();
    decimator.process(&input, &mut output_simd);
    
    // Manually verify equivalence on sliding windows using compute_scalar
    let hist_len = taps.len() - 1;
    let mut idx = 0;
    let total_len = hist_len + input.len();
    
    let history = vec![Complex::new(0.0, 0.0); hist_len];
    let mut output_scalar = Vec::new();
    
    // Boundary phase scalar
    while idx < hist_len && idx + taps.len() <= total_len {
        let mut sum = Complex::new(0.0, 0.0);
        for n in 0..taps.len() {
            let sample_idx = idx + n;
            let sample = if sample_idx < hist_len {
                history[sample_idx]
            } else {
                input[sample_idx - hist_len]
            };
            sum += sample * taps[n];
        }
        output_scalar.push(sum);
        idx += 4;
    }
    
    // Main phase scalar
    while idx + taps.len() <= total_len {
        let input_offset = idx - hist_len;
        let window = &input[input_offset..input_offset + taps.len()];
        let val = decimator.compute_scalar(window);
        output_scalar.push(val);
        idx += 4;
    }
    
    assert_eq!(output_simd.len(), output_scalar.len(), "Outputs length mismatch");
    
    for i in 0..output_simd.len() {
        let diff_re = (output_simd[i].re - output_scalar[i].re).abs();
        let diff_im = (output_simd[i].im - output_scalar[i].im).abs();
        
        assert!(
            diff_re < 1e-5,
            "Real mismatch at index {}: SIMD {} vs Scalar {}, diff = {}",
            i, output_simd[i].re, output_scalar[i].re, diff_re
        );
        assert!(
            diff_im < 1e-5,
            "Imag mismatch at index {}: SIMD {} vs Scalar {}, diff = {}",
            i, output_simd[i].im, output_scalar[i].im, diff_im
        );
    }
}

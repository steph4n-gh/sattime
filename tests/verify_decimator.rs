use num_complex::Complex;
use sattime::dsp::FirDecimator;

#[test]
fn test_decimator_continuous_streaming() {
    let taps = vec![0.1, 0.2, 0.4, 0.2, 0.1];
    let decimation_factor = 3;

    // Generate random-like input
    let mut input = Vec::new();
    for i in 0..1000 {
        let re = (i as f32 * 0.1).sin();
        let im = (i as f32 * 0.15).cos();
        input.push(Complex::new(re, im));
    }

    // Process in a single block
    let mut single_dec = FirDecimator::new(taps.clone(), decimation_factor);
    let mut single_output = Vec::new();
    single_dec.process(&input, &mut single_output);

    // Process in chunks of different sizes
    let chunk_sizes = vec![1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377];
    let mut chunk_dec = FirDecimator::new(taps.clone(), decimation_factor);
    let mut chunk_output = Vec::new();

    let mut start = 0;
    let mut chunk_idx = 0;
    while start < input.len() {
        let size = chunk_sizes[chunk_idx % chunk_sizes.len()];
        let end = (start + size).min(input.len());
        let chunk = &input[start..end];

        let mut temp_out = Vec::new();
        chunk_dec.process(chunk, &mut temp_out);
        chunk_output.extend(temp_out);

        start = end;
        chunk_idx += 1;
    }

    assert_eq!(
        single_output.len(),
        chunk_output.len(),
        "Output lengths mismatch"
    );
    for i in 0..single_output.len() {
        let diff_re = (single_output[i].re - chunk_output[i].re).abs();
        let diff_im = (single_output[i].im - chunk_output[i].im).abs();
        assert!(
            diff_re < 1e-5,
            "Mismatch at index {} re: {} vs {}",
            i,
            single_output[i].re,
            chunk_output[i].re
        );
        assert!(
            diff_im < 1e-5,
            "Mismatch at index {} im: {} vs {}",
            i,
            single_output[i].im,
            chunk_output[i].im
        );
    }
}

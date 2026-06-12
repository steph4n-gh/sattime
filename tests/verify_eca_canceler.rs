use num_complex::Complex;
use sattime::dsp::EcaCanceler;

#[test]
fn test_eca_clutter_suppression() {
    let mut canceler = EcaCanceler::new();
    let mut input = vec![Complex::new(0.0, 0.0); 1000];

    // Create a strong direct path carrier wave (sine/cosine wave)
    for i in 0..1000 {
        let val = Complex::new((i as f32 * 0.02).cos(), (i as f32 * 0.02).sin());
        input[i] = val * 1000.0; // Massive direct path leakage
    }

    let mut clean_output = vec![Complex::new(0.0, 0.0); 1000];
    canceler.process_block(&input, &mut clean_output);

    // Verify that after the initial filter transient, the massive direct path is suppressed
    for i in 50..1000 {
        let mag = (clean_output[i].re * clean_output[i].re + clean_output[i].im * clean_output[i].im).sqrt();
        assert!(
            mag < 0.5,
            "Direct path not suppressed at index {}: got magnitude {}",
            i, mag
        );
    }
}

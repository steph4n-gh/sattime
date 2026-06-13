use num_complex::Complex;
use sattime::dsp::{EcaCanceler, clean_ambiguity_map};

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

#[test]
fn test_clean_algorithm_omp() {
    // Generate a simple 10x10 ambiguity map with a large target and a small target
    let mut map = vec![vec![0.0f32; 10]; 10];
    map[3][4] = 100.0; // Large airliner target
    map[6][7] = 45.0;  // Small drone target

    // Add some sidelobes from airliner using Gaussian spread
    for r in 0..10 {
        let dr = (r as f32 - 3.0).powi(2);
        for c in 0..10 {
            let dc = (c as f32 - 4.0).powi(2);
            map[r][c] += 100.0 * (-dr/8.0 - dc/8.0).exp();
        }
    }
    // Set exact peak values again
    map[3][4] = 100.0;
    map[6][7] = 45.0;

    let components = clean_ambiguity_map(&mut map, 2, 0.8);
    
    assert_eq!(components.len(), 2);
    // First component should be airliner at (3, 4)
    assert_eq!(components[0].0, 3);
    assert_eq!(components[0].1, 4);
    assert!(components[0].2 > 90.0);

    // Second component should be drone at (6, 7)
    assert_eq!(components[1].0, 6);
    assert_eq!(components[1].1, 7);
}

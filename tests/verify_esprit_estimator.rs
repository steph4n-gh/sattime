use num_complex::Complex;
use sattime::dsp::estimate_frequency_esprit;

#[test]
fn test_esprit_estimator_accuracy() {
    let sample_rate = 50000.0;
    let f0 = 1234.5;
    let n = 500;
    let m = 10;
    
    // Simple LCG for deterministic noise generation
    let mut state: u64 = 1337;
    let mut next_noise = || -> f32 {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        // map to [-1.0, 1.0]
        let val = (state >> 32) as f32 / (u32::MAX as f32);
        val * 2.0 - 1.0
    };

    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let t = (i as f64) / sample_rate;
        let angle = 2.0 * std::f64::consts::PI * f0 * t;
        let signal = Complex::new(angle.cos() as f32, angle.sin() as f32);
        
        // Add slight noise: SNR is high, so noise amplitude is small, say 0.01
        let noise = Complex::new(next_noise() * 0.01, next_noise() * 0.01);
        
        samples.push(signal + noise);
    }
    
    let estimated = estimate_frequency_esprit(&samples, sample_rate, m);
    println!("Estimated frequency: {}, Target: {}", estimated, f0);
    
    assert!(
        (estimated.abs() - f0).abs() < 1.0,
        "Estimated frequency magnitude {} is not within 1.0 Hz of target {}",
        estimated.abs(),
        f0
    );
}

#[test]
fn test_esprit_estimator_signed() {
    let sample_rate = 50000.0;
    let m = 10;
    let n = 500;

    for &f0 in &[1500.0, -1500.0] {
        let mut samples = Vec::with_capacity(n);
        for i in 0..n {
            let t = (i as f64) / sample_rate;
            let angle = 2.0 * std::f64::consts::PI * f0 * t;
            let signal = Complex::new(angle.cos() as f32, angle.sin() as f32);
            samples.push(signal);
        }
        let estimated = estimate_frequency_esprit(&samples, sample_rate, m);
        println!("Estimated signed frequency: {}, Target: {}", estimated, f0);
        assert!(
            (estimated - f0).abs() < 1.0,
            "Estimated frequency {} is not within 1.0 Hz of target {}",
            estimated,
            f0
        );
    }
}

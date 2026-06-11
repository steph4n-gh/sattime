use num_complex::Complex;
use sattime::dsp::{DemodChannel, Modulation};

#[test]
fn test_esprit_bootstrap_carrier() {
    let nominal_freq = 1000.0;
    let sample_rate = 10000.0;
    // 200 Hz frequency offset from nominal frequency
    let target_freq = 1200.0;
    
    let mut channel = DemodChannel::new(nominal_freq, sample_rate, "SAT_TEST_CARRIER".to_string());
    
    // Generate samples
    let n = 2000;
    let mut input = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / sample_rate;
        let angle = 2.0 * std::f64::consts::PI * target_freq * t;
        input.push(Complex::new(angle.cos() as f32, angle.sin() as f32));
    }
    
    // Process block
    channel.process_block(&input);
    
    // Verify ESPRIT bootstrap initialized the EKF frequency estimate
    let x_freq = channel.pll_tracker.x[1];
    let est_offset_hz = x_freq / (2.0 * std::f64::consts::PI);
    
    println!("Carrier test: target offset = 200 Hz, ESPRIT initialized EKF offset = {} Hz", est_offset_hz);
    
    assert!(
        (est_offset_hz - 200.0).abs() < 10.0,
        "Carrier ESPRIT bootstrap frequency offset {} is not within 10 Hz of target 200 Hz",
        est_offset_hz
    );
}

#[test]
fn test_esprit_bootstrap_bpsk() {
    let nominal_freq = 1000.0;
    let sample_rate = 10000.0;
    // 200 Hz frequency offset from nominal frequency
    let target_freq = 1200.0;
    
    // For BPSK, we create the channel using new_prod to set Modulation::Bpsk
    let taps = sattime::dsp::design_lowpass_filter(0.4 * sample_rate, sample_rate, 127);
    let mut channel = DemodChannel::new_prod(
        1,
        sample_rate,
        1000.0, // sym_rate
        Modulation::Bpsk,
        false,
        false,
        false,
        true, // no_multihypothesis
        3.0,
        5.0,
        taps,
        4, // decimate_factor
    );
    channel.nominal_freq = nominal_freq;
    channel.sample_rate = sample_rate;
    channel.sat_name = "SAT_TEST_BPSK".to_string();
    channel.frequency = nominal_freq;
    channel.status = sattime::dsp::ChannelStatus::Acquisition;
    
    // Generate BPSK samples
    let n = 2000;
    let mut input = Vec::with_capacity(n);
    
    // Alternating/random symbols every 10 samples
    let mut current_symbol = 1.0f32;
    for i in 0..n {
        if i % 10 == 0 {
            current_symbol = if (i / 10) % 2 == 0 { 1.0 } else { -1.0 };
        }
        let t = i as f64 / sample_rate;
        let angle = 2.0 * std::f64::consts::PI * target_freq * t;
        input.push(Complex::new(
            current_symbol * angle.cos() as f32,
            current_symbol * angle.sin() as f32,
        ));
    }
    
    // Process block
    channel.process_block(&input);
    
    // Verify ESPRIT bootstrap initialized the EKF frequency estimate to twice the offset
    let x_freq = channel.pll_tracker.x[1];
    let est_offset_hz = x_freq / (2.0 * std::f64::consts::PI);
    
    println!("BPSK test: target offset (doubled) = 400 Hz, ESPRIT initialized EKF offset = {} Hz", est_offset_hz);
    
    assert!(
        (est_offset_hz - 400.0).abs() < 20.0,
        "BPSK ESPRIT bootstrap frequency offset {} is not within 20 Hz of target 400 Hz",
        est_offset_hz
    );
}

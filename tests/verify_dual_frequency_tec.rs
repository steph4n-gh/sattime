use num_complex::Complex;
use sattime::dsp::{DemodChannel, Modulation, ChannelStatus};

#[test]
fn test_dual_frequency_tec_tracking() {
    let nominal_freq = 1000.0;
    let frequency2 = 1500.0;
    let sample_rate = 10000.0;
    let decimate_factor = 2;
    let sym_rate = sample_rate / decimate_factor as f64; // 5000 Hz decimated rate

    // Target frequencies with Doppler shift (50 Hz and 75 Hz)
    let target_freq1 = 1050.0;
    let target_freq2 = 1575.0;

    // Design a narrowband filter with 150 Hz cutoff to isolate f1 (+50 Hz) from f2 (+575 Hz) and vice versa.
    let taps = sattime::dsp::design_lowpass_filter(150.0, sample_rate, 127);
    let mut channel = DemodChannel::new_prod(
        0,
        sample_rate,
        sym_rate,
        Modulation::Carrier,
        false, // no_adaptive_ekf
        false, // no_dual_lock
        true,  // no_gardner = true
        true,  // no_multihypothesis = true
        1.0,   // min_snr
        5.0,   // fade_timeout
        taps,
        decimate_factor,
        false, // eca_enabled
    );

    // Set channel status to Acquisition so it is processed
    channel.status = ChannelStatus::Acquisition;

    // Configure dual frequency
    channel.is_dual = true;
    channel.nominal_freq = nominal_freq;
    channel.frequency2 = frequency2;
    channel.target_freq = nominal_freq;
    channel.target_freq2 = frequency2;
    channel.pll_tracker.set_frequency_ratio(frequency2 / nominal_freq);

    // Reset tracker with initial carrier 1 Doppler guess of 50.0 Hz.
    // The EKF reset will automatically scale carrier 2 Doppler to 50.0 * 1.5 = 75.0 Hz.
    channel.pll_tracker.reset(0.0, 50.0, 0.0);

    // Generate dual-frequency samples (we sum them together with 0.5 amplitude to prevent clipping)
    let n = 6000;
    let mut raw_iq = Vec::with_capacity(n);
    
    // Simulating some dispersive phase shift (e.g. 0.5 rad on f1, 0.2 rad on f2)
    let phase_offset1 = 0.5;
    let phase_offset2 = 0.2;

    for i in 0..n {
        let t = i as f64 / sample_rate;
        let angle1 = 2.0 * std::f64::consts::PI * target_freq1 * t + phase_offset1;
        let angle2 = 2.0 * std::f64::consts::PI * target_freq2 * t + phase_offset2;
        // Combine the two carriers with scaling
        let s = Complex::new(
            0.5 * (angle1.cos() as f32 + angle2.cos() as f32),
            0.5 * (angle1.sin() as f32 + angle2.sin() as f32),
        );
        raw_iq.push(s);
    }

    // Process the block of samples in chunks
    for (chunk_idx, chunk) in raw_iq.chunks(1000).enumerate() {
        channel.process_block(chunk);
        let f1_est = channel.pll_tracker.x[1] / (2.0 * std::f64::consts::PI);
        let f2_est = channel.pll_tracker.x[4] / (2.0 * std::f64::consts::PI);
        println!(
            "Chunk {}: is_locked={}, lock_metric={:.4}, x[1] (f1_est)={:.2} Hz, x[4] (f2_est)={:.2} Hz, current_tec={:.6e}",
            chunk_idx,
            channel.is_locked,
            channel.pll_tracker.lock_metric,
            f1_est,
            f2_est,
            channel.current_tec
        );
    }

    let est_freq = channel.frequency;
    let tec = channel.current_tec;

    println!("Dual frequency test results:");
    println!("  Tracked Ionosphere-Free Frequency: {:.2} Hz", est_freq);
    println!("  Estimated TEC: {:.6e} TECU", tec);

    // Check if tracking loop locked
    assert!(channel.is_locked, "DemodChannel failed to lock in dual-frequency mode");

    // Under the cancelation formula:
    // f_free = (f1_sq * target_freq1 - f2_sq * target_freq2) / (f1_sq - f2_sq)
    // For f1=1000, f2=1500, target_freq1=1050, target_freq2=1575:
    // f_free = (1000^2 * 1050 - 1500^2 * 1575) / (1000^2 - 1500^2) = 1995.0
    let expected_f_free = 1995.0;
    assert!(
        (est_freq - expected_f_free).abs() < 2.0,
        "Tracked frequency {:.2} is not within 2 Hz of expected {:.2} Hz",
        est_freq,
        expected_f_free
    );

    // Verify TEC is non-zero and reasonable (positive)
    assert!(tec > 0.0, "Estimated TEC should be positive");
}

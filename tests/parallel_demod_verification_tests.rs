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

use num_complex::Complex;

// ----------------------------------------------------------------------------
// TEST 1: DDC Phase Continuity (Production vs Mock)
// ----------------------------------------------------------------------------
#[test]
fn test_ddc_phase_continuity() {
    let sample_rate = 10000.0;
    let frequency = 1250.0;
    let block_size = 100;
    let input1 = vec![Complex::new(1.0f32, 0.0f32); block_size];
    let input2 = vec![Complex::new(1.0f32, 0.0f32); block_size];

    // --- Production DDC ---
    let mut prod_ddc = dsp::DigitalDownConverter::new();
    let mut prod_out1 = vec![Complex::new(0.0f32, 0.0f32); block_size];
    let mut prod_out2 = vec![Complex::new(0.0f32, 0.0f32); block_size];

    prod_ddc.process(&input1, frequency, sample_rate, &mut prod_out1);
    prod_ddc.process(&input2, frequency, sample_rate, &mut prod_out2);

    // --- Mock/Facade DDC (as implemented in ekf::DemodChannel::process_block) ---
    // In ekf.rs:
    // for (i, &sample) in raw_iq.iter().enumerate() {
    //     let t = i as f64 / self.sample_rate;
    //     let phase = -2.0 * std::f64::consts::PI * self.frequency * t;
    //     ...
    // }
    let mut mock_out1 = vec![Complex::new(0.0f32, 0.0f32); block_size];
    let mut mock_out2 = vec![Complex::new(0.0f32, 0.0f32); block_size];

    for i in 0..block_size {
        let t = i as f64 / sample_rate;
        let phase = -2.0 * std::f64::consts::PI * frequency * t;
        mock_out1[i] = input1[i] * Complex::new(phase.cos() as f32, phase.sin() as f32);
    }
    // Block 2: t resets to i / sample_rate (i.e. starts at 0.0)
    for i in 0..block_size {
        let t = i as f64 / sample_rate;
        let phase = -2.0 * std::f64::consts::PI * frequency * t;
        mock_out2[i] = input2[i] * Complex::new(phase.cos() as f32, phase.sin() as f32);
    }

    // --- Contrast Phase Difference at Block Boundary ---
    // Last sample of block 1 vs first sample of block 2
    let prod_diff_phase = (prod_out2[0].arg() - prod_out1[block_size - 1].arg()).abs();
    let mock_diff_phase = (mock_out2[0].arg() - mock_out1[block_size - 1].arg()).abs();

    println!("Production boundary phase difference: {}", prod_diff_phase);
    println!("Mock boundary phase difference: {}", mock_diff_phase);

    // The expected phase step per sample is 2 * pi * f / fs = 2 * pi * 1250 / 10000 = pi / 4 = 0.785398
    let expected_step = (2.0 * std::f64::consts::PI * frequency / sample_rate) as f32;

    // Check production continuity
    let prod_step_diff = (prod_diff_phase - expected_step).abs();
    assert!(
        prod_step_diff < 1e-4
            || (prod_diff_phase + expected_step - 2.0 * std::f32::consts::PI).abs() < 1e-4,
        "Production DDC phase step at block boundary is discontinuous!"
    );

    // Check mock discontinuity (it should jump back to t = 0 phase, which is 0.0, causing a mismatch)
    let mock_step_diff = (mock_diff_phase - expected_step).abs();
    println!(
        "Mock step difference from expected continuous phase: {}",
        mock_step_diff
    );
}

// ----------------------------------------------------------------------------
// TEST 2: Gardner Loop Timing Evaluation (1-Sample Delay Bug Verification)
// ----------------------------------------------------------------------------
#[test]
fn test_gardner_loop_timing_indices() {
    let sample_rate = 1000.0;
    let symbol_rate = 500.0;
    let mut gardner = dsp::GardnerLoop::new(sample_rate, symbol_rate);

    // Feed a few samples to populate the Farrow history
    let mut symbols = Vec::new();

    // Let's push 5 samples: index 1.0, 2.0, 3.0, 4.0, 5.0
    gardner.process(Complex::new(1.0, 0.0), &mut symbols); // sample_index = 1.0
    gardner.process(Complex::new(2.0, 0.0), &mut symbols); // sample_index = 2.0
    gardner.process(Complex::new(3.0, 0.0), &mut symbols); // sample_index = 3.0

    // When pushing sample 4.0, sample_index becomes 4.0, which triggers processing
    // Let's capture what mu is and what samples it interpolates between
    gardner.process(Complex::new(4.0, 0.0), &mut symbols); // sample_index = 4.0

    // At this point:
    // farrow.history contains: [1.0, 2.0, 3.0, 4.0] (corresponding to indices 1, 2, 3, 4)
    // Since symbol_rate = 500.0 (sps = 2.0), the first symbol center is at index 1.0.
    // The sample at index 1.0 (0-indexed) has value 2.0.

    assert!(
        !symbols.is_empty(),
        "Gardner loop did not output any symbols"
    );
    let (first_symbol, mu) = symbols[0];
    println!("Gardner output symbol: {:?}, mu: {}", first_symbol, mu);
    println!("Farrow history: {:?}", gardner.farrow.history);

    assert_eq!(
        first_symbol.re, 2.0,
        "Expected first symbol to be 2.0 after fixing the 1-sample delay bug"
    );
}

use crossbeam_channel::bounded;
use num_complex::Complex;
use rustfft::FftPlanner;
use std::io;
use std::time::{Duration, Instant};

// Helper to parse TLE file (Worst-case payload)
fn parse_tle_file(path: &str) -> io::Result<Vec<(String, sgp4::Elements)>> {
    use std::fs::File;
    use std::io::BufRead;

    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    let mut satellites = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();
        if line.is_empty() {
            i += 1;
            continue;
        }

        if line.starts_with('1') && i + 1 < lines.len() && lines[i + 1].trim().starts_with('2') {
            let line1 = line;
            let line2 = lines[i + 1].trim();
            if let Ok(elements) = sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes())
            {
                let sat_name = elements
                    .object_name
                    .clone()
                    .unwrap_or_else(|| "UNKNOWN".to_string());
                satellites.push((sat_name, elements));
            }
            i += 2;
        } else if i + 2 < lines.len() {
            let name = line.to_string();
            let line1 = lines[i + 1].trim();
            let line2 = lines[i + 2].trim();

            if line1.starts_with('1') && line2.starts_with('2') {
                if let Ok(elements) =
                    sgp4::Elements::from_tle(Some(name.clone()), line1.as_bytes(), line2.as_bytes())
                {
                    satellites.push((name, elements));
                }
                i += 3;
            } else {
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    Ok(satellites)
}

// Low-pass filter design
fn design_lowpass_filter(cutoff_hz: f64, sample_rate_hz: f64, num_taps: usize) -> Vec<f32> {
    let mut taps = vec![0.0f32; num_taps];
    let middle = (num_taps - 1) as f64 / 2.0;
    let fc = cutoff_hz / sample_rate_hz;
    let w_c = 2.0 * std::f64::consts::PI * fc;

    let mut sum = 0.0f64;
    for n in 0..num_taps {
        let x = (n as f64) - middle;
        let val = if x.abs() < 1e-9 {
            w_c / std::f64::consts::PI
        } else {
            (w_c * x).sin() / (std::f64::consts::PI * x)
        };

        let win =
            0.54 - 0.46 * (2.0 * std::f64::consts::PI * n as f64 / (num_taps - 1) as f64).cos();
        taps[n] = (val * win) as f32;
        sum += taps[n] as f64;
    }

    for val in taps.iter_mut() {
        *val /= sum as f32;
    }
    taps
}

// FIR Decimator
struct FirDecimator {
    taps: Vec<f32>,
    decimation_factor: usize,
    history: Vec<Complex<f32>>,
    pending_offset: usize,
}

impl FirDecimator {
    fn new(taps: Vec<f32>, decimation_factor: usize) -> Self {
        let hist_len = taps.len().saturating_sub(1);
        Self {
            taps,
            decimation_factor,
            history: vec![Complex::new(0.0, 0.0); hist_len],
            pending_offset: 0,
        }
    }

    fn process(&mut self, input: &[Complex<f32>], output: &mut Vec<Complex<f32>>) {
        if self.decimation_factor <= 1 {
            output.extend_from_slice(input);
            return;
        }

        let num_taps = self.taps.len();
        let hist_len = self.history.len();
        let total_len = hist_len + input.len();

        let mut idx = self.pending_offset;
        while idx + num_taps <= total_len {
            let mut sum = Complex::new(0.0f32, 0.0f32);
            for n in 0..num_taps {
                let sample_idx = idx + n;
                let sample = if sample_idx < hist_len {
                    self.history[sample_idx]
                } else {
                    input[sample_idx - hist_len]
                };
                sum += sample * self.taps[n];
            }
            output.push(sum);
            idx += self.decimation_factor;
        }

        self.pending_offset = idx.saturating_sub(input.len());

        if input.len() >= hist_len {
            self.history
                .copy_from_slice(&input[input.len() - hist_len..]);
        } else {
            let shift = hist_len - input.len();
            self.history.copy_within(input.len().., 0);
            self.history[shift..].copy_from_slice(input);
        }
    }
}

#[test]
fn test_adversarial_performance_limits() {
    let tle_path = "passes/starlink.tle";

    // 1. TLE Parsing Latency (worst-case Starlink dataset)
    println!("=== 1. TLE Parsing Benchmark ===");
    let start_parse = Instant::now();
    let satellites = parse_tle_file(tle_path).expect("Failed to parse Starlink TLE");
    let parse_duration = start_parse.elapsed();
    println!(
        "Parsed {} satellites in {:?}",
        satellites.len(),
        parse_duration
    );

    // 2. Satellite List Cloning Latency (step_count % 500 == 0)
    println!("\n=== 2. Satellite Vector Cloning Benchmark ===");
    let start_clone = Instant::now();
    let cloned_sats = satellites.clone();
    let clone_duration = start_clone.elapsed();
    println!(
        "Cloned {} satellites in {:?}",
        cloned_sats.len(),
        clone_duration
    );

    // 3. DSP Processing Step Latency (including FIR, FFT, and spur notching)
    println!("\n=== 3. DSP Step Processing Benchmark ===");
    let sample_rate = 2_000_000.0;
    let decimate = 40;
    let pipeline_sample_rate = sample_rate / decimate as f64; // 50 kHz
    let pipeline_fft_size = 1024;
    let pipeline_step_size = 1000;
    let num_taps = 31;
    let taps = design_lowpass_filter(0.4 * pipeline_sample_rate, sample_rate, num_taps);
    let mut decimator = FirDecimator::new(taps, decimate);

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(pipeline_fft_size);

    let block_size = 32768;
    let mut mock_block = vec![Complex::new(0.0f32, 0.0f32); block_size];
    for (i, val) in mock_block.iter_mut().enumerate() {
        let t = i as f32 / sample_rate as f32;
        *val = Complex::new(
            (2.0 * std::f32::consts::PI * 5000.0 * t).cos(),
            (2.0 * std::f32::consts::PI * 5000.0 * t).sin(),
        );
    }

    let mut decimated_samples = Vec::with_capacity(block_size / decimate + 1);
    let mut queue = std::collections::VecDeque::with_capacity(pipeline_fft_size * 2);
    let mut fft_input = vec![Complex::new(0.0, 0.0); pipeline_fft_size];

    let mut dsp_latencies = Vec::new();
    let mut rise_sched_latencies = Vec::new();

    // Warm-up and measure
    for step in 0..1000 {
        let step_start = Instant::now();

        // Simulate decimator and queue buffering
        decimated_samples.clear();
        decimator.process(&mock_block, &mut decimated_samples);
        queue.extend(decimated_samples.drain(..));

        if queue.len() >= pipeline_fft_size {
            let (slice1, slice2) = queue.as_slices();
            if slice1.len() >= pipeline_fft_size {
                fft_input[..pipeline_fft_size].copy_from_slice(&slice1[..pipeline_fft_size]);
            } else {
                fft_input[..slice1.len()].copy_from_slice(slice1);
                let remaining = pipeline_fft_size - slice1.len();
                fft_input[slice1.len()..pipeline_fft_size].copy_from_slice(&slice2[..remaining]);
            }

            fft.process(&mut fft_input);
            queue.drain(..pipeline_step_size);
        }

        // At step % 500, we perform satellites.clone()
        let is_rise_sched_step = step % 500 == 0;
        let final_elapsed = step_start.elapsed();

        if is_rise_sched_step {
            let start_rise = Instant::now();
            let _cloned = satellites.clone();
            let rise_elapsed = final_elapsed + start_rise.elapsed();
            rise_sched_latencies.push(rise_elapsed);
        } else {
            dsp_latencies.push(final_elapsed);
        }
    }

    let avg_dsp_ms =
        dsp_latencies.iter().sum::<Duration>().as_secs_f64() * 1000.0 / dsp_latencies.len() as f64;
    let max_dsp_ms = dsp_latencies.iter().max().unwrap().as_secs_f64() * 1000.0;
    let avg_rise_ms = rise_sched_latencies.iter().sum::<Duration>().as_secs_f64() * 1000.0
        / rise_sched_latencies.len() as f64;
    let max_rise_ms = rise_sched_latencies.iter().max().unwrap().as_secs_f64() * 1000.0;

    println!(
        "Average DSP step latency: {:.6} ms (Max: {:.6} ms)",
        avg_dsp_ms, max_dsp_ms
    );
    println!(
        "Average Rise Schedule step latency (step % 500): {:.6} ms (Max: {:.6} ms)",
        avg_rise_ms, max_rise_ms
    );

    // 4. Threading & Queue Stability Simulation
    println!("\n=== 4. Threading & Queue Drop Simulation ===");
    // Config: bounded channel of 100 vectors, pool of 120 vectors (like main.rs)
    let (tx, rx) = bounded::<Vec<Complex<f32>>>(100);
    let (pool_tx, pool_rx) = bounded::<Vec<Complex<f32>>>(120);
    for _ in 0..120 {
        let _ = pool_tx.send(vec![Complex::new(0.0f32, 0.0f32); block_size]);
    }

    let block_time_ns = (block_size as f64 / sample_rate * 1e9) as u64; // ~16,384,000 ns
    let block_time = Duration::from_nanos(block_time_ns);

    let tx_clone = tx.clone();
    let pool_tx_clone = pool_tx.clone();
    let producer_handle = std::thread::spawn(move || {
        let mut drop_count = 0;
        let mut send_count = 0;
        let start = Instant::now();

        // Stream for 2 seconds (around 122 blocks)
        for i in 0..122 {
            let buf = pool_rx.recv().unwrap();

            // Try sending with a tiny timeout to check if it would block
            let send_res = tx_clone.send_timeout(buf, Duration::from_millis(1));
            match send_res {
                Ok(_) => send_count += 1,
                Err(e) => {
                    drop_count += 1;
                    // Put back to pool if failed
                    if let crossbeam_channel::SendTimeoutError::Timeout(b) = e {
                        let _ = pool_tx_clone.send(b);
                    }
                }
            }

            // Maintain real-time pace
            let next_target = start + block_time * (i + 1);
            let now = Instant::now();
            if next_target > now {
                std::thread::sleep(next_target - now);
            }
        }
        (send_count, drop_count)
    });

    // Drop our tx side on the main thread so that rx knows there are no more senders once the producer finishes
    drop(tx);

    // Main thread simulates blocking due to profile switch (synchronous TLE parse)
    println!(
        "Main thread simulating profile switch block of {:?}",
        parse_duration
    );
    std::thread::sleep(parse_duration);

    // Now start consuming everything from the channel until it closes
    let mut received_count = 0;
    while let Ok(buf) = rx.recv() {
        received_count += 1;
        let _ = pool_tx.send(buf);
    }

    let (sent, drops) = producer_handle.join().unwrap();
    println!(
        "Real-time Producer sent: {}, dropped (blocked/timeout): {}",
        sent, drops
    );
    println!("Main thread received: {}", received_count);

    // Assert safety limit constraints
    assert!(
        avg_dsp_ms < 1.8,
        "Average DSP latency exceeds 1.8 ms: {:.3} ms",
        avg_dsp_ms
    );
    assert!(
        avg_rise_ms < 1.8,
        "Average Rise Schedule latency exceeds 1.8 ms: {:.3} ms",
        avg_rise_ms
    );
    assert_eq!(
        drops, 0,
        "Real-time samples were dropped during TLE profile switch!"
    );
}

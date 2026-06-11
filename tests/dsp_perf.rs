use num_complex::Complex;
use rustfft::FftPlanner;
use std::time::Instant;

// Duplicate low-pass filter design
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

// Duplicate FirDecimator implementation
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
fn test_dsp_pipeline_performance() {
    let sample_rate = 2_000_000.0;
    let decimate = 40;
    let pipeline_sample_rate = sample_rate / decimate as f64; // 50,000 Hz
    let pipeline_fft_size = 1024;
    let pipeline_step_size = 1000;
    let num_taps = 31;
    let notch_spurs = true;

    // Design filter
    let taps = design_lowpass_filter(0.4 * pipeline_sample_rate, sample_rate, num_taps);
    let mut decimator = FirDecimator::new(taps, decimate);

    // Setup FFT planner
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(pipeline_fft_size);

    // Generate mock input (10 seconds of data = 20,000,000 samples)
    // We process in blocks of 32768 samples (typical SDR block size)
    let block_size = 32768;
    let total_samples = 20_000_000;
    let num_blocks = total_samples / block_size;

    let mut mock_block = vec![Complex::new(0.1f32, 0.2f32); block_size];
    // Add a carrier + some noise
    for (i, val) in mock_block.iter_mut().enumerate() {
        let t = i as f32 / sample_rate as f32;
        let re = (2.0 * std::f32::consts::PI * 5000.0 * t).cos() + 0.1;
        let im = (2.0 * std::f32::consts::PI * 5000.0 * t).sin() - 0.1;
        *val = Complex::new(re, im);
    }

    let mut decimated_samples = Vec::with_capacity(block_size / decimate + 1);
    let mut queue = std::collections::VecDeque::with_capacity(pipeline_fft_size * 2);
    let mut fft_input = vec![Complex::new(0.0, 0.0); pipeline_fft_size];

    // Search indices (+-20 kHz skipped DC)
    let max_k = ((20000.0 / pipeline_sample_rate) * pipeline_fft_size as f64).round() as usize;
    let mut search_indices = Vec::new();
    for k in 0..pipeline_fft_size {
        let is_in_band = k <= max_k || k >= pipeline_fft_size - max_k;
        if !is_in_band {
            continue;
        }
        let is_dc_region = k <= 5 || k >= pipeline_fft_size - 5;
        if is_dc_region {
            continue;
        }
        search_indices.push(k);
    }

    let mut fft_mag_ema = vec![0.0f32; pipeline_fft_size];
    let mut ema_initialized = false;
    let mut spurs = vec![false; pipeline_fft_size];
    let mut spur_consecutive_counts = vec![0u8; pipeline_fft_size];
    let mut step_count = 0;

    println!("Starting DSP pipeline performance benchmark (simulating 10s of 2 MSPS input)...");
    let start_time = Instant::now();

    for _ in 0..num_blocks {
        decimated_samples.clear();
        decimator.process(&mock_block, &mut decimated_samples);

        queue.extend(decimated_samples.drain(..));

        while queue.len() >= pipeline_fft_size {
            let (slice1, slice2) = queue.as_slices();
            if slice1.len() >= pipeline_fft_size {
                fft_input[..pipeline_fft_size].copy_from_slice(&slice1[..pipeline_fft_size]);
            } else {
                fft_input[..slice1.len()].copy_from_slice(slice1);
                let remaining = pipeline_fft_size - slice1.len();
                fft_input[slice1.len()..pipeline_fft_size].copy_from_slice(&slice2[..remaining]);
            }

            fft.process(&mut fft_input);

            if notch_spurs {
                let alpha = 0.001f32;
                if !ema_initialized {
                    for &k in &search_indices {
                        fft_mag_ema[k] = fft_input[k].norm_sqr();
                    }
                    ema_initialized = true;
                } else {
                    for &k in &search_indices {
                        fft_mag_ema[k] =
                            (1.0 - alpha) * fft_mag_ema[k] + alpha * fft_input[k].norm_sqr();
                    }
                }

                if step_count > 0 && step_count % 1000 == 0 {
                    let mut sorted = Vec::with_capacity(search_indices.len());
                    for &k in &search_indices {
                        sorted.push(fft_mag_ema[k]);
                    }
                    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    let median = sorted[sorted.len() / 2];
                    let threshold = 16.0 * median.max(0.0225f32);
                    for &k in &search_indices {
                        if fft_mag_ema[k] > threshold {
                            spur_consecutive_counts[k] =
                                spur_consecutive_counts[k].saturating_add(1);
                        } else {
                            spur_consecutive_counts[k] = 0;
                        }

                        if spur_consecutive_counts[k] >= 10 {
                            spurs[k] = true;
                        } else {
                            spurs[k] = false;
                        }
                    }
                }
            }

            let mut max_mag_sqr = -1.0f32;
            let mut k_max = 0;
            let mut power_sum = 0.0f32;
            let mut count_sum = 0;

            for &k in &search_indices {
                let val = fft_input[k];
                let mag_sqr = val.norm_sqr();
                if notch_spurs && spurs[k] {
                    continue;
                }
                power_sum += mag_sqr;
                count_sum += 1;

                if mag_sqr > max_mag_sqr {
                    max_mag_sqr = mag_sqr;
                    k_max = k;
                }
            }

            if max_mag_sqr >= 0.0 {
                let peak_power = max_mag_sqr;
                let noise_power = if count_sum > 1 {
                    (power_sum - peak_power) / (count_sum - 1) as f32
                } else {
                    1.0
                };
                let _snr_db = if noise_power > 1e-10 {
                    10.0 * (peak_power / noise_power).log10()
                } else {
                    0.0
                };

                let k_prev = (k_max + pipeline_fft_size - 1) % pipeline_fft_size;
                let k_next = (k_max + 1) % pipeline_fft_size;

                let y0 = fft_input[k_max].norm();
                let y_prev = fft_input[k_prev].norm();
                let y_next = fft_input[k_next].norm();

                let y0_log = (y0 + 1e-10).ln();
                let y_prev_log = (y_prev + 1e-10).ln();
                let y_next_log = (y_next + 1e-10).ln();

                let denom = y_prev_log - 2.0 * y0_log + y_next_log;
                let delta = if denom.abs() > 1e-6 {
                    (y_prev_log - y_next_log) / (2.0 * denom)
                } else {
                    0.0
                };
                let delta = delta.clamp(-0.5, 0.5);
                let _k_interp = (k_max as f64) + (delta as f64);
            }

            step_count += 1;
            // Advance the queue
            queue.drain(..pipeline_step_size);
        }
    }

    let elapsed = start_time.elapsed();
    let elapsed_sec = elapsed.as_secs_f64();
    let cpu_utilization = (elapsed_sec / 10.0) * 100.0;
    println!("Benchmark completed in {:.6} seconds.", elapsed_sec);
    println!("Simulated 10.0 seconds of IQ streaming at 2 MSPS.");
    println!(
        "Computed CPU utilization: {:.3}% on a single core.",
        cpu_utilization
    );
    assert!(
        cpu_utilization < 10.0,
        "CPU utilization is too high: {:.3}%",
        cpu_utilization
    );
}

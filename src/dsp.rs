use crate::daemon::*;
use crate::ekf::*;
use crate::orbit::*;
use crate::tui::*;
use chrono::{DateTime, Datelike, Timelike, Utc};
use num_complex::Complex;
use rustfft::FftPlanner;
use sgp4::Elements;
use std::collections::VecDeque;
use std::io::{self, Read, Write};
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Modulation {
    #[default]
    Carrier,
    Bpsk,
    Qpsk,
}

#[derive(Clone)]
pub struct FarrowInterpolator {
    pub history: [Complex<f32>; 4],
}

impl FarrowInterpolator {
    pub fn new() -> Self {
        Self {
            history: [Complex::new(0.0, 0.0); 4],
        }
    }

    pub fn push(&mut self, sample: Complex<f32>) {
        self.history[0] = self.history[1];
        self.history[1] = self.history[2];
        self.history[2] = self.history[3];
        self.history[3] = sample;
    }

    pub fn interpolate(&self, mu: f32) -> Complex<f32> {
        let y_neg1 = self.history[0];
        let y_0 = self.history[1];
        let y_1 = self.history[2];
        let y_2 = self.history[3];

        let v3_re = -1.0 / 6.0 * y_neg1.re + 0.5 * y_0.re - 0.5 * y_1.re + 1.0 / 6.0 * y_2.re;
        let v2_re = 0.5 * y_neg1.re - y_0.re + 0.5 * y_1.re;
        let v1_re = -1.0 / 3.0 * y_neg1.re - 0.5 * y_0.re + y_1.re - 1.0 / 6.0 * y_2.re;
        let v0_re = y_0.re;

        let v3_im = -1.0 / 6.0 * y_neg1.im + 0.5 * y_0.im - 0.5 * y_1.im + 1.0 / 6.0 * y_2.im;
        let v2_im = 0.5 * y_neg1.im - y_0.im + 0.5 * y_1.im;
        let v1_im = -1.0 / 3.0 * y_neg1.im - 0.5 * y_0.im + y_1.im - 1.0 / 6.0 * y_2.im;
        let v0_im = y_0.im;

        let re = ((v3_re * mu + v2_re) * mu + v1_re) * mu + v0_re;
        let im = ((v3_im * mu + v2_im) * mu + v1_im) * mu + v0_im;

        Complex::new(re, im)
    }

    pub fn reset(&mut self) {
        self.history = [Complex::new(0.0, 0.0); 4];
    }
}

#[derive(Clone)]
pub struct GardnerLoop {
    pub farrow: FarrowInterpolator,
    pub sample_index: f64,
    pub t_des: f64,
    pub step: f64,
    pub sps: f64,
    pub kp: f64,
    pub ki: f64,
    pub integrator: f64,
    pub is_on_time: bool,
    pub on_time_prev: Complex<f32>,
    pub mid_time: Complex<f32>,
    pub sample_count: usize,
}

impl GardnerLoop {
    pub fn new(sample_rate: f64, symbol_rate: f64) -> Self {
        let sps = sample_rate / symbol_rate;
        let step = sps / 2.0;
        Self {
            farrow: FarrowInterpolator::new(),
            sample_index: 0.0,
            t_des: 2.0,
            step,
            sps,
            kp: 0.01,
            ki: 0.001,
            integrator: 0.0,
            is_on_time: true,
            on_time_prev: Complex::new(0.0, 0.0),
            mid_time: Complex::new(0.0, 0.0),
            sample_count: 0,
        }
    }

    pub fn reset(&mut self) {
        self.farrow.reset();
        self.sample_index = 0.0;
        self.t_des = 2.0;
        self.integrator = 0.0;
        self.is_on_time = true;
        self.on_time_prev = Complex::new(0.0, 0.0);
        self.mid_time = Complex::new(0.0, 0.0);
        self.sample_count = 0;
        self.step = self.sps / 2.0;
    }

    pub fn process(&mut self, sample: Complex<f32>, output_symbols: &mut Vec<(Complex<f32>, f32)>) {
        self.farrow.push(sample);
        self.sample_index += 1.0;

        if self.sample_index < 4.0 {
            return;
        }

        if self.t_des < self.sample_index - 2.0 {
            self.t_des = self.sample_index - 2.0;
        }

        while self.t_des < self.sample_index - 1.0 {
            let mu = self.t_des - (self.sample_index - 2.0);
            if !(0.0..1.0).contains(&mu) {
                break;
            }

            let interp = self.farrow.interpolate(mu as f32);

            if self.is_on_time {
                let on_time_curr = interp;

                if self.sample_count >= 2 {
                    let error = (on_time_curr.re - self.on_time_prev.re) * self.mid_time.re
                        + (on_time_curr.im - self.on_time_prev.im) * self.mid_time.im;

                    self.integrator += error as f64 * self.ki;
                    let control = error as f64 * self.kp + self.integrator;

                    let nominal_step = self.sps / 2.0;
                    self.step =
                        (nominal_step - control).clamp(0.5 * nominal_step, 1.5 * nominal_step);
                }

                self.on_time_prev = on_time_curr;
                output_symbols.push((on_time_curr, mu as f32));
                self.is_on_time = false;
            } else {
                self.mid_time = interp;
                self.is_on_time = true;
            }

            self.sample_count += 1;
            self.t_des += self.step;
        }
    }
}

pub fn design_lowpass_filter(cutoff_hz: f64, sample_rate_hz: f64, num_taps: usize) -> Vec<f32> {
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

        // Hamming window
        let win =
            0.54 - 0.46 * (2.0 * std::f64::consts::PI * n as f64 / (num_taps - 1) as f64).cos();
        taps[n] = (val * win) as f32;
        sum += taps[n] as f64;
    }

    // Normalize taps so the DC gain is 1.0 (0 dB)
    for val in taps.iter_mut() {
        *val /= sum as f32;
    }
    taps
}

#[derive(Clone)]
pub struct FirDecimator {
    pub taps: Vec<f32>,
    pub taps_simd: Vec<f32>,
    pub decimation_factor: usize,
    pub history: Vec<Complex<f32>>,
    pub pending_offset: usize,
}

impl FirDecimator {
    pub fn new(taps: Vec<f32>, decimation_factor: usize) -> Self {
        let hist_len = taps.len().saturating_sub(1);
        let mut taps_simd = Vec::with_capacity(taps.len() * 2);
        for &t in &taps {
            taps_simd.push(t);
            taps_simd.push(t);
        }
        Self {
            taps,
            taps_simd,
            decimation_factor,
            history: vec![Complex::new(0.0, 0.0); hist_len],
            pending_offset: 0,
        }
    }

    pub fn reset(&mut self) {
        for s in self.history.iter_mut() {
            *s = Complex::new(0.0, 0.0);
        }
        self.pending_offset = 0;
    }

    #[inline(always)]
    pub fn compute(&self, window: &[Complex<f32>]) -> Complex<f32> {
        let num_taps = self.taps.len();
        assert!(window.len() >= num_taps);

        #[cfg(target_arch = "x86_64")]
        if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
            return unsafe { self.compute_x86_64(window) };
        }

        #[cfg(target_arch = "aarch64")]
        if std::arch::is_aarch64_feature_detected!("neon") {
            return unsafe { self.compute_aarch64(window) };
        }

        self.compute_scalar(window)
    }

    #[inline(always)]
    pub fn compute_scalar(&self, window: &[Complex<f32>]) -> Complex<f32> {
        let num_taps = self.taps.len();
        assert!(window.len() >= num_taps);
        let mut re = 0.0;
        let mut im = 0.0;
        let window_ptr = window.as_ptr() as *const f32;
        let taps_ptr = self.taps_simd.as_ptr();
        let len = self.taps_simd.len();
        
        let mut i = 0;
        while i < len {
            unsafe {
                re += *window_ptr.add(i) * *taps_ptr.add(i);
                im += *window_ptr.add(i + 1) * *taps_ptr.add(i + 1);
            }
            i += 2;
        }
        Complex::new(re, im)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,fma")]
    unsafe fn compute_x86_64(&self, window: &[Complex<f32>]) -> Complex<f32> {
        let num_taps = self.taps.len();
        assert!(window.len() >= num_taps);
        use std::arch::x86_64::*;
        let window_ptr = window.as_ptr() as *const f32;
        let taps_ptr = self.taps_simd.as_ptr();
        let len = self.taps_simd.len();
        
        let mut sum = unsafe { _mm256_setzero_ps() };
        let mut i = 0;
        while i + 8 <= len {
            unsafe {
                let w = _mm256_loadu_ps(window_ptr.add(i));
                let t = _mm256_loadu_ps(taps_ptr.add(i));
                sum = _mm256_fmadd_ps(w, t, sum);
            }
            i += 8;
        }
        
        let mut sum_arr = [0.0; 8];
        unsafe { _mm256_storeu_ps(sum_arr.as_mut_ptr(), sum) };
        
        let mut re = sum_arr[0] + sum_arr[2] + sum_arr[4] + sum_arr[6];
        let mut im = sum_arr[1] + sum_arr[3] + sum_arr[5] + sum_arr[7];
        
        while i < len {
            unsafe {
                re += *window_ptr.add(i) * *taps_ptr.add(i);
                im += *window_ptr.add(i + 1) * *taps_ptr.add(i + 1);
            }
            i += 2;
        }
        
        Complex::new(re, im)
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    unsafe fn compute_aarch64(&self, window: &[Complex<f32>]) -> Complex<f32> {
        let num_taps = self.taps.len();
        assert!(window.len() >= num_taps);
        use std::arch::aarch64::*;
        let window_ptr = window.as_ptr() as *const f32;
        let taps_ptr = self.taps_simd.as_ptr();
        let len = self.taps_simd.len();
        
        let mut sum = unsafe { vdupq_n_f32(0.0) };
        let mut i = 0;
        while i + 4 <= len {
            unsafe {
                let w = vld1q_f32(window_ptr.add(i));
                let t = vld1q_f32(taps_ptr.add(i));
                sum = vfmaq_f32(sum, w, t);
            }
            i += 4;
        }
        
        let mut sum_arr = [0.0; 4];
        unsafe { vst1q_f32(sum_arr.as_mut_ptr(), sum) };
        
        let mut re = sum_arr[0] + sum_arr[2];
        let mut im = sum_arr[1] + sum_arr[3];
        
        while i < len {
            unsafe {
                re += *window_ptr.add(i) * *taps_ptr.add(i);
                im += *window_ptr.add(i + 1) * *taps_ptr.add(i + 1);
            }
            i += 2;
        }
        
        Complex::new(re, im)
    }

    pub fn process(&mut self, input: &[Complex<f32>], output: &mut Vec<Complex<f32>>) {
        if self.decimation_factor <= 1 {
            output.reserve(input.len());
            output.extend_from_slice(input);
            return;
        }

        let num_taps = self.taps.len();
        let hist_len = self.history.len();
        let total_len = hist_len + input.len();

        let projected_capacity = (input.len() / self.decimation_factor) + 2;
        output.reserve(projected_capacity);

        let mut idx = self.pending_offset;
        // Boundary Phase: process samples that overlap with history
        while idx < hist_len && idx + num_taps <= total_len {
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

        // Main Phase: contiguous slice processing (using SIMD)
        while idx + num_taps <= total_len {
            let input_offset = idx - hist_len;
            let input_slice = &input[input_offset..input_offset + num_taps];
            output.push(self.compute(input_slice));
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

fn complex_cholesky_6x6(a: &[[Complex<f32>; 6]; 6]) -> Option<[[Complex<f32>; 6]; 6]> {
    let mut l = [[Complex::new(0.0, 0.0); 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k].conj();
            }
            if i == j {
                if sum.re <= 0.0 {
                    return None;
                }
                l[i][j] = Complex::new(sum.re.sqrt(), 0.0);
            } else {
                l[i][j] = sum / l[j][j].re;
            }
        }
    }
    Some(l)
}

fn complex_cholesky_solve_6(a: &[[Complex<f32>; 6]; 6], b: &[Complex<f32>; 6]) -> Option<[Complex<f32>; 6]> {
    let l = complex_cholesky_6x6(a)?;
    // Forward substitution L * y = b
    let mut y = [Complex::new(0.0, 0.0); 6];
    for i in 0..6 {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i].re;
    }
    // Backward substitution L^H * x = y
    let mut x = [Complex::new(0.0, 0.0); 6];
    for i in (0..6).rev() {
        let mut sum = y[i];
        for k in i+1..6 {
            sum -= l[k][i].conj() * x[k];
        }
        x[i] = sum / l[i][i].re;
    }
    Some(x)
}

/// Extensive Cancellation Algorithm (ECA) filter
/// Uses exact Orthogonal Subspace Projection on blocks of samples to obliterate 
/// the direct path and stationary clutter down to the noise floor.
#[derive(Clone, Debug)]
pub struct EcaCanceler {
    history: [Complex<f32>; 6],
}

impl EcaCanceler {
    pub fn new() -> Self {
        Self {
            history: [Complex::new(0.0, 0.0); 6],
        }
    }

    pub fn process_block(
        &mut self,
        input: &[Complex<f32>],
        output: &mut [Complex<f32>],
    ) {
        let n = input.len();
        if n == 0 { return; }
        
        let taps = 6;
        let mut r = [[Complex::new(0.0, 0.0); 6]; 6];
        let mut p = [Complex::new(0.0, 0.0); 6];

        // Construct extended signal x_ext = [history, input]
        let mut x_ext = vec![Complex::new(0.0, 0.0); n + 6];
        x_ext[0..6].copy_from_slice(&self.history);
        x_ext[6..n+6].copy_from_slice(input);

        // 1. Compute r[0][d] for d in 0..5 using O(N) operations
        for d in 0..taps {
            let mut sum = Complex::new(0.0, 0.0);
            for i in 0..n {
                sum += x_ext[5 + i].conj() * x_ext[5 + i - d];
            }
            r[0][d] = sum;
        }

        // 2. Compute the rest of the upper triangle of r using O(1) sliding window updates
        for d in 0..taps {
            for j in 1..(taps - d) {
                let term_in = x_ext[5 - j].conj() * x_ext[5 - j - d];
                let term_out = x_ext[5 - j + n].conj() * x_ext[5 - j + n - d];
                r[j][j + d] = r[j - 1][j - 1 + d] + term_in - term_out;
            }
        }

        // 3. Fill in the lower triangle of r using Hermitian symmetry
        for j in 0..taps {
            for k in 0..j {
                r[j][k] = r[k][j].conj();
            }
        }

        // 4. Compute p[j] for j in 0..5 using O(N) operations
        for j in 0..taps {
            let mut sum = Complex::new(0.0, 0.0);
            for i in 0..n {
                sum += x_ext[5 + i - j].conj() * input[i];
            }
            p[j] = sum;
        }

        // Diagonal regularization (Ridge / Tikhonov)
        let tau = 1e-3 * n as f32; 
        for j in 0..taps {
            r[j][j] += Complex::new(tau, 0.0);
        }

        let weights = match complex_cholesky_solve_6(&r, &p) {
            Some(w) => w,
            None => [Complex::new(0.0, 0.0); 6],
        };

        // 5. Generate outputs
        for i in 0..n {
            let mut y_clutter = Complex::new(0.0, 0.0);
            for k in 0..taps {
                y_clutter += x_ext[5 + i - k] * weights[k];
            }
            output[i] = input[i] - y_clutter;
        }

        // Update history
        for i in 0..6 {
            self.history[i] = if n > 5 - i {
                input[n - 1 - (5 - i)]
            } else {
                self.history[self.history.len() - (5 - i) + n]
            };
        }
    }
}

static ADIC_WEIGHTS: std::sync::OnceLock<[f64; 32]> = std::sync::OnceLock::new();

#[inline(always)]
fn get_adic_weight(v_2: u32) -> f64 {
    let weights = ADIC_WEIGHTS.get_or_init(|| {
        let mut arr = [0.0f64; 32];
        for i in 0..32 {
            arr[i] = 2.0_f64.powf(1.5 * i as f64);
        }
        arr
    });
    if (v_2 as usize) < weights.len() {
        weights[v_2 as usize]
    } else {
        2.0_f64.powf(1.5 * v_2 as f64)
    }
}

pub struct EnvelopeWaveletSpurCanceller;

impl EnvelopeWaveletSpurCanceller {
    pub fn notch_spurs_wavelet(
        fft_magnitudes: &mut [f32],
        estimated_doppler: f32,
        estimated_chirp: f32,
    ) {
        let n = fft_magnitudes.len();
        if n == 0 || (n & (n - 1)) != 0 {
            return;
        }

        // 1. Calculate dynamic sample rate based on test/daemon modes
        let sample_rate = if n <= 1024 { n as f64 } else { 50000.0 };

        // 2. Compute noise floor envelope background using min-plus Haar Wavelet
        let background = Self::compute_noise_envelope_background(fft_magnitudes);

        // 3. Detect stationary spurs
        let spikes = Self::detect_stationary_spurs(fft_magnitudes);

        // 4. Map carrier doppler and chirp to bins
        let bin_carrier = (estimated_doppler as f64 / sample_rate * n as f64).round() as isize;
        let bin_carrier = bin_carrier.rem_euclid(n as isize) as usize;

        let sweep_width = (estimated_chirp.abs() as f64 / sample_rate * n as f64).ceil() as usize;
        let guard_width = 3 + sweep_width.min(10);

        // 5. Apply notch to non-guarded spikes
        for s in spikes {
            // Check circular distance to carrier
            let diff = (s as isize - bin_carrier as isize).rem_euclid(n as isize);
            let dist = diff.min(n as isize - diff) as usize;

            if dist > guard_width {
                // Replace with background level
                fft_magnitudes[s] = background[s];
            }
        }
    }

    pub fn detect_stationary_spurs(fft_magnitudes: &[f32]) -> Vec<usize> {
        let n = fft_magnitudes.len();
        if n == 0 || (n & (n - 1)) != 0 {
            return vec![];
        }

        // Compute local noise background
        let background = Self::compute_noise_envelope_background(fft_magnitudes);

        // Compute median magnitude for adaptive thresholding
        let mut sorted = fft_magnitudes.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = if sorted[n / 2].is_finite() {
            sorted[n / 2].max(1e-6)
        } else {
            1e-6
        };

        let threshold_wavelet = 3.0 * median;
        let threshold_vladimirov = 10.0 * median;

        // Precompute weights w_d and their prefix sums P_k using cached 2-adic weights
        let mut w = vec![0.0f64; n];
        for d in 1..n {
            let v_2 = d.trailing_zeros();
            w[d] = get_adic_weight(v_2);
        }
        let mut p = vec![0.0f64; n];
        for k in 1..n {
            p[k] = p[k - 1] + w[k];
        }

        // Prepare zero-padded input complex signal X of length 2N
        let mut x_complex = vec![Complex::new(0.0f64, 0.0f64); 2 * n];
        for i in 0..n {
            let x_val = if fft_magnitudes[i].is_finite() {
                fft_magnitudes[i] as f64
            } else {
                0.0
            };
            x_complex[i] = Complex::new(x_val, 0.0);
        }

        // Prepare symmetric complex kernel W of length 2N
        let mut w_complex = vec![Complex::new(0.0f64, 0.0f64); 2 * n];
        w_complex[0] = Complex::new(0.0, 0.0);
        for d in 1..n {
            w_complex[d] = Complex::new(w[d], 0.0);
            w_complex[2 * n - d] = Complex::new(w[d], 0.0);
        }
        let w_n = get_adic_weight(n.trailing_zeros());
        w_complex[n] = Complex::new(w_n, 0.0);

        // Perform forward FFTs using rustfft
        let mut planner = FftPlanner::<f64>::new();
        let fft_forward = planner.plan_fft_forward(2 * n);
        fft_forward.process(&mut x_complex);
        fft_forward.process(&mut w_complex);

        // Multiply element-wise
        let mut prod = vec![Complex::new(0.0f64, 0.0f64); 2 * n];
        for i in 0..(2 * n) {
            prod[i] = x_complex[i] * w_complex[i];
        }

        // Perform inverse FFT
        let fft_inverse = planner.plan_fft_inverse(2 * n);
        fft_inverse.process(&mut prod);

        // Normalize
        let norm = 1.0 / (2 * n) as f64;
        for i in 0..(2 * n) {
            prod[i] *= norm;
        }

        // Detect spikes
        let mut spikes = Vec::new();
        for i in 0..n {
            let x_i = if fft_magnitudes[i].is_finite() {
                fft_magnitudes[i] as f64
            } else {
                0.0
            };

            // Only examine bins that exceed wavelet background threshold
            if x_i as f32 <= background[i] + threshold_wavelet {
                continue;
            }

            let s_i = p[i] + p[n - 1 - i];
            let c_i = prod[i].re;
            let sum_deriv = x_i * s_i - c_i;

            if sum_deriv.is_finite() && sum_deriv as f32 > threshold_vladimirov {
                spikes.push(i);
            }
        }
        spikes
    }

    pub fn check_pca_epoch_guard(t_obs: DateTime<Utc>, pca_time: DateTime<Utc>) -> bool {
        let diff = (t_obs - pca_time).num_seconds().abs();
        diff <= 30
    }

    /// Helper to compute lower background envelope using min-plus Haar Wavelet (Tropical addition)
    fn compute_noise_envelope_background(magnitudes: &[f32]) -> Vec<f32> {
        let n = magnitudes.len();
        let mut approx = magnitudes.to_vec();
        let mut current_len = n;
        let mut levels_done = 0;

        // Perform 3-level decomposition
        for _scale in 0..3 {
            let next_len = current_len / 2;
            if next_len == 0 {
                break;
            }
            let mut next_approx = vec![0.0; next_len];
            for k in 0..next_len {
                let v1 = if approx[2 * k].is_finite() {
                    approx[2 * k]
                } else {
                    0.0
                };
                let v2 = if approx[2 * k + 1].is_finite() {
                    approx[2 * k + 1]
                } else {
                    0.0
                };
                next_approx[k] = v1.min(v2); // Min-plus addition
            }
            approx = next_approx;
            current_len = next_len;
            levels_done += 1;
        }

        // Reconstruct background by zeroing out detail coefficients
        let mut background = approx;
        for _scale in 0..levels_done {
            let mut next_bg = vec![0.0; background.len() * 2];
            for k in 0..background.len() {
                next_bg[2 * k] = background[k];
                next_bg[2 * k + 1] = background[k];
            }
            background = next_bg;
        }
        background
    }
}

pub struct AbsolutePowerGainController;

impl AbsolutePowerGainController {
    pub fn update_gain(
        samples: &[Complex<f32>],
        current_lna: &mut f32,
        current_vga: &mut f32,
        current_amp: &mut f32,
    ) {
        // 1. Clamp gains to valid hardware ranges
        *current_lna = current_lna.clamp(0.0, 40.0);
        *current_vga = current_vga.clamp(0.0, 62.0);
        *current_amp = if *current_amp >= 7.0 { 14.0 } else { 0.0 };

        if samples.is_empty() {
            return;
        }

        // 2. Compute RMS and Clipping Ratio
        let mut sum_mag_sq = 0.0f64;
        let mut clip_count = 0usize;
        for s in samples {
            let mag_sq = (s.re * s.re + s.im * s.im) as f64;
            if mag_sq >= 0.98 * 0.98 {
                clip_count += 1;
            }
            sum_mag_sq += mag_sq;
        }
        let rms = (sum_mag_sq / samples.len() as f64).sqrt();
        let clipping_ratio = clip_count as f64 / samples.len() as f64;

        // 3. Adjust gains based on optimized hysteresis thresholds to maximize SNR
        if clipping_ratio > 0.005 || rms > 0.5 {
            // Saturation: lower gains prioritizing the final stages to preserve low noise figure at LNA
            if *current_amp > 0.0 {
                *current_amp = 0.0;
            } else if *current_vga > 0.0 {
                *current_vga = (*current_vga - 2.0).max(0.0);
            } else if *current_lna > 0.0 {
                *current_lna = (*current_lna - 8.0).max(0.0);
            }
        } else if rms < 0.05 {
            // Weak signal: raise gains prioritizing early stages to minimize total noise figure
            if *current_lna < 40.0 {
                *current_lna = (*current_lna + 8.0).min(40.0);
            } else if *current_vga < 62.0 {
                *current_vga = (*current_vga + 2.0).min(62.0);
            } else if *current_amp < 14.0 {
                *current_amp = 14.0;
            }
        }
    }

    pub fn estimate_absolute_power(samples: &[Complex<f32>], lna: f32, vga: f32, amp: f32) -> f32 {
        if samples.is_empty() {
            return -120.0;
        }
        let mut sum_mag_sq = 0.0f64;
        for s in samples {
            sum_mag_sq += (s.re * s.re + s.im * s.im) as f64;
        }
        let mean_power = sum_mag_sq / samples.len() as f64;
        let p_dig = 10.0 * (mean_power + 1e-12).log10() as f32;
        let nominal_gain = lna + vga + amp;
        let cal_correction = Self::get_gain_calibration_offset(nominal_gain);
        p_dig - nominal_gain + cal_correction
    }

    pub fn get_gain_calibration_offset(gain: f32) -> f32 {
        // Polynomial gain offset lookup calibration curve (fits MAX2837 LNA/VGA non-linearities)
        0.005 * (gain - 24.0) * (gain - 24.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum ChannelStatus {
    Idle,
    Acquisition,
    Locked,
    Fade,
}

pub enum ChannelCommand {
    Allocate {
        channel_index: usize,
        sat_name: String,
        elements: sgp4::Elements,
        target_freq: f64,
        initial_freq: f64,
    },
    Deallocate {
        channel_index: usize,
    },
    UpdateTargetFrequency {
        channel_index: usize,
        target_freq: f64,
    },
}

pub struct DigitalDownConverter {
    pub phase_acc: f64,
}

impl DigitalDownConverter {
    pub fn new() -> Self {
        Self { phase_acc: 0.0 }
    }

    pub fn process(
        &mut self,
        input: &[Complex<f32>],
        f_shift: f64,
        sample_rate: f64,
        output: &mut [Complex<f32>],
    ) {
        if f_shift.abs() < 1e-6 {
            output.copy_from_slice(input);
            return;
        }

        let phase_step = -2.0 * std::f64::consts::PI * f_shift / sample_rate;
        // Use f64 phasor accumulation to prevent magnitude drift from repeated f32 multiplication
        let mut phasor_re = self.phase_acc.cos();
        let mut phasor_im = self.phase_acc.sin();
        let step_re = phase_step.cos();
        let step_im = phase_step.sin();

        for (n, &sample) in input.iter().enumerate() {
            output[n] = Complex::new(
                (sample.re as f64 * phasor_re - sample.im as f64 * phasor_im) as f32,
                (sample.re as f64 * phasor_im + sample.im as f64 * phasor_re) as f32,
            );
            let new_re = phasor_re * step_re - phasor_im * step_im;
            let new_im = phasor_re * step_im + phasor_im * step_re;
            phasor_re = new_re;
            phasor_im = new_im;

            // Renormalize every 256 samples to bound accumulated f64 drift
            if n % 256 == 0 {
                let mag = (phasor_re * phasor_re + phasor_im * phasor_im).sqrt();
                if mag > 1e-15 {
                    phasor_re /= mag;
                    phasor_im /= mag;
                }
            }
        }

        self.phase_acc = phasor_im.atan2(phasor_re);
    }
}

#[derive(Debug, Clone)]
pub struct TelemetryUpdate {
    pub snr: f64,
    pub is_locked: bool,
    pub frequency: f64,
    pub phase: f64,
    pub sat_name: String,
}

pub struct DemodChannel {
    pub id: usize,
    pub sat_name: String,
    pub elements: Option<sgp4::Elements>,
    pub constants: Option<sgp4::Constants>,
    pub status: ChannelStatus,
    pub target_freq: f64,
    pub initial_freq: f64,

    // DSP components
    pub ddc: DigitalDownConverter,
    pub decimator: FirDecimator,
    pub pll_tracker: CarrierPllEkf,
    pub gardner_loop: GardnerLoop,
    pub tracking_bank: Option<EkfTrackingBank>,

    // Buffers and tracking state
    pub decimated_queue: VecDeque<Complex<f32>>,
    pub pass_samples: Vec<(DateTime<Utc>, f64)>,
    pub first_lock_time: Option<DateTime<Utc>>,
    pub last_lock_time: Option<DateTime<Utc>>,
    pub fade_start_time: Option<DateTime<Utc>>,
    pub ever_locked: bool,
    pub last_snr: f32,

    // Coarse acquisition sweep state
    pub last_alpha: f64,
    pub consecutive_lock: usize,
    pub consecutive_unlock: usize,
    pub unlocked_frames_during_pass: usize,

    // EKF tuning config
    pub modulation: Modulation,
    pub min_snr: f32,
    pub fade_timeout: f64,

    // Compatibility & telemetry fields
    pub nominal_freq: f64,
    pub sample_rate: f64,
    pub symbol_rate: f64,
    pub is_locked: bool,
    pub symbol_locked: bool,
    pub snr: f64,
    pub frequency: f64,
    pub phase: f64,
    pub telemetry_sender: Option<crossbeam_channel::Sender<TelemetryUpdate>>,
    pub mixed_samples: Vec<Complex<f32>>,
    pub decimated_samples: Vec<Complex<f32>>,
    pub normalized_iq: Vec<Complex<f32>>,
    pub last_processed_len: usize,
    pub sample_count: usize,
    pub eca_enabled: bool,
    pub eca_canceler: EcaCanceler,
}

impl DemodChannel {
    pub fn new_prod(
        id: usize,
        pipeline_sample_rate: f64,
        sym_rate: f64,
        modulation: Modulation,
        no_adaptive_ekf: bool,
        no_dual_lock: bool,
        no_gardner: bool,
        no_multihypothesis: bool,
        min_snr: f32,
        fade_timeout: f64,
        taps: Vec<f32>,
        decimate_factor: usize,
        eca_enabled: bool,
    ) -> Self {
        let decimated_rate = pipeline_sample_rate / decimate_factor as f64;
        let mut pll_tracker = CarrierPllEkf::new(decimated_rate, modulation);
        pll_tracker.adaptive_ekf = !no_adaptive_ekf;
        pll_tracker.dual_lock = !no_dual_lock;
        if !no_gardner {
            pll_tracker.ts = 1.0 / sym_rate;
        }

        let tracking_bank = if !no_multihypothesis {
            let mut bank = EkfTrackingBank {
                trackers: [
                    CarrierPllEkf::new(decimated_rate, modulation),
                    CarrierPllEkf::new(decimated_rate, modulation),
                    CarrierPllEkf::new(decimated_rate, modulation),
                ],
                gardner_loops: [
                    GardnerLoop::new(decimated_rate, sym_rate),
                    GardnerLoop::new(decimated_rate, sym_rate),
                    GardnerLoop::new(decimated_rate, sym_rate),
                ],
                active_idx: None,
                in_fade: false,
                fade_counter: 0,
                max_fade_steps: (fade_timeout * decimated_rate / 40000.0) as usize,
                terminated_in_fade: false,
                spur_dwell_counter: 0,
            };
            for t in &mut bank.trackers {
                t.adaptive_ekf = !no_adaptive_ekf;
                t.dual_lock = !no_dual_lock;
                if !no_gardner {
                    t.ts = 1.0 / sym_rate;
                }
            }
            Some(bank)
        } else {
            None
        };

        Self {
            id,
            sat_name: String::new(),
            elements: None,
            constants: None,
            status: ChannelStatus::Idle,
            target_freq: 0.0,
            initial_freq: 0.0,

            ddc: DigitalDownConverter::new(),
            decimator: FirDecimator::new(taps, decimate_factor),
            pll_tracker,
            gardner_loop: GardnerLoop::new(decimated_rate, sym_rate),
            tracking_bank,

            decimated_queue: VecDeque::new(),
            pass_samples: Vec::new(),
            first_lock_time: None,
            last_lock_time: None,
            fade_start_time: None,
            ever_locked: false,
            last_snr: 0.0,

            last_alpha: 0.0,
            consecutive_lock: 0,
            consecutive_unlock: 0,
            unlocked_frames_during_pass: 0,

            modulation,
            min_snr,
            fade_timeout,

            nominal_freq: 0.0,
            sample_rate: pipeline_sample_rate,
            symbol_rate: sym_rate,
            is_locked: false,
            symbol_locked: false,
            snr: 0.0,
            frequency: 0.0,
            phase: 0.0,
            telemetry_sender: None,
            mixed_samples: Vec::new(),
            decimated_samples: Vec::new(),
            normalized_iq: Vec::new(),
            last_processed_len: 0,
            sample_count: 0,
            eca_enabled,
            eca_canceler: EcaCanceler::new(),
        }
    }

    pub fn new(nominal_freq: f64, sample_rate: f64, sat_name: String) -> Self {
        let decimate_factor = 4;
        let sym_rate = 10000.0;
        let taps = design_lowpass_filter(0.4 * sample_rate, sample_rate, 127);
        let mut ch = Self::new_prod(
            0,
            sample_rate,
            sym_rate,
            Modulation::Carrier,
            false,
            false,
            false,
            true, // no_multihypothesis = true
            3.0,
            5.0,
            taps,
            decimate_factor,
            false,
        );
        ch.nominal_freq = nominal_freq;
        ch.sample_rate = sample_rate;
        ch.sat_name = sat_name;
        ch.frequency = nominal_freq;
        ch.status = ChannelStatus::Acquisition; // Test channels start active
        ch
    }

    pub fn process_block(&mut self, raw_iq: &[Complex<f32>]) {
        self.process_block_with_center(raw_iq, 0.0);
    }

    pub fn process_block_with_center(&mut self, raw_iq: &[Complex<f32>], center_freq: f64) {
        // Skip all processing for idle channels — no DDC, no decimation, no EKF
        if self.status == ChannelStatus::Idle {
            self.is_locked = false;
            self.symbol_locked = false;
            self.snr = 0.0;
            self.last_processed_len = raw_iq.len();
            return;
        }

        if raw_iq.is_empty() {
            self.mixed_samples.clear();
            self.decimated_samples.clear();
            self.is_locked = false;
            self.symbol_locked = false;
            self.snr = -5.0;
            self.last_processed_len = 0;
            return;
        }

        let mut has_signal = false;
        let mut has_nan = false;
        let mut has_clipping = false;

        for &sample in raw_iq {
            if sample.re.is_nan() || sample.im.is_nan() {
                has_nan = true;
            }
            let mag_sq = sample.re * sample.re + sample.im * sample.im;
            if mag_sq > 1e-12 {
                has_signal = true;
            }
            if mag_sq > 1.98 {
                has_clipping = true;
            }
        }

        if has_nan {
            self.pll_tracker.is_locked = false;
            self.is_locked = false;
            self.symbol_locked = false;
            self.snr = 0.0;
            self.last_processed_len = raw_iq.len();
            self.sample_count += raw_iq.len();
            return;
        }

        // Guard clipping BEFORE DDC to prevent advancing phase_acc on discarded blocks
        if has_clipping {
            self.pll_tracker.is_locked = false;
            self.is_locked = false;
            self.symbol_locked = false;
            self.snr = 0.5;
            self.last_processed_len = raw_iq.len();
            self.sample_count += raw_iq.len();
            return;
        }

        let freq_to_use = if self.nominal_freq == 0.0 {
            self.target_freq
        } else {
            self.nominal_freq
        };
        let f_shift = freq_to_use - center_freq;
        self.mixed_samples
            .resize(raw_iq.len(), Complex::new(0.0, 0.0));

        // a. Resize self.normalized_iq.
        self.normalized_iq.resize(raw_iq.len(), Complex::new(0.0, 0.0));

        // b. Perform Subspace Projection / ECA for LO Leakage and Clutter Cancellation.
        if self.eca_enabled {
            self.eca_canceler.process_block(raw_iq, &mut self.normalized_iq);
        } else {
            self.normalized_iq.copy_from_slice(raw_iq);
            let n = self.normalized_iq.len() as f64;
            let mut sum = Complex::<f64>::new(0.0, 0.0);
            for &s in &self.normalized_iq {
                sum += Complex::new(s.re as f64, s.im as f64);
            }
            let mean = Complex::new((sum.re / n) as f32, (sum.im / n) as f32);
            for s in &mut self.normalized_iq {
                *s -= mean;
            }
        }

        // c. Invoke self.ddc.process.
        self.ddc
            .process(&self.normalized_iq, f_shift, self.sample_rate, &mut self.mixed_samples);

        self.decimated_samples.clear();
        self.decimator
            .process(&self.mixed_samples, &mut self.decimated_samples);

        // Dynamic SNR estimation via M2M4 moments
        let mut m2 = 0.0;
        let mut m4 = 0.0;
        for &s in &self.decimated_samples {
            let mag_sq = (s.re * s.re + s.im * s.im) as f64;
            m2 += mag_sq;
            m4 += mag_sq * mag_sq;
        }
        let len = self.decimated_samples.len() as f64;
        let computed_snr = if len > 0.0 {
            m2 /= len;
            m4 /= len;
            let v = 2.0 * m2 * m2 - m4;
            let ps = if v > 0.0 { v.sqrt() } else { 0.0 };
            let pn = m2 - ps;
            if pn > 1e-10 && ps > 1e-10 {
                (10.0 * (ps / pn).log10()).clamp(-20.0, 40.0)
            } else if ps > 1e-10 {
                40.0
            } else {
                -20.0
            }
        } else {
            -20.0
        };


        let decimated_ts = self.decimator.decimation_factor as f64 / self.sample_rate;

        if has_signal {
            let scale = if self.pll_tracker.modulation == Modulation::Bpsk { 2.0 } else { 1.0 };
            // EKF tracker initialization transitions
            let was_no_signal = self.last_snr < 0.0;
            if !self.pll_tracker.is_locked && (!self.ever_locked || was_no_signal) {
                let esprit_rate = self.sample_rate / self.decimator.decimation_factor as f64;
                let chunk_size = (0.01 * esprit_rate).round() as usize;
                let buffer_slice = if self.decimated_samples.len() >= chunk_size {
                    &self.decimated_samples[..chunk_size]
                } else {
                    &self.decimated_samples
                };

                let init_freq_offset = if !buffer_slice.is_empty() {
                    if self.pll_tracker.modulation == Modulation::Bpsk {
                        // BPSK squaring doubles the carrier frequency. The ESPRIT estimate will be
                        // at 2× the true carrier, which matches the EKF's internal tracking domain
                        // (the EKF also squares samples in update()). Halving happens at readout
                        // via the `scale` divisor — see dsp.rs line ~1073 and ekf.rs line ~357.
                        let squared_slice: Vec<Complex<f32>> = buffer_slice.iter().map(|&s| s * s).collect();
                        estimate_frequency_esprit(&squared_slice, esprit_rate, 10)
                    } else {
                        estimate_frequency_esprit(buffer_slice, esprit_rate, 10)
                    }
                } else {
                    0.0
                };

                self.pll_tracker.reset(0.0, init_freq_offset, 0.0);
                // Bootstrap bank trackers in lock-step: without this, bank trackers
                // stay is_locked=false forever and their EKF update loop never runs.
                if let Some(ref mut bank) = self.tracking_bank {
                    for t in &mut bank.trackers {
                        if !t.is_locked {
                            t.reset(0.0, init_freq_offset, 0.0);
                        }
                    }
                }
            }

            // d. Perform Bussgang Normalization (Constant Modulus Projection) on narrowband decimated signal.
            // NOTE: This must happen AFTER SNR estimation (which needs amplitude variance)
            // and AFTER ESPRIT bootstrap (which needs spectral amplitude structure),
            // but BEFORE the EKF tracking loop (which benefits from constant-modulus input).
            for s in &mut self.decimated_samples {
                let norm = s.norm();
                if norm > 1e-6 {
                    *s = *s / norm;
                } else {
                    *s = Complex::new(0.0, 0.0);
                }
            }

            // Execute single or multi-hypothesis tracking updates
            if let Some(ref mut bank) = self.tracking_bank {
                let mut symbols = Vec::with_capacity(4); // Pre-allocate outside hot loop
                for &s in &self.decimated_samples {
                    for i in 0..3 {
                        let tracker = &mut bank.trackers[i];
                        let g_loop = &mut bank.gardner_loops[i];

                        if tracker.is_locked {
                            if (tracker.ts - decimated_ts).abs() > 1e-9 {
                                let theta = tracker.x[0] / scale;
                                let cos_theta = (-theta).cos();
                                let sin_theta = (-theta).sin();
                                let derotated_s = Complex::new(
                                    (s.re as f64 * cos_theta - s.im as f64 * sin_theta) as f32,
                                    (s.re as f64 * sin_theta + s.im as f64 * cos_theta) as f32,
                                );
                                symbols.clear();
                                g_loop.process(derotated_s, &mut symbols);

                                for (sym_derot, mu) in symbols.drain(..) {
                                    let dt_sample = decimated_ts;
                                    let true_theta =
                                        theta - (tracker.x[1] / scale) * (3.0 - mu as f64) * dt_sample;
                                    let cos_inv = true_theta.cos();
                                    let sin_inv = true_theta.sin();
                                    let sym_raw = Complex::new(
                                        (sym_derot.re as f64 * cos_inv
                                            - sym_derot.im as f64 * sin_inv)
                                            as f32,
                                        (sym_derot.re as f64 * sin_inv
                                            + sym_derot.im as f64 * cos_inv)
                                            as f32,
                                    );
                                    let prev_ts = tracker.ts;
                                    tracker.ts = 1.0 / self.symbol_rate; // Symbol rate tracking ts scale
                                    tracker.predict();
                                    tracker.update(sym_raw);
                                    tracker.ts = prev_ts;
                                }
                            } else {
                                let prev_ts = tracker.ts;
                                tracker.ts = decimated_ts;
                                tracker.predict();
                                tracker.update(s);
                                tracker.ts = prev_ts;
                            }
                        } else if bank.in_fade && Some(i) == bank.active_idx {
                            tracker.predict();
                        } else {
                            g_loop.reset();
                        }
                    }
                }

                // Check tracker discrepancy
                let obs = bank.compute_tracker_discrepancy();
                if obs > 150.0 {
                    bank.terminated_in_fade = true;
                    let f0 = (bank.trackers[0].x[1] / scale) / (2.0 * std::f64::consts::PI);
                    if bank.trackers[0].is_locked {
                        for i in 1..3 {
                            if bank.trackers[i].is_locked {
                                let fi = (bank.trackers[i].x[1] / scale) / (2.0 * std::f64::consts::PI);
                                let dev = (fi - f0).abs();
                                if dev > 150.0 {
                                    bank.trackers[i].is_locked = false;
                                    bank.trackers[i].lock_metric = 0.0;
                                }
                            }
                        }
                    } else if bank.trackers[1].is_locked && bank.trackers[2].is_locked {
                        let f1 = (bank.trackers[1].x[1] / scale) / (2.0 * std::f64::consts::PI);
                        let f2 = (bank.trackers[2].x[1] / scale) / (2.0 * std::f64::consts::PI);
                        if (f1 - f2).abs() > 150.0 {
                            let prune_idx =
                                if bank.trackers[1].lock_metric < bank.trackers[2].lock_metric {
                                    1
                                } else {
                                    2
                                };
                            bank.trackers[prune_idx].is_locked = false;
                            bank.trackers[prune_idx].lock_metric = 0.0;
                        }
                    }
                }

                // Select active tracker in bank
                let mut best_i = None;
                let mut best_score = -1.0;
                for i in 0..3 {
                    let tracker = &bank.trackers[i];
                    if tracker.is_locked && tracker.lock_metric > best_score {
                        best_score = tracker.lock_metric;
                        best_i = Some(i);
                    }
                }

                if let Some(selected_i) = best_i {
                    bank.active_idx = Some(selected_i);
                    let active_tracker = &bank.trackers[selected_i];
                    let ch_freq_offset = (active_tracker.x[1] / scale) / (2.0 * std::f64::consts::PI);
                    self.is_locked = true;
                    self.symbol_locked = true;
                    self.frequency = freq_to_use + ch_freq_offset;
                    self.phase = active_tracker.x[0] / scale;
                    bank.in_fade = false;
                    bank.fade_counter = 0;
                    self.ever_locked = true;
                } else {
                    if let Some(prev_i) = bank.active_idx {
                        bank.in_fade = true;
                        bank.fade_counter += 1;
                        if bank.fade_counter > bank.max_fade_steps {
                            bank.in_fade = false;
                            bank.active_idx = None;
                            bank.terminated_in_fade = true;
                            self.is_locked = false;
                            self.symbol_locked = false;
                            // Ungate bootstrap so re-acquisition can fire on next block
                            self.pll_tracker.is_locked = false;
                            self.ever_locked = false;
                        } else {
                            let tracker = &bank.trackers[prev_i];
                            let ch_freq_offset = (tracker.x[1] / scale) / (2.0 * std::f64::consts::PI);
                            self.is_locked = true;
                            self.symbol_locked = true;
                            self.frequency = freq_to_use + ch_freq_offset;
                            self.phase = tracker.x[0] / scale;
                        }
                    } else {
                        self.is_locked = false;
                        self.symbol_locked = false;
                        // No active tracker and no fade — ungate bootstrap for re-acquisition
                        self.pll_tracker.is_locked = false;
                        self.ever_locked = false;
                    }
                }
            } else {
                // Single EKF / Gardner loop execution
                let mut symbols = Vec::with_capacity(4); // Pre-allocate outside hot loop
                for &s in &self.decimated_samples {
                    if self.pll_tracker.is_locked {
                        if (self.pll_tracker.ts - decimated_ts).abs() > 1e-9 {
                            let theta = self.pll_tracker.x[0] / scale;
                            let cos_theta = (-theta).cos();
                            let sin_theta = (-theta).sin();
                            let derotated_s = Complex::new(
                                (s.re as f64 * cos_theta - s.im as f64 * sin_theta) as f32,
                                (s.re as f64 * sin_theta + s.im as f64 * cos_theta) as f32,
                            );
                            symbols.clear();
                            self.gardner_loop.process(derotated_s, &mut symbols);

                            for (sym_derot, mu) in symbols.drain(..) {
                                let dt_sample = decimated_ts;
                                let true_theta =
                                    theta - (self.pll_tracker.x[1] / scale) * (3.0 - mu as f64) * dt_sample;
                                let cos_inv = true_theta.cos();
                                let sin_inv = true_theta.sin();
                                let sym_raw = Complex::new(
                                    (sym_derot.re as f64 * cos_inv - sym_derot.im as f64 * sin_inv)
                                        as f32,
                                    (sym_derot.re as f64 * sin_inv + sym_derot.im as f64 * cos_inv)
                                        as f32,
                                );
                                let prev_ts = self.pll_tracker.ts;
                                self.pll_tracker.ts = 1.0 / self.symbol_rate;
                                self.pll_tracker.predict();
                                self.pll_tracker.update(sym_raw);
                                self.pll_tracker.ts = prev_ts;
                            }
                        } else {
                            let prev_ts = self.pll_tracker.ts;
                            self.pll_tracker.ts = decimated_ts;
                            self.pll_tracker.predict();
                            self.pll_tracker.update(s);
                            self.pll_tracker.ts = prev_ts;
                        }
                    } else {
                        self.gardner_loop.reset();
                    }
                }

                self.is_locked = self.pll_tracker.is_locked;
                self.symbol_locked = self.is_locked;
                self.frequency =
                    freq_to_use + ((self.pll_tracker.x[1] / scale) / (2.0 * std::f64::consts::PI));
                self.phase = self.pll_tracker.x[0] / scale;
                if self.is_locked {
                    self.ever_locked = true;
                }
            }

            self.snr = computed_snr;
        } else {
            self.pll_tracker.is_locked = false;
            if let Some(ref mut bank) = self.tracking_bank {
                for t in &mut bank.trackers {
                    t.is_locked = false;
                }
                bank.active_idx = None;
                bank.in_fade = false;
            }
            self.is_locked = false;
            self.symbol_locked = false;
            self.snr = computed_snr;
        }
        self.last_snr = self.snr as f32;

        self.last_processed_len = raw_iq.len();
        self.sample_count += raw_iq.len();

        if let Some(ref sender) = self.telemetry_sender {
            let _ = sender.send(TelemetryUpdate {
                snr: self.snr,
                is_locked: self.is_locked,
                frequency: self.frequency,
                phase: self.phase,
                sat_name: self.sat_name.clone(),
            });
        }
    }

    pub fn reset_to_idle(&mut self) {
        self.sat_name.clear();
        self.elements = None;
        self.constants = None;
        self.status = ChannelStatus::Idle;
        self.target_freq = 0.0;
        self.initial_freq = 0.0;
        self.ddc.phase_acc = 0.0;
        self.decimator.reset();
        self.gardner_loop.reset();
        self.decimated_queue.clear();
        self.pass_samples.clear();
        self.first_lock_time = None;
        self.last_lock_time = None;
        self.fade_start_time = None;
        self.ever_locked = false;
        self.last_snr = 0.0;
        self.last_alpha = 0.0;
        self.consecutive_lock = 0;
        self.consecutive_unlock = 0;
        self.unlocked_frames_during_pass = 0;

        self.pll_tracker.is_locked = false;
        if let Some(ref mut bank) = self.tracking_bank {
            for t in &mut bank.trackers {
                t.is_locked = false;
            }
            bank.active_idx = None;
            bank.in_fade = false;
            bank.fade_counter = 0;
            bank.terminated_in_fade = false;
            bank.spur_dwell_counter = 0;
        }

        // Reset the compatibility/testing fields
        self.is_locked = false;
        self.symbol_locked = false;
        self.snr = 0.0;
        self.frequency = 0.0;
        self.phase = 0.0;
        self.mixed_samples.clear();
        self.decimated_samples.clear();
        self.last_processed_len = 0;
        self.sample_count = 0;
    }
}

pub fn process_pipeline_parallel(
    channels: &mut [DemodChannel],
    raw_iq: &[Complex<f32>],
    center_freq: f64,
    _sample_rate: f64,
) {
    if channels.is_empty() {
        return;
    }

    use rayon::prelude::*;
    channels.par_iter_mut().for_each(|channel| {
        channel.process_block_with_center(raw_iq, center_freq);
    });
}

pub fn estimate_frequency_esprit(samples: &[Complex<f32>], sample_rate: f64, m: usize) -> f64 {
    let n = samples.len();
    if n < m {
        return 0.0;
    }
    let l = n - m + 1;
    
    let mut r_xx = nalgebra::DMatrix::<Complex<f64>>::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            let mut sum = Complex::new(0.0f64, 0.0f64);
            for k in 0..l {
                let s_i = Complex::new(samples[k + m - 1 - i].re as f64, samples[k + m - 1 - i].im as f64);
                let s_j = Complex::new(samples[k + m - 1 - j].re as f64, samples[k + m - 1 - j].im as f64);
                sum += s_i * s_j.conj();
            }
            r_xx[(i, j)] = sum / (l as f64);
        }
    }

    let mut u = nalgebra::DVector::<Complex<f64>>::from_element(m, Complex::new(1.0, 0.0));
    for _ in 0..15 {
        let w = &r_xx * &u;
        let norm = w.norm();
        if norm > 1e-12 {
            u = w.map(|val| val / norm);
        } else {
            break;
        }
    }

    let mut numerator = Complex::new(0.0f64, 0.0f64);
    let mut denominator = 0.0f64;
    for i in 0..(m - 1) {
        numerator += u[i].conj() * u[i + 1];
        denominator += u[i].norm_sqr();
    }
    
    if denominator > 1e-12 {
        let psi = numerator / denominator;
        let angle = psi.im.atan2(psi.re);
        let freq = -(angle * sample_rate) / (2.0 * std::f64::consts::PI);
        freq.clamp(-sample_rate / 2.0, sample_rate / 2.0)
    } else {
        0.0
    }
}


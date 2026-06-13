use nalgebra::{SMatrix, SVector};

/// Speed of light in m/s
pub const C: f64 = 299792458.0;

#[derive(Clone, Debug)]
pub struct MasterNavEkf {
    /// State vector: [X, Y, Z, Vx, Vy, Vz, dt_clk, df_clk]^T
    /// Coordinates are in ECEF meters and m/s.
    /// dt_clk is clock bias in seconds.
    /// df_clk is clock drift in seconds/second.
    pub x: SVector<f64, 8>,
    
    /// Error covariance matrix
    pub p: SMatrix<f64, 8, 8>,
    
    /// Process noise spectral densities
    pub q_acc: f64,       // acceleration random walk spectral density (m^2/s^3)
    pub q_clk_bias: f64,  // clock bias random walk spectral density (s^2/s)
    pub q_clk_drift: f64, // clock drift random walk spectral density (s^2/s^3)
    
    /// Measurement noise floor values
    pub r_phase_base: f64, // baseline carrier phase/range measurement variance (m^2)
    pub r_freq_base: f64,  // baseline Doppler/range rate measurement variance (m^2/s^2)
}

impl MasterNavEkf {
    pub fn new(init_pos: [f64; 3]) -> Self {
        let mut x = SVector::<f64, 8>::zeros();
        x[0] = init_pos[0];
        x[1] = init_pos[1];
        x[2] = init_pos[2];
        
        let mut p = SMatrix::<f64, 8, 8>::zeros();
        // 100 meters uncertainty in initial position
        p[(0, 0)] = 1e4;
        p[(1, 1)] = 1e4;
        p[(2, 2)] = 1e4;
        // 10 m/s uncertainty in initial velocity
        p[(3, 3)] = 1e2;
        p[(4, 4)] = 1e2;
        p[(5, 5)] = 1e2;
        // 1 second uncertainty in initial clock bias
        p[(6, 6)] = 1.0;
        // 10 PPM (1e-5) uncertainty in initial clock drift
        p[(7, 7)] = 1e-10;

        Self {
            x,
            p,
            q_acc: 1.0,           // 1 m^2/s^3 velocity walk
            q_clk_bias: 1e-12,    // TCXO phase walk
            q_clk_drift: 1e-14,   // TCXO drift walk
            r_phase_base: 0.0025, // (0.05 m)^2
            r_freq_base: 0.04,    // (0.2 m/s)^2
        }
    }

    /// Propagate state covariance and state vector by time step dt (seconds)
    pub fn predict(&mut self, dt: f64) {
        if dt <= 0.0 {
            return;
        }

        // 1. Propagate state using Newtonian kinematics
        let mut f = SMatrix::<f64, 8, 8>::identity();
        f[(0, 3)] = dt;
        f[(1, 4)] = dt;
        f[(2, 5)] = dt;
        f[(6, 7)] = dt;

        self.x = f * self.x;

        // 2. Compute process noise matrix Q
        let mut q = SMatrix::<f64, 8, 8>::zeros();
        let dt2 = dt * dt;
        let dt3_3 = dt2 * dt / 3.0;
        let dt2_2 = dt2 / 2.0;

        // Kinematic state process noise from velocity random walk
        for i in 0..3 {
            q[(i, i)] = self.q_acc * dt3_3;
            q[(i, i + 3)] = self.q_acc * dt2_2;
            q[(i + 3, i)] = self.q_acc * dt2_2;
            q[(i + 3, i + 3)] = self.q_acc * dt;
        }

        // Clock state process noise
        q[(6, 6)] = self.q_clk_bias * dt + self.q_clk_drift * dt3_3;
        q[(6, 7)] = self.q_clk_drift * dt2_2;
        q[(7, 6)] = self.q_clk_drift * dt2_2;
        q[(7, 7)] = self.q_clk_drift * dt;

        // Propagate covariance
        self.p = f * self.p * f.transpose() + q;

        // Force symmetry to maintain numerical stability
        self.p = (self.p + self.p.transpose()) * 0.5;
    }

    /// Perform sequential 1D measurement update for a single channel's residuals
    pub fn update_channel(
        &mut self,
        sat_pos: [f64; 3],
        _sat_vel: [f64; 3],
        phase_residual_rad: f64,
        freq_residual_hz: f64,
        carrier_freq_hz: f64,
        noise_scale: f64,
    ) {
        if carrier_freq_hz <= 0.0 {
            return;
        }

        let lambda = C / carrier_freq_hz;

        // 1. Geometry calculations
        let dx = sat_pos[0] - self.x[0];
        let dy = sat_pos[1] - self.x[1];
        let dz = sat_pos[2] - self.x[2];
        let rho = (dx * dx + dy * dy + dz * dz).sqrt();
        if rho < 1.0 {
            return;
        }

        // Unit line-of-sight vector from receiver to satellite
        let ux = dx / rho;
        let uy = dy / rho;
        let uz = dz / rho;

        // 2. Carrier phase (range) measurement update
        // Wavelength conversion: 1 rad of phase = lambda / (2*pi) meters of range
        let y_range = - (lambda / (2.0 * std::f64::consts::PI)) * phase_residual_rad;
        let mut h_range = SVector::<f64, 8>::zeros();
        h_range[0] = ux;
        h_range[1] = uy;
        h_range[2] = uz;
        h_range[6] = C; // clock bias sensitivity (scaled by speed of light)

        let r_range = self.r_phase_base * noise_scale;
        self.update_1d(y_range, &h_range, r_range);

        // 3. Doppler frequency (range rate) measurement update
        // Wavelength conversion: 1 Hz of frequency = lambda meters/sec of range rate
        let y_rate = - lambda * freq_residual_hz;
        let mut h_rate = SVector::<f64, 8>::zeros();
        h_rate[3] = ux;
        h_rate[4] = uy;
        h_rate[5] = uz;
        h_rate[7] = C; // clock drift sensitivity (scaled by speed of light)

        let r_rate = self.r_freq_base * noise_scale;
        self.update_1d(y_rate, &h_rate, r_rate);
    }

    /// Internal 1-dimensional EKF update using Joseph form covariance propagation
    fn update_1d(&mut self, y: f64, h: &SVector<f64, 8>, r: f64) {
        let s = (h.transpose() * self.p * h)[0] + r;
        if s.abs() >= 1e-12 {
            let k = (self.p * h) / s;
            self.x += k * y;
            let i = SMatrix::<f64, 8, 8>::identity();
            let a = i - k * h.transpose();
            self.p = a * self.p * a.transpose() + k * r * k.transpose();
            
            // Force symmetry and enforce covariance floor
            self.p = (self.p + self.p.transpose()) * 0.5;
            for idx in 0..8 {
                self.p[(idx, idx)] = self.p[(idx, idx)].max(1e-15);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_master_nav_ekf_predict_and_update() {
        let init_pos = [0.0, 0.0, 6378137.0];
        let mut ekf = MasterNavEkf::new(init_pos);

        // Verify initial state
        assert_eq!(ekf.x[0], 0.0);
        assert_eq!(ekf.x[1], 0.0);
        assert_eq!(ekf.x[2], 6378137.0);
        assert_eq!(ekf.x[3], 0.0);
        assert_eq!(ekf.x[7], 0.0);

        // Record initial diagonal values of P
        let init_cov_pos = ekf.p[(2, 2)];

        // Set a small clock bias uncertainty in test to isolate position updates
        ekf.p[(6, 6)] = 1e-16;

        // Predict
        ekf.predict(1.0);
        
        // P variance should increase due to process noise
        assert!(ekf.p[(2, 2)] > init_cov_pos);

        // Simulated satellite overhead at [0.0, 0.0, 7000000.0]
        let sat_pos = [0.0, 0.0, 7000000.0];
        let sat_vel = [0.0, 2000.0, 0.0]; // Moving in Y direction at 2 km/s

        // Perform an update step
        // phase_residual_rad = 1.5 rad, freq_residual_hz = 10.0 Hz
        ekf.update_channel(
            sat_pos,
            sat_vel,
            1.5,
            10.0,
            150800000.0,
            1.0,
        );

        // Verify that covariance has decreased (converged) after the measurement update
        assert!(ekf.p[(2, 2)] < init_cov_pos);

        // Check that state vector is updated
        assert!(ekf.x[2] != 6378137.0);
    }
}

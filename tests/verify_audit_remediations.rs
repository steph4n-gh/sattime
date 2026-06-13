#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use num_complex::Complex;
use chrono::TimeZone;
use sattime::daemon::{ConsensusSteeringEngine, CompletedPassData};
use sattime::orbit::{
    datetime_to_jd, apply_sagnac_correction, saastamoinen_tropospheric_delay,
    wgs84_to_ecef
};
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use sattime::dsp::DigitalDownConverter;

#[test]
fn test_bounded_consensus_engine_history() {
    let mut engine = ConsensusSteeringEngine::new();
    assert_eq!(engine.passes.len(), 0);

    for i in 0..60 {
        engine.add_pass_result(CompletedPassData {
            sat_name: format!("SAT_{}", i),
            timestamp: chrono::Utc::now(),
            offset_seconds: 0.002,
            freq_drift_ppm: 0.01,
            snr: 15.0,
            max_elevation: 60.0,
            fit_rmse: 0.5,
        });
    }

    assert_eq!(engine.passes.len(), 50);
    // Should contain SAT_10 through SAT_59 (since SAT_0 to SAT_9 were popped)
    assert_eq!(engine.passes[0].sat_name, "SAT_10");
    assert_eq!(engine.passes[49].sat_name, "SAT_59");
}

#[test]
fn test_two_part_julian_date_precision() {
    // 2000-01-01T12:00:00Z is exactly Julian Date 2451545.0
    let dt = chrono::Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
    let jd = datetime_to_jd(dt);
    
    // Day fraction at noon is 0.5
    assert_eq!(jd.1, 0.5);
    // Integer base should be 2451544.5 (since base + fraction = 2451545.0)
    assert_eq!(jd.0, 2451544.5);
    assert_eq!(jd.0 + jd.1, 2451545.0);

    // 2000-01-01T18:00:00Z should have day fraction 0.75
    let dt2 = chrono::Utc.with_ymd_and_hms(2000, 1, 1, 18, 0, 0).unwrap();
    let jd2 = datetime_to_jd(dt2);
    assert_eq!(jd2.1, 0.75);
    assert_eq!(jd2.0, 2451544.5);
    assert_eq!(jd2.0 + jd2.1, 2451545.25);
}

#[test]
fn test_sagnac_correction_correctness() {
    let pos_sat = [7000e3, 0.0, 0.0];
    let vel_sat = [0.0, 7500.0, 0.0];
    let pos_obs = [6378e3, 0.0, 0.0];

    // Sat to obs distance: 622 km. Time of flight: ~2.07 ms.
    // Sagnac angle should be roughly -1.5e-7 rad.
    let (pos_corr, vel_corr) = apply_sagnac_correction(pos_sat, vel_sat, pos_obs);
    
    // Rotation should affect y (since it rotates around Z)
    assert_ne!(pos_corr[0], pos_sat[0]);
    assert_ne!(pos_corr[1], pos_sat[1]);
    assert_eq!(pos_corr[2], pos_sat[2]);

    assert_ne!(vel_corr[0], vel_sat[0]);
    assert_ne!(vel_corr[1], vel_sat[1]);
    assert_eq!(vel_corr[2], vel_sat[2]);
}

#[test]
fn test_saastamoinen_delay_output_values() {
    let obs_ecef = wgs84_to_ecef(45.0, -75.0, 100.0);
    // Place satellite exactly overhead
    let sat_ecef = [obs_ecef[0] * 1.1, obs_ecef[1] * 1.1, obs_ecef[2] * 1.1];

    let delay = saastamoinen_tropospheric_delay(sat_ecef, obs_ecef);
    // At zenith (elevation = 90 deg), delay is roughly 2.3 / (1.0 + 0.00143) = 2.296 meters
    assert!((delay - 2.296).abs() < 0.01);
}

#[test]
fn test_avx2_ddc_equivalence() {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            let mut ddc_scalar = DigitalDownConverter::new();
            let mut ddc_avx2 = DigitalDownConverter::new();

            // Generate test input
            let mut input = vec![Complex::new(0.0f32, 0.0f32); 1024];
            for i in 0..1024 {
                let phi = (i as f32) * 0.05;
                input[i] = Complex::new(phi.cos(), phi.sin());
            }

            let mut out_scalar = vec![Complex::new(0.0f32, 0.0f32); 1024];
            let mut out_avx2 = vec![Complex::new(0.0f32, 0.0f32); 1024];

            let f_shift = 15000.0;
            let sample_rate = 2e6;

            // Run process on both
            ddc_scalar.process(&input, f_shift, sample_rate, &mut out_scalar);
            
            // Run avx2 directly
            unsafe {
                ddc_avx2.process_avx2(&input, f_shift, sample_rate, &mut out_avx2);
            }

            // Compare outputs
            for i in 0..1024 {
                let diff_re = (out_scalar[i].re - out_avx2[i].re).abs();
                let diff_im = (out_scalar[i].im - out_avx2[i].im).abs();
                assert!(diff_re < 1e-4, "Mismatch at index {} re: {} vs {}", i, out_scalar[i].re, out_avx2[i].re);
                assert!(diff_im < 1e-4, "Mismatch at index {} im: {} vs {}", i, out_scalar[i].im, out_avx2[i].im);
            }

            assert!((ddc_scalar.phase_acc - ddc_avx2.phase_acc).abs() < 1e-4);
        }
    }
}

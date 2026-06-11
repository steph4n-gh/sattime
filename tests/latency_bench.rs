use crossbeam_channel::bounded;
use num_complex::Complex;
use std::io;
use std::time::{Duration, Instant};

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

        // Handle 2-line format
        if line.starts_with('1') && i + 1 < lines.len() && lines[i + 1].trim().starts_with('2') {
            let line1 = line;
            let line2 = lines[i + 1].trim();
            match sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes()) {
                Ok(elements) => {
                    let sat_name = elements
                        .object_name
                        .clone()
                        .unwrap_or_else(|| "UNKNOWN".to_string());
                    satellites.push((sat_name, elements));
                }
                Err(e) => {
                    eprintln!("Error parsing 2-line TLE at line {}: {:?}", i + 1, e);
                }
            }
            i += 2;
        } else if i + 2 < lines.len() {
            let name = line.to_string();
            let line1 = lines[i + 1].trim();
            let line2 = lines[i + 2].trim();

            if line1.starts_with('1') && line2.starts_with('2') {
                match sgp4::Elements::from_tle(
                    Some(name.clone()),
                    line1.as_bytes(),
                    line2.as_bytes(),
                ) {
                    Ok(elements) => {
                        satellites.push((name, elements));
                    }
                    Err(e) => {
                        eprintln!("Error parsing TLE for '{}': {:?}", name, e);
                    }
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

#[test]
fn test_latency_measurements() {
    let tle_path = "passes/starlink.tle";

    // 1. Measure TLE parsing latency (Profile Switching Scenario)
    println!("--- 1. TLE Parsing Latency (Profile Switch) ---");
    let start_parse = Instant::now();
    let satellites = match parse_tle_file(tle_path) {
        Ok(sats) => sats,
        Err(e) => {
            panic!("Failed to parse TLE file {}: {:?}", tle_path, e);
        }
    };
    let parse_duration = start_parse.elapsed();
    println!(
        "Parsed {} satellites from {} in {:?}",
        satellites.len(),
        tle_path,
        parse_duration
    );

    // 2. Measure Satellites Cloning Latency (Rise Schedule Calculation Step, step_count % 500 == 0)
    println!("\n--- 2. Satellites Cloning Latency (step_count % 500 == 0) ---");
    let start_clone = Instant::now();
    let satellites_cloned = satellites.clone();
    let clone_duration = start_clone.elapsed();
    println!(
        "Cloned {} satellites in {:?}",
        satellites_cloned.len(),
        clone_duration
    );

    // 3. Channel behavior and backpressure verification
    println!("\n--- 3. Channel Buffer & Backpressure Simulation ---");
    let (tx, rx) = bounded::<Vec<Complex<f32>>>(4);

    // Simulate main thread blocking (e.g. during TLE parsing on profile switch)
    let main_blocked_time = parse_duration;

    // Spawn SDR thread simulating real-time streaming at 2 MSPS
    // Block size: 32768 samples.
    // Time per block: 32768 / 2,000,000 = 16.384 ms.
    let block_size = 32768;
    let sample_rate = 2_000_000.0;
    let block_time = Duration::from_secs_f64(block_size as f64 / sample_rate);

    let tx_clone = tx.clone();
    let producer = std::thread::spawn(move || {
        let mut sent_count = 0;
        let mut drop_count = 0;
        let start = Instant::now();

        for _ in 0..10 {
            // Simulate the SDR thread generating samples and trying to send them
            let samples = vec![Complex::new(0.0f32, 0.0f32); block_size];

            // In the real code, it does a blocking send:
            // if tx_clone.send(samples).is_err() { break; }
            // Let's measure if it blocks
            let send_start = Instant::now();
            let res = tx_clone.send_timeout(samples, Duration::from_millis(1));
            let send_elapsed = send_start.elapsed();

            if res.is_err() {
                drop_count += 1;
            } else {
                sent_count += 1;
            }

            // Sleep to simulate real-time interval
            let elapsed_since_start = start.elapsed();
            let target_time = block_time * (sent_count + drop_count);
            if target_time > elapsed_since_start {
                std::thread::sleep(target_time - elapsed_since_start);
            }
        }
        (sent_count, drop_count)
    });

    // Simulate main thread taking `main_blocked_time` to process / parse
    println!("Main thread simulating block of {:?}", main_blocked_time);
    std::thread::sleep(main_blocked_time);

    // Now main thread starts reading
    let mut received = 0;
    while let Ok(_) = rx.try_recv() {
        received += 1;
    }

    let (sent, drops) = producer.join().unwrap();
    println!(
        "SDR Producer sent: {}, dropped (or timeout/blocked): {}",
        sent, drops
    );
    println!("Main Consumer received: {}", received);

    // Assertions and limits check
    println!("\n--- Latency and Safety Verdict ---");
    let limit_ms = 1.8;
    let clone_ms = clone_duration.as_secs_f64() * 1000.0;
    let parse_ms = parse_duration.as_secs_f64() * 1000.0;

    println!(
        "Satellite list clone time: {:.3} ms (Limit: {:.3} ms)",
        clone_ms, limit_ms
    );
    println!("TLE parse time: {:.3} ms (Blocks main thread)", parse_ms);

    // Let's print the status explicitly
    if clone_ms >= limit_ms {
        println!(
            "[FAIL] Satellites cloning duration ({:.3} ms) exceeds the step latency budget of {:.3} ms!",
            clone_ms, limit_ms
        );
    } else {
        println!(
            "[PASS] Satellites cloning duration ({:.3} ms) is within the step latency budget of {:.3} ms.",
            clone_ms, limit_ms
        );
    }

    if parse_ms > 65.5 {
        println!(
            "[FAIL] Profile switching TLE parse duration ({:.3} ms) exceeds the SDR thread buffer depth (65.5 ms at 2 MSPS), causing channel blockage / dropped samples!",
            parse_ms
        );
    } else {
        println!(
            "[PASS] Profile switching TLE parse duration ({:.3} ms) is within the SDR thread buffer depth.",
            parse_ms
        );
    }
}

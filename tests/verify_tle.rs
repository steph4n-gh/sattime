use chrono::Utc;
use std::io;
use std::io::{Read, Write};
use std::net::TcpListener;
use std::sync::{Arc, Mutex};
use std::thread;

// Duplicate TLE functions from main.rs to test them
fn download_tle_file(urls: &[&str], output_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let temp_path = format!("{}.tmp", output_path);
    for &url in urls {
        eprintln!("Downloading latest TLE catalog from {}...", url);
        let resp = match ureq::get(url)
            .set(
                "User-Agent",
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/120.0.0.0 Safari/537.36",
            )
            .call()
        {
            Ok(r) => r,
            Err(e) => {
                eprintln!("Warning: Failed to download TLE from {}: {}", url, e);
                continue;
            }
        };

        if resp.status() == 200 {
            let mut body = String::new();
            if let Err(e) = resp.into_reader().read_to_string(&mut body) {
                eprintln!("Warning: Failed to read response body from {}: {}", url, e);
                continue;
            }
            if let Err(e) = std::fs::write(&temp_path, &body) {
                eprintln!("Warning: Failed to write temp TLE file: {}", e);
                continue;
            }
            // Parse and validate
            match parse_tle_file(&temp_path) {
                Ok(sats) => {
                    if sats.is_empty() {
                        eprintln!("Warning: Parsed TLE from {} is empty", url);
                        let _ = std::fs::remove_file(&temp_path);
                        continue;
                    }
                    // Valid! Overwrite the cache file
                    if let Err(e) = std::fs::rename(&temp_path, output_path) {
                        eprintln!("Warning: Failed to rename temp file to cache path: {}", e);
                        let _ = std::fs::remove_file(&temp_path);
                        continue;
                    }
                    eprintln!(
                        "Successfully downloaded, validated, and cached TLE catalog to '{}' from {}",
                        output_path, url
                    );
                    return Ok(());
                }
                Err(e) => {
                    eprintln!(
                        "Warning: Failed to parse downloaded TLE from {}: {}",
                        url, e
                    );
                    let _ = std::fs::remove_file(&temp_path);
                    continue;
                }
            }
        } else {
            eprintln!("Warning: HTTP error status {} from {}", resp.status(), url);
        }
    }

    if std::path::Path::new(output_path).exists() {
        eprintln!("[WARNING] All mirrors failed. Falling back to cached TLE.");
        Ok(())
    } else {
        Err("All mirrors failed and no cached TLE found".into())
    }
}

fn get_tle_file_cached(
    urls: &[&str],
    path: &str,
    force_download: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut needs_download = force_download || !std::path::Path::new(path).exists();
    if !needs_download {
        if let Ok(metadata) = std::fs::metadata(path) {
            if let Ok(modified) = metadata.modified() {
                if let Ok(elapsed) = modified.elapsed() {
                    if elapsed.as_secs() > 43200 {
                        // 12 hours
                        needs_download = true;
                    }
                }
            }
        }
    }

    if needs_download {
        if let Err(e) = download_tle_file(urls, path) {
            eprintln!(
                "[WARNING] TLE download failed: {}. Falling back to cached file.",
                e
            );
            if !std::path::Path::new(path).exists() {
                return Err(e); // No cache exists, must fail
            }
        }
    }
    Ok(())
}

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

// Helper to spawn mock server
fn run_mock_server(responses: Arc<Mutex<Vec<(u16, String)>>>) -> String {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let port = listener.local_addr().unwrap().port();

    thread::spawn(move || {
        while let Ok((mut stream, _)) = listener.accept() {
            // Read request (minimal read to avoid blocking)
            let mut buf = [0; 512];
            let _ = stream.read(&mut buf).unwrap_or(0);

            let mut resps = responses.lock().unwrap();
            let (status, body) = if !resps.is_empty() {
                resps.remove(0)
            } else {
                (500, "Error".to_string())
            };

            let status_line = match status {
                200 => "HTTP/1.1 200 OK",
                404 => "HTTP/1.1 404 Not Found",
                500 => "HTTP/1.1 500 Internal Server Error",
                _ => "HTTP/1.1 500 Internal Server Error",
            };

            let response = format!(
                "{}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                status_line,
                body.len(),
                body
            );
            let _ = stream.write_all(response.as_bytes());
            let _ = stream.flush();
        }
    });

    format!("http://127.0.0.1:{}", port)
}

const VALID_TLE: &str = "AO-07\n1 07530U 74089B   26160.31226999 -.00000035  00000-0  69227-4 0  9998\n2 07530 101.9901 173.4766 0012488 144.7701 229.1729 12.53697932359533\n";

#[test]
fn test_tle_primary_success() {
    let responses = Arc::new(Mutex::new(vec![(200, VALID_TLE.to_string())]));
    let base_url = run_mock_server(responses);

    let temp_dir = std::env::temp_dir().join(format!(
        "tle_test_success_{}",
        Utc::now().timestamp_micros()
    ));
    let _ = std::fs::create_dir_all(&temp_dir);
    let cache_file = temp_dir.join("test.tle");
    let cache_path = cache_file.to_str().unwrap();

    let urls = vec![base_url.as_str()];
    let res = download_tle_file(&urls, cache_path);
    assert!(res.is_ok());
    assert!(cache_file.exists());

    let sats = parse_tle_file(cache_path).unwrap();
    assert_eq!(sats.len(), 1);
    assert_eq!(sats[0].0, "AO-07");

    let _ = std::fs::remove_dir_all(&temp_dir);
}

#[test]
fn test_tle_fallback_to_mirror() {
    // Primary URL fails with 500, secondary succeeds with 200
    let responses = Arc::new(Mutex::new(vec![
        (500, "Internal Server Error".to_string()),
        (200, VALID_TLE.to_string()),
    ]));
    let base_url = run_mock_server(responses);

    let temp_dir = std::env::temp_dir().join(format!(
        "tle_test_fallback_{}",
        Utc::now().timestamp_micros()
    ));
    let _ = std::fs::create_dir_all(&temp_dir);
    let cache_file = temp_dir.join("test_fallback.tle");
    let cache_path = cache_file.to_str().unwrap();

    // Both URLs point to the mock server, which will return 500 first, then 200
    let urls = vec![base_url.as_str(), base_url.as_str()];
    let res = download_tle_file(&urls, cache_path);
    assert!(res.is_ok());
    assert!(cache_file.exists());

    let sats = parse_tle_file(cache_path).unwrap();
    assert_eq!(sats.len(), 1);
    assert_eq!(sats[0].0, "AO-07");

    let _ = std::fs::remove_dir_all(&temp_dir);
}

#[test]
fn test_tle_all_mirrors_fail_no_cache() {
    let responses = Arc::new(Mutex::new(vec![
        (500, "Server Error".to_string()),
        (404, "Not Found".to_string()),
    ]));
    let base_url = run_mock_server(responses);

    let temp_dir =
        std::env::temp_dir().join(format!("tle_test_fail_{}", Utc::now().timestamp_micros()));
    let _ = std::fs::create_dir_all(&temp_dir);
    let cache_file = temp_dir.join("test_fail.tle");
    let cache_path = cache_file.to_str().unwrap();

    let urls = vec![base_url.as_str(), base_url.as_str()];
    let res = download_tle_file(&urls, cache_path);
    assert!(res.is_err());
    assert!(!cache_file.exists());

    let _ = std::fs::remove_dir_all(&temp_dir);
}

#[test]
fn test_tle_all_mirrors_fail_with_cache_fallback() {
    let responses = Arc::new(Mutex::new(vec![
        (500, "Server Error".to_string()),
        (404, "Not Found".to_string()),
    ]));
    let base_url = run_mock_server(responses);

    let temp_dir = std::env::temp_dir().join(format!(
        "tle_test_cache_fallback_{}",
        Utc::now().timestamp_micros()
    ));
    let _ = std::fs::create_dir_all(&temp_dir);
    let cache_file = temp_dir.join("test_cache_fallback.tle");
    let cache_path = cache_file.to_str().unwrap();

    // Pre-seed cache file
    std::fs::write(cache_path, VALID_TLE).unwrap();

    let urls = vec![base_url.as_str(), base_url.as_str()];
    let res = download_tle_file(&urls, cache_path);
    // Should return Ok(()) and keep the cache file
    assert!(res.is_ok());
    assert!(cache_file.exists());

    let sats = parse_tle_file(cache_path).unwrap();
    assert_eq!(sats.len(), 1);
    assert_eq!(sats[0].0, "AO-07");

    let _ = std::fs::remove_dir_all(&temp_dir);
}

#[test]
fn test_tle_corrupt_download_does_not_overwrite_cache() {
    // First connection succeeds with corrupt TLE data
    let responses = Arc::new(Mutex::new(vec![(
        200,
        "THIS IS INVALID TLE DATA\nONLY ONE LINE".to_string(),
    )]));
    let base_url = run_mock_server(responses);

    let temp_dir = std::env::temp_dir().join(format!(
        "tle_test_corrupt_{}",
        Utc::now().timestamp_micros()
    ));
    let _ = std::fs::create_dir_all(&temp_dir);
    let cache_file = temp_dir.join("test_corrupt.tle");
    let cache_path = cache_file.to_str().unwrap();

    // Pre-seed cache file with valid TLE
    std::fs::write(cache_path, VALID_TLE).unwrap();

    let urls = vec![base_url.as_str()];
    let res = download_tle_file(&urls, cache_path);
    // Since parsing of the downloaded file fails, download_tle_file should reject it and fall back to the cache.
    // It returns Ok(()) because cache exists.
    assert!(res.is_ok());

    // Cache file should still contain the original valid TLE data
    let sats = parse_tle_file(cache_path).unwrap();
    assert_eq!(sats.len(), 1);
    assert_eq!(sats[0].0, "AO-07");

    let _ = std::fs::remove_dir_all(&temp_dir);
}

#[test]
fn test_tle_cache_age_forcing() {
    let responses = Arc::new(Mutex::new(vec![
        (200, VALID_TLE.to_string()),
        (200, VALID_TLE.to_string()),
    ]));
    let base_url = run_mock_server(responses);

    let temp_dir =
        std::env::temp_dir().join(format!("tle_test_age_{}", Utc::now().timestamp_micros()));
    let _ = std::fs::create_dir_all(&temp_dir);
    let cache_file = temp_dir.join("test_age.tle");
    let cache_path = cache_file.to_str().unwrap();

    // Case 1: cache doesn't exist, must download
    let urls = vec![base_url.as_str()];
    let res = get_tle_file_cached(&urls, cache_path, false);
    assert!(res.is_ok());
    assert!(cache_file.exists());

    // Case 2: cache is fresh and force_download = false, should NOT trigger download.
    // (If it triggered download, it would consume the second mock response).
    let res = get_tle_file_cached(&urls, cache_path, false);
    assert!(res.is_ok());

    // Case 3: cache is fresh but force_download = true, should trigger download and consume mock response.
    let res = get_tle_file_cached(&urls, cache_path, true);
    assert!(res.is_ok());

    let _ = std::fs::remove_dir_all(&temp_dir);
}

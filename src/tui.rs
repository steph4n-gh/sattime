use crossterm::{
    event::{DisableMouseCapture, EnableMouseCapture},
    execute,
    terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode},
};
use ratatui::{
    Terminal,
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout},
    style::{Color, Modifier, Style},
    symbols,
    text::{Line, Span},
    widgets::{
        Axis, Block, Borders, Cell, Chart, Dataset, Paragraph, Row, Scrollbar,
        ScrollbarOrientation, ScrollbarState, Sparkline, Table,
    },
};

use crate::daemon::*;
use crate::orbit::*;
use crate::{Args, get_process_rss_mb};
use crossbeam_channel::Receiver;
use std::collections::VecDeque;
use std::io::{self, Write};

const CHANNEL_COLORS: [Color; 8] = [
    Color::Yellow,
    Color::LightMagenta,
    Color::Cyan,
    Color::LightGreen,
    Color::LightRed,
    Color::LightBlue,
    Color::LightYellow,
    Color::White,
];

pub struct TuiState {
    pub history_offsets: VecDeque<(f64, f64)>,
    pub channel_histories: Vec<VecDeque<(f64, f64)>>,
    pub spectrum: Vec<f32>,
    pub was_locked: bool,
    pub active_sat_trail: VecDeque<(f64, f64)>,
    pub active_sat_name: Option<String>,
    pub frames_since_active: usize,
    pub telemetry_scroll_offset: usize,
}

pub fn format_slider(value: f64, max_val: f64, width: usize) -> String {
    let filled = if max_val > 0.0 {
        (((value / max_val) * width as f64).round() as usize).min(width)
    } else {
        0
    };
    let mut bar = String::new();
    bar.push('[');
    for i in 0..width {
        if i < filled {
            bar.push('█');
        } else {
            bar.push('░');
        }
    }
    bar.push(']');
    format!("{} {:.0} dB", bar, value)
}

pub struct TerminalGuard {
    #[cfg(unix)]
    pub original_stderr: Option<libc::c_int>,
}

impl Drop for TerminalGuard {
    fn drop(&mut self) {
        let _ = disable_raw_mode();
        let mut stdout = io::stdout();
        let _ = execute!(stdout, LeaveAlternateScreen, DisableMouseCapture);
        let _ = execute!(stdout, crossterm::cursor::Show);
        #[cfg(unix)]
        if let Some(fd) = self.original_stderr
            && fd >= 0
        {
            unsafe {
                libc::dup2(fd, libc::STDERR_FILENO);
                libc::close(fd);
            }
        }
    }
}

pub static VISIBLE_SATS: std::sync::OnceLock<std::sync::Mutex<Vec<VisibleSat>>> =
    std::sync::OnceLock::new();

#[derive(Clone, Debug)]
pub struct VisibleSat {
    pub name: String,
    pub az: f64,
    pub el: f64,
    pub freq_expected: f64,
    pub freq_expected2: f64,
    pub range: f64,
    pub pass_progress: f64,
}

pub fn get_visible_sats() -> &'static std::sync::Mutex<Vec<VisibleSat>> {
    VISIBLE_SATS.get_or_init(|| std::sync::Mutex::new(Vec::new()))
}

#[derive(Clone, Debug)]
pub struct ChannelTelemetry {
    pub id: usize,
    pub sat_name: String,
    pub status: String,
    pub target_freq: f64,
    pub freq_offset: f64,
    pub doppler_rate: f64,
    pub snr_db: f32,
    pub tec: f64,
    pub is_dual: bool,
}

pub struct TuiManager {
    pub terminal: Terminal<CrosstermBackend<io::Stdout>>,
    pub _guard: TerminalGuard,
    pub state: TuiState,
    pub log_rx: Option<Receiver<String>>,
    pub event_logs: VecDeque<String>,
}

impl TuiManager {
    pub fn new(
        fft_size: usize,
        log_rx: Option<Receiver<String>>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        enable_raw_mode()?;
        let mut stdout = io::stdout();
        execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
        let backend = CrosstermBackend::new(stdout);
        let mut terminal = Terminal::new(backend)?;
        terminal.clear()?;
        terminal.hide_cursor()?;

        #[cfg(unix)]
        let original_stderr_fd = unsafe { libc::dup(libc::STDERR_FILENO) };
        #[cfg(unix)]
        if original_stderr_fd >= 0 {
            let dev_null = std::fs::OpenOptions::new().write(true).open("/dev/null");
            if let Ok(file) = dev_null {
                use std::os::unix::io::AsRawFd;
                unsafe {
                    libc::dup2(file.as_raw_fd(), libc::STDERR_FILENO);
                }
            }
        }

        let default_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            let _ = disable_raw_mode();
            let mut stdout = io::stdout();
            let _ = execute!(stdout, LeaveAlternateScreen, DisableMouseCapture);
            let _ = execute!(stdout, crossterm::cursor::Show);
            #[cfg(unix)]
            if original_stderr_fd >= 0 {
                unsafe {
                    libc::dup2(original_stderr_fd, libc::STDERR_FILENO);
                }
            }
            default_hook(info);
        }));

        Ok(Self {
            terminal,
            _guard: TerminalGuard {
                #[cfg(unix)]
                original_stderr: if original_stderr_fd >= 0 {
                    Some(original_stderr_fd)
                } else {
                    None
                },
            },
            state: TuiState {
                history_offsets: VecDeque::with_capacity(350),
                channel_histories: Vec::new(),
                spectrum: vec![0.0f32; fft_size],
                was_locked: false,
                active_sat_trail: VecDeque::with_capacity(20),
                active_sat_name: None,
                frames_since_active: 0,
                telemetry_scroll_offset: 0,
            },
            log_rx,
            event_logs: VecDeque::with_capacity(100),
        })
    }

    pub fn suspend(&mut self) {
        let _ = disable_raw_mode();
        let mut stdout = io::stdout();
        let _ = execute!(stdout, LeaveAlternateScreen, DisableMouseCapture);
        let _ = execute!(stdout, crossterm::cursor::Show);
        #[cfg(unix)]
        if let Some(fd) = self._guard.original_stderr {
            unsafe {
                libc::dup2(fd, libc::STDERR_FILENO);
            }
        }
    }

    pub fn resume(&mut self) {
        let _ = enable_raw_mode();
        let mut stdout = io::stdout();
        let _ = execute!(stdout, EnterAlternateScreen, EnableMouseCapture);
        let _ = self.terminal.clear();
        let _ = self.terminal.hide_cursor();
        #[cfg(unix)]
        if self._guard.original_stderr.is_some() {
            let dev_null = std::fs::OpenOptions::new().write(true).open("/dev/null");
            if let Ok(file) = dev_null {
                use std::os::unix::io::AsRawFd;
                unsafe {
                    libc::dup2(file.as_raw_fd(), libc::STDERR_FILENO);
                }
            }
        }
    }

    pub fn draw(
        &mut self,
        state_name: &str,
        elapsed_seconds: f64,
        tracked_peak: f64,
        freq_offset: f64,
        doppler_rate: f64,
        snr_db: f32,
        _raw_snr_db: f32,
        _overhead_count: &str,
        current_lna_gain: f64,
        current_vga_gain: f64,
        amp_gain: f64,
        args: &Args,
        current_freq: f64,
        active_profile_name: &str,
        _active_guided_window: f64,
        _satellites: &[(String, crate::orbit::OrbitModel)],
        _buffer_len: usize,
        _last_action: &str,
        agc: bool,
        notch_spurs: bool,
        decimate: usize,
        active_min_snr: f32,
        next_pass_info: &str,
        captured_duration: Option<f64>,
        channels: &[ChannelTelemetry],
    ) {
        let is_locked = state_name == "LOCKED";
        if is_locked && !self.state.was_locked {
            print!("\x07");
            let _ = std::io::stdout().flush();
        }
        self.state.was_locked = is_locked;

        // Update historical trajectories for all channels
        if self.state.channel_histories.len() < channels.len() {
            self.state.channel_histories.resize_with(channels.len(), || VecDeque::with_capacity(350));
        }
        for ch in channels {
            let ch_idx = ch.id;
            if ch_idx >= self.state.channel_histories.len() {
                self.state.channel_histories.resize_with(ch_idx + 1, || VecDeque::with_capacity(350));
            }
            let history = &mut self.state.channel_histories[ch_idx];
            if ch.status != "IDLE" {
                history.push_back((elapsed_seconds, ch.freq_offset));
            } else {
                history.clear();
            }
            if history.len() > 300 {
                history.pop_front();
            }
        }

        if let Some(rx) = &self.log_rx {
            while let Ok(msg) = rx.try_recv() {
                let trimmed = msg.trim_end_matches(|c| c == '\r' || c == '\n').to_string();
                if !trimmed.is_empty() {
                    let is_duplicate = self
                        .event_logs
                        .iter()
                        .any(|existing| trimmed.contains(existing));
                    if !is_duplicate {
                        self.event_logs.push_back(trimmed);
                        if self.event_logs.len() > 50 {
                            self.event_logs.pop_front();
                        }
                    }
                }
            }
        }

        let state = &mut self.state;
        let _ = self.terminal.draw(|rect| {
            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .margin(1)
                .constraints([
                    Constraint::Length(3), // Status Header
                    Constraint::Min(10),   // Main Diagnostics & Charts
                    Constraint::Length(4), // Interactive Help Menu
                ].as_ref())
                .split(rect.size());

            let sat_info = if state_name == "LOCKED" {
                " | ACTIVE SATELLITE: LOCKED".to_string()
            } else {
                "".to_string()
            };
            let header = Paragraph::new(format!(
                " CHRONOS ORBITAL TIME SERVER | FREQ: {:.3} MHz | STATE: {}{}",
                current_freq / 1e6,
                state_name,
                sat_info
            ))
            .block(Block::default().borders(Borders::ALL).title(" System Status "));
            rect.render_widget(header, chunks[0]);

            let main_chunks = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Percentage(30),
                    Constraint::Percentage(70),
                ].as_ref())
                .split(chunks[1]);

            // Split diagnostics sidebar vertically into 3 panels
            let sidebar_chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Percentage(45), // General Diagnostics Table
                    Constraint::Percentage(25), // NTP Steering Cockpit
                    Constraint::Percentage(30), // Geolocation & Reverse-GPS
                ].as_ref())
                .split(main_chunks[0]);

            // Panel 1: General Diagnostics
            let rows = vec![
                Row::new(vec![Cell::from("Profile"), Cell::from(active_profile_name)]).style(Style::default().fg(Color::Cyan)),
                Row::new(vec![Cell::from("Target Freq"), Cell::from(format!("{:.3} MHz", current_freq / 1e6))]),
                Row::new(vec![Cell::from("Tracked Peak"), Cell::from(format!("{:.3} MHz", tracked_peak / 1e6))]),
                Row::new(vec![Cell::from("Freq Offset"), Cell::from(format!("{:+.1} Hz", freq_offset))]),
                Row::new(vec![Cell::from("Doppler Rate"), Cell::from(format!("{:+.1} Hz/s", doppler_rate))]),
                Row::new(vec![Cell::from("SNR (Lock)"), Cell::from(format!("{:.1} dB ({:.1} dB)", snr_db, active_min_snr))])
                    .style(Style::default().fg(if snr_db >= active_min_snr { Color::Green } else { Color::Yellow })),
                Row::new(vec![Cell::from("LNA Gain"), Cell::from(format_slider(current_lna_gain, 40.0, 10))]),
                Row::new(vec![Cell::from("VGA Gain"), Cell::from(format_slider(current_vga_gain, 62.0, 10))]),
                Row::new(vec![Cell::from("AMP Gain"), Cell::from(format_slider(amp_gain, 14.0, 10))]),
                Row::new(vec![Cell::from("AGC"), Cell::from(if agc { "ENABLED" } else { "DISABLED" })]).style(Style::default().fg(if agc { Color::Green } else { Color::Red })),
                Row::new(vec![Cell::from("Spur Notch"), Cell::from(if notch_spurs { "ON" } else { "OFF" })]),
                Row::new(vec![Cell::from("Decimation"), Cell::from(format!("{}x", decimate))]),
                Row::new(vec![Cell::from("RSS Memory"), Cell::from(format!("{:.2} MB", get_process_rss_mb()))]).style(Style::default().fg(Color::Magenta)),
            ];
            let sidebar = Table::new(rows, [Constraint::Percentage(45), Constraint::Percentage(55)])
                .block(Block::default().borders(Borders::ALL).title(" Diagnostics "))
                .header(Row::new(vec![
                    Cell::from("Metric").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                    Cell::from("Value").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                ]));
            rect.render_widget(sidebar, sidebar_chunks[0]);

            // Panel 2: LEODO NTP Steering Status
            let leodo_stats = {
                let lock = get_leodo_loop().lock().unwrap();
                (lock.last_offset, lock.last_freq_err_ppm, lock.last_target_adjustment, lock.last_status.clone())
            };

            let loop_state = if !args.leodo {
                "FREE RUN (DRY)"
            } else if leodo_stats.3.contains("SUCCESS") {
                "LOCKED (SLEWING)"
            } else if leodo_stats.3.contains("EPERM") {
                "EPERM (BLOCKED)"
            } else {
                "WARM UP"
            };

            let ntp_rows = vec![
                Row::new(vec![Cell::from("Loop State"), Cell::from(loop_state)]).style(Style::default().fg(if args.leodo { Color::Green } else { Color::Yellow })),
                Row::new(vec![Cell::from("Offset (dt)"), Cell::from(format!("{:+.6} s", leodo_stats.0))]),
                Row::new(vec![Cell::from("Drift (df)"), Cell::from(format!("{:+.3} PPM", leodo_stats.1))]),
                Row::new(vec![Cell::from("Slew Sched"), Cell::from(format!("{:+.6} s", leodo_stats.2))]),
                Row::new(vec![Cell::from("Last Status"), Cell::from(leodo_stats.3)]),
            ];
            let ntp_panel = Table::new(ntp_rows, [Constraint::Percentage(45), Constraint::Percentage(55)])
                .block(Block::default().borders(Borders::ALL).title(" LEODO NTP Steering "));
            rect.render_widget(ntp_panel, sidebar_chunks[1]);

            // Panel 3: Reverse-GPS Geolocation Status
            let geo_stats = {
                let lock = get_geolocation_result().lock().unwrap();
                (lock.lat, lock.lon, lock.alt, lock.rmse, lock.converged, lock.num_passes, lock.gdop, lock.uncertainty_km)
            };

            let mut captured_passes = 0;
            if let Ok(entries) = std::fs::read_dir(&args.output_dir) {
                for entry in entries.flatten() {
                    let path = entry.path();
                    if path.is_file() && path.extension().is_some_and(|ext| ext == "csv") {
                        captured_passes += 1;
                    }
                }
            }

            let progress_bar = match captured_passes {
                0 => "[░░░░░] 0/3",
                1 => "[██░░░] 1/3",
                2 => "[████░] 2/3",
                3 => "[█████] 3/3",
                4 => "[█████+] 4/5",
                _ => "[█████++] 5/5",
            };

            let convergence_state = if geo_stats.4 {
                "CONVERGED"
            } else if captured_passes >= 3 {
                "READY TO SOLVE"
            } else if captured_passes > 0 {
                "ACQUIRING"
            } else {
                "UNSOLVED"
            };

            let geo_rows = vec![
                Row::new(vec![Cell::from("GPS Lock"), Cell::from(convergence_state)]).style(Style::default().fg(if geo_stats.4 { Color::Green } else { Color::Yellow })),
                Row::new(vec![Cell::from("Passes Saved"), Cell::from(progress_bar)]),
                Row::new(vec![Cell::from("Solved Lat"), Cell::from(if geo_stats.4 { format!("{:.5}°", geo_stats.0) } else { "N/A".to_string() })]),
                Row::new(vec![Cell::from("Solved Lon"), Cell::from(if geo_stats.4 { format!("{:.5}°", geo_stats.1) } else { "N/A".to_string() })]),
                Row::new(vec![Cell::from("Fit RMSE"), Cell::from(if geo_stats.4 { format!("{:.2} Hz", geo_stats.3) } else { "N/A".to_string() })]),
                Row::new(vec![Cell::from("GDOP Quality"), Cell::from(if geo_stats.4 { format!("{:.2}", geo_stats.6) } else { "N/A".to_string() })]),
                Row::new(vec![Cell::from("Uncertainty"), Cell::from(if geo_stats.4 { format!("{:.3} km", geo_stats.7) } else { "N/A".to_string() })]),
            ];
            let geo_panel = Table::new(geo_rows, [Constraint::Percentage(45), Constraint::Percentage(55)])
                .block(Block::default().borders(Borders::ALL).title(" Reverse-GPS Geolocation "));
            rect.render_widget(geo_panel, sidebar_chunks[2]);

            let visuals_chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Percentage(20), // Spectrum Analyzer
                    Constraint::Percentage(30), // Doppler S-Curve Plot
                    Constraint::Percentage(35), // Sky-Track Radar Zenith Plot
                    Constraint::Percentage(15), // Console Logs
                ].as_ref())
                .split(main_chunks[1]);

            let fft_len = state.spectrum.len();
            let display_width = (visuals_chunks[0].width as usize).saturating_sub(2).max(10);
            let mut spark_data = vec![0u64; display_width];
            if fft_len > 0 {
                let mut sum = 0.0f32;
                for &v in &state.spectrum {
                    sum += v;
                }
                let avg = sum / (fft_len as f32);
                let noise_floor = if avg > 1e-6 { avg } else { 1.0f32 };

                let half = display_width / 2;
                // Left half: negative frequencies (bins N/2 to N-1). Group and take max to preserve peak.
                let neg_bins_count = fft_len / 2;
                for i in 0..half {
                    let bin_start = neg_bins_count + i * neg_bins_count / half;
                    let bin_end = neg_bins_count + (i + 1) * neg_bins_count / half;

                    let mut max_mag = 0.0f32;
                    for b in bin_start..bin_end.min(fft_len) {
                        let m = state.spectrum[b];
                        if m > max_mag { max_mag = m; }
                    }

                    let ratio = (max_mag / noise_floor).max(1.0) as f64;
                    let db = 20.0 * ratio.log10();
                    let val = (db * 0.35).clamp(0.0, 7.0) as u64;
                    spark_data[i] = val;
                }
                // Right half: positive frequencies (bins 0 to N/2-1). Group and take max.
                let pos_bins_count = fft_len / 2;
                let remaining = display_width - half;
                for i in 0..remaining {
                    let bin_start = i * pos_bins_count / remaining;
                    let bin_end = (i + 1) * pos_bins_count / remaining;

                    let mut max_mag = 0.0f32;
                    for b in bin_start..bin_end.min(pos_bins_count) {
                        let m = state.spectrum[b];
                        if m > max_mag { max_mag = m; }
                    }

                    let ratio = (max_mag / noise_floor).max(1.0) as f64;
                    let db = 20.0 * ratio.log10();
                    let val = (db * 0.35).clamp(0.0, 7.0) as u64;
                    spark_data[half + i] = val;
                }
            }


            let sparkline = Sparkline::default()
                .block(Block::default().borders(Borders::ALL).title(" Real-Time Spectrum Analyzer (Full Bandwidth) "))
                .data(&spark_data)
                .style(Style::default().fg(Color::Green));
            rect.render_widget(sparkline, visuals_chunks[0]);

            let chart_layout = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Percentage(40), // Doppler Chart
                    Constraint::Percentage(60), // Parallel Channels Status
                ].as_ref())
                .split(visuals_chunks[1]);

            // Convert channel histories to vectors of points for rendering
            let mut channel_points: Vec<Vec<(f64, f64)>> = Vec::new();
            for ch_idx in 0..state.channel_histories.len() {
                let history = &state.channel_histories[ch_idx];
                channel_points.push(history.iter().copied().collect());
            }

            let mut x_min = f64::MAX;
            let mut x_max = f64::MIN;
            let mut y_min = -100.0;
            let mut y_max = 100.0;
            let mut has_points = false;
            let mut min_val = f64::MAX;
            let mut max_val = f64::MIN;

            for pts in &channel_points {
                if !pts.is_empty() {
                    has_points = true;
                    for &(x, val) in pts {
                        if x < x_min { x_min = x; }
                        if x > x_max { x_max = x; }
                        if val < min_val { min_val = val; }
                        if val > max_val { max_val = val; }
                    }
                }
            }

            if has_points {
                let range = max_val - min_val;
                let margin = (range * 0.15).max(20.0); // 15% margin or at least 20 Hz
                y_min = min_val - margin;
                y_max = max_val + margin;
                if x_min >= x_max {
                    x_max = x_min + 1.0;
                }
            } else {
                // Fallback to legacy single-channel history_offsets if no multi-channel lock exists
                let points: Vec<(f64, f64)> = state.history_offsets.iter().copied().collect();
                if !points.is_empty() {
                    has_points = true;
                    x_max = points.last().unwrap().0;
                    x_min = if points.len() >= 300 { points[0].0 } else { 0.0 };
                    for &(_, val) in &points {
                        if val < min_val { min_val = val; }
                        if val > max_val { max_val = val; }
                    }
                    let range = max_val - min_val;
                    let margin = (range * 0.15).max(20.0);
                    y_min = min_val - margin;
                    y_max = max_val + margin;
                } else {
                    x_min = 0.0;
                    x_max = 1.0;
                }
            }

            let mut datasets = Vec::new();
            if has_points && !channel_points.is_empty() {
                for ch in channels {
                    let ch_idx = ch.id;
                    if ch_idx < channel_points.len() && !channel_points[ch_idx].is_empty() {
                        let color = CHANNEL_COLORS[ch_idx % CHANNEL_COLORS.len()];
                        let dataset = Dataset::default()
                            .name(format!("Ch {} ({})", ch.id, if ch.sat_name.is_empty() { "---" } else { &ch.sat_name }))
                            .marker(symbols::Marker::Braille)
                            .style(Style::default().fg(color))
                            .data(&channel_points[ch_idx]);
                        datasets.push(dataset);
                    }
                }
            }

            // Fallback dataset if empty
            let fallback_points: Vec<(f64, f64)> = state.history_offsets.iter().copied().collect();
            if datasets.is_empty() {
                let dataset = Dataset::default()
                    .name("Doppler Drift")
                    .marker(symbols::Marker::Braille)
                    .style(Style::default().fg(Color::Yellow))
                    .data(&fallback_points);
                datasets.push(dataset);
            }

            let chart = Chart::new(datasets)
                .block(Block::default().borders(Borders::ALL).title(" Doppler S-Curve Plot (Hz Offset over Time) "))
                .x_axis(Axis::default()
                    .title("Time (seconds)")
                    .style(Style::default().fg(Color::Gray))
                    .bounds([x_min, x_max.max(1.0)]))
                .y_axis(Axis::default()
                    .title("Offset (Hz)")
                    .style(Style::default().fg(Color::Gray))
                    .bounds([y_min, y_max])
                    .labels(vec![
                        Span::raw(format!("{:.0} Hz", y_min)),
                        Span::raw(format!("{:.0} Hz", (y_min + y_max) / 2.0)),
                        Span::raw(format!("{:.0} Hz", y_max)),
                    ]));
            rect.render_widget(chart, chart_layout[0]);

            let mut channel_rows = Vec::new();
            for ch in channels {
                let status_style = match ch.status.as_str() {
                    "LOCKED" => Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
                    "ACQUISITION" => Style::default().fg(Color::Yellow),
                    "FADE" => Style::default().fg(Color::Red),
                    _ => Style::default().fg(Color::DarkGray),
                };

                let ch_color = CHANNEL_COLORS[ch.id % CHANNEL_COLORS.len()];
                let ch_style = Style::default().fg(ch_color).add_modifier(Modifier::BOLD);

                let row = Row::new(vec![
                    Cell::from(format!("Ch {}", ch.id)).style(ch_style),
                    Cell::from(if ch.sat_name.is_empty() { "---".to_string() } else { ch.sat_name.clone() }),
                    Cell::from(ch.status.clone()).style(status_style),
                    Cell::from(if ch.status == "IDLE" { "---".to_string() } else { format!("{:.3} MHz", ch.target_freq / 1e6) }),
                    Cell::from(if ch.status == "IDLE" { "---".to_string() } else { format!("{:+.1} Hz", ch.freq_offset) }),
                    Cell::from(if ch.status == "IDLE" { "---".to_string() } else { format!("{:+.1} Hz/s", ch.doppler_rate) }),
                    Cell::from(if ch.status == "IDLE" { "---".to_string() } else { format!("{:.1} dB", ch.snr_db) }),
                    Cell::from(if ch.status == "IDLE" || !ch.is_dual { "---".to_string() } else { format!("{:.2} TECU", ch.tec) }),
                ]);
                channel_rows.push(row);
            }

            let ch_header = Row::new(vec![
                Cell::from("Channel").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Satellite").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Status").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Target Freq").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Offset").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Doppler Rate").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("SNR").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("TEC").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
            ]);

            let ch_table = Table::new(
                channel_rows,
                [
                    Constraint::Percentage(10),
                    Constraint::Percentage(20),
                    Constraint::Percentage(12),
                    Constraint::Percentage(14),
                    Constraint::Percentage(12),
                    Constraint::Percentage(12),
                    Constraint::Percentage(8),
                    Constraint::Percentage(12),
                ]
            )
            .header(ch_header)
            .block(Block::default().borders(Borders::ALL).title(" Parallel Demodulation Channels Status "));

            rect.render_widget(ch_table, chart_layout[1]);

            // 1. Initialize empty 15x31 radar grid
            let mut radar_grid = vec![vec![' '; 31]; 15];

            // 2. Draw crosshair guidelines
            for y in 1..14 {
                if y != 7 {
                    radar_grid[y][15] = ':';
                }
            }
            for x in 2..29 {
                if x != 15 {
                    radar_grid[7][x] = '-';
                }
            }

            // 3. Draw concentric elevation rings
            for i in 0..120 {
                let angle = (i as f64) * 2.0 * std::f64::consts::PI / 120.0;

                // Horizon (el = 0.0)
                let x0 = 15.0 + 14.0 * angle.sin();
                let y0 = 7.0 - 7.0 * angle.cos();
                let px0 = x0.round() as usize;
                let py0 = y0.round() as usize;
                if px0 < 31 && py0 < 15 {
                    radar_grid[py0][px0] = '.';
                }

                // 30 deg elevation ring (radius = 2/3 of horizon)
                let x30 = 15.0 + (14.0 * 2.0 / 3.0) * angle.sin();
                let y30 = 7.0 - (7.0 * 2.0 / 3.0) * angle.cos();
                let px30 = x30.round() as usize;
                let py30 = y30.round() as usize;
                if px30 < 31 && py30 < 15 {
                    radar_grid[py30][px30] = '.';
                }

                // 60 deg elevation ring (radius = 1/3 of horizon)
                let x60 = 15.0 + (14.0 / 3.0) * angle.sin();
                let y60 = 7.0 - (7.0 / 3.0) * angle.cos();
                let px60 = x60.round() as usize;
                let py60 = y60.round() as usize;
                if px60 < 31 && py60 < 15 {
                    radar_grid[py60][px60] = '.';
                }
            }

            // 4. Draw Center (Zenith) & Headings
            radar_grid[7][15] = '+';
            radar_grid[0][15] = 'N';
            radar_grid[14][15] = 'S';
            radar_grid[7][0] = 'W';
            radar_grid[7][30] = 'E';



            // Retrieve visible satellites from the background propagator thread cache
            let mut visible_sats = Vec::new();
            {
                let lock = get_visible_sats().lock().unwrap();
                for sat in lock.iter() {
                    let freq_diff = (sat.freq_expected - tracked_peak).abs();
                    visible_sats.push((sat.name.clone(), sat.az, sat.el, freq_diff));
                }
            }

            let closest_idx = if state_name == "LOCKED" && !visible_sats.is_empty() {
                let mut min_diff = f64::MAX;
                let mut best_idx = 0;
                for (idx, &(_, _, _, diff)) in visible_sats.iter().enumerate() {
                    if diff < min_diff {
                        min_diff = diff;
                        best_idx = idx;
                    }
                }
                Some(best_idx)
            } else {
                None
            };

            // Update active satellite trail in TuiState
            let mut active_name = None;
            let mut active_coords = None;
            if let Some(idx) = closest_idx {
                active_name = Some(visible_sats[idx].0.clone());
                active_coords = Some((visible_sats[idx].1, visible_sats[idx].2));
                state.frames_since_active = 0;
            } else {
                state.frames_since_active += 1;
            }

            let should_clear_trail = match (&active_name, &state.active_sat_name) {
                (Some(new_name), Some(old_name)) => new_name != old_name,
                (None, Some(_)) => state.frames_since_active > 75,
                _ => false,
            };

            if should_clear_trail {
                state.active_sat_trail.clear();
                state.active_sat_name = active_name.clone();
            } else if active_name.is_some() {
                state.active_sat_name = active_name.clone();
            }

            if let Some((az, el)) = active_coords {
                let should_push = match state.active_sat_trail.back() {
                    Some(&(last_az, last_el)) => {
                        let d_az = (az - last_az).abs();
                        let d_el = (el - last_el).abs();
                        d_az >= 0.01 || d_el >= 0.01
                    }
                    None => true,
                };
                if should_push {
                    state.active_sat_trail.push_back((az, el));
                    if state.active_sat_trail.len() > 6 {
                        state.active_sat_trail.pop_front();
                    }
                }
            } else if state.frames_since_active > 75 {
                state.active_sat_trail.clear();
                state.active_sat_name = None;
            }

            // Projection helper function
            let project_coord = |az: f64, el: f64| -> (usize, usize) {
                let r_y = 7.0 * (std::f64::consts::FRAC_PI_2 - el) / std::f64::consts::FRAC_PI_2;
                let r_x = 14.0 * (std::f64::consts::FRAC_PI_2 - el) / std::f64::consts::FRAC_PI_2;

                let x = 15.0 + r_x * az.sin();
                let y = 7.0 - r_y * az.cos();

                let px = x.round().clamp(0.0, 30.0) as usize;
                let py = y.round().clamp(0.0, 14.0) as usize;
                (px, py)
            };

            // Plot active satellite trail dots
            if state.active_sat_trail.len() > 1 {
                for i in 0..(state.active_sat_trail.len() - 1) {
                    let (az, el) = state.active_sat_trail[i];
                    let (px, py) = project_coord(az, el);
                    radar_grid[py][px] = '·';
                }
            }

            // Plot candidate satellites
            for (idx, &(_, az, el, _)) in visible_sats.iter().enumerate() {
                if Some(idx) != closest_idx {
                    let (px, py) = project_coord(az, el);
                    radar_grid[py][px] = 's';
                }
            }

            // Plot actively tracked satellite
            if let Some(idx) = closest_idx {
                let (_, az, el, _) = visible_sats[idx];
                let (px, py) = project_coord(az, el);
                radar_grid[py][px] = 'S';
            }

            let radar_layout = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([
                    Constraint::Length(34), // Radar Grid
                    Constraint::Length(28), // Legend Panel
                    Constraint::Min(20),    // Live Pass Telemetry Table
                ].as_ref())
                .split(visuals_chunks[2]);

            // Construct colored span lines for the radar grid
            let mut spans_lines = Vec::new();
            for (y, row) in radar_grid.iter().enumerate() {
                let mut line_spans = Vec::new();
                line_spans.push(Span::raw(" "));
                for (x, &ch) in row.iter().enumerate() {
                    let span = match ch {
                        'N' | 'E' | 'W' => {
                            Span::styled(ch.to_string(), Style::default().fg(Color::White).add_modifier(Modifier::BOLD))
                        }
                        'S' => {
                            if y == 14 && x == 15 {
                                Span::styled(ch.to_string(), Style::default().fg(Color::White).add_modifier(Modifier::BOLD))
                            } else {
                                Span::styled(ch.to_string(), Style::default().fg(Color::Green).add_modifier(Modifier::BOLD))
                            }
                        }
                        '+' | '-' | ':' => {
                            Span::styled(ch.to_string(), Style::default().fg(Color::DarkGray))
                        }
                        '.' => {
                            Span::styled(ch.to_string(), Style::default().fg(Color::DarkGray))
                        }
                        's' => {
                            Span::styled(ch.to_string(), Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD))
                        }
                        '·' => {
                            Span::styled(ch.to_string(), Style::default().fg(Color::Green))
                        }
                        ' ' => {
                            Span::raw(" ")
                        }
                        ch_val => {
                            Span::raw(ch_val.to_string())
                        }
                    };
                    line_spans.push(span);
                }
                spans_lines.push(Line::from(line_spans));
            }

            let map_title = if geo_stats.4 {
                " Radar [S:Active|s:Pass|·:Trail] (SOLVED) ".to_string()
            } else {
                " Radar [S:Active|s:Pass|·:Trail] (UNSOLVED) ".to_string()
            };

            let map_panel = Paragraph::new(spans_lines)
                .block(Block::default().borders(Borders::ALL).title(map_title));

            // Build the live pass telemetry table
            let mut active_row = None;
            let mut candidate_rows = Vec::new();
            {
                let lock = get_visible_sats().lock().unwrap();
                for (idx, sat) in lock.iter().enumerate() {
                    let is_active = Some(idx) == closest_idx;
                    let doppler_hz = sat.freq_expected - current_freq;

                    let (name_style, row_style, status_text, status_style) = if is_active {
                        let name_st = Style::default().fg(Color::Green).add_modifier(Modifier::BOLD);

                        let (row_st, stat_txt, stat_st) = if state_name == "LOCKED" {
                            if let Some(dur) = captured_duration {
                                if dur >= 60.0 {
                                    (
                                        Style::default().fg(Color::LightGreen).add_modifier(Modifier::BOLD),
                                        format!("LOCKED (DONE: {:.0}s)", dur),
                                        Style::default().fg(Color::LightGreen).add_modifier(Modifier::BOLD),
                                    )
                                } else {
                                    (
                                        Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
                                        format!("LOCKED (CAP: {:.0}s)", dur),
                                        Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
                                    )
                                }
                            } else {
                                (
                                    Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
                                    format!("LOCKED ({:.1} dB)", snr_db),
                                    Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
                                )
                            }
                        } else {
                            (
                                Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
                                "SEARCHING".to_string(),
                                Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
                            )
                        };

                        (name_st, row_st, stat_txt, stat_st)
                    } else {
                        let name_st = Style::default().fg(Color::Yellow);
                        let row_st = Style::default().fg(Color::Gray);
                        let stat_txt = format!("In View ({:.0}%)", sat.pass_progress * 100.0);
                        let stat_st = Style::default().fg(Color::Gray);

                        (name_st, row_st, stat_txt, stat_st)
                    };

                    let row = Row::new(vec![
                        Cell::from(sat.name.clone()).style(name_style),
                        Cell::from(format!("{:.1}°", sat.az.to_degrees())),
                        Cell::from(format!("{:.1}°", sat.el.to_degrees())),
                        Cell::from(format!("{:.1} km", sat.range / 1000.0)),
                        Cell::from(format!("{:+.2} kHz", doppler_hz / 1000.0)),
                        Cell::from(status_text).style(status_style),
                    ]).style(row_style);

                    if is_active {
                        active_row = Some(row);
                    } else {
                        candidate_rows.push(row);
                    }
                }
            }

            let mut sat_rows = Vec::new();
            if let Some(row) = active_row {
                sat_rows.push(row);
            }
            sat_rows.extend(candidate_rows);

            let total_sats = sat_rows.len();
            let max_visible_rows = (radar_layout[2].height as usize).saturating_sub(3).max(1);

            let scroll_offset = state.telemetry_scroll_offset.min(total_sats.saturating_sub(max_visible_rows));
            state.telemetry_scroll_offset = scroll_offset;

            let showing_str = if total_sats > max_visible_rows {
                format!(" (Showing {}-{}/{})", scroll_offset + 1, (scroll_offset + max_visible_rows).min(total_sats), total_sats)
            } else {
                "".to_string()
            };

            let mut visible_rows = sat_rows.iter().skip(scroll_offset).take(max_visible_rows).cloned().collect::<Vec<_>>();

            if visible_rows.is_empty() {
                visible_rows.push(Row::new(vec![
                    Cell::from("No satellites in range").style(Style::default().fg(Color::DarkGray)),
                    Cell::from(""),
                    Cell::from(""),
                    Cell::from(""),
                    Cell::from(""),
                    Cell::from(""),
                ]));
            }

            let sat_header = Row::new(vec![
                Cell::from("Satellite").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Azimuth").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Elevation").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Range").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Doppler").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
                Cell::from("Status").style(Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan)),
            ]);

            let sat_table_title = format!(" Live Pass Telemetry{} | {} ", showing_str, next_pass_info);

            let sat_table = Table::new(
                visible_rows,
                [
                    Constraint::Percentage(24),
                    Constraint::Percentage(13),
                    Constraint::Percentage(13),
                    Constraint::Percentage(15),
                    Constraint::Percentage(15),
                    Constraint::Percentage(20),
                ]
            )
            .header(sat_header)
            .block(Block::default().borders(Borders::ALL).title(sat_table_title));

            let legend_lines = vec![
                Line::from(vec![
                    Span::styled(" S ", Style::default().fg(Color::Green).add_modifier(Modifier::BOLD)),
                    Span::raw("Active Tracked Sat"),
                ]),
                Line::from(vec![
                    Span::styled(" s ", Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)),
                    Span::raw("Candidate (Visible)"),
                ]),
                Line::from(vec![
                    Span::styled(" · ", Style::default().fg(Color::Green)),
                    Span::raw("Active Sat Trail"),
                ]),
                Line::from(vec![
                    Span::styled(" + ", Style::default().fg(Color::DarkGray)),
                    Span::raw("Zenith (Overhead)"),
                ]),
                Line::from(vec![
                    Span::styled(" . ", Style::default().fg(Color::DarkGray)),
                    Span::raw("Elev Rings (30/60)"),
                ]),
                Line::from(vec![
                    Span::styled(" ─/│", Style::default().fg(Color::DarkGray)),
                    Span::raw("Compass Axes"),
                ]),
                Line::from(""),
                Line::from(Span::styled(" * s vs S explanation:", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))),
                Line::from(" S is the actively tracked"),
                Line::from(" satellite. s represents"),
                Line::from(" other candidates in sky"),
                Line::from(" but not locked on."),
            ];

            let legend_panel = Paragraph::new(legend_lines)
                .block(Block::default().borders(Borders::ALL).title(" Radar Legend "));

            rect.render_widget(map_panel, radar_layout[0]);
            rect.render_widget(legend_panel, radar_layout[1]);
            rect.render_widget(sat_table, radar_layout[2]);

            if total_sats > max_visible_rows {
                let mut scrollbar_state = ScrollbarState::default()
                    .content_length(total_sats)
                    .position(scroll_offset)
                    .viewport_content_length(max_visible_rows);
                rect.render_stateful_widget(
                    Scrollbar::default()
                        .orientation(ScrollbarOrientation::VerticalRight)
                        .begin_symbol(Some("▲"))
                        .end_symbol(Some("▼"))
                        .track_symbol(Some("░"))
                        .thumb_symbol("█"),
                    radar_layout[2].inner(&ratatui::layout::Margin { vertical: 1, horizontal: 0 }),
                    &mut scrollbar_state,
                );
            }

            let log_entries = {
                let max_lines = (visuals_chunks[3].height as usize).saturating_sub(2).max(1);
                let start_idx = self.event_logs.len().saturating_sub(max_lines);
                self.event_logs.iter().skip(start_idx).cloned().collect::<Vec<String>>().join("\n")
            };
            let logs_panel = Paragraph::new(log_entries)
                .block(Block::default().borders(Borders::ALL).title(" Console Event Logs "))
                .style(Style::default().fg(Color::Gray));
            rect.render_widget(logs_panel, visuals_chunks[3]);

            let help_menu = Paragraph::new(
                " [q/Esc] Quit | [1-5] Switch Profile | [a] Toggle AGC | [n] Toggle Notching \n [Up/Down] Adj LNA (+/-8dB) | [Left/Right] Adj VGA (+/-2dB) | [g] Toggle AMP (0/14dB) | [j/k] Scroll Tel"
            )
            .block(Block::default().borders(Borders::ALL).title(" Shortcuts & Keybindings "));
            rect.render_widget(help_menu, chunks[2]);
        });
    }
}

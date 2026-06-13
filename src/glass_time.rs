use std::net::{TcpListener, TcpStream};
use std::sync::{Arc, Mutex};
use std::thread;
use std::io::{Read, Write};
use crossbeam_channel::{Sender, unbounded};
use serde::Serialize;

/// SHA-1 implementation to calculate the WebSocket Accept Key without external dependencies
fn sha1(input: &str) -> [u8; 20] {
    let mut h0: u32 = 0x67452301;
    let mut h1: u32 = 0xEFCDAB89;
    let mut h2: u32 = 0x98BADCFE;
    let mut h3: u32 = 0x10325476;
    let mut h4: u32 = 0xC3D2E1F0;

    let bytes = input.as_bytes();
    let bit_len = (bytes.len() as u64) * 8;
    
    let mut padded = bytes.to_vec();
    padded.push(0x80);
    while (padded.len() + 8) % 64 != 0 {
        padded.push(0x00);
    }
    
    for shift in (0..8).rev() {
        padded.push(((bit_len >> (shift * 8)) & 0xFF) as u8);
    }

    for chunk in padded.chunks_exact(64) {
        let mut w = [0u32; 80];
        for i in 0..16 {
            w[i] = ((chunk[i * 4] as u32) << 24)
                | ((chunk[i * 4 + 1] as u32) << 16)
                | ((chunk[i * 4 + 2] as u32) << 8)
                | (chunk[i * 4 + 3] as u32);
        }
        for i in 16..80 {
            w[i] = (w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16]).rotate_left(1);
        }

        let mut a = h0;
        let mut b = h1;
        let mut c = h2;
        let mut d = h3;
        let mut e = h4;

        for i in 0..80 {
            let (f, k) = if i < 20 {
                ((b & c) | (!b & d), 0x5A827999)
            } else if i < 40 {
                (b ^ c ^ d, 0x6ED9EBA1)
            } else if i < 60 {
                ((b & c) | (b & d) | (c & d), 0x8F1BBCDC)
            } else {
                (b ^ c ^ d, 0xCA62C1D6)
            };

            let temp = a.rotate_left(5)
                .wrapping_add(f)
                .wrapping_add(e)
                .wrapping_add(k)
                .wrapping_add(w[i]);
            e = d;
            d = c;
            c = b.rotate_left(30);
            b = a;
            a = temp;
        }

        h0 = h0.wrapping_add(a);
        h1 = h1.wrapping_add(b);
        h2 = h2.wrapping_add(c);
        h3 = h3.wrapping_add(d);
        h4 = h4.wrapping_add(e);
    }

    let mut result = [0u8; 20];
    for (i, &val) in [h0, h1, h2, h3, h4].iter().enumerate() {
        result[i * 4] = ((val >> 24) & 0xFF) as u8;
        result[i * 4 + 1] = ((val >> 16) & 0xFF) as u8;
        result[i * 4 + 2] = ((val >> 8) & 0xFF) as u8;
        result[i * 4 + 3] = (val & 0xFF) as u8;
    }
    result
}

/// Base64 encoder to format the WebSocket Accept Key
fn base64_encode(input: &[u8]) -> String {
    const CHARSET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut result = String::with_capacity((input.len() + 2) / 3 * 4);
    for chunk in input.chunks(3) {
        match chunk.len() {
            3 => {
                let b0 = chunk[0] as usize;
                let b1 = chunk[1] as usize;
                let b2 = chunk[2] as usize;
                result.push(CHARSET[b0 >> 2] as char);
                result.push(CHARSET[((b0 & 0x03) << 4) | (b1 >> 4)] as char);
                result.push(CHARSET[((b1 & 0x0F) << 2) | (b2 >> 6)] as char);
                result.push(CHARSET[b2 & 0x3F] as char);
            }
            2 => {
                let b0 = chunk[0] as usize;
                let b1 = chunk[1] as usize;
                result.push(CHARSET[b0 >> 2] as char);
                result.push(CHARSET[((b0 & 0x03) << 4) | (b1 >> 4)] as char);
                result.push(CHARSET[(b1 & 0x0F) << 2] as char);
                result.push('=');
            }
            1 => {
                let b0 = chunk[0] as usize;
                result.push(CHARSET[b0 >> 2] as char);
                result.push(CHARSET[(b0 & 0x03) << 4] as char);
                result.push('=');
                result.push('=');
            }
            _ => unreachable!(),
        }
    }
    result
}

/// Encodes a WebSocket text frame from a payload string
fn encode_ws_frame(payload: &str) -> Vec<u8> {
    let bytes = payload.as_bytes();
    let len = bytes.len();
    let mut frame = Vec::new();
    
    // Fin=1, Opcode=1 (Text)
    frame.push(0x81);
    
    if len <= 125 {
        frame.push(len as u8);
    } else if len <= 65535 {
        frame.push(126);
        frame.push(((len >> 8) & 0xFF) as u8);
        frame.push((len & 0xFF) as u8);
    } else {
        frame.push(127);
        for shift in (0..8).rev() {
            frame.push(((len >> (shift * 8)) & 0xFF) as u8);
        }
    }
    
    frame.extend_from_slice(bytes);
    frame
}

/// Handles the WebSocket handshake with a TCP client
fn handle_handshake(stream: &mut TcpStream) -> Result<(), std::io::Error> {
    let mut buf = [0u8; 2048];
    let n = stream.read(&mut buf)?;
    let request = String::from_utf8_lossy(&buf[..n]);
    
    let key_header = "Sec-WebSocket-Key: ";
    if let Some(pos) = request.find(key_header) {
        let start = pos + key_header.len();
        if let Some(end) = request[start..].find("\r\n") {
            let key = request[start..start+end].trim();
            let accept_val = base64_encode(&sha1(&format!("{}{}", key, "258EAFA5-E914-47DA-95CA-C5AB0DC85B11")));
            
            let response = format!(
                "HTTP/1.1 101 Switching Protocols\r\n\
                 Upgrade: websocket\r\n\
                 Connection: Upgrade\r\n\
                 Sec-WebSocket-Accept: {}\r\n\r\n",
                accept_val
            );
            stream.write_all(response.as_bytes())?;
            return Ok(());
        }
    }
    Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid WebSocket Handshake"))
}

#[derive(Serialize, Clone, Debug)]
pub struct TelemetryFrame {
    pub timestamp: String,
    pub rx_position: [f64; 3],
    pub rx_velocity: [f64; 3],
    pub clock_bias: f64,
    pub clock_drift: f64,
    pub active_channels: Vec<ChannelTelem>,
}

#[derive(Serialize, Clone, Debug)]
pub struct ChannelTelem {
    pub id: usize,
    pub sat_name: String,
    pub status: String,
    pub target_freq: f64,
    pub freq_offset: f64,
    pub snr: f64,
    pub tec: f64,
    pub s4: f64,
    pub sigma_phi: f64,
    pub is_dual: bool,
    pub sat_position: [f64; 3],
}

pub struct GlassTimeServer {
    pub telemetry_tx: Sender<TelemetryFrame>,
}

impl GlassTimeServer {
    pub fn start(port: u16) -> Self {
        let (telemetry_tx, telemetry_rx) = unbounded::<TelemetryFrame>();
        let clients = Arc::new(Mutex::new(Vec::new()));
        let clients_clone = clients.clone();

        // Spawn a thread to accept incoming TCP WebSocket connections
        thread::spawn(move || {
            let listener = match TcpListener::bind(format!("0.0.0.0:{}", port)) {
                Ok(l) => l,
                Err(e) => {
                    tracing::warn!("Failed to bind WebSocket server on port {}: {}", port, e);
                    return;
                }
            };
            listener.set_nonblocking(false).ok();

            for stream in listener.incoming() {
                if let Ok(mut stream) = stream {
                    let clients_list = clients_clone.clone();
                    thread::spawn(move || {
                        if handle_handshake(&mut stream).is_ok() {
                            stream.set_nonblocking(true).ok();
                            if let Ok(mut list) = clients_list.lock() {
                                list.push(stream);
                            }
                        }
                    });
                }
            }
        });

        // Spawn a thread to broadcast serialized telemetry frames to all clients
        thread::spawn(move || {
            for frame in telemetry_rx {
                if let Ok(payload) = serde_json::to_string(&frame) {
                    let frame_bytes = encode_ws_frame(&payload);
                    if let Ok(mut list) = clients.lock() {
                        let mut active_clients = Vec::new();
                        for mut client in list.drain(..) {
                            if client.write_all(&frame_bytes).is_ok() {
                                active_clients.push(client);
                            }
                        }
                        *list = active_clients;
                    }
                }
            }
        });

        Self { telemetry_tx }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sha1_and_base64() {
        let hash = sha1("dGhlIHNhbXBsZSBub25jZQ==258EAFA5-E914-47DA-95CA-C5AB0DC85B11");
        let b64 = base64_encode(&hash);
        assert_eq!(b64, "s3pPLMBiTxaQ9kYGzzhZRbK+xOo=");
    }

    #[test]
    fn test_websocket_server_e2e() {
        let server = GlassTimeServer::start(19876);
        std::thread::sleep(std::time::Duration::from_millis(100));

        let mut client = TcpStream::connect("127.0.0.1:19876").expect("Failed to connect client");

        let handshake = "GET / HTTP/1.1\r\n\
                         Host: 127.0.0.1:19876\r\n\
                         Upgrade: websocket\r\n\
                         Connection: Upgrade\r\n\
                         Sec-WebSocket-Key: dGhlIHNhbXBsZSBub25jZQ==\r\n\
                         Sec-WebSocket-Version: 13\r\n\r\n";
        client.write_all(handshake.as_bytes()).expect("Failed to send handshake");

        let mut buf = [0u8; 1024];
        let n = client.read(&mut buf).expect("Failed to read response");
        let resp = String::from_utf8_lossy(&buf[..n]);
        assert!(resp.contains("HTTP/1.1 101 Switching Protocols"));
        assert!(resp.contains("Sec-WebSocket-Accept: s3pPLMBiTxaQ9kYGzzhZRbK+xOo="));

        let frame = TelemetryFrame {
            timestamp: "2026-06-13T00:00:00Z".to_string(),
            rx_position: [1.0, 2.0, 3.0],
            rx_velocity: [4.0, 5.0, 6.0],
            clock_bias: 1e-6,
            clock_drift: 1.2,
            active_channels: vec![],
        };
        server.telemetry_tx.send(frame).expect("Failed to send telemetry frame");

        std::thread::sleep(std::time::Duration::from_millis(100));
        let n2 = client.read(&mut buf).expect("Failed to read frame");
        
        assert_eq!(buf[0], 0x81);
        let len_byte = buf[1];
        assert!(len_byte > 0);
        
        let payload_start = if len_byte <= 125 {
            2
        } else if len_byte == 126 {
            4
        } else {
            10
        };
        
        let payload = String::from_utf8_lossy(&buf[payload_start..n2]);
        assert!(payload.contains("2026-06-13T00:00:00Z"));
        assert!(payload.contains("clock_drift"));
    }
}

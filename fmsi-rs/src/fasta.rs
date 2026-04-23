use std::fs::File;
use std::io::{self, Cursor, Read};
use std::process::{Command, Stdio};

#[derive(Clone, Debug)]
pub struct SeqRecord {
    pub name: String,
    pub seq: Vec<u8>,
}

pub fn read_masked_superstring(path: &str) -> Result<String, String> {
    let records = read_records(path)?;
    if records.is_empty() {
        return Err("the FASTA/FASTQ file is empty".to_string());
    }
    if records.len() > 1 {
        eprintln!("Warning: the fasta file contains more than one entry. Only the first entry will be used.");
    }
    String::from_utf8(records[0].seq.clone()).map_err(|e| e.to_string())
}

pub fn read_records(path: &str) -> Result<Vec<SeqRecord>, String> {
    let mut reader = open_reader(path)?;
    let mut data = String::new();
    reader.read_to_string(&mut data).map_err(|e| e.to_string())?;
    parse_fasta_fastq(&data)
}

fn open_reader(path: &str) -> Result<Cursor<Vec<u8>>, String> {
    if path == "-" {
        let mut buf = Vec::new();
        io::stdin().read_to_end(&mut buf).map_err(|e| e.to_string())?;
        return Ok(Cursor::new(buf));
    }
    if path.ends_with(".gz") {
        let out = Command::new("gzip")
            .arg("-cd")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|e| format!("failed to spawn gzip -cd for {path}: {e}"))?
            .wait_with_output()
            .map_err(|e| e.to_string())?;
        if !out.status.success() {
            return Err(format!("gzip -cd failed for {path}"));
        }
        return Ok(Cursor::new(out.stdout));
    }
    let mut file = File::open(path).map_err(|e| e.to_string())?;
    let mut buf = Vec::new();
    file.read_to_end(&mut buf).map_err(|e| e.to_string())?;
    Ok(Cursor::new(buf))
}

fn parse_fasta_fastq(data: &str) -> Result<Vec<SeqRecord>, String> {
    let bytes = data.as_bytes();
    let mut i = 0usize;
    let mut records = Vec::new();
    while i < bytes.len() {
        while i < bytes.len() && (bytes[i] == b'\n' || bytes[i] == b'\r') {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        match bytes[i] {
            b'>' => {
                i += 1;
                let name_start = i;
                while i < bytes.len() && bytes[i] != b'\n' && bytes[i] != b'\r' {
                    i += 1;
                }
                let name = data[name_start..i].trim().to_string();
                while i < bytes.len() && (bytes[i] == b'\n' || bytes[i] == b'\r') {
                    i += 1;
                }
                let mut seq = Vec::new();
                while i < bytes.len() && bytes[i] != b'>' && bytes[i] != b'@' {
                    if bytes[i] != b'\n' && bytes[i] != b'\r' {
                        seq.push(bytes[i]);
                    }
                    i += 1;
                }
                records.push(SeqRecord { name, seq });
            }
            b'@' => {
                i += 1;
                let name_start = i;
                while i < bytes.len() && bytes[i] != b'\n' && bytes[i] != b'\r' {
                    i += 1;
                }
                let name = data[name_start..i].trim().to_string();
                while i < bytes.len() && (bytes[i] == b'\n' || bytes[i] == b'\r') {
                    i += 1;
                }
                let mut seq = Vec::new();
                while i < bytes.len() && bytes[i] != b'+' {
                    if bytes[i] != b'\n' && bytes[i] != b'\r' {
                        seq.push(bytes[i]);
                    }
                    i += 1;
                }
                if i >= bytes.len() || bytes[i] != b'+' {
                    return Err("invalid FASTQ: missing '+' separator".to_string());
                }
                while i < bytes.len() && bytes[i] != b'\n' && bytes[i] != b'\r' {
                    i += 1;
                }
                while i < bytes.len() && (bytes[i] == b'\n' || bytes[i] == b'\r') {
                    i += 1;
                }
                let mut qual_consumed = 0usize;
                while i < bytes.len() && qual_consumed < seq.len() {
                    if bytes[i] != b'\n' && bytes[i] != b'\r' {
                        qual_consumed += 1;
                    }
                    i += 1;
                }
                records.push(SeqRecord { name, seq });
            }
            other => {
                return Err(format!("unsupported input format starting with byte '{}'", other as char));
            }
        }
    }
    Ok(records)
}

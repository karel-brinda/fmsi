pub const INVALID_NUCLEOTIDE: u8 = 4;

#[inline]
pub fn nucleotide_to_int(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => INVALID_NUCLEOTIDE,
    }
}

#[inline]
pub fn complementary_nucleotide(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

#[inline]
pub fn is_upper(c: u8) -> bool {
    c.is_ascii_uppercase()
}

pub fn reverse_complement_string(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| complementary_nucleotide(b))
        .collect()
}

pub fn are_strings_equal(a: &[u8], b: &[u8]) -> bool {
    a == b
}

pub fn infer_k(ms: &str) -> Result<usize, String> {
    let bytes = ms.as_bytes();
    if bytes.is_empty() {
        return Err("masked superstring is empty".to_string());
    }
    let mut k = 1usize;
    while k <= bytes.len() && !is_upper(bytes[bytes.len() - k]) {
        k += 1;
    }
    if k > bytes.len() {
        return Err("failed to infer k from trailing lowercase mask convention".to_string());
    }
    Ok(k)
}

pub fn next_invalid_character_or_end(sequence: &[u8]) -> usize {
    for (i, &c) in sequence.iter().enumerate() {
        if nucleotide_to_int(c) == INVALID_NUCLEOTIDE {
            return i;
        }
    }
    sequence.len()
}

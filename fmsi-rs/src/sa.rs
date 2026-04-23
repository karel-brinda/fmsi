/// Suffix-array backend.
///
/// This is intentionally split so the first faithful Rust translation can keep the
/// MBWT/FMS logic unchanged while swapping in a faster low-level constructor later.
///
/// For now the default backend is a naive lexicographic suffix sort. It is correct,
/// self-contained, and suitable for tiny regression cases. Before serious benchmarks,
/// replace this with libsais or divsufsort.

pub fn build_suffix_array_with_sentinel(text: &[u8]) -> Vec<usize> {
    let n = text.len();
    let mut sa: Vec<usize> = (0..=n).collect();
    sa.sort_by(|&a, &b| compare_suffixes(text, a, b));
    sa
}

fn compare_suffixes(text: &[u8], a: usize, b: usize) -> std::cmp::Ordering {
    let n = text.len();
    if a == b {
        return std::cmp::Ordering::Equal;
    }
    let mut i = a;
    let mut j = b;
    loop {
        let ac = if i == n { None } else { Some(text[i]) };
        let bc = if j == n { None } else { Some(text[j]) };
        match (ac, bc) {
            (None, None) => return std::cmp::Ordering::Equal,
            (None, Some(_)) => return std::cmp::Ordering::Less,
            (Some(_), None) => return std::cmp::Ordering::Greater,
            (Some(x), Some(y)) => {
                if x != y {
                    return x.cmp(&y);
                }
                i += 1;
                j += 1;
            }
        }
    }
}

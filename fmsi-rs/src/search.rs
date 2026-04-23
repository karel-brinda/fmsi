use crate::dna::nucleotide_to_int;
use crate::index::FmsIndex;

#[inline]
pub fn rank(index: &FmsIndex, i: usize, c: u8) -> usize {
    let gt_position = index.ac_gt.rank1(i);
    if c >= 2 {
        let t_position = index.gt.rank1(gt_position);
        if c == 2 {
            gt_position - t_position
        } else {
            t_position
        }
    } else {
        let c_position = index.ac.rank1(i - gt_position);
        if c == 0 {
            i - gt_position - c_position - usize::from(i >= index.dollar_position + 1)
        } else {
            c_position
        }
    }
}

#[inline]
pub fn access(index: &FmsIndex, i: usize) -> u8 {
    let gt_position = index.ac_gt.rank1(i);
    if index.ac_gt.get(i) {
        2 + u8::from(index.gt.get(gt_position))
    } else {
        u8::from(index.ac.get(i - gt_position))
    }
}

#[inline]
pub fn update_range(index: &FmsIndex, i: &mut usize, j: &mut usize, c: u8) {
    if *j == *i {
        return;
    }
    let count = index.counts[c as usize];
    *i = count + rank(index, *i, c);
    *j = count + rank(index, *j, c);
}

#[inline]
pub fn extend_range_with_klcp(index: &FmsIndex, i: &mut usize, j: &mut usize) {
    let Some(klcp) = &index.klcp else {
        return;
    };
    while *j > 0 && *j < klcp.len() && klcp.get(*j - 1) {
        *j += 1;
    }
    while *i > 0 && klcp.get(*i - 1) {
        *i -= 1;
    }
}

pub fn get_range_with_pattern(index: &FmsIndex, pattern: &[u8]) -> (usize, usize) {
    let mut sa_start = 0usize;
    let mut sa_end = index.sa_transformed_mask.len();
    for &c in pattern.iter().rev() {
        if sa_start == sa_end {
            break;
        }
        update_range(index, &mut sa_start, &mut sa_end, nucleotide_to_int(c));
    }
    (sa_start, sa_end)
}

pub fn infer_presence_or(index: &FmsIndex, sa_start: usize, sa_end: usize) -> i32 {
    if sa_start == sa_end {
        return -1;
    }
    for i in sa_start..sa_end {
        if index.sa_transformed_mask.get(i) {
            return 1;
        }
    }
    0
}

pub fn infer_presence_all(index: &FmsIndex, sa_start: usize, sa_end: usize) -> i32 {
    if sa_start != sa_end {
        i32::from(index.sa_transformed_mask.get(sa_start))
    } else {
        -1
    }
}

pub fn kmer_order(index: &FmsIndex, sa_start: usize) -> i64 {
    index.sa_transformed_mask.rank1(sa_start) as i64
}

pub fn kmer_order_if_present(index: &FmsIndex, sa_start: usize, sa_end: usize) -> i64 {
    if infer_presence_or(index, sa_start, sa_end) == 1 {
        kmer_order(index, sa_start)
    } else {
        -1
    }
}

pub fn single_query_or(index: &FmsIndex, pattern: &[u8]) -> i32 {
    let (sa_start, sa_end) = get_range_with_pattern(index, pattern);
    infer_presence_or(index, sa_start, sa_end)
}

pub fn single_query_all(index: &FmsIndex, pattern: &[u8]) -> i32 {
    let (sa_start, sa_end) = get_range_with_pattern(index, pattern);
    infer_presence_all(index, sa_start, sa_end)
}

pub fn single_query_order(index: &FmsIndex, pattern: &[u8]) -> i64 {
    let (sa_start, sa_end) = get_range_with_pattern(index, pattern);
    kmer_order_if_present(index, sa_start, sa_end)
}

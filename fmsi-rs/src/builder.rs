use crate::bitvec::BitVecRank;
use crate::dna::{infer_k, is_upper, nucleotide_to_int, INVALID_NUCLEOTIDE};
use crate::index::FmsIndex;
use crate::sa::build_suffix_array_with_sentinel;

pub fn construct(ms: &str, k: usize, use_klcp: bool) -> Result<FmsIndex, String> {
    if k == 0 || k > 32 {
        return Err("this Rust v1 keeps the agreed scope k <= 32".to_string());
    }
    let text = convert_superstring(ms)?;
    let sa = build_suffix_array_with_sentinel(&text);
    let n = text.len();

    let klcp = if use_klcp {
        Some(construct_klcp(&sa, &text, k - 1))
    } else {
        None
    };

    let mut dollar_position = 0usize;
    let mut bwt = vec![0u8; n + 1];
    let mut mask = BitVecRank::with_len(n + 1);
    let ms_bytes = ms.as_bytes();
    for (i, &sa_i) in sa.iter().enumerate() {
        if sa_i == 0 {
            dollar_position = i;
            bwt[i] = 0;
        } else {
            bwt[i] = text[sa_i - 1];
        }
        if sa_i != n {
            mask.set(i, is_upper(ms_bytes[sa_i]));
        }
    }
    mask.rebuild_rank();

    let mut ac_gt = BitVecRank::with_len(bwt.len());
    let mut gt_count = 0usize;
    for (i, &x) in bwt.iter().enumerate() {
        let is_gt = x >= 2;
        if is_gt {
            gt_count += 1;
        }
        ac_gt.set(i, is_gt);
    }
    ac_gt.rebuild_rank();

    let ac_count = bwt.len() - gt_count;
    let mut ac = BitVecRank::with_len(ac_count);
    let mut gt = BitVecRank::with_len(gt_count);
    let mut a_count = 0usize;
    let mut g_count = 0usize;
    let mut ac_idx = 0usize;
    let mut gt_idx = 0usize;
    for (i, &x) in bwt.iter().enumerate() {
        let is_one = (x & 1) != 0;
        if !ac_gt.get(i) {
            ac.set(ac_idx, is_one);
            if !is_one {
                a_count += 1;
            }
            ac_idx += 1;
        } else {
            gt.set(gt_idx, is_one);
            if !is_one {
                g_count += 1;
            }
            gt_idx += 1;
        }
    }
    ac.rebuild_rank();
    gt.rebuild_rank();

    Ok(FmsIndex {
        ac_gt,
        ac,
        gt,
        sa_transformed_mask: mask,
        counts: [1, a_count, ac_count, ac_count + g_count],
        dollar_position,
        klcp,
        k,
        predictor: Default::default(),
    })
}

pub fn infer_or_validate_k(ms: &str, k: Option<usize>) -> Result<usize, String> {
    let inferred = infer_k(ms)?;
    match k {
        Some(user_k) => Ok(user_k),
        None => Ok(inferred),
    }
}

fn convert_superstring(ms: &str) -> Result<Vec<u8>, String> {
    let mut ret = Vec::with_capacity(ms.len());
    for &c in ms.as_bytes() {
        let x = nucleotide_to_int(c);
        if x == INVALID_NUCLEOTIDE {
            return Err(format!("invalid nucleotide in masked superstring: '{}'", c as char));
        }
        ret.push(x);
    }
    Ok(ret)
}

fn construct_klcp(sa: &[usize], text: &[u8], k_minus_1: usize) -> BitVecRank {
    let n = text.len();
    let mut klcp = BitVecRank::with_len(n + 1);
    if k_minus_1 == 0 {
        klcp.rebuild_rank();
        return klcp;
    }
    for i in 0..n {
        let sa_i = sa[i];
        let sa_i1 = sa[i + 1];
        if sa_i > n - k_minus_1 || sa_i1 > n - k_minus_1 {
            continue;
        }
        if text[sa_i..sa_i + k_minus_1] == text[sa_i1..sa_i1 + k_minus_1] {
            klcp.set(i, true);
        }
    }
    klcp.rebuild_rank();
    klcp
}

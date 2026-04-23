use crate::dna::{are_strings_equal, next_invalid_character_or_end, reverse_complement_string};
use crate::index::FmsIndex;
use crate::search::{
    extend_range_with_klcp, get_range_with_pattern, infer_presence_all, infer_presence_or,
    kmer_order_if_present, single_query_all, single_query_or, single_query_order, update_range,
};

#[derive(Copy, Clone, Eq, PartialEq)]
pub enum QueryMode {
    Or,
    All,
}

pub fn query_record(
    index: &mut FmsIndex,
    sequence: &[u8],
    k: usize,
    use_klcp: bool,
    mode: QueryMode,
    output_orders: bool,
) -> String {
    let mut out = String::new();
    if sequence.len() < k {
        return out;
    }

    let mut offset = 0usize;
    let mut output_comma = false;
    let mut remaining = sequence;
    let max_chunk_len = {
        let base = (sequence.len() as f64).sqrt() as usize;
        k + 10.max(400.min(2 * base))
    };

    while !remaining.is_empty() {
        let current_length = next_invalid_character_or_end(remaining);
        let mut local = current_length;
        while local >= k {
            if output_orders && output_comma {
                out.push(',');
            }
            output_comma = true;
            let chunk_length = local.min(max_chunk_len);
            let chunk = &remaining[..chunk_length];
            let chunk_string = if use_klcp {
                query_kmers_streaming(index, chunk, k, mode, output_orders)
            } else {
                query_kmers_single(index, chunk, k, mode, output_orders)
            };
            out.push_str(&chunk_string);
            let consumed = chunk_length - k + 1;
            remaining = &remaining[consumed..];
            local -= consumed;
            offset += consumed;
        }

        if remaining.len() <= current_length {
            break;
        }
        let zeros = k.min(current_length + 1);
        for _ in 0..zeros {
            if output_orders {
                if output_comma {
                    out.push(',');
                }
                output_comma = true;
                out.push_str("-1");
            } else {
                out.push('0');
            }
        }
        remaining = &remaining[current_length + 1..];
        offset += current_length + 1;
    }

    let _ = offset;
    out
}

fn query_kmers_streaming(
    index: &mut FmsIndex,
    sequence: &[u8],
    k: usize,
    mode: QueryMode,
    output_orders: bool,
) -> String {
    let mut result = vec![-1i64; sequence.len() - k + 1];
    let rc_storage = reverse_complement_string(sequence);
    let mut forward_seq = sequence;
    let mut reverse_seq = rc_storage.as_slice();

    let should_swap = index.predictor.predict_swap();
    if should_swap {
        std::mem::swap(&mut forward_seq, &mut reverse_seq);
    }

    let mut forward_predictor_result = 0i32;
    let mut backward_predictor_result = 0i32;

    let mut sa_start = usize::MAX;
    let mut sa_end = usize::MAX;
    for i in 0..=sequence.len() - k {
        let i_back = sequence.len() - k - i;
        if sa_start == sa_end {
            (sa_start, sa_end) = get_range_with_pattern(index, &forward_seq[i_back..i_back + k]);
        } else {
            extend_range_with_klcp(index, &mut sa_start, &mut sa_end);
            update_range(
                index,
                &mut sa_start,
                &mut sa_end,
                crate::dna::nucleotide_to_int(forward_seq[i_back]),
            );
        }
        if output_orders {
            result[i_back] = kmer_order_if_present(index, sa_start, sa_end);
            if result[i_back] >= 0 {
                forward_predictor_result += 1;
            } else {
                forward_predictor_result -= 1;
            }
        } else {
            result[i_back] = match mode {
                QueryMode::Or => infer_presence_or(index, sa_start, sa_end) as i64,
                QueryMode::All => infer_presence_all(index, sa_start, sa_end) as i64,
            };
            forward_predictor_result += result[i_back] as i32;
        }
    }

    sa_start = usize::MAX;
    sa_end = usize::MAX;
    for i in 0..=sequence.len() - k {
        if (output_orders && result[i] >= 0)
            || (!output_orders && result[i] == 1)
            || (!output_orders && mode == QueryMode::All && result[i] == 0)
        {
            sa_start = usize::MAX;
            sa_end = usize::MAX;
            continue;
        }
        let i_back = sequence.len() - k - i;
        if sa_start == sa_end {
            (sa_start, sa_end) = get_range_with_pattern(index, &reverse_seq[i_back..i_back + k]);
        } else {
            extend_range_with_klcp(index, &mut sa_start, &mut sa_end);
            update_range(
                index,
                &mut sa_start,
                &mut sa_end,
                crate::dna::nucleotide_to_int(reverse_seq[i_back]),
            );
        }
        let res = if output_orders {
            let x = kmer_order_if_present(index, sa_start, sa_end);
            if x >= 0 {
                backward_predictor_result += 1;
            } else {
                backward_predictor_result -= 1;
            }
            x
        } else {
            let x = match mode {
                QueryMode::Or => infer_presence_or(index, sa_start, sa_end) as i64,
                QueryMode::All => infer_presence_all(index, sa_start, sa_end) as i64,
            };
            backward_predictor_result += x as i32;
            x
        };
        result[i] = result[i].max(res);
    }

    if should_swap {
        result.reverse();
        std::mem::swap(&mut forward_predictor_result, &mut backward_predictor_result);
    }
    index
        .predictor
        .log_result(forward_predictor_result, backward_predictor_result);

    render_result(&result, output_orders)
}

fn query_kmers_single(
    index: &mut FmsIndex,
    sequence: &[u8],
    k: usize,
    mode: QueryMode,
    output_orders: bool,
) -> String {
    let rc_sequence = reverse_complement_string(sequence);
    let mut rendered = String::new();
    for i in 0..=sequence.len() - k {
        let mut kmer = &sequence[i..i + k];
        let rc_start = sequence.len() - k - i;
        let mut rc_kmer = &rc_sequence[rc_start..rc_start + k];
        let should_swap = index.predictor.predict_swap();
        if should_swap {
            std::mem::swap(&mut kmer, &mut rc_kmer);
        }
        let mut forward_predictor_result = 0i32;
        let mut backward_predictor_result = 0i32;

        let got = if output_orders {
            let mut got = single_query_order(index, kmer);
            forward_predictor_result = if got >= 0 { 1 } else { -1 };
            if got < 0 {
                got = single_query_order(index, rc_kmer);
                backward_predictor_result = if got >= 0 { 1 } else { -1 };
            }
            got
        } else {
            let mut got = match mode {
                QueryMode::Or => single_query_or(index, kmer) as i64,
                QueryMode::All => single_query_all(index, kmer) as i64,
            };
            forward_predictor_result = got as i32;
            match mode {
                QueryMode::Or => {
                    if got != 1 {
                        got = single_query_or(index, rc_kmer) as i64;
                        backward_predictor_result = got as i32;
                    }
                }
                QueryMode::All => {
                    if got == -1 {
                        got = single_query_all(index, rc_kmer) as i64;
                        backward_predictor_result = got as i32;
                    }
                }
            }
            got
        };

        if should_swap {
            std::mem::swap(&mut forward_predictor_result, &mut backward_predictor_result);
        }
        index
            .predictor
            .log_result(forward_predictor_result, backward_predictor_result);

        if output_orders {
            if i > 0 {
                rendered.push(',');
            }
            rendered.push_str(&got.to_string());
        } else if got == 1 {
            rendered.push('1');
        } else {
            rendered.push('0');
        }

        let _ = are_strings_equal(kmer, rc_kmer);
    }
    rendered
}

fn render_result(result: &[i64], output_orders: bool) -> String {
    if output_orders {
        result
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    } else {
        result
            .iter()
            .map(|&x| if x == 1 { '1' } else { '0' })
            .collect()
    }
}

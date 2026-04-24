use std::collections::HashSet;

use crate::models::GuideInfo;
use crate::sequence::reverse_complement;

#[derive(Debug, Clone)]
pub struct GuideMatchResult {
    pub sg_rna_sequences: Vec<String>,
    pub sg_rna_intervals: Vec<(usize, usize)>,
    pub sg_rna_cut_points: Vec<usize>,
    pub sg_rna_mismatches: Vec<Vec<usize>>,
    pub sg_rna_names: Vec<String>,
    pub sg_rna_include_idxs: Vec<Vec<usize>>,
    pub include_idxs: Vec<usize>,
    pub exclude_idxs: Vec<usize>,
}

pub fn get_amplicon_info_for_guides(
    ref_seq: &str,
    guides: &[GuideInfo],
    quant_window_coordinates: Option<&str>,
    exclude_bp_from_left: usize,
    exclude_bp_from_right: usize,
    plot_window_size: usize,
) -> GuideMatchResult {
    let ref_seq_upper = ref_seq.to_uppercase();
    let ref_seq_length = ref_seq_upper.len();

    let mut sg_rna_sequences = Vec::new();
    let mut sg_rna_intervals = Vec::new();
    let mut sg_rna_cut_points = Vec::new();
    let mut sg_rna_mismatches = Vec::new();
    let mut sg_rna_names = Vec::new();
    let mut sg_rna_include_idxs = Vec::new();
    let mut all_include_idxs = Vec::new();
    let mut exclude_idxs = Vec::new();

    exclude_idxs.extend(0..exclude_bp_from_left.min(ref_seq_length));
    if exclude_bp_from_right > 0 && exclude_bp_from_right <= ref_seq_length {
        exclude_idxs.extend((ref_seq_length - exclude_bp_from_right)..ref_seq_length);
    }

    let _window_around_cut = plot_window_size.max(1);
    let mut seen_names = HashSet::new();

    for guide in guides {
        if guide.sequence.is_empty() {
            continue;
        }

        let guide_upper = guide.sequence.to_uppercase();
        let guide_rc = reverse_complement(&guide_upper);
        let offset_fw = guide.qw_center + guide_upper.len() as i32 - 1;
        let offset_rc = -guide.qw_center - 1;

        let fw_matches = find_all_matches(&ref_seq_upper, &guide_upper);
        let rv_matches = find_all_matches(&ref_seq_upper, &guide_rc);
        let match_count = fw_matches.len() + rv_matches.len();

        for start in fw_matches {
            let cut_p = start as i32 + offset_fw;
            if cut_p < 0 || cut_p as usize >= ref_seq_length {
                continue;
            }
            let cut_p = cut_p as usize;
            sg_rna_cut_points.push(cut_p);
            sg_rna_intervals.push((start, start + guide_upper.len() - 1));
            sg_rna_mismatches.push(guide.mismatches.clone());
            sg_rna_include_idxs.push(window_from_cut(
                cut_p,
                guide.qw_size,
                ref_seq_length,
                &mut all_include_idxs,
            ));
            sg_rna_names.push(resolve_guide_name(
                &guide.name,
                &guide_upper,
                match_count,
                start,
                &mut seen_names,
            ));
            sg_rna_sequences.push(guide_upper.clone());
        }

        for start in rv_matches {
            let cut_p = start as i32 + offset_rc;
            if cut_p < 0 || cut_p as usize >= ref_seq_length {
                continue;
            }
            let cut_p = cut_p as usize;
            sg_rna_cut_points.push(cut_p);
            sg_rna_intervals.push((start, start + guide_upper.len() - 1));
            sg_rna_mismatches.push(
                guide
                    .mismatches
                    .iter()
                    .map(|&x| guide_upper.len() - (x + 1))
                    .collect(),
            );
            sg_rna_include_idxs.push(window_from_cut(
                cut_p,
                guide.qw_size,
                ref_seq_length,
                &mut all_include_idxs,
            ));
            sg_rna_names.push(resolve_guide_name(
                &guide.name,
                &guide_upper,
                match_count,
                start,
                &mut seen_names,
            ));
            sg_rna_sequences.push(guide_upper.clone());
        }
    }

    let mut include_idxs = if let Some(coords) = quant_window_coordinates {
        if coords != "0" && !coords.is_empty() {
            parse_quant_window_coordinates(coords, ref_seq_length)
        } else if !sg_rna_cut_points.is_empty() && !all_include_idxs.is_empty() {
            dedup_sorted(all_include_idxs)
        } else {
            (0..ref_seq_length).collect()
        }
    } else if !sg_rna_cut_points.is_empty() && !all_include_idxs.is_empty() {
        dedup_sorted(all_include_idxs)
    } else {
        (0..ref_seq_length).collect()
    };

    let exclude_set: HashSet<usize> = exclude_idxs.iter().copied().collect();
    include_idxs.retain(|idx| !exclude_set.contains(idx));

    GuideMatchResult {
        sg_rna_sequences,
        sg_rna_intervals,
        sg_rna_cut_points,
        sg_rna_mismatches,
        sg_rna_names,
        sg_rna_include_idxs,
        include_idxs,
        exclude_idxs: dedup_sorted(exclude_idxs),
    }
}

pub fn build_gap_incentive(seq_len: usize, cut_points: &[usize], incentive_value: i64) -> Vec<i64> {
    let mut gap_incentive = vec![0_i64; seq_len + 1];
    for &cut_point in cut_points {
        if cut_point + 1 < gap_incentive.len() {
            gap_incentive[cut_point + 1] = incentive_value;
        }
    }
    gap_incentive
}

pub fn generate_seeds(
    seq: &str,
    exclude_left: usize,
    exclude_right: usize,
    seed_len: usize,
    seed_step: usize,
) -> (Vec<String>, Vec<String>) {
    let seq_upper = seq.to_uppercase();
    let seq_rc = reverse_complement(&seq_upper);
    let seq_len = seq_upper.len();
    let mut seeds = Vec::new();
    let mut rc_seeds = Vec::new();

    if seed_len == 0 || seed_len > seq_len || seed_step == 0 {
        return (seeds, rc_seeds);
    }

    let stop = seq_len.saturating_sub(exclude_right + seed_len);
    let mut seed_start = exclude_left;
    while seed_start < stop {
        let mut attempts_to_find_seed = 0usize;
        let mut this_seed_start = seed_start;
        let mut potential_seed = &seq_upper[this_seed_start..this_seed_start + seed_len];

        while seq_rc.contains(potential_seed) || seeds.iter().any(|s| s == potential_seed) {
            attempts_to_find_seed += 1;
            if attempts_to_find_seed > 100 {
                break;
            }
            if this_seed_start > seq_len.saturating_sub(seed_len) {
                this_seed_start = 0;
            }
            this_seed_start += 1;
            if this_seed_start + seed_len > seq_len {
                break;
            }
            potential_seed = &seq_upper[this_seed_start..this_seed_start + seed_len];
        }

        if potential_seed.len() == seed_len {
            let seed_rc = reverse_complement(potential_seed);
            if !seq_upper.contains(&seed_rc) && !seq_rc.contains(potential_seed) {
                seeds.push(potential_seed.to_string());
                rc_seeds.push(seed_rc);
            }
        }

        seed_start += seed_step;
    }

    (seeds, rc_seeds)
}

pub fn get_quant_window_ranges_from_include_idxs(include_idxs: &[usize]) -> Vec<(usize, usize)> {
    if include_idxs.is_empty() {
        return Vec::new();
    }
    let mut ranges = Vec::new();
    let mut start = include_idxs[0];
    let mut last = include_idxs[0];
    for &idx in &include_idxs[1..] {
        if idx == last + 1 {
            last = idx;
        } else {
            ranges.push((start, last));
            start = idx;
            last = idx;
        }
    }
    ranges.push((start, last));
    ranges
}

fn find_all_matches(haystack: &str, needle: &str) -> Vec<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return Vec::new();
    }
    let mut starts = Vec::new();
    let last_start = haystack.len() - needle.len();
    for start in 0..=last_start {
        if &haystack[start..start + needle.len()] == needle {
            starts.push(start);
        }
    }
    starts
}

fn parse_quant_window_coordinates(coords: &str, ref_len: usize) -> Vec<usize> {
    let mut result = Vec::new();
    for part in coords.split('_') {
        let mut pieces = part.split('-');
        let start = pieces.next().and_then(|v| v.parse::<usize>().ok());
        let end = pieces.next().and_then(|v| v.parse::<usize>().ok());
        if let (Some(start), Some(end)) = (start, end) {
            for idx in start..=end {
                if idx < ref_len {
                    result.push(idx);
                }
            }
        }
    }
    dedup_sorted(result)
}

fn window_from_cut(
    cut_p: usize,
    qw_size: i32,
    ref_len: usize,
    all_include: &mut Vec<usize>,
) -> Vec<usize> {
    if qw_size <= 0 || ref_len == 0 {
        return Vec::new();
    }
    let st = (cut_p as i32 - qw_size + 1).max(0) as usize;
    let en_exclusive = (cut_p + qw_size as usize + 1).min(ref_len.saturating_sub(1));
    if st >= en_exclusive {
        return Vec::new();
    }
    let idxs: Vec<usize> = (st..en_exclusive).collect();
    all_include.extend(idxs.iter().copied());
    idxs
}

fn resolve_guide_name(
    base_name: &str,
    guide_seq: &str,
    match_count: usize,
    position: usize,
    seen_names: &mut HashSet<String>,
) -> String {
    if match_count == 1 {
        return base_name.to_string();
    }

    let base = if base_name.is_empty() {
        guide_seq.to_string()
    } else {
        base_name.to_string()
    };

    let mut candidate = format!("{}_{}", base, position);
    let mut idx = 1usize;
    while seen_names.contains(&candidate) {
        candidate = format!("{}_{}_{}", base, position, idx);
        idx += 1;
    }
    seen_names.insert(candidate.clone());
    candidate
}

fn dedup_sorted(mut values: Vec<usize>) -> Vec<usize> {
    values.sort_unstable();
    values.dedup();
    values
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_quant_window_coordinates() {
        let result = parse_quant_window_coordinates("2-5_8-10", 20);
        assert_eq!(result, vec![2, 3, 4, 5, 8, 9, 10]);
    }

    #[test]
    fn test_build_gap_incentive() {
        let result = build_gap_incentive(10, &[3, 7], 1);
        assert_eq!(result.len(), 11);
        assert_eq!(result[4], 1);
        assert_eq!(result[8], 1);
        assert_eq!(result[0], 0);
    }

    #[test]
    fn test_generate_seeds_basic() {
        let seq = "ATCGATCGATCGATCG";
        let (seeds, rc_seeds) = generate_seeds(seq, 0, 0, 5, 3);
        assert_eq!(seeds.len(), rc_seeds.len());
    }

    #[test]
    fn test_get_quant_window_ranges_from_include_idxs() {
        assert_eq!(
            get_quant_window_ranges_from_include_idxs(&[0, 1, 2, 10, 11, 12]),
            vec![(0, 2), (10, 12)]
        );
    }
}

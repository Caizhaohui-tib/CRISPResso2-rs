use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EditPayload {
    pub all_insertion_positions: Vec<usize>,
    pub all_insertion_left_positions: Vec<usize>,
    pub insertion_positions: Vec<usize>,
    pub insertion_coordinates: Vec<(usize, usize)>,
    pub insertion_sizes: Vec<usize>,
    pub insertion_n: usize,
    pub all_deletion_positions: Vec<usize>,
    pub all_deletion_coordinates: Vec<(usize, usize)>,
    pub deletion_positions: Vec<usize>,
    pub deletion_coordinates: Vec<(usize, usize)>,
    pub deletion_sizes: Vec<usize>,
    pub deletion_n: usize,
    pub all_substitution_positions: Vec<usize>,
    pub substitution_positions: Vec<usize>,
    pub all_substitution_values: Vec<char>,
    pub substitution_values: Vec<char>,
    pub substitution_n: usize,
    pub ref_positions: Vec<isize>,
}

pub fn find_indels_substitutions(
    read_seq_al: &str,
    ref_seq_al: &str,
    include_idx: &[usize],
) -> EditPayload {
    let include_idx_set: HashSet<usize> = include_idx.iter().copied().collect();
    let read_chars: Vec<char> = read_seq_al.chars().collect();
    let ref_chars: Vec<char> = ref_seq_al.chars().collect();

    let mut ref_positions = Vec::with_capacity(ref_chars.len());
    let mut all_substitution_positions = Vec::new();
    let mut substitution_positions = Vec::new();
    let mut all_substitution_values = Vec::new();
    let mut substitution_values = Vec::new();

    let mut all_deletion_positions = Vec::new();
    let mut all_deletion_coordinates = Vec::new();
    let mut deletion_positions = Vec::new();
    let mut deletion_coordinates = Vec::new();
    let mut deletion_sizes = Vec::new();
    let mut start_deletion: Option<usize> = None;

    let mut all_insertion_positions = Vec::new();
    let mut all_insertion_left_positions = Vec::new();
    let mut insertion_positions = Vec::new();
    let mut insertion_coordinates = Vec::new();
    let mut insertion_sizes = Vec::new();
    let mut start_insertion: Option<usize> = None;
    let mut current_insertion_size = 0_usize;

    let mut idx = 0_usize;
    for (idx_c, &ref_char) in ref_chars.iter().enumerate() {
        if ref_char != '-' {
            ref_positions.push(idx as isize);
            let read_char = read_chars[idx_c];
            if ref_char != read_char && read_char != '-' && read_char != 'N' {
                all_substitution_positions.push(idx);
                all_substitution_values.push(read_char);
                if include_idx_set.contains(&idx) {
                    substitution_positions.push(idx);
                    substitution_values.push(read_char);
                }
            }
            if let Some(start) = start_insertion {
                all_insertion_left_positions.push(start);
                all_insertion_positions.push(start);
                all_insertion_positions.push(idx);
                if include_idx_set.contains(&start) && include_idx_set.contains(&idx) {
                    insertion_coordinates.push((start, idx));
                    insertion_positions.push(start);
                    insertion_positions.push(idx);
                    insertion_sizes.push(current_insertion_size);
                }
                start_insertion = None;
            }
            current_insertion_size = 0;
            idx += 1;
        } else {
            if idx == 0 {
                ref_positions.push(-1);
            } else {
                ref_positions.push(-(idx as isize));
            }
            if idx > 0 && start_insertion.is_none() {
                start_insertion = Some(idx - 1);
            }
            current_insertion_size += 1;
        }

        if read_chars[idx_c] == '-' && start_deletion.is_none() {
            start_deletion = Some(if idx_c >= 1 {
                ref_positions[idx_c] as usize
            } else {
                0
            });
        } else if read_chars[idx_c] != '-' {
            if let Some(start) = start_deletion {
                let end = ref_positions[idx_c] as usize;
                all_deletion_positions.extend(start..end);
                all_deletion_coordinates.push((start, end));
                if (start..end).any(|pos| include_idx_set.contains(&pos)) {
                    deletion_positions.extend(start..end);
                    deletion_coordinates.push((start, end));
                    deletion_sizes.push(end - start);
                }
                start_deletion = None;
            }
        }
    }

    if let Some(start) = start_deletion {
        let end = ref_positions[ref_positions.len() - 1] as usize + 1;
        all_deletion_positions.extend(start..end);
        all_deletion_coordinates.push((start, end));
        if (start..end).any(|pos| include_idx_set.contains(&pos)) {
            deletion_positions.extend(start..end);
            deletion_coordinates.push((start, end));
            deletion_sizes.push(end - start);
        }
    }

    let _substitution_n = substitution_positions.len();
    let deletion_n = deletion_sizes.iter().sum();
    let insertion_n = insertion_sizes.iter().sum();

    let substitution_n = substitution_values.len();

    EditPayload {
        all_insertion_positions,
        all_insertion_left_positions,
        insertion_positions,
        insertion_coordinates,
        insertion_sizes,
        insertion_n,
        all_deletion_positions,
        all_deletion_coordinates,
        deletion_positions,
        deletion_coordinates,
        deletion_sizes,
        deletion_n,
        all_substitution_positions,
        substitution_positions,
        all_substitution_values,
        substitution_values,
        substitution_n,
        ref_positions,
    }
}

pub fn find_indels_substitutions_legacy(
    read_seq_al: &str,
    ref_seq_al: &str,
    include_idx: &[usize],
) -> EditPayload {
    let include_idx_set: HashSet<usize> = include_idx.iter().copied().collect();
    let read_chars: Vec<char> = read_seq_al.chars().collect();
    let ref_chars: Vec<char> = ref_seq_al.chars().collect();

    let mut ref_positions = Vec::with_capacity(ref_chars.len());
    let mut all_substitution_positions = Vec::new();
    let mut substitution_positions = Vec::new();
    let mut all_substitution_values = Vec::new();
    let mut substitution_values = Vec::new();

    let mut idx = 0usize;
    for (idx_c, &ref_char) in ref_chars.iter().enumerate() {
        if ref_char != '-' {
            ref_positions.push(idx as isize);
            let read_char = read_chars[idx_c];
            if ref_char != read_char && read_char != '-' && read_char != 'N' {
                all_substitution_positions.push(idx);
                all_substitution_values.push(read_char);
                if include_idx_set.contains(&idx) {
                    substitution_positions.push(idx);
                    substitution_values.push(read_char);
                }
            }
            idx += 1;
        } else if idx == 0 {
            ref_positions.push(-1);
        } else {
            ref_positions.push(-(idx as isize));
        }
    }

    let mut all_deletion_positions = Vec::new();
    let mut all_deletion_coordinates = Vec::new();
    let mut deletion_positions = Vec::new();
    let mut deletion_coordinates = Vec::new();
    let mut deletion_sizes = Vec::new();
    let mut idx_c = 0usize;
    while idx_c < read_chars.len() {
        if read_chars[idx_c] != '-' {
            idx_c += 1;
            continue;
        }
        let start = idx_c;
        while idx_c < read_chars.len() && read_chars[idx_c] == '-' {
            idx_c += 1;
        }
        let end = idx_c;
        let ref_st = if start > 1 {
            ref_positions[start] as usize
        } else {
            0
        };
        let ref_en = if end < ref_positions.len() {
            ref_positions[end] as usize
        } else {
            idx.saturating_sub(1)
        };
        all_deletion_positions.extend(ref_st..ref_en);
        all_deletion_coordinates.push((ref_st, ref_en));
        if (ref_st..ref_en).any(|pos| include_idx_set.contains(&pos)) {
            deletion_positions.extend(ref_st..ref_en);
            deletion_coordinates.push((ref_st, ref_en));
            deletion_sizes.push(end - start);
        }
    }

    let mut all_insertion_positions = Vec::new();
    let mut all_insertion_left_positions = Vec::new();
    let mut insertion_positions = Vec::new();
    let mut insertion_coordinates = Vec::new();
    let mut insertion_sizes = Vec::new();
    idx_c = 0;
    while idx_c < ref_chars.len() {
        if ref_chars[idx_c] != '-' {
            idx_c += 1;
            continue;
        }
        let start = idx_c;
        while idx_c < ref_chars.len() && ref_chars[idx_c] == '-' {
            idx_c += 1;
        }
        let end = idx_c;
        if start == 0 || end == ref_chars.len() {
            continue;
        }
        let ref_st = ref_positions[start - 1] as usize;
        let ref_en = ref_positions[end] as usize;
        all_insertion_left_positions.push(ref_st);
        all_insertion_positions.push(ref_st);
        all_insertion_positions.push(ref_en);
        if include_idx_set.contains(&ref_st) || include_idx_set.contains(&ref_en) {
            insertion_coordinates.push((ref_st, ref_en));
            insertion_positions.push(ref_st);
            insertion_positions.push(ref_en);
            insertion_sizes.push(end - start);
        }
    }

    let substitution_n = substitution_values.len();

    EditPayload {
        all_insertion_positions,
        all_insertion_left_positions,
        insertion_positions,
        insertion_coordinates,
        insertion_sizes: insertion_sizes.clone(),
        insertion_n: insertion_sizes.iter().sum(),
        all_deletion_positions,
        all_deletion_coordinates,
        deletion_positions,
        deletion_coordinates,
        deletion_sizes: deletion_sizes.clone(),
        deletion_n: deletion_sizes.iter().sum(),
        all_substitution_positions,
        substitution_positions,
        all_substitution_values,
        substitution_values,
        substitution_n,
        ref_positions,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_single_deletion() {
        let payload = find_indels_substitutions("A-T", "ACT", &[0, 1, 2]);
        assert_eq!(payload.ref_positions, vec![0, 1, 2]);
        assert_eq!(payload.all_deletion_positions, vec![1]);
        assert_eq!(payload.all_deletion_coordinates, vec![(1, 2)]);
        assert_eq!(payload.deletion_positions, vec![1]);
        assert_eq!(payload.deletion_sizes, vec![1]);
        assert_eq!(payload.deletion_n, 1);
    }

    #[test]
    fn detects_single_insertion_between_included_anchors() {
        let payload = find_indels_substitutions("ATGT", "AT-T", &[0, 1, 2]);
        assert_eq!(payload.ref_positions, vec![0, 1, -2, 2]);
        assert_eq!(payload.all_insertion_left_positions, vec![1]);
        assert_eq!(payload.all_insertion_positions, vec![1, 2]);
        assert_eq!(payload.insertion_coordinates, vec![(1, 2)]);
        assert_eq!(payload.insertion_sizes, vec![1]);
        assert_eq!(payload.insertion_n, 1);
    }

    #[test]
    fn insertion_requires_both_anchors_in_quantification_window() {
        let payload = find_indels_substitutions("ATGT", "AT-T", &[1]);
        assert_eq!(payload.all_insertion_positions, vec![1, 2]);
        assert!(payload.insertion_positions.is_empty());
        assert_eq!(payload.insertion_n, 0);
    }

    #[test]
    fn detects_substitutions_and_ignores_read_n() {
        let payload = find_indels_substitutions("AGNA", "ATTA", &[0, 1, 2, 3]);
        assert_eq!(payload.all_substitution_positions, vec![1]);
        assert_eq!(payload.substitution_positions, vec![1]);
        assert_eq!(payload.all_substitution_values, vec!['G']);
        assert_eq!(payload.substitution_n, 1);
    }

    #[test]
    fn deletion_outside_window_is_counted_separately() {
        let payload = find_indels_substitutions("A-T", "ACT", &[0]);
        assert_eq!(payload.all_deletion_coordinates, vec![(1, 2)]);
        assert!(payload.deletion_coordinates.is_empty());
        assert_eq!(payload.deletion_n, 0);
    }

    #[test]
    fn trailing_deletion_matches_python_end_coordinate() {
        let payload = find_indels_substitutions("AT--", "ATCG", &[0, 1, 2, 3]);
        assert_eq!(payload.all_deletion_positions, vec![2, 3]);
        assert_eq!(payload.all_deletion_coordinates, vec![(2, 4)]);
        assert_eq!(payload.deletion_n, 2);
    }

    #[test]
    fn first_position_deletion_matches_python_unit_test() {
        let include: Vec<usize> = (0..8).collect();
        let payload = find_indels_substitutions("-TGCGTAC", "ATGCGTAC", &include);
        assert_eq!(payload.deletion_coordinates, vec![(0, 1)]);
    }

    #[test]
    fn last_position_deletion_matches_python_unit_test() {
        let include: Vec<usize> = (0..8).collect();
        let payload = find_indels_substitutions("ATGCGTA-", "ATGCGTAC", &include);
        assert_eq!(payload.deletion_n, 1);
        assert_eq!(payload.deletion_positions, vec![7]);
        assert_eq!(payload.all_deletion_positions, vec![7]);
        assert_eq!(payload.deletion_coordinates, vec![(7, 8)]);
        assert_eq!(payload.all_deletion_coordinates, vec![(7, 8)]);
    }

    #[test]
    fn second_and_third_position_deletions_match_python_unit_tests() {
        let include: Vec<usize> = (0..8).collect();
        let second = find_indels_substitutions("A-GCGTAC", "ATGCGTAC", &include);
        assert_eq!(second.deletion_coordinates, vec![(1, 2)]);
        let third = find_indels_substitutions("AT-CGTAC", "ATGCGTAC", &include);
        assert_eq!(third.deletion_coordinates, vec![(2, 3)]);
    }

    #[test]
    fn end_deletion_spans_match_python_unit_tests() {
        let include: Vec<usize> = (0..8).collect();
        let two = find_indels_substitutions("ATGCGT--", "ATGCGTAC", &include);
        assert_eq!(two.deletion_n, 2);
        assert_eq!(two.deletion_coordinates, vec![(6, 8)]);
        assert_eq!(two.deletion_positions, vec![6, 7]);

        let three = find_indels_substitutions("ATGCG---", "ATGCGTAC", &include);
        assert_eq!(three.deletion_n, 3);
        assert_eq!(three.deletion_coordinates, vec![(5, 8)]);
        assert_eq!(three.deletion_positions, vec![5, 6, 7]);
    }

    #[test]
    fn complex_payload_matches_python_nonlegacy_case() {
        let aln_seq = "TAATCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAAGAGGGCGGCTTTGGGCGGGGTC-CAGTTCCGGGATTA--GCGAACTTAGAGCAC-----ACGTCTGAACTCCAGTCACCGATGTATATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAACTTACTCTCACTTAACTCTTGCTTCCCTCCTGACGCCGATG";
        let ref_seq = "----CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGC------------ACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCA-AGCACTACCTACGTCAGCACCTGGGACCCCGCCAC------CGTGCGCCGGGC----CTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGG--CGCTTTGGTCGG-----";
        let payload = find_indels_substitutions(aln_seq, ref_seq, &[91, 92]);
        assert_eq!(
            payload.all_insertion_positions,
            vec![91, 92, 126, 127, 161, 162, 173, 174, 210, 211]
        );
        assert_eq!(
            payload.all_insertion_left_positions,
            vec![91, 126, 161, 173, 210]
        );
        assert_eq!(payload.insertion_positions, vec![91, 92]);
        assert_eq!(payload.insertion_coordinates, vec![(91, 92)]);
        assert_eq!(payload.insertion_sizes, vec![12]);
        assert_eq!(payload.insertion_n, 12);
        assert_eq!(
            payload.all_deletion_positions,
            vec![101, 116, 117, 132, 133, 134, 135, 136]
        );
        assert_eq!(
            payload.all_deletion_coordinates,
            vec![(101, 102), (116, 118), (132, 137)]
        );
        assert!(payload.deletion_positions.is_empty());
        assert_eq!(payload.deletion_n, 0);
        assert_eq!(payload.substitution_positions, vec![91, 92]);
        assert_eq!(payload.substitution_n, 2);
        assert_eq!(payload.ref_positions.len(), ref_seq.len());
    }

    #[test]
    fn complex_payload_matches_python_legacy_case() {
        let aln_seq = "TAATCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAAGAGGGCGGCTTTGGGCGGGGTC-CAGTTCCGGGATTA--GCGAACTTAGAGCAC-----ACGTCTGAACTCCAGTCACCGATGTATATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAACTTACTCTCACTTAACTCTTGCTTCCCTCCTGACGCCGATG";
        let ref_seq = "----CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGC------------ACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCA-AGCACTACCTACGTCAGCACCTGGGACCCCGCCAC------CGTGCGCCGGGC----CTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGG--CGCTTTGGTCGG-----";
        let payload = find_indels_substitutions_legacy(aln_seq, ref_seq, &[91, 92]);
        assert_eq!(
            payload.all_insertion_positions,
            vec![91, 92, 126, 127, 161, 162, 173, 174, 210, 211]
        );
        assert_eq!(payload.insertion_positions, vec![91, 92]);
        assert_eq!(payload.insertion_sizes, vec![12]);
        assert_eq!(payload.insertion_n, 12);
        assert!(payload.deletion_positions.is_empty());
        assert_eq!(payload.substitution_positions, vec![91, 92]);
        assert_eq!(payload.substitution_n, 2);
        assert_eq!(payload.ref_positions.len(), ref_seq.len());
    }

    #[test]
    fn legacy_insertion_quantification_counts_single_anchor_overlap() {
        let nonlegacy = find_indels_substitutions("ATGT", "AT-T", &[1]);
        assert_eq!(nonlegacy.insertion_n, 0);

        let legacy = find_indels_substitutions_legacy("ATGT", "AT-T", &[1]);
        assert_eq!(legacy.insertion_coordinates, vec![(1, 2)]);
        assert_eq!(legacy.insertion_positions, vec![1, 2]);
        assert_eq!(legacy.insertion_sizes, vec![1]);
        assert_eq!(legacy.insertion_n, 1);
    }
}

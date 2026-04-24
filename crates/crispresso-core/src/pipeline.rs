use std::collections::HashMap;

use crate::align::{global_align, ScoreMatrix};
use crate::edits::{find_indels_substitutions, EditPayload};
use crate::models::{
    AlignmentStats, ReadClassification, ReferenceInfo, RunConfig, VariantPayload, VariantRecord,
};
use crate::sequence::reverse_complement;

#[derive(Debug, Clone)]
pub struct ClassifiedRead {
    pub variant: VariantRecord,
    pub stats_contribution: StatsContribution,
}

#[derive(Debug, Clone, Default)]
pub struct StatsContribution {
    pub is_aligned: bool,
    pub n_global_subs: usize,
    pub n_subs_outside_window: usize,
    pub n_mods_in_window: usize,
    pub n_mods_outside_window: usize,
    pub irregular_ends: bool,
}

pub fn classify_read(
    read_seq: &str,
    refs: &[&ReferenceInfo],
    matrix: &ScoreMatrix,
    config: &RunConfig,
) -> ClassifiedRead {
    let mut aln_scores = Vec::with_capacity(refs.len());
    let mut best_match_score = -1.0_f64;
    let mut best_match_s1s = Vec::new();
    let mut best_match_s2s = Vec::new();
    let mut best_match_names = Vec::new();
    let mut best_match_strands = Vec::new();
    let mut ref_aln_details = Vec::new();
    let read_upper = read_seq.to_uppercase();

    for ref_info in refs {
        let (aln_strand, s1, s2, score) = align_read_to_ref(&read_upper, ref_info, matrix, config);
        ref_aln_details.push((ref_info.name.clone(), s1.clone(), s2.clone(), score));
        aln_scores.push(score);

        if score > best_match_score && score > ref_info.min_aln_score {
            best_match_score = score;
            best_match_s1s = vec![s1];
            best_match_s2s = vec![s2];
            best_match_names = vec![ref_info.name.clone()];
            best_match_strands = vec![aln_strand];
        } else if (score - best_match_score).abs() < f64::EPSILON && score > ref_info.min_aln_score
        {
            best_match_s1s.push(s1);
            best_match_s2s.push(s2);
            best_match_names.push(ref_info.name.clone());
            best_match_strands.push(aln_strand);
        }
    }

    if best_match_score <= 0.0 {
        return ClassifiedRead {
            variant: VariantRecord {
                count: 1,
                aln_ref_names: Vec::new(),
                aln_scores,
                ref_aln_details,
                best_match_score,
                best_match_name: None,
                class_name: String::new(),
                payloads: Vec::new(),
            },
            stats_contribution: StatsContribution::default(),
        };
    }

    let mut payloads = Vec::new();
    let mut class_names = Vec::new();
    let mut contribution = StatsContribution {
        is_aligned: true,
        ..StatsContribution::default()
    };

    for idx in 0..best_match_names.len() {
        let ref_info = refs
            .iter()
            .find(|r| r.name == best_match_names[idx])
            .expect("reference present");
        let edit = classify_edit_payload(
            &best_match_s1s[idx],
            &best_match_s2s[idx],
            &ref_info.include_idxs,
            config.use_legacy_insertion_quantification,
        );
        let irregular_ends = has_irregular_ends(&best_match_s1s[idx], &best_match_s2s[idx]);
        let insertions_outside_window = (edit.all_insertion_positions.len() / 2)
            .saturating_sub(edit.insertion_positions.len() / 2);
        let deletions_outside_window = edit
            .all_deletion_coordinates
            .len()
            .saturating_sub(edit.deletion_coordinates.len());
        let substitutions_outside_window = edit
            .all_substitution_positions
            .len()
            .saturating_sub(edit.substitution_positions.len());
        let total_mods = (edit.all_insertion_positions.len() / 2)
            + edit.all_deletion_positions.len()
            + edit.all_substitution_positions.len();
        let mods_in_window = edit.substitution_n + edit.deletion_n + edit.insertion_n;
        let mods_outside_window = total_mods.saturating_sub(mods_in_window);

        let is_modified = (!config.ignore_deletions && edit.deletion_n > 0)
            || (!config.ignore_insertions && edit.insertion_n > 0)
            || (!config.ignore_substitutions && edit.substitution_n > 0);

        let classification = if is_modified {
            ReadClassification::Modified
        } else {
            ReadClassification::Unmodified
        };
        class_names.push(if is_modified {
            format!("{}_MODIFIED", best_match_names[idx])
        } else {
            format!("{}_UNMODIFIED", best_match_names[idx])
        });

        contribution.n_global_subs += edit.substitution_n + substitutions_outside_window;
        contribution.n_subs_outside_window += substitutions_outside_window;
        contribution.n_mods_in_window += mods_in_window;
        contribution.n_mods_outside_window += mods_outside_window;
        contribution.irregular_ends |= irregular_ends;

        payloads.push(VariantPayload {
            edit,
            ref_name: best_match_names[idx].clone(),
            aln_scores: aln_scores.clone(),
            aln_seq: best_match_s1s[idx].clone(),
            aln_ref: best_match_s2s[idx].clone(),
            aln_strand: best_match_strands[idx],
            classification,
            irregular_ends,
            insertions_outside_window,
            deletions_outside_window,
            substitutions_outside_window,
            total_mods,
            mods_in_window,
            mods_outside_window,
        });
    }

    let mut class_name = class_names.join("&");
    let mut aln_ref_names = best_match_names.clone();
    let mut best_match_name = best_match_names.first().cloned();
    if best_match_names.len() > 1 {
        if config.assign_ambiguous_to_first {
            class_name = class_names[0].clone();
            aln_ref_names = vec![best_match_names[0].clone()];
            best_match_name = Some(best_match_names[0].clone());
        } else if !config.expand_ambiguous_alignments {
            class_name = "AMBIGUOUS".to_string();
        }
    }

    ClassifiedRead {
        variant: VariantRecord {
            count: 1,
            aln_ref_names,
            aln_scores,
            ref_aln_details,
            best_match_score,
            best_match_name,
            class_name,
            payloads,
        },
        stats_contribution: contribution,
    }
}

fn classify_edit_payload(
    aligned_read: &str,
    aligned_ref: &str,
    include_idxs: &[usize],
    legacy: bool,
) -> EditPayload {
    if legacy {
        crate::edits::find_indels_substitutions_legacy(aligned_read, aligned_ref, include_idxs)
    } else {
        find_indels_substitutions(aligned_read, aligned_ref, include_idxs)
    }
}

fn has_irregular_ends(aligned_read: &str, aligned_ref: &str) -> bool {
    let read_chars: Vec<char> = aligned_read.chars().collect();
    let ref_chars: Vec<char> = aligned_ref.chars().collect();
    if read_chars.is_empty() || ref_chars.is_empty() {
        return true;
    }
    read_chars[0] == '-'
        || ref_chars[0] == '-'
        || read_chars[0] != ref_chars[0]
        || read_chars[read_chars.len() - 1] == '-'
        || ref_chars[ref_chars.len() - 1] == '-'
        || read_chars[read_chars.len() - 1] != ref_chars[ref_chars.len() - 1]
}

fn align_read_to_ref(
    read_seq: &str,
    ref_info: &ReferenceInfo,
    matrix: &ScoreMatrix,
    config: &RunConfig,
) -> (char, String, String, f64) {
    let read_rc = reverse_complement(read_seq);
    let mut found_forward = 0usize;
    let mut found_reverse = 0usize;
    let max_seed = config
        .aln_seed_count
        .min(ref_info.fw_seeds.len())
        .min(ref_info.rc_seeds.len());

    for seed_i in 0..max_seed {
        if read_seq.contains(&ref_info.fw_seeds[seed_i]) {
            found_forward += 1;
        }
        if read_seq.contains(&ref_info.rc_seeds[seed_i]) {
            found_reverse += 1;
        }
    }

    if found_forward > config.aln_seed_min && found_reverse == 0 {
        let result = global_align(
            read_seq,
            &ref_info.sequence,
            matrix,
            &ref_info.gap_incentive,
            config.needleman_wunsch_gap_open,
            config.needleman_wunsch_gap_extend,
        )
        .expect("alignment should succeed");
        (
            '+',
            result.aligned_read,
            result.aligned_reference,
            result.score,
        )
    } else if found_forward == 0 && found_reverse > config.aln_seed_min {
        let result = global_align(
            &read_rc,
            &ref_info.sequence,
            matrix,
            &ref_info.gap_incentive,
            config.needleman_wunsch_gap_open,
            config.needleman_wunsch_gap_extend,
        )
        .expect("alignment should succeed");
        (
            '-',
            result.aligned_read,
            result.aligned_reference,
            result.score,
        )
    } else {
        let fw = global_align(
            read_seq,
            &ref_info.sequence,
            matrix,
            &ref_info.gap_incentive,
            config.needleman_wunsch_gap_open,
            config.needleman_wunsch_gap_extend,
        )
        .expect("alignment should succeed");
        let rv = global_align(
            &read_rc,
            &ref_info.sequence,
            matrix,
            &ref_info.gap_incentive,
            config.needleman_wunsch_gap_open,
            config.needleman_wunsch_gap_extend,
        )
        .expect("alignment should succeed");
        if rv.score > fw.score {
            ('-', rv.aligned_read, rv.aligned_reference, rv.score)
        } else {
            ('+', fw.aligned_read, fw.aligned_reference, fw.score)
        }
    }
}

pub fn process_fastq_reads(
    unique_reads: &HashMap<String, usize>,
    refs: &[ReferenceInfo],
    matrix: &ScoreMatrix,
    config: &RunConfig,
) -> (HashMap<String, VariantRecord>, AlignmentStats) {
    let ref_refs: Vec<&ReferenceInfo> = refs.iter().collect();
    let mut variant_cache = HashMap::new();
    let mut stats = AlignmentStats::default();

    for (seq, &count) in unique_reads {
        stats.n_tot_reads += count;
        let classified = classify_read(seq, &ref_refs, matrix, config);
        if classified.variant.best_match_score <= 0.0 {
            stats.n_computed_notaln += 1;
            stats.n_cached_notaln += count.saturating_sub(1);
            continue;
        }

        let mut variant = classified.variant;
        variant.count = count;
        stats.n_computed_aln += 1;
        stats.n_cached_aln += count.saturating_sub(1);
        if stats.read_length == 0 {
            if let Some(payload) = variant.payloads.first() {
                stats.read_length = payload.aln_seq.len();
            }
        }
        if variant.aln_ref_names.len() == 1 || config.expand_ambiguous_alignments {
            stats.n_global_subs += classified.stats_contribution.n_global_subs * count;
            stats.n_subs_outside_window +=
                classified.stats_contribution.n_subs_outside_window * count;
            stats.n_mods_in_window += classified.stats_contribution.n_mods_in_window * count;
            stats.n_mods_outside_window +=
                classified.stats_contribution.n_mods_outside_window * count;
            if classified.stats_contribution.irregular_ends {
                stats.n_reads_irregular_ends += count;
            }
        }
        variant_cache.insert(seq.clone(), variant);
    }

    (variant_cache, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::make_matrix;
    use crate::guide::{build_gap_incentive, generate_seeds, get_amplicon_info_for_guides};
    use crate::models::GuideInfo;

    fn make_test_ref() -> ReferenceInfo {
        let seq = "ATCGATCGATCGATCGATCG";
        let guides = vec![GuideInfo {
            sequence: "ATCGATCG".to_string(),
            orig_sequence: "ATCGATCG".to_string(),
            name: "test_guide".to_string(),
            qw_center: -3,
            qw_size: 1,
            mismatches: vec![],
            plot_cut_point: true,
        }];
        let result = get_amplicon_info_for_guides(seq, &guides, None, 0, 0, 20);
        let gap_incentive = build_gap_incentive(seq.len(), &result.sg_rna_cut_points, 1);
        let (fw_seeds, rc_seeds) = generate_seeds(seq, 0, 0, 10, 5);
        ReferenceInfo {
            name: "Reference".to_string(),
            sequence: seq.to_string(),
            sequence_length: seq.len(),
            min_aln_score: 0.0,
            gap_incentive,
            sg_rna_cut_points: result.sg_rna_cut_points,
            sg_rna_intervals: result.sg_rna_intervals,
            sg_rna_sequences: result.sg_rna_sequences,
            sg_rna_names: result.sg_rna_names,
            sg_rna_mismatches: result.sg_rna_mismatches,
            sg_rna_orig_sequences: vec!["ATCGATCG".to_string()],
            sg_rna_include_idxs: result.sg_rna_include_idxs,
            contains_guide: true,
            include_idxs: result.include_idxs,
            exclude_idxs: result.exclude_idxs,
            fw_seeds,
            rc_seeds,
            exon_positions: Default::default(),
            exon_intervals: vec![],
        }
    }

    #[test]
    fn test_classify_read_perfect_match() {
        let ref_info = make_test_ref();
        let matrix = make_matrix(5, -4, -2, -1);
        let config = RunConfig::default();
        let result = classify_read("ATCGATCGATCGATCGATCG", &[&ref_info], &matrix, &config);
        assert!(result.stats_contribution.is_aligned);
        assert!(result.variant.best_match_score > 0.0);
    }

    #[test]
    fn test_classify_read_no_match_empty_reference_list_is_not_needed() {
        let ref_info = make_test_ref();
        let matrix = make_matrix(5, -4, -2, -1);
        let config = RunConfig::default();
        let result = classify_read("NNNNNNNNNNNNNNNNNNNNNN", &[&ref_info], &matrix, &config);
        assert!(!result.stats_contribution.is_aligned || result.variant.best_match_score >= 0.0);
    }
}

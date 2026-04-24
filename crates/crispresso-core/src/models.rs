use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct ReferenceInfo {
    pub name: String,
    pub sequence: String,
    pub sequence_length: usize,
    pub min_aln_score: f64,
    pub gap_incentive: Vec<i64>,
    pub sg_rna_cut_points: Vec<usize>,
    pub sg_rna_intervals: Vec<(usize, usize)>,
    pub sg_rna_sequences: Vec<String>,
    pub sg_rna_names: Vec<String>,
    pub sg_rna_mismatches: Vec<Vec<usize>>,
    pub sg_rna_orig_sequences: Vec<String>,
    pub sg_rna_include_idxs: Vec<Vec<usize>>,
    pub contains_guide: bool,
    pub include_idxs: Vec<usize>,
    pub exclude_idxs: Vec<usize>,
    pub fw_seeds: Vec<String>,
    pub rc_seeds: Vec<String>,
    pub exon_positions: HashSet<usize>,
    pub exon_intervals: Vec<(usize, usize)>,
}

#[derive(Debug, Clone)]
pub struct GuideInfo {
    pub sequence: String,
    pub orig_sequence: String,
    pub name: String,
    pub qw_center: i32,
    pub qw_size: i32,
    pub mismatches: Vec<usize>,
    pub plot_cut_point: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadClassification {
    Modified,
    Unmodified,
    Ambiguous,
}

#[derive(Debug, Clone)]
pub struct VariantPayload {
    pub edit: crate::edits::EditPayload,
    pub ref_name: String,
    pub aln_scores: Vec<f64>,
    pub aln_seq: String,
    pub aln_ref: String,
    pub aln_strand: char,
    pub classification: ReadClassification,
    pub irregular_ends: bool,
    pub insertions_outside_window: usize,
    pub deletions_outside_window: usize,
    pub substitutions_outside_window: usize,
    pub total_mods: usize,
    pub mods_in_window: usize,
    pub mods_outside_window: usize,
}

#[derive(Debug, Clone)]
pub struct VariantRecord {
    pub count: usize,
    pub aln_ref_names: Vec<String>,
    pub aln_scores: Vec<f64>,
    pub ref_aln_details: Vec<(String, String, String, f64)>,
    pub best_match_score: f64,
    pub best_match_name: Option<String>,
    pub class_name: String,
    pub payloads: Vec<VariantPayload>,
}

#[derive(Debug, Clone, Default)]
pub struct AlignmentStats {
    pub n_tot_reads: usize,
    pub n_cached_aln: usize,
    pub n_cached_notaln: usize,
    pub n_computed_aln: usize,
    pub n_computed_notaln: usize,
    pub n_global_subs: usize,
    pub n_subs_outside_window: usize,
    pub n_mods_in_window: usize,
    pub n_mods_outside_window: usize,
    pub n_reads_irregular_ends: usize,
    pub read_length: usize,
}

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub needleman_wunsch_gap_open: i32,
    pub needleman_wunsch_gap_extend: i32,
    pub needleman_wunsch_gap_incentive: i64,
    pub aln_seed_len: usize,
    pub aln_seed_count: usize,
    pub aln_seed_min: usize,
    pub default_min_aln_score: f64,
    pub exclude_bp_from_left: usize,
    pub exclude_bp_from_right: usize,
    pub plot_window_size: usize,
    pub quantification_window_center: i32,
    pub quantification_window_size: i32,
    pub ignore_deletions: bool,
    pub ignore_insertions: bool,
    pub ignore_substitutions: bool,
    pub use_legacy_insertion_quantification: bool,
    pub expand_ambiguous_alignments: bool,
    pub assign_ambiguous_to_first: bool,
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            needleman_wunsch_gap_open: -20,
            needleman_wunsch_gap_extend: -2,
            needleman_wunsch_gap_incentive: 1,
            aln_seed_len: 10,
            aln_seed_count: 5,
            aln_seed_min: 3,
            default_min_aln_score: 0.0,
            exclude_bp_from_left: 0,
            exclude_bp_from_right: 0,
            plot_window_size: 20,
            quantification_window_center: -3,
            quantification_window_size: 1,
            ignore_deletions: false,
            ignore_insertions: false,
            ignore_substitutions: false,
            use_legacy_insertion_quantification: false,
            expand_ambiguous_alignments: false,
            assign_ambiguous_to_first: false,
        }
    }
}

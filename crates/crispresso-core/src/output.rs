use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use crate::models::{AlignmentStats, ReadClassification, ReferenceInfo, VariantRecord};

pub fn write_allele_frequency_table(
    path: &Path,
    variant_cache: &HashMap<String, VariantRecord>,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "Aligned_Sequence\tReference_Sequence\t#Reads\tReads_Percentage\tReference_Name\tClassification\tRead_Status"
    )?;

    let total_reads: usize = variant_cache
        .values()
        .filter(|variant| variant.best_match_score > 0.0)
        .map(|variant| variant.count)
        .sum();

    let mut entries: Vec<_> = variant_cache
        .values()
        .filter(|variant| variant.best_match_score > 0.0)
        .collect();
    entries.sort_by(|a, b| b.count.cmp(&a.count));

    for variant in entries {
        for payload in &variant.payloads {
            let pct = if total_reads > 0 {
                100.0 * variant.count as f64 / total_reads as f64
            } else {
                0.0
            };
            let status = match payload.classification {
                ReadClassification::Modified => "MODIFIED",
                ReadClassification::Unmodified => "UNMODIFIED",
                ReadClassification::Ambiguous => "AMBIGUOUS",
            };
            writeln!(
                file,
                "{}\t{}\t{}\t{:.6}\t{}\t{}\t{}",
                payload.aln_seq,
                payload.aln_ref,
                variant.count,
                pct,
                payload.ref_name,
                variant.class_name,
                status,
            )?;
        }
    }

    Ok(())
}

pub fn write_quantification_of_editing_frequency(
    path: &Path,
    refs: &[ReferenceInfo],
    variant_cache: &HashMap<String, VariantRecord>,
    reads_in_input: usize,
    reads_aligned_all_amplicons: usize,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "Amplicon\tUnmodified%\tModified%\tReads_in_input\tReads_aligned_all_amplicons\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions"
    )?;

    for ref_info in refs {
        let mut reads_aligned = 0usize;
        let mut unmodified = 0usize;
        let mut modified = 0usize;
        let mut insertions = 0usize;
        let mut deletions = 0usize;
        let mut substitutions = 0usize;
        let mut only_insertions = 0usize;
        let mut only_deletions = 0usize;
        let mut only_substitutions = 0usize;
        let mut ins_del = 0usize;
        let mut ins_sub = 0usize;
        let mut del_sub = 0usize;
        let mut ins_del_sub = 0usize;

        for variant in variant_cache.values() {
            for payload in &variant.payloads {
                if payload.ref_name != ref_info.name {
                    continue;
                }
                reads_aligned += variant.count;
                let has_ins = payload.edit.insertion_n > 0;
                let has_del = payload.edit.deletion_n > 0;
                let has_sub = payload.edit.substitution_n > 0;
                match payload.classification {
                    ReadClassification::Modified => modified += variant.count,
                    ReadClassification::Unmodified => unmodified += variant.count,
                    ReadClassification::Ambiguous => {}
                }
                if has_ins {
                    insertions += variant.count;
                }
                if has_del {
                    deletions += variant.count;
                }
                if has_sub {
                    substitutions += variant.count;
                }
                match (has_ins, has_del, has_sub) {
                    (true, false, false) => only_insertions += variant.count,
                    (false, true, false) => only_deletions += variant.count,
                    (false, false, true) => only_substitutions += variant.count,
                    (true, true, false) => ins_del += variant.count,
                    (true, false, true) => ins_sub += variant.count,
                    (false, true, true) => del_sub += variant.count,
                    (true, true, true) => ins_del_sub += variant.count,
                    (false, false, false) => {}
                }
            }
        }

        let round8 = |v: f64| (v * 1e8).round() / 1e8;
        let modified_pct = if reads_aligned > 0 {
            round8(100.0 * modified as f64 / reads_aligned as f64)
        } else {
            0.0
        };
        let unmodified_pct = if reads_aligned > 0 {
            round8(100.0 * unmodified as f64 / reads_aligned as f64)
        } else {
            0.0
        };

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            ref_info.name,
            unmodified_pct,
            modified_pct,
            reads_in_input,
            reads_aligned_all_amplicons,
            reads_aligned,
            unmodified,
            modified,
            insertions,
            deletions,
            substitutions,
            only_insertions,
            only_deletions,
            only_substitutions,
            ins_del,
            ins_sub,
            del_sub,
            ins_del_sub,
        )?;
    }

    Ok(())
}

pub fn write_nucleotide_frequency_table(
    path: &Path,
    variant_cache: &HashMap<String, VariantRecord>,
    ref_info: &ReferenceInfo,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    write!(file, "\t")?;
    for (idx, ch) in ref_info.sequence.chars().enumerate() {
        write!(file, "{}", ch)?;
        if idx + 1 < ref_info.sequence_length {
            write!(file, "\t")?;
        }
    }
    writeln!(file)?;

    let bases = ['A', 'C', 'G', 'T', 'N', '-'];
    let mut counts: HashMap<char, Vec<usize>> = HashMap::new();
    for base in bases {
        counts.insert(base, vec![0; ref_info.sequence_length]);
    }

    for variant in variant_cache.values() {
        for payload in &variant.payloads {
            if payload.ref_name != ref_info.name {
                continue;
            }
            let ref_positions = &payload.edit.ref_positions;
            let aln_chars: Vec<char> = payload.aln_seq.chars().collect();
            for (aln_idx, &ref_pos) in ref_positions.iter().enumerate() {
                if ref_pos < 0 {
                    continue;
                }
                let pos = ref_pos as usize;
                if pos >= ref_info.sequence_length || aln_idx >= aln_chars.len() {
                    continue;
                }
                let base = aln_chars[aln_idx].to_ascii_uppercase();
                let key = if counts.contains_key(&base) {
                    base
                } else {
                    'N'
                };
                if let Some(vector) = counts.get_mut(&key) {
                    vector[pos] += variant.count;
                }
            }
        }
    }

    for base in bases {
        write!(file, "{}", base)?;
        if let Some(vector) = counts.get(&base) {
            for value in vector {
                write!(file, "\t{:.1}", *value as f64)?;
            }
        }
        writeln!(file)?;
    }

    Ok(())
}

pub fn write_run_info_json(
    path: &Path,
    refs: &[ReferenceInfo],
    stats: &AlignmentStats,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "{{")?;
    writeln!(file, "  \"running_info\": {{")?;
    writeln!(file, "    \"version\": \"0.1.0\"")?;
    writeln!(file, "  }},")?;
    writeln!(file, "  \"results\": {{")?;
    writeln!(
        file,
        "    \"ref_names\": [{}],",
        join_quoted(refs.iter().map(|r| r.name.as_str()).collect())
    )?;
    writeln!(file, "    \"alignment_stats\": {{")?;
    writeln!(file, "      \"N_TOT_READS\": {},", stats.n_tot_reads)?;
    writeln!(file, "      \"N_CACHED_ALN\": {},", stats.n_cached_aln)?;
    writeln!(
        file,
        "      \"N_CACHED_NOTALN\": {},",
        stats.n_cached_notaln
    )?;
    writeln!(file, "      \"N_COMPUTED_ALN\": {},", stats.n_computed_aln)?;
    writeln!(
        file,
        "      \"N_COMPUTED_NOTALN\": {},",
        stats.n_computed_notaln
    )?;
    writeln!(file, "      \"N_GLOBAL_SUBS\": {},", stats.n_global_subs)?;
    writeln!(
        file,
        "      \"N_SUBS_OUTSIDE_WINDOW\": {},",
        stats.n_subs_outside_window
    )?;
    writeln!(
        file,
        "      \"N_MODS_IN_WINDOW\": {},",
        stats.n_mods_in_window
    )?;
    writeln!(
        file,
        "      \"N_MODS_OUTSIDE_WINDOW\": {},",
        stats.n_mods_outside_window
    )?;
    writeln!(
        file,
        "      \"N_READS_IRREGULAR_ENDS\": {},",
        stats.n_reads_irregular_ends
    )?;
    writeln!(file, "      \"READ_LENGTH\": {}", stats.read_length)?;
    writeln!(file, "    }},")?;
    writeln!(file, "    \"refs\": {{")?;
    for (idx, ref_info) in refs.iter().enumerate() {
        writeln!(
            file,
            "      \"{}\": {{ \"sequence\": \"{}\", \"include_idxs\": [{}] }}{}",
            escape_json(&ref_info.name),
            escape_json(&ref_info.sequence),
            join_numbers(&ref_info.include_idxs),
            if idx + 1 == refs.len() { "" } else { "," }
        )?;
    }
    writeln!(file, "    }}")?;
    writeln!(file, "  }}")?;
    writeln!(file, "}}")?;
    Ok(())
}

fn join_quoted(values: Vec<&str>) -> String {
    values
        .into_iter()
        .map(|v| format!("\"{}\"", escape_json(v)))
        .collect::<Vec<_>>()
        .join(", ")
}

fn join_numbers(values: &[usize]) -> String {
    values
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

fn escape_json(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edits::EditPayload;
    use crate::models::VariantPayload;

    fn make_test_variant(seq: &str, count: usize, modified: bool) -> VariantRecord {
        let classification = if modified {
            ReadClassification::Modified
        } else {
            ReadClassification::Unmodified
        };
        let empty_edit = EditPayload {
            all_insertion_positions: vec![],
            all_insertion_left_positions: vec![],
            insertion_positions: vec![],
            insertion_coordinates: vec![],
            insertion_sizes: vec![],
            insertion_n: 0,
            all_deletion_positions: vec![],
            all_deletion_coordinates: vec![],
            deletion_positions: vec![],
            deletion_coordinates: vec![],
            deletion_sizes: vec![],
            deletion_n: 0,
            all_substitution_positions: vec![],
            substitution_positions: vec![],
            all_substitution_values: vec![],
            substitution_values: vec![],
            substitution_n: 0,
            ref_positions: vec![0, 1, 2, 3],
        };
        VariantRecord {
            count,
            aln_ref_names: vec!["Reference".to_string()],
            aln_scores: vec![95.0],
            ref_aln_details: vec![],
            best_match_score: 95.0,
            best_match_name: Some("Reference".to_string()),
            class_name: if modified {
                "Reference_MODIFIED".to_string()
            } else {
                "Reference_UNMODIFIED".to_string()
            },
            payloads: vec![VariantPayload {
                edit: empty_edit,
                ref_name: "Reference".to_string(),
                aln_scores: vec![95.0],
                aln_seq: seq.to_string(),
                aln_ref: seq.to_string(),
                aln_strand: '+',
                classification,
                irregular_ends: false,
                insertions_outside_window: 0,
                deletions_outside_window: 0,
                substitutions_outside_window: 0,
                total_mods: 0,
                mods_in_window: 0,
                mods_outside_window: 0,
            }],
        }
    }

    #[test]
    fn test_write_allele_frequency_table() {
        let dir = std::env::temp_dir().join("crispresso_test_output");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("alleles.tsv");
        let mut cache = HashMap::new();
        cache.insert("ATCG".to_string(), make_test_variant("ATCG", 100, false));
        cache.insert("ATAG".to_string(), make_test_variant("ATAG", 50, true));
        write_allele_frequency_table(&path, &cache).unwrap();
        assert!(path.exists());
        std::fs::remove_dir_all(&dir).ok();
    }
}

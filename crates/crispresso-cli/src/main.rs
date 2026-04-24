use std::env;
use std::error::Error;
use std::path::PathBuf;

use crispresso_core::*;

fn main() -> Result<(), Box<dyn Error>> {
    let args = parse_args(env::args().skip(1).collect())?;
    run(args)
}

fn run(args: Args) -> Result<(), Box<dyn Error>> {
    if args.fastq_r2.is_some() {
        return Err("--fastq_r2 paired-end processing is not implemented yet".into());
    }
    if args.bam_input.is_some() {
        return Err("--bam_input processing is not implemented yet".into());
    }
    if args.vcf_output {
        return Err("--vcf_output is not implemented yet".into());
    }
    let fastq_path = args.fastq_r1.as_ref().ok_or("Please provide --fastq_r1")?;
    let amplicon_seq = args
        .amplicon_seq
        .as_ref()
        .ok_or("Please provide --amplicon_seq")?
        .to_uppercase();

    let guides = parse_guides(
        args.guide_seq.as_deref(),
        args.guide_name.as_deref(),
        args.qw_center,
        args.qw_size,
    );

    let run_name = args.name.clone().unwrap_or_else(|| {
        fastq_path
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "crispresso_run".to_string())
    });
    let output_dir = args
        .output_folder
        .clone()
        .unwrap_or_else(|| PathBuf::from("."))
        .join(format!("CRISPResso_on_{}", run_name));
    std::fs::create_dir_all(&output_dir)?;

    let config = RunConfig {
        needleman_wunsch_gap_open: args.gap_open,
        needleman_wunsch_gap_extend: args.gap_extend,
        needleman_wunsch_gap_incentive: args.gap_incentive,
        aln_seed_len: args.seed_len,
        aln_seed_count: args.seed_count,
        default_min_aln_score: args.min_aln_score,
        exclude_bp_from_left: args.exclude_left,
        exclude_bp_from_right: args.exclude_right,
        quantification_window_center: args.qw_center,
        quantification_window_size: args.qw_size,
        ignore_deletions: args.ignore_deletions,
        ignore_insertions: args.ignore_insertions,
        ignore_substitutions: args.ignore_substitutions,
        use_legacy_insertion_quantification: args.use_legacy_insertion_quantification,
        ..RunConfig::default()
    };

    let matrix = make_matrix(5, -4, -2, -1);
    let guide_match = get_amplicon_info_for_guides(
        &amplicon_seq,
        &guides,
        None,
        config.exclude_bp_from_left,
        config.exclude_bp_from_right,
        config.plot_window_size,
    );
    let gap_incentive = build_gap_incentive(
        amplicon_seq.len(),
        &guide_match.sg_rna_cut_points,
        config.needleman_wunsch_gap_incentive,
    );
    let (fw_seeds, rc_seeds) = generate_seeds(
        &amplicon_seq,
        config.exclude_bp_from_left,
        config.exclude_bp_from_right,
        config.aln_seed_len,
        config.aln_seed_count,
    );

    let ref_info = ReferenceInfo {
        name: "Reference".to_string(),
        sequence: amplicon_seq.clone(),
        sequence_length: amplicon_seq.len(),
        min_aln_score: config.default_min_aln_score,
        gap_incentive,
        sg_rna_cut_points: guide_match.sg_rna_cut_points,
        sg_rna_intervals: guide_match.sg_rna_intervals,
        sg_rna_sequences: guide_match.sg_rna_sequences,
        sg_rna_names: guide_match.sg_rna_names,
        sg_rna_mismatches: guide_match.sg_rna_mismatches,
        sg_rna_orig_sequences: guides.iter().map(|g| g.orig_sequence.clone()).collect(),
        sg_rna_include_idxs: guide_match.sg_rna_include_idxs,
        contains_guide: !guides.is_empty(),
        include_idxs: guide_match.include_idxs,
        exclude_idxs: guide_match.exclude_idxs,
        fw_seeds,
        rc_seeds,
        exon_positions: Default::default(),
        exon_intervals: vec![],
    };

    let refs = vec![ref_info];
    let unique_reads = count_unique_reads(fastq_path)?;
    let reads_in_input: usize = unique_reads.values().sum();
    let (variant_cache, stats) = process_fastq_reads(&unique_reads, &refs, &matrix, &config);
    let reads_aligned_all_amplicons: usize = variant_cache.values().map(|v| v.count).sum();

    write_allele_frequency_table(
        &output_dir.join("Alleles_frequency_table.txt"),
        &variant_cache,
    )?;
    write_quantification_of_editing_frequency(
        &output_dir.join("CRISPResso_quantification_of_editing_frequency.txt"),
        &refs,
        &variant_cache,
        reads_in_input,
        reads_aligned_all_amplicons,
    )?;
    write_nucleotide_frequency_table(
        &output_dir.join("Nucleotide_frequency_table.txt"),
        &variant_cache,
        &refs[0],
    )?;
    write_run_info_json(&output_dir.join("CRISPResso2_info.json"), &refs, &stats)?;

    Ok(())
}

#[derive(Debug, Default)]
struct Args {
    fastq_r1: Option<PathBuf>,
    fastq_r2: Option<PathBuf>,
    bam_input: Option<PathBuf>,
    amplicon_seq: Option<String>,
    guide_seq: Option<String>,
    guide_name: Option<String>,
    name: Option<String>,
    output_folder: Option<PathBuf>,
    gap_open: i32,
    gap_extend: i32,
    gap_incentive: i64,
    qw_center: i32,
    qw_size: i32,
    min_aln_score: f64,
    exclude_left: usize,
    exclude_right: usize,
    seed_len: usize,
    seed_count: usize,
    ignore_deletions: bool,
    ignore_insertions: bool,
    ignore_substitutions: bool,
    use_legacy_insertion_quantification: bool,
    vcf_output: bool,
    amplicon_coordinates: Option<String>,
}

fn parse_args(argv: Vec<String>) -> Result<Args, Box<dyn Error>> {
    let mut args = Args {
        gap_open: -20,
        gap_extend: -2,
        gap_incentive: 1,
        qw_center: -3,
        qw_size: 1,
        min_aln_score: 60.0,
        exclude_left: 0,
        exclude_right: 0,
        seed_len: 10,
        seed_count: 5,
        ..Args::default()
    };

    let mut i = 0usize;
    while i < argv.len() {
        match argv[i].as_str() {
            "-r" | "--fastq_r1" => {
                i += 1;
                args.fastq_r1 = argv.get(i).map(PathBuf::from);
            }
            "-r2" | "--fastq_r2" => {
                i += 1;
                args.fastq_r2 = argv.get(i).map(PathBuf::from);
            }
            "--bam_input" => {
                i += 1;
                args.bam_input = argv.get(i).map(PathBuf::from);
            }
            "-a" | "--amplicon_seq" => {
                i += 1;
                args.amplicon_seq = argv.get(i).cloned();
            }
            "-g" | "--guide_seq" => {
                i += 1;
                args.guide_seq = argv.get(i).cloned();
            }
            "--guide_name" => {
                i += 1;
                args.guide_name = argv.get(i).cloned();
            }
            "-n" | "--name" => {
                i += 1;
                args.name = argv.get(i).cloned();
            }
            "--output_folder" => {
                i += 1;
                args.output_folder = argv.get(i).map(PathBuf::from);
            }
            "--needleman_wunsch_gap_open" => {
                i += 1;
                args.gap_open = argv.get(i).ok_or("missing gap_open")?.parse()?;
            }
            "--needleman_wunsch_gap_extend" => {
                i += 1;
                args.gap_extend = argv.get(i).ok_or("missing gap_extend")?.parse()?;
            }
            "--needleman_wunsch_gap_incentive" => {
                i += 1;
                args.gap_incentive = argv.get(i).ok_or("missing gap_incentive")?.parse()?;
            }
            "--quantification_window_center" => {
                i += 1;
                args.qw_center = argv.get(i).ok_or("missing qw_center")?.parse()?;
            }
            "--quantification_window_size" => {
                i += 1;
                args.qw_size = argv.get(i).ok_or("missing qw_size")?.parse()?;
            }
            "--default_min_aln_score" => {
                i += 1;
                args.min_aln_score = argv.get(i).ok_or("missing min score")?.parse()?;
            }
            "--exclude_bp_from_left" => {
                i += 1;
                args.exclude_left = argv.get(i).ok_or("missing exclude_left")?.parse()?;
            }
            "--exclude_bp_from_right" => {
                i += 1;
                args.exclude_right = argv.get(i).ok_or("missing exclude_right")?.parse()?;
            }
            "--aln_seed_len" => {
                i += 1;
                args.seed_len = argv.get(i).ok_or("missing seed_len")?.parse()?;
            }
            "--aln_seed_count" => {
                i += 1;
                args.seed_count = argv.get(i).ok_or("missing seed_count")?.parse()?;
            }
            "--ignore_deletions" => args.ignore_deletions = true,
            "--ignore_insertions" => args.ignore_insertions = true,
            "--ignore_substitutions" => args.ignore_substitutions = true,
            "--use_legacy_insertion_quantification" => {
                args.use_legacy_insertion_quantification = true
            }
            "--vcf_output" => args.vcf_output = true,
            "--amplicon_coordinates" => {
                i += 1;
                args.amplicon_coordinates = argv.get(i).cloned();
            }
            _ => {}
        }
        i += 1;
    }

    Ok(args)
}

fn parse_guides(
    guide_seq: Option<&str>,
    guide_name: Option<&str>,
    qw_center: i32,
    qw_size: i32,
) -> Vec<GuideInfo> {
    let seqs: Vec<&str> = guide_seq
        .map(|s| s.split(',').collect())
        .unwrap_or_default();
    let names: Vec<&str> = guide_name
        .map(|s| s.split(',').collect())
        .unwrap_or_default();

    seqs.into_iter()
        .enumerate()
        .map(|(i, seq)| GuideInfo {
            sequence: seq.to_uppercase(),
            orig_sequence: seq.to_uppercase(),
            name: names.get(i).map(|s| s.to_string()).unwrap_or_default(),
            qw_center,
            qw_size,
            mismatches: vec![],
            plot_cut_point: true,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn parse_minimal_args() {
        let args = parse_args(vec![
            "-r".to_string(),
            "input.fastq".to_string(),
            "-a".to_string(),
            "ATCG".to_string(),
        ])
        .unwrap();
        assert_eq!(args.fastq_r1, Some(PathBuf::from("input.fastq")));
        assert_eq!(args.amplicon_seq.as_deref(), Some("ATCG"));
    }

    #[test]
    fn parse_guides_basic() {
        let guides = parse_guides(Some("AAAA,CCCC"), Some("g1,g2"), -3, 1);
        assert_eq!(guides.len(), 2);
        assert_eq!(guides[0].name, "g1");
        assert_eq!(guides[1].sequence, "CCCC");
    }

    #[test]
    fn parse_extended_parameter_surface() {
        let args = parse_args(vec![
            "--fastq_r1".to_string(),
            "input.fastq".to_string(),
            "--fastq_r2".to_string(),
            "input_r2.fastq".to_string(),
            "--bam_input".to_string(),
            "input.bam".to_string(),
            "--use_legacy_insertion_quantification".to_string(),
            "--vcf_output".to_string(),
            "--amplicon_coordinates".to_string(),
            "chr1:100".to_string(),
        ])
        .unwrap();
        assert_eq!(args.fastq_r2, Some(PathBuf::from("input_r2.fastq")));
        assert_eq!(args.bam_input, Some(PathBuf::from("input.bam")));
        assert!(args.use_legacy_insertion_quantification);
        assert!(args.vcf_output);
        assert_eq!(args.amplicon_coordinates.as_deref(), Some("chr1:100"));
    }

    #[test]
    fn fanc_cas9_golden_output_matches_expected() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .to_path_buf();
        let out_dir =
            std::env::temp_dir().join(format!("crispresso_fanc_golden_{}", std::process::id()));
        let _ = fs::remove_dir_all(&out_dir);
        fs::create_dir_all(&out_dir).unwrap();

        run(
            parse_args(vec![
                "--fastq_r1".to_string(),
                repo_root
                    .join("CRISPResso2-master/tests/FANC.Cas9.fastq")
                    .display()
                    .to_string(),
                "--amplicon_seq".to_string(),
                "CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG".to_string(),
                "--guide_seq".to_string(),
                "GGAATCCCTTCTGCAGCACC".to_string(),
                "--output_folder".to_string(),
                out_dir.display().to_string(),
            ])
            .unwrap(),
        )
        .unwrap();

        let actual_dir = out_dir.join("CRISPResso_on_FANC.Cas9");
        let expected_dir =
            repo_root.join("CRISPResso2-master/tests/expectedResults/CRISPResso_on_FANC.Cas9");
        let quant_actual = fs::read_to_string(
            actual_dir.join("CRISPResso_quantification_of_editing_frequency.txt"),
        )
        .unwrap();
        let quant_expected = fs::read_to_string(
            expected_dir.join("CRISPResso_quantification_of_editing_frequency.txt"),
        )
        .unwrap();
        assert_eq!(quant_actual, quant_expected);

        let nuc_actual =
            fs::read_to_string(actual_dir.join("Nucleotide_frequency_table.txt")).unwrap();
        let nuc_expected =
            fs::read_to_string(expected_dir.join("Nucleotide_frequency_table.txt")).unwrap();
        assert_eq!(nuc_actual, nuc_expected);

        let _ = fs::remove_dir_all(&out_dir);
    }

    #[test]
    fn run_rejects_unimplemented_paired_end_and_bam_and_vcf_modes() {
        let paired_err = run(parse_args(vec![
            "--fastq_r1".to_string(),
            "input.fastq".to_string(),
            "--fastq_r2".to_string(),
            "input_r2.fastq".to_string(),
            "--amplicon_seq".to_string(),
            "ATCG".to_string(),
        ])
        .unwrap())
        .unwrap_err()
        .to_string();
        assert!(paired_err.contains("--fastq_r2"));

        let bam_err = run(parse_args(vec![
            "--bam_input".to_string(),
            "input.bam".to_string(),
            "--amplicon_seq".to_string(),
            "ATCG".to_string(),
        ])
        .unwrap())
        .unwrap_err()
        .to_string();
        assert!(bam_err.contains("--bam_input"));

        let vcf_err = run(parse_args(vec![
            "--fastq_r1".to_string(),
            "input.fastq".to_string(),
            "--amplicon_seq".to_string(),
            "ATCG".to_string(),
            "--vcf_output".to_string(),
        ])
        .unwrap())
        .unwrap_err()
        .to_string();
        assert!(vcf_err.contains("--vcf_output"));
    }
}

use std::{error::Error, fmt, fs, path::Path};

#[derive(Debug)]
pub enum AlignmentError {
    EmptyMatrixFile,
    InvalidMatrixHeader,
    InvalidMatrixValue(String),
    GapIncentiveLength { expected: usize, actual: usize },
    EmptySequence,
    MatrixOutOfBounds(u8, u8),
    Io(std::io::Error),
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMatrixFile => write!(f, "alignment matrix file is empty"),
            Self::InvalidMatrixHeader => write!(f, "invalid matrix header"),
            Self::InvalidMatrixValue(value) => write!(f, "invalid matrix value '{value}'"),
            Self::GapIncentiveLength { expected, actual } => {
                write!(
                    f,
                    "gap_incentive length mismatch: expected {expected}, got {actual}"
                )
            }
            Self::EmptySequence => write!(f, "global alignment requires non-empty sequences"),
            Self::MatrixOutOfBounds(row, col) => {
                write!(
                    f,
                    "matrix does not contain a score for byte values {row} and {col}"
                )
            }
            Self::Io(err) => write!(f, "{err}"),
        }
    }
}

impl Error for AlignmentError {}

impl From<std::io::Error> for AlignmentError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentResult {
    pub aligned_read: String,
    pub aligned_reference: String,
    pub score: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScoreMatrix {
    scores: Vec<Vec<i64>>,
}

impl ScoreMatrix {
    pub fn new(size: usize) -> Self {
        Self {
            scores: vec![vec![0; size]; size],
        }
    }

    pub fn set(&mut self, row: u8, col: u8, value: i64) {
        let row = row as usize;
        let col = col as usize;
        if row >= self.scores.len() || col >= self.scores.len() {
            let new_size = row.max(col) + 1;
            for existing_row in &mut self.scores {
                existing_row.resize(new_size, 0);
            }
            self.scores.resize(new_size, vec![0; new_size]);
        }
        self.scores[row][col] = value;
    }

    fn get(&self, row: u8, col: u8) -> Result<i32, AlignmentError> {
        self.scores
            .get(row as usize)
            .and_then(|r| r.get(col as usize))
            .copied()
            .map(|v| v as i32)
            .ok_or(AlignmentError::MatrixOutOfBounds(row, col))
    }
}

pub fn make_matrix(
    match_score: i64,
    mismatch_score: i64,
    n_mismatch_score: i64,
    n_match_score: i64,
) -> ScoreMatrix {
    let letters = [b'A', b'T', b'C', b'G', b'N'];
    let max_header = *letters.iter().max().unwrap() as usize;
    let mut matrix = ScoreMatrix::new(max_header + 1);
    let nucs = [b'A', b'T', b'C', b'G'];

    for &nuc in &nucs {
        for &other in &nucs {
            matrix.set(
                nuc,
                other,
                if nuc == other {
                    match_score
                } else {
                    mismatch_score
                },
            );
        }
        matrix.set(nuc, b'N', n_mismatch_score);
        matrix.set(b'N', nuc, n_mismatch_score);
    }
    matrix.set(b'N', b'N', n_match_score);
    matrix
}

pub fn read_matrix(path: impl AsRef<Path>) -> Result<ScoreMatrix, AlignmentError> {
    let contents = fs::read_to_string(path)?;
    let mut lines = contents.lines();
    let headers = loop {
        let line = lines.next().ok_or(AlignmentError::EmptyMatrixFile)?.trim();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let headers = line
            .split(' ')
            .filter(|value| !value.is_empty())
            .map(|value| {
                value
                    .as_bytes()
                    .first()
                    .copied()
                    .ok_or(AlignmentError::InvalidMatrixHeader)
            })
            .collect::<Result<Vec<_>, _>>()?;
        break headers;
    };

    let max_header = headers
        .iter()
        .copied()
        .max()
        .ok_or(AlignmentError::InvalidMatrixHeader)? as usize;
    let mut matrix = ScoreMatrix::new(max_header + 1);

    for (row_idx, line) in lines.enumerate() {
        let values = line
            .split(' ')
            .skip(1)
            .filter(|value| !value.is_empty())
            .map(|value| {
                value
                    .parse::<i64>()
                    .map_err(|_| AlignmentError::InvalidMatrixValue(value.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if row_idx >= headers.len() {
            break;
        }
        for (&col_header, value) in headers.iter().zip(values) {
            matrix.set(headers[row_idx], col_header, value);
        }
    }

    Ok(matrix)
}

pub fn global_align(
    read: &str,
    reference: &str,
    matrix: &ScoreMatrix,
    gap_incentive: &[i64],
    gap_open: i32,
    gap_extend: i32,
) -> Result<AlignmentResult, AlignmentError> {
    let max_j = read.len();
    let max_i = reference.len();
    if max_i == 0 || max_j == 0 {
        return Err(AlignmentError::EmptySequence);
    }
    if gap_incentive.len() != max_i + 1 {
        return Err(AlignmentError::GapIncentiveLength {
            expected: max_i + 1,
            actual: gap_incentive.len(),
        });
    }

    let seqj = read.as_bytes();
    let seqi = reference.as_bytes();
    let min_score = gap_open * max_j as i32 * max_i as i32;

    let mut m_score = vec![vec![0_i32; max_j + 1]; max_i + 1];
    let mut i_score = vec![vec![0_i32; max_j + 1]; max_i + 1];
    let mut j_score = vec![vec![0_i32; max_j + 1]; max_i + 1];
    let mut m_pointer = vec![vec![0_u8; max_j + 1]; max_i + 1];
    let mut i_pointer = vec![vec![0_u8; max_j + 1]; max_i + 1];
    let mut j_pointer = vec![vec![0_u8; max_j + 1]; max_i + 1];

    const MARRAY: u8 = 1;
    const IARRAY: u8 = 2;
    const JARRAY: u8 = 3;

    for j in 1..=max_j {
        m_score[0][j] = min_score;
        m_pointer[0][j] = IARRAY;
    }
    for i in 1..=max_i {
        m_score[i][0] = min_score;
        m_pointer[i][0] = JARRAY;
    }

    for j in 1..=max_j {
        i_score[0][j] = gap_extend * j as i32 + gap_incentive[0] as i32;
        i_pointer[0][j] = IARRAY;
    }
    for row in i_score.iter_mut().take(max_i + 1) {
        row[0] = min_score;
    }

    for i in 1..=max_i {
        j_score[i][0] = gap_extend * i as i32 + gap_incentive[0] as i32;
        j_pointer[i][0] = JARRAY;
    }
    for j in 0..=max_j {
        j_score[0][j] = min_score;
    }

    for i in 1..max_i {
        let ci = seqi[i - 1];
        for j in 1..max_j {
            let cj = seqj[j - 1];
            fill_cell(
                i,
                j,
                ci,
                cj,
                matrix,
                gap_incentive,
                gap_open,
                gap_extend,
                &mut m_score,
                &mut i_score,
                &mut j_score,
                &mut m_pointer,
                &mut i_pointer,
                &mut j_pointer,
            )?;
        }
    }

    let j = max_j;
    let cj = seqj[j - 1];
    for i in 1..max_i {
        let ci = seqi[i - 1];
        fill_cell(
            i,
            j,
            ci,
            cj,
            matrix,
            gap_incentive,
            gap_extend,
            gap_extend,
            &mut m_score,
            &mut i_score,
            &mut j_score,
            &mut m_pointer,
            &mut i_pointer,
            &mut j_pointer,
        )?;
    }

    let i = max_i;
    let ci = seqi[i - 1];
    for j in 1..=max_j {
        let cj = seqj[j - 1];
        fill_cell(
            i,
            j,
            ci,
            cj,
            matrix,
            gap_incentive,
            gap_extend,
            gap_extend,
            &mut m_score,
            &mut i_score,
            &mut j_score,
            &mut m_pointer,
            &mut i_pointer,
            &mut j_pointer,
        )?;
    }

    let mut i = max_i;
    let mut j = max_j;
    let mut ci = seqi[i - 1];
    let mut cj = seqj[j - 1];
    let mut curr_matrix = if m_score[i][j] > j_score[i][j] {
        if m_score[i][j] > i_score[i][j] {
            MARRAY
        } else {
            IARRAY
        }
    } else if j_score[i][j] > i_score[i][j] {
        JARRAY
    } else {
        IARRAY
    };

    let mut aligned_read = Vec::with_capacity(max_i + max_j);
    let mut aligned_reference = Vec::with_capacity(max_i + max_j);
    let mut match_count = 0_usize;

    while i > 0 || j > 0 {
        if curr_matrix == MARRAY {
            curr_matrix = m_pointer[i][j];
            aligned_read.push(cj);
            aligned_reference.push(ci);
            if cj == ci {
                match_count += 1;
            }
            if i > 1 {
                i -= 1;
                ci = seqi[i - 1];
            } else {
                i = 0;
                ci = seqi[0];
            }
            if j > 1 {
                j -= 1;
                cj = seqj[j - 1];
            } else {
                j = 0;
                cj = seqj[0];
            }
        } else if curr_matrix == JARRAY {
            curr_matrix = j_pointer[i][j];
            aligned_read.push(b'-');
            aligned_reference.push(ci);
            if i > 1 {
                i -= 1;
                ci = seqi[i - 1];
            } else {
                i = 0;
                ci = seqi[0];
            }
        } else if curr_matrix == IARRAY {
            curr_matrix = i_pointer[i][j];
            aligned_read.push(cj);
            aligned_reference.push(b'-');
            if j > 1 {
                j -= 1;
                cj = seqj[j - 1];
            } else {
                j = 0;
                cj = seqj[0];
            }
        } else {
            return Err(AlignmentError::InvalidMatrixHeader);
        }
    }

    aligned_read.reverse();
    aligned_reference.reverse();
    let align_len = aligned_read.len();
    let raw_score = 100.0 * match_count as f64 / align_len as f64;

    Ok(AlignmentResult {
        aligned_read: String::from_utf8(aligned_read).expect("ASCII alignment"),
        aligned_reference: String::from_utf8(aligned_reference).expect("ASCII alignment"),
        score: (raw_score * 1000.0).round() / 1000.0,
    })
}

#[allow(clippy::too_many_arguments)]
fn fill_cell(
    i: usize,
    j: usize,
    ci: u8,
    cj: u8,
    matrix: &ScoreMatrix,
    gap_incentive: &[i64],
    gap_open_for_new: i32,
    gap_extend: i32,
    m_score: &mut [Vec<i32>],
    i_score: &mut [Vec<i32>],
    j_score: &mut [Vec<i32>],
    m_pointer: &mut [Vec<u8>],
    i_pointer: &mut [Vec<u8>],
    j_pointer: &mut [Vec<u8>],
) -> Result<(), AlignmentError> {
    const MARRAY: u8 = 1;
    const IARRAY: u8 = 2;
    const JARRAY: u8 = 3;

    let i_from_m = gap_open_for_new + m_score[i][j - 1] + gap_incentive[i] as i32;
    let i_extend = gap_extend + i_score[i][j - 1] + gap_incentive[i] as i32;
    if i_from_m > i_extend {
        i_score[i][j] = i_from_m;
        i_pointer[i][j] = MARRAY;
    } else {
        i_score[i][j] = i_extend;
        i_pointer[i][j] = IARRAY;
    }

    let j_from_m = gap_open_for_new + m_score[i - 1][j] + gap_incentive[i - 1] as i32;
    let j_extend = gap_extend + j_score[i - 1][j];
    if j_from_m > j_extend {
        j_score[i][j] = j_from_m;
        j_pointer[i][j] = MARRAY;
    } else {
        j_score[i][j] = j_extend;
        j_pointer[i][j] = JARRAY;
    }

    let subst = matrix.get(ci, cj)?;
    let m_val = m_score[i - 1][j - 1] + subst;
    let i_val = i_score[i - 1][j - 1] + subst;
    let j_val = j_score[i - 1][j - 1] + subst;

    if m_val > j_val {
        if m_val > i_val {
            m_score[i][j] = m_val;
            m_pointer[i][j] = MARRAY;
        } else {
            m_score[i][j] = i_val;
            m_pointer[i][j] = IARRAY;
        }
    } else if j_val > i_val {
        m_score[i][j] = j_val;
        m_pointer[i][j] = JARRAY;
    } else {
        m_score[i][j] = i_val;
        m_pointer[i][j] = IARRAY;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln_matrix() -> ScoreMatrix {
        make_matrix(5, -4, -2, -1)
    }

    #[test]
    fn global_align_perfect_match() {
        let result =
            global_align("ATTA", "ATTA", &aln_matrix(), &[0, 0, 0, 0, 0], -20, -2).unwrap();
        assert_eq!(result.aligned_read, "ATTA");
        assert_eq!(result.aligned_reference, "ATTA");
        assert_eq!(result.score, 100.0);
    }

    #[test]
    fn global_align_with_n_matches_python_behavior() {
        let result =
            global_align("ATNA", "ATTA", &aln_matrix(), &[0, 0, 0, 0, 0], -20, -2).unwrap();
        assert_eq!(result.aligned_read, "ATNA");
        assert_eq!(result.aligned_reference, "ATTA");
        assert_eq!(result.score, 75.0);
    }

    // NOTE: ANNG vs ATCG with make_matrix(5,-4,-2,-1): both produce same score for
    // ANNG/ATCG (score=6, mScore wins) so no gap is introduced. The Python unit test
    // test_global_align_with_n uses EDNAFULL matrix (not make_matrix), which changes
    // the tie-breaking and produces A-NNG/ATC-G. This test uses make_matrix behavior.
    #[test]
    fn global_align_with_n_no_gap_make_matrix_behavior() {
        let result =
            global_align("ANNG", "ATCG", &aln_matrix(), &[0, 0, 0, 0, 0], -20, -2).unwrap();
        assert_eq!(result.aligned_read, "ANNG");
        assert_eq!(result.aligned_reference, "ATCG");
        assert_eq!(result.score, 50.0);
    }

    // gap_incentive test: with make_matrix, ATTTA vs ATTA ends up at ATTA- (trailing gap).
    // The Python test with EDNAFULL matrix + gap_incentive=[0,0,1,0,0] gives AT-TA.
    #[test]
    fn global_align_gap_incentive_prefers_cut_site() {
        let result =
            global_align("ATTTA", "ATTA", &aln_matrix(), &[0, 0, 1, 0, 0], -20, -2).unwrap();
        assert_eq!(result.aligned_read, "ATTTA");
        assert_eq!(result.aligned_reference, "ATTA-");
        assert_eq!(result.score, 60.0);
    }

    #[test]
    fn global_align_completely_different_matches_python_behavior() {
        let result =
            global_align("AAAA", "TTTT", &aln_matrix(), &[0, 0, 0, 0, 0], -20, -2).unwrap();
        assert_eq!(result.aligned_read, "---AAAA");
        assert_eq!(result.aligned_reference, "TTTT---");
        assert_eq!(result.score, 0.0);
    }

    #[test]
    fn matrix_rejects_bad_gap_incentive_len() {
        let err = global_align("ATTA", "ATTA", &aln_matrix(), &[0], -20, -2).unwrap_err();
        assert!(matches!(err, AlignmentError::GapIncentiveLength { .. }));
    }

    // -----------------------------------------------------------------------
    // EDNAFULL-based tests mirroring Python test_CRISPResso2Align.py exactly
    // (default gap_open=-1, gap_extend=-1 as in Python global_align signature)
    // -----------------------------------------------------------------------

    fn ednafull() -> ScoreMatrix {
        read_matrix(concat!(env!("CARGO_MANIFEST_DIR"), "/EDNAFULL")).expect("EDNAFULL matrix")
    }

    #[test]
    fn ednafull_perfect_match() {
        let m = ednafull();
        let r = global_align("ATTA", "ATTA", &m, &[0, 0, 0, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "ATTA");
        assert_eq!(r.aligned_reference, "ATTA");
        assert_eq!(r.score, 100.0);
    }

    #[test]
    fn ednafull_with_n_produces_gap() {
        // Python test_global_align_with_n: uses EDNAFULL, gap_open=-1, gap_extend=-1
        let m = ednafull();
        let r = global_align("ANNG", "ATCG", &m, &[0, 0, 0, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "A-NNG");
        assert_eq!(r.aligned_reference, "ATC-G");
        assert_eq!(r.score, 40.0);
    }

    #[test]
    fn ednafull_gap_incentive_pos1() {
        // Python test_global_align_gap_incentive_s1: gap_incentive=[0,1,0,0,0]
        let m = ednafull();
        let r = global_align("ATTTA", "ATTA", &m, &[0, 1, 0, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "ATTTA");
        assert_eq!(r.aligned_reference, "A-TTA");
        assert_eq!(r.score, (100.0 * 4.0 / 5.0_f64 * 1000.0).round() / 1000.0);
    }

    #[test]
    fn ednafull_gap_incentive_pos2() {
        // Python: gap_incentive=[0,0,1,0,0] => AT-TA
        let m = ednafull();
        let r = global_align("ATTTA", "ATTA", &m, &[0, 0, 1, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "ATTTA");
        assert_eq!(r.aligned_reference, "AT-TA");
        assert_eq!(r.score, (100.0 * 4.0 / 5.0_f64 * 1000.0).round() / 1000.0);
    }

    #[test]
    fn ednafull_homopolymer_no_incentive() {
        // Python: TTTTT vs TTTT, no incentive => trailing gap
        let m = ednafull();
        let r = global_align("TTTTT", "TTTT", &m, &[0, 0, 0, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "TTTTT");
        assert_eq!(r.aligned_reference, "TTTT-");
        assert_eq!(r.score, (100.0 * 4.0 / 5.0_f64 * 1000.0).round() / 1000.0);
    }

    #[test]
    fn ednafull_homopolymer_incentive_pos0() {
        // Python: TTTTT vs TTTT, gap_incentive=[1,0,0,0,0] => leading gap
        let m = ednafull();
        let r = global_align("TTTTT", "TTTT", &m, &[1, 0, 0, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "TTTTT");
        assert_eq!(r.aligned_reference, "-TTTT");
        assert_eq!(r.score, (100.0 * 4.0 / 5.0_f64 * 1000.0).round() / 1000.0);
    }

    #[test]
    fn ednafull_homopolymer_incentive_pos2() {
        // Python: TTTTT vs TTTT, gap_incentive=[0,0,1,0,0] => TT-TT
        let m = ednafull();
        let r = global_align("TTTTT", "TTTT", &m, &[0, 0, 1, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "TTTTT");
        assert_eq!(r.aligned_reference, "TT-TT");
        assert_eq!(r.score, (100.0 * 4.0 / 5.0_f64 * 1000.0).round() / 1000.0);
    }

    #[test]
    fn ednafull_completely_different() {
        let m = ednafull();
        let r = global_align("AAAA", "TTTT", &m, &[0, 0, 0, 0, 0], -1, -1).unwrap();
        assert_eq!(r.aligned_read, "---AAAA");
        assert_eq!(r.aligned_reference, "TTTT---");
        assert_eq!(r.score, 0.0);
    }
}

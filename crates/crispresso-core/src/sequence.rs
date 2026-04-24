pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        'a' => 't',
        't' => 'a',
        'c' => 'g',
        'g' => 'c',
        'n' => 'n',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_preserves_crispresso_basic_behavior() {
        assert_eq!(reverse_complement("ATCGN"), "NCGAT");
        assert_eq!(reverse_complement("atcgn"), "ncgat");
    }
}

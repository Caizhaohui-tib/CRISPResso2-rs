use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub id: String,
    pub seq: String,
    pub qual: String,
}

pub fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>, std::io::Error> {
    let file = File::open(path)?;
    Ok(Box::new(BufReader::new(file)))
}

pub struct FastqReader {
    reader: Box<dyn BufRead>,
}

impl FastqReader {
    pub fn from_path(path: &Path) -> Result<Self, std::io::Error> {
        Ok(Self {
            reader: open_fastq(path)?,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<FastqRecord>, std::io::Error> {
        let mut id_line = String::new();
        let mut seq_line = String::new();
        let mut plus_line = String::new();
        let mut qual_line = String::new();

        loop {
            id_line.clear();
            if self.reader.read_line(&mut id_line)? == 0 {
                return Ok(None);
            }
            if !id_line.trim().is_empty() {
                break;
            }
        }

        if self.reader.read_line(&mut seq_line)? == 0 {
            return Ok(None);
        }
        if self.reader.read_line(&mut plus_line)? == 0 {
            return Ok(None);
        }
        if self.reader.read_line(&mut qual_line)? == 0 {
            return Ok(None);
        }

        Ok(Some(FastqRecord {
            id: id_line.trim_end().to_string(),
            seq: seq_line.trim_end().to_string(),
            qual: qual_line.trim_end().to_string(),
        }))
    }
}

pub fn count_unique_reads(path: &Path) -> Result<HashMap<String, usize>, std::io::Error> {
    let mut reader = FastqReader::from_path(path)?;
    let mut cache = HashMap::new();
    while let Some(record) = reader.read_record()? {
        *cache.entry(record.seq).or_insert(0) += 1;
    }
    Ok(cache)
}

pub fn count_total_reads(path: &Path) -> Result<usize, std::io::Error> {
    let mut reader = FastqReader::from_path(path)?;
    let mut count = 0usize;
    while reader.read_record()?.is_some() {
        count += 1;
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_test_fastq(path: &Path, records: &[(&str, &str, &str)], extra_newlines: &str) {
        let mut f = File::create(path).unwrap();
        for (id, seq, qual) in records {
            writeln!(f, "{}", id).unwrap();
            writeln!(f, "{}", seq).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", qual).unwrap();
        }
        write!(f, "{}", extra_newlines).unwrap();
    }

    #[test]
    fn test_fastq_counts() {
        let dir = std::env::temp_dir().join("crispresso_test_fastq_counts");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.fastq");

        write_test_fastq(
            &path,
            &[
                ("@r1", "ATCG", "IIII"),
                ("@r2", "GCTA", "IIII"),
                ("@r3", "ATCG", "IIII"),
            ],
            "\n\n",
        );

        let counts = count_unique_reads(&path).unwrap();
        assert_eq!(counts.get("ATCG"), Some(&2));
        assert_eq!(counts.get("GCTA"), Some(&1));
        assert_eq!(count_total_reads(&path).unwrap(), 3);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_empty_fastq() {
        let dir = std::env::temp_dir().join("crispresso_test_fastq_empty");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("empty.fastq");
        File::create(&path).unwrap();
        assert_eq!(count_total_reads(&path).unwrap(), 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_fastq_counts_no_trailing_newline() {
        let dir = std::env::temp_dir().join("crispresso_test_fastq_no_newline");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.fastq");
        let mut f = File::create(&path).unwrap();
        write!(f, "@r1\nATCG\n+\nIIII").unwrap();
        drop(f);

        assert_eq!(count_total_reads(&path).unwrap(), 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_fastq_counts_many_extra_newlines() {
        let dir = std::env::temp_dir().join("crispresso_test_fastq_newlines");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.fastq");
        write_test_fastq(&path, &[("@r1", "ATCG", "IIII")], "\n\n\n\n\n\n\n\n");

        assert_eq!(count_total_reads(&path).unwrap(), 1);
        std::fs::remove_dir_all(&dir).ok();
    }
}

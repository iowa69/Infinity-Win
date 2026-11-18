//! FASTA reading (section 6.2.1 and section 9).
//!
//! Plain and gzip input. One file is one genome; each record is a contig.
//! Sequences are upper-cased and every character that is not A, C, G or T
//! becomes N before anything else looks at them.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::seq::normalize_byte;

/// One contig of a genome, already normalised.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub seq: Vec<u8>,
}

impl Contig {
    pub fn len(&self) -> u64 {
        self.seq.len() as u64
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// A whole genome in memory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeData {
    /// Genome name: the file stem, with a trailing `.fasta`/`.fa`/`.fna` removed.
    pub name: String,
    pub contigs: Vec<Contig>,
}

impl GenomeData {
    pub fn total_length(&self) -> u64 {
        self.contigs.iter().map(|c| c.len()).sum()
    }
}

/// The gzip magic bytes, used to detect compression from content rather than
/// from the file extension: real download pipelines produce `.fasta.gz`,
/// `.fa.gz` and plain `.gz` interchangeably.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Read a genome from a path, detecting gzip by magic bytes.
pub fn read_genome(path: &Path) -> io::Result<GenomeData> {
    let file = File::open(path)?;
    let mut peek = BufReader::new(file);
    let head = peek.fill_buf()?;
    let gz = head.len() >= 2 && head[0] == GZIP_MAGIC[0] && head[1] == GZIP_MAGIC[1];
    let name = genome_name(path);
    if gz {
        let dec = MultiGzDecoder::new(peek);
        parse_fasta(dec, &name)
    } else {
        parse_fasta(peek, &name)
    }
}

/// Derive a genome name from a file name: strip the directory and every
/// sequence extension, including a doubled one such as `.fasta.gz`.
pub fn genome_name(path: &Path) -> String {
    let file = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned());
    let mut stem = file.as_str();
    loop {
        let lower = stem.to_ascii_lowercase();
        let cut = [".gz", ".fa", ".fasta", ".fna", ".fas", ".seq"]
            .iter()
            .find_map(|ext| {
                if lower.ends_with(ext) && stem.len() > ext.len() {
                    Some(stem.len() - ext.len())
                } else {
                    None
                }
            });
        match cut {
            Some(at) => stem = &stem[..at],
            None => break,
        }
    }
    if stem.is_empty() {
        file
    } else {
        stem.to_string()
    }
}

/// Parse FASTA from any byte stream. Lines starting with `>` begin a record;
/// the record name is the first whitespace-delimited token after `>`.
fn parse_fasta<R: Read>(reader: R, genome_name: &str) -> io::Result<GenomeData> {
    let mut buf = BufReader::new(reader);
    let mut contigs: Vec<Contig> = Vec::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    let mut line = Vec::new();

    loop {
        line.clear();
        let got = buf.read_until(b'\n', &mut line)?;
        if got == 0 {
            break;
        }
        let trimmed = trim_eol(&line);
        if trimmed.is_empty() {
            continue;
        }
        if trimmed[0] == b';' {
            continue; // legacy FASTA comment
        }
        if trimmed[0] == b'>' {
            if let Some(n) = name.take() {
                contigs.push(Contig {
                    name: n,
                    seq: std::mem::take(&mut seq),
                });
            }
            name = Some(parse_record_name(trimmed));
            continue;
        }
        if name.is_none() {
            // Sequence before any header: treat the file as a single unnamed
            // contig rather than failing, so odd but usable inputs still work.
            name = Some(format!("{genome_name}_contig1"));
        }
        for &b in trimmed {
            if b.is_ascii_whitespace() {
                continue;
            }
            seq.push(normalize_byte(b));
        }
    }
    if let Some(n) = name.take() {
        contigs.push(Contig {
            name: n,
            seq: std::mem::take(&mut seq),
        });
    }
    // An empty record is kept (it costs one contig and no bases) so that the
    // contig numbering in the output still matches the input file.
    Ok(GenomeData {
        name: genome_name.to_string(),
        contigs,
    })
}

fn trim_eol(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}

fn parse_record_name(header: &[u8]) -> String {
    let body = &header[1..];
    let end = body
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(body.len());
    let n = String::from_utf8_lossy(&body[..end]).into_owned();
    if n.is_empty() {
        "unnamed".to_string()
    } else {
        n
    }
}

/// The file extensions accepted as input (section 9).
pub const INPUT_EXTENSIONS: &[&str] = &["fa", "fasta", "fna", "gz"];

/// Expand the command-line inputs into a sorted, de-duplicated list of genome
/// files. A directory means every accepted file directly inside it; recursion
/// is one level deep only, so a folder of folders does not silently pull in an
/// entire project tree.
pub fn collect_inputs(inputs: &[std::path::PathBuf]) -> io::Result<Vec<std::path::PathBuf>> {
    let mut out: Vec<std::path::PathBuf> = Vec::new();
    for input in inputs {
        if input.is_dir() {
            for entry in std::fs::read_dir(input)? {
                let p = entry?.path();
                if p.is_file() && is_sequence_file(&p) {
                    out.push(p);
                }
            }
        } else if input.is_file() {
            out.push(input.clone());
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("input not found: {}", input.display()),
            ));
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
}

fn is_sequence_file(path: &Path) -> bool {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".gz") {
        // Accept <something>.fa.gz / .fasta.gz / .fna.gz only.
        let stem = &lower[..lower.len() - 3];
        return stem.ends_with(".fa") || stem.ends_with(".fasta") || stem.ends_with(".fna");
    }
    lower.ends_with(".fa") || lower.ends_with(".fasta") || lower.ends_with(".fna")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("mtato-fasta-test-{}-{}", std::process::id(), name));
        p
    }

    #[test]
    fn parses_multi_record_fasta_and_normalises() {
        let p = tmp("plain.fasta");
        {
            let mut f = File::create(&p).unwrap();
            writeln!(f, ">ctg1 some description here").unwrap();
            writeln!(f, "acgtacgt").unwrap();
            writeln!(f, "RYNN").unwrap();
            writeln!(f, ">ctg2").unwrap();
            writeln!(f, "AAAA").unwrap();
        }
        let g = read_genome(&p).unwrap();
        assert_eq!(g.contigs.len(), 2);
        assert_eq!(g.contigs[0].name, "ctg1");
        assert_eq!(g.contigs[0].seq, b"ACGTACGTNNNN");
        assert_eq!(g.contigs[1].seq, b"AAAA");
        assert_eq!(g.total_length(), 16);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn gzip_is_detected_from_magic_not_extension() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let p = tmp("weird_name.dat");
        {
            let f = File::create(&p).unwrap();
            let mut e = GzEncoder::new(f, Compression::fast());
            e.write_all(b">c1\nACGTACGT\n").unwrap();
            e.finish().unwrap();
        }
        let g = read_genome(&p).unwrap();
        assert_eq!(g.contigs[0].seq, b"ACGTACGT");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn sequence_without_a_header_is_still_read() {
        let p = tmp("noheader.fa");
        std::fs::write(&p, b"ACGTACGT\n").unwrap();
        let g = read_genome(&p).unwrap();
        assert_eq!(g.contigs.len(), 1);
        assert_eq!(g.contigs[0].seq, b"ACGTACGT");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn genome_name_strips_every_sequence_extension() {
        assert_eq!(genome_name(Path::new("/x/GCF_00001.fa")), "GCF_00001");
        assert_eq!(genome_name(Path::new("/x/GCF_00001.fasta.gz")), "GCF_00001");
        assert_eq!(genome_name(Path::new("/x/GCF_00001.fna.gz")), "GCF_00001");
        assert_eq!(
            genome_name(Path::new("SRR1__NCBI_SKESA.fasta")),
            "SRR1__NCBI_SKESA"
        );
    }

    #[test]
    fn directory_input_takes_only_sequence_files() {
        let dir = tmp("indir");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("a.fasta"), b">c\nACGT\n").unwrap();
        std::fs::write(dir.join("b.fna.gz"), b"not really gz").unwrap();
        std::fs::write(dir.join("notes.txt"), b"ignore me").unwrap();
        std::fs::write(dir.join("c.bam"), b"ignore me too").unwrap();
        let got = collect_inputs(std::slice::from_ref(&dir)).unwrap();
        assert_eq!(got.len(), 2);
        assert!(got.iter().all(|p| is_sequence_file(p)));
        std::fs::remove_dir_all(&dir).ok();
    }
}

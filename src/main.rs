use std::{cmp::max, io::BufRead, vec};

use aho_corasick::AhoCorasick;
use clap::{Arg, Command};
use jseqio::writer::SeqRecordWriter;

// Local alignment of needle against the haystack.
// Returns the match score of the best match. The maximum possible score
// is the length of the shorter sequence. 
fn smith_waterman(needle: &[u8], haystack: &[u8]) -> usize {
    let m = needle.len();
    let n = haystack.len();

    let mut score_matrix = vec![vec![0_isize; n + 1]; m + 1];
    let mut max_score = 0;

    // Fill the score matrix
    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch_score = if needle[i - 1] == haystack[j - 1] { 1 } else { 0 } as isize;

            score_matrix[i][j] = max(0, max(
                score_matrix[i - 1][j - 1] + match_mismatch_score,
                max(
                    score_matrix[i - 1][j] - 1,
                    score_matrix[i][j - 1] - 1,
                ),
            ));

            max_score = max(max_score, score_matrix[i][j]);
        }
    }

    max_score as usize
}

fn read_barcodes(filename: &str) -> Vec<Vec<u8>> {
    let file = std::io::BufReader::new(std::fs::File::open(filename).unwrap());
    file.lines().map(|x| x.unwrap().as_bytes().to_owned()).collect()
}

fn get_reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&x| match x {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => panic!("Invalid character in barcodes: ASCII value {}", x),
    }).collect()
}

// Score should be between 0 and 1
fn identify_barcodes_smith_waterman(barcodes: &[&[u8]], read: &[u8], score_threshold: f64) -> Vec<usize> {
    assert!(score_threshold >= 0.0 && score_threshold <= 1.0);
    let mut found_barcodes = Vec::<usize>::new();
    for (i, &barcode) in barcodes.iter().enumerate() {
        let score = smith_waterman(barcode, read);
        let fscore = score as f64 / barcode.len() as f64;
        if fscore >= score_threshold {
            found_barcodes.push(i);
        }
    }
    found_barcodes
}

// The automaton should have the barcodes, and then their reverse complements
fn identify_barcodes_aho_corasick(automaton: AhoCorasick, read: &[u8]) -> Vec<usize> {
    let mut found_barcodes = Vec::<usize>::new();
    for mat in automaton.find_iter(read) {
        let n_barcodes = automaton.patterns_len() / 2;
        let id = mat.pattern().as_usize() % n_barcodes; // Reverse complements are at the second half, hence modulo
        found_barcodes.push(id);
    }

    // We might have the reverse complement as well as forward, so we need to deduplicate
    found_barcodes.sort();
    found_barcodes.dedup();
    found_barcodes
}

fn main() {

    let cli = Command::new("barcode-binner")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg_required_else_help(true)
        .about("Reads fastq data from standard input, writes to fastq files by barcode")
        .arg(Arg::new("out-prefix")
            .long("out-prefix")
            .short('o')
            .help("Prefix for the output files") 
            .value_parser(clap::value_parser!(std::path::PathBuf))
        )
        .arg(Arg::new("barcodes")
            .long("barcodes")
            .short('b')
            .help("File with barcode sequences, one per line") 
            .value_parser(clap::value_parser!(std::path::PathBuf))
        );

    let matches = cli.get_matches();
    let barcode_filepath = matches.get_one::<std::path::PathBuf>("barcodes").unwrap();
    let out_prefix = matches.get_one::<std::path::PathBuf>("out-prefix").unwrap();

    let barcodes = read_barcodes(barcode_filepath.to_str().unwrap());
    let rc_barcodes = barcodes.iter().map(|x| get_reverse_complement(x)).collect::<Vec<_>>(); 
    let barcodes_and_rc = barcodes.iter().chain(rc_barcodes.iter()).collect::<Vec<_>>();

    let n_barcodes = barcodes.len();
    let aho_corasick = AhoCorasick::new(barcodes_and_rc).unwrap();

    let mut reader = jseqio::reader::DynamicFastXReader::from_stdin().unwrap();
    
    // Create writers for each barcode
    let mut writers = vec![];
    for barcode_id in 0..n_barcodes {
        let filename = format!("{}-barcode{}.fastq", out_prefix.to_str().unwrap(), barcode_id + 1); // 1-based indexing
        writers.push(jseqio::writer::DynamicFastXWriter::new_to_file(&filename).unwrap());
    }
    let mixed_filename = format!("{}-mixed.fastq", out_prefix.to_str().unwrap());
    writers.push(jseqio::writer::DynamicFastXWriter::new_to_file(&mixed_filename).unwrap());
    let none_filename = format!("{}-none.fastq", out_prefix.to_str().unwrap());
    writers.push(jseqio::writer::DynamicFastXWriter::new_to_file(&none_filename).unwrap());
    let mut written_counts = vec![0_usize; writers.len()];
    let mut hit_counts = vec![0_usize; n_barcodes];

    while let Some(rec) = reader.read_next().unwrap() {
        let mut barcode_id: Option<usize> = None;
        let mut have_multiple = false;
        for mat in aho_corasick.find_iter(rec.seq) {
            let id = mat.pattern().as_usize() % n_barcodes; // Reverse complements are at the second half, hence modulo
            hit_counts[id] += 1;
            if barcode_id.is_some_and(|x| x != id) {
                have_multiple = true;
            }
            barcode_id = Some(id);
        }
        let writer_idx = match barcode_id {
            Some(id) => {
                if have_multiple {
                    n_barcodes // Mixed writer
                } else {
                    id
                }
            },
            None => n_barcodes+1 // None writer
        };
        writers[writer_idx].write_ref_record(&rec).unwrap();
        written_counts[writer_idx] += 1;
    }
    
    for (idx, count) in written_counts.iter().enumerate() {
        #[allow(clippy::comparison_chain)]
        if idx < n_barcodes {
            eprintln!("Found {} reads with barcode {}", count, idx + 1); // 1-based indexing
        } else if idx == n_barcodes {
            eprintln!("Found {} reads with multiple barcodes", count);
        } else {
            eprintln!("Found {} reads with no barcodes", count);
        }
    }

    for (idx, count) in hit_counts.iter().enumerate() {
        eprintln!("{} total occurrences of barcode {}", count, idx+1); // 1-based indexing
    }

}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_smith_waterman(){
        let s1 = b"TAGATACGTACGTACGTGAAG";
        let s2 =      b"ACGTAAGTACGT"; // 1 substitution
        let s3 =      b"ACTACGTACXXGT"; // 1 dels, 2 inserts

        assert_eq!(smith_waterman(s2, s1), s2.len()-1);
        assert_eq!(smith_waterman(s3, s1), 12-1-2); // 12 matches, 1 del, 2 inserts
    }
}
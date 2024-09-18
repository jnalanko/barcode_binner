use std::{io::BufRead, vec};

use aho_corasick::AhoCorasick;
use clap::{Arg, Command};
use jseqio::writer::SeqRecordWriter;

fn read_barcodes(filename: &str) -> Vec<Vec<u8>> {
    let file = std::io::BufReader::new(std::fs::File::open(filename).unwrap());
    file.lines().map(|x| x.unwrap().as_bytes().to_owned()).collect()
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
    let n_barcodes = barcodes.len();
    let ac = AhoCorasick::new(barcodes).unwrap();

    let mut reader = jseqio::reader::DynamicFastXReader::from_stdin().unwrap();
    
    // Create writers for each barcode
    let mut writers = vec![];
    for barcode_id in 0..n_barcodes {
        let filename = format!("{}-barcode{}.fastq", out_prefix.to_str().unwrap(), barcode_id);
        writers.push(jseqio::writer::DynamicFastXWriter::new_to_file(&filename).unwrap());
    }
    let mixed_filename = format!("{}-mixed.fastq", out_prefix.to_str().unwrap());
    writers.push(jseqio::writer::DynamicFastXWriter::new_to_file(&mixed_filename).unwrap());
    let none_filename = format!("{}-none.fastq", out_prefix.to_str().unwrap());
    writers.push(jseqio::writer::DynamicFastXWriter::new_to_file(&none_filename).unwrap());
    let mut written_counts = vec![0; writers.len()];

    while let Some(rec) = reader.read_next().unwrap() {
        let mut barcode_id: Option<usize> = None;
        let mut have_multiple = false;
        for mat in ac.find_iter(rec.seq) {
            let id = mat.pattern().as_usize();
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
            eprintln!("Found {} reads with barcode {}", count, idx);
        } else if idx == n_barcodes {
            eprintln!("Found {} reads with multiple barcodes", count);
        } else {
            eprintln!("Found {} reads with no barcodes", count);
        }
    }

}

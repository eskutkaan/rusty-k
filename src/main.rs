use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::{Arc, Mutex};
use std::thread;
use xxhash_rust::xxh3::xxh3_64;

// Compute the canonical representation of a k-mer
fn canonical_kmer(kmer: &str) -> String {
    let rev_comp = revcomp(kmer);
    if kmer < rev_comp.as_str() {
        kmer.to_string()
    } else {
        rev_comp
    }
}

// Compute the reverse complement of a DNA sequence
fn revcomp(sequence: &str) -> String {
    let mut rev_comp = String::with_capacity(sequence.len());
    for c in sequence.chars().rev() {
        rev_comp.push(complement(c));
    }
    rev_comp
}

// Complement a DNA base
fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => base,
    }
}

// Hash a k-mer using xxHash
fn hash_kmer(kmer: &str) -> u64 {
    xxh3_64(kmer.as_bytes())
}

// Count canonical k-mers in a DNA sequence
fn count_canonical_kmers(sequence: &str, k: usize) -> HashMap<u64, u32> {
    let mut kmer_counts = HashMap::new();

    for i in 0..(sequence.len() - k + 1) {
        let kmer = &sequence[i..i+k];
        let canonical_kmer = canonical_kmer(kmer);
        let hash = hash_kmer(&canonical_kmer);
        *kmer_counts.entry(hash).or_insert(0) += 1;
    }

    kmer_counts
}

fn count_kmers_from_file_threaded(
    file_path: &str,
    k: usize,
    num_threads: usize,
) -> HashMap<u64, u32> {
    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut sequence = String::new();
    let mut is_fastq = false;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !sequence.is_empty() {
                sequences.push(sequence.clone());
                sequence.clear();
            }
        } else if line.starts_with('@') {
            if !sequence.is_empty() {
                sequences.push(sequence.clone());
                sequence.clear();
            }
            is_fastq = true;
        } else if is_fastq && (line.starts_with('+') || line.starts_with('#')) {
            // Skip quality lines in FASTQ
        } else {
            sequence.push_str(&line);
        }
    }

    if !sequence.is_empty() {
        sequences.push(sequence);
    }

    let kmer_counts = Arc::new(Mutex::new(HashMap::new()));
    let mut handles = Vec::with_capacity(num_threads);

    let chunk_size = (sequences.len() + num_threads - 1) / num_threads;
    let mut start = 0;

    for _ in 0..num_threads {
        let end = std::cmp::min(start + chunk_size, sequences.len());
        let sequences_chunk = sequences[start..end].to_vec();
        let kmer_counts = Arc::clone(&kmer_counts);

        let handle = thread::spawn(move || {
            let mut thread_counts = HashMap::new();

            for sequence in sequences_chunk {
                let counts = count_canonical_kmers(&sequence, k);
                for (kmer_hash, count) in counts {
                    *thread_counts.entry(kmer_hash).or_insert(0) += count;
                }
            }

            let mut global_counts = kmer_counts.lock().unwrap();
            for (kmer_hash, count) in thread_counts {
                *global_counts.entry(kmer_hash).or_insert(0) += count;
            }
        });

        handles.push(handle);
        start = end;
    }

    for handle in handles {
        handle.join().expect("Failed to join thread");
    }

    Arc::try_unwrap(kmer_counts)
        .expect("Failed to unwrap Arc")
        .into_inner()
        .expect("Failed to acquire mutex")
}

// Convert a k-mer hash back to the DNA sequence
fn hash_to_kmer(hash: u64, k: usize) -> String {
    let mut kmer = String::with_capacity(k);
    let mut remaining_hash = hash;

    for _ in 0..k {
        let base = match remaining_hash & 0x3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        };
        kmer.push(base);
        remaining_hash >>= 2;
    }

    kmer
}

fn main() {
    let mut args = env::args().skip(1); // Skip the program name
    let mut output_file = None;
    let mut k: Option<usize> = None;
    let mut input_file = None;
    let mut num_threads = 1;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-o" | "--output" => {
                output_file = args.next();
            }
            "-t" | "--threads" => {
                num_threads = args
                    .next()
                    .expect("Missing thread count")
                    .parse()
                    .expect("Invalid thread count");
            }
            _ if k.is_none() => {
                k = Some(arg.parse().expect("Invalid k-mer size"));
            }
            _ => {
                input_file = Some(arg);
            }
        }
    }

    let k = k.expect("K-mer size not provided");
    let input_file = input_file.expect("Input file not provided");

    let kmer_counts = count_kmers_from_file_threaded(&input_file, k, num_threads);

    if let Some(output_file) = output_file {
        let mut output = File::create(output_file).expect("Failed to create output file");
        for (kmer_hash, count) in &kmer_counts {
            let kmer = hash_to_kmer(*kmer_hash, k);
            writeln!(output, "{}	{}", kmer, count).expect("Failed to write to output file");
        }
    } else {
        for (kmer_hash, count) in &kmer_counts {
            let kmer = hash_to_kmer(*kmer_hash, k);
            println!("{} {}", kmer, count);
        }
    }
}


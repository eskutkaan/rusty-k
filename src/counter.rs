use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::thread;
use rusty_k::canonical_kmer;

pub fn count_canonical_kmers(sequence: &str, k: usize) -> HashMap<String, u64> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=(sequence.len() - k) {
        let kmer = &sequence[i..i + k];
        let canonical = canonical_kmer(kmer);
        *kmer_counts.entry(canonical).or_insert(0) += 1;
    }

    kmer_counts
}

pub fn count_kmers_from_file_threaded(
    file_path: &str,
    k: usize,
    num_threads: usize,
) -> anyhow::Result<HashMap<String, u64>> {
    let file = File::open(file_path)
        .map_err(|e| anyhow::anyhow!("Failed to open file '{}': {}", file_path, e))?;
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut sequence = String::new();

    for line in reader.lines() {
        let line = line.map_err(|e| anyhow::anyhow!("Failed to read line: {}", e))?;
        if line.starts_with('>') || line.starts_with('@') {
            if !sequence.is_empty() {
                sequences.push(sequence.clone());
                sequence.clear();
            }
        } else if line.starts_with('+') || line.starts_with('#') {
            // Skip quality lines in FASTQ
        } else {
            sequence.push_str(&line);
        }
    }

    if !sequence.is_empty() {
        sequences.push(sequence);
    }

    if sequences.is_empty() {
        return Ok(HashMap::new());
    }

    let kmer_counts = Arc::new(Mutex::new(HashMap::new()));
    let mut handles = Vec::with_capacity(num_threads);

    let chunk_size = sequences.len().div_ceil(num_threads);
    let mut start = 0;

    for _ in 0..num_threads {
        let end = std::cmp::min(start + chunk_size, sequences.len());
        if start >= sequences.len() {
            break;
        }
        let sequences_chunk = sequences[start..end].to_vec();
        let kmer_counts = Arc::clone(&kmer_counts);

        let handle = thread::spawn(move || {
            let mut thread_counts = HashMap::new();

            for sequence in sequences_chunk {
                let counts = count_canonical_kmers(&sequence, k);
                for (kmer, count) in counts {
                    *thread_counts.entry(kmer).or_insert(0) += count;
                }
            }

            let mut global_counts = kmer_counts.lock().unwrap();
            for (kmer, count) in thread_counts {
                *global_counts.entry(kmer).or_insert(0) += count;
            }
        });

        handles.push(handle);
        start = end;
    }

    for handle in handles {
        handle.join().map_err(|_| anyhow::anyhow!("Thread panicked"))?;
    }

    let counts = Arc::try_unwrap(kmer_counts)
        .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc"))?
        .into_inner()
        .map_err(|_| anyhow::anyhow!("Failed to acquire mutex"))?;

    Ok(counts)
}

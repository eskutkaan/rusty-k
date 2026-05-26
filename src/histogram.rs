use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::thread;

pub fn read_kmer_counts_threaded(
    file_path: &str,
    num_threads: usize,
) -> anyhow::Result<HashMap<String, u64>> {
    let file = File::open(file_path)
        .map_err(|e| anyhow::anyhow!("Failed to open file '{}': {}", file_path, e))?;
    let reader = BufReader::new(file);

    let mut lines = Vec::new();
    for line in reader.lines() {
        lines.push(line.map_err(|e| anyhow::anyhow!("Failed to read line: {}", e))?);
    }

    if lines.is_empty() {
        return Ok(HashMap::new());
    }

    let kmer_counts = Arc::new(Mutex::new(HashMap::new()));
    let mut handles = Vec::with_capacity(num_threads);

    let chunk_size = lines.len().div_ceil(num_threads);
    let mut start = 0;

    for _ in 0..num_threads {
        let end = std::cmp::min(start + chunk_size, lines.len());
        if start >= lines.len() {
            break;
        }
        let lines_chunk = lines[start..end].to_vec();
        let kmer_counts = Arc::clone(&kmer_counts);

        let handle = thread::spawn(move || {
            let mut thread_counts = HashMap::new();

            for line in lines_chunk {
                let mut parts = line.split_whitespace();
                if let (Some(kmer), Some(count_str)) = (parts.next(), parts.next()) {
                    if let Ok(count) = count_str.parse::<u64>() {
                        *thread_counts.entry(kmer.to_string()).or_insert(0) += count;
                    }
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

pub fn create_histogram(kmer_counts: &HashMap<String, u64>) -> HashMap<u64, u64> {
    let mut histogram = HashMap::new();

    for count in kmer_counts.values() {
        *histogram.entry(*count).or_insert(0) += 1;
    }

    histogram
}

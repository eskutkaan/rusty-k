use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::{Arc, Mutex};
use std::thread;

fn main() {
    let mut args = env::args().skip(1); // Skip the program name
    let mut output_file = None;
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
            _ => {
                input_file = Some(arg);
            }
        }
    }

    let input_file = input_file.expect("Input file not provided");
    let kmer_counts = read_kmer_counts_threaded(&input_file, num_threads);
    let mut histogram: Vec<_> = create_histogram(&kmer_counts).into_iter().collect();

    // Sort the histogram numerically by the first column (count)
    histogram.sort_by(|a, b| a.0.cmp(&b.0));

    if let Some(output_file) = output_file {
        let mut output = File::create(output_file).expect("Failed to create output file");
        for (count, frequency) in histogram {
            writeln!(output, "{}	{}", count, frequency).expect("Failed to write to output file");
        }
    } else {
        for (count, frequency) in histogram {
            println!("{} {}", count, frequency);
        }
    }
}

fn read_kmer_counts_threaded(file_path: &str, num_threads: usize) -> HashMap<String, u32> {
    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut lines = Vec::new();
    for line in reader.lines() {
        lines.push(line.expect("Failed to read line"));
    }

    let kmer_counts = Arc::new(Mutex::new(HashMap::new()));
    let mut handles = Vec::with_capacity(num_threads);

    let chunk_size = (lines.len() + num_threads - 1) / num_threads;
    let mut start = 0;

    for _ in 0..num_threads {
        let end = std::cmp::min(start + chunk_size, lines.len());
        let lines_chunk = lines[start..end].to_vec();
        let kmer_counts = Arc::clone(&kmer_counts);

        let handle = thread::spawn(move || {
            let mut thread_counts = HashMap::new();

            for line in lines_chunk {
                let mut parts = line.split_whitespace();
                let kmer = parts.next().expect("Invalid line format").to_string();
                let count: u32 = parts.next().expect("Invalid line format").parse().expect("Invalid count");
                *thread_counts.entry(kmer).or_insert(0) += count;
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
        handle.join().expect("Failed to join thread");
    }

    Arc::try_unwrap(kmer_counts)
        .expect("Failed to unwrap Arc")
        .into_inner()
        .expect("Failed to acquire mutex")
}

fn create_histogram(kmer_counts: &HashMap<String, u32>) -> HashMap<u32, u32> {
    let mut histogram = HashMap::new();

    for count in kmer_counts.values() {
        *histogram.entry(*count).or_insert(0) += 1;
    }

    histogram
}


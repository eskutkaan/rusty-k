use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::thread;

use rusty_k::canonical_kmer;
use xxhash_rust::xxh3::xxh3_64;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatRegion {
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub supporting_hits: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TandemRepeat {
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub period: usize,
    pub copies: usize,
    pub motif: String,
    pub consensus: String,
    pub mismatches: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TandemRepeatSummary {
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub period: usize,
    pub copies: usize,
    pub mismatches: usize,
    pub supporting_calls: usize,
    pub motif: String,
    pub consensus: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RepeatBenchmarkRow {
    pub kmer_size: usize,
    pub region_count: usize,
    pub total_repeat_bp: usize,
    pub repeat_fraction: f64,
    pub repeat_bp_per_mb: f64,
    pub region_density_per_mb: f64,
    pub average_region_len: f64,
    pub supporting_hits: u64,
}

pub fn recommend_benchmark_row<'a>(
    rows: &'a [RepeatBenchmarkRow],
    target_repeat_fraction: f64,
) -> Option<&'a RepeatBenchmarkRow> {
    if rows.is_empty() {
        return None;
    }

    let target_repeat_fraction = if target_repeat_fraction.is_finite() && target_repeat_fraction > 0.0 {
        target_repeat_fraction
    } else {
        0.05
    };

    if let Some(row) = rows
        .iter()
        .find(|row| row.repeat_fraction > 0.0 && row.repeat_fraction <= target_repeat_fraction)
    {
        return Some(row);
    }

    rows.iter().min_by(|left, right| {
        left.repeat_fraction
            .partial_cmp(&right.repeat_fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| right.region_count.cmp(&left.region_count))
            .then_with(|| left.kmer_size.cmp(&right.kmer_size))
    })
}

pub fn summarize_tandem_repeats(mut repeats: Vec<TandemRepeat>) -> Vec<TandemRepeatSummary> {
    if repeats.is_empty() {
        return Vec::new();
    }

    repeats.sort_by(|left, right| {
        left.contig
            .cmp(&right.contig)
            .then_with(|| left.start.cmp(&right.start))
            .then_with(|| left.end.cmp(&right.end))
            .then_with(|| left.period.cmp(&right.period))
    });

    let mut summaries = Vec::new();
    let mut current_cluster: Vec<TandemRepeat> = Vec::new();

    let flush_cluster = |cluster: &mut Vec<TandemRepeat>, summaries: &mut Vec<TandemRepeatSummary>| {
        if cluster.is_empty() {
            return;
        }

        cluster.sort_by(|left, right| {
            left.mismatches
                .cmp(&right.mismatches)
                .then_with(|| (right.end - right.start).cmp(&(left.end - left.start)))
                .then_with(|| right.copies.cmp(&left.copies))
                .then_with(|| left.period.cmp(&right.period))
        });

        let representative = cluster[0].clone();
        let start = cluster.iter().map(|repeat| repeat.start).min().unwrap();
        let end = cluster.iter().map(|repeat| repeat.end).max().unwrap();
        let supporting_calls = cluster.len();

        summaries.push(TandemRepeatSummary {
            contig: representative.contig,
            start,
            end,
            period: representative.period,
            copies: representative.copies,
            mismatches: representative.mismatches,
            supporting_calls,
            motif: representative.motif,
            consensus: representative.consensus,
        });

        cluster.clear();
    };

    for repeat in repeats {
        let should_start_new_cluster = current_cluster
            .last()
            .map(|previous| previous.contig != repeat.contig || repeat.start > previous.end)
            .unwrap_or(true);

        if should_start_new_cluster {
            flush_cluster(&mut current_cluster, &mut summaries);
        }

        current_cluster.push(repeat);
    }

    flush_cluster(&mut current_cluster, &mut summaries);
    summaries
}

fn canonical_kmer_hash(kmer: &str) -> u64 {
    let canonical = canonical_kmer(kmer);
    xxh3_64(canonical.as_bytes())
}

fn is_valid_dna_kmer(kmer: &str) -> bool {
    kmer.chars()
        .all(|base| matches!(base, 'A' | 'C' | 'G' | 'T'))
}

fn is_valid_dna_sequence(sequence: &str) -> bool {
    sequence
        .chars()
        .all(|base| matches!(base, 'A' | 'C' | 'G' | 'T'))
}

fn hamming_distance(left: &str, right: &str) -> usize {
    left.chars()
        .zip(right.chars())
        .filter(|(left_base, right_base)| left_base != right_base)
        .count()
}

fn consensus_sequence(copies: &[&str], period: usize) -> String {
    let mut consensus = String::with_capacity(period);

    for column in 0..period {
        let mut counts = [0usize; 4];

        for copy in copies {
            if let Some(base) = copy.as_bytes().get(column) {
                match base {
                    b'A' => counts[0] += 1,
                    b'C' => counts[1] += 1,
                    b'G' => counts[2] += 1,
                    b'T' => counts[3] += 1,
                    _ => {}
                }
            }
        }

        let consensus_base = match counts
            .iter()
            .enumerate()
            .max_by_key(|(_, count)| **count)
            .map(|(index, _)| index)
            .unwrap_or(0)
        {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };

        consensus.push(consensus_base);
    }

    consensus
}

pub fn read_fasta_records(file_path: &str) -> anyhow::Result<Vec<FastaRecord>> {
    let file = File::open(file_path)
        .map_err(|e| anyhow::anyhow!("Failed to open file '{}': {}", file_path, e))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut sequence = String::new();

    for line in reader.lines() {
        let line = line.map_err(|e| anyhow::anyhow!("Failed to read line: {}", e))?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if !sequence.is_empty() {
                    records.push(FastaRecord {
                        id,
                        sequence: std::mem::take(&mut sequence),
                    });
                }
            }

            current_id = Some(
                header
                    .split_whitespace()
                    .next()
                    .unwrap_or(header)
                    .to_string(),
            );
        } else if line.starts_with(';') {
            continue;
        } else if current_id.is_some() {
            sequence.push_str(&line.to_ascii_uppercase());
        } else {
            return Err(anyhow::anyhow!(
                "Encountered sequence data before a FASTA header in '{}'",
                file_path
            ));
        }
    }

    if let Some(id) = current_id.take() {
        if !sequence.is_empty() {
            records.push(FastaRecord { id, sequence });
        }
    }

    Ok(records)
}

fn chunk_ranges(len: usize, num_threads: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return Vec::new();
    }

    let num_threads = num_threads.max(1);
    let chunk_size = len.div_ceil(num_threads);
    let mut ranges = Vec::new();
    let mut start = 0;

    while start < len {
        let end = std::cmp::min(start + chunk_size, len);
        ranges.push((start, end));
        start = end;
    }

    ranges
}

fn count_canonical_kmer_hashes_in_range(
    sequence: &str,
    k: usize,
    start: usize,
    core_end: usize,
) -> HashMap<u64, u64> {
    let mut kmer_counts = HashMap::new();

    if start >= core_end || sequence.len() < k {
        return kmer_counts;
    }

    let scan_end = std::cmp::min(core_end + k - 1, sequence.len());

    for pos in start..core_end {
        if pos + k > scan_end {
            break;
        }

        let kmer = &sequence[pos..pos + k];
        if !is_valid_dna_kmer(kmer) {
            continue;
        }

        *kmer_counts.entry(canonical_kmer_hash(kmer)).or_insert(0) += 1;
    }

    kmer_counts
}

fn collect_repeat_hits_in_range(
    sequence: &str,
    k: usize,
    start: usize,
    core_end: usize,
    repeated_kmers: &HashSet<u64>,
) -> Vec<(usize, usize)> {
    let mut hits = Vec::new();

    if start >= core_end || sequence.len() < k {
        return hits;
    }

    let scan_end = std::cmp::min(core_end + k - 1, sequence.len());

    for pos in start..core_end {
        if pos + k > scan_end {
            break;
        }

        let kmer = &sequence[pos..pos + k];
        if !is_valid_dna_kmer(kmer) {
            continue;
        }

        if repeated_kmers.contains(&canonical_kmer_hash(kmer)) {
            hits.push((pos, pos + k));
        }
    }

    hits
}

fn merge_repeat_hits(contig: &str, mut hits: Vec<(usize, usize)>) -> Vec<RepeatRegion> {
    if hits.is_empty() {
        return Vec::new();
    }

    hits.sort_by_key(|hit| hit.0);

    let mut regions = Vec::new();
    let mut current_start = hits[0].0;
    let mut current_end = hits[0].1;
    let mut supporting_hits = 1u64;

    for (start, end) in hits.into_iter().skip(1) {
        if start <= current_end {
            current_end = current_end.max(end);
            supporting_hits += 1;
        } else {
            regions.push(RepeatRegion {
                contig: contig.to_string(),
                start: current_start,
                end: current_end,
                supporting_hits,
            });

            current_start = start;
            current_end = end;
            supporting_hits = 1;
        }
    }

    regions.push(RepeatRegion {
        contig: contig.to_string(),
        start: current_start,
        end: current_end,
        supporting_hits,
    });

    regions
}

fn count_approximate_tandem_repeat_copies(
    sequence: &str,
    start: usize,
    period: usize,
    max_mismatches_per_copy: usize,
) -> usize {
    let motif = &sequence[start..start + period];
    let mut copies = 1usize;

    while start + (copies + 1) * period <= sequence.len() {
        let next_start = start + copies * period;
        let next_end = next_start + period;
        let next_copy = &sequence[next_start..next_end];

        if hamming_distance(motif, next_copy) <= max_mismatches_per_copy {
            copies += 1;
        } else {
            break;
        }
    }

    copies
}

fn approximate_repeat_mismatches(sequence: &str, start: usize, period: usize, copies: usize) -> usize {
    let mut total_mismatches = 0usize;

    if copies == 0 {
        return total_mismatches;
    }

    let copy_slices: Vec<&str> = (0..copies)
        .map(|copy_index| {
            let copy_start = start + copy_index * period;
            &sequence[copy_start..copy_start + period]
        })
        .collect();

    let consensus = consensus_sequence(&copy_slices, period);

    for copy in copy_slices {
        total_mismatches += hamming_distance(copy, &consensus);
    }

    total_mismatches
}

fn find_tandem_repeats_in_range(
    sequence: &str,
    contig: &str,
    start: usize,
    core_end: usize,
    min_period: usize,
    max_period: usize,
    min_copies: usize,
    min_length: usize,
    max_mismatches_per_copy: usize,
) -> Vec<TandemRepeat> {
    let mut repeats = Vec::new();

    if start >= core_end || sequence.len() < min_period.saturating_mul(min_copies) {
        return repeats;
    }

    for pos in start..core_end {
        for period in min_period..=max_period {
            if pos + period * min_copies > sequence.len() {
                break;
            }

            if pos >= period && &sequence[pos - period..pos] == &sequence[pos..pos + period] {
                continue;
            }

            let motif = &sequence[pos..pos + period];
            if !is_valid_dna_sequence(motif) {
                continue;
            }

            let copies = count_approximate_tandem_repeat_copies(
                sequence,
                pos,
                period,
                max_mismatches_per_copy,
            );
            let repeat_length = period * copies;
            if copies < min_copies || repeat_length < min_length {
                continue;
            }

            let mismatches = approximate_repeat_mismatches(sequence, pos, period, copies);
            let copy_slices: Vec<&str> = (0..copies)
                .map(|copy_index| {
                    let copy_start = pos + copy_index * period;
                    &sequence[copy_start..copy_start + period]
                })
                .collect();
            let consensus = consensus_sequence(&copy_slices, period);

            repeats.push(TandemRepeat {
                contig: contig.to_string(),
                start: pos,
                end: pos + repeat_length,
                period,
                copies,
                motif: motif.to_string(),
                consensus,
                mismatches,
            });
        }
    }

    repeats
}

fn dedupe_tandem_repeats(mut repeats: Vec<TandemRepeat>) -> Vec<TandemRepeat> {
    if repeats.is_empty() {
        return repeats;
    }

    repeats.sort_by(|left, right| {
        left.contig
            .cmp(&right.contig)
            .then_with(|| left.start.cmp(&right.start))
            .then_with(|| left.end.cmp(&right.end))
            .then_with(|| left.period.cmp(&right.period))
            .then_with(|| right.copies.cmp(&left.copies))
    });

    let mut deduped: Vec<TandemRepeat> = Vec::new();

    for repeat in repeats {
        match deduped.last_mut() {
            Some(previous)
                if previous.contig == repeat.contig && repeat.start < previous.end =>
            {
                let previous_span = previous.end.saturating_sub(previous.start);
                let repeat_span = repeat.end.saturating_sub(repeat.start);

                if repeat_span > previous_span
                    || (repeat_span == previous_span && repeat.copies > previous.copies)
                    || (repeat_span == previous_span
                        && repeat.copies == previous.copies
                        && repeat.period < previous.period)
                {
                    *previous = repeat;
                }
            }
            _ => deduped.push(repeat),
        }
    }

    deduped
}

fn find_tandem_repeats_in_sequence_threaded(
    sequence: &str,
    contig: &str,
    min_period: usize,
    max_period: usize,
    min_copies: usize,
    min_length: usize,
    max_mismatches_per_copy: usize,
    num_threads: usize,
) -> Vec<TandemRepeat> {
    let ranges = chunk_ranges(sequence.len(), num_threads);
    if ranges.is_empty() {
        return Vec::new();
    }

    let sequence = Arc::new(sequence.to_string());
    let contig = Arc::new(contig.to_string());
    let repeats = Arc::new(Mutex::new(Vec::new()));
    let mut handles = Vec::with_capacity(ranges.len());

    for (start, core_end) in ranges {
        let sequence = Arc::clone(&sequence);
        let contig = Arc::clone(&contig);
        let repeats = Arc::clone(&repeats);

        let handle = thread::spawn(move || {
            let local_repeats = find_tandem_repeats_in_range(
                &sequence,
                &contig,
                start,
                core_end,
                min_period,
                max_period,
                min_copies,
                min_length,
                max_mismatches_per_copy,
            );

            let mut global_repeats = repeats.lock().unwrap();
            global_repeats.extend(local_repeats);
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    let repeats = Arc::try_unwrap(repeats)
        .unwrap()
        .into_inner()
        .unwrap();
    dedupe_tandem_repeats(repeats)
}

pub fn find_tandem_repeats_from_records_threaded(
    records: &[FastaRecord],
    min_period: usize,
    max_period: usize,
    min_copies: usize,
    min_length: usize,
    max_mismatches_per_copy: usize,
    num_threads: usize,
) -> Vec<TandemRepeat> {
    let mut tandem_repeats = Vec::new();

    for record in records {
        if record.sequence.len() < min_period.saturating_mul(min_copies) {
            continue;
        }

        tandem_repeats.extend(find_tandem_repeats_in_sequence_threaded(
            &record.sequence,
            &record.id,
            min_period,
            max_period,
            min_copies,
            min_length,
            max_mismatches_per_copy,
            num_threads,
        ));
    }

    tandem_repeats
}

pub fn find_tandem_repeats_from_fasta_threaded(
    file_path: &str,
    min_period: usize,
    max_period: usize,
    min_copies: usize,
    min_length: usize,
    max_mismatches_per_copy: usize,
    num_threads: usize,
) -> anyhow::Result<Vec<TandemRepeat>> {
    let records = read_fasta_records(file_path)?;
    Ok(find_tandem_repeats_from_records_threaded(
        &records,
        min_period,
        max_period,
        min_copies,
        min_length,
        max_mismatches_per_copy,
        num_threads,
    ))
}

fn count_repeated_kmers_in_sequence(
    sequence: &str,
    k: usize,
    num_threads: usize,
) -> HashMap<u64, u64> {
    let ranges = chunk_ranges(sequence.len(), num_threads);
    if ranges.is_empty() {
        return HashMap::new();
    }

    let counts = Arc::new(Mutex::new(HashMap::new()));
    let mut handles = Vec::with_capacity(ranges.len());

    for (start, core_end) in ranges {
        let sequence = sequence.to_string();
        let counts = Arc::clone(&counts);

        let handle = thread::spawn(move || {
            let local_counts = count_canonical_kmer_hashes_in_range(&sequence, k, start, core_end);
            let mut global_counts = counts.lock().unwrap();

            for (kmer_hash, count) in local_counts {
                *global_counts.entry(kmer_hash).or_insert(0) += count;
            }
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    Arc::try_unwrap(counts)
        .unwrap()
        .into_inner()
        .unwrap()
}

fn collect_repeat_hits_in_sequence(
    sequence: &str,
    k: usize,
    num_threads: usize,
    repeated_kmers: &HashSet<u64>,
) -> Vec<(usize, usize)> {
    let ranges = chunk_ranges(sequence.len(), num_threads);
    if ranges.is_empty() {
        return Vec::new();
    }

    let repeated_kmers = Arc::new(repeated_kmers.clone());
    let hits = Arc::new(Mutex::new(Vec::new()));
    let mut handles = Vec::with_capacity(ranges.len());

    for (start, core_end) in ranges {
        let sequence = sequence.to_string();
        let repeated_kmers = Arc::clone(&repeated_kmers);
        let hits = Arc::clone(&hits);

        let handle = thread::spawn(move || {
            let local_hits = collect_repeat_hits_in_range(&sequence, k, start, core_end, &repeated_kmers);
            let mut global_hits = hits.lock().unwrap();
            global_hits.extend(local_hits);
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    Arc::try_unwrap(hits).unwrap().into_inner().unwrap()
}

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

pub fn find_repeat_regions_from_records_threaded(
    records: &[FastaRecord],
    k: usize,
    num_threads: usize,
    min_count: u64,
) -> Vec<RepeatRegion> {
    let mut repeat_regions = Vec::new();

    for record in records {
        if record.sequence.len() < k {
            continue;
        }

        let counts = count_repeated_kmers_in_sequence(&record.sequence, k, num_threads);
        let repeated_kmers: HashSet<u64> = counts
            .into_iter()
            .filter_map(|(kmer_hash, count)| (count >= min_count).then_some(kmer_hash))
            .collect();

        if repeated_kmers.is_empty() {
            continue;
        }

        let hits = collect_repeat_hits_in_sequence(
            &record.sequence,
            k,
            num_threads,
            &repeated_kmers,
        );

        repeat_regions.extend(merge_repeat_hits(&record.id, hits));
    }

    repeat_regions
}

pub fn find_repeat_regions_from_fasta_threaded(
    file_path: &str,
    k: usize,
    num_threads: usize,
    min_count: u64,
) -> anyhow::Result<Vec<RepeatRegion>> {
    let records = read_fasta_records(file_path)?;
    Ok(find_repeat_regions_from_records_threaded(
        &records,
        k,
        num_threads,
        min_count,
    ))
}

fn summarize_repeat_regions(
    regions: &[RepeatRegion],
    assembly_length: usize,
    kmer_size: usize,
) -> RepeatBenchmarkRow {
    let region_count = regions.len();
    let total_repeat_bp = regions
        .iter()
        .map(|region| region.end.saturating_sub(region.start))
        .sum();
    let supporting_hits = regions.iter().map(|region| region.supporting_hits).sum();
    let average_region_len = if region_count == 0 {
        0.0
    } else {
        total_repeat_bp as f64 / region_count as f64
    };
    let assembly_mb = if assembly_length == 0 {
        0.0
    } else {
        assembly_length as f64 / 1_000_000.0
    };
    let repeat_fraction = if assembly_length == 0 {
        0.0
    } else {
        total_repeat_bp as f64 / assembly_length as f64
    };
    let repeat_bp_per_mb = if assembly_mb == 0.0 {
        0.0
    } else {
        total_repeat_bp as f64 / assembly_mb
    };
    let region_density_per_mb = if assembly_mb == 0.0 {
        0.0
    } else {
        region_count as f64 / assembly_mb
    };

    RepeatBenchmarkRow {
        kmer_size,
        region_count,
        total_repeat_bp,
        repeat_fraction,
        repeat_bp_per_mb,
        region_density_per_mb,
        average_region_len,
        supporting_hits,
    }
}

pub fn benchmark_repeat_k_sizes_from_records(
    records: &[FastaRecord],
    min_k: usize,
    max_k: usize,
    step: usize,
    num_threads: usize,
    min_count: u64,
) -> anyhow::Result<Vec<RepeatBenchmarkRow>> {
    if step == 0 {
        anyhow::bail!("step must be greater than 0");
    }
    if min_k > max_k {
        anyhow::bail!("min_k must be less than or equal to max_k");
    }

    let assembly_length: usize = records.iter().map(|record| record.sequence.len()).sum();
    let mut rows = Vec::new();

    let mut k = min_k;
    while k <= max_k {
        let regions = find_repeat_regions_from_records_threaded(records, k, num_threads, min_count);
        rows.push(summarize_repeat_regions(&regions, assembly_length, k));

        match k.checked_add(step) {
            Some(next_k) => k = next_k,
            None => break,
        }
    }

    Ok(rows)
}

pub fn benchmark_repeat_k_sizes_from_fasta(
    file_path: &str,
    min_k: usize,
    max_k: usize,
    step: usize,
    num_threads: usize,
    min_count: u64,
) -> anyhow::Result<Vec<RepeatBenchmarkRow>> {
    let records = read_fasta_records(file_path)?;
    benchmark_repeat_k_sizes_from_records(&records, min_k, max_k, step, num_threads, min_count)
}

#[cfg(test)]
mod tests {
    use super::{
        benchmark_repeat_k_sizes_from_records, find_repeat_regions_from_records_threaded,
        find_tandem_repeats_from_records_threaded, merge_repeat_hits,
        recommend_benchmark_row, summarize_tandem_repeats, FastaRecord, RepeatBenchmarkRow,
        TandemRepeat,
    };

    #[test]
    fn merges_overlapping_repeat_hits() {
        let regions = merge_repeat_hits("chr1", vec![(0, 4), (2, 6), (10, 14)]);

        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].contig, "chr1");
        assert_eq!(regions[0].start, 0);
        assert_eq!(regions[0].end, 6);
        assert_eq!(regions[0].supporting_hits, 2);
        assert_eq!(regions[1].start, 10);
        assert_eq!(regions[1].end, 14);
        assert_eq!(regions[1].supporting_hits, 1);
    }

    #[test]
    fn finds_repeat_regions_in_simple_sequence() {
        let records = vec![FastaRecord {
            id: "chr1".to_string(),
            sequence: "AAAAAA".to_string(),
        }];

        let regions = find_repeat_regions_from_records_threaded(&records, 3, 2, 2);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].contig, "chr1");
        assert_eq!(regions[0].start, 0);
        assert_eq!(regions[0].end, 6);
        assert_eq!(regions[0].supporting_hits, 4);
    }

    #[test]
    fn benchmarks_multiple_k_sizes() {
        let records = vec![FastaRecord {
            id: "chr1".to_string(),
            sequence: "AAAAAA".to_string(),
        }];

        let rows = benchmark_repeat_k_sizes_from_records(&records, 3, 5, 1, 2, 2).unwrap();

        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].kmer_size, 3);
        assert_eq!(rows[0].region_count, 1);
        assert!(rows[0].repeat_fraction > 0.0);
        assert!(rows[0].repeat_bp_per_mb > 0.0);
        assert!(rows[0].average_region_len > 0.0);
    }

    #[test]
    fn finds_tandem_repeats_in_simple_sequence() {
        let records = vec![FastaRecord {
            id: "chr1".to_string(),
            sequence: "ATATATGG".to_string(),
        }];

        let repeats = find_tandem_repeats_from_records_threaded(&records, 2, 4, 2, 6, 1, 2);

        assert_eq!(repeats.len(), 1);
        assert_eq!(repeats[0].contig, "chr1");
        assert_eq!(repeats[0].start, 0);
        assert_eq!(repeats[0].end, 6);
        assert_eq!(repeats[0].period, 2);
        assert_eq!(repeats[0].copies, 3);
        assert_eq!(repeats[0].motif, "AT");
        assert_eq!(repeats[0].consensus, "AT");
        assert_eq!(repeats[0].mismatches, 0);
    }

    #[test]
    fn tandem_repeat_dedupes_shifted_phases() {
        let repeats = find_tandem_repeats_from_records_threaded(
            &[FastaRecord {
                id: "chr1".to_string(),
                sequence: "ATATATGG".to_string(),
            }],
            2,
            4,
            2,
            6,
            1,
            2,
        );

        assert_eq!(repeats.len(), 1);
        assert_eq!(repeats[0].start, 0);
        assert_eq!(repeats[0].end, 6);
        assert_eq!(repeats[0].period, 2);
        assert_eq!(repeats[0].copies, 3);
        assert_eq!(repeats[0].consensus, "AT");
    }

    #[test]
    fn tandem_repeat_min_length_filters_short_calls() {
        let repeats = find_tandem_repeats_from_records_threaded(
            &[FastaRecord {
                id: "chr1".to_string(),
                sequence: "ATATGG".to_string(),
            }],
            2,
            4,
            2,
            6,
            1,
            2,
        );

        assert!(repeats.is_empty());
    }

    #[test]
    fn approximate_tandem_repeats_allow_one_mismatch_per_copy() {
        let repeats = find_tandem_repeats_from_records_threaded(
            &[FastaRecord {
                id: "chr1".to_string(),
                sequence: "ATATAG".to_string(),
            }],
            2,
            4,
            3,
            6,
            1,
            2,
        );

        assert_eq!(repeats[0].mismatches, 1);
        assert_eq!(repeats[0].period, 2);
        assert_eq!(repeats[0].copies, 3);
        assert_eq!(repeats[0].consensus, "AT");
        assert!(repeats[0].mismatches <= 2);
    }

    #[test]
    fn benchmark_recommendation_prefers_first_informative_k() {
        let rows = vec![
            RepeatBenchmarkRow {
                kmer_size: 15,
                region_count: 100,
                total_repeat_bp: 1000,
                repeat_fraction: 0.29,
                repeat_bp_per_mb: 290000.0,
                region_density_per_mb: 10000.0,
                average_region_len: 10.0,
                supporting_hits: 1000,
            },
            RepeatBenchmarkRow {
                kmer_size: 17,
                region_count: 60,
                total_repeat_bp: 400,
                repeat_fraction: 0.08,
                repeat_bp_per_mb: 80000.0,
                region_density_per_mb: 6000.0,
                average_region_len: 6.0,
                supporting_hits: 500,
            },
            RepeatBenchmarkRow {
                kmer_size: 19,
                region_count: 30,
                total_repeat_bp: 200,
                repeat_fraction: 0.04,
                repeat_bp_per_mb: 40000.0,
                region_density_per_mb: 3000.0,
                average_region_len: 7.0,
                supporting_hits: 200,
            },
            RepeatBenchmarkRow {
                kmer_size: 21,
                region_count: 20,
                total_repeat_bp: 150,
                repeat_fraction: 0.03,
                repeat_bp_per_mb: 30000.0,
                region_density_per_mb: 2000.0,
                average_region_len: 7.5,
                supporting_hits: 150,
            },
        ];

        let recommended = recommend_benchmark_row(&rows, 0.05).unwrap();
        assert_eq!(recommended.kmer_size, 19);
    }

    #[test]
    fn summarize_tandem_repeats_groups_overlapping_calls() {
        let summaries = summarize_tandem_repeats(vec![
            TandemRepeat {
                contig: "chr1".to_string(),
                start: 0,
                end: 6,
                period: 2,
                copies: 3,
                motif: "AT".to_string(),
                consensus: "AT".to_string(),
                mismatches: 1,
            },
            TandemRepeat {
                contig: "chr1".to_string(),
                start: 2,
                end: 8,
                period: 2,
                copies: 3,
                motif: "TA".to_string(),
                consensus: "TA".to_string(),
                mismatches: 2,
            },
            TandemRepeat {
                contig: "chr1".to_string(),
                start: 20,
                end: 30,
                period: 5,
                copies: 2,
                motif: "AACGT".to_string(),
                consensus: "AACGT".to_string(),
                mismatches: 0,
            },
        ]);

        assert_eq!(summaries.len(), 2);
        assert_eq!(summaries[0].contig, "chr1");
        assert_eq!(summaries[0].start, 0);
        assert_eq!(summaries[0].end, 8);
        assert_eq!(summaries[0].supporting_calls, 2);
        assert_eq!(summaries[1].start, 20);
        assert_eq!(summaries[1].supporting_calls, 1);
    }
}

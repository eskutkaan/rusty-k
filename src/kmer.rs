//! Canonical k-mer encoding (k ≤ 32) and parallel counting.

use crate::error::{Error, Result};
use crossbeam_channel::bounded;
use fxhash::FxHashMap;
use needletail::parse_fastx_file;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

/// 2-bit encoding: A=00, C=01, G=10, T=11.  Ambiguous bases are skipped.
#[inline]
fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Reverse-complement a 2-bit encoded k-mer of length k.
#[inline]
pub fn revcomp(kmer: u64, k: u8) -> u64 {
    let mut x = !kmer;
    // reverse the 2-bit chunks
    x = ((x >> 2) & 0x3333_3333_3333_3333) | ((x & 0x3333_3333_3333_3333) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    x = ((x >> 8) & 0x00FF_00FF_00FF_00FF) | ((x & 0x00FF_00FF_00FF_00FF) << 8);
    x = ((x >> 16) & 0x0000_FFFF_0000_FFFF) | ((x & 0x0000_FFFF_0000_FFFF) << 16);
    x = (x >> 32) | (x << 32);
    x >> (64 - 2 * k as u64)
}

/// Decode a 2-bit k-mer back to an ASCII string (for output).
pub fn decode_kmer(kmer: u64, k: u8) -> String {
    const LUT: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut s = vec![0u8; k as usize];
    for i in 0..k {
        let shift = 2 * (k - 1 - i) as u64;
        s[i as usize] = LUT[((kmer >> shift) & 0b11) as usize];
    }
    unsafe { String::from_utf8_unchecked(s) }
}

/// Extract all valid k-mers from a sequence, optionally canonicalised.
#[allow(dead_code)]
pub fn extract_kmers(seq: &[u8], k: u8, canonical: bool) -> Vec<u64> {
    if seq.len() < k as usize {
        return Vec::new();
    }
    let mut kmers = Vec::with_capacity(seq.len() - k as usize + 1);
    let mut current: u64 = 0;
    let mut valid = 0u8;
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    for &b in seq {
        match encode_base(b) {
            Some(bits) => {
                current = ((current << 2) | bits) & mask;
                valid = valid.saturating_add(1);
                if valid >= k {
                    let kmer = if canonical {
                        let rc = revcomp(current, k);
                        current.min(rc)
                    } else {
                        current
                    };
                    kmers.push(kmer);
                }
            }
            None => {
                valid = 0;
                current = 0;
            }
        }
    }
    kmers
}

/// Count k-mers across an entire FASTA/FASTQ file (parallel over records).
pub fn count_kmers(path: &Path, k: u8, canonical: bool) -> Result<FxHashMap<u64, u64>> {
    count_kmers_with_threads(path, k, canonical, 1)
}

/// Count k-mers with bounded parallel record processing.
pub fn count_kmers_with_threads(
    path: &Path,
    k: u8,
    canonical: bool,
    threads: usize,
) -> Result<FxHashMap<u64, u64>> {
    if k == 0 || k > 32 {
        return Err(Error::InvalidK(k));
    }

    let workers = threads.max(1);
    if workers == 1 {
        return count_kmers_sequential(path, k, canonical);
    }

    let (sender, receiver) = bounded::<CountTask>(workers * 2);
    let (result_sender, result_receiver) = bounded::<FxHashMap<u64, u64>>(workers);
    let mut parse_error = None;

    std::thread::scope(|scope| {
        for _ in 0..workers {
            let receiver = receiver.clone();
            let result_sender = result_sender.clone();
            scope.spawn(move || {
                while let Ok(task) = receiver.recv() {
                    let counts = count_sequence_chunk(&task, k, canonical);
                    if result_sender.send(counts).is_err() {
                        break;
                    }
                }
            });
        }
        drop(receiver);
        drop(result_sender);

        let mut reader = match parse_fastx_file(path).map_err(|e| Error::Other(e.to_string())) {
            Ok(reader) => reader,
            Err(error) => {
                parse_error = Some(error);
                drop(sender);
                return;
            }
        };
        while let Some(rec) = reader.next() {
            match rec {
                Ok(rec) => {
                    let seq = rec.seq();
                    let chunk_size = 4 * 1024 * 1024;
                    let mut chunk_start = 0usize;
                    while chunk_start < seq.len() {
                        let chunk_end = (chunk_start + chunk_size).min(seq.len());
                        let left = chunk_start.saturating_sub(k as usize - 1);
                        let right = (chunk_end + k as usize - 1).min(seq.len());
                        let task = CountTask {
                            sequence: seq[left..right].to_vec(),
                            global_offset: left,
                            owned_start: chunk_start,
                            owned_end: chunk_end,
                        };
                        if sender.send(task).is_err() {
                            break;
                        }
                        chunk_start = chunk_end;
                    }
                }
                Err(error) => {
                    parse_error = Some(Error::Other(error.to_string()));
                    break;
                }
            }
        }
        drop(sender);
    });

    if let Some(error) = parse_error {
        return Err(error);
    }
    let mut global = FxHashMap::default();
    for local in result_receiver {
        merge_counts(&mut global, local);
    }
    Ok(global)
}

fn count_kmers_sequential(path: &Path, k: u8, canonical: bool) -> Result<FxHashMap<u64, u64>> {
    let mut reader = parse_fastx_file(path).map_err(|e| Error::Other(e.to_string()))?;
    let mut global = FxHashMap::default();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| Error::Other(e.to_string()))?;
        merge_counts(
            &mut global,
            count_sequence(rec.seq().as_ref(), k, canonical),
        );
    }
    Ok(global)
}

fn count_sequence(seq: &[u8], k: u8, canonical: bool) -> FxHashMap<u64, u64> {
    let mut counts = FxHashMap::default();
    for_each_kmer_position(seq, k, canonical, |_, value| {
        *counts.entry(value).or_insert(0) += 1;
    });
    counts
}

struct CountTask {
    sequence: Vec<u8>,
    global_offset: usize,
    owned_start: usize,
    owned_end: usize,
}

fn count_sequence_chunk(task: &CountTask, k: u8, canonical: bool) -> FxHashMap<u64, u64> {
    let mut counts = FxHashMap::default();
    for_each_kmer_position(&task.sequence, k, canonical, |position, value| {
        let global_position = task.global_offset + position;
        if global_position >= task.owned_start && global_position < task.owned_end {
            *counts.entry(value).or_insert(0) += 1;
        }
    });
    counts
}

fn merge_counts(target: &mut FxHashMap<u64, u64>, source: FxHashMap<u64, u64>) {
    for (value, count) in source {
        *target.entry(value).or_insert(0) += count;
    }
}

/// Count only k-mers meeting `min_count`, using temporary hash shards on disk.
///
/// This keeps the in-memory table bounded approximately by `max_memory_mb` per
/// shard. It is intended for repeat calling, where low-abundance k-mers are
/// irrelevant and retaining every distinct k-mer is prohibitively expensive.
pub fn count_kmers_at_least(
    path: &Path,
    k: u8,
    canonical: bool,
    min_count: u64,
    max_memory_mb: usize,
) -> Result<FxHashMap<u64, u64>> {
    count_kmers_at_least_with_threads(path, k, canonical, min_count, max_memory_mb, 1)
}

pub fn count_kmers_at_least_with_threads(
    path: &Path,
    k: u8,
    canonical: bool,
    min_count: u64,
    max_memory_mb: usize,
    threads: usize,
) -> Result<FxHashMap<u64, u64>> {
    if k == 0 || k > 32 {
        return Err(Error::InvalidK(k));
    }
    if min_count <= 1 {
        return count_kmers_with_threads(path, k, canonical, threads);
    }

    let bytes = fs::metadata(path)?.len().max(1);
    let memory_bytes = (max_memory_mb.max(1) as u64) * 1024 * 1024;
    let estimated_entry_bytes = 40u64;
    let target_entries = (memory_bytes / estimated_entry_bytes).max(1);
    let shard_count = ((bytes / target_entries).max(16) as usize).min(4096);
    let temp_dir = temporary_shard_dir()?;
    let result = count_kmers_from_shards(
        path,
        k,
        canonical,
        min_count,
        shard_count,
        &temp_dir,
        threads,
    );
    let _ = fs::remove_dir_all(&temp_dir);
    result
}

fn count_kmers_from_shards(
    path: &Path,
    k: u8,
    canonical: bool,
    min_count: u64,
    shard_count: usize,
    temp_dir: &Path,
    threads: usize,
) -> Result<FxHashMap<u64, u64>> {
    let workers = threads.max(1);
    let mut writers = Vec::with_capacity(shard_count);
    for shard in 0..shard_count {
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(temp_dir.join(format!("shard-{shard:04}.bin")))?;
        writers.push(BufWriter::new(file));
    }

    let (task_sender, task_receiver) = bounded::<CountTask>(workers * 2);
    let (map_sender, map_receiver) = bounded::<FxHashMap<u64, u64>>(workers * 2);
    let mut parse_error = None;
    let mut staging_error = None;

    std::thread::scope(|scope| {
        let merger = scope.spawn(move || -> Result<()> {
            let mut writers = writers;
            for local_counts in map_receiver {
                for (value, count) in local_counts {
                    let shard = shard_index(value, shard_count);
                    writers[shard].write_all(&value.to_le_bytes())?;
                    writers[shard].write_all(&count.to_le_bytes())?;
                }
            }
            for writer in &mut writers {
                writer.flush()?;
            }
            Ok(())
        });

        for _ in 0..workers {
            let task_receiver = task_receiver.clone();
            let map_sender = map_sender.clone();
            scope.spawn(move || {
                while let Ok(task) = task_receiver.recv() {
                    if map_sender
                        .send(count_sequence_chunk(&task, k, canonical))
                        .is_err()
                    {
                        break;
                    }
                }
            });
        }
        drop(task_receiver);
        drop(map_sender);

        let mut reader = match parse_fastx_file(path).map_err(|e| Error::Other(e.to_string())) {
            Ok(reader) => reader,
            Err(error) => {
                parse_error = Some(error);
                drop(task_sender);
                return;
            }
        };
        while let Some(rec) = reader.next() {
            match rec {
                Ok(rec) => {
                    let seq = rec.seq();
                    let chunk_size = 4 * 1024 * 1024;
                    let mut chunk_start = 0usize;
                    while chunk_start < seq.len() {
                        let chunk_end = (chunk_start + chunk_size).min(seq.len());
                        let left = chunk_start.saturating_sub(k as usize - 1);
                        let right = (chunk_end + k as usize - 1).min(seq.len());
                        let task = CountTask {
                            sequence: seq[left..right].to_vec(),
                            global_offset: left,
                            owned_start: chunk_start,
                            owned_end: chunk_end,
                        };
                        if task_sender.send(task).is_err() {
                            break;
                        }
                        chunk_start = chunk_end;
                    }
                }
                Err(error) => {
                    parse_error = Some(Error::Other(error.to_string()));
                    break;
                }
            }
        }
        drop(task_sender);
        if let Err(error) = merger.join().expect("shard merger thread panicked") {
            staging_error = Some(error);
        }
    });

    if let Some(error) = parse_error {
        return Err(error);
    }
    if let Some(error) = staging_error {
        return Err(error);
    }

    let mut retained = FxHashMap::default();
    for shard in 0..shard_count {
        let path = temp_dir.join(format!("shard-{shard:04}.bin"));
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut counts = FxHashMap::default();
        let mut bytes = [0u8; 16];
        loop {
            match reader.read_exact(&mut bytes) {
                Ok(()) => {
                    let value = u64::from_le_bytes(bytes[..8].try_into().unwrap());
                    let count = u64::from_le_bytes(bytes[8..].try_into().unwrap());
                    *counts.entry(value).or_insert(0) += count;
                }
                Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(error) => return Err(error.into()),
            }
        }
        for (value, count) in counts {
            if count >= min_count {
                retained.insert(value, count);
            }
        }
    }
    Ok(retained)
}

fn shard_index(value: u64, shard_count: usize) -> usize {
    let mixed = value ^ value.rotate_right(29).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    (mixed as usize) % shard_count
}

fn temporary_shard_dir() -> Result<PathBuf> {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| Error::Other(e.to_string()))?
        .as_nanos();
    let path = std::env::temp_dir().join(format!("rusty-k-{}-{stamp}", std::process::id()));
    fs::create_dir(&path)?;
    Ok(path)
}

/// Sliding-window k-mer iterator that also yields the start position.
#[allow(dead_code)]
pub fn kmer_positions(seq: &[u8], k: u8, canonical: bool) -> Vec<(usize, u64)> {
    let mut out = Vec::new();
    for_each_kmer_position(seq, k, canonical, |pos, value| out.push((pos, value)));
    out
}

/// Visit each valid k-mer and its start position without allocating a result vector.
pub fn for_each_kmer_position<F>(seq: &[u8], k: u8, canonical: bool, mut visit: F)
where
    F: FnMut(usize, u64),
{
    if seq.len() < k as usize {
        return;
    }
    let mut current: u64 = 0;
    let mut valid = 0u8;
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    for (i, &b) in seq.iter().enumerate() {
        match encode_base(b) {
            Some(bits) => {
                current = ((current << 2) | bits) & mask;
                valid = valid.saturating_add(1);
                if valid >= k {
                    let pos = i + 1 - k as usize;
                    let kmer = if canonical {
                        let rc = revcomp(current, k);
                        current.min(rc)
                    } else {
                        current
                    };
                    visit(pos, kmer);
                }
            }
            None => {
                valid = 0;
                current = 0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn skips_ambiguous_bases() {
        let kmers = extract_kmers(b"ACGTNACG", 3, false);
        let decoded: Vec<_> = kmers.iter().map(|&km| decode_kmer(km, 3)).collect();
        assert_eq!(decoded, ["ACG", "CGT", "ACG"]);
    }

    #[test]
    fn canonical_count_collapses_reverse_complements() {
        let forward = extract_kmers(b"ACG", 3, true);
        let reverse = extract_kmers(b"CGT", 3, true);
        assert_eq!(forward, reverse);
    }

    #[test]
    fn supports_maximum_k() {
        let kmers = extract_kmers(&[b'A'; 32], 32, false);
        assert_eq!(kmers.len(), 1);
        assert_eq!(decode_kmer(kmers[0], 32), "A".repeat(32));
    }

    #[test]
    fn disk_counter_matches_thresholded_memory_counter() {
        let path = std::path::Path::new("test/ecoli_MG1655.fna");
        let expected: FxHashMap<_, _> = count_kmers(path, 5, true)
            .unwrap()
            .into_iter()
            .filter(|(_, count)| *count >= 100)
            .collect();
        let actual = count_kmers_at_least(path, 5, true, 100, 1).unwrap();
        assert_eq!(actual, expected);
    }
}

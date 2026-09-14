//! Repetitive-region detection from high-abundance k-mers.
//!
//! 1. Count all k-mers (canonical).
//! 2. For every position in the assembly, look up the count of the k-mer
//!    that starts there.
//! 3. Positions whose k-mer count ≥ min_count are marked “repetitive”.
//! 4. Consecutive marked positions are collapsed into intervals,
//!    gaps ≤ merge_gap are closed, and intervals shorter than min_len
//!    are discarded.

use crate::error::Result;
use crate::io::BedInterval;
use crate::kmer::{count_kmers, count_kmers_at_least, for_each_kmer_position};
use fxhash::FxHashMap;
use needletail::parse_fastx_file;
use std::collections::BinaryHeap;
use std::path::Path;

/// Returns (BED intervals, optional per-contig coverage vectors).
#[allow(dead_code)]
pub fn detect_repeats(
    path: &Path,
    k: u8,
    min_count: u64,
    min_len: u32,
    merge_gap: u32,
    canonical: bool,
    collect_coverage: bool,
) -> Result<(Vec<BedInterval>, Vec<(String, Vec<u64>)>)> {
    let counts = count_kmers(path, k, canonical)?;
    detect_repeats_with_counts(
        path,
        &counts,
        k,
        min_count,
        min_len,
        merge_gap,
        canonical,
        collect_coverage,
    )
}

/// Call repeats with a disk-backed threshold count table.
pub fn detect_repeats_with_memory_budget(
    path: &Path,
    k: u8,
    min_count: u64,
    min_len: u32,
    merge_gap: u32,
    canonical: bool,
    collect_coverage: bool,
    max_memory_mb: usize,
) -> Result<(Vec<BedInterval>, Vec<(String, Vec<u64>)>)> {
    let counts = count_kmers_at_least(path, k, canonical, min_count, max_memory_mb)?;
    detect_repeats_with_counts(
        path,
        &counts,
        k,
        min_count,
        min_len,
        merge_gap,
        canonical,
        collect_coverage,
    )
}

/// Call repetitive regions using already computed global k-mer counts.
pub fn detect_repeats_with_counts(
    path: &Path,
    counts: &FxHashMap<u64, u64>,
    k: u8,
    min_count: u64,
    min_len: u32,
    merge_gap: u32,
    canonical: bool,
    collect_coverage: bool,
) -> Result<(Vec<BedInterval>, Vec<(String, Vec<u64>)>)> {
    let mut reader =
        parse_fastx_file(path).map_err(|e| crate::error::Error::Other(e.to_string()))?;
    let mut all_ivs = Vec::new();
    let mut all_cov = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| crate::error::Error::Other(e.to_string()))?;
        let id = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq = rec.seq();
        if collect_coverage {
            let mut coverage = vec![0u64; seq.len()];
            stream_max_coverage(seq.as_ref(), k, canonical, counts, |pos, value| {
                coverage[pos] = value;
            });
            all_ivs.extend(coverage_to_intervals(
                &id, &coverage, min_count, min_len, merge_gap,
            ));
            all_cov.push((id, coverage));
        } else {
            all_ivs.extend(intervals_from_stream(
                &id,
                seq.as_ref(),
                k,
                canonical,
                counts,
                min_count,
                min_len,
                merge_gap,
            ));
        }
    }
    all_ivs.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    Ok((all_ivs, all_cov))
}

fn stream_max_coverage<F>(
    seq: &[u8],
    k: u8,
    canonical: bool,
    counts: &FxHashMap<u64, u64>,
    mut emit: F,
) where
    F: FnMut(usize, u64),
{
    let mut active = BinaryHeap::new();
    let mut cursor = 0usize;
    for_each_kmer_position(seq, k, canonical, |start, value| {
        for position in cursor..start {
            while active.peek().is_some_and(|&(_, end)| end <= position) {
                active.pop();
            }
            emit(position, active.peek().map_or(0, |&(count, _)| count));
        }
        while active.peek().is_some_and(|&(_, end)| end <= start) {
            active.pop();
        }
        active.push((*counts.get(&value).unwrap_or(&0), start + k as usize));
        emit(start, active.peek().map_or(0, |&(count, _)| count));
        cursor = start + 1;
    });
    for position in cursor..seq.len() {
        while active.peek().is_some_and(|&(_, end)| end <= position) {
            active.pop();
        }
        emit(position, active.peek().map_or(0, |&(count, _)| count));
    }
}

fn intervals_from_stream(
    chrom: &str,
    seq: &[u8],
    k: u8,
    canonical: bool,
    counts: &FxHashMap<u64, u64>,
    min_count: u64,
    min_len: u32,
    merge_gap: u32,
) -> Vec<BedInterval> {
    let mut intervals = Vec::new();
    let mut current_start = None;
    let mut current_end = 0usize;
    let mut current_score = 0u64;
    let mut last_high = None;
    let mut emit = |position: usize, value: u64| {
        if value >= min_count {
            if current_start.is_none() {
                current_start = Some(position);
            }
            current_end = position + 1;
            current_score = current_score.max(value);
            last_high = Some(position);
        } else if let (Some(start), Some(last)) = (current_start, last_high) {
            if position.saturating_sub(last) > merge_gap as usize {
                if current_end - start >= min_len as usize {
                    intervals.push(BedInterval {
                        chrom: chrom.to_string(),
                        start: start as u64,
                        end: current_end as u64,
                        name: "repetitive".into(),
                        score: current_score,
                        strand: ".",
                    });
                }
                current_start = None;
                current_score = 0;
            }
        }
    };
    stream_max_coverage(seq, k, canonical, counts, &mut emit);
    if let Some(start) = current_start {
        if current_end - start >= min_len as usize {
            intervals.push(BedInterval {
                chrom: chrom.to_string(),
                start: start as u64,
                end: current_end as u64,
                name: "repetitive".into(),
                score: current_score,
                strand: ".",
            });
        }
    }
    intervals
}

fn coverage_to_intervals(
    chrom: &str,
    cov: &[u64],
    min_count: u64,
    min_len: u32,
    merge_gap: u32,
) -> Vec<BedInterval> {
    let mut raw = Vec::new();
    let mut i = 0usize;
    while i < cov.len() {
        if cov[i] >= min_count {
            let start = i;
            while i < cov.len() && cov[i] >= min_count {
                i += 1;
            }
            raw.push((start, i)); // [start, end)
        } else {
            i += 1;
        }
    }

    // Merge gaps ≤ merge_gap
    if raw.is_empty() {
        return Vec::new();
    }
    let mut merged = Vec::new();
    let (mut s, mut e) = raw[0];
    for &(ns, ne) in raw.iter().skip(1) {
        if ns <= e + merge_gap as usize {
            e = e.max(ne);
        } else {
            if e - s >= min_len as usize {
                merged.push(BedInterval {
                    chrom: chrom.to_string(),
                    start: s as u64,
                    end: e as u64,
                    name: "repetitive".into(),
                    score: *cov[s..e].iter().max().unwrap_or(&0),
                    strand: ".",
                });
            }
            s = ns;
            e = ne;
        }
    }
    if e - s >= min_len as usize {
        merged.push(BedInterval {
            chrom: chrom.to_string(),
            start: s as u64,
            end: e as u64,
            name: "repetitive".into(),
            score: *cov[s..e].iter().max().unwrap_or(&0),
            strand: ".",
        });
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_short_gaps_and_keeps_peak_score() {
        let coverage = [5, 5, 0, 5, 5];
        let intervals = coverage_to_intervals("ctg", &coverage, 5, 4, 1);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 0);
        assert_eq!(intervals[0].end, 5);
        assert_eq!(intervals[0].score, 5);
    }

    #[test]
    fn filters_short_regions() {
        let coverage = [5, 5, 0, 5, 5];
        let intervals = coverage_to_intervals("ctg", &coverage, 5, 3, 0);
        assert!(intervals.is_empty());
    }
}

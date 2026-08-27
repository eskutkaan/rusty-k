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
use crate::kmer::{count_kmers, kmer_positions};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::path::Path;

/// Returns (BED intervals, optional per-contig coverage vectors).
pub fn detect_repeats(
    path: &Path,
    k: u8,
    min_count: u64,
    min_len: u32,
    merge_gap: u32,
    canonical: bool,
) -> Result<(Vec<BedInterval>, Vec<(String, Vec<u64>)>)> {
    // Global k-mer counts
    let counts = count_kmers(path, k, canonical)?;

    // Load sequences
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = parse_fastx_file(path).map_err(|e| crate::error::Error::Other(e.to_string()))?;
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| crate::error::Error::Other(e.to_string()))?;
        let id = String::from_utf8_lossy(rec.id()).split_whitespace().next().unwrap_or("").to_string();
        records.push((id, rec.seq().to_vec()));
    }

    let results: Vec<(Vec<BedInterval>, (String, Vec<u64>))> = records
        .par_iter()
        .map(|(id, seq)| {
            let positions = kmer_positions(seq, k, canonical);
            let mut coverage = vec![0u64; seq.len()];

            for &(pos, km) in &positions {
                let c = *counts.get(&km).unwrap_or(&0);
                // annotate the whole k-mer span
                for p in pos..pos + k as usize {
                    if p < coverage.len() {
                        coverage[p] = coverage[p].max(c);
                    }
                }
            }

            let intervals = coverage_to_intervals(id, &coverage, min_count, min_len, merge_gap);
            (intervals, (id.clone(), coverage))
        })
        .collect();

    let mut all_ivs = Vec::new();
    let mut all_cov = Vec::new();
    for (ivs, cov) in results {
        all_ivs.extend(ivs);
        all_cov.push(cov);
    }
    Ok((all_ivs, all_cov))
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

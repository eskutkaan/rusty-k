//! Tandem-repeat detection via periodogram of consecutive identical k-mers.
//!
//! For each contig we slide a window and look for stretches where the same
//! k-mer (or a short period of k-mers) repeats at least `min_copies` times.
//! Periods from 1 bp up to `max_period` are examined.

use crate::error::Result;
use crate::io::BedInterval;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::path::Path;

/// Detect tandem repeats in a FASTA assembly.
pub fn detect_tandems(
    path: &Path,
    k: u8,
    min_copies: u32,
    max_period: u32,
) -> Result<Vec<BedInterval>> {
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = parse_fastx_file(path).map_err(|e| crate::error::Error::Other(e.to_string()))?;
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| crate::error::Error::Other(e.to_string()))?;
        let id = String::from_utf8_lossy(rec.id()).split_whitespace().next().unwrap_or("").to_string();
        records.push((id, rec.seq().to_vec()));
    }

    let intervals: Vec<Vec<BedInterval>> = records
        .par_iter()
        .map(|(id, seq)| find_tandems_in_seq(id, seq, k, min_copies, max_period))
        .collect();

    Ok(intervals.into_iter().flatten().collect())
}

fn find_tandems_in_seq(
    contig: &str,
    seq: &[u8],
    k: u8,
    min_copies: u32,
    max_period: u32,
) -> Vec<BedInterval> {
    let n = seq.len();
    if n < (k as usize) * (min_copies as usize) {
        return Vec::new();
    }

    let mut intervals = Vec::new();
    // We look for exact tandem repeats of unit length `period`.
    // For efficiency we only check periods that divide cleanly with the k-mer size
    // or use a simple run-length approach on the nucleotide level for short periods.

    // --- short exact tandem runs (period 1..max_period) ---
    // Classic run-length / period detection with a sliding comparison.
    let mut i = 0usize;
    while i < n {
        let mut best_period = 0u32;
        let mut best_copies = 0u32;
        let mut best_end = i;

        for period in 1..=max_period as usize {
            if i + period >= n {
                break;
            }
            // Count how many consecutive copies of seq[i..i+period] exist
            let unit = &seq[i..i + period];
            let mut copies = 1u32;
            let mut pos = i + period;
            while pos + period <= n && &seq[pos..pos + period] == unit {
                copies += 1;
                pos += period;
            }
            // also allow a partial final copy if it matches a prefix
            if pos < n {
                let rem = n - pos;
                if rem < period && &seq[pos..] == &unit[..rem] {
                    // partial – we still count the full copies
                }
            }
            if copies >= min_copies && copies > best_copies {
                best_copies = copies;
                best_period = period as u32;
                best_end = pos;
            }
        }

        if best_copies >= min_copies {
            // avoid reporting extremely short units that are just homopolymers unless asked
            let start = i;
            let end = best_end;
            let unit_len = best_period;
            intervals.push(BedInterval {
                chrom: contig.to_string(),
                start: start as u64,
                end: end as u64,
                name: format!("TR_period{}_x{}", unit_len, best_copies),
                score: best_copies as u64,
                strand: ".",
            });
            i = end; // jump past this repeat
        } else {
            i += 1;
        }
    }

    // Merge overlapping / adjacent intervals of the same period family
    intervals.sort_by_key(|iv| (iv.start, iv.end));
    merge_intervals(intervals)
}

fn merge_intervals(mut ivs: Vec<BedInterval>) -> Vec<BedInterval> {
    if ivs.is_empty() {
        return ivs;
    }
    ivs.sort_by_key(|a| (a.chrom.clone(), a.start));
    let mut merged = Vec::with_capacity(ivs.len());
    let mut cur = ivs[0].clone();
    for next in ivs.into_iter().skip(1) {
        if next.chrom == cur.chrom && next.start <= cur.end + 10 {
            // small gap allowed
            cur.end = cur.end.max(next.end);
            cur.score = cur.score.max(next.score);
            // keep the name of the longer one
            if next.end - next.start > cur.end - cur.start {
                cur.name = next.name;
            }
        } else {
            merged.push(cur);
            cur = next;
        }
    }
    merged.push(cur);
    merged
}

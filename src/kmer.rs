//! Canonical k-mer encoding (k ≤ 32) and parallel counting.

use crate::error::{Error, Result};
use fxhash::FxHashMap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::path::Path;

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
    if k == 0 || k > 32 {
        return Err(Error::InvalidK(k));
    }

    // Collect sequences first (assemblies are usually modest in size)
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = parse_fastx_file(path).map_err(|e| Error::Other(e.to_string()))?;
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| Error::Other(e.to_string()))?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let seq = rec.seq().to_vec();
        records.push((id, seq));
    }

    // Parallel extraction + local maps, then merge
    let partials: Vec<FxHashMap<u64, u64>> = records
        .par_iter()
        .map(|(_, seq)| {
            let mut map = FxHashMap::default();
            for km in extract_kmers(seq, k, canonical) {
                *map.entry(km).or_insert(0) += 1;
            }
            map
        })
        .collect();

    let mut global = FxHashMap::default();
    for m in partials {
        for (k, v) in m {
            *global.entry(k).or_insert(0) += v;
        }
    }
    Ok(global)
}

/// Sliding-window k-mer iterator that also yields the start position.
pub fn kmer_positions(seq: &[u8], k: u8, canonical: bool) -> Vec<(usize, u64)> {
    if seq.len() < k as usize {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(seq.len() - k as usize + 1);
    let mut current: u64 = 0;
    let mut valid = 0u8;
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

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
                    out.push((pos, kmer));
                }
            }
            None => {
                valid = 0;
                current = 0;
            }
        }
    }
    out
}

//! Input / output helpers (BED, TSV, JSON, coverage).

use crate::error::Result;
use crate::kmer::decode_kmer;
use fxhash::FxHashMap;
use serde::Serialize;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug, Serialize)]
pub struct BedInterval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: u64,
    pub strand: &'static str,
}

pub fn write_bed(path: &Path, intervals: &[BedInterval]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for iv in intervals {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            iv.chrom, iv.start, iv.end, iv.name, iv.score, iv.strand
        )?;
    }
    Ok(())
}

pub fn write_kmer_counts(
    path: &Path,
    counts: &FxHashMap<u64, u64>,
    min_count: u64,
    as_json: bool,
) -> Result<()> {
    // We need the k value to decode; we infer a reasonable k from the highest bit used.
    // Caller should guarantee consistent k; we store k in the first line as a comment for TSV.
    let mut pairs: Vec<(u64, u64)> = counts
        .iter()
        .filter(|(_, &c)| c >= min_count)
        .map(|(&k, &c)| (k, c))
        .collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    if as_json {
        // decode on the fly – we guess k from the bit length of the largest kmer
        let max_kmer = pairs.iter().map(|(k, _)| *k).max().unwrap_or(0);
        let k = if max_kmer == 0 {
            21
        } else {
            ((64 - max_kmer.leading_zeros() + 1) / 2).max(1) as u8
        };
        let objs: Vec<_> = pairs
            .iter()
            .map(|(km, c)| {
                serde_json::json!({
                    "kmer": decode_kmer(*km, k),
                    "count": c
                })
            })
            .collect();
        let f = File::create(path)?;
        serde_json::to_writer_pretty(f, &objs)?;
    } else {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "kmer\tcount")?;
        // Infer k the same way
        let max_kmer = pairs.iter().map(|(k, _)| *k).max().unwrap_or(0);
        let k = if max_kmer == 0 {
            21
        } else {
            ((64 - max_kmer.leading_zeros() + 1) / 2).max(1) as u8
        };
        for (km, c) in pairs {
            writeln!(w, "{}\t{}", decode_kmer(km, k), c)?;
        }
    }
    Ok(())
}

pub fn write_coverage(path: &Path, cov: &[(String, Vec<u64>)]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "contig\tposition\tcoverage")?;
    for (id, v) in cov {
        for (i, &c) in v.iter().enumerate() {
            if c > 0 {
                writeln!(w, "{}\t{}\t{}", id, i, c)?;
            }
        }
    }
    Ok(())
}

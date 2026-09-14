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
    k: u8,
    min_count: u64,
    as_json: bool,
) -> Result<()> {
    let mut pairs: Vec<(u64, u64)> = counts
        .iter()
        .filter(|(_, &c)| c >= min_count)
        .map(|(&k, &c)| (k, c))
        .collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    if as_json {
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

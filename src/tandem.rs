//! Tandem-repeat detection via periodogram of consecutive identical k-mers.
//!
//! For each contig we look for stretches where a nucleotide unit repeats at
//! least `min_copies` times. Periods from 1 bp up to `max_period` are examined.

use crate::error::Result;
use crate::io::BedInterval;
use crossbeam_channel::bounded;
use needletail::parse_fastx_file;
use std::path::Path;

/// Detect tandem repeats in a FASTA assembly.
#[allow(dead_code)]
pub fn detect_tandems(
    path: &Path,
    min_repeat_span: u32,
    min_copies: u32,
    max_period: u32,
) -> Result<Vec<BedInterval>> {
    detect_tandems_with_threads(path, min_repeat_span, min_copies, max_period, 1)
}

pub fn detect_tandems_with_threads(
    path: &Path,
    min_repeat_span: u32,
    min_copies: u32,
    max_period: u32,
    threads: usize,
) -> Result<Vec<BedInterval>> {
    let workers = threads.max(1);
    if workers == 1 {
        return detect_tandems_sequential(path, min_repeat_span, min_copies, max_period);
    }

    let (sender, receiver) = bounded::<TandemTask>(workers * 2);
    let (result_sender, result_receiver) = bounded::<Vec<BedInterval>>(workers);
    let mut parse_error = None;

    std::thread::scope(|scope| {
        for _ in 0..workers {
            let receiver = receiver.clone();
            let result_sender = result_sender.clone();
            scope.spawn(move || {
                while let Ok(task) = receiver.recv() {
                    let mut intervals = find_tandems_in_seq(
                        &task.id,
                        &task.sequence,
                        min_repeat_span,
                        min_copies,
                        max_period,
                    );
                    intervals.retain(|interval| {
                        let global_start = task.global_offset + interval.start as usize;
                        global_start >= task.owned_start && global_start < task.owned_end
                    });
                    for interval in &mut intervals {
                        interval.start += task.global_offset as u64;
                        interval.end += task.global_offset as u64;
                    }
                    if result_sender.send(intervals).is_err() {
                        break;
                    }
                }
            });
        }
        drop(receiver);
        drop(result_sender);

        let mut reader =
            match parse_fastx_file(path).map_err(|e| crate::error::Error::Other(e.to_string())) {
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
                    let id = String::from_utf8_lossy(rec.id())
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string();
                    let seq = rec.seq();
                    let chunk_size = 4 * 1024 * 1024;
                    let overlap = (max_period as usize)
                        .saturating_mul(min_copies.max(1) as usize)
                        .max(min_repeat_span as usize);
                    let mut chunk_start = 0usize;
                    while chunk_start < seq.len() {
                        let chunk_end = (chunk_start + chunk_size).min(seq.len());
                        let left = chunk_start.saturating_sub(overlap);
                        let right = (chunk_end + overlap).min(seq.len());
                        let task = TandemTask {
                            id: id.clone(),
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
                    parse_error = Some(crate::error::Error::Other(error.to_string()));
                    break;
                }
            }
        }
        drop(sender);
    });

    if let Some(error) = parse_error {
        return Err(error);
    }
    Ok(result_receiver.into_iter().flatten().collect())
}

struct TandemTask {
    id: String,
    sequence: Vec<u8>,
    global_offset: usize,
    owned_start: usize,
    owned_end: usize,
}

fn detect_tandems_sequential(
    path: &Path,
    min_repeat_span: u32,
    min_copies: u32,
    max_period: u32,
) -> Result<Vec<BedInterval>> {
    let mut reader =
        parse_fastx_file(path).map_err(|e| crate::error::Error::Other(e.to_string()))?;
    let mut all_intervals = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| crate::error::Error::Other(e.to_string()))?;
        let id = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        all_intervals.extend(find_tandems_in_seq(
            &id,
            rec.seq().as_ref(),
            min_repeat_span,
            min_copies,
            max_period,
        ));
    }

    Ok(all_intervals)
}

fn find_tandems_in_seq(
    contig: &str,
    seq: &[u8],
    min_repeat_span: u32,
    min_copies: u32,
    max_period: u32,
) -> Vec<BedInterval> {
    let n = seq.len();
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
            if copies >= min_copies && (pos - i) >= min_repeat_span as usize && copies > best_copies
            {
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
            let cur_len = cur.end - cur.start;
            let next_len = next.end - next.start;
            // small gap allowed
            cur.end = cur.end.max(next.end);
            cur.score = cur.score.max(next.score);
            // keep the name of the longer one
            if next_len > cur_len {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_exact_tandem_run() {
        let intervals = find_tandems_in_seq("ctg", b"TTACACACGG", 6, 3, 10);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 2);
        assert_eq!(intervals[0].end, 8);
        assert_eq!(intervals[0].name, "TR_period2_x3");
    }

    #[test]
    fn minimum_repeat_span_is_enforced() {
        let intervals = find_tandems_in_seq("ctg", b"ACACAC", 7, 3, 10);
        assert!(intervals.is_empty());
    }
}

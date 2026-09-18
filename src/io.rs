//! Streaming k-mer count output helpers.

use crate::error::Result;
use crate::kmer::decode_kmer;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

pub struct KmerWriter {
    writer: BufWriter<File>,
    output_path: PathBuf,
    temporary_path: PathBuf,
    sorter: Option<SortRuns>,
    k: u8,
    as_json: bool,
    first: bool,
}

impl KmerWriter {
    pub fn create(
        path: &Path,
        k: u8,
        as_json: bool,
        sort: bool,
        sort_temp_dir: &Path,
    ) -> Result<Self> {
        let temporary_path = temporary_output_path(path)?;
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary_path)?;
        let mut writer = BufWriter::new(file);
        if !sort {
            write_header(&mut writer, as_json)?;
        }
        Ok(Self {
            writer,
            output_path: path.to_path_buf(),
            temporary_path,
            sorter: sort.then(|| SortRuns::new(sort_temp_dir.to_path_buf())),
            k,
            as_json,
            first: true,
        })
    }

    pub fn write_count(&mut self, kmer: u64, count: u64) -> Result<()> {
        if let Some(sorter) = &mut self.sorter {
            sorter.push(kmer, count)?;
        } else {
            write_count(
                &mut self.writer,
                self.k,
                self.as_json,
                &mut self.first,
                kmer,
                count,
            )?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        if let Some(sorter) = &mut self.sorter {
            write_header(&mut self.writer, self.as_json)?;
            sorter.finish(&mut self.writer, self.k, self.as_json, &mut self.first)?;
        } else if self.as_json {
            writeln!(self.writer, "]")?;
        }
        self.writer.flush()?;
        fs::rename(&self.temporary_path, &self.output_path)?;
        Ok(())
    }
}

impl Drop for KmerWriter {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.temporary_path);
        if let Some(sorter) = &mut self.sorter {
            sorter.cleanup();
        }
    }
}

const SORT_BUFFER_RECORDS: usize = 65_536;

struct SortRuns {
    temp_dir: PathBuf,
    buffer: Vec<(u64, u64)>,
    runs: Vec<PathBuf>,
}

impl SortRuns {
    fn new(temp_dir: PathBuf) -> Self {
        Self {
            temp_dir,
            buffer: Vec::with_capacity(SORT_BUFFER_RECORDS),
            runs: Vec::new(),
        }
    }

    fn push(&mut self, kmer: u64, count: u64) -> Result<()> {
        self.buffer.push((kmer, count));
        if self.buffer.len() >= SORT_BUFFER_RECORDS {
            self.flush_run()?;
        }
        Ok(())
    }

    fn flush_run(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer.sort_unstable_by_key(|&(kmer, _)| kmer);
        let path = self.temp_dir.join(format!(
            "rusty-k-sort-{}-{}-{}.bin",
            std::process::id(),
            self.runs.len(),
            timestamp_nanos()?
        ));
        let mut writer = BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)?,
        );
        for &(kmer, count) in &self.buffer {
            writer.write_all(&kmer.to_le_bytes())?;
            writer.write_all(&count.to_le_bytes())?;
        }
        writer.flush()?;
        self.runs.push(path);
        self.buffer.clear();
        Ok(())
    }

    fn finish(
        &mut self,
        writer: &mut BufWriter<File>,
        k: u8,
        as_json: bool,
        first: &mut bool,
    ) -> Result<()> {
        self.flush_run()?;
        let mut readers = Vec::with_capacity(self.runs.len());
        let mut heap = BinaryHeap::new();
        for (run, path) in self.runs.iter().enumerate() {
            let mut reader = BufReader::new(File::open(path)?);
            if let Some((kmer, count)) = read_record(&mut reader)? {
                heap.push(HeapItem { kmer, count, run });
            }
            readers.push(reader);
        }
        while let Some(item) = heap.pop() {
            write_count(writer, k, as_json, first, item.kmer, item.count)?;
            if let Some((kmer, count)) = read_record(&mut readers[item.run])? {
                heap.push(HeapItem {
                    kmer,
                    count,
                    run: item.run,
                });
            }
        }
        Ok(())
    }

    fn cleanup(&mut self) {
        for path in self.runs.drain(..) {
            let _ = fs::remove_file(path);
        }
    }
}

impl Drop for SortRuns {
    fn drop(&mut self) {
        self.cleanup();
    }
}

struct HeapItem {
    kmer: u64,
    count: u64,
    run: usize,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .kmer
            .cmp(&self.kmer)
            .then_with(|| other.run.cmp(&self.run))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer && self.run == other.run
    }
}

impl Eq for HeapItem {}

fn read_record(reader: &mut BufReader<File>) -> Result<Option<(u64, u64)>> {
    let mut kmer_bytes = [0u8; 8];
    match reader.read_exact(&mut kmer_bytes) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(error) => return Err(error.into()),
    }
    let mut count_bytes = [0u8; 8];
    reader.read_exact(&mut count_bytes)?;
    Ok(Some((
        u64::from_le_bytes(kmer_bytes),
        u64::from_le_bytes(count_bytes),
    )))
}

fn write_header(writer: &mut BufWriter<File>, as_json: bool) -> Result<()> {
    if as_json {
        write!(writer, "[")?;
    } else {
        writeln!(writer, "kmer\tcount")?;
    }
    Ok(())
}

fn write_count(
    writer: &mut BufWriter<File>,
    k: u8,
    as_json: bool,
    first: &mut bool,
    kmer: u64,
    count: u64,
) -> Result<()> {
    if as_json {
        if !*first {
            write!(writer, ",")?;
        }
        write!(
            writer,
            "{{\"kmer\":\"{}\",\"count\":{}}}",
            decode_kmer(kmer, k),
            count
        )?;
    } else {
        writeln!(writer, "{}\t{}", decode_kmer(kmer, k), count)?;
    }
    *first = false;
    Ok(())
}

fn timestamp_nanos() -> Result<u128> {
    Ok(SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| crate::error::Error::Other(error.to_string()))?
        .as_nanos())
}

fn temporary_output_path(path: &Path) -> Result<PathBuf> {
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    let name = path
        .file_name()
        .ok_or_else(|| crate::error::Error::Other("output path has no file name".into()))?
        .to_string_lossy();
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| crate::error::Error::Other(error.to_string()))?
        .as_nanos();
    Ok(parent.join(format!(".{name}.rusty-k-{stamp}-{}", std::process::id())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorted_writer_orders_kmers_lexicographically() {
        let root = std::env::temp_dir().join(format!("rusty-k-writer-{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let output = root.join("output.tsv");
        let mut writer = KmerWriter::create(&output, 1, false, true, &root).unwrap();
        writer.write_count(3, 1).unwrap();
        writer.write_count(0, 2).unwrap();
        writer.write_count(1, 3).unwrap();
        writer.finish().unwrap();

        let contents = fs::read_to_string(&output).unwrap();
        assert_eq!(contents, "kmer\tcount\nA\t2\nC\t3\nT\t1\n");
        fs::remove_dir_all(root).unwrap();
    }
}

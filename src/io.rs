//! Streaming k-mer count output helpers.

use crate::error::Result;
use crate::kmer::decode_kmer;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct KmerWriter {
    writer: BufWriter<File>,
    k: u8,
    as_json: bool,
    first: bool,
}

impl KmerWriter {
    pub fn create(path: &Path, k: u8, as_json: bool) -> Result<Self> {
        let mut writer = BufWriter::new(File::create(path)?);
        if as_json {
            write!(writer, "[")?;
        } else {
            writeln!(writer, "kmer\tcount")?;
        }
        Ok(Self {
            writer,
            k,
            as_json,
            first: true,
        })
    }

    pub fn write_count(&mut self, kmer: u64, count: u64) -> Result<()> {
        if self.as_json {
            if !self.first {
                write!(self.writer, ",")?;
            }
            write!(
                self.writer,
                "{{\"kmer\":\"{}\",\"count\":{}}}",
                decode_kmer(kmer, self.k),
                count
            )?;
        } else {
            writeln!(self.writer, "{}\t{}", decode_kmer(kmer, self.k), count)?;
        }
        self.first = false;
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        if self.as_json {
            writeln!(self.writer, "]")?;
        }
        Ok(())
    }
}

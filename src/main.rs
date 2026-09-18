//! Fast k-mer counting for genome assemblies.

mod error;
mod io;
mod kmer;

use clap::Parser;
use log::LevelFilter;
use std::fs;
use std::path::{Path, PathBuf};

use crate::error::Result;

#[derive(Parser, Debug)]
#[command(name = "rusty-k", version, about = "Fast k-mer counter")]
struct Cli {
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
    #[arg(short, long)]
    input: PathBuf,
    #[arg(short, long, default_value_t = 21)]
    k: u8,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long)]
    json: bool,
    /// Sort output k-mers lexicographically using bounded external sorting.
    #[arg(long)]
    sort: bool,
    #[arg(long, default_value_t = 1)]
    min_count: u64,
    /// Approximate RAM budget per disk-counting shard in MiB.
    #[arg(long, default_value_t = 4096)]
    max_memory_mb: usize,
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    canonical: bool,
    /// Number of counting workers; 0 uses all available CPUs.
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
    /// Directory for temporary shard files. Defaults to tmp beside the executable.
    #[arg(long, value_name = "DIR")]
    tmp_dir: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if paths_refer_to_same_file(&cli.input, &cli.output)? {
        return Err(crate::error::Error::Other(
            "input and output must be different files".into(),
        ));
    }
    let level = match cli.verbose {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    env_logger::Builder::new()
        .filter_level(level)
        .format_timestamp_secs()
        .init();

    let threads = if cli.threads == 0 {
        num_cpus::get()
    } else {
        cli.threads
    };
    let temp_dir = temporary_parent(cli.tmp_dir.as_deref())?;
    fs::create_dir_all(&temp_dir)?;
    let mut output = io::KmerWriter::create(&cli.output, cli.k, cli.json, cli.sort, &temp_dir)?;
    let written = kmer::count_kmers_streaming_with_temp_dir(
        &cli.input,
        cli.k,
        cli.canonical,
        cli.min_count,
        cli.max_memory_mb,
        threads,
        Some(&temp_dir),
        |kmer, count| output.write_count(kmer, count),
    )?;
    output.finish()?;
    log::info!("Wrote {} distinct k-mers to {:?}", written, cli.output);
    Ok(())
}

fn temporary_parent(path: Option<&Path>) -> Result<PathBuf> {
    match path {
        Some(path) => Ok(path.to_path_buf()),
        None => Ok(std::env::current_exe()?
            .parent()
            .ok_or_else(|| crate::error::Error::Other("executable has no parent directory".into()))?
            .join("tmp")),
    }
}

fn paths_refer_to_same_file(input: &std::path::Path, output: &std::path::Path) -> Result<bool> {
    let input = fs::canonicalize(input)?;
    let output = if output.exists() {
        fs::canonicalize(output)?
    } else {
        let parent = output.parent().unwrap_or_else(|| std::path::Path::new("."));
        fs::canonicalize(parent)?.join(
            output
                .file_name()
                .ok_or_else(|| crate::error::Error::Other("output path has no file name".into()))?,
        )
    };
    Ok(input == output)
}

//! Fast k-mer counting for genome assemblies.

mod error;
mod io;
mod kmer;

use clap::Parser;
use log::LevelFilter;
use std::path::PathBuf;

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
}

fn main() -> Result<()> {
    let cli = Cli::parse();
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
    let mut output = io::KmerWriter::create(&cli.output, cli.k, cli.json)?;
    let written = kmer::count_kmers_streaming(
        &cli.input,
        cli.k,
        cli.canonical,
        cli.min_count,
        cli.max_memory_mb,
        threads,
        |kmer, count| output.write_count(kmer, count),
    )?;
    output.finish()?;
    log::info!("Wrote {} distinct k-mers to {:?}", written, cli.output);
    Ok(())
}

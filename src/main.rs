//! rusty-k – efficient k-mer counting, tandem-repeat detection,
//! and repetitive-region calling for genome assemblies.

mod error;
mod io;
mod kmer;
mod repeats;
mod tandem;

use clap::{Parser, Subcommand};
use log::{info, LevelFilter};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

use crate::error::Result;

#[derive(Parser, Debug)]
#[command(name = "rusty-k", version, about = "K-mer counter with tandem-repeat and repetitive-region detection", long_about = None)]
struct Cli {
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Count {
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
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        canonical: bool,
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },
    Tandem {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value_t = 11)]
        k: u32,
        #[arg(long, default_value_t = 3)]
        min_copies: u32,
        #[arg(long, default_value_t = 100)]
        max_period: u32,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },
    Repeats {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value_t = 21)]
        k: u8,
        #[arg(long, default_value_t = 5)]
        min_count: u64,
        #[arg(long, default_value_t = 100)]
        min_len: u32,
        #[arg(long, default_value_t = 50)]
        merge_gap: u32,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        coverage: Option<PathBuf>,
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        canonical: bool,
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },
    All {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value_t = 21)]
        k: u8,
        #[arg(long, default_value_t = 11)]
        tandem_k: u32,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value_t = 5)]
        min_count: u64,
        #[arg(long, default_value_t = 3)]
        tandem_min_copies: u32,
        #[arg(long, default_value_t = 100)]
        max_period: u32,
        #[arg(long, default_value_t = 100)]
        min_len: u32,
        #[arg(long, default_value_t = 50)]
        merge_gap: u32,
        #[arg(long)]
        coverage: Option<PathBuf>,
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        canonical: bool,
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },
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

    match cli.command {
        Commands::Count {
            input,
            k,
            output,
            json,
            min_count,
            canonical,
            threads,
        } => {
            set_threads(threads);
            let counts = kmer::count_kmers(&input, k, canonical)?;
            io::write_kmer_counts(&output, &counts, k, min_count, json)?;
            info!("Wrote {} distinct k-mers to {:?}", counts.len(), output);
        }
        Commands::Tandem {
            input,
            k,
            min_copies,
            max_period,
            output,
            threads,
        } => {
            set_threads(threads);
            let intervals = tandem::detect_tandems(&input, k, min_copies, max_period)?;
            io::write_bed(&output, &intervals)?;
            info!("Found {} tandem-repeat intervals", intervals.len());
        }
        Commands::Repeats {
            input,
            k,
            min_count,
            min_len,
            merge_gap,
            output,
            coverage,
            canonical,
            threads,
        } => {
            set_threads(threads);
            let (intervals, cov) =
                repeats::detect_repeats(&input, k, min_count, min_len, merge_gap, canonical)?;
            io::write_bed(&output, &intervals)?;
            if let Some(path) = coverage {
                io::write_coverage(&path, &cov)?;
            }
            info!("Found {} repetitive intervals", intervals.len());
        }
        Commands::All {
            input,
            k,
            tandem_k,
            output,
            min_count,
            tandem_min_copies,
            max_period,
            min_len,
            merge_gap,
            coverage,
            canonical,
            threads,
        } => {
            set_threads(threads);
            std::fs::create_dir_all(&output)?;
            let counts = kmer::count_kmers(&input, k, canonical)?;
            io::write_kmer_counts(&output.join("kmers.tsv"), &counts, k, 1, false)?;
            let tandems = tandem::detect_tandems(&input, tandem_k, tandem_min_copies, max_period)?;
            io::write_bed(&output.join("tandems.bed"), &tandems)?;
            let (reps, cov) = repeats::detect_repeats_with_counts(
                &input, &counts, k, min_count, min_len, merge_gap, canonical,
            )?;
            io::write_bed(&output.join("repeats.bed"), &reps)?;
            if let Some(path) = coverage {
                io::write_coverage(&path, &cov)?;
            }
            let summary = serde_json::json!({
                "schema_version": 1,
                "software": { "name": env!("CARGO_PKG_NAME"), "version": env!("CARGO_PKG_VERSION") },
                "input": input_provenance(&input)?,
                "threads": if threads == 0 { num_cpus::get() } else { threads },
                "k": k,
                "tandem_k": tandem_k,
                "canonical": canonical,
                "tandem_min_copies": tandem_min_copies,
                "max_period": max_period,
                "min_count": min_count,
                "min_len": min_len,
                "merge_gap": merge_gap,
                "distinct_kmers": counts.len(),
                "tandem_intervals": tandems.len(),
                "repetitive_intervals": reps.len(),
            });
            std::fs::write(
                output.join("summary.json"),
                serde_json::to_string_pretty(&summary)?,
            )?;
            info!("Analysis outputs written to {:?}", output);
        }
    }
    Ok(())
}

fn set_threads(n: usize) {
    let n = if n == 0 { num_cpus::get() } else { n };
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .ok();
    info!("Using {} threads", n);
}

fn input_provenance(path: &std::path::Path) -> Result<serde_json::Value> {
    let metadata = std::fs::metadata(path)?;
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let bytes_read = file.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }
    Ok(serde_json::json!({
        "path": path.display().to_string(),
        "size_bytes": metadata.len(),
        "sha256": format!("{:x}", hasher.finalize()),
    }))
}

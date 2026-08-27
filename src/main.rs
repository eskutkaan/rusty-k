//! kmer_tool – efficient k-mer counting, tandem-repeat detection,
//! and repetitive-region calling for genome assemblies.

mod kmer;
mod tandem;
mod repeats;
mod io;
mod error;

use clap::{Parser, Subcommand};
use log::{info, LevelFilter};
use std::path::PathBuf;

use crate::error::Result;

#[derive(Parser, Debug)]
#[command(
    name = "kmer_tool",
    version,
    about = "K-mer counter with tandem-repeat and repetitive-region detection",
    long_about = None
)]
struct Cli {
    /// Increase verbosity (-v, -vv)
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Count k-mers in a FASTA/FASTQ assembly or read set
    Count {
        /// Input FASTA/FASTQ (supports .gz)
        #[arg(short, long)]
        input: PathBuf,

        /// k-mer length (1–32 recommended for 64-bit encoding)
        #[arg(short, long, default_value_t = 21)]
        k: u8,

        /// Output TSV of k-mer counts (or JSON with --json)
        #[arg(short, long)]
        output: PathBuf,

        /// Write JSON instead of TSV
        #[arg(long)]
        json: bool,

        /// Minimum count to report
        #[arg(long, default_value_t = 1)]
        min_count: u64,

        /// Canonicalise k-mers (min of forward / reverse-complement)
        #[arg(long, default_value_t = true)]
        canonical: bool,

        /// Number of threads (0 = all)
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },

    /// Detect tandem repeats by looking for periodic k-mer runs
    Tandem {
        /// Input FASTA assembly
        #[arg(short, long)]
        input: PathBuf,

        /// k-mer length used for periodicity detection
        #[arg(short, long, default_value_t = 11)]
        k: u8,

        /// Minimum number of tandem copies
        #[arg(long, default_value_t = 3)]
        min_copies: u32,

        /// Maximum period (in bp) to consider
        #[arg(long, default_value_t = 100)]
        max_period: u32,

        /// Output BED file of tandem-repeat intervals
        #[arg(short, long)]
        output: PathBuf,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },

    /// Detect repetitive regions from high-frequency k-mers
    Repeats {
        /// Input FASTA assembly
        #[arg(short, long)]
        input: PathBuf,

        /// k-mer length
        #[arg(short, long, default_value_t = 21)]
        k: u8,

        /// Minimum k-mer count to consider a position “repetitive”
        #[arg(long, default_value_t = 5)]
        min_count: u64,

        /// Minimum length (bp) of a repetitive region
        #[arg(long, default_value_t = 100)]
        min_len: u32,

        /// Merge gaps shorter than this many bases
        #[arg(long, default_value_t = 50)]
        merge_gap: u32,

        /// Output BED file of repetitive intervals
        #[arg(short, long)]
        output: PathBuf,

        /// Also write a TSV of per-position coverage
        #[arg(long)]
        coverage: Option<PathBuf>,

        /// Canonicalise k-mers
        #[arg(long, default_value_t = true)]
        canonical: bool,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
    },

    /// Run count + tandem + repeats in one pass and write a summary directory
    All {
        /// Input FASTA assembly
        #[arg(short, long)]
        input: PathBuf,

        /// k-mer length for counting / repeats
        #[arg(short, long, default_value_t = 21)]
        k: u8,

        /// k-mer length for tandem detection
        #[arg(long, default_value_t = 11)]
        tandem_k: u8,

        /// Output directory
        #[arg(short, long)]
        output: PathBuf,

        /// Minimum count for repeat calling
        #[arg(long, default_value_t = 5)]
        min_count: u64,

        /// Number of threads
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
            info!("Counting {}-mers from {:?}", k, input);
            let counts = kmer::count_kmers(&input, k, canonical)?;
            io::write_kmer_counts(&output, &counts, min_count, json)?;
            info!("Wrote {} distinct k-mers (min_count={}) to {:?}", counts.len(), min_count, output);
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
            info!("Detecting tandem repeats (k={}, min_copies={}, max_period={})", k, min_copies, max_period);
            let intervals = tandem::detect_tandems(&input, k, min_copies, max_period)?;
            io::write_bed(&output, &intervals)?;
            info!("Found {} tandem-repeat intervals → {:?}", intervals.len(), output);
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
            info!("Calling repetitive regions (k={}, min_count={}, min_len={})", k, min_count, min_len);
            let (intervals, cov) = repeats::detect_repeats(
                &input, k, min_count, min_len, merge_gap, canonical,
            )?;
            io::write_bed(&output, &intervals)?;
            if let Some(cov_path) = coverage {
                io::write_coverage(&cov_path, &cov)?;
            }
            info!("Found {} repetitive intervals → {:?}", intervals.len(), output);
        }

        Commands::All {
            input,
            k,
            tandem_k,
            output,
            min_count,
            threads,
        } => {
            set_threads(threads);
            std::fs::create_dir_all(&output)?;
            info!("Running full analysis pipeline on {:?}", input);

            // 1. k-mer counts
            let counts = kmer::count_kmers(&input, k, true)?;
            let count_path = output.join("kmers.tsv");
            io::write_kmer_counts(&count_path, &counts, 1, false)?;
            info!("K-mer counts → {:?}", count_path);

            // 2. tandem repeats
            let tandems = tandem::detect_tandems(&input, tandem_k, 3, 100)?;
            let tandem_path = output.join("tandems.bed");
            io::write_bed(&tandem_path, &tandems)?;
            info!("Tandem repeats → {:?}", tandem_path);

            // 3. repetitive regions
            let (reps, _cov) = repeats::detect_repeats(&input, k, min_count, 100, 50, true)?;
            let rep_path = output.join("repeats.bed");
            io::write_bed(&rep_path, &reps)?;
            info!("Repetitive regions → {:?}", rep_path);

            // summary
            let summary = serde_json::json!({
                "k": k,
                "tandem_k": tandem_k,
                "distinct_kmers": counts.len(),
                "tandem_intervals": tandems.len(),
                "repetitive_intervals": reps.len(),
            });
            let summary_path = output.join("summary.json");
            std::fs::write(&summary_path, serde_json::to_string_pretty(&summary)?)?;
            info!("Summary → {:?}", summary_path);
        }
    }

    Ok(())
}

fn set_threads(n: usize) {
    let n = if n == 0 { num_cpus::get() } else { n };
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .ok(); // ignore if already initialised
    info!("Using {} threads", n);
}

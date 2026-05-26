mod counter;
mod histogram;

use clap::{Parser, Subcommand};
use std::fs::File;
use std::io::{BufWriter, Write};
use rusty_k::validate_k_size;

#[derive(Parser)]
#[command(name = "rusty-k")]
#[command(about = "A fast k-mer counter for DNA sequences", long_about = None)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Count k-mers in DNA sequences (FASTA/FASTQ format)
    Count {
        /// Path to input FASTA or FASTQ file
        #[arg(value_name = "FILE")]
        input: String,

        /// K-mer size
        #[arg(short, long, value_name = "SIZE")]
        kmer_size: usize,

        /// Output file (default: stdout)
        #[arg(short, long, value_name = "FILE")]
        output: Option<String>,

        /// Number of threads (default: 1)
        #[arg(short, long, value_name = "NUM", default_value = "1")]
        threads: usize,
    },

    /// Generate histogram of k-mer frequencies
    Histogram {
        /// Path to k-mer counts file
        #[arg(value_name = "FILE")]
        input: String,

        /// Output file (default: stdout)
        #[arg(short, long, value_name = "FILE")]
        output: Option<String>,

        /// Number of threads (default: 1)
        #[arg(short, long, value_name = "NUM", default_value = "1")]
        threads: usize,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Count {
            input,
            kmer_size,
            output,
            threads,
        } => {
            validate_k_size(kmer_size)?;
            let counts = counter::count_kmers_from_file_threaded(&input, kmer_size, threads)?;

            let mut sorted: Vec<_> = counts.into_iter().collect();
            sorted.sort_by_key(|a| a.0.clone());

            if let Some(output_path) = output {
                let file = File::create(&output_path)
                    .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", output_path, e))?;
                let mut writer = BufWriter::new(file);
                for (kmer, count) in sorted {
                    writeln!(writer, "{}\t{}", kmer, count)
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                }
            } else {
                for (kmer, count) in sorted {
                    println!("{}\t{}", kmer, count);
                }
            }
        }

        Commands::Histogram { input, output, threads } => {
            let counts = histogram::read_kmer_counts_threaded(&input, threads)?;
            let hist = histogram::create_histogram(&counts);

            let mut sorted: Vec<_> = hist.into_iter().collect();
            sorted.sort_by_key(|a| a.0);

            if let Some(output_path) = output {
                let file = File::create(&output_path)
                    .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", output_path, e))?;
                let mut writer = BufWriter::new(file);
                for (count, frequency) in sorted {
                    writeln!(writer, "{}\t{}", count, frequency)
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                }
            } else {
                for (count, frequency) in sorted {
                    println!("{}\t{}", count, frequency);
                }
            }
        }
    }

    Ok(())
}


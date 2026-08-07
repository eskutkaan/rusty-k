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

    /// Find repeat candidates in genome assemblies
    #[command(name = "repeatcandidates", alias = "repeats")]
    RepeatCandidates {
        /// Path to input genome FASTA file
        #[arg(value_name = "FILE")]
        input: String,

        /// K-mer size
        #[arg(short, long, value_name = "SIZE")]
        kmer_size: usize,

        /// Minimum count to treat a k-mer as repetitive
        #[arg(short = 'm', long = "min-count", value_name = "COUNT", default_value = "2")]
        min_count: u64,

        /// Output file (default: stdout)
        #[arg(short, long, value_name = "FILE")]
        output: Option<String>,

        /// Number of threads (default: 1)
        #[arg(short, long, value_name = "NUM", default_value = "1")]
        threads: usize,
    },

    /// Find tandem repeats in genome assemblies
    Tandem {
        /// Path to input genome FASTA file
        #[arg(value_name = "FILE")]
        input: String,

        /// Minimum motif size to test
        #[arg(long, value_name = "SIZE", default_value = "4")]
        min_period: usize,

        /// Maximum motif size to test
        #[arg(long, value_name = "SIZE", default_value = "50")]
        max_period: usize,

        /// Minimum number of adjacent motif copies required
        #[arg(short = 'c', long = "min-copies", value_name = "COUNT", default_value = "3")]
        min_copies: usize,

        /// Minimum total repeat length required
        #[arg(long = "min-length", value_name = "BP", default_value = "12")]
        min_length: usize,

        /// Maximum mismatches allowed per repeat copy
        #[arg(long = "max-mismatches", value_name = "COUNT", default_value = "1")]
        max_mismatches: usize,

        /// Emit the full per-call tandem repeat list instead of grouped locus summaries
        #[arg(long)]
        raw: bool,

        /// Output file (default: stdout)
        #[arg(short, long, value_name = "FILE")]
        output: Option<String>,

        /// Number of threads (default: 1)
        #[arg(short, long, value_name = "NUM", default_value = "1")]
        threads: usize,
    },

    /// Benchmark repeat finding across a range of k-mer sizes
    Benchmark {
        /// Path to input genome FASTA file
        #[arg(value_name = "FILE")]
        input: String,

        /// Minimum k-mer size to test
        #[arg(long, value_name = "SIZE")]
        min_k: usize,

        /// Maximum k-mer size to test
        #[arg(long, value_name = "SIZE")]
        max_k: usize,

        /// Step between tested k-mer sizes
        #[arg(long, value_name = "SIZE", default_value = "1")]
        step: usize,

        /// Target repeat fraction used to recommend a k-mer size
        #[arg(long = "target-repeat-fraction", value_name = "FRACTION", default_value = "0.05")]
        target_repeat_fraction: f64,

        /// Minimum count to treat a k-mer as repetitive
        #[arg(short = 'm', long = "min-count", value_name = "COUNT", default_value = "2")]
        min_count: u64,

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

        Commands::RepeatCandidates {
            input,
            kmer_size,
            min_count,
            output,
            threads,
        } => {
            validate_k_size(kmer_size)?;
            let repeats = counter::find_repeat_regions_from_fasta_threaded(
                &input,
                kmer_size,
                threads,
                min_count,
            )?;

            let mut sorted: Vec<_> = repeats.into_iter().collect();
            sorted.sort_by(|a, b| a.contig.cmp(&b.contig).then_with(|| a.start.cmp(&b.start)));

            if let Some(output_path) = output {
                let file = File::create(&output_path)
                    .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", output_path, e))?;
                let mut writer = BufWriter::new(file);
                for region in sorted {
                    writeln!(writer, "{}\t{}\t{}\t{}", region.contig, region.start, region.end, region.supporting_hits)
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                }
            } else {
                for region in sorted {
                    println!("{}\t{}\t{}\t{}", region.contig, region.start, region.end, region.supporting_hits);
                }
            }
        }

        Commands::Tandem {
            input,
            min_period,
            max_period,
            min_copies,
            min_length,
            max_mismatches,
            raw,
            output,
            threads,
        } => {
            if min_period == 0 {
                anyhow::bail!("min_period must be greater than 0");
            }
            if max_period < min_period {
                anyhow::bail!("max_period must be greater than or equal to min_period");
            }
            if min_copies < 2 {
                anyhow::bail!("min_copies must be at least 2");
            }
            if min_length == 0 {
                anyhow::bail!("min_length must be greater than 0");
            }
            if max_mismatches > min_period {
                anyhow::bail!("max_mismatches should not exceed min_period");
            }

            let repeats = counter::find_tandem_repeats_from_fasta_threaded(
                &input,
                min_period,
                max_period,
                min_copies,
                min_length,
                max_mismatches,
                threads,
            )?;

            if raw {
                let mut rows = repeats;
                rows.sort_by(|left, right| {
                    left.contig
                        .cmp(&right.contig)
                        .then_with(|| left.start.cmp(&right.start))
                        .then_with(|| left.end.cmp(&right.end))
                        .then_with(|| left.period.cmp(&right.period))
                });

                if let Some(output_path) = output {
                    let file = File::create(&output_path)
                        .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", output_path, e))?;
                    let mut writer = BufWriter::new(file);
                    writeln!(writer, "contig\tstart\tend\tperiod\tcopies\tmismatches\tmotif\tconsensus")
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                    for repeat in rows {
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            repeat.contig,
                            repeat.start,
                            repeat.end,
                            repeat.period,
                            repeat.copies,
                            repeat.mismatches,
                            repeat.motif,
                            repeat.consensus
                        )
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                    }
                } else {
                    println!("contig\tstart\tend\tperiod\tcopies\tmismatches\tmotif\tconsensus");
                    for repeat in rows {
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            repeat.contig,
                            repeat.start,
                            repeat.end,
                            repeat.period,
                            repeat.copies,
                            repeat.mismatches,
                            repeat.motif,
                            repeat.consensus
                        );
                    }
                }
            } else {
                let mut rows = counter::summarize_tandem_repeats(repeats);
                rows.sort_by(|left, right| {
                    left.contig
                        .cmp(&right.contig)
                        .then_with(|| left.start.cmp(&right.start))
                        .then_with(|| left.end.cmp(&right.end))
                        .then_with(|| left.period.cmp(&right.period))
                });

                if let Some(output_path) = output {
                    let file = File::create(&output_path)
                        .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", output_path, e))?;
                    let mut writer = BufWriter::new(file);
                    writeln!(writer, "contig\tstart\tend\tperiod\tcopies\tmismatches\tsupporting_calls\tconsensus")
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                    for repeat in rows {
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            repeat.contig,
                            repeat.start,
                            repeat.end,
                            repeat.period,
                            repeat.copies,
                            repeat.mismatches,
                            repeat.supporting_calls,
                            repeat.consensus
                        )
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                    }
                } else {
                    println!("contig\tstart\tend\tperiod\tcopies\tmismatches\tsupporting_calls\tconsensus");
                    for repeat in rows {
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            repeat.contig,
                            repeat.start,
                            repeat.end,
                            repeat.period,
                            repeat.copies,
                            repeat.mismatches,
                            repeat.supporting_calls,
                            repeat.consensus
                        );
                    }
                }
            }
        }

        Commands::Benchmark {
            input,
            min_k,
            max_k,
            step,
            target_repeat_fraction,
            min_count,
            output,
            threads,
        } => {
            validate_k_size(min_k)?;
            validate_k_size(max_k)?;
            let rows = counter::benchmark_repeat_k_sizes_from_fasta(
                &input,
                min_k,
                max_k,
                step,
                threads,
                min_count,
            )?;

            let best_row = counter::recommend_benchmark_row(&rows, target_repeat_fraction).cloned();

            if let Some(output_path) = output {
                let file = File::create(&output_path)
                    .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", output_path, e))?;
                let mut writer = BufWriter::new(file);
                writeln!(writer, "kmer_size\tregion_count\ttotal_repeat_bp\trepeat_fraction\trepeat_bp_per_mb\tregion_density_per_mb\taverage_region_len\tsupporting_hits")
                    .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                for row in &rows {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{:.6}\t{:.3}\t{:.3}\t{:.3}\t{}",
                        row.kmer_size,
                        row.region_count,
                        row.total_repeat_bp,
                        row.repeat_fraction,
                        row.repeat_bp_per_mb,
                        row.region_density_per_mb,
                        row.average_region_len,
                        row.supporting_hits
                    )
                    .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                }
                if let Some(best_row) = best_row {
                    writeln!(writer, "# target_repeat_fraction\t{:.3}", target_repeat_fraction)
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                    writeln!(writer, "# recommended_k\t{}", best_row.kmer_size)
                        .map_err(|e| anyhow::anyhow!("Failed to write to file: {}", e))?;
                }
            } else {
                println!("kmer_size\tregion_count\ttotal_repeat_bp\trepeat_fraction\trepeat_bp_per_mb\tregion_density_per_mb\taverage_region_len\tsupporting_hits");
                for row in &rows {
                    println!(
                        "{}\t{}\t{}\t{:.6}\t{:.3}\t{:.3}\t{:.3}\t{}",
                        row.kmer_size,
                        row.region_count,
                        row.total_repeat_bp,
                        row.repeat_fraction,
                        row.repeat_bp_per_mb,
                        row.region_density_per_mb,
                        row.average_region_len,
                        row.supporting_hits
                    );
                }
                if let Some(best_row) = best_row {
                    println!("target_repeat_fraction\t{:.3}", target_repeat_fraction);
                    println!("recommended_k\t{}", best_row.kmer_size);
                }
            }
        }
    }

    Ok(())
}


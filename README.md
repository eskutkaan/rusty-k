# rusty-k

A fast, parallel k-mer counter for DNA sequences written in Rust.

## Features

- Counts canonical k-mers in FASTA/FASTQ files
- Multi-threaded processing for large datasets
- Generates frequency histograms
- Proper error handling and validation
- Unified CLI with subcommands

## Installation

```bash
cargo build --release
```

The binary will be at `target/release/rusty-k`.

## Usage

### Count K-mers

Count canonical k-mers in a DNA sequence file:

```bash
rusty-k count <INPUT> --kmer-size <K> [-o OUTPUT] [-t THREADS]
```

**Arguments:**
- `INPUT`: Path to FASTA or FASTQ file

**Options:**
- `-k, --kmer-size <SIZE>`: K-mer size (1-32, required)
- `-o, --output <FILE>`: Output file (default: stdout)
- `-t, --threads <NUM>`: Number of threads (default: 1)

**Example:**
```bash
rusty-k count input.fasta --kmer-size 21 -o counts.txt -t 4
```

### Generate Histogram

Create a frequency histogram from k-mer counts:

```bash
rusty-k histogram <INPUT> [-o OUTPUT] [-t THREADS]
```

**Arguments:**
- `INPUT`: Path to k-mer counts file (output from `count` command)

**Options:**
- `-o, --output <FILE>`: Output file (default: stdout)
- `-t, --threads <NUM>`: Number of threads (default: 1)

**Example:**
```bash
rusty-k histogram counts.txt -o histogram.txt -t 4
```

## Complete Workflow

```bash
# Count k-mers
rusty-k count input.fasta --kmer-size 21 -o counts.txt -t 4

# Generate histogram
rusty-k histogram counts.txt -o histogram.txt
```

## Help

```bash
rusty-k --help
rusty-k count --help
rusty-k histogram --help
```


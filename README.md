# rusty-k

`rusty-k` is a fast Rust tool for k-mer counting and repeat analysis on genome assemblies.

## Install

```bash
cargo build --release
```

Binary: `target/release/rusty-k`

## Commands

### Count

```bash
rusty-k count input.fasta --kmer-size 21 -o counts.tsv -t 4
```

Counts canonical k-mers from FASTA or FASTQ.

### Histogram

```bash
rusty-k histogram counts.tsv -o histogram.tsv -t 4
```

Builds a frequency histogram from k-mer counts.

### Repeat Candidates

```bash
rusty-k repeatcandidates assembly.fasta --kmer-size 21 -m 2 -o repeats.bed -t 4
```

Reports repeat candidates from an assembly in BED-like format.

### Tandem

```bash
rusty-k tandem assembly.fasta --min-period 4 --max-period 50 -c 3 --min-length 12 --max-mismatches 1 -o tandem.tsv -t 4
```

Detects approximate tandem repeats. Default output is grouped by locus. Use `--raw` for call-level output.

Grouped output columns: `contig`, `start`, `end`, `period`, `copies`, `mismatches`, `supporting_calls`, `consensus`

Raw output columns: `contig`, `start`, `end`, `period`, `copies`, `mismatches`, `motif`, `consensus`

### Benchmark

```bash
rusty-k benchmark assembly.fasta --min-k 15 --max-k 31 --step 2 --target-repeat-fraction 0.05 -m 2 -o benchmark.tsv -t 4
```

Sweeps k-mer sizes and recommends the first k whose repeat fraction falls at or below the target fraction.

## Help

```bash
rusty-k --help
rusty-k tandem --help
```


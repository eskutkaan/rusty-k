# rusty-k

A fast, parallel Rust CLI for:

* **k-mer counting** (canonical, 2-bit packed, k ≤ 32)
* **tandem-repeat detection** (exact periods 1–N bp)
* **repetitive-region calling** from high-abundance k-mers

Designed for genome assemblies, including large contig and scaffold sets.

## Getting started

Build the optimized binary before running analyses:

```bash
cargo build --release
./target/release/rusty-k --help
```

Use `target/release/rusty-k` for production runs and performance comparisons.

## Large genomes

For very large genomes, use a count threshold and memory budget:

```bash
rusty-k repeats -i salamander.fa -o repeats.bed \
    --k 31 --min-count 5 --max-memory-mb 4096

rusty-k count -i salamander.fa -o frequent-kmers.tsv \
    --k 31 --min-count 5 --max-memory-mb 4096

rusty-k all -i salamander.fa -o results/ \
    --k 31 --kmer-min-count 5 --min-count 5 --max-memory-mb 4096 \
    --skip-tandem
```

`--max-memory-mb` sets the approximate memory budget for temporary counting
shards. Temporary files can be large, so ensure the system temporary directory
has sufficient free space. `--threads` controls parallel counting and tandem
detection. Use `--skip-tandem` for a fast first pass on chromosome-scale data.

## Sub-commands

### `count` – k-mer frequencies

```bash
rusty-k count -k 21 -i assembly.fa -o kmers.tsv
rusty-k count -k 21 -i assembly.fa -o kmers.json --json --min-count 2
```

Useful options:

```text
--k <N>                 k-mer length, 1–32
--min-count <N>         report only k-mers at or above this count
--max-memory-mb <N>     memory budget for thresholded counting
--canonical <true|false> count reverse complements together or separately
--threads <N>           number of worker threads; 0 uses available CPUs
```

### `tandem` – tandem repeats

```bash
rusty-k tandem -k 11 -i assembly.fa --min-copies 3 --max-period 50 -o tandems.bed
```

For `tandem`, `-k` is the minimum total repeat span in bases.

Output is BED6: `chrom  start  end  name  score  strand`  
where `name` encodes the period and copy number (e.g. `TR_period6_x12`).

### `repeats` – repetitive regions

```bash
rusty-k repeats -k 21 -i assembly.fa --min-count 5 --min-len 100 --merge-gap 50 \
                  -o repeats.bed --coverage coverage.tsv
```

A position is marked repetitive when an overlapping k-mer has genome-wide count
at least `min_count`. Neighboring high-coverage stretches separated by at most
`merge_gap` bases are merged.

### `all` – run everything

```bash
rusty-k all -k 21 --tandem-k 11 -i assembly.fa -o results/
```

The `all` command accepts the same key thresholds used by the individual
analyses, so a pipeline can be tuned without running separate commands:

```bash
rusty-k all -k 21 --tandem-k 11 --tandem-min-copies 4 --max-period 50 \
    --min-count 8 --min-len 150 --merge-gap 25 --coverage results/coverage.tsv \
    -i assembly.fa -o results/
```

Produces:

```
results/
├── kmers.tsv
├── tandems.bed
├── repeats.bed
└── summary.json
```

When `--coverage` is provided, the requested TSV is written in addition to
these files. `summary.json` records the selected parameters and input-file
provenance, including size and SHA-256 checksum.

## Output formats

- `kmers.tsv`: `kmer` and `count` columns, sorted by count.
- `tandems.bed`: BED6 intervals; names include repeat period and copy count.
- `repeats.bed`: BED6 repetitive intervals; coordinates are zero-based and
    half-open.
- `coverage.tsv`: optional non-zero per-position repeat coverage.
- `summary.json`: run parameters, input provenance, and result counts.

## Test data

*E. coli* K-12 MG1655 complete genome (NC_000913.3, 4.64 Mb) is included under `test/`:

```bash
./rusty-k all -k 21 --tandem-k 11 \
    -i test/ecoli_MG1655.fna -o test/results
```

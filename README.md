# rusty-k

A fast Rust CLI for counting DNA/RNA k-mers in FASTA and FASTQ files.

## Build

```bash
cargo build --release
```

## Usage

```bash
./target/release/rusty-k \
    --input assembly.fa \
    --output kmers.tsv \
    --k 31 \
    --threads 0
```

Temporary shard files are written to `tmp` beside the executable by default.
Use `--tmp-dir PATH` to choose another location:

```bash
./target/release/rusty-k \
    --input assembly.fa \
    --output kmers.tsv \
    --tmp-dir /scratch/rusty-k
```

Output follows shard order by default. Pass `--sort` to sort the final k-mer
table lexicographically by k-mer. Sorting uses bounded external runs in the
temporary directory, so it does not load the complete result table into RAM.

By default, reverse complements are counted together. Use `--canonical false`
to count each strand separately. `--threads 0` uses all available CPUs; set a
positive value to choose a worker count.

Counting is disk-backed for every run, so RAM does not grow with the number of
distinct k-mers. Plain FASTA is fed in bounded sequence chunks; FASTQ and
compressed inputs use the parser’s record buffering. `--min-count N` filters
results while reducing each shard. Adjust the approximate per-shard memory
budget with `--max-memory-mb`. Temporary shards require free disk space and are
removed when counting finishes. Results are written to a temporary file and
renamed into place only after counting completes successfully.

Output is TSV with `kmer` and `count` columns. Pass `--json` for a JSON array
instead. Sorting by count is not performed; `--sort` sorts by k-mer sequence.

K-mer lengths from 1 through 32 are supported. Ambiguous bases split the input
sequence and are not included in any k-mer.

## Test

```bash
cargo test
```

## Benchmark

Build first, then measure an end-to-end run. `/usr/bin/time` reports wall time
and peak resident memory on macOS:

```bash
/usr/bin/time -l ./target/release/rusty-k \
    --input test/ecoli_MG1655.fna \
    --output /tmp/kmers.tsv \
    --k 31 \
    --threads 0 \
    --max-memory-mb 1024
```

For meaningful comparisons, use the same input, thread count, memory budget,
canonical setting, storage device, and output destination. Report elapsed time,
peak memory, temporary-disk usage, and distinct k-mers written.

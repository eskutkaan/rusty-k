# rusty-k

A fast, parallel Rust CLI for:

* **k-mer counting** (canonical, 2-bit packed, k ≤ 32)
* **tandem-repeat detection** (exact periods 1–N bp)
* **repetitive-region calling** from high-abundance k-mers

Designed for genome *assemblies* (contigs / scaffolds). 

## Build

```bash
cargo build --release
# binary → target/release/rusty-k
```

## Sub-commands

### `count` – k-mer frequencies

```bash
rusty-k count -k 21 -i assembly.fa -o kmers.tsv
rusty-k count -k 21 -i assembly.fa -o kmers.json --json --min-count 2
```

### `tandem` – tandem repeats

```bash
rusty-k tandem -k 11 -i assembly.fa --min-copies 3 --max-period 50 -o tandems.bed
```

Output is BED6: `chrom  start  end  name  score  strand`  
where `name` encodes the period and copy number (e.g. `TR_period6_x12`).

### `repeats` – repetitive regions

```bash
rusty-k repeats -k 21 -i assembly.fa --min-count 5 --min-len 100 --merge-gap 50 \
                  -o repeats.bed --coverage coverage.tsv
```

A position is marked repetitive when the k-mer that starts there (or any overlapping k-mer) has genome-wide count ≥ `min_count`.  Neighbouring high-coverage stretches separated by ≤ `merge_gap` bases are merged.

### `all` – run everything

```bash
rusty-k all -k 21 --tandem-k 11 -i assembly.fa -o results/
```

Produces:

```
results/
├── kmers.tsv
├── tandems.bed
├── repeats.bed
└── summary.json
```

## Algorithm notes

| Component | Method |
|-----------|--------|
| Encoding | 2-bit DNA, 64-bit integers, k ≤ 32 |
| Canonical | `min(forward, reverse-complement)` |
| Counting | Parallel per-contig hash maps (FxHash) → merge |
| Tandem | Exact unit matching for every period 1…max_period; run-length collapse |
| Repeats | Global k-mer counts → per-base coverage track → interval collapsing |

## Test data

*E. coli* K-12 MG1655 complete genome (NC_000913.3, 4.64 Mb) is included under `test/`:

```bash
./rusty-k all -k 21 --tandem-k 11 \
    -i test/ecoli_MG1655.fna -o test/results
```

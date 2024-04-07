# rusty-k
a work-in-progress k-mer counter written in Rust

Usage
rusty-k <k-mer size> <input file> [OPTIONS]

Options

    -o, --output: Specify the output file to write the k-mer counts. If not provided, the output will be printed to the console.
    -t, --threads: Specify the number of threads to use for parallel processing. Default is 1.

Arguments

    <k-mer size>: The size of the k-mers to count (e.g., 21).
    <input file>: The path to the input FASTA or FASTQ file containing the DNA sequences.

Example

rusty-k 21 input.fasta/fastq -o output.txt -t 4

This command will count the canonical k-mers of size 21 in the input.fasta/fastq file using 4 threads and write the results to output.txt.

Histogram

The histogram script takes the output of the k-mer counting program and generates a histogram of k-mer frequencies, sorted numerically by the k-mer count.

Usage
bin/histogram [OPTIONS] <input file>

Options

    -o, --output: Specify the output file to write the histogram. If not provided, the histogram will be printed to the console.
    -t, --threads: Specify the number of threads to use for parallel processing. Default is 1.

Arguments

    <input file>: The path to the input file containing the k-mer counts (output from the main program).

Example

bin/histogram input.txt -o histogram.txt -t 4

This command will read the k-mer counts from input.txt using 4 threads, create a histogram, sort it numerically, and write the sorted histogram to histogram.txt.


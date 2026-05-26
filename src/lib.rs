use std::io::Write;

pub fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => base,
    }
}

pub fn revcomp(sequence: &str) -> String {
    let mut rev_comp = String::with_capacity(sequence.len());
    for c in sequence.chars().rev() {
        rev_comp.push(complement(c));
    }
    rev_comp
}

pub fn canonical_kmer(kmer: &str) -> String {
    let rev_comp = revcomp(kmer);
    if kmer < rev_comp.as_str() {
        kmer.to_string()
    } else {
        rev_comp
    }
}

pub fn write_output<W: Write>(
    writer: &mut W,
    items: &[(String, u64)],
    format_fn: fn(&mut W, &str, u64) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    for (key, count) in items {
        format_fn(writer, key, *count)?;
    }
    Ok(())
}

pub fn validate_k_size(k: usize) -> anyhow::Result<()> {
    if !(1..=32).contains(&k) {
        anyhow::bail!("k-mer size must be between 1 and 32, got {}", k);
    }
    Ok(())
}

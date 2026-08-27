use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid k-mer length: {0} (must be 1–32)")]
    InvalidK(u8),

    #[error("Sequence too short for k={0}")]
    SeqTooShort(u8),

    #[error("Parse error: {0}")]
    Parse(String),

    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),

    #[error("{0}")]
    Other(String),
}

use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid k-mer length: {0} (must be 1–32)")]
    InvalidK(u8),

    #[error("{0}")]
    Other(String),
}

//! Output dispatch for SAM, BAM, FASTA, and FASTQ formats.
//!
//! Manages the output pipeline including optional gzip/bzip2 compression
//! and multi-threaded BGZF for BAM output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use bzip2::write::BzEncoder;
use flate2::write::GzEncoder;

/// Default buffer size for the output writer (256 KB).
const OUTPUT_BUF_SIZE: usize = 256 * 1024;

/// Compression mode for output.
pub enum CompressionMode {
    None,
    Gzip,
    Bzip2,
}

/// Abstraction over output destinations (stdout or file, with optional compression).
pub struct OutputWriter {
    inner: BufWriter<Box<dyn Write>>,
}

impl OutputWriter {
    /// Create a writer to stdout with the given compression.
    pub fn stdout_with_compression(mode: CompressionMode) -> Self {
        let raw: Box<dyn Write> = match mode {
            CompressionMode::None => Box::new(io::stdout().lock()),
            CompressionMode::Gzip => {
                Box::new(GzEncoder::new(io::stdout().lock(), flate2::Compression::default()))
            }
            CompressionMode::Bzip2 => {
                Box::new(BzEncoder::new(io::stdout().lock(), bzip2::Compression::default()))
            }
        };
        Self { inner: BufWriter::with_capacity(OUTPUT_BUF_SIZE, raw) }
    }

    /// Create a writer to a file with the given compression.
    pub fn from_path_with_compression(path: &Path, mode: CompressionMode) -> Result<Self> {
        let file = File::create(path)?;
        let raw: Box<dyn Write> = match mode {
            CompressionMode::None => Box::new(file),
            CompressionMode::Gzip => Box::new(GzEncoder::new(file, flate2::Compression::default())),
            CompressionMode::Bzip2 => Box::new(BzEncoder::new(file, bzip2::Compression::default())),
        };
        Ok(Self { inner: BufWriter::with_capacity(OUTPUT_BUF_SIZE, raw) })
    }

    /// Write a complete SAM header (already formatted with newlines).
    pub fn write_header(&mut self, header: &str) -> Result<()> {
        self.inner.write_all(header.as_bytes())?;
        Ok(())
    }

    /// Write pre-formatted bytes to the output.
    ///
    /// The bytes may contain one or more SAM lines, or any other pre-formatted
    /// content. No framing or newlines are added — callers are responsible for
    /// including record delimiters.
    pub fn write_bytes(&mut self, data: &[u8]) -> Result<()> {
        self.inner.write_all(data)?;
        Ok(())
    }

    /// Flush and finalize the output.
    pub fn finish(mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}

impl io::Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

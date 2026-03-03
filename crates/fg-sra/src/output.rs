//! Output dispatch for SAM, BAM, FASTA, and FASTQ formats.
//!
//! Manages the output pipeline including optional gzip/bzip2 compression
//! and BGZF for BAM output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use bzip2::write::BzEncoder;
use flate2::write::GzEncoder;

/// Default buffer size for the text output writer (256 KB).
const OUTPUT_BUF_SIZE: usize = 256 * 1024;

/// Compression mode for text output.
pub enum CompressionMode {
    None,
    Gzip,
    Bzip2,
}

/// Abstraction over output destinations (stdout or file, with optional compression).
///
/// For SAM/FASTA/FASTQ text output, uses a `BufWriter` with optional gzip/bzip2
/// compression. For BAM output, uses a BGZF writer that handles block compression.
pub struct OutputWriter {
    inner: WriterInner,
}

/// Internal writer variant: text (buffered) or BAM (BGZF-compressed).
enum WriterInner {
    Text(BufWriter<Box<dyn Write>>),
    Bgzf(noodles_bgzf::io::Writer<Box<dyn Write>>),
}

impl OutputWriter {
    /// Create a text writer to stdout with the given compression.
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
        Self { inner: WriterInner::Text(BufWriter::with_capacity(OUTPUT_BUF_SIZE, raw)) }
    }

    /// Create a text writer to a file with the given compression.
    pub fn from_path_with_compression(path: &Path, mode: CompressionMode) -> Result<Self> {
        let file = File::create(path)?;
        let raw: Box<dyn Write> = match mode {
            CompressionMode::None => Box::new(file),
            CompressionMode::Gzip => Box::new(GzEncoder::new(file, flate2::Compression::default())),
            CompressionMode::Bzip2 => Box::new(BzEncoder::new(file, bzip2::Compression::default())),
        };
        Ok(Self { inner: WriterInner::Text(BufWriter::with_capacity(OUTPUT_BUF_SIZE, raw)) })
    }

    /// Create a BAM writer to stdout (BGZF-compressed).
    pub fn bam_stdout() -> Self {
        let raw: Box<dyn Write> = Box::new(io::stdout().lock());
        Self { inner: WriterInner::Bgzf(noodles_bgzf::io::Writer::new(raw)) }
    }

    /// Create a BAM writer to a file (BGZF-compressed).
    pub fn bam_from_path(path: &Path) -> Result<Self> {
        let file = File::create(path)?;
        let raw: Box<dyn Write> = Box::new(file);
        Ok(Self { inner: WriterInner::Bgzf(noodles_bgzf::io::Writer::new(raw)) })
    }

    /// Write a SAM header (text mode) or BAM header (BAM mode).
    ///
    /// For text mode, writes the header string verbatim.
    /// For BAM mode, encodes the BAM header: magic bytes, SAM header text,
    /// and reference sequence dictionary parsed from `@SQ` lines.
    pub fn write_header(&mut self, header: &str) -> Result<()> {
        match &mut self.inner {
            WriterInner::Text(w) => {
                w.write_all(header.as_bytes())?;
            }
            WriterInner::Bgzf(w) => {
                write_bam_header(w, header)?;
            }
        }
        Ok(())
    }

    /// Write pre-formatted bytes to the output.
    pub fn write_bytes(&mut self, data: &[u8]) -> Result<()> {
        match &mut self.inner {
            WriterInner::Text(w) => w.write_all(data)?,
            WriterInner::Bgzf(w) => w.write_all(data)?,
        }
        Ok(())
    }

    /// Flush and finalize the output.
    ///
    /// For BAM mode, writes the BGZF EOF marker.
    pub fn finish(self) -> Result<()> {
        match self.inner {
            WriterInner::Text(mut w) => {
                w.flush()?;
            }
            WriterInner::Bgzf(w) => {
                w.finish().context("failed to finalize BAM output")?;
            }
        }
        Ok(())
    }
}

impl io::Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            WriterInner::Text(w) => w.write(buf),
            WriterInner::Bgzf(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            WriterInner::Text(w) => w.flush(),
            WriterInner::Bgzf(w) => w.flush(),
        }
    }
}

/// Write the BAM header section through a BGZF writer.
///
/// Encodes: magic (`BAM\1`), SAM header text, and reference sequence
/// dictionary (from `@SQ` lines in the header text).
fn write_bam_header(writer: &mut impl Write, header_text: &str) -> Result<()> {
    // BAM magic.
    writer.write_all(b"BAM\x01")?;

    // SAM header text.
    let text_bytes = header_text.as_bytes();
    writer.write_all(&(text_bytes.len() as i32).to_le_bytes())?;
    writer.write_all(text_bytes)?;

    // Parse @SQ lines for the reference dictionary.
    let refs = parse_sq_lines(header_text);

    writer.write_all(&(refs.len() as i32).to_le_bytes())?;
    for (name, length) in &refs {
        let name_bytes = name.as_bytes();
        writer.write_all(&((name_bytes.len() + 1) as i32).to_le_bytes())?;
        writer.write_all(name_bytes)?;
        writer.write_all(&[0])?; // null terminator
        writer.write_all(&length.to_le_bytes())?;
    }

    Ok(())
}

/// Parse `@SQ` lines from a SAM header, returning `(name, length)` pairs.
fn parse_sq_lines(header: &str) -> Vec<(&str, i32)> {
    let mut refs = Vec::new();
    for line in header.lines() {
        if !line.starts_with("@SQ\t") {
            continue;
        }
        let mut name = "";
        let mut length: i32 = 0;
        for field in line.split('\t').skip(1) {
            if let Some(val) = field.strip_prefix("SN:") {
                name = val;
            } else if let Some(val) = field.strip_prefix("LN:") {
                length = val.parse().unwrap_or(0);
            }
        }
        if !name.is_empty() {
            refs.push((name, length));
        }
    }
    refs
}

/// Build a map from reference sequence name to BAM reference ID (0-based index).
///
/// Parses `@SQ` lines from the SAM header to establish the mapping.
pub fn build_ref_name_to_id(header: &str) -> std::collections::HashMap<String, i32> {
    parse_sq_lines(header)
        .into_iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_owned(), i as i32))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sq_lines() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n";
        let refs = parse_sq_lines(header);
        assert_eq!(refs, vec![("chr1", 248956422), ("chr2", 242193529)]);
    }

    #[test]
    fn test_parse_sq_lines_empty() {
        let header = "@HD\tVN:1.6\n@CO\tsome comment\n";
        let refs = parse_sq_lines(header);
        assert!(refs.is_empty());
    }

    #[test]
    fn test_build_ref_name_to_id() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n";
        let map = build_ref_name_to_id(header);
        assert_eq!(map.get("chr1"), Some(&0));
        assert_eq!(map.get("chr2"), Some(&1));
        assert_eq!(map.get("chr3"), None);
    }

    #[test]
    fn test_write_bam_header() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
        let mut buf = Vec::new();
        write_bam_header(&mut buf, header).unwrap();

        // Check magic.
        assert_eq!(&buf[0..4], b"BAM\x01");

        // Check header text length.
        let text_len = i32::from_le_bytes(buf[4..8].try_into().unwrap());
        assert_eq!(text_len, header.len() as i32);

        // Check header text.
        let text_end = 8 + text_len as usize;
        assert_eq!(&buf[8..text_end], header.as_bytes());

        // Check number of references.
        let n_ref = i32::from_le_bytes(buf[text_end..text_end + 4].try_into().unwrap());
        assert_eq!(n_ref, 1);

        // Check first reference name.
        let name_len = i32::from_le_bytes(buf[text_end + 4..text_end + 8].try_into().unwrap());
        assert_eq!(name_len, 5); // "chr1\0"
        let name_end = text_end + 8 + name_len as usize;
        assert_eq!(&buf[text_end + 8..name_end - 1], b"chr1");
        assert_eq!(buf[name_end - 1], 0); // null terminator

        // Check reference length.
        let ref_len = i32::from_le_bytes(buf[name_end..name_end + 4].try_into().unwrap());
        assert_eq!(ref_len, 100);
    }
}

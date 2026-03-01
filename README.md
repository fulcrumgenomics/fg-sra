[![Build](https://github.com/fulcrumgenomics/fg-sratosam/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/fg-sratosam/actions/workflows/ci.yml)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fg-sratosam/blob/main/LICENSE)

# fg-sratosam

High-performance SRA-to-SAM/BAM converter, replacing NCBI's `sam-dump` with
multi-threaded processing for significantly higher throughput.

<p>
<a href="https://fulcrumgenomics.com"><img src="https://raw.githubusercontent.com/fulcrumgenomics/fg-sratosam/main/.github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Overview

fg-sratosam converts NCBI SRA archives to SAM or BAM format. It uses FFI
bindings to the NCBI VDB C library (`libncbi-vdb`) for reading SRA data and
processes references in parallel for high throughput.

Key features:
- **Multi-threaded** reference processing with ordered output
- **SAM and BAM** output (BAM via multi-threaded BGZF compression)
- **gzip/bzip2** compression for SAM output
- **Region filtering** by genomic coordinates
- **RNA splice detection** with XS tag generation
- **MD tag computation**
- **Quality quantization**
- **Mate cache** for proper SAM flag and mate-pair information

## Installation

### Building from source

```bash
git clone --recurse-submodules https://github.com/fulcrumgenomics/fg-sratosam
cd fg-sratosam
cargo build --release
```

#### Prerequisites

- Rust (stable toolchain)
- CMake (for building the vendored ncbi-vdb C library)

The vendored ncbi-vdb library is built automatically during `cargo build`.

#### Pre-built VDB

To use a pre-installed VDB library instead of building from source, set:

```bash
export VDB_INCDIR=/path/to/ncbi-vdb/interfaces
export VDB_LIBDIR=/path/to/lib/containing/libncbi-vdb.a
cargo build --release
```

## Usage

```bash
# Convert an SRA accession to SAM
fg-sratosam SRR390728

# Convert to BAM
fg-sratosam --output-format bam --output-file output.bam SRR390728

# Primary alignments only, with unaligned reads
fg-sratosam -1 -u SRR390728

# Filter by region
fg-sratosam --aligned-region chr1:1000000-2000000 SRR390728

# Multi-threaded with 8 threads
fg-sratosam -t 8 SRR390728
```

For full usage, run:

```bash
fg-sratosam --help
```

## Workspace Structure

```
fg-sratosam/
├── crates/
│   ├── fg-sra-vdb-sys/    # Raw FFI bindings to libncbi-vdb
│   ├── fg-sra-vdb/        # Safe Rust wrappers over VDB
│   └── fg-sratosam/       # Binary crate (the converter)
└── vendor/
    └── ncbi-vdb/           # Vendored VDB library (git submodule)
```

## Resources

- [Issues](https://github.com/fulcrumgenomics/fg-sratosam/issues): Report a bug or request a feature
- [Pull requests](https://github.com/fulcrumgenomics/fg-sratosam/pulls): Submit a patch or new feature
- [Contributors guide](https://github.com/fulcrumgenomics/fg-sratosam/blob/main/CONTRIBUTING.md)
- [License](https://github.com/fulcrumgenomics/fg-sratosam/blob/main/LICENSE): Released under the MIT license

## Authors

- [Nils Homer](https://github.com/nh13)

## Sponsors

Development of fg-sratosam is supported by [Fulcrum Genomics](https://www.fulcrumgenomics.com).

[Become a sponsor](https://github.com/sponsors/fulcrumgenomics)

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](https://github.com/fulcrumgenomics/fg-sratosam/blob/main/LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/fg-sratosam/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/fg-sratosam/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.

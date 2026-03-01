//! Integration tests for VDB wrappers.
//!
//! These tests require a working VDB installation and network access to
//! resolve SRA accessions. Run with: `cargo test -- --ignored`

use fg_sra_vdb::manager::VdbManager;

/// A small test accession with aligned reads.
const TEST_ACCESSION: &str = "SRR390728";

#[test]
#[ignore = "requires VDB and network access to resolve SRA accessions"]
fn test_manager_make_read() {
    let mgr = VdbManager::make_read().expect("failed to create VDB manager");
    mgr.disable_pagemap_thread()
        .expect("failed to disable pagemap thread");
}

#[test]
#[ignore = "requires VDB and network access to resolve SRA accessions"]
fn test_open_database() {
    let mgr = VdbManager::make_read().unwrap();
    let db = mgr
        .open_db_read(TEST_ACCESSION)
        .expect("failed to open test accession");
    let tables = db.list_tables().expect("failed to list tables");
    assert!(
        tables
            .iter()
            .any(|t| t == "PRIMARY_ALIGNMENT" || t == "SEQUENCE"),
        "expected PRIMARY_ALIGNMENT or SEQUENCE table, found: {tables:?}"
    );
}

#[test]
#[ignore = "requires VDB and network access to resolve SRA accessions"]
fn test_open_table_and_cursor() {
    let mgr = VdbManager::make_read().unwrap();
    let db = mgr.open_db_read(TEST_ACCESSION).unwrap();
    let tbl = db
        .open_table_read("SEQUENCE")
        .expect("failed to open SEQUENCE table");
    let cursor = tbl.create_cursor_read().expect("failed to create cursor");
    let col_idx = cursor
        .add_column("(ascii)READ")
        .expect("failed to add READ column");
    cursor.open().expect("failed to open cursor");
    let (first, count) = cursor.id_range(col_idx).expect("failed to get id range");
    assert!(count > 0, "expected non-empty table, got count={count}");

    // Read the first row.
    let seq = cursor
        .read_str(first, col_idx)
        .expect("failed to read sequence");
    assert!(
        !seq.is_empty(),
        "expected non-empty sequence for row {first}"
    );
}

#[test]
#[ignore = "requires VDB and network access to resolve SRA accessions"]
fn test_reference_list() {
    let mgr = VdbManager::make_read().unwrap();
    let db = mgr.open_db_read(TEST_ACCESSION).unwrap();
    let reflist = fg_sra_vdb::reference::ReferenceList::make_database(&db, 0, 0)
        .expect("failed to create reference list");
    let count = reflist.count().expect("failed to count references");
    assert!(count > 0, "expected at least one reference");

    let ref_obj = reflist.get(0).expect("failed to get first reference");
    let name = ref_obj.name().expect("failed to get reference name");
    assert!(!name.is_empty(), "expected non-empty reference name");

    let seq_len = ref_obj
        .seq_length()
        .expect("failed to get reference length");
    assert!(seq_len > 0, "expected non-zero reference length");
}

#[test]
#[ignore = "requires VDB and network access to resolve SRA accessions"]
fn test_metadata_read() {
    let mgr = VdbManager::make_read().unwrap();
    let db = mgr.open_db_read(TEST_ACCESSION).unwrap();
    let meta = db.open_metadata_read().expect("failed to open metadata");

    // BAM_HEADER may or may not exist, but the metadata API should work.
    match meta.open_node_read("BAM_HEADER") {
        Ok(node) => {
            let header = node.read_all().expect("failed to read BAM_HEADER node");
            assert!(
                header.starts_with("@HD") || header.starts_with("@SQ"),
                "expected SAM header content"
            );
        }
        Err(_) => {
            // Not all accessions have BAM_HEADER — that's OK.
        }
    }
}

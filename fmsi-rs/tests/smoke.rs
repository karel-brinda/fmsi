use fmsi_rust::builder::construct;
use fmsi_rust::query::{query_record, QueryMode};

#[test]
fn roundtrip_export_small() {
    let ms = "ACGTac";
    let index = construct(ms, 3, true).unwrap();
    assert_eq!(index.export_ms(), ms);
}

#[test]
fn query_small_presence() {
    let ms = "ACGTac";
    let mut index = construct(ms, 3, false).unwrap();
    let got = query_record(&mut index, b"ACGT", 3, false, QueryMode::Or, false);
    assert_eq!(got.len(), 2);
}

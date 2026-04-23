# fmsi-rust

A faithful stable-core Rust port scaffold of FMSI.

What is included:
- `index`
- `query`
- `lookup`
- `export`
- MBWT/FMS logical layout preserved
- optional `klcp`
- `k <= 32`
- Rust-native serialization while keeping the familiar `.fmsi.*` component naming

What is intentionally not included:
- all experimental set-operation code
- file-format compatibility with SDSL
- production-grade suffix-array construction

## Important

The current `sa` module uses a naive suffix-array backend so the project stays self-contained.
For any serious benchmark, swap `src/sa.rs` to a `libsais` or `divsufsort` backend and keep the
rest of the crate unchanged.

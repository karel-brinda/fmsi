use crate::builder::{construct, infer_or_validate_k};
use crate::fasta::{read_masked_superstring, read_records};
use crate::index::FmsIndex;
use crate::query::{query_record, QueryMode};

pub fn run(args: &[String]) -> Result<i32, String> {
    if args.len() < 2 {
        print_usage();
        return Ok(1);
    }
    match args[1].as_str() {
        "index" => run_index(&args[2..]),
        "query" => run_query(&args[2..], false),
        "lookup" => run_query(&args[2..], true),
        "export" => run_export(&args[2..]),
        "-h" | "--help" => {
            print_usage();
            Ok(0)
        }
        "-v" | "--version" => {
            println!("0.1.0");
            Ok(0)
        }
        other => Err(format!("unrecognized command: {other}")),
    }
}

fn run_index(args: &[String]) -> Result<i32, String> {
    let mut k: Option<usize> = None;
    let mut no_streaming = false;
    let mut path: Option<String> = None;
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" => {
                print_usage_index();
                return Ok(0);
            }
            "-k" => {
                i += 1;
                k = Some(args.get(i).ok_or("missing value after -k")?.parse::<usize>().map_err(|e| e.to_string())?);
            }
            "-x" => no_streaming = true,
            s if !s.starts_with('-') => path = Some(s.to_string()),
            other => return Err(format!("unrecognized index option: {other}")),
        }
        i += 1;
    }
    let path = path.ok_or("path to masked superstring is required")?;
    let ms = read_masked_superstring(&path)?;
    let inferred = crate::dna::infer_k(&ms)?;
    let k = infer_or_validate_k(&ms, k)?;
    if k != inferred {
        eprintln!("WARNING: provided k ({k}) does not match the k inferred from the mask convention ({inferred}).");
    }
    let use_klcp = !no_streaming;
    let index = construct(&ms, k, use_klcp)?;
    index.dump_to_prefix(&path)?;
    Ok(0)
}

fn run_query(args: &[String], output_orders: bool) -> Result<i32, String> {
    let mut k: Option<usize> = None;
    let mut query_path = "-".to_string();
    let mut use_klcp = false;
    let mut mode = QueryMode::Or;
    let mut prefix: Option<String> = None;
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" => {
                if output_orders {
                    print_usage_lookup();
                } else {
                    print_usage_query();
                }
                return Ok(0);
            }
            "-q" => {
                i += 1;
                query_path = args.get(i).ok_or("missing value after -q")?.clone();
            }
            "-k" => {
                i += 1;
                k = Some(args.get(i).ok_or("missing value after -k")?.parse::<usize>().map_err(|e| e.to_string())?);
            }
            "-S" => use_klcp = true,
            "-O" if !output_orders => mode = QueryMode::All,
            s if !s.starts_with('-') => prefix = Some(s.to_string()),
            other => return Err(format!("unrecognized query option: {other}")),
        }
        i += 1;
    }
    let prefix = prefix.ok_or("index prefix is required")?;
    let mut index = FmsIndex::load_from_prefix(&prefix, use_klcp)?;
    if use_klcp && index.klcp.is_none() {
        return Err("kLCP array was not constructed for this index".to_string());
    }
    let k = match k {
        Some(value) if value != index.k => {
            return Err(format!("provided k ({value}) does not match the k of the index ({})", index.k))
        }
        Some(value) => value,
        None => index.k,
    };
    let records = read_records(&query_path)?;
    for record in records {
        let result = query_record(&mut index, &record.seq, k, use_klcp, mode, output_orders);
        println!("{}\t{}", record.name, result);
    }
    Ok(0)
}

fn run_export(args: &[String]) -> Result<i32, String> {
    let mut prefix: Option<String> = None;
    for arg in args {
        match arg.as_str() {
            "-h" => {
                print_usage_export();
                return Ok(0);
            }
            s if !s.starts_with('-') => prefix = Some(s.to_string()),
            other => return Err(format!("unrecognized export option: {other}")),
        }
    }
    let prefix = prefix.ok_or("index prefix is required")?;
    let index = FmsIndex::load_from_prefix(&prefix, false)?;
    println!(">exported masked superstring");
    println!("{}", index.export_ms());
    Ok(0)
}

fn print_usage() {
    eprintln!("Usage: fmsi <command> [options]");
    eprintln!("Commands:");
    eprintln!("  index   Create an MBWT/FMS index from a masked superstring");
    eprintln!("  query   Query k-mers against the index");
    eprintln!("  lookup  Return k-mer order identifiers for present k-mers");
    eprintln!("  export  Print the underlying masked superstring");
}

fn print_usage_index() {
    eprintln!("Usage: fmsi index [options] <masked-superstring-input>");
    eprintln!("  -k INT  size of k-mers [default: infer from mask trailing zeros - 1]");
    eprintln!("  -x      do not compute the kLCP array used for faster streaming queries");
}

fn print_usage_query() {
    eprintln!("Usage: fmsi query [options] <index-prefix>");
    eprintln!("  -q FILE path to FASTA/FASTQ with queries [default: stdin]");
    eprintln!("  -k INT  size of k-mers [default: infer from index]");
    eprintln!("  -S      use the kLCP array for streamed queries");
    eprintln!("  -O      use max-one masked-superstring semantics");
}

fn print_usage_lookup() {
    eprintln!("Usage: fmsi lookup [options] <index-prefix>");
    eprintln!("  -q FILE path to FASTA/FASTQ with queries [default: stdin]");
    eprintln!("  -k INT  size of k-mers [default: infer from index]");
    eprintln!("  -S      use the kLCP array for streamed queries");
}

fn print_usage_export() {
    eprintln!("Usage: fmsi export <index-prefix>");
}

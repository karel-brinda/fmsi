fn main() {
    let args: Vec<String> = std::env::args().collect();
    match fmsi_rust::cli::run(&args) {
        Ok(code) => std::process::exit(code),
        Err(err) => {
            eprintln!("ERROR: {err}");
            std::process::exit(1);
        }
    }
}

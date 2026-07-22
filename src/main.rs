use anyhow::{anyhow, bail, Context, Result};
use std::io::{self, Write};
use std::time::Instant;

mod binencode;
mod certificate;
mod checker;
mod numerics;
mod postprocess;

use binencode::write_certificate_bin;
use certificate::{certificate_to_string, Certificate};
use checker::{check_certificate, check_generated_certificate};
use postprocess::{
    build_full_dim_certificate, build_item_neighbors_and_lifts,
    build_simplex_graph, certificate_inequalities, choose_default_k0,
    parse_lrs_ext_items, parse_lrs_hrep, root_certificate,
};

#[derive(Debug, Clone)]
enum Command {
    Postprocess(PostprocessArgs),
    Check(CheckArgs),
}

#[derive(Debug, Clone)]
struct PostprocessArgs {
    ine_path: String,
    ext_path: String,
    pretty: bool,
    bin: bool,
    k0: Option<usize>,
}

#[derive(Debug, Clone)]
struct CheckArgs {
    ine_path: String,
    certificate_path: String,
}

fn print_help(program: &str) {
    eprintln!(
        "Usage:
  {program} postprocess input.ine output.ext [OPTIONS]
  {program} check input.ine certificate.json

Postprocess options:
  --pretty          Pretty-print JSON output
  --bin             Write binary output for coq-binreader
  --k0 <INDEX>      Root vertex item index, 0-based.
                    Defaults to a vertex with the smallest number of simplices.

General options:
  -h, --help        Show this help message

Examples:
  {program} postprocess cross3.ine cross3.ext --pretty
  {program} postprocess cross3.ine cross3.ext --bin > cross3-cert.bin
  {program} postprocess cross3.ine cross3.ext --pretty --k0 0
  {program} check cross3.ine certificate.json

Backward-compatible shorthand:
  {program} input.ine output.ext [OPTIONS]
"
    );
}

fn parse_usize_flag(value: Option<String>, flag: &str) -> Result<usize> {
    let value = value.ok_or_else(|| anyhow!("{flag} must be followed by a 0-based integer"))?;

    value
        .parse::<usize>()
        .with_context(|| format!("invalid value for {flag}: `{value}`"))
}

fn parse_postprocess_args<I>(
    program: &str,
    mut it: I,
    first_positional: Option<String>,
) -> Result<Command>
where
    I: Iterator<Item = String>,
{
    let mut positional: Vec<String> = Vec::new();

    if let Some(arg) = first_positional {
        positional.push(arg);
    }

    let mut pretty = false;
    let mut bin = false;
    let mut k0: Option<usize> = None;

    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_help(program);
                std::process::exit(0);
            }

            "--pretty" => {
                pretty = true;
            }

            "--bin" => {
                bin = true;
            }

            "--k0" => {
                k0 = Some(parse_usize_flag(it.next(), "--k0")?);
            }

            _ if arg.starts_with('-') => {
                bail!("unknown option `{arg}`; use --help for usage");
            }

            _ => {
                positional.push(arg);
            }
        }
    }

    if pretty && bin {
        bail!("options --pretty and --bin are incompatible");
    }

    if positional.len() != 2 {
        print_help(program);
        bail!(
            "postprocess expects exactly 2 positional arguments: input.ine and output.ext; got {}",
            positional.len()
        );
    }

    Ok(Command::Postprocess(PostprocessArgs {
        ine_path: positional[0].clone(),
        ext_path: positional[1].clone(),
        pretty,
        bin,
        k0,
    }))
}

fn parse_check_args<I>(program: &str, mut it: I) -> Result<Command>
where
    I: Iterator<Item = String>,
{
    let mut positional: Vec<String> = Vec::new();

    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_help(program);
                std::process::exit(0);
            }

            _ if arg.starts_with('-') => {
                bail!("unknown option `{arg}`; use --help for usage");
            }

            _ => {
                positional.push(arg);
            }
        }
    }

    if positional.len() != 2 {
        print_help(program);
        bail!(
            "check expects exactly 2 positional arguments: input.ine and certificate.json; got {}",
            positional.len()
        );
    }

    Ok(Command::Check(CheckArgs {
        ine_path: positional[0].clone(),
        certificate_path: positional[1].clone(),
    }))
}

fn parse_args() -> Result<Command> {
    let mut it = std::env::args();
    let program = it.next().unwrap_or_else(|| "lrs-postprocess".to_string());

    let Some(first) = it.next() else {
        print_help(&program);
        bail!("missing command or positional arguments");
    };

    match first.as_str() {
        "-h" | "--help" => {
            print_help(&program);
            std::process::exit(0);
        }

        "postprocess" => parse_postprocess_args(&program, it, None),

        "check" => parse_check_args(&program, it),

        // Backward-compatible shorthand:
        //   program input.ine output.ext [OPTIONS]
        _ => parse_postprocess_args(&program, it, Some(first)),
    }
}

fn run_postprocess(args: PostprocessArgs) -> Result<()> {
    // Report each construction phase immediately. This makes long-running
    // certificate generation observable while preserving an aggregate timer.
    let certificate_start = Instant::now();

    let start = Instant::now();
    let h = parse_lrs_hrep(&args.ine_path).context("failed to parse H-representation")?;
    let read_ine_elapsed = start.elapsed();
    eprintln!(
        "Certificate: read and parse .ine file: {:.6} s",
        read_ine_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let (items, labels) = parse_lrs_ext_items(&args.ext_path, h.d, h.a.len())
        .context("failed to parse lrs ext output")?;
    let read_ext_elapsed = start.elapsed();
    eprintln!(
        "Certificate: stream .ext and build vertices/facets: {:.6} s",
        read_ext_elapsed.as_secs_f64()
    );

    if labels.is_empty() {
        bail!("no vertex/cobasis records found in lrs ext output");
    }

    let start = Instant::now();
    let graph = build_simplex_graph(labels).context("failed to build simplex graph")?;
    let build_graph_elapsed = start.elapsed();
    eprintln!(
        "Certificate: build facet graph: {:.6} s",
        build_graph_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let k0 = match args.k0 {
        Some(k0) => {
            if k0 >= items.len() {
                bail!("--k0={k0} out of range; there are {} vertex items", items.len());
            }

            if !graph.lbl.iter().any(|label| label.owner == k0) {
                bail!("--k0={k0} owns no simplex");
            }

            k0
        }

        None => choose_default_k0(&items, &graph)?,
    };
    let choose_root_elapsed = start.elapsed();
    eprintln!(
        "Certificate: choose root vertex: {:.6} s",
        choose_root_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let root = root_certificate(&h, &items, &graph, k0)
        .context("failed to build root certificate")?;
    let root_elapsed = start.elapsed();
    eprintln!(
        "Certificate: build root certificate: {:.6} s",
        root_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let inequalities = certificate_inequalities(&h)?;
    let inequalities_elapsed = start.elapsed();
    eprintln!(
        "Certificate: convert inequalities: {:.6} s",
        inequalities_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let (neighbors, geom_edge_lifts) =
        build_item_neighbors_and_lifts(&graph, items.len())?;
    let geom_graph_elapsed = start.elapsed();
    eprintln!(
        "Certificate: build geometric graph and lifts: {:.6} s",
        geom_graph_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let root_owner = graph
        .lbl
        .get(root.simplex_id)
        .ok_or_else(|| {
            anyhow!(
                "root simplex id {} is out of bounds",
                root.simplex_id
            )
        })?
        .owner;
    let full_dim = build_full_dim_certificate(&items, &neighbors, root_owner, h.d)?;
    let full_dim_elapsed = start.elapsed();
    eprintln!(
        "Certificate: build full-dimensionality certificate: {:.6} s",
        full_dim_elapsed.as_secs_f64()
    );

    let start = Instant::now();
    let cert = Certificate {
        n_inequalities: h.a.len(),
        dimension: h.d,
        inequalities,
        items,
        graph,
        neighbors,
        geom_edge_lifts,
        full_dim,
        root,
    };
    let assemble_elapsed = start.elapsed();
    eprintln!(
        "Certificate: assemble value: {:.6} s",
        assemble_elapsed.as_secs_f64()
    );
    let certificate_elapsed = certificate_start.elapsed();
    eprintln!(
        "Create certificate in memory: {:.6} s",
        certificate_elapsed.as_secs_f64()
    );

    // Check the generated value before any serialization or file I/O. This is
    // the production path used by both JSON and binary certificate generation.
    let check_start = Instant::now();
    check_generated_certificate(&h, &cert)?;
    eprintln!(
        "Check certificate: {:.6} s",
        check_start.elapsed().as_secs_f64()
    );
    eprintln!("Generated certificate accepted");

    if args.bin {
        let stdout = io::stdout();
        let mut out = io::BufWriter::new(stdout.lock());

        let start = Instant::now();
        write_certificate_bin(&mut out, &cert)?;
        out.flush()?;
        eprintln!(
            "Generate and write binary certificate: {:.6} s",
            start.elapsed().as_secs_f64()
        );
    } else {
        let start = Instant::now();
        let text = certificate_to_string(&cert, args.pretty)?;
        eprintln!(
            "Serialize JSON certificate: {:.6} s",
            start.elapsed().as_secs_f64()
        );

        let stdout = io::stdout();
        let mut out = io::BufWriter::new(stdout.lock());
        let start = Instant::now();
        out.write_all(text.as_bytes())?;
        out.write_all(b"\n")?;
        out.flush()?;
        eprintln!(
            "Write JSON certificate: {:.6} s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

fn run_check(args: CheckArgs) -> Result<()> {
    check_certificate(&args.ine_path, &args.certificate_path)?;
    println!("certificate accepted");
    Ok(())
}

fn main() -> Result<()> {
    match parse_args()? {
        Command::Postprocess(args) => run_postprocess(args),
        Command::Check(args) => run_check(args),
    }
}

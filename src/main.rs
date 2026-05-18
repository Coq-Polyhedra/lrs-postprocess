use anyhow::{anyhow, bail, Context, Result};

mod certificate;
mod checker;
mod numerics;
mod postprocess;

use certificate::certificate_to_string;
use checker::check_certificate;
use postprocess::{
    build_items, choose_default_k0, parse_lrs_ext_records, parse_lrs_hrep, root_certificate,
    to_certificate,
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
  --k0 <INDEX>      Root vertex item index, 0-based.
                    Defaults to a vertex with the smallest number of simplices.

General options:
  -h, --help        Show this help message

Examples:
  {program} postprocess cross3.ine cross3.ext --pretty
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
    let h = parse_lrs_hrep(&args.ine_path).context("failed to parse H-representation")?;

    let records = parse_lrs_ext_records(&args.ext_path, h.d, h.a.len())
        .context("failed to parse lrs ext output")?;

    if records.is_empty() {
        bail!("no vertex/cobasis records found in lrs ext output");
    }

    let items = build_items(records).context("failed to build vertex items")?;

    let k0 = match args.k0 {
        Some(k0) => {
            if k0 >= items.len() {
                bail!("--k0={k0} out of range; there are {} vertex items", items.len());
            }

            if items[k0].simplices.is_empty() {
                bail!("--k0={k0} has no simplices");
            }

            k0
        }

        None => choose_default_k0(&items)?,
    };

    let root = Some(root_certificate(&h, &items, k0).context("failed to build root certificate")?);

    let cert = to_certificate(items, root);
    println!("{}", certificate_to_string(&cert, args.pretty)?);

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
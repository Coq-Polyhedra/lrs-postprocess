use anyhow::{anyhow, bail, Context, Result};

mod certificate;
mod numerics;
mod postprocess;

use certificate::certificate_to_string;
use postprocess::{
    build_items, choose_default_k0, parse_lrs_ext_records, parse_lrs_hrep, root_certificate,
    to_certificate,
};

#[derive(Debug, Clone)]
struct Args {
    ine_path: String,
    ext_path: String,
    pretty: bool,
    k0: Option<usize>,
}

fn print_help(program: &str) {
    eprintln!(
        "Usage:
  {program} input.ine output.ext [OPTIONS]

Options:
  --pretty          Pretty-print JSON output
  --k0 <INDEX>      Root vertex item index, 0-based.
                    Defaults to a vertex with the smallest number of simplices.
  -h, --help        Show this help message

Examples:
  {program} cross3.ine cross3.ext --pretty
  {program} cross3.ine cross3.ext --pretty --k0 0
"
    );
}

fn parse_args() -> Result<Args> {
    let mut it = std::env::args();
    let program = it.next().unwrap_or_else(|| "lrs-postprocess".to_string());

    let mut positional: Vec<String> = Vec::new();
    let mut pretty = false;
    let mut k0: Option<usize> = None;

    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_help(&program);
                std::process::exit(0);
            }

            "--pretty" => {
                pretty = true;
            }

            "--k0" => {
                let value = it
                    .next()
                    .ok_or_else(|| anyhow!("--k0 must be followed by a 0-based integer"))?;

                k0 = Some(
                    value
                        .parse::<usize>()
                        .with_context(|| format!("invalid value for --k0: `{value}`"))?,
                );
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
        print_help(&program);
        bail!(
            "expected exactly 2 positional arguments: input.ine and output.ext; got {}",
            positional.len()
        );
    }

    Ok(Args {
        ine_path: positional[0].clone(),
        ext_path: positional[1].clone(),
        pretty,
        k0,
    })
}

fn main() -> Result<()> {
    let args = parse_args()?;

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
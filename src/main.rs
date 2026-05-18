use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;

#[derive(Debug, Clone)]
struct LrsRecord {
    cobasis: Vec<usize>,
    incident: Vec<usize>,
    vertex: Vec<String>,
}

#[derive(Debug, Serialize)]
struct VertexItemOut {
    vertex: Vec<String>,
    incident: Vec<usize>,
    simplices: Vec<Vec<usize>>,
}

fn strip_star(tok: &str) -> &str {
    tok.trim_end_matches('*')
}

fn parse_lrs_index(tok: &str) -> Result<usize> {
    let t = strip_star(tok);
    let idx_1based: usize = t
        .parse()
        .with_context(|| format!("bad lrs index token `{tok}`"))?;
    if idx_1based == 0 {
        bail!("lrs indices are expected to be 1-based, got 0 in `{tok}`");
    }
    Ok(idx_1based - 1)
}

fn is_index_token(tok: &str) -> bool {
    let t = strip_star(tok);
    !t.is_empty() && t.chars().all(|c| c.is_ascii_digit())
}

fn is_output_vector_line(line: &str) -> bool {
    let toks: Vec<_> = line.split_whitespace().collect();
    if toks.is_empty() {
        return false;
    }
    toks[0] == "1" || toks[0] == "0"
}

/// Parse one lrs header line with `facets ... : ... I#...`.
///
/// For H-representation output with `printcobasis 1` and `incidence`:
/// - entries after `facets` and before `:` are the cobasis facets;
/// - entries after `:` and before `I#...`, `det=...`, etc. are additional incident/tight inequalities;
/// - incident set = cobasis union additional tight inequalities.
fn parse_header(line: &str) -> Result<Option<(Vec<usize>, Vec<usize>)>> {
    let toks: Vec<&str> = line.split_whitespace().collect();

    // We focus on H-representation vertex/ray output lines.
    if !line.starts_with("V#") || !toks.iter().any(|&t| t == "facets") {
        return Ok(None);
    }

    let facets_pos = toks
        .iter()
        .position(|&t| t == "facets")
        .expect("checked above");

    let mut cobasis = Vec::new();
    let mut extra_incident = Vec::new();

    let mut pos = facets_pos + 1;
    let mut after_colon = false;

    while pos < toks.len() {
        let tok = toks[pos];
        if tok == ":" {
            after_colon = true;
            pos += 1;
            continue;
        }

        if tok.starts_with("I#") || tok.starts_with("det=") || tok.starts_with("in_det=") {
            break;
        }

        // Some lrs lines can have `det= 8`, so stop at the `det=` token above;
        // standalone numbers before this point are indices.
        if is_index_token(tok) {
            let idx = parse_lrs_index(tok)?;
            if after_colon {
                extra_incident.push(idx);
            } else {
                cobasis.push(idx);
            }
            pos += 1;
        } else {
            // Unknown token: stop the facet list rather than guess.
            break;
        }
    }

    if cobasis.is_empty() {
        bail!("found `facets` line with empty cobasis: `{line}`");
    }

    let mut incident = cobasis.clone();
    incident.extend(extra_incident);
    incident.sort_unstable();
    incident.dedup();

    Ok(Some((cobasis, incident)))
}

/// Parse lrs output produced with options like:
/// allbases
/// incidence
/// printcobasis 1
///
/// We trust lrs. No rational parsing or feasibility check is done.
fn parse_lrs_records(path: &str) -> Result<Vec<LrsRecord>> {
    let text = fs::read_to_string(path)?;
    let lines: Vec<&str> = text.lines().collect();
    let mut records = Vec::new();
    let mut i = 0usize;

    while i < lines.len() {
        let line = lines[i].trim();
        if let Some((cobasis, incident)) = parse_header(line)? {
            // Find next output vector line. For vertex enumeration, keep only lines
            // beginning with homogeneous coordinate `1`; ignore rays beginning with `0`.
            i += 1;
            while i < lines.len() {
                let vline = lines[i].trim();
                if vline.is_empty()
                    || vline.starts_with('*')
                    || vline.eq_ignore_ascii_case("begin")
                    || vline.eq_ignore_ascii_case("end")
                {
                    i += 1;
                    continue;
                }

                if is_output_vector_line(vline) {
                    let vtoks: Vec<&str> = vline.split_whitespace().collect();
                    if vtoks[0] == "1" {
                        let vertex = vtoks[1..].iter().map(|s| s.to_string()).collect();
                        records.push(LrsRecord {
                            cobasis,
                            incident,
                            vertex,
                        });
                    }
                    break;
                }

                // If we hit another header before a vector, stop looking.
                if vline.starts_with("V#") || vline.starts_with("F#") {
                    i = i.saturating_sub(1);
                    break;
                }

                i += 1;
            }
        }
        i += 1;
    }

    Ok(records)
}

fn cobasis_to_local_simplex(cobasis: &[usize], incident: &[usize]) -> Result<Vec<usize>> {
    let pos: BTreeMap<usize, usize> = incident
        .iter()
        .copied()
        .enumerate()
        .map(|(local, global)| (global, local))
        .collect();

    let mut local = Vec::with_capacity(cobasis.len());
    for &j in cobasis {
        let p = pos
            .get(&j)
            .copied()
            .with_context(|| format!("cobasis row {j} is missing from incident list"))?;
        local.push(p);
    }
    local.sort_unstable();
    Ok(local)
}

fn build_items(records: Vec<LrsRecord>, mod2: bool) -> Result<Vec<VertexItemOut>> {
    // Group records by exact vertex coordinate strings.
    let mut by_vertex: BTreeMap<Vec<String>, Vec<LrsRecord>> = BTreeMap::new();
    for rec in records {
        by_vertex.entry(rec.vertex.clone()).or_default().push(rec);
    }

    let mut out = Vec::new();

    for (vertex, recs) in by_vertex {
        // Incident set for the vertex: union of all incident lists reported by lrs.
        let mut incident_set = BTreeSet::new();
        for rec in &recs {
            for &j in &rec.incident {
                incident_set.insert(j);
            }
        }
        let incident: Vec<usize> = incident_set.into_iter().collect();

        // Each lrs cobasis gives one simplex, expressed locally in the union incident list.
        let mut simplex_parity: BTreeMap<Vec<usize>, bool> = BTreeMap::new();
        let mut simplices = Vec::new();

        for rec in &recs {
            let local = cobasis_to_local_simplex(&rec.cobasis, &incident)?;
            if mod2 {
                let e = simplex_parity.entry(local).or_insert(false);
                *e = !*e;
            } else {
                simplices.push(local);
            }
        }

        if mod2 {
            simplices = simplex_parity
                .into_iter()
                .filter_map(|(s, parity)| if parity { Some(s) } else { None })
                .collect();
        }

        out.push(VertexItemOut {
            vertex,
            incident,
            simplices,
        });
    }

    Ok(out)
}

fn main() -> Result<()> {
    let mut path: Option<String> = None;
    let mut pretty = false;
    let mut mod2 = false;

    for arg in env::args().skip(1) {
        match arg.as_str() {
            "--pretty" => pretty = true,
            "--mod2" => mod2 = true,
            "-h" | "--help" => {
                print_help();
                return Ok(());
            }
            _ => {
                if path.is_some() {
                    bail!("unexpected extra argument `{arg}`");
                }
                path = Some(arg);
            }
        }
    }

    let Some(path) = path else {
        print_help();
        bail!("missing lrs output file");
    };

    let records = parse_lrs_records(&path)?;
    if records.is_empty() {
        bail!("no vertex records with `facets` and homogeneous leading coordinate 1 found");
    }

    let items = build_items(records, mod2)?;

    if pretty {
        println!("{}", serde_json::to_string_pretty(&items)?);
    } else {
        println!("{}", serde_json::to_string(&items)?);
    }

    Ok(())
}

fn print_help() {
    eprintln!("lrs-postprocess");
    eprintln!();
    eprintln!("Usage:");
    eprintln!("  lrs-postprocess <lrs-output.txt> [--pretty] [--mod2]");
    eprintln!();
    eprintln!("Expected lrs options:");
    eprintln!("  allbases");
    eprintln!("  incidence");
    eprintln!("  printcobasis 1");
    eprintln!();
    eprintln!("Semantics for H-representation output:");
    eprintln!("  V#... facets a b c : d e I#... ");
    eprintln!("  cobasis simplex = [a,b,c] converted from 1-based to 0-based");
    eprintln!("  incident list    = [a,b,c,d,e] converted from 1-based to 0-based");
    eprintln!("  starred indices, e.g. 9*, are stripped before conversion");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  JSON list of {{ vertex: Vec<String>, incident: Vec<usize>, simplices: Vec<Vec<usize>> }}");
    eprintln!("  simplex entries are local indices into the vertex's incident list");
}

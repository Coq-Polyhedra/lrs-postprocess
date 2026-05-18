use anyhow::{anyhow, bail, Context, Result};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;

type Q = BigRational;

#[derive(Debug, Clone)]
struct Args {
    ine_path: String,
    lrs_output_path: String,
    pretty: bool,
    k0: Option<usize>,
    s0: Option<usize>,
}

#[derive(Debug, Clone)]
struct HRep {
    /// Matrix rows for A x <= b.
    a: Vec<Vec<Q>>,
    b: Vec<Q>,
    d: usize,
}

#[derive(Debug, Clone)]
struct LrsRecord {
    /// Vertex coordinates as strings, without the leading homogeneous coordinate.
    vertex_tokens: Vec<String>,

    /// Cobasis facets before `:` in the lrs incidence/cobasis line.
    /// Converted to 0-based row indices.
    cobasis: Vec<usize>,

    /// Additional incident inequalities after `:`.
    /// Converted to 0-based row indices.
    additional_incident: Vec<usize>,
}

#[derive(Debug, Clone)]
struct VertexItem {
    vertex_tokens: Vec<String>,
    incident: Vec<usize>,
    simplices: Vec<Vec<usize>>,
}

#[derive(Debug, Serialize)]
struct Output {
    items: Vec<VertexItemOut>,
    special: Option<SpecialOut>,
}

#[derive(Debug, Serialize)]
struct VertexItemOut {
    vertex: Vec<String>,
    incident: Vec<usize>,
    simplices: Vec<Vec<usize>>,
}

#[derive(Debug, Serialize)]
struct SpecialOut {
    k0: usize,
    s0: usize,

    /// Global row indices of the special simplex.
    special_rows: Vec<usize>,

    /// c_star = sum of the special row normals.
    c_star: Vec<String>,

    /// Inverse of the special basis matrix, provided as row vectors.
    ///
    /// If special rows are j_1,...,j_d and R has columns A_{j_r}^T,
    /// then inverse_rows[a] dot A_{j_b} = delta_{ab}.
    inverse_rows: Vec<Vec<String>>,

    /// For every other simplex of vertex k0, an inverse-row index a such that
    /// -inverse_rows[a] separates c_star from that cone.
    same_label_separators: Vec<SameLabelSeparatorOut>,
}

#[derive(Debug, Serialize)]
struct SameLabelSeparatorOut {
    simplex: usize,
    inverse_row: usize,
}

fn print_help(program: &str) {
    eprintln!(
        "Usage:
  {program} input.ine output.ext [OPTIONS]

Options:
  --pretty          Pretty-print JSON output
  --k0 <INDEX>      Special vertex item index, 0-based
  --s0 <INDEX>      Special simplex index inside items[k0].simplices, 0-based
  -h, --help        Show this help message

Examples:
  {program} cross3.ine cross3.out --pretty
  {program} cross3.ine cross3.out --pretty --k0 0 --s0 0
"
    );
}

fn parse_args() -> Result<Args> {
    let mut it = std::env::args();
    let program = it.next().unwrap_or_else(|| "lrs-postprocess".to_string());

    let mut positional: Vec<String> = Vec::new();
    let mut pretty = false;
    let mut k0: Option<usize> = None;
    let mut s0: Option<usize> = None;

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

            "--s0" => {
                let value = it
                    .next()
                    .ok_or_else(|| anyhow!("--s0 must be followed by a 0-based integer"))?;
                s0 = Some(
                    value
                        .parse::<usize>()
                        .with_context(|| format!("invalid value for --s0: `{value}`"))?,
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

    match (k0, s0) {
        (Some(_), Some(_)) | (None, None) => {}
        (Some(_), None) => bail!("--k0 was provided but --s0 is missing"),
        (None, Some(_)) => bail!("--s0 was provided but --k0 is missing"),
    }

    Ok(Args {
        ine_path: positional[0].clone(),
        lrs_output_path: positional[1].clone(),
        pretty,
        k0,
        s0,
    })
}

fn parse_q(s: &str) -> Result<Q> {
    let s = s.trim();

    if let Some((num, den)) = s.split_once('/') {
        let n: BigInt = num
            .parse()
            .with_context(|| format!("bad numerator `{num}`"))?;
        let d: BigInt = den
            .parse()
            .with_context(|| format!("bad denominator `{den}`"))?;
        if d.is_zero() {
            bail!("zero denominator in rational `{s}`");
        }
        Ok(BigRational::new(n, d))
    } else {
        let n: BigInt = s
            .parse()
            .with_context(|| format!("bad integer rational `{s}`"))?;
        Ok(BigRational::from_integer(n))
    }
}

fn q_to_string(x: &Q) -> String {
    if x.denom().is_one() {
        x.numer().to_string()
    } else {
        format!("{}/{}", x.numer(), x.denom())
    }
}

fn dot(a: &[Q], x: &[Q]) -> Q {
    assert_eq!(a.len(), x.len());
    a.iter()
        .zip(x)
        .fold(Q::zero(), |acc, (ai, xi)| acc + ai * xi)
}

fn mat_mul(a: &[Vec<Q>], b: &[Vec<Q>]) -> Vec<Vec<Q>> {
    let n = a.len();
    let k = if n == 0 { 0 } else { a[0].len() };
    let m = if b.is_empty() { 0 } else { b[0].len() };

    assert_eq!(b.len(), k);

    let mut c = vec![vec![Q::zero(); m]; n];

    for i in 0..n {
        for r in 0..k {
            for j in 0..m {
                c[i][j] += &a[i][r] * &b[r][j];
            }
        }
    }

    c
}

fn identity(n: usize) -> Vec<Vec<Q>> {
    let mut id = vec![vec![Q::zero(); n]; n];
    for i in 0..n {
        id[i][i] = Q::one();
    }
    id
}

fn invert_matrix(a: &[Vec<Q>]) -> Result<Vec<Vec<Q>>> {
    let n = a.len();

    if n == 0 {
        bail!("cannot invert empty matrix");
    }

    if a.iter().any(|row| row.len() != n) {
        bail!("matrix is not square");
    }

    let mut aug = vec![vec![Q::zero(); 2 * n]; n];

    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j].clone();
        }
        aug[i][n + i] = Q::one();
    }

    for col in 0..n {
        let pivot = (col..n).find(|&r| !aug[r][col].is_zero());

        let Some(pivot) = pivot else {
            bail!("matrix is singular");
        };

        if pivot != col {
            aug.swap(pivot, col);
        }

        let pivot_val = aug[col][col].clone();

        for j in 0..2 * n {
            aug[col][j] /= pivot_val.clone();
        }

        for r in 0..n {
            if r == col {
                continue;
            }

            if aug[r][col].is_zero() {
                continue;
            }

            let factor = aug[r][col].clone();

            for j in 0..2 * n {
                let sub = factor.clone() * aug[col][j].clone();
                aug[r][j] -= sub;
            }
        }
    }

    let inv = aug
        .into_iter()
        .map(|row| row[n..].to_vec())
        .collect::<Vec<_>>();

    Ok(inv)
}

/// Parse an lrs-style matrix between begin/end.
fn parse_lrs_matrix(path: &str) -> Result<Vec<Vec<Q>>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('*'));

    while let Some(line) = lines.next() {
        if line.eq_ignore_ascii_case("begin") {
            let header = lines.next().context("missing matrix header after begin")?;
            let h: Vec<_> = header.split_whitespace().collect();

            if h.len() < 2 {
                bail!("bad matrix header `{header}`");
            }

            let rows: usize = h[0].parse()?;
            let cols: usize = h[1].parse()?;

            let mut mat = Vec::with_capacity(rows);

            for r in 0..rows {
                let row_line = lines.next().with_context(|| format!("missing matrix row {r}"))?;
                let toks: Vec<_> = row_line.split_whitespace().collect();

                if toks.len() != cols {
                    bail!("row {r}: expected {cols} columns, got {}", toks.len());
                }

                let row = toks
                    .iter()
                    .map(|t| parse_q(t))
                    .collect::<Result<Vec<_>>>()?;
                mat.push(row);
            }

            return Ok(mat);
        }
    }

    bail!("no begin/end matrix found in `{path}`")
}

/// lrs H-rep row: beta + alpha*x >= 0.
/// Convert to A x <= b by A = -alpha, b = beta.
fn parse_lrs_hrep(path: &str) -> Result<HRep> {
    let mat = parse_lrs_matrix(path)?;

    if mat.is_empty() {
        bail!("empty H-representation");
    }

    let cols = mat[0].len();

    if cols < 2 {
        bail!("H-representation must have at least 2 columns");
    }

    let d = cols - 1;

    let mut a = Vec::with_capacity(mat.len());
    let mut b = Vec::with_capacity(mat.len());

    for (r, row) in mat.into_iter().enumerate() {
        if row.len() != cols {
            bail!("inconsistent row length at H row {r}");
        }

        let beta = row[0].clone();
        let normal: Vec<Q> = row[1..].iter().map(|x| -x.clone()).collect();

        b.push(beta);
        a.push(normal);
    }

    Ok(HRep { a, b, d })
}

fn strip_star(tok: &str) -> &str {
    tok.trim_end_matches('*')
}

fn parse_lrs_index_1based(tok: &str) -> Result<usize> {
    let t = strip_star(tok);
    let value = t
        .parse::<usize>()
        .with_context(|| format!("bad lrs index `{tok}`"))?;

    if value == 0 {
        bail!("lrs indices are expected to be 1-based, got 0");
    }

    Ok(value - 1)
}

fn is_index_token(tok: &str) -> bool {
    let t = strip_star(tok);
    !t.is_empty() && t.chars().all(|c| c.is_ascii_digit())
}

fn is_vertex_or_ray_row(line: &str, d: usize) -> bool {
    let toks: Vec<&str> = line.split_whitespace().collect();

    if toks.len() != d + 1 {
        return false;
    }

    toks.first() == Some(&"1") || toks.first() == Some(&"0")
}

/// Parse lrs output produced with H-representation input and:
///
/// allbases
/// incidence
/// printcobasis 1
///
/// For a line like:
///
/// V#1 R#0 B#1 h=0 facets  12 14 15 16 : 9 10 11 13 I#8 det= 8
///  1  0  0  0  1
///
/// - 12 14 15 16 are the cobasis facets and define a simplex.
/// - 9 10 11 13 are additional incident inequalities.
/// - all are converted from 1-based to 0-based.
fn parse_lrs_output_records(path: &str, d: usize, m_ineq: usize) -> Result<Vec<LrsRecord>> {
    let text = fs::read_to_string(path)?;
    let lines: Vec<&str> = text.lines().collect();

    let mut records = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        if line.starts_with("V#") && line.contains(" facets ") {
            let toks: Vec<&str> = line.split_whitespace().collect();

            let facets_pos = toks
                .iter()
                .position(|&t| t == "facets")
                .ok_or_else(|| anyhow!("line has `facets` text but could not locate token: `{line}`"))?;

            let mut cobasis = Vec::new();
            let mut additional = Vec::new();

            let mut p = facets_pos + 1;
            let mut after_colon = false;

            while p < toks.len() {
                let tok = toks[p];

                if tok == ":" {
                    after_colon = true;
                    p += 1;
                    continue;
                }

                if tok.starts_with("I#")
                    || tok.starts_with("det=")
                    || tok.starts_with("in_det=")
                    || tok.starts_with("h=")
                    || tok.starts_with("V#")
                    || tok.starts_with("R#")
                    || tok.starts_with("B#")
                {
                    break;
                }

                if is_index_token(tok) {
                    let idx = parse_lrs_index_1based(tok)?;

                    if idx >= m_ineq {
                        bail!(
                            "facet index {} out of range 1..={} in line `{}`",
                            idx + 1,
                            m_ineq,
                            line
                        );
                    }

                    if after_colon {
                        additional.push(idx);
                    } else {
                        cobasis.push(idx);
                    }

                    p += 1;
                    continue;
                }

                break;
            }

            if cobasis.len() != d {
                bail!(
                    "cobasis has {} facets, expected d={d}, in line `{}`",
                    cobasis.len(),
                    line
                );
            }

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

                if is_vertex_or_ray_row(vline, d) {
                    let vtoks: Vec<&str> = vline.split_whitespace().collect();

                    // Keep only vertices. Ignore rays.
                    if vtoks[0] == "1" {
                        records.push(LrsRecord {
                            vertex_tokens: vtoks[1..].iter().map(|s| (*s).to_string()).collect(),
                            cobasis,
                            additional_incident: additional,
                        });
                    }

                    break;
                }

                i += 1;
            }
        }

        i += 1;
    }

    Ok(records)
}

fn sorted_unique(mut v: Vec<usize>) -> Vec<usize> {
    v.sort_unstable();
    v.dedup();
    v
}

fn build_items(records: Vec<LrsRecord>) -> Result<Vec<VertexItem>> {
    let mut vertex_to_id: BTreeMap<Vec<String>, usize> = BTreeMap::new();
    let mut items: Vec<VertexItem> = Vec::new();

    for rec in records {
        let vid = if let Some(&id) = vertex_to_id.get(&rec.vertex_tokens) {
            id
        } else {
            let id = items.len();
            vertex_to_id.insert(rec.vertex_tokens.clone(), id);
            items.push(VertexItem {
                vertex_tokens: rec.vertex_tokens.clone(),
                incident: Vec::new(),
                simplices: Vec::new(),
            });
            id
        };

        let item = &mut items[vid];

        let mut incident = item.incident.clone();
        incident.extend(rec.cobasis.iter().copied());
        incident.extend(rec.additional_incident.iter().copied());
        item.incident = sorted_unique(incident);

        let mut simplex_global = rec.cobasis;
        simplex_global.sort_unstable();

        let pos: HashMap<usize, usize> = item
            .incident
            .iter()
            .copied()
            .enumerate()
            .map(|(local, global)| (global, local))
            .collect();

        let simplex_local = simplex_global
            .iter()
            .map(|j| {
                pos.get(j)
                    .copied()
                    .ok_or_else(|| anyhow!("internal error: simplex row missing from incident set"))
            })
            .collect::<Result<Vec<_>>>()?;

        item.simplices.push(simplex_local);
    }

    Ok(items)
}

fn row_matrix_columns_from_global_rows(h: &HRep, rows: &[usize]) -> Vec<Vec<Q>> {
    let d = h.d;
    let mut m = vec![vec![Q::zero(); d]; d];

    for (col, &j) in rows.iter().enumerate() {
        for row in 0..d {
            m[row][col] = h.a[j][row].clone();
        }
    }

    m
}

fn decode_simplex_global(item: &VertexItem, simplex_index: usize) -> Result<Vec<usize>> {
    let simplex = item
        .simplices
        .get(simplex_index)
        .ok_or_else(|| anyhow!("invalid simplex index {simplex_index}"))?;

    simplex
        .iter()
        .map(|&local| {
            item.incident
                .get(local)
                .copied()
                .ok_or_else(|| anyhow!("local simplex index {local} out of range"))
        })
        .collect()
}

fn compute_c_star(h: &HRep, rows: &[usize]) -> Vec<Q> {
    let mut c = vec![Q::zero(); h.d];

    for &j in rows {
        for q in 0..h.d {
            c[q] += h.a[j][q].clone();
        }
    }

    c
}

fn special_certificate(h: &HRep, items: &[VertexItem], k0: usize, s0: usize) -> Result<SpecialOut> {
    let item = items
        .get(k0)
        .ok_or_else(|| anyhow!("invalid k0={k0}; only {} items", items.len()))?;

    let special_rows = decode_simplex_global(item, s0)
        .with_context(|| format!("cannot decode special simplex s0={s0} for k0={k0}"))?;

    if special_rows.len() != h.d {
        bail!(
            "special simplex has {} rows, expected d={}",
            special_rows.len(),
            h.d
        );
    }

    let rmat = row_matrix_columns_from_global_rows(h, &special_rows);
    let inv = invert_matrix(&rmat).context("special basis matrix is singular")?;

    let prod = mat_mul(&inv, &rmat);
    if prod != identity(h.d) {
        bail!("internal error: computed inverse does not multiply to identity");
    }

    let c_star = compute_c_star(h, &special_rows);

    let mut same_label_separators = Vec::new();

    for sidx in 0..item.simplices.len() {
        if sidx == s0 {
            continue;
        }

        let rows = decode_simplex_global(item, sidx)?;

        let mut found = None;

        for a in 0..h.d {
            let beta = &inv[a];

            let ok = rows.iter().all(|&j| dot(beta, &h.a[j]) <= Q::zero());

            if ok {
                found = Some(a);
                break;
            }
        }

        let Some(inverse_row) = found else {
            bail!(
                "could not find inverse-row separator for simplex {sidx} of vertex item {k0}"
            );
        };

        same_label_separators.push(SameLabelSeparatorOut {
            simplex: sidx,
            inverse_row,
        });
    }

    Ok(SpecialOut {
        k0,
        s0,
        special_rows,
        c_star: c_star.iter().map(q_to_string).collect(),
        inverse_rows: inv
            .iter()
            .map(|row| row.iter().map(q_to_string).collect())
            .collect(),
        same_label_separators,
    })
}

fn to_output(items: Vec<VertexItem>, special: Option<SpecialOut>) -> Output {
    Output {
        items: items
            .into_iter()
            .map(|item| VertexItemOut {
                vertex: item.vertex_tokens,
                incident: item.incident,
                simplices: item.simplices,
            })
            .collect(),
        special,
    }
}

fn main() -> Result<()> {
    let args = parse_args()?;

    let h = parse_lrs_hrep(&args.ine_path).context("failed to parse H-representation")?;

    let records = parse_lrs_output_records(&args.lrs_output_path, h.d, h.a.len())
        .context("failed to parse lrs output")?;

    if records.is_empty() {
        bail!("no vertex/cobasis records found in lrs output");
    }

    let items = build_items(records).context("failed to build vertex items")?;

    let special = match (args.k0, args.s0) {
        (Some(k0), Some(s0)) => Some(
            special_certificate(&h, &items, k0, s0)
                .context("failed to build special-direction certificate")?,
        ),
        _ => None,
    };

    let output = to_output(items, special);

    if args.pretty {
        println!("{}", serde_json::to_string_pretty(&output)?);
    } else {
        println!("{}", serde_json::to_string(&output)?);
    }

    Ok(())
}
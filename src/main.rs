use anyhow::{anyhow, bail, Context, Result};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero, Signed};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::env;
use std::fs;

type Q = BigRational;

#[derive(Debug, Clone)]
struct HRep {
    /// Rows A_i of A x <= b.
    a: Vec<Vec<Q>>,
    d: usize,
}

#[derive(Debug, Clone)]
struct LrsRecord {
    cobasis_0based: Vec<usize>,
    incident_0based: Vec<usize>,
    vertex: Vec<String>,
}

#[derive(Debug, Serialize)]
struct VertexItemOut {
    vertex: Vec<String>,
    /// Global 0-based inequality indices active at this vertex.
    incident: Vec<usize>,
    /// Simplices are local indices into `incident`.
    simplices: Vec<Vec<usize>>,
}

#[derive(Debug, Serialize)]
struct SameLabelSeparatorOut {
    /// Index in items[k0].simplices.
    simplex: usize,
    /// Row index of inverse_rows used as beta_a.
    /// The separator is - beta_a.
    inverse_row: usize,
}

#[derive(Debug, Serialize)]
struct SpecialCertificateOut {
    k0: usize,
    s0: usize,
    /// Global row indices of the special simplex.
    special_global_rows: Vec<usize>,
    /// c* = sum of the row normals in special_global_rows.
    c_star: Vec<String>,
    /// Rows of R^{-1}, where R has columns A_j^T for j in special_global_rows.
    /// Verification can be done by inner products <inverse_rows[a], A_j_b> = delta_ab.
    inverse_rows: Vec<Vec<String>>,
    /// For every other simplex at k0, an inverse row beta_a such that
    /// <beta_a, A_h> <= 0 for every generator h of that simplex.
    same_label_separators: Vec<SameLabelSeparatorOut>,
}

#[derive(Debug, Serialize)]
struct Output {
    dimension: usize,
    items: Vec<VertexItemOut>,
    special: SpecialCertificateOut,
}

#[derive(Debug, Clone)]
struct Config {
    ine_path: String,
    lrs_out_path: String,
    pretty: bool,
    mod2: bool,
    k0: Option<usize>,
    s0: Option<usize>,
}

fn parse_args() -> Result<Config> {
    let mut args = env::args().skip(1);
    let ine_path = args.next().ok_or_else(|| anyhow!("missing input .ine file"))?;
    let lrs_out_path = args.next().ok_or_else(|| anyhow!("missing lrs output file"))?;

    let mut cfg = Config {
        ine_path,
        lrs_out_path,
        pretty: false,
        mod2: false,
        k0: None,
        s0: None,
    };

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--pretty" => cfg.pretty = true,
            "--mod2" => cfg.mod2 = true,
            "--k0" => {
                let v = args.next().ok_or_else(|| anyhow!("--k0 requires an argument"))?;
                cfg.k0 = Some(v.parse()?);
            }
            "--s0" => {
                let v = args.next().ok_or_else(|| anyhow!("--s0 requires an argument"))?;
                cfg.s0 = Some(v.parse()?);
            }
            "-h" | "--help" => {
                print_help();
                std::process::exit(0);
            }
            _ => bail!("unknown argument `{arg}`"),
        }
    }

    Ok(cfg)
}

fn print_help() {
    eprintln!("Usage:");
    eprintln!("  lrs-postprocess input.ine lrs-output.txt [--pretty] [--mod2] [--k0 K --s0 S]");
    eprintln!();
    eprintln!("Input assumptions:");
    eprintln!("  - input.ine is an lrs H-representation with rows b + alpha*x >= 0.");
    eprintln!("  - lrs output was produced with allbases, incidence, printcobasis 1.");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  JSON with vertex triples (vertex, incident inequalities, simplices)");
    eprintln!("  and a compact special certificate (k0, s0, inverse rows, separator-row indices).");
}

fn parse_q(s: &str) -> Result<Q> {
    let s = s.trim();
    if let Some((num, den)) = s.split_once('/') {
        let n: BigInt = num.parse().with_context(|| format!("bad numerator `{num}`"))?;
        let d: BigInt = den.parse().with_context(|| format!("bad denominator `{den}`"))?;
        if d.is_zero() {
            bail!("zero denominator `{s}`");
        }
        Ok(BigRational::new(n, d))
    } else {
        let n: BigInt = s.parse().with_context(|| format!("bad integer `{s}`"))?;
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
                mat.push(toks.iter().map(|t| parse_q(t)).collect::<Result<Vec<_>>>()?);
            }
            return Ok(mat);
        }
    }
    bail!("no begin/end matrix found in `{path}`")
}

/// lrs H-row: beta + alpha*x >= 0. Convert to A x <= b via A = -alpha.
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
    for row in mat {
        if row.len() != cols {
            bail!("inconsistent H row length");
        }
        a.push(row[1..].iter().map(|x| -x.clone()).collect());
    }
    Ok(HRep { a, d })
}

fn strip_star(tok: &str) -> &str {
    tok.trim_end_matches('*')
}

fn parse_lrs_index(tok: &str) -> Option<usize> {
    let t = strip_star(tok);
    if t.is_empty() || !t.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    t.parse::<usize>().ok()
}

fn is_record_header(line: &str) -> bool {
    (line.starts_with("V#") || line.starts_with("F#")) && (line.contains(" facets ") || line.contains(" vertices/rays "))
}

/// Parse only H-representation vertex records:
/// V#... facets c1 c2 ... : extra1 extra2 ... I#... det=...
/// next output vector line begins with homogeneous coordinate 1 for vertices, 0 for rays.
fn parse_lrs_output_records(path: &str, d: usize, m_ineq: usize) -> Result<Vec<LrsRecord>> {
    let text = fs::read_to_string(path)?;
    let lines: Vec<&str> = text.lines().collect();
    let mut records = Vec::new();
    let mut i = 0usize;

    while i < lines.len() {
        let line = lines[i].trim();
        if line.starts_with("V#") && line.contains(" facets ") {
            let toks: Vec<&str> = line.split_whitespace().collect();
            let facets_pos = toks
                .iter()
                .position(|&t| t == "facets")
                .ok_or_else(|| anyhow!("no facets token in `{line}`"))?;

            let mut before_colon = Vec::new();
            let mut after_colon = Vec::new();
            let mut in_after = false;
            let mut p = facets_pos + 1;
            while p < toks.len() {
                let t = toks[p];
                if t == ":" {
                    in_after = true;
                    p += 1;
                    continue;
                }
                if t.starts_with("I#") || t.starts_with("det=") || t.starts_with("in_det=") || t == "det" || t == "in_det" {
                    break;
                }
                if let Some(idx1) = parse_lrs_index(t) {
                    if idx1 == 0 || idx1 > m_ineq {
                        bail!("lrs inequality index {idx1} out of range 1..={m_ineq} in `{line}`");
                    }
                    if in_after {
                        after_colon.push(idx1 - 1);
                    } else {
                        before_colon.push(idx1 - 1);
                    }
                } else {
                    // Stop at the first non-index token we do not recognize.
                    break;
                }
                p += 1;
            }

            if before_colon.len() != d {
                bail!(
                    "cobasis has {} facets but expected d={}: `{}`",
                    before_colon.len(),
                    d,
                    line
                );
            }

            let mut incident_set: BTreeSet<usize> = before_colon.iter().copied().collect();
            incident_set.extend(after_colon.iter().copied());
            let incident_0based: Vec<usize> = incident_set.into_iter().collect();

            // Find the following vertex/ray row. Keep only leading homogeneous coordinate 1.
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
                if is_record_header(vline) {
                    // Malformed/no vector row; resume from this header.
                    i -= 1;
                    break;
                }

                let vtoks: Vec<&str> = vline.split_whitespace().collect();
                if vtoks.len() == d + 1 {
                    if vtoks[0] == "1" {
                        records.push(LrsRecord {
                            cobasis_0based: sorted(before_colon),
                            incident_0based,
                            vertex: vtoks[1..].iter().map(|s| s.to_string()).collect(),
                        });
                    }
                    // If leading 0, this is a ray; ignore.
                    break;
                }
                i += 1;
            }
        }
        i += 1;
    }

    Ok(records)
}

fn sorted(mut v: Vec<usize>) -> Vec<usize> {
    v.sort_unstable();
    v
}

fn dot(a: &[Q], b: &[Q]) -> Q {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).fold(Q::zero(), |acc, (x, y)| acc + x * y)
}

fn invert_square(mut a: Vec<Vec<Q>>) -> Result<Vec<Vec<Q>>> {
    let n = a.len();
    if n == 0 || a.iter().any(|r| r.len() != n) {
        bail!("matrix to invert must be nonempty square");
    }

    let mut inv = vec![vec![Q::zero(); n]; n];
    for i in 0..n {
        inv[i][i] = Q::one();
    }

    for col in 0..n {
        let pivot = (col..n).find(|&r| !a[r][col].is_zero());
        let Some(piv) = pivot else {
            bail!("special basis matrix is singular");
        };
        if piv != col {
            a.swap(piv, col);
            inv.swap(piv, col);
        }

        let pv = a[col][col].clone();
        for j in 0..n {
            a[col][j] /= pv.clone();
            inv[col][j] /= pv.clone();
        }

        for r in 0..n {
            if r == col || a[r][col].is_zero() {
                continue;
            }
            let factor = a[r][col].clone();
            for j in 0..n {
                let av = factor.clone() * a[col][j].clone();
                let iv = factor.clone() * inv[col][j].clone();
                a[r][j] -= av;
                inv[r][j] -= iv;
            }
        }
    }

    Ok(inv)
}

fn matrix_from_rows_as_columns(h: &HRep, rows: &[usize]) -> Vec<Vec<Q>> {
    let d = h.d;
    let mut m = vec![vec![Q::zero(); d]; d];
    for (col, &row_idx) in rows.iter().enumerate() {
        for r in 0..d {
            m[r][col] = h.a[row_idx][r].clone();
        }
    }
    m
}

fn sum_normals(h: &HRep, rows: &[usize]) -> Vec<Q> {
    let mut c = vec![Q::zero(); h.d];
    for &j in rows {
        for r in 0..h.d {
            c[r] += h.a[j][r].clone();
        }
    }
    c
}

fn build_items(records: Vec<LrsRecord>, mod2: bool) -> Result<Vec<VertexItemOut>> {
    #[derive(Default)]
    struct Accum {
        incident: BTreeSet<usize>,
        bases: BTreeMap<Vec<usize>, bool>,
    }

    let mut vertex_order: Vec<Vec<String>> = Vec::new();
    let mut vertex_to_id: BTreeMap<Vec<String>, usize> = BTreeMap::new();
    let mut accums: Vec<Accum> = Vec::new();

    for rec in records {
        let vid = if let Some(&id) = vertex_to_id.get(&rec.vertex) {
            id
        } else {
            let id = vertex_order.len();
            vertex_to_id.insert(rec.vertex.clone(), id);
            vertex_order.push(rec.vertex.clone());
            accums.push(Accum::default());
            id
        };

        accums[vid].incident.extend(rec.incident_0based.iter().copied());
        accums[vid].incident.extend(rec.cobasis_0based.iter().copied());

        let basis = sorted(rec.cobasis_0based);
        if mod2 {
            let e = accums[vid].bases.entry(basis).or_insert(false);
            *e = !*e;
        } else {
            accums[vid].bases.insert(basis, true);
        }
    }

    let mut items = Vec::new();
    for (vid, acc) in accums.into_iter().enumerate() {
        let incident: Vec<usize> = acc.incident.into_iter().collect();
        let pos: HashMap<usize, usize> = incident.iter().copied().enumerate().map(|(i, r)| (r, i)).collect();

        let mut simplices = Vec::new();
        for (basis, keep) in acc.bases {
            if !keep {
                continue;
            }
            let mut local = Vec::with_capacity(basis.len());
            for j in basis {
                let p = *pos
                    .get(&j)
                    .ok_or_else(|| anyhow!("vertex {vid}: basis row {j} absent from incident set"))?;
                local.push(p);
            }
            simplices.push(local);
        }

        items.push(VertexItemOut {
            vertex: vertex_order[vid].clone(),
            incident,
            simplices,
        });
    }
    Ok(items)
}

fn decode_local_simplex(item: &VertexItemOut, s: usize) -> Result<Vec<usize>> {
    let loc = item
        .simplices
        .get(s)
        .ok_or_else(|| anyhow!("invalid simplex index {s}"))?;
    let mut rows = Vec::with_capacity(loc.len());
    for &i in loc {
        let j = *item
            .incident
            .get(i)
            .ok_or_else(|| anyhow!("local simplex index {i} out of bounds"))?;
        rows.push(j);
    }
    Ok(rows)
}

fn build_special(h: &HRep, items: &[VertexItemOut], k0_opt: Option<usize>, s0_opt: Option<usize>) -> Result<SpecialCertificateOut> {
    let k0 = k0_opt.unwrap_or_else(|| {
        items
            .iter()
            .position(|it| !it.simplices.is_empty())
            .unwrap_or(0)
    });
    if k0 >= items.len() {
        bail!("k0={k0} out of range");
    }
    if items[k0].simplices.is_empty() {
        bail!("chosen k0={k0} has no simplices");
    }
    let s0 = s0_opt.unwrap_or(0);
    if s0 >= items[k0].simplices.len() {
        bail!("s0={s0} out of range for k0={k0}");
    }

    let special_rows = decode_local_simplex(&items[k0], s0)?;
    if special_rows.len() != h.d {
        bail!("special simplex has size {}, expected d={}", special_rows.len(), h.d);
    }

    let rmat = matrix_from_rows_as_columns(h, &special_rows);
    let inv = invert_square(rmat).context("failed to invert special basis matrix")?;
    let cstar = sum_normals(h, &special_rows);

    let mut same_label_separators = Vec::new();
    for s in 0..items[k0].simplices.len() {
        if s == s0 {
            continue;
        }
        let rows = decode_local_simplex(&items[k0], s)?;
        let mut found = None;
        'row_loop: for (a, beta) in inv.iter().enumerate() {
            for &j in &rows {
                if dot(beta, &h.a[j]).is_positive() {
                    continue 'row_loop;
                }
            }
            found = Some(a);
            break;
        }
        let Some(a) = found else {
            bail!(
                "no inverse-row separator found for simplex {s} at k0={k0}; choose another special simplex or extend certificate format"
            );
        };
        same_label_separators.push(SameLabelSeparatorOut { simplex: s, inverse_row: a });
    }

    Ok(SpecialCertificateOut {
        k0,
        s0,
        special_global_rows: special_rows,
        c_star: cstar.iter().map(q_to_string).collect(),
        inverse_rows: inv
            .iter()
            .map(|r| r.iter().map(q_to_string).collect())
            .collect(),
        same_label_separators,
    })
}

fn main() -> Result<()> {
    let cfg = parse_args()?;
    let h = parse_lrs_hrep(&cfg.ine_path).context("failed to parse .ine H-representation")?;
    let records = parse_lrs_output_records(&cfg.lrs_out_path, h.d, h.a.len())
        .context("failed to parse lrs output")?;
    if records.is_empty() {
        bail!("no vertex records found in lrs output");
    }

    let items = build_items(records, cfg.mod2)?;
    let special = build_special(&h, &items, cfg.k0, cfg.s0)?;

    let out = Output {
        dimension: h.d,
        items,
        special,
    };

    if cfg.pretty {
        println!("{}", serde_json::to_string_pretty(&out)?);
    } else {
        println!("{}", serde_json::to_string(&out)?);
    }
    Ok(())
}

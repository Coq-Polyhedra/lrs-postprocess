use anyhow::{anyhow, bail, Context, Result};
use num_traits::Zero;
use std::collections::BTreeMap;
use std::fs;

use crate::certificate::{Certificate, Root, VertexItem};
use crate::numerics::{dot, identity, invert_matrix, mat_mul, parse_q, q_to_string, Q};

#[derive(Debug, Clone)]
pub struct HRep {
    /// Matrix rows for A x <= b.
    pub a: Vec<Vec<Q>>,
    pub b: Vec<Q>,
    pub d: usize,
}

#[derive(Debug, Clone)]
pub struct LrsRecord {
    /// Vertex coordinates as strings, without the leading homogeneous coordinate.
    pub vertex: Vec<String>,

    /// Cobasis facets before `:` in the lrs incidence/cobasis line.
    /// Converted to 0-based row indices.
    pub cobasis: Vec<usize>,

    /// Additional incident inequalities after `:`.
    /// Converted to 0-based row indices.
    pub additional_incident: Vec<usize>,
}

fn sorted_unique(mut v: Vec<usize>) -> Vec<usize> {
    v.sort_unstable();
    v.dedup();
    v
}

/// Parse an lrs-style matrix between begin/end.
pub fn parse_lrs_matrix(path: &str) -> Result<Vec<Vec<Q>>> {
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
pub fn parse_lrs_hrep(path: &str) -> Result<HRep> {
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
pub fn parse_lrs_ext_records(path: &str, d: usize, m_ineq: usize) -> Result<Vec<LrsRecord>> {
    let text = fs::read_to_string(path)?;
    let lines: Vec<&str> = text.lines().collect();

    let mut records = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        if !(line.starts_with("V#") && line.contains(" facets ")) {
            i += 1;
            continue;
        }

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

        // Pair this cobasis/incidence line with the next output row.
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

            if !is_vertex_or_ray_row(vline, d) {
                bail!(
                    "expected vertex/ray row immediately after cobasis line, got `{}` after `{}`",
                    vline,
                    line
                );
            }

            let vtoks: Vec<&str> = vline.split_whitespace().collect();

            // Keep only vertices. Ignore rays.
            if vtoks[0] == "1" {
                records.push(LrsRecord {
                    vertex: vtoks[1..].iter().map(|s| (*s).to_string()).collect(),
                    cobasis,
                    additional_incident: additional,
                });
            }

            i += 1;
            break;
        }
    }

    Ok(records)
}

pub fn build_items(records: Vec<LrsRecord>) -> Result<Vec<VertexItem>> {
    let mut items: Vec<VertexItem> = Vec::new();
    let mut vertex_to_id: BTreeMap<Vec<String>, usize> = BTreeMap::new();

    for rec in records {
        let mut incident = rec.cobasis.clone();
        incident.extend(rec.additional_incident.iter().copied());
        incident = sorted_unique(incident);

        let id = match vertex_to_id.get(&rec.vertex) {
            Some(&id) => {
                if items[id].incident != incident {
                    bail!(
                        "inconsistent incident list for vertex {:?}: first={:?}, new={:?}",
                        rec.vertex,
                        items[id].incident,
                        incident
                    );
                }
                id
            }
            None => {
                let id = items.len();

                vertex_to_id.insert(rec.vertex.clone(), id);

                items.push(VertexItem {
                    vertex: rec.vertex.clone(),
                    incident: incident.clone(),
                    simplices: Vec::new(),
                });

                id
            }
        };

        let simplex = rec
            .cobasis
            .iter()
            .map(|j| {
                items[id]
                    .incident
                    .iter()
                    .position(|x| x == j)
                    .ok_or_else(|| anyhow!("cobasis row {j} is not in incident list"))
            })
            .collect::<Result<Vec<_>>>()?;

        items[id].simplices.push(simplex);
    }

    Ok(items)
}

pub fn choose_default_k0(items: &[VertexItem]) -> Result<usize> {
    items
        .iter()
        .enumerate()
        .filter(|(_, item)| !item.simplices.is_empty())
        .min_by_key(|(_, item)| item.simplices.len())
        .map(|(k, _)| k)
        .ok_or_else(|| anyhow!("cannot choose k0: no vertex item has simplices"))
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

pub fn decode_simplex_global(item: &VertexItem, simplex_index: usize) -> Result<Vec<usize>> {
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

pub fn root_certificate(h: &HRep, items: &[VertexItem], k0: usize) -> Result<Root> {
    let item = items
        .get(k0)
        .ok_or_else(|| anyhow!("invalid k0={k0}; only {} items", items.len()))?;

    if item.simplices.is_empty() {
        bail!("root vertex item k0={k0} has no simplices");
    }

    let rows = decode_simplex_global(item, 0)
        .with_context(|| format!("cannot decode root simplex 0 for k0={k0}"))?;

    if rows.len() != h.d {
        bail!("root simplex has {} rows, expected d={}", rows.len(), h.d);
    }

    let rmat = row_matrix_columns_from_global_rows(h, &rows);
    let inv = invert_matrix(&rmat).context("root basis matrix is singular")?;

    let prod = mat_mul(&inv, &rmat);
    if prod != identity(h.d) {
        bail!("internal error: computed inverse does not multiply to identity");
    }

    let mut same_label_separators = Vec::new();

    for sidx in 1..item.simplices.len() {
        let other_rows = decode_simplex_global(item, sidx)?;

        let inverse_row = (0..h.d)
            .find(|&a| {
                let beta = &inv[a];
                other_rows.iter().all(|&j| dot(beta, &h.a[j]) <= Q::zero())
            })
            .ok_or_else(|| {
                anyhow!(
                    "could not find inverse-row separator for simplex {sidx} of vertex item {k0}"
                )
            })?;

        same_label_separators.push(inverse_row);
    }

    Ok(Root {
        k0,
        rows,
        inverse_rows: inv
            .iter()
            .map(|row| row.iter().map(q_to_string).collect())
            .collect(),
        same_label_separators,
    })
}

pub fn to_certificate(items: Vec<VertexItem>, root: Option<Root>) -> Certificate {
    Certificate { items, root }
}
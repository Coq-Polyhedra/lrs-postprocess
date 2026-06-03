use anyhow::{anyhow, bail, Context, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, HashMap};
use std::fs;

use crate::certificate::{Certificate, GraphLabel, Inequality, ItemSimplex, LocalSimplexRef, Root, SimplexGraph, VertexCoords, VertexItem};
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


fn parse_rational_parts(s: &str) -> Result<(BigInt, BigInt)> {
    let s = s.trim();

    let (mut num, mut den) = if let Some((a, b)) = s.split_once('/') {
        let num = BigInt::parse_bytes(a.trim().as_bytes(), 10)
            .ok_or_else(|| anyhow!("invalid rational numerator `{}`", a.trim()))?;
        let den = BigInt::parse_bytes(b.trim().as_bytes(), 10)
            .ok_or_else(|| anyhow!("invalid rational denominator `{}`", b.trim()))?;
        (num, den)
    } else {
        let num = BigInt::parse_bytes(s.as_bytes(), 10)
            .ok_or_else(|| anyhow!("invalid integer rational `{s}`"))?;
        (num, BigInt::one())
    };

    if den.is_zero() {
        bail!("invalid rational `{s}` with zero denominator");
    }

    if den.sign() == Sign::Minus {
        num = -num;
        den = -den;
    }

    let g = num.abs().gcd(&den);
    num /= &g;
    den /= &g;

    Ok((num, den))
}

fn vertex_coords_from_strings(v: &[String]) -> Result<VertexCoords> {
    let mut nums = Vec::<BigInt>::with_capacity(v.len());
    let mut dens = Vec::<BigInt>::with_capacity(v.len());

    for s in v {
        let (num, den) = parse_rational_parts(s)
            .with_context(|| format!("failed to parse vertex coordinate `{s}`"))?;
        nums.push(num);
        dens.push(den);
    }

    let lcm = dens
        .iter()
        .fold(BigInt::one(), |acc, den| acc.lcm(den));

    let cleared = nums
        .iter()
        .zip(&dens)
        .map(|(num, den)| (num * (&lcm / den)).to_string())
        .collect();

    Ok(VertexCoords {
        num: cleared,
        den: lcm.to_string(),
    })
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
    let mut vertex_to_id: BTreeMap<VertexCoords, usize> = BTreeMap::new();
    for rec in records {
        let vertex = vertex_coords_from_strings(&rec.vertex)?;

        let mut incident = rec.cobasis.clone();
        incident.extend(rec.additional_incident.iter().copied());
        incident = sorted_unique(incident);

        let id = match vertex_to_id.get(&vertex) {
            Some(&id) => {
                if items[id].incident != incident {
                    bail!(
                        "inconsistent incident list for vertex {:?}: first={:?}, new={:?}",
                        vertex,
                        items[id].incident,
                        incident
                    );
                }
                id
            }
            None => {
                let id = items.len();

                vertex_to_id.insert(vertex.clone(), id);

                items.push(VertexItem {
                    vertex: vertex.clone(),
                    incident: incident.clone(),
                    simplices: Vec::new(),
                });

                id
            }
        };

        let indices = rec
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

        items[id].simplices.push(ItemSimplex {
            indices,
            node: 0,
        });
    }

    Ok(items)
}

fn build_simplex_graph(items: &mut [VertexItem]) -> SimplexGraph {
    let node_count = items.iter().map(|item| item.simplices.len()).sum::<usize>();

    let mut entries = Vec::<(Vec<usize>, usize, usize)>::with_capacity(node_count);

    for (item_index, item) in items.iter().enumerate() {
        for (simplex_index, simplex) in item.simplices.iter().enumerate() {
            let rows = simplex
                .indices
                .iter()
                .filter_map(|&local| item.incident.get(local).copied())
                .collect::<Vec<_>>();

            entries.push((rows, item_index, simplex_index));
        }
    }

    // Canonicalize graph node numbering by lexicographic order of global labels.
    entries.sort_by(|a, b| a.0.cmp(&b.0));

    let mut lbl = Vec::<GraphLabel>::with_capacity(node_count);

    for (node, (rows, item_index, simplex_index)) in entries.into_iter().enumerate() {
        items[item_index].simplices[simplex_index].node = node;
        lbl.push(GraphLabel {
            simplex: rows,
            owner: LocalSimplexRef {
                item: item_index,
                simplex: simplex_index,
            },
        });
    }

    let mut ridge_map: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();

    for (node, label) in lbl.iter().enumerate() {
        let sigma = &label.simplex;
        for r in 0..sigma.len() {
            let mut ridge = sigma.clone();
            ridge.remove(r);
            ridge_map.entry(ridge).or_default().push(node);
        }
    }

    let mut g = vec![Vec::<usize>::new(); node_count];

    for nodes in ridge_map.into_values() {
        if nodes.len() == 2 {
            let a = nodes[0];
            let b = nodes[1];
            if a < node_count && b < node_count && a != b {
                g[a].push(b);
                g[b].push(a);
            }
        }
    }

    for adj in &mut g {
        adj.sort_unstable();
        adj.dedup();
    }

    SimplexGraph { g, lbl }
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
        .indices
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

fn clear_rationals_to_integer_row(values: &[String]) -> Result<Vec<String>> {
    let mut nums = Vec::<BigInt>::with_capacity(values.len());
    let mut dens = Vec::<BigInt>::with_capacity(values.len());

    for s in values {
        let (num, den) = parse_rational_parts(s)
            .with_context(|| format!("failed to parse rational `{s}` while clearing denominators"))?;
        nums.push(num);
        dens.push(den);
    }

    let lcm = dens.iter().fold(BigInt::one(), |acc, den| acc.lcm(den));

    Ok(nums
        .iter()
        .zip(&dens)
        .map(|(num, den)| (num * (&lcm / den)).to_string())
        .collect())
}

fn integer_inequality_from_q_row(a: &[Q], b: &Q) -> Result<Inequality> {
    let mut values = a.iter().map(q_to_string).collect::<Vec<_>>();
    values.push(q_to_string(b));

    let cleared = clear_rationals_to_integer_row(&values)?;
    let rhs = cleared
        .last()
        .cloned()
        .ok_or_else(|| anyhow!("empty inequality row after clearing denominators"))?;
    let coeffs = cleared[..cleared.len() - 1].to_vec();

    Ok(Inequality { a: coeffs, b: rhs })
}

fn certificate_inequalities(h: &HRep) -> Result<Vec<Inequality>> {
    h.a.iter()
        .zip(&h.b)
        .enumerate()
        .map(|(i, (a, b))| {
            integer_inequality_from_q_row(a, b)
                .with_context(|| format!("failed to clear denominators of inequality {i}"))
        })
        .collect()
}

pub fn to_certificate(h: &HRep, mut items: Vec<VertexItem>, root: Root) -> Result<Certificate> {
    let inequalities = certificate_inequalities(h)?;
    let graph = build_simplex_graph(&mut items);
    Ok(Certificate {
        inequalities,
        items,
        graph,
        root,
    })
}

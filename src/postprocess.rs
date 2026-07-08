use anyhow::{anyhow, bail, Context, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;

use crate::certificate::{Certificate, GraphLabel, Inequality, Root, SimplexGraph, VertexCoords, VertexItem};
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
fn sorted_subset(a: &[usize], b: &[usize]) -> bool {
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            i += 1;
            j += 1;
        } else if a[i] > b[j] {
            j += 1;
        } else {
            return false;
        }
    }

    i == a.len()
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

pub fn build_items(records: Vec<LrsRecord>) -> Result<(Vec<VertexItem>, Vec<GraphLabel>)> {
    let mut items: Vec<VertexItem> = Vec::new();
    let mut vertex_to_id: BTreeMap<VertexCoords, usize> = BTreeMap::new();
    let mut pending_labels: Vec<(Vec<usize>, VertexCoords)> = Vec::new();

    for rec in records {
        let vertex = vertex_coords_from_strings(&rec.vertex)?;

        let mut incident = rec.cobasis.clone();
        incident.extend(rec.additional_incident.iter().copied());
        incident = sorted_unique(incident);

        match vertex_to_id.get(&vertex) {
            Some(&id) => {
                if items[id].incident != incident {
                    bail!(
                        "inconsistent incident list for vertex {:?}: first={:?}, new={:?}",
                        vertex,
                        items[id].incident,
                        incident
                    );
                }
            }
            None => {
                let id = items.len();
                vertex_to_id.insert(vertex.clone(), id);
                items.push(VertexItem {
                    incident: incident.clone(),
                    vertex: vertex.clone(),
                });
            }
        }

        let rows = sorted_unique(rec.cobasis.clone());
        if rows.len() != rec.cobasis.len() {
            bail!("cobasis contains duplicate row indices: {:?}", rec.cobasis);
        }

        // The owner is stored by vertex coordinates for now.  We sort items below,
        // then convert owners to the final item indices.
        pending_labels.push((rows, vertex));
    }

    // Canonicalize item ordering for the certificate: items are sorted
    // lexicographically by their incident inequality lists. This makes the
    // Coq-side item array searchable by its first component.
    items.sort_by(|a, b| a.incident.cmp(&b.incident));

    let mut coord_to_sorted_id = BTreeMap::<VertexCoords, usize>::new();
    for (id, item) in items.iter().enumerate() {
        coord_to_sorted_id.insert(item.vertex.clone(), id);
    }

    let mut lbl = Vec::<GraphLabel>::with_capacity(pending_labels.len());
    for (simplex, vertex) in pending_labels {
        let owner = *coord_to_sorted_id
            .get(&vertex)
            .ok_or_else(|| anyhow!("internal error: lost vertex while remapping simplex owners"))?;

        if !sorted_subset(&simplex, &items[owner].incident) {
            bail!(
                "simplex {:?} is not contained in incident list of owner item {}: {:?}",
                simplex,
                owner,
                items[owner].incident
            );
        }

        lbl.push(GraphLabel { simplex, owner });
    }

    // Canonicalize graph node numbering by lexicographic order of global simplices.
    lbl.sort_by(|a, b| a.simplex.cmp(&b.simplex).then(a.owner.cmp(&b.owner)));

    Ok((items, lbl))
}

pub fn build_simplex_graph(lbl: Vec<GraphLabel>) -> Result<SimplexGraph> {
    // For every ridge occurrence (node, r), store the unique adjacent node
    // through the ridge obtained by deleting lbl[node].simplex[r].
    let mut ridge_map: HashMap<Vec<usize>, Vec<(usize, usize)>> = HashMap::new();

    for (node, label) in lbl.iter().enumerate() {
        let sigma = &label.simplex;
        for r in 0..sigma.len() {
            let mut ridge = sigma.clone();
            ridge.remove(r);
            ridge_map.entry(ridge).or_default().push((node, r));
        }
    }

    let mut g = lbl
        .iter()
        .map(|label| vec![usize::MAX; label.simplex.len()])
        .collect::<Vec<_>>();

    for (ridge, occs) in ridge_map {
        if occs.len() != 2 {
            bail!(
                "ridge {:?} has {} incident simplex occurrences, expected exactly 2",
                ridge,
                occs.len()
            );
        }

        let (a, ra) = occs[0];
        let (b, rb) = occs[1];

        if a == b {
            bail!(
                "ridge {:?} is paired with two occurrences of the same graph node {a}",
                ridge
            );
        }

        g[a][ra] = b;
        g[b][rb] = a;
    }

    for (node, adj) in g.iter().enumerate() {
        for (r, &other) in adj.iter().enumerate() {
            if other == usize::MAX {
                bail!("missing graph neighbor for node {node}, ridge position {r}");
            }
        }
    }

    Ok(SimplexGraph { g, lbl })
}


fn build_item_neighbors(graph: &SimplexGraph, item_count: usize) -> Result<Vec<Vec<usize>>> {
    let mut sets = (0..item_count)
        .map(|_| BTreeSet::<usize>::new())
        .collect::<Vec<_>>();

    if graph.g.len() != graph.lbl.len() {
        bail!(
            "internal error: graph.g has length {}, graph.lbl has length {}",
            graph.g.len(),
            graph.lbl.len()
        );
    }

    for (node, adj) in graph.g.iter().enumerate() {
        let v = graph.lbl[node].owner;
        if v >= item_count {
            bail!(
                "internal error: graph node {node} has owner item {v}, but there are only {item_count} items"
            );
        }

        for &other in adj {
            if other >= graph.lbl.len() {
                bail!(
                    "internal error: graph adjacency {node}->{other} is out of range 0..{}",
                    graph.lbl.len()
                );
            }

            let w = graph.lbl[other].owner;
            if w >= item_count {
                bail!(
                    "internal error: graph node {other} has owner item {w}, but there are only {item_count} items"
                );
            }

            if v != w {
                sets[v].insert(w);
            }
        }
    }

    Ok(sets
        .into_iter()
        .map(|set| set.into_iter().collect::<Vec<_>>())
        .collect())
}


pub fn choose_default_k0(items: &[VertexItem], graph: &SimplexGraph) -> Result<usize> {
    let mut counts = vec![0usize; items.len()];

    for (node, label) in graph.lbl.iter().enumerate() {
        if label.owner >= items.len() {
            bail!(
                "graph.lbl[{node}].owner={} out of range 0..{}",
                label.owner,
                items.len()
            );
        }
        counts[label.owner] += 1;
    }

    counts
        .iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
        .min_by_key(|(_, &count)| count)
        .map(|(k, _)| k)
        .ok_or_else(|| anyhow!("cannot choose k0: no vertex item owns a simplex"))
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

pub fn root_certificate(h: &HRep, items: &[VertexItem], graph: &SimplexGraph, k0: usize) -> Result<Root> {
    let item = items
        .get(k0)
        .ok_or_else(|| anyhow!("invalid k0={k0}; only {} items", items.len()))?;

    let rows = graph
        .lbl
        .iter()
        .find(|label| label.owner == k0)
        .map(|label| label.simplex.clone())
        .ok_or_else(|| anyhow!("root vertex item k0={k0} owns no simplex"))?;

    if !sorted_subset(&rows, &item.incident) {
        bail!("root.rows {:?} are not contained in I(k0)={:?}", rows, item.incident);
    }

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

    for label in graph.lbl.iter().filter(|label| label.owner == k0) {
        if label.simplex == rows {
            continue;
        }

        let other_rows = &label.simplex;
        let inverse_row = (0..h.d)
            .find(|&a| {
                let beta = &inv[a];
                other_rows.iter().all(|&j| dot(beta, &h.a[j]) <= Q::zero())
            })
            .ok_or_else(|| {
                anyhow!(
                    "could not find inverse-row separator for same-label simplex {:?} of vertex item {k0}",
                    other_rows
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

pub fn to_certificate(h: &HRep, items: Vec<VertexItem>, graph: SimplexGraph, root: Root) -> Result<Certificate> {
    let inequalities = certificate_inequalities(h)?;
    let neighbors = build_item_neighbors(&graph, items.len())?;
    Ok(Certificate {
        inequalities,
        items,
        graph,
        neighbors,
        root,
    })
}

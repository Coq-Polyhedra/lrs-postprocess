use anyhow::{anyhow, bail, Context, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;

use crate::certificate::{Certificate, FullDimCertificate, GraphLabel, Inequality, Root, SimplexGraph, VertexCoords, VertexItem};
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


fn build_item_neighbors_and_lifts(
    graph: &SimplexGraph,
    item_count: usize,
) -> Result<(Vec<Vec<usize>>, Vec<Vec<(usize, usize)>>)> {
    if graph.g.len() != graph.lbl.len() {
        bail!(
            "simplex graph adjacency has length {}, but labels have length {}",
            graph.g.len(),
            graph.lbl.len()
        );
    }

    // For each directed owner edge v -> w, retain the lexicographically
    // smallest oriented simplex-graph edge (s,t) mapping to it.
    let mut lifts = (0..item_count)
        .map(|_| BTreeMap::<usize, (usize, usize)>::new())
        .collect::<Vec<_>>();

    for (s, adj) in graph.g.iter().enumerate() {
        let v = graph.lbl[s].owner;
        if v >= item_count {
            bail!("graph label {s} has owner {v}, out of range 0..{item_count}");
        }

        for &t in adj {
            if t >= graph.lbl.len() {
                bail!(
                    "graph adjacency[{s}] contains node {t}, out of range 0..{}",
                    graph.lbl.len()
                );
            }
            let w = graph.lbl[t].owner;
            if w >= item_count {
                bail!("graph label {t} has owner {w}, out of range 0..{item_count}");
            }
            if v == w {
                continue;
            }

            lifts[v]
                .entry(w)
                .and_modify(|old| {
                    if (s, t) < *old {
                        *old = (s, t);
                    }
                })
                .or_insert((s, t));
        }
    }

    let mut neighbors = Vec::with_capacity(item_count);
    let mut geom_edge_lifts = Vec::with_capacity(item_count);
    for map in lifts {
        let mut ns = Vec::with_capacity(map.len());
        let mut es = Vec::with_capacity(map.len());
        for (w, edge) in map {
            ns.push(w);
            es.push(edge);
        }
        neighbors.push(ns);
        geom_edge_lifts.push(es);
    }

    Ok((neighbors, geom_edge_lifts))
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



fn q_from_bigint(x: &BigInt) -> Q {
    Q::from_integer(x.clone())
}

fn dot_bigint(a: &[BigInt], x: &[BigInt]) -> BigInt {
    a.iter()
        .zip(x)
        .fold(BigInt::zero(), |acc, (ai, xi)| acc + ai * xi)
}

fn clear_q_vec_to_bigints(xs: &[Q]) -> Vec<BigInt> {
    let lcm = xs
        .iter()
        .fold(BigInt::one(), |acc, x| acc.lcm(x.denom()));

    xs.iter()
        .map(|x| x.numer() * (&lcm / x.denom()))
        .collect()
}

fn integer_inequality_coefficients(h: &HRep) -> Result<Vec<Vec<BigInt>>> {
    h.a.iter()
        .zip(&h.b)
        .enumerate()
        .map(|(i, (a, b))| {
            let mut values = a.iter().map(q_to_string).collect::<Vec<_>>();
            values.push(q_to_string(b));
            let cleared = clear_rationals_to_integer_row(&values)
                .with_context(|| format!("failed to clear denominators of inequality {i}"))?;
            cleared[..cleared.len() - 1]
                .iter()
                .enumerate()
                .map(|(j, s)| {
                    BigInt::parse_bytes(s.as_bytes(), 10).ok_or_else(|| {
                        anyhow!("internal error: cleared coefficient {j} of inequality {i} is not an integer: {s}")
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect()
}

fn root_integer_basis_vectors(a_int: &[Vec<BigInt>], rows: &[usize], d: usize) -> Result<Vec<Vec<BigInt>>> {
    let mut root_matrix = vec![vec![Q::zero(); d]; d];

    for (r, &row_id) in rows.iter().enumerate() {
        for c in 0..d {
            root_matrix[r][c] = q_from_bigint(&a_int[row_id][c]);
        }
    }

    let inv = invert_matrix(&root_matrix).context("root integer basis matrix is singular")?;

    let mut basis = Vec::with_capacity(d);
    for j in 0..d {
        let col = (0..d).map(|r| inv[r][j].clone()).collect::<Vec<_>>();
        let f_j = clear_q_vec_to_bigints(&col);
        basis.push(f_j);
    }

    Ok(basis)
}

fn build_root_m_matrix(a_int: &[Vec<BigInt>], incident: &[usize], basis: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    incident
        .iter()
        .map(|&row_id| {
            basis
                .iter()
                .map(|f_j| dot_bigint(&a_int[row_id], f_j))
                .collect()
        })
        .collect()
}

fn combinations_until<T, F>(n: usize, k: usize, f: &mut F) -> Option<T>
where
    F: FnMut(&[usize]) -> Option<T>,
{
    fn rec<T, F>(n: usize, k: usize, start: usize, cur: &mut Vec<usize>, f: &mut F) -> Option<T>
    where
        F: FnMut(&[usize]) -> Option<T>,
    {
        if cur.len() == k {
            return f(cur);
        }

        let need = k - cur.len();
        for x in start..=n - need {
            cur.push(x);
            if let Some(ans) = rec(n, k, x + 1, cur, f) {
                return Some(ans);
            }
            cur.pop();
        }
        None
    }

    if k > n {
        return None;
    }
    let mut cur = Vec::with_capacity(k);
    rec(n, k, 0, &mut cur, f)
}

fn solve_square_q(a: &[Vec<Q>], b: &[Q]) -> Option<Vec<Q>> {
    let inv = invert_matrix(a).ok()?;
    let n = b.len();
    let mut x = vec![Q::zero(); n];

    for i in 0..n {
        for j in 0..n {
            x[i] += &inv[i][j] * &b[j];
        }
    }

    Some(x)
}

fn q_dot_bigint_row(row: &[BigInt], x: &[Q], support: &[usize]) -> Q {
    support
        .iter()
        .zip(x)
        .fold(Q::zero(), |acc, (&j, xj)| acc + q_from_bigint(&row[j]) * xj)
}

fn sparse_nonnegative_certificate_for_simplex(
    m_matrix: &[Vec<BigInt>],
    incident: &[usize],
    simplex: &[usize],
    d: usize,
) -> Result<Vec<(usize, String)>> {
    let local_rows = simplex
        .iter()
        .map(|&row| {
            incident.binary_search(&row).map_err(|_| {
                anyhow!("same-owner simplex row {row} is not contained in the owner's incident set")
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let constraint_rows = local_rows
        .iter()
        .map(|&p| m_matrix[p].clone())
        .collect::<Vec<_>>();

    for support_size in 1..=d {
        if let Some(answer) = combinations_until(d, support_size, &mut |support| {
            let active_size = support_size - 1;
            combinations_until(constraint_rows.len(), active_size, &mut |active| {
                let mut eqs = Vec::<Vec<Q>>::with_capacity(support_size);
                eqs.push(vec![Q::one(); support_size]);
                for &row_idx in active {
                    eqs.push(
                        support
                            .iter()
                            .map(|&j| q_from_bigint(&constraint_rows[row_idx][j]))
                            .collect(),
                    );
                }

                let mut rhs = vec![Q::zero(); support_size];
                rhs[0] = Q::one();

                let x = solve_square_q(&eqs, &rhs)?;

                if x.iter().any(|v| v < &Q::zero()) {
                    return None;
                }

                if x.iter().all(|v| v.is_zero()) {
                    return None;
                }

                if constraint_rows
                    .iter()
                    .any(|row| q_dot_bigint_row(row, &x, support) > Q::zero())
                {
                    return None;
                }

                let coeffs = clear_q_vec_to_bigints(&x);
                let sparse = support
                    .iter()
                    .zip(coeffs.iter())
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(&j, c)| (j, c.to_string()))
                    .collect::<Vec<_>>();

                if sparse.is_empty() {
                    None
                } else {
                    Some(sparse)
                }
            })
        }) {
            return Ok(answer);
        }
    }

    bail!(
        "could not find a nonnegative sparse root certificate vector for same-owner simplex {:?}",
        simplex
    )
}

pub fn root_certificate(h: &HRep, items: &[VertexItem], graph: &SimplexGraph, k0: usize) -> Result<Root> {
    let item = items
        .get(k0)
        .ok_or_else(|| anyhow!("invalid k0={k0}; only {} items", items.len()))?;

    let simplex_id = graph
        .lbl
        .iter()
        .position(|label| label.owner == k0)
        .ok_or_else(|| anyhow!("root vertex item k0={k0} owns no simplex"))?;

    let root_label = &graph.lbl[simplex_id];
    let rows = &root_label.simplex;

    if !sorted_subset(rows, &item.incident) {
        bail!("root simplex {:?} is not contained in I(k0)={:?}", rows, item.incident);
    }

    if rows.len() != h.d {
        bail!("root simplex has {} rows, expected d={}", rows.len(), h.d);
    }

    let a_int = integer_inequality_coefficients(h)?;
    let basis = root_integer_basis_vectors(&a_int, rows, h.d)?;
    let m_matrix = build_root_m_matrix(&a_int, &item.incident, &basis);

    for (root_pos, &row) in rows.iter().enumerate() {
        let local_pos = item.incident.binary_search(&row).map_err(|_| {
            anyhow!("root row {row} is not contained in I(k0)={:?}", item.incident)
        })?;

        for j in 0..h.d {
            let val = &m_matrix[local_pos][j];
            if j == root_pos {
                if val <= &BigInt::zero() {
                    bail!(
                        "root diagonal condition failed for row {row}, column {j}: expected > 0, got {val}"
                    );
                }
            } else if !val.is_zero() {
                bail!(
                    "root off-diagonal condition failed for row {row}, column {j}: expected 0, got {val}"
                );
            }
        }
    }

    let mut q_vectors = Vec::new();
    for (label_id, label) in graph.lbl.iter().enumerate() {
        if label.owner != k0 || label_id == simplex_id {
            continue;
        }

        let q = sparse_nonnegative_certificate_for_simplex(
            &m_matrix,
            &item.incident,
            &label.simplex,
            h.d,
        )
        .with_context(|| {
            format!(
                "failed to build sparse root Q-vector for same-owner simplex graph node {label_id} with simplex {:?}",
                label.simplex
            )
        })?;
        q_vectors.push(q);
    }

    let mut inverse_incident_map = vec![h.a.len(); h.a.len()];
    for (local_pos, &row_id) in item.incident.iter().enumerate() {
        if row_id >= h.a.len() {
            bail!(
                "root owner incident row {row_id} out of range 0..{}",
                h.a.len()
            );
        }
        inverse_incident_map[row_id] = local_pos;
    }

    Ok(Root {
        simplex_id,
        inverse_incident_map,
        basis_vectors: basis
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect(),
        m_matrix: m_matrix
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect(),
        q_vectors,
    })
}



fn vertex_coords_to_q(vertex: &VertexCoords, d: usize) -> Result<Vec<Q>> {
    if vertex.num.len() != d {
        bail!(
            "vertex has {} coordinates, expected dimension {d}",
            vertex.num.len()
        );
    }
    let den = BigInt::parse_bytes(vertex.den.trim().as_bytes(), 10)
        .ok_or_else(|| anyhow!("invalid vertex denominator `{}`", vertex.den))?;
    if den <= BigInt::zero() {
        bail!("vertex denominator must be positive, got {den}");
    }
    vertex
        .num
        .iter()
        .enumerate()
        .map(|(j, x)| {
            let num = BigInt::parse_bytes(x.trim().as_bytes(), 10)
                .ok_or_else(|| anyhow!("invalid vertex numerator at coordinate {j}: `{x}`"))?;
            Ok(Q::new(num, den.clone()))
        })
        .collect()
}

fn pivot_columns_by_row_echelon(matrix: &mut [Vec<Q>]) -> Vec<usize> {
    if matrix.is_empty() {
        return Vec::new();
    }

    let nrows = matrix.len();
    let ncols = matrix[0].len();
    let mut pivot_row = 0usize;
    let mut pivot_columns = Vec::with_capacity(nrows);

    for col in 0..ncols {
        if pivot_row == nrows {
            break;
        }

        let Some(pivot) = (pivot_row..nrows).find(|&r| !matrix[r][col].is_zero()) else {
            continue;
        };

        matrix.swap(pivot_row, pivot);
        let pivot_value = matrix[pivot_row][col].clone();

        // It is enough to eliminate below the pivot.  We do not normalize the
        // pivot row, which avoids many rational divisions and gcd reductions.
        for r in (pivot_row + 1)..nrows {
            if matrix[r][col].is_zero() {
                continue;
            }

            let factor = matrix[r][col].clone() / pivot_value.clone();
            matrix[r][col] = Q::zero();
            for c in (col + 1)..ncols {
                let correction = factor.clone() * matrix[pivot_row][c].clone();
                matrix[r][c] -= correction;
            }
        }

        pivot_columns.push(col);
        pivot_row += 1;
    }

    pivot_columns
}

fn build_full_dim_certificate(
    items: &[VertexItem],
    neighbors: &[Vec<usize>],
    root_owner: usize,
    d: usize,
) -> Result<FullDimCertificate> {
    let root_item = items.get(root_owner).ok_or_else(|| {
        anyhow!(
            "root owner {root_owner} is out of bounds for {} items",
            items.len()
        )
    })?;
    let root_neighbors = neighbors.get(root_owner).ok_or_else(|| {
        anyhow!(
            "missing geometric-neighbor list for root owner {root_owner}"
        )
    })?;

    if root_neighbors.len() < d {
        bail!(
            "root vertex {root_owner} has only {} geometric neighbors, fewer than dimension {d}",
            root_neighbors.len()
        );
    }

    let x0 = vertex_coords_to_q(&root_item.vertex, d)?;

    // Build once the d x deg(v*) rational matrix whose columns are w - v*.
    let mut edge_matrix = vec![vec![Q::zero(); root_neighbors.len()]; d];
    for (col, &neighbor) in root_neighbors.iter().enumerate() {
        let w_item = items.get(neighbor).ok_or_else(|| {
            anyhow!(
                "neighbors[{root_owner}] contains item {neighbor}, but there are only {} items",
                items.len()
            )
        })?;
        let w = vertex_coords_to_q(&w_item.vertex, d)?;
        for row in 0..d {
            edge_matrix[row][col] = w[row].clone() - x0[row].clone();
        }
    }

    // A single row-echelon decomposition gives pivot columns of the original
    // matrix, hence d linearly independent neighbor directions.
    let pivot_columns = pivot_columns_by_row_echelon(&mut edge_matrix);
    if pivot_columns.len() < d {
        bail!(
            "neighbor directions at root vertex {root_owner} have rank {}, expected {d}",
            pivot_columns.len()
        );
    }

    let selected_neighbors = pivot_columns
        .into_iter()
        .take(d)
        .map(|col| root_neighbors[col])
        .collect::<Vec<_>>();

    // Convert only the selected d neighbors to rational coordinates, then
    // clear denominators once for v* and those neighbors.
    let selected_points = selected_neighbors
        .iter()
        .map(|&neighbor| vertex_coords_to_q(&items[neighbor].vertex, d))
        .collect::<Result<Vec<_>>>()?;

    let mut q = BigInt::one();
    for coord in &x0 {
        q = q.lcm(coord.denom());
    }
    for x in &selected_points {
        for coord in x {
            q = q.lcm(coord.denom());
        }
    }
    if q <= BigInt::zero() {
        bail!("internal error: common denominator for full-dimensionality certificate is {q}");
    }

    let p = x0
        .iter()
        .map(|x| x.numer() * (&q / x.denom()))
        .collect::<Vec<_>>();

    // Store R by columns: r[j] is the integer direction numerator
    // q (w^j - v*).
    let mut r = vec![vec![BigInt::zero(); d]; d];
    for j in 0..d {
        for k in 0..d {
            let w_num = selected_points[j][k].numer()
                * (&q / selected_points[j][k].denom());
            r[j][k] = w_num - &p[k];
        }
    }

    // Convert the column representation to the usual row representation only
    // for the exact inversion. The matrix entry R_{k,j} is r[j][k].
    let r_q = (0..d)
        .map(|k| {
            (0..d)
                .map(|j| Q::from_integer(r[j][k].clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let inv = invert_matrix(&r_q).context(
        "pivot columns unexpectedly produced a singular integer direction matrix",
    )?;

    let mut scale = BigInt::one();
    for row in &inv {
        for x in row {
            scale = scale.lcm(x.denom());
        }
    }
    if scale <= BigInt::zero() {
        bail!("internal error: inverse-matrix denominator scale is {scale}");
    }

    let u = inv
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| x.numer() * (&scale / x.denom()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    Ok(FullDimCertificate {
        denominator: q.to_string(),
        point: p.iter().map(ToString::to_string).collect(),
        directions: r
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect(),
        left_inverse: u
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect(),
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
    let (neighbors, geom_edge_lifts) =
        build_item_neighbors_and_lifts(&graph, items.len())?;
    let root_owner = graph
        .lbl
        .get(root.simplex_id)
        .ok_or_else(|| anyhow!("root simplex id {} is out of bounds", root.simplex_id))?
        .owner;
    let full_dim = build_full_dim_certificate(&items, &neighbors, root_owner, h.d)?;
    Ok(Certificate {
        n_inequalities: h.a.len(),
        dimension: h.d,
        inequalities,
        items,
        graph,
        neighbors,
        geom_edge_lifts,
        full_dim,
        root,
    })
}

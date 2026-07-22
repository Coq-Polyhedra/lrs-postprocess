use anyhow::{anyhow, bail, Context, Result};
use rug::{Complete, Integer};
use std::collections::{hash_map::Entry, BTreeMap, HashMap};
use std::fs::{self, File};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, BufReader};

use crate::certificate::{FullDimCertificate, GraphLabel, Inequality, Root, SimplexGraph, VertexCoords, VertexItem};
use crate::numerics::{invert_matrix, parse_q, q_to_string, Q};

#[derive(Debug, Clone)]
pub struct HRep {
    /// Matrix rows for A x <= b.
    pub a: Vec<Vec<Q>>,
    pub b: Vec<Q>,
    pub d: usize,
}

fn parse_rational_parts(s: &str) -> Result<(Integer, Integer)> {
    let value = parse_q(s)?;
    Ok((value.numer().clone(), value.denom().clone()))
}

fn vertex_coords_from_strings(v: &[String]) -> Result<VertexCoords> {
    let mut nums = Vec::<Integer>::with_capacity(v.len());
    let mut dens = Vec::<Integer>::with_capacity(v.len());

    for s in v {
        let (num, den) = parse_rational_parts(s)
            .with_context(|| format!("failed to parse vertex coordinate `{s}`"))?;
        nums.push(num);
        dens.push(den);
    }

    let lcm = dens
        .iter()
        .fold(Integer::from(1), |acc, den| acc.lcm_ref(den).complete());

    let cleared = nums
        .iter()
        .zip(&dens)
        .map(|(num, den)| {
            let mut value = num.clone();
            value *= (&lcm / den).complete();
            value.to_string()
        })
        .collect();

    Ok(VertexCoords {
        num: cleared,
        den: lcm.to_string(),
    })
}

fn vertex_coords_from_common_denominator(
    denominator: &str,
    numerators: &[&str],
) -> Result<VertexCoords> {
    let denominator = Integer::parse(denominator.trim())
        .ok()
        .map(Integer::from)
        .context("invalid common vertex denominator")?;

    if denominator <= 0 {
        bail!("common vertex denominator must be positive, got {denominator}");
    }

    let numerators = numerators
        .iter()
        .enumerate()
        .map(|(j, numerator)| {
            Integer::parse(numerator.trim())
                .ok()
                .map(Integer::from)
                .map(|value| value.to_string())
                .with_context(|| format!("invalid common vertex numerator {j}: `{numerator}`"))
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(VertexCoords {
        num: numerators,
        den: denominator.to_string(),
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

/// Parse lrs output produced with H-representation input and directly build
/// the vertex items and global simplex labels.
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
///
/// This is deliberately a streaming operation.  In particular, it does not
/// retain the `.ext` text or a record containing the (usually very large)
/// complete incidence list for every basis.  Only the d-element simplex and a
/// provisional owner id survive each iteration.
pub fn parse_lrs_ext_items(
    path: &str,
    d: usize,
    m_ineq: usize,
) -> Result<(Vec<VertexItem>, Vec<GraphLabel>)> {
    let file = File::open(path).with_context(|| format!("failed to open `{path}`"))?;
    let mut lines = BufReader::new(file).lines();

    let mut items = Vec::<VertexItem>::new();
    let mut incident_to_id = BTreeMap::<Vec<usize>, usize>::new();
    let mut labels = Vec::<GraphLabel>::new();

    while let Some(line) = lines.next() {
        let line = line?;
        let line = line.trim();

        if !(line.starts_with("V#") && line.contains(" facets ")) {
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
        cobasis.sort_unstable();
        if cobasis.windows(2).any(|pair| pair[0] == pair[1]) {
            bail!("cobasis contains duplicate row indices: {cobasis:?}");
        }

        // Pair this cobasis/incidence line with the next output row. A
        // `*vertexcommon` line is data, despite using lrs's usual comment
        // prefix, so recognize it before skipping other `*` lines.
        let (vline, common_denominator) = loop {
            let vline = lines
                .next()
                .transpose()?
                .with_context(|| format!("missing vertex/ray row after `{line}`"))?;
            let trimmed = vline.trim();
            let vtoks = trimmed.split_whitespace().collect::<Vec<_>>();

            if vtoks.first() == Some(&"*vertexcommon") {
                if vtoks.len() != d + 2 {
                    bail!(
                        "common-denominator vertex row has {} values after its tag, expected {} (denominator plus {d} numerators): `{}`",
                        vtoks.len().saturating_sub(1),
                        d + 1,
                        trimmed
                    );
                }
                break (vline, true);
            }

            if trimmed.is_empty()
                || trimmed.starts_with('*')
                || trimmed.eq_ignore_ascii_case("begin")
                || trimmed.eq_ignore_ascii_case("end")
            {
                continue;
            }

            if !is_vertex_or_ray_row(trimmed, d) {
                bail!(
                    "expected vertex/ray row immediately after cobasis line, got `{}` after `{}`",
                    trimmed,
                    line
                );
            }
            break (vline, false);
        };

        let vtoks = vline.split_whitespace().collect::<Vec<_>>();
        if !common_denominator && vtoks[0] == "0" {
            continue;
        }

        let mut incident = cobasis.clone();
        incident.extend(additional.iter().copied());
        incident = sorted_unique(incident);

        let owner = if let Some(&id) = incident_to_id.get(&incident) {
            id
        } else {
            let vertex = if common_denominator {
                vertex_coords_from_common_denominator(vtoks[1], &vtoks[2..])?
            } else {
                let coordinates = vtoks[1..]
                    .iter()
                    .map(|value| (*value).to_string())
                    .collect::<Vec<_>>();
                vertex_coords_from_strings(&coordinates)?
            };
            let id = items.len();
            items.push(VertexItem {
                incident: incident.clone(),
                vertex,
            });
            incident_to_id.insert(incident, id);
            id
        };

        labels.push(GraphLabel {
            simplex: cobasis,
            owner,
        });
    }

    // Canonicalize item ordering for the certificate: items are sorted
    // lexicographically by their incident inequality lists. This makes the
    // Coq-side item array searchable by its first component.
    let mut indexed_items = items.into_iter().enumerate().collect::<Vec<_>>();
    indexed_items.sort_unstable_by(|(_, a), (_, b)| a.incident.cmp(&b.incident));

    let mut old_to_new = vec![usize::MAX; indexed_items.len()];
    let mut items = Vec::with_capacity(indexed_items.len());
    for (new_id, (old_id, item)) in indexed_items.into_iter().enumerate() {
        old_to_new[old_id] = new_id;
        items.push(item);
    }

    for label in &mut labels {
        label.owner = old_to_new[label.owner];
        if !sorted_subset(&label.simplex, &items[label.owner].incident) {
            bail!(
                "simplex {:?} is not contained in incident list of owner item {}: {:?}",
                label.simplex,
                label.owner,
                items[label.owner].incident
            );
        }
    }

    // Canonicalize graph node numbering by lexicographic order of global simplices.
    labels.sort_unstable_by(|a, b| a.simplex.cmp(&b.simplex).then(a.owner.cmp(&b.owner)));

    Ok((items, labels))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct RidgeFingerprint {
    xor: u64,
    sum: u64,
}

impl RidgeFingerprint {
    fn remove(self, value_fingerprint: (u64, u64)) -> Self {
        let (x, y) = value_fingerprint;
        Self {
            xor: self.xor ^ x,
            sum: self.sum.wrapping_sub(y),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct RidgeRef {
    node: usize,
    removed_pos: usize,
}

#[derive(Debug)]
struct RidgeGroup {
    first: RidgeRef,
    paired: bool,
}

/// Normally a fingerprint identifies a single ridge, so keep that group
/// inline. The vector is used only if distinct ridges have the same two-word
/// fingerprint. Exact comparison below makes such collisions harmless.
#[derive(Debug)]
struct RidgeBucket {
    first: RidgeGroup,
    collisions: Vec<RidgeGroup>,
}

/// Ridge fingerprints have already been mixed. Avoid applying SipHash to both
/// words on every incidence; exact key equality and `same_ridge` still handle
/// all collisions.
#[derive(Default)]
struct RidgeKeyHasher(u64);

impl Hasher for RidgeKeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = (self.0 ^ value)
            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
            .rotate_left(27);
    }
}

type RidgeMap = HashMap<
    RidgeFingerprint,
    RidgeBucket,
    BuildHasherDefault<RidgeKeyHasher>,
>;

fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

fn fingerprint_value(value: usize) -> (u64, u64) {
    let value = value as u64;
    (
        mix64(value),
        mix64(value ^ 0xd6e8_feb8_6659_fd93),
    )
}

fn simplex_fingerprint(
    simplex: &[usize],
    value_fingerprints: &[(u64, u64)],
) -> RidgeFingerprint {
    let mut fingerprint = RidgeFingerprint { xor: 0, sum: 0 };
    for &value in simplex {
        let (x, y) = value_fingerprints[value];
        fingerprint.xor ^= x;
        fingerprint.sum = fingerprint.sum.wrapping_add(y);
    }
    fingerprint
}

fn same_ridge(lbl: &[GraphLabel], a: RidgeRef, b: RidgeRef) -> bool {
    let sa = &lbl[a.node].simplex;
    let sb = &lbl[b.node].simplex;

    if sa.len() != sb.len()
        || a.removed_pos >= sa.len()
        || b.removed_pos >= sb.len()
    {
        return false;
    }

    // Both simplices are sorted. Compare the two views without allocating the
    // vectors obtained by deleting the indicated positions.
    for k in 0..sa.len().saturating_sub(1) {
        let ia = k + if k >= a.removed_pos { 1 } else { 0 };
        let ib = k + if k >= b.removed_pos { 1 } else { 0 };
        if sa[ia] != sb[ib] {
            return false;
        }
    }
    true
}

fn materialize_ridge(lbl: &[GraphLabel], ridge: RidgeRef) -> Vec<usize> {
    lbl[ridge.node]
        .simplex
        .iter()
        .enumerate()
        .filter_map(|(pos, &value)| (pos != ridge.removed_pos).then_some(value))
        .collect()
}

pub fn build_simplex_graph(lbl: Vec<GraphLabel>) -> Result<SimplexGraph> {
    // A ridge occurrence is represented only by its simplex and the position
    // to skip. Its two-word fingerprint is obtained in O(1) after computing
    // one fingerprint per simplex. Matching views are compared exactly, so
    // fingerprint collisions can affect performance but never correctness.
    let ridge_incidence_count = lbl
        .iter()
        .fold(0usize, |count, label| count.saturating_add(label.simplex.len()));
    let mut ridge_map = RidgeMap::with_capacity_and_hasher(
        ridge_incidence_count / 2,
        BuildHasherDefault::default(),
    );

    // Row indices are dense in lrs certificates. Hash every possible row only
    // once instead of running the mixing function for every occurrence.
    let value_fingerprints = lbl
        .iter()
        .flat_map(|label| label.simplex.iter().copied())
        .max()
        .map(|max_value| (0..=max_value).map(fingerprint_value).collect::<Vec<_>>())
        .unwrap_or_default();

    let mut g = lbl
        .iter()
        .map(|label| vec![usize::MAX; label.simplex.len()])
        .collect::<Vec<_>>();

    for (node, label) in lbl.iter().enumerate() {
        let sigma = &label.simplex;
        let simplex_fingerprint = simplex_fingerprint(sigma, &value_fingerprints);

        for removed_pos in 0..sigma.len() {
            let current = RidgeRef { node, removed_pos };
            let fingerprint =
                simplex_fingerprint.remove(value_fingerprints[sigma[removed_pos]]);

            let bucket = match ridge_map.entry(fingerprint) {
                Entry::Vacant(entry) => {
                    entry.insert(RidgeBucket {
                        first: RidgeGroup {
                            first: current,
                            paired: false,
                        },
                        collisions: Vec::new(),
                    });
                    continue;
                }
                Entry::Occupied(entry) => entry.into_mut(),
            };

            let group = if same_ridge(&lbl, bucket.first.first, current) {
                &mut bucket.first
            } else if let Some(pos) = bucket
                .collisions
                .iter()
                .position(|group| same_ridge(&lbl, group.first, current))
            {
                &mut bucket.collisions[pos]
            } else {
                bucket.collisions.push(RidgeGroup {
                    first: current,
                    paired: false,
                });
                continue;
            };

            if group.paired {
                bail!(
                    "ridge {:?} has at least 3 incident simplex occurrences, expected exactly 2",
                    materialize_ridge(&lbl, current)
                );
            }

            let first = group.first;
            if first.node == current.node {
                bail!(
                    "ridge {:?} is paired with two occurrences of the same graph node {}",
                    materialize_ridge(&lbl, current),
                    current.node
                );
            }

            group.paired = true;
            g[first.node][first.removed_pos] = current.node;
            g[current.node][current.removed_pos] = first.node;
        }
    }

    for bucket in ridge_map.values() {
        for group in std::iter::once(&bucket.first).chain(bucket.collisions.iter()) {
            if !group.paired {
                bail!(
                    "ridge {:?} has 1 incident simplex occurrence, expected exactly 2",
                    materialize_ridge(&lbl, group.first)
                );
            }
        }
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


pub fn build_item_neighbors_and_lifts(
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



fn q_from_integer(x: &Integer) -> Q {
    Q::from(x.clone())
}

fn dot_integer(a: &[Integer], x: &[Integer]) -> Integer {
    let mut result = Integer::new();
    for (ai, xi) in a.iter().zip(x) {
        result += ai * xi;
    }
    result
}

fn clear_q_vec_to_integers(xs: &[Q]) -> Vec<Integer> {
    let lcm = xs
        .iter()
        .fold(Integer::from(1), |acc, x| acc.lcm_ref(x.denom()).complete());

    xs.iter()
        .map(|x| {
            let mut value = x.numer().clone();
            value *= (&lcm / x.denom()).complete();
            value
        })
        .collect()
}

fn integer_inequality_coefficients(h: &HRep) -> Result<Vec<Vec<Integer>>> {
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
                    Integer::parse(s)
                        .ok()
                        .map(Integer::from)
                        .ok_or_else(|| {
                        anyhow!("internal error: cleared coefficient {j} of inequality {i} is not an integer: {s}")
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect()
}

fn root_integer_basis_vectors(a_int: &[Vec<Integer>], rows: &[usize], d: usize) -> Result<Vec<Vec<Integer>>> {
    let mut root_matrix = vec![vec![Q::new(); d]; d];

    for (r, &row_id) in rows.iter().enumerate() {
        for c in 0..d {
            root_matrix[r][c] = q_from_integer(&a_int[row_id][c]);
        }
    }

    let inv = invert_matrix(&root_matrix).context("root integer basis matrix is singular")?;

    let mut basis = Vec::with_capacity(d);
    for j in 0..d {
        let col = (0..d).map(|r| inv[r][j].clone()).collect::<Vec<_>>();
        let f_j = clear_q_vec_to_integers(&col);
        basis.push(f_j);
    }

    Ok(basis)
}

fn build_root_m_matrix(a_int: &[Vec<Integer>], incident: &[usize], basis: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    incident
        .iter()
        .map(|&row_id| {
            basis
                .iter()
                .map(|f_j| dot_integer(&a_int[row_id], f_j))
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
    let mut x = vec![Q::new(); n];

    for i in 0..n {
        for j in 0..n {
            let term: Q = (&inv[i][j] * &b[j]).complete();
            x[i] += &term;
        }
    }

    Some(x)
}

fn q_dot_integer_row(row: &[Integer], x: &[Q], support: &[usize]) -> Q {
    let mut result = Q::new();
    for (&j, xj) in support.iter().zip(x) {
        result += q_from_integer(&row[j]) * xj;
    }
    result
}

fn sparse_nonnegative_certificate_for_simplex(
    m_matrix: &[Vec<Integer>],
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
                eqs.push(vec![Q::from(1); support_size]);
                for &row_idx in active {
                    eqs.push(
                        support
                            .iter()
                            .map(|&j| q_from_integer(&constraint_rows[row_idx][j]))
                            .collect(),
                    );
                }

                let mut rhs = vec![Q::new(); support_size];
                rhs[0] = Q::from(1);

                let x = solve_square_q(&eqs, &rhs)?;

                if x.iter().any(|v| v < &0) {
                    return None;
                }

                if x.iter().all(|v| v == &0) {
                    return None;
                }

                if constraint_rows
                    .iter()
                    .any(|row| q_dot_integer_row(row, &x, support) > 0)
                {
                    return None;
                }

                let coeffs = clear_q_vec_to_integers(&x);
                let sparse = support
                    .iter()
                    .zip(coeffs.iter())
                    .filter(|(_, c)| *c != &0)
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
                if val <= &0 {
                    bail!(
                        "root diagonal condition failed for row {row}, column {j}: expected > 0, got {val}"
                    );
                }
            } else if val != &0 {
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
    let den = Integer::parse(vertex.den.trim())
        .ok()
        .map(Integer::from)
        .ok_or_else(|| anyhow!("invalid vertex denominator `{}`", vertex.den))?;
    if den <= 0 {
        bail!("vertex denominator must be positive, got {den}");
    }
    vertex
        .num
        .iter()
        .enumerate()
        .map(|(j, x)| {
            let num = Integer::parse(x.trim())
                .ok()
                .map(Integer::from)
                .ok_or_else(|| anyhow!("invalid vertex numerator at coordinate {j}: `{x}`"))?;
            Ok(Q::from((num, den.clone())))
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

        let Some(pivot) = (pivot_row..nrows).find(|&r| matrix[r][col] != 0) else {
            continue;
        };

        matrix.swap(pivot_row, pivot);
        let pivot_value = matrix[pivot_row][col].clone();

        // It is enough to eliminate below the pivot.  We do not normalize the
        // pivot row, which avoids many rational divisions and gcd reductions.
        for r in (pivot_row + 1)..nrows {
            if matrix[r][col] == 0 {
                continue;
            }

            let factor = matrix[r][col].clone() / pivot_value.clone();
            matrix[r][col] = Q::new();
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

pub fn build_full_dim_certificate(
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
    let mut edge_matrix = vec![vec![Q::new(); root_neighbors.len()]; d];
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

    let mut q = Integer::from(1);
    for coord in &x0 {
        q = q.lcm_ref(coord.denom()).complete();
    }
    for x in &selected_points {
        for coord in x {
            q = q.lcm_ref(coord.denom()).complete();
        }
    }
    if q <= 0 {
        bail!("internal error: common denominator for full-dimensionality certificate is {q}");
    }

    let p = x0
        .iter()
        .map(|x| {
            let mut value = x.numer().clone();
            value *= (&q / x.denom()).complete();
            value
        })
        .collect::<Vec<_>>();

    // Store R by columns: r[j] is the integer direction numerator
    // q (w^j - v*).
    let mut r = vec![vec![Integer::new(); d]; d];
    for j in 0..d {
        for k in 0..d {
            let mut w_num = selected_points[j][k].numer().clone();
            w_num *= (&q / selected_points[j][k].denom()).complete();
            r[j][k] = w_num - &p[k];
        }
    }

    // Convert the column representation to the usual row representation only
    // for the exact inversion. The matrix entry R_{k,j} is r[j][k].
    let r_q = (0..d)
        .map(|k| {
            (0..d)
                .map(|j| Q::from(r[j][k].clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let inv = invert_matrix(&r_q).context(
        "pivot columns unexpectedly produced a singular integer direction matrix",
    )?;

    let mut scale = Integer::from(1);
    for row in &inv {
        for x in row {
            scale = scale.lcm_ref(x.denom()).complete();
        }
    }
    if scale <= 0 {
        bail!("internal error: inverse-matrix denominator scale is {scale}");
    }

    let u = inv
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| {
                    let mut value = x.numer().clone();
                    value *= (&scale / x.denom()).complete();
                    value
                })
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn temporary_ext_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "lrs-postprocess-{name}-{}-{}.ext",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ))
    }

    #[test]
    fn parses_consecutive_common_denominator_vertices() {
        let path = temporary_ext_path("vertexcommon");
        let ext = "\
V#1 R#0 B#1 h=0 facets  1 2 I#2 det= 1
*vertexcommon 3 1 -2
V#2 R#0 B#2 h=1 facets  3 4 I#2 det= 1
*vertexcommon 5 6 7
";
        fs::write(&path, ext).unwrap();

        let (items, labels) = parse_lrs_ext_items(path.to_str().unwrap(), 2, 4).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(items.len(), 2);
        assert_eq!(labels.len(), 2);
        assert_eq!(items[0].vertex.den, "3");
        assert_eq!(items[0].vertex.num, ["1", "-2"]);
        assert_eq!(items[1].vertex.den, "5");
        assert_eq!(items[1].vertex.num, ["6", "7"]);
    }
}

fn clear_rationals_to_integer_row(values: &[String]) -> Result<Vec<String>> {
    let mut nums = Vec::<Integer>::with_capacity(values.len());
    let mut dens = Vec::<Integer>::with_capacity(values.len());

    for s in values {
        let (num, den) = parse_rational_parts(s)
            .with_context(|| format!("failed to parse rational `{s}` while clearing denominators"))?;
        nums.push(num);
        dens.push(den);
    }

    let lcm = dens
        .iter()
        .fold(Integer::from(1), |acc, den| acc.lcm_ref(den).complete());

    Ok(nums
        .iter()
        .zip(&dens)
        .map(|(num, den)| {
            let mut value = num.clone();
            value *= (&lcm / den).complete();
            value.to_string()
        })
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

pub fn certificate_inequalities(h: &HRep) -> Result<Vec<Inequality>> {
    h.a.iter()
        .zip(&h.b)
        .enumerate()
        .map(|(i, (a, b))| {
            integer_inequality_from_q_row(a, b)
                .with_context(|| format!("failed to clear denominators of inequality {i}"))
        })
        .collect()
}

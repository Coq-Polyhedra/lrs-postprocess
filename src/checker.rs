use anyhow::{anyhow, bail, Context, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::time::Instant;

use crate::certificate::{read_certificate, Certificate, Root, VertexCoords};
use crate::numerics::{dot, mat_mul, parse_q, q_to_string, Q};
use crate::postprocess::{parse_lrs_hrep, HRep};


fn parse_bigint_decimal(s: &str, what: &str) -> Result<BigInt> {
    BigInt::parse_bytes(s.trim().as_bytes(), 10)
        .ok_or_else(|| anyhow!("invalid integer `{}` for {what}", s.trim()))
}

fn check_vertex_coords(v: &VertexCoords, expected_dim: usize, what: &str) -> Result<()> {
    if v.num.len() != expected_dim {
        bail!(
            "{what}: vertex numerator vector has dimension {}, expected {}",
            v.num.len(),
            expected_dim
        );
    }

    for (i, x) in v.num.iter().enumerate() {
        let _ = parse_bigint_decimal(x, &format!("{what}: vertex numerator {i}"))?;
    }

    let den = parse_bigint_decimal(&v.den, &format!("{what}: vertex denominator"))?;
    if den.sign() != Sign::Plus {
        bail!("{what}: vertex denominator must be positive, got {}", v.den);
    }

    Ok(())
}

fn parse_q_vec(v: &[String]) -> Result<Vec<Q>> {
    v.iter().map(|s| parse_q(s)).collect()
}

fn parse_q_matrix(m: &[Vec<String>]) -> Result<Vec<Vec<Q>>> {
    m.iter().map(|row| parse_q_vec(row)).collect()
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

fn clear_rationals_to_integer_row(values: &[String]) -> Result<Vec<BigInt>> {
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
        .map(|(num, den)| num * (&lcm / den))
        .collect())
}

fn integer_inequality_from_hrep(h: &HRep, i: usize) -> Result<(Vec<BigInt>, BigInt)> {
    let mut values = h.a[i]
        .iter()
        .map(q_to_string)
        .collect::<Vec<_>>();
    values.push(q_to_string(&h.b[i]));

    let cleared = clear_rationals_to_integer_row(&values)?;
    let b = cleared
        .last()
        .cloned()
        .ok_or_else(|| anyhow!("empty inequality row after clearing denominators"))?;
    let a = cleared[..cleared.len() - 1].to_vec();
    Ok((a, b))
}

fn identity(n: usize) -> Vec<Vec<Q>> {
    let mut id = vec![vec![Q::zero(); n]; n];

    for i in 0..n {
        id[i][i] = Q::one();
    }

    id
}

fn strictly_sorted_len(v: &[usize]) -> Result<usize> {
    for i in 1..v.len() {
        if v[i - 1] >= v[i] {
            bail!(
                "expected strictly sorted list, but found {} followed by {} at positions {} and {}",
                v[i - 1],
                v[i],
                i - 1,
                i
            );
        }
    }

    Ok(v.len())
}

fn contains_sorted(v: &[usize], x: usize) -> bool {
    v.binary_search(&x).is_ok()
}

fn lex_less(a: &[usize], b: &[usize]) -> bool {
    a < b
}

fn check_strictly_lex_sorted(xs: &[Vec<usize>], what: &str) -> Result<()> {
    for i in 1..xs.len() {
        if !lex_less(&xs[i - 1], &xs[i]) {
            bail!(
                "{what} is not strictly lexicographically sorted at positions {} and {}: {:?} then {:?}",
                i - 1,
                i,
                xs[i - 1],
                xs[i]
            );
        }
    }

    Ok(())
}

fn row_matrix_columns_from_global_rows(h: &HRep, rows: &[usize]) -> Result<Vec<Vec<Q>>> {
    if strictly_sorted_len(rows)? != h.d {
        bail!(
            "expected {} strictly sorted global rows for a full-dimensional simplex, got {}",
            h.d,
            rows.len()
        );
    }

    let mut m = vec![vec![Q::zero(); h.d]; h.d];

    for (col, &j) in rows.iter().enumerate() {
        if j >= h.a.len() {
            bail!("global row index {j} out of range 0..{}", h.a.len());
        }

        for row in 0..h.d {
            m[row][col] = h.a[j][row].clone();
        }
    }

    Ok(m)
}


fn check_explicit_sizes(h: &HRep, cert: &Certificate) -> Result<()> {
    if cert.n_inequalities != h.a.len() {
        bail!(
            "certificate declares {} inequalities, but input has {}",
            cert.n_inequalities,
            h.a.len()
        );
    }
    if cert.dimension != h.d {
        bail!(
            "certificate declares dimension {}, but input has dimension {}",
            cert.dimension,
            h.d
        );
    }
    if cert.inequalities.len() != cert.n_inequalities {
        bail!(
            "certificate declares {} inequalities, but stores {} inequality rows",
            cert.n_inequalities,
            cert.inequalities.len()
        );
    }
    Ok(())
}

fn check_certificate_inequalities(h: &HRep, cert: &Certificate) -> Result<()> {
    if cert.inequalities.len() != h.a.len() {
        bail!(
            "certificate contains {} inequalities, but input has {}",
            cert.inequalities.len(),
            h.a.len()
        );
    }

    for (i, ineq) in cert.inequalities.iter().enumerate() {
        if ineq.a.len() != h.d {
            bail!(
                "certificate inequality {i} has {} coefficients, expected d={}",
                ineq.a.len(),
                h.d
            );
        }

        let a = ineq
            .a
            .iter()
            .enumerate()
            .map(|(j, s)| parse_bigint_decimal(s, &format!("inequality {i} coefficient {j}")))
            .collect::<Result<Vec<_>>>()?;
        let b = parse_bigint_decimal(&ineq.b, &format!("inequality {i} rhs"))?;

        let (expected_a, expected_b) = integer_inequality_from_hrep(h, i)
            .with_context(|| format!("failed to clear input inequality {i}"))?;

        if a != expected_a {
            bail!(
                "certificate inequality {i} integer coefficients do not match input H-representation after clearing denominators"
            );
        }

        if b != expected_b {
            bail!(
                "certificate inequality {i} integer rhs does not match input H-representation after clearing denominators"
            );
        }
    }

    Ok(())
}

fn check_items_basic(h: &HRep, cert: &Certificate) -> Result<()> {
    for (k, item) in cert.items.iter().enumerate() {
        check_vertex_coords(&item.vertex, h.d, &format!("item {k}"))
            .with_context(|| format!("item {k}: invalid vertex coordinates"))?;

        strictly_sorted_len(&item.incident)
            .with_context(|| format!("item {k}: incident list is not strictly sorted"))?;

        for (pos, &j) in item.incident.iter().enumerate() {
            if j >= h.a.len() {
                bail!(
                    "item {k}: incident[{pos}]={j} out of range 0..{}",
                    h.a.len()
                );
            }
        }
    }

    for k in 1..cert.items.len() {
        if cert.items[k - 1].incident >= cert.items[k].incident {
            bail!(
                "items are not strictly lexicographically sorted by incident sets at positions {} and {}: {:?} then {:?}",
                k - 1,
                k,
                cert.items[k - 1].incident,
                cert.items[k].incident
            );
        }
    }

    Ok(())
}


fn sorted_subset(a: &[usize], b: &[usize]) -> bool {
    // Assumes both are strictly sorted.
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

fn check_incident_sets_antichain(cert: &Certificate) -> Result<()> {
    let mut sizes = Vec::with_capacity(cert.items.len());

    for (k, item) in cert.items.iter().enumerate() {
        let size = strictly_sorted_len(&item.incident)
            .with_context(|| format!("item {k}: incident list is not strictly sorted"))?;
        sizes.push(size);
    }

    for i in 0..cert.items.len() {
        for j in (i + 1)..cert.items.len() {
            let si = sizes[i];
            let sj = sizes[j];

            if si < sj {
                if sorted_subset(&cert.items[i].incident, &cert.items[j].incident) {
                    bail!(
                        "incident set of item {i} is strictly contained in incident set of item {j}"
                    );
                }
            } else if sj < si {
                if sorted_subset(&cert.items[j].incident, &cert.items[i].incident) {
                    bail!(
                        "incident set of item {j} is strictly contained in incident set of item {i}"
                    );
                }
            }
            // If si == sj, strict inclusion is impossible.
        }
    }

    Ok(())
}

fn check_graph_labels(h: &HRep, cert: &Certificate) -> Result<()> {
    if cert.graph.g.len() != cert.graph.lbl.len() {
        bail!(
            "graph.g has length {}, but graph.lbl has length {}",
            cert.graph.g.len(),
            cert.graph.lbl.len()
        );
    }

    let labels = cert
        .graph
        .lbl
        .iter()
        .map(|label| label.simplex.clone())
        .collect::<Vec<_>>();

    check_strictly_lex_sorted(&labels, "graph.lbl[*].simplex")?;

    for (node, label) in cert.graph.lbl.iter().enumerate() {
        let lbl = &label.simplex;

        let len = strictly_sorted_len(lbl)
            .with_context(|| format!("graph.lbl[{node}].simplex is not strictly sorted"))?;

        if len != h.d {
            bail!(
                "graph.lbl[{node}].simplex has length {}, expected d={}",
                len,
                h.d
            );
        }

        for (pos, &j) in lbl.iter().enumerate() {
            if j >= h.a.len() {
                bail!(
                    "graph.lbl[{node}].simplex[{pos}]={j} out of range 0..{}",
                    h.a.len()
                );
            }
        }

        let owner = label.owner;
        let item = cert.items.get(owner).ok_or_else(|| {
            anyhow!(
                "graph.lbl[{node}].owner={owner} out of range 0..{}",
                cert.items.len()
            )
        })?;

        if !sorted_subset(lbl, &item.incident) {
            bail!(
                "graph.lbl[{node}].simplex={:?} is not contained in I(owner={})={:?}",
                lbl,
                owner,
                item.incident
            );
        }
    }

    Ok(())
}

fn diff_pos(a: &[usize], b: &[usize]) -> Result<usize> {
    let mut missing = Vec::new();

    for (pos, &x) in a.iter().enumerate() {
        if !contains_sorted(b, x) {
            missing.push(pos);
        }
    }

    if missing.len() != 1 {
        bail!(
            "expected labels to differ by exactly one element, got {} missing elements from {:?} to {:?}",
            missing.len(),
            a,
            b
        );
    }

    Ok(missing[0])
}

fn check_simplex_graph(h: &HRep, cert: &Certificate) -> Result<()> {
    check_graph_labels(h, cert)?;

    let m = cert.graph.lbl.len();

    for (node, adj) in cert.graph.g.iter().enumerate() {
        if adj.len() != h.d {
            bail!(
                "graph.g[{node}] has length {}, expected d={} (one neighbor per ridge position)",
                adj.len(),
                h.d
            );
        }

        let sigma = &cert.graph.lbl[node].simplex;

        for (r, &other) in adj.iter().enumerate() {
            if other >= m {
                bail!(
                    "graph.g[{node}][{r}]={other} out of range 0..{}",
                    m
                );
            }

            if other == node {
                bail!("graph.g[{node}][{r}] is a self-loop");
            }

            let other_sigma = &cert.graph.lbl[other].simplex;

            let missing_from_node = diff_pos(sigma, other_sigma).with_context(|| {
                format!("graph edge {node}->{other} is not a ridge adjacency")
            })?;

            if missing_from_node != r {
                bail!(
                    "graph.g[{node}][{r}]={other} does not correspond to deleting position {r}; \
                     the unique element of graph.lbl[{node}] missing from graph.lbl[{other}] is at position {missing_from_node}"
                );
            }

            let missing_from_other = diff_pos(other_sigma, sigma).with_context(|| {
                format!("graph edge {other}->{node} is not a ridge adjacency")
            })?;

            if cert.graph.g[other].len() != h.d {
                bail!(
                    "graph.g[{other}] has length {}, expected d={} (one neighbor per ridge position)",
                    cert.graph.g[other].len(),
                    h.d
                );
            }

            if cert.graph.g[other][missing_from_other] != node {
                bail!(
                    "graph adjacency is not reciprocal through the corresponding ridge: \
                     graph.g[{node}][{r}]={other}, but graph.g[{other}][{missing_from_other}]={}",
                    cert.graph.g[other][missing_from_other]
                );
            }
        }
    }

    Ok(())
}

fn sorted_difference_values(a: &[usize], b: &[usize]) -> Vec<usize> {
    // Return a \ b, assuming both lists are strictly sorted.
    let mut i = 0;
    let mut j = 0;
    let mut out = Vec::new();

    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            out.push(a[i]);
            i += 1;
        } else {
            j += 1;
        }
    }

    out.extend_from_slice(&a[i..]);
    out
}

fn sorted_not_subset_witness(a: &[usize], b: &[usize]) -> Option<usize> {
    // Return x in a \ b, if one exists. Both lists are assumed strictly sorted.
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            return Some(a[i]);
        } else {
            j += 1;
        }
    }

    a.get(i).copied()
}

fn graph_item_neighbors(cert: &Certificate) -> Result<Vec<Vec<usize>>> {
    let item_count = cert.items.len();
    let mut sets = (0..item_count)
        .map(|_| std::collections::BTreeSet::<usize>::new())
        .collect::<Vec<_>>();

    if cert.graph.g.len() != cert.graph.lbl.len() {
        bail!(
            "graph.g has length {}, but graph.lbl has length {}",
            cert.graph.g.len(),
            cert.graph.lbl.len()
        );
    }

    for (node, adj) in cert.graph.g.iter().enumerate() {
        let v = cert.graph.lbl[node].owner;
        if v >= item_count {
            bail!("graph.lbl[{node}].owner={v} out of range 0..{item_count}");
        }

        for &other in adj {
            if other >= cert.graph.lbl.len() {
                bail!(
                    "graph.g[{node}] contains node {other}, out of range 0..{}",
                    cert.graph.lbl.len()
                );
            }

            let w = cert.graph.lbl[other].owner;
            if w >= item_count {
                bail!("graph.lbl[{other}].owner={w} out of range 0..{item_count}");
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

fn check_neighbor_lists(cert: &Certificate) -> Result<()> {
    let n = cert.items.len();

    if cert.neighbors.len() != n {
        bail!(
            "cert.neighbors has length {}, expected one list for each of the {n} items",
            cert.neighbors.len()
        );
    }

    for (v, ns) in cert.neighbors.iter().enumerate() {
        strictly_sorted_len(ns)
            .with_context(|| format!("neighbors[{v}] is not strictly sorted"))?;

        for &w in ns {
            if w >= n {
                bail!("neighbors[{v}] contains item {w}, out of range 0..{n}");
            }
            if w == v {
                bail!("neighbors[{v}] contains a self-neighbor");
            }
            if !contains_sorted(&cert.neighbors[w], v) {
                bail!(
                    "neighbor relation is not symmetric: {w} occurs in neighbors[{v}], but {v} does not occur in neighbors[{w}]"
                );
            }
        }
    }

    let expected = graph_item_neighbors(cert)
        .context("failed to derive item-neighbor lists from simplex graph")?;
    if cert.neighbors != expected {
        bail!(
            "cert.neighbors does not match cross-label adjacencies of the simplex graph: certificate={:?}, graph-derived={:?}",
            cert.neighbors,
            expected
        );
    }

    Ok(())
}

fn check_local_edge_test(cert: &Certificate) -> Result<()> {
    // For fixed v, let D_w = I(v) \ I(w) for w in neighbors[v].
    // The local edge test, for both orders of every pair of distinct
    // neighbors, is exactly pairwise incomparability of the D_w.
    for (v, ns) in cert.neighbors.iter().enumerate() {
        let iv = &cert.items[v].incident;
        let diffs = ns
            .iter()
            .map(|&w| sorted_difference_values(iv, &cert.items[w].incident))
            .collect::<Vec<_>>();

        for j in 0..diffs.len() {
            for k in (j + 1)..diffs.len() {
                if sorted_not_subset_witness(&diffs[j], &diffs[k]).is_none() {
                    bail!(
                        "local edge test failed at item {v}: D_{{{}}}=I({v})\\I({}) is contained in D_{{{}}}=I({v})\\I({}); D_{{{}}}={:?}, D_{{{}}}={:?}",
                        ns[j],
                        ns[j],
                        ns[k],
                        ns[k],
                        ns[j],
                        diffs[j],
                        ns[k],
                        diffs[k]
                    );
                }
                if sorted_not_subset_witness(&diffs[k], &diffs[j]).is_none() {
                    bail!(
                        "local edge test failed at item {v}: D_{{{}}}=I({v})\\I({}) is contained in D_{{{}}}=I({v})\\I({}); D_{{{}}}={:?}, D_{{{}}}={:?}",
                        ns[k],
                        ns[k],
                        ns[j],
                        ns[j],
                        ns[k],
                        diffs[k],
                        ns[j],
                        diffs[j]
                    );
                }
            }
        }
    }

    Ok(())
}


fn parse_bigint_vec(xs: &[String], what: &str) -> Result<Vec<BigInt>> {
    xs.iter()
        .enumerate()
        .map(|(i, s)| parse_bigint_decimal(s, &format!("{what}[{i}]")))
        .collect()
}

fn parse_bigint_matrix(m: &[Vec<String>], rows: usize, cols: usize, what: &str) -> Result<Vec<Vec<BigInt>>> {
    if m.len() != rows {
        bail!("{what} has {} rows, expected {rows}", m.len());
    }

    m.iter()
        .enumerate()
        .map(|(i, row)| {
            if row.len() != cols {
                bail!("{what}[{i}] has {} columns, expected {cols}", row.len());
            }
            parse_bigint_vec(row, &format!("{what}[{i}]"))
        })
        .collect()
}

fn dot_bigint(a: &[BigInt], x: &[BigInt]) -> BigInt {
    a.iter()
        .zip(x)
        .fold(BigInt::zero(), |acc, (ai, xi)| acc + ai * xi)
}

fn parse_certificate_integer_rows(cert: &Certificate, d: usize) -> Result<Vec<Vec<BigInt>>> {
    cert.inequalities
        .iter()
        .enumerate()
        .map(|(i, ineq)| {
            if ineq.a.len() != d {
                bail!(
                    "certificate inequality {i} has {} coefficients, expected d={d}",
                    ineq.a.len()
                );
            }
            parse_bigint_vec(&ineq.a, &format!("inequality {i}.a"))
        })
        .collect()
}

fn check_root(h: &HRep, cert: &Certificate) -> Result<()> {
    let root = &cert.root;

    if root.simplex_id >= cert.graph.lbl.len() {
        bail!(
            "root.simplex_id={} out of range 0..{}",
            root.simplex_id,
            cert.graph.lbl.len()
        );
    }

    let root_label = &cert.graph.lbl[root.simplex_id];
    let k0 = root_label.owner;
    let item = cert
        .items
        .get(k0)
        .ok_or_else(|| anyhow!("root owner k0={k0} out of range 0..{}", cert.items.len()))?;

    if root_label.simplex.len() != h.d {
        bail!(
            "root simplex graph node {} has length {}, expected d={}",
            root.simplex_id,
            root_label.simplex.len(),
            h.d
        );
    }

    if !sorted_subset(&root_label.simplex, &item.incident) {
        bail!(
            "root simplex {:?} is not contained in I(root owner {})={:?}",
            root_label.simplex,
            k0,
            item.incident
        );
    }

    if root.inverse_incident_map.len() != cert.n_inequalities {
        bail!(
            "root.inverse_incident_map has length {}, expected n={}",
            root.inverse_incident_map.len(),
            cert.n_inequalities
        );
    }

    let sentinel = cert.n_inequalities;
    for (global_row, &local_pos) in root.inverse_incident_map.iter().enumerate() {
        if local_pos == sentinel {
            if item.incident.binary_search(&global_row).is_ok() {
                bail!(
                    "root.inverse_incident_map[{global_row}] is the sentinel {sentinel}, but row {global_row} belongs to I(root owner {k0})"
                );
            }
        } else {
            if local_pos >= item.incident.len() {
                bail!(
                    "root.inverse_incident_map[{global_row}]={local_pos}, expected a local position below {} or sentinel {sentinel}",
                    item.incident.len()
                );
            }
            if item.incident[local_pos] != global_row {
                bail!(
                    "root.inverse_incident_map[{global_row}]={local_pos}, but I(root owner {k0})[{local_pos}]={}",
                    item.incident[local_pos]
                );
            }
        }
    }
    for (local_pos, &global_row) in item.incident.iter().enumerate() {
        if root.inverse_incident_map[global_row] != local_pos {
            bail!(
                "root inverse map mismatch: I(root owner {k0})[{local_pos}]={global_row}, but inverse_incident_map[{global_row}]={}",
                root.inverse_incident_map[global_row]
            );
        }
    }

    let basis = parse_bigint_matrix(&root.basis_vectors, h.d, h.d, "root.basis_vectors")?;
    let m_matrix = parse_bigint_matrix(
        &root.m_matrix,
        item.incident.len(),
        h.d,
        "root.m_matrix",
    )?;
    let a_int = parse_certificate_integer_rows(cert, h.d)?;

    for (p, &row_id) in item.incident.iter().enumerate() {
        let a = a_int.get(row_id).ok_or_else(|| {
            anyhow!(
                "item {k0}: incident row {row_id} out of range 0..{}",
                a_int.len()
            )
        })?;

        for j in 0..h.d {
            let expected = dot_bigint(a, &basis[j]);
            if m_matrix[p][j] != expected {
                bail!(
                    "root M mismatch at active row position {p} (global row {row_id}), column {j}: got {}, expected {}",
                    m_matrix[p][j],
                    expected
                );
            }
        }
    }

    for (root_pos, &row_id) in root_label.simplex.iter().enumerate() {
        let local_pos = root.inverse_incident_map[row_id];
        if local_pos == sentinel {
            bail!(
                "root row {row_id} is not contained in I(root owner {})={:?}",
                k0,
                item.incident
            );
        }

        for j in 0..h.d {
            let val = &m_matrix[local_pos][j];
            if j == root_pos {
                if val <= &BigInt::zero() {
                    bail!(
                        "root diagonal condition failed for row {row_id}, column {j}: expected > 0, got {val}"
                    );
                }
            } else if !val.is_zero() {
                bail!(
                    "root off-diagonal condition failed for row {row_id}, column {j}: expected 0, got {val}"
                );
            }
        }
    }

    let same_owner_nonroot = cert
        .graph
        .lbl
        .iter()
        .enumerate()
        .filter(|(label_id, label)| label.owner == k0 && *label_id != root.simplex_id)
        .collect::<Vec<_>>();

    if root.q_vectors.len() != same_owner_nonroot.len() {
        bail!(
            "root.q_vectors has length {}, expected {} same-owner non-root simplices",
            root.q_vectors.len(),
            same_owner_nonroot.len()
        );
    }

    for (q_idx, ((label_id, label), sparse)) in same_owner_nonroot
        .iter()
        .zip(root.q_vectors.iter())
        .enumerate()
    {
        let mut last_index = None::<usize>;
        let mut nonzero = false;
        let mut parsed = Vec::<(usize, BigInt)>::with_capacity(sparse.len());

        for (entry_pos, (coord, value_s)) in sparse.iter().enumerate() {
            if *coord >= h.d {
                bail!(
                    "root.q_vectors[{q_idx}][{entry_pos}] has coordinate {}, expected in 0..{}",
                    coord,
                    h.d
                );
            }
            if let Some(prev) = last_index {
                if prev >= *coord {
                    bail!(
                        "root.q_vectors[{q_idx}] is not strictly sorted by coordinate: {prev} then {coord}"
                    );
                }
            }
            last_index = Some(*coord);

            let value = parse_bigint_decimal(value_s, &format!("root.q_vectors[{q_idx}][{entry_pos}].value"))?;
            if value < BigInt::zero() {
                bail!(
                    "root.q_vectors[{q_idx}][{entry_pos}] has negative value {value}"
                );
            }
            if !value.is_zero() {
                nonzero = true;
            }
            parsed.push((*coord, value));
        }

        if !nonzero {
            bail!("root.q_vectors[{q_idx}] for graph node {label_id} is zero");
        }

        for &row_id in &label.simplex {
            if row_id >= root.inverse_incident_map.len() {
                bail!(
                    "same-owner simplex graph node {label_id} contains row {row_id} out of range 0..{}",
                    root.inverse_incident_map.len()
                );
            }
            let local_pos = root.inverse_incident_map[row_id];
            if local_pos == sentinel {
                bail!(
                    "same-owner simplex graph node {label_id} contains row {row_id}, not in I(root owner {})={:?}",
                    k0,
                    item.incident
                );
            }

            let mut sum = BigInt::zero();
            for (coord, value) in &parsed {
                sum += value * &m_matrix[local_pos][*coord];
            }

            if sum > BigInt::zero() {
                bail!(
                    "root.q_vectors[{q_idx}] fails on same-owner graph node {label_id}, row {row_id}: dot = {sum} > 0"
                );
            }
        }
    }

    Ok(())
}


fn time_check<T, F>(name: &str, f: F) -> Result<T>
where
    F: FnOnce() -> Result<T>,
{
    let t0 = Instant::now();
    let res = f();
    let dt = t0.elapsed();
    eprintln!("{name}: {:.6} s", dt.as_secs_f64());
    res
}

pub fn check_certificate(ine_path: &str, certificate_path: &str) -> Result<()> {
    let h = time_check("parse H-representation", || {
        parse_lrs_hrep(ine_path).context("failed to parse H-representation")
    })?;

    let cert = time_check("read certificate", || read_certificate(certificate_path))?;

    time_check("explicit size check", || check_explicit_sizes(&h, &cert))?;

    time_check("certificate inequality check", || {
        check_certificate_inequalities(&h, &cert).context("certificate inequality check failed")
    })?;

    time_check("basic item consistency check", || {
        check_items_basic(&h, &cert).context("basic item consistency check failed")
    })?;

    time_check("simplex graph check", || {
        check_simplex_graph(&h, &cert).context("simplex graph check failed")
    })?;

    time_check("root check", || {
        check_root(&h, &cert).context("root check failed")
    })?;

    time_check("neighbor-list check", || {
        check_neighbor_lists(&cert).context("neighbor-list check failed")
    })?;

    time_check("local edge test", || {
        check_local_edge_test(&cert).context("local edge test failed")
    })?;

    Ok(())
}

use anyhow::{anyhow, bail, Context, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use crate::certificate::{read_certificate, Certificate, LocalSimplexRef, Root, VertexCoords, VertexItem};
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

fn decode_simplex_global(item: &VertexItem, simplex_index: usize) -> Result<Vec<usize>> {
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

        for (sidx, simplex) in item.simplices.iter().enumerate() {
            let mut global_rows = Vec::with_capacity(simplex.indices.len());

            for &local in &simplex.indices {
                if local >= item.incident.len() {
                    bail!(
                        "item {k}, simplex {sidx}: local index {local} out of range 0..{}",
                        item.incident.len()
                    );
                }

                global_rows.push(item.incident[local]);
            }

            let len = strictly_sorted_len(&global_rows).with_context(|| {
                format!("item {k}, simplex {sidx}: decoded global rows are not strictly sorted")
            })?;

            if len != h.d {
                bail!(
                    "item {k}, simplex {sidx}: decoded simplex has length {}, expected d={}",
                    len,
                    h.d
                );
            }
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

/// Checks the two inverse maps between local simplices and graph nodes.
///
/// Forward map f: local simplex -> graph node is stored in `ItemSimplex::node`.
/// Reverse map g: graph node -> local simplex is stored in `graph.lbl[node].owner`.
///
/// The checker verifies both identities: g(f(s)) = s and f(g(k)) = k.
fn check_node_bijection(cert: &Certificate) -> Result<Vec<LocalSimplexRef>> {
    let m = cert.graph.lbl.len();

    if cert.graph.g.len() != m {
        bail!(
            "graph.g has length {}, but graph.lbl has length {}",
            cert.graph.g.len(),
            m
        );
    }

    for (k, item) in cert.items.iter().enumerate() {
        for (l, simplex) in item.simplices.iter().enumerate() {
            let node = simplex.node;

            if node >= m {
                bail!("item {k}, simplex {l}: graph node {node} out of range 0..{m}");
            }

            let owner = cert.graph.lbl[node].owner;
            if owner.item != k || owner.simplex != l {
                bail!(
                    "node map identities fail: item {k}, simplex {l} points to node {node}, \
                     but graph.lbl[{node}].owner = ({},{})",
                    owner.item,
                    owner.simplex
                );
            }
        }
    }

    for (node, label) in cert.graph.lbl.iter().enumerate() {
        let owner = label.owner;
        let item = cert.items.get(owner.item).ok_or_else(|| {
            anyhow!(
                "graph.lbl[{node}].owner refers to item {}, but there are only {} items",
                owner.item,
                cert.items.len()
            )
        })?;

        let simplex = item.simplices.get(owner.simplex).ok_or_else(|| {
            anyhow!(
                "graph.lbl[{node}].owner refers to item {}, simplex {}, but item has only {} simplices",
                owner.item,
                owner.simplex,
                item.simplices.len()
            )
        })?;

        if simplex.node != node {
            bail!(
                "node map identities fail: graph.lbl[{node}].owner = ({},{}), but that simplex points to node {}",
                owner.item,
                owner.simplex,
                simplex.node
            );
        }
    }

    Ok(cert.graph.lbl.iter().map(|label| label.owner).collect())
}

fn check_graph_labels(h: &HRep, cert: &Certificate, node_to_simplex: &[LocalSimplexRef]) -> Result<()> {
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

        let owner = node_to_simplex[node];
        let item_index = owner.item;
        let simplex_index = owner.simplex;
        let decoded = decode_simplex_global(&cert.items[item_index], simplex_index)
            .with_context(|| format!("cannot decode item {item_index}, simplex {simplex_index}"))?;

        if decoded != *lbl {
            bail!(
                "graph.lbl[{node}].simplex does not match decoded item {item_index}, simplex {simplex_index}: graph={:?}, decoded={:?}",
                lbl,
                decoded
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
    let node_to_simplex = check_node_bijection(cert)?;
    check_graph_labels(h, cert, &node_to_simplex)?;

    let m = cert.graph.lbl.len();

    for (node, adj) in cert.graph.g.iter().enumerate() {
        strictly_sorted_len(adj)
            .with_context(|| format!("graph.g[{node}] is not strictly sorted"))?;

        if adj.len() != h.d {
            bail!(
                "graph.g[{node}] has length {}, expected d={} (one neighbor per ridge)",
                adj.len(),
                h.d
            );
        }

        let sigma = &cert.graph.lbl[node].simplex;
        let mut seen_missing_positions = vec![false; h.d];

        for (apos, &other) in adj.iter().enumerate() {
            if other >= m {
                bail!(
                    "graph.g[{node}][{apos}]={other} out of range 0..{}",
                    m
                );
            }

            if other == node {
                bail!("graph.g[{node}] contains a self-loop");
            }

            if !contains_sorted(&cert.graph.g[other], node) {
                bail!(
                    "directed graph is not symmetric: {node} lists {other}, but {other} does not list {node}"
                );
            }

            let other_sigma = &cert.graph.lbl[other].simplex;

            let missing_from_node = diff_pos(sigma, other_sigma).with_context(|| {
                format!("graph edge {node}->{other} is not a ridge adjacency")
            })?;

            let _missing_from_other = diff_pos(other_sigma, sigma).with_context(|| {
                format!("graph edge {other}->{node} is not a ridge adjacency")
            })?;

            if seen_missing_positions[missing_from_node] {
                bail!(
                    "graph.g[{node}] has two neighbors through the same ridge position {missing_from_node}"
                );
            }
            seen_missing_positions[missing_from_node] = true;
        }
    }

    Ok(())
}

fn check_root_inverse(h: &HRep, cert: &Certificate, root: &Root) -> Result<Vec<Vec<Q>>> {
    let inv = parse_q_matrix(&root.inverse_rows).context("failed to parse root.inverse_rows")?;

    if inv.len() != h.d || inv.iter().any(|row| row.len() != h.d) {
        bail!("root.inverse_rows must be a {} x {} matrix", h.d, h.d);
    }

    let rmat =
        row_matrix_columns_from_global_rows(h, &root.rows).context("failed to build root matrix")?;

    let prod = mat_mul(&inv, &rmat);

    if prod != identity(h.d) {
        bail!("root.inverse_rows do not invert the root basis matrix");
    }

    let item = cert
        .items
        .get(root.k0)
        .ok_or_else(|| anyhow!("root.k0={} out of range", root.k0))?;

    if item.simplices.is_empty() {
        bail!("root item k0={} has no simplices", root.k0);
    }

    let decoded = decode_simplex_global(item, 0)
        .with_context(|| format!("cannot decode root simplex 0 for k0={}", root.k0))?;

    if decoded != root.rows {
        bail!("root.rows do not match decoded simplex items[k0].simplices[0]");
    }

    Ok(inv)
}

fn check_same_label_separators(
    h: &HRep,
    cert: &Certificate,
    root: &Root,
    inverse_rows: &[Vec<Q>],
) -> Result<()> {
    let item = cert
        .items
        .get(root.k0)
        .ok_or_else(|| anyhow!("root.k0={} out of range", root.k0))?;

    let expected_len = item.simplices.len().saturating_sub(1);

    if root.same_label_separators.len() != expected_len {
        bail!(
            "root.same_label_separators has length {}, expected {}",
            root.same_label_separators.len(),
            expected_len
        );
    }

    for (offset, &inverse_row) in root.same_label_separators.iter().enumerate() {
        let sidx = offset + 1;

        if inverse_row >= h.d {
            bail!(
                "separator for simplex {sidx} has inverse_row={}, expected in 0..{}",
                inverse_row,
                h.d
            );
        }

        let beta = &inverse_rows[inverse_row];

        let rows = decode_simplex_global(item, sidx)
            .with_context(|| format!("cannot decode same-label simplex {sidx}"))?;

        strictly_sorted_len(&rows).with_context(|| {
            format!("same-label simplex {sidx}: decoded global rows are not strictly sorted")
        })?;

        for &j in &rows {
            let val = dot(beta, &h.a[j]);

            if val > Q::zero() {
                bail!(
                    "invalid separator for same-label simplex {sidx}: inverse row {inverse_row} has positive value {val} on row {j}"
                );
            }
        }
    }

    Ok(())
}

fn check_root(h: &HRep, cert: &Certificate) -> Result<()> {
    let root = &cert.root;

    let inverse_rows = check_root_inverse(h, cert, root)?;

    check_same_label_separators(h, cert, root, &inverse_rows)
        .context("same-label separator check failed")?;

    Ok(())
}

pub fn check_certificate(ine_path: &str, certificate_path: &str) -> Result<()> {
    let h = parse_lrs_hrep(ine_path).context("failed to parse H-representation")?;
    let cert = read_certificate(certificate_path)?;

    check_certificate_inequalities(&h, &cert).context("certificate inequality check failed")?;

    check_items_basic(&h, &cert).context("basic item consistency check failed")?;

    check_simplex_graph(&h, &cert).context("simplex graph check failed")?;

    check_root(&h, &cert).context("root check failed")?;

    check_incident_sets_antichain(&cert).context("incident-set antichain check failed")?;

    Ok(())
}

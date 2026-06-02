use anyhow::{anyhow, bail, Context, Result};
use num_traits::{One, Zero};

use crate::certificate::{read_certificate, AdjacentSimplex, Certificate, Root, VertexItem};
use crate::numerics::{dot, mat_mul, parse_q, Q};
use crate::postprocess::{parse_lrs_hrep, HRep};

fn parse_q_vec(v: &[String]) -> Result<Vec<Q>> {
    v.iter().map(|s| parse_q(s)).collect()
}

fn parse_q_matrix(m: &[Vec<String>]) -> Result<Vec<Vec<Q>>> {
    m.iter().map(|row| parse_q_vec(row)).collect()
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

fn ridge_from_simplex(sigma: &[usize], missing_pos: usize) -> Result<Vec<usize>> {
    if missing_pos >= sigma.len() {
        bail!(
            "missing position {missing_pos} out of range 0..{}",
            sigma.len()
        );
    }

    let mut ridge = Vec::with_capacity(sigma.len().saturating_sub(1));

    for (i, &x) in sigma.iter().enumerate() {
        if i != missing_pos {
            ridge.push(x);
        }
    }

    Ok(ridge)
}

fn find_missing_pos_for_ridge(sigma: &[usize], ridge: &[usize]) -> Option<usize> {
    if sigma.len() != ridge.len() + 1 {
        return None;
    }

    (0..sigma.len()).find(|&r| ridge_from_simplex(sigma, r).ok().as_deref() == Some(ridge))
}

fn check_items_basic(h: &HRep, cert: &Certificate) -> Result<()> {
    for (k, item) in cert.items.iter().enumerate() {
        if item.vertex.len() != h.d {
            bail!(
                "item {k}: vertex has dimension {}, expected {}",
                item.vertex.len(),
                h.d
            );
        }

        let _x = parse_q_vec(&item.vertex)
            .with_context(|| format!("item {k}: failed to parse vertex coordinates"))?;

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
            if simplex.adj.len() != h.d {
                bail!(
                    "item {k}, simplex {sidx}: adjacency list has length {}, expected d={}",
                    simplex.adj.len(),
                    h.d
                );
            }

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

            for (r, adj) in simplex.adj.iter().enumerate() {
                if adj.item >= cert.items.len() {
                    bail!(
                        "item {k}, simplex {sidx}, ridge {r}: adjacent item {} out of range 0..{}",
                        adj.item,
                        cert.items.len()
                    );
                }

                if adj.simplex >= cert.items[adj.item].simplices.len() {
                    bail!(
                        "item {k}, simplex {sidx}, ridge {r}: adjacent simplex {} out of range 0..{} for item {}",
                        adj.simplex,
                        cert.items[adj.item].simplices.len(),
                        adj.item
                    );
                }
            }
        }
    }

    Ok(())
}

/// Checks the explicit adjacency pointers stored in every simplex.
///
/// If simplex `(k,s)` has global rows `sigma`, then `sigma \ {sigma[r]}`
/// is the ridge represented by occurrence `(k,s,r)`.  The certificate stores
/// `simplices[s].adj[r] = (k',s')`.  The checker verifies that `(k',s')`
/// really has the same ridge, and that its corresponding adjacency pointer
/// points back to `(k,s)`.
fn check_ridge_adjacencies(h: &HRep, cert: &Certificate) -> Result<()> {
    let d = h.d;

    for (k, item) in cert.items.iter().enumerate() {
        for (sidx, simplex) in item.simplices.iter().enumerate() {
            let sigma = decode_simplex_global(item, sidx)
                .with_context(|| format!("failed decoding item {k}, simplex {sidx}"))?;

            if strictly_sorted_len(&sigma)? != d {
                bail!(
                    "item {k}, simplex {sidx}: decoded simplex has length {}, expected {d}",
                    sigma.len()
                );
            }

            for r in 0..d {
                let adj = simplex.adj[r];

                if adj.item == k && adj.simplex == sidx {
                    bail!("item {k}, simplex {sidx}, ridge {r}: self-adjacency is invalid");
                }

                let ridge = ridge_from_simplex(&sigma, r)?;

                let adj_item = &cert.items[adj.item];
                let adj_sigma = decode_simplex_global(adj_item, adj.simplex).with_context(|| {
                    format!(
                        "failed decoding adjacent simplex ({},{}) from ({},{}) ridge {}",
                        adj.item, adj.simplex, k, sidx, r
                    )
                })?;

                if strictly_sorted_len(&adj_sigma)? != d {
                    bail!(
                        "adjacent simplex ({},{}) has length {}, expected {d}",
                        adj.item,
                        adj.simplex,
                        adj_sigma.len()
                    );
                }

                let Some(adj_missing_pos) = find_missing_pos_for_ridge(&adj_sigma, &ridge) else {
                    bail!(
                        "item {k}, simplex {sidx}, ridge {r} = {:?} is not a ridge of adjacent simplex ({},{}) = {:?}",
                        ridge,
                        adj.item,
                        adj.simplex,
                        adj_sigma
                    );
                };

                let back = cert.items[adj.item].simplices[adj.simplex].adj[adj_missing_pos];
                if back != (AdjacentSimplex { item: k, simplex: sidx }) {
                    bail!(
                        "adjacency is not reciprocal: ({k},{sidx}) ridge {r} points to ({},{}), whose matching ridge {} points to ({},{})",
                        adj.item,
                        adj.simplex,
                        adj_missing_pos,
                        back.item,
                        back.simplex
                    );
                }

                if k != adj.item {
                    let missing_row = sigma[r];
                    let adj_missing_row = adj_sigma[adj_missing_pos];

                    if contains_sorted(&cert.items[adj.item].incident, missing_row) {
                        bail!(
                            "local ridge test failed: simplex ({k},{sidx}) has exchanged row {missing_row}, but this row belongs to incident set of item {}",
                            adj.item
                        );
                    }

                    if contains_sorted(&cert.items[k].incident, adj_missing_row) {
                        bail!(
                            "local ridge test failed: simplex ({},{}) has exchanged row {}, but this row belongs to incident set of item {k}",
                            adj.item,
                            adj.simplex,
                            adj_missing_row
                        );
                    }
                }
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
    let Some(root) = &cert.root else {
        bail!("certificate has no `root` field; cannot check B3");
    };

    let inverse_rows = check_root_inverse(h, cert, root)?;

    check_same_label_separators(h, cert, root, &inverse_rows)
        .context("same-label separator check failed")?;

    Ok(())
}

pub fn check_certificate(ine_path: &str, certificate_path: &str) -> Result<()> {
    let h = parse_lrs_hrep(ine_path).context("failed to parse H-representation")?;
    let cert = read_certificate(certificate_path)?;

    check_items_basic(&h, &cert).context("basic item consistency check failed")?;

    check_ridge_adjacencies(&h, &cert).context("ridge adjacency check failed")?;

    check_root(&h, &cert).context("root check failed")?;

    check_incident_sets_antichain(&cert).context("incident-set antichain check failed")?;

    Ok(())
}

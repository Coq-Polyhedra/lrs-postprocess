use anyhow::{anyhow, bail, Context, Result};
use num_traits::{One, Zero};
use std::collections::HashMap;

use crate::certificate::{read_certificate, Certificate, Root, VertexItem};
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
        .iter()
        .map(|&local| {
            item.incident
                .get(local)
                .copied()
                .ok_or_else(|| anyhow!("local simplex index {local} out of range"))
        })
        .collect()
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
            let mut global_rows = Vec::with_capacity(simplex.len());

            for &local in simplex {
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

fn check_mod2_cycle(h: &HRep, cert: &Certificate) -> Result<()> {
    let d = h.d;
    let mut ridge_parity: HashMap<Vec<usize>, bool> = HashMap::new();

    for (k, item) in cert.items.iter().enumerate() {
        for sidx in 0..item.simplices.len() {
            let sigma = decode_simplex_global(item, sidx)
                .with_context(|| format!("failed decoding item {k}, simplex {sidx}"))?;

            let len = strictly_sorted_len(&sigma).with_context(|| {
                format!("item {k}, simplex {sidx}: decoded global rows are not strictly sorted")
            })?;

            if len != d {
                bail!(
                    "item {k}, simplex {sidx}: decoded simplex has length {}, expected {d}",
                    len
                );
            }

            for r in 0..d {
                let mut ridge = sigma.clone();
                ridge.remove(r);

                let entry = ridge_parity.entry(ridge).or_insert(false);
                *entry = !*entry;
            }
        }
    }

    let bad = ridge_parity
        .into_iter()
        .find(|(_, odd)| *odd)
        .map(|(ridge, _)| ridge);

    if let Some(ridge) = bad {
        bail!("mod-2 cycle check failed: odd ridge {:?}", ridge);
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

    /* check_root_vertex_strictly_maximizes(cert, root, &c_star)
        .context("root direction does not uniquely select the root vertex")?; */

    check_same_label_separators(h, cert, root, &inverse_rows)
        .context("same-label separator check failed")?;

    Ok(())
}

fn is_subset_sorted(a: &[usize], b: &[usize]) -> bool {
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

fn check_incident_antichain(cert: &Certificate) -> Result<()> {
    for k in 0..cert.items.len() {
        for l in 0..cert.items.len() {
            if k == l {
                continue;
            }

            let lk = &cert.items[k].incident;
            let ll = &cert.items[l].incident;

            if is_subset_sorted(lk, ll) {
                bail!(
                    "incident-set antichain check failed: incident set of item {k} is contained in incident set of item {l}"
                );
            }
        }
    }

    Ok(())
}

pub fn check_certificate(ine_path: &str, certificate_path: &str) -> Result<()> {
    let h = parse_lrs_hrep(ine_path).context("failed to parse H-representation")?;
    let cert = read_certificate(certificate_path)?;

    check_items_basic(&h, &cert).context("basic item consistency check failed")?;
    check_mod2_cycle(&h, &cert).context("mod-2 cycle check failed")?;
    check_root(&h, &cert).context("root check failed")?;
    check_incident_antichain(&cert).context("incident-set antichain check failed")?;

    Ok(())
}
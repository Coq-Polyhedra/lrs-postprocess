use anyhow::{bail, Result};
use rug::{Complete, Integer};
use std::io::Write;

use crate::certificate::{
    Certificate, FullDimCertificate, GraphLabel, Inequality, Root, SimplexGraph, VertexCoords,
    VertexItem,
};

#[derive(Clone, Debug)]
enum Descr {
    Int63,
    BigN,
    BigZ,
    Pair(Box<Descr>, Box<Descr>),
    Array(Box<Descr>),
}

fn pair(a: Descr, b: Descr) -> Descr {
    Descr::Pair(Box::new(a), Box::new(b))
}

fn array(a: Descr) -> Descr {
    Descr::Array(Box::new(a))
}

fn write_i63<W: Write>(w: &mut W, x: u64) -> Result<()> {
    if x >= (1u64 << 63) {
        bail!("integer {x} does not fit in Coq int63");
    }
    w.write_all(&x.to_le_bytes())?;
    Ok(())
}

fn write_descr<W: Write>(w: &mut W, d: &Descr) -> Result<()> {
    match d {
        Descr::Int63 => write_i63(w, 0x00),
        Descr::BigN => write_i63(w, 0x01),
        Descr::BigZ => write_i63(w, 0x02),
        Descr::Pair(a, b) => {
            write_i63(w, 0x04)?;
            write_descr(w, a)?;
            write_descr(w, b)
        }
        Descr::Array(a) => {
            write_i63(w, 0x05)?;
            write_descr(w, a)
        }
    }
}

fn write_default<W: Write>(w: &mut W, d: &Descr) -> Result<()> {
    match d {
        Descr::Int63 => write_i63(w, 0),
        Descr::BigN => write_bign(w, &Integer::new()),
        Descr::BigZ => write_bigz(w, &Integer::new()),
        Descr::Pair(a, b) => {
            write_default(w, a)?;
            write_default(w, b)
        }
        Descr::Array(a) => {
            write_i63(w, 0)?; // length
            write_default(w, a) // default element
        }
    }
}

fn write_int63_usize<W: Write>(w: &mut W, x: usize) -> Result<()> {
    write_i63(w, x as u64)
}

fn write_bign<W: Write>(w: &mut W, n: &Integer) -> Result<()> {
    if n < &0 {
        bail!("cannot encode negative integer {n} as BigN");
    }
    if n == &0 {
        write_i63(w, 0)?;
        return Ok(());
    }

    let mask = Integer::from((1u64 << 63) - 1);

    let mut limbs = Vec::<u64>::new();
    let mut x = n.clone();

    while x != 0 {
        let limb_value: Integer = (&x & &mask).complete();
        let limb = limb_value
            .to_u64()
            .ok_or_else(|| anyhow::anyhow!("internal error while extracting 63-bit limb"))?;
        limbs.push(limb);
        x >>= 63usize;
    }

    write_i63(w, limbs.len() as u64)?;
    for limb in limbs {
        write_i63(w, limb)?;
    }

    Ok(())
}

fn write_bigz<W: Write>(w: &mut W, z: &Integer) -> Result<()> {
    let nonnegative = z >= &0;
    write_i63(w, if nonnegative { 1 } else { 0 })?;

    let mut abs = z.clone();
    abs.abs_mut();
    write_bign(w, &abs)
}

fn write_array<W, T, F>(w: &mut W, elem_descr: &Descr, xs: &[T], mut write_elem: F) -> Result<()>
where
    W: Write,
    F: FnMut(&mut W, &T) -> Result<()>,
{
    write_i63(w, xs.len() as u64)?;
    write_default(w, elem_descr)?;
    for x in xs {
        write_elem(w, x)?;
    }
    Ok(())
}

fn inequality_descr() -> Descr {
    // inequality := integer coefficients * integer rhs
    pair(array(Descr::BigZ), Descr::BigZ)
}

fn vertex_coords_descr() -> Descr {
    // vertex := integer numerators * common denominator
    // The denominator is positive, hence encoded as BigN.
    pair(array(Descr::BigZ), Descr::BigN)
}

fn item_descr() -> Descr {
    // flag := local_inequalities * witness_items
    let flag = pair(array(Descr::Int63), array(Descr::Int63));

    // item := incident * (vertex * flag)
    pair(array(Descr::Int63), pair(vertex_coords_descr(), flag))
}

fn graph_label_descr() -> Descr {
    // graph_label := simplex * owner_item
    pair(array(Descr::Int63), Descr::Int63)
}

fn graph_descr() -> Descr {
    // graph := g * lbl
    pair(array(array(Descr::Int63)), array(graph_label_descr()))
}

fn sparse_entry_descr() -> Descr {
    // sparse_entry := coordinate * coefficient
    pair(Descr::Int63, Descr::BigZ)
}

fn sparse_vector_descr() -> Descr {
    // sparse_vector := array sparse_entry
    array(sparse_entry_descr())
}

fn full_dim_descr() -> Descr {
    // full_dim := (point * denominator) * (direction_columns * left_inverse)
    pair(
        pair(array(Descr::BigZ), Descr::BigN),
        pair(array(array(Descr::BigZ)), array(array(Descr::BigZ))),
    )
}

fn root_descr() -> Descr {
    // root := simplex_id * (inverse_incident_map * (basis_vectors * (m_matrix * q_vectors)))
    pair(
        Descr::Int63,
        pair(
            array(Descr::Int63),
            pair(
                array(array(Descr::BigZ)),
                pair(array(array(Descr::BigZ)), array(sparse_vector_descr())),
            ),
        ),
    )
}

fn certificate_descr() -> Descr {
    // The geometric graph and its two parallel lift maps are grouped as:
    //   geom := neighbors * (geom_edge_sources * geom_edge_local_targets)
    // where
    //   geom_edge_sources[v][j] = source simplex node, and
    //   geom_edge_local_targets[v][j] = local index of the target in graph[source].
    //
    // certificate := n_inequalities *
    //   (dimension *
    //    (inequalities *
    //     (items *
    //      (graph *
    //       (geom *
    //        (full_dim * root))))))
    let int_matrix = || array(array(Descr::Int63));
    let geom_descr = pair(int_matrix(), pair(int_matrix(), int_matrix()));

    pair(
        Descr::Int63,
        pair(
            Descr::Int63,
            pair(
                array(inequality_descr()),
                pair(
                    array(item_descr()),
                    pair(
                        graph_descr(),
                        pair(geom_descr, pair(full_dim_descr(), root_descr())),
                    ),
                ),
            ),
        ),
    )
}

fn write_usize_array<W: Write>(w: &mut W, xs: &[usize]) -> Result<()> {
    write_array(w, &Descr::Int63, xs, |w, x| write_int63_usize(w, *x))
}

fn write_bigz_matrix<W: Write>(w: &mut W, m: &[Vec<Integer>]) -> Result<()> {
    let row_descr = array(Descr::BigZ);
    write_array(w, &row_descr, m, |w, row| write_bigz_array(w, row))
}

fn write_usize_matrix<W: Write>(w: &mut W, m: &[Vec<usize>]) -> Result<()> {
    let row_descr = array(Descr::Int63);
    write_array(w, &row_descr, m, |w, row| write_usize_array(w, row))
}

fn write_inequality<W: Write>(w: &mut W, ineq: &Inequality) -> Result<()> {
    // inequality := integer coefficients * integer rhs
    write_bigz_array(w, &ineq.a)?;
    write_bigz(w, &ineq.b)?;
    Ok(())
}

fn write_positive_bign<W: Write>(w: &mut W, n: &Integer) -> Result<()> {
    if n <= &0 {
        bail!("BigN denominator must be positive, got {n}");
    }

    write_bign(w, n)
}

fn write_bigz_array<W: Write>(w: &mut W, xs: &[Integer]) -> Result<()> {
    write_array(w, &Descr::BigZ, xs, write_bigz)
}

fn write_vertex_coords<W: Write>(w: &mut W, vertex: &VertexCoords) -> Result<()> {
    // vertex := integer numerators * common positive denominator
    write_bigz_array(w, &vertex.num)?;
    write_positive_bign(w, &vertex.den)?;
    Ok(())
}

fn write_item<W: Write>(w: &mut W, item: &VertexItem) -> Result<()> {
    // item := incident * (vertex * (local_inequalities * witness_items))
    write_usize_array(w, &item.incident)?;
    write_vertex_coords(w, &item.vertex)?;
    write_usize_array(w, &item.flag.0)?;
    write_usize_array(w, &item.flag.1)?;
    Ok(())
}

fn write_graph_label<W: Write>(w: &mut W, label: &GraphLabel) -> Result<()> {
    // graph_label := simplex * owner_item
    write_usize_array(w, &label.simplex)?;
    write_int63_usize(w, label.owner)?;
    Ok(())
}

fn write_graph<W: Write>(w: &mut W, graph: &SimplexGraph) -> Result<()> {
    // graph := g * lbl
    write_usize_matrix(w, &graph.g)?;

    let label_d = graph_label_descr();
    write_array(w, &label_d, &graph.lbl, |w, label| {
        write_graph_label(w, label)
    })?;

    Ok(())
}

fn write_sparse_entry<W: Write>(w: &mut W, entry: &(usize, Integer)) -> Result<()> {
    write_int63_usize(w, entry.0)?;
    write_bigz(w, &entry.1)?;
    Ok(())
}

fn write_sparse_vector<W: Write>(w: &mut W, sparse: &[(usize, Integer)]) -> Result<()> {
    let entry_d = sparse_entry_descr();
    write_array(w, &entry_d, sparse, |w, entry| write_sparse_entry(w, entry))
}

fn split_geom_edge_lifts(cert: &Certificate) -> Result<(Vec<Vec<usize>>, Vec<Vec<usize>>)> {
    if cert.geom_edge_lifts.len() != cert.neighbors.len() {
        bail!(
            "geom_edge_lifts has {} rows, but neighbors has {} rows",
            cert.geom_edge_lifts.len(),
            cert.neighbors.len()
        );
    }

    let mut sources = Vec::with_capacity(cert.geom_edge_lifts.len());
    let mut local_targets = Vec::with_capacity(cert.geom_edge_lifts.len());

    for (v, (lift_row, neighbor_row)) in cert
        .geom_edge_lifts
        .iter()
        .zip(cert.neighbors.iter())
        .enumerate()
    {
        if lift_row.len() != neighbor_row.len() {
            bail!(
                "geom_edge_lifts[{v}] has length {}, but neighbors[{v}] has length {}",
                lift_row.len(),
                neighbor_row.len()
            );
        }

        let mut source_row = Vec::with_capacity(lift_row.len());
        let mut target_row = Vec::with_capacity(lift_row.len());

        for (j, (&(source, target), &neighbor)) in
            lift_row.iter().zip(neighbor_row.iter()).enumerate()
        {
            let source_label = cert.graph.lbl.get(source).ok_or_else(|| {
                anyhow::anyhow!("geom_edge_lifts[{v}][{j}] has invalid source simplex {source}")
            })?;

            if source_label.owner != v {
                bail!(
                    "geom_edge_lifts[{v}][{j}] has source simplex {source} owned by {}, expected {v}",
                    source_label.owner
                );
            }

            let target_label = cert.graph.lbl.get(target).ok_or_else(|| {
                anyhow::anyhow!("geom_edge_lifts[{v}][{j}] has invalid target simplex {target}")
            })?;

            if target_label.owner != neighbor {
                bail!(
                    "geom_edge_lifts[{v}][{j}] targets simplex {target} owned by {}, expected geometric neighbor {neighbor}",
                    target_label.owner
                );
            }

            let source_neighbors = cert.graph.g.get(source).ok_or_else(|| {
                anyhow::anyhow!(
                    "geom_edge_lifts[{v}][{j}] has source simplex {source}, but the graph has only {} rows",
                    cert.graph.g.len()
                )
            })?;

            let target_local_index = source_neighbors
                .iter()
                .position(|&u| u == target)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "geom_edge_lifts[{v}][{j}] contains ({source}, {target}), but {target} is not a neighbor of {source}"
                    )
                })?;

            source_row.push(source);
            target_row.push(target_local_index);
        }

        sources.push(source_row);
        local_targets.push(target_row);
    }

    Ok((sources, local_targets))
}

fn write_full_dim<W: Write>(w: &mut W, full_dim: &FullDimCertificate) -> Result<()> {
    // full_dim := (point * denominator) * (direction_columns * left_inverse)
    write_bigz_array(w, &full_dim.point)?;
    write_positive_bign(w, &full_dim.denominator)?;
    write_bigz_matrix(w, &full_dim.directions)?;
    write_bigz_matrix(w, &full_dim.left_inverse)?;
    Ok(())
}

fn write_root<W: Write>(w: &mut W, root: &Root) -> Result<()> {
    // root := simplex_id * (inverse_incident_map * (basis_vectors * (m_matrix * q_vectors)))
    write_int63_usize(w, root.simplex_id)?;
    write_usize_array(w, &root.inverse_incident_map)?;
    write_bigz_matrix(w, &root.basis_vectors)?;
    write_bigz_matrix(w, &root.m_matrix)?;

    let sparse_d = sparse_vector_descr();
    write_array(w, &sparse_d, &root.q_vectors, |w, sparse| {
        write_sparse_vector(w, sparse)
    })?;
    Ok(())
}

fn write_certificate_value<W: Write>(w: &mut W, cert: &Certificate) -> Result<()> {
    write_int63_usize(w, cert.n_inequalities)?;
    write_int63_usize(w, cert.dimension)?;

    let inequality_d = inequality_descr();
    write_array(w, &inequality_d, &cert.inequalities, |w, ineq| {
        write_inequality(w, ineq)
    })?;

    let item_d = item_descr();
    write_array(w, &item_d, &cert.items, |w, item| write_item(w, item))?;

    write_graph(w, &cert.graph)?;

    // geom := neighbors * (geom_edge_sources * geom_edge_local_targets)
    write_usize_matrix(w, &cert.neighbors)?;
    let (geom_edge_sources, geom_edge_local_targets) = split_geom_edge_lifts(cert)?;
    write_usize_matrix(w, &geom_edge_sources)?;
    write_usize_matrix(w, &geom_edge_local_targets)?;

    write_full_dim(w, &cert.full_dim)?;
    write_root(w, &cert.root)?;

    Ok(())
}

pub fn write_certificate_bin<W: Write>(
    w: &mut W,
    cert: &crate::certificate::Certificate,
) -> anyhow::Result<()> {
    let d = certificate_descr();
    write_descr(w, &d)?;
    write_certificate_value(w, cert)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(bytes: &[u8]) -> Vec<u64> {
        bytes
            .chunks_exact(8)
            .map(|chunk| u64::from_le_bytes(chunk.try_into().unwrap()))
            .collect()
    }

    #[test]
    fn bign_uses_63_bit_little_endian_limbs() {
        let n = (Integer::from(1) << 126usize) + (Integer::from(7) << 63usize) + 5;
        let mut output = Vec::new();
        write_bign(&mut output, &n).unwrap();
        assert_eq!(words(&output), vec![3, 5, 7, 1]);
    }

    #[test]
    fn bigz_preserves_the_sign_encoding() {
        let mut output = Vec::new();
        write_bigz(&mut output, &Integer::from(-5)).unwrap();
        assert_eq!(words(&output), vec![0, 1, 5]);
    }

    #[test]
    fn item_encoding_appends_the_complete_flag_pair() {
        let mut descriptor = Vec::new();
        write_descr(&mut descriptor, &item_descr()).unwrap();
        assert_eq!(
            words(&descriptor),
            vec![4, 5, 0, 4, 4, 5, 2, 1, 4, 5, 0, 5, 0]
        );

        let item = VertexItem {
            incident: vec![2, 5],
            vertex: VertexCoords {
                num: vec![Integer::from(7), Integer::from(-3)],
                den: Integer::from(11),
            },
            flag: (vec![0, 1], vec![4, 1]),
        };
        let mut value = Vec::new();
        write_item(&mut value, &item).unwrap();
        assert_eq!(
            words(&value),
            vec![
                2, 0, 2, 5, // incident
                2, 1, 0, 1, 1, 7, 0, 1, 3, 1, 11, // vertex
                2, 0, 0, 1, // local flag inequalities
                2, 0, 4, 1, // global witness items
            ]
        );
    }
}

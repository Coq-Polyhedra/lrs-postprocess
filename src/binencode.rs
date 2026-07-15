use anyhow::{bail, Context, Result};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::certificate::{read_certificate, Certificate, FullDimCertificate, GraphLabel, Inequality, Root, SimplexGraph, VertexCoords, VertexItem};

#[derive(Clone, Debug)]
enum Descr {
    Int63,
    BigN,
    BigZ,
    BigQ,
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
        Descr::BigQ => write_i63(w, 0x03),
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
        Descr::BigN => write_bign(w, &BigUint::zero()),
        Descr::BigZ => write_bigz(w, &BigInt::zero()),
        Descr::BigQ => write_bigq_parts(w, &BigInt::zero(), &BigUint::one()),
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

fn write_bign<W: Write>(w: &mut W, n: &BigUint) -> Result<()> {
    if n.is_zero() {
        write_i63(w, 0)?;
        return Ok(());
    }

    let base = BigUint::one() << 63usize;
    let mask = &base - BigUint::one();

    let mut limbs = Vec::<u64>::new();
    let mut x = n.clone();

    while !x.is_zero() {
        let limb = (&x & &mask)
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

fn write_bigz<W: Write>(w: &mut W, z: &BigInt) -> Result<()> {
    let nonnegative = z.sign() != Sign::Minus;
    write_i63(w, if nonnegative { 1 } else { 0 })?;

    let abs = z.abs().to_biguint().unwrap();
    write_bign(w, &abs)
}

fn write_bigq_parts<W: Write>(w: &mut W, num: &BigInt, den: &BigUint) -> Result<()> {
    if den.is_zero() {
        bail!("BigQ denominator is zero");
    }

    write_bigz(w, num)?;
    write_bign(w, den)
}

fn parse_bigq_string(s: &str) -> Result<(BigInt, BigUint)> {
    let s = s.trim();

    let (mut num, mut den) = if let Some((a, b)) = s.split_once('/') {
        let num = BigInt::parse_bytes(a.trim().as_bytes(), 10)
            .ok_or_else(|| anyhow::anyhow!("invalid rational numerator `{}`", a.trim()))?;
        let den = BigInt::parse_bytes(b.trim().as_bytes(), 10)
            .ok_or_else(|| anyhow::anyhow!("invalid rational denominator `{}`", b.trim()))?;
        (num, den)
    } else {
        let num = BigInt::parse_bytes(s.as_bytes(), 10)
            .ok_or_else(|| anyhow::anyhow!("invalid integer rational `{s}`"))?;
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

    Ok((num, den.to_biguint().unwrap()))
}

fn write_bigq_str<W: Write>(w: &mut W, s: &str) -> Result<()> {
    let (num, den) = parse_bigq_string(s)?;
    write_bigq_parts(w, &num, &den)
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
    // item := incident * vertex
    pair(array(Descr::Int63), vertex_coords_descr())
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
        pair(
            array(array(Descr::BigZ)),
            array(array(Descr::BigZ)),
        ),
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
    let geom_descr = pair(
        int_matrix(),
        pair(int_matrix(), int_matrix()),
    );

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
                        pair(
                            geom_descr,
                            pair(full_dim_descr(), root_descr()),
                        ),
                    ),
                ),
            ),
        ),
    )
}

fn write_usize_array<W: Write>(w: &mut W, xs: &[usize]) -> Result<()> {
    write_array(w, &Descr::Int63, xs, |w, x| write_int63_usize(w, *x))
}

fn write_bigq_string_array<W: Write>(w: &mut W, xs: &[String]) -> Result<()> {
    write_array(w, &Descr::BigQ, xs, |w, x| write_bigq_str(w, x))
}

fn write_bigq_string_matrix<W: Write>(w: &mut W, m: &[Vec<String>]) -> Result<()> {
    let row_descr = array(Descr::BigQ);
    write_array(w, &row_descr, m, |w, row| write_bigq_string_array(w, row))
}

fn write_bigz_string_matrix<W: Write>(w: &mut W, m: &[Vec<String>]) -> Result<()> {
    let row_descr = array(Descr::BigZ);
    write_array(w, &row_descr, m, |w, row| write_bigz_string_array(w, row))
}

fn write_usize_matrix<W: Write>(w: &mut W, m: &[Vec<usize>]) -> Result<()> {
    let row_descr = array(Descr::Int63);
    write_array(w, &row_descr, m, |w, row| write_usize_array(w, row))
}

fn write_inequality<W: Write>(w: &mut W, ineq: &Inequality) -> Result<()> {
    // inequality := integer coefficients * integer rhs
    write_bigz_string_array(w, &ineq.a)?;
    write_bigz_str(w, &ineq.b)?;
    Ok(())
}


fn write_bigz_str<W: Write>(w: &mut W, s: &str) -> Result<()> {
    let z = BigInt::parse_bytes(s.trim().as_bytes(), 10)
        .ok_or_else(|| anyhow::anyhow!("invalid integer `{}`", s.trim()))?;
    write_bigz(w, &z)
}

fn write_bign_str<W: Write>(w: &mut W, s: &str) -> Result<()> {
    let n = BigUint::parse_bytes(s.trim().as_bytes(), 10)
        .ok_or_else(|| anyhow::anyhow!("invalid natural integer `{}`", s.trim()))?;

    if n.is_zero() {
        bail!("BigN denominator must be positive, got 0");
    }

    write_bign(w, &n)
}

fn write_bigz_string_array<W: Write>(w: &mut W, xs: &[String]) -> Result<()> {
    write_array(w, &Descr::BigZ, xs, |w, x| write_bigz_str(w, x))
}

fn write_vertex_coords<W: Write>(w: &mut W, vertex: &VertexCoords) -> Result<()> {
    // vertex := integer numerators * common positive denominator
    write_bigz_string_array(w, &vertex.num)?;
    write_bign_str(w, &vertex.den)?;
    Ok(())
}

fn write_item<W: Write>(w: &mut W, item: &VertexItem) -> Result<()> {
    // item := incident * vertex
    write_usize_array(w, &item.incident)?;
    write_vertex_coords(w, &item.vertex)?;
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
    write_array(w, &label_d, &graph.lbl, |w, label| write_graph_label(w, label))?;

    Ok(())
}

fn write_sparse_entry<W: Write>(w: &mut W, entry: &(usize, String)) -> Result<()> {
    write_int63_usize(w, entry.0)?;
    write_bigz_str(w, &entry.1)?;
    Ok(())
}

fn write_sparse_vector<W: Write>(w: &mut W, sparse: &[(usize, String)]) -> Result<()> {
    let entry_d = sparse_entry_descr();
    write_array(w, &entry_d, sparse, |w, entry| write_sparse_entry(w, entry))
}

fn split_geom_edge_lifts(
    cert: &Certificate,
) -> Result<(Vec<Vec<usize>>, Vec<Vec<usize>>)> {
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
                anyhow::anyhow!(
                    "geom_edge_lifts[{v}][{j}] has invalid source simplex {source}"
                )
            })?;

            if source_label.owner != v {
                bail!(
                    "geom_edge_lifts[{v}][{j}] has source simplex {source} owned by {}, expected {v}",
                    source_label.owner
                );
            }

            let target_label = cert.graph.lbl.get(target).ok_or_else(|| {
                anyhow::anyhow!(
                    "geom_edge_lifts[{v}][{j}] has invalid target simplex {target}"
                )
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
    write_bigz_string_array(w, &full_dim.point)?;
    write_bign_str(w, &full_dim.denominator)?;
    write_bigz_string_matrix(w, &full_dim.directions)?;
    write_bigz_string_matrix(w, &full_dim.left_inverse)?;
    Ok(())
}

fn write_root<W: Write>(w: &mut W, root: &Root) -> Result<()> {
    // root := simplex_id * (inverse_incident_map * (basis_vectors * (m_matrix * q_vectors)))
    write_int63_usize(w, root.simplex_id)?;
    write_usize_array(w, &root.inverse_incident_map)?;
    write_bigz_string_matrix(w, &root.basis_vectors)?;
    write_bigz_string_matrix(w, &root.m_matrix)?;

    let sparse_d = sparse_vector_descr();
    write_array(w, &sparse_d, &root.q_vectors, |w, sparse| write_sparse_vector(w, sparse))?;
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
    let (geom_edge_sources, geom_edge_local_targets) =
        split_geom_edge_lifts(cert)?;
    write_usize_matrix(w, &geom_edge_sources)?;
    write_usize_matrix(w, &geom_edge_local_targets)?;

    write_full_dim(w, &cert.full_dim)?;
    write_root(w, &cert.root)?;

    Ok(())
}

pub fn write_certificate_bin<W: std::io::Write>(
    w: &mut W,
    cert: &crate::certificate::Certificate,
) -> anyhow::Result<()> {
    let d = certificate_descr();
    write_descr(w, &d)?;
    write_certificate_value(w, cert)?;
    Ok(())
}

pub fn convert_certificate_json_to_bin<P: AsRef<std::path::Path>, Q: AsRef<std::path::Path>>(
    json_path: P,
    bin_path: Q,
) -> anyhow::Result<()> {
    let json_path_ref = json_path.as_ref();
    let json_path_str = json_path_ref
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("certificate path is not valid UTF-8"))?;

    let cert = read_certificate(json_path_str)
        .with_context(|| format!("failed to read certificate `{}`", json_path_ref.display()))?;

    let file = std::fs::File::create(bin_path.as_ref())
        .with_context(|| format!("failed to create `{}`", bin_path.as_ref().display()))?;
    let mut out = std::io::BufWriter::new(file);

    write_certificate_bin(&mut out, &cert)?;
    out.flush()?;

    Ok(())
}
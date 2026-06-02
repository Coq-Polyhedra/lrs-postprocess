use anyhow::{bail, Context, Result};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::certificate::{read_certificate, AdjacentSimplex, Certificate, Root, Simplex, VertexItem};

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

fn adj_descr() -> Descr {
    // adjacency pointer := item * simplex
    pair(Descr::Int63, Descr::Int63)
}

fn simplex_descr() -> Descr {
    // simplex := indices * adj
    pair(array(Descr::Int63), array(adj_descr()))
}

fn item_descr() -> Descr {
    // item := vertex * (incident * simplices)
    pair(
        array(Descr::BigQ),
        pair(
            array(Descr::Int63),
            array(simplex_descr()),
        ),
    )
}

fn root_descr() -> Descr {
    // root := k0 * (rows * (inverse_rows * same_label_separators))
    pair(
        Descr::Int63,
        pair(
            array(Descr::Int63),
            pair(
                array(array(Descr::BigQ)),
                array(Descr::Int63),
            ),
        ),
    )
}

fn certificate_descr() -> Descr {
    // certificate := items * roots
    pair(array(item_descr()), array(root_descr()))
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

fn write_adj<W: Write>(w: &mut W, adj: &AdjacentSimplex) -> Result<()> {
    write_int63_usize(w, adj.item)?;
    write_int63_usize(w, adj.simplex)?;
    Ok(())
}

fn write_adj_array<W: Write>(w: &mut W, xs: &[AdjacentSimplex]) -> Result<()> {
    let d = adj_descr();
    write_array(w, &d, xs, |w, x| write_adj(w, x))
}

fn write_simplex<W: Write>(w: &mut W, simplex: &Simplex) -> Result<()> {
    // simplex := indices * adj
    write_usize_array(w, &simplex.indices)?;
    write_adj_array(w, &simplex.adj)?;
    Ok(())
}

fn write_simplex_array<W: Write>(w: &mut W, simplices: &[Simplex]) -> Result<()> {
    let d = simplex_descr();
    write_array(w, &d, simplices, |w, simplex| write_simplex(w, simplex))
}

fn write_item<W: Write>(w: &mut W, item: &VertexItem) -> Result<()> {
    // item := vertex * (incident * simplices)
    write_bigq_string_array(w, &item.vertex)?;
    write_usize_array(w, &item.incident)?;
    write_simplex_array(w, &item.simplices)?;
    Ok(())
}

fn write_root<W: Write>(w: &mut W, root: &Root) -> Result<()> {
    // root := k0 * (rows * (inverse_rows * same_label_separators))
    write_int63_usize(w, root.k0)?;
    write_usize_array(w, &root.rows)?;
    write_bigq_string_matrix(w, &root.inverse_rows)?;
    write_usize_array(w, &root.same_label_separators)?;
    Ok(())
}

fn write_certificate_value<W: Write>(w: &mut W, cert: &Certificate) -> Result<()> {
    let item_d = item_descr();
    write_array(w, &item_d, &cert.items, |w, item| write_item(w, item))?;

    let root_d = root_descr();
    match &cert.root {
        None => {
            let empty: Vec<Root> = Vec::new();
            write_array(w, &root_d, &empty, |w, root| write_root(w, root))?;
        }
        Some(root) => {
            write_array(w, &root_d, std::slice::from_ref(root), |w, root| {
                write_root(w, root)
            })?;
        }
    }

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
use anyhow::{bail, Context, Result};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

pub type Q = BigRational;

pub fn parse_q(s: &str) -> Result<Q> {
    let s = s.trim();

    if let Some((num, den)) = s.split_once('/') {
        let n: BigInt = num
            .parse()
            .with_context(|| format!("bad numerator `{num}`"))?;
        let d: BigInt = den
            .parse()
            .with_context(|| format!("bad denominator `{den}`"))?;

        if d.is_zero() {
            bail!("zero denominator in rational `{s}`");
        }

        Ok(BigRational::new(n, d))
    } else {
        let n: BigInt = s
            .parse()
            .with_context(|| format!("bad integer rational `{s}`"))?;
        Ok(BigRational::from_integer(n))
    }
}

pub fn q_to_string(x: &Q) -> String {
    if x.denom().is_one() {
        x.numer().to_string()
    } else {
        format!("{}/{}", x.numer(), x.denom())
    }
}

pub fn dot(a: &[Q], x: &[Q]) -> Q {
    assert_eq!(a.len(), x.len());

    a.iter()
        .zip(x)
        .fold(Q::zero(), |acc, (ai, xi)| acc + ai * xi)
}

pub fn mat_mul(a: &[Vec<Q>], b: &[Vec<Q>]) -> Vec<Vec<Q>> {
    let n = a.len();
    let k = if n == 0 { 0 } else { a[0].len() };
    let m = if b.is_empty() { 0 } else { b[0].len() };

    assert_eq!(b.len(), k);

    let mut c = vec![vec![Q::zero(); m]; n];

    for i in 0..n {
        for r in 0..k {
            for j in 0..m {
                c[i][j] += &a[i][r] * &b[r][j];
            }
        }
    }

    c
}

pub fn identity(n: usize) -> Vec<Vec<Q>> {
    let mut id = vec![vec![Q::zero(); n]; n];

    for i in 0..n {
        id[i][i] = Q::one();
    }

    id
}

pub fn invert_matrix(a: &[Vec<Q>]) -> Result<Vec<Vec<Q>>> {
    let n = a.len();

    if n == 0 {
        bail!("cannot invert empty matrix");
    }

    if a.iter().any(|row| row.len() != n) {
        bail!("matrix is not square");
    }

    let mut aug = vec![vec![Q::zero(); 2 * n]; n];

    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j].clone();
        }
        aug[i][n + i] = Q::one();
    }

    for col in 0..n {
        let pivot = (col..n).find(|&r| !aug[r][col].is_zero());

        let Some(pivot) = pivot else {
            bail!("matrix is singular");
        };

        if pivot != col {
            aug.swap(pivot, col);
        }

        let pivot_val = aug[col][col].clone();

        for j in 0..2 * n {
            aug[col][j] /= pivot_val.clone();
        }

        for r in 0..n {
            if r == col || aug[r][col].is_zero() {
                continue;
            }

            let factor = aug[r][col].clone();

            for j in 0..2 * n {
                let sub = factor.clone() * aug[col][j].clone();
                aug[r][j] -= sub;
            }
        }
    }

    Ok(aug
        .into_iter()
        .map(|row| row[n..].to_vec())
        .collect())
}
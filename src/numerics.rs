use anyhow::{bail, Context, Result};
use rug::{Integer, Rational};

pub type Q = Rational;

pub fn parse_q(s: &str) -> Result<Q> {
    let s = s.trim();

    let (n, d) = if let Some((num, den)) = s.split_once('/') {
        let n = Integer::parse(num.trim())
            .map(Integer::from)
            .with_context(|| format!("bad numerator `{num}`"))?;
        let d = Integer::parse(den.trim())
            .map(Integer::from)
            .with_context(|| format!("bad denominator `{den}`"))?;
        (n, d)
    } else {
        let n = Integer::parse(s)
            .map(Integer::from)
            .with_context(|| format!("bad integer rational `{s}`"))?;
        (n, Integer::from(1))
    };

    if d == 0 {
        bail!("zero denominator in rational `{s}`");
    }

    Ok(Q::from((n, d)))
}

pub fn q_to_string(x: &Q) -> String {
    if x.denom() == &1 {
        x.numer().to_string()
    } else {
        format!("{}/{}", x.numer(), x.denom())
    }
}

pub fn invert_matrix(a: &[Vec<Q>]) -> Result<Vec<Vec<Q>>> {
    let n = a.len();

    if n == 0 {
        bail!("cannot invert empty matrix");
    }

    if a.iter().any(|row| row.len() != n) {
        bail!("matrix is not square");
    }

    let mut aug = vec![vec![Q::new(); 2 * n]; n];

    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j].clone();
        }
        aug[i][n + i] = Q::from(1);
    }

    for col in 0..n {
        let pivot = (col..n).find(|&r| aug[r][col] != 0);

        let Some(pivot) = pivot else {
            bail!("matrix is singular");
        };

        if pivot != col {
            aug.swap(pivot, col);
        }

        let pivot_val = aug[col][col].clone();

        for j in 0..2 * n {
            aug[col][j] /= &pivot_val;
        }

        for r in 0..n {
            if r == col || aug[r][col] == 0 {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_parsing_is_canonical() {
        assert_eq!(q_to_string(&parse_q("-6/-8").unwrap()), "3/4");
        assert_eq!(q_to_string(&parse_q("10/5").unwrap()), "2");
        assert!(parse_q("1/0").is_err());
    }

    #[test]
    fn matrix_inverse_is_exact() {
        let a = vec![
            vec![Q::from(2), Q::from(1)],
            vec![Q::from(1), Q::from(1)],
        ];
        let expected = vec![
            vec![Q::from(1), Q::from(-1)],
            vec![Q::from(-1), Q::from(2)],
        ];
        assert_eq!(invert_matrix(&a).unwrap(), expected);
    }
}

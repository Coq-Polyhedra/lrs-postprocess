use anyhow::{Context, Result};
use rug::Integer;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Certificate {
    pub n_inequalities: usize,
    pub dimension: usize,
    pub inequalities: Vec<Inequality>,
    pub items: Vec<VertexItem>,
    pub graph: SimplexGraph,
    pub neighbors: Vec<Vec<usize>>,
    pub geom_edge_lifts: Vec<Vec<(usize, usize)>>,
    pub full_dim: FullDimCertificate,
    pub root: Root,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Inequality {
    /// Integer coefficients after clearing denominators in this row.
    #[serde(with = "decimal_vec")]
    pub a: Vec<Integer>,

    /// Integer right-hand side after clearing denominators in this row.
    #[serde(with = "decimal")]
    pub b: Integer,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct VertexCoords {
    /// Integer numerators over the common denominator.
    #[serde(with = "decimal_vec")]
    pub num: Vec<Integer>,

    /// Common positive denominator, encoded as BigN in the binary format.
    #[serde(with = "decimal")]
    pub den: Integer,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VertexItem {
    pub incident: Vec<usize>,
    pub vertex: VertexCoords,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphLabel {
    pub simplex: Vec<usize>,
    pub owner: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimplexGraph {
    pub g: Vec<Vec<usize>>,
    pub lbl: Vec<GraphLabel>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FullDimCertificate {
    #[serde(with = "decimal")]
    pub denominator: Integer,

    #[serde(with = "decimal_vec")]
    pub point: Vec<Integer>,

    /// Columns of the integer direction matrix.
    #[serde(with = "decimal_matrix")]
    pub directions: Vec<Vec<Integer>>,

    /// Rows of an integer left inverse up to a nonzero diagonal scaling.
    #[serde(with = "decimal_matrix")]
    pub left_inverse: Vec<Vec<Integer>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Root {
    pub simplex_id: usize,
    pub inverse_incident_map: Vec<usize>,

    #[serde(with = "decimal_matrix")]
    pub basis_vectors: Vec<Vec<Integer>>,

    #[serde(with = "decimal_matrix")]
    pub m_matrix: Vec<Vec<Integer>>,

    #[serde(with = "sparse_decimal_vectors")]
    pub q_vectors: Vec<Vec<(usize, Integer)>>,
}

fn parse_decimal<E: serde::de::Error>(value: &str) -> std::result::Result<Integer, E> {
    Integer::parse(value.trim())
        .map(Integer::from)
        .map_err(|_| E::custom("invalid decimal integer"))
}

mod decimal {
    use super::*;

    pub fn serialize<S>(value: &Integer, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&value.to_string())
    }

    pub fn deserialize<'de, D>(deserializer: D) -> std::result::Result<Integer, D::Error>
    where
        D: Deserializer<'de>,
    {
        parse_decimal(&String::deserialize(deserializer)?)
    }
}

mod decimal_vec {
    use super::*;

    pub fn serialize<S>(values: &[Integer], serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        values
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> std::result::Result<Vec<Integer>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Vec::<String>::deserialize(deserializer)?
            .into_iter()
            .map(|value| parse_decimal(&value))
            .collect()
    }
}

mod decimal_matrix {
    use super::*;

    pub fn serialize<S>(
        values: &[Vec<Integer>],
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        values
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect::<Vec<_>>())
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D>(
        deserializer: D,
    ) -> std::result::Result<Vec<Vec<Integer>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Vec::<Vec<String>>::deserialize(deserializer)?
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|value| parse_decimal(&value))
                    .collect::<std::result::Result<Vec<_>, D::Error>>()
            })
            .collect()
    }
}

mod sparse_decimal_vectors {
    use super::*;

    pub fn serialize<S>(
        values: &[Vec<(usize, Integer)>],
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        values
            .iter()
            .map(|row| {
                row.iter()
                    .map(|(index, value)| (*index, value.to_string()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D>(
        deserializer: D,
    ) -> std::result::Result<Vec<Vec<(usize, Integer)>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Vec::<Vec<(usize, String)>>::deserialize(deserializer)?
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|(index, value)| Ok((index, parse_decimal::<D::Error>(&value)?)))
                    .collect::<std::result::Result<Vec<_>, D::Error>>()
            })
            .collect()
    }
}

/// Reads JSON while decoding decimal strings directly into GMP integers.
pub fn read_certificate(path: &str) -> Result<Certificate> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read certificate `{path}`"))?;

    serde_json::from_str(&text)
        .with_context(|| format!("failed to parse certificate JSON and integers `{path}`"))
}

pub fn certificate_to_string(cert: &Certificate, pretty: bool) -> Result<String> {
    if pretty {
        Ok(serde_json::to_string_pretty(cert)?)
    } else {
        Ok(serde_json::to_string(cert)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_fields_round_trip_as_decimal_strings() {
        let coords = VertexCoords {
            num: vec![Integer::from(-7), Integer::from(1) << 200usize],
            den: Integer::from(11),
        };
        let json = serde_json::to_string(&coords).unwrap();
        assert!(json.contains("\"-7\""));
        let large = Integer::from(1) << 200usize;
        assert!(json.contains(&format!("\"{large}\"")));
        assert_eq!(serde_json::from_str::<VertexCoords>(&json).unwrap(), coords);
    }
}

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Certificate {
    /// Inequalities A x <= b used by the certificate.
    pub inequalities: Vec<Inequality>,

    pub items: Vec<VertexItem>,
    pub graph: SimplexGraph,
    pub root: Root,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Inequality {
    /// Integer coefficients of a_i after clearing denominators in this row.
    pub a: Vec<String>,

    /// Integer right-hand side b_i after clearing denominators in this row.
    pub b: String,
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct VertexCoords {
    /// Integer numerators after clearing denominators.
    pub num: Vec<String>,

    /// Common positive denominator, encoded as BigN in the binary format.
    pub den: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VertexItem {
    /// Vertex coordinates represented as num / den, componentwise.
    pub vertex: VertexCoords,
    pub incident: Vec<usize>,
    pub simplices: Vec<ItemSimplex>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ItemSimplex {
    /// Local indices into `incident`.
    pub indices: Vec<usize>,

    /// Node of the global simplex adjacency graph.
    pub node: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct LocalSimplexRef {
    pub item: usize,
    pub simplex: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphLabel {
    /// Global row indices of the simplex. This list must be sorted and have length d.
    pub simplex: Vec<usize>,

    /// Inverse map from this graph node to the corresponding local simplex.
    pub owner: LocalSimplexRef,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimplexGraph {
    /// Directed adjacency lists. The checker verifies that the directed graph is symmetric.
    pub g: Vec<Vec<usize>>,

    /// Labels of graph nodes. Each entry contains the global simplex and its local owner.
    pub lbl: Vec<GraphLabel>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Root {
    /// Root vertex item.
    pub k0: usize,

    /// Global row indices of the root simplex, namely items[k0].simplices[0].
    pub rows: Vec<usize>,

    /// Inverse of the root basis matrix, provided as row vectors.
    ///
    /// If rows are j_1,...,j_d and R has columns A_{j_r}^T,
    /// then inverse_rows[a] dot A_{j_b} = delta_{ab}.
    pub inverse_rows: Vec<Vec<String>>,

    /// For each simplex of items[k0].simplices except simplex 0, in increasing
    /// simplex order, this gives the row of inverse_rows used as separator.
    pub same_label_separators: Vec<usize>,
}

pub fn write_certificate(path: &str, cert: &Certificate, pretty: bool) -> Result<()> {
    let text = certificate_to_string(cert, pretty)?;
    fs::write(path, text).with_context(|| format!("failed to write certificate `{path}`"))?;
    Ok(())
}

pub fn read_certificate(path: &str) -> Result<Certificate> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read certificate `{path}`"))?;

    serde_json::from_str(&text)
        .with_context(|| format!("failed to parse certificate JSON `{path}`"))
}

pub fn certificate_to_string(cert: &Certificate, pretty: bool) -> Result<String> {
    if pretty {
        Ok(serde_json::to_string_pretty(cert)?)
    } else {
        Ok(serde_json::to_string(cert)?)
    }
}

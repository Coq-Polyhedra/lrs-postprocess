use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Certificate {
    /// Inequalities A x <= b used by the certificate.
    pub inequalities: Vec<Inequality>,

    /// Candidate vertices/items.  The simplices are no longer stored locally
    /// here: all simplices live globally in `graph.lbl`.
    pub items: Vec<VertexItem>,

    pub graph: SimplexGraph,

    /// Candidate edge-neighbors for each item.
    ///
    /// `neighbors[v]` is a strictly sorted list of item indices `w`.
    /// The postprocessor derives it from cross-label ridge adjacencies
    /// in the simplex graph. The checker verifies symmetry and exact
    /// agreement with the graph-derived neighbor relation.
    pub neighbors: Vec<Vec<usize>>,

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
    /// Incident inequality indices.
    ///
    /// This field is intentionally first in the binary encoding of items.
    /// The item array is sorted lexicographically by this list.
    pub incident: Vec<usize>,

    /// Vertex coordinates represented as num / den, componentwise.
    pub vertex: VertexCoords,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphLabel {
    /// Global row indices of the simplex. This list must be sorted and have length d.
    pub simplex: Vec<usize>,

    /// Owner item/vertex of this global simplex.
    pub owner: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimplexGraph {
    /// Position-indexed directed adjacency lists.
    ///
    /// If `lbl[k].simplex = sigma`, then `g[k][j]` is the node adjacent
    /// to `k` through the ridge `sigma \ {sigma[j]}`.
    /// The checker verifies the reciprocal condition.
    pub g: Vec<Vec<usize>>,

    /// Labels of graph nodes. Each entry contains the global simplex and its owner item.
    pub lbl: Vec<GraphLabel>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Root {
    /// Root vertex item.
    pub k0: usize,

    /// Global row indices of the root simplex.
    pub rows: Vec<usize>,

    /// Inverse of the root basis matrix, provided as row vectors.
    ///
    /// If rows are j_1,...,j_d and R has columns A_{j_r}^T,
    /// then inverse_rows[a] dot A_{j_b} = delta_{ab}.
    pub inverse_rows: Vec<Vec<String>>,

    /// For each other global simplex whose owner is k0, in graph-label order,
    /// this gives the row of inverse_rows used as separator.
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

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Certificate {
    /// Explicit number of inequalities.
    pub n_inequalities: usize,

    /// Explicit ambient dimension.
    pub dimension: usize,

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

    /// For every directed geometric edge `v -> neighbors[v][j]`, stores
    /// one oriented simplex-graph edge `(s, t)` with owner(s)=v and
    /// owner(t)=neighbors[v][j]. This array is positionally parallel to
    /// `neighbors`.
    pub geom_edge_lifts: Vec<Vec<(usize, usize)>>,

    /// Integer certificate that the polyhedron is full-dimensional.
    pub full_dim: FullDimCertificate,

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
pub struct FullDimCertificate {
    /// Common positive denominator q.
    pub denominator: String,

    /// Integer numerator p of the feasible base point x0 = p / q.
    pub point: Vec<String>,

    /// The d integer direction numerators r^j. Each outer entry is one
    /// column of R, so directions[j][k] = R_{k,j}.
    pub directions: Vec<Vec<String>>,

    /// Integer d x d matrix U, stored by rows, such that U R is diagonal
    /// with nonzero diagonal entries.
    pub left_inverse: Vec<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Root {
    /// Graph-node index of the distinguished root simplex sigma*.
    pub simplex_id: usize,

    /// Inverse map to the root owner incident list.
    ///
    /// This has length `n_inequalities`. If the root owner has incident list I,
    /// then `inverse_incident_map[i]` is the position of `i` in I when `i in I`,
    /// and is `n_inequalities` otherwise.
    pub inverse_incident_map: Vec<usize>,

    /// Integer vectors f_1, ..., f_d.
    ///
    /// The outer array is indexed by j, and `basis_vectors[j]` is the
    /// coordinate vector f_j in the ambient dimension d.
    pub basis_vectors: Vec<Vec<String>>,

    /// Matrix M indexed by local active-inequality position and basis-vector index.
    ///
    /// If sigma* is owned by v and I_v = items[v].incident, then
    /// `m_matrix[p][j] = a_{I_v[p]}^T f_j`, where the integer row a_i is
    /// the numerator row stored in `inequalities[i].a`.
    pub m_matrix: Vec<Vec<String>>,

    /// Sparse nonnegative certificates for same-owner simplices distinct from sigma*.
    ///
    /// Entries are listed in graph-label order over labels with the same owner as
    /// sigma*, excluding sigma* itself. A sparse vector is encoded as sorted pairs
    /// `(coordinate, positive_integer_value)`.
    pub q_vectors: Vec<Vec<(usize, String)>>,
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

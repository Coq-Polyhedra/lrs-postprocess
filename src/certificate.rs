use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Certificate {
    pub items: Vec<VertexItem>,
    pub root: Option<Root>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VertexItem {
    pub vertex: Vec<String>,
    pub incident: Vec<usize>,
    pub simplices: Vec<Vec<usize>>,
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
use anyhow::{anyhow, Context, Result};
use rug::{Complete, Integer};
use std::time::Instant;

use crate::certificate::{read_certificate, Certificate, Inequality, VertexItem};
use crate::postprocess::{parse_lrs_hrep, HRep};

fn has_matrix_shape<T>(matrix: &[Vec<T>], rows: usize, columns: usize) -> bool {
    matrix.len() == rows && matrix.iter().all(|row| row.len() == columns)
}

fn dot(x: &[Integer], y: &[Integer]) -> Option<Integer> {
    if x.len() != y.len() {
        return None;
    }
    let mut result = Integer::new();
    for (xi, yi) in x.iter().zip(y) {
        // `xi * yi` is a rug incomplete computation. AddAssign consumes it
        // through GMP's fused add-multiply path without allocating an
        // intermediate Integer for the product.
        result += xi * yi;
    }
    Some(result)
}

fn sparse_dot(weight: &[(usize, Integer)], x: &[Integer]) -> Option<Integer> {
    let mut result = Integer::new();
    for (i, value) in weight {
        result += value * x.get(*i)?;
    }
    Some(result)
}

fn strictly_sorted(xs: &[usize]) -> bool {
    xs.windows(2).all(|w| w[0] < w[1])
}

fn flag_is_well_formed(item: &VertexItem, dimension: usize, item_count: usize) -> bool {
    let (local_inequalities, witnesses) = &item.flag;
    local_inequalities.len() == dimension
        && witnesses.len() == dimension
        && local_inequalities
            .iter()
            .all(|&local_index| local_index < item.incident.len())
        && witnesses.iter().all(|&witness| witness < item_count)
}

fn complete_flag_check(items: &[VertexItem], item: &VertexItem, dimension: usize) -> bool {
    if !flag_is_well_formed(item, dimension, items.len()) {
        return false;
    }

    let (local_inequalities, witnesses) = &item.flag;
    for k in 0..dimension {
        let Some(&global_inequality) = local_inequalities
            .get(k)
            .and_then(|&local_index| item.incident.get(local_index))
        else {
            return false;
        };
        let Some(witness) = witnesses
            .get(k)
            .and_then(|&witness_index| items.get(witness_index))
        else {
            return false;
        };

        // The kth inequality must cut the preceding face strictly.
        if witness.incident.binary_search(&global_inequality).is_ok() {
            return false;
        }

        // The witness belongs to the face defined by the preceding prefix.
        for &previous_local_index in &local_inequalities[..k] {
            let Some(&previous_global_inequality) = item.incident.get(previous_local_index) else {
                return false;
            };
            if witness
                .incident
                .binary_search(&previous_global_inequality)
                .is_err()
            {
                return false;
            }
        }
    }

    true
}

fn sorted_subset(xs: &[usize], ys: &[usize]) -> bool {
    let mut i = 0;
    let mut j = 0;
    while i < xs.len() && j < ys.len() {
        if xs[i] == ys[j] {
            i += 1;
            j += 1;
        } else if xs[i] > ys[j] {
            j += 1;
        } else {
            return false;
        }
    }
    i == xs.len()
}

fn sorted_difference(xs: &[usize], ys: &[usize]) -> Vec<usize> {
    let mut i = 0;
    let mut j = 0;
    let mut result = Vec::new();
    while i < xs.len() && j < ys.len() {
        if xs[i] == ys[j] {
            i += 1;
            j += 1;
        } else if xs[i] < ys[j] {
            result.push(xs[i]);
            i += 1;
        } else {
            j += 1;
        }
    }
    result.extend_from_slice(&xs[i..]);
    result
}

fn is_undirected(graph: &[Vec<usize>]) -> bool {
    graph.iter().enumerate().all(|(i, neighbors)| {
        neighbors.iter().all(|&j| {
            graph
                .get(j)
                .is_some_and(|reverse_neighbors| reverse_neighbors.contains(&i))
        })
    })
}

fn integer_inequality_from_hrep(h: &HRep, i: usize) -> Option<(Vec<Integer>, Integer)> {
    let a = h.a.get(i)?;
    let b = h.b.get(i)?;
    let lcm = a
        .iter()
        .chain(std::iter::once(b))
        .fold(Integer::from(1), |acc, value| {
            acc.lcm_ref(value.denom()).complete()
        });
    let clear = |value: &rug::Rational| {
        let mut integer = value.numer().clone();
        integer *= (&lcm / value.denom()).complete();
        integer
    };
    Some((a.iter().map(|value| clear(value)).collect(), clear(b)))
}

// Rust-specific check: the Rocq checker receives the integer inequalities
// directly and therefore has no external H-representation to compare against.
fn check_ine_file(h: &HRep, cert: &Certificate) -> bool {
    cert.n_inequalities == h.a.len()
        && cert.dimension == h.d
        && cert.inequalities.len() == h.a.len()
        && cert.inequalities.iter().enumerate().all(|(i, inequality)| {
            let Some((expected_normal, expected_bound)) =
                integer_inequality_from_hrep(h, i)
            else {
                return false;
            };
            inequality.a == expected_normal && inequality.b == expected_bound
        })
}

mod rocq {
    #![allow(non_snake_case)]

    use super::*;

    pub fn areInequalitiesWellFormed(cert: &Certificate) -> bool {
        cert.inequalities.len() == cert.n_inequalities
            && cert.inequalities.iter().all(|inequality| {
                inequality.a.len() == cert.dimension
                    && inequality.a.iter().any(|x| x != &0)
            })
    }

    pub fn arePointsWellFormed(cert: &Certificate) -> bool {
        cert.items.iter().all(|item| {
            item.vertex.den > 0 && item.vertex.num.len() == cert.dimension
        })
    }

    pub fn areActiveSetsWellFormed(cert: &Certificate) -> bool {
        cert.items.iter().all(|item| {
            strictly_sorted(&item.incident)
                && item
                    .incident
                    .iter()
                    .all(|&i| i < cert.n_inequalities)
        })
    }

    pub fn areFlagsWellFormed(cert: &Certificate) -> bool {
        cert.items
            .iter()
            .all(|item| flag_is_well_formed(item, cert.dimension, cert.items.len()))
    }

    pub fn areVerticesWellFormed(cert: &Certificate) -> bool {
        arePointsWellFormed(cert) && areActiveSetsWellFormed(cert) && areFlagsWellFormed(cert)
    }

    pub fn isGraphWellFormed(cert: &Certificate) -> bool {
        let graph = &cert.graph.g;
        let number_of_facets = cert.graph.lbl.len();
        graph.len() == number_of_facets
            && graph
                .iter()
                .flatten()
                .all(|&facet| facet < number_of_facets)
            && graph
                .iter()
                .enumerate()
                .all(|(i, neighbors)| !neighbors.contains(&i))
            && is_undirected(graph)
    }

    pub fn areDescriptionsWellFormed(cert: &Certificate) -> bool {
        cert.graph.lbl.iter().all(|facet| {
            facet.simplex.len() == cert.dimension
                && strictly_sorted(&facet.simplex)
                && facet
                    .simplex
                    .iter()
                    .all(|&i| i < cert.n_inequalities)
        })
    }

    pub fn isMappingWellFormed(cert: &Certificate) -> bool {
        cert.graph
            .lbl
            .iter()
            .all(|facet| facet.owner < cert.items.len())
    }

    pub fn areFacetsWellFormed(cert: &Certificate) -> bool {
        areDescriptionsWellFormed(cert) && isMappingWellFormed(cert)
    }

    pub fn isGeomGraphWellFormed(cert: &Certificate) -> bool {
        let graph = &cert.neighbors;
        let number_of_vertices = cert.items.len();
        graph.len() == number_of_vertices
            && graph.iter().all(|neighbors| {
                strictly_sorted(neighbors)
                    && neighbors.iter().all(|&vertex| vertex < number_of_vertices)
            })
            && graph
                .iter()
                .enumerate()
                .all(|(i, neighbors)| !neighbors.contains(&i))
            && is_undirected(graph)
    }

    // A Rust lift stores (source facet, target facet), whereas Rocq stores
    // (source facet, local position of the target in the source adjacency row).
    pub fn areGeomEdgeSourcesWellFormed(cert: &Certificate) -> bool {
        let number_of_vertices = cert.items.len();
        let number_of_facets = cert.graph.lbl.len();
        cert.geom_edge_lifts.len() == number_of_vertices
            && cert
                .geom_edge_lifts
                .iter()
                .enumerate()
                .all(|(i, lifts)| {
                    cert.neighbors.get(i).is_some_and(|neighbors| {
                        lifts.len() == neighbors.len()
                            && lifts.iter().all(|&(source, _)| source < number_of_facets)
                    })
                })
    }

    pub fn areGeomEdgeLocalTargetsWellFormed(cert: &Certificate) -> bool {
        let number_of_vertices = cert.items.len();
        let number_of_facets = cert.graph.lbl.len();
        cert.geom_edge_lifts.len() == number_of_vertices
            && cert
                .geom_edge_lifts
                .iter()
                .enumerate()
                .all(|(i, lifts)| {
                    cert.neighbors.get(i).is_some_and(|neighbors| {
                        lifts.len() == neighbors.len()
                            && lifts.iter().all(|&(source, target)| {
                                target < number_of_facets
                                    && cert
                                        .graph
                                        .g
                                        .get(source)
                                        .is_some_and(|adjacency| adjacency.contains(&target))
                            })
                    })
                })
    }

    pub fn isFullDimPointWellFormed(cert: &Certificate) -> bool {
        cert.full_dim.denominator > 0
            && cert.full_dim.point.len() == cert.dimension
    }

    pub fn isFullDimDirWellFormed(cert: &Certificate) -> bool {
        has_matrix_shape(
            &cert.full_dim.directions,
            cert.dimension,
            cert.dimension,
        )
    }

    pub fn isFullDimInverseWellFormed(cert: &Certificate) -> bool {
        has_matrix_shape(
            &cert.full_dim.left_inverse,
            cert.dimension,
            cert.dimension,
        )
    }

    pub fn isFullDimWellFormed(cert: &Certificate) -> bool {
        isFullDimPointWellFormed(cert)
            && isFullDimDirWellFormed(cert)
            && isFullDimInverseWellFormed(cert)
    }

    pub fn isSimplexIndexWellFormed(cert: &Certificate) -> bool {
        cert.root.simplex_id < cert.graph.lbl.len()
    }

    pub fn isActiveInverseWellFormed(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        let inverse = &cert.root.inverse_incident_map;
        let active = &root_vertex.incident;
        let sentinel = cert.n_inequalities;

        inverse.len() == cert.n_inequalities
            && active
                .iter()
                .enumerate()
                .all(|(i, &row)| inverse.get(row) == Some(&i))
            && inverse.iter().enumerate().all(|(row, &value)| {
                active.binary_search(&row).is_ok() || value == sentinel
            })
    }

    pub fn areWitnessesWellFormed(cert: &Certificate) -> bool {
        has_matrix_shape(
            &cert.root.basis_vectors,
            cert.dimension,
            cert.dimension,
        )
    }

    pub fn areScalarProductsWellFormed(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        has_matrix_shape(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        )
    }

    pub fn isSparseVectorWellFormed(dimension: usize, weight: &[(usize, Integer)]) -> bool {
        weight
            .iter()
            .all(|(i, value)| *i < dimension && value != &0)
            && weight.windows(2).all(|w| w[0].0 < w[1].0)
    }

    pub fn areWeightsWellFormed(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let expected = cert
            .graph
            .lbl
            .iter()
            .enumerate()
            .filter(|(i, facet)| {
                *i != cert.root.simplex_id && facet.owner == root_facet.owner
            })
            .count();
        cert.root.q_vectors.len() == expected
            && cert
                .root
                .q_vectors
                .iter()
                .all(|weight| isSparseVectorWellFormed(cert.dimension, weight))
    }

    pub fn isRootWellFormed(cert: &Certificate) -> bool {
        isSimplexIndexWellFormed(cert)
            && isActiveInverseWellFormed(cert)
            && areWitnessesWellFormed(cert)
            && areScalarProductsWellFormed(cert)
            && areWeightsWellFormed(cert)
    }

    pub fn areActiveSetsUnique(cert: &Certificate) -> bool {
        cert.items
            .windows(2)
            .all(|items| items[0].incident < items[1].incident)
    }

    pub fn areFacetsUnique(cert: &Certificate) -> bool {
        cert.graph
            .lbl
            .windows(2)
            .all(|facets| facets[0].simplex < facets[1].simplex)
    }

    pub fn check_ineqs(
        inequalities: &[Inequality],
        active_set: &[usize],
        numerators: &[Integer],
        denominator: &Integer,
    ) -> bool {
        let mut active_position = 0;

        for (i, inequality) in inequalities.iter().enumerate() {
            let Some(lhs) = dot(&inequality.a, numerators) else {
                return false;
            };
            let rhs = (&inequality.b * denominator).complete();

            if active_set.get(active_position) == Some(&i) {
                if lhs != rhs {
                    return false;
                }
                active_position += 1;
            } else if lhs >= rhs {
                return false;
            }
        }

        active_position == active_set.len()
    }

    pub fn feasibility_check(cert: &Certificate) -> bool {
        cert.items.iter().all(|item| {
            check_ineqs(
                &cert.inequalities,
                &item.incident,
                &item.vertex.num,
                &item.vertex.den,
            )
        })
    }

    pub fn vertexhood_check(cert: &Certificate) -> bool {
        cert.items
            .iter()
            .all(|item| complete_flag_check(&cert.items, item, cert.dimension))
    }

    pub fn isRidgeInFacet(facet1: &[usize], facet2: &[usize], value: usize) -> bool {
        let difference = sorted_difference(facet1, facet2);
        difference.len() == 1 && difference[0] == value
    }

    pub fn graph_check(cert: &Certificate) -> bool {
        cert.graph.g.iter().enumerate().all(|(i, neighbors)| {
            neighbors.len() <= cert.dimension
                && neighbors.iter().enumerate().all(|(j, &neighbor)| {
                    let Some(facet) = cert.graph.lbl.get(i) else {
                        return false;
                    };
                    let Some(other_facet) = cert.graph.lbl.get(neighbor) else {
                        return false;
                    };
                    facet.simplex.get(j).is_some_and(|&value| {
                        isRidgeInFacet(&facet.simplex, &other_facet.simplex, value)
                    })
                })
        })
    }

    pub fn mapping_check(cert: &Certificate) -> bool {
        cert.graph.lbl.iter().all(|facet| {
            cert.items
                .get(facet.owner)
                .is_some_and(|vertex| sorted_subset(&facet.simplex, &vertex.incident))
        })
    }

    pub fn scalarProducts_check(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        if !has_matrix_shape(
            &cert.root.basis_vectors,
            cert.dimension,
            cert.dimension,
        ) {
            return false;
        }
        if !has_matrix_shape(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        ) {
            return false;
        }

        root_vertex
            .incident
            .iter()
            .enumerate()
            .all(|(i, &row)| {
                let Some(inequality) = cert.inequalities.get(row) else {
                    return false;
                };
                (0..cert.dimension).all(|j| {
                    dot(&inequality.a, &cert.root.basis_vectors[j])
                        .is_some_and(|expected| cert.root.m_matrix[i][j] == expected)
                })
            })
    }

    pub fn inversibility_check(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        if !has_matrix_shape(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        ) {
            return false;
        }

        root_facet
            .simplex
            .iter()
            .enumerate()
            .all(|(i, &row)| {
                let Some(&active_position) = cert.root.inverse_incident_map.get(row) else {
                    return false;
                };
                let Some(product_row) = cert.root.m_matrix.get(active_position) else {
                    return false;
                };
                (0..cert.dimension).all(|j| {
                    if i == j {
                        product_row[j] > 0
                    } else {
                        product_row[j] == 0
                    }
                })
            })
    }

    pub fn isSparseVectorPositive(weight: &[(usize, Integer)]) -> bool {
        !weight.is_empty()
            && weight
                .iter()
                .all(|(_, value)| value >= &0)
    }

    pub fn separability_check(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        if !has_matrix_shape(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        ) {
            return false;
        }

        let same_owner_nonroot = cert
            .graph
            .lbl
            .iter()
            .enumerate()
            .filter(|(i, facet)| {
                *i != cert.root.simplex_id && facet.owner == root_facet.owner
            })
            .map(|(_, facet)| facet)
            .collect::<Vec<_>>();

        cert.root.q_vectors.len() == same_owner_nonroot.len()
            && same_owner_nonroot
                .iter()
                .zip(&cert.root.q_vectors)
                .all(|(facet, weight)| {
                    isSparseVectorPositive(weight)
                        && facet.simplex.iter().all(|&row| {
                            let Some(&active_position) =
                                cert.root.inverse_incident_map.get(row)
                            else {
                                return false;
                            };
                            cert.root.m_matrix.get(active_position).is_some_and(|product_row| {
                                sparse_dot(weight, product_row)
                                    .is_some_and(|value| value <= 0)
                            })
                        })
                })
    }

    pub fn ridge_graph_image_check(cert: &Certificate) -> bool {
        let graph_edges_are_mapped =
            cert.graph.g.iter().enumerate().all(|(source_facet, neighbors)| {
                neighbors.iter().all(|&target_facet| {
                    let Some(source) = cert.graph.lbl.get(source_facet) else {
                        return false;
                    };
                    let Some(target) = cert.graph.lbl.get(target_facet) else {
                        return false;
                    };
                    source.owner == target.owner
                        || cert
                            .neighbors
                            .get(source.owner)
                            .is_some_and(|neighbors| {
                                neighbors.binary_search(&target.owner).is_ok()
                            })
                })
            });

        let geom_edges_are_images =
            cert.neighbors.iter().enumerate().all(|(source_vertex, neighbors)| {
                let Some(lifts) = cert.geom_edge_lifts.get(source_vertex) else {
                    return false;
                };
                lifts.len() == neighbors.len()
                    && neighbors.iter().zip(lifts).all(
                        |(&target_vertex, &(source_facet, target_facet))| {
                            cert.graph
                                .lbl
                                .get(source_facet)
                                .zip(cert.graph.lbl.get(target_facet))
                                .is_some_and(|(source, target)| {
                                    source.owner == source_vertex
                                        && target.owner == target_vertex
                                        && cert
                                            .graph
                                            .g
                                            .get(source_facet)
                                            .is_some_and(|adjacency| {
                                                adjacency.contains(&target_facet)
                                            })
                                })
                        },
                    )
            });

        graph_edges_are_mapped && geom_edges_are_images
    }

    pub fn full_dim_feasibility_check(cert: &Certificate) -> bool {
        if !has_matrix_shape(
            &cert.full_dim.directions,
            cert.dimension,
            cert.dimension,
        ) {
            return false;
        }

        cert.inequalities.iter().all(|inequality| {
            let Some(base) = dot(&inequality.a, &cert.full_dim.point) else {
                return false;
            };
            let rhs = (&inequality.b * &cert.full_dim.denominator).complete();
            base <= rhs
                && cert.full_dim.directions.iter().all(|direction| {
                    dot(&inequality.a, direction).is_some_and(|increment| {
                        let mut value = base.clone();
                        value += increment;
                        value <= rhs
                    })
                })
        })
    }

    pub fn full_dim_inverse_check(cert: &Certificate) -> bool {
        if !has_matrix_shape(
            &cert.full_dim.directions,
            cert.dimension,
            cert.dimension,
        ) {
            return false;
        }
        if !has_matrix_shape(
            &cert.full_dim.left_inverse,
            cert.dimension,
            cert.dimension,
        ) {
            return false;
        }

        (0..cert.dimension).all(|i| {
            (0..cert.dimension).all(|j| {
                dot(
                    &cert.full_dim.directions[i],
                    &cert.full_dim.left_inverse[j],
                )
                .is_some_and(|value| {
                    if i == j {
                        value != 0
                    } else {
                        value == 0
                    }
                })
            })
        })
    }

    pub fn full_dim_check(cert: &Certificate) -> bool {
        full_dim_feasibility_check(cert) && full_dim_inverse_check(cert)
    }

    pub fn well_formedness_check(cert: &Certificate) -> bool {
        areInequalitiesWellFormed(cert)
            && areVerticesWellFormed(cert)
            && isGraphWellFormed(cert)
            && areFacetsWellFormed(cert)
            && isGeomGraphWellFormed(cert)
            && areGeomEdgeSourcesWellFormed(cert)
            && areGeomEdgeLocalTargetsWellFormed(cert)
            && isFullDimWellFormed(cert)
            && isRootWellFormed(cert)
    }

    pub fn uniqueness_check(cert: &Certificate) -> bool {
        areActiveSetsUnique(cert) && areFacetsUnique(cert)
    }

    pub fn root_check(cert: &Certificate) -> bool {
        scalarProducts_check(cert)
            && inversibility_check(cert)
            && separability_check(cert)
    }

}

fn timed_bool(name: &str, check: impl FnOnce() -> bool) -> bool {
    let start = Instant::now();
    let result = check();
    eprintln!("{name}: {:.6} s", start.elapsed().as_secs_f64());
    result
}

fn check_in_memory_certificate(h: &HRep, cert: &Certificate) -> Result<()> {
    let inequality_file_check =
        timed_bool("Inequality file check", || check_ine_file(h, cert));
    let well_formedness =
        timed_bool("Well-formedness check", || rocq::well_formedness_check(cert));
    let uniqueness = timed_bool("Uniqueness check", || rocq::uniqueness_check(cert));
    let feasibility = timed_bool("Feasibility check", || rocq::feasibility_check(cert));
    let vertexhood = timed_bool("Vertexhood flag check", || rocq::vertexhood_check(cert));
    let graph = timed_bool("Graph check", || rocq::graph_check(cert));
    let mapping = timed_bool("Mapping check", || rocq::mapping_check(cert));
    let root = timed_bool("Root check", || rocq::root_check(cert));
    let ridge_graph_image = timed_bool("Ridge graph image check", || {
        rocq::ridge_graph_image_check(cert)
    });
    let full_dimension =
        timed_bool("Full dimension check", || rocq::full_dim_check(cert));

    let accepted = inequality_file_check
        && well_formedness
        && uniqueness
        && feasibility
        && vertexhood
        && graph
        && mapping
        && root
        && ridge_graph_image
        && full_dimension;

    if accepted {
        return Ok(());
    }

    let failed = [
        ("inequality-file", inequality_file_check),
        ("well-formedness", well_formedness),
        ("uniqueness", uniqueness),
        ("feasibility", feasibility),
        ("vertexhood flag", vertexhood),
        ("graph", graph),
        ("mapping", mapping),
        ("root", root),
        ("ridge graph image", ridge_graph_image),
        ("full dimension", full_dimension),
    ]
    .into_iter()
    .filter_map(|(name, passed)| (!passed).then_some(name))
    .collect::<Vec<_>>()
    .join(", ");

    Err(anyhow!("certificate rejected; failed checks: {failed}"))
}

/// Checks a certificate immediately after generation, without a JSON
/// serialization/write/read/parse round trip.
pub fn check_generated_certificate(h: &HRep, cert: &Certificate) -> Result<()> {
    check_in_memory_certificate(h, cert)
}

pub fn check_certificate(ine_path: &str, certificate_path: &str) -> Result<()> {
    let start = Instant::now();
    let h = parse_lrs_hrep(ine_path).context("failed to parse H-representation")?;
    eprintln!(
        "Parse H-representation: {:.6} s",
        start.elapsed().as_secs_f64()
    );

    let start = Instant::now();
    let cert = read_certificate(certificate_path)?;
    eprintln!("Read certificate: {:.6} s", start.elapsed().as_secs_f64());

    check_in_memory_certificate(&h, &cert)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flag_item(
        incident: Vec<usize>,
        local_inequalities: Vec<usize>,
        witnesses: Vec<usize>,
    ) -> VertexItem {
        VertexItem {
            incident,
            vertex: crate::certificate::VertexCoords {
                num: Vec::new(),
                den: Integer::from(1),
            },
            flag: (local_inequalities, witnesses),
        }
    }

    fn graph_checks_certificate() -> Certificate {
        let items = vec![
            flag_item(vec![0, 1, 2], Vec::new(), Vec::new()),
            flag_item(vec![1, 2, 3], Vec::new(), Vec::new()),
            flag_item(vec![0, 2, 4], Vec::new(), Vec::new()),
        ];
        let graph = crate::certificate::SimplexGraph {
            g: vec![vec![1, 2], vec![0], vec![0]],
            lbl: (0..3)
                .map(|owner| crate::certificate::GraphLabel {
                    simplex: Vec::new(),
                    owner,
                })
                .collect(),
        };

        Certificate {
            n_inequalities: 5,
            dimension: 0,
            inequalities: Vec::new(),
            items,
            graph,
            neighbors: vec![vec![1, 2], vec![0], vec![0]],
            geom_edge_lifts: vec![vec![(0, 1), (0, 2)], vec![(1, 0)], vec![(2, 0)]],
            full_dim: crate::certificate::FullDimCertificate {
                denominator: Integer::from(1),
                point: Vec::new(),
                directions: Vec::new(),
                left_inverse: Vec::new(),
            },
            root: crate::certificate::Root {
                simplex_id: 0,
                inverse_incident_map: Vec::new(),
                basis_vectors: Vec::new(),
                m_matrix: Vec::new(),
                q_vectors: Vec::new(),
            },
        }
    }

    #[test]
    fn rug_dot_is_exact_for_large_signed_values() {
        let large = Integer::from(1) << 200usize;
        let x = vec![large.clone(), -large.clone(), Integer::from(7)];
        let y = vec![Integer::from(3), Integer::from(5), Integer::from(-11)];
        let mut expected = -(large << 1usize);
        expected -= 77;
        assert_eq!(dot(&x, &y), Some(expected));
    }

    #[test]
    fn rug_dot_rejects_different_lengths() {
        assert_eq!(dot(&[Integer::from(1)], &[]), None);
    }

    #[test]
    fn complete_flag_resolves_inequalities_locally() {
        // The flag [1, 0] denotes global inequalities [5, 2] in the target's
        // active set. Witness 1 avoids 5; witness 2 saturates 5 and avoids 2.
        let items = vec![
            flag_item(vec![2, 5], vec![1, 0], vec![1, 2]),
            flag_item(vec![2], Vec::new(), Vec::new()),
            flag_item(vec![5], Vec::new(), Vec::new()),
        ];

        assert!(flag_is_well_formed(&items[0], 2, items.len()));
        assert!(complete_flag_check(&items, &items[0], 2));
    }

    #[test]
    fn complete_flag_rejects_bad_indices() {
        let bad_local = flag_item(vec![2, 5], vec![2, 0], vec![1, 2]);
        let bad_witness = flag_item(vec![2, 5], vec![1, 0], vec![1, 3]);

        assert!(!flag_is_well_formed(&bad_local, 2, 3));
        assert!(!flag_is_well_formed(&bad_witness, 2, 3));
    }

    #[test]
    fn complete_flag_rejects_missing_prefix_incidence() {
        let items = vec![
            flag_item(vec![2, 5], vec![1, 0], vec![1, 2]),
            flag_item(vec![2], Vec::new(), Vec::new()),
            flag_item(Vec::new(), Vec::new(), Vec::new()),
        ];

        assert!(!complete_flag_check(&items, &items[0], 2));
    }

    #[test]
    fn complete_flag_rejects_nonstrict_cut() {
        let items = vec![
            flag_item(vec![2, 5], vec![1, 0], vec![1, 2]),
            flag_item(vec![2, 5], Vec::new(), Vec::new()),
            flag_item(vec![5], Vec::new(), Vec::new()),
        ];

        assert!(!complete_flag_check(&items, &items[0], 2));
    }

    #[test]
    fn ridge_image_check_rejects_invalid_lift() {
        let mut bad_ridge_image = graph_checks_certificate();
        assert!(rocq::ridge_graph_image_check(&bad_ridge_image));
        bad_ridge_image.geom_edge_lifts[0][0] = (0, 2);
        assert!(!rocq::ridge_graph_image_check(&bad_ridge_image));
    }
}

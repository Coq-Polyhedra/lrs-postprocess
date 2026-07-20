use anyhow::{anyhow, Context, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::VecDeque;
use std::time::Instant;

use crate::certificate::{read_certificate, Certificate};
use crate::numerics::q_to_string;
use crate::postprocess::{parse_lrs_hrep, HRep};

fn parse_bigint(s: &str) -> Option<BigInt> {
    BigInt::parse_bytes(s.trim().as_bytes(), 10)
}

fn parse_bigint_vec(xs: &[String]) -> Option<Vec<BigInt>> {
    xs.iter().map(|x| parse_bigint(x)).collect()
}

fn parse_bigint_matrix(
    matrix: &[Vec<String>],
    rows: usize,
    columns: usize,
) -> Option<Vec<Vec<BigInt>>> {
    if matrix.len() != rows || matrix.iter().any(|row| row.len() != columns) {
        return None;
    }
    matrix.iter().map(|row| parse_bigint_vec(row)).collect()
}

type ParsedInequality = (Vec<BigInt>, BigInt);

fn parse_inequalities(cert: &Certificate) -> Option<Vec<ParsedInequality>> {
    cert.inequalities
        .iter()
        .map(|inequality| {
            Some((
                parse_bigint_vec(&inequality.a)?,
                parse_bigint(&inequality.b)?,
            ))
        })
        .collect()
}

fn dot(x: &[BigInt], y: &[BigInt]) -> Option<BigInt> {
    if x.len() != y.len() {
        return None;
    }
    Some(
        x.iter()
            .zip(y)
            .fold(BigInt::zero(), |acc, (xi, yi)| acc + xi * yi),
    )
}

fn sparse_dot(weight: &[(usize, String)], x: &[BigInt]) -> Option<BigInt> {
    let mut result = BigInt::zero();
    for (i, value) in weight {
        let coefficient = parse_bigint(value)?;
        result += coefficient * x.get(*i)?;
    }
    Some(result)
}

fn strictly_sorted(xs: &[usize]) -> bool {
    xs.windows(2).all(|w| w[0] < w[1])
}

fn strictly_lexicographically_sorted(xs: &[Vec<usize>]) -> bool {
    xs.windows(2).all(|w| w[0] < w[1])
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

fn parse_rational_parts(s: &str) -> Option<(BigInt, BigInt)> {
    let s = s.trim();
    let (mut numerator, mut denominator) = if let Some((a, b)) = s.split_once('/') {
        (parse_bigint(a)?, parse_bigint(b)?)
    } else {
        (parse_bigint(s)?, BigInt::one())
    };

    if denominator.is_zero() {
        return None;
    }
    if denominator.sign() == Sign::Minus {
        numerator = -numerator;
        denominator = -denominator;
    }

    let gcd = numerator.abs().gcd(&denominator);
    Some((numerator / &gcd, denominator / gcd))
}

fn clear_rationals_to_integer_row(values: &[String]) -> Option<Vec<BigInt>> {
    let parts = values
        .iter()
        .map(|value| parse_rational_parts(value))
        .collect::<Option<Vec<_>>>()?;
    let lcm = parts
        .iter()
        .fold(BigInt::one(), |acc, (_, denominator)| {
            acc.lcm(denominator)
        });
    Some(
        parts
            .into_iter()
            .map(|(numerator, denominator)| numerator * (&lcm / denominator))
            .collect(),
    )
}

fn integer_inequality_from_hrep(h: &HRep, i: usize) -> Option<(Vec<BigInt>, BigInt)> {
    let mut values = h.a.get(i)?.iter().map(q_to_string).collect::<Vec<_>>();
    values.push(q_to_string(h.b.get(i)?));
    let cleared = clear_rationals_to_integer_row(&values)?;
    let bound = cleared.last()?.clone();
    Some((cleared[..cleared.len() - 1].to_vec(), bound))
}

// Rust-specific check: the Rocq checker receives the integer inequalities
// directly and therefore has no external H-representation to compare against.
fn check_ine_file(h: &HRep, cert: &Certificate) -> bool {
    cert.n_inequalities == h.a.len()
        && cert.dimension == h.d
        && cert.inequalities.len() == h.a.len()
        && cert.inequalities.iter().enumerate().all(|(i, inequality)| {
            let Some(normal) = parse_bigint_vec(&inequality.a) else {
                return false;
            };
            let Some(bound) = parse_bigint(&inequality.b) else {
                return false;
            };
            let Some((expected_normal, expected_bound)) =
                integer_inequality_from_hrep(h, i)
            else {
                return false;
            };
            normal == expected_normal && bound == expected_bound
        })
}

mod rocq {
    #![allow(non_snake_case)]

    use super::*;

    pub fn areInequalitiesWellFormed(cert: &Certificate) -> bool {
        cert.inequalities.len() == cert.n_inequalities
            && cert.inequalities.iter().all(|inequality| {
                let Some(normal) = parse_bigint_vec(&inequality.a) else {
                    return false;
                };
                parse_bigint(&inequality.b).is_some()
                    && normal.len() == cert.dimension
                    && normal.iter().any(|x| !x.is_zero())
            })
    }

    pub fn arePointsWellFormed(cert: &Certificate) -> bool {
        cert.items.iter().all(|item| {
            let Some(numerators) = parse_bigint_vec(&item.vertex.num) else {
                return false;
            };
            let Some(denominator) = parse_bigint(&item.vertex.den) else {
                return false;
            };
            denominator > BigInt::zero() && numerators.len() == cert.dimension
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

    pub fn areVerticesWellFormed(cert: &Certificate) -> bool {
        arePointsWellFormed(cert) && areActiveSetsWellFormed(cert)
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
        let Some(point) = parse_bigint_vec(&cert.full_dim.point) else {
            return false;
        };
        let Some(denominator) = parse_bigint(&cert.full_dim.denominator) else {
            return false;
        };
        denominator > BigInt::zero() && point.len() == cert.dimension
    }

    pub fn isFullDimDirWellFormed(cert: &Certificate) -> bool {
        parse_bigint_matrix(
            &cert.full_dim.directions,
            cert.dimension,
            cert.dimension,
        )
        .is_some()
    }

    pub fn isFullDimInverseWellFormed(cert: &Certificate) -> bool {
        parse_bigint_matrix(
            &cert.full_dim.left_inverse,
            cert.dimension,
            cert.dimension,
        )
        .is_some()
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
        parse_bigint_matrix(
            &cert.root.basis_vectors,
            cert.dimension,
            cert.dimension,
        )
        .is_some()
    }

    pub fn areScalarProductsWellFormed(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        parse_bigint_matrix(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        )
        .is_some()
    }

    pub fn isSparseVectorWellFormed(dimension: usize, weight: &[(usize, String)]) -> bool {
        weight.iter().all(|(i, value)| {
            *i < dimension && parse_bigint(value).is_some_and(|x| !x.is_zero())
        }) && weight.windows(2).all(|w| w[0].0 < w[1].0)
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
        let active_sets = cert
            .items
            .iter()
            .map(|item| item.incident.clone())
            .collect::<Vec<_>>();
        strictly_lexicographically_sorted(&active_sets)
    }

    pub fn areFacetsUnique(cert: &Certificate) -> bool {
        let descriptions = cert
            .graph
            .lbl
            .iter()
            .map(|facet| facet.simplex.clone())
            .collect::<Vec<_>>();
        strictly_lexicographically_sorted(&descriptions)
    }

    pub fn check_ineqs(
        inequalities: &[ParsedInequality],
        active_set: &[usize],
        numerators: &[BigInt],
        denominator: &BigInt,
    ) -> bool {
        let mut active_position = 0;

        for (i, (normal, bound)) in inequalities.iter().enumerate() {
            let Some(lhs) = dot(normal, numerators) else {
                return false;
            };
            let rhs = bound * denominator;

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
        let Some(inequalities) = parse_inequalities(cert) else {
            return false;
        };

        cert.items.iter().all(|item| {
            let Some(numerators) = parse_bigint_vec(&item.vertex.num) else {
                return false;
            };
            let Some(denominator) = parse_bigint(&item.vertex.den) else {
                return false;
            };

            check_ineqs(
                &inequalities,
                &item.incident,
                &numerators,
                &denominator,
            )
        })
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
        let Some(witnesses) = parse_bigint_matrix(
            &cert.root.basis_vectors,
            cert.dimension,
            cert.dimension,
        ) else {
            return false;
        };
        let Some(products) = parse_bigint_matrix(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        ) else {
            return false;
        };

        root_vertex
            .incident
            .iter()
            .enumerate()
            .all(|(i, &row)| {
                let Some(inequality) = cert.inequalities.get(row) else {
                    return false;
                };
                let Some(normal) = parse_bigint_vec(&inequality.a) else {
                    return false;
                };
                (0..cert.dimension).all(|j| {
                    dot(&normal, &witnesses[j])
                        .is_some_and(|expected| products[i][j] == expected)
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
        let Some(products) = parse_bigint_matrix(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        ) else {
            return false;
        };

        root_facet
            .simplex
            .iter()
            .enumerate()
            .all(|(i, &row)| {
                let Some(&active_position) = cert.root.inverse_incident_map.get(row) else {
                    return false;
                };
                let Some(product_row) = products.get(active_position) else {
                    return false;
                };
                (0..cert.dimension).all(|j| {
                    if i == j {
                        product_row[j] > BigInt::zero()
                    } else {
                        product_row[j].is_zero()
                    }
                })
            })
    }

    pub fn isSparseVectorPositive(weight: &[(usize, String)]) -> bool {
        !weight.is_empty()
            && weight.iter().all(|(_, value)| {
                parse_bigint(value).is_some_and(|x| x >= BigInt::zero())
            })
    }

    pub fn separability_check(cert: &Certificate) -> bool {
        let Some(root_facet) = cert.graph.lbl.get(cert.root.simplex_id) else {
            return false;
        };
        let Some(root_vertex) = cert.items.get(root_facet.owner) else {
            return false;
        };
        let Some(products) = parse_bigint_matrix(
            &cert.root.m_matrix,
            root_vertex.incident.len(),
            cert.dimension,
        ) else {
            return false;
        };

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
                            products.get(active_position).is_some_and(|product_row| {
                                sparse_dot(weight, product_row)
                                    .is_some_and(|value| value <= BigInt::zero())
                            })
                        })
                })
    }

    pub fn graph_image_check(cert: &Certificate) -> bool {
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

    fn not_subset(xs: &[usize], ys: &[usize]) -> bool {
        !sorted_subset(xs, ys)
    }

    fn incomparable(xs: &[usize], ys: &[usize]) -> bool {
        not_subset(xs, ys) && not_subset(ys, xs)
    }

    pub fn geom_edge_pairwise_check(cert: &Certificate) -> bool {
        cert.items.iter().enumerate().all(|(i, vertex)| {
            let Some(neighbors) = cert.neighbors.get(i) else {
                return false;
            };
            let Some(differences) = neighbors
                .iter()
                .map(|&neighbor| {
                    cert.items.get(neighbor).map(|other_vertex| {
                        sorted_difference(&vertex.incident, &other_vertex.incident)
                    })
                })
                .collect::<Option<Vec<_>>>()
            else {
                return false;
            };

            differences.iter().all(|difference| !difference.is_empty())
                && (0..differences.len()).all(|j| {
                    ((j + 1)..differences.len())
                        .all(|k| incomparable(&differences[j], &differences[k]))
                })
        })
    }

    pub fn connectivity_check(cert: &Certificate) -> bool {
        let graph = &cert.neighbors;
        if graph.is_empty() {
            return true;
        }

        let mut visited = vec![false; graph.len()];
        let mut queue = VecDeque::new();
        visited[0] = true;
        queue.push_back(0usize);
        let mut count = 1usize;

        while let Some(vertex) = queue.pop_front() {
            let Some(neighbors) = graph.get(vertex) else {
                return false;
            };
            for &neighbor in neighbors {
                if neighbor >= graph.len() {
                    return false;
                }
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    count += 1;
                    queue.push_back(neighbor);
                }
            }
        }

        count == graph.len()
    }

    pub fn full_dim_feasibility_check(cert: &Certificate) -> bool {
        let Some(point) = parse_bigint_vec(&cert.full_dim.point) else {
            return false;
        };
        let Some(denominator) = parse_bigint(&cert.full_dim.denominator) else {
            return false;
        };
        let Some(directions) = parse_bigint_matrix(
            &cert.full_dim.directions,
            cert.dimension,
            cert.dimension,
        ) else {
            return false;
        };

        cert.inequalities.iter().all(|inequality| {
            let Some(normal) = parse_bigint_vec(&inequality.a) else {
                return false;
            };
            let Some(bound) = parse_bigint(&inequality.b) else {
                return false;
            };
            let Some(base) = dot(&normal, &point) else {
                return false;
            };
            let rhs = bound * &denominator;
            base <= rhs
                && directions.iter().all(|direction| {
                    dot(&normal, direction)
                        .is_some_and(|increment| &base + increment <= rhs)
                })
        })
    }

    pub fn full_dim_inverse_check(cert: &Certificate) -> bool {
        let Some(directions) = parse_bigint_matrix(
            &cert.full_dim.directions,
            cert.dimension,
            cert.dimension,
        ) else {
            return false;
        };
        let Some(inverse) = parse_bigint_matrix(
            &cert.full_dim.left_inverse,
            cert.dimension,
            cert.dimension,
        ) else {
            return false;
        };

        (0..cert.dimension).all(|i| {
            (0..cert.dimension).all(|j| {
                dot(&directions[i], &inverse[j]).is_some_and(|value| {
                    if i == j {
                        !value.is_zero()
                    } else {
                        value.is_zero()
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

    pub fn geom_graph_check(cert: &Certificate) -> bool {
        graph_image_check(cert)
            && geom_edge_pairwise_check(cert)
            && connectivity_check(cert)
    }

    // This follows the literal Rocq definition. As in the Rocq benchmark,
    // full_dim_check is evaluated separately by the file-level entry point.
    pub fn check_certificate(cert: &Certificate) -> bool {
        well_formedness_check(cert)
            && uniqueness_check(cert)
            && feasibility_check(cert)
            && graph_check(cert)
            && mapping_check(cert)
            && root_check(cert)
            && geom_graph_check(cert)
    }
}

fn timed_bool(name: &str, check: impl FnOnce() -> bool) -> bool {
    let start = Instant::now();
    let result = check();
    eprintln!("{name}: {:.6} s", start.elapsed().as_secs_f64());
    result
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

    let inequality_file_check =
        timed_bool("Inequality file check", || check_ine_file(&h, &cert));
    let well_formedness =
        timed_bool("Well-formedness check", || rocq::well_formedness_check(&cert));
    let uniqueness = timed_bool("Uniqueness check", || rocq::uniqueness_check(&cert));
    let feasibility = timed_bool("Feasibility check", || rocq::feasibility_check(&cert));
    let graph = timed_bool("Graph check", || rocq::graph_check(&cert));
    let mapping = timed_bool("Mapping check", || rocq::mapping_check(&cert));
    let root = timed_bool("Root check", || rocq::root_check(&cert));
    let geometric_graph =
        timed_bool("Geometric graph check", || rocq::geom_graph_check(&cert));
    let full_dimension =
        timed_bool("Full dimension check", || rocq::full_dim_check(&cert));

    let accepted = inequality_file_check
        && well_formedness
        && uniqueness
        && feasibility
        && graph
        && mapping
        && root
        && geometric_graph
        && full_dimension;

    accepted
        .then_some(())
        .ok_or_else(|| anyhow!("certificate rejected"))
}

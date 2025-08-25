//Jacobi polynomials enum (database)
//det Jacobi database

//calculations

//first calculate mesh quality

//then number of optimal gaussian quadrature points

use std::fs::File;                          
use std::io::{self};    // For input/output operations 

use crate::database::*;                // Import mesh data structures and error types from database module



#[derive(Debug, Clone, Copy, PartialEq)]
enum Topology {
    Line,
    Triangle,
    Quad,
    Tetra,
    Pyramid,
    Wedge,
    Hexahedron,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Family { //serendipity and lagrange elements are the same for the first order elements
    Lagrange,
    Serendipity,
    Hierarchical,
    Mixed,
}

#[derive(Debug, Clone, Copy)]
struct ElementSignature {
    topology: Topology,
    family: Family,
    order: usize,
}

#[derive(Debug, Clone)]
struct ShapeFunction {
    values: Vec<f64>,
    derivatives: Vec<Vec<f64>>, // derivatives[node][dimension]
}

#[derive(Debug, Clone)]
pub struct Jacobian {
    pub matrix: Vec<Vec<f64>>, // Jacobian matrix J[i][j] = dx_i/dxi_j
    pub determinant: f64,
}

impl ElementType {
    fn from_element_type_to_signature(element_type: ElementType) -> Result<ElementSignature, ParseError> {
        match element_type {
            ElementType::Line => Ok(ElementSignature { 
                topology: Topology::Line,
                family: Family::Lagrange,  //Lagrange??
                order: 1,
            }),
            ElementType::QuadraticEdge => Ok(ElementSignature {
                topology: Topology::Line,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::Triangle => Ok(ElementSignature {
                topology: Topology::Triangle,
                family: Family::Lagrange,
                order: 1,
            }),
            ElementType::QuadraticTriangle => Ok(ElementSignature {
                topology: Topology::Triangle,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::Quad => Ok(ElementSignature {
                topology: Topology::Quad,
                family: Family::Lagrange,
                order: 1,
            }),
            ElementType::QuadraticQuad => Ok(ElementSignature {
                topology: Topology::Quad,
                family: Family::Serendipity,
                order: 2,
            }),
            ElementType::BiquadraticQuad => Ok(ElementSignature {
                topology: Topology::Quad,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::Tetra => Ok(ElementSignature {
                topology: Topology::Tetra,
                family: Family::Lagrange,
                order: 1,
            }),
            ElementType::QuadraticTetra => Ok(ElementSignature {
                topology: Topology::Tetra,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::Pyramid => Ok(ElementSignature {
                topology: Topology::Pyramid,
                family: Family::Lagrange,
                order: 1,
            }),
            ElementType::QuadraticPyramid => Ok(ElementSignature {
                topology: Topology::Pyramid,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::Wedge => Ok(ElementSignature {
                topology: Topology::Wedge,
                family: Family::Lagrange,
                order: 1,
            }),
            ElementType::QuadraticWedge => Ok(ElementSignature {
                topology: Topology::Wedge,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::BiquadraticQuadraticWedge => Ok(ElementSignature {
                topology: Topology::Wedge,
                family: Family::Mixed,
                order: 2,
            }),
            ElementType::Hexahedron => Ok(ElementSignature {
                topology: Topology::Hexahedron,
                family: Family::Lagrange,
                order: 1,
            }),
            ElementType::QuadraticHexahedron => Ok(ElementSignature {
                topology: Topology::Hexahedron,
                family: Family::Serendipity,
                order: 2,
            }),
            ElementType::BiquadraticQuadraticHexahedron => Ok(ElementSignature {
                topology: Topology::Hexahedron,
                family: Family::Mixed,
                order: 2,
            }),
            ElementType::TriquadraticHexahedron => Ok(ElementSignature {
                topology: Topology::Hexahedron,
                family: Family::Lagrange,
                order: 2,
            }),
            ElementType::Vertex => Err(ParseError::ElementError("Vertex element - skipping Jacobian calculation".to_string())),
        }
    }

    fn get_shape_functions(
        signature: ElementSignature,
        xi: &[f64], // Natural coordinates
    ) -> Result<ShapeFunction, ParseError> {
        match signature.topology {
            Topology::Line => {
                get_line_shape_functions(signature.family, signature.order, xi)
            }
            Topology::Triangle => {
                get_triangle_shape_functions(signature.family, signature.order, xi)
            }
            Topology::Quad => {
                get_quad_shape_functions(signature.family, signature.order, xi)
            }
            Topology::Tetra => {
                get_tetra_shape_functions(signature.family, signature.order, xi)
            }
            _ => Err(ParseError::FormatError(format!(
                "No shape function implementation for {:?} {:?} order {}",
                signature.topology, signature.family, signature.order
            ))),
        }
    }

    /// Compute Jacobian matrix for an element at given natural coordinates
    pub fn get_jacobian(
        element: &Element,
        nodes: &[Node],
        element_type: ElementType,
        xi: &[f64], // Natural coordinates
    ) -> Result<Jacobian, ParseError> {
        // Handle vertex elements gracefully
        if matches!(element_type, ElementType::Vertex) {
            return Err(ParseError::ElementError("Vertex element - skipping Jacobian calculation".to_string()));
        }

        // Get element signature and shape functions
        let signature = Self::from_element_type_to_signature(element_type)?;
        let shape_fn = Self::get_shape_functions(signature, xi)?;
        
        // Determine spatial dimension from nodes
        if nodes.is_empty() {
            return Err(ParseError::FormatError("No nodes provided for Jacobian calculation".to_string()));
        }
        
        let spatial_dim = nodes[0].coordinates.len();
        let natural_dim = xi.len();
        
        // Initialize Jacobian matrix
        let mut jacobian_matrix = vec![vec![0.0; natural_dim]; spatial_dim];
        
        // Compute Jacobian: J[i][j] = sum(N_k,j * x_k_i) where N_k,j is derivative of shape function k w.r.t. natural coordinate j
        for i in 0..spatial_dim {
            for j in 0..natural_dim {
                for (node_idx, &global_node_id) in element.nodes.iter().enumerate() {
                    // Find the node in the nodes array
                    if let Some(node) = nodes.iter().find(|n| n.id == global_node_id) {
                        if node_idx < shape_fn.derivatives.len() && j < shape_fn.derivatives[node_idx].len() {
                            jacobian_matrix[i][j] += shape_fn.derivatives[node_idx][j] * node.coordinates[i];
                        }
                    }
                }
            }
        }
        
        // Calculate determinant
        let det = Self::calculate_determinant(&jacobian_matrix)?;
        
        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant: det,
        })
    }

    /// Calculate determinant of a matrix
    fn calculate_determinant(matrix: &[Vec<f64>]) -> Result<f64, ParseError> {
        let n = matrix.len();
        if n == 0 {
            return Err(ParseError::FormatError("Cannot calculate determinant of empty matrix".to_string()));
        }
        
        // Check if matrix is square
        for row in matrix {
            if row.len() != matrix[0].len() {
                return Err(ParseError::FormatError("Matrix must be rectangular for determinant calculation".to_string()));
            }
        }
        
        let m = matrix[0].len();
        
        match (n, m) {
            (1, 1) => Ok(matrix[0][0]),
            (2, 2) => Ok(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]),
            (3, 3) => {
                Ok(matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))
            }
            _ => {
                // For non-square matrices or higher dimensions, use a more general approach
                // For now, return an error for unsupported cases
                Err(ParseError::FormatError(format!(
                    "Determinant calculation not implemented for {}x{} matrices", n, m
                )))
            }
        }
    }
}

fn get_line_shape_functions(
    family: Family,
    order: usize,
    xi: &[f64],
) -> Result<ShapeFunction, ParseError> {
    if xi.len() != 1 {
        return Err(ParseError::FormatError("Line elements require 1 natural coordinate".to_string()));
    }
    
    let xi = xi[0];
    
    match (family, order) {
        (Family::Lagrange, 1) => {
            // Linear Lagrange line element
            let values = vec![
                0.5 * (1.0 - xi), // N1
                0.5 * (1.0 + xi), // N2
            ];
            let derivatives = vec![
                vec![-0.5], // dN1/dxi
                vec![0.5],  // dN2/dxi
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        (Family::Lagrange, 2) => {
            // Quadratic Lagrange line element
            let values = vec![
                0.5 * xi * (xi - 1.0),     // N1
                1.0 - xi * xi,             // N2 (middle node)
                0.5 * xi * (xi + 1.0),     // N3
            ];
            let derivatives = vec![
                vec![xi - 0.5],            // dN1/dxi
                vec![-2.0 * xi],           // dN2/dxi
                vec![xi + 0.5],            // dN3/dxi
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        (Family::Lagrange, order) if order > 2 => {
            // General Lagrange interpolation for higher orders
            let n_nodes = order + 1;
            let mut values = vec![0.0; n_nodes];
            let mut derivatives = vec![vec![0.0; 1]; n_nodes];
            
            // Node positions in natural coordinates
            let mut xi_nodes = vec![0.0; n_nodes];
            for i in 0..n_nodes {
                xi_nodes[i] = -1.0 + 2.0 * (i as f64) / (order as f64);
            }
            
            // Compute Lagrange shape functions
            for i in 0..n_nodes {
                let mut li = 1.0;
                let mut dli_dxi = 0.0;
                
                for j in 0..n_nodes {
                    if i != j {
                        li *= (xi - xi_nodes[j]) / (xi_nodes[i] - xi_nodes[j]);
                        
                        // Derivative using product rule
                        let mut temp = 1.0 / (xi_nodes[i] - xi_nodes[j]);
                        for k in 0..n_nodes {
                            if k != i && k != j {
                                temp *= (xi - xi_nodes[k]) / (xi_nodes[i] - xi_nodes[k]);
                            }
                        }
                        dli_dxi += temp;
                    }
                }
                
                values[i] = li;
                derivatives[i][0] = dli_dxi;
            }
            
            Ok(ShapeFunction { values, derivatives })
        }
        _ => Err(ParseError::FormatError(format!(
            "No shape function implementation for Line {:?} order {}",
            family, order
        ))),
    }
}

fn get_triangle_shape_functions(
    family: Family,
    order: usize,
    xi: &[f64],
) -> Result<ShapeFunction, ParseError> {
    if xi.len() != 2 {
        return Err(ParseError::FormatError("Triangle elements require 2 natural coordinates".to_string()));
    }
    
    let xi = xi[0];
    let eta = xi[1];
    let zeta = 1.0 - xi - eta;
    
    match (family, order) {
        (Family::Lagrange, 1) => {
            // Linear triangle
            let values = vec![
                zeta,  // N1
                xi,    // N2
                eta,   // N3
            ];
            let derivatives = vec![
                vec![-1.0, -1.0], // dN1/dxi, dN1/deta
                vec![1.0, 0.0],   // dN2/dxi, dN2/deta
                vec![0.0, 1.0],   // dN3/dxi, dN3/deta
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        (Family::Lagrange, 2) => {
            // Quadratic triangle (6 nodes)
            let values = vec![
                zeta * (2.0 * zeta - 1.0),           // N1 (corner)
                xi * (2.0 * xi - 1.0),               // N2 (corner)
                eta * (2.0 * eta - 1.0),             // N3 (corner)
                4.0 * xi * zeta,                     // N4 (edge 1-2)
                4.0 * xi * eta,                      // N5 (edge 2-3)
                4.0 * eta * zeta,                    // N6 (edge 3-1)
            ];
            let derivatives = vec![
                vec![4.0 * xi + 4.0 * eta - 3.0, 4.0 * xi + 4.0 * eta - 3.0], // dN1
                vec![4.0 * xi - 1.0, 0.0],                                      // dN2
                vec![0.0, 4.0 * eta - 1.0],                                     // dN3
                vec![4.0 * (1.0 - 2.0 * xi - eta), -4.0 * xi],                 // dN4
                vec![4.0 * eta, 4.0 * xi],                                      // dN5
                vec![-4.0 * eta, 4.0 * (1.0 - xi - 2.0 * eta)],                // dN6
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        _ => Err(ParseError::FormatError(format!(
            "No shape function implementation for Triangle {:?} order {}",
            family, order
        ))),
    }
}

fn get_quad_shape_functions(
    family: Family,
    order: usize,
    xi: &[f64],
) -> Result<ShapeFunction, ParseError> {
    if xi.len() != 2 {
        return Err(ParseError::FormatError("Quad elements require 2 natural coordinates".to_string()));
    }
    
    let xi = xi[0];
    let eta = xi[1];
    
    match (family, order) {
        (Family::Lagrange, 1) => {
            // Bilinear quad
            let values = vec![
                0.25 * (1.0 - xi) * (1.0 - eta), // N1
                0.25 * (1.0 + xi) * (1.0 - eta), // N2
                0.25 * (1.0 + xi) * (1.0 + eta), // N3
                0.25 * (1.0 - xi) * (1.0 + eta), // N4
            ];
            let derivatives = vec![
                vec![-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)], // dN1
                vec![0.25 * (1.0 - eta), -0.25 * (1.0 + xi)],  // dN2
                vec![0.25 * (1.0 + eta), 0.25 * (1.0 + xi)],   // dN3
                vec![-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)],  // dN4
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        (Family::Serendipity, 2) => {
            // 8-node serendipity quad
            let values = vec![
                0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0), // N1 (corner)
                0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),  // N2 (corner)
                0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),  // N3 (corner)
                0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0), // N4 (corner)
                0.5 * (1.0 - xi * xi) * (1.0 - eta),                 // N5 (edge)
                0.5 * (1.0 + xi) * (1.0 - eta * eta),                // N6 (edge)
                0.5 * (1.0 - xi * xi) * (1.0 + eta),                 // N7 (edge)
                0.5 * (1.0 - xi) * (1.0 - eta * eta),                // N8 (edge)
            ];
            
            // Calculate derivatives for serendipity elements
            let derivatives = vec![
                vec![-0.25 * (1.0 - eta) * (-xi - eta - 1.0) - 0.25 * (1.0 - xi) * (1.0 - eta),
                     -0.25 * (1.0 - xi) * (-xi - eta - 1.0) - 0.25 * (1.0 - xi) * (1.0 - eta)], // dN1
                vec![0.25 * (1.0 - eta) * (xi - eta - 1.0) + 0.25 * (1.0 + xi) * (1.0 - eta),
                     -0.25 * (1.0 + xi) * (xi - eta - 1.0) - 0.25 * (1.0 + xi) * (1.0 - eta)], // dN2
                vec![0.25 * (1.0 + eta) * (xi + eta - 1.0) + 0.25 * (1.0 + xi) * (1.0 + eta),
                     0.25 * (1.0 + xi) * (xi + eta - 1.0) + 0.25 * (1.0 + xi) * (1.0 + eta)], // dN3
                vec![-0.25 * (1.0 + eta) * (-xi + eta - 1.0) - 0.25 * (1.0 - xi) * (1.0 + eta),
                     0.25 * (1.0 - xi) * (-xi + eta - 1.0) + 0.25 * (1.0 - xi) * (1.0 + eta)], // dN4
                vec![-xi * (1.0 - eta), -0.5 * (1.0 - xi * xi)], // dN5
                vec![0.5 * (1.0 - eta * eta), -(1.0 + xi) * eta], // dN6
                vec![-xi * (1.0 + eta), 0.5 * (1.0 - xi * xi)], // dN7
                vec![-0.5 * (1.0 - eta * eta), -(1.0 - xi) * eta], // dN8
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        _ => Err(ParseError::FormatError(format!(
            "No shape function implementation for Quad {:?} order {}",
            family, order
        ))),
    }
}

fn get_tetra_shape_functions(
    family: Family,
    order: usize,
    xi: &[f64],
) -> Result<ShapeFunction, ParseError> {
    if xi.len() != 3 {
        return Err(ParseError::FormatError("Tetra elements require 3 natural coordinates".to_string()));
    }
    
    let xi = xi[0];
    let eta = xi[1];
    let zeta = xi[2];
    let lambda = 1.0 - xi - eta - zeta;
    
    match (family, order) {
        (Family::Lagrange, 1) => {
            // Linear tetrahedron
            let values = vec![
                lambda, // N1
                xi,     // N2
                eta,    // N3
                zeta,   // N4
            ];
            let derivatives = vec![
                vec![-1.0, -1.0, -1.0], // dN1
                vec![1.0, 0.0, 0.0],    // dN2
                vec![0.0, 1.0, 0.0],    // dN3
                vec![0.0, 0.0, 1.0],    // dN4
            ];
            Ok(ShapeFunction { values, derivatives })
        }
        _ => Err(ParseError::FormatError(format!(
            "No shape function implementation for Tetra {:?} order {}",
            family, order
        ))),
    }
}
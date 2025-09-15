//use std::fs::File;                          
use std::io::{self};    // For input/output operations 
use std::f64;            // Import f64 constants like INFINITY

//use crate::database::*;                     // Import mesh data structures and error types from database module

use crate::lib::*;                  // Import mesh data structures and error types from lib module
use crate::error::*;                     // Import mesh data structures and error types from error module

pub struct GeometricAnalysis;  // Defines a structure to represent the geometric analysis for the mesh elements

impl GeometricAnalysis {

    /// Main method to analyze quality of all elements in the mesh using reference elements
    pub fn analyse_mesh_quality(mesh_data: &MeshData) -> Result<MeshQualityReport, ElementError> {
        let mut element_qualities = Vec::new();  // Vector to store quality metrics for each element
        let mut processed_elements = 0;         // Counter for successfully processed elements
        
        // Process elements by type groups using the element_type_info from parsed data
        for type_info in &mesh_data.element_type_info {
            // Skip vertex elements as they don't have meaningful quality metrics
            if matches!(type_info.element_type, ElementType::Vertex) {
                continue;  // Move to next element type
            }
            
            // Calculate the range of elements for this type
            let start_idx = type_info.start_index;                    // First element index of this type
            let end_idx = start_idx + type_info.num_elements;         // One past last element index of this type
            
            // Process each element of this specific type
            for element_idx in start_idx..end_idx {
                // Check if element index is valid (within bounds of elements vector)
                if element_idx < mesh_data.elements.len() {
                    let element = &mesh_data.elements[element_idx];  // Get reference to current element
                    
                    // Attempt to calculate quality metrics for this element
                    match Self::calculate_element_quality(element, &type_info.element_type, &mesh_data.nodes) {
                        Ok(quality) => {
                            // Success: store the quality metrics and increment counter
                            element_qualities.push(quality);
                            processed_elements += 1;
                        },
                        Err(e) => {
                            // Failure: print warning but continue processing other elements
                            println!("Warning: Failed to analyze element {}: {:?}", element.id, e);
                        }
                    }
                }
            }
        }
        
        // Check if we successfully analyzed any elements
        if element_qualities.is_empty() {
            return Err(ElementError::GeometryError("No elements could be analyzed".to_string()));
        }
        
        // Calculate overall statistics from individual element qualities
        let statistics = Self::calculate_statistics(&element_qualities);
        
        
        // Return the complete quality report
        Ok(MeshQualityReport {
            total_elements: processed_elements,
            element_qualities,
            statistics,
        })
    }

    /// Calculate element quality using Jacobian analysis
    fn calculate_element_quality(
        element: &Element,                    // The element to analyze
        element_type: &ElementType,
        nodes: &[Node]                       // All nodes in the mesh
    ) -> Result<ElementQuality, ElementError> {
        // Get the actual node coordinates for this element
        let element_nodes = Self::get_element_nodes(element, nodes)?; //from MeshData
        
        // Get shape functions and their derivatives for the element type
        let shape_function = Self::get_shape_functions(element_type)?;
        let shape_derivatives = shape_function.derivatives;

        // Calculate Jacobian
        let jacobian = Self::calculate_jacobian(&element_nodes, &shape_derivatives)?;
        
        // Return all quality metrics
        Ok(ElementQuality {
            element_id: element.id,
            det_jacobian: jacobian.determinant,
        })
    }

    /// Create shape functions for reference elements
    pub fn get_shape_functions(element_type: &ElementType) -> Result<ShapeFunction, ParseError> {  //ShapeFunctionError
        match element_type {
            ElementType::Line => {
                let reference_coords = vec![
                    vec![0.0],    // Node 0 at xi = 0
                    vec![1.0],     // Node 1 at xi = 1
                ];
                let num_nodes = 2;

                let xi = 0.2113 ;
                
                let values = vec![
                    1.0 - xi,   // N0: node at xi = 0
                    xi,         // N1: node at xi = 1
                ];
                let derivatives = vec![
                    vec![-1.0], // dN0/dxi
                    vec![1.0],  // dN1/dxi
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::QuadraticEdge => {
                let reference_coords = vec![
                    vec![0.0],      // Node 0 at xi = 0
                    vec![1.0],      // Node 1 at xi = 1
                    vec![0.5],      // Node 2 at xi = 0.5 (middle)
                ];
                let num_nodes = 3;

                let xi = 0.1127;
                
                let values = vec![
                    2.0 * (xi - 0.5) * (xi - 1.0),      // N0: node at xi = 0
                    2.0 * xi * (xi - 0.5),              // N1: node at xi = 1
                    4.0 * xi * (1.0 - xi),              // N2: node at xi = 0.5
                ];
                let derivatives = vec![
                    vec![4.0 * xi - 3.0],               // dN0/dxi
                    vec![4.0 * xi - 1.0],               // dN1/dxi
                    vec![4.0 - xi * 8.0],               // dN2/dxi
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::Triangle => {
                let reference_coords = vec![
                    vec![0.0, 0.0],    // Node 0 at (0,0)
                    vec![1.0, 0.0],    // Node 1 at (1,0)
                    vec![0.0, 1.0],    // Node 2 at (0,1)
                ];
                let num_nodes = 3;

                let xi = 1.0/3.0;
                let eta = 1.0/3.0;

                let values = vec![
                    1.0 - xi - eta,               // N0: node at (0,0)
                    xi,                              // N1: node at (1,0)
                    eta,                             // N2: node at (0,1)
                ];
                let derivatives = vec![
                    vec![-1.0, -1.0],       // dN0/dxi, dN0/deta
                    vec![1.0, 0.0],         // dN1/dxi, dN1/deta
                    vec![0.0, 1.0],         // dN2/dxi, dN2/deta
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::QuadraticTriangle => {
                let reference_coords = vec![
                    vec![0.0, 0.0],    // Node 0 at (0,0)
                    vec![1.0, 0.0],    // Node 1 at (1,0)
                    vec![0.0, 1.0],    // Node 2 at (0,1)
                    vec![0.5, 0.0],    // Node 3: mid-edge between 0-1
                    vec![0.5, 0.5],    // Node 4: mid-edge between 1-2
                    vec![0.0, 0.5],    // Node 5: mid-edge between 2-0
                ];
                let num_nodes = 6;

                let xi = 2.0 / 3.0;
                let eta = 1.0 / 6.0;

                let lambda = 1.0 - xi - eta;
                
                let values = vec![
                    lambda * (2.0 * lambda - 1.0),      // N0: corner at 
                    xi * (2.0 * xi - 1.0),              // N1: corner at 
                    eta * (2.0 * eta - 1.0),            // N2: corner at 
                    4.0 * xi * lambda,                  // N3: mid-edge at 
                    4.0 * xi * eta,                     // N4: mid-edge at 
                    4.0 * eta * lambda,                 // N5: mid-edge at 
                ];
                let derivatives = vec![
                    vec![4.0 * xi + 4.0 * eta - 3.0, 4.0 * xi + 4.0 * eta - 3.0],       // dN0
                    vec![4.0 * xi - 1.0, 0.0],                                          // dN1
                    vec![0.0, 4.0 * eta - 1.0],                                         // dN2
                    vec![4.0 - 8.0 * xi - 4.0 * eta, -4.0 * xi],                        // dN3
                    vec![4.0 * eta,  4.0 * xi],                                         // dN4
                    vec![-4.0 * eta, 4.0 - 8.0 * eta - 4.0 * xi],                        // dN5
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::Quad => {
                let reference_coords = vec![
                    vec![0.0, 0.0],    // Node 0: bottom-left corner
                    vec![1.0, 0.0],    // Node 1: bottom-right corner
                    vec![1.0, 1.0],    // Node 2: top-right corner
                    vec![0.0, 1.0],    // Node 3: top-left corner
                ];
                let num_nodes = 4;

                let xi = 0.2113;
                let eta = 0.7887;

                let xim = 1.0 - xi;
                let etam = 1.0 - eta;
                
                let values = vec![
                    xim * etam,         // N0: 
                    xi * etam,          // N1: 
                    xi * eta,           // N2: 
                    xim * eta,          // N3: 
                ];
                let derivatives = vec![
                    vec![-etam, -xim],      // dN0
                    vec![etam, -xi],        // dN1
                    vec![eta, xi],          // dN2
                    vec![-eta, xim],        // dN3
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::QuadraticQuad => {
                let reference_coords = vec![
                    vec![0.0, 0.0],    // Node 0: corner bottom-left
                    vec![1.0, 0.0],     // Node 1: corner bottom-right
                    vec![1.0, 1.0],      // Node 2: corner top-right
                    vec![0.0, 1.0],     // Node 3: corner top-left
                    vec![0.5, 0.0],     // Node 4: mid-edge bottom
                    vec![1.0, 0.5],      // Node 5: mid-edge right
                    vec![0.5, 1.0],      // Node 6: mid-edge top
                    vec![0.0, 0.5],     // Node 7: mid-edge left
                ];
                let num_nodes = 8;

                let xi = 0.1127;
                let eta = 0.7887;

                let values = vec![
                    (1.0 - xi) * (1.0 - eta) - 0.5 * (4.0 * xi * (1.0 - xi) * (1.0 - eta) + 4.0 * (1.0 - xi) * (1.0 - eta) * eta),          // N0: corner
                    xi * (1.0 - eta) - 0.5 * (4.0 * xi * (1.0 - xi) * (1.0 - eta) + 4.0 * xi * (1.0 - eta) * eta),                          // N1: corner
                    xi * eta - 0.5 * (4.0 * xi * (1.0 - eta) * eta + 4.0 * xi * (1.0 - xi) * eta),                                          // N2: corner
                    (1.0 - xi) * eta - 0.5 * (4.0 * xi * (1.0 - xi) * eta + 4.0 * (1.0 - xi) * (1.0 - eta) * eta),                          // N3: corner
                    4.0 * xi * (1.0 - xi) * (1.0 - eta),                                                                                  // N4: mid-edge bottom
                    4.0 * xi * (1.0 - eta) * eta,                                                                                         // N5: mid-edge right
                    4.0 * xi * (1.0 - xi) * eta,                                                                                          // N6: mid-edge top
                    4.0 * (1.0 - xi) * (1.0 - eta) * eta,                                                                                 // N7: mid-edge left
                ];
                let derivatives = vec![
                    vec![-(1.0 - eta) - 0.5 * (4.0 * (1.0 - eta) * (1.0 - 2.0 * xi) - 4.0 * (1.0 - eta) * eta), -(1.0 - xi) - 0.5 * (-4.0 * xi * (1.0 - xi) +  4.0 * (1.0 - xi) * (1.0 - 2.0 * eta))],      // dN0
                    vec![(1.0 - eta) - 0.5 * (4.0 * (1.0 - eta) * (1.0 - 2.0 * xi) + 4.0 * (1.0 - eta) * eta), - xi - 0.5 * (-4.0 * xi * (1.0 - xi) +  4.0 * xi * (1.0 - 2.0 * eta))],               // dN1
                    vec![eta - 0.5 * (4.0 * (1.0 - eta) * eta +  4.0 * eta * (1.0 - 2.0 * xi)), xi - 0.5 * (4.0 * xi * (1.0 - 2.0 * eta) + 4.0 * xi * (1.0 - xi))],                                         // dN2
                    vec![-eta - 0.5 * (4.0 * eta * (1.0 - 2.0 * xi) + (-4.0) * (1.0 - eta) * eta), (1.0 - xi) - 0.5 * (4.0 * xi * (1.0 - xi) +  4.0 * (1.0 - xi) * (1.0 - 2.0 * eta))],                // dN3
                    vec![4.0 * (1.0 - eta) * (1.0 - 2.0 * xi), -4.0 * xi * (1.0 - xi)],                                                                                                             // dN4
                    vec![4.0 * (1.0 - eta) * eta, 4.0 * xi * (1.0 - 2.0 * eta)],                                                                                                                            // dN5
                    vec![4.0 * eta * (1.0 - 2.0 * xi), 4.0 * xi * (1.0 - xi)],                                                                                                                      // dN6
                    vec![-4.0 * (1.0 - eta) * eta, 4.0 * (1.0 - xi) * (1.0 - 2.0 * eta)],                                                                                                                   // dN7
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::BiquadraticQuad => {
                let reference_coords = vec![
                    vec![0.0, 0.0],     // Node 0: corner bottom-left
                    vec![1.0, 0.0],     // Node 1: corner bottom-right
                    vec![1.0, 1.0],     // Node 2: corner top-right
                    vec![0.0, 1.0],     // Node 3: corner top-left
                    vec![0.5, 0.0],     // Node 4: mid-edge bottom
                    vec![1.0, 0.5],     // Node 5: mid-edge right
                    vec![0.5, 1.0],     // Node 6: mid-edge top
                    vec![0.0, 0.5],     // Node 7: mid-edge left
                    vec![0.5, 0.5],     // Node 8: center node
                ];
                let num_nodes = 9;

                let xi = 0.1127;
                let eta = 0.7887;

                let values = vec![
                    4.0 * (1.0 - xi) * (xi - 0.5) * (1.0 - eta) * (eta - 0.5),      // N0:  
                    -4.0 * xi * (xi - 0.5) * (1.0 - eta) * (eta - 0.5),             // N1:  
                    4.0 * xi * (xi - 0.5) * eta * (eta - 0.5),                      // N2:  
                    -4.0 * (1.0 - xi) * (xi - 0.5) * eta * (eta - 0.5),             // N3:  
                    8.0 * xi * (1.0 - xi) * (1.0 - eta) * (0.5 - eta),              // N4:  
                    -8.0 * xi * (0.5 - xi) * (1.0 - eta) * eta,                     // N5:  
                    -8.0 * xi * (1.0 - xi) * eta * (0.5 - eta),                     // N6:  
                    8.0 * (1.0 - xi) * (0.5 - xi) * (1.0 - eta) * eta,              // N7:  
                    16.0 * xi * (1.0 - xi) * (1.0 - eta) * eta,                     // N8: 
                ];
                let derivatives = vec![
                    vec![4.0 * (1.5 - 2.0 * xi) * (1.0 - eta) * (eta - 0.5), 4.0 * (1.0 - xi) * (xi - 0.5) * (1.5 - 2.0 * eta)],    // dN0
                    vec![-4.0 * (2.0 * xi - 0.5) * (1.0 - eta) * (eta - 0.5), -4.0 * xi * (xi - 0.5) * (1.5 - 2.0 * eta)],          // dN1
                    vec![4.0 * (2.0 * xi - 0.5) * eta * (eta - 0.5), 4.0 * xi * (xi - 0.5) * (2.0 * eta - 0.5)],                    // dN2
                    vec![-4.0 * (1.5 - 2.0 * xi) * eta * (eta - 0.5), -4.0 * (1.0 - xi) * (xi - 0.5) * (2.0 * eta - 0.5)],          // dN3
                    vec![8.0 * (1.0 - 2.0 * xi) * (1.0 - eta) * (0.5 - eta), 8.0 * xi * (1.0 - xi) * (2.0 * eta - 1.5)],            // dN4
                    vec![-8.0 * (0.5 - 2.0 * xi) * (1.0 - eta) * eta, -8.0 * xi * (0.5 - xi) * (1.0 - 2.0 * eta)],                  // dN5
                    vec![-8.0 * (1.0 - 2.0 * xi) * eta * (0.5 - eta), -8.0 * xi * (1.0 - xi) * (0.5 - 2.0 * eta)],                  // dN6
                    vec![8.0 * (2.0 * xi - 1.5) * (1.0 - eta) * eta, 8.0 * (1.0 - xi) * (0.5 - xi) * (1.0 - 2.0 * eta)],            // dN7
                    vec![16.0 * (1.0 - 2.0 * xi) * (1.0 - eta) * eta, 16.0 * xi * (1.0 - xi) * (1.0 - 2.0 * eta)],                  // dN8
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::Tetra => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],     // Node 0: at origin
                    vec![1.0, 0.0, 0.0],     // Node 1: on xxx axis
                    vec![0.0, 1.0, 0.0],     // Node 2: on xxx axis
                    vec![0.0, 0.0, 1.0],     // Node 3: on xxx axis
                ];
                let num_nodes = 4;

                let xi = 0.25;
                let eta = 0.25;
                let psi = 0.25;

                let values = vec![
                    1.0 - xi - eta - psi,       // N0: node at 
                    xi,                         // N1: node at 
                    eta,                        // N2: node at 
                    psi,                        // N3: node at
                ];
                let derivatives = vec![
                    vec![-1.0, -1.0, -1.0],     // dN0
                    vec![1.0, 0.0, 0.0],        // dN1
                    vec![0.0, 1.0, 0.0],        // dN2
                    vec![0.0, 0.0, 1.0],        // dN3
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::QuadraticTetra => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],     // Node 0: at origin
                    vec![1.0, 0.0, 0.0],     // Node 1: on xxx axis
                    vec![0.0, 1.0, 0.0],     // Node 2: on xxx axis
                    vec![0.0, 0.0, 1.0],     // Node 3: on xxx axis
                    vec![0.5, 0.0, 0.0],     // Node 4: mid-edge 0-1
                    vec![0.5, 0.5, 0.0],     // Node 5: mid-edge 1-2
                    vec![0.0, 0.5, 0.0],     // Node 6: mid-edge 0-2
                    vec![0.0, 0.0, 0.5],     // Node 7: mid-edge 0-3
                    vec![0.5, 0.0, 0.5],     // Node 8: mid-edge 1-3
                    vec![0.0, 0.5, 0.5],     // Node 9: mid-edge 2-3
                ];
                let num_nodes = 10;

                let xi = 0.5854;
                let eta = 0.1382;
                let psi = 0.1382;

                let lambda = 1.0 - xi - eta - psi;

                let values = vec![
                    lambda * (2.0 * lambda - 1.0),          // N0: corner
                    xi * (2.0 * xi - 1.0),                  // N1: corner
                    eta * (2.0 * eta - 1.0),                // N2: corner
                    psi * (2.0 * psi - 1.0),                // N3: corner
                    4.0 * lambda * xi,                      // N4: mid-edge
                    4.0 * xi * eta,                         // N5: mid-edge
                    4.0 * eta * lambda,                     // N6: mid-edge
                    4.0 * lambda * psi,                     // N7: mid-edge
                    4.0 * xi * psi,                         // N8: mid-edge
                    4.0 * eta * psi,                        // N9: mid-edge
                ];
                let derivatives = vec![
                    vec![4.0 * (xi + eta + psi) - 3.0, 4.0 * (xi + eta + psi) - 3.0, 4.0 * (xi + eta + psi) - 3.0], // dN0
                    vec![4.0 * xi - 1.0, 0.0, 0.0],                                                                 // dN1
                    vec![0.0, 4.0 * eta - 1.0, 0.0],                                                                // dN2
                    vec![0.0, 0.0, 4.0 * psi - 1.0],                                                                // dN3
                    vec![4.0 - 8.0 * xi - 4.0 * eta - 4.0 * psi, -4.0 * xi, -4.0 * xi],                             // dN4
                    vec![4.0 * eta, 4.0 * xi, 0.0],                                                                 // dN5
                    vec![-4.0 * eta, 4.0 - 4.0 * xi - 8.0 * eta - 4.0 * psi, -4.0 * eta],                           // dN6
                    vec![-4.0 * psi, -4.0 * psi, 4.0 - 4.0 * xi - 4.0 * eta - 8.0 * psi],                           // dN7
                    vec![4.0 * psi, 0.0, 4.0 * xi],                                                                 // dN8
                    vec![0.0, 4.0 * psi, 4.0 * eta],                                                                // dN9
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::Pyramid => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],    // Node 0: bottom-left corner at base
                    vec![1.0, 0.0, 0.0],     // Node 1: bottom-right corner at base
                    vec![1.0, 1.0, 0.0],      // Node 2: top-right corner at base
                    vec![0.0, 1.0, 0.0],     // Node 3: top-left corner at base
                    vec![0.0, 0.0, 1.0],       // Node 4: at apex
                ];
                let num_nodes = 5;

                let xi = 1.0/3.0;
                let eta = 1.0/3.0;
                let psi = 0.2113;

                let xim = 1.0 - xi;
                let etam = 1.0 - eta;
                let psim = 1.0 - psi;

                let values = vec![
                    xim * etam * psim,      // N0: base 
                    xi * etam * psim,       // N1: base 
                    xi * eta * psim,        // N2: base 
                    xim * eta * psim,       // N3: base 
                    psi,                    // N4: apex 
                ];
                let derivatives = vec![
                    vec![-etam * psim, -xim * psim, -xim * etam],   // dN0
                    vec![etam * psim, -xi * psim, -xi * etam],      // dN1
                    vec![eta * psim, xi * psim, -xi * eta],         // dN2
                    vec![-eta * psim, xim * psim, -xim * eta],      // dN3
                    vec![0.0, 0.0, 1.0],                            // dN4
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            /* 
            ElementType::Pyramid => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],    // Node 0: bottom-left corner at base
                ];
                let num_nodes = 5;

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let values = vec![
                    xim[0] * etam[0] * psim[0],      // N0: base 
                ];
                let derivatives = vec![
                    vec![],   // dN0
                ];

                Ok(ShapeFunction {values, derivatives})
            }, 
            */

            ElementType::Wedge => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],       // Node 0: bottom face, xxx
                    vec![1.0, 0.0, 0.0],       // Node 1: bottom face, xxx
                    vec![0.0, 1.0, 0.0],       // Node 2: bottom face, xxx
                    vec![0.0, 0.0, 1.0],        // Node 3: top face, xxx
                    vec![1.0, 0.0, 1.0],        // Node 4: top face, xxx
                    vec![0.0, 1.0, 1.0],        // Node 5: top face, xxx
                ];
                let num_nodes = 6;

                let xi = 1.0/3.0;
                let eta = 1.0/3.0;
                let psi = 0.2113;

                let values = vec![
                    (1.0 - xi - eta) * (1.0 - psi),     // N0: bottom
                    xi * (1.0 - psi),                   // N1: bottom
                    eta * (1.0 - psi),                  // N2: bottom
                    (1.0 - xi - eta) * psi,             // N3: top
                    xi * psi,                           // N4: top
                    eta * psi,                          // N5: top
                ];
                let derivatives = vec![
                    vec![-1.0 + psi, -1.0 + psi, -1.0 + xi + eta],  // dN0
                    vec![1.0 - psi, 0.0, -xi],                      // dN1
                    vec![0.0, 1.0 - psi, -eta],                     // dN2
                    vec![-psi, -psi, 1.0 - xi - eta],               // dN3
                    vec![psi, 0.0, xi],                             // dN4
                    vec![0.0, psi, eta],                            // dN5
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::QuadraticWedge => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],       // Node 0: bottom face, xxx
                    vec![1.0, 0.0, 0.0],       // Node 1: bottom face, xxx
                    vec![0.0, 1.0, 0.0],       // Node 2: bottom face, xxx
                    vec![0.0, 0.0, 1.0],       // Node 3: top face, xxx
                    vec![1.0, 0.0, 1.0],       // Node 4: top face, xxx
                    vec![0.0, 1.0, 1.0],       // Node 5: top face, xxx
                    // mid-edge nodes
                    vec![0.5, 0.0, 0.0],       // Node 6: bottom face, mid-edge 0-1
                    vec![0.5, 0.5, 0.0],       // Node 7: bottom face, mid-edge 1-2
                    vec![0.0, 0.5, 0.0],       // Node 8: bottom face, mid-edge 0-2
                    vec![0.5, 0.0, 1.0],       // Node 9: top face, mid-edge 3-4
                    vec![0.5, 0.5, 1.0],       // Node 10: top face, mid-edge 4-5
                    vec![0.0, 0.5, 1.0],       // Node 11: top face, mid-edge 3-5
                    vec![0.0, 0.0, 0.5],       // Node 12: mid-edge 0-3
                    vec![1.0, 0.0, 0.5],       // Node 13: mid-edge 1-4
                    vec![0.0, 1.0, 0.5],       // Node 14: mid-edge 2-5
                ];
                let num_nodes = 15;

                let xi = 0.1127;
                let eta = 0.5;
                let psi = 0.8873;

                let values = vec![
                    // Corner nodes
                    2.0 * (1.0 - xi - eta) * (1.0 - psi) * (0.5 - xi - eta - psi),    // N0: bottom corner
                    2.0 * xi * (1.0 - psi) * (xi - psi - 0.5),                      // N1: bottom corner
                    2.0 * eta * (1.0 - psi) * (eta - psi - 0.5),                    // N2: bottom corner
                    2.0 * (1.0 - xi - eta) * psi * (psi - xi - eta - 0.5),          // N3: top corner
                    2.0 * xi * psi * (xi + psi - 1.5),                            // N4: top corner
                    2.0 * eta * psi * (eta + psi - 1.5),                          // N5: top corner
                    // Mid-edge nodes on faces
                    4.0 * xi * (1.0 - xi - eta) * (1.0 - psi),                        // N6: bottom face
                    4.0 * xi * eta * (1.0 - psi),                                   // N7: bottom face
                    4.0 * (1.0 - xi - eta) * eta * (1.0 - psi),                       // N8: bottom face
                    4.0 * xi * (1.0 - xi - eta) * psi,                              // N9: top face
                    4.0 * xi * eta * psi,                                         // N10: top face
                    4.0 * (1.0 - xi - eta) * eta * psi,                             // N11: top face
                    // Vertical mid-edge nodes
                    4.0 * psi * (1.0 - xi - eta) * (1.0 - psi),                       // N12: vertical
                    4.0 * psi * xi * (1.0 - psi),                                   // N13: vertical
                    4.0 * psi * eta * (1.0 - psi),                                  // N14: vertical
                ];
                let derivatives = vec![
                    // Corner nodes
                    vec![2.0 * (1.0 - psi) * (-1.5 + 2.0 * xi + 2.0 * eta + psi), 2.0 * (1.0 - psi) * (-1.5 + 2.0 * xi + 2.0 * eta + psi), 2.0 * (1.0 - xi - eta) * (-1.5 + xi + eta + 2.0 * psi)], // dN0
                    vec![2.0 * (1.0 - psi) * (-0.5 + 2.0 * xi + 2.0 * eta - psi), 0.0, 2.0 * xi * (-0.5 - xi + 2.0 * psi)], // dN1
                    vec![0.0, 2.0 * (1.0 - psi) * (-0.5 + 2.0 * eta - psi), 2.0 * eta * (-0.5 - eta + 2.0 * psi)], // dN2
                    vec![2.0 * psi * (-0.5 + 2.0 * xi + 2.0 * eta - psi), 2.0 * psi * (-0.5 + 2.0 * xi + 2.0 * eta - psi), 2.0 * (1.0 - xi - eta) * (-0.5 - xi - eta + 2.0 * psi)], // dN3
                    vec![2.0 * psi * (-1.5 + 2.0 * xi + psi), 0.0, 2.0 * xi * (-1.5 + xi + 2.0 * psi)], // dN4
                    vec![0.0, 2.0 * psi * (-1.5 + 2.0 * eta + psi), 2.0 * eta * (-1.5 + eta + 2.0 * psi)], // dN5
                    // Mid-edge derivatives
                    vec![4.0 * (1.0 - psi) * (1.0 - 2.0 * xi - eta), -4.0 * (1.0 - psi) * xi, -4.0 * xi * (1.0 - xi - eta)], // dN6
                    vec![4.0 * (1.0 - psi) * eta, 4.0 * (1.0 - psi) * xi, -4.0 * xi * eta], // dN7
                    vec![-4.0 * (1.0 - psi) * eta, 4.0 * (1.0 - psi) * (1.0 - xi - 2.0 * eta), -4.0 * eta * (1.0 - xi - eta)], // dN8
                    vec![4.0 * psi * (1.0 - 2.0 * xi - eta), -4.0 * xi * psi, 4.0 * xi * (1.0 - xi - eta)], // dN9
                    vec![4.0 * eta * psi, 4.0 * xi * psi, 4.0 * xi * eta], // dN10
                    vec![-4.0 * eta * psi, 4.0 * psi * (1.0 - xi - 2.0 * eta), 4.0 * eta * (1.0 - xi - eta)], // dN11
                    // Vertical mid-edge derivatives
                    vec![-4.0 * psi * (1.0 - psi), -4.0 * psi * (1.0 - psi), 4.0 * (1.0 - 2.0 * psi) * (1.0 - xi - eta)], // dN12
                    vec![4.0 * psi * (1.0 - psi), 0.0, 4.0 * (1.0 - 2.0 * psi) * xi], // dN13
                    vec![0.0, 4.0 * psi * (1.0 - psi), 4.0 * (1.0 - 2.0 * psi) * eta], // dN14
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::BiquadraticQuadraticWedge => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],        // Node 0: bottom face, xxx
                    vec![1.0, 0.0, 0.0],        // Node 1: bottom face, xxx
                    vec![0.0, 1.0, 0.0],        // Node 2: bottom face, xxx
                    vec![0.0, 0.0, 1.0],        // Node 3: top face, xxx
                    vec![1.0, 0.0, 1.0],        // Node 4: top face, xxx
                    vec![0.0, 1.0, 1.0],        // Node 5: top face, xxx
                    // mid-edge nodes
                    vec![0.5, 0.0, 0.0],        // Node 6: bottom face, mid-edge 0-1
                    vec![0.5, 0.5, 0.0],        // Node 7: bottom face, mid-edge 1-2
                    vec![0.0, 0.5, 0.0],        // Node 8: bottom face, mid-edge 0-2
                    vec![0.5, 0.0, 1.0],        // Node 9: top face, mid-edge 3-4
                    vec![0.5, 0.5, 1.0],        // Node 10: top face, mid-edge 4-5
                    vec![0.0, 0.5, 1.0],        // Node 11: top face, mid-edge 3-5
                    vec![0.0, 0.0, 0.5],        // Node 12: mid-edge 0-3
                    vec![1.0, 0.0, 0.5],        // Node 13: mid-edge 1-4
                    vec![0.0, 1.0, 0.5],        // Node 14: mid-edge 2-5
                    // mid-face nodes
                    vec![0.5, 0.0, 0.5],        // Node 15: mid-edge 12-13
                    vec![0.5, 0.5, 0.5],        // Node 16: mid-edge 13-14
                    vec![0.0, 0.5, 0.5],        // Node 17: mid-edge 12-14
                ];
                let num_nodes = 18;

                let xi = 0.1127;
                let eta = 0.5;
                let psi = 0.8873;

                let xim = 2.0 * (xi - 0.5);
                let etam = 2.0 * (eta - 0.5);
                let psim = 2.0 * (psi - 0.5);

                let values = vec![
                    // Corner nodes
                    -0.25 * (xim + etam) * (xim + etam + 1.0) * psim * (1.0 - psim),    // N0: bottom corner
                    -0.25 * xim * (xim + 1.0) * psim * (1.0 - psim),                    // N1: bottom corner
                    -0.25 * etam * (1.0 + etam) * psim * (1.0 - psim),                  // N2: bottom corner
                    0.25 * (xim + etam) * (xim + etam + 1.0) * psim * (1.0 + psim),     // N3: top corner
                    0.25 * xim * (xim + 1.0) * psim * (1.0 + psim),                     // N4: top corner
                    0.25 * etam * (1.0 + etam) * psim * (1.0 + psim),                   // N5: top corner
                    // Mid-edge nodes on faces
                    (xim + 1.0) * (xim + etam) * 0.5 * psim * (1.0 - psim),                 // N6: bottom face
                    -(xim + 1.0) * (etam + 1.0) * 0.5 * psim * (1.0 - psim),                  // N7: bottom face
                    (etam + 1.0) * (xim + etam) * 0.5 * psim * (1.0 - psim),                // N8: bottom face
                    -(xim + 1.0) * (xim + etam) * 0.5 * psim * (1.0 + psim),                // N9: top face
                    (xim + 1.0) * (etam + 1.0) * 0.5 * psim * (1.0 + psim),                   // N10: top face
                    -(etam + 1.0) * (xim + etam) * 0.5 * psim * (1.0 + psim),               // N11: top face
                    // Vertical mid-edge nodes
                    0.5 * (xim + etam) * (xim + etam + 1.0) * (1.0 + psim) * (1.0 - psim),    // N12: vertical
                    0.5 * xim * (xim + 1.0) * (1.0 + psim) * (1.0 - psim),                    // N13: vertical
                    0.5 * etam * (1.0 + etam) * (1.0 + psim) * (1.0 - psim),                  // N14: vertical
                    // Mid-face nodes
                    -(xim + 1.0)*(xim + etam) * (1.0 + psim) * (1.0 - psim),                  // N15: mid-face
                    (xim + 1.0)*(etam + 1.0) * (1.0 + psim) * (1.0 - psim),                   // N16: mid-face
                    -(etam + 1.0)*(xim + etam) * (1.0 + psim) * (1.0 - psim),                 // N17: mid-face
                ];
                let derivatives = vec![
                    // we compute derivatives in [-1; 1] but we need them in [ 0; 1]
                        // for(int i = 0; i < 54; i++)
                            // derivs[i] *= 2;

                    // Corner derivatives
                    vec![2.0 * (-0.25 * (2.0 * xim + 2.0 * etam + 1.0) * psim * (1.0 - psim)), 
                        2.0 * (-0.25 * (2.0 * etam + 2.0 * xim + 1.0) * psim * (1.0 - psim)), 
                        2.0 * (-0.25 * (xim + etam) * (xim + etam + 1.0) * (1.0 - 2.0 * psim))],   // dN0
                    vec![2.0 * (-0.25 * (2.0 * xim + 1.0) * psim * (1.0 - psim)), 
                        2.0 * 0.0, 
                        2.0 * (-0.25 * xim * (xim + 1.0) * (1.0 - 2.0 * psim))],                   // dN1
                    vec![2.0 * 0.0, 
                        2.0 * (-0.25 * (2.0 * etam + 1.0) * psim * (1.0 - psim)), 
                        2.0 * (-0.25 * etam * (1.0 + etam) * (1.0 - 2.0 * psim))],                 // dN2
                    vec![2.0 * (0.25 * (2.0 * xim + 2.0 * etam + 1.0) * psim * (1.0 + psim)), 
                        2.0 * (0.25 * (2.0 * etam + 2.0 * xim + 1.0) * psim * (1.0 + psim)), 
                        2.0 * (0.25 * (xim + etam) * (xim + etam + 1.0) * (1.0 + 2.0 * psim))],    // dN3
                    vec![2.0 * (0.25 * (2.0 * xim + 1.0) * psim * (1.0 + psim)), 
                        2.0 * 0.0, 
                        2.0 * (0.25 * xim * (xim + 1.0) * (1.0 + 2.0 * psim))],                    // dN4
                    vec![2.0 * 0.0, 
                        2.0 * (0.25 * (2.0 * etam + 1.0) * psim * (1.0 + psim)), 
                        2.0 * (0.25 * etam * (1.0 + etam) * (1.0 + 2.0 * psim))],                  // dN5
                    // Face mid-edge derivatives
                    vec![2.0 * ((2.0 * xim + etam + 1.0) * 0.5 * psim * (1.0 - psim)), 
                        2.0 * ((xim + 1.0) * 0.5 * psim * (1.0 - psim)), 
                        2.0 * ((xim + 1.0) * (xim + etam) *  0.5 * (1.0 - 2.0 * psim))],           // dN6
                    vec![2.0 * (-(etam + 1.0) * 0.5 * psim * (1.0 - psim)), 
                        2.0 * (-(xim + 1.0) * 0.5 * psim * (1.0 - psim)), 
                        2.0 * (-(xim + 1.0) * (etam + 1.0) *  0.5 * (1.0 - 2.0 * psim))],          // dN7
                    vec![2.0 * ((etam + 1.0) * 0.5 * psim * (1.0 - psim)), 
                        2.0 * ((2.0 * etam + xim + 1.0) * 0.5 * psim * (1.0 - psim)), 
                        2.0 * ((etam + 1.0) * (xim + etam) *  0.5 * (1.0 - 2.0 * psim))],          // dN8
                    vec![2.0 * (-(2.0 * xim + etam + 1.0) * 0.5 * psim * (1.0 + psim)), 
                        2.0 * (-(xim + 1.0) * 0.5 * psim * (1.0 + psim)), 
                        2.0 * (-(xim + 1.0) * (xim + etam) *  0.5 * (1.0 + 2.0 * psim))],          // dN9
                    vec![2.0 * ((etam + 1.0) * 0.5 * psim * (1.0 + psim)), 
                        2.0 * ((xim + 1.0) * 0.5 * psim * (1.0 + psim)), 
                        2.0 * ((xim + 1.0) * (etam + 1.0) *  0.5 * (1.0 + 2.0 * psim))],           // dN10
                    vec![2.0 * (-(etam + 1.0) * 0.5 * psim * (1.0 + psim)), 
                        2.0 * (-(2.0 * etam + xim + 1.0) * 0.5 * psim * (1.0 + psim)), 
                        2.0 * (-(etam + 1.0) * (xim + etam) *  0.5 * (1.0 + 2.0 * psim))],         // dN11
                    // Vertical mid-edge derivatives
                    vec![2.0 * (0.5 * (2.0 * xim + 2.0 * etam + 1.0) * (1.0 + psim) * (1.0 - psim)),  
                        2.0 * (0.5 * (2.0 * etam + 2.0 * xim + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * (0.5 * (xim + etam) * (xim + etam + 1.0) * (-2.0 * psim))],          // dN12
                    vec![2.0 * (0.5 * (2.0 * xim + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * 0.0, 
                        2.0 * (0.5 * xim * (xim + 1.0) * (-2.0 * psim))],                          // dN13
                    vec![2.0 * 0.0, 
                        2.0 * (0.5 * (2.0 * etam + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * (0.5 * etam * (1.0 + etam) * (-2.0 * psim))],                        // dN14
                    // Mid-face derivatives
                    vec![2.0 * (-(2.0 * xim + etam + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * (-(xim + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * (-(xim + 1.0) * (xim + etam) * (-2.0 * psim))],                      // dN15
                    vec![2.0 * ((etam + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * ((xim + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * ((xim + 1.0) * (etam + 1.0) * (-2.0 * psim))],                       // dN16
                    vec![2.0 * (-(etam + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * (-(2.0 * etam + xim + 1.0) * (1.0 + psim) * (1.0 - psim)), 
                        2.0 * (-(etam + 1.0) * (xim + etam) * (-2.0 * psim))],                     // dN17
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::Hexahedron => {
                let reference_coords = vec![
                    vec![0.0, 0.0, 0.0],        // Node 0: bottom face, top-left
                    vec![1.0, 0.0, 0.0],        // Node 1: bottom face, bottom-left
                    vec![1.0, 1.0, 0.0],        // Node 2: bottom face, bottom-right
                    vec![0.0, 1.0, 0.0],        // Node 3: bottom face, top-right
                    vec![0.0, 0.0, 1.0],        // Node 4: top face, top-left
                    vec![1.0, 0.0, 1.0],        // Node 5: top face, bottom-left
                    vec![1.0, 1.0, 1.0],        // Node 6: top face, bottom-right
                    vec![0.0, 1.0, 1.0],        // Node 7: top face, top-right
                ];
                let num_nodes = 8;

                let xi = 0.2113;
                let eta = 0.2113;
                let psi = 0.7887;

                let xim = 1.0 - xi;
                let etam = 1.0 - eta;
                let psim = 1.0 - psi;

                let values = vec![
                    xim * etam * psim,      // N0: 
                    xi * etam * psim,       // N1: 
                    xi * eta * psim,        // N2: 
                    xim * eta * psim,       // N3: 
                    xim * etam * psi,       // N4: 
                    xi * etam * psi,        // N5: 
                    xi * eta * psi,         // N6: 
                    xim * eta * psi,       // N7: 
                ];
                let derivatives = vec![
                    vec![-etam * psim, -xim * psim, -xim * etam], // dN0
                    vec![etam * psim, -xi * psim, -xi * etam], // dN1
                    vec![eta * psim, xi * psim, -xi * eta], // dN2
                    vec![-eta * psim, xim * psim, -xim * eta], // dN3
                    vec![-etam * psi, -xim * psi, xim * etam], // dN4
                    vec![etam * psi, -xi * psi, xi * etam], // dN5
                    vec![eta * psi, xi * psi, xi * eta], // dN6
                    vec![-eta * psi, xim * psi, xim * eta], // dN7
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::QuadraticHexahedron => {
                let reference_coords = vec![
                    // corner nodes
                    vec![0.0, 0.0, 0.0],        // Node 0: bottom face, bottom-left
                    vec![1.0, 0.0, 0.0],        // Node 1: bottom face, bottom-right
                    vec![1.0, 1.0, 0.0],        // Node 2: bottom face, top-right
                    vec![0.0, 1.0, 0.0],        // Node 3: bottom face, top-left
                    vec![0.0, 0.0, 1.0],        // Node 4: top face, bottom-left
                    vec![1.0, 0.0, 1.0],        // Node 5: top face, bottom-right
                    vec![1.0, 1.0, 1.0],        // Node 6: top face, top-right
                    vec![0.0, 1.0, 1.0],        // Node 7: top face, top-left
                    // mid-edge nodes
                    vec![0.5, 0.0, 0.0],   // Node 8: bottom face, left edge
                    vec![1.0, 0.5, 0.0],    // Node 9: bottom face, front edge
                    vec![0.5, 1.0, 0.0],    // Node 10: bottom face, right edge
                    vec![0.0, 0.5, 0.0],   // Node 11: bottom face, back edge
                    vec![0.5, 0.0, 1.0],   // Node 12: top face, left edge
                    vec![1.0, 0.5, 1.0],     // Node 13: top face, front edge
                    vec![0.5, 1.0, 1.0],     // Node 14: top face, right edge
                    vec![0.0, 0.5, 1.0],    // Node 15: top face, back edge
                    vec![0.0, 0.0, 0.5],   // Node 16: left face, left edge
                    vec![1.0, 0.0, 0.5],    // Node 17: left face, right edge
                    vec![1.0, 1.0, 0.5],     // Node 18: right face, left edge
                    vec![0.0, 1.0, 0.5],    // Node 19: right face, right edge
                ];
                let num_nodes = 20;

                let xi = 0.1127;
                let eta = 0.5;
                let psi = 0.8873;

                let xik = 2.0 * (xi - 0.5);
                let etak = 2.0 * (eta - 0.5);
                let psik = 2.0 * (psi - 0.5);

                let xim = 1.0 - xik;
                let xip = 1.0 + xik;
                let etam = 1.0 - etak;
                let etap = 1.0 + etak;
                let psim = 1.0 - psik;
                let psip = 1.0 + psik;

                let xi2 = 1.0 - xik * xik;
                let eta2 = 1.0 - etak * etak;
                let psi2 = 1.0 - psik * psik;

                let values = vec![
                    // Corner nodes
                    0.125 * xim * etam * psim * (-xik - etak - psik - 2.0),     // N0: 
                    0.125 * xip * etam * psim * (xik - etak - psik - 2.0),      // N1: 
                    0.125 * xip * etap * psim * (xik + etak - psik - 2.0),      // N2: 
                    0.125 * xim * etap * psim * (-xik + etak - psik - 2.0),     // N3: 
                    0.125 * xim * etam * psip * (-xik - etak + psik - 2.0),     // N4: 
                    0.125 * xip * etam * psip * (xik - etak + psik - 2.0),      // N5: 
                    0.125 * xip * etap * psip * (xik + etak + psik - 2.0),      // N6: 
                    0.125 * xim * etap * psip * (-xik + etak + psik - 2.0),     // N7: 
                    // Mid-edge nodes
                    0.25 * xi2 * etam * psim,       // N8:  
                    0.25 * eta2 * xip * psim,       // N9:
                    0.25 * xi2 * etap * psim,       // N10:
                    0.25 * eta2 * xim * psim,       // N11:
                    0.25 * xi2 * etam * psip,       // N12:
                    0.25 * eta2 * xip * psip,       // N13:
                    0.25 * xi2 * etap * psip,       // N14:
                    0.25 * eta2 * xim * psip,       // N15:
                    0.25 * psi2 * xim * etam,       // N16:
                    0.25 * psi2 * xip * etam,       // N17:
                    0.25 * psi2 * xip * etap,       // N18:
                    0.25 * psi2 * xim * etap,       // N19:
                ];
                let derivatives = vec![
                    vec![2.0 * (-0.125 * (etam * psim - 2.0 * xik * etam * psim - etak * etam * psim - psik * etam * psim - 2.0 * etam * psim)), 
                        2.0 * (-0.125 * (xim * psim - 2.0 * etak * xim * psim - xik * xim * psim - psik * xim * psim - 2.0 * xim * psim)), 
                        2.0 * (-0.125 * (xim * etam - 2.0 * psik * xim * etam - xik * xim * etam - etak * xim * etam - 2.0 * xim * etam))],            // dN0
                    vec![2.0 * (0.125 * (etam * psim + 2.0 * xik * etam * psim - etak * etam * psim - psik * etam * psim - 2.0 * etam * psim)), 
                        2.0 * (-0.125 * (xip * psim - 2.0 * etak * xip * psim + xik * xip * psim - psik * xip * psim - 2.0 * xip * psim)), 
                        2.0 * (-0.125 * (xip * etam - 2.0 * psik * xip * etam + xik * xip * etam - etak * xip * etam - 2.0 * xip * etam))],            // dN1
                    vec![2.0 * (0.125 * (etap * psim + 2.0 * xik * etap * psim + etak * etap * psim - psik * etap * psim - 2.0 * etap * psim)), 
                        2.0 * (0.125 * (xip * psim + 2.0 * etak * xip * psim + xik * xip * psim - psik * xip * psim - 2.0 * xip * psim)), 
                        2.0 * (-0.125 * (xip * etap - 2.0 * psik * xip * etap + xik * xip * etap + etak * xip * etap - 2.0 * xip * etap))],            // dN2
                    vec![2.0 * (-0.125 * (etap * psim - 2.0 * xik * etap * psim + etak * etap * psim - psik * etap * psim - 2.0 * etap * psim)), 
                        2.0 * (0.125 * (xim * psim + 2.0 * etak * xim * psim - xik * xim * psim - psik * xim * psim - 2.0 * xim * psim)), 
                        2.0 * (-0.125 * (xim * etap - 2.0 * psik * xim * etap - xik * xim * etap + etak * xim * etap - 2.0 * xim * etap))],            // dN3
                    vec![2.0 * (-0.125 * (etam * psip - 2.0 * xik * etam * psip - etak * etam * psip + psik * etam * psip - 2.0 * etam * psip)), 
                        2.0 * (-0.125 * (xim * psip - 2.0 * etak * xim * psip - xik * xim * psip + psik * xim * psip - 2.0 * xim * psip)), 
                        2.0 * (0.125 * (xim * etam + 2.0 * psik * xim * etam - xik * xim * etam - etak * xim * etam - 2.0 * xim * etam))],             // dN4
                    vec![2.0 * (0.125 * (etam * psip + 2.0 * xik * etam * psip - etak * etam * psip + psik * etam * psip - 2.0 * etam * psip)), 
                        2.0 * (-0.125 * (xip * psip - 2.0 * etak * xip * psip + xik * xip * psip - psik * xip * psip - 2.0 * xip * psip)), 
                        2.0 * (0.125 * (xip * etam + 2.0 * psik * xip * etam + xik * xip * etam - etak * xip * etam - 2.0 * xip * etam))],             // dN5
                    vec![2.0 * (0.125 * (etap * psip + 2.0 * xik * etap * psip + etak * etap * psip + psik * etap * psip - 2.0 * etap * psip)), 
                        2.0 * (0.125 * (xip * psip + 2.0 * etak * xip * psip + xik * xip * psip + psik * xip * psip - 2.0 * xip * psip)), 
                        2.0 * (0.125 * (xip * etap + 2.0 * psik * xip * etap + xik * xip * etap + etak * xip * etap - 2.0 * xip * etap))],             // dN6
                    vec![2.0 * (-0.125 * (etap * psip - 2.0 * xik * etap * psip + etak * etap * psip + psik * etap * psip - 2.0 * etap * psip)), 
                        2.0 * (0.125 * (xim * psip + 2.0 * etak * xim * psip - xik * xim * psip - psik * xim * psip - 2.0 * xim * psip)), 
                        2.0 * (0.125 * (xim * etap + 2.0 * psik * xim * etap - xik * xim * etap + etak * xim * etap - 2.0 * xim * etap))],             // dN7
                    vec![2.0 * (-0.5 * xik * etam * psim), 2.0 * (-0.25 * (psim - xik * xik * psim)), 2.0 * (-0.25 * (etam - xik * xik * etam))],       // dN8
                    vec![2.0 * (0.25 * (psim - etak * etak * psim)), 2.0 * (-0.5 * etak * xip * psim), 2.0 * (-0.25 * (xip - etak * etak * xip))],      // dN9
                    vec![2.0 * (-0.5 * xik * etap * psim), 2.0 * (0.25 * (psim - xik * xik * psim)), 2.0 * (-0.25 * (etap - xik * xik * etap))],        // dN10
                    vec![2.0 * (-0.25 * (psim - etak * etak * psim)), 2.0 * (-0.5 * etak * xim * psim), 2.0 * (-0.25 * (xim - etak * etak * xim))],     // dN11
                    vec![2.0 * (-0.5 * xik * etam * psip), 2.0 * (-0.25 * (psip - xik * xik * psip)), 2.0 * (0.25 * (etam - xik * xik * etam))],        // dN12
                    vec![2.0 * (0.25 * (psip - etak * etak * psip)), 2.0 * (-0.5 * etak * xip * psip), 2.0 * (0.25 * (xip - etak * etak * xip))],       // dN13
                    vec![2.0 * (-0.5 * xik * etap * psip), 2.0 * (0.25 * (psip - xik * xik * psip)), 2.0 * (0.25 * (etap - xik * xik * etap))],         // dN14
                    vec![2.0 * (-0.25 * (psip - etak * etak * psip)), 2.0 * (-0.5 * etak * xim * psip), 2.0 * (-0.25 * (xim - etak * etak * xim))],      // dN15
                    vec![2.0 * (-0.25 * (etam - psik * psik * etam)), 2.0 * (-0.25 * (xim - psik * psik * xim)), 2.0 * (-0.5 * psik * xim * etam)],     // dN16
                    vec![2.0 * (0.25 * (etam - psik * psik * etam)), 2.0 * (-0.25 * (xip - psik * psik * xip)), 2.0 * (-0.5 * psik * xip * etam)],      // dN17
                    vec![2.0 * (0.25 * (etap - psik * psik * etap)), 2.0 * (0.25 * (xip - psik * psik * xip)), 2.0 * (-0.5 * psik * xip * etap)],       // dN18
                    vec![2.0 * (-0.25 * (etap - psik * psik * etap)), 2.0 * (0.25 * (xim - psik * psik * xim)), 2.0 * (-0.5 * psik * xim * etap)],      // dN19
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::BiquadraticQuadraticHexahedron => {
                let reference_coords = vec![
                    // corner nodes
                    vec![0.0, 0.0, 0.0],        // Node 0: bottom face, bottom-left
                    vec![1.0, 0.0, 0.0],        // Node 1: bottom face, bottom-right
                    vec![1.0, 1.0, 0.0],        // Node 2: bottom face, top-right
                    vec![0.0, 1.0, 0.0],        // Node 3: bottom face, top-left
                    vec![0.0, 0.0, 1.0],        // Node 4: top face, bottom-left
                    vec![1.0, 0.0, 1.0],        // Node 5: top face, bottom-right
                    vec![1.0, 1.0, 1.0],        // Node 6: top face, top-right
                    vec![0.0, 1.0, 1.0],        // Node 7: top face, top-left
                    // mid-edge nodes
                    vec![0.5, 0.0, 0.0],   // Node 8: bottom face, left edge
                    vec![1.0, 0.5, 0.0],    // Node 9: bottom face, front edge
                    vec![0.5, 1.0, 0.0],    // Node 10: bottom face, right edge
                    vec![0.0, 0.5, 0.0],   // Node 11: bottom face, back edge
                    vec![0.5, 0.0, 1.0],   // Node 12: top face, left edge
                    vec![1.0, 0.5, 1.0],     // Node 13: top face, front edge
                    vec![0.5, 1.0, 1.0],     // Node 14: top face, right edge
                    vec![0.0, 0.5, 1.0],    // Node 15: top face, back edge
                    vec![0.0, 0.0, 0.5],   // Node 16: left face, left edge
                    vec![1.0, 0.0, 0.5],    // Node 17: left face, right edge
                    vec![1.0, 1.0, 0.5],     // Node 18: right face, left edge
                    vec![0.0, 1.0, 0.5],    // Node 19: right face, right edge
                    // mid-face nodes
                    vec![0.0, 0.5, 0.5],    // Node 20: back face, mid
                    vec![1.0, 0.5, 0.5],     // Node 21: front face, mid
                    vec![0.5, 0.0, 0.5],    // Node 22: left face, mid
                    vec![0.5, 1.0, 0.5],     // Node 23: right face, mid
                ];
                let num_nodes = 24;

                let xi = 0.1127;
                let eta = 0.5;
                let psi = 0.8873;

                let xik = 2.0 * (xi - 0.5);
                let etak = 2.0 * (eta - 0.5);
                let psik = 2.0 * (psi - 0.5);

                let values = vec![
                    // Corner nodes
                    (0.25 * (xik * (1.0 - xik)) * (etak * (1.0 - etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (-0.5 * psik * (1.0 - psik)),       // N0:
                    (-0.25 * (xik * (1.0 + xik)) * (etak * (1.0 - etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (-0.5 * psik * (1.0 - psik)),      // N1:
                    (0.25 * (xik * (1.0 + xik)) * (etak * (1.0 + etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (-0.5 * psik * (1.0 - psik)),       // N2:
                    (-0.25 * (xik * (1.0 - xik)) * (etak * (1.0 + etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (-0.5 * psik * (1.0 - psik)),      // N3:
                    (0.25 * (xik * (1.0 - xik)) * (etak * (1.0 - etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (0.5 * psik * (1.0 + psik)),        // N4:
                    (-0.25 * (xik * (1.0 + xik)) * (etak * (1.0 - etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (0.5 * psik * (1.0 + psik)),       // N5:
                    (0.25 * (xik * (1.0 + xik)) * (etak * (1.0 + etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (0.5 * psik * (1.0 + psik)),        // N6:
                    (-0.25 * (xik * (1.0 - xik)) * (etak * (1.0 + etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * (0.5 * psik * (1.0 + psik)),       // N7:
                    // Mid-edge nodes
                    0.5 * ((1.0 + xik) * (1.0 - xik)) * (1.0 - etak) * (-0.5 * psik * (1.0 - psik)),        // N8:
                    0.5 * ((1.0 + etak) * (1.0 - etak)) * (1.0 + xik) * (-0.5 * psik * (1.0 - psik)),       // N9:
                    0.5 * ((1.0 + xik) * (1.0 - xik)) * (1.0 + etak) * (-0.5 * psik * (1.0 - psik)),        // N10:
                    0.5 * ((1.0 + etak) * (1.0 - etak)) * (1.0 - xik) * (-0.5 * psik * (1.0 - psik)),       // N11:
                    0.5 * ((1.0 + xik) * (1.0 - xik)) * (1.0 - etak) * (0.5 * psik * (1.0 + psik)),         // N12:
                    0.5 * ((1.0 + etak) * (1.0 - etak)) * (1.0 + xik) * (0.5 * psik * (1.0 + psik)),        // N13:
                    0.5 * ((1.0 + xik) * (1.0 - xik)) * (1.0 + etak) * (0.5 * psik * (1.0 + psik)),         // N14:
                    0.5 * ((1.0 + etak) * (1.0 - etak)) * (1.0 - xik) * (0.5 * psik * (1.0 + psik)),        // N15:
                    (0.25 * (xik * (1.0 - xik)) * (etak * (1.0 - etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * ((1.0 + psik) * (1.0 - psik)),        // N16:
                    (-0.25 * (xik * (1.0 + xik)) * (etak * (1.0 - etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * ((1.0 + psik) * (1.0 - psik)),       // N17:
                    (0.25 * (xik * (1.0 + xik)) * (etak * (1.0 + etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * ((1.0 + psik) * (1.0 - psik)),        // N18:
                    (-0.25 * (xik * (1.0 - xik)) * (etak * (1.0 + etak)) - 0.25 * (1.0 + xik) * (1.0 - xik) * (1.0 + etak) * (1.0 - etak)) * ((1.0 + psik) * (1.0 - psik)),       // N19:
                    // Mid-face nodes
                    0.5 * ((1.0 + etak) * (1.0 - etak)) * (1.0 - xik) * ((1.0 + psik) * (1.0 - psik)),        // N20:
                    0.5 * ((1.0 + etak) * (1.0 - etak)) * (1.0 + xik) * ((1.0 + psik) * (1.0 - psik)),        // N21:
                    0.5 * ((1.0 + xik) * (1.0 - xik)) * (1.0 - etak) * ((1.0 + psik) * (1.0 - psik)),         // N22:
                    0.5 * ((1.0 + xik) * (1.0 - xik)) * (1.0 + etak) * ((1.0 + psik) * (1.0 - psik)),         // N23:
                ];
                let derivatives = vec![
                    vec![2.0 * (-((etak * etak + (2.0 * xik - 1.0) * etak - 2.0 * xik) * psik * psik + (-etak * etak + (1.0 - 2.0 * xik) * etak + 2.0 * xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak + xik * xik - xik) * psik * psik + ((2.0 - 2.0 * xik) * etak - xik * xik + xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak * etak + (2.0 * xik * xik - 2.0 * xik) * etak - 2.0 * xik * xik + 2.0) * psik + (1.0 - xik) * etak * etak + (xik - xik * xik) * etak + xik * xik - 1.0) / 8.0)],     // dN0
                    vec![2.0 * (((etak * etak + (-2.0 * xik - 1.0) * etak + 2.0 * xik) * psik * psik + (-etak * etak + (2.0 * xik + 1.0) * etak - 2.0 * xik) * psik) / 8.0), 
                        2.0 * ((((2.0 * xik + 2.0) * etak - xik * xik - xik) * psik * psik + ((-2.0 * xik - 2.0) * etak + xik * xik + xik) * psik) / 8.0),
                        2.0 * ((((2.0 * xik + 2.0) * etak * etak + (-2.0 * xik * xik - 2.0 * xik) * etak + 2.0 * xik * xik - 2.0) * psik + (-xik - 1.0) * etak * etak + (xik * xik + xik) * etak - xik * xik + 1.0) / 8.0)],      // dN1
                    vec![2.0 * (((etak * etak + (2.0 * xik + 1.0) * etak + 2.0 * xik) * psik * psik + (-etak * etak + (-2.0 * xik - 1.0) * etak - 2.0 * xik) * psik) / 8.0), 
                        2.0 * ((((2.0 * xik + 2.0) * etak + xik * xik + xik) * psik * psik + ((-2.0 * xik - 2.0) * etak + xik * xik + xik) * psik) / 8.0),
                        2.0 * ((((2.0 * xik + 2.0) * etak * etak + (2.0 * xik * xik + 2.0 * xik) * etak + 2.0 * xik * xik - 2.0) * psik + (-xik - 1.0) * etak * etak + (-xik * xik - xik) * etak - xik * xik + 1.0) / 8.0)],       // dN2
                    vec![2.0 * (-((etak * etak + (1.0 - 2.0 * xik) * etak - 2.0 * xik) * psik * psik + (-etak * etak + (2.0 * xik - 1.0) * etak + 2.0 * xik) * psik) / 8.0),
                        2.0 * (-(((2.0 * xik - 2.0) * etak - xik * xik + xik) * psik * psik + ((2.0 - 2.0 * xik) * etak + xik * xik - xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak * etak + (2.0 * xik - 2.0 * xik * xik) * etak - 2.0 * xik * xik + 2.0) * psik + (1.0 - xik) * etak * etak + (xik * xik - xik) * etak + xik * xik - 1.0) / 8.0)],      // dN3
                    vec![2.0 * (-((etak * etak + (2.0 * xik - 1.0) * etak - 2.0 * xik) * psik * psik + (etak * etak + (2.0 * xik - 1.0) * etak - 2.0 * xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak + xik * xik - xik) * psik * psik + ((2.0 * xik - 2.0) * etak + xik * xik - xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak * etak + (2.0 * xik * xik - 2.0 * xik) * etak - 2.0 * xik * xik + 2.0) * psik + (xik - 1.0) * etak * etak + (xik * xik - xik) * etak - xik * xik + 1.0) / 8.0)],      // dN4
                    vec![2.0 * (((etak * etak + (-2.0 * xik - 1.0) * etak + 2.0 * xik) * psik * psik + (etak * etak + (-2.0 * xik - 1.0) * etak + 2.0 * xik) * psik) / 8.0), 
                        2.0 * ((((2.0 * xik + 2.0) * etak - xik * xik - xik) * psik * psik + ((2.0 * xik + 2.0) * etak - xik * xik - xik) * psik) / 8.0), 
                        2.0 * ((((2.0 * xik + 2.0) * etak * etak + (-2.0 * xik * xik - 2.0 * xik) * etak + 2.0 * xik * xik - 2.0) * psik + (xik + 1.0) * etak * etak + (-xik * xik - xik) * etak + xik * xik - 1.0) / 8.0)],       // dN5
                    vec![2.0 * (((etak * etak + (2.0 * xik + 1.0) * etak + 2.0 * xik) * psik * psik + (etak * etak + (2.0 * xik + 1.0) * etak + 2.0 * xik) * psik) / 8.0), 
                        2.0 * ((((2.0 * xik + 2.0) * etak + xik * xik + xik) * psik * psik + ((2.0 * xik + 2.0) * etak + xik * xik + xik) * psik) / 8.0), 
                        2.0 * ((((2.0 * xik + 2.0) * etak * etak + (2.0 * xik * xik + 2.0 * xik) * etak + 2.0 * xik * xik - 2.0) * psik + (xik + 1.0) * etak * etak + (xik * xik + xik) * etak + xik * xik - 1.0) / 8.0)],        // dN6
                    vec![2.0 * (-((etak * etak + (1.0 - 2.0 * xik) * etak - 2.0 * xik) * psik * psik + (etak * etak + (1.0 - 2.0 * xik) * etak - 2.0 * xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak - xik * xik + xik) * psik * psik + ((2.0 * xik - 2.0) * etak - xik * xik + xik) * psik) / 8.0), 
                        2.0 * (-(((2.0 * xik - 2.0) * etak * etak + (2.0 * xik - 2.0 * xik * xik) * etak - 2.0 * xik * xik + 2.0) * psik + (1.0 - xik) * etak * etak + (xik - xik * xik) * etak - xik * xik + 1.0) / 8.0)],       // dN7
                    vec![2.0 * (((xik * etak - xik) * psik * psik + (xik - xik * etak) * psik) / 2.0), 
                        2.0 * (((xik * xik - 1.0) * psik * psik + (1.0 - xik * xik) * psik) / 4.0), 
                        2.0 * ((((2.0 * xik * xik - 2.0) * etak - 2.0 * xik * xik + 2.0) * psik + (1.0 - xik * xik) * etak + xik * xik - 1.0) / 4.0)],       // dN8
                    vec![2.0 * (-((etak * etak - 1.0) * psik * psik + (1.0 - etak * etak) * psik) / 4.0), 
                        2.0 * (-((xik + 1.0) * etak * psik * psik + (-xik - 1.0) * etak * psik) / 2.0), 
                        2.0 * (-(((2.0 * xik + 2.0) * etak * etak - 2.0 * xik - 2.0) * psik + (-xik - 1.0) * etak * etak + xik + 1.0) / 4.0)],      // dN9
                    vec![2.0 * (-((xik * etak + xik) * psik * psik + (-xik * etak - xik) * psik) / 2.0), 
                        2.0 * (-((xik * xik - 1.0) * psik * psik + (1.0 - xik * xik) * psik) / 4.0), 
                        2.0 * (-(((2.0 * xik * xik - 2.0) * etak + 2.0 * xik * xik - 2.0) * psik + (1.0 - xik * xik) * etak - xik * xik + 1.0) / 4.0)],        // dN10
                    vec![2.0 * (((etak * etak - 1.0) * psik * psik + (1.0 - etak * etak) * psik) / 4.0), 
                        2.0 * (((xik - 1.0) * etak * psik * psik + (1.0 - xik) * etak * psik) / 2.0), 
                        2.0 * ((((2.0 * xik - 2.0) * etak * etak - 2.0 * xik + 2.0) * psik + (1.0 - xik) * etak * etak + xik - 1.0) / 4.0)],     // dN11
                    vec![2.0 * (((xik * etak - xik) * psik * psik + (xik * etak - xik) * psik) / 2.0), 
                        2.0 * (((xik * xik - 1.0) * psik * psik + (xik * xik - 1.0) * psik) / 4.0), 
                        2.0 * ((((2.0 * xik * xik - 2.0) * etak - 2.0 * xik * xik + 2.0) * psik + (xik * xik - 1.0) * etak - xik * xik + 1.0) / 4.0)],        // dN12
                    vec![2.0 * (-((etak * etak - 1.0) * psik * psik + (etak * etak - 1.0) * psik) / 4.0), 
                        2.0 * (-((xik + 1.0) * etak * psik * psik + (xik + 1.0) * etak * psik) / 2.0), 
                        2.0 * (-(((2.0 * xik + 2.0) * etak * etak - 2.0 * xik - 2.0) * psik + (xik + 1.0) * etak * etak - xik - 1.0) / 4.0)],       // dN13
                    vec![2.0 * (-((xik * etak + xik) * psik * psik + (xik * etak + xik) * psik) / 2.0), 
                        2.0 * (-((xik * xik - 1.0) * psik * psik + (xik * xik - 1.0) * psik) / 4.0), 
                        2.0 * (-(((2.0 * xik * xik - 2.0) * etak + 2.0 * xik * xik - 2.0) * psik + (xik * xik - 1.0) * etak + xik * xik - 1.0) / 4.0)],         // dN14
                    vec![2.0 * (((etak * etak + (2.0 * xik - 1.0) * etak - 2.0 * xik) * psik * psik - etak * etak + (1.0 - 2.0 * xik) * etak + 2.0 * xik) / 4.0), 
                        2.0 * ((((2.0 * xik - 2.0) * etak + xik * xik - xik) * psik * psik + (2.0 - 2.0 * xik) * etak - xik * xik + xik) / 4.0), 
                        2.0 * (((xik - 1.0) * etak * etak + (xik * xik - xik) * etak - xik * xik + 1.0) * psik / 2.0)],     // dN15
                    vec![2.0 * (-((etak * etak + (-2.0 * xik - 1.0) * etak + 2.0 * xik) * psik * psik - etak * etak + (2.0 * xik + 1.0) * etak - 2.0 * xik) / 4.0), 
                        2.0 * (-(((2.0 * xik + 2.0) * etak - xik * xik - xik) * psik * psik + (-2.0 * xik - 2.0) * etak + xik * xik + xik) / 4.0), 
                        2.0 * (-((xik + 1.0) * etak * etak + (-xik * xik - xik) * etak + xik * xik - 1.0) * psik / 2.0)],      // dN16
                    vec![2.0 * (-((etak * etak + (2.0 * xik + 1.0) * etak + 2.0 * xik) * psik * psik - etak * etak + (-2.0 * xik - 1.0) * etak - 2.0 * xik) / 4.0), 
                        2.0 * (-(((2.0 * xik + 2.0) * etak + xik * xik + xik) * psik * psik + (-2.0 * xik - 2.0) * etak - xik * xik - xik) / 4.0), 
                        2.0 * (-((xik + 1.0) * etak * etak + (xik * xik + xik) * etak + xik * xik - 1.0) * psik / 2.0)],       // dN17
                    vec![2.0 * (((etak * etak + (1.0 - 2.0 * xik) * etak - 2.0 * xik) * psik * psik - etak * etak + (2.0 * xik - 1.0) * etak + 2.0 * xik) / 4.0), 
                        2.0 * ((((2.0 * xik - 2.0) * etak - xik * xik + xik) * psik * psik + (2.0 - 2.0 * xik) * etak + xik * xik - xik) / 4.0), 
                        2.0 * (((xik - 1.0) * etak * etak + (xik - xik * xik) * etak - xik * xik + 1.0) * psik / 2.0)],      // dN18
                    vec![2.0 * (-((etak * etak - 1.0) * psik * psik - etak * etak + 1.0) / 2.0), 
                        2.0 * ((1.0 - xik) * etak * psik * psik + (xik - 1.0) * etak), 
                        2.0 * (((1.0 - xik) * etak * etak + xik - 1.0) * psik)],       // dN19
                    vec![2.0 * (((etak * etak - 1.0) * psik * psik - etak * etak + 1.0) / 2.0), 
                        2.0 * ((xik + 1.0) * etak * psik * psik + (-xik - 1.0) * etak), 
                        2.0 * (((xik + 1.0) * etak * etak - xik - 1.0) * psik)],      // dN20
                    vec![2.0 * ((xik - xik * etak) * psik * psik + xik * etak - xik), 
                        2.0 * (-((xik * xik - 1.0) * psik * psik - xik * xik + 1.0) / 2.0), 
                        2.0 * (((1.0 - xik * xik) * etak + xik * xik - 1.0) * psik)],       // dN21
                    vec![2.0 * ((xik * etak + xik) * psik * psik - xik * etak - xik), 
                        2.0 * (((xik * xik - 1.0) * psik * psik - xik * xik + 1.0) / 2.0), 
                        2.0 * (((xik * xik - 1.0) * etak + xik * xik - 1.0) * psik)],      // dN22
                    vec![2.0 * ((xik * etak + xik) * psik * psik - xik * etak - xik), 
                        2.0 * (((xik * xik - 1.0) * psik * psik - xik * xik + 1.0) / 2.0), 
                        2.0 * (((xik * xik - 1.0) * etak + xik * xik - 1.0) * psik)],      // dN23
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            ElementType::TriquadraticHexahedron => {
                let reference_coords = vec![
                    // corner nodes
                    vec![0.0, 0.0, 0.0],        // Node 0: bottom face, bottom-left
                    vec![1.0, 0.0, 0.0],        // Node 1: bottom face, bottom-right
                    vec![1.0, 1.0, 0.0],        // Node 2: bottom face, top-right
                    vec![0.0, 1.0, 0.0],        // Node 3: bottom face, top-left
                    vec![0.0, 0.0, 1.0],        // Node 4: top face, bottom-left
                    vec![1.0, 0.0, 1.0],        // Node 5: top face, bottom-right
                    vec![1.0, 1.0, 1.0],        // Node 6: top face, top-right
                    vec![0.0, 1.0, 1.0],        // Node 7: top face, top-left
                    // mid-edge nodes
                    vec![0.5, 0.0, 0.0],   // Node 8: bottom face, left edge
                    vec![1.0, 0.5, 0.0],    // Node 9: bottom face, front edge
                    vec![0.5, 1.0, 0.0],    // Node 10: bottom face, right edge
                    vec![0.0, 0.5, 0.0],   // Node 11: bottom face, back edge
                    vec![0.5, 0.0, 1.0],   // Node 12: top face, left edge
                    vec![1.0, 0.5, 1.0],     // Node 13: top face, front edge
                    vec![0.5, 1.0, 1.0],     // Node 14: top face, right edge
                    vec![0.0, 0.5, 1.0],    // Node 15: top face, back edge
                    vec![0.0, 0.0, 0.5],   // Node 16: left face, left edge
                    vec![1.0, 0.0, 0.5],    // Node 17: left face, right edge
                    vec![1.0, 1.0, 0.5],     // Node 18: right face, left edge
                    vec![0.0, 1.0, 0.5],    // Node 19: right face, right edge
                    // mid-face nodes
                    vec![0.0, 0.5, 0.5],    // Node 20: back face, mid
                    vec![1.0, 0.5, 0.5],     // Node 21: front face, mid
                    vec![0.5, 0.0, 0.5],    // Node 22: left face, mid
                    vec![0.5, 1.0, 0.5],     // Node 23: right face, mid
                    vec![0.5, 0.5, 0.0],    // Node 24: bottom face, mid
                    vec![0.5, 0.5, 1.0],     // Node 25: top face, mid
                    // center node
                    vec![0.5, 0.5, 0.5],     // Node 26: center 
                ];
                let num_nodes = 27;

                // VTK needs parametric coordinates to be between (0,1). Isoparametric
                // shape functions are formulated between (-1,1). Here we do a
                // coordinate system conversion from (0,1) to (-1,1).

                let xi = 0.1127;
                let eta = 0.5;
                let psi = 0.8873;

                let xik = 2.0 * (xi - 0.5);
                let etak = 2.0 * (eta - 0.5);
                let psik = 2.0 * (psi - 0.5);

                let xi1 = -0.5 * xik * (1.0 - xik);
                let eta1 = -0.5 * etak * (1.0 - etak);
                let psi1 = -0.5 * psik * (1.0 - psik);

                let xi2 = (1.0 + xik) * (1.0 - xik);
                let eta2 = (1.0 + etak) * (1.0 - etak);
                let psi2 = (1.0 + psik) * (1.0 - psik);

                let xi3 = 0.5 * xik * (1.0 + xik);
                let eta3 = 0.5 * etak * (1.0 + etak);
                let psi3 = 0.5 * psik * (1.0 + psik);

                let xi1_xi = xik - 0.5;
                let eta1_eta = etak - 0.5;
                let psi1_psi = psik - 0.5;

                let xi2_xi = -2.0 * xik;
                let eta2_eta = -2.0 * etak;
                let psi2_psi = -2.0 * psik;

                let xi3_xi = 0.5 * xik + 0.5;
                let eta3_eta = 0.5 * etak + 0.5;
                let psi3_psi = 0.5 * psik + 0.5;

                let values = vec![
                    // Corner nodes
                    xi1 * eta1 * psi1,      // N0:
                    xi3 * eta1 * psi1,      // N1:
                    xi3 * eta3 * psi1,      // N2:
                    xi1 * eta3 * psi1,      // N3:
                    xi1 * eta1 * psi3,      // N4:
                    xi3 * eta1 * psi3,      // N5:
                    xi3 * eta3 * psi3,      // N6:
                    xi1 * eta3 * psi3,      // N7:
                    // Mid-edge nodes
                    xi2 * eta1 * psi1,      // N8:
                    xi3 * eta2 * psi1,      // N9:
                    xi2 * eta3 * psi1,  // N10:
                    xi1 * eta2 * psi1,  // N11:
                    xi2 * eta1 * psi3,  // N12:
                    xi3 * eta2 * psi3,  // N13:
                    xi2 * eta3 * psi3,  // N14:
                    xi1 * eta2 * psi3,  // N15:
                    xi1 * eta1 * psi2,  // N16:
                    xi3 * eta1 * psi2,  // N17:
                    xi3 * eta3 * psi2,  // N18:
                    xi1 * eta3 * psi2,  // N19:
                    // Mid-face nodes
                    xi2 * eta1 * psi2,  // N20:
                    xi3 * eta2 * psi2,  // N21:
                    xi2 * eta3 * psi2,  // N22:
                    xi1 * eta2 * psi2,  // N23:
                    xi2 * eta2 * psi1,  // N24:
                    xi2 * eta2 * psi3,  // N25:
                    // Center node
                    xi2 * eta2 * psi2,  // N26:
                ];
                let derivatives = vec![
                    vec![2.0 * (xi1_xi * eta1 * psi1), 2.0 * (xi1 * eta1_eta * psi1), 2.0 * (xi1 * eta1 * psi1_psi)],     // dN0
                    vec![2.0 * (xi3_xi * eta1 * psi1), 2.0 * (xi3 * eta1_eta * psi1), 2.0 * (xi3 * eta1 * psi1_psi)],      // dN1
                    vec![2.0 * (xi3_xi * eta3 * psi1), 2.0 * (xi3 * eta3_eta * psi1), 2.0 * (xi3 * eta3 * psi1_psi)],       // dN2
                    vec![2.0 * (xi1_xi * eta3 * psi1), 2.0 * (xi1 * eta3_eta * psi1), 2.0 * (xi1 * eta3 * psi1_psi)],      // dN3
                    vec![2.0 * (xi1_xi * eta1 * psi3), 2.0 * (xi1 * eta1_eta * psi3), 2.0 * (xi1 * eta1 * psi3_psi)],      // dN4
                    vec![2.0 * (xi3_xi * eta1 * psi3), 2.0 * (xi3 * eta1_eta * psi3), 2.0 * (xi3 * eta1 * psi3_psi)],       // dN5
                    vec![2.0 * (xi3_xi * eta3 * psi3), 2.0 * (xi3 * eta3_eta * psi3), 2.0 * (xi3 * eta3 * psi3_psi)],        // dN6
                    vec![2.0 * (xi1_xi * eta3 * psi3), 2.0 * (xi1 * eta3_eta * psi3), 2.0 * (xi1 * eta3 * psi3_psi)],       // dN7
                    vec![2.0 * (xi2_xi * eta1 * psi1), 2.0 * (xi2 * eta1_eta * psi1), 2.0 * (xi2 * eta1 * psi1_psi)],       // dN8
                    vec![2.0 * (xi3_xi * eta2 * psi1), 2.0 * (xi3 * eta2_eta * psi1), 2.0 * (xi3 * eta2 * psi1_psi)],      // dN9
                    vec![2.0 * (xi2_xi * eta3 * psi1), 2.0 * (xi2 * eta3_eta * psi1), 2.0 * (xi2 * eta3 * psi1_psi)],        // dN10
                    vec![2.0 * (xi1_xi * eta2 * psi1), 2.0 * (xi1 * eta2_eta * psi1), 2.0 * (xi1 * eta2 * psi1_psi)],     // dN11
                    vec![2.0 * (xi2_xi * eta1 * psi3), 2.0 * (xi2 * eta1_eta * psi3), 2.0 * (xi2 * eta1 * psi3_psi)],        // dN12
                    vec![2.0 * (xi3_xi * eta2 * psi3), 2.0 * (xi3 * eta2_eta * psi3), 2.0 * (xi3 * eta2 * psi3_psi)],       // dN13
                    vec![2.0 * (xi2_xi * eta3 * psi3), 2.0 * (xi2 * eta3_eta * psi3), 2.0 * (xi2 * eta3 * psi3_psi)],         // dN14
                    vec![2.0 * (xi1_xi * eta2 * psi3), 2.0 * (xi1 * eta2_eta * psi3), 2.0 * (xi1 * eta2 * psi3_psi)],      // dN15
                    vec![2.0 * (xi1_xi * eta1 * psi2), 2.0 * (xi1 * eta1_eta * psi2), 2.0 * (xi1 * eta1 * psi2_psi)],     // dN16
                    vec![2.0 * (xi3_xi * eta1 * psi2), 2.0 * (xi3 * eta1_eta * psi2), 2.0 * (xi3 * eta1 * psi2_psi)],      // dN17
                    vec![2.0 * (xi3_xi * eta3 * psi2), 2.0 * (xi3 * eta3_eta * psi2), 2.0 * (xi3 * eta3 * psi2_psi)],       // dN18
                    vec![2.0 * (xi1_xi * eta3 * psi2), 2.0 * (xi1 * eta3_eta * psi2), 2.0 * (xi1 * eta3 * psi2_psi)],      // dN19
                    vec![2.0 * (xi1_xi * eta2 * psi2), 2.0 * (xi1 * eta2_eta * psi2), 2.0 * (xi1 * eta2 * psi2_psi)],       // dN20
                    vec![2.0 * (xi3_xi * eta2 * psi2), 2.0 * (xi3 * eta2_eta * psi2), 2.0 * (xi3 * eta2 * psi2_psi)],      // dN21
                    vec![2.0 * (xi2_xi * eta1 * psi2), 2.0 * (xi2 * eta1_eta * psi2), 2.0 * (xi2 * eta1 * psi2_psi)],       // dN22
                    vec![2.0 * (xi2_xi * eta3 * psi2), 2.0 * (xi2 * eta3_eta * psi2), 2.0 * (xi2 * eta3 * psi2_psi)],      // dN23
                    vec![2.0 * (xi2_xi * eta2 * psi1), 2.0 * (xi2 * eta2_eta * psi1), 2.0 * (xi2 * eta2 * psi1_psi)],       // dN24
                    vec![2.0 * (xi2_xi * eta2 * psi3), 2.0 * (xi2 * eta2_eta * psi3), 2.0 * (xi2 * eta2 * psi3_psi)],      // dN25
                    vec![2.0 * (xi2_xi * eta2 * psi2), 2.0 * (xi2 * eta2_eta * psi2), 2.0 * (xi2 * eta2 * psi2_psi)],       // dN26
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            _ => Err(ParseError::ElementError(format!("Shape functions for element type {:?} are not implemented", element_type))),
        }
    }

    /// Get the nodes associated with an element
    pub fn get_element_nodes(element: &Element, nodes: &[Node]) -> Result<Vec<Node>, ElementError> {
        let mut element_nodes = Vec::new();
        
        for &node_id in &element.nodes {
            if let Some(node) = nodes.iter().find(|n| n.id == node_id) {
                element_nodes.push(node.clone());
            } else {
                return Err(ElementError::InvalidElement(format!(
                    "Node {} not found for element {}", node_id, element.id
                )));
            }
        }
        
        if element_nodes.is_empty() {
            return Err(ElementError::InvalidElement(format!(
                "No valid nodes found for element {}", element.id
            )));
        }
        
        Ok(element_nodes)
    }

    /// Calculate Jacobian matrix: J[i][j] = sum(dN_k/dxi_j * x_k_i)
    pub fn calculate_jacobian(
        element_nodes: &[Node],
        shape_derivatives: &[Vec<f64>], // shape_derivatives[node][natural_dim]
    ) -> Result<Jacobian, ElementError> {
        if element_nodes.is_empty() {
            return Err(ElementError::GeometryError("No nodes provided".to_string()));
        }
        
        if element_nodes.len() != shape_derivatives.len() {
            return Err(ElementError::GeometryError(format!(
                "Node count ({}) doesn't match shape derivative count ({})",
                element_nodes.len(), shape_derivatives.len()
            )));
        }
        
        let spatial_dim = element_nodes[0].coordinates.len();
        let natural_dim = shape_derivatives[0].len();
        
        // Initialize Jacobian matrix
        let mut jacobian_matrix = vec![vec![0.0; natural_dim]; spatial_dim];
        
        // Calculate Jacobian: J[i][j] = sum_k(dN_k/dxi_j * x_k_i)
        for i in 0..spatial_dim {           // spatial coordinate index
            for j in 0..natural_dim {       // natural coordinate index  
                for (k, node) in element_nodes.iter().enumerate() { // node index
                    if j < shape_derivatives[k].len() && i < node.coordinates.len() {
                        jacobian_matrix[i][j] += shape_derivatives[k][j] * node.coordinates[i];
                    }
                }
            }
        }
        
        // Calculate determinant or sqrt(det(J^T * J)) for non-square matrices
        let determinant = Self::calculate_jacobian_determinant(&jacobian_matrix)?;
        
        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant,
        })
    }

    /// Calculate determinant of Jacobian matrix with proper handling of dimension mismatch
    /// Returns signed determinant for square matrices, signed magnitude for non-square matrices
    fn calculate_jacobian_determinant(matrix: &[Vec<f64>]) -> Result<f64, ParseError> {
        let spatial_dim = matrix.len();  // Number of rows = spatial dimensions of the mesh
        if spatial_dim == 0 {
            return Err(ParseError::FormatError("Cannot calculate determinant of empty matrix".to_string()));
        }
        
        // Verify matrix is rectangular (all rows have same length)
        for row in matrix {
            if row.len() != matrix[0].len() {
                return Err(ParseError::FormatError("Matrix must be rectangular for determinant calculation".to_string()));
            }
        }
        
        let natural_dim = matrix[0].len();  // Number of columns = natural/parametric dimensions of element
        
        match (spatial_dim, natural_dim) {
            // Case 1: Square matrices - element dimension matches mesh dimension
            // Use standard determinant calculation det(J)
            (1, 1) => {
                // 1D element in 1D mesh: det(J) = J[0,0] = dx/dxi
                Ok(matrix[0][0])
            },
            
            (2, 2) => {
                // 2D element in 2D mesh: det(J) = J[0,0]*J[1,1] - J[0,1]*J[1,0]
                // This gives the area scaling factor (can be negative if element is inverted)
                Ok(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0])
            },
            
            (3, 3) => {
                // 3D element in 3D mesh: standard 3x3 determinant
                // This gives the volume scaling factor (can be negative if element is inverted)
                Ok(matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))
            }
            
            // Case 2: Non-square matrices - element dimension != mesh dimension
            // Use signed sqrt(det(G)) where G = J^T * J for Jacobian magnitude
            (spatial, natural) if spatial != natural => {
                // Step 1: Calculate G = J^T * J (Gram matrix)
                // G[i,j] = sum_k(J[k,i] * J[k,j]) for all spatial dimensions k
                let mut g = vec![vec![0.0; natural]; natural];
                
                for i in 0..natural {           // Loop over natural coordinate i
                    for j in 0..natural {       // Loop over natural coordinate j
                        for k in 0..spatial {   // Sum over all spatial dimensions k
                            // G[i,j] += J[k,i] * J[k,j]
                            g[i][j] += matrix[k][i] * matrix[k][j];
                        }
                    }
                }
                
                // Step 2: Calculate determinant of G
                let det_g = match natural {
                    1 => {
                        // 1D element: G is 1x1, det(G) = G[0,0] = ||dx/dxi||²
                        g[0][0]
                    },
                    
                    2 => {
                        // 2D element: G is 2x2, standard 2x2 determinant
                        g[0][0] * g[1][1] - g[0][1] * g[1][0]
                    },
                    
                    3 => {
                        // 3D element: G is 3x3, standard 3x3 determinant  
                        g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
                        - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
                        + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0])
                    }
                    
                    _ => return Err(ParseError::FormatError(format!(
                        "Determinant calculation not implemented for {}x{} G matrices", natural, natural
                    ))),
                };
                
                // Step 3: Handle sign for element inversion detection
                // Since G = J^T * J is always positive semi-definite, det(G) ≥ 0
                // We need to determine the sign from the orientation of the element
                let sign = Self::calculate_element_orientation_sign(matrix, spatial, natural)?;
                
                // Step 4: Return signed sqrt(det(G)) as the Jacobian magnitude
                if det_g < 0.0 {
                    // This mathematically shouldn't happen for G = J^T*J, but handle gracefully
                    Ok(0.0)
                } else {
                    Ok(sign * det_g.sqrt())
                }
            }
            
            // Case 3: Unsupported square matrix sizes
            _ => {
                Err(ParseError::FormatError(format!(
                    "Determinant calculation not implemented for {}x{} matrices", spatial_dim, natural_dim
                )))
            }
        }
    }

    /// Calculate the orientation sign for non-square Jacobian matrices
    /// This determines if the element is inverted (negative orientation)
    fn calculate_element_orientation_sign(
        matrix: &[Vec<f64>], 
        spatial_dim: usize, 
        natural_dim: usize
    ) -> Result<f64, ParseError> {
        match (spatial_dim, natural_dim) {
            // Case: 1D element in 2D or 3D mesh
            (2, 1) | (3, 1) => {
                // For 1D elements (lines), orientation is determined by the direction vector
                // J = [dx/dxi, dy/dxi] or [dx/dxi, dy/dxi, dz/dxi]
                // Sign is positive if element follows expected parametric direction
                // For simplicity, we can use the sign of the first non-zero component
                for i in 0..spatial_dim {
                    if matrix[i][0].abs() > 1e-12 {  // Found first significant component
                        return Ok(if matrix[i][0] > 0.0 { 1.0 } else { -1.0 });
                    }
                }
                // All components are zero - degenerate element
                Ok(1.0)  // Default to positive
            },
                
            // Case: 2D element in 3D mesh  
            (3, 2) => {
                // For 2D elements in 3D (surfaces), orientation is determined by cross product
                // J = [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
                // Cross product: (dx/dxi, dy/dxi, dz/dxi) × (dx/deta, dy/deta, dz/deta)
                let u = [matrix[0][0], matrix[1][0], matrix[2][0]];  // First column: ∂r/∂xi
                let v = [matrix[0][1], matrix[1][1], matrix[2][1]];  // Second column: ∂r/∂eta
                    
                // Cross product w = u × v
                let cross_product = [
                    u[1] * v[2] - u[2] * v[1],  // i component
                    u[2] * v[0] - u[0] * v[2],  // j component  
                    u[0] * v[1] - u[1] * v[0],  // k component
                ];
                    
                // The sign is determined by the z-component of the normal vector
                // (assuming standard right-hand rule orientation)
                // For more robust orientation, we could use the largest component
                let max_component = cross_product.iter()
                    .max_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(&0.0);
                    
                Ok(if *max_component >= 0.0 { 1.0 } else { -1.0 })
            },
                
            // Other cases not yet implemented
            _ => {
                // Default to positive orientation for unhandled cases
                Ok(1.0)
            }
        }
    }

    /// Calculate statistics for all element qualities
    fn calculate_statistics(qualities: &[ElementQuality]) -> QualityStatistics {
        let jacobians: Vec<f64> = qualities.iter().map(|q| q.det_jacobian).collect();
        
        let min_jacobian = jacobians.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_jacobian = jacobians.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let avg_jacobian = jacobians.iter().sum::<f64>() / jacobians.len() as f64;
        let negative_jacobian_count = jacobians.iter().filter(|&&j| j < 0.0).count();

        QualityStatistics {
            min_jacobian,
            max_jacobian,
            avg_jacobian,
            negative_jacobian_count,
        }
    }
}
            
               
 



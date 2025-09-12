//use std::fs::File;                          
use std::io::{self};    // For input/output operations 
use std::f64;            // Import f64 constants like INFINITY

//use crate::database::*;                     // Import mesh data structures and error types from database module

use crate::lib::*;                  // Import mesh data structures and error types from lib module
use crate::error::*;                     // Import mesh data structures and error types from error module


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

#[derive(Debug, Clone)]
pub struct ElementQuality {
    pub element_id: usize,      // ID of the element being analyzed
    pub det_jacobian: f64,      // determinant of the Jacobian matrix

    // more quality metrics can be added here
}

// Structure to hold the complete mesh quality analysis results
#[derive(Debug, Clone)]
pub struct MeshQualityReport {
    pub total_elements: usize,                    // Total number of elements that were successfully analyzed
    pub element_qualities: Vec<ElementQuality>,   // Quality metrics for each individual element
    pub statistics: QualityStatistics,           // Overall statistical summary of mesh quality
}
 
// Structure to hold statistical summary of mesh quality metrics
#[derive(Debug, Clone)]
pub struct QualityStatistics {
    pub min_jacobian: f64,              // Minimum Jacobian determinant in the mesh
    pub max_jacobian: f64,              // Maximum Jacobian determinant in the mesh
    pub avg_jacobian: f64,              // Average Jacobian determinant across all elements
    pub negative_jacobian_count: usize, // Number of elements with negative Jacobian (invalid elements)
}

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
    fn get_shape_functions(element_type: &ElementType) -> Result<ShapeFunction, ParseError> {  //ShapeFunctionError
        match element_type {
            ElementType::Line => {
                let reference_coords = vec![
                    vec![0.0],    // Node 0 at xi = 0
                    vec![1.0],     // Node 1 at xi = 1
                ];
                let num_nodes = 2;

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                
                let values = vec![
                    1.0 - xi[0],   // N0: node at xi = 0
                    xi[1],         // N1: node at xi = 1
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                
                let values = vec![
                    2.0 * (xi[0] - 0.5) * (xi[0] - 1.0),      // N0: node at xi = 0
                    2.0 * xi[1] * (xi[1] - 0.5),              // N1: node at xi = 1
                    4.0 * xi[2] * (1.0 - xi[2]),              // N2: node at xi = 0.5
                ];
                let derivatives = vec![
                    vec![4.0 * xi[0] - 3.0],               // dN0/dxi
                    vec![4.0 * xi[1] - 1.0],               // dN1/dxi
                    vec![4.0 - xi[2] * 8.0],               // dN2/dxi
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                
                let values = vec![
                    1.0 - xi[0] - eta[0],               // N0: node at (0,0)
                    xi[1],                              // N1: node at (1,0)
                    eta[2],                             // N2: node at (0,1)
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();

                let lambda: Vec<f64> = xi.iter()
                    .zip(eta.iter())
                    .map(|(&x, &e)| 1.0 - x - e)
                    .collect();
                
                let values = vec![
                    lambda[0] * (2.0 * lambda[0] - 1.0),      // N0: corner at 
                    xi[1] * (2.0 * xi[1] - 1.0),              // N1: corner at 
                    eta[2] * (2.0 * eta[2] - 1.0),            // N2: corner at 
                    4.0 * xi[3] * lambda[3],                  // N3: mid-edge at 
                    4.0 * xi[4] * eta[4],                     // N4: mid-edge at 
                    4.0 * eta[5] * lambda[5],                 // N5: mid-edge at 
                ];
                let derivatives = vec![
                    vec![4.0 * xi[0] + 4.0 * eta[0] - 3.0, 4.0 * xi[0] + 4.0 * eta[0] - 3.0],       // dN0
                    vec![4.0 * xi[1] - 1.0, 0.0],                                          // dN1
                    vec![0.0, 4.0 * eta[2] - 1.0],                                         // dN2
                    vec![4.0 - 8.0 * xi[3] - 4.0 * eta[3], -4.0 * xi[3]],                        // dN3
                    vec![4.0 * eta[4],  4.0 * xi[4]],                                         // dN4
                    vec![-4.0 * xi[5], 4.0 - 8.0 * eta[5] - 4.0 * xi[5]],                        // dN5
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();

                let xim: Vec<f64> = xi.iter()
                    .map(|&x| 1.0 - x)
                    .collect();
                let etam: Vec<f64> = eta.iter()
                    .map(|&e| 1.0 - e)
                    .collect();
                
                let values = vec![
                    xim[0] * etam[0],         // N0: 
                    xi[1] * etam[1],          // N1: 
                    xi[2] * eta[2],           // N2: 
                    xim[3] * eta[3],          // N3: 
                ];
                let derivatives = vec![
                    vec![-etam[0], -xim[0]],      // dN0
                    vec![etam[1], -xi[1]],        // dN1
                    vec![eta[2], xi[2]],          // dN2
                    vec![-eta[3], xim[3]],        // dN3
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();

                let values = vec![
                    (1.0 - xi[0]) * (1.0 - eta[0]) - 0.5 * (4.0 * xi[0] * (1.0 - xi[0]) * (1.0 - eta[0]) + 4.0 * (1.0 - xi[0]) * (1.0 - eta[0]) * eta[0]),          // N0: corner
                    xi[1] * (1.0 - eta[1]) - 0.5 * (4.0 * xi[1] * (1.0 - xi[1]) * (1.0 - eta[1]) + 4.0 * xi[1] * (1.0 - eta[1]) * eta[1]),                          // N1: corner
                    xi[2] * eta[2] - 0.5 * (4.0 * xi[2] * (1.0 - eta[2]) * eta[2] + 4.0 * xi[2] * (1.0 - xi[2]) * eta[2]),                                          // N2: corner
                    (1.0 - xi[3]) * eta[3] - 0.5 * (4.0 * xi[3] * (1.0 - xi[3]) * eta[3] + 4.0 * (1.0 - xi[3]) * (1.0 - eta[3]) * eta[3]),                          // N3: corner
                    4.0 * xi[4] * (1.0 - xi[4]) * (1.0 - eta[4]),                                                                                  // N4: mid-edge bottom
                    4.0 * xi[5] * (1.0 - eta[5]) * eta[5],                                                                                         // N5: mid-edge right
                    4.0 * xi[6] * (1.0 - xi[6]) * eta[6],                                                                                          // N6: mid-edge top
                    4.0 * (1.0 - xi[7]) * (1.0 - eta[7]) * eta[7],                                                                                 // N7: mid-edge left
                ];
                let derivatives = vec![
                    vec![-(1.0 - eta[0]) - 0.5 * (4.0 * (1.0 - eta[0]) * (1.0 - 2.0 * xi[0]) - 4.0 * (1.0 - eta[0]) * eta[0]), -(1.0 - xi[0]) - 0.5 * (-4.0 * xi[0] * (1.0 - xi[0]) +  4.0 * (1.0 - xi[0]) * (1.0 - 2.0 * eta[0]))],      // dN0
                    vec![(1.0 - eta[1]) - 0.5 * (4.0 * (1.0 - eta[1]) * (1.0 - 2.0 * xi[1]) + 4.0 * (1.0 - eta[1]) * eta[1]), -xi[1] - 0.5 * (-4.0 * xi[1] * (1.0 - xi[1]) +  4.0 * (1.0 - xi[1]) * (1.0 - 2.0 * eta[1]))],               // dN1
                    vec![eta[2] - 0.5 * (4.0 * (1.0 - eta[2]) * eta[2] +  4.0 * eta[2] * (1.0 - 2.0 * xi[2])), xi[2] - 0.5 * (4.0 * xi[2] * (1.0 - 2.0 * eta[2]) + 4.0 * xi[2] * (1.0 - xi[2]))],                                         // dN2
                    vec![-eta[3] - 0.5 * (4.0 * (1.0 - eta[3]) * (1.0 - 2.0 * xi[3]) - 4.0 * (1.0 - eta[3]) * eta[3]), (1.0 - xi[3]) - 0.5 * (4.0 * xi[3] * (1.0 - xi[3]) +  4.0 * (1.0 - xi[3]) * (1.0 - 2.0 * eta[3]))],                // dN3
                    vec![4.0 * (1.0 - eta[4]) * (1.0 - 2.0 * xi[4]), -4.0 * (1.0 - xi[4]) * (1.0 - xi[4])],                                                                                                             // dN4
                    vec![4.0 * (1.0 - eta[5]) * eta[5], 4.0 * xi[5] * (1.0 - 2.0 * eta[5])],                                                                                                                            // dN5
                    vec![4.0 * (1.0 - eta[6]) * (1.0 - 2.0 * xi[6]), 4.0 * xi[6] * (1.0 - xi[6])],                                                                                                                      // dN6
                    vec![-4.0 * (1.0 - eta[7]) * eta[7], 4.0 * (1.0 - xi[7]) * (1.0 - 2.0 * eta[7])],                                                                                                                   // dN7
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();

                let values = vec![
                    4.0 * (1.0 - xi[0]) * (xi[0] - 0.5) * (1.0 - eta[0]) * (eta[0] - 0.5),      // N0:  
                    -4.0 * xi[1] * (xi[1] - 0.5) * (1.0 - eta[1]) * (eta[1] - 0.5),             // N1:  
                    4.0 * xi[2] * (xi[2] - 0.5) * eta[2] * (eta[2] - 0.5),                      // N2:  
                    -4.0 * (1.0 - xi[3]) * (xi[3] - 0.5) * eta[3] * (eta[3] - 0.5),             // N3:  
                    8.0 * xi[4] * (1.0 - xi[4]) * (1.0 - eta[4]) * (0.5 - eta[4]),              // N4:  
                    -8.0 * xi[5] * (0.5 - xi[5]) * (1.0 - eta[5]) * eta[5],                     // N5:  
                    -8.0 * xi[6] * (1.0 - xi[6]) * eta[6] * (0.5 - eta[6]),                     // N6:  
                    8.0 * (1.0 - xi[7]) * (0.5 - xi[7]) * (1.0 - eta[7]) * eta[7],              // N7:  
                    16.0 * xi[8] * (1.0 - xi[8]) * (1.0 - eta[8]) * eta[8],                     // N8: 
                ];
                let derivatives = vec![
                    vec![4.0 * (1.5 - 2.0 * xi[0]) * (1.0 - eta[0]) * (eta[0] - 0.5), 4.0 * (1.0 - xi[0]) * (xi[0] - 0.5) * (1.5 - 2.0 * eta[0])],    // dN0
                    vec![-4.0 * (2.0 * xi[1] - 0.5) * (1.0 - eta[1]) * (eta[1] - 0.5), -4.0 * xi[1] * (xi[1] - 0.5) * (1.5 - 2.0 * eta[1])],          // dN1
                    vec![4.0 * (2.0 * xi[2] - 0.5) * eta[2] * (eta[2] - 0.5), 4.0 * xi[2] * (xi[2] - 0.5) * (2.0 * eta[2] - 0.5)],                    // dN2
                    vec![-4.0 * (1.5 - 2.0 * xi[3]) * eta[3] * (eta[3] - 0.5), -4.0 * (1.0 - xi[3]) * (xi[3] - 0.5) * (2.0 * eta[3] - 0.5)],          // dN3
                    vec![8.0 * (1.0 - 2.0 * xi[4]) * (1.0 - eta[4]) * (0.5 - eta[4]), 8.0 * xi[4] * (1.0 - xi[4]) * (2.0 * eta[4] - 1.5)],            // dN4
                    vec![-8.0 * (0.5 - 2.0 * xi[5]) * (1.0 - eta[5]) * eta[5], -8.0 * xi[5] * (0.5 - xi[5]) * (1.0 - 2.0 * eta[5])],                  // dN5
                    vec![-8.0 * (1.0 - 2.0 * xi[6]) * eta[6] * (0.5 - eta[6]), -8.0 * xi[6] * (1.0 - xi[6]) * (0.5 - 2.0 * eta[6])],                  // dN6
                    vec![8.0 * (2.0 * xi[7] - 1.5) * (1.0 - eta[7]) * eta[7], 8.0 * (1.0 - xi[7]) * (0.5 - xi[7]) * (1.0 - 2.0 * eta[7])],            // dN7
                    vec![16.0 * (1.0 - 2.0 * xi[8]) * (1.0 - eta[8]) * eta[8], 16.0 * xi[8] * (1.0 - xi[8]) * (1.0 - 2.0 * eta[8])],                  // dN8
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
                    1.0 - xi[0] - eta[0] - psi[0],       // N0: node at 
                    xi[1],                         // N1: node at 
                    eta[2],                        // N2: node at 
                    psi[3],                        // N3: node at
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let lambda: Vec<f64> = xi.iter()
                    .zip(eta.iter())
                    .zip(psi.iter())
                    .map(|((&x, &e), &p)| 1.0 - x - e - p)
                    .collect();

                let values = vec![
                    lambda[0] * (2.0 * lambda[0] - 1.0),          // N0: corner
                    xi[1] * (2.0 * xi[1] - 1.0),                  // N1: corner  
                    eta[2] * (2.0 * eta[2] - 1.0),                // N2: corner 
                    psi[3] * (2.0 * psi[3] - 1.0),                // N3: corner 
                    4.0 * lambda[4] * xi[4],                      // N4: mid-edge 
                    4.0 * xi[5] * eta[5],                         // N5: mid-edge 
                    4.0 * eta[6] * lambda[6],                     // N6: mid-edge 
                    4.0 * lambda[7] * psi[7],                     // N7: mid-edge 
                    4.0 * xi[8] * psi[8],                         // N8: mid-edge 
                    4.0 * eta[9] * psi[9],                        // N9: mid-edge 
                ];
                let derivatives = vec![
                    vec![4.0 * (xi[0] + eta[0] + psi[0]) - 3.0, 4.0 * (xi[0] + eta[0] + psi[0]) - 3.0, 4.0 * (xi[0] + eta[0] + psi[0]) - 3.0], // dN0
                    vec![4.0 * xi[1] - 1.0, 0.0, 0.0],                                                                 // dN1
                    vec![0.0, 4.0 * eta[2] - 1.0, 0.0],                                                                // dN2
                    vec![0.0, 0.0, 4.0 * psi[3] - 1.0],                                                                // dN3
                    vec![4.0 - 8.0 * xi[4] - 4.0 * eta[4] - 4.0 * psi[4], -4.0 * xi[4], -4.0 * xi[4]],                             // dN4
                    vec![4.0 * eta[5], 4.0 * xi[5], 0.0],                                                                 // dN5
                    vec![-4.0 * eta[6], 4.0 - 4.0 * xi[6] - 8.0 * eta[6] - 4.0 * psi[6], -4.0 * eta[6]],                           // dN6
                    vec![-4.0 * psi[7], -4.0 * psi[7], 4.0 - 4.0 * xi[7] - 4.0 * eta[7] - 8.0 * psi[7]],                           // dN7
                    vec![4.0 * psi[8], 0.0, 4.0 * xi[8]],                                                                 // dN8
                    vec![0.0, 4.0 * psi[9], 4.0 * eta[9]],                                                                // dN9
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let xim: Vec<f64> = xi.iter()
                    .map(|&x| 1.0 - x)
                    .collect();
                let etam: Vec<f64> = eta.iter()
                    .map(|&e| 1.0 - e)
                    .collect();
                let psim: Vec<f64> = psi.iter()
                    .map(|&p| 1.0 - p)
                    .collect();

                let values = vec![
                    xim[0] * etam[0] * psim[0],      // N0: base 
                    xi[1] * etam[1] * psim[1],       // N1: base 
                    xi[2] * eta[2] * psim[2],        // N2: base 
                    xim[3] * eta[3] * psim[3],       // N3: base 
                    psi[4],                    // N4: apex 
                ];
                let derivatives = vec![
                    vec![-etam[0] * psim[0], -xim[0] * psim[0], -xim[0] * etam[0]],   // dN0
                    vec![etam[1] * psim[1], -xi[1] * psim[1], -xi[1] * etam[1]],      // dN1
                    vec![eta[2] * psim[2], xi[2] * psim[2], -xi[2] * eta[2]],         // dN2
                    vec![-eta[3] * psim[3], xim[3] * psim[3], -xim[3] * eta[3]],      // dN3
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
                    (1.0 - xi[0] - eta[0]) * (1.0 - psi[0]),     // N0: bottom
                    xi[1] * (1.0 - psi[1]),                   // N1: bottom
                    eta[2] * (1.0 - psi[2]),                  // N2: bottom
                    (1.0 - xi[3] - eta[3]) * psi[3],             // N3: top
                    xi[4] * psi[4],                           // N4: top
                    eta[5] * psi[5],                          // N5: top
                ];
                let derivatives = vec![
                    vec![-1.0 + psi[0], -1.0 + psi[0], -1.0 + xi[0] + eta[0]],  // dN0
                    vec![1.0 - psi[1], 0.0, -xi[1]],                      // dN1
                    vec![0.0, 1.0 - psi[2], -eta[2]],                     // dN2
                    vec![-psi[3], -psi[3], 1.0 - xi[3] - eta[3]],               // dN3
                    vec![psi[4], 0.0, xi[4]],                             // dN4
                    vec![0.0, psi[5], eta[5]],                            // dN5
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
                    // Corner nodes
                    2.0 * (1.0 - xi[0] - eta[0]) * (1.0 - psi[0]) * (0.5 - xi[0] - eta[0] - psi[0]),    // N0: bottom corner
                    2.0 * xi[1] * (1.0 - psi[1]) * (xi[1] - psi[1] - 0.5),                      // N1: bottom corner
                    2.0 * eta[2] * (1.0 - psi[2]) * (eta[2] - psi[2] - 0.5),                    // N2: bottom corner
                    2.0 * (1.0 - xi[3] - eta[3]) * psi[3] * (psi[3] - xi[3] - eta[3] - 0.5),          // N3: top corner
                    2.0 * xi[4] * psi[4] * (xi[4] + psi[4] - 1.5),                            // N4: top corner
                    2.0 * eta[5] * psi[5] * (eta[5] + psi[5] - 1.5),                          // N5: top corner
                    // Mid-edge nodes on faces
                    4.0 * xi[6] * (1.0 - xi[6] - eta[6]) * (1.0 - psi[6]),                        // N6: bottom face
                    4.0 * xi[7] * eta[7] * (1.0 - psi[7]),                                   // N7: bottom face
                    4.0 * (1.0 - xi[8] - eta[8]) * eta[8] * (1.0 - psi[8]),                       // N8: bottom face
                    4.0 * xi[9] * (1.0 - xi[9] - eta[9]) * psi[9],                              // N9: top face
                    4.0 * xi[10] * eta[10] * psi[10],                                         // N10: top face
                    4.0 * (1.0 - xi[11] - eta[11]) * eta[11] * psi[11],                             // N11: top face
                    // Vertical mid-edge nodes
                    4.0 * psi[12] * (1.0 - xi[12] - eta[12]) * (1.0 - psi[12]),                       // N12: vertical
                    4.0 * psi[13] * xi[13] * (1.0 - psi[13]),                                   // N13: vertical
                    4.0 * psi[14] * eta[14] * (1.0 - psi[14]),                                  // N14: vertical
                ];
                let derivatives = vec![
                    // Corner nodes
                    vec![2.0 * (1.0 - psi[0]) * (-1.5 + 2.0 * xi[0] + 2.0 * eta[0] + psi[0]), 2.0 * (1.0 - psi[0]) * (-1.5 + 2.0 * xi[0] + 2.0 * eta[0] + psi[0]), 2.0 * (1.0 - xi[0] - eta[0]) * (-1.5 + xi[0] + eta[0] + 2.0 * psi[0])], // dN0
                    vec![2.0 * (1.0 - psi[1]) * (-0.5 + 2.0 * xi[1] - psi[1]), 0.0, 2.0 * xi[1] * (-0.5 - xi[1] + 2.0 * psi[1])], // dN1
                    vec![0.0, 2.0 * (1.0 - psi[2]) * (-0.5 + 2.0 * eta[2] - psi[2]), 2.0 * eta[2] * (-0.5 - eta[2] + 2.0 * psi[2])], // dN2
                    vec![2.0 * psi[3] * (-0.5 + 2.0 * xi[3] + 2.0 * eta[3] - psi[3]), 2.0 * psi[3] * (-0.5 + 2.0 * xi[3] + 2.0 * eta[3] - psi[3]), 2.0 * (1.0 - xi[3] - eta[3]) * (-0.5 - xi[3] - eta[3] + 2.0 * psi[3])], // dN3
                    vec![2.0 * psi[4] * (-1.5 + 2.0 * xi[4] + psi[4]), 0.0, 2.0 * xi[4] * (-1.5 + xi[4] + 2.0 * psi[4])], // dN4
                    vec![0.0, 2.0 * psi[5] * (-1.5 + 2.0 * eta[5] + psi[5]), 2.0 * eta[5] * (-1.5 + eta[5] + 2.0 * psi[5])], // dN5
                    // Mid-edge derivatives
                    vec![4.0 * (1.0 - psi[6]) * (1.0 - 2.0 * xi[6] - eta[6]), -4.0 * (1.0 - psi[6]) * xi[6], -4.0 * xi[6] * (1.0 - xi[6] - eta[6])], // dN6
                    vec![4.0 * (1.0 - psi[7]) * eta[7], 4.0 * (1.0 - psi[7]) * xi[7], -4.0 * xi[7] * eta[7]], // dN7
                    vec![-4.0 * (1.0 - psi[8]) * eta[8], 4.0 * (1.0 - psi[8]) * (1.0 - xi[8] - 2.0 * eta[8]), -4.0 * eta[8] * (1.0 - xi[8] - eta[8])], // dN8
                    vec![4.0 * psi[9] * (1.0 - 2.0 * xi[9] - eta[9]), -4.0 * xi[9] * psi[9], 4.0 * xi[9] * (1.0 - xi[9] - eta[9])], // dN9
                    vec![4.0 * eta[10] * psi[10], 4.0 * xi[10] * psi[10], 4.0 * xi[10] * eta[10]], // dN10
                    vec![-4.0 * eta[11] * psi[11], 4.0 * psi[11] * (1.0 - xi[11] - 2.0 * eta[11]), 4.0 * eta[11] * (1.0 - xi[11] - eta[11])], // dN11
                    // Vertical mid-edge derivatives
                    vec![-4.0 * psi[12] * (1.0 - psi[12]), -4.0 * psi[12] * (1.0 - psi[12]), 4.0 * (1.0 - 2.0 * psi[12]) * (1.0 - xi[12] - eta[12])], // dN12
                    vec![4.0 * psi[13] * (1.0 - psi[13]), 0.0, 4.0 * (1.0 - 2.0 * psi[13]) * xi[13]], // dN13
                    vec![0.0, 4.0 * psi[14] * (1.0 - psi[14]), 4.0 * (1.0 - 2.0 * psi[14]) * eta[14]], // dN14
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let xim: Vec<f64> = xi.iter()
                    .map(|&x| 2.0 * (x - 0.5))
                    .collect();
                let etam: Vec<f64> = eta.iter()
                    .map(|&e| 2.0 * (e - 0.5))
                    .collect();
                let psim: Vec<f64> = psi.iter()
                    .map(|&p| 2.0 * (p - 0.5))
                    .collect();

                let values = vec![
                    // Corner nodes
                    -0.25 * (xim[0] + etam[0]) * (xim[0] + etam[0] + 1.0) * psim[0] * (1.0 - psim[0]),    // N0: bottom corner
                    -0.25 * xim[1] * (xim[1] + 1.0) * psim[1] * (1.0 - psim[1]),                    // N1: bottom corner
                    -0.25 * etam[2] * (1.0 + etam[2]) * psim[2] * (1.0 - psim[2]),                  // N2: bottom corner
                    0.25 * (xim[3] + etam[3]) * (xim[3] + etam[3] + 1.0) * psim[3] * (1.0 + psim[3]),     // N3: top corner
                    0.25 * xim[4] * (xim[4] + 1.0) * psim[4] * (1.0 + psim[4]),                     // N4: top corner
                    0.25 * etam[5] * (1.0 + etam[5]) * psim[5] * (1.0 + psim[5]),                   // N5: top corner
                    // Mid-edge nodes on faces
                    (xim[6] + 1.0) * (xim[6] + etam[6]) * 0.5 * psim[6] * (1.0 - psim[6]),                 // N6: bottom face
                    -(xim[7] + 1.0) * (etam[7] + 1.0) * 0.5 * psim[7] * (1.0 - psim[7]),                  // N7: bottom face
                    (etam[8] + 1.0) * (xim[8] + etam[8]) * 0.5 * psim[8] * (1.0 - psim[8]),                // N8: bottom face
                    -(xim[9] + 1.0) * (xim[9] + etam[9]) * 0.5 * psim[9] * (1.0 + psim[9]),                // N9: top face
                    (xim[10] + 1.0) * (etam[10] + 1.0) * 0.5 * psim[10] * (1.0 + psim[10]),                   // N10: top face
                    -(etam[11] + 1.0) * (xim[11] + etam[11]) * 0.5 * psim[11] * (1.0 + psim[11]),               // N11: top face
                    // Vertical mid-edge nodes
                    0.5 * (xim[12] + etam[12]) * (xim[12] + etam[12] + 1.0) * (1.0 + psim[12]) * (1.0 - psim[12]),    // N12: vertical
                    0.5 * xim[13] * (xim[13] + 1.0) * (1.0 + psim[13]) * (1.0 - psim[13]),                    // N13: vertical
                    0.5 * etam[14] * (1.0 + etam[14]) * (1.0 + psim[14]) * (1.0 - psim[14]),                  // N14: vertical
                    // Mid-face nodes
                    -(xim[15] + 1.0)*(xim[15] + etam[15]) * (1.0 + psim[15]) * (1.0 - psim[15]),                  // N15: mid-face
                    (xim[16] + 1.0)*(etam[16] + 1.0) * (1.0 + psim[16]) * (1.0 - psim[16]),                   // N16: mid-face
                    -(etam[17] + 1.0)*(xim[17] + etam[17]) * (1.0 + psim[17]) * (1.0 - psim[17]),                 // N17: mid-face
                ];
                let derivatives = vec![
                    // we compute derivatives in [-1; 1] but we need them in [ 0; 1]
                        // for(int i = 0; i < 54; i++)
                            // derivs[i] *= 2;

                    // Corner derivatives
                    vec![2.0 * (-0.25 * (2.0 * xim[0] + 2.0 * etam[0] + 1.0) * psim[0] * (1.0 - psim[0])), 
                        2.0 * (-0.25 * (2.0 * etam[0] + 2.0 * xim[0] + 1.0) * psim[0] * (1.0 - psim[0])), 
                        2.0 * (-0.25 * (xim[0] + etam[0]) * (xim[0] + etam[0] + 1.0) * (1.0 - 2.0 * psim[0]))],   // dN0
                    vec![2.0 * (-0.25 * (2.0 * xim[1] + 1.0) * psim[1] * (1.0 - psim[1])), 
                        2.0 * 0.0, 
                        2.0 * (-0.25 * xim[1] * (xim[1] + 1.0) * (1.0 - 2.0 * psim[1]))],                   // dN1
                    vec![2.0 * 0.0, 
                        2.0 * (-0.25 * (2.0 * etam[2] + 1.0) * psim[2] * (1.0 - psim[2])), 
                        2.0 * (-0.25 * etam[2] * (1.0 + etam[2]) * (1.0 - 2.0 * psim[2]))],                 // dN2
                    vec![2.0 * (0.25 * (2.0 * xim[3] + 2.0 * etam[3] + 1.0) * psim[3] * (1.0 + psim[3])), 
                        2.0 * (0.25 * (2.0 * etam[3] + 2.0 * xim[3] + 1.0) * psim[3] * (1.0 + psim[3])), 
                        2.0 * (0.25 * (xim[3] + etam[3]) * (xim[3] + etam[3] + 1.0) * (1.0 + 2.0 * psim[3]))],    // dN3
                    vec![2.0 * (0.25 * (2.0 * xim[4] + 1.0) * psim[4] * (1.0 + psim[4])), 
                        2.0 * 0.0, 
                        2.0 * (0.25 * xim[4] * (xim[4] + 1.0) * (1.0 + 2.0 * psim[4]))],                    // dN4
                    vec![2.0 * 0.0, 
                        2.0 * (0.25 * (2.0 * etam[5] + 1.0) * psim[5] * (1.0 + psim[5])), 
                        2.0 * (0.25 * etam[5] * (1.0 + etam[5]) * (1.0 + 2.0 * psim[5]))],                  // dN5
                    // Face mid-edge derivatives
                    vec![2.0 * ((2.0 * xim[6] + etam[6] + 1.0) * 0.5 * psim[6] * (1.0 - psim[6])), 
                        2.0 * ((xim[6] + 1.0) * 0.5 * psim[6] * (1.0 - psim[6])), 
                        2.0 * ((xim[6] + 1.0) * (xim[6] + etam[6]) *  0.5 * (1.0 - 2.0 * psim[6]))],           // dN6
                    vec![2.0 * (-(etam[7] + 1.0) * 0.5 * psim[7] * (1.0 - psim[7])), 
                        2.0 * (-(xim[7] + 1.0) * 0.5 * psim[7] * (1.0 - psim[7])), 
                        2.0 * (-(xim[7] + 1.0) * (etam[7] + 1.0) *  0.5 * (1.0 - 2.0 * psim[7]))],          // dN7
                    vec![2.0 * ((etam[8] + 1.0) * 0.5 * psim[8] * (1.0 - psim[8])), 
                        2.0 * ((2.0 * etam[8] + xim[8] + 1.0) * 0.5 * psim[8] * (1.0 - psim[8])), 
                        2.0 * ((etam[8] + 1.0) * (xim[8] + etam[8]) *  0.5 * (1.0 - 2.0 * psim[8]))],          // dN8
                    vec![2.0 * (-(2.0 * xim[9] + etam[9] + 1.0) * 0.5 * psim[9] * (1.0 + psim[9])), 
                        2.0 * (-(xim[9] + 1.0) * 0.5 * psim[9] * (1.0 + psim[9])), 
                        2.0 * (-(xim[9] + 1.0) * (xim[9] + etam[9]) *  0.5 * (1.0 + 2.0 * psim[9]))],          // dN9
                    vec![2.0 * ((etam[10] + 1.0) * 0.5 * psim[10] * (1.0 + psim[10])), 
                        2.0 * ((xim[10] + 1.0) * 0.5 * psim[10] * (1.0 + psim[10])), 
                        2.0 * ((xim[10] + 1.0) * (etam[10] + 1.0) *  0.5 * (1.0 + 2.0 * psim[10]))],           // dN10
                    vec![2.0 * (-(etam[11] + 1.0) * 0.5 * psim[11] * (1.0 + psim[11])), 
                        2.0 * (-(2.0 * etam[11] + xim[11] + 1.0) * 0.5 * psim[11] * (1.0 + psim[11])), 
                        2.0 * (-(etam[11] + 1.0) * (xim[11] + etam[11]) *  0.5 * (1.0 + 2.0 * psim[11]))],         // dN11
                    // Vertical mid-edge derivatives
                    vec![2.0 * (0.5 * (2.0 * xim[12] + 2.0 * etam[12] + 1.0) * (1.0 + psim[12]) * (1.0 - psim[12])),  
                        2.0 * (0.5 * (2.0 * etam[12] + 2.0 * xim[12] + 1.0) * (1.0 + psim[12]) * (1.0 - psim[12])), 
                        2.0 * (0.5 * (xim[12] + etam[12]) * (xim[12] + etam[12] + 1.0) * (-2.0 * psim[12]))],          // dN12
                    vec![2.0 * (0.5 * (2.0 * xim[13] + 1.0) * (1.0 + psim[13]) * (1.0 - psim[13])), 
                        2.0 * 0.0, 
                        2.0 * (0.5 * xim[13] * (xim[13] + 1.0) * (-2.0 * psim[13]))],                          // dN13
                    vec![2.0 * 0.0, 
                        2.0 * (0.5 * (2.0 * etam[14] + 1.0) * (1.0 + psim[14]) * (1.0 - psim[14])), 
                        2.0 * (0.5 * etam[14] * (1.0 + etam[14]) * (-2.0 * psim[14]))],                        // dN14
                    // Mid-face derivatives
                    vec![2.0 * (-(2.0 * xim[15] + etam[15] + 1.0) * (1.0 + psim[15]) * (1.0 - psim[15])), 
                        2.0 * (-(xim[15] + 1.0) * (1.0 + psim[15]) * (1.0 - psim[15])), 
                        2.0 * (-(xim[15] + 1.0) * (xim[15] + etam[15]) * (-2.0 * psim[15]))],                      // dN15
                    vec![2.0 * ((etam[16] + 1.0) * (1.0 + psim[16]) * (1.0 - psim[16])), 
                        2.0 * ((xim[16] + 1.0) * (1.0 + psim[16]) * (1.0 - psim[16])), 
                        2.0 * ((xim[16] + 1.0) * (etam[16] + 1.0) * (-2.0 * psim[16]))],                       // dN16
                    vec![2.0 * (-(etam[17] + 1.0) * (1.0 + psim[17]) * (1.0 - psim[17])), 
                        2.0 * (-(2.0 * etam[17] + xim[17] + 1.0) * (1.0 + psim[17]) * (1.0 - psim[17])), 
                        2.0 * (-(etam[17] + 1.0) * (xim[17] + etam[17]) * (-2.0 * psim[17]))],                     // dN17
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let xim: Vec<f64> = xi.iter()
                    .map(|&x| 1.0 - x)
                    .collect();
                let etam: Vec<f64> = eta.iter()
                    .map(|&e| 1.0 - e)
                    .collect();
                let psim: Vec<f64> = psi.iter()
                    .map(|&p| 1.0 - p)
                    .collect();

                let values = vec![
                    xim[0] * etam[0] * psim[0],      // N0: 
                    xi[1] * etam[1] * psim[1],       // N1: 
                    xi[2] * eta[2] * psim[2],        // N2: 
                    xim[3] * eta[3] * psim[3],       // N3: 
                    xim[4] * etam[4] * psi[4],       // N4: 
                    xi[5] * etam[5] * psi[5],        // N5: 
                    xi[6] * eta[6] * psi[6],         // N6: 
                    xim[7] * eta[7] * psi[7],        // N7: 
                ];
                let derivatives = vec![
                    vec![-etam[0] * psim[0], -xim[0] * psim[0], -xim[0] * etam[0]], // dN0
                    vec![etam[1] * psim[1], -xi[1] * psim[1], -xi[1] * etam[1]], // dN1
                    vec![eta[2] * psim[2], xi[2] * psim[2], -xi[2] * eta[2]], // dN2
                    vec![-eta[3] * psim[3], xim[3] * psim[3], -xim[3] * eta[3]], // dN3
                    vec![-etam[4] * psi[4], -xim[4] * psi[4], xim[4] * etam[4]], // dN4
                    vec![etam[5] * psi[5], -xi[5] * psi[5], xi[5] * etam[5]], // dN5
                    vec![eta[6] * psi[6], xi[6] * psi[6], xi[6] * eta[6]], // dN6
                    vec![-eta[7] * psi[7], xim[7] * psi[7], xim[7] * eta[7]], // dN7
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let xik: Vec<f64> = xi.iter()
                    .map(|&x| 2.0 * (x - 0.5))
                    .collect();
                let etak: Vec<f64> = eta.iter()
                    .map(|&e| 2.0 * (e - 0.5))
                    .collect();
                let psik: Vec<f64> = psi.iter()
                    .map(|&p| 2.0 * (p - 0.5))
                    .collect();

                let xim: Vec<f64> = xik.iter()
                    .map(|&x2| 1.0 - x2)
                    .collect();
                let xip: Vec<f64> = xik.iter()
                    .map(|&x2| 1.0 + x2)
                    .collect();
                let etam: Vec<f64> = etak.iter()
                    .map(|&e2| 1.0 - e2)
                    .collect();
                let etap: Vec<f64> = etak.iter()
                    .map(|&e2| 1.0 + e2)
                    .collect();
                let psim: Vec<f64> = psik.iter()
                    .map(|&p2| 1.0 - p2)
                    .collect();
                let psip: Vec<f64> = psik.iter()
                    .map(|&p2| 1.0 + p2)
                    .collect();

                let xi2: Vec<f64> = xik.iter()
                    .map(|&x2| 1.0 - x2 * x2)
                    .collect();
                let eta2: Vec<f64> = etak.iter()
                    .map(|&e2| 1.0 - e2 * e2)
                    .collect();
                let psi2: Vec<f64> = psik.iter()
                    .map(|&p2| 1.0 - p2 * p2)
                    .collect();

                let values = vec![
                    // Corner nodes
                    0.125 * xim[0] * etam[0] * psim[0] * (-xik[0] - etak[0] - psik[0] - 2.0),     // N0: 
                    0.125 * xip[1] * etam[1] * psim[1] * (xik[1] - etak[1] - psik[1] - 2.0),      // N1: 
                    0.125 * xip[2] * etap[2] * psim[2] * (xik[2] + etak[2] - psik[2] - 2.0),      // N2: 
                    0.125 * xim[3] * etap[3] * psim[3] * (-xik[3] + etak[3] - psik[3] - 2.0),     // N3: 
                    0.125 * xim[4] * etam[4] * psip[4] * (-xik[4] - etak[4] + psik[4] - 2.0),     // N4: 
                    0.125 * xip[5] * etam[5] * psip[5] * (xik[5] - etak[5] + psik[5] - 2.0),      // N5: 
                    0.125 * xip[6] * etap[6] * psip[6] * (xik[6] + etak[6] + psik[6] - 2.0),      // N6: 
                    0.125 * xim[7] * etap[7] * psip[7] * (-xik[7] + etak[7] + psik[7] - 2.0),     // N7: 
                    // Mid-edge nodes
                    0.25 * xi2[8] * etam[8] * psim[8],       // N8:  
                    0.25 * eta2[9] * xip[9] * psim[9],       // N9:
                    0.25 * xi2[10] * etap[10] * psim[10],       // N10:
                    0.25 * eta2[11] * xim[11] * psim[11],       // N11:
                    0.25 * xi2[12] * etam[12] * psip[12],       // N12:
                    0.25 * eta2[13] * xip[13] * psip[13],       // N13:
                    0.25 * xi2[14] * etap[14] * psip[14],       // N14:
                    0.25 * eta2[15] * xim[15] * psip[15],       // N15:
                    0.25 * psi2[16] * xim[16] * etam[16],       // N16:
                    0.25 * psi2[17] * xip[17] * etam[17],       // N17:
                    0.25 * psi2[18] * xip[18] * etap[18],       // N18:
                    0.25 * psi2[19] * xim[19] * etap[19],       // N19:
                ];
                let derivatives = vec![
                    vec![2.0 * (-0.125 * (etam[0] * psim[0] - 2.0 * xik[0] * etam[0] * psim[0] - etak[0] * etam[0] * psim[0] - psik[0] * etam[0] * psim[0] - 2.0 * etam[0] * psim[0])), 
                        2.0 * (-0.125 * (xim[0] * psim[0] - 2.0 * etak[0] * xim[0] * psim[0] - xik[0] * xim[0] * psim[0] - psik[0] * xim[0] * psim[0] - 2.0 * xim[0] * psim[0])), 
                        2.0 * (-0.125 * (xim[0] * etam[0] - 2.0 * psik[0] * xim[0] * etam[0] - xik[0] * xim[0] * etam[0] - etak[0] * xim[0] * etam[0] - 2.0 * xim[0] * etam[0]))],            // dN0
                    vec![2.0 * (0.125 * (etam[1] * psim[1] + 2.0 * xik[1] * etam[1] * psim[1] - etak[1] * etam[1] * psim[1] - psik[1] * etam[1] * psim[1] - 2.0 * etam[1] * psim[1])), 
                        2.0 * (-0.125 * (xip[1] * psim[1] - 2.0 * etak[1] * xip[1] * psim[1] + xik[1] * xip[1] * psim[1] - psik[1] * xip[1] * psim[1] - 2.0 * xip[1] * psim[1])), 
                        2.0 * (-0.125 * (xip[1] * etam[1] - 2.0 * psik[1] * xip[1] * etam[1] + xik[1] * xip[1] * etam[1] - etak[1] * xip[1] * etam[1] - 2.0 * xip[1] * etam[1]))],            // dN1
                    vec![2.0 * (0.125 * (etap[2] * psim[2] + 2.0 * xik[2] * etap[2] * psim[2] + etak[2] * etap[2] * psim[2] - psik[2] * etap[2] * psim[2] - 2.0 * etap[2] * psim[2])), 
                        2.0 * (0.125 * (xip[2] * psim[2] + 2.0 * etak[2] * xip[2] * psim[2] + xik[2] * xip[2] * psim[2] - psik[2] * xip[2] * psim[2] - 2.0 * xip[2] * psim[2])), 
                        2.0 * (-0.125 * (xip[2] * etap[2] - 2.0 * psik[2] * xip[2] * etap[2] + xik[2] * xip[2] * etap[2] + etak[2] * xip[2] * etap[2] - 2.0 * xip[2] * etap[2]))],            // dN2
                    vec![2.0 * (-0.125 * (etap[3] * psim[3] - 2.0 * xik[3] * etap[3] * psim[3] + etak[3] * etap[3] * psim[3] - psik[3] * etap[3] * psim[3] - 2.0 * etap[3] * psim[3])), 
                        2.0 * (0.125 * (xim[3] * psim[3] + 2.0 * etak[3] * xim[3] * psim[3] - xik[3] * xim[3] * psim[3] - psik[3] * xim[3] * psim[3] - 2.0 * xim[3] * psim[3])), 
                        2.0 * (-0.125 * (xim[3] * etap[3] - 2.0 * psik[3] * xim[3] * etap[3] - xik[3] * xim[3] * etap[3] + etak[3] * xim[3] * etap[3] - 2.0 * xim[3] * etap[3]))],            // dN3
                    vec![2.0 * (-0.125 * (etam[4] * psip[4] - 2.0 * xik[4] * etam[4] * psip[4] - etak[4] * etam[4] * psip[4] + psik[4] * etam[4] * psip[4] - 2.0 * etam[4] * psip[4])), 
                        2.0 * (-0.125 * (xim[4] * psip[4] - 2.0 * etak[4] * xim[4] * psip[4] - xik[4] * xim[4] * psip[4] + psik[4] * xim[4] * psip[4] - 2.0 * xim[4] * psip[4])), 
                        2.0 * (0.125 * (xim[4] * etam[4] + 2.0 * psik[4] * xim[4] * etam[4] - xik[4] * xim[4] * etam[4] - etak[4] * xim[4] * etam[4] - 2.0 * xim[4] * etam[4]))],             // dN4
                    vec![2.0 * (0.125 * (etam[5] * psip[5] + 2.0 * xik[5] * etam[5] * psip[5] - etak[5] * etam[5] * psip[5] + psik[5] * etam[5] * psip[5] - 2.0 * etam[5] * psip[5])), 
                        2.0 * (-0.125 * (xip[5] * psip[5] - 2.0 * etak[5] * xip[5] * psip[5] + xik[5] * xip[5] * psip[5] + psik[5] * xip[5] * psip[5] - 2.0 * xip[5] * psip[5])), 
                        2.0 * (0.125 * (xip[5] * etam[5] + 2.0 * psik[5] * xip[5] * etam[5] + xik[5] * xip[5] * etam[5] - etak[5] * xip[5] * etam[5] - 2.0 * xip[5] * etam[5]))],             // dN5
                    vec![2.0 * (0.125 * (etap[6] * psip[6] + 2.0 * xik[6] * etap[6] * psip[6] + etak[6] * etap[6] * psip[6] + psik[6] * etap[6] * psip[6] - 2.0 * etap[6] * psip[6])), 
                        2.0 * (0.125 * (xip[6] * psip[6] + 2.0 * etak[6] * xip[6] * psip[6] + xik[6] * xip[6] * psip[6] + psik[6] * xip[6] * psip[6] - 2.0 * xip[6] * psip[6])), 
                        2.0 * (0.125 * (xip[6] * etap[6] + 2.0 * psik[6] * xip[6] * etap[6] + xik[6] * xip[6] * etap[6] + etak[6] * xip[6] * etap[6] - 2.0 * xip[6] * etap[6]))],             // dN6
                    vec![2.0 * (-0.125 * (etap[7] * psip[7] - 2.0 * xik[7] * etap[7] * psip[7] + etak[7] * etap[7] * psip[7] + psik[7] * etap[7] * psip[7] - 2.0 * etap[7] * psip[7])), 
                        2.0 * (0.125 * (xim[7] * psip[7] + 2.0 * etak[7] * xim[7] * psip[7] - xik[7] * xim[7] * psip[7] + psik[7] * xim[7] * psip[7] - 2.0 * xim[7] * psip[7])), 
                        2.0 * (0.125 * (xim[7] * etap[7] + 2.0 * psik[7] * xim[7] * etap[7] - xik[7] * xim[7] * etap[7] + etak[7] * xim[7] * etap[7] - 2.0 * xim[7] * etap[7]))],             // dN7
                    vec![2.0 * (-0.5 * xik[8] * etam[8] * psim[8]), 2.0 * (-0.25 * (psim[8] - xik[8] * xik[8] * psim[8])), 2.0 * (-0.25 * (etam[8] - xik[8] * xik[8] * etam[8]))],       // dN8
                    vec![2.0 * (0.25 * (psim[9] - etak[9] * etak[9] * psim[9])), 2.0 * (-0.5 * etak[9] * xip[9] * psim[9]), 2.0 * (-0.25 * (xip[9] - etak[9] * etak[9] * xip[9]))],      // dN9
                    vec![2.0 * (-0.5 * xik[10] * etap[10] * psim[10]), 2.0 * (0.25 * (psim[10] - xik[10] * xik[10] * psim[10])), 2.0 * (-0.25 * (etap[10] - xik[10] * xik[10] * etap[10]))],        // dN10
                    vec![2.0 * (-0.25 * (psim[11] - etak[11] * etak[11] * psim[11])), 2.0 * (-0.5 * etak[11] * xim[11] * psim[11]), 2.0 * (-0.25 * (xim[11] - etak[11] * etak[11] * xim[11]))],     // dN11
                    vec![2.0 * (-0.5 * xik[12] * etam[12] * psip[12]), 2.0 * (-0.25 * (psip[12] - xik[12] * xik[12] * psip[12])), 2.0 * (0.25 * (etam[12] - xik[12] * xik[12] * etam[12]))],        // dN12
                    vec![2.0 * (0.25 * (psip[13] - etak[13] * etak[13] * psip[13])), 2.0 * (-0.5 * etak[13] * xip[13] * psip[13]), 2.0 * (0.25 * (xip[13] - etak[13] * etak[13] * xip[13]))],       // dN13
                    vec![2.0 * (-0.5 * xik[14] * etap[14] * psip[14]), 2.0 * (0.25 * (psip[14] - xik[14] * xik[14] * psip[14])), 2.0 * (0.25 * (etap[14] - xik[14] * xik[14] * etap[14]))],         // dN14
                    vec![2.0 * (-0.25 * (psip[15] - etak[15] * etak[15] * psip[15])), 2.0 * (-0.5 * etak[15] * xim[15] * psip[15]), 2.0 * (0.25 * (xim[15] - etak[15] * etak[15] * xim[15]))],      // dN15
                    vec![2.0 * (-0.25 * (etam[16] - psik[16] * psik[16] * etam[16])), 2.0 * (-0.25 * (xim[16] - psik[16] * psik[16] * xim[16])), 2.0 * (-0.5 * psik[16] * xim[16] * etam[16])],     // dN16
                    vec![2.0 * (0.25 * (etam[17] - psik[17] * psik[17] * etam[17])), 2.0 * (-0.25 * (xip[17] - psik[17] * psik[17] * xip[17])), 2.0 * (-0.5 * psik[17] * xip[17] * etam[17])],      // dN17
                    vec![2.0 * (0.25 * (etap[18] - psik[18] * psik[18] * etap[18])), 2.0 * (0.25 * (xip[18] - psik[18] * psik[18] * xip[18])), 2.0 * (-0.5 * psik[18] * xip[18] * etap[18])],       // dN18
                    vec![2.0 * (-0.25 * (etap[19] - psik[19] * psik[19] * etap[19])), 2.0 * (0.25 * (xim[19] - psik[19] * psik[19] * xim[19])), 2.0 * (-0.5 * psik[19] * xim[19] * etap[19])],      // dN19
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let xik: Vec<f64> = xi.iter()
                    .map(|&x| 2.0 * (x - 0.5))
                    .collect();
                let etak: Vec<f64> = eta.iter()
                    .map(|&e| 2.0 * (e - 0.5))
                    .collect();
                let psik: Vec<f64> = psi.iter()
                    .map(|&p| 2.0 * (p - 0.5))
                    .collect();

                let values = vec![
                    // Corner nodes
                    (0.25 * (xik[0] * (1.0 - xik[0])) * (etak[0] * (1.0 - etak[0])) - 0.25 * (1.0 + xik[0]) * (1.0 - xik[0]) * (1.0 + etak[0]) * (1.0 - etak[0])) * (-0.5 * psik[0] * (1.0 - psik[0])),       // N0:
                    (-0.25 * (xik[1] * (1.0 + xik[1])) * (etak[1] * (1.0 - etak[1])) - 0.25 * (1.0 + xik[1]) * (1.0 - xik[1]) * (1.0 + etak[1]) * (1.0 - etak[1])) * (-0.5 * psik[1] * (1.0 - psik[1])),      // N1:
                    (0.25 * (xik[2] * (1.0 + xik[2])) * (etak[2] * (1.0 + etak[2])) - 0.25 * (1.0 + xik[2]) * (1.0 - xik[2]) * (1.0 + etak[2]) * (1.0 - etak[2])) * (-0.5 * psik[2] * (1.0 - psik[2])),       // N2:
                    (-0.25 * (xik[3] * (1.0 - xik[3])) * (etak[3] * (1.0 + etak[3])) - 0.25 * (1.0 + xik[3]) * (1.0 - xik[3]) * (1.0 + etak[3]) * (1.0 - etak[3])) * (-0.5 * psik[3] * (1.0 - psik[3])),      // N3:
                    (0.25 * (xik[4] * (1.0 - xik[4])) * (etak[4] * (1.0 - etak[4])) - 0.25 * (1.0 + xik[4]) * (1.0 - xik[4]) * (1.0 + etak[4]) * (1.0 - etak[4])) * (0.5 * psik[4] * (1.0 + psik[4])),        // N4: 
                    (-0.25 * (xik[5] * (1.0 + xik[5])) * (etak[5] * (1.0 - etak[5])) - 0.25 * (1.0 + xik[5]) * (1.0 - xik[5]) * (1.0 + etak[5]) * (1.0 - etak[5])) * (0.5 * psik[5] * (1.0 + psik[5])),       // N5: 
                    (0.25 * (xik[6] * (1.0 + xik[6])) * (etak[6] * (1.0 + etak[6])) - 0.25 * (1.0 + xik[6]) * (1.0 - xik[6]) * (1.0 + etak[6]) * (1.0 - etak[6])) * (0.5 * psik[6] * (1.0 + psik[6])),        // N6:
                    (-0.25 * (xik[7] * (1.0 - xik[7])) * (etak[7] * (1.0 + etak[7])) - 0.25 * (1.0 + xik[7]) * (1.0 - xik[7]) * (1.0 + etak[7]) * (1.0 - etak[7])) * (0.5 * psik[7] * (1.0 + psik[7])),       // N7:
                    // Mid-edge nodes
                    0.5 * ((1.0 + xik[8]) * (1.0 - xik[8])) * (1.0 - etak[8]) * (-0.5 * psik[8] * (1.0 - psik[8])),        // N8:
                    0.5 * ((1.0 + etak[9]) * (1.0 - etak[9])) * (1.0 + xik[9]) * (-0.5 * psik[9] * (1.0 - psik[9])),       // N9:
                    0.5 * ((1.0 + xik[10]) * (1.0 - xik[10])) * (1.0 + etak[10]) * (-0.5 * psik[10] * (1.0 - psik[10])),        // N10:
                    0.5 * ((1.0 + etak[11]) * (1.0 - etak[11])) * (1.0 - xik[11]) * (-0.5 * psik[11] * (1.0 - psik[11])),       // N11:
                    0.5 * ((1.0 + xik[12]) * (1.0 - xik[12])) * (1.0 - etak[12]) * (0.5 * psik[12] * (1.0 + psik[12])),         // N12:
                    0.5 * ((1.0 + etak[13]) * (1.0 - etak[13])) * (1.0 + xik[13]) * (0.5 * psik[13] * (1.0 + psik[13])),        // N13:
                    0.5 * ((1.0 + xik[14]) * (1.0 - xik[14])) * (1.0 + etak[14]) * (0.5 * psik[14] * (1.0 + psik[14])),         // N14:
                    0.5 * ((1.0 + etak[15]) * (1.0 - etak[15])) * (1.0 - xik[15]) * (0.5 * psik[15] * (1.0 + psik[15])),        // N15:
                    (0.25 * (xik[16] * (1.0 - xik[16])) * (etak[16] * (1.0 - etak[16])) - 0.25 * (1.0 + xik[16]) * (1.0 - xik[16]) * (1.0 + etak[16]) * (1.0 - etak[16])) * ((1.0 + psik[16]) * (1.0 - psik[16])),        // N16:
                    (-0.25 * (xik[17] * (1.0 + xik[17])) * (etak[17] * (1.0 - etak[17])) - 0.25 * (1.0 + xik[17]) * (1.0 - xik[17]) * (1.0 + etak[17]) * (1.0 - etak[17])) * ((1.0 + psik[17]) * (1.0 - psik[17])),       // N17:
                    (0.25 * (xik[18] * (1.0 + xik[18])) * (etak[18] * (1.0 + etak[18])) - 0.25 * (1.0 + xik[18]) * (1.0 - xik[18]) * (1.0 + etak[18]) * (1.0 - etak[18])) * ((1.0 + psik[18]) * (1.0 - psik[18])),        // N18:
                    (-0.25 * (xik[19] * (1.0 - xik[19])) * (etak[19] * (1.0 + etak[19])) - 0.25 * (1.0 + xik[19]) * (1.0 - xik[19]) * (1.0 + etak[19]) * (1.0 - etak[19])) * ((1.0 + psik[19]) * (1.0 - psik[19])),       // N19:
                    // Mid-face nodes
                    0.5 * ((1.0 + etak[20]) * (1.0 - etak[20])) * (1.0 - xik[20]) * ((1.0 + psik[20]) * (1.0 - psik[20])),        // N20:
                    0.5 * ((1.0 + etak[21]) * (1.0 - etak[21])) * (1.0 + xik[21]) * ((1.0 + psik[21]) * (1.0 - psik[21])),        // N21:
                    0.5 * ((1.0 + xik[22]) * (1.0 - xik[22])) * (1.0 - etak[22]) * ((1.0 + psik[22]) * (1.0 - psik[22])),         // N22:
                    0.5 * ((1.0 + xik[23]) * (1.0 - xik[23])) * (1.0 + etak[23]) * ((1.0 + psik[23]) * (1.0 - psik[23])),         // N23:
                ];
                let derivatives = vec![
                    vec![2.0 * (-((etak[0] * etak[0] + (2.0 * xik[0] - 1.0) * etak[0] - 2.0 * xik[0]) * psik[0] * psik[0] + (-etak[0] * etak[0] + (1.0 - 2.0 * xik[0]) * etak[0] + 2.0 * xik[0]) * psik[0]) / 8.0), 
                        2.0 * (-(((2.0 * xik[0] - 2.0) * etak[0] + xik[0] * xik[0] - xik[0]) * psik[0] * psik[0] + ((2.0 - 2.0 * xik[0]) * etak[0] - xik[0] * xik[0] + xik[0]) * psik[0]) / 8.0), 
                        2.0 * (-(((2.0 * xik[0] - 2.0) * etak[0] * etak[0] + (2.0 * xik[0] * xik[0] - 2.0 * xik[0]) * etak[0] - 2.0 * xik[0] * xik[0] + 2.0) * psik[0] + (1.0 - xik[0]) * etak[0] * etak[0] + (xik[0] - xik[0] * xik[0]) * etak[0] + xik[0] * xik[0] - 1.0) / 8.0)],     // dN0
                    vec![2.0 * (((etak[1] * etak[1] + (-2.0 * xik[1] - 1.0) * etak[1] + 2.0 * xik[1]) * psik[1] * psik[1] + (-etak[1] * etak[1] + (2.0 * xik[1] + 1.0) * etak[1] - 2.0 * xik[1]) * psik[1]) / 8.0), 
                        2.0 * ((((2.0 * xik[1] + 2.0) * etak[1] - xik[1] * xik[1] - xik[1]) * psik[1] * psik[1] + ((-2.0 * xik[1] - 2.0) * etak[1] + xik[1] * xik[1] + xik[1]) * psik[1]) / 8.0),
                        2.0 * ((((2.0 * xik[1] + 2.0) * etak[1] * etak[1] + (-2.0 * xik[1] * xik[1] - 2.0 * xik[1]) * etak[1] + 2.0 * xik[1] * xik[1] - 2.0) * psik[1] + (-xik[1] - 1.0) * etak[1] * etak[1] + (xik[1] * xik[1] + xik[1]) * etak[1] - xik[1] * xik[1] + 1.0) / 8.0)],      // dN1
                    vec![2.0 * (((etak[2] * etak[2] + (2.0 * xik[2] + 1.0) * etak[2] + 2.0 * xik[2]) * psik[2] * psik[2] + (-etak[2] * etak[2] + (-2.0 * xik[2] - 1.0) * etak[2] - 2.0 * xik[2]) * psik[2]) / 8.0), 
                        2.0 * ((((2.0 * xik[2] + 2.0) * etak[2] + xik[2] * xik[2] + xik[2]) * psik[2] * psik[2] + ((-2.0 * xik[2] - 2.0) * etak[2] - xik[2] * xik[2] - xik[2]) * psik[2]) / 8.0),
                        2.0 * ((((2.0 * xik[2] + 2.0) * etak[2] * etak[2] + (2.0 * xik[2] * xik[2] + 2.0 * xik[2]) * etak[2] + 2.0 * xik[2] * xik[2] - 2.0) * psik[2] + (-xik[2] - 1.0) * etak[2] * etak[2] + (-xik[2] * xik[2] - xik[2]) * etak[2] - xik[2] * xik[2] + 1.0) / 8.0)],       // dN2
                    vec![2.0 * (-((etak[3] * etak[3] + (1.0 - 2.0 * xik[3]) * etak[3] - 2.0 * xik[3]) * psik[3] * psik[3] + (-etak[3] * etak[3] + (2.0 * xik[3] - 1.0) * etak[3] + 2.0 * xik[3]) * psik[3]) / 8.0),
                        2.0 * (-(((2.0 * xik[3] - 2.0) * etak[3] - xik[3] * xik[3] + xik[3]) * psik[3] * psik[3] + ((2.0 - 2.0 * xik[3]) * etak[3] + xik[3] * xik[3] - xik[3]) * psik[3]) / 8.0), 
                        2.0 * (-(((2.0 * xik[3] - 2.0) * etak[3] * etak[3] + (2.0 * xik[3] - 2.0 * xik[3] * xik[3]) * etak[3] - 2.0 * xik[3] * xik[3] + 2.0) * psik[3] + (1.0 - xik[3]) * etak[3] * etak[3] + (xik[3] * xik[3] - xik[3]) * etak[3] + xik[3] * xik[3] - 1.0) / 8.0)],      // dN3
                    vec![2.0 * (-((etak[4] * etak[4] + (2.0 * xik[4] - 1.0) * etak[4] - 2.0 * xik[4]) * psik[4] * psik[4] + (etak[4] * etak[4] + (2.0 * xik[4] - 1.0) * etak[4] - 2.0 * xik[4]) * psik[4]) / 8.0), 
                        2.0 * (-(((2.0 * xik[4] - 2.0) * etak[4] + xik[4] * xik[4] - xik[4]) * psik[4] * psik[4] + ((2.0 * xik[4] - 2.0) * etak[4] + xik[4] * xik[4] - xik[4]) * psik[4]) / 8.0), 
                        2.0 * (-(((2.0 * xik[4] - 2.0) * etak[4] * etak[4] + (2.0 * xik[4] * xik[4] - 2.0 * xik[4]) * etak[4] - 2.0 * xik[4] * xik[4] + 2.0) * psik[4] + (xik[4] - 1.0) * etak[4] * etak[4] + (xik[4] * xik[4] - xik[4]) * etak[4] - xik[4] * xik[4] + 1.0) / 8.0)],      // dN4
                    vec![2.0 * (((etak[5] * etak[5] + (-2.0 * xik[5] - 1.0) * etak[5] + 2.0 * xik[5]) * psik[5] * psik[5] + (etak[5] * etak[5] + (-2.0 * xik[5] - 1.0) * etak[5] + 2.0 * xik[5]) * psik[5]) / 8.0), 
                        2.0 * ((((2.0 * xik[5] + 2.0) * etak[5] - xik[5] * xik[5] - xik[5]) * psik[5] * psik[5] + ((2.0 * xik[5] + 2.0) * etak[5] - xik[5] * xik[5] - xik[5]) * psik[5]) / 8.0), 
                        2.0 * ((((2.0 * xik[5] + 2.0) * etak[5] * etak[5] + (-2.0 * xik[5] * xik[5] - 2.0 * xik[5]) * etak[5] + 2.0 * xik[5] * xik[5] - 2.0) * psik[5] + (xik[5] + 1.0) * etak[5] * etak[5] + (-xik[5] * xik[5] - xik[5]) * etak[5] + xik[5] * xik[5] - 1.0) / 8.0)],       // dN5
                    vec![2.0 * (((etak[6] * etak[6] + (2.0 * xik[6] + 1.0) * etak[6] + 2.0 * xik[6]) * psik[6] * psik[6] + (etak[6] * etak[6] + (2.0 * xik[6] + 1.0) * etak[6] + 2.0 * xik[6]) * psik[6]) / 8.0), 
                        2.0 * ((((2.0 * xik[6] + 2.0) * etak[6] + xik[6] * xik[6] + xik[6]) * psik[6] * psik[6] + ((2.0 * xik[6] + 2.0) * etak[6] + xik[6] * xik[6] + xik[6]) * psik[6]) / 8.0), 
                        2.0 * ((((2.0 * xik[6] + 2.0) * etak[6] * etak[6] + (2.0 * xik[6] * xik[6] + 2.0 * xik[6]) * etak[6] + 2.0 * xik[6] * xik[6] - 2.0) * psik[6] + (xik[6] + 1.0) * etak[6] * etak[6] + (xik[6] * xik[6] + xik[6]) * etak[6] + xik[6] * xik[6] - 1.0) / 8.0)],        // dN6
                    vec![2.0 * (-((etak[7] * etak[7] + (1.0 - 2.0 * xik[7]) * etak[7] - 2.0 * xik[7]) * psik[7] * psik[7] + (etak[7] * etak[7] + (1.0 - 2.0 * xik[7]) * etak[7] - 2.0 * xik[7]) * psik[7]) / 8.0), 
                        2.0 * (-(((2.0 * xik[7] - 2.0) * etak[7] - xik[7] * xik[7] + xik[7]) * psik[7] * psik[7] + ((2.0 * xik[7] - 2.0) * etak[7] - xik[7] * xik[7] + xik[7]) * psik[7]) / 8.0), 
                        2.0 * (-(((2.0 * xik[7] - 2.0) * etak[7] * etak[7] + (2.0 * xik[7] - 2.0 * xik[7] * xik[7]) * etak[7] - 2.0 * xik[7] * xik[7] + 2.0) * psik[7] + (1.0 - xik[7]) * etak[7] * etak[7] + (xik[7] - xik[7] * xik[7]) * etak[7] - xik[7] * xik[7] + 1.0) / 8.0)],       // dN7
                    vec![2.0 * (((xik[8] * etak[8] - xik[8]) * psik[8] * psik[8] + (xik[8] - xik[8] * etak[8]) * psik[8]) / 2.0), 
                        2.0 * (((xik[8] * xik[8] - 1.0) * psik[8] * psik[8] + (1.0 - xik[8] * xik[8]) * psik[8]) / 4.0), 
                        2.0 * ((((2.0 * xik[8] * xik[8] - 2.0) * etak[8] - 2.0 * xik[8] * xik[8] + 2.0) * psik[8] + (1.0 - xik[8] * xik[8]) * etak[8] + xik[8] * xik[8] - 1.0) / 4.0)],       // dN8
                    vec![2.0 * (-((etak[9] * etak[9] - 1.0) * psik[9] * psik[9] + (1.0 - etak[9] * etak[9]) * psik[9]) / 4.0), 
                        2.0 * (-((xik[9] + 1.0) * etak[9] * psik[9] * psik[9] + (-xik[9] - 1.0) * etak[9] * psik[9]) / 2.0), 
                        2.0 * (-(((2.0 * xik[9] + 2.0) * etak[9] * etak[9] - 2.0 * xik[9] - 2.0) * psik[9] + (-xik[9] - 1.0) * etak[9] * etak[9] + xik[9] + 1.0) / 4.0)],      // dN9
                    vec![2.0 * (-((xik[10] * etak[10] + xik[10]) * psik[10] * psik[10] + (-xik[10] * etak[10] - xik[10]) * psik[10]) / 2.0), 
                        2.0 * (-((xik[10] * xik[10] - 1.0) * psik[10] * psik[10] + (1.0 - xik[10] * xik[10]) * psik[10]) / 4.0), 
                        2.0 * (-(((2.0 * xik[10] * xik[10] - 2.0) * etak[10] + 2.0 * xik[10] * xik[10] - 2.0) * psik[10] + (1.0 - xik[10] * xik[10]) * etak[10] - xik[10] * xik[10] + 1.0) / 4.0)],        // dN10
                    vec![2.0 * (((etak[11] * etak[11] - 1.0) * psik[11] * psik[11] + (1.0 - etak[11] * etak[11]) * psik[11]) / 4.0), 
                        2.0 * (((xik[11] - 1.0) * etak[11] * psik[11] * psik[11] + (1.0 - xik[11]) * etak[11] * psik[11]) / 2.0), 
                        2.0 * ((((2.0 * xik[11] - 2.0) * etak[11] * etak[11] - 2.0 * xik[11] + 2.0) * psik[11] + (1.0 - xik[11]) * etak[11] * etak[11] + xik[11] - 1.0) / 4.0)],     // dN11
                    vec![2.0 * (((xik[12] * etak[12] - xik[12]) * psik[12] * psik[12] + (xik[12] * etak[12] - xik[12]) * psik[12]) / 2.0), 
                        2.0 * (((xik[12] * xik[12] - 1.0) * psik[12] * psik[12] + (xik[12] * xik[12] - 1.0) * psik[12]) / 4.0), 
                        2.0 * ((((2.0 * xik[12] * xik[12] - 2.0) * etak[12] - 2.0 * xik[12] * xik[12] + 2.0) * psik[12] + (xik[12] * xik[12] - 1.0) * etak[12] - xik[12] * xik[12] + 1.0) / 4.0)],        // dN12
                    vec![2.0 * (-((etak[13] * etak[13] - 1.0) * psik[13] * psik[13] + (etak[13] * etak[13] - 1.0) * psik[13]) / 4.0), 
                        2.0 * (-((xik[13] + 1.0) * etak[13] * psik[13] * psik[13] + (xik[13] + 1.0) * etak[13] * psik[13]) / 2.0), 
                        2.0 * (-(((2.0 * xik[13] + 2.0) * etak[13] * etak[13] - 2.0 * xik[13] - 2.0) * psik[13] + (xik[13] + 1.0) * etak[13] * etak[13] - xik[13] - 1.0) / 4.0)],       // dN13
                    vec![2.0 * (-((xik[14] * etak[14] + xik[14]) * psik[14] * psik[14] + (xik[14] * etak[14] + xik[14]) * psik[14]) / 2.0), 
                        2.0 * (-((xik[14] * xik[14] - 1.0) * psik[14] * psik[14] + (xik[14] * xik[14] - 1.0) * psik[14]) / 4.0), 
                        2.0 * (-(((2.0 * xik[14] * xik[14] - 2.0) * etak[14] + 2.0 * xik[14] * xik[14] - 2.0) * psik[14] + (xik[14] * xik[14] - 1.0) * etak[14] + xik[14] * xik[14] - 1.0) / 4.0)],         // dN14
                    vec![2.0 * (((etak[15] * etak[15] + (2.0 * xik[15] - 1.0) * etak[15] - 2.0 * xik[15]) * psik[15] * psik[15] - etak[15] * etak[15] + (1.0 - 2.0 * xik[15]) * etak[15] + 2.0 * xik[15]) / 4.0), 
                        2.0 * ((((2.0 * xik[15] - 2.0) * etak[15] + xik[15] * xik[15] - xik[15]) * psik[15] * psik[15] + (2.0 - 2.0 * xik[15]) * etak[15] - xik[15] * xik[15] + xik[15]) / 4.0), 
                        2.0 * (((xik[15] - 1.0) * etak[15] * etak[15] + (xik[15] * xik[15] - xik[15]) * etak[15] - xik[15] * xik[15] + 1.0) * psik[15] / 2.0)],     // dN15
                    vec![2.0 * (-((etak[16] * etak[16] + (-2.0 * xik[16] - 1.0) * etak[16] + 2.0 * xik[16]) * psik[16] * psik[16] - etak[16] * etak[16] + (2.0 * xik[16] + 1.0) * etak[16] - 2.0 * xik[16]) / 4.0), 
                        2.0 * (-(((2.0 * xik[16] + 2.0) * etak[16] - xik[16] * xik[16] - xik[16]) * psik[16] * psik[16] + (-2.0 * xik[16] - 2.0) * etak[16] + xik[16] * xik[16] + xik[16]) / 4.0), 
                        2.0 * (-((xik[16] + 1.0) * etak[16] * etak[16] + (-xik[16] * xik[16] - xik[16]) * etak[16] + xik[16] * xik[16] - 1.0) * psik[16] / 2.0)],      // dN16
                    vec![2.0 * (-((etak[17] * etak[17] + (2.0 * xik[17] + 1.0) * etak[17] + 2.0 * xik[17]) * psik[17] * psik[17] - etak[17] * etak[17] + (-2.0 * xik[17] - 1.0) * etak[17] - 2.0 * xik[17]) / 4.0), 
                        2.0 * (-(((2.0 * xik[17] + 2.0) * etak[17] + xik[17] * xik[17] + xik[17]) * psik[17] * psik[17] + (-2.0 * xik[17] - 2.0) * etak[17] - xik[17] * xik[17] - xik[17]) / 4.0), 
                        2.0 * (-((xik[17] + 1.0) * etak[17] * etak[17] + (xik[17] * xik[17] + xik[17]) * etak[17] + xik[17] * xik[17] - 1.0) * psik[17] / 2.0)],       // dN17
                    vec![2.0 * (((etak[18] * etak[18] + (1.0 - 2.0 * xik[18]) * etak[18] - 2.0 * xik[18]) * psik[18] * psik[18] - etak[18] * etak[18] + (2.0 * xik[18] - 1.0) * etak[18] + 2.0 * xik[18]) / 4.0), 
                        2.0 * ((((2.0 * xik[18] - 2.0) * etak[18] - xik[18] * xik[18] + xik[18]) * psik[18] * psik[18] + (2.0 - 2.0 * xik[18]) * etak[18] + xik[18] * xik[18] - xik[18]) / 4.0), 
                        2.0 * (((xik[18] - 1.0) * etak[18] * etak[18] + (xik[18] - xik[18] * xik[18]) * etak[18] - xik[18] * xik[18] + 1.0) * psik[18] / 2.0)],      // dN18
                    vec![2.0 * (-((etak[19] * etak[19] - 1.0) * psik[19] * psik[19] - etak[19] * etak[19] + 1.0) / 2.0), 
                        2.0 * ((1.0 - xik[19]) * etak[19] * psik[19] * psik[19] + (xik[19] - 1.0) * etak[19]), 
                        2.0 * (((1.0 - xik[19]) * etak[19] * etak[19] + xik[19] - 1.0) * psik[19])],       // dN19
                    vec![2.0 * (((etak[20] * etak[20] - 1.0) * psik[20] * psik[20] - etak[20] * etak[20] + 1.0) / 2.0), 
                        2.0 * ((xik[20] + 1.0) * etak[20] * psik[20] * psik[20] + (-xik[20] - 1.0) * etak[20]), 
                        2.0 * (((xik[20] + 1.0) * etak[20] * etak[20] - xik[20] - 1.0) * psik[20])],      // dN20
                    vec![2.0 * ((xik[21] - xik[21] * etak[21]) * psik[21] * psik[21] + xik[21] * etak[21] - xik[21]), 
                        2.0 * (-((xik[21] * xik[21] - 1.0) * psik[21] * psik[21] - xik[21] * xik[21] + 1.0) / 2.0), 
                        2.0 * (((1.0 - xik[21] * xik[21]) * etak[21] + xik[21] * xik[21] - 1.0) * psik[21])],       // dN21
                    vec![2.0 * ((xik[22] * etak[22] + xik[22]) * psik[22] * psik[22] - xik[22] * etak[22] - xik[22]), 
                        2.0 * (((xik[22] * xik[22] - 1.0) * psik[22] * psik[22] - xik[22] * xik[22] + 1.0) / 2.0), 
                        2.0 * (((xik[22] * xik[22] - 1.0) * etak[22] + xik[22] * xik[22] - 1.0) * psik[22])],      // dN22
                    vec![2.0 * ((xik[23] * etak[23] + xik[23]) * psik[23] * psik[23] - xik[23] * etak[23] - xik[23]), 
                        2.0 * (((xik[23] * xik[23] - 1.0) * psik[23] * psik[23] - xik[23] * xik[23] + 1.0) / 2.0), 
                        2.0 * (((xik[23] * xik[23] - 1.0) * etak[23] + xik[23] * xik[23] - 1.0) * psik[23])],      // dN23
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

                let xi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[0])
                    .collect();
                let eta: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[1])
                    .collect();
                let psi: Vec<f64> = reference_coords.iter()
                    .map(|coord| coord[2])
                    .collect();

                let xik: Vec<f64> = xi.iter()
                    .map(|&x| 2.0 * (x - 0.5))
                    .collect();
                let etak: Vec<f64> = eta.iter()
                    .map(|&e| 2.0 * (e - 0.5))
                    .collect();
                let psik: Vec<f64> = psi.iter()
                    .map(|&p| 2.0 * (p - 0.5))
                    .collect();

                let xi1: Vec<f64> = xik.iter()
                    .map(|&x2| -0.5 * x2 * (1.0 - x2))
                    .collect();
                let eta1: Vec<f64> = etak.iter()
                    .map(|&e2| -0.5 * e2 * (1.0 - e2))
                    .collect();
                let psi1: Vec<f64> = psik.iter()
                    .map(|&p2| -0.5 * p2 * (1.0 - p2))
                    .collect();

                let xi2: Vec<f64> = xik.iter()
                    .map(|&x2| (1.0 + x2) * (1.0 - x2))
                    .collect();
                let eta2: Vec<f64> = etak.iter()
                    .map(|&e2| (1.0 + e2) * (1.0 - e2))
                    .collect();
                let psi2: Vec<f64> = psik.iter()
                    .map(|&p2| (1.0 + p2) * (1.0 - p2))
                    .collect();

                let xi3: Vec<f64> = xik.iter()
                    .map(|&x2| 0.5 * x2 * (1.0 + x2))
                    .collect();
                let eta3: Vec<f64> = etak.iter()
                    .map(|&e2| 0.5 * e2 * (1.0 + e2))
                    .collect();
                let psi3: Vec<f64> = psik.iter()
                    .map(|&p2| 0.5 * p2 * (1.0 + p2))
                    .collect();

                let xi1_xi: Vec<f64> = xik.iter()
                    .map(|&x2| x2 - 0.5)
                    .collect();
                let eta1_eta: Vec<f64> = etak.iter()
                    .map(|&e2| e2 - 0.5)
                    .collect();
                let psi1_psi: Vec<f64> = psik.iter()
                    .map(|&p2| p2 - 0.5)
                    .collect();

                let xi2_xi: Vec<f64> = xik.iter()
                    .map(|&x2| -2.0 * x2)
                    .collect();
                let eta2_eta: Vec<f64> = etak.iter()
                    .map(|&e2| -2.0 * e2)
                    .collect();
                let psi2_psi: Vec<f64> = psik.iter()
                    .map(|&p2| -2.0 * p2)
                    .collect();

                let xi3_xi: Vec<f64> = xik.iter()
                    .map(|&x2| x2 + 0.5)
                    .collect();
                let eta3_eta: Vec<f64> = etak.iter()
                    .map(|&e2| e2 + 0.5)
                    .collect();
                let psi3_psi: Vec<f64> = psik.iter()
                    .map(|&p2| p2 + 0.5)
                    .collect();

                let values = vec![
                    // Corner nodes
                    xi1[0] * eta1[0] * psi1[0],      // N0:
                    xi3[1] * eta1[1] * psi1[1],      // N1:
                    xi3[2] * eta3[2] * psi1[2],      // N2:
                    xi1[3] * eta3[3] * psi1[3],      // N3:
                    xi1[4] * eta1[4] * psi3[4],      // N4:
                    xi3[5] * eta1[5] * psi3[5],      // N5:
                    xi3[6] * eta3[6] * psi3[6],      // N6:
                    xi1[7] * eta3[7] * psi3[7],      // N7:
                    // Mid-edge nodes
                    xi2[8] * eta1[8] * psi1[8],      // N8:
                    xi3[9] * eta2[9] * psi1[9],      // N9:
                    xi2[10] * eta3[10] * psi1[10],  // N10:
                    xi1[11] * eta2[11] * psi1[11],  // N11:
                    xi2[12] * eta1[12] * psi3[12],  // N12:
                    xi3[13] * eta2[13] * psi3[13],  // N13:
                    xi2[14] * eta3[14] * psi3[14],  // N14:
                    xi1[15] * eta2[15] * psi3[15],  // N15:
                    xi1[16] * eta1[16] * psi2[16],  // N16:
                    xi3[17] * eta1[17] * psi2[17],  // N17:
                    xi3[18] * eta3[18] * psi2[18],  // N18:
                    xi1[19] * eta3[19] * psi2[19],  // N19:
                    // Mid-face nodes
                    xi2[20] * eta1[20] * psi2[20],  // N20:
                    xi3[21] * eta2[21] * psi2[21],  // N21:
                    xi2[22] * eta3[22] * psi2[22],  // N22:
                    xi1[23] * eta2[23] * psi2[23],  // N23:
                    xi2[24] * eta2[24] * psi1[24],  // N24:
                    xi2[25] * eta2[25] * psi3[25],  // N25:
                    // Center node
                    xi2[26] * eta2[26] * psi2[26],  // N26:
                ];
                let derivatives = vec![
                    vec![2.0 * (xi1_xi[0] * eta1[0] * psi1[0]), 2.0 * (xi1[0] * eta1_eta[0] * psi1[0]), 2.0 * (xi1[0] * eta1[0] * psi1_psi[0])],     // dN0
                    vec![2.0 * (xi3_xi[1] * eta1[1] * psi1[1]), 2.0 * (xi3[1] * eta1_eta[1] * psi1[1]), 2.0 * (xi3[1] * eta1[1] * psi1_psi[1])],      // dN1
                    vec![2.0 * (xi3_xi[2] * eta3[2] * psi1[2]), 2.0 * (xi3[2] * eta3_eta[2] * psi1[2]), 2.0 * (xi3[2] * eta3[2] * psi1_psi[2])],       // dN2
                    vec![2.0 * (xi1_xi[3] * eta3[3] * psi1[3]), 2.0 * (xi1[3] * eta3_eta[3] * psi1[3]), 2.0 * (xi1[3] * eta3[3] * psi1_psi[3])],      // dN3
                    vec![2.0 * (xi1_xi[4] * eta1[4] * psi3[4]), 2.0 * (xi1[4] * eta1_eta[4] * psi3[4]), 2.0 * (xi1[4] * eta1[4] * psi3_psi[4])],      // dN4
                    vec![2.0 * (xi3_xi[5] * eta1[5] * psi3[5]), 2.0 * (xi3[5] * eta1_eta[5] * psi3[5]), 2.0 * (xi3[5] * eta1[5] * psi3_psi[5])],       // dN5
                    vec![2.0 * (xi3_xi[6] * eta3[6] * psi3[6]), 2.0 * (xi3[6] * eta3_eta[6] * psi3[6]), 2.0 * (xi3[6] * eta3[6] * psi3_psi[6])],        // dN6
                    vec![2.0 * (xi1_xi[7] * eta3[7] * psi3[7]), 2.0 * (xi1[7] * eta3_eta[7] * psi3[7]), 2.0 * (xi1[7] * eta3[7] * psi3_psi[7])],       // dN7
                    vec![2.0 * (xi2_xi[8] * eta1[8] * psi1[8]), 2.0 * (xi2[8] * eta1_eta[8] * psi1[8]), 2.0 * (xi2[8] * eta1[8] * psi1_psi[8])],       // dN8
                    vec![2.0 * (xi3_xi[9] * eta2[9] * psi1[9]), 2.0 * (xi3[9] * eta2_eta[9] * psi1[9]), 2.0 * (xi3[9] * eta2[9] * psi1_psi[9])],      // dN9
                    vec![2.0 * (xi2_xi[10] * eta3[10] * psi1[10]), 2.0 * (xi2[10] * eta3_eta[10] * psi1[10]), 2.0 * (xi2[10] * eta3[10] * psi1_psi[10])],        // dN10
                    vec![2.0 * (xi1_xi[11] * eta2[11] * psi1[11]), 2.0 * (xi1[11] * eta2_eta[11] * psi1[11]), 2.0 * (xi1[11] * eta2[11] * psi1_psi[11])],     // dN11
                    vec![2.0 * (xi2_xi[12] * eta1[12] * psi3[12]), 2.0 * (xi2[12] * eta1_eta[12] * psi3[12]), 2.0 * (xi2[12] * eta1[12] * psi3_psi[12])],        // dN12
                    vec![2.0 * (xi3_xi[13] * eta2[13] * psi3[13]), 2.0 * (xi3[13] * eta2_eta[13] * psi3[13]), 2.0 * (xi3[13] * eta2[13] * psi3_psi[13])],       // dN13
                    vec![2.0 * (xi2_xi[14] * eta3[14] * psi3[14]), 2.0 * (xi2[14] * eta3_eta[14] * psi3[14]), 2.0 * (xi2[14] * eta3[14] * psi3_psi[14])],         // dN14
                    vec![2.0 * (xi1_xi[15] * eta2[15] * psi3[15]), 2.0 * (xi1[15] * eta2_eta[15] * psi3[15]), 2.0 * (xi1[15] * eta2[15] * psi3_psi[15])],      // dN15
                    vec![2.0 * (xi1_xi[16] * eta1[16] * psi2[16]), 2.0 * (xi1[16] * eta1_eta[16] * psi2[16]), 2.0 * (xi1[16] * eta1[16] * psi2_psi[16])],     // dN16
                    vec![2.0 * (xi3_xi[17] * eta1[17] * psi2[17]), 2.0 * (xi3[17] * eta1_eta[17] * psi2[17]), 2.0 * (xi3[17] * eta1[17] * psi2_psi[17])],      // dN17
                    vec![2.0 * (xi3_xi[18] * eta3[18] * psi2[18]), 2.0 * (xi3[18] * eta3_eta[18] * psi2[18]), 2.0 * (xi3[18] * eta3[18] * psi2_psi[18])],       // dN18
                    vec![2.0 * (xi1_xi[19] * eta3[19] * psi2[19]), 2.0 * (xi1[19] * eta3_eta[19] * psi2[19]), 2.0 * (xi1[19] * eta3[19] * psi2_psi[19])],      // dN19
                    vec![2.0 * (xi1_xi[20] * eta2[20] * psi2[20]), 2.0 * (xi1[20] * eta2_eta[20] * psi2[20]), 2.0 * (xi1[20] * eta2[20] * psi2_psi[20])],       // dN20
                    vec![2.0 * (xi3_xi[21] * eta2[21] * psi2[21]), 2.0 * (xi3[21] * eta2_eta[21] * psi2[21]), 2.0 * (xi3[21] * eta2[21] * psi2_psi[21])],      // dN21
                    vec![2.0 * (xi2_xi[22] * eta1[22] * psi2[22]), 2.0 * (xi2[22] * eta1_eta[22] * psi2[22]), 2.0 * (xi2[22] * eta1[22] * psi2_psi[22])],       // dN22
                    vec![2.0 * (xi2_xi[23] * eta3[23] * psi2[23]), 2.0 * (xi2[23] * eta3_eta[23] * psi2[23]), 2.0 * (xi2[23] * eta3[23] * psi2_psi[23])],      // dN23
                    vec![2.0 * (xi2_xi[24] * eta2[24] * psi1[24]), 2.0 * (xi2[24] * eta2_eta[24] * psi1[24]), 2.0 * (xi2[24] * eta2[24] * psi1_psi[24])],       // dN24
                    vec![2.0 * (xi2_xi[25] * eta2[25] * psi3[25]), 2.0 * (xi2[25] * eta2_eta[25] * psi3[25]), 2.0 * (xi2[25] * eta2[25] * psi3_psi[25])],      // dN25
                    vec![2.0 * (xi2_xi[26] * eta2[26] * psi2[26]), 2.0 * (xi2[26] * eta2_eta[26] * psi2[26]), 2.0 * (xi2[26] * eta2[26] * psi2_psi[26])],       // dN26
                ];

                Ok(ShapeFunction {values, derivatives})
            },

            _ => Err(ParseError::ElementError(format!("Shape functions for element type {:?} are not implemented", element_type))),
        }
    }

    /// Get the nodes associated with an element
    fn get_element_nodes(element: &Element, nodes: &[Node]) -> Result<Vec<Node>, ElementError> {
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
    fn calculate_jacobian(
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
            
               
 



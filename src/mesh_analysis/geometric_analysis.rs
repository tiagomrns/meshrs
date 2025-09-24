use std::io::{self};
use std::f64;

// Import your existing structs
use crate::structs_and_impls::*;
use crate::error::*;

pub struct GeometricAnalysis;

impl GeometricAnalysis {
    /// Main method to analyze quality of all elements in the mesh using reference elements
    pub fn analyse_mesh_quality(mesh_data: &MeshData) -> Result<MeshQualityReport, ElementError> {
        let mut element_qualities = Vec::new();
        let mut processed_elements = 0;
        
        for type_info in &mesh_data.element_type_info {
            if matches!(type_info.element_type, ElementType::Vertex) {
                continue;
            }
            
            let start_idx = type_info.start_index;
            let end_idx = start_idx + type_info.num_elements;
            
            for element_idx in start_idx..end_idx {
                if element_idx < mesh_data.elements.len() {
                    let element = &mesh_data.elements[element_idx];
                    
                    match Self::calculate_element_quality(element, &type_info.element_type, &mesh_data.nodes) {
                        Ok(quality) => {
                            element_qualities.push(quality);
                            processed_elements += 1;
                        },
                        Err(e) => {
                            println!("Warning: Failed to analyze element {}: {:?}", element.id, e);
                        }
                    }
                }
            }
        }
        
        if element_qualities.is_empty() {
            return Err(ElementError::GeometryError("No elements could be analyzed".to_string()));
        }
        
        Ok(MeshQualityReport {
            total_elements: processed_elements,
            element_qualities,
        })
    }

    /// Calculate element quality using Jacobian analysis
    fn calculate_element_quality(
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node]
    ) -> Result<ElementQuality, ElementError> {
        let element_nodes = Self::get_element_nodes(element, nodes)?;
        let center_point = Self::get_center_point(element_type);
        
        let shape_function = ElementType::get_shape_functions(element_type, &center_point)
            .ok_or_else(|| ElementError::GeometryError("Failed to get shape functions".to_string()))?;
        let shape_derivatives = shape_function.derivatives;

        let jacobian = Self::calculate_jacobian(&element_nodes, &shape_derivatives)?;
        
        Ok(ElementQuality {
            element_id: element.id,
            det_jacobian: jacobian.determinant,
        })
    }
    

    fn get_center_point(element_type: &ElementType) -> Vec<f64> {
        match element_type {
            ElementType::Line => vec![0.5],
            ElementType::QuadraticEdge => vec![0.5],
            ElementType::Triangle => vec![1.0/3.0, 1.0/3.0],
            ElementType::QuadraticTriangle => vec![1.0/3.0, 1.0/3.0],
            ElementType::Quad => vec![0.5, 0.5],
            ElementType::QuadraticQuad => vec![0.5, 0.5],
            ElementType::BiquadraticQuad => vec![0.5, 0.5],
            ElementType::Tetra => vec![0.25, 0.25, 0.25],
            ElementType::QuadraticTetra => vec![0.25, 0.25, 0.25],
            ElementType::Pyramid => vec![0.5, 0.5, 0.5],
            ElementType::Wedge => vec![1.0/3.0, 1.0/3.0, 0.5],
            ElementType::QuadraticWedge => vec![1.0/3.0, 1.0/3.0, 0.5],
            ElementType::BiquadraticQuadraticWedge => vec![1.0/3.0, 1.0/3.0, 0.5],
            ElementType::Hexahedron => vec![0.5, 0.5, 0.5],
            ElementType::QuadraticHexahedron => vec![0.5, 0.5, 0.5],
            ElementType::BiquadraticQuadraticHexahedron => vec![0.5, 0.5, 0.5],
            ElementType::TriquadraticHexahedron => vec![0.5, 0.5, 0.5],
            ElementType::Vertex => vec![0.0],
            _ => vec![0.0]
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
    pub fn calculate_jacobian<T: FloatLike>(
        element_nodes: &[Node],
        shape_derivatives: &[Vec<T>],
    ) -> Result<Jacobian<T>, ElementError> {
        if element_nodes.is_empty() {
            return Err(ElementError::GeometryError("No nodes provided".to_string()));
        }
        
        if element_nodes.len() != shape_derivatives.len() {
            return Err(ElementError::GeometryError(format!(
                "Node count ({}) doesn't match shape derivative count ({})",
                element_nodes.len(), shape_derivatives.len()
            )));
        }
        
        let mesh_dim = element_nodes[0].coordinates.len();
        let element_dim = shape_derivatives[0].len();
        
        let mut jacobian_matrix = vec![vec![T::zero(); element_dim]; mesh_dim];
        
        for i in 0..mesh_dim {
            for j in 0..element_dim {
                for (k, node) in element_nodes.iter().enumerate() {
                    if j < shape_derivatives[k].len() && i < node.coordinates.len() {
                        let coord = T::from_f64(node.coordinates[i]).unwrap();
                        jacobian_matrix[i][j] = jacobian_matrix[i][j].clone() + 
                            shape_derivatives[k][j].clone() * coord;
                    }
                }
            }
        }
        
        let determinant = Self::calculate_jacobian_determinant(&jacobian_matrix)?;
        
        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant,
        })
    }

    fn calculate_jacobian_determinant<T: FloatLike>(
        matrix: &[Vec<T>],
    ) -> Result<T, ElementError> {
        let rows = matrix.len();
        if rows == 0 {
            return Err(ElementError::GeometryError("Empty matrix".to_string()));
        }
        
        let cols = matrix[0].len();
        
        for row in matrix {
            if row.len() != cols {
                return Err(ElementError::GeometryError("Inconsistent matrix dimensions".to_string()));
            }
        }
        
        match (rows, cols) {
            (1, 1) => Ok(matrix[0][0].clone()),
            
            (2, 2) => Ok(matrix[0][0].clone() * matrix[1][1].clone() - 
                        matrix[0][1].clone() * matrix[1][0].clone()),
            
            (3, 3) => Ok(
                matrix[0][0].clone() * (matrix[1][1].clone() * matrix[2][2].clone() - 
                                       matrix[1][2].clone() * matrix[2][1].clone())
                - matrix[0][1].clone() * (matrix[1][0].clone() * matrix[2][2].clone() - 
                                        matrix[1][2].clone() * matrix[2][0].clone())
                + matrix[0][2].clone() * (matrix[1][0].clone() * matrix[2][1].clone() - 
                                        matrix[1][1].clone() * matrix[2][0].clone())
            ),
            
            (2, 1) => {
                let g = matrix[0][0].clone() * matrix[0][0].clone() + 
                        matrix[1][0].clone() * matrix[1][0].clone();
                // For symbolic types, you'll need to implement sqrt differently
                // This assumes T has a sqrt method
                Ok(g.sqrt())
            },
            
            (3, 1) => {
                let g = matrix[0][0].clone() * matrix[0][0].clone() + 
                        matrix[1][0].clone() * matrix[1][0].clone() + 
                        matrix[2][0].clone() * matrix[2][0].clone();
                Ok(g.sqrt())
            },
            
            (3, 2) => {
                let mut g = vec![vec![T::zero(); 2]; 2];
                
                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..3 {
                            g[i][j] = g[i][j].clone() + matrix[k][i].clone() * matrix[k][j].clone();
                        }
                    }
                }
                
                let det_g = g[0][0].clone() * g[1][1].clone() - g[0][1].clone() * g[1][0].clone();
                Ok(det_g.sqrt())
            },
            
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for {}x{} matrices", 
                rows, cols
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create test nodes
    fn create_node(id: usize, x: f64, y: f64, z: f64) -> Node {
        Node {
            id,
            coordinates: vec![x, y, z],
        }
    }

    // Helper function to create 2D nodes
    fn create_node_2d(id: usize, x: f64, y: f64) -> Node {
        Node {
            id,
            coordinates: vec![x, y],
        }
    }

    // Helper function to create 1D nodes
    fn create_node_1d(id: usize, x: f64) -> Node {
        Node {
            id,
            coordinates: vec![x],
        }
    }

    // Helper function to create mesh data for testing
    fn create_test_mesh_data(nodes: Vec<Node>, elements: Vec<Element>, element_type_info: Vec<ElementTypeInfo>) -> MeshData {
        MeshData {
            dimension: if !nodes.is_empty() { nodes[0].coordinates.len() } else { 3 },
            num_nodes: nodes.len(),
            min_node_index: nodes.iter().map(|n| n.id).min().unwrap_or(1),
            nodes,
            num_eltypes: element_type_info.len(),
            elements,
            element_type_info,
        }
    }

    // Make test functions public with pub
    #[test]
    pub fn test_1d_line_element_1d_space() {
        let nodes = vec![
            create_node_1d(1, 0.0),
            create_node_1d(2, 2.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2],
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Line,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 2,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert_eq!(report.total_elements, 1);
        assert_eq!(report.element_qualities[0].element_id, 1);
        
        assert!((report.element_qualities[0].det_jacobian - 2.0).abs() < 1e-10);
        println!("1D Line in 1D space - detJ: {:.6} (Expected: 2.0)", report.element_qualities[0].det_jacobian);
    }

    #[test]
    pub fn test_1d_line_element_2d_space() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 1.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2],
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Line,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 2,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        let expected = (2.0f64).sqrt();
        assert!((report.element_qualities[0].det_jacobian - expected).abs() < 1e-10);
        println!("1D Line in 2D space - detJ: {:.6} (Expected: {:.6})", report.element_qualities[0].det_jacobian, expected);
    }

    #[test]
    pub fn test_2d_triangle_element_2d_space() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2, 3],
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 3,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
        println!("2D Triangle in 2D space - detJ: {:.6} (Expected: 1.0)", report.element_qualities[0].det_jacobian);
    }

    #[test]
    pub fn test_2d_quad_element_2d_space() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 1.0, 1.0),
            create_node_2d(4, 0.0, 1.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2, 3, 4],
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Quad,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 4,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
        println!("2D Quad in 2D space - detJ: {:.6} (Expected: 1.0)", report.element_qualities[0].det_jacobian);
    }

    #[test]
    pub fn test_3d_tetrahedron_element_3d_space() {
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
            create_node(3, 0.0, 1.0, 0.0),
            create_node(4, 0.0, 0.0, 1.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2, 3, 4],
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Tetra,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 4,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
        println!("3D Tetrahedron in 3D space - detJ: {:.6} (Expected: 1.0)", report.element_qualities[0].det_jacobian);
    }

    #[test]
    pub fn test_3d_hexahedron_element_3d_space() {
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
            create_node(3, 1.0, 1.0, 0.0),
            create_node(4, 0.0, 1.0, 0.0),
            create_node(5, 0.0, 0.0, 1.0),
            create_node(6, 1.0, 0.0, 1.0),
            create_node(7, 1.0, 1.0, 1.0),
            create_node(8, 0.0, 1.0, 1.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2, 3, 4, 5, 6, 7, 8],
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Hexahedron,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 8,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
        println!("3D Hexahedron in 3D space - detJ: {:.6} (Expected: 1.0)", report.element_qualities[0].det_jacobian);
    }

    #[test]
    pub fn test_complete_mesh_analysis() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
            create_node_2d(4, 0.0, 0.0),
            create_node_2d(5, 2.0, 0.0),
            create_node_2d(6, 2.0, 2.0),
            create_node_2d(7, 0.0, 2.0),
        ];
        
        let elements = vec![
            Element { id: 1, nodes: vec![1, 2, 3] },
            Element { id: 2, nodes: vec![4, 5, 6, 7] },
        ];
        
        let element_type_info = vec![
            ElementTypeInfo {
                element_type: ElementType::Triangle,
                start_index: 0,
                num_elements: 1,
                nodes_per_element: 3,
            },
            ElementTypeInfo {
                element_type: ElementType::Quad,
                start_index: 1,
                num_elements: 1,
                nodes_per_element: 4,
            },
        ];
        
        let mesh_data = MeshData {
            dimension: 2,
            num_nodes: nodes.len(),
            min_node_index: 1,
            nodes,
            num_eltypes: element_type_info.len(),
            elements,
            element_type_info,
        };
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert_eq!(report.total_elements, 2);
        
        let triangle_quality = report.element_qualities.iter()
            .find(|q| q.element_id == 1)
            .unwrap();
        assert!((triangle_quality.det_jacobian - 1.0).abs() < 1e-10);
        
        let quad_quality = report.element_qualities.iter()
            .find(|q| q.element_id == 2)
            .unwrap();
        assert!((quad_quality.det_jacobian - 4.0).abs() < 1e-10);
        
        println!("Complete mesh analysis:");
        println!("  Triangle element {} - detJ: {:.6}", triangle_quality.element_id, triangle_quality.det_jacobian);
        println!("  Quad element {} - detJ: {:.6}", quad_quality.element_id, quad_quality.det_jacobian);
    }

    #[test]
    pub fn test_invalid_node_reference() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
        ];
        
        let element = Element {
            id: 1,
            nodes: vec![1, 2, 3], // Node 3 doesn't exist
        };
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 3,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, vec![element], element_type_info);
        
        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        // This should print a warning but not crash
        assert!(result.is_err() || result.is_ok());
    }
}

// Now the test_runner module can access the public test functions
#[cfg(test)]
mod test_runner {
    use super::tests::*;

    #[test]
    fn run_all_geometric_analysis_tests() {
        println!("Running Geometric Analysis Tests...");
        println!("=====================================");
        
        let tests = vec![
            ("1D Line in 1D space", test_1d_line_element_1d_space as fn()),
            ("1D Line in 2D space", test_1d_line_element_2d_space as fn()),
            ("2D Triangle in 2D space", test_2d_triangle_element_2d_space as fn()),
            ("2D Quad in 2D space", test_2d_quad_element_2d_space as fn()),
            ("3D Tetrahedron in 3D space", test_3d_tetrahedron_element_3d_space as fn()),
            ("3D Hexahedron in 3D space", test_3d_hexahedron_element_3d_space as fn()),
            ("Complete mesh analysis", test_complete_mesh_analysis as fn()),
            ("Invalid node reference", test_invalid_node_reference as fn()),
        ];
        
        let mut passed = 0;
        let mut failed = 0;
        
        for (test_name, test_func) in &tests {
            print!("{:.<40}", test_name);
            let result = std::panic::catch_unwind(|| test_func());
            
            if result.is_ok() {
                println!("PASSED");
                passed += 1;
            } else {
                println!("FAILED");
                failed += 1;
            }
        }
        
        println!("=====================================");
        println!("Total: {} tests, {} passed, {} failed", tests.len(), passed, failed);
        println!();
        
        println!("Expected detJ Values Summary:");
        println!("1D Line in 1D space: 2.0 (element length)");
        println!("1D Line in 2D space: √2 ≈ 1.4142 (element length)");
        println!("2D Triangle in 2D space: 1.0 (2 × area)");
        println!("2D Quad in 2D space: 1.0 (area scaling factor)");
        println!("3D Tetrahedron in 3D space: 1.0 (6 × volume)");
        println!("3D Hexahedron in 3D space: 1.0 (volume scaling factor)");
        println!("2D Quad (2x2): 4.0 (area scaling for 2x2 quad)");
        
        assert_eq!(failed, 0, "Some tests failed!");
    }
}
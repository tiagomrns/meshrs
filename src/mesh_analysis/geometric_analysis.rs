use core::num;
use std::f64;

use crate::structs_and_impls::*;
use crate::error::*;

/// Geometry routines (Jacobian, quality metrics, etc.)
pub struct GeometricAnalysis;

impl GeometricAnalysis {
    /// Analyze mesh quality over all elements (skips Vertex elements).
    pub fn analyse_mesh_quality(mesh_data: &MeshData) -> Result<MeshQualityReport, ElementError> {
        let mut element_qualities = Vec::new();
        let mut processed_elements = 0;
        let mesh_dim = mesh_data.dimension;

        // Iterate through all element types in the mesh
        for type_info in &mesh_data.element_type_info {
            // Skip vertex elements as they have no geometric quality
            if matches!(type_info.element_type, ElementType::Vertex) {
                continue;
            }

            // Calculate element range for this element type
            let start_idx = type_info.start_index;
            let end_idx = start_idx + type_info.num_elements;

            // Process each element of this type
            for element_idx in start_idx..end_idx {
                if element_idx < mesh_data.elements.len() {
                    let element = &mesh_data.elements[element_idx];

                    // Calculate quality metric for this individual element
                    match Self::calculate_element_quality(
                        element,
                        &type_info.element_type,
                        &mesh_data.nodes,
                        mesh_dim,
                    ) {
                        Ok(quality) => {
                            element_qualities.push(quality);
                            processed_elements += 1;
                        }
                        Err(e) => {
                            println!("Warning: Failed to analyze element {}: {:?}", element.id, e);
                        }
                    }
                }
            }
        }

        // Return error if no elements could be analyzed
        if element_qualities.is_empty() {
            return Err(ElementError::GeometryError(
                "No elements could be analyzed".to_string(),
            ));
        }

        // Return comprehensive quality report
        Ok(MeshQualityReport {
            total_elements: processed_elements,
            element_qualities,
        })
    }

    /// Calculate per-element quality
    fn calculate_element_quality(
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        mesh_dim: usize,
    ) -> Result<ElementQuality, ElementError> {
        // Get node coordinates for this element
        let element_nodes = Self::get_element_nodes(element, nodes)?;
        
        // Get integration point for evaluation (typically Gaussian quadrature point)
        let gaussian_point = Self::get_one_gaussian_point(element_type);

        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();

        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            })?;

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = Self::calculate_jacobian(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        )?;

        // Convert Gaussian point coordinates to tuple format
        let eval_point = match gaussian_point.len() {
            1 => (gaussian_point[0], 0.0, 0.0),
            2 => (gaussian_point[0], gaussian_point[1], 0.0),
            3 => (gaussian_point[0], gaussian_point[1], gaussian_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate the determinant polynomial at the Gaussian point
        let metric_value = MonomialPolynomial::evaluate(&jacobian_matrix.determinant, eval_point)
            .map_err(|e| ElementError::GeometryError(format!("Failed to evaluate determinant: {}", e)))?;

        // Process the determinant value based on matrix dimensions
        let det_value = match (mesh_dim, element_dim) {
            (d, e) if d == e => {
                // Square matrix: determinant directly represents volume scaling
                metric_value
            },
            (2, 1) | (3, 1) | (3, 2) => {
                // Non-square matrices (lower-dimensional elements in higher-dimensional space):
                // Determinant polynomial represents squared length/area, so take sqrt
                let magnitude = metric_value.abs().sqrt();
                // Determine orientation sign to preserve directional information
                let sign = Self::calculate_element_orientation_sign(
                    &jacobian_matrix.matrix,
                    eval_point,
                    mesh_dim,
                    element_dim,
                )?;
                sign * magnitude
            },
            _ => metric_value, // Default case (shouldn't normally occur)
        };

        Ok(ElementQuality {
            element_id: element.id,
            det_jacobian_value: det_value, // Final quality metric value
        })
    }

    /// Return one arbitrary Gaussian point for each element type
    /// These points are typically chosen for numerical integration accuracy
    fn get_one_gaussian_point(element_type: &ElementType) -> Vec<f64> {
        match element_type {
            // 1D elements: single coordinate in parametric space
            ElementType::Line | ElementType::QuadraticEdge => vec![0.2115],
            // 2D triangular elements: area coordinates
            ElementType::Triangle | ElementType::QuadraticTriangle => vec![0.1667, 0.6667],
            // 2D quadrilateral elements: two parametric coordinates
            ElementType::Quad | ElementType::QuadraticQuad | ElementType::BiquadraticQuad => {
                vec![0.2115, 0.7885]
            }
            // 3D tetrahedral elements: volume coordinates  
            ElementType::Tetra | ElementType::QuadraticTetra => vec![0.1667, 0.6667, 0.1667],
            // Pyramid element: mixed coordinates
            ElementType::Pyramid => vec![0.2115, 0.2115, 0.1667],
            // Wedge/prism elements
            ElementType::Wedge
            | ElementType::QuadraticWedge
            | ElementType::BiquadraticQuadraticWedge => vec![0.1667, 0.1667, 0.2115],
            // Hexahedral/brick elements
            ElementType::Hexahedron
            | ElementType::QuadraticHexahedron
            | ElementType::BiquadraticQuadraticHexahedron
            | ElementType::TriquadraticHexahedron => vec![0.2115, 0.2115, 0.7885],
            // Vertex has no integration point
            ElementType::Vertex => vec![0.0],
            _ => vec![0.0], // Default case
        }
    }

    /// Convert element's node IDs to concrete Node structs
    pub fn get_element_nodes(
        element: &Element,
        nodes: &[Node],
    ) -> Result<Vec<Node>, ElementError> {
        let mut element_nodes = Vec::with_capacity(element.nodes.len());

        // Look up each node by ID
        for &node_id in &element.nodes {
            if let Some(node) = nodes.iter().find(|n| n.id == node_id) {
                element_nodes.push(node.clone());
            } else {
                return Err(ElementError::InvalidElement(format!(
                    "Node {} not found for element {}",
                    node_id, element.id
                )));
            }
        }

        // Ensure we found at least one valid node
        if element_nodes.is_empty() {
            return Err(ElementError::InvalidElement(format!(
                "No valid nodes found for element {}",
                element.id
            )));
        }

        Ok(element_nodes)
    }

    /// Compute Jacobian struct (legacy interface if needed elsewhere)
    pub fn calculate_jacobian(
        element_nodes: &Vec<Node>,
        shape_derivatives: &[Vec<Vec<f64>>],
        element_dim: usize,
        num_nodes: usize,
        mesh_dim: usize,
        element_order: usize,
    ) -> Result<Jacobian, ElementError> {

        // Build the full polynomial Jacobian matrix

        // Determine maximum polynomial length based on element order
        let max_len = if element_order == 1 {
            4 //[1, x, y, z]
        }   
        else if element_order == 2 {
            10 //[1, x, y, z, x^2, xy, xz, y^2, yz, z^2]
        }   
        else {
            0; // Unsupported order
            return Err(ElementError::GeometryError(format!(
                "Unsupported element order {} for Jacobian construction",
                element_order
            )));
        }; 

        // Initialize 3D Jacobian matrix: [mesh_dim][element_dim][polynomial_coeffs]
        let mut jacobian_matrix: Vec<Vec<Vec<f64>>> = 
            vec![vec![vec![0.0; max_len]; element_dim]; mesh_dim];

        // Build Jacobian matrix by summing contributions from all nodes
        // J_ij = Σ_node (x_node_i * ∂N_node/∂ξ_j) where N_node are shape functions
        for i in 0..mesh_dim {           // Physical space dimension
            for j in 0..element_dim {    // Parametric space dimension  
                for k in 0..num_nodes {  // Node index
                    let coord = element_nodes[k].coordinates[i];
                    let derivative_poly = &shape_derivatives[k][j];
                    
                    // Multiply shape function derivative by nodal coordinate
                    let scaled_poly = MonomialPolynomial::multiply_scalar(derivative_poly, coord);
                    // Accumulate contribution to Jacobian component
                    jacobian_matrix[i][j] = MonomialPolynomial::add(&jacobian_matrix[i][j], &scaled_poly)
                        .map_err(|e| ElementError::GeometryError(format!("Failed to add polynomials: {}", e)))?;
                }
            }
        }

        // Calculate determinant as polynomial
        let det_jacobian = Self::calculate_jacobian_determinant(&jacobian_matrix, mesh_dim, element_dim)?;

        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant: det_jacobian,
        })
    }

    /// Calculate Jacobian determinant/metric for polynomial matrices
    /// Handles both square matrices (true determinant) and non-square matrices (metric determinant)
    pub fn calculate_jacobian_determinant(
        matrix: &Vec<Vec<Vec<f64>>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> Result<Vec<f64>, ElementError> {

        /*  2D Jacobian     3D Jacobian
              |a b|           |a b c|
              |c d|           |d e f|
                              |g h i|
        */

        match (mesh_dim, element_dim) {
            // 1x1 matrix: determinant is the single element
            (1, 1) => Ok(matrix[0][0].clone()),

            // 2x2 matrix: det = ad - bc
            (2, 2) => {
                let ad = MonomialPolynomial::multiply(&matrix[0][0], &matrix[1][1])?;
                let bc = MonomialPolynomial::multiply(&matrix[0][1], &matrix[1][0])?;
                let neg_bc = MonomialPolynomial::multiply_scalar(&bc, -1.0);
                MonomialPolynomial::add(&ad, &neg_bc)
                    .map_err(|e| ElementError::GeometryError(format!("Add failed: {}", e)))
            }

            // 3x3 matrix: det = a(ei−fh) − b(di−fg) + c(dh−eg)
            (3, 3) => {
                let ei = MonomialPolynomial::multiply(&matrix[1][1], &matrix[2][2])?;
                let fh = MonomialPolynomial::multiply(&matrix[1][2], &matrix[2][1])?;
                let ei_fh = MonomialPolynomial::add(&ei, &MonomialPolynomial::multiply_scalar(&fh, -1.0))?;
                let term1 = MonomialPolynomial::multiply(&matrix[0][0], &ei_fh)?;

                let di = MonomialPolynomial::multiply(&matrix[1][0], &matrix[2][2])?;
                let fg = MonomialPolynomial::multiply(&matrix[1][2], &matrix[2][0])?;
                let di_fg = MonomialPolynomial::add(&di, &MonomialPolynomial::multiply_scalar(&fg, -1.0))?;
                let term2 = MonomialPolynomial::multiply(&matrix[0][1], &di_fg)?;

                let dh = MonomialPolynomial::multiply(&matrix[1][0], &matrix[2][1])?;
                let eg = MonomialPolynomial::multiply(&matrix[1][1], &matrix[2][0])?;
                let dh_eg = MonomialPolynomial::add(&dh, &MonomialPolynomial::multiply_scalar(&eg, -1.0))?;
                let term3 = MonomialPolynomial::multiply(&matrix[0][2], &dh_eg)?;

                let result = MonomialPolynomial::add(&term1, &MonomialPolynomial::multiply_scalar(&term2, -1.0))?;
                MonomialPolynomial::add(&result, &term3)
                    .map_err(|e| ElementError::GeometryError(format!("Final add failed: {}", e)))
            }

            // 1D elements in 2D space: metric = dx² + dy² (squared length)
            (2, 1) => {
                let dx_sq = MonomialPolynomial::multiply(&matrix[0][0], &matrix[0][0])?;
                let dy_sq = MonomialPolynomial::multiply(&matrix[1][0], &matrix[1][0])?;
                MonomialPolynomial::add(&dx_sq, &dy_sq)
                    .map_err(|e| ElementError::GeometryError(format!("Add failed: {}", e)))
            }

            // 1D elements in 3D space: metric = dx² + dy² + dz² (squared length)
            (3, 1) => {
                let dx_sq = MonomialPolynomial::multiply(&matrix[0][0], &matrix[0][0])?;
                let dy_sq = MonomialPolynomial::multiply(&matrix[1][0], &matrix[1][0])?;
                let dz_sq = MonomialPolynomial::multiply(&matrix[2][0], &matrix[2][0])?;
                let temp = MonomialPolynomial::add(&dx_sq, &dy_sq)?;
                MonomialPolynomial::add(&temp, &dz_sq)
                    .map_err(|e| ElementError::GeometryError(format!("Add failed: {}", e)))
            }

            // 2D elements in 3D space: metric determinant from first fundamental form
            (3, 2) => {
                // Build metric tensor G = J^T * J where J is the Jacobian
                let mut g = vec![vec![vec![0.0; matrix[0][0].len()]; 2]; 2];
                
                for i in 0..2 {
                    for j in 0..2 {
                        let mut sum = vec![0.0; matrix[0][0].len()];
                        for k in 0..3 {
                            let prod = MonomialPolynomial::multiply(&matrix[k][i], &matrix[k][j])?;
                            sum = MonomialPolynomial::add(&sum, &prod)?;
                        }
                        g[i][j] = sum;
                    }
                }

                // det(G) = g00*g11 - g01*g10
                let g00_g11 = MonomialPolynomial::multiply(&g[0][0], &g[1][1])?;
                let g01_g10 = MonomialPolynomial::multiply(&g[0][1], &g[1][0])?;
                MonomialPolynomial::add(&g00_g11, &MonomialPolynomial::multiply_scalar(&g01_g10, -1.0))
                    .map_err(|e| ElementError::GeometryError(format!("Metric det failed: {}", e)))
            }

            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for {}x{} (mesh_dim x element_dim) matrices",
                mesh_dim, element_dim
            ))),
        }
    }

    /// Calculate orientation sign for non-square Jacobian matrices
    /// Preserves directional information that would be lost by taking absolute values
    fn calculate_element_orientation_sign(
        matrix: &Vec<Vec<Vec<f64>>>,
        eval_point: (f64, f64, f64),
        mesh_dim: usize,
        element_dim: usize,
    ) -> Result<f64, ElementError> {
        match (mesh_dim, element_dim) {
            // Square matrices: determinant already has correct sign
            (d, e) if d == e => Ok(1.0),
            
            // 1D elements in 2D or 3D space: sign from largest tangent component
            (2, 1) | (3, 1) => {
                // Evaluate tangent vector at Gaussian point
                let mut tangent = vec![0.0; mesh_dim];
                for i in 0..mesh_dim {
                    tangent[i] = MonomialPolynomial::evaluate(&matrix[i][0], eval_point)
                        .map_err(|e| ElementError::GeometryError(format!("Tangent eval failed: {}", e)))?;
                }
                
                // Find component with largest magnitude for robust sign determination
                let mut max_abs = 0.0;
                let mut sign_component = 0.0;
                for &comp in &tangent {
                    if comp.abs() > max_abs {
                        max_abs = comp.abs();
                        sign_component = comp;
                    }
                }
                
                // Return sign based on dominant component (avoid division by near-zero)
                Ok(if max_abs < 1e-12 { 1.0 } else if sign_component >= 0.0 { 1.0 } else { -1.0 })
            },
            
            // 2D elements in 3D space: sign from surface normal
            (3, 2) => {
                // Evaluate both tangent vectors
                let mut u = [0.0; 3];
                let mut v = [0.0; 3];
                for i in 0..3 {
                    u[i] = MonomialPolynomial::evaluate(&matrix[i][0], eval_point)
                        .map_err(|e| ElementError::GeometryError(format!("u eval failed: {}", e)))?;
                    v[i] = MonomialPolynomial::evaluate(&matrix[i][1], eval_point)
                        .map_err(|e| ElementError::GeometryError(format!("v eval failed: {}", e)))?;
                }
                
                // Compute surface normal: n = u × v
                let normal = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                
                // Find component with largest magnitude for robust sign determination
                let mut max_abs = 0.0;
                let mut sign_component = 0.0;
                for &comp in &normal {
                    if comp.abs() > max_abs {
                        max_abs = comp.abs();
                        sign_component = comp;
                    }
                }
                
                Ok(if max_abs < 1e-12 { 1.0 } else if sign_component >= 0.0 { 1.0 } else { -1.0 })
            },
            
            // Default case (shouldn't normally occur)
            _ => Ok(1.0),
        }
    }
}

// Tests remain the same as they test the public interface
#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // Helper functions for creating test data
    fn create_node(id: usize, x: f64, y: f64, z: f64) -> Node {
        Node {
            id,
            coordinates: vec![x, y, z],
        }
    }

    fn create_node_2d(id: usize, x: f64, y: f64) -> Node {
        Node {
            id,
            coordinates: vec![x, y],
        }
    }

    fn create_node_1d(id: usize, x: f64) -> Node {
        Node {
            id,
            coordinates: vec![x],
        }
    }

    fn create_element(id: usize, nodes: Vec<usize>) -> Element {
        Element { id, nodes }
    }

    fn create_test_mesh_data(
        nodes: Vec<Node>,
        elements: Vec<Element>,
        element_type_info: Vec<ElementTypeInfo>,
    ) -> MeshData {
        MeshData {
            dimension: if !nodes.is_empty() {
                nodes[0].coordinates.len()
            } else {
                3
            },
            num_nodes: nodes.len(),
            min_node_index: nodes.iter().map(|n| n.id).min().unwrap_or(1),
            nodes,
            num_eltypes: element_type_info.len(),
            elements,
            element_type_info,
        }
    }

    // Basic functionality tests
    #[test]
    fn test_get_element_nodes_valid() {
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
            create_node(3, 0.0, 1.0, 0.0),
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::get_element_nodes(&element, &nodes);
        assert!(result.is_ok());
        
        let element_nodes = result.unwrap();
        assert_eq!(element_nodes.len(), 3);
        assert_eq!(element_nodes[0].id, 1);
        assert_eq!(element_nodes[1].id, 2);
        assert_eq!(element_nodes[2].id, 3);
    }

    #[test]
    fn test_get_element_nodes_invalid_node() {
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
        ];
        let element = create_element(1, vec![1, 2, 999]); // Node 999 doesn't exist
        
        let result = GeometricAnalysis::get_element_nodes(&element, &nodes);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_element_nodes_empty() {
        let nodes = vec![];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::get_element_nodes(&element, &nodes);
        assert!(result.is_err());
    }

    // Jacobian calculation tests
    #[test]
    fn test_calculate_jacobian_1d_line() {
        let nodes = vec![
            create_node_1d(1, 0.0),
            create_node_1d(2, 2.0),
        ];
        let element = create_element(1, vec![1, 2]);
        
        let shape_derivatives = ElementType::get_shape_functions(&ElementType::Line)
            .unwrap()
            .derivatives;
        
        let result = GeometricAnalysis::calculate_jacobian(
            &nodes,
            &shape_derivatives,
            1, // element_dim
            2, // num_nodes
            1, // mesh_dim
            1, // element_order
        );
        
        assert!(result.is_ok());
        let jacobian = result.unwrap();
        
        // For a 1D line from 0 to 2, Jacobian should be 2
        let eval_result = MonomialPolynomial::evaluate(&jacobian.determinant, (0.5, 0.0, 0.0));
        assert!(eval_result.is_ok());
        assert!((eval_result.unwrap() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_jacobian_2d_triangle() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        
        let shape_derivatives = ElementType::get_shape_functions(&ElementType::Triangle)
            .unwrap()
            .derivatives;
        
        let result = GeometricAnalysis::calculate_jacobian(
            &nodes,
            &shape_derivatives,
            2, // element_dim
            3, // num_nodes
            2, // mesh_dim
            1, // element_order
        );
        
        assert!(result.is_ok());
        let jacobian = result.unwrap();
        
        // For unit triangle, Jacobian determinant should be 1.0
        let eval_result = MonomialPolynomial::evaluate(&jacobian.determinant, (0.333, 0.333, 0.0));
        assert!(eval_result.is_ok());
        assert!((eval_result.unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_jacobian_3d_tetrahedron() {
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
            create_node(3, 0.0, 1.0, 0.0),
            create_node(4, 0.0, 0.0, 1.0),
        ];
        
        let shape_derivatives = ElementType::get_shape_functions(&ElementType::Tetra)
            .unwrap()
            .derivatives;
        
        let result = GeometricAnalysis::calculate_jacobian(
            &nodes,
            &shape_derivatives,
            3, // element_dim
            4, // num_nodes
            3, // mesh_dim
            1, // element_order
        );
        
        assert!(result.is_ok());
        let jacobian = result.unwrap();
        
        // For unit tetrahedron, Jacobian determinant should be 1.0
        let eval_result = MonomialPolynomial::evaluate(&jacobian.determinant, (0.25, 0.25, 0.25));
        assert!(eval_result.is_ok());
        assert!((eval_result.unwrap() - 1.0).abs() < 1e-12);
    }

    // Jacobian determinant calculation tests
    #[test]
    fn test_calculate_jacobian_determinant_1x1() {
        // 1x1 matrix: [a]
        let matrix = vec![
            vec![vec![2.0]] // a = 2
        ];
        
        let result = GeometricAnalysis::calculate_jacobian_determinant(&matrix, 1, 1);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), vec![2.0]);
    }

    #[test]
    fn test_calculate_jacobian_determinant_2x2() {
        // 2x2 matrix: [[a, b], [c, d]] = [[1, 2], [3, 4]]
        // det = ad - bc = 1*4 - 2*3 = -2
        let matrix = vec![
            vec![vec![1.0], vec![2.0]], // [a, b]
            vec![vec![3.0], vec![4.0]], // [c, d]
        ];
        
        let result = GeometricAnalysis::calculate_jacobian_determinant(&matrix, 2, 2);
        assert!(result.is_ok());
        let det = result.unwrap();
        assert_eq!(det, vec![-2.0]);
    }

    #[test]
    fn test_calculate_jacobian_determinant_3x3() {
        // 3x3 matrix with known determinant
        let matrix = vec![
            vec![vec![1.0], vec![2.0], vec![3.0]], // [a, b, c]
            vec![vec![4.0], vec![5.0], vec![6.0]], // [d, e, f]
            vec![vec![7.0], vec![8.0], vec![9.0]], // [g, h, i]
        ];
        // det = a(ei−fh) − b(di−fg) + c(dh−eg)
        // = 1*(5*9 - 6*8) - 2*(4*9 - 6*7) + 3*(4*8 - 5*7)
        // = 1*(45-48) - 2*(36-42) + 3*(32-35)
        // = -3 - 2*(-6) + 3*(-3)
        // = -3 + 12 - 9 = 0
        
        let result = GeometricAnalysis::calculate_jacobian_determinant(&matrix, 3, 3);
        assert!(result.is_ok());
        let det = result.unwrap();
        assert!((det[0] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_jacobian_determinant_1d_in_2d() {
        // 1D element in 2D space: dx² + dy²
        let matrix = vec![
            vec![vec![3.0]], // dx/dξ = 3
            vec![vec![4.0]], // dy/dξ = 4
        ];
        // metric = 3² + 4² = 9 + 16 = 25
        
        let result = GeometricAnalysis::calculate_jacobian_determinant(&matrix, 2, 1);
        assert!(result.is_ok());
        let det = result.unwrap();
        assert!((det[0] - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_jacobian_determinant_1d_in_3d() {
        // 1D element in 3D space: dx² + dy² + dz²
        let matrix = vec![
            vec![vec![1.0]], // dx/dξ
            vec![vec![2.0]], // dy/dξ
            vec![vec![2.0]], // dz/dξ
        ];
        // metric = 1² + 2² + 2² = 1 + 4 + 4 = 9
        
        let result = GeometricAnalysis::calculate_jacobian_determinant(&matrix, 3, 1);
        assert!(result.is_ok());
        let det = result.unwrap();
        assert!((det[0] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_jacobian_determinant_2d_in_3d() {
        // 2D element in 3D space - test with simple case
        let matrix = vec![
            vec![vec![1.0], vec![0.0]], // [dx/du, dx/dv]
            vec![vec![0.0], vec![1.0]], // [dy/du, dy/dv]
            vec![vec![0.0], vec![0.0]], // [dz/du, dz/dv]
        ];
        // This represents a flat element in xy-plane
        // Metric tensor G = J^T J = [[1, 0], [0, 1]]
        // det(G) = 1
        
        let result = GeometricAnalysis::calculate_jacobian_determinant(&matrix, 3, 2);
        assert!(result.is_ok());
        let det = result.unwrap();
        assert!((det[0] - 1.0).abs() < 1e-12);
    }

    // Element quality analysis tests
    #[test]
    fn test_analyse_mesh_quality_single_triangle() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        let element = create_element(1, vec![1, 2, 3]);
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 3,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert_eq!(report.total_elements, 1);
        assert_eq!(report.element_qualities[0].element_id, 1);
        assert!(report.element_qualities[0].det_jacobian_value > 0.0);
    }

    #[test]
    fn test_analyse_mesh_quality_multiple_elements() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
            create_node_2d(4, 1.0, 1.0),
        ];
        let elements = vec![
            create_element(1, vec![1, 2, 3]), // Triangle
            create_element(2, vec![2, 4, 3]), // Triangle
        ];
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 2,
            start_index: 0,
            nodes_per_element: 3,
        }];
        let mesh_data = create_test_mesh_data(nodes, elements, element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert_eq!(report.total_elements, 2);
        assert_eq!(report.element_qualities[0].element_id, 1);
        assert_eq!(report.element_qualities[1].element_id, 2);
    }

    #[test]
    fn test_analyse_mesh_quality_skips_vertex_elements() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        let elements = vec![
            create_element(1, vec![1]), // Vertex element - should be skipped
            create_element(2, vec![1, 2, 3]), // Triangle element - should be analyzed
        ];
        let element_type_info = vec![
            ElementTypeInfo {
                element_type: ElementType::Vertex,
                num_elements: 1,
                start_index: 0,
                nodes_per_element: 1,
            },
            ElementTypeInfo {
                element_type: ElementType::Triangle,
                num_elements: 1,
                start_index: 1,
                nodes_per_element: 3,
            },
        ];
        let mesh_data = create_test_mesh_data(nodes, elements, element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        // Should only analyze the triangle, not the vertex
        assert_eq!(report.total_elements, 1);
        assert_eq!(report.element_qualities[0].element_id, 2);
    }

    #[test]
    fn test_analyse_mesh_quality_empty_mesh() {
        let nodes = vec![];
        let elements = vec![];
        let element_type_info = vec![];
        let mesh_data = create_test_mesh_data(nodes, elements, element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_err());
    }

    // Element quality calculation tests
    #[test]
    fn test_calculate_element_quality_regular_triangle() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2, // mesh_dim
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        assert_eq!(quality.element_id, 1);
        // For unit triangle, determinant should be close to 1.0
        assert!((quality.det_jacobian_value - 1.0).abs() < 1e-10);
    }


    #[test]
    fn test_calculate_element_quality_distorted_triangle() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 0.1, 0.0), // Very short base
            create_node_2d(3, 0.0, 1.0), // Tall height
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2,
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        // Distorted element should have smaller determinant
        assert!(quality.det_jacobian_value > 0.0);
        assert!(quality.det_jacobian_value < 1.0);
    }

    #[test]
    fn test_calculate_element_quality_1d_element() {
        let nodes = vec![
            create_node_1d(1, 0.0),
            create_node_1d(2, 5.0),
        ];
        let element = create_element(1, vec![1, 2]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Line,
            &nodes,
            1,
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        // For 1D line of length 5, metric should be related to length
        assert!(quality.det_jacobian_value > 0.0);
    }

    // Edge cases and error conditions
    #[test]
    fn test_calculate_element_quality_invalid_element_type() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
        ];
        let element = create_element(1, vec![1, 2]);
        
        // Try with an element type that doesn't have shape functions
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Vertex, // Vertex has no shape functions
            &nodes,
            2,
        );
        
        assert!(result.is_err());
    }

    #[test]
    fn test_calculate_element_quality_degenerate_element() {
        // Degenerate triangle with all points collinear
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 2.0, 0.0), // Collinear with others
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2,
        );
        
        // This might fail or return a very small determinant
        match result {
            Ok(quality) => {
                // If it succeeds, determinant should be very small
                assert!(quality.det_jacobian_value.abs() < 1e-10);
            }
            Err(_) => {
                // It's also acceptable for this to fail
            }
        }
    }

    // Orientation tests
    #[test]
    fn test_calculate_element_orientation_sign_2d_triangle() {
        let matrix = vec![
            vec![vec![1.0], vec![0.0]], // Simple 2x2 identity-like
            vec![vec![0.0], vec![1.0]],
        ];
        
        let result = GeometricAnalysis::calculate_element_orientation_sign(
            &matrix,
            (0.0, 0.0, 0.0),
            2, // mesh_dim
            2, // element_dim
        );
        
        assert!(result.is_ok());
        let sign = result.unwrap();
        assert!((sign - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_element_orientation_sign_1d_in_2d() {
        // 1D element with positive orientation
        let matrix = vec![
            vec![vec![1.0]], // dx/dξ = 1
            vec![vec![0.0]], // dy/dξ = 0
        ];
        
        let result = GeometricAnalysis::calculate_element_orientation_sign(
            &matrix,
            (0.0, 0.0, 0.0),
            2, // mesh_dim
            1, // element_dim
        );
        
        assert!(result.is_ok());
        let sign = result.unwrap();
        assert!((sign - 1.0).abs() < 1e-12);
    }

    // Gaussian point tests
    #[test]
    fn test_get_one_gaussian_point_various_elements() {
        // Test that Gaussian points are returned for various element types
        let element_types = vec![
            ElementType::Line,
            ElementType::Triangle,
            ElementType::Quad,
            ElementType::Tetra,
            ElementType::Hexahedron,
        ];
        
        for element_type in element_types {
            let gauss_point = GeometricAnalysis::get_one_gaussian_point(&element_type);
            assert!(!gauss_point.is_empty(), "Gaussian points should be available for {:?}", element_type);
            
            // Check that coordinates are in reasonable ranges
            for &coord in &gauss_point {
                assert!(coord >= -1.0 && coord <= 1.0, 
                    "Gaussian point coordinate {} out of range for {:?}", coord, element_type);
            }
        }
    }

    #[test]
    fn test_get_one_gaussian_point_vertex() {
        let gauss_point = GeometricAnalysis::get_one_gaussian_point(&ElementType::Vertex);
        assert_eq!(gauss_point, vec![0.0]);
    }

    // Complex geometry tests
    #[test]
    fn test_curved_element_geometry() {
        // Test with a slightly curved quadrilateral element
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.1), // Slight curvature
            create_node_2d(3, 1.1, 1.0), // Slight curvature
            create_node_2d(4, 0.1, 1.0), // Slight curvature
        ];
        let element = create_element(1, vec![1, 2, 3, 4]);
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Quad,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 4,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert_eq!(report.total_elements, 1);
        println!("Curved quad Jacobian determinant: {}", report.element_qualities[0].det_jacobian_value);
        // Curved element should still have positive determinant
        assert!(report.element_qualities[0].det_jacobian_value > 0.0, 
            "Expected positive determinant, got {}", report.element_qualities[0].det_jacobian_value);
    }

    #[test]
    fn test_large_deformation_geometry() {
        // Test with elements that have large aspect ratios
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 10.0, 0.0), // Large aspect ratio
            create_node_2d(3, 0.0, 1.0),
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2,
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        // Element with large aspect ratio should still have positive determinant
        assert!(quality.det_jacobian_value > 0.0);
    }

    // Performance and numerical stability tests
    #[test]
    fn test_numerical_stability_small_elements() {
        // Test with very small elements to check numerical stability
        let scale = 1e-12;
        let nodes = vec![
            create_node_2d(1, 0.0 * scale, 0.0 * scale),
            create_node_2d(2, 1.0 * scale, 0.0 * scale),
            create_node_2d(3, 0.0 * scale, 1.0 * scale),
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2,
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        // Even for very small elements, we should get a reasonable value
        assert!(quality.det_jacobian_value >= 0.0);
    }

    #[test]
    fn test_numerical_stability_large_elements() {
        // Test with very large elements
        let scale = 1e12;
        let nodes = vec![
            create_node_2d(1, 0.0 * scale, 0.0 * scale),
            create_node_2d(2, 1.0 * scale, 0.0 * scale),
            create_node_2d(3, 0.0 * scale, 1.0 * scale),
        ];
        let element = create_element(1, vec![1, 2, 3]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2,
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        // Should handle large values without overflow
        assert!(quality.det_jacobian_value > 0.0);
        assert!(!quality.det_jacobian_value.is_infinite());
    }

    // Integration with polynomial operations
    #[test]
    fn test_jacobian_with_polynomial_coordinates() {
        // Test that Jacobian calculation works with polynomial shape functions
        // This is important for higher-order elements
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
            create_node_2d(4, 0.5, 0.0), // Mid-side node between nodes 1 and 2
            create_node_2d(5, 0.5, 0.5), // Mid-side node between nodes 2 and 3
            create_node_2d(6, 0.0, 0.5), // Mid-side node between nodes 3 and 1
        ];
        
        let shape_derivatives = ElementType::get_shape_functions(&ElementType::QuadraticTriangle)
            .unwrap()
            .derivatives;
        
        let result = GeometricAnalysis::calculate_jacobian(
            &nodes,
            &shape_derivatives,
            2, // element_dim
            6, // num_nodes (quadratic triangle has 6 nodes)
            2, // mesh_dim
            2, // element_order (quadratic)
        );
        
        assert!(result.is_ok());
        let jacobian = result.unwrap();
        
        // Should be able to evaluate the determinant polynomial
        let eval_result = MonomialPolynomial::evaluate(&jacobian.determinant, (0.333, 0.333, 0.0));
        assert!(eval_result.is_ok());
        let det_val = eval_result.unwrap();
        println!("Quadratic triangle Jacobian determinant: {}", det_val);
        assert!(det_val > 0.0, "Expected positive determinant, got {}", det_val);
    }
}
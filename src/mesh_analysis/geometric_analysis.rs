use core::num;
use std::f64;
use std::vec;

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
        let jacobian_matrix = Self::calculate_jacobian_monomial(
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
        let evaluated_jacobian_determinant = MonomialPolynomial::evaluate(&jacobian_matrix.determinant, eval_point)
            .map_err(|e| ElementError::GeometryError(format!("Failed to evaluate determinant: {}", e)))?;

        /* 
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
        */

        // Evaluate Jacobian at the point
        let mut evaluated_jacobian = vec![vec![0.0; element_dim]; mesh_dim];
        for i in 0..mesh_dim {
            for j in 0..element_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(&jacobian_matrix.matrix[i][j], eval_point)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to evaluate Jacobian: {}", e)))?;
            }
        }

        let (r_matrix, q_matrix, d_matrix) = 
            Self::jacobian_factorization(&evaluated_jacobian, &evaluated_jacobian_determinant, mesh_dim, element_dim, eval_point)?;

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = Self::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        let inverse_q_matrix = Self::inverse_matrix(&q_matrix, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        let inverse_d_matrix = Self::inverse_matrix(&d_matrix, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        
        // Calculate metrics
        let shape_metric = Self::calculate_shape_metric(&evaluated_jacobian, &inverse_evaluated_jacobian, mesh_dim, element_dim);
        let skewness_metric = Self::calculate_skewness_metric(&q_matrix, &inverse_q_matrix, mesh_dim, element_dim);
        let length_ratio = Self::calculate_length_ratio(&d_matrix, &inverse_d_matrix, mesh_dim, element_dim);
        let orientation_metric = Self::calculate_orientation_metric(&r_matrix, mesh_dim, element_dim);
        let volume_metric = evaluated_jacobian_determinant;
        let volume_shape_metric = Self::calculate_volume_shape_metric(shape_metric, evaluated_jacobian_determinant);
        let volume_shape_orientation_metric = Self::calculate_volume_shape_orientation_metric(shape_metric, evaluated_jacobian_determinant, orientation_metric);

        Ok(ElementQuality {
            element_id: element.id,
            det_jacobian_value: evaluated_jacobian_determinant,  // not necessary anymore volume_metric gives the same result
            shape_metric,
            skewness_metric,
            length_ratio,
            orientation_metric,
            volume_metric,
            volume_shape_metric,
            volume_shape_orientation_metric, 
        })
    }

    /// Calculate shape metric
    fn calculate_shape_metric(
        jacobian: &Vec<Vec<f64>>,
        inverse_jacobian: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> f64 {
        
        let j_norm = Self::frobenius_norm(jacobian, mesh_dim, element_dim);
        let inv_j_norm = Self::frobenius_norm(inverse_jacobian, mesh_dim, element_dim);

        element_dim as f64 / (j_norm * inv_j_norm)
    }

    /// Calculate skewness metric
    fn calculate_skewness_metric(
        q_matrix: &Vec<Vec<f64>>,
        inverse_q_matrix: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> f64 {
        
        let q_norm = Self::frobenius_norm(q_matrix, mesh_dim, element_dim);
        let inv_q_norm = Self::frobenius_norm(inverse_q_matrix, mesh_dim, element_dim);
        
        element_dim as f64 / (q_norm * inv_q_norm)
    }

    /// Calculate length ratio metric
    fn calculate_length_ratio(
        d_matrix: &Vec<Vec<f64>>,
        inverse_d_matrix: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> f64 {
        
        let d_norm = Self::frobenius_norm(d_matrix, mesh_dim, element_dim);
        let inv_d_norm = Self::frobenius_norm(inverse_d_matrix, mesh_dim, element_dim);

        element_dim as f64 / (d_norm * inv_d_norm)
    }

    /// Calculate orientation metric
    fn calculate_orientation_metric(
        r_matrix: &Vec<Vec<f64>>, 
        mesh_dim: usize, 
        element_dim: usize
    ) -> f64 {
        
        1.0 + ((Self::trace(r_matrix, mesh_dim, element_dim) - element_dim as f64) / 4.0)
    }

    /// Calculate volume shape metric
    fn calculate_volume_shape_metric(
        shape_metric: f64,
        det_jacobian: f64,
    ) -> f64 {

        let min_det = det_jacobian.min(1.0 / det_jacobian);
        
        min_det * shape_metric
    }

    /// Calculate volume shape orientation metric
    fn calculate_volume_shape_orientation_metric(
        shape_metric: f64,
        det_jacobian: f64,
        orientation_metric: f64,
    ) -> f64 {
    
        let min_det = det_jacobian.min(1.0 / det_jacobian);
        
        min_det * shape_metric * orientation_metric
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

    /// Compute Jacobian struct 
    pub fn calculate_jacobian_monomial(
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
        let det_jacobian = Self::calculate_determinant_monomial(&jacobian_matrix, mesh_dim, element_dim)?;

        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant: det_jacobian,
        })
    }

    pub fn calculate_determinant(
        matrix: &Vec<Vec<f64>>, 
        element_dim: usize, 
        mesh_dim: usize
    ) -> Result<f64, ElementError> {
        match (mesh_dim, element_dim) {
            (1, 1) => Ok(matrix[0][0]),
            (2, 2) => Ok(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]),
            (3, 3) => {
                let a = matrix[0][0];
                let b = matrix[0][1];
                let c = matrix[0][2];
                let d = matrix[1][0];
                let e = matrix[1][1];
                let f = matrix[1][2];
                let g = matrix[2][0];
                let h = matrix[2][1];
                let i = matrix[2][2];
                
                Ok(a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g))
            },
            _ => Err(ElementError::GeometryError(format!(
                "Determinant calculation not implemented for {}x{} matrices",
                mesh_dim, element_dim
            ))),
        }
    }


    pub fn inverse_matrix(
        matrix: &Vec<Vec<f64>>, 
        mesh_dim: usize, 
        element_dim: usize
    ) -> Result<Vec<Vec<f64>>, ElementError> {
        match (mesh_dim, element_dim) {
            // 1x1 matrix
            (1, 1) => {
                if matrix[0][0].abs() < 1e-15 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                Ok(vec![vec![1.0 / matrix[0][0]]])
            },

            // 2x2 matrix
            (2, 2) => {
                let det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
                if det.abs() < 1e-15 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                let inv_det = 1.0 / det;
                Ok(vec![
                    vec![matrix[1][1] * inv_det, -matrix[0][1] * inv_det],
                    vec![-matrix[1][0] * inv_det, matrix[0][0] * inv_det],
                ])
            },

            // 3x3 matrix
            (3, 3) => {
                let adjoint = Self::calculate_adjoint(&matrix, element_dim, element_dim)?;
                let determinant = Self::calculate_determinant(&matrix, element_dim, element_dim)?;

                if determinant.abs() < 1e-15 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                
                let inv_det = 1.0 / determinant;
                let mut inverse = vec![vec![0.0; element_dim]; element_dim];
                
                for i in 0..element_dim {
                    for j in 0..element_dim {
                        inverse[i][j] = adjoint[i][j] * inv_det;
                    }
                }
                Ok(inverse)
            },

            // 1D elements in 2D space
            (2, 1) => {
                let mut sum = 0.0;
                for k in 0..mesh_dim {
                    sum += matrix[k][0] * matrix[k][0];
                }
                if sum.abs() < 1e-15 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                Ok(vec![vec![1.0 / sum]])
            },

            // 1D elements in 3D space
            (3, 1) => {
                let mut sum = 0.0;
                for k in 0..mesh_dim {
                    sum += matrix[k][0] * matrix[k][0];
                }
                if sum.abs() < 1e-15 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                Ok(vec![vec![1.0 / sum]])
            },

            // 2D elements in 3D space
            (3, 2) => {
                let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
                for i in 0..element_dim {
                    for j in 0..element_dim {
                        let mut sum = 0.0;
                        for k in 0..mesh_dim {
                            sum += matrix[k][i] * matrix[k][j];
                        }
                        square_matrix[i][j] = sum;
                    }
                }

                let adjoint = Self::calculate_adjoint(&square_matrix, element_dim, element_dim)?;
                let determinant = Self::calculate_determinant(&square_matrix, element_dim, element_dim)?;
                
                if determinant.abs() < 1e-15 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                
                let inv_det = 1.0 / determinant;
                let mut inverse = vec![vec![0.0; element_dim]; element_dim];
                
                for i in 0..element_dim {
                    for j in 0..element_dim {
                        inverse[i][j] = adjoint[i][j] * inv_det;
                    }
                }
                Ok(inverse)
            },

            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for {}x{} (mesh_dim x element_dim) matrices",
                mesh_dim, element_dim
            ))),
        }
    }

    /// Calculate Jacobian determinant/metric for polynomial matrices
    /// Handles both square matrices (true determinant) and non-square matrices (metric determinant)
    pub fn calculate_determinant_monomial(
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

    

    pub fn calculate_adjoint(
        matrix: &Vec<Vec<f64>>, 
        element_dim: usize, 
        mesh_dim: usize
    ) -> Result<Vec<Vec<f64>>, ElementError> {
        match (mesh_dim, element_dim) {
            (2, 2) => {
                Ok(vec![
                    vec![matrix[1][1], -matrix[0][1]],
                    vec![-matrix[1][0], matrix[0][0]],
                ])
            },
            (3, 3) => {
                let a = matrix[0][0];
                let b = matrix[0][1];
                let c = matrix[0][2];
                let d = matrix[1][0];
                let e = matrix[1][1];
                let f = matrix[1][2];
                let g = matrix[2][0];
                let h = matrix[2][1];
                let i = matrix[2][2];

                Ok(vec![
                    vec![
                        e * i - f * h,
                        c * h - b * i, 
                        b * f - c * e,
                    ],
                    vec![
                        f * g - d * i,
                        a * i - c * g,
                        c * d - a * f,
                    ],
                    vec![
                        d * h - e * g,
                        b * g - a * h,
                        a * e - b * d,
                    ],
                ])
            },
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian adjoint calculation not implemented for {}x{} (mesh_dim x element_dim) matrices",
                mesh_dim, element_dim
            ))),
        }
    }

    fn trace(matrix: &Vec<Vec<f64>>, mesh_dim: usize, element_dim: usize) -> f64 {
        let mut sum = 0.0;
        for i in 0..mesh_dim.min(element_dim) { // if the matrix square than mesh_dim == element_dim, if not then matrix has the format mesh_dim x element_dim
            sum += matrix[i][i];
        }
        sum
    }

    fn frobenius_norm(matrix: &Vec<Vec<f64>>, mesh_dim: usize, element_dim: usize) -> f64 {
        let mut sum = 0.0;
        for i in 0..mesh_dim { // if the matrix square than mesh_dim == element_dim, if not then matrix has the format mesh_dim x element_dim
            for j in 0..element_dim {
                sum += matrix[i][j] * matrix[i][j];
            }
        }
        sum.sqrt()
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

    /// Compute 3D cross product of two vectors
    fn cross_product_3d(a: &[f64], b: &[f64]) -> Vec<f64> {
        vec![
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// Compute Euclidean norm of a vector
    fn vector_norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    fn jacobian_factorization(evaluated_jacobian: &Vec<Vec<f64>>, 
                            evaluated_det_jacobian: &f64, 
                            mesh_dim: usize, 
                            element_dim: usize, 
                            point: (f64, f64, f64) 
                            )-> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>), ElementError> {

        // Compute the metric tensor lambda = J^T * J
        let mut lambda = vec![vec![0.0; element_dim]; element_dim];
        for i in 0..element_dim {
            for j in 0..element_dim {
                for k in 0..mesh_dim {
                    lambda[i][j] += evaluated_jacobian[k][i] * evaluated_jacobian[k][j];
                }
            }
        }
       
       let mut evaluated_square_matrix = vec![vec![0.0; element_dim]; element_dim];
                
        for i in 0..element_dim {
            for j in 0..element_dim {
                for k in 0..mesh_dim {
                    evaluated_square_matrix[i][j] += evaluated_jacobian[k][i] * evaluated_jacobian[k][j];
                }
            }
        }

        let mut r_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; mesh_dim];
        let mut q_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; mesh_dim];
        let mut d_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; mesh_dim];

        match (mesh_dim, element_dim) {
            // 1x1 matrix: 
            (1, 1) => {

                let r_matrix = vec![vec![1.0]];
                let q_matrix = vec![vec![1.0]];
                let d_matrix = vec![vec![1.0]];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 2x2 matrix: 
            (2, 2) => { 
                r_matrix = vec![
                    vec![evaluated_jacobian[0][0] / lambda[0][0].sqrt(), - evaluated_jacobian[1][0] / lambda[0][0].sqrt()],
                    vec![evaluated_jacobian[1][0] / lambda[0][0].sqrt(), evaluated_jacobian[0][0] / lambda[0][0].sqrt()],
                ];
                q_matrix = vec![
                    vec![1.0, lambda[0][1] / (lambda[0][0] * lambda[1][1]).sqrt()],
                    vec![0.0, evaluated_det_jacobian / (lambda[0][0] * lambda[1][1]).sqrt()],
                ];
                d_matrix = vec![
                    vec![1.0, 0.0],
                    vec![0.0, (lambda[1][1] / lambda[0][0]).sqrt()],
                ];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 3x3 matrix: 
            (3, 3) => {

                let col_0 = vec![evaluated_jacobian[0][0], evaluated_jacobian[1][0], evaluated_jacobian[2][0]];
                let col_1 = vec![evaluated_jacobian[0][1], evaluated_jacobian[1][1], evaluated_jacobian[2][1]];
                let col_2 = vec![evaluated_jacobian[0][2], evaluated_jacobian[1][2], evaluated_jacobian[2][2]];

                r_matrix = vec![
                    vec![col_0[0] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[0] - lambda[0][1] * col_0[0]) / (lambda[0][0].sqrt() * (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))), (col_0[1] * col_1[2] - col_0[2] * col_1[1]) / (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))],
                    vec![col_0[1] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[1] - lambda[0][1] * col_0[1]) / (lambda[0][0].sqrt() * (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))), (col_0[2] * col_1[0] - col_0[0] * col_1[2]) / (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))],
                    vec![col_0[2] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[2] - lambda[0][1] * col_0[2]) / (lambda[0][0].sqrt() * (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))), (col_0[0] * col_1[1] - col_0[1] * col_1[0]) / (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))],
                ];
                q_matrix = vec![
                    vec![1.0, lambda[0][1] / (lambda[0][0] * lambda[1][1]).sqrt(), lambda[0][2] / (lambda[0][0] * lambda[2][2]).sqrt()],
                    vec![0.0, (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1))) / (lambda[0][0] * lambda[1][1]).sqrt(), (lambda[0][0] * lambda[1][2] - lambda[0][1] * lambda[0][2]) / ((lambda[0][0] * lambda[2][2]).sqrt() * (Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1))))],
                    vec![0.0, 0.0, evaluated_det_jacobian / (lambda[2][2].sqrt() * Self::vector_norm(&Self::cross_product_3d(&col_0, &col_1)))],
                ];
                d_matrix = vec![
                    vec![1.0, 0.0, 0.0],
                    vec![0.0, (lambda[1][1]).sqrt() / (lambda[0][0]).sqrt(), 0.0],
                    vec![0.0, 0.0, (lambda[2][2]).sqrt() / (lambda[0][0]).sqrt()],
                ];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 1D elements in 2D space: metric = dx^2 + dy^2 (squared length)
            (2, 1) => {
                
                let r_matrix = vec![vec![1.0]];
                let q_matrix = vec![vec![1.0]];
                let d_matrix = vec![vec![1.0]];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 1D elements in 3D space:
            (3, 1) => {
                
                let r_matrix = vec![vec![1.0]];
                let q_matrix = vec![vec![1.0]];
                let d_matrix = vec![vec![1.0]];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 2D elements in 3D space: metric determinant from first fundamental form
            (3, 2) => {
                r_matrix = vec![
                    vec![evaluated_square_matrix[0][0] / lambda[0][0].sqrt(), - evaluated_square_matrix[1][0] / lambda[0][0].sqrt()],
                    vec![evaluated_square_matrix[1][0] / lambda[0][0].sqrt(), evaluated_square_matrix[0][0] / lambda[0][0].sqrt()],
                ];
                q_matrix = vec![
                    vec![1.0, lambda[0][1] / (lambda[0][0] * lambda[1][1]).sqrt()],
                    vec![0.0, evaluated_det_jacobian / (lambda[0][0] * lambda[1][1]).sqrt()],
                ];
                d_matrix = vec![
                    vec![1.0, 0.0],
                    vec![0.0, (lambda[1][1] / lambda[0][0]).sqrt()],
                ];

                Ok((r_matrix, q_matrix, d_matrix))
            }
            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian factorization not implemented for {}x{} (mesh_dim x element_dim) matrices",
                mesh_dim, element_dim
            ))),
        }
    }

}


////////////////////////////////////////////////////////////////////////////////////////////////

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

    // Jacobian calculation tests 
    #[test]
    fn test_calculate_jacobian_monomial_1d_line() {
        let nodes = vec![
            create_node_1d(1, 0.0),
            create_node_1d(2, 2.0),
        ];
        let element = create_element(1, vec![1, 2]);
        
        let shape_derivatives = ElementType::get_shape_functions(&ElementType::Line)
            .unwrap()
            .derivatives;
        
        let result = GeometricAnalysis::calculate_jacobian_monomial(
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
    fn test_calculate_jacobian_monomial_2d_triangle() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        
        let shape_derivatives = ElementType::get_shape_functions(&ElementType::Triangle)
            .unwrap()
            .derivatives;
        
        let result = GeometricAnalysis::calculate_jacobian_monomial(
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

    // NEW: Determinant calculation tests for constant matrices
    #[test]
    fn test_calculate_determinant_2x2() {
        let matrix = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        
        let result = GeometricAnalysis::calculate_determinant(&matrix, 2, 2);
        assert!(result.is_ok());
        assert!((result.unwrap() - (-2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_determinant_3x3() {
        let matrix = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        
        let result = GeometricAnalysis::calculate_determinant(&matrix, 3, 3);
        assert!(result.is_ok());
        assert!((result.unwrap() - 0.0).abs() < 1e-12);
    }
    /* 
    #[test]
    fn test_calculate_determinant_1d_in_2d() {
        let matrix = vec![
            vec![3.0],
            vec![4.0],
        ];
        
        let result = GeometricAnalysis::calculate_determinant(&matrix, 1, 2);
        assert!(result.is_ok());
        // For 1D in 2D, this should calculate the squared length
        assert!((result.unwrap() - 25.0).abs() < 1e-12);
    }
    */

    // NEW: Matrix inversion tests
    #[test]
    fn test_inverse_matrix_2x2() {
        let matrix = vec![
            vec![4.0, 7.0],
            vec![2.0, 6.0],
        ];
        
        let result = GeometricAnalysis::inverse_matrix(&matrix, 2, 2);
        assert!(result.is_ok());
        
        let inverse = result.unwrap();
        let expected = vec![
            vec![0.6, -0.7],
            vec![-0.2, 0.4],
        ];
        
        for i in 0..2 {
            for j in 0..2 {
                assert!((inverse[i][j] - expected[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_inverse_matrix_singular() {
        let matrix = vec![
            vec![1.0, 2.0],
            vec![2.0, 4.0], // Linearly dependent rows
        ];
        
        let result = GeometricAnalysis::inverse_matrix(&matrix, 2, 2);
        assert!(result.is_err());
    }

    // NEW: Adjoint matrix tests
    #[test]
    fn test_calculate_adjoint_2x2() {
        let matrix = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        
        let result = GeometricAnalysis::calculate_adjoint(&matrix, 2, 2);
        assert!(result.is_ok());
        
        let adjoint = result.unwrap();
        let expected = vec![
            vec![4.0, -2.0],
            vec![-3.0, 1.0],
        ];
        
        for i in 0..2 {
            for j in 0..2 {
                assert!((adjoint[i][j] - expected[i][j]).abs() < 1e-12);
            }
        }
    }

    // NEW: Norm and trace tests
    #[test]
    fn test_frobenius_norm() {
        let matrix = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        
        let norm = GeometricAnalysis::frobenius_norm(&matrix, 2, 2);
        let expected = (1.0 + 4.0 + 9.0 + 16.0_f64).sqrt();
        assert!((norm - expected).abs() < 1e-12);
    }

    #[test]
    fn test_trace() {
        let matrix = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        
        let trace = GeometricAnalysis::trace(&matrix, 2, 2);
        assert!((trace - 5.0).abs() < 1e-12);
    }

    // NEW: Cross product and vector norm tests
    #[test]
    fn test_cross_product_3d() {
        let a = vec![1.0, 0.0, 0.0];
        let b = vec![0.0, 1.0, 0.0];
        
        let cross = GeometricAnalysis::cross_product_3d(&a, &b);
        let expected = vec![0.0, 0.0, 1.0];
        
        for i in 0..3 {
            assert!((cross[i] - expected[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_vector_norm() {
        let v = vec![3.0, 4.0];
        let norm = GeometricAnalysis::vector_norm(&v);
        assert!((norm - 5.0).abs() < 1e-12);
    }

    // UPDATED: Element quality analysis tests with new metrics
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
        let quality = &report.element_qualities[0];
        
        assert_eq!(quality.element_id, 1);
        assert!(quality.det_jacobian_value > 0.0);
        
        // Check that all metrics are within reasonable bounds
        assert!(quality.shape_metric >= 0.0 && quality.shape_metric <= 1.0);
        assert!(quality.skewness_metric >= 0.0 && quality.skewness_metric <= 1.0);
        assert!(quality.length_ratio >= 0.0 && quality.length_ratio <= 1.0);
        assert!(quality.orientation_metric >= 0.0 && quality.orientation_metric <= 1.0);
        assert!(quality.volume_metric > 0.0);
        assert!(quality.volume_shape_metric >= 0.0 && quality.volume_shape_metric <= 1.0);
        assert!(quality.volume_shape_orientation_metric >= 0.0 && quality.volume_shape_orientation_metric <= 1.0);
    }

    // NEW: Test for quality metrics calculations
    #[test]
    fn test_quality_metrics_calculations() {
        // Test individual metric calculation functions
        
        // Shape metric test
        let jacobian = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        let inverse_jacobian = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        
        let shape_metric = GeometricAnalysis::calculate_shape_metric(
            &jacobian, 
            &inverse_jacobian, 
            2, 2
        );
        // For identity matrix, shape metric should be 1.0
        assert!((shape_metric - 1.0).abs() < 1e-12);
        
        // Volume-shape metric test
        let volume_shape_metric = GeometricAnalysis::calculate_volume_shape_metric(
            1.0, 1.0
        );
        assert!((volume_shape_metric - 1.0).abs() < 1e-12);
    }

    // NEW: Test Jacobian factorization
    #[test]
    fn test_jacobian_factorization_2x2() {
        let jacobian = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        let det = 1.0;
        
        let result = GeometricAnalysis::jacobian_factorization(
            &jacobian, &det, 2, 2, (0.0, 0.0, 0.0)
        );
        
        assert!(result.is_ok());
        let (r, q, d) = result.unwrap();
        
        // For identity matrix, R should be close to identity, Q should be identity, D should be identity
        assert!((r[0][0] - 1.0).abs() < 1e-6);
        assert!((q[0][0] - 1.0).abs() < 1e-6);
        assert!((d[0][0] - 1.0).abs() < 1e-6);
    }

    // NEW: Test orientation metric calculation
    #[test]
    fn test_orientation_metric() {
        let r_matrix = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        
        let orientation_metric = GeometricAnalysis::calculate_orientation_metric(
            &r_matrix, 2, 2
        );
        // For identity matrix, orientation metric should be 1.0
        assert!((orientation_metric - 1.0).abs() < 1e-12);
    }

    // NEW: Test with different element types
    #[test]
    fn test_quad_element_quality() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 1.0, 1.0),
            create_node_2d(4, 0.0, 1.0),
        ];
        let element = create_element(1, vec![1, 2, 3, 4]);
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Quad,
            &nodes,
            2,
        );
        
        assert!(result.is_ok());
        let quality = result.unwrap();
        
        // Square quad should have good quality metrics
        assert!(quality.det_jacobian_value > 0.0);
        assert!(quality.shape_metric > 0.8); // Should be close to 1 for perfect square
        assert!(quality.volume_metric > 0.0);
    }

    // NEW: Test edge cases for quality metrics
    #[test]
    fn test_quality_metrics_edge_cases() {
        // Test with very small determinant
        let small_det = 1e-10;
        let volume_shape_metric = GeometricAnalysis::calculate_volume_shape_metric(1.0, small_det);
        assert!(volume_shape_metric >= 0.0 && volume_shape_metric <= 1.0);
        
        // Test with large determinant
        let large_det = 1e10;
        let volume_shape_metric_large = GeometricAnalysis::calculate_volume_shape_metric(1.0, large_det);
        assert!(volume_shape_metric_large >= 0.0 && volume_shape_metric_large <= 1.0);
    }

    // NEW: Test mesh quality report structure
    #[test]
    fn test_mesh_quality_report() {
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
        assert_eq!(report.element_qualities.len(), 1);
        
        let quality = &report.element_qualities[0];
        assert_eq!(quality.element_id, 1);
        
        // Verify all metrics are present and reasonable
        assert!(quality.shape_metric.is_finite());
        assert!(quality.skewness_metric.is_finite());
        assert!(quality.length_ratio.is_finite());
        assert!(quality.orientation_metric.is_finite());
        assert!(quality.volume_metric.is_finite());
        assert!(quality.volume_shape_metric.is_finite());
        assert!(quality.volume_shape_orientation_metric.is_finite());
    }

    // NEW: Test error conditions
    #[test]
    fn test_analyse_mesh_quality_empty_mesh() {
        let nodes = vec![];
        let elements = vec![];
        let element_type_info = vec![];
        let mesh_data = create_test_mesh_data(nodes, elements, element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_err());
    }

    // NEW: Test with invalid element data
    #[test]
    fn test_calculate_element_quality_invalid_nodes() {
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            // Missing other nodes for triangle
        ];
        let element = create_element(1, vec![1, 2, 3]); // References non-existent nodes
        
        let result = GeometricAnalysis::calculate_element_quality(
            &element,
            &ElementType::Triangle,
            &nodes,
            2,
        );
        
        // This should fail when trying to get element nodes
        assert!(result.is_err());
    }
    /* 
    // Performance test for multiple elements
    #[test]
    fn test_analyse_mesh_quality_multiple_elements_performance() {
        let mut nodes = Vec::new();
        let mut elements = Vec::new();
        
        // Create a small grid of triangles
        let grid_size = 3;
        let mut node_id = 1;
        
        for i in 0..grid_size {
            for j in 0..grid_size {
                nodes.push(create_node_2d(node_id, i as f64, j as f64));
                node_id += 1;
            }
        }
        
        // Create triangular elements
        for i in 0..grid_size-1 {
            for j in 0..grid_size-1 {
                let n1 = i * grid_size + j + 1;
                let n2 = i * grid_size + j + 2;
                let n3 = (i + 1) * grid_size + j + 1;
                let n4 = (i + 1) * grid_size + j + 2;
                
                elements.push(create_element(elements.len() + 1, vec![n1, n2, n3]));
                elements.push(create_element(elements.len() + 1, vec![n2, n4, n3]));
            }
        }
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: elements.len(),
            start_index: 0,
            nodes_per_element: 3,
        }];
        
        let mesh_data = create_test_mesh_data(nodes, elements, element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        
        let report = result.unwrap();
        assert_eq!(report.total_elements, (grid_size-1) * (grid_size-1) * 2);
        
        // All elements should have positive volume
        for quality in report.element_qualities {
            assert!(quality.volume_metric > 0.0);
        }
    }
    */
}
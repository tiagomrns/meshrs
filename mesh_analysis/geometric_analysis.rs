use core::num;
use std::f64;
use std::vec;

use crate::mesh_analysis::gaussian_quadrature::GaussianQuadrature;
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
        
        // Get multiple evaluation points
        let evaluation_points = Self::get_evaluation_points(element_type);

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
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        )?;

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_zero_detj) = Self::find_min_max_jacobian_determinant(
            &jacobian_matrix.determinant,
            &evaluation_points,
        )?;

        // EARLY RETURN: If ANY point has zero Jacobian, return all metrics as 0
        if has_zero_detj {
            return Ok(ElementQuality {
                element_id: element.id,
                shape_metric: 0.0,
                skewness_metric: 0.0,
                length_ratio: 0.0,
                orientation_metric: 0.0,
                signed_volume_metric: 0.0,
                jacobian_ratio: 0.0,
            });
        }
        
        // Calculate Jacobian ratio (max/min)
        let jacobian_ratio = if (min_detj.abs() < 1e-12) || (max_detj.abs() < 1e-12) {
            0.0 // either min detJ or max detJ is zero means invalid element
        } else {
            max_detj / min_detj
        };

        // Convert worst point to tuple format for evaluation (we'll use min point for detailed analysis)
        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };
        
        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; element_dim]; mesh_dim];
        for i in 0..mesh_dim {
            for j in 0..element_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix.matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                ))?;
            }
        }

        // Calculate orientation sign only for non-square matrices (-1.0 or 1.0)
        let orientation_sign = if mesh_dim != element_dim {
            Self::calculate_element_orientation_sign(
                &evaluated_jacobian,
                mesh_dim,
                element_dim,
            )?
        } else {
            1.0 // Default to 1.0 for square matrices
        };

        let (r_matrix, q_matrix, d_matrix) = 
            Self::jacobian_factorization(&evaluated_jacobian, &min_detj, mesh_dim, element_dim)?;

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = Self::inverse_matrix(&evaluated_jacobian, mesh_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        let inverse_q_matrix = Self::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        let inverse_d_matrix = Self::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        

        // Calculate metrics
        let shape_metric = Self::calculate_shape_metric(&evaluated_jacobian, &inverse_evaluated_jacobian, mesh_dim, element_dim);
        let skewness_metric = Self::calculate_skewness_metric(&q_matrix, &inverse_q_matrix, mesh_dim, element_dim);
        let length_ratio = Self::calculate_length_ratio(&d_matrix, &inverse_d_matrix, mesh_dim, element_dim);
        let orientation_metric = Self::calculate_orientation_metric(&r_matrix, mesh_dim, element_dim);
        let volume_metric = min_detj;
        let signed_volume_metric = orientation_sign * volume_metric;
        // more metrics can be added here if needed

        
        Ok(ElementQuality {
            element_id: element.id,
            shape_metric,
            skewness_metric,
            length_ratio,
            orientation_metric,
            signed_volume_metric,   // Signed volume metric
            jacobian_ratio,         // max(detJ) / min(detJ)
        })
    }

    // Find both the minimum and maximum Jacobian determinant and check if there is detJ = 0
    fn find_min_max_jacobian_determinant(
        determinant_poly: &Vec<f64>,
        evaluation_points: &[Vec<f64>],
    ) -> Result<(Vec<f64>, f64, f64, bool), ElementError> {
        if evaluation_points.is_empty() {
            return Err(ElementError::GeometryError(
                "No evaluation points provided".to_string()
            ));
        }

        let mut min_point = &evaluation_points[0];
        let mut min_det = f64::MAX;
        let mut max_det = f64::MIN;
        let mut has_zero = false;
        
        for point in evaluation_points {
            let eval_point = match point.len() {
                1 => (point[0], 0.0, 0.0),
                2 => (point[0], point[1], 0.0),
                3 => (point[0], point[1], point[2]),
                _ => (0.0, 0.0, 0.0),
            };
            
            let det_value = MonomialPolynomial::evaluate(determinant_poly, eval_point)
                .map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate determinant at point {:?}: {}", point, e)
                ))?;
            
            // Check for zero determinant
            if det_value.abs() < 1e-12 {
                has_zero = true;
                // We can break early since we found a zero
                break;
            }
            
            // Find minimum determinant (worst case)
            if det_value < min_det {
                min_det = det_value;
                min_point = point;
            }
            
            // Find maximum determinant
            if det_value > max_det {
                max_det = det_value;
            }
        }

        Ok((min_point.clone(), min_det, max_det, has_zero))
    }


    /// Calculate orientation sign for non-square Jacobian matrices
    /// Preserves directional information that would be lost by taking absolute values
    fn calculate_element_orientation_sign(
        jacobian: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> Result<f64, ElementError> {
        // This function is now only called for non-square matrices
        match (mesh_dim, element_dim) {
            // 1D elements in 2D or 3D space: sign from largest tangent component
            (2, 1) | (3, 1) => {
                // Extract tangent vector from Jacobian columns
                let mut tangent = vec![0.0; mesh_dim];
                for i in 0..mesh_dim {
                    tangent[i] = jacobian[i][0];
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
                
                // Return sign based on dominant component
                Ok(if max_abs < 1e-12 { 1.0 } else if sign_component >= 0.0 { 1.0 } else { -1.0 })
            },
            
            // 2D elements in 3D space: sign from surface normal
            (3, 2) => {
                // Extract both tangent vectors from Jacobian columns
                let u = [jacobian[0][0], jacobian[1][0], jacobian[2][0]];
                let v = [jacobian[0][1], jacobian[1][1], jacobian[2][1]];
                
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
            
            // Default case (shouldn't normally occur for non-square)
            _ => Ok(1.0),
        }
    }

    /// Calculate shape metric
    fn calculate_shape_metric(
        jacobian: &Vec<Vec<f64>>,
        inverse_jacobian: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> f64 {
        
        let j_norm = Self::frobenius_norm(jacobian, element_dim);
        let inv_j_norm = Self::frobenius_norm(inverse_jacobian, element_dim);

        element_dim as f64 / (j_norm * inv_j_norm)
    }

    /// Calculate skewness metric
    fn calculate_skewness_metric(
        q_matrix: &Vec<Vec<f64>>,
        inverse_q_matrix: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> f64 {
        
        let q_norm = Self::frobenius_norm(q_matrix, element_dim);
        let inv_q_norm = Self::frobenius_norm(inverse_q_matrix, element_dim);
        
        element_dim as f64 / (q_norm * inv_q_norm)
    }

    /// Calculate length ratio metric
    fn calculate_length_ratio(
        d_matrix: &Vec<Vec<f64>>,
        inverse_d_matrix: &Vec<Vec<f64>>,
        mesh_dim: usize,
        element_dim: usize,
    ) -> f64 {

        let d_norm = Self::frobenius_norm(d_matrix, element_dim);
        let inv_d_norm = Self::frobenius_norm(inverse_d_matrix, element_dim);

        element_dim as f64 / (d_norm * inv_d_norm)
    }

    /// Calculate orientation metric
    fn calculate_orientation_metric(
        r_matrix: &Vec<Vec<f64>>, 
        mesh_dim: usize, 
        element_dim: usize
    ) -> f64 {

        1.0 + ((Self::trace(r_matrix, element_dim) - element_dim as f64) / 4.0)
    }

    // Return evaluation points for each element type to find min/max Jacobian determinant
    // corner and center points of each element
    fn get_evaluation_points(element_type: &ElementType) -> Vec<Vec<f64>> {
        match element_type {
            // 1D elements - check endpoints and center
            ElementType::Line | ElementType::QuadraticEdge => {
                vec![
                    vec![0.0],
                    vec![1.0],
                    vec![0.5], // Center
                ]
            }
            // 2D triangular elements - check vertices and centroid
            ElementType::Triangle | ElementType::QuadraticTriangle => {
                vec![
                    vec![0.0, 0.0],
                    vec![1.0, 0.0],
                    vec![0.0, 1.0],
                    vec![1.0/3.0, 1.0/3.0], // Centroid
                ]
            }
            // 2D quadrilateral elements - check corners and center
            ElementType::Quad | ElementType::QuadraticQuad | ElementType::BiquadraticQuad => {
                vec![
                    vec![0.0, 0.0],    
                    vec![1.0, 0.0],
                    vec![1.0, 1.0],
                    vec![0.0, 1.0],
                    vec![0.5, 0.5], // Center
                ]
            }
            // 3D tetrahedral elements - vertices and centroid  
            ElementType::Tetra | ElementType::QuadraticTetra => {
                vec![
                    vec![0.0, 0.0, 0.0],
                    vec![1.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 1.0],
                    vec![0.25, 0.25, 0.25],  // Centroid
                ]
            }
            // Pyramid element - key points including apex
            ElementType::Pyramid => {
                vec![
                    vec![0.0, 0.0, 0.0],    
                    vec![1.0, 0.0, 0.0],     
                    vec![1.0, 1.0, 0.0],      
                    vec![0.0, 1.0, 0.0],   
                    vec![0.0, 0.0, 1.0],  
                ]
            }
            // Wedge/prism elements - key locations
            ElementType::Wedge | ElementType::QuadraticWedge | ElementType::BiquadraticQuadraticWedge => {
                vec![
                    vec![0.0, 0.0, 0.0],       
                    vec![1.0, 0.0, 0.0],       
                    vec![0.0, 1.0, 0.0],       
                    vec![0.0, 0.0, 1.0],       
                    vec![1.0, 0.0, 1.0],       
                    vec![0.0, 1.0, 1.0],       
                    vec![0.5, 0.5, 0.5],       // Centroid
                ]
            }
            // 3D hexahedral elements - corners and center
            ElementType::Hexahedron | ElementType::QuadraticHexahedron | ElementType::BiquadraticQuadraticHexahedron | ElementType::TriquadraticHexahedron => {
                vec![
                    vec![0.0, 0.0, 0.0],        
                    vec![1.0, 0.0, 0.0],
                    vec![1.0, 1.0, 0.0],
                    vec![0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 1.0],
                    vec![1.0, 0.0, 1.0],
                    vec![1.0, 1.0, 1.0],
                    vec![0.0, 1.0, 1.0],
                    vec![0.5, 0.5, 0.5], // Center
                ]
            }
            // Vertex has no integration point
            ElementType::Vertex => vec![vec![0.0]],
            _ => vec![vec![0.0]], // Default case
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
                
                Ok(a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e))
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
                if matrix[0][0].abs() < 1e-12 {
                    return Err(ElementError::GeometryError("Matrix is singular".to_string()));
                }
                Ok(vec![vec![1.0 / matrix[0][0]]])
            },

            // 2x2 and 3x3 matrices
            (2, 2) | (3, 3) => {
                let adjoint = Self::calculate_adjoint(&matrix, element_dim, element_dim)?;
                let determinant = Self::calculate_determinant(&matrix, element_dim, element_dim)?;

                if determinant.abs() < 1e-12 {
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

            // 1D elements in 2D and 3D spaces
            (2, 1) | (3, 1) => {
                let mut sum = 0.0;
                for k in 0..mesh_dim {
                    sum += matrix[k][0] * matrix[k][0];
                }
                if sum.abs() < 1e-12 {
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
                
                if determinant.abs() < 1e-12 {
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

    fn trace(matrix: &Vec<Vec<f64>>, element_dim: usize) -> f64 {
        let mut sum = 0.0;
        for i in 0..element_dim {
            sum += matrix[i][i];
        }
        sum
    }

    fn frobenius_norm(matrix: &Vec<Vec<f64>>, element_dim: usize) -> f64 {
        let mut sum = 0.0;
        for i in 0..element_dim {
            for j in 0..element_dim {
                sum += matrix[i][j] * matrix[i][j];
            }
        }
        sum.sqrt()
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
       
       let evaluated_square_matrix = lambda.clone();
                

        let mut r_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; element_dim];
        let mut q_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; element_dim];
        let mut d_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; element_dim];

        match (mesh_dim, element_dim) {
            // 1x1 matrix: 
            (1, 1) => {

                let r_matrix = vec![vec![1.0]];
                let q_matrix = vec![vec![1.0]];
                let d_matrix = vec![vec![evaluated_jacobian[0][0]]];

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
                // let col_2 = vec![evaluated_jacobian[0][2], evaluated_jacobian[1][2], evaluated_jacobian[2][2]];

                let cross_product_01 = Self::cross_product_3d(&col_0, &col_1);
                let norm_cross_product_01 = Self::vector_norm(&cross_product_01);

                r_matrix = vec![
                    vec![col_0[0] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[0] - lambda[0][1] * col_0[0]) / (lambda[0][0].sqrt() * norm_cross_product_01), cross_product_01[0] / norm_cross_product_01],
                    vec![col_0[1] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[1] - lambda[0][1] * col_0[1]) / (lambda[0][0].sqrt() * norm_cross_product_01), cross_product_01[1] / norm_cross_product_01],
                    vec![col_0[2] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[2] - lambda[0][1] * col_0[2]) / (lambda[0][0].sqrt() * norm_cross_product_01), cross_product_01[2] / norm_cross_product_01],
                ];
                q_matrix = vec![
                    vec![1.0, lambda[0][1] / (lambda[0][0] * lambda[1][1]).sqrt(), lambda[0][2] / (lambda[0][0] * lambda[2][2]).sqrt()],
                    vec![0.0, norm_cross_product_01 / (lambda[0][0] * lambda[1][1]).sqrt(), (lambda[0][0] * lambda[1][2] - lambda[0][1] * lambda[0][2]) / ((lambda[0][0] * lambda[2][2]).sqrt() * norm_cross_product_01)],
                    vec![0.0, 0.0, evaluated_det_jacobian / (lambda[2][2].sqrt() * norm_cross_product_01)],
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
                let d_matrix = vec![vec![lambda[0][0].sqrt()]];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 1D elements in 3D space:
            (3, 1) => {
                
                let r_matrix = vec![vec![1.0]];
                let q_matrix = vec![vec![1.0]];
                let d_matrix = vec![vec![lambda[0][0].sqrt()]];

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

    // Determinant calculation tests for constant matrices
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
        
        let norm = GeometricAnalysis::frobenius_norm(&matrix,  2);
        let expected = (1.0 + 4.0 + 9.0 + 16.0_f64).sqrt();
        assert!((norm - expected).abs() < 1e-12);
    }

    #[test]
    fn test_trace() {
        let matrix = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        
        let trace = GeometricAnalysis::trace(&matrix, 2);
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

    // NEW: Test orientation sign calculation
    #[test]
    fn test_calculate_element_orientation_sign() {
        // Test with identity matrix (2D)
        let jacobian_2d = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        
        let result = GeometricAnalysis::calculate_element_orientation_sign(&jacobian_2d, 2, 2);
        assert!(result.is_ok());
        assert!((result.unwrap() - 1.0).abs() < 1e-12);
        
        // Test with negative determinant (2D)
        let jacobian_2d_neg = vec![
            vec![1.0, 0.0],
            vec![0.0, -1.0],
        ];
        
        let result = GeometricAnalysis::calculate_element_orientation_sign(&jacobian_2d_neg, 2, 2);
        assert!(result.is_ok());
        assert!((result.unwrap() - (-1.0)).abs() < 1e-12);
        
        // Test 1D element in 2D space
        let jacobian_1d = vec![
            vec![1.0],
            vec![2.0],
        ];
        
        let result = GeometricAnalysis::calculate_element_orientation_sign(&jacobian_1d, 2, 1);
        assert!(result.is_ok());
        assert!((result.unwrap() - 1.0).abs() < 1e-12);
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
        //assert!(quality.det_jacobian_value > 0.0);
        
        // Check that all metrics are within reasonable bounds
        assert!(quality.shape_metric >= 0.0 && quality.shape_metric <= 1.0);
        assert!(quality.skewness_metric >= 0.0 && quality.skewness_metric <= 1.0);
        assert!(quality.length_ratio >= 0.0 && quality.length_ratio <= 1.0);
        assert!(quality.orientation_metric >= 0.0 && quality.orientation_metric <= 1.0);
        //assert!(quality.volume_metric > 0.0);
        //assert!(quality.volume_shape_metric >= 0.0 && quality.volume_shape_metric <= 1.0);
        //assert!(quality.volume_shape_orientation_metric >= 0.0 && quality.volume_shape_orientation_metric <= 1.0);
        
        // NEW: Test Jacobian ratio and scaled Jacobian
        assert!(quality.jacobian_ratio >= 1.0); // Should be at least 1.0
        //assert!(quality.scaled_jacobian >= 0.0 && quality.scaled_jacobian <= 1.0);
        //assert!(!quality.min_detj_point.is_empty());
        //assert!(!quality.max_detj_point.is_empty());
        //assert!(quality.max_detj_value > 0.0);
        //assert!(quality.orientation_sign == 1.0 || quality.orientation_sign == -1.0);
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
        //let volume_shape_metric = GeometricAnalysis::calculate_volume_shape_metric(
        //    1.0, 1.0
        //);
        //assert!((volume_shape_metric - 1.0).abs() < 1e-12);
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
            &jacobian, &det, 2, 2
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

    // NEW: Test scaled Jacobian calculation
    #[test]
    fn test_scaled_jacobian() {
        // Test with identity matrix
        let jacobian = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        
        //let result = GeometricAnalysis::calculate_scaled_jacobian(&jacobian, 2, 2);
        //assert!(result.is_ok());
        //assert!((result.unwrap() - 1.0).abs() < 1e-12);
    }

    // NEW: Test min-max Jacobian point finding
    #[test]
    fn test_find_min_max_jacobian_points() {
        // Create a simple determinant polynomial (constant)
        let determinant_poly = vec![2.0]; // Constant polynomial: detJ = 2.0 everywhere
        let evaluation_points = vec![
            vec![0.0],
            vec![0.5],
            vec![1.0],
        ];
        
        let result = GeometricAnalysis::find_min_max_jacobian_determinant(
            &determinant_poly,
            &evaluation_points
        );
        
        assert!(result.is_ok());
        let (min_point, min_det, max_det) = result.unwrap();
        
        // All points should have the same determinant value
        assert!((min_det - 2.0).abs() < 1e-12);
        assert!((max_det - 2.0).abs() < 1e-12);
        assert!(!min_point.is_empty());
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
        //assert!(quality.det_jacobian_value > 0.0);
        assert!(quality.shape_metric > 0.8); // Should be close to 1 for perfect square
        //assert!(quality.volume_metric > 0.0);
        assert!(quality.jacobian_ratio >= 1.0);
        //assert!(quality.scaled_jacobian > 0.8); // Should be close to 1 for perfect square
        //assert!(quality.orientation_sign == 1.0 || quality.orientation_sign == -1.0);
    }

    // NEW: Test edge cases for quality metrics
    #[test]
    fn test_quality_metrics_edge_cases() {
        // Test with very small determinant
        let small_det = 1e-10;
        //let volume_shape_metric = GeometricAnalysis::calculate_volume_shape_metric(1.0, small_det);
        //assert!(volume_shape_metric >= 0.0 && volume_shape_metric <= 1.0);
        
        // Test with large determinant
        let large_det = 1e10;
        //let volume_shape_metric_large = GeometricAnalysis::calculate_volume_shape_metric(1.0, large_det);
        //assert!(volume_shape_metric_large >= 0.0 && volume_shape_metric_large <= 1.0);
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
        //assert!(quality.volume_metric.is_finite());
        //assert!(quality.volume_shape_metric.is_finite());
        //assert!(quality.volume_shape_orientation_metric.is_finite());
        assert!(quality.jacobian_ratio.is_finite());
        //assert!(quality.scaled_jacobian.is_finite());
        //assert!(!quality.min_detj_point.is_empty());
        //assert!(!quality.max_detj_point.is_empty());
        //assert!(quality.max_detj_value.is_finite());
        //assert!(quality.orientation_sign == 1.0 || quality.orientation_sign == -1.0);
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
}
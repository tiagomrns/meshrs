// Standard library imports for floating-point and vector operations
use std::f64;
use std::vec;

// Import Gaussian quadrature utilities for numerical integration
use crate::mesh_analysis::gaussian_quadrature::GaussianQuadrature;
// Import mesh data structures (MeshData, Element, Node, etc.)
use crate::structs_and_impls::*;
// Import custom error types for element-related errors
use crate::error::*;

/// GeometricAnalysis provides methods for computing mesh quality metrics
/// based on Knupp's quality metrics (shape, skewness, length ratio, orientation, volume, jacobian ratio).
/// These metrics assess the geometric quality of finite element mesh elements.
pub struct GeometricAnalysis;

impl GeometricAnalysis {
    /// Analyze mesh quality over all elements (skips Vertex elements).
    /// 
    /// Computes quality metrics for each element in the mesh and returns a comprehensive report.
    /// Vertex elements are skipped as they have no geometric extent.
    /// 
    /// # Arguments
    /// * `mesh_data` - Reference to the mesh containing elements, nodes, and element type information
    /// 
    /// # Returns
    /// * `Ok(MeshQualityReport)` - Report containing quality metrics for all successfully analyzed elements
    /// * `Err(ElementError)` - If no elements could be analyzed
    pub fn analyse_mesh_quality(mesh_data: &MeshData) -> Result<MeshQualityReport, ElementError> {
        // Initialize storage for quality metrics of all elements
        let mut element_qualities = Vec::new();
        let mut processed_elements = 0;
        let mut failed_elements = 0;
        let mesh_dim = mesh_data.dimension;

        println!("Starting mesh quality analysis for {} elements", mesh_data.elements.len());

        // Iterate through all element types in the mesh
        for type_info in &mesh_data.element_type_info {
            // Skip vertex elements as they have no geometric quality and quadratic pyramids (due to missing node issue)
            if matches!(type_info.element_type, ElementType::Vertex | ElementType::QuadraticPyramid) {
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
                            failed_elements += 1;
                        }
                    }
                }
            }
        }

        println!("Mesh quality analysis completed: {} succeeded, {} failed", 
                 processed_elements, failed_elements);

        // Return error if no elements could be analyzed
        if element_qualities.is_empty() {
            return Err(ElementError::GeometryError(
                format!("No elements could be analyzed ({} failed)", failed_elements)
            ));
        }

        // Return comprehensive quality report
        Ok(MeshQualityReport {
            total_elements: processed_elements,
            element_qualities,
        })
    }

    /// Calculate quality metrics for a single element.
    /// 
    /// This is the core function that computes all geometric quality metrics for an individual element
    /// using Knupp's quality metrics framework. The metrics include:
    /// - Shape metric: Measures deviation from ideal shape
    /// - Skewness metric: Measures angular deviation from orthogonality
    /// - Length ratio: Measures aspect ratio of element edges
    /// - Orientation metric: Measures element orientation relative to coordinate axes
    /// - Volume metric: The minimum Jacobian determinant (volume/area scaling factor)
    /// - Jacobian ratio: Ratio of minimum to maximum Jacobian determinant
    /// 
    /// # Arguments
    /// * `element` - The element to analyze
    /// * `element_type` - The type of the element (Triangle, Quad, Tetra, Hex, etc.)
    /// * `nodes` - All nodes in the mesh (used to look up element nodes)
    /// * `mesh_dim` - Dimension of the mesh space (1D, 2D, or 3D)
    /// 
    /// # Returns
    /// * `Ok(ElementQuality)` - Quality metrics for the element
    /// * `Err(ElementError)` - If the element cannot be analyzed
    fn calculate_element_quality(
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        mesh_dim: usize,
    ) -> Result<ElementQuality, ElementError> {
        // Get node coordinates for this element by looking up node IDs
        let element_nodes = Self::get_element_nodes(element, nodes)?;
        
        // Get evaluation points where Jacobian will be sampled (corners, center, etc.)
        // These points are in the reference element coordinate system
        let evaluation_points = Self::get_evaluation_points(element_type);

        // Get element dimensionality (1D for edges, 2D for faces, 3D for volumes)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();

        // Get element order (1 for linear, 2 for quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        // Shape functions N_i map reference coordinates to physical coordinates: x = Σ N_i * x_i
        let shape_function = ElementType::get_shape_functions(element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            })?;

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix J = ∂x/∂ξ as a polynomial in reference coordinates
        // The Jacobian maps differential lengths/areas/volumes from reference to physical space
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        )?;

        // Find minimum and maximum Jacobian determinant values across evaluation points
        // det(J) represents the local volume/area/length scaling factor
        // Positive det(J) indicates proper element orientation; zero or negative indicates degenerate/inverted element
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = Self::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        )?;

        // EARLY RETURN: If ANY evaluation point has non-positive Jacobian determinant,
        // the element is degenerate or inverted. Return all metrics as 0 to indicate complete failure.
        if has_nonpositive_detj {
            return Ok(ElementQuality {
                element_id: element.id,
                shape_metric: 0.0,
                skewness_metric: 0.0,
                length_ratio: 0.0,
                orientation_metric: 0.0,
                volume_metric: 0.0,
                jacobian_ratio: 0.0,
            });
        }
        
        // Calculate Jacobian ratio: min(det J) / max(det J)
        // Ideal value is 1.0 (uniform scaling throughout element)
        // Values close to 0 indicate severe distortion
        let jacobian_ratio = if (min_detj.abs() < 1e-12) || (max_detj.abs() < 1e-12) {
            0.0 // Either min or max det(J) is zero, indicating invalid element
        } else {
            min_detj / max_detj
        };

        // Convert worst point (minimum Jacobian) to tuple format for polynomial evaluation
        // We analyze quality at the worst point since that determines element usability
        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),        // 1D: (ξ, 0, 0)
            2 => (min_point[0], min_point[1], 0.0), // 2D: (ξ, η, 0)
            3 => (min_point[0], min_point[1], min_point[2]), // 3D: (ξ, η, ζ)
            _ => (0.0, 0.0, 0.0),
        };
        
        // Evaluate Jacobian matrix at the worst (minimum Jacobian) point
        // This gives us the actual numeric Jacobian matrix J at that point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                // Evaluate each polynomial entry J[i,j](ξ) at the evaluation point
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                ))?;
            }
        }

        // Construct metric Jacobian matrix for Knupp's quality metrics
        // Knupp's metrics require a square matrix, so we handle non-square cases
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices (element embedded in same dimensional space):
            // Use Jacobian as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices (element embedded in higher dimensional space):
            // Compute metric tensor J * J^T to get a square matrix [element_dim × element_dim]
            // This represents the induced metric on the element manifold
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        // Perform Knupp's Jacobian factorization: J = R * Q * D
        // R = rotation/reflection matrix (measures orientation)
        // Q = skew matrix (measures angular distortion/non-orthogonality)
        // D = diagonal scaling matrix (measures aspect ratio)
        let (r_matrix, q_matrix, d_matrix) = 
            Self::jacobian_factorization(&metric_jacobian, element_dim)?;

        // Calculate matrix inverses needed for quality metrics
        // These inverses are used in computing normalized condition numbers
        let inverse_metric_jacobian = Self::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        let inverse_q_matrix = Self::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        let inverse_d_matrix = Self::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e)))?;
        

        // Calculate Knupp's quality metrics (all range from 0 to 1, where 1 is ideal)
        
        // Shape metric: n / (||J|| * ||J^-1||)
        // Measures overall shape quality combining size and distortion
        let shape_metric = Self::calculate_shape_metric(&metric_jacobian, &inverse_metric_jacobian, element_dim);
        
        // Skewness metric: n / (||Q|| * ||Q^-1||)
        // Measures angular distortion (deviation from orthogonality)
        let skewness_metric = Self::calculate_skewness_metric(&q_matrix, &inverse_q_matrix, element_dim);
        
        // Length ratio metric: n / (||D|| * ||D^-1||)
        // Measures aspect ratio (ratio of longest to shortest edge)
        let length_ratio = Self::calculate_length_ratio(&d_matrix, &inverse_d_matrix, element_dim);
        
        // Orientation metric: 1 + (trace(R) - n) / 4
        // Measures element orientation relative to coordinate axes
        let orientation_metric = Self::calculate_orientation_metric(&r_matrix, element_dim);
        
        // Volume metric: minimum det(J)
        // The minimum Jacobian determinant (local volume/area scaling factor)
        let volume_metric = min_detj;

        
        // Return all quality metrics for this element
        Ok(ElementQuality {
            element_id: element.id,
            shape_metric,
            skewness_metric,
            length_ratio,
            orientation_metric,
            volume_metric,          // Minimum Jacobian determinant
            jacobian_ratio,         // Ratio: min(det J) / max(det J)
        })
    }

    /// Find both the minimum and maximum Jacobian determinant across evaluation points.
    /// 
    /// Evaluates the Jacobian determinant at multiple points (corners, center, etc.) to find
    /// the range of scaling factors across the element. Also checks for non-positive determinants
    /// which indicate degenerate or inverted elements.
    /// 
    /// # Arguments
    /// * `jacobian_matrix` - Jacobian matrix as polynomial coefficients [element_dim][mesh_dim][coefficients]
    /// * `evaluation_points` - Points in reference coordinates where to evaluate det(J)
    /// * `element_dim` - Dimension of the element (1, 2, or 3)
    /// * `mesh_dim` - Dimension of the mesh space (1, 2, or 3)
    /// 
    /// # Returns
    /// * Tuple: (point_of_min, min_det, max_det, has_nonpositive)
    ///   - point_of_min: Reference coordinates where minimum occurs
    ///   - min_det: Minimum Jacobian determinant value
    ///   - max_det: Maximum Jacobian determinant value
    ///   - has_nonpositive: True if any det(J) ≤ 0 (degenerate/inverted element)
    fn find_min_max_jacobian_determinant(
        jacobian_matrix: &Vec<Vec<Vec<f64>>>,  // Jacobian matrix as polynomial [element_dim][mesh_dim][coefficients]
        evaluation_points: &[Vec<f64>],
        element_dim: usize,
        mesh_dim: usize,
    ) -> Result<(Vec<f64>, f64, f64, bool), ElementError> {
        // Validate that we have evaluation points
        if evaluation_points.is_empty() {
            return Err(ElementError::GeometryError("No evaluation points provided".to_string()));
        }

        // Validate Jacobian matrix dimensions
        if jacobian_matrix.len() != element_dim {
            return Err(ElementError::GeometryError(format!(
                "Jacobian matrix has {} rows, expected {} for element dimension",
                jacobian_matrix.len(), element_dim
            )));
        }
        for row in jacobian_matrix {
            if row.len() != mesh_dim {
                return Err(ElementError::GeometryError(format!(
                    "Jacobian matrix row has {} columns, expected {} for mesh dimension",
                    row.len(), mesh_dim
                )));
            }
        }

        // Initialize tracking variables
        let mut min_point = &evaluation_points[0];
        let mut min_det = f64::MAX;  // Start with maximum possible value
        let mut max_det = f64::MIN;  // Start with minimum possible value
        let mut has_nonpositive = false;
        
        // Evaluate Jacobian determinant at each sampling point
        for point in evaluation_points {
            // Convert point vector to tuple format required by polynomial evaluator
            let eval_point = match point.len() {
                1 => (point[0], 0.0, 0.0),               // 1D point
                2 => (point[0], point[1], 0.0),          // 2D point
                3 => (point[0], point[1], point[2]),     // 3D point
                _ => (0.0, 0.0, 0.0),
            };
            
            // Evaluate Jacobian matrix at this point: J(ξ) → J_numeric
            let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];
            
            for i in 0..element_dim {
                for j in 0..mesh_dim {
                    // Evaluate polynomial J[i,j](ξ) at this point
                    evaluated_jacobian[i][j] = MonomialPolynomial::evaluate_safe(
                        &jacobian_matrix[i][j], eval_point
                    );
                }
            }
            
            // Calculate determinant of the evaluated Jacobian matrix
            // det(J) represents the local volume/area/length scaling factor
            let det_value = Self::calculate_determinant(
                &evaluated_jacobian, 
                element_dim, 
                mesh_dim
            )?;
            
            // Check for non-positive determinant (indicates degenerate or inverted element)
            if det_value <= 1e-12 {
                has_nonpositive = true;
                min_point = point;
                min_det = 0.0;
                max_det = 0.0;
                // Break early since we found a degenerate point
                break;
            }
            
            // Track minimum determinant (worst quality point)
            if det_value < min_det {
                min_det = det_value;
                min_point = point;
            }
            
            // Track maximum determinant
            if det_value > max_det {
                max_det = det_value;
            }
        }

        Ok((min_point.clone(), min_det, max_det, has_nonpositive))
    }

    /// Calculate shape metric using Knupp's formula.
    /// 
    /// Shape metric = n / (||J|| * ||J^-1||)
    /// where n is the element dimension and || || is the Frobenius norm.
    /// 
    /// This metric combines size and distortion effects. For an ideal element (perfect reference shape),
    /// the value is 1.0. Distorted elements have values < 1.0, approaching 0 for severely distorted elements.
    /// 
    /// # Arguments
    /// * `jacobian` - The Jacobian matrix
    /// * `inverse_jacobian` - The inverse of the Jacobian matrix
    /// * `element_dim` - Dimension of the element (1, 2, or 3)
    /// 
    /// # Returns
    /// * Shape quality metric in range [0, 1], where 1 is ideal
    fn calculate_shape_metric(
        jacobian: &Vec<Vec<f64>>,
        inverse_jacobian: &Vec<Vec<f64>>,
        element_dim: usize,
    ) -> f64 {
        // Compute Frobenius norm of J: ||J|| = sqrt(Σ J[i,j]²)
        let j_norm = Self::frobenius_norm(jacobian, element_dim);
        // Compute Frobenius norm of J^-1
        let inv_j_norm = Self::frobenius_norm(inverse_jacobian, element_dim);

        // Shape metric: n / (||J|| * ||J^-1||)
        // Product ||J|| * ||J^-1|| is the condition number of J
        element_dim as f64 / (j_norm * inv_j_norm)
    }

    /// Calculate skewness metric using Knupp's formula.
    /// 
    /// Skewness metric = n / (||Q|| * ||Q^-1||)
    /// where Q is the skew matrix from J = R * Q * D factorization.
    /// 
    /// This metric measures angular distortion (deviation from orthogonality).
    /// Value of 1.0 indicates perfectly orthogonal element edges.
    /// Values < 1.0 indicate skewed angles between element edges.
    /// 
    /// # Arguments
    /// * `q_matrix` - The Q (skew) matrix from Jacobian factorization
    /// * `inverse_q_matrix` - The inverse of Q
    /// * `element_dim` - Dimension of the element
    /// 
    /// # Returns
    /// * Skewness quality metric in range [0, 1], where 1 is ideal (orthogonal)
    fn calculate_skewness_metric(
        q_matrix: &Vec<Vec<f64>>,
        inverse_q_matrix: &Vec<Vec<f64>>,
        element_dim: usize,
    ) -> f64 {
        let q_norm = Self::frobenius_norm(q_matrix, element_dim);
        let inv_q_norm = Self::frobenius_norm(inverse_q_matrix, element_dim);
        
        // Skewness metric: n / (||Q|| * ||Q^-1||)
        element_dim as f64 / (q_norm * inv_q_norm)
    }

    /// Calculate length ratio metric using Knupp's formula.
    /// 
    /// Length ratio = n / (||D|| * ||D^-1||)
    /// where D is the diagonal scaling matrix from J = R * Q * D factorization.
    /// 
    /// This metric measures the aspect ratio (ratio of longest to shortest edge).
    /// Value of 1.0 indicates all edges have similar lengths.
    /// Values < 1.0 indicate high aspect ratio (long thin elements).
    /// 
    /// # Arguments
    /// * `d_matrix` - The D (diagonal) matrix from Jacobian factorization
    /// * `inverse_d_matrix` - The inverse of D
    /// * `element_dim` - Dimension of the element
    /// 
    /// # Returns
    /// * Length ratio quality metric in range [0, 1], where 1 is ideal (uniform edge lengths)
    fn calculate_length_ratio(
        d_matrix: &Vec<Vec<f64>>,
        inverse_d_matrix: &Vec<Vec<f64>>,
        element_dim: usize,
    ) -> f64 {
        let d_norm = Self::frobenius_norm(d_matrix, element_dim);
        let inv_d_norm = Self::frobenius_norm(inverse_d_matrix, element_dim);

        // Length ratio: n / (||D|| * ||D^-1||)
        element_dim as f64 / (d_norm * inv_d_norm)
    }

    /// Calculate orientation metric using Knupp's formula.
    /// 
    /// Orientation metric = 1 + (trace(R) - n) / 4
    /// where R is the rotation matrix from J = R * Q * D factorization.
    /// 
    /// This metric measures element orientation relative to coordinate axes.
    /// Value of 1.0 indicates optimal orientation (aligned with axes).
    /// Values < 1.0 indicate rotation from ideal orientation.
    /// 
    /// # Arguments
    /// * `r_matrix` - The R (rotation/reflection) matrix from Jacobian factorization
    /// * `element_dim` - Dimension of the element
    /// 
    /// # Returns
    /// * Orientation quality metric in range [0, 1], where 1 is ideal
    fn calculate_orientation_metric(
        r_matrix: &Vec<Vec<f64>>,
        element_dim: usize
    ) -> f64 {
        // Orientation metric: 1 + (trace(R) - n) / 4
        // trace(R) for identity matrix = n, giving metric = 1.0
        1.0 + ((Self::trace(r_matrix, element_dim) - element_dim as f64) / 4.0)
    }

    /// Get evaluation points in reference coordinates for quality assessment.
    /// 
    /// Returns strategic sampling points (corners, center, etc.) where Jacobian determinant
    /// will be evaluated to find minimum and maximum values. These points are in the
    /// reference element coordinate system (ξ, η, ζ).
    /// 
    /// # Arguments
    /// * `element_type` - The type of element (Line, Triangle, Quad, Tetra, Hex, etc.)
    /// 
    /// # Returns
    /// * Vector of points in reference coordinates, each point is a vector of coordinates
    fn get_evaluation_points(element_type: &ElementType) -> Vec<Vec<f64>> {
        match element_type {
            // 1D elements - check endpoints and center
            ElementType::Line | ElementType::QuadraticEdge => {
                vec![
                    vec![0.0],    // Left endpoint
                    vec![1.0],    // Right endpoint
                    vec![0.5],    // Center point
                ]
            }
            // 2D triangular elements - check vertices and centroid
            ElementType::Triangle | ElementType::QuadraticTriangle => {
                vec![
                    vec![0.0, 0.0],         // Vertex 1
                    vec![1.0, 0.0],         // Vertex 2
                    vec![0.0, 1.0],         // Vertex 3
                    vec![1.0/3.0, 1.0/3.0], // Centroid (geometric center)
                ]
            }
            // 2D quadrilateral elements - check four corners and center
            ElementType::Quad | ElementType::QuadraticQuad | ElementType::BiquadraticQuad => {
                vec![
                    vec![0.0, 0.0],    // Bottom-left corner
                    vec![1.0, 0.0],    // Bottom-right corner
                    vec![1.0, 1.0],    // Top-right corner
                    vec![0.0, 1.0],    // Top-left corner
                    vec![0.5, 0.5],    // Center point
                ]
            }
            // 3D tetrahedral elements - check four vertices and centroid
            ElementType::Tetra | ElementType::QuadraticTetra => {
                vec![
                    vec![0.0, 0.0, 0.0],        // Vertex 1
                    vec![1.0, 0.0, 0.0],        // Vertex 2
                    vec![0.0, 1.0, 0.0],        // Vertex 3
                    vec![0.0, 0.0, 1.0],        // Vertex 4 (apex)
                    vec![0.25, 0.25, 0.25],     // Centroid
                ]
            }
            // Pyramid element - check base corners, apex, and interior
            ElementType::Pyramid => {
                vec![
                    vec![0.0, 0.0, 0.0],    // Base corner 1
                    vec![1.0, 0.0, 0.0],    // Base corner 2
                    vec![1.0, 1.0, 0.0],    // Base corner 3
                    vec![0.0, 1.0, 0.0],    // Base corner 4
                    vec![0.0, 0.0, 1.0],    // Apex
                ]
            }
            // Wedge/prism elements (triangular prism) - check key vertices
            ElementType::Wedge | ElementType::QuadraticWedge | ElementType::BiquadraticQuadraticWedge => {
                vec![
                    vec![0.0, 0.0, 0.0],       // Bottom triangle vertex 1
                    vec![1.0, 0.0, 0.0],       // Bottom triangle vertex 2
                    vec![0.0, 1.0, 0.0],       // Bottom triangle vertex 3
                    vec![0.0, 0.0, 1.0],       // Top triangle vertex 1
                    vec![1.0, 0.0, 1.0],       // Top triangle vertex 2
                    vec![0.0, 1.0, 1.0],       // Top triangle vertex 3
                    vec![0.5, 0.5, 0.5],       // Centroid
                ]
            }
            // 3D hexahedral elements - check eight corners and center
            ElementType::Hexahedron | ElementType::QuadraticHexahedron | ElementType::BiquadraticQuadraticHexahedron | ElementType::TriquadraticHexahedron => {
                vec![
                    vec![0.0, 0.0, 0.0],        // Bottom face corners
                    vec![1.0, 0.0, 0.0],
                    vec![1.0, 1.0, 0.0],
                    vec![0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 1.0],        // Top face corners
                    vec![1.0, 0.0, 1.0],
                    vec![1.0, 1.0, 1.0],
                    vec![0.0, 1.0, 1.0],
                    vec![0.5, 0.5, 0.5],        // Center point
                ]
            }
            // Vertex has no geometric extent, return dummy point
            ElementType::Vertex => vec![vec![0.0]],
            // Default case for any unhandled element types
            _ => vec![vec![0.0]],
        }
    }

    /// Convert element's node IDs to concrete Node structs.
    /// 
    /// Elements store references to nodes by ID. This function looks up each node ID
    /// in the mesh's node list and returns the actual Node structs with coordinates.
    /// 
    /// # Arguments
    /// * `element` - The element whose nodes to retrieve
    /// * `nodes` - The complete list of nodes in the mesh
    /// 
    /// # Returns
    /// * `Ok(Vec<Node>)` - Vector of Node structs corresponding to the element's node IDs
    /// * `Err(ElementError)` - If any node ID is not found or element has no nodes
    pub fn get_element_nodes(
        element: &Element,
        nodes: &[Node],
    ) -> Result<Vec<Node>, ElementError> {
        // Pre-allocate vector with exact capacity needed
        let mut element_nodes = Vec::with_capacity(element.nodes.len());

        // Look up each node by ID and clone it into our result vector
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

    /// Calculate the determinant of a matrix.
    /// 
    /// Supports both square and non-square matrices. For non-square matrices (element embedded
    /// in higher dimensional space), computes det(J * J^T)^(1/2) which represents the volume/area
    /// scaling factor for the element manifold.
    /// 
    /// # Arguments
    /// * `matrix` - The matrix (Jacobian) to compute determinant of
    /// * `element_dim` - Row dimension (element dimension)
    /// * `mesh_dim` - Column dimension (mesh space dimension)
    /// 
    /// # Returns
    /// * `Ok(f64)` - The determinant value
    /// * `Err(ElementError)` - If matrix dimensions are unsupported
    pub fn calculate_determinant(
        matrix: &Vec<Vec<f64>>, 
        element_dim: usize, 
        mesh_dim: usize
    ) -> Result<f64, ElementError> {
        match (element_dim, mesh_dim) {
            // 1×1 matrix: determinant is the single element
            (1, 1) => Ok(matrix[0][0]),
            
            // 2×2 matrix: det = ad - bc
            (2, 2) => Ok(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]),
            
            // 3×3 matrix: Sarrus rule / cofactor expansion
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
                
                // det = a(ei - fh) - b(di - fg) + c(dh - eg)
                Ok(a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e))
            },
            
            // 1D element in 2D or 3D space (curve in plane/space)
            // Compute metric tensor G = J·J^T, then det(J) = sqrt(det(G)) = ||J||
            (1, 2) | (1, 3) => {
                // Metric tensor is 1×1: G = [J·J^T] = [Σ J[0,k]²]
                let mut metric_tensor = vec![vec![0.0; 1]; 1];
                
                for k in 0..mesh_dim {
                    metric_tensor[0][0] += matrix[0][k] * matrix[0][k];
                }
                    
                // For 1×1 metric tensor, det(J) = sqrt(G[0,0]) = curve length scaling
                Ok((metric_tensor[0][0]).sqrt())
            },
            
            // 2D element in 3D space (surface in space)
            // Compute metric tensor G = J·J^T, then det(J) = sqrt(det(G))
            (2, 3) => {
                // Metric tensor is 2×2: G = J·J^T
                let mut metric_tensor = vec![vec![0.0; 2]; 2];
                
                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..mesh_dim {
                            metric_tensor[i][j] += matrix[i][k] * matrix[j][k];
                        }
                    }
                }
                
                // det(J) = sqrt(det(G)) = sqrt(g_00*g_11 - g_01*g_10) = surface area scaling
                Ok((metric_tensor[0][0] * metric_tensor[1][1] - metric_tensor[0][1] * metric_tensor[1][0]).sqrt())
            },
            
            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Determinant calculation not implemented for {}×{} (element_dim × mesh_dim) matrices",
                element_dim, mesh_dim
            ))),
        }
    }

    /// Calculate the inverse of a matrix.
    /// 
    /// Supports both square and non-square matrices. For non-square matrices, computes
    /// the Moore-Penrose pseudoinverse using the metric tensor approach.
    /// 
    /// # Arguments
    /// * `matrix` - The matrix to invert
    /// * `element_dim` - Row dimension
    /// * `mesh_dim` - Column dimension
    /// 
    /// # Returns
    /// * `Ok(Vec<Vec<f64>>)` - The inverse (or pseudoinverse) matrix
    /// * `Err(ElementError)` - If matrix is singular or dimensions are unsupported
    pub fn inverse_matrix(
        matrix: &Vec<Vec<f64>>, 
        element_dim: usize, 
        mesh_dim: usize
    ) -> Result<Vec<Vec<f64>>, ElementError> {
        match (element_dim, mesh_dim) {
            // 1×1 matrix: inverse is simply 1/a
            (1, 1) => {
                if matrix[0][0].abs() < 1e-12 {
                    return Err(ElementError::GeometryError("Matrix is singular (determinant ≈ 0)".to_string()));
                }
                Ok(vec![vec![1.0 / matrix[0][0]]])
            },

            // 2×2 and 3×3 matrices: use adjoint method (A^-1 = adj(A) / det(A))
            (2, 2) | (3, 3) => {
                // Calculate adjoint (adjugate) matrix
                let adjoint = Self::calculate_adjoint(&matrix, element_dim, mesh_dim)?;
                let determinant = Self::calculate_determinant(&matrix, element_dim, mesh_dim)?;

                // Check for singular matrix
                if determinant.abs() < 1e-12 {
                    return Err(ElementError::GeometryError("Matrix is singular (determinant ≈ 0)".to_string()));
                }
                
                // Compute inverse: A^-1 = adj(A) / det(A)
                let inv_det = 1.0 / determinant;
                let mut inverse = vec![vec![0.0; mesh_dim]; element_dim];
                
                for i in 0..element_dim {
                    for j in 0..mesh_dim {
                        inverse[i][j] = adjoint[i][j] * inv_det;
                    }
                }
                Ok(inverse)
            },

            // 1D element in 2D or 3D space: pseudoinverse is J^T / (J·J^T)
            (1, 2) | (1, 3) => {
                // Compute norm squared: ||J||² = J·J^T
                let mut sum = 0.0;
                for k in 0..mesh_dim {
                    sum += matrix[0][k] * matrix[0][k];
                }
                if sum.abs() < 1e-12 {
                    return Err(ElementError::GeometryError("Matrix is singular (norm ≈ 0)".to_string()));
                }
                // Pseudoinverse: J^+ = J^T / ||J||²
                let inv_sum = 1.0 / sum;
                let mut inverse = vec![vec![0.0; element_dim]; mesh_dim];
                for j in 0..mesh_dim {
                    inverse[j][0] = matrix[0][j] * inv_sum;
                }
                Ok(inverse)
            },

            // 2D element in 3D space: pseudoinverse using metric tensor
            // J^+ = J^T * (J * J^T)^-1
            (2, 3) => {
                // Compute metric tensor: G = J * J^T (2×2 matrix)
                let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
                for i in 0..element_dim {
                    for j in 0..element_dim {
                        let mut sum = 0.0;
                        for k in 0..mesh_dim {
                            sum += matrix[i][k] * matrix[j][k];
                        }
                        square_matrix[i][j] = sum;
                    }
                }

                // Compute inverse of metric tensor: G^-1
                let adjoint = Self::calculate_adjoint(&square_matrix, element_dim, element_dim)?;
                let determinant = Self::calculate_determinant(&square_matrix, element_dim, element_dim)?;
                
                if determinant.abs() < 1e-12 {
                    return Err(ElementError::GeometryError("Matrix is singular (metric tensor determinant ≈ 0)".to_string()));
                }
                
                // Compute G^-1
                let inv_det = 1.0 / determinant;
                let mut inv_square = vec![vec![0.0; element_dim]; element_dim];
                
                for i in 0..element_dim {
                    for j in 0..element_dim {
                        inv_square[i][j] = adjoint[i][j] * inv_det;
                    }
                }
                
                // Compute pseudoinverse: J^+ = J^T * (J * J^T)^-1
                // Result is [mesh_dim × element_dim] = [3 × 2]
                let mut inverse = vec![vec![0.0; element_dim]; mesh_dim];
                for i in 0..mesh_dim {
                    for j in 0..element_dim {
                        let mut sum = 0.0;
                        for k in 0..element_dim {
                            // J^T[i][k] = J[k][i]
                            sum += matrix[k][i] * inv_square[j][k];
                        }
                        inverse[i][j] = sum;
                    }
                }
                Ok(inverse)
            },

            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian inverse calculation not implemented for {}×{} (element_dim × mesh_dim) matrices",
                element_dim, mesh_dim
            ))),
        }
    }

    /// Calculate the adjoint (adjugate) matrix.
    /// 
    /// The adjoint matrix is the transpose of the cofactor matrix. Used in computing
    /// matrix inverse: A^-1 = adj(A) / det(A).
    /// 
    /// # Arguments
    /// * `matrix` - The matrix to compute adjoint of
    /// * `element_dim` - Row dimension
    /// * `mesh_dim` - Column dimension
    /// 
    /// # Returns
    /// * `Ok(Vec<Vec<f64>>)` - The adjoint matrix
    /// * `Err(ElementError)` - If dimensions are unsupported
    pub fn calculate_adjoint(
        matrix: &Vec<Vec<f64>>, 
        element_dim: usize, 
        mesh_dim: usize
    ) -> Result<Vec<Vec<f64>>, ElementError> {
        match (element_dim, mesh_dim) {
            // 2×2 matrix: adj = [[d, -b], [-c, a]] for matrix [[a, b], [c, d]]
            (2, 2) => {
                Ok(vec![
                    vec![matrix[1][1], -matrix[0][1]],
                    vec![-matrix[1][0], matrix[0][0]],
                ])
            },
            // 3×3 matrix: compute cofactor matrix and transpose
            (3, 3) => {
                // Extract matrix elements for readability
                let a = matrix[0][0];
                let b = matrix[0][1];
                let c = matrix[0][2];
                let d = matrix[1][0];
                let e = matrix[1][1];
                let f = matrix[1][2];
                let g = matrix[2][0];
                let h = matrix[2][1];
                let i = matrix[2][2];

                // Compute cofactor matrix entries and transpose to get adjoint
                Ok(vec![
                    vec![
                        e * i - f * h,  // Cofactor(0,0)
                        c * h - b * i,  // Cofactor(1,0) transposed
                        b * f - c * e,  // Cofactor(2,0) transposed
                    ],
                    vec![
                        f * g - d * i,  // Cofactor(0,1) transposed
                        a * i - c * g,  // Cofactor(1,1)
                        c * d - a * f,  // Cofactor(2,1) transposed
                    ],
                    vec![
                        d * h - e * g,  // Cofactor(0,2) transposed
                        b * g - a * h,  // Cofactor(1,2) transposed
                        a * e - b * d,  // Cofactor(2,2)
                    ],
                ])
            },
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian adjoint calculation not implemented for {}×{} (element_dim × mesh_dim) matrices",
                element_dim, mesh_dim
            ))),
        }
    }

    /// Calculate the trace of a square matrix.
    /// 
    /// The trace is the sum of diagonal elements: tr(A) = Σ A[i,i].
    /// Used in the orientation metric calculation.
    /// 
    /// # Arguments
    /// * `matrix` - The square matrix
    /// * `element_dim` - Dimension of the matrix
    /// 
    /// # Returns
    /// * The trace value
    fn trace(matrix: &Vec<Vec<f64>>, element_dim: usize) -> f64 {
        let mut sum = 0.0;
        // Sum diagonal elements
        for i in 0..element_dim {
            sum += matrix[i][i];
        }
        sum
    }

    /// Calculate the Frobenius norm of a matrix.
    /// 
    /// Frobenius norm: ||A||_F = sqrt(Σ_ij A[i,j]²).
    /// This is the matrix equivalent of the Euclidean (L2) norm for vectors.
    /// Used extensively in Knupp's quality metrics.
    /// 
    /// # Arguments
    /// * `matrix` - The matrix
    /// * `element_dim` - Dimension of the square matrix
    /// 
    /// # Returns
    /// * The Frobenius norm
    fn frobenius_norm(matrix: &Vec<Vec<f64>>, element_dim: usize) -> f64 {
        let mut sum = 0.0;
        // Sum of squares of all matrix elements
        for i in 0..element_dim {
            for j in 0..element_dim {
                sum += matrix[i][j] * matrix[i][j];
            }
        }
        sum.sqrt()
    }

    /// Compute 3D cross product of two vectors.
    /// 
    /// Cross product: a × b = [a_y*b_z - a_z*b_y, a_z*b_x - a_x*b_z, a_x*b_y - a_y*b_x]
    /// Used in Jacobian factorization for 3D elements.
    /// 
    /// # Arguments
    /// * `a` - First vector (length 3)
    /// * `b` - Second vector (length 3)
    /// 
    /// # Returns
    /// * Cross product vector perpendicular to both inputs
    fn cross_product_3d(a: &[f64], b: &[f64]) -> Vec<f64> {
        vec![
            a[1] * b[2] - a[2] * b[1],  // x-component
            a[2] * b[0] - a[0] * b[2],  // y-component
            a[0] * b[1] - a[1] * b[0],  // z-component
        ]
    }

    /// Compute Euclidean (L2) norm of a vector.
    /// 
    /// ||v|| = sqrt(Σ v_i²)
    /// 
    /// # Arguments
    /// * `v` - The vector
    /// 
    /// # Returns
    /// * The Euclidean norm (magnitude) of the vector
    fn vector_norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    /// Perform Knupp's Jacobian factorization: J = R * Q * D.
    /// 
    /// This decomposes the Jacobian matrix into three components that isolate different
    /// geometric distortions:
    /// - R: Rotation/reflection matrix (measures orientation)
    /// - Q: Upper triangular matrix with unit diagonal (measures skewness/angles)
    /// - D: Diagonal scaling matrix (measures aspect ratio/edge lengths)
    /// 
    /// The factorization allows us to separately measure different types of element quality:
    /// orientation, skewness, and length ratio.
    /// 
    /// # Arguments
    /// * `metric_jacobian` - The (square) Jacobian or metric tensor matrix
    /// * `element_dim` - Dimension of the element (1, 2, or 3)
    /// 
    /// # Returns
    /// * `Ok((R, Q, D))` - The three factorization matrices
    /// * `Err(ElementError)` - If factorization fails or dimensions are unsupported
    fn jacobian_factorization(metric_jacobian: &Vec<Vec<f64>>, 
                            element_dim: usize,  
                            )-> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>), ElementError> {

        // Compute the metric tensor λ = J * J^T
        // This is used in the factorization formulas
        let mut lambda = vec![vec![0.0; element_dim]; element_dim];
        for i in 0..element_dim {
            for j in 0..element_dim {
                for k in 0..element_dim {
                    lambda[i][j] += metric_jacobian[i][k] * metric_jacobian[j][k];
                }
            }
        }

        // Calculate determinant for use in factorization
        let det_metric_jacobian = Self::calculate_determinant(&metric_jacobian, element_dim, element_dim)?;
                

        // Initialize result matrices
        let mut r_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; element_dim];
        let mut q_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; element_dim];
        let mut d_matrix: Vec<Vec<f64>> = vec![vec![0.0; element_dim]; element_dim];

        match (element_dim, element_dim) {
            // 1×1 matrix: trivial factorization
            (1, 1) => {
                // For 1D, R and Q are identity, D contains the Jacobian value
                let r_matrix = vec![vec![1.0]];
                let q_matrix = vec![vec![1.0]];
                let d_matrix = vec![vec![metric_jacobian[0][0]]];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 2×2 matrix: 2D factorization formulas
            (2, 2) => { 
                // R matrix: rotation component (orthogonal matrix)
                r_matrix = vec![
                    vec![metric_jacobian[0][0] / lambda[0][0].sqrt(), - metric_jacobian[1][0] / lambda[0][0].sqrt()],
                    vec![metric_jacobian[1][0] / lambda[0][0].sqrt(), metric_jacobian[0][0] / lambda[0][0].sqrt()],
                ];
                // Q matrix: skew component (upper triangular with unit diagonal)
                q_matrix = vec![
                    vec![1.0, lambda[0][1] / (lambda[0][0] * lambda[1][1]).sqrt()],
                    vec![0.0, det_metric_jacobian / (lambda[0][0] * lambda[1][1]).sqrt()],
                ];
                // D matrix: diagonal scaling component
                d_matrix = vec![
                    vec![1.0, 0.0],
                    vec![0.0, (lambda[1][1] / lambda[0][0]).sqrt()],
                ];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // 3×3 matrix: 3D factorization formulas
            (3, 3) => {
                // Extract Jacobian columns for cross product computation
                let col_0 = vec![metric_jacobian[0][0], metric_jacobian[1][0], metric_jacobian[2][0]];
                let col_1 = vec![metric_jacobian[0][1], metric_jacobian[1][1], metric_jacobian[2][1]];

                // Compute cross product of first two columns to get orthogonal direction
                let cross_product_01 = Self::cross_product_3d(&col_0, &col_1);
                let norm_cross_product_01 = Self::vector_norm(&cross_product_01);

                // R matrix: 3D rotation component
                // Each column is normalized and orthogonalized
                r_matrix = vec![
                    vec![col_0[0] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[0] - lambda[0][1] * col_0[0]) / (lambda[0][0].sqrt() * norm_cross_product_01), cross_product_01[0] / norm_cross_product_01],
                    vec![col_0[1] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[1] - lambda[0][1] * col_0[1]) / (lambda[0][0].sqrt() * norm_cross_product_01), cross_product_01[1] / norm_cross_product_01],
                    vec![col_0[2] / lambda[0][0].sqrt(), (lambda[0][0] * col_1[2] - lambda[0][1] * col_0[2]) / (lambda[0][0].sqrt() * norm_cross_product_01), cross_product_01[2] / norm_cross_product_01],
                ];
                // Q matrix: 3D skew component (upper triangular)
                q_matrix = vec![
                    vec![1.0, lambda[0][1] / (lambda[0][0] * lambda[1][1]).sqrt(), lambda[0][2] / (lambda[0][0] * lambda[2][2]).sqrt()],
                    vec![0.0, norm_cross_product_01 / (lambda[0][0] * lambda[1][1]).sqrt(), (lambda[0][0] * lambda[1][2] - lambda[0][1] * lambda[0][2]) / ((lambda[0][0] * lambda[2][2]).sqrt() * norm_cross_product_01)],
                    vec![0.0, 0.0, det_metric_jacobian / (lambda[2][2].sqrt() * norm_cross_product_01)],
                ];
                // D matrix: 3D diagonal scaling component
                d_matrix = vec![
                    vec![1.0, 0.0, 0.0],
                    vec![0.0, (lambda[1][1]).sqrt() / (lambda[0][0]).sqrt(), 0.0],
                    vec![0.0, 0.0, (lambda[2][2]).sqrt() / (lambda[0][0]).sqrt()],
                ];

                Ok((r_matrix, q_matrix, d_matrix))
            }

            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian factorization not implemented for {}×{} (element_dim × element_dim) matrices",
                element_dim, element_dim
            ))),
        }
    }
}

// ============================================================================
// UNIT TESTS
// ============================================================================
// The following tests verify geometric quality metrics for various element types
// in different distortion scenarios: perfect reference, bad shape, bad skewness,
// bad length ratio, bad orientation, bad Jacobian ratio, and inverted elements.
 
#[cfg(test)]
mod geometric_quality_tests {
    use super::*;

    #[test]
    /// Test 1D quadratic edge element in 1D mesh with perfect reference configuration
    fn test_1d_quadratic_edge_in_1d_mesh_perfect_reference() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![1.0] },  
            Node { id: 2, coordinates: vec![0.5] },  
            
        ];
        let element_type = ElementType::QuadraticEdge;
        let mesh_dim = 1;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 1.0).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10); 
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 1.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 1.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 2: 1D quadratic edge element in 1D mesh - BAD Jacobian Ratio
    fn test_1d_quadratic_edge_in_1d_mesh_bad_jacobianratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![1.0] },  
            Node { id: 2, coordinates: vec![0.26] },  
            
        ];
        let element_type = ElementType::QuadraticEdge;
        let mesh_dim = 1;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.04).abs() < 1e-10);
        assert!((max_detj - 1.96).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 0.04).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 0.04).abs() < 1e-10);  

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 25.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((inverse_d_matrix[0][0] - 25.0).abs() < 1e-10); 
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 1.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.04).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.0204081632653061).abs() < 1e-10);
    }

    #[test]
    // Test 3: 1D quadratic edge element in 1D mesh - Inverted element
    fn test_1d_quadratic_edge_in_1d_mesh_inverted() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![4.0] },  
            Node { id: 2, coordinates: vec![1.0] },  
            
        ];
        let element_type = ElementType::QuadraticEdge;
        let mesh_dim = 1;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!(has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.0).abs() < 1e-10);
    }

    #[test]
    // Test 4: Perfect 2D quad element
    fn test_2d_quad_perfect_reference() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0] },  
            Node { id: 2, coordinates: vec![1.0, 1.0] },  
            Node { id: 3, coordinates: vec![0.0, 1.0] },
            
        ];
        let element_type = ElementType::Quad;
        let mesh_dim = 2;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 1.0).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 1.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 1.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 5: 2D quad element - BAD Skewness
    fn test_2d_quad_bad_shape_skewness() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0] },  
            Node { id: 2, coordinates: vec![3.0, 0.5] },  
            Node { id: 3, coordinates: vec![2.0, 1.0] },
            
        ];
        let element_type = ElementType::Quad;
        let mesh_dim = 2;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.5).abs() < 1e-10);
        assert!((max_detj - 2.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }
        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][0] - 2.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][1] - 0.5).abs() < 1e-10);
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - (-2.0)).abs() < 1e-10);
        assert!((r_matrix[1][0] - 2.0).abs() < 1e-10);  assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.9701425001453319).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 0.242535625036333).abs() < 1e-10);
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 2.06155281280883).abs() < 1e-10);

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][0] - (-4.0)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][1] - 2.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - (-4.0)).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 4.123105625617661).abs() < 1e-10);
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 0.4850712500726659).abs() < 1e-10);
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.1904761904761905).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.2425356250363330).abs() < 1e-10);
        assert!((quality.length_ratio - 0.7853534524986019).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.5).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.25).abs() < 1e-10);
    }

    #[test]
    // Test 6: 2D quad element - BAD Length Ratio
    fn test_2d_quad_bad_shape_lengthratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0] },
            Node { id: 1, coordinates: vec![10.0, 0.0] },  
            Node { id: 2, coordinates: vec![10.0, 1.0] },  
            Node { id: 3, coordinates: vec![0.0, 1.0] },
            
        ];
        let element_type = ElementType::Quad;
        let mesh_dim = 2;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 10.0).abs() < 1e-10);
        assert!((max_detj - 10.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }
        assert!((evaluated_jacobian[0][0] - 10.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 0.1).abs() < 1e-10);

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 0.1).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 10.0).abs() < 1e-10);
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.198019801980198).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.9999999999999998).abs() < 1e-10);
        assert!((quality.length_ratio - 0.198019801980198).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 10.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 7: 2D quad element - BAD Orientation
    fn test_2d_quad_bad_orientation() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0] },
            Node { id: 1, coordinates: vec![-0.707, 0.707] },  
            Node { id: 2, coordinates: vec![-1.414, 0.0] },  
            Node { id: 3, coordinates: vec![-0.707, -0.707] },
            
        ];
        let element_type = ElementType::Quad;
        let mesh_dim = 2;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.999698).abs() < 1e-10);
        assert!((max_detj - 0.999698).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }
        assert!((evaluated_jacobian[0][0] - (-0.707)).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.707).abs() < 1e-10);
        assert!((evaluated_jacobian[1][0] - (-0.707)).abs() < 1e-10);  assert!((evaluated_jacobian[1][1] - (-0.707)).abs() < 1e-10);
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - (-0.7071067811865475)).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.7071067811865475).abs() < 1e-10);
        assert!((r_matrix[1][0] - (-0.7071067811865475)).abs() < 1e-10);  assert!((r_matrix[1][1] - (-0.7071067811865475)).abs() < 1e-10);
        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - (-0.7072135785007072)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - (-0.7072135785007072)).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][0] - 0.7072135785007072).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][1] - (-0.7072135785007072)).abs() < 1e-10);
        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.9999999999999998).abs() < 1e-10);
        assert!((quality.length_ratio - 0.9999999999999998).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.1464466094067263).abs() < 1e-10);
        assert!((quality.volume_metric - 0.9996979999999999).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 8: 2D quad element - BAD Jacobian Ratio
    fn test_2d_quad_bad_jacobianratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-1.5, 0.0] },
            Node { id: 1, coordinates: vec![1.5, 0.0] },  
            Node { id: 2, coordinates: vec![0.05, 1.0] },  
            Node { id: 3, coordinates: vec![-0.04, 1.0] },
            
        ];
        let element_type = ElementType::Quad;
        let mesh_dim = 2;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.09).abs() < 1e-10);
        assert!((max_detj - 3.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }
        assert!((evaluated_jacobian[0][0] - 0.09).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][0] - (-1.45)).abs() < 1e-10);  assert!((evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 16.11111111111111).abs() < 1e-10);
        assert!((r_matrix[1][0] - (-16.11111111111111)).abs() < 1e-10);  assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);
        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - (-0.8232127859153063)).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 0.5677329558036595).abs() < 1e-10);
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 19.57101661534282).abs() < 1e-10);

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 11.11111111111111).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][0] - 16.11111111111111).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 1.45).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 1.761391495380854).abs() < 1e-10);
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 0.05109596602232936).abs() < 1e-10);
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0578666495209927).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.5677329558036593).abs() < 1e-10);
        assert!((quality.length_ratio - 0.1019258243324611).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.09).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.03).abs() < 1e-10);
    }

    #[test]
    // Test 9: 2D quad element - Inverted element
    fn test_2d_quad_inverted() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0] },  
            Node { id: 2, coordinates: vec![0.0, 1.0] },  
            Node { id: 3, coordinates: vec![1.0, 1.0] },
            
        ];
        let element_type = ElementType::Quad;
        let mesh_dim = 2;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!(has_nonpositive_detj);

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.0).abs() < 1e-10);
    }

    #[test]
    // Test 10: Perfect 3D wedge element
    fn test_3d_wedge_perfect_reference() { 
        // Perfect square element
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![0.0, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![0.0, 0.0, 1.0] },
            Node { id: 4, coordinates: vec![1.0, 0.0, 1.0] },
            Node { id: 5, coordinates: vec![0.0, 1.0, 1.0] },
            
        ];
        let element_type = ElementType::Wedge;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        println!("min_point: {:?}", min_point);
        assert!((min_detj - 1.0).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][0] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((r_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][2] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][2] - 1.0).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][2] - 1.0).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][2] - 1.0).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][2] - 1.0).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 1.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 1.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 11: 3D wedge element - BAD Skewness
    fn test_3d_wedge_bad_shape_skewness() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![0.0, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![3.0, 0.0, 1.0] },
            Node { id: 4, coordinates: vec![4.0, 0.0, 1.0] },
            Node { id: 5, coordinates: vec![3.0, 1.0, 1.0] },
            
        ];
        let element_type = ElementType::Wedge;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 1.0).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][0] - 3.0).abs() < 1e-10);  assert!((evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[0][2] - (-0.9486832980505138)).abs() < 1e-10);  assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][1] - 0.3162277660168379).abs() < 1e-10);  assert!((r_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][0] - 3.0).abs() < 1e-10);  assert!((r_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][2] - 0.3162277660168379).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[0][2] - 0.9486832980505138).abs() < 1e-10);  assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][1] - 3.162277660168379).abs() < 1e-10);  assert!((q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][2] - 0.1).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][2] - 3.162277660168379).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][0] - (-3.0)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][2] - (-9.486832980505138)).abs() < 1e-10);  assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][1] - 0.3162277660168379).abs() < 1e-10);  assert!((inverse_q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][2] - 10.0).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][2] - 0.3162277660168379).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.25).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0628832574081254).abs() < 1e-10);
        assert!((quality.length_ratio - 0.5976143046671967).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.658113883008419).abs() < 1e-10);
        assert!((quality.volume_metric - 1.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 12: 3D wedge element - BAD Length Ratio
    fn test_3d_wedge_bad_shape_lengthratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![0.0, 10.0, 0.0] },  
            Node { id: 3, coordinates: vec![0.0, 0.0, 1.0] },
            Node { id: 4, coordinates: vec![1.0, 0.0, 1.0] },
            Node { id: 5, coordinates: vec![0.0, 10.0, 1.0] },
            
        ];
        let element_type = ElementType::Wedge;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 10.0).abs() < 1e-10);
        assert!((max_detj - 10.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][1] - 10.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][0] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((r_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][2] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][2] - 1.0).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][1] - 10.0).abs() < 1e-10);  assert!((d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][2] - 1.0).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][1] - 0.1).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][2] - 1.0).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][1] - 0.1).abs() < 1e-10);  assert!((inverse_d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][2] - 1.0).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.2095188685475281).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.2095188685475281).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 10.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 13: 3D wedge element - BAD Orientation
    fn test_3d_wedge_bad_orientation() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![-0.5, 0.707, -0.5] },  
            Node { id: 2, coordinates: vec![-0.5, -0.707, -0.5] },  
            Node { id: 3, coordinates: vec![-0.707, 0.0, 0.707] },
            Node { id: 4, coordinates: vec![-1.207, 0.707, 0.207] },
            Node { id: 5, coordinates: vec![-1.207, -0.707, 0.207] },
            
        ];
        let element_type = ElementType::Wedge;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.999698).abs() < 1e-10);
        assert!((max_detj - 0.999698).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - (-0.5)).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.707).abs() < 1e-10);
        assert!((evaluated_jacobian[0][2] - (-0.5)).abs() < 1e-10);  assert!((evaluated_jacobian[1][0] - (-0.5)).abs() < 1e-10);
        assert!((evaluated_jacobian[1][1] - (-0.707)).abs() < 1e-10);  assert!((evaluated_jacobian[1][2] - (-0.5)).abs() < 1e-10);
        assert!((evaluated_jacobian[2][0] - (-0.707)).abs() < 1e-10);  assert!((evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][2] - 0.707).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - (-0.5000377542757255)).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.7071823039935743).abs() < 1e-10);
        assert!((r_matrix[0][2] - (-0.4999622428732884)).abs() < 1e-10);  assert!((r_matrix[1][0] - (-0.5000377542757255)).abs() < 1e-10);
        assert!((r_matrix[1][1] - (-0.7070312583795209)).abs() < 1e-10);  assert!((r_matrix[1][2] - (-0.4999622428732884)).abs() < 1e-10);
        assert!((r_matrix[2][0] - (-0.7070533845458759)).abs() < 1e-10);  assert!((r_matrix[2][1] - 0.0001067892491358361).abs() < 1e-10);
        assert!((r_matrix[2][2] - 0.7071601737953161).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0001510228044435387).abs() < 1e-10);
        assert!((q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][1] - 0.999924485746577).abs() < 1e-10);  assert!((q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][2] - 1.000075508551451).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][2] - 0.999924485746577).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - (-0.5)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - (-0.5)).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[0][2] - (-0.7072135785007072)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][0] - 0.7072135785007072).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][1] - (-0.7072135785007072)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][0] - (-0.5)).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[2][1] - (-0.5)).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][2] - 0.7072135785007072).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - (-0.000151034209679124)).abs() < 1e-10);
        assert!((inverse_q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][1] - 1.000075519956256).abs() < 1e-10);  assert!((inverse_q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][2] - 0.9999244971496598).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][2] - 1.000075519956256).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.9999999898631611).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.9999999847941676).abs() < 1e-10);
        assert!((quality.length_ratio - 0.9999999974654075).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.1250227902850174).abs() < 1e-10);
        assert!((quality.volume_metric - 0.9996979999999999).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 14: 3D wedge element - BAD Jacobian Ratio
    fn test_3d_wedge_bad_jacobianratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![0.0, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![0.0, 0.0, 1.0] },
            Node { id: 4, coordinates: vec![0.032, 0.0, 1.0] },
            Node { id: 5, coordinates: vec![0.0, 0.032, 1.0] },
            
        ];
        let element_type = ElementType::Wedge;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.001024).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 0.032).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][1] - 0.032).abs() < 1e-10);  assert!((evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][0] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  
        
        let metric_jacobian = evaluated_jacobian.clone();

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((r_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[2][2] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[2][2] - 1.0).abs() < 1e-10);  

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[2][2] - 31.25).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_evaluated_jacobian = GeometricAnalysis::inverse_matrix(&evaluated_jacobian, element_dim, mesh_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_evaluated_jacobian[0][0] - 31.25).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[1][1] - 31.25).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_evaluated_jacobian[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_evaluated_jacobian[2][2] - 1.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[2][2] - 1.0).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[0][2] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][2] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[2][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[2][2] - 0.032).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0677954927313605).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0677954927313605).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.001024).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.001024).abs() < 1e-10);
    }

    #[test]
    // Test 15: 3D wedge element - Inverted element
    fn test_3d_wedge_inverted() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![0.0, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![0.0, 0.0, 1.0] },
            Node { id: 4, coordinates: vec![0.0, 0.0, 1.0] },
            Node { id: 5, coordinates: vec![1.0, 1.0, 1.0] },
            
        ];
        let element_type = ElementType::Wedge;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        println!("min_point: {:?}", min_point);
        assert!(has_nonpositive_detj); 
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.0).abs() < 1e-10);
    }

    #[test]
    // Test 16: Perfect 2D quadratic triangle element in 3D mesh
    fn test_2d_quadratic_triangle_in_3d_mesh_perfect_reference() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-0.5, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![0.5, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![-0.5, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 4, coordinates: vec![0.0, 0.5, 0.0] },
            Node { id: 5, coordinates: vec![-0.5, 0.5, 0.0] },
            
        ];
        let element_type = ElementType::QuadraticTriangle;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 1.0).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        assert!((evaluated_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[0][2] - 0.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][0] - 0.0).abs() < 1e-10);
        assert!((evaluated_jacobian[1][1] - 1.0).abs() < 1e-10);  assert!((evaluated_jacobian[1][2] - 0.0).abs() < 1e-10);  
        
        // Construct metric Jacobian matrix only for knupps metrics
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices, use as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices, compute J * J^T to get a square matrix [element_dim × element_dim]
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        assert!((metric_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((metric_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((metric_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((metric_jacobian[1][1] - 1.0).abs() < 1e-10);

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);  
          

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_metric_jacobian = GeometricAnalysis::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_metric_jacobian[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_metric_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_metric_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_metric_jacobian[1][1] - 1.0).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.9999999999999998).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.9999999999999998).abs() < 1e-10);
        assert!((quality.length_ratio - 0.9999999999999998).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 1.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 17: 2D quadratic triangle element in 3D mesh - BAD Skewness 
    fn test_2d_quadratic_triangle_in_3d_mesh_bad_shape_skewness() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.000, 0.000, 0.000] },
            Node { id: 1, coordinates: vec![0.866, 0.000, -0.500] },  
            Node { id: 2, coordinates: vec![1.732, 0.612, -0.854] },  
            Node { id: 3, coordinates: vec![0.433, 0.000, -0.250] },
            Node { id: 4, coordinates: vec![1.299, 0.306, -0.677] },
            Node { id: 5, coordinates: vec![0.866, 0.306, -0.427] },
            
        ];
        let element_type = ElementType::QuadraticTriangle;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.6249108593711586).abs() < 1e-10);
        assert!((max_detj - 0.6249108593711586).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        } 
        
        // Construct metric Jacobian matrix only for knupps metrics
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices, use as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices, compute J * J^T to get a square matrix [element_dim × element_dim]
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        assert!((metric_jacobian[0][0] - 0.999956).abs() < 1e-10);  assert!((metric_jacobian[0][1] - 1.926912).abs() < 1e-10);
        assert!((metric_jacobian[1][0] - 1.926912).abs() < 1e-10);  assert!((metric_jacobian[1][1] - 4.103684).abs() < 1e-10);

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 0.4606135319997494).abs() < 1e-10);  assert!((r_matrix[0][1] - (-0.8876007966077519)).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.8876007966077519).abs() < 1e-10);  assert!((r_matrix[1][1] - 0.4606135319997494).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.9992125096812894).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 0.03967821186015499).abs() < 1e-10);  
          

        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 2.088313304399514).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_metric_jacobian = GeometricAnalysis::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_metric_jacobian[0][0] - 10.50842835555627).abs() < 1e-10);  assert!((inverse_metric_jacobian[0][1] - (-4.934302129370009)).abs() < 1e-10);
        assert!((inverse_metric_jacobian[1][0] - (-4.934302129370009)).abs() < 1e-10);  assert!((inverse_metric_jacobian[1][1] - 2.560617724149477).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - (-25.18290171953798)).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 25.20274864009696).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 0.478855350819855).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0309120413880750).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0396782118601550).abs() < 1e-10);
        assert!((quality.length_ratio - 0.7790684090559280).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.7303067659998748).abs() < 1e-10);
        assert!((quality.volume_metric - 0.6249108593711586).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }
    
    #[test]
    // Test 18: 2D quadratic triangle element in 3D mesh - BAD Length Ratio 
    fn test_2d_quadratic_triangle_in_3d_mesh_bad_lengthratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.000, 0.000, 0.000] },
            Node { id: 1, coordinates: vec![8.660, 0.000, -5.000] },  
            Node { id: 2, coordinates: vec![0.433, 0.612, 0.354] },  
            Node { id: 3, coordinates: vec![4.330, 0.000, -2.500] },
            Node { id: 4, coordinates: vec![4.547, 0.306, -2.323] },
            Node { id: 5, coordinates: vec![0.217, 0.306, 0.177] },
            
        ];
        let element_type = ElementType::QuadraticTriangle;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 8.0441179762606669).abs() < 1e-10);
        assert!((max_detj - 8.0571123621307414).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        
        // Construct metric Jacobian matrix only for knupps metrics
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices, use as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices, compute J * J^T to get a square matrix [element_dim × element_dim]
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        assert!((metric_jacobian[0][0] - 99.9956).abs() < 1e-10);  assert!((metric_jacobian[0][1] - 1.96246).abs() < 1e-10);
        assert!((metric_jacobian[1][0] - 1.96246).abs() < 1e-10);  assert!((metric_jacobian[1][1] - 0.685621).abs() < 1e-10);

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 0.9998074762031988).abs() < 1e-10);  assert!((r_matrix[0][1] - (-0.01962168515164397)).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.01962168515164397).abs() < 1e-10);  assert!((r_matrix[1][1] - 0.9998074762031988).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.95033405197289).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 0.3112317298425536).abs() < 1e-10);  
          
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 0.02078471026206253).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_metric_jacobian = GeometricAnalysis::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_metric_jacobian[0][0] - 0.0105956413226638).abs() < 1e-10);  assert!((inverse_metric_jacobian[0][1] - (-0.03032801251722862)).abs() < 1e-10);
        assert!((inverse_metric_jacobian[1][0] - (-0.03032801251722862)).abs() < 1e-10);  assert!((inverse_metric_jacobian[1][1] - 1.545339934810282).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - (-3.053461330737861)).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 3.213040008825198).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 48.11228962980825).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0129321359214928).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.3112317298425537).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0415514701153219).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.9999037381015994).abs() < 1e-10);
        assert!((quality.volume_metric - 8.0441179762606669).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.9983872155077358).abs() < 1e-10);
    }

    #[test]
    // Test 19: 2D quadratic triangle element in 3D mesh - BAD Orientation (we loose orientation quality because of JJ' ??)
    fn test_2d_quadratic_triangle_in_3d_mesh_bad_orientation() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.000, 0.000, 0.000] },
            Node { id: 1, coordinates: vec![8.660, 0.000, -5.000] },  
            Node { id: 2, coordinates: vec![0.433, 0.612, 0.354] },  
            Node { id: 3, coordinates: vec![4.330, 0.000, -2.500] },
            Node { id: 4, coordinates: vec![4.547, 0.306, -2.323] },
            Node { id: 5, coordinates: vec![0.217, 0.306, 0.177] },
            
        ];
        let element_type = ElementType::QuadraticTriangle;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 8.0441179762606669).abs() < 1e-10);
        assert!((max_detj - 8.0571123621307414).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        
        // Construct metric Jacobian matrix only for knupps metrics
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices, use as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices, compute J * J^T to get a square matrix [element_dim × element_dim]
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        assert!((metric_jacobian[0][0] - 99.9956).abs() < 1e-10);  assert!((metric_jacobian[0][1] - 1.96246).abs() < 1e-10);
        assert!((metric_jacobian[1][0] - 1.96246).abs() < 1e-10);  assert!((metric_jacobian[1][1] - 0.685621).abs() < 1e-10);

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 0.9998074762031988).abs() < 1e-10);  assert!((r_matrix[0][1] - (-0.01962168515164397)).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.01962168515164397).abs() < 1e-10);  assert!((r_matrix[1][1] - 0.9998074762031988).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.95033405197289).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 0.3112317298425536).abs() < 1e-10);  
          
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 0.02078471026206253).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_metric_jacobian = GeometricAnalysis::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_metric_jacobian[0][0] - 0.0105956413226638).abs() < 1e-10);  assert!((inverse_metric_jacobian[0][1] - (-0.03032801251722862)).abs() < 1e-10);
        assert!((inverse_metric_jacobian[1][0] - (-0.03032801251722862)).abs() < 1e-10);  assert!((inverse_metric_jacobian[1][1] - 1.545339934810282).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - (-3.053461330737861)).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 3.213040008825198).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 48.11228962980825).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0129321359214928).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.3112317298425537).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0415514701153219).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.9999037381015994).abs() < 1e-10);
        assert!((quality.volume_metric - 8.0441179762606669).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.9983872155077358).abs() < 1e-10);
    }

    #[test]
    // Test 20: 2D quadratic triangle element in 3D mesh - BAD Jacobian Ratio 
    fn test_2d_quadratic_triangle_in_3d_mesh_bad_jacobianratio() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-0.5, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![0.5, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![-0.5, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![-0.242, 0.0, 0.0] },
            Node { id: 4, coordinates: vec![0.0, 0.5, 0.0] },
            Node { id: 5, coordinates: vec![-0.5, 0.258, 0.0] },
            
        ];
        let element_type = ElementType::QuadraticTriangle;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 0.001024).abs() < 1e-10);
        assert!((max_detj - 1.968).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }

        
        // Construct metric Jacobian matrix only for knupps metrics
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices, use as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices, compute J * J^T to get a square matrix [element_dim × element_dim]
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        assert!((metric_jacobian[0][0] - 0.001024).abs() < 1e-10);  assert!((metric_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((metric_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((metric_jacobian[1][1] - 0.001024).abs() < 1e-10);

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((r_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((r_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((r_matrix[1][1] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((q_matrix[1][1] - 1.0).abs() < 1e-10);  
          
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((d_matrix[1][1] - 1.0).abs() < 1e-10);  
        

        // Calculate inverses of matrices
        let inverse_metric_jacobian = GeometricAnalysis::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_metric_jacobian[0][0] - 976.5625).abs() < 1e-10);  assert!((inverse_metric_jacobian[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_metric_jacobian[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_metric_jacobian[1][1] - 976.5625).abs() < 1e-10);  

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_q_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_q_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_q_matrix[1][1] - 1.0).abs() < 1e-10);  
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10);  assert!((inverse_d_matrix[0][1] - 0.0).abs() < 1e-10);
        assert!((inverse_d_matrix[1][0] - 0.0).abs() < 1e-10);  assert!((inverse_d_matrix[1][1] - 1.0).abs() < 1e-10);  
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 1.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.001024).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.000520325203252).abs() < 1e-10);
    }

    #[test]
    // Test 21: 2D quadratic triangle element in 3D mesh - Inverted element
    fn test_2d_quadratic_triangle_in_3d_mesh_inverted() { 
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4, 5] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-0.5, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![0.5, 0.0, 0.0] },  
            Node { id: 2, coordinates: vec![-0.5, 1.0, 0.0] },  
            Node { id: 3, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 4, coordinates: vec![0.0, 0.5, 0.0] },
            Node { id: 5, coordinates: vec![-0.5, 0.25, 0.0] },
            
        ];
        let element_type = ElementType::QuadraticTriangle;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!(has_nonpositive_detj); 
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.0).abs() < 1e-10);
    }

    #[test]
    // Test 22: Perfect 1D linear line element in 3D mesh
    fn test_1d_linear_line_in_3d_mesh_perfect_reference() { 
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![1.0, 2.0, 3.0] },
            Node { id: 1, coordinates: vec![2.0, 2.0, 3.0] },  

        ];
        let element_type = ElementType::Line;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!((min_detj - 1.0).abs() < 1e-10);
        assert!((max_detj - 1.0).abs() < 1e-10);
        assert!(!has_nonpositive_detj);

        let eval_point_tuple = match min_point.len() {
            1 => (min_point[0], 0.0, 0.0),
            2 => (min_point[0], min_point[1], 0.0),
            3 => (min_point[0], min_point[1], min_point[2]),
            _ => (0.0, 0.0, 0.0),
        };

        // Evaluate Jacobian matrix at the worst (minimum) point
        let mut evaluated_jacobian = vec![vec![0.0; mesh_dim]; element_dim];  // [element_dim][mesh_dim]
        for i in 0..element_dim {
            for j in 0..mesh_dim {
                evaluated_jacobian[i][j] = MonomialPolynomial::evaluate(
                    &jacobian_matrix[i][j], eval_point_tuple
                ).map_err(|e| ElementError::GeometryError(
                    format!("Failed to evaluate Jacobian: {}", e)
                )).unwrap();
            }
        }
  
        
        // Construct metric Jacobian matrix only for knupps metrics
        let metric_jacobian = if mesh_dim == element_dim {
            // For square matrices, use as-is: [element_dim × mesh_dim] = [n × n]
            evaluated_jacobian.clone()
        } else {
            // For non-square matrices, compute J * J^T to get a square matrix [element_dim × element_dim]
            let mut square_matrix = vec![vec![0.0; element_dim]; element_dim];
            for i in 0..element_dim {
                for j in 0..element_dim {
                    for k in 0..mesh_dim {
                        // (J * J^T)[i][j] = Σ_k J[i][k] * J[j][k]
                        square_matrix[i][j] += evaluated_jacobian[i][k] * evaluated_jacobian[j][k];
                    }
                }
            }
            square_matrix
        };

        assert!((metric_jacobian[0][0] - 1.0).abs() < 1e-10); 

        let (r_matrix, q_matrix, d_matrix) = 
            GeometricAnalysis::jacobian_factorization(&metric_jacobian, element_dim).unwrap();

        assert!((r_matrix[0][0] - 1.0).abs() < 1e-10);  

        assert!((q_matrix[0][0] - 1.0).abs() < 1e-10);   
          
        assert!((d_matrix[0][0] - 1.0).abs() < 1e-10);    
        

        // Calculate inverses of matrices
        let inverse_metric_jacobian = GeometricAnalysis::inverse_matrix(&metric_jacobian, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_q_matrix = GeometricAnalysis::inverse_matrix(&q_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        let inverse_d_matrix = GeometricAnalysis::inverse_matrix(&d_matrix, element_dim, element_dim)
            .map_err(|e| ElementError::GeometryError(format!("Failed to calculate inverse: {:?}", e))).unwrap();
        
        assert!((inverse_metric_jacobian[0][0] - 1.0).abs() < 1e-10);   

        assert!((inverse_q_matrix[0][0] - 1.0).abs() < 1e-10);   
        
        assert!((inverse_d_matrix[0][0] - 1.0).abs() < 1e-10); 
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 1.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 1.0).abs() < 1e-10);
        assert!((quality.length_ratio - 1.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 1.0).abs() < 1e-10);
        assert!((quality.volume_metric - 1.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    // Test 23: 1D linear line element in 3D mesh - Inverted element
    fn test_1d_linear_line_in_3d_mesh_inverted() { 
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![2.0, 2.0, 3.0] },
            Node { id: 1, coordinates: vec![2.0, 2.0, 3.0] },  

        ];
        let element_type = ElementType::Line;
        let mesh_dim = 3;

        // Get node coordinates for this element
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        // Get multiple evaluation points
        let evaluation_points = GeometricAnalysis::get_evaluation_points(&element_type);
        // Get element dimensionality (1D, 2D, 3D)
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        // Get element order (linear, quadratic, etc.)
        let element_order = ElementType::get_element_order(&element_type).unwrap();

        // Retrieve shape functions and their derivatives for this element type
        let shape_function = ElementType::get_shape_functions(&element_type)
            .ok_or_else(|| {
                ElementError::GeometryError("Failed to get shape functions".to_string())
            }).unwrap();

        let num_nodes = shape_function.num_nodes;

        // Calculate Jacobian matrix and its determinant as polynomial
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        // Find min, max Jacobian determinant AND check for zeros
        let (min_point, min_detj, max_detj, has_nonpositive_detj) = GeometricAnalysis::find_min_max_jacobian_determinant(
            &jacobian_matrix,
            &evaluation_points,
            element_dim,
            mesh_dim,
        ).unwrap();
        assert!(has_nonpositive_detj); 
        

        let quality = GeometricAnalysis::calculate_element_quality(
            &element,
            &element_type,
            &nodes,
            mesh_dim,
        ).unwrap();

        println!("Perfect Square Test - shape: {:.3}, skew: {:.3}, length: {:.3}, orient: {:.3}, volume: {:.3}, jac_ratio: {:.3}",
                 quality.shape_metric, quality.skewness_metric, quality.length_ratio, 
                 quality.orientation_metric, quality.volume_metric, quality.jacobian_ratio);

        assert!((quality.shape_metric - 0.0).abs() < 1e-10);
        assert!((quality.skewness_metric - 0.0).abs() < 1e-10);
        assert!((quality.length_ratio - 0.0).abs() < 1e-10);
        assert!((quality.orientation_metric - 0.0).abs() < 1e-10);
        assert!((quality.volume_metric - 0.0).abs() < 1e-10);
        assert!((quality.jacobian_ratio - 0.0).abs() < 1e-10);
    }
}
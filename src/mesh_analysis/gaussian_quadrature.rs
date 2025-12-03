use core::num;
use std::f64;
use std::vec;
use crate::structs_and_impls::*;
use crate::error::*;
use super::geometric_analysis::GeometricAnalysis;
use once_cell::sync::Lazy;
use rayon::prelude::*;
use std::time::Instant;


// Precompute factorials for efficient error calculation
static FACT: Lazy<[f64; 100]> = Lazy::new(|| {
    let mut v = [1.0_f64; 100];
    for i in 1..v.len() {
        v[i] = v[i - 1] * (i as f64);
    }
    v
});

/// Calculate factorial with fallback to Stirling's approximation for large n
fn factorial(n: usize) -> f64 {
    if n < FACT.len() {
        FACT[n]
    } else {
        // Stirling with first correction
        let nf = n as f64;
        (2.0 * std::f64::consts::PI * nf).sqrt()
            * (nf / std::f64::consts::E).powi(n as i32)
            * (1.0 + 1.0 / (12.0 * nf))
    }
}

pub struct GaussianQuadrature;

impl GaussianQuadrature {

    /// Detect polynomial orders in 1D coefficient vector
    /// Returns the maximum degree based on non-zero coefficients
    pub fn detect_polynomial_order_1d(coeffs: &[f64]) -> usize {
        if coeffs.is_empty() {
            return 0;
        }

        let mut max_deg = 0;
        for (i, &coeff) in coeffs.iter().enumerate() {
            if coeff.abs() > 1e-12 {
                max_deg = i;
            }
        }
        max_deg
    }

    /// Detect polynomial orders in 2D coefficient matrix
    /// Returns (max_degree_x, max_degree_y) based on non-zero coefficients
    pub fn detect_polynomial_orders_2d(coeff_matrix: &[Vec<f64>]) -> (usize, usize) {
        if coeff_matrix.is_empty() {
            return (0, 0);
        }

        let (mut max_i, mut max_j) = (0, 0);
        // Scan through all coefficients to find maximum non-zero indices
        for i in 0..coeff_matrix.len() {
            for j in 0..coeff_matrix[i].len() {
                if coeff_matrix[i][j].abs() > 1e-12 { // Tolerance for "non-zero"
                    max_i = max_i.max(i);
                    max_j = max_j.max(j);
                }
            }
        }
        (max_i, max_j)
    }

    /// Detect polynomial orders in 3D coefficient tensor
    /// Returns (max_degree_x, max_degree_y, max_degree_z)
    pub fn detect_polynomial_orders_3d(coeff_tensor: &[Vec<Vec<f64>>]) -> (usize, usize, usize) {
        if coeff_tensor.is_empty() {
            return (0, 0, 0);
        }

        let (mut max_i, mut max_j, mut max_k) = (0, 0, 0);
        // Three-dimensional scan for non-zero coefficients
        for i in 0..coeff_tensor.len() {
            for j in 0..coeff_tensor[i].len() {
                for k in 0..coeff_tensor[i][j].len() {
                    if coeff_tensor[i][j][k].abs() > 1e-12 {
                        max_i = max_i.max(i);
                        max_j = max_j.max(j);
                        max_k = max_k.max(k);
                    }
                }
            }
        }
        (max_i, max_j, max_k)
    }

    /// Check if material matrix C is symmetric (in physical space)
    fn is_symmetric_material_matrix(matrix: &Vec<Vec<Vec<f64>>>) -> bool {
        let n = matrix.len();
        
        // Check dimensions are square
        for row in matrix {
            if row.len() != n {
                return false;
            }
        }
        
        // Check symmetry: C[i][j] == C[j][i]
        for i in 0..n {
            for j in (i+1)..n {
                // Compare polynomial coefficients
                let poly_ij = &matrix[i][j];
                let poly_ji = &matrix[j][i];
                
                // Handle different lengths
                let max_len = poly_ij.len().max(poly_ji.len());
                for k in 0..max_len {
                    let coeff_ij = if k < poly_ij.len() { poly_ij[k] } else { 0.0 };
                    let coeff_ji = if k < poly_ji.len() { poly_ji[k] } else { 0.0 };
                    
                    if (coeff_ij - coeff_ji).abs() > 1e-12 {
                        return false;
                    }
                }
            }
        }
        
        true
    }

    /// Check if a polynomial component is effectively zero
    fn is_zero_polynomial(poly: &[f64]) -> bool {
        poly.iter().all(|&c| c.abs() < 1e-12)
    }

    /// Find optimal Gauss points for entire mesh - PARALLEL VERSION WITH FALLBACKS
    pub fn find_optimal_gauss_points_number_mesh(
        mesh_data: &MeshData,
        int_type: IntegrationType,
        tolerance: f64,
        material_property: &MaterialProperty, 
    ) -> Result<GaussianPointNumberReport, GaussError> {

        let mesh_dim = mesh_data.dimension;
        let total_start_time = Instant::now();

        // Collect all elements that need processing
        let collection_start = Instant::now();
        let all_elements: Vec<(usize, &Element, &ElementType)> = mesh_data.element_type_info
            .iter()
            .filter(|type_info| !matches!(type_info.element_type, ElementType::Vertex))
            .flat_map(|type_info| {
                let start_idx = type_info.start_index;
                let end_idx = start_idx + type_info.num_elements;
                
                (start_idx..end_idx)
                    .filter(|&element_idx| element_idx < mesh_data.elements.len())
                    .map(move |element_idx| {
                        let element = &mesh_data.elements[element_idx];
                        (element_idx, element, &type_info.element_type)
                    })
            })
            .collect();
        let collection_time = collection_start.elapsed();

        println!("Element collection took: {:.2?}", collection_time);
        println!("Processing {} elements in parallel...", all_elements.len());

        let parallel_start = Instant::now();

        // Process elements in parallel
        let gauss_point_numbers: Vec<GaussianPointNumber> = all_elements
            .par_iter()
            .map(|(element_idx, element, element_type)| {
                let result = match Self::find_optimal_gauss_points(
                    int_type.clone(),
                    element,
                    element_type,
                    &mesh_data.nodes,
                    mesh_dim,
                    tolerance,
                    material_property,
                ) {
                    Ok(gauss_point_number) => {
                        if element_idx % 100 == 0 {
                            println!("Processed element {} successfully", element_idx);
                        }
                        gauss_point_number
                    }
                    Err(e) => {
                        println!("Element {} failed: {:?}, using fallback", element_idx, e);
                        // Fallback with reasonable defaults
                        GaussianPointNumber {
                            element_id: element.id,
                            theoretical_number: 2,
                            optimal_number: 2,
                        }
                    }
                };
                result
            })
            .collect();

        let parallel_time = parallel_start.elapsed();
        let total_time = total_start_time.elapsed();

        // Calculate basic statistics
        let elements_per_second = if parallel_time.as_secs_f64() > 0.0 {
            gauss_point_numbers.len() as f64 / parallel_time.as_secs_f64()
        } else {
            0.0
        };

        let avg_element_time = if !gauss_point_numbers.is_empty() {
            parallel_time / gauss_point_numbers.len() as u32
        } else {
            std::time::Duration::ZERO
        };

        // Print timing report
        println!("\n{}", "=".repeat(50));
        println!("GAUSSIAN QUADRATURE TIMING");
        println!("{}", "=".repeat(50));
        println!("Total elements:      {}", gauss_point_numbers.len());
        println!("Collection:          {:.2?}", collection_time);
        println!("Parallel processing: {:.2?}", parallel_time);
        println!("Total:               {:.2?}", total_time);
        println!("Average per element: {:.2?}", avg_element_time);
        println!("Elements per second: {:.1}", elements_per_second);
        println!("{}", "=".repeat(50));

        Ok(GaussianPointNumberReport {
            total_elements: gauss_point_numbers.len(),
            gauss_point_numbers,
        })
    }

    /// Find optimal Gaussian points that minimize error under tolerance
    pub fn find_optimal_gauss_points(
        int_type: IntegrationType,
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        mesh_dim: usize,
        tolerance: f64,
        material_property: &MaterialProperty,
    ) -> Result<GaussianPointNumber, GaussError> {

        // Get element dimension
        let element_dim = ElementType::get_element_dimension(element_type)
            .ok_or_else(|| GaussError::UnsupportedDimension(0))?;

        // Get element nodes
        let element_nodes: Vec<Node> = GeometricAnalysis::get_element_nodes(element, nodes)
            .map_err(|e| e)?;

        // Calculate the integrand
        let integrand = Self::calculate_integrand(
            &int_type,
            &element_type,
            &element_nodes,
            material_property,
            mesh_dim,
        )?;

        // Calculate theoretical maximum points based on polynomial degree
        let max_total_degree = Self::find_max_degree_in_matrix(&integrand);
        let max_num_gp = ((max_total_degree as f64 + 1.0) / 2.0).ceil() as usize;
        let theoretical_points = max_num_gp;

        // Find optimal points
        let mut optimal_points = 1;
        let mut min_error = 100.0; // Start with a large error

        for n in 1..=max_num_gp {
            let error_matrix = match Self::calculate_error(&integrand, n, element_dim) {
                Ok(matrix) => matrix,
                Err(_) => continue,
            };

            // Find maximum error in the matrix
            let max_error = Self::find_max_in_matrix(&error_matrix);

            if max_error < min_error {
                min_error = max_error;
                optimal_points = n;
            }

            // If we meet tolerance, we can stop early
            if max_error <= tolerance {
                break;
            }
        }

        Ok(GaussianPointNumber {
            element_id: element.id,
            theoretical_number: theoretical_points,
            optimal_number: optimal_points,
        })
    }

    /// Helper function to find maximum degree in integrand matrix
    fn find_max_degree_in_matrix(matrix: &[Vec<Vec<f64>>]) -> u32 {
        let mut max_deg = 0;
        for row in matrix {
            for poly in row {
                let deg = MonomialPolynomial::total_degree_polynomial(poly);
                max_deg = max_deg.max(deg);
            }
        }
        max_deg
    }

    /// Helper function to find maximum value in error matrix
    fn find_max_in_matrix(matrix: &[Vec<f64>]) -> f64 {
        let mut max_val = 0.0;
        for row in matrix {
            for &val in row {
                let abs_val = val.abs();
                if abs_val > max_val {
                    max_val = abs_val;
                }
            }
        }
        max_val
    }

    /// Check if a polynomial is constant (all coefficients beyond first are zero)
    fn is_constant_polynomial(poly: &[f64]) -> bool {
        poly.iter().enumerate().all(|(i, &c)| i == 0 || c.abs() < 1e-12)
    }

    /// Check if a material property matrix has constant components
    fn is_constant_material_matrix(matrix: &Vec<Vec<Vec<f64>>>) -> bool {
        for row in matrix {
            for poly in row {
                if !Self::is_constant_polynomial(poly) {
                    return false;
                }
            }
        }
        true
    }

    // Calculate integrand for given element and integration type
    pub fn calculate_integrand(
        int_type: &IntegrationType,
        element_type: &ElementType,
        element_nodes: &Vec<Node>,
        material_property: &MaterialProperty,
        mesh_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {

        // Get element properties
        let element_dim = ElementType::get_element_dimension(element_type)
            .ok_or(GaussError::UnsupportedDimension(0))?;
        let element_order = ElementType::get_element_order(element_type)
            .ok_or(GaussError::InvalidElement("Element order not found".to_string()))?;
        
        // Get shape functions and their derivatives
        let shape_function = ElementType::get_shape_functions(element_type)
            .ok_or_else(|| GaussError::InvalidElement("Shape functions not found".to_string()))?;
        
        let num_nodes = shape_function.num_nodes;
        
        // Calculate Jacobian with correct dimensions [element_dim × mesh_dim]
        // This maps from parametric space (xi,eta,psi) to physical space (x,y,z)
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,    // Number of parametric dimensions (rows)
            num_nodes,
            mesh_dim,       // Number of physical dimensions (columns)
            element_order,
        )?;
        
        // Extract determinant (or metric) for integration measure
        let det_j = GaussianQuadrature::calculate_determinant_monomial(
            &jacobian_matrix,
            mesh_dim,
            element_dim,
        )?;
        
        // Validate determinant
        if det_j.is_empty() {
            return Err(GaussError::GeometryError("Jacobian determinant is empty".to_string()));
        }

        // Initialize integrand matrix (size depends on mesh dimension and number of nodes)
        let matrix_size = mesh_dim * num_nodes;
        let mut integrand = vec![vec![vec![0.0]; matrix_size]; matrix_size];

        match int_type {
            // Mass matrix: M = ∫ ρ(x) N_i(x) N_j(x) det(J) dξ
            // Represents inertial properties of the element
            IntegrationType::Mass => {
                // Get density material property
                let rho_physical = material_property.as_scalar()?;

                if rho_physical.is_empty() {
                    return Err(GaussError::InvalidMaterialProperty(
                        "Density polynomial is empty".to_string()
                    ));
                }

                // Transform density from physical to isoparametric coordinates if non-constant
                let rho_isoparametric = if MonomialPolynomial::is_constant_polynomial(rho_physical) {
                    rho_physical.to_vec()  // No transformation needed for constant
                } else {
                    Self::transform_material_property_scalar(
                        rho_physical,
                        element_nodes,
                        &shape_function.values,
                        mesh_dim,
                    )?
                };

                // Calculate ρ * det(J) (integrand without shape functions)
                let rho_detj = MonomialPolynomial::multiply_optimized(&rho_isoparametric, &det_j);

                // OPTIMIZATION: Mass matrix is symmetric, compute only upper triangular part
                let mut mass_entries = vec![vec![vec![0.0]; num_nodes]; num_nodes];
                
                for i in 0..num_nodes {
                    for j in i..num_nodes {  // Only compute upper triangular
                        // Calculate N_i * N_j
                        let ni_nj = MonomialPolynomial::multiply_optimized(
                            &shape_function.values[i],
                            &shape_function.values[j]
                        );
                        // Multiply by ρ * det(J)
                        mass_entries[i][j] = MonomialPolynomial::multiply_optimized(&rho_detj, &ni_nj);
                        
                        // Mirror to lower triangular (symmetry)
                        if i != j {
                            mass_entries[j][i] = mass_entries[i][j].clone();
                        }
                    }
                }
                
                // Distribute to all spatial dimensions (block diagonal structure)
                // Each DOF (u, v, w) at node i couples with same DOF at node j
                for i in 0..num_nodes {
                    for j in 0..num_nodes {
                        for dim in 0..mesh_dim {
                            let row = mesh_dim * i + dim;
                            let col = mesh_dim * j + dim;
                            integrand[row][col] = mass_entries[i][j].clone();
                        }
                    }
                }
                
                Ok(integrand)
            }
            
            // Stiffness matrix: K = ∫ B^T * C * B * det(J) dξ
            // Represents elastic properties of the element
            IntegrationType::Stiffness => {
                // FIXED: Calculate inverse Jacobian with correct parameter order
                // J^(-1) transforms derivatives from parametric to physical space
                // Dimensions: [mesh_dim × element_dim]
                let inv_jacobian = Self::inverse_matrix_monomial(
                    &jacobian_matrix,
                    element_dim,  // FIXED: rows of J
                    mesh_dim,     // FIXED: columns of J
                ).map_err(|e| {
                    GaussError::GeometryError(format!("Failed to calculate inverse Jacobian: {:?}", e))
                })?;
                
                // Build B matrix (strain-displacement matrix)
                // B relates nodal displacements to element strains
                let b_matrix = Self::build_b_matrix(
                    &shape_function.derivatives,
                    &inv_jacobian,
                    num_nodes,
                    mesh_dim,
                    element_dim,
                )?;

                // Get elasticity matrix (constitutive relationship)
                let c_matrix_physical = material_property.as_matrix()?;

                // Validate material matrix dimensions
                let expected_c_size = match mesh_dim {
                    1 => 1,  // 1D: σ = E * ε
                    2 => 3,  // 2D: [σ_xx, σ_yy, τ_xy]^T = C * [ε_xx, ε_yy, γ_xy]^T
                    3 => 6,  // 3D: 6 stress components
                    _ => return Err(GaussError::UnsupportedDimension(mesh_dim)),
                };

                if c_matrix_physical.len() != expected_c_size {
                    return Err(GaussError::InvalidMaterialProperty(format!(
                        "Material matrix has {} rows, expected {} for {}D elements",
                        c_matrix_physical.len(), expected_c_size, mesh_dim
                    )));
                }

                // Validate each row has correct number of columns
                for (i, row) in c_matrix_physical.iter().enumerate() {
                    if row.len() != expected_c_size {
                        return Err(GaussError::InvalidMaterialProperty(format!(
                            "Material matrix row {} has {} columns, expected {}",
                            i, row.len(), expected_c_size
                        )));
                    }
                }

                // Validate each component is a non-empty polynomial
                for (i, row) in c_matrix_physical.iter().enumerate() {
                    for (j, poly) in row.iter().enumerate() {
                        if poly.is_empty() {
                            return Err(GaussError::InvalidMaterialProperty(format!(
                                "Material matrix component ({},{}) is empty",
                                i, j
                            )));
                        }
                    }
                }

                // Detect symmetry and constancy BEFORE transformation
                let is_symmetric = Self::is_symmetric_material_matrix(&c_matrix_physical);
                let is_constant = Self::is_constant_material_matrix(&c_matrix_physical);
                
                // Use optimized transformation for symmetric matrices
                let c_matrix_isoparametric = if is_constant {
                    // No transformation needed for constant
                    c_matrix_physical.clone()
                } else if is_symmetric {
                    // Transform only upper triangle, then mirror
                    Self::transform_material_property_tensor_symmetric(
                        c_matrix_physical,
                        element_nodes,
                        &shape_function.values,
                        mesh_dim,
                    )?
                } else {
                    // Full transformation for non-symmetric
                    Self::transform_material_property_tensor(
                        c_matrix_physical,
                        element_nodes,
                        &shape_function.values,
                        mesh_dim,
                    )?
                };

                let voigt_size = b_matrix.len();

                if voigt_size != expected_c_size {
                    return Err(GaussError::GeometryError(format!(
                        "B matrix has {} rows but material matrix expects {}",
                        voigt_size, expected_c_size
                    )));
                }
                
                //Build sparsity map for C matrix (track zero entries)
                let mut c_is_zero = vec![vec![false; voigt_size]; voigt_size];
                for i in 0..voigt_size {
                    for j in 0..voigt_size {
                        c_is_zero[i][j] = Self::is_zero_polynomial(&c_matrix_isoparametric[i][j]);
                    }
                }
                
                // Step 1: Compute B^T × C, exploiting C sparsity and symmetry
                // Result: [matrix_size × voigt_size] matrix
                let mut bt_times_c = vec![vec![vec![0.0]; voigt_size]; matrix_size];
                
                if is_symmetric {
                    // Exploit C symmetry for better cache locality
                    for i in 0..matrix_size {
                        for j in 0..voigt_size {
                            let mut sum = vec![0.0];
                            
                            // B^T[i][k] = B[k][i]
                            for k in 0..voigt_size {
                                if c_is_zero[k][j] {
                                    continue;  // Skip zero entries
                                }
                                
                                // Compute B[k][i] × C[k][j]
                                let term = MonomialPolynomial::multiply_optimized(
                                    &b_matrix[k][i], 
                                    &c_matrix_isoparametric[k][j]
                                );
                                sum = MonomialPolynomial::add_optimized(&sum, &term);
                            }
                            
                            bt_times_c[i][j] = sum;
                        }
                    }
                } else {
                    // Non-symmetric: standard computation with zero-skipping
                    for i in 0..matrix_size {
                        for j in 0..voigt_size {
                            let mut sum = vec![0.0];
                            
                            for k in 0..voigt_size {
                                if c_is_zero[k][j] {
                                    continue;
                                }
                                
                                let term = MonomialPolynomial::multiply_optimized(
                                    &b_matrix[k][i], 
                                    &c_matrix_isoparametric[k][j]
                                );
                                sum = MonomialPolynomial::add_optimized(&sum, &term);
                            }
                            
                            bt_times_c[i][j] = sum;
                        }
                    }
                }
                
                // Step 2: Compute K = (B^T × C) × B × det(J)
                if is_symmetric {
                    // K is also symmetric, compute only upper triangle
                    for i in 0..matrix_size {
                        for j in i..matrix_size {  // Only upper triangle
                            let mut sum = vec![0.0];
                            for k in 0..voigt_size {
                                // (B^T × C)[i][k] × B[k][j]
                                let term = MonomialPolynomial::multiply_optimized(
                                    &bt_times_c[i][k], 
                                    &b_matrix[k][j]
                                );
                                sum = MonomialPolynomial::add_optimized(&sum, &term);
                            }
                            integrand[i][j] = MonomialPolynomial::multiply_optimized(&sum, &det_j);
                            
                            // Mirror to lower triangle
                            if i != j {
                                integrand[j][i] = integrand[i][j].clone();
                            }
                        }
                    }
                } else {
                    // Non-symmetric: compute full matrix
                    for i in 0..matrix_size {
                        for j in 0..matrix_size {
                            let mut sum = vec![0.0];
                            for k in 0..voigt_size {
                                let term = MonomialPolynomial::multiply_optimized(
                                    &bt_times_c[i][k], 
                                    &b_matrix[k][j]
                                );
                                sum = MonomialPolynomial::add_optimized(&sum, &term);
                            }
                            integrand[i][j] = MonomialPolynomial::multiply_optimized(&sum, &det_j);
                        }
                    }
                }
                
                Ok(integrand)
            }
        }
    }

    // Build B matrix with updated shape derivatives format and optimized operations
    fn build_b_matrix(
        shape_derivatives: &Vec<Vec<Vec<f64>>>, // [param_dim][node][coefficients]
        inv_jacobian: &Vec<Vec<Vec<f64>>>,      // [mesh_dim][element_dim]
        num_nodes: usize,
        mesh_dim: usize,
        element_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {
        
        // Step 1: Compute dN/dx = J^(-1) * dN/dxi for all nodes
        // This transforms shape function derivatives from parametric to physical space
        // Result: [num_nodes × mesh_dim] matrix
        let mut dn_dx = vec![vec![vec![0.0]; mesh_dim]; num_nodes];

        for node in 0..num_nodes {
            for phys_dim in 0..mesh_dim {  // Loop over x, y, z
                let mut sum = vec![0.0];
                for param_dim in 0..element_dim {  // Loop over xi, eta, psi
                    // Compute: dN_node/dx = Σ_param J^(-1)[phys][param] * dN_node/dxi
                    // inv_jacobian[phys_dim][param_dim] since J^(-1) is [mesh_dim × element_dim]
                    let term = MonomialPolynomial::multiply_optimized(
                        &inv_jacobian[phys_dim][param_dim],      // J^(-1) component
                        &shape_derivatives[param_dim][node]      // dN/dxi component
                    );
                    sum = MonomialPolynomial::add_optimized(&sum, &term);
                }
                // Store result: dn_dx[node][phys_dim] = dN_node/dx
                dn_dx[node][phys_dim] = MonomialPolynomial::trim_to_valid_length_fallback(&sum);
            }
        }

        // Step 2: Build B matrix according to standard Voigt notation
        // B matrix maps nodal displacements to strain components
        match mesh_dim {
            // 1D: Only axial strain ε_xx = du/dx
            // B = [dN1/dx, dN2/dx, ...] with size [1 × num_nodes]
            1 => {
                let mut b_matrix = vec![vec![vec![0.0]; num_nodes]; 1];
                for node in 0..num_nodes {
                    b_matrix[0][node] = dn_dx[node][0].clone();  // ε_xx row
                }
                Ok(b_matrix)
            }
            
            // 2D: Plane stress/strain with 3 strain components
            // Voigt notation: ε = [ε_xx, ε_yy, γ_xy]^T
            // B matrix size: [3 × (2*num_nodes)]
            2 => {
                let mut b_matrix = vec![vec![vec![0.0]; 2 * num_nodes]; 3];
                
                for node in 0..num_nodes {
                    let u_col = 2 * node;      // u displacement column
                    let v_col = 2 * node + 1;  // v displacement column
                    
                    // Row 0: ε_xx = du/dx
                    b_matrix[0][u_col] = dn_dx[node][0].clone();
                    
                    // Row 1: ε_yy = dv/dy
                    b_matrix[1][v_col] = dn_dx[node][1].clone();
                    
                    // Row 2: γ_xy = du/dy + dv/dx (engineering shear strain)
                    b_matrix[2][u_col] = dn_dx[node][1].clone();  // du/dy
                    b_matrix[2][v_col] = dn_dx[node][0].clone();  // dv/dx
                }
                Ok(b_matrix)
            }
            
            // 3D: Full 3D stress state with 6 strain components
            // Voigt notation: ε = [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]^T
            // B matrix size: [6 × (3*num_nodes)]
            3 => {
                let mut b_matrix = vec![vec![vec![0.0]; 3 * num_nodes]; 6];
                
                for node in 0..num_nodes {
                    let u_col = 3 * node;      // u displacement column
                    let v_col = 3 * node + 1;  // v displacement column
                    let w_col = 3 * node + 2;  // w displacement column

                    // Ensure no empty polynomials (safety check)
                    if dn_dx[node][0].is_empty() { dn_dx[node][0] = vec![0.0]; }
                    if dn_dx[node][1].is_empty() { dn_dx[node][1] = vec![0.0]; }
                    if dn_dx[node][2].is_empty() { dn_dx[node][2] = vec![0.0]; }
                    
                    // Normal strains (rows 0-2)
                    b_matrix[0][u_col] = dn_dx[node][0].clone();  // ε_xx = du/dx
                    b_matrix[1][v_col] = dn_dx[node][1].clone();  // ε_yy = dv/dy
                    b_matrix[2][w_col] = dn_dx[node][2].clone();  // ε_zz = dw/dz
                    
                    // Shear strains (rows 3-5)
                    // γ_xy = du/dy + dv/dx
                    b_matrix[3][u_col] = dn_dx[node][1].clone();  // du/dy
                    b_matrix[3][v_col] = dn_dx[node][0].clone();  // dv/dx
                    
                    // γ_yz = dv/dz + dw/dy
                    b_matrix[4][v_col] = dn_dx[node][2].clone();  // dv/dz
                    b_matrix[4][w_col] = dn_dx[node][1].clone();  // dw/dy
                    
                    // γ_zx = dw/dx + du/dz
                    b_matrix[5][w_col] = dn_dx[node][0].clone();  // dw/dx
                    b_matrix[5][u_col] = dn_dx[node][2].clone();  // du/dz
                }

                // Final safety check: ensure no empty polynomials in B matrix
                for i in 0..b_matrix.len() {
                    for j in 0..b_matrix[i].len() {
                        if b_matrix[i][j].is_empty() {
                            b_matrix[i][j] = vec![0.0];
                        }
                    }
                }
                
                Ok(b_matrix)
            }
            _ => Err(GaussError::UnsupportedDimension(mesh_dim)),
        }
    }

    // Transform scalar material property (density ρ) from physical to isoparametric coordinates 
    pub fn transform_material_property_scalar(
        material_property_physical: &[f64],
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>],  
        mesh_dim: usize,
    ) -> Result<Vec<f64>, GaussError> {

        // Use safe trimming
        let material_property_physical = MonomialPolynomial::trim_to_valid_length(material_property_physical)
            .map_err(|e| GaussError::PolynomialError(format!("Failed to trim material property: {}", e)))?;

        // Build coordinate mappings: x(xi,eta,psi), y(xi,eta,psi), z(xi,eta,psi)
        let coordinate_maps = Self::build_coordinate_maps(
            element_nodes,
            shape_function_values,
            mesh_dim,
        )?;

        if coordinate_maps.iter().any(|m| m.is_empty()) {
            return Err(GaussError::GeometryError(
                "Empty coordinate mapping detected".to_string()
            ));
        }

        // Substitute physical coordinates with isoparametric coordinate polynomials
        Self::substitute_polynomial(&material_property_physical, &coordinate_maps)
    }

    // Transform tensor material property (elasticity matrix C) from physical to isoparametric
    pub fn transform_material_property_tensor(
        material_property_physical: &Vec<Vec<Vec<f64>>>,
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>], 
        mesh_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {
            
        // Validate input tensor
        for i in 0..material_property_physical.len() {
            for j in 0..material_property_physical[i].len() {
                if material_property_physical[i][j].is_empty() {
                    return Err(GaussError::InvalidMaterialProperty(
                        format!("Empty polynomial in material tensor at ({},{})", i, j)
                    ));
                }
            }
        }

        // Build coordinate maps once for all components
        let coordinate_maps = Self::build_coordinate_maps(
            element_nodes,
            shape_function_values,
            mesh_dim,
        )?;

        // Transform each component of the C matrix independently
        let matrix_size = material_property_physical.len();
        let mut material_property_isoparametric = vec![vec![vec![0.0]; matrix_size]; matrix_size];
        
        for i in 0..matrix_size {
            for j in 0..matrix_size {
                let component = &material_property_physical[i][j];
                
                // Skip transformation for constant components
                if MonomialPolynomial::is_constant_polynomial(component) {
                    material_property_isoparametric[i][j] = component.clone();
                } else {
                    // Only transform non-constant components
                    let trimmed_input = MonomialPolynomial::trim_to_valid_length(component)
                        .map_err(|e| {
                            GaussError::PolynomialError(format!("Failed to trim component ({},{}): {}", i, j, e))
                        })?;
                        
                    material_property_isoparametric[i][j] = Self::substitute_polynomial(
                        &trimmed_input,
                        &coordinate_maps,
                    )?;
                }
            }
        }

        Ok(material_property_isoparametric)
    }

    /// Transform symmetric tensor material property - only upper triangle
    pub fn transform_material_property_tensor_symmetric(
        material_property_physical: &Vec<Vec<Vec<f64>>>,
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>], 
        mesh_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {
        
        let matrix_size = material_property_physical.len();
        
        // Validate input tensor
        for i in 0..matrix_size {
            if material_property_physical[i].len() != matrix_size {
                return Err(GaussError::InvalidMaterialProperty(
                    format!("Material tensor row {} has wrong size", i)
                ));
            }
            for j in 0..matrix_size {
                if material_property_physical[i][j].is_empty() {
                    return Err(GaussError::InvalidMaterialProperty(
                        format!("Empty polynomial in material tensor at ({},{})", i, j)
                    ));
                }
            }
        }

        // Build coordinate maps once for all components
        let coordinate_maps = Self::build_coordinate_maps(
            element_nodes,
            shape_function_values,
            mesh_dim,
        )?;

        let mut material_property_isoparametric = vec![vec![vec![0.0]; matrix_size]; matrix_size];
        
        // OPTIMIZATION: Only transform upper triangle, then mirror
        for i in 0..matrix_size {
            for j in i..matrix_size {  // Only upper triangle (including diagonal)
                let component = &material_property_physical[i][j];
                
                // Skip transformation for constant components
                if MonomialPolynomial::is_constant_polynomial(component) {
                    material_property_isoparametric[i][j] = component.clone();
                } else {
                    let trimmed_input = MonomialPolynomial::trim_to_valid_length(component)
                        .map_err(|e| {
                            GaussError::PolynomialError(format!("Failed to trim component ({},{}): {}", i, j, e))
                        })?;
                        
                    material_property_isoparametric[i][j] = Self::substitute_polynomial(
                        &trimmed_input,
                        &coordinate_maps,
                    )?;
                }
                
                // Mirror to lower triangle
                if i != j {
                    material_property_isoparametric[j][i] = material_property_isoparametric[i][j].clone();
                }
            }
        }

        Ok(material_property_isoparametric)
    }

    
    // Build coordinate mappings for isoparametric formulation
    fn build_coordinate_maps(
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>],
        mesh_dim: usize,
    ) -> Result<Vec<Vec<f64>>, GaussError> {
        
        let num_nodes = element_nodes.len();
        let mut coordinate_maps = Vec::with_capacity(mesh_dim);
        
        // For each spatial dimension (x, y, z)
        for i in 0..mesh_dim {
            let mut coord_poly = vec![0.0];
            
            // Build: x_i(xi,eta,psi) = Σ_k N_k(xi,eta,psi) * x_ik
            for j in 0..num_nodes {
                // Get physical coordinate value for this node
                let coord_value = element_nodes[j].coordinates[i];  // x_ik

                // Skip if coordinate is zero (optimization)
                if coord_value.abs() < 1e-12 {
                    continue;
                }
                
                // N_k(xi,eta,psi) * x_ik
                let term = MonomialPolynomial::multiply_scalar(
                    &shape_function_values[j],
                    coord_value,
                );
                
                // Use safe addition
                coord_poly = MonomialPolynomial::add_optimized(&coord_poly, &term);
            }
            
            // If all nodes had zero coordinate, set to zero polynomial
            if coord_poly.is_empty() {
                coord_poly = vec![0.0];
            }
            
            // Use safe trimming
            let trimmed = MonomialPolynomial::trim_to_valid_length(&coord_poly)
                .map_err(|e| GaussError::PolynomialError(format!("Failed to trim coordinate map: {}", e)))?;
            coordinate_maps.push(trimmed);
        }
        
        Ok(coordinate_maps)
    }
    
    // Substitute physical coordinates with parametric coordinate polynomials
    fn substitute_polynomial(
        physical_poly: &[f64],
        coordinate_maps: &Vec<Vec<f64>>,
    ) -> Result<Vec<f64>, GaussError> {

        // Check for empty inputs
        if physical_poly.is_empty() {
            return Ok(vec![0.0]);
        }
        
        if coordinate_maps.iter().any(|m| m.is_empty()) {
            return Err(GaussError::GeometryError(
                "Empty coordinate map in substitution".to_string()
            ));
        }

        // Trim input polynomial to valid length
        let physical_poly = MonomialPolynomial::trim_to_valid_length(physical_poly)
            .map_err(|e| GaussError::PolynomialError(format!("Failed to trim physical poly: {}", e)))?;

        // Get the degree of the physical polynomial
        let physical_degree = MonomialPolynomial::infer_max_degree(physical_poly.len())
            .map_err(|e| GaussError::PolynomialError(format!("Failed to infer degree: {}", e)))?;
        
        // Get the degrees of the coordinate mappings
        let mut coord_degrees = Vec::new();
        for coord_map in coordinate_maps {
            let deg = MonomialPolynomial::infer_max_degree(coord_map.len())
                .map_err(|e| GaussError::PolynomialError(format!("Failed to infer coord degree: {}", e)))?;
            coord_degrees.push(deg);
        }
        
        // For safety, use a conservative estimate: physical_degree + max_coord_degree * physical_degree
        let max_coord_degree = coord_degrees.iter().max().unwrap_or(&0);
        let conservative_max_degree = physical_degree + max_coord_degree * physical_degree;
        
        // Generate basis for the conservative maximum degree
        let basis = MonomialPolynomial::generate_basis(conservative_max_degree);
        
        let mut result = vec![0.0];
        
        // Process each term: aᵢⱼₖ × xⁱ × yʲ × zᵏ
        for (idx, &(i, j, k)) in basis.iter().enumerate() {
            if idx >= physical_poly.len() {
                break;
            }
            
            let coefficient = physical_poly[idx];
            
            // Skip negligible terms for numerical stability
            if coefficient.abs() < 1e-12 {
                continue;
            }
            
            // Start with coefficient as constant polynomial: [aᵢⱼₖ, 0, 0, ...]
            let mut term_poly = vec![coefficient];
            
            // Multiply by x(xi,eta,psi)^i
            if i > 0 && !coordinate_maps.is_empty() {
                let x_power = Self::polynomial_power_safe(&coordinate_maps[0], i as usize);
                term_poly = MonomialPolynomial::multiply_optimized(&term_poly, &x_power);
            }
            
            // Multiply by y(xi,eta,psi)^j (if y dimension exists)
            if j > 0 && coordinate_maps.len() > 1 {
                let y_power = Self::polynomial_power_safe(&coordinate_maps[1], j as usize);
                term_poly = MonomialPolynomial::multiply_optimized(&term_poly, &y_power);
            }
            
            // Multiply by z(xi,eta,psi)^k (if z dimension exists)
            if k > 0 && coordinate_maps.len() > 2 {
                let z_power = Self::polynomial_power_safe(&coordinate_maps[2], k as usize);
                term_poly = MonomialPolynomial::multiply_optimized(&term_poly, &z_power);
            }
            
            // Accumulate: add this term to the result
            result = MonomialPolynomial::add_optimized(&result, &term_poly);
        }
        
        // If all terms were zero, return zero polynomial
        if result.is_empty() {
            result = vec![0.0];
        }
        
        // Trim output polynomial to valid length
        MonomialPolynomial::trim_to_valid_length(&result)
            .map_err(|e| GaussError::PolynomialError(format!("Failed to trim result: {}", e)))
    }

    /// Compute polynomial raised to integer power: poly^n
    fn polynomial_power_safe(poly: &Vec<f64>, power: usize) -> Vec<f64> {
        match power {
            0 => vec![1.0],
            1 => MonomialPolynomial::trim_to_valid_length_fallback(poly),
            _ => {
                let poly_trimmed = MonomialPolynomial::trim_to_valid_length_fallback(poly);
                let mut result = poly_trimmed.clone();
                for _ in 1..power {
                    result = MonomialPolynomial::multiply_optimized(&result, &poly_trimmed);
                }
                result
            }
        }
    }

    
    // J[param_dim][phys_dim] = [element_dim × mesh_dim]
    // 
    // Standard Jacobian matrix:
    //     J = [dx/dxi   dy/dxi   dz/dxi  ]
    //         [dx/deta   dy/deta   dz/deta  ]
    //         [dx/dpsi   dy/dpsi   dz/dpsi  ]
    //
    // Dimensions: [element_dim × mesh_dim]
    // Rows: parametric coordinates (xi, eta, psi)
    // Columns: physical coordinates (x, y, z)

    // Calculate Jacobian matrix with standard indexing
    pub fn calculate_jacobian_monomial(
        element_nodes: &Vec<Node>,
        shape_derivatives: &Vec<Vec<Vec<f64>>>, // [param_dim][node][coefficients]
        element_dim: usize,
        num_nodes: usize,
        mesh_dim: usize,
        element_order: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, ElementError> {

        // Determine maximum polynomial length based on element order
        let max_len = match element_order {
            1 => 4,  // Linear: [1, x, y, z]
            2 => 10, // Quadratic: [1, x, y, z, x^2, xy, xz, y^2, yz, z^2]
            3 => 20, // Cubic: [1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, ...]
            _ => return Err(ElementError::GeometryError(format!(
                "Unsupported element order {}", element_order
            ))),
        };

        // Initialize Jacobian matrix 
        // J[parametric_dim][physical_dim] = [element_dim × mesh_dim]
        // Each entry is a polynomial represented as Vec<f64>
        let mut jacobian_matrix: Vec<Vec<Vec<f64>>> = 
            vec![vec![vec![0.0; max_len]; mesh_dim]; element_dim];

        // Build Jacobian 
        // J[param_dim][phys_dim] = Σ_node (x_node_phys * dN_node/dxi)
        // This computes dx/dxi for each entry
        for param_dim in 0..element_dim {   // Loop over rows: xi, eta, psi (parametric dims)
            for phys_dim in 0..mesh_dim {   // Loop over columns: x, y, z (physical dims)
                for node in 0..num_nodes {  // Sum over all nodes in element
                    // Get physical coordinate of this node (e.g., x_node, y_node, z_node)
                    let coord = element_nodes[node].coordinates[phys_dim];
                    
                    // Get dN_node/dxi (shape function derivative for this node)
                    let derivative_poly = &shape_derivatives[param_dim][node];
                    
                    // Multiply: x_node * dN_node/dxi
                    let scaled_poly = MonomialPolynomial::multiply_scalar(derivative_poly, coord);
                    
                    // Accumulate: J[param_dim][phys_dim] += x_node * dN_node/dxi
                    jacobian_matrix[param_dim][phys_dim] = MonomialPolynomial::add_optimized(
                        &jacobian_matrix[param_dim][phys_dim], 
                        &scaled_poly
                    );
                }
                // Trim polynomial to valid length (remove trailing zeros)
                jacobian_matrix[param_dim][phys_dim] = 
                    MonomialPolynomial::trim_to_valid_length(&jacobian_matrix[param_dim][phys_dim]).unwrap();
            }
        }

        Ok(jacobian_matrix)
    }
    
    // Calculate inverse of Jacobian matrix
    // For square matrices: J^(-1) has same dimensions as J
    // For non-square matrices: pseudo-inverse J⁺ = J^T(JJ^T)^(-1)
    pub fn inverse_matrix_monomial(
        matrix: &Vec<Vec<Vec<f64>>>, 
        element_dim: usize,  // Number of rows of J (parametric dimensions)
        mesh_dim: usize      // Number of columns of J (physical dimensions)
    ) -> Result<Vec<Vec<Vec<f64>>>, ElementError> {
        
        match (element_dim, mesh_dim) {
            // 1x1 matrix: inverse is simply 1/det(J)
            (1, 1) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0]; element_dim]; mesh_dim];
                let det_matrix = &matrix[0][0];
                let inv_determinant = Self::calculate_inverse_monomial(&det_matrix)?;
                inverse_matrix[0][0] = inv_determinant;
                Ok(inverse_matrix)
            },
            
            // 2x2 and 3x3 square matrices: J^(-1) = adj(J) / det(J)
            // Result has same dimensions [n×n] as input
            (2, 2) | (3, 3) => {
                // Initialize inverse matrix with same dimensions as J
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0]; element_dim]; mesh_dim];
                
                // Calculate adjoint matrix (transpose of cofactor matrix)
                let adjoint = Self::calculate_adjoint_monomial(matrix, element_dim)?;
                
                // Calculate determinant
                let det_matrix = Self::calculate_determinant_monomial(matrix, mesh_dim, element_dim)?;
                
                // Calculate 1/det(J)
                let inv_determinant = Self::calculate_inverse_monomial(&det_matrix)?;
                
                // J^(-1) = adj(J) / det(J)
                for i in 0..mesh_dim {
                    for j in 0..element_dim {
                        inverse_matrix[i][j] = MonomialPolynomial::multiply_optimized(&adjoint[i][j], &inv_determinant);
                    }
                }
                Ok(inverse_matrix)
            },
            
            // 1D element in 2D space
            // J is [1×2] = [dx/dxi, dy/dxi]
            // Pseudo-inverse J⁺ = J^T(JJ^T)^(-1) has dimensions [2×1]
            (1, 2) => {
                // Initialize pseudo-inverse [2×1]
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0]; element_dim]; mesh_dim];
                
                // Calculate G = J*J^T = (dx/dxi)^2 + (dy/dxi)^2 (scalar, [1×1])
                let mut g_scalar = vec![0.0];
                for k in 0..mesh_dim {
                    let prod = MonomialPolynomial::multiply_optimized(&matrix[0][k], &matrix[0][k]);
                    g_scalar = MonomialPolynomial::add_optimized(&g_scalar, &prod);
                }
                
                // Calculate G^(-1) = 1/G
                let inv_g = Self::calculate_inverse_monomial(&g_scalar)?;
                
                // Calculate J⁺ = J^T * G^(-1)
                // J^T is [2×1], G^(-1) is scalar, result is [2×1]
                for j in 0..mesh_dim {
                    inverse_matrix[j][0] = MonomialPolynomial::multiply_optimized(&inv_g, &matrix[0][j]);
                }
                
                Ok(inverse_matrix)
            },
            
            // 1D element in 3D space
            // J is [1×3] = [dx/dxi, dy/dxi, dz/dxi]
            // Pseudo-inverse J⁺ = J^T(JJ^T)^(-1) has dimensions [3×1]
            (1, 3) => {
                // Initialize pseudo-inverse [3×1]
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0]; element_dim]; mesh_dim];
                
                // Calculate G = J*J^T = (dx/dxi)^2 + (dy/dxi)^2 + (dz/dxi)^2 (scalar, [1×1])
                let mut g_scalar = vec![0.0];
                for k in 0..mesh_dim {
                    let prod = MonomialPolynomial::multiply_optimized(&matrix[0][k], &matrix[0][k]);
                    g_scalar = MonomialPolynomial::add_optimized(&g_scalar, &prod);
                }
                
                // Calculate G^(-1) = 1/G
                let inv_g = Self::calculate_inverse_monomial(&g_scalar)?;
                
                // Calculate J⁺ = J^T * G^(-1)
                // J^T is [3×1], G^(-1) is scalar, result is [3×1]
                for j in 0..mesh_dim {
                    inverse_matrix[j][0] = MonomialPolynomial::multiply_optimized(&inv_g, &matrix[0][j]);
                }
                
                Ok(inverse_matrix)
            },
            
            // 2D element in 3D space
            // J is [2×3] = |dx/dxi  dy/dxi  dz/dxi|
            //              |dx/deta  dy/deta  dz/deta|
            // Pseudo-inverse J⁺ = J^T(JJ^T)^(-1) has dimensions [3×2]
            (2, 3) => {
                // Initialize pseudo-inverse [3×2]
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0]; element_dim]; mesh_dim];
                
                // Calculate G = J*J^T, which is [2×2]
                // G represents the metric tensor (first fundamental form)
                let mut g_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0]; 2]; 2];
                for i in 0..2 {      // Loop over rows of G
                    for j in 0..2 {  // Loop over columns of G
                        let mut sum = vec![0.0];
                        for k in 0..3 {  // Sum over physical dimensions
                            // G[i][j] = Σ_k J[i][k] * J[j][k]
                            let prod = MonomialPolynomial::multiply_optimized(&matrix[i][k], &matrix[j][k]);
                            sum = MonomialPolynomial::add_optimized(&sum, &prod);
                        }
                        g_matrix[i][j] = sum;
                    }
                }
                
                // Calculate G^(-1) using 2×2 matrix inverse formula
                let adj_g = Self::calculate_adjoint_monomial(&g_matrix, 2)?;
                let det_g = Self::calculate_determinant_monomial(&g_matrix, 2, 2)?;
                let inv_det_g = Self::calculate_inverse_monomial(&det_g)?;
                
                // Build G^(-1) = adj(G) / det(G)
                let mut inv_g = vec![vec![vec![0.0]; 2]; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        inv_g[i][j] = MonomialPolynomial::multiply_optimized(&adj_g[i][j], &inv_det_g);
                    }
                }
                
                // Calculate J⁺ = J^T * G^(-1)
                // J^T is [3×2], G^(-1) is [2×2], result is [3×2]
                for i in 0..mesh_dim {      // Loop over rows (physical dims)
                    for j in 0..element_dim { // Loop over columns (parametric dims)
                        let mut sum = vec![0.0];
                        for k in 0..2 {  // Matrix multiplication
                            // J^T[i][k] = J[k][i] (transpose)
                            let term = MonomialPolynomial::multiply_optimized(&matrix[k][i], &inv_g[k][j]);
                            sum = MonomialPolynomial::add_optimized(&sum, &term);
                        }
                        inverse_matrix[i][j] = sum;
                    }
                }
                
                Ok(inverse_matrix)
            },
            
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian inverse calculation not implemented for [{}×{}] matrices",
                element_dim, mesh_dim
            ))),
        }
    }


    // Adjoint of a matrix with monomial polynomial entries
    pub fn calculate_adjoint_monomial(matrix: &Vec<Vec<Vec<f64>>>, element_dim: usize) -> Result<Vec<Vec<Vec<f64>>>, ElementError> {
        // Adjoint calculation only implemented for 2x2 and 3x3 matrices

        match element_dim {
            2 => {
                let adj = vec![
                    vec![matrix[1][1].clone(), MonomialPolynomial::multiply_scalar(&matrix[0][1], -1.0)],
                    vec![MonomialPolynomial::multiply_scalar(&matrix[1][0], -1.0), matrix[0][0].clone()],
                ];
                Ok(adj)
            },
            3 => {
                let a0 = &matrix[0][0];
                let a1 = &matrix[0][1];
                let a2 = &matrix[0][2];
                let b0 = &matrix[1][0];
                let b1 = &matrix[1][1];
                let b2 = &matrix[1][2];
                let c0 = &matrix[2][0];
                let c1 = &matrix[2][1];
                let c2 = &matrix[2][2];

                let adj = vec![
                    vec![
                        MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(b1, c2), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(b2, c1), -1.0)),
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(b0, c2), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(b2, c0), -1.0)), -1.0),
                        MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(b0, c1), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(b1, c0), -1.0)),
                    ],
                    vec![
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(a1, c2), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(a2, c1), -1.0)), -1.0),
                        MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(a0, c2), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(a2, c0), -1.0)),
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(a0, c1), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(a1, c0), -1.0)), -1.0),
                    ],
                    vec![
                        MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(a1, b2), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(a2, b1), -1.0)),
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(a0, b2), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(a2, b0), -1.0)), -1.0),
                        MonomialPolynomial::add_optimized(&MonomialPolynomial::multiply_optimized(a0, b1), &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply_optimized(a1, b0), -1.0)),
                    ],
                ];
                Ok(adj)
            },
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian adjoint calculation not implemented for {}x{} (element_dim x element_dim) matrices",
                element_dim, element_dim
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

        match (element_dim, mesh_dim) {
            // 1x1 matrix: determinant is the single element
            (1, 1) => Ok(matrix[0][0].clone()),

            // 2x2 square matrix: det = ad - bc
            (2, 2) => {
                let ad = MonomialPolynomial::multiply_optimized(&matrix[0][0], &matrix[1][1]);
                let bc = MonomialPolynomial::multiply_optimized(&matrix[0][1], &matrix[1][0]);
                let neg_bc = MonomialPolynomial::multiply_scalar(&bc, -1.0);
                Ok(MonomialPolynomial::add_optimized(&ad, &neg_bc))
            }

            // 3x3 square matrix: det = a(ei−fh) − b(di−fg) + c(dh−eg)
            (3, 3) => {
                let ei = MonomialPolynomial::multiply_optimized(&matrix[1][1], &matrix[2][2]);
                let fh = MonomialPolynomial::multiply_optimized(&matrix[1][2], &matrix[2][1]);
                let ei_fh = MonomialPolynomial::add_optimized(&ei, &MonomialPolynomial::multiply_scalar(&fh, -1.0));
                let term1 = MonomialPolynomial::multiply_optimized(&matrix[0][0], &ei_fh);

                let di = MonomialPolynomial::multiply_optimized(&matrix[1][0], &matrix[2][2]);
                let fg = MonomialPolynomial::multiply_optimized(&matrix[1][2], &matrix[2][0]);
                let di_fg = MonomialPolynomial::add_optimized(&di, &MonomialPolynomial::multiply_scalar(&fg, -1.0));
                let term2 = MonomialPolynomial::multiply_optimized(&matrix[0][1], &di_fg);

                let dh = MonomialPolynomial::multiply_optimized(&matrix[1][0], &matrix[2][1]);
                let eg = MonomialPolynomial::multiply_optimized(&matrix[1][1], &matrix[2][0]);
                let dh_eg = MonomialPolynomial::add_optimized(&dh, &MonomialPolynomial::multiply_scalar(&eg, -1.0));
                let term3 = MonomialPolynomial::multiply_optimized(&matrix[0][2], &dh_eg);

                let result = MonomialPolynomial::add_optimized(&term1, &MonomialPolynomial::multiply_scalar(&term2, -1.0));
                Ok(MonomialPolynomial::add_optimized(&result, &term3))
            }

            // 1D elements in 2D space: J is [1×2], metric = ||J||
            (1, 2) => { 
                let dx_sq = MonomialPolynomial::multiply_optimized(&matrix[0][0], &matrix[0][0]);
                let dy_sq = MonomialPolynomial::multiply_optimized(&matrix[0][1], &matrix[0][1]);
                let squared_length = MonomialPolynomial::add_optimized(&dx_sq, &dy_sq);
                Self::calculate_sqrt_monomial(&squared_length)
            }

            // 1D elements in 3D space: J is [1×3], metric = ||J||
            (1, 3) => {
                let dx_sq = MonomialPolynomial::multiply_optimized(&matrix[0][0], &matrix[0][0]);
                let dy_sq = MonomialPolynomial::multiply_optimized(&matrix[0][1], &matrix[0][1]);
                let dz_sq = MonomialPolynomial::multiply_optimized(&matrix[0][2], &matrix[0][2]);
                let temp = MonomialPolynomial::add_optimized(&dx_sq, &dy_sq);
                let squared_length = MonomialPolynomial::add_optimized(&temp, &dz_sq);
                Self::calculate_sqrt_monomial(&squared_length)
            }

            // 2D elements in 3D space: J is [2×3], metric = sqrt(det(J^T*J))
            (2, 3) => {
                // Build metric tensor G = J^T * J where J is [2×3]
                // G will be [3×3] but we only need [2×2] part
                let mut g = vec![vec![vec![0.0; matrix[0][0].len()]; 2]; 2];
                
                for i in 0..2 {
                    for j in 0..2 {
                        let mut sum = vec![0.0; matrix[0][0].len()];
                        for k in 0..3 {
                            // matrix[i][k] and matrix[j][k] since J is [2×3]
                            let prod = MonomialPolynomial::multiply_optimized(&matrix[i][k], &matrix[j][k]);
                            sum = MonomialPolynomial::add_optimized(&sum, &prod);
                        }
                        g[i][j] = sum;
                    }
                }

                // det(G) = g00*g11 - g01*g10
                let g00_g11 = MonomialPolynomial::multiply_optimized(&g[0][0], &g[1][1]);
                let g01_g10 = MonomialPolynomial::multiply_optimized(&g[0][1], &g[1][0]);
                let squared_area = MonomialPolynomial::add_optimized(&g00_g11, &MonomialPolynomial::multiply_scalar(&g01_g10, -1.0));
                
                Self::calculate_sqrt_monomial(&squared_area)
            }

            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for [{}×{}] matrices",
                element_dim, mesh_dim
            ))),
        }
    }

    /// Calculate square root of Jacobian determinant/metric for polynomial matrices
    pub fn calculate_sqrt_monomial(
        poly: &Vec<f64>,
    ) -> Result<Vec<f64>, ElementError> {
        let wished_n = 5; // 5-term binomial expansion
        
        // Handle zero determinant
        if poly[0].abs() < 1e-12 {
            return Err(ElementError::GeometryError(
                "Cannot take square root of zero determinant".to_string()
            ));
        }
        
        // For constant determinant, return simple square root
        if poly.len() == 1 {
            return Ok(vec![poly[0].sqrt()]);
        }
        
        // Normalize: p(x) = a₀ * (1 + q(x)) where q(x) = (p(x)/a₀ - 1)
        let a0 = poly[0];
        let sqrt_a0 = a0.sqrt();
        let normalized_poly = MonomialPolynomial::multiply_scalar(&poly, 1.0/a0);

        // Create q(x) = normalized_poly - 1
        let mut q = normalized_poly.clone();
        q[0] = 0.0; // Now q(x) = p(x)/a₀ - 1
        
        // Binomial coefficients for sqrt(1+x): [1, 1/2, -1/8, 1/16, -5/128]
        let binomial_coeffs = [1.0, 0.5, -0.125, 0.0625, -0.0390625, 0.02734375, -0.0205078125, 0.01611328125, -0.013092041015625, 0.0109100341796875, -0.009273529052734375, 0.008008956909179688];  // for 12 terms
        
        let mut term = vec![vec![0.0]; wished_n];
        term[0] = vec![1.0]; // Constant term is 1
        
        // Build terms: term[k] = q^k
        term[1] = q.clone();
        for i in 2..wished_n {
            term[i] = MonomialPolynomial::multiply_optimized(&term[i-1], &q);
        }
        
        // Sum the series: result = sum_{k=0}^{4} [binomial_coeffs[k] * term[k]]
        let mut result = vec![0.0];
        for j in 0..wished_n {
            let scaled_term = MonomialPolynomial::multiply_scalar(&term[j], binomial_coeffs[j]);
            result = MonomialPolynomial::add_optimized(&result, &scaled_term);
        }
        
        // Multiply by sqrt(a₀)
        Ok(MonomialPolynomial::multiply_scalar(&result, sqrt_a0))
    }

    /// Calculate inverse of Jacobian determinant/metric for polynomial matrices
    pub fn calculate_inverse_monomial( 
        poly: &Vec<f64>, 
    ) -> Result<Vec<f64>, ElementError> {
        // Handle constant polynomial case directly
        if poly.len() == 1 {
            if poly[0].abs() < 1e-12 {
                return Err(ElementError::GeometryError(
                    "Cannot invert zero polynomial".to_string()
                ));
            }
            return Ok(vec![1.0 / poly[0]]);
        }

        let wished_n = 5; // Can be parameterized later

        // Check if leading coefficient is zero
        if poly[0].abs() < 1e-12 {
            return Err(ElementError::GeometryError(
                "Leading coefficient is zero, cannot normalize polynomial".to_string()
            ));
        }

        let a0 = poly[0];
        let inv_a0 = 1.0 / a0;

        let nominal_poly = MonomialPolynomial::multiply_scalar(&poly, inv_a0);
        let mut term = vec![vec![0.0]; wished_n];
        term[0] = vec![1.0]; // Constant term is 1
        term[1] = nominal_poly;
        for i in 2..(wished_n-1) {
            term[i] = MonomialPolynomial::multiply_optimized(&term[i-1], &term[1]);
        }
        let mut result = vec![0.0];
        for j in 0..(wished_n-1) {
            let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
            let signed_term = MonomialPolynomial::multiply_scalar(&term[j], sign);
            result = MonomialPolynomial::add_optimized(&result, &signed_term);
        }
        
        // Multiply by 1/poly[0]
        let result = MonomialPolynomial::multiply_scalar(&result, inv_a0);
        
        Ok(result)
    }


    /// Calculate error based on element type and dispatch to appropriate quadrature rule
    pub fn calculate_error(
        integrand: &Vec<Vec<Vec<f64>>>,
        n: usize,
        element_dim: usize,
    ) -> Result<Vec<Vec<f64>>, GaussError> {

        let matrix_size = integrand.len();
        let mut errors = vec![vec![0.0; matrix_size]; matrix_size];
        
        for i in 0..matrix_size {
            for j in 0..matrix_size {
                match Self::gauss_legendre_error(&integrand[i][j], n, element_dim) {
                    Ok(error) => {
                        errors[i][j] = error;
                    }
                    Err(_) => {
                        errors[i][j] = 0.0; // Use zero error as fallback
                    }
                }
            }
        }
        
        Ok(errors)
    }

    // Gauss-Legendre quadrature error estimation
    fn gauss_legendre_error(
        integrand: &[f64],
        n: usize,
        element_dim: usize,
    ) -> Result<f64, GaussError> {

        // Handle empty polynomials immediately
        if integrand.is_empty() {
            return Ok(0.0);
        }

        let c = factorial(n).powi(4) / ((2 * n + 1) as f64 * factorial(2 * n).powi(3));

        match element_dim {
            1 => {  //1D Gauss-Legendre error function
                let coeff_vector = match MonomialPolynomial::get_coefficients_1d(&integrand) {
                    Ok(vec) => vec,
                    Err(_) => return Ok(0.0),
                };

                let poly_degree = Self::detect_polynomial_order_1d(&coeff_vector);
                let mut s = 0.0;
                if poly_degree >= 2 * n && coeff_vector.len() > 2 * n {
                    let upper = poly_degree.min(coeff_vector.len() - 1);
                    for k in (2 * n)..=upper {
                        let ak = coeff_vector[k].abs();
                        if ak > 1e-12 {
                            s += ak * (factorial(k) / factorial(k - 2 * n));
                        }
                    }
                }
                
                Ok(c * s)
            }
            2 => {  //2D Gauss-Legendre error function
                let coeff_matrix = match MonomialPolynomial::get_coefficients_2d(&integrand) {
                    Ok(matrix) => matrix,
                    Err(_) => return Ok(0.0),
                };

                let (poly_degree_x, poly_degree_y) = Self::detect_polynomial_orders_2d(&coeff_matrix);
                let mut term1 = 0.0;
                for i in (2 * n)..=poly_degree_x.min(coeff_matrix.len().saturating_sub(1)) {
                    for j in 0..=poly_degree_y.min(coeff_matrix[i].len().saturating_sub(1)) {
                        let aij = coeff_matrix[i][j].abs();
                        if aij > 1e-12 {
                            term1 += aij * (factorial(i) / factorial(i - 2 * n));
                        }
                    }
                }

                let mut term2 = 0.0;
                for i in 0..=poly_degree_x.min(coeff_matrix.len().saturating_sub(1)) {
                    for j in (2 * n)..=poly_degree_y.min(coeff_matrix[i].len().saturating_sub(1)) {
                        let aij = coeff_matrix[i][j].abs();
                        if aij > 1e-12 {
                            term2 += aij * (factorial(j) / factorial(j - 2 * n));
                        }
                    }
                }

                Ok(c * (term1 + term2))
            }
            3 => {  //3D Gauss-Legendre error function
                let coeff_tensor = match MonomialPolynomial::get_coefficients_3d(&integrand) {
                    Ok(tensor) => tensor,
                    Err(_) => return Ok(0.0),
                };

                let (poly_degree_x, poly_degree_y, poly_degree_z) = Self::detect_polynomial_orders_3d(&coeff_tensor);
                let mut term1 = 0.0;
                for i in (2 * n)..=poly_degree_x.min(coeff_tensor.len().saturating_sub(1)) {
                    for j in 0..=poly_degree_y.min(coeff_tensor[i].len().saturating_sub(1)) {
                        for k in 0..=poly_degree_z.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                            let aijk = coeff_tensor[i][j][k].abs();
                            if aijk > 1e-12 {
                                term1 += aijk * (factorial(i) / factorial(i - 2 * n));
                            }
                        }
                    }
                }

                let mut term2 = 0.0;
                for i in 0..=poly_degree_x.min(coeff_tensor.len().saturating_sub(1)) {
                    for j in (2 * n)..=poly_degree_y.min(coeff_tensor[i].len().saturating_sub(1)) {
                        for k in 0..=poly_degree_z.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                            let aijk = coeff_tensor[i][j][k].abs();
                            if aijk > 1e-12 {
                                term2 += aijk * (factorial(j) / factorial(j - 2 * n));
                            }
                        }
                    }
                }

                let mut term3 = 0.0;
                for i in 0..=poly_degree_x.min(coeff_tensor.len().saturating_sub(1)) {
                    for j in 0..=poly_degree_y.min(coeff_tensor[i].len().saturating_sub(1)) {
                        for k in (2 * n)..=poly_degree_z.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                            let aijk = coeff_tensor[i][j][k].abs();
                            if aijk > 1e-12 {
                                term3 += aijk * (factorial(k) / factorial(k - 2 * n));
                            }
                        }
                    }
                }

                Ok(c * (term1 + term2 + term3))
            }
            _ => Err(GaussError::GeometryError(format!(
                "Unsupported element dimension {} for Gauss-Legendre error estimation",
                element_dim
            ))),
        }
    }
}



#[cfg(test)]
mod minimal_tests {
    use super::*;
    use crate::structs_and_impls::*;

    // Helper macro for polynomial assertions
    macro_rules! assert_polynomial {
        ($actual:expr, $expected:expr) => {
            for (i, &expected_val) in $expected.iter().enumerate() {
                assert!(
                    ($actual[i] - expected_val).abs() < 1e-10, // Tolerance for floating-point comparison
                    "Coefficient {}: expected {}, got {}",
                    i,
                    expected_val,
                    $actual[i]
                );
            }
        };
    }

    // TEST 1: 1D Linear in 1D mesh - CONSTANT material properties - MASS matrix
    #[test]
    fn test_1d_linear_1d_constant_material() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        let element_order = ElementType::get_element_order(&element_type).unwrap();
        
        // CONSTANT density
        let material_property = MaterialProperty::Scalar(vec![1.0]);

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();

        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_functions.derivatives, 
            element_dim,
            shape_functions.num_nodes,
            1,
            element_order,
        ).unwrap();

        let detj = GaussianQuadrature::calculate_determinant_monomial(
            &jacobian_matrix,
            1,
            element_dim,
        ).unwrap();

        assert!((jacobian_matrix[0][0][0] - 2.0).abs() < 1e-10);
        assert!((detj[0] - 2.0).abs() < 1e-10);
        
        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &element_nodes,
            &material_property,
            1,
        ).unwrap();
        
        let m00 = &integrand[0][0];     // 2 - 4*xi + 2*xi^2
        let m01 = &integrand[0][1];     // 2*xi - 2*xi^2
        let m11 = &integrand[1][1];     // 2*xi^2

        let m10 = &integrand[1][0];     // 2*xi - 2*xi^2
        

        assert_polynomial!(m00, &[2.0, -4.0, 0.0, 0.0, 2.0]);
        assert_polynomial!(m01, &[0.0, 2.0, 0.0, 0.0, -2.0]);
        assert_polynomial!(m11, &[0.0, 0.0, 0.0, 0.0, 2.0]);
        
        // Test symmetry 
        assert_polynomial!(m01, m10);
    }

    // TEST 2: 1D Linear in 1D mesh - VARIABLE material properties - MASS matrix
    #[test]
    fn test_1d_linear_1d_variable_material() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        let mesh_dim = 1;
        
        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();

        // VARIABLE density: rho(x) = 1 + x
        // In physical space: x(xi) = 2*xi, so rho(xi) = 1 + 2*xi
        let material_property_physical = MaterialProperty::Scalar(vec![1.0, 1.0, 0.0, 0.0]); // 1 + x
        
        let coordinate_maps = GaussianQuadrature::build_coordinate_maps(&element_nodes, &shape_functions.values, mesh_dim).unwrap(); // x(xi) = 2*xi
        assert_polynomial!(coordinate_maps[0], &[0.0, 2.0, 0.0, 0.0]);

        let rho_isoparametric = GaussianQuadrature::substitute_polynomial(material_property_physical.as_scalar().unwrap(), &coordinate_maps).unwrap(); // rho(xi) = 1 + 2*xi
        assert_polynomial!(rho_isoparametric, &[1.0, 2.0, 0.0, 0.0]);

        let material_property_isoparametric = GaussianQuadrature::transform_material_property_scalar(&material_property_physical.as_scalar().unwrap(), &element_nodes, &shape_functions.values, mesh_dim).unwrap();
        assert_polynomial!(material_property_isoparametric, &[1.0, 2.0, 0.0, 0.0]);
        

        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &element_nodes,
            &material_property_physical,
            mesh_dim,
        ).unwrap();

        
        let m00 = &integrand[0][0];     // 2 - 6*xi^2 + 4*xi^3 
        let m01 = &integrand[0][1];     // 2*xi + 2*xi^2 - 4*xi^3
        let m11 = &integrand[1][1];     // 2*xi^2 + 4*xi^3

        let m10 = &integrand[1][0];     // 2*xi + 2*xi^2 - 4*xi^3
        
        assert_polynomial!(m00, &[2.0, 0.0, 0.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0]);
        assert_polynomial!(m01, &[0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0]);
        assert_polynomial!(m11, &[0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0]);
        assert_polynomial!(m01, m10);
    } 

    // TEST 3: 1D Linear in 1D mesh - CONSTANT material - STIFFNESS matrix
    #[test]
    fn test_1d_stiffness_constant_material() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        let element_order = ElementType::get_element_order(&element_type).unwrap();
        let mesh_dim = 1;
        

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();
        
        // CONSTANT elasticity: C = 1.0
        let material_property = MaterialProperty::Matrix(vec![vec![vec![1.0]]]);
        

        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_functions.derivatives, 
            element_dim,
            shape_functions.num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();
        assert!((jacobian_matrix[0][0][0] - 2.0).abs() < 1e-10);
        
        let inv_jacobian = GaussianQuadrature::inverse_matrix_monomial(
            &jacobian_matrix, 
            element_dim, 
            mesh_dim,
        ).unwrap();
        assert!((inv_jacobian[0][0][0] - 0.5).abs() < 1e-10);

        let b_matrix = GaussianQuadrature::build_b_matrix(&shape_functions.derivatives, &inv_jacobian, shape_functions.num_nodes, mesh_dim, element_dim).unwrap();
        assert!((b_matrix[0][0][0] - (-0.5)).abs() < 1e-10);
        assert!((b_matrix[0][1][0] - 0.5).abs() < 1e-10);

        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Stiffness,
            &element_type,
            &element_nodes,
            &material_property,
            1,
        ).unwrap();
        
        let k00 = &integrand[0][0];     // 0.5
        let k01 = &integrand[0][1];     // -0.5
        let k11 = &integrand[1][1];     // 0.5

        let k10 = &integrand[1][0];     // -0.5
        
        assert!((k00[0] - 0.5).abs() < 1e-10);
        assert!((k01[0] - (-0.5)).abs() < 1e-10);
        assert!((k11[0] - 0.5).abs() < 1e-10);

        assert!((k01[0] - k10[0]).abs() < 1e-10);
        

    }

    // TEST 4: 1D Linear in 1D mesh - VARIABLE material - STIFFNESS matrix
    #[test]
    fn test_1d_stiffness_variable_material() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        let mesh_dim = 1;

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();

        // VARIABLE elasticity: C(x) = 5 - x
        // In physical space: x(xi) = 2*xi, so C(xi) = 5 - 2*xi
        let material_property_physical = MaterialProperty::Matrix(vec![vec![vec![5.0, -1.0, 0.0, 0.0]]]); // 5 - x
        
        let coordinate_maps = GaussianQuadrature::build_coordinate_maps(&element_nodes, &shape_functions.values, mesh_dim).unwrap(); // x(xi) = 2*xi

        let c_isoparametric = GaussianQuadrature::substitute_polynomial(&material_property_physical.as_matrix().unwrap()[0][0], &coordinate_maps).unwrap(); // C(xi) = 5 - 2*xi
        assert_polynomial!(c_isoparametric, &[5.0, -2.0, 0.0, 0.0]);

        let material_property_isoparametric = GaussianQuadrature::transform_material_property_tensor(&material_property_physical.as_matrix().unwrap(), &element_nodes, &shape_functions.values, mesh_dim).unwrap();
        assert_polynomial!(material_property_isoparametric[0][0], &[5.0, -2.0, 0.0, 0.0]);
        
        
        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Stiffness,
            &element_type,
            &element_nodes,
            &material_property_physical,
            mesh_dim,
        ).unwrap();

        
        let k00 = &integrand[0][0];     // 2.5 - xi
        let k01 = &integrand[0][1];     // -2.5 + xi
        let k11 = &integrand[1][1];     // 2.5 - xi

        let k10 = &integrand[1][0];     // -2.5 + xi
        
        
        assert_polynomial!(k00, &[2.5, -1.0, 0.0, 0.0]);
        assert_polynomial!(k01, &[-2.5, 1.0, 0.0, 0.0]);
        assert_polynomial!(k11, &[2.5, -1.0, 0.0, 0.0]);

        assert_polynomial!(k01, k10);
    }

    // TEST 5: 2D Triangle in 2D mesh - CONSTANT material - MASS matrix
    #[test]
    fn test_2d_triangle_2d_constant_material() {
        let element = Element { id: 0, nodes: vec![0, 1, 2] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-1.0, -1.0] },
            Node { id: 1, coordinates: vec![2.0, -1.0] },
            Node { id: 2, coordinates: vec![-1.0, 3.0] },
        ];
        let element_type = ElementType::Triangle;
        let mesh_dim = 2;

        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        let element_order = ElementType::get_element_order(&element_type).unwrap();
        
        // CONSTANT density
        let material_property = MaterialProperty::Scalar(vec![2.0]);

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();

        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_functions.derivatives, 
            element_dim,
            shape_functions.num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        let detj = GaussianQuadrature::calculate_determinant_monomial(
            &jacobian_matrix,
            mesh_dim,
            element_dim,
        ).unwrap();

        assert!((jacobian_matrix[0][0][0] - 3.0).abs() < 1e-10);
        assert!((jacobian_matrix[0][1][0] - 0.0).abs() < 1e-10);
        assert!((jacobian_matrix[1][0][0] - 0.0).abs() < 1e-10);
        assert!((jacobian_matrix[1][1][0] - 4.0).abs() < 1e-10);

        assert!((detj[0] - 12.0).abs() < 1e-10);
        
        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &element_nodes,
            &material_property,
            mesh_dim,
        ).unwrap();
             
        
        // Test non-zero components 
        assert_polynomial!(integrand[0][0], &[24.0, -48.0, -48.0, 0.0, 24.0, 48.0, 0.0, 24.0, 0.0, 0.0]); // 24 - 48*xi - 48*eta + 24*xi^2 + 48*xi*eta + 24*eta^2
        assert_polynomial!(integrand[0][2], &[0.0, 24.0, 0.0, 0.0, -24.0, -24.0, 0.0, 0.0, 0.0, 0.0]); // 24*xi - 24*xi^2 - 24*xi*eta
        assert_polynomial!(integrand[0][4], &[0.0, 0.0, 24.0, 0.0, 0.0, -24.0, 0.0, -24.0, 0.0, 0.0]); // 24*eta - 24*xi*eta - 24*eta^2
        assert_polynomial!(integrand[2][2], &[0.0, 0.0, 0.0, 0.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0]); // 24*xi^2
        assert_polynomial!(integrand[2][4], &[0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 0.0, 0.0, 0.0, 0.0]); // 24*xi*eta
        assert_polynomial!(integrand[4][4], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 0.0, 0.0]); // 24*eta^2
        
        // Test zero components 
        let zero_indices = [(0,1), (0,3), (0,5), (1,2), (1,4), (1,5), (2,3), (2,5), (3,4), (4,5)];
        for (i, j) in zero_indices {
            assert_polynomial!(integrand[i][j], &[0.0]);
        }
        
        // Test symmetry 
        assert_polynomial!(integrand[1][1], integrand[0][0]); 
        assert_polynomial!(integrand[3][3], integrand[2][2]);
        assert_polynomial!(integrand[5][5], integrand[4][4]);
        
        let symmetric_pairs = [
        (1,0), (2,0), (2,1), (3,0), (3,1), (3,2), 
        (4,0), (4,1), (4,2), (4,3), (5,0), (5,1), 
        (5,2), (5,3), (5,4)
        ];
        for (i, j) in symmetric_pairs {
            assert_polynomial!(integrand[i][j], integrand[j][i]);
        }
    }

    // TEST 6: 2D Triangle in 2D mesh - VARIABLE material - MASS matrix
    #[test]
    fn test_2d_triangle_2d_variable_material() {
        let element = Element { id: 0, nodes: vec![0, 1, 2] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-1.0, -1.0] },
            Node { id: 1, coordinates: vec![2.0, -1.0] },
            Node { id: 2, coordinates: vec![-1.0, 3.0] },
        ];
        let element_type = ElementType::Triangle;
        let mesh_dim = 2;

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();

        // VARIABLE density: rho(x) = 10 - 2*x + y
        // In physical space: x(xi,eta) = -1 + 3*xi and y(xi,eta) = -1 + 4*eta, so rho(xi,eta) = 11 - 6*xi + 4*eta
        let material_property_physical = MaterialProperty::Scalar(vec![10.0, -2.0, 1.0, 0.0]); // 10 - 2*x + y
        
        let coordinate_maps = GaussianQuadrature::build_coordinate_maps(&element_nodes, &shape_functions.values, mesh_dim).unwrap(); // x(xi) = 2*xi
        assert_polynomial!(coordinate_maps[0], &[-1.0, 3.0, 0.0, 0.0]);
        assert_polynomial!(coordinate_maps[1], &[-1.0, 0.0, 4.0, 0.0]);

        let rho_isoparametric = GaussianQuadrature::substitute_polynomial(material_property_physical.as_scalar().unwrap(), &coordinate_maps).unwrap(); // ρ(xi) = 1 + 2*xi
        assert_polynomial!(rho_isoparametric, &[11.0, -6.0, 4.0, 0.0]);

        let material_property_isoparametric = GaussianQuadrature::transform_material_property_scalar(&material_property_physical.as_scalar().unwrap(), &element_nodes, &shape_functions.values, mesh_dim).unwrap();
        assert_polynomial!(material_property_isoparametric, &[11.0, -6.0, 4.0, 0.0]);
        
        
        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &element_nodes,
            &material_property_physical,
            mesh_dim,
        ).unwrap();


        // Test non-zero components 

        // m00 = 48*eta^3 + 24*eta^2 *xi + 36*eta^2 - 96*eta*xi^2 + 312*eta*xi - 216*eta - 72*xi^3 + 276*xi^2 - 336*xi + 132
        assert_polynomial!(integrand[0][0], &[132.0, -336.0, -216.0, 0.0, 276.0, 312.0, 0.0, 36.0, 0.0, 0.0, -72.0, -96.0, 0.0, 24.0, 0.0, 0.0, 48.0, 0.0, 0.0, 0.0]);
        // m02 = -48*eta^2 *xi + 24*eta*xi^2 - 84*eta*xi + 72*xi^3 - 204*xi^2 + 132*xi
        assert_polynomial!(integrand[0][2], &[0.0, 132.0, 0.0, 0.0, -204.0, -84.0, 0.0, 0.0, 0.0, 0.0, 72.0, 24.0, 0.0, -48.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // m04 = -48*eta^3 + 24*eta^2 *xi - 84*eta^2 + 72*eta*xi^2 - 204*eta*xi + 132*eta
        assert_polynomial!(integrand[0][4], &[0.0, 0.0, 132.0, 0.0, 0.0, -204.0, 0.0, -84.0, 0.0, 0.0, 0.0, 72.0, 0.0, 24.0, 0.0, 0.0, -48.0, 0.0, 0.0, 0.0]);
        // m22 = 48*eta*xi^2 - 72*xi^3 + 132*xi^2
        assert_polynomial!(integrand[2][2], &[0.0, 0.0, 0.0, 0.0, 132.0, 0.0, 0.0, 0.0, 0.0, 0.0, -72.0, 48.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // m24 = 48*eta^2 *xi - 72*eta*xi^2 + 132*eta*xi
        assert_polynomial!(integrand[2][4], &[0.0, 0.0, 0.0, 0.0, 0.0, 132.0, 0.0, 0.0, 0.0, 0.0, 0.0, -72.0, 0.0, 48.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // m44 = 48*eta^3 - 72*eta^2 *xi + 132*eta^2
        assert_polynomial!(integrand[4][4], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 132.0, 0.0, 0.0, 0.0, 0.0, 0.0, -72.0, 0.0, 0.0, 48.0, 0.0, 0.0, 0.0]);
        
        // Test zero components 
        let zero_indices = [(0,1), (0,3), (0,5), (1,2), (1,4), (1,5), (2,3), (2,5), (3,4), (4,5)];
        for (i, j) in zero_indices {
            assert_polynomial!(integrand[i][j], &[0.0]);
        }
        
        // Test symmetry patterns
        assert_polynomial!(integrand[1][1], integrand[0][0]);
        assert_polynomial!(integrand[3][3], integrand[2][2]);
        assert_polynomial!(integrand[5][5], integrand[4][4]);
        
        // Test symmetry for remaining pairs
        let symmetric_pairs = [
            (1,0), (2,0), (2,1), (3,0), (3,1), (3,2), 
            (4,0), (4,1), (4,2), (4,3), (5,0), (5,1), 
            (5,2), (5,3), (5,4)
        ];
        for (i, j) in symmetric_pairs {
            assert_polynomial!(integrand[i][j], integrand[j][i]);
        }
    }
 
    // TEST 7: 2D Triangle in 2D mesh - VARIABLE material - STIFFNESS matrix
    #[test]
    fn test_2d_stiffness_variable_anisotropic() {
        let element = Element { id: 0, nodes: vec![0, 1, 2] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-1.0, -1.0] },
            Node { id: 1, coordinates: vec![2.0, -1.0] },
            Node { id: 2, coordinates: vec![-1.0, 3.0] },
        ];
        let element_type = ElementType::Triangle;
        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        let element_order = ElementType::get_element_order(&element_type).unwrap();
        let mesh_dim = 2;

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();
        
        // VARIABLE anisotropic material: C(x,y) 
        let c_matrix = vec![
            vec![vec![1.0, 1.0, 0.0, 0.0], vec![0.5, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 0.0]], // C11 = 1+x, C12 = 0.5, C13 = 0
            vec![vec![0.5, 0.0, 1.0, 0.0], vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 0.0]], // C21 = 0.5+y, C22 = 1.0, C23 = 0
            vec![vec![0.0, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 0.0], vec![0.5, 1.0, 1.0, 0.0]], // C31 = 0, C32 = 0, C33 = 0.5+x+y
        ];

        let material_property_physical = MaterialProperty::Matrix(c_matrix);


        let material_property_isoparametric = GaussianQuadrature::transform_material_property_tensor(&material_property_physical.as_matrix().unwrap(), &element_nodes, &shape_functions.values, mesh_dim).unwrap();
        assert_polynomial!(material_property_isoparametric[0][0], &[0.0, 3.0, 0.0, 0.0]);
        assert_polynomial!(material_property_isoparametric[0][1], &[0.5]);
        assert_polynomial!(material_property_isoparametric[0][2], &[0.0]);
        assert_polynomial!(material_property_isoparametric[1][0], &[-0.5, 0.0, 4.0, 0.0]);
        assert_polynomial!(material_property_isoparametric[1][1], &[1.0]);
        assert_polynomial!(material_property_isoparametric[1][2], &[0.0]);
        assert_polynomial!(material_property_isoparametric[2][0], &[0.0]);
        assert_polynomial!(material_property_isoparametric[2][1], &[0.0]);
        assert_polynomial!(material_property_isoparametric[2][2], &[-1.5, 3.0, 4.0, 0.0]);


        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_functions.derivatives, 
            element_dim,
            shape_functions.num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        let adjoint_jacobian = GaussianQuadrature::calculate_adjoint_monomial(
            &jacobian_matrix, 
            element_dim
        ).unwrap();
        assert_polynomial!(adjoint_jacobian[0][0], &[4.0]);
        assert_polynomial!(adjoint_jacobian[0][1], &[0.0]);
        assert_polynomial!(adjoint_jacobian[1][0], &[0.0]);
        assert_polynomial!(adjoint_jacobian[1][1], &[3.0]);
        
        
        let inv_jacobian = GaussianQuadrature::inverse_matrix_monomial(
            &jacobian_matrix, 
            mesh_dim, 
            element_dim,
        ).unwrap();
        assert_polynomial!(inv_jacobian[0][0], &[0.3333333333]);
        assert_polynomial!(inv_jacobian[0][1], &[0.0]);
        assert_polynomial!(inv_jacobian[1][0], &[0.0]);
        assert_polynomial!(inv_jacobian[1][1], &[0.25]);


        let b_matrix = GaussianQuadrature::build_b_matrix(&shape_functions.derivatives, &inv_jacobian, shape_functions.num_nodes, mesh_dim, element_dim).unwrap();
        assert_polynomial!(b_matrix[0][0], &[-0.3333333333]);
        assert_polynomial!(b_matrix[0][1], &[0.0]);
        assert_polynomial!(b_matrix[0][2], &[0.3333333333]);
        assert_polynomial!(b_matrix[0][3], &[0.0]);
        assert_polynomial!(b_matrix[0][4], &[0.0]);
        assert_polynomial!(b_matrix[0][5], &[0.0]);

        assert_polynomial!(b_matrix[1][0], &[0.0]);
        assert_polynomial!(b_matrix[1][1], &[-0.25]);
        assert_polynomial!(b_matrix[1][2], &[0.0]);
        assert_polynomial!(b_matrix[1][3], &[0.0]);
        assert_polynomial!(b_matrix[1][4], &[0.0]);
        assert_polynomial!(b_matrix[1][5], &[0.25]);

        assert_polynomial!(b_matrix[2][0], &[-0.25]);
        assert_polynomial!(b_matrix[2][1], &[-0.3333333333]);
        assert_polynomial!(b_matrix[2][2], &[0.0]);
        assert_polynomial!(b_matrix[2][3], &[0.3333333333]);
        assert_polynomial!(b_matrix[2][4], &[0.25]);
        assert_polynomial!(b_matrix[2][5], &[0.0]);

        
        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Stiffness,
            &element_type,
            &nodes,
            &material_property_physical,
            mesh_dim,
        ).unwrap();

        
        // Test stiffness matrix components
        assert_polynomial!(integrand[0][0], &[-1.125, 6.25, 3.0, 0.0]); // k00 = 12 * [ eta/4 + (25/48)*xi - 3/32 ] (12 is detJ) 
        assert_polynomial!(integrand[0][1], &[-1.0, 3.0, 4.0, 0.0]);    // k01 = 12 * [ eta/3 + xi/4 - 1/12 ]
        assert_polynomial!(integrand[0][2], &[0.0, -4.0, 0.0, 0.0]);    // k02 = 12 * [ -xi/3 ]
        assert_polynomial!(integrand[0][3], &[1.5, -3.0, -4.0, 0.0]);   // k03 = 12 * [ -eta/3 - xi/4 + 1/8 ]
        assert_polynomial!(integrand[0][4], &[1.125, -2.25, -3.0, 0.0]);    // k04 = 12 * [ -eta/4 - (3/16)*xi + 3/32 ]
        assert_polynomial!(integrand[0][5], &[-0.5]);   // k05 = 12 * [ -1/24 ]

        assert_polynomial!(integrand[1][0], &[-2.0, 3.0, 8.0, 0.0]);    // k10 = 12 * [ (2/3)*eta + xi/4 - 1/6 ]
        assert_polynomial!(integrand[1][1], &[-1.25, 4.0, 5.3333333333, 0.0]);  // k11 = 12 * [ (4/9)*eta + xi/3 - 5/48 ]
        assert_polynomial!(integrand[1][2], &[0.5, 0.0, -4.0, 0.0]);    // k12 = 12 * [ 1/24 - eta/3 ]
        assert_polynomial!(integrand[1][3], &[2.0, -4.0, -5.3333333333, 0.0]);  // k13 = 12 * [ -(4/9)*eta - xi/3 + 1/6 ]
        assert_polynomial!(integrand[1][4], &[1.5, -3.0, -4.0, 0.0]);   // k14 = 12 * [ -eta/3 - xi/4 + 1/8 ]
        assert_polynomial!(integrand[1][5], &[-0.75]);  // k15 = 12 * [ -1/16 ]

        assert_polynomial!(integrand[2][0], &[0.0, -4.0, 0.0, 0.0]);    // k20 = 12 * [ -xi/3 ]
        assert_polynomial!(integrand[2][1], &[-0.5]);   // k21 = 12 * [ -1/24 ]
        assert_polynomial!(integrand[2][2], &[0.0, 4.0, 0.0, 0.0]); // k22 = 12 * [ xi/3 ]
        assert_polynomial!(integrand[2][3], &[0.0]);    // k23 = 0
        assert_polynomial!(integrand[2][4], &[0.0]);    // k24 = 0
        assert_polynomial!(integrand[2][5], &[0.5]);    // k25 = 12 * [ 1/24 ]
        
        assert_polynomial!(integrand[3][0], &[1.5, -3.0, -4.0, 0.0]);   // k30 = 12 * [ -eta/3 - xi/4 + 1/8 ]
        assert_polynomial!(integrand[3][1], &[2.0, -4.0, -5.3333333333, 0.0]);  // k31 = 12 * [ -(4/9)*eta - xi/3 + 1/6 ]
        assert_polynomial!(integrand[3][2], &[0.0]);    // k32 = 0  
        assert_polynomial!(integrand[3][3], &[-2.0, 4.0, 5.3333333333, 0.0]);   // k33 = 12 * [ (4/9)*eta + xi/3 - 1/6 ]
        assert_polynomial!(integrand[3][4], &[-1.5, 3.0, 4.0, 0.0]);    // k34 = 12 * [ eta/3 + xi/4 - 1/8 ]  
        assert_polynomial!(integrand[3][5], &[0.0]);    // k35 = 0

        assert_polynomial!(integrand[4][0], &[1.125, -2.25, -3.0, 0.0]);    // k40 = 12 * [ -eta/4 - (3/16)*xi + 3/32 ]
        assert_polynomial!(integrand[4][1], &[1.5, -3.0, -4.0, 0.0]);   // k41 = 12 * [ -eta/3 - xi/4 + 1/8 ] 
        assert_polynomial!(integrand[4][2], &[0.0]);    // k42 = 0
        assert_polynomial!(integrand[4][3], &[-1.5, 3.0, 4.0, 0.0]);    // k43 = 12 * [ eta/3 + xi/4 - 1/8 ]  
        assert_polynomial!(integrand[4][4], &[-1.125, 2.25, 3.0, 0.0]); // k44 = 12 * [ eta/4 + (3/16)*xi - 3/32 ]  
        assert_polynomial!(integrand[4][5], &[0.0]);    // k45 = 0

        assert_polynomial!(integrand[5][0], &[0.5, 0.0, -4.0, 0.0]);    // k50 = 12 * [ 1/24 - eta/3] 
        assert_polynomial!(integrand[5][1], &[-0.75]);  // k51 = 12 * [ -1/16 ] 
        assert_polynomial!(integrand[5][2], &[-0.5, 0.0, 4.0, 0.0]);    // k52 = 12 * [ -1/24 + eta/3 ]  
        assert_polynomial!(integrand[5][3], &[0.0]);    // k53 = 0 
        assert_polynomial!(integrand[5][4], &[0.0]);    // k54 = 0 
        assert_polynomial!(integrand[5][5], &[0.75]);   // k55 = 12 * [ 1/16 ] 
    } 

    // TEST 8: 3D Pyramid in 3D mesh - CONSTANT material - MASS matrix
    #[test]
    fn test_3d_pyramid_constant_mass() {
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-2.0, -1.0, 0.0] },  // base corner 1
            Node { id: 1, coordinates: vec![1.0, -1.5, 0.0] },   // base corner 2  
            Node { id: 2, coordinates: vec![1.0, 1.0, 0.0] },    // base corner 3
            Node { id: 3, coordinates: vec![-1.0, 1.0, 0.0] },   // base corner 4
            Node { id: 4, coordinates: vec![0.5, 0.5, 2.0] },    // apex
        ];
        let element_type = ElementType::Pyramid;
        let mesh_dim = 3;

        let element_dim = ElementType::get_element_dimension(&element_type).unwrap();
        let element_order = ElementType::get_element_order(&element_type).unwrap();
        
        // CONSTANT density
        let material_property = MaterialProperty::Scalar(vec![2.0]);

        let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
        let shape_functions = ElementType::get_shape_functions(
            &element_type,
        ).unwrap();
        println!("Shape functions: {:?}", shape_functions.values);
        println!("Derivative of shape functions: {:?}", shape_functions.derivatives);

        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_functions.derivatives, 
            element_dim,
            shape_functions.num_nodes,
            mesh_dim,
            element_order,
        ).unwrap();

        let detj = GaussianQuadrature::calculate_determinant_monomial(
            &jacobian_matrix,
            mesh_dim,
            element_dim,
        ).unwrap();
        
        println!("The Jacobian matrix: {:?}", jacobian_matrix);
        println!("Determinant of the Jacobian matrix: {:?}", detj);

        // Test Jacobian matrix components
        assert_polynomial!(jacobian_matrix[0][0], &[3.0, 0.0, -1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
        assert_polynomial!(jacobian_matrix[0][1], &[-0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0]);
        assert_polynomial!(jacobian_matrix[0][2], &[0.0]);

        assert_polynomial!(jacobian_matrix[1][0], &[1.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
        assert_polynomial!(jacobian_matrix[1][1], &[2.0, 0.5, 0.0, -2.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0]);
        assert_polynomial!(jacobian_matrix[1][2], &[0.0]);

        assert_polynomial!(jacobian_matrix[2][0], &[2.5, -3.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        assert_polynomial!(jacobian_matrix[2][1], &[1.5, 0.5, -2.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0]);
        assert_polynomial!(jacobian_matrix[2][2], &[2.0]);

        // Test Jacobian determinant
        assert_polynomial!(detj, &[13.0, 2.0, -5.0, -26.0, 0.0, 0.0, -4.0, 0.0, 10.0, 13.0,
                                                    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, -5.0, 0.0]);

        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &element_nodes,
            &material_property,
            mesh_dim,
        ).unwrap();

        print!("Integrand00: {:?}", integrand[0][0]);
        print!("Integrand01: {:?}", integrand[0][1]);
        print!("Integrand02: {:?}", integrand[0][2]);
        print!("Integrand03: {:?}", integrand[0][3]);
        // Test non-zero components 
        assert_polynomial!(integrand[0][0], &[26.0, -48.0, -62.0, -104.0, 18.0, 116.0, 192.0, 46.0, 248.0, 156.0, 4.0, -46.0, -72.0, -88.0, -464.0, -288.0, 
                                              -10.0, -184.0, -372.0, -104.0, 0.0, -8.0, -16.0, 38.0, 184.0, 108.0, 20.0, 352.0, 696.0, 192.0, 0.0, 40.0, 
                                              276.0, 248.0, 26.0, 0.0, 0.0, 0.0, 4.0, 32.0, 24.0, -10.0, -152.0, -276.0, -72.0, 0.0, -80.0, -528.0, -464.0, 
                                              -48.0, 0.0, 0.0, -60.0, -184.0, -62.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -48.0, -16.0, 0.0, 40.0, 228.0, 
                                              184.0, 18.0, 0.0, 0.0, 120.0, 352.0, 116.0, 0.0, 0.0, 0.0, 40.0, 46.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 24.0, 32.0, 4.0, 0.0, 0.0, -60.0, -152.0, -46.0, 0.0, 0.0, 0.0, -80.0, -88.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -8.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 40.0, 38.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); //20.0*eta^3*psi^4*xi - 10.0*eta^3*psi^4*xi^2 - 10.0*eta^3*psi^4 
                                              // + 40.0*eta^3*psi^3*xi^2 - 80.0*eta^3*psi^3*xi + 40.0*eta^3*psi^3 - 60.0*eta^3*psi^2*xi^2 + 120.0*eta^3*psi^2*xi 
                                              // - 60.0*eta^3*psi^2 + 40.0*eta^3*psi*xi^2 - 80.0*eta^3*psi*xi + 40.0*eta^3*psi - 10.0*eta^3*xi^2 + 20.0*eta^3*xi 
                                              // - 10.0*eta^3 + 4.0*eta^2*psi^4*xi^3 + 38.0*eta^2*psi^4*xi^2 - 88.0*eta^2*psi^4*xi + 46.0*eta^2*psi^4 
                                              // - 16.0*eta^2*psi^3*xi^3 - 152.0*eta^2*psi^3*xi^2 + 352.0*eta^2*psi^3*xi - 184.0*eta^2*psi^3 + 24.0*eta^2*psi^2*xi^3 
                                              // + 228.0*eta^2*psi^2*xi^2 - 528.0*eta^2*psi^2*xi + 276.0*eta^2*psi^2 - 16.0*eta^2*psi*xi^3 - 152.0*eta^2*psi*xi^2 
                                              // + 352.0*eta^2*psi*xi - 184.0*eta^2*psi + 4.0*eta^2*xi^3 + 38.0*eta^2*xi^2 - 88.0*eta^2*xi + 46.0*eta^2 
                                              // - 8.0*eta*psi^4*xi^3 - 46.0*eta*psi^4*xi^2 + 116.0*eta*psi^4*xi - 62.0*eta*psi^4 + 32.0*eta*psi^3*xi^3 
                                              // + 184.0*eta*psi^3*xi^2 - 464.0*eta*psi^3*xi + 248.0*eta*psi^3 - 48.0*eta*psi^2*xi^3 - 276.0*eta*psi^2*xi^2 
                                              // + 696.0*eta*psi^2*xi - 372.0*eta*psi^2 + 32.0*eta*psi*xi^3 + 184.0*eta*psi*xi^2 - 464.0*eta*psi*xi + 248.0*eta*psi 
                                              // - 8.0*eta*xi^3 - 46.0*eta*xi^2 + 116.0*eta*xi - 62.0*eta + 4.0*psi^4*xi^3 + 18.0*psi^4*xi^2 - 48.0*psi^4*xi 
                                              // + 26.0*psi^4 - 16.0*psi^3*xi^3 - 72.0*psi^3*xi^2 + 192.0*psi^3*xi - 104.0*psi^3 + 24.0*psi^2*xi^3 + 108.0*psi^2*xi^2 
                                              // - 288.0*psi^2*xi + 156.0*psi^2 - 16.0*psi*xi^3 - 72.0*psi*xi^2 + 192.0*psi*xi - 104.0*psi + 4.0*xi^3 + 18.0*xi^2 
                                              // - 48.0*xi + 26.0
  
        assert_polynomial!(integrand[0][3], &[0.0, 26.0, 0.0, 0.0, -22.0, -62.0, -104.0, 0.0, 248.0, 156.0, -4.0, 54.0, 88.0, 46.0, -216.0, -132.0, 0.0, 
                                              168.0, 324.0, 88.0, 0.0, 8.0, 16.0, -42.0, -216.0, -132.0, -10.0, -184.0, -372.0, -104.0, 0.0, -32.0, -252.0, 
                                              -184.0, 26.0, 0.0, 0.0, 0.0, -4.0, 16.0, 24.0, 10.0, 168.0, 324.0, 88.0, 0.0, 40.0, 276.0, 248.0, 26.0, 0.0, 
                                              0.0, 60.0, 168.0, 54.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 48.0, 16.0, 0.0, -40.0, -252.0, -216.0, 
                                              -22.0, 0.0, 0.0, -60.0, -184.0, -62.0, 0.0, 0.0, 0.0, -40.0, -42.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, -24.0, -32.0, -4.0, 0.0, 0.0, 60.0, 168.0, 54.0, 0.0, 0.0, 0.0, 40.0, 46.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 
                                              8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -40.0, -42.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); // 10.0*eta^3*psi^4*xi^2 
                                              // - 10.0*eta^3*psi^4*xi - 40.0*eta^3*psi^3*xi^2 + 40.0*eta^3*psi^3*xi + 60.0*eta^3*psi^2*xi^2 
                                              // - 60.0*eta^3*psi^2*xi - 40.0*eta^3*psi*xi^2 + 40.0*eta^3*psi*xi + 10.0*eta^3*xi^2 - 10.0*eta^3*xi 
                                              // - 4.0*eta^2*psi^4*xi^3 - 42.0*eta^2*psi^4*xi^2 + 46.0*eta^2*psi^4*xi + 16.0*eta^2*psi^3*xi^3 
                                              // + 168.0*eta^2*psi^3*xi^2 - 184.0*eta^2*psi^3*xi - 24.0*eta^2*psi^2*xi^3 - 252.0*eta^2*psi^2*xi^2 
                                              // + 276.0*eta^2*psi^2*xi + 16.0*eta^2*psi*xi^3 + 168.0*eta^2*psi*xi^2 - 184.0*eta^2*psi*xi - 4.0*eta^2*xi^3 
                                              // - 42.0*eta^2*xi^2 + 46.0*eta^2*xi + 8.0*eta*psi^4*xi^3 + 54.0*eta*psi^4*xi^2 - 62.0*eta*psi^4*xi 
                                              // - 32.0*eta*psi^3*xi^3 - 216.0*eta*psi^3*xi^2 + 248.0*eta*psi^3*xi + 48.0*eta*psi^2*xi^3 + 324.0*eta*psi^2*xi^2 
                                              // - 372.0*eta*psi^2*xi - 32.0*eta*psi*xi^3 - 216.0*eta*psi*xi^2 + 248.0*eta*psi*xi + 8.0*eta*xi^3 + 54.0*eta*xi^2 
                                              // - 62.0*eta*xi - 4.0*psi^4*xi^3 - 22.0*psi^4*xi^2 + 26.0*psi^4*xi + 16.0*psi^3*xi^3 + 88.0*psi^3*xi^2 
                                              // - 104.0*psi^3*xi - 24.0*psi^2*xi^3 - 132.0*psi^2*xi^2 + 156.0*psi^2*xi + 16.0*psi*xi^3 + 88.0*psi*xi^2 
                                              // - 104.0*psi*xi - 4.0*xi^3 - 22.0*xi^2 + 26.0*xi
 
        assert_polynomial!(integrand[0][6], &[0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, -104.0, 0.0, 0.0, -22.0, 0.0, -36.0, 88.0, 0.0, 10.0, 144.0, -132.0, 
                                              0.0, 0.0, -4.0, 0.0, 32.0, 88.0, 0.0, -10.0, -128.0, 156.0, 0.0, 0.0, 16.0, 192.0, 144.0, 26.0, 0.0, 0.0, 
                                              0.0, 4.0, -16.0, -24.0, 0.0, -128.0, -132.0, 0.0, 0.0, -40.0, -216.0, -104.0, 0.0, 0.0, 0.0, -60.0, -128.0, 
                                              -22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -24.0, 0.0, 0.0, 40.0, 192.0, 88.0, 0.0, 0.0, 0.0, 60.0, 
                                              144.0, 26.0, 0.0, 0.0, 0.0, 40.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 16.0, 
                                              -4.0, 0.0, 0.0, -60.0, -128.0, -22.0, 0.0, 0.0, 0.0, -40.0, -36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              40.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); // 10.0*eta^3*psi^4*xi - 10.0*eta^3*psi^4*xi^2 + 40.0*eta^3*psi^3*xi^2 
                                              // - 40.0*eta^3*psi^3*xi - 60.0*eta^3*psi^2*xi^2 + 60.0*eta^3*psi^2*xi + 40.0*eta^3*psi*xi^2 - 40.0*eta^3*psi*xi 
                                              // - 10.0*eta^3*xi^2 + 10.0*eta^3*xi + 4.0*eta^2*psi^4*xi^3 + 32.0*eta^2*psi^4*xi^2 - 36.0*eta^2*psi^4*xi 
                                              // - 16.0*eta^2*psi^3*xi^3 - 128.0*eta^2*psi^3*xi^2 + 144.0*eta^2*psi^3*xi + 24.0*eta^2*psi^2*xi^3 
                                              // + 192.0*eta^2*psi^2*xi^2 - 216.0*eta^2*psi^2*xi - 16.0*eta^2*psi*xi^3 - 128.0*eta^2*psi*xi^2 + 144.0*eta^2*psi*xi 
                                              // + 4.0*eta^2*xi^3 + 32.0*eta^2*xi^2 - 36.0*eta^2*xi - 4.0*eta*psi^4*xi^3 - 22.0*eta*psi^4*xi^2 + 26.0*eta*psi^4*xi 
                                              // + 16.0*eta*psi^3*xi^3 + 88.0*eta*psi^3*xi^2 - 104.0*eta*psi^3*xi - 24.0*eta*psi^2*xi^3 - 132.0*eta*psi^2*xi^2 
                                              // + 156.0*eta*psi^2*xi + 16.0*eta*psi*xi^3 + 88.0*eta*psi*xi^2 - 104.0*eta*psi*xi - 4.0*eta*xi^3 - 22.0*eta*xi^2 
                                              // + 26.0*eta*xi
 
        assert_polynomial!(integrand[0][9], &[0.0, 0.0, 26.0, 0.0, 0.0, -48.0, 0.0, -36.0, -104.0, 0.0, 0.0, 18.0, 0.0, 68.0, 192.0, 0.0, 10.0, -272.0, -288.0, 
                                              0.0, 0.0, 4.0, 0.0, -28.0, -72.0, 0.0, -10.0, 112.0, 192.0, 0.0, 0.0, -16.0, -168.0, -272.0, 26.0, 0.0, 0.0, 0.0, 
                                              -4.0, 16.0, 24.0, 0.0, 112.0, 108.0, 0.0, 0.0, 80.0, 408.0, 192.0, 0.0, 0.0, 0.0, 120.0, 112.0, 18.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 24.0, 0.0, 0.0, -40.0, -168.0, -72.0, 0.0, 0.0, 0.0, -60.0, -272.0, -48.0, 0.0, 
                                              0.0, 0.0, -40.0, -28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -24.0, -16.0, 4.0, 0.0, 0.0, 
                                              60.0, 112.0, 18.0, 0.0, 0.0, 0.0, 40.0, 68.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -40.0, -28.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                              // 10.0*eta^3*psi^4*xi^2 - 20.0*eta^3*psi^4*xi + 10.0*eta^3*psi^4 - 40.0*eta^3*psi^3*xi^2 + 80.0*eta^3*psi^3*xi 
                                              // - 40.0*eta^3*psi^3 + 60.0*eta^3*psi^2*xi^2 - 120.0*eta^3*psi^2*xi + 60.0*eta^3*psi^2 - 40.0*eta^3*psi*xi^2 
                                              // + 80.0*eta^3*psi*xi - 40.0*eta^3*psi + 10.0*eta^3*xi^2 - 20.0*eta^3*xi + 10.0*eta^3 - 4.0*eta^2*psi^4*xi^3 
                                              // - 28.0*eta^2*psi^4*xi^2 + 68.0*eta^2*psi^4*xi - 36.0*eta^2*psi^4 + 16.0*eta^2*psi^3*xi^3 + 112.0*eta^2*psi^3*xi^2 
                                              // - 272.0*eta^2*psi^3*xi + 144.0*eta^2*psi^3 - 24.0*eta^2*psi^2*xi^3 - 168.0*eta^2*psi^2*xi^2 + 408.0*eta^2*psi^2*xi 
                                              // - 216.0*eta^2*psi^2 + 16.0*eta^2*psi*xi^3 + 112.0*eta^2*psi*xi^2 - 272.0*eta^2*psi*xi + 144.0*eta^2*psi - 4.0*eta^2*xi^3 
                                              // - 28.0*eta^2*xi^2 + 68.0*eta^2*xi - 36.0*eta^2 + 4.0*eta*psi^4*xi^3 + 18.0*eta*psi^4*xi^2 - 48.0*eta*psi^4*xi 
                                              // + 26.0*eta*psi^4 - 16.0*eta*psi^3*xi^3 - 72.0*eta*psi^3*xi^2 + 192.0*eta*psi^3*xi - 104.0*eta*psi^3 + 24.0*eta*psi^2*xi^3 
                                              // + 108.0*eta*psi^2*xi^2 - 288.0*eta*psi^2*xi + 156.0*eta*psi^2 - 16.0*eta*psi*xi^3 - 72.0*eta*psi*xi^2 + 192.0*eta*psi*xi 
                                              // - 104.0*eta*psi + 4.0*eta*xi^3 + 18.0*eta*xi^2 - 48.0*eta*xi + 26.0*eta
 
        assert_polynomial!(integrand[0][12], &[0.0, 0.0, 0.0, 26.0, 0.0, 0.0, -22.0, 0.0, -36.0, -78.0, 0.0, 0.0, 22.0, 0.0, 32.0, 66.0, 0.0, 10.0, -96.0, 66.0, 0.0, 
                                               0.0, -4.0, 0.0, -32.0, -66.0, -10.0, 96.0, -96.0, -78.0, 0.0, 0.0, 12.0, 0.0, 36.0, 0.0, 0.0, 0.0, 4.0, 0.0, 12.0, 0.0, 
                                               0.0, 96.0, 66.0, 0.0, -10.0, 0.0, 32.0, -26.0, 0.0, 0.0, 30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 
                                               -4.0, 0.0, 30.0, 0.0, -32.0, 22.0, 0.0, 0.0, -30.0, 0.0, 36.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 4.0, 0.0, 0.0, -30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                               // 10.0*eta^2*psi^4*xi - 10.0*eta^2*psi^4 - 30.0*eta^2*psi^3*xi + 30.0*eta^2*psi^3 + 30.0*eta^2*psi^2*xi - 30.0*eta^2*psi^2 
                                               // - 10.0*eta^2*psi*xi + 10.0*eta^2*psi - 4.0*eta*psi^4*xi^2 - 32.0*eta*psi^4*xi + 36.0*eta*psi^4 + 12.0*eta*psi^3*xi^2 
                                               // + 96.0*eta*psi^3*xi - 108.0*eta*psi^3 - 12.0*eta*psi^2*xi^2 - 96.0*eta*psi^2*xi + 108.0*eta*psi^2 + 4.0*eta*psi*xi^2 
                                               // + 32.0*eta*psi*xi - 36.0*eta*psi + 4.0*psi^4*xi^2 + 22.0*psi^4*xi - 26.0*psi^4 - 12.0*psi^3*xi^2 - 66.0*psi^3*xi 
                                               // + 78.0*psi^3 + 12.0*psi^2*xi^2 + 66.0*psi^2*xi - 78.0*psi^2 - 4.0*psi*xi^2 - 22.0*psi*xi + 26.0*psi
 
        assert_polynomial!(integrand[3][3], &[0.0, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, 0.0, 0.0, 156.0, 4.0, 0.0, 0.0, 46.0, 0.0, 0.0, -10.0, -184.0, -372.0, -104.0, 0.0, 
                                              -8.0, -16.0, 46.0, 248.0, 156.0, 0.0, -184.0, 276.0, 0.0, 0.0, 32.0, 276.0, 248.0, 26.0, 0.0, 0.0, 0.0, 4.0, -16.0, 24.0, 
                                              0.0, -184.0, -372.0, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, -60.0, -184.0, -62.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              -16.0, -48.0, -16.0, 0.0, 40.0, 276.0, 248.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 40.0, 46.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 32.0, -8.0, 0.0, 0.0, -60.0, -184.0, -62.0, 0.0, 0.0, 0.0, 40.0, 46.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 
                                              -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 40.0, 46.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); // 40.0*eta^3*psi^3*xi^2 - 10.0*eta^3*psi^4*xi^2 - 60.0*eta^3*psi^2*xi^2 
                                              // + 40.0*eta^3*psi*xi^2 - 10.0*eta^3*xi^2 + 4.0*eta^2*psi^4*xi^3 + 46.0*eta^2*psi^4*xi^2 - 16.0*eta^2*psi^3*xi^3 
                                              // - 184.0*eta^2*psi^3*xi^2 + 24.0*eta^2*psi^2*xi^3 + 276.0*eta^2*psi^2*xi^2 - 16.0*eta^2*psi*xi^3 - 184.0*eta^2*psi*xi^2 
                                              // + 4.0*eta^2*xi^3 + 46.0*eta^2*xi^2 - 8.0*eta*psi^4*xi^3 - 62.0*eta*psi^4*xi^2 + 32.0*eta*psi^3*xi^3 + 248.0*eta*psi^3*xi^2 
                                              // - 48.0*eta*psi^2*xi^3 - 372.0*eta*psi^2*xi^2 + 32.0*eta*psi*xi^3 + 248.0*eta*psi*xi^2 - 8.0*eta*xi^3 - 62.0*eta*xi^2 
                                              // + 4.0*psi^4*xi^3 + 26.0*psi^4*xi^2 - 16.0*psi^3*xi^3 - 104.0*psi^3*xi^2 + 24.0*psi^2*xi^3 + 156.0*psi^2*xi^2 - 16.0*psi*xi^3 
                                              // - 104.0*psi*xi^2 + 4.0*xi^3 + 26.0*xi^2
 
        assert_polynomial!(integrand[3][6], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 156.0, 4.0, 0.0, 0.0, -36.0, 0.0, 0.0, 10.0, 144.0, -216.0, -104.0, 0.0, 4.0, 
                                              -16.0, -36.0, 144.0, 156.0, 0.0, 144.0, -216.0, 0.0, 0.0, -16.0, -216.0, 144.0, 26.0, 0.0, 0.0, 0.0, -4.0, 16.0, -24.0, 0.0, 
                                              144.0, 156.0, 0.0, 0.0, -40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60.0, 144.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 24.0, 
                                              -16.0, 0.0, -40.0, -216.0, 144.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -40.0, -36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -24.0, -16.0, 4.0, 0.0, 0.0, 60.0, 144.0, 26.0, 0.0, 0.0, 0.0, -40.0, -36.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              -40.0, -36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                              // 10.0*eta^3*psi^4*xi^2 - 40.0*eta^3*psi^3*xi^2 + 60.0*eta^3*psi^2*xi^2 - 40.0*eta^3*psi*xi^2 + 10.0*eta^3*xi^2 
                                              // - 4.0*eta^2*psi^4*xi^3 - 36.0*eta^2*psi^4*xi^2 + 16.0*eta^2*psi^3*xi^3 + 144.0*eta^2*psi^3*xi^2 - 24.0*eta^2*psi^2*xi^3 
                                              // - 216.0*eta^2*psi^2*xi^2 + 16.0*eta^2*psi*xi^3 + 144.0*eta^2*psi*xi^2 - 4.0*eta^2*xi^3 - 36.0*eta^2*xi^2 + 4.0*eta*psi^4*xi^3 
                                              // + 26.0*eta*psi^4*xi^2 - 16.0*eta*psi^3*xi^3 - 104.0*eta*psi^3*xi^2 + 24.0*eta*psi^2*xi^3 + 156.0*eta*psi^2*xi^2 
                                              // - 16.0*eta*psi*xi^3 - 104.0*eta*psi*xi^2 + 4.0*eta*xi^3 + 26.0*eta*xi^2
 
        assert_polynomial!(integrand[3][9], &[0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, -104.0, 0.0, 0.0, -22.0, 0.0, -36.0, 88.0, 0.0, 10.0, 144.0, -132.0, 0.0, 0.0, -4.0, 
                                              0.0, 32.0, 88.0, 0.0, -10.0, -128.0, 156.0, 0.0, 0.0, 16.0, 192.0, 144.0, 26.0, 0.0, 0.0, 0.0, 4.0, -16.0, -24.0, 0.0, -128.0, 
                                              -132.0, 0.0, 0.0, -40.0, -216.0, -104.0, 0.0, 0.0, 0.0, -60.0, -128.0, -22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -24.0, 
                                              0.0, 0.0, 40.0, 192.0, 88.0, 0.0, 0.0, 0.0, 60.0, 144.0, 26.0, 0.0, 0.0, 0.0, 40.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 24.0, 16.0, -4.0, 0.0, 0.0, -60.0, -128.0, -22.0, 0.0, 0.0, 0.0, -40.0, -36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 40.0, 
                                              32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                              // 10.0*eta^3*psi^4*xi - 10.0*eta^3*psi^4*xi^2 + 40.0*eta^3*psi^3*xi^2 - 40.0*eta^3*psi^3*xi - 60.0*eta^3*psi^2*xi^2 
                                              // + 60.0*eta^3*psi^2*xi + 40.0*eta^3*psi*xi^2 - 40.0*eta^3*psi*xi - 10.0*eta^3*xi^2 + 10.0*eta^3*xi + 4.0*eta^2*psi^4*xi^3 
                                              // + 32.0*eta^2*psi^4*xi^2 - 36.0*eta^2*psi^4*xi - 16.0*eta^2*psi^3*xi^3 - 128.0*eta^2*psi^3*xi^2 + 144.0*eta^2*psi^3*xi 
                                              // + 24.0*eta^2*psi^2*xi^3 + 192.0*eta^2*psi^2*xi^2 - 216.0*eta^2*psi^2*xi - 16.0*eta^2*psi*xi^3 - 128.0*eta^2*psi*xi^2 
                                              // + 144.0*eta^2*psi*xi + 4.0*eta^2*xi^3 + 32.0*eta^2*xi^2 - 36.0*eta^2*xi - 4.0*eta*psi^4*xi^3 - 22.0*eta*psi^4*xi^2 
                                              // + 26.0*eta*psi^4*xi + 16.0*eta*psi^3*xi^3 + 88.0*eta*psi^3*xi^2 - 104.0*eta*psi^3*xi - 24.0*eta*psi^2*xi^3 
                                              // - 132.0*eta*psi^2*xi^2 + 156.0*eta*psi^2*xi + 16.0*eta*psi*xi^3 + 88.0*eta*psi*xi^2 - 104.0*eta*psi*xi - 4.0*eta*xi^3 
                                              // - 22.0*eta*xi^2 + 26.0*eta*xi
 
        assert_polynomial!(integrand[3][12], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 0.0, -36.0, -78.0, 0.0, 0.0, -26.0, 0.0, 36.0, 78.0, 0.0, 0.0, -108.0, 108.0, 0.0, 0.0, 
                                               4.0, 0.0, 36.0, -108.0, 0.0, 108.0, -108.0, -78.0, 0.0, 0.0, -12.0, 0.0, 36.0, 0.0, 0.0, 0.0, -4.0, 0.0, 12.0, 0.0, 0.0, 
                                               108.0, 78.0, 0.0, 10.0, 0.0, 36.0, 0.0, 0.0, 0.0, -30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -4.0, 0.0, 
                                               30.0, 0.0, 36.0, -26.0, 0.0, 0.0, -30.0, 0.0, 36.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, -12.0, 0.0, 4.0, 0.0, 0.0, 30.0, 0.0, 36.0, 0.0, 0.0, 0.0, -30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                               // 30.0*eta^2*psi^3*xi - 10.0*eta^2*psi^4*xi - 30.0*eta^2*psi^2*xi + 10.0*eta^2*psi*xi + 4.0*eta*psi^4*xi^2 + 36.0*eta*psi^4*xi 
                                               // - 12.0*eta*psi^3*xi^2 - 108.0*eta*psi^3*xi + 12.0*eta*psi^2*xi^2 + 108.0*eta*psi^2*xi - 4.0*eta*psi*xi^2 - 36.0*eta*psi*xi 
                                               // - 4.0*psi^4*xi^2 - 26.0*psi^4*xi + 12.0*psi^3*xi^2 + 78.0*psi^3*xi - 12.0*psi^2*xi^2 - 78.0*psi^2*xi + 4.0*psi*xi^2 
                                               // + 26.0*psi*xi
 
        assert_polynomial!(integrand[6][6], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 156.0, 4.0, 0.0, 0.0, 26.0, 0.0, 0.0, -10.0, -104.0, -372.0, -104.0, 0.0, 4.0, 
                                              -16.0, 26.0, 248.0, 156.0, 0.0, -104.0, 156.0, 0.0, 0.0, -16.0, 156.0, 248.0, 26.0, 0.0, 0.0, 0.0, 4.0, -16.0, 24.0, 0.0, 
                                              -104.0, -372.0, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, -60.0, -104.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -48.0, 
                                              -16.0, 0.0, 40.0, 156.0, 248.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 40.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 32.0, 4.0, 0.0, 0.0, -60.0, -104.0, 26.0, 0.0, 0.0, 0.0, 40.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 40.0, 26.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                              // 40.0*eta^3*psi^3*xi^2 - 10.0*eta^3*psi^4*xi^2 - 60.0*eta^3*psi^2*xi^2 + 40.0*eta^3*psi*xi^2 - 10.0*eta^3*xi^2 
                                              // + 4.0*eta^2*psi^4*xi^3 + 26.0*eta^2*psi^4*xi^2 - 16.0*eta^2*psi^3*xi^3 - 104.0*eta^2*psi^3*xi^2 + 24.0*eta^2*psi^2*xi^3 
                                              // + 156.0*eta^2*psi^2*xi^2 - 16.0*eta^2*psi*xi^3 - 104.0*eta^2*psi*xi^2 + 4.0*eta^2*xi^3 + 26.0*eta^2*xi^2
 
        assert_polynomial!(integrand[6][9], &[0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, -104.0, 0.0, -4.0, -22.0, 0.0, 26.0, 88.0, 0.0, -10.0, 88.0, -132.0, 0.0, 0.0, -4.0, 
                                              16.0, -22.0, 88.0, 0.0, 10.0, -132.0, 156.0, 0.0, 0.0, 16.0, 88.0, 88.0, 26.0, 0.0, 0.0, 0.0, 4.0, -24.0, -24.0, 0.0, 88.0, 
                                              -132.0, 0.0, 0.0, -40.0, 156.0, -104.0, 0.0, 0.0, 0.0, 60.0, 88.0, -22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -24.0, 16.0, 
                                              0.0, 40.0, 88.0, 88.0, 26.0, 0.0, 0.0, -60.0, 88.0, 26.0, 0.0, 0.0, 0.0, -40.0, -22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, -24.0, 16.0, -4.0, 0.0, 0.0, 60.0, 88.0, -22.0, 0.0, 0.0, 0.0, -40.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -24.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 40.0, 
                                              -22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                              // 10.0*eta^3*psi^4*xi^2 - 10.0*eta^3*psi^4*xi - 40.0*eta^3*psi^3*xi^2 + 40.0*eta^3*psi^3*xi + 60.0*eta^3*psi^2*xi^2 
                                              // - 60.0*eta^3*psi^2*xi - 40.0*eta^3*psi*xi^2 + 40.0*eta^3*psi*xi + 10.0*eta^3*xi^2 - 10.0*eta^3*xi - 4.0*eta^2*psi^4*xi^3 
                                              // - 22.0*eta^2*psi^4*xi^2 + 26.0*eta^2*psi^4*xi + 16.0*eta^2*psi^3*xi^3 + 88.0*eta^2*psi^3*xi^2 - 104.0*eta^2*psi^3*xi 
                                              // - 24.0*eta^2*psi^2*xi^3 - 132.0*eta^2*psi^2*xi^2 + 156.0*eta^2*psi^2*xi + 16.0*eta^2*psi*xi^3 + 88.0*eta^2*psi*xi^2 
                                              // - 104.0*eta^2*psi*xi - 4.0*eta^2*xi^3 - 22.0*eta^2*xi^2 + 26.0*eta^2*xi
 
        assert_polynomial!(integrand[6][12], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 0.0, -36.0, -78.0, 0.0, 0.0, -26.0, 0.0, 36.0, 78.0, 0.0, 10.0, -108.0, 108.0, 0.0, 0.0, 
                                               4.0, 0.0, 36.0, -108.0, 0.0, 108.0, -108.0, -78.0, 0.0, 0.0, -12.0, 0.0, 36.0, 0.0, 0.0, 0.0, -4.0, 0.0, 12.0, 0.0, 0.0, 108.0, 
                                               78.0, 0.0, -10.0, 0.0, 36.0, -26.0, 0.0, 0.0, 30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -4.0, 0.0, -30.0, 
                                               0.0, 36.0, -26.0, 0.0, 0.0, 30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, -12.0, 0.0, 4.0, 0.0, 0.0, -30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                               // 10.0*eta^2*psi^4*xi - 30.0*eta^2*psi^3*xi + 30.0*eta^2*psi^2*xi - 10.0*eta^2*psi*xi - 4.0*eta*psi^4*xi^2 - 26.0*eta*psi^4*xi 
                                               // + 12.0*eta*psi^3*xi^2 + 78.0*eta*psi^3*xi - 12.0*eta*psi^2*xi^2 - 78.0*eta*psi^2*xi + 4.0*eta*psi*xi^2 + 26.0*eta*psi*xi
 
        assert_polynomial!(integrand[9][9], &[26.0, 0.0, 26.0, -104.0, 18.0, 0.0, 0.0, 26.0, -104.0, 156.0, 4.0, 18.0, 0.0, -48.0, 0.0, 0.0, -10.0, -72.0, 108.0, 0.0, 0.0, 
                                              4.0, -16.0, 18.0, 0.0, 108.0, -10.0, -72.0, 0.0, -104.0, 0.0, 0.0, 24.0, -16.0, 26.0, 0.0, 0.0, 0.0, 4.0, -16.0, 24.0, 0.0, 
                                              -72.0, 108.0, 0.0, 0.0, 20.0, 192.0, 0.0, 26.0, 0.0, 0.0, -60.0, -72.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -72.0, 
                                              -16.0, 0.0, -80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, -72.0, -48.0, 0.0, 0.0, 0.0, 40.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 192.0, 4.0, 0.0, 0.0, -60.0, -72.0, 18.0, 0.0, 0.0, 0.0, 40.0, -48.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -80.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              40.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                              // 20.0*eta^3*psi^4*xi - 10.0*eta^3*psi^4*xi^2 - 10.0*eta^3*psi^4 + 40.0*eta^3*psi^3*xi^2 - 80.0*eta^3*psi^3*xi + 40.0*eta^3*psi^3 
                                              // - 60.0*eta^3*psi^2*xi^2 + 120.0*eta^3*psi^2*xi - 60.0*eta^3*psi^2 + 40.0*eta^3*psi*xi^2 - 80.0*eta^3*psi*xi + 40.0*eta^3*psi 
                                              // - 10.0*eta^3*xi^2 + 20.0*eta^3*xi - 10.0*eta^3 + 4.0*eta^2*psi^4*xi^3 + 18.0*eta^2*psi^4*xi^2 - 48.0*eta^2*psi^4*xi 
                                              // + 26.0*eta^2*psi^4 - 16.0*eta^2*psi^3*xi^3 - 72.0*eta^2*psi^3*xi^2 + 192.0*eta^2*psi^3*xi - 104.0*eta^2*psi^3 
                                              // + 24.0*eta^2*psi^2*xi^3 + 108.0*eta^2*psi^2*xi^2 - 288.0*eta^2*psi^2*xi + 156.0*eta^2*psi^2 - 16.0*eta^2*psi*xi^3 
                                              // - 72.0*eta^2*psi*xi^2 + 192.0*eta^2*psi*xi - 104.0*eta^2*psi + 4.0*eta^2*xi^3 + 18.0*eta^2*xi^2 - 48.0*eta^2*xi + 26.0*eta^2
 
        assert_polynomial!(integrand[9][12], &[0.0, 0.0, 0.0, 26.0, 0.0, 0.0, -22.0, 0.0, -36.0, -78.0, 0.0, 0.0, 22.0, 0.0, 32.0, 66.0, 0.0, 10.0, -96.0, 66.0, 0.0, 0.0, -4.0, 
                                               0.0, -32.0, -66.0, -10.0, 96.0, -96.0, -78.0, 0.0, 0.0, 12.0, 0.0, 36.0, 0.0, 0.0, 0.0, 4.0, 0.0, 12.0, 0.0, 0.0, 96.0, 66.0, 0.0, 
                                               10.0, 0.0, 32.0, -26.0, 0.0, 0.0, -30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, -4.0, 0.0, -30.0, 0.0, -32.0, 
                                               22.0, 0.0, 0.0, 30.0, 0.0, 36.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 
                                               4.0, 0.0, 0.0, 30.0, 0.0, 36.0, 0.0, 0.0, 0.0, -30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); // 10.0*eta^2*psi^4 - 10.0*eta^2*psi^4*xi + 30.0*eta^2*psi^3*xi - 30.0*eta^2*psi^3 
                                               // - 30.0*eta^2*psi^2*xi + 30.0*eta^2*psi^2 + 10.0*eta^2*psi*xi - 10.0*eta^2*psi + 4.0*eta*psi^4*xi^2 + 22.0*eta*psi^4*xi 
                                               // - 26.0*eta*psi^4 - 12.0*eta*psi^3*xi^2 - 66.0*eta*psi^3*xi + 78.0*eta*psi^3 + 12.0*eta*psi^2*xi^2 + 66.0*eta*psi^2*xi 
                                               // - 78.0*eta*psi^2 - 4.0*eta*psi*xi^2 - 22.0*eta*psi*xi + 26.0*eta*psi
 
        assert_polynomial!(integrand[12][12], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                0.0, 0.0, 0.0, 20.0, -52.0, 0.0, 4.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 
                                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]); 
                                                // 4.0*xi*psi^4 - 8.0*xi*psi^3 + 4.0*xi*psi^2 - 10.0*eta*psi^4 + 20.0*eta*psi^3 - 10.0*eta*psi^2 + 26.0*psi^4 - 52.0*psi^3 + 26.0*psi^2
 

        // Test zero components 
        let zero_indices = [(0,1), (0,2), (0,4), (0,5), (0,7), (0,8), (0,10), (0,11), (0,13), (0,14), (1,2), (1,3), (1,5), (1,6), (1,8), (1,9), (1,11), (1,12), 
                                                  (1,14), (2,3), (2,4), (2,6), (2,7), (2,9), (2,10), (2,12), (2,13), (3,4), (3,5), (3,7), (3,8), (3,10), (3,11), (3,13), (3,14), 
                                                  (4,5), (4,6), (4,8), (4,9), (4,11), (4,12), (4,14), (5,6), (5,7), (5,9), (5,10), (5,12), (5,13), (6,7), (6,8), (6,10), (6,11), 
                                                  (6,13), (6,14), (7,8), (7,9), (7,11), (7,12), (7,14), (8,9), (8,10), (8,12), (8,13), (9,10), (9,11), (9,13), (9,14), (10,11), 
                                                  (10,12), (10,14), (11,12), (11,13), (12,13), (12,14), (13,14)];
        for (i, j) in zero_indices {
            assert_polynomial!(integrand[i][j], &[0.0]); // 0
        }
        
        // Test symmetry 
        assert_polynomial!(integrand[1][1], integrand[0][0]); // Same as (0,0)
        assert_polynomial!(integrand[2][2], integrand[0][0]); // Same as (0,0)
        assert_polynomial!(integrand[4][4], integrand[3][3]); // Same as (3,3)
        assert_polynomial!(integrand[5][5], integrand[3][3]); // Same as (3,3)
        assert_polynomial!(integrand[7][7], integrand[6][6]); // Same as (6,6)
        assert_polynomial!(integrand[8][8], integrand[6][6]); // Same as (6,6)
        assert_polynomial!(integrand[10][10], integrand[9][9]); // Same as (9,9)
        assert_polynomial!(integrand[11][11], integrand[9][9]); // Same as (9,9)
        assert_polynomial!(integrand[13][13], integrand[12][12]); // Same as (12,12)
        assert_polynomial!(integrand[14][14], integrand[12][12]); // Same as (12,12)
        
        let symmetric_pairs = [
            (1,0), (2,0), (2,1), (3,0), (3,1), (3,2), (4,0), (4,1), (4,2), (4,3), (5,0), 
            (5,1), (5,2), (5,3), (5,4), (6,0), (6,1), (6,2), (6,3), (6,4), (6,5), (7,0), (7,1), (7,2), (7,3), 
            (7,4), (7,5), (7,6), (8,0), (8,1), (8,2), (8,3), (8,4), (8,5), (8,6), (8,7), 
            (9,0), (9,1), (9,2), (9,3), (9,4), (9,5), (9,6), (9,7), (9,8), (10,0), (10,1), 
            (10,2), (10,3), (10,4), (10,5), (10,6), (10,7), (10,8), (10,9), (11,0), (11,1), 
            (11,2), (11,3), (11,4), (11,5), (11,6), (11,7), (11,8), (11,9), (11,10), (12,0), 
            (12,1), (12,2), (12,3), (12,4), (12,5), (12,6), (12,7), (12,8), (12,9), (12,10), 
            (12,11), (13,0), (13,1), (13,2), (13,3), (13,4), (13,5), (13,6), (13,7), (13,8), 
            (13,9), (13,10), (13,11), (13,12), (14,0), (14,1), (14,2), (14,3), (14,4), 
            (14,5),(14,6),(14,7),(14,8),(14,9),(14,10),(14,11),(14,12),(14,13)
        ];
        for (i, j) in symmetric_pairs {
            assert_polynomial!(integrand[i][j], integrand[j][i]); // Symmetry check
        }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    }

    // TEST 9: 3D Pyramid in 3D mesh - VARIABLE material - MASS matrix
    #[test]
    fn test_3d_pyramid_variable_mass() {
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3, 4] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![-1.0, -1.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, -1.0, 0.0] },
            Node { id: 2, coordinates: vec![1.0, 1.0, 0.0] },
            Node { id: 3, coordinates: vec![-1.0, 1.0, 0.0] },
            Node { id: 4, coordinates: vec![0.0, 0.0, 2.0] },
        ];
        let element_type = ElementType::Pyramid;
        let mesh_dim = 3;

        // VARIABLE anisotropic material for 3D (6x6 elasticity matrix)
        let c_matrix = vec![
            vec![
                vec![1.0, 0.5, 0.0, 0.0, 0.0, 0.0], // C11 = 1 + 0.5x
                vec![0.3, 0.0, 0.0, 0.0, 0.0, 0.0], // C12 = 0.3
                vec![0.2, 0.0, 0.0, 0.0, 0.0, 0.0], // C13 = 0.2
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C14 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C15 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C16 = 0
            ],
            vec![
                vec![0.3, 0.0, 0.0, 0.0, 0.0, 0.0], // C21 = 0.3
                vec![2.0, 0.0, 0.3, 0.0, 0.0, 0.0], // C22 = 2 + 0.3y
                vec![0.4, 0.0, 0.0, 0.0, 0.0, 0.0], // C23 = 0.4
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C24 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C25 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C26 = 0
            ],
            vec![
                vec![0.2, 0.0, 0.0, 0.0, 0.0, 0.0], // C31 = 0.2
                vec![0.4, 0.0, 0.0, 0.0, 0.0, 0.0], // C32 = 0.4
                vec![3.0, 0.0, 0.0, 0.4, 0.0, 0.0], // C33 = 3 + 0.4z
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C34 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C35 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C36 = 0
            ],
            vec![
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C41 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C42 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C43 = 0
                vec![0.5, 0.2, 0.0, 0.0, 0.0, 0.0], // C44 = 0.5 + 0.2x
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C45 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C46 = 0
            ],
            vec![
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C51 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C52 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C53 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C54 = 0
                vec![0.6, 0.0, 0.3, 0.0, 0.0, 0.0], // C55 = 0.6 + 0.3y
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C56 = 0
            ],
            vec![
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C61 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C62 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C63 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C64 = 0
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // C65 = 0
                vec![0.7, 0.0, 0.0, 0.5, 0.0, 0.0], // C66 = 0.7 + 0.5z
            ],
        ];

        let material_property_physical = MaterialProperty::Matrix(c_matrix);

        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Stiffness,
            &element_type,
            &nodes,
            &material_property_physical,
            mesh_dim,
        ).unwrap();

        let matrix_size = mesh_dim * 5; // 5 nodes × 3 DOF each = 15
        
        // Check symmetry of stiffness matrix
        for i in 0..matrix_size {
            for j in 0..matrix_size {
                let diff: f64 = integrand[i][j].iter().zip(&integrand[j][i]).map(|(a, b)| (a - b).abs()).sum();
                assert!(diff < 1e-10, "Stiffness matrix not symmetric at ({}, {})", i, j);
            }
        }
    }

    // TEST 9: Optimal Gauss points - CONSTANT material
    #[test]
    fn test_optimal_gauss_points_constant_material() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        
        // CONSTANT material
        let material_property = MaterialProperty::Scalar(vec![1.0]);
        
        let result = GaussianQuadrature::find_optimal_gauss_points(
            IntegrationType::Mass,
            &element,
            &element_type,
            &nodes,
            1,
            1e-10,
            &material_property,
        ).unwrap();

        // For linear element with constant density: integrand is quadratic
        // Need n ≥ 2 for exact integration (2n-1 ≥ 2 ⇒ n ≥ 1.5 ⇒ n=2)
        assert_eq!(result.optimal_number, 2);
    }

    // TEST 10: Optimal Gauss points - VARIABLE material
    #[test]
    fn test_optimal_gauss_points_variable_material() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        
        // VARIABLE material: ρ(x) = 1 + x (linear)
        let material_property = MaterialProperty::Scalar(vec![1.0, 1.0, 0.0, 0.0]);
        
        let result = GaussianQuadrature::find_optimal_gauss_points(
            IntegrationType::Mass,
            &element,
            &element_type,
            &nodes,
            1,
            1e-10,
            &material_property,
        ).unwrap();

        // For linear element with linear density: integrand is cubic
        // Need n ≥ 2 for exact integration (2n-1 ≥ 3 ⇒ n ≥ 2 ⇒ n=2)
        // Actually cubic requires n=2 since 2×2-1=3
        assert_eq!(result.optimal_number, 2);
    }

    // TEST 11: Material property transformation accuracy
    #[test]
    fn test_material_transformation_accuracy() {
        let element = Element { id: 0, nodes: vec![0, 1] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0] },
            Node { id: 1, coordinates: vec![2.0] },
        ];
        let element_type = ElementType::Line;
        let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
        
        // Test transformation: ρ(x) = 1 + x in physical space
        // Should become ρ(xi) = 1 + 2xi in parametric space
        let material_physical = vec![1.0, 1.0, 0.0, 0.0]; // 1 + x
        
        let transformed = GaussianQuadrature::transform_material_property_scalar(
            &material_physical,
            &nodes,
            &shape_function.values,
            1,
        ).unwrap();

        // Expected: ρ(xi) = 1 + 2xi
        let expected = vec![1.0, 2.0, 0.0, 0.0];
        
        assert_eq!(transformed.len(), expected.len());
        for (i, (&actual, &expected)) in transformed.iter().zip(expected.iter()).enumerate() {
            assert!((actual - expected).abs() < 1e-10, 
                   "Coefficient {}: expected {}, got {}", i, expected, actual);
        }
    }

    // TEST 12: Complex 3D case with variable material
    #[test]
    fn test_3d_tetra_variable_material() {
        let element = Element { id: 0, nodes: vec![0, 1, 2, 3] };
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0, 0.0] },
            Node { id: 2, coordinates: vec![0.0, 1.0, 0.0] },
            Node { id: 3, coordinates: vec![0.0, 0.0, 1.0] },
        ];
        let element_type = ElementType::Tetra;
        
        // VARIABLE density in 3D: ρ(x,y,z) = 1 + x + y + z
        let material_property = MaterialProperty::Scalar(vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        
        let integrand = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &nodes,
            &material_property,
            3,
        ).unwrap();

        // Just test it runs and produces results
        assert_eq!(integrand.len(), 12); // 4 nodes × 3 DOF = 12×12
        
        // Check that higher order terms are present due to variable material
        let max_degree = GaussianQuadrature::find_max_degree_in_matrix(&integrand);
        assert!(max_degree >= 2, "Variable material should increase polynomial degree");
    }
}

/* 
#[cfg(test)]
mod tests {

    use super::*;
    use crate::structs_and_impls::*;

    // Test helper functions
    mod helper_tests {
        use super::*;

        #[test]
        fn test_factorial_small_values() {
            // Test factorial for small values that are precomputed
            assert_eq!(factorial(0), 1.0);
            assert_eq!(factorial(1), 1.0);
            assert_eq!(factorial(2), 2.0);
            assert_eq!(factorial(3), 6.0);
            assert_eq!(factorial(4), 24.0);
            assert_eq!(factorial(5), 120.0);
        }

        #[test]
        fn test_factorial_large_values() {
            // Test factorial for larger values using Stirling's approximation
            // We'll check that it's reasonably close to expected values
            let fact_10 = factorial(10);
            let expected_10 = 3628800.0;
            assert!((fact_10 - expected_10).abs() / expected_10 < 0.01, 
                   "Factorial(10) should be close to 3628800, got {}", fact_10);

            let fact_20 = factorial(20);
            let expected_20 = 2432902008176640000.0;
            assert!((fact_20 - expected_20).abs() / expected_20 < 0.001,
                   "Factorial(20) should be reasonably accurate");
        }

        #[test]
        fn test_detect_polynomial_order_1d() {
            // Test constant polynomial
            let coeffs_constant = vec![1.0, 0.0, 0.0, 0.0];
            assert_eq!(GaussianQuadrature::detect_polynomial_order_1d(&coeffs_constant), 0);

            // Test linear polynomial
            let coeffs_linear = vec![1.0, 2.0, 0.0, 0.0];
            assert_eq!(GaussianQuadrature::detect_polynomial_order_1d(&coeffs_linear), 1);

            // Test quadratic polynomial
            let coeffs_quadratic = vec![1.0, 2.0, 3.0, 0.0];
            assert_eq!(GaussianQuadrature::detect_polynomial_order_1d(&coeffs_quadratic), 2);

            // Test cubic polynomial
            let coeffs_cubic = vec![1.0, 2.0, 3.0, 4.0];
            assert_eq!(GaussianQuadrature::detect_polynomial_order_1d(&coeffs_cubic), 3);

            // Test with tolerance - small coefficients should be ignored
            let coeffs_with_noise = vec![1.0, 1e-13, 1e-13, 2.0];
            assert_eq!(GaussianQuadrature::detect_polynomial_order_1d(&coeffs_with_noise), 0);
        }

        #[test]
        fn test_detect_polynomial_orders_2d() {
            // Test constant polynomial in 2D
            let coeffs_constant = vec![
                vec![1.0, 0.0],
                vec![0.0, 0.0],
            ];
            assert_eq!(GaussianQuadrature::detect_polynomial_orders_2d(&coeffs_constant), (0, 0));

            // Test mixed orders
            let coeffs_mixed = vec![
                vec![1.0, 2.0, 0.0],
                vec![3.0, 4.0, 5.0],
                vec![0.0, 6.0, 0.0],
            ];
            assert_eq!(GaussianQuadrature::detect_polynomial_orders_2d(&coeffs_mixed), (2, 2));
        }

        #[test]
        fn test_detect_polynomial_orders_3d() {
            // Simple 3D test case
            let coeffs = vec![
                vec![
                    vec![1.0, 0.0],
                    vec![0.0, 2.0],
                ],
                vec![
                    vec![0.0, 0.0],
                    vec![3.0, 0.0],
                ],
            ];
            assert_eq!(GaussianQuadrature::detect_polynomial_orders_3d(&coeffs), (1, 1, 1));
        }
    }

    // Test monomial polynomial operations
    mod monomial_tests {
        use super::*;

        #[test]
        fn test_total_degree_polynomial() {
            // Test constant polynomial
            assert_eq!(MonomialPolynomial::total_degree_polynomial(&vec![1.0]), 0);
            
            // Test linear polynomial
            assert_eq!(MonomialPolynomial::total_degree_polynomial(&vec![1.0, 2.0]), 1);
            
            // Test quadratic polynomial
            assert_eq!(MonomialPolynomial::total_degree_polynomial(&vec![1.0, 2.0, 3.0]), 2);
        }
    }

    // Main test for 1D line element mass integration
    mod line_element_1d_tests {
        use super::*;

        pub fn create_1d_line_element() -> (Element, Vec<Node>, ElementType, MaterialProperty, MeshData) {
            // Create a simple 1D line element from x=0 to x=1 with 0-based indexing
            let element = Element {
                id: 0,
                nodes: vec![0, 1],  // 0-based node indices
            };

            let nodes = vec![
                Node {
                    id: 0,
                    coordinates: vec![-1.0, 0.0, 0.0],
                },
                Node {
                    id: 1,
                    coordinates: vec![1.0, 0.0, 0.0],
                },
            ];

            // Constant density material property
            let material_property = MaterialProperty::Scalar(vec![1.0]); // ρ = 1.0

            // Create minimal mesh data for testing
            let mesh_data = MeshData {
                dimension: 3,
                num_nodes: 2,
                min_node_index: 0,
                nodes: nodes.clone(),
                num_eltypes: 1,
                elements: vec![element.clone()],
                element_type_info: vec![ElementTypeInfo {
                    element_type: ElementType::Line,
                    num_elements: 1,
                    start_index: 0,
                    nodes_per_element: 2,
                }],
            };

            (element, nodes, ElementType::Line, material_property, mesh_data)
        }

        #[test]
        fn test_1d_line_shape_functions() {
            let (_, _, element_type, _, _) = create_1d_line_element();
            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            // Line2 element should have 2 nodes
            assert_eq!(shape_function.num_nodes, 2);
            
            // Shape function values should be for parametric coordinate xi ∈ [-1, 1]
            // N1 = (1 - xi)/2, N2 = (1 + xi)/2
            let n1 = &shape_function.values[0]; // N1 coefficients
            let n2 = &shape_function.values[1]; // N2 coefficients
            
            // N1 = 0.5 - 0.5xi → coefficients: [0.5, -0.5] for basis [1, xi]
            assert_eq!(n1.len(), 2);
            assert!((n1[0] - 0.5).abs() < 1e-12, "N1 constant term should be 0.5");
            assert!((n1[1] - (-0.5)).abs() < 1e-12, "N1 linear term should be -0.5");
            
            // N2 = 0.5 + 0.5xi → coefficients: [0.5, 0.5] for basis [1, xi]
            assert!((n2[0] - 0.5).abs() < 1e-12, "N2 constant term should be 0.5");
            assert!((n2[1] - 0.5).abs() < 1e-12, "N2 linear term should be 0.5");
        }

        #[test]
        fn test_1d_line_shape_function_derivatives() {
            let (_, _, element_type, _, _) = create_1d_line_element();
            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            // Derivatives: dN1/dxi = -0.5, dN2/dxi = 0.5
            let dn1_dxi = &shape_function.derivatives[0][0]; // dN1/dxi
            let dn2_dxi = &shape_function.derivatives[1][0]; // dN2/dxi
            
            // Derivatives should be constant polynomials
            assert_eq!(dn1_dxi.len(), 1);
            assert!((dn1_dxi[0] - (-0.5)).abs() < 1e-12, "dN1/dxi should be -0.5");
            
            assert_eq!(dn2_dxi.len(), 1);
            assert!((dn2_dxi[0] - 0.5).abs() < 1e-12, "dN2/dxi should be 0.5");
        }

        #[test]
        fn test_1d_line_jacobian() {
            let (_, nodes, element_type, _, _) = create_1d_line_element();
            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            let jacobian = GaussianQuadrature::calculate_jacobian_monomial(
                &nodes,
                &shape_function.derivatives,
                1, // element_dim
                2, // num_nodes  
                3, // mesh_dim (3D space)
                1, // element_order
            ).unwrap();

            // For a straight line from (0,0,0) to (1,0,0):
            // x(xi) = (1-xi)/2 * 0 + (1+xi)/2 * 1 = (1+xi)/2
            // dx/dxi = 0.5
            // Similarly, y and z coordinates are zero, so their derivatives are zero
            
            // Jacobian matrix should be [dx/dxi, dy/dxi, dz/dxi] = [0.5, 0, 0]
            let j_matrix = &jacobian.matrix;
            
            // Check dx/dxi = 0.5
            assert_eq!(j_matrix[0][0].len(), 1); // Constant polynomial
            assert!((j_matrix[0][0][0] - 0.5).abs() < 1e-12, "dx/dxi should be 0.5");
            
            // Check dy/dxi = 0
            assert!((j_matrix[1][0][0]).abs() < 1e-12, "dy/dxi should be 0");
            
            // Check dz/dxi = 0
            assert!((j_matrix[2][0][0]).abs() < 1e-12, "dz/dxi should be 0");
        }

        #[test]
        fn test_1d_line_jacobian_determinant() {
            let (_, nodes, element_type, _, _) = create_1d_line_element();
            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            let jacobian = GaussianQuadrature::calculate_jacobian_monomial(
                &nodes,
                &shape_function.derivatives,
                1, // element_dim
                2, // num_nodes  
                3, // mesh_dim (3D space)
                1, // element_order
            ).unwrap();

            // For 1D element in 3D space, determinant is metric: dx^2 + dy^2 + dz^2
            // dx/dxi = 0.5, dy/dxi = 0, dz/dxi = 0 → det = (0.5)^2 = 0.25
            let det_j = &jacobian.determinant;
            
            assert_eq!(det_j.len(), 1); // Should be constant
            assert!((det_j[0] - 0.25).abs() < 1e-12, "detJ should be 0.25");
        }

        #[test]
        fn test_1d_line_inverse_jacobian() {
            let (_, nodes, element_type, _, _) = create_1d_line_element();
            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            let jacobian = GaussianQuadrature::calculate_jacobian_monomial(
                &nodes,
                &shape_function.derivatives,
                1, // element_dim
                2, // num_nodes  
                3, // mesh_dim (3D space)
                1, // element_order
            ).unwrap();

            // For 1D in 3D, we compute pseudo-inverse
            let inv_j = GaussianQuadrature::inverse_matrix_monomial(
                &jacobian.matrix,
                3, // mesh_dim
                1, // element_dim
            ).unwrap();

            // For J = [0.5, 0, 0]^T, J⁺ = [2, 0, 0] (pseudo-inverse)
            // Dimensions: [element_dim × mesh_dim] = [1 × 3]
            assert_eq!(inv_j.len(), 1); // 1 row (element_dim)
            assert_eq!(inv_j[0].len(), 3); // 3 columns (mesh_dim)
            
            // Check inv_j[0][0] = 2.0
            assert!((inv_j[0][0][0] - 2.0).abs() < 1e-12, "invJ[0][0] should be 2.0");
            assert!((inv_j[0][1][0]).abs() < 1e-12, "invJ[0][1] should be 0");
            assert!((inv_j[0][2][0]).abs() < 1e-12, "invJ[0][2] should be 0");
        }

        #[test]
        fn test_1d_line_mass_integrand() {
            let (element, nodes, element_type, material_property, _) = create_1d_line_element();
            
            let integrand = GaussianQuadrature::calculate_integrand(
                &IntegrationType::Mass,
                &element_type,
                &nodes,
                &material_property,
                3, // mesh_dim
            ).unwrap();

            // Mass matrix integrand: ρ * N_i * N_j * detJ
            // ρ = 1.0, detJ = 0.25
            // N1 = 0.5 - 0.5xi, N2 = 0.5 + 0.5xi
            
            // For 2 nodes in 3D, mass matrix is 6x6 (2 nodes × 3 DOF each)
            assert_eq!(integrand.len(), 6);
            assert_eq!(integrand[0].len(), 6);

            // Check M11 component: ρ * N1 * N1 * detJ = 1.0 * (0.5 - 0.5xi)^2 * 0.25
            // = 0.25 * (0.25 - 0.5xi + 0.25xi^2) = 0.0625 - 0.125xi + 0.0625xi^2
            let m11 = &integrand[0][0]; // First DOF of first node
            let expected_m11 = vec![0.0625, -0.125, 0.0625]; // [constant, xi, xi^2] coefficients
            
            assert_eq!(m11.len(), expected_m11.len());
            for (actual, expected) in m11.iter().zip(expected_m11.iter()) {
                assert!((actual - expected).abs() < 1e-12, 
                       "M11 coefficient: expected {}, got {}", expected, actual);
            }

            // Check M12 component: ρ * N1 * N2 * detJ = 1.0 * (0.5 - 0.5xi)(0.5 + 0.5xi) * 0.25
            // = 0.25 * (0.25 - 0.25xi^2) = 0.0625 - 0.0625xi^2
            let m12 = &integrand[0][3]; // First DOF of first node × first DOF of second node
            let expected_m12 = vec![0.0625, 0.0, -0.0625];
            
            assert_eq!(m12.len(), expected_m12.len());
            for (actual, expected) in m12.iter().zip(expected_m12.iter()) {
                assert!((actual - expected).abs() < 1e-12,
                       "M12 coefficient: expected {}, got {}", expected, actual);
            }
        }

        #[test]
        fn test_1d_line_gauss_error_estimation() {
            let (element, nodes, element_type, material_property, _) = create_1d_line_element();
            
            // Calculate integrand first
            let integrand = GaussianQuadrature::calculate_integrand(
                &IntegrationType::Mass,
                &element_type,
                &nodes,
                &material_property,
                3, // mesh_dim
            ).unwrap();

            // Test error calculation for different numbers of Gauss points
            let error_1pt = GaussianQuadrature::calculate_error(&integrand, 1, 1).unwrap();
            let error_2pt = GaussianQuadrature::calculate_error(&integrand, 2, 1).unwrap();

            // For mass matrix with quadratic polynomials, 2-point Gauss should be exact
            // So error for 2 points should be very small
            let max_error_2pt = GaussianQuadrature::find_max_in_matrix(&error_2pt);
            assert!(max_error_2pt < 1e-12, "2-point Gauss should be exact for quadratic polynomials, error: {}", max_error_2pt);

            // 1-point Gauss should have some error for quadratic polynomials
            let max_error_1pt = GaussianQuadrature::find_max_in_matrix(&error_1pt);
            assert!(max_error_1pt > 1e-12, "1-point Gauss should have error for quadratic polynomials");
        }

        #[test]
        fn test_1d_line_optimal_gauss_points() {
            let (element, nodes, element_type, material_property, _) = create_1d_line_element();
            
            let tolerance = 1e-10;
            
            let result = GaussianQuadrature::find_optimal_gauss_points(
                IntegrationType::Mass,
                &element,
                &element_type,
                &nodes,
                3, // mesh_dim
                tolerance,
                &material_property,
            ).unwrap();

            // For quadratic polynomials (max degree 2), theoretical points = ceil((2+1)/2) = 2
            assert_eq!(result.theoretical_number, 2);
            
            // Optimal should also be 2 since it's exact
            assert_eq!(result.optimal_number, 2);
        }

        #[test]
        fn test_1d_line_polynomial_order_detection() {
            let (_, nodes, element_type, material_property, _) = create_1d_line_element();
            
            let integrand = GaussianQuadrature::calculate_integrand(
                &IntegrationType::Mass,
                &element_type,
                &nodes,
                &material_property,
                3, // mesh_dim
            ).unwrap();

            // Find maximum degree in the integrand matrix
            let max_degree = GaussianQuadrature::find_max_degree_in_matrix(&integrand);
            
            // Mass matrix integrand for line element should have max degree 2
            // (from N_i * N_j terms where both are linear)
            assert_eq!(max_degree, 2);
        }

        #[test]
        fn test_1d_line_gauss_legendre_error_function() {
            // Test the error function directly with known polynomials
            
            // Quadratic polynomial: x^2 coefficients in graded lex order for 1D
            // For 1D, basis is [1, x, x^2, x³, ...]
            let quadratic_poly = vec![0.0, 0.0, 1.0]; // x^2
            
            // Error for 1-point Gauss on quadratic polynomial should be non-zero
            let error_1pt = GaussianQuadrature::gauss_legendre_error(&quadratic_poly, 1, 1).unwrap();
            assert!(error_1pt > 0.0, "1-point Gauss should have error for quadratic polynomial");
            
            // Error for 2-point Gauss on quadratic polynomial should be zero (exact)
            let error_2pt = GaussianQuadrature::gauss_legendre_error(&quadratic_poly, 2, 1).unwrap();
            assert!(error_2pt.abs() < 1e-12, "2-point Gauss should be exact for quadratic polynomial");
        }

        #[test]
        fn test_1d_line_mesh_level_analysis() {
            let (_, _, _, material_property, mesh_data) = create_1d_line_element();
            
            let tolerance = 1e-10;
            
            let result = GaussianQuadrature::find_optimal_gauss_points_number_mesh(
                &mesh_data,
                IntegrationType::Mass,
                tolerance,
                &material_property,
            ).unwrap();

            // Should analyze 1 element
            assert_eq!(result.total_elements, 1);
            assert_eq!(result.gauss_point_numbers.len(), 1);
            
            let element_report = &result.gauss_point_numbers[0];
            assert_eq!(element_report.element_id, 0); // 0-based element ID
            assert_eq!(element_report.theoretical_number, 2);
            assert_eq!(element_report.optimal_number, 2);
        }
    }

    // Additional tests for edge cases
    mod edge_case_tests {
        use super::*;

        fn create_zero_length_element() -> (Element, Vec<Node>, ElementType, MaterialProperty) {
            // Test with degenerate element (both nodes at same location)
            let element = Element {
                id: 0,
                nodes: vec![0, 1],
            };

            let nodes = vec![
                Node {
                    id: 0,
                    coordinates: vec![0.0, 0.0, 0.0],
                },
                Node {
                    id: 1,
                    coordinates: vec![0.0, 0.0, 0.0], // Same as node 0
                },
            ];

            let material_property = MaterialProperty::Scalar(vec![1.0]);

            (element, nodes, ElementType::Line, material_property)
        }

        #[test]
        fn test_zero_length_element() {
            let (element, nodes, element_type, material_property) = create_zero_length_element();

            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            // Jacobian should be zero
            let jacobian = GaussianQuadrature::calculate_jacobian_monomial(
                &nodes,
                &shape_function.derivatives,
                1, 2, 3, 1,
            ).unwrap();

            // All Jacobian components should be zero
            for i in 0..3 {
                assert!((jacobian.matrix[i][0][0]).abs() < 1e-12, 
                       "Jacobian component should be zero for degenerate element");
            }

            // Determinant should be zero
            assert!((jacobian.determinant[0]).abs() < 1e-12, 
                   "detJ should be zero for degenerate element");
        }

        #[test]
        fn test_constant_material_property() {
            // Test with non-constant material property
            let element = Element {
                id: 0,
                nodes: vec![0, 1],
            };

            let nodes = vec![
                Node {
                    id: 0,
                    coordinates: vec![0.0, 0.0, 0.0],
                },
                Node {
                    id: 1,
                    coordinates: vec![1.0, 0.0, 0.0],
                },
            ];

            // Linear density: ρ(x) = 1 + x
            // In monomial basis for 1D: [constant, x] = [1.0, 1.0]
            let material_property = MaterialProperty::Scalar(vec![1.0, 1.0]);

            let result = GaussianQuadrature::find_optimal_gauss_points(
                IntegrationType::Mass,
                &element,
                &ElementType::Line,
                &nodes,
                3,
                1e-10,
                &material_property,
            );

            // Should succeed and find appropriate Gauss points
            assert!(result.is_ok());
            let report = result.unwrap();
            
            // With linear density and quadratic shape functions, integrand becomes cubic
            // Theoretical points = ceil((3+1)/2) = 2, but might need more for tolerance
            assert!(report.theoretical_number >= 2);
        }

        #[test]
        fn test_element_node_retrieval() {
            let (element, nodes, _, _, _) = super::line_element_1d_tests::create_1d_line_element();
            
            // Test that GeometricAnalysis::get_element_nodes works correctly with 0-based indexing
            let element_nodes = GeometricAnalysis::get_element_nodes(&element, &nodes).unwrap();
            
            assert_eq!(element_nodes.len(), 2);
            assert_eq!(element_nodes[0].id, 0);
            assert_eq!(element_nodes[1].id, 1);
            assert_eq!(element_nodes[0].coordinates, vec![0.0, 0.0, 0.0]);
            assert_eq!(element_nodes[1].coordinates, vec![1.0, 0.0, 0.0]);
        }
    }
}

*/
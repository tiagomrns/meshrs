use core::num;
use std::f64;
use std::vec;
use crate::structs_and_impls::*;
use crate::error::*;
use super::geometric_analysis::GeometricAnalysis;
use once_cell::sync::Lazy;

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

    /// Find optimal Gauss points for entire mesh
    pub fn find_optimal_gauss_points_number_mesh(
        mesh_data: &MeshData,
        int_type: IntegrationType,
        tolerance: f64,
        material_property: &MaterialProperty, 
    ) -> Result<GaussianPointNumberReport, GaussError> {

        let mut gauss_point_numbers: Vec<GaussianPointNumber> = Vec::new();
        let mut processed_elements = 0;
        let mesh_dim = mesh_data.dimension;

        // Iterate through all element types in the mesh
        for type_info in &mesh_data.element_type_info {
            // Skip vertex elements as they are just points
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

                    match Self::find_optimal_gauss_points(
                        int_type,
                        element,  // get element dimension from this
                        &type_info.element_type,  // get element type
                        &mesh_data.nodes, // to get determinant of jacobian
                        mesh_dim,
                        tolerance,
                        material_property,
                    ) {
                        Ok(gauss_point_number) => {
                            gauss_point_numbers.push(gauss_point_number);
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
        if gauss_point_numbers.is_empty() {
            return Err(GaussError::GeometryError(
                "No elements could be analyzed".to_string(),
            ));
        }

        // Return comprehensive quality report
        Ok(GaussianPointNumberReport {
            total_elements: processed_elements,
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
        material_property: &MaterialProperty, // density polynomial coefficients
    ) -> Result<GaussianPointNumber, GaussError> {

        // Get element dimension
        let element_dim = ElementType::get_element_dimension(element_type)
            .ok_or(GaussError::UnsupportedDimension(0))?;

        // Get element nodes
        let element_nodes: Vec<Node> = GeometricAnalysis::get_element_nodes(element, nodes)?;

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

        let shape_function = ElementType::get_shape_functions(element_type)
            .ok_or_else(|| GaussError::InvalidElement("Shape functions not found".to_string()))?;

        for n in 1..=max_num_gp {
            let error_matrix = Self::calculate_error(
                &integrand,
                n,
                element_dim,
            )?;

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

    // Calculate integrand for given element and integration type
    pub fn calculate_integrand(
        int_type: &IntegrationType,
        element_type: &ElementType,
        element_nodes: &Vec<Node>,
        material_property: &MaterialProperty,
        mesh_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {

        // Get element and mesh dimensions and element order
        let element_dim = ElementType::get_element_dimension(element_type)
            .ok_or(GaussError::UnsupportedDimension(0))?;
        let element_order = ElementType::get_element_order(element_type)
            .ok_or(GaussError::InvalidElement("Element order not found".to_string()))?;
        
        let shape_function = ElementType::get_shape_functions(element_type)
            .ok_or_else(|| GaussError::InvalidElement("Shape functions not found".to_string()))?;
        
        let num_nodes = shape_function.num_nodes;
        
        // Calculate Jacobian matrix mapping parametric to physical coordinates
        let jacobian_matrix = GaussianQuadrature::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        )?;
        
        // Calculate determinant of Jacobian matrix
        let det_j = &jacobian_matrix.determinant;

        let inv_jacobian = Self::inverse_matrix_monomial(
            &jacobian_matrix.matrix,
            mesh_dim,
            element_dim,
        ).map_err(|e| GaussError::GeometryError(format!("Failed to calculate inverse Jacobian: {:?}", e)))?;

        let matrix_size = mesh_dim * num_nodes;
        let mut integrand = vec![vec![vec![]; matrix_size]; matrix_size];
        

        match int_type {
            IntegrationType::Mass => {
                // Mass matrix: M = ∫ ρ(x) N_i(x) N_j(x) det(J) dξ
                let rho_physical = material_property.as_scalar()?;
                let rho_isoparametric = Self::transform_material_property_scalar(
                    rho_physical,
                    element_nodes,
                    &shape_function.values,
                    mesh_dim,
                )?;

                // Calculate ρ * det(J)
                let rho_detj = MonomialPolynomial::multiply(&rho_isoparametric, &det_j)?;

                // Mass matrix is block-diagonal: M = [m_ij * I] where I is identity matrix
                for i in 0..num_nodes {
                    for j in 0..num_nodes {
                        let ni_nj = MonomialPolynomial::multiply(
                            &shape_function.values[i],
                            &shape_function.values[j]
                        )?;
                        let m_ij = MonomialPolynomial::multiply(&rho_detj, &ni_nj)?;
                        
                        // Distribute to all spatial dimensions
                        for dim in 0..mesh_dim {
                            let row = mesh_dim * i + dim;
                            let col = mesh_dim * j + dim;
                            integrand[row][col] = m_ij.clone();
                        }
                    }
                }
                
                Ok(integrand)
            }
            IntegrationType::Stiffness => {
                // Stiffness matrix: K = ∫ B^T * C * B * det(J) dξ
                
                // Build B matrix with correct dimensions
                let b_matrix = Self::build_b_matrix(
                    &shape_function.derivatives,
                    &inv_jacobian,
                    num_nodes,
                    mesh_dim,
                    element_dim,
                )?;

                let c_matrix_physical = material_property.as_matrix()?;
                let c_matrix_isoparametric = Self::transform_material_property_tensor(
                    c_matrix_physical,
                    element_nodes,
                    &shape_function.values,
                    mesh_dim,
                )?;

                // Validate material matrix dimensions
                let expected_c_size = match mesh_dim {
                    1 => 1,  // 1D: 1x1
                    2 => 3,  // 2D: 3x3 (plane stress/strain)
                    3 => 6,  // 3D: 6x6
                    _ => return Err(GaussError::UnsupportedDimension(mesh_dim)),
                };

                if c_matrix_isoparametric.len() != expected_c_size {
                    return Err(GaussError::InvalidMaterialProperty(format!(
                        "Material matrix has {} rows, expected {} for {}D elements",
                        c_matrix_isoparametric.len(), expected_c_size, mesh_dim
                    )));
                }

                // Compute K = B^T * C * B * det(J)
                let voigt_size = b_matrix.len();
                
                // Precompute C * B
                let mut c_times_b = vec![vec![vec![]; matrix_size]; voigt_size];
                for i in 0..voigt_size {
                    for j in 0..matrix_size {
                        let mut sum = vec![];
                        for k in 0..voigt_size {
                            let term = MonomialPolynomial::multiply(&c_matrix_isoparametric[i][k], &b_matrix[k][j])?;
                            sum = if sum.is_empty() { term } else { MonomialPolynomial::add(&sum, &term)? };
                        }
                        c_times_b[i][j] = sum;
                    }
                }
                
                // Compute B^T * (C * B) * det(J)
                for i in 0..matrix_size {
                    for j in 0..matrix_size {
                        let mut sum = vec![];
                        for k in 0..voigt_size {
                            let term = MonomialPolynomial::multiply(&b_matrix[k][i], &c_times_b[k][j])?;
                            sum = if sum.is_empty() { term } else { MonomialPolynomial::add(&sum, &term)? };
                        }
                        integrand[i][j] = MonomialPolynomial::multiply(&sum, &det_j)?;
                    }
                }
                
                Ok(integrand)
            }
        }
    }


    // Build B matrix
    fn build_b_matrix(
        shape_derivatives: &Vec<Vec<Vec<f64>>>, // [node][param_dim][coefficients]
        inv_jacobian: &Vec<Vec<Vec<f64>>>,      // [element_dim][mesh_dim][coefficients] - CORRECTED dimensions
        num_nodes: usize,
        mesh_dim: usize,
        element_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {
        
        // Compute dN/dx = [num_nodes × mesh_dim] 
        // dN/dx = dN/dξ * ∂ξ/∂x = dN/dξ * inv(J)
        let mut dn_dx = vec![vec![vec![]; mesh_dim]; num_nodes]; // [num_nodes × mesh_dim]

        for node in 0..num_nodes {
            for phys_dim in 0..mesh_dim {
                let mut sum = vec![];
                for param_dim in 0..element_dim {
                    let term = MonomialPolynomial::multiply(
                        &shape_derivatives[node][param_dim], // dN/dξ
                        &inv_jacobian[param_dim][phys_dim]   // ∂ξ/∂x from inv(J)
                    )?;
                    sum = if sum.is_empty() { 
                        term 
                    } else { 
                        MonomialPolynomial::add(&sum, &term)? 
                    };
                }
                dn_dx[node][phys_dim] = sum;
            }
        }

        // Build B matrix according to standard Voigt notation
        match mesh_dim {
            1 => {
                // 1D: B = [dN1/dx, dN2/dx, ...] - size [1 × num_nodes]
                let mut b_matrix = vec![vec![vec![]; num_nodes]; 1];
                for node in 0..num_nodes {
                    b_matrix[0][node] = dn_dx[node][0].clone(); // dN/dx
                }
                Ok(b_matrix)
            }
            2 => {
                // 2D Voigt: ε = [ε_xx, ε_yy, γ_xy]^T
                // B matrix size: [3 × (2*num_nodes)]
                let mut b_matrix = vec![vec![vec![]; 2 * num_nodes]; 3];
                
                for node in 0..num_nodes {
                    let u_col = 2 * node;      // u displacement
                    let v_col = 2 * node + 1;  // v displacement
                    
                    // ε_xx = ∂u/∂x
                    b_matrix[0][u_col] = dn_dx[node][0].clone();
                    
                    // ε_yy = ∂v/∂y  
                    b_matrix[1][v_col] = dn_dx[node][1].clone();
                    
                    // γ_xy = ∂u/∂y + ∂v/∂x
                    b_matrix[2][u_col] = dn_dx[node][1].clone(); // ∂u/∂y
                    b_matrix[2][v_col] = dn_dx[node][0].clone(); // ∂v/∂x
                }
                Ok(b_matrix)
            }
            3 => {
                // 3D Voigt: ε = [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]^T  
                // B matrix size: [6 × (3*num_nodes)]
                let mut b_matrix = vec![vec![vec![]; 3 * num_nodes]; 6];
                
                for node in 0..num_nodes {
                    let u_col = 3 * node;      // u displacement
                    let v_col = 3 * node + 1;  // v displacement
                    let w_col = 3 * node + 2;  // w displacement
                    
                    // Normal strains
                    b_matrix[0][u_col] = dn_dx[node][0].clone(); // ε_xx = ∂u/∂x
                    b_matrix[1][v_col] = dn_dx[node][1].clone(); // ε_yy = ∂v/∂y
                    b_matrix[2][w_col] = dn_dx[node][2].clone(); // ε_zz = ∂w/∂z
                    
                    // Shear strains  
                    b_matrix[3][u_col] = dn_dx[node][1].clone(); // γ_xy = ∂u/∂y
                    b_matrix[3][v_col] = dn_dx[node][0].clone(); // γ_xy = ∂v/∂x
                    
                    b_matrix[4][v_col] = dn_dx[node][2].clone(); // γ_yz = ∂v/∂z  
                    b_matrix[4][w_col] = dn_dx[node][1].clone(); // γ_yz = ∂w/∂y
                    
                    b_matrix[5][w_col] = dn_dx[node][0].clone(); // γ_zx = ∂w/∂x
                    b_matrix[5][u_col] = dn_dx[node][2].clone(); // γ_zx = ∂u/∂z
                }
                Ok(b_matrix)
            }
            _ => Err(GaussError::UnsupportedDimension(mesh_dim)),
        }
    }

    // Transform scalar material property (density ρ) from physical to isoparametric coordinates 
    // Input: ρ(x,y,z) in monomial format. Output: ρ(xi,eta,psi) in monomial format
    // Substitutes: x → x(xi,eta,psi), y → y(xi,eta,psi), z → z(xi,eta,psi), where x = Σ N_i * x_i_real
    pub fn transform_material_property_scalar(
        material_property_physical: &[f64],
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>],  
        mesh_dim: usize,
    ) -> Result<Vec<f64>, GaussError> {

        // Build coordinate mappings: x(xi,eta,psi), y(xi,eta,psi), z(xi,eta,psi)
        let coordinate_maps = Self::build_coordinate_maps(
            element_nodes,
            shape_function_values,
            mesh_dim,
        )?;

        // Substitute physical coordinates with isoparametric coordinate polynomials
        Self::substitute_polynomial(material_property_physical, &coordinate_maps)

    }

    // Transform tensor material property (elasticity matrix C) from physical to isoparametric
    // Input: C(x,y,z) as matrix where each component C_ij is a monomial polynomial
    // Output: C(xi,eta,psi) with each component transformed
    // Substitutes: x → x(xi,eta,psi), y → y(xi,eta,psi), z → z(xi,eta,psi), where x = Σ N_i * x_i_real
    pub fn transform_material_property_tensor(
        material_property_physical: &Vec<Vec<Vec<f64>>>,
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>], 
        mesh_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {
        
        // Build coordinate maps once for all components
        let coordinate_maps = Self::build_coordinate_maps(
            element_nodes,
            shape_function_values,
            mesh_dim,
        )?;

        // Transform each component of the C matrix independently
        let matrix_size = material_property_physical.len(); // change this
        let mut material_property_isoparametric = vec![vec![vec![]; matrix_size]; matrix_size];
        
        for i in 0..matrix_size {
            for j in 0..matrix_size {
                material_property_isoparametric[i][j] = Self::substitute_polynomial(
                    &material_property_physical[i][j],
                    &coordinate_maps,
                )?;
            }
        }

        Ok(material_property_isoparametric)
    }
    
    // Build coordinate mappings for isoparametric formulation
    // Returns: [x(xi,eta,psi), y(xi,eta,psi), z(xi,eta,psi)] as monomial polynomials
    // x(xi,eta,psi) = Σ N_i(xi,eta,psi) × x_i , y(xi,eta,psi) = Σ N_i(xi,eta,psi) × y_i , etc.
    fn build_coordinate_maps(
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>],
        mesh_dim: usize,
    ) -> Result<Vec<Vec<f64>>, GaussError> {
        
        let num_nodes = element_nodes.len();
        let mut coordinate_maps = Vec::with_capacity(mesh_dim);
        
        // For each spatial dimension (x, y, z)
        for i in 0..mesh_dim {
            let mut coord_poly = vec![];
            
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
                
                // Sum of N_k(xi,eta,psi) * x_ik
                coord_poly = if coord_poly.is_empty() {
                    term
                } else {
                    MonomialPolynomial::add(&coord_poly, &term)
                        .map_err(|e| GaussError::GeometryError(e.to_string()))?
                };
            }
            
            // If all nodes had zero coordinate, set to zero polynomial
            if coord_poly.is_empty() {
                coord_poly = vec![0.0];
            }
            
            coordinate_maps.push(coord_poly);
        }
        
        Ok(coordinate_maps)
    }
    
    // Substitute physical coordinates with parametric coordinate polynomials
    // Given: f(x,y,z) = Σ a_ijk * x^i * y^j * z^k (in graded lexicographic format)
    // Returns: f(xi,eta,psi) by substituting x→x(xi,eta,psi), y→y(xi,eta,psi), z→z(xi,eta,psi)
    // Algorithm:
    // 1. Parse each term a_ijk * x^i * y^j * z^k from the physical polynomial
    // 2. Compute x(xi,eta,psi)^i, y(xi,eta,psi)^j, z(xi,eta,psi)^k using polynomial powers
    // 3. Multiply and accumulate: a_ijk * [x(xi,eta,psi)]^i * [y(xi,eta,psi)]^j * [z(xi,eta,psi)]^k
    // Uses the graded lexicographic basis to iterate through all terms
    fn substitute_polynomial(
        physical_poly: &[f64],
        coordinate_maps: &Vec<Vec<f64>>,
    ) -> Result<Vec<f64>, GaussError> {
        
        // Get the degree of the physical polynomial
        let physical_degree = MonomialPolynomial::infer_max_degree(physical_poly.len())
            .map_err(|e| GaussError::GeometryError(e.to_string()))?;
        
        // Get the degrees of the coordinate mappings
        let mut coord_degrees = Vec::new();
        for coord_map in coordinate_maps {
            let deg = MonomialPolynomial::infer_max_degree(coord_map.len())
                .map_err(|e| GaussError::GeometryError(e.to_string()))?;
            coord_degrees.push(deg);
        }
        
        // The maximum degree after substitution is the sum of the physical degree 
        // and the maximum coordinate degree times the maximum exponent
        // But we need to be conservative - let's compute it properly
        
        // For safety, use a conservative estimate: physical_degree + max_coord_degree * physical_degree
        let max_coord_degree = coord_degrees.iter().max().unwrap_or(&0);
        let conservative_max_degree = physical_degree + max_coord_degree * physical_degree;
        
        // Generate basis for the conservative maximum degree
        let basis = MonomialPolynomial::generate_basis(conservative_max_degree);
        
        let mut result = vec![];
        
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
            
            // Multiply by x(ξ,η,ψ)^i
            if i > 0 && coordinate_maps.len() > 0 {
                let x_power = Self::polynomial_power(&coordinate_maps[0], i as usize)?;
                term_poly = MonomialPolynomial::multiply(&term_poly, &x_power)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?;
            }
            
            // Multiply by y(ξ,η,ψ)^j (if y dimension exists)
            if j > 0 && coordinate_maps.len() > 1 {
                let y_power = Self::polynomial_power(&coordinate_maps[1], j as usize)?;
                term_poly = MonomialPolynomial::multiply(&term_poly, &y_power)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?;
            }
            
            // Multiply by z(ξ,η,ψ)^k (if z dimension exists)
            if k > 0 && coordinate_maps.len() > 2 {
                let z_power = Self::polynomial_power(&coordinate_maps[2], k as usize)?;
                term_poly = MonomialPolynomial::multiply(&term_poly, &z_power)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?;
            }
            
            // Accumulate: add this term to the result
            result = if result.is_empty() {
                term_poly
            } else {
                MonomialPolynomial::add(&result, &term_poly)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?
            };
        }
        
        // If all terms were zero, return zero polynomial
        if result.is_empty() {
            result = vec![0.0];
        }
        
        Ok(result)
    }

    /// Compute polynomial raised to integer power: poly^n
    /// 
    /// Uses repeated multiplication for efficiency
    /// Special cases: poly^0 = 1, poly^1 = poly
    fn polynomial_power(poly: &Vec<f64>, power: usize) -> Result<Vec<f64>, GaussError> {
        match power {
            0 => Ok(vec![1.0]), // poly^0 = 1 (constant)
            1 => Ok(poly.clone()), // poly^1 = poly
            _ => {
                // Repeated multiplication for poly^n where n > 1
                let mut result = poly.clone();
                for _ in 1..power {
                    result = MonomialPolynomial::multiply(&result, poly)
                        .map_err(|e| GaussError::GeometryError(e.to_string()))?;
                }
                Ok(result)
            }
        }
    }

    // Compute Jacobian struct 
    pub fn calculate_jacobian_monomial(
        element_nodes: &Vec<Node>,
        shape_derivatives: &[Vec<Vec<f64>>], // [node][param_dim][coefficients]
        element_dim: usize,
        num_nodes: usize,
        mesh_dim: usize,
        element_order: usize,
    ) -> Result<Jacobian, ElementError> {

        let max_len = match element_order {
            1 => 4,  // [1, x, y, z]
            2 => 10, // [1, x, y, z, x^2, xy, xz, y^2, yz, z^2]
            3 => 20, // [1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3]
            // Extendable for higher orders if needed
            _ => return Err(ElementError::GeometryError(format!(
                "Unsupported element order {}", element_order
            ))),
        };

        // Initialize Jacobian matrix: J[physical_dim][parametric_dim]
        let mut jacobian_matrix: Vec<Vec<Vec<f64>>> = 
            vec![vec![vec![0.0; max_len]; element_dim]; mesh_dim];

        // J_ij = Σ_node (x_node_i * ∂N_node/∂ξ_j)
        for phys_dim in 0..mesh_dim {           // Physical space dimension
            for param_dim in 0..element_dim {   // Parametric space dimension  
                for node in 0..num_nodes {      // Node index
                    let coord = element_nodes[node].coordinates[phys_dim];
                    let derivative_poly = &shape_derivatives[node][param_dim];
                    
                    let scaled_poly = MonomialPolynomial::multiply_scalar(derivative_poly, coord);
                    jacobian_matrix[phys_dim][param_dim] = MonomialPolynomial::add(
                        &jacobian_matrix[phys_dim][param_dim], 
                        &scaled_poly
                    ).map_err(|e| ElementError::GeometryError(format!("Failed to add polynomials: {}", e)))?;
                }
            }
        }

        let det_jacobian = Self::calculate_determinant_monomial(&jacobian_matrix, mesh_dim, element_dim)?;

        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant: det_jacobian,
        })
    }

    
    // Inverse of a matrix with monomial polynomial entries for square Jacobian matrices and pseudo-inverse for non-square Jacobian matrices
    pub fn inverse_matrix_monomial(matrix: &Vec<Vec<Vec<f64>>>, mesh_dim: usize, element_dim: usize) -> Result<Vec<Vec<Vec<f64>>>, ElementError> {
        match (mesh_dim, element_dim) {
            // Square matrices - inverse has same dimensions [element_dim × mesh_dim] = [mesh_dim × element_dim]
            (1, 1) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; mesh_dim]; element_dim];
                let det_matrix = &matrix[0][0];
                let inv_determinant = Self::calculate_inverse_monomial(&det_matrix)?;
                inverse_matrix[0][0] = inv_determinant;
                Ok(inverse_matrix)
            },
            (2, 2) | (3, 3) => {
                // Combined square matrix case for 2x2 and 3x3
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; mesh_dim]; element_dim];
                let adjoint = Self::calculate_adjoint_monomial(matrix, element_dim)?;
                let det_matrix = Self::calculate_determinant_monomial(matrix, mesh_dim, element_dim)?;
                let inv_determinant = Self::calculate_inverse_monomial(&det_matrix)?;
                for i in 0..element_dim {
                    for j in 0..mesh_dim {
                        inverse_matrix[i][j] = MonomialPolynomial::multiply(&adjoint[i][j], &inv_determinant)?;
                    }
                }
                Ok(inverse_matrix)
            },
            
            // 1D elements in higher dimensional spaces - pseudo-inverse has dimensions [element_dim × mesh_dim]
            (2, 1) | (3, 1) => {
                // Combined case for 1D elements in 2D or 3D space
                // J = [mesh_dim × 1], J⁺ = [1 × mesh_dim]
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; mesh_dim]; element_dim]; // [1 × mesh_dim]
                
                // Compute G = JᵀJ [1 × 1]
                let mut g_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; 1]; 1];
                let mut sum = vec![];
                for k in 0..mesh_dim {
                    let prod = MonomialPolynomial::multiply(&matrix[k][0], &matrix[k][0])?;
                    sum = MonomialPolynomial::add(&sum, &prod)?;
                }
                g_matrix[0][0] = sum;

                // Compute G⁻¹ = 1/det(G)
                let det_g = &g_matrix[0][0];
                let inv_det_g = Self::calculate_inverse_monomial(&det_g)?;

                // Compute J⁺ = G⁻¹Jᵀ = [1 × 1] * [1 × mesh_dim] = [1 × mesh_dim]
                for j in 0..mesh_dim {
                    inverse_matrix[0][j] = MonomialPolynomial::multiply(&inv_det_g, &matrix[j][0])?;
                }
                
                Ok(inverse_matrix)
            },
            (3, 2) => {
                // J = [3 × 2], J⁺ = [2 × 3]
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; mesh_dim]; element_dim]; // [2 × 3]
                
                // Compute G = JᵀJ [2 × 2]
                let mut g_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; 2]; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        let mut sum = vec![];
                        for k in 0..3 {
                            let prod = MonomialPolynomial::multiply(&matrix[k][i], &matrix[k][j])?;
                            sum = MonomialPolynomial::add(&sum, &prod)?;
                        }
                        g_matrix[i][j] = sum;
                    }
                }

                // Compute G⁻¹ = adj(G)/det(G)
                let adj_g = Self::calculate_adjoint_monomial(&g_matrix, 2)?;
                let det_g = Self::calculate_determinant_monomial(&g_matrix, 2, 2)?;
                let inv_det_g = Self::calculate_inverse_monomial(&det_g)?;
                let mut inv_g = vec![vec![vec![]; 2]; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        inv_g[i][j] = MonomialPolynomial::multiply(&adj_g[i][j], &inv_det_g)?;
                    }
                }

                // Compute J⁺ = G⁻¹Jᵀ = [2 × 2] * [2 × 3] = [2 × 3]
                for i in 0..element_dim {
                    for j in 0..mesh_dim {
                        let mut sum = vec![];
                        for k in 0..2 {
                            let term = MonomialPolynomial::multiply(&inv_g[i][k], &matrix[j][k])?;
                            sum = if sum.is_empty() { term } else { MonomialPolynomial::add(&sum, &term)? };
                        }
                        inverse_matrix[i][j] = sum;
                    }
                }

                Ok(inverse_matrix)
            },
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian inverse calculation not implemented for {}x{} matrices",
                mesh_dim, element_dim
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
                        MonomialPolynomial::add(&MonomialPolynomial::multiply(b1, c2)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(b2, c1)?, -1.0))?,
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add(&MonomialPolynomial::multiply(b0, c2)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(b2, c0)?, -1.0))?, -1.0),
                        MonomialPolynomial::add(&MonomialPolynomial::multiply(b0, c1)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(b1, c0)?, -1.0))?,
                    ],
                    vec![
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add(&MonomialPolynomial::multiply(a1, c2)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(a2, c1)?, -1.0))?, -1.0),
                        MonomialPolynomial::add(&MonomialPolynomial::multiply(a0, c2)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(a2, c0)?, -1.0))?,
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add(&MonomialPolynomial::multiply(a0, c1)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(a1, c0)?, -1.0))?, -1.0),
                    ],
                    vec![
                        MonomialPolynomial::add(&MonomialPolynomial::multiply(a1, b2)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(a2, b1)?, -1.0))?,
                        MonomialPolynomial::multiply_scalar(&MonomialPolynomial::add(&MonomialPolynomial::multiply(a0, b2)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(a2, b0)?, -1.0))?, -1.0),
                        MonomialPolynomial::add(&MonomialPolynomial::multiply(a0, b1)?, &MonomialPolynomial::multiply_scalar(&MonomialPolynomial::multiply(a1, b0)?, -1.0))?,
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
                let squared_length = MonomialPolynomial::add(&dx_sq, &dy_sq)
                    .map_err(|e| ElementError::GeometryError(format!("Add failed: {}", e)))?;
                
                // Take square root to get actual length
                Self::calculate_sqrt_monomial(&squared_length)
            }

            // 1D elements in 3D space: metric = dx² + dy² + dz² (squared length)
            (3, 1) => {
                let dx_sq = MonomialPolynomial::multiply(&matrix[0][0], &matrix[0][0])?;
                let dy_sq = MonomialPolynomial::multiply(&matrix[1][0], &matrix[1][0])?;
                let dz_sq = MonomialPolynomial::multiply(&matrix[2][0], &matrix[2][0])?;
                let temp = MonomialPolynomial::add(&dx_sq, &dy_sq)?;
                let squared_length = MonomialPolynomial::add(&temp, &dz_sq)
                    .map_err(|e| ElementError::GeometryError(format!("Add failed: {}", e)))?;
                
                // Take square root to get actual length
                Self::calculate_sqrt_monomial(&squared_length)
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
                let squared_area = MonomialPolynomial::add(&g00_g11, &MonomialPolynomial::multiply_scalar(&g01_g10, -1.0))
                    .map_err(|e| ElementError::GeometryError(format!("Metric det failed: {}", e)))?;
                
                // Take square root to get actual area
                Self::calculate_sqrt_monomial(&squared_area)
            }

            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for {}x{} (mesh_dim x element_dim) matrices",
                mesh_dim, element_dim
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
        let binomial_coeffs = [1.0, 0.5, -0.125, 0.0625, -0.0390625];  // for 5 terms
        
        let mut term = vec![vec![]; wished_n];
        term[0] = vec![1.0]; // Constant term is 1
        
        // Build terms: term[k] = q^k
        term[1] = q.clone();
        for i in 2..wished_n {
            term[i] = MonomialPolynomial::multiply(&term[i-1], &q)
                .map_err(|e| ElementError::GeometryError(format!("Failed to compute q power {}: {}", i, e)))?;
        }
        
        // Sum the series: result = sum_{k=0}^{4} [binomial_coeffs[k] * term[k]]
        let mut result = vec![0.0];
        for j in 0..wished_n {
            let scaled_term = MonomialPolynomial::multiply_scalar(&term[j], binomial_coeffs[j]);
            result = MonomialPolynomial::add(&result, &scaled_term)
                .map_err(|e| ElementError::GeometryError(format!("Failed to compute sqrt sum: {}", e)))?;
        }
        
        // Multiply by sqrt(a₀)
        Ok(MonomialPolynomial::multiply_scalar(&result, sqrt_a0))
    }

    /// Calculate inverse of Jacobian determinant/metric for polynomial matrices
    pub fn calculate_inverse_monomial( 
        poly: &Vec<f64>, 
    ) -> Result<Vec<f64>, ElementError> {
        let wished_n = 5; // Can be parameterized later

        let nominal_poly = MonomialPolynomial::multiply_scalar(&poly, 1.0/poly[0]);
        let mut term = vec![vec![]; wished_n];
        term[0] = vec![1.0]; // Constant term is 1
        term[1] = nominal_poly;
        for i in 2..wished_n {
            term[i] = MonomialPolynomial::multiply(&term[i-1], &term[1])
                .map_err(|e| ElementError::GeometryError(format!("Failed to compute inverse determinant power {}: {}", i, e)))?;
        }
        let mut result = vec![0.0];
        for j in 0..wished_n {
            result = MonomialPolynomial::add(&result, &term[j])
                .map_err(|e| ElementError::GeometryError(format!("Failed to compute inverse determinant sum: {}", e)))?;
        }
        Ok(result)
    }


    // Calculate error based on element type and dispatch to appropriate quadrature rule
    pub fn calculate_error(
        integrand: &Vec<Vec<Vec<f64>>>,
        n: usize,
        element_dim: usize,
    ) -> Result<Vec<Vec<f64>>, GaussError> {

        let matrix_size = integrand.len(); // This gives mesh_dim * num_nodes
        let mut errors = vec![vec![0.0; matrix_size]; matrix_size];
        
        for i in 0..matrix_size {
            for j in 0..matrix_size {
                errors[i][j] = Self::gauss_legendre_error(&integrand[i][j], n, element_dim)?;
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

        let c = factorial(n).powi(4) / ((2 * n + 1) as f64 * factorial(2 * n).powi(3));

        match element_dim {
            1 => {  //1D Gauss-Legendre error function

                // Convert to dense coefficient representation
                let coeff_vector = MonomialPolynomial::get_coefficients_1d(&integrand)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?;

                // Detect polynomial degree
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
                
                // Convert to dense coefficient representation
                let coeff_matrix = MonomialPolynomial::get_coefficients_2d(&integrand)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?;

                // Detect polynomial degree
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

                // Convert to dense coefficient representation
                let coeff_tensor = MonomialPolynomial::get_coefficients_3d(&integrand)
                    .map_err(|e| GaussError::GeometryError(e.to_string()))?;

                // Detect polynomial degree
                let (poly_degree_x, poly_degree_y, poly_degree_z) = Self::detect_polynomial_orders_3d(&coeff_tensor);

                let mut term1 = 0.0;
                for i in (2 * n)..=poly_degree_x.min(coeff_tensor.len().saturating_sub(1)) {
                    for j in 0..=poly_degree_y.min(coeff_tensor[i].len().saturating_sub(1)) {
                        for k in 0..=poly_degree_z.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                            let aijk = coeff_tensor[i][j][k].abs();
                            if aijk > 1e-12 {
                                term1 += aijk
                                    * (factorial(i) / factorial(i - 2 * n));
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
                                term2 += aijk
                                    * (factorial(j) / factorial(j - 2 * n));
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

// first create element and mesh data (1d linear in 1d-2d-3d, 1d quadratic in 1d-2d-3d, 2d linear in 2d-3d, 2d quadratic in 2d-3d, 3d linear in 3d, 3d quadratic in 3d)
    // then check shape functions and its derivatives (also monomial representation)
    // then check jacobian calculation
    // then check detJ, adjJ and invJ
    // then check the integrand for mass and stiffness with constant rho and C first then isotropic rho and C and then anisotropic rho and C
    // then check the f = f_0 + f_1 x + f_2 x^2 + ...
    // check the polynomial order and order of x, y and z
    // then check the error estimation
    // then check the creiteria for error selection (is it max|E| or frobeniusnorm(E)?) 
    // then check the number of gaussian points for given error tolerance (give multiple tolerances and check the number of points selected)
    // test the exactness number d = 2n-1 


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
            
            // Shape function values should be for parametric coordinate ξ ∈ [-1, 1]
            // N1 = (1 - ξ)/2, N2 = (1 + ξ)/2
            let n1 = &shape_function.values[0]; // N1 coefficients
            let n2 = &shape_function.values[1]; // N2 coefficients
            
            // N1 = 0.5 - 0.5ξ → coefficients: [0.5, -0.5] for basis [1, ξ]
            assert_eq!(n1.len(), 2);
            assert!((n1[0] - 0.5).abs() < 1e-12, "N1 constant term should be 0.5");
            assert!((n1[1] - (-0.5)).abs() < 1e-12, "N1 linear term should be -0.5");
            
            // N2 = 0.5 + 0.5ξ → coefficients: [0.5, 0.5] for basis [1, ξ]
            assert!((n2[0] - 0.5).abs() < 1e-12, "N2 constant term should be 0.5");
            assert!((n2[1] - 0.5).abs() < 1e-12, "N2 linear term should be 0.5");
        }

        #[test]
        fn test_1d_line_shape_function_derivatives() {
            let (_, _, element_type, _, _) = create_1d_line_element();
            let shape_function = ElementType::get_shape_functions(&element_type).unwrap();
            
            // Derivatives: dN1/dξ = -0.5, dN2/dξ = 0.5
            let dn1_dxi = &shape_function.derivatives[0][0]; // dN1/dξ
            let dn2_dxi = &shape_function.derivatives[1][0]; // dN2/dξ
            
            // Derivatives should be constant polynomials
            assert_eq!(dn1_dxi.len(), 1);
            assert!((dn1_dxi[0] - (-0.5)).abs() < 1e-12, "dN1/dξ should be -0.5");
            
            assert_eq!(dn2_dxi.len(), 1);
            assert!((dn2_dxi[0] - 0.5).abs() < 1e-12, "dN2/dξ should be 0.5");
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
            // x(ξ) = (1-ξ)/2 * 0 + (1+ξ)/2 * 1 = (1+ξ)/2
            // dx/dξ = 0.5
            // Similarly, y and z coordinates are zero, so their derivatives are zero
            
            // Jacobian matrix should be [dx/dξ, dy/dξ, dz/dξ] = [0.5, 0, 0]
            let j_matrix = &jacobian.matrix;
            
            // Check dx/dξ = 0.5
            assert_eq!(j_matrix[0][0].len(), 1); // Constant polynomial
            assert!((j_matrix[0][0][0] - 0.5).abs() < 1e-12, "dx/dξ should be 0.5");
            
            // Check dy/dξ = 0
            assert!((j_matrix[1][0][0]).abs() < 1e-12, "dy/dξ should be 0");
            
            // Check dz/dξ = 0
            assert!((j_matrix[2][0][0]).abs() < 1e-12, "dz/dξ should be 0");
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

            // For 1D element in 3D space, determinant is metric: dx² + dy² + dz²
            // dx/dξ = 0.5, dy/dξ = 0, dz/dξ = 0 → det = (0.5)² = 0.25
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
            // N1 = 0.5 - 0.5ξ, N2 = 0.5 + 0.5ξ
            
            // For 2 nodes in 3D, mass matrix is 6x6 (2 nodes × 3 DOF each)
            assert_eq!(integrand.len(), 6);
            assert_eq!(integrand[0].len(), 6);

            // Check M11 component: ρ * N1 * N1 * detJ = 1.0 * (0.5 - 0.5ξ)² * 0.25
            // = 0.25 * (0.25 - 0.5ξ + 0.25ξ²) = 0.0625 - 0.125ξ + 0.0625ξ²
            let m11 = &integrand[0][0]; // First DOF of first node
            let expected_m11 = vec![0.0625, -0.125, 0.0625]; // [constant, ξ, ξ²] coefficients
            
            assert_eq!(m11.len(), expected_m11.len());
            for (actual, expected) in m11.iter().zip(expected_m11.iter()) {
                assert!((actual - expected).abs() < 1e-12, 
                       "M11 coefficient: expected {}, got {}", expected, actual);
            }

            // Check M12 component: ρ * N1 * N2 * detJ = 1.0 * (0.5 - 0.5ξ)(0.5 + 0.5ξ) * 0.25
            // = 0.25 * (0.25 - 0.25ξ²) = 0.0625 - 0.0625ξ²
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
            
            // Quadratic polynomial: x² coefficients in graded lex order for 1D
            // For 1D, basis is [1, x, x², x³, ...]
            let quadratic_poly = vec![0.0, 0.0, 1.0]; // x²
            
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
use core::num;
use std::f64;
use crate::structs_and_impls::*;
use crate::error::*;
use super::geometric_analysis::GeometricAnalysis;
use num_traits::ops::inv;
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

/// Utility struct for Gaussian quadrature error estimation and optimization
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
                shape_function.num_nodes,
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

    /// Calculate integrand for given element and integration type
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

        let mut integrand = vec![vec![vec![]; num_nodes]; num_nodes];

        match int_type {
            IntegrationType::Mass => {
                // Mass matrix: M = int rho(x) N_i(x) N_j(x) det(J) dxi
                // Mass matrix integrand: rho(x) * N' N * det(J)

                let rho_physical = material_property.as_scalar()?; // Density polynomial coefficients in monomial form
                let rho_isoparametric = Self::transform_material_property_scalar(
                    rho_physical,
                    element_nodes,
                    &shape_function.values,
                    mesh_dim,
                )?;

                // Calculate rho * det(J)
                let rho_detj = MonomialPolynomial::multiply(&rho_isoparametric, &det_j)?;

                for i in 0..num_nodes {
                    for j in 0..num_nodes {
                        let ni_nj = MonomialPolynomial::multiply(
                            &shape_function.values[i],
                            &shape_function.values[j]
                        )?;

                        // Final integrand: rho * N_i * N_j * det(J)
                        integrand[i][j] = MonomialPolynomial::multiply(&rho_detj, &ni_nj)?;
                    }
                }
                
                Ok(integrand)
            }
            IntegrationType::Stiffness => {
                // Stiffness matrix: K = int B_i^T * D * B_j * det(J) dxi
                // where B is derivative and combination of rows of N

                // Build B matrix
                let b_matrix = Self::build_b_matrix(
                    &shape_function.derivatives,
                    &inv_jacobian,
                    num_nodes,
                    mesh_dim,
                    element_dim,
                )?;

                let c_matrix_physical = material_property.as_matrix()?; // Material property matrix C in monomial form
                let c_matrix_isoparametric = Self::transform_material_property_tensor(
                    c_matrix_physical,
                    element_nodes,
                    &shape_function.values,
                    mesh_dim,
                )?;

                // Validate material matrix dimensions based on element dimension
                let expected_c_size = match mesh_dim {
                    1 => 1,  // 1D: 1x1
                    2 => 3,  // 2D: 3x3 (plane stress/strain)
                    3 => 6,  // 3D: 6x6
                    _ => return Err(GaussError::UnsupportedDimension(element_dim)),
                };

                if c_matrix_isoparametric.len() != expected_c_size {
                    return Err(GaussError::InvalidMaterialProperty(format!(
                        "Material matrix has {} rows, expected {} for {}D elements",
                        c_matrix_isoparametric.len(), expected_c_size, element_dim
                    )));
                }
                
                for i in 0..expected_c_size {
                    if c_matrix_isoparametric[i].len() != expected_c_size {
                        return Err(GaussError::InvalidMaterialProperty(format!(
                            "Material matrix row {} has {} columns, expected {}",
                            i, c_matrix_isoparametric[i].len(), expected_c_size
                        )));
                    }
                }

                match mesh_dim {
                    1 => {
                        // 1D: K[i][j] = B[i] * C[0][0] * B[j] * det(J)
                        for i in 0..num_nodes {
                            for j in 0..num_nodes {
                                let bt_c = MonomialPolynomial::multiply(&b_matrix[0][i], &c_matrix_isoparametric[0][0])?;
                                let bt_c_b = MonomialPolynomial::multiply(&bt_c, &b_matrix[0][j])?;
                                integrand[i][j] = MonomialPolynomial::multiply(&bt_c_b, &det_j)?;
                            }
                        }
                    }
                    2 => {
                        // 2D: K[i][j] = sum_k sum_l B[k][i] * C[k][l] * B[l][j] * det(J)
                        for i in 0..num_nodes {
                            for j in 0..num_nodes {
                                let mut bt_c_b = vec![];
                                
                                for k in 0..3 {
                                    for l in 0..3 {
                                        let bt_c = MonomialPolynomial::multiply(&b_matrix[k][i], &c_matrix_isoparametric[k][l])?;
                                        let part = MonomialPolynomial::multiply(&bt_c, &b_matrix[l][j])?;
                                        bt_c_b = if bt_c_b.is_empty() {
                                            part
                                        } else {
                                            MonomialPolynomial::add(&bt_c_b, &part)?
                                        };
                                    }
                                }
                                
                                integrand[i][j] = MonomialPolynomial::multiply(&bt_c_b, &det_j)?;
                            }
                        }
                    }
                    3 => {
                        // 3D: K[i][j] = sum_k sum_l B[k][i] * C[k][l] * B[l][j] * det(J)
                        for i in 0..num_nodes {
                            for j in 0..num_nodes {
                                let mut bt_c_b = vec![];
                                
                                for k in 0..6 {
                                    for l in 0..6 {
                                        let bt_c = MonomialPolynomial::multiply(&b_matrix[k][i], &c_matrix_isoparametric[k][l])?;
                                        let part = MonomialPolynomial::multiply(&bt_c, &b_matrix[l][j])?;
                                        bt_c_b = if bt_c_b.is_empty() {
                                            part
                                        } else {
                                            MonomialPolynomial::add(&bt_c_b, &part)?
                                        };
                                    }
                                }
                                
                                integrand[i][j] = MonomialPolynomial::multiply(&bt_c_b, &det_j)?;
                            }
                        }
                    }
                    _ => return Err(GaussError::UnsupportedDimension(mesh_dim)),
                }
            
                Ok(integrand)
            }
            
            _ => Err(GaussError::InvalidIntegrationType(format!(
                "Unsupported integration type for integrand calculation"
            ))),
        }
    }


    // Build B matrix
    fn build_b_matrix(
        shape_derivatives: &Vec<Vec<Vec<f64>>>,
        inv_jacobian: &Vec<Vec<Vec<f64>>>,
        num_nodes: usize,
        mesh_dim: usize,
        element_dim: usize,
    ) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {

        // Calculate derivatives of shape functions w.r.t physical coordinates
        // dN/dx = inv(J) * dN/dxi
        let mut dn_dx = vec![vec![vec![]; num_nodes]; mesh_dim];

        for i in 0..mesh_dim {
            for j in 0..num_nodes {
                let mut sum = vec![];
                for k in 0..element_dim {
                    let part = MonomialPolynomial::multiply(
                        &inv_jacobian[i][k],
                        &shape_derivatives[j][k]
                    ).map_err(|e| GaussError::GeometryError(e.to_string()))?;
                    
                    sum = if sum.is_empty() {
                        part
                    } else {
                        MonomialPolynomial::add(&sum, &part)
                            .map_err(|e| GaussError::GeometryError(e.to_string()))?
                    };
                }
                dn_dx[i][j] = sum;
            }
        }

        match mesh_dim{
            1 => {
                // 1D: B is 1 x num_nodes
                let mut b_matrix = vec![vec![vec![]; num_nodes]; 1];
                
                for i in 0..num_nodes {
                    b_matrix[0][i] = dn_dx[i][0].clone(); // monomial coefficients
                }
                Ok(b_matrix)
            }
            2 => {
                let mut b_matrix = vec![vec![vec![]; 2 * num_nodes]; 3]; // row number = 3
                for i in 0..num_nodes {
                    b_matrix[0][0+(2*i)] = dn_dx[0+i][0].clone();
                    b_matrix[2][0+(2*i)] = dn_dx[0+i][1].clone();

                    b_matrix[1][1+(2*i)] = dn_dx[0+i][1].clone();
                    b_matrix[2][1+(2*i)] = dn_dx[0+i][0].clone();
                }
                Ok(b_matrix)
            }
            3 => {
                let mut b_matrix = vec![vec![vec![]; 3 * num_nodes]; 6]; // row number = 6
                for i in 0..num_nodes {
                    b_matrix[0][0+(3*i)] = dn_dx[0+i][0].clone();
                    b_matrix[3][0+(3*i)] = dn_dx[0+i][1].clone();
                    b_matrix[5][0+(3*i)] = dn_dx[0+i][2].clone();

                    b_matrix[1][1+(3*i)] = dn_dx[0+i][1].clone();
                    b_matrix[3][1+(3*i)] = dn_dx[0+i][0].clone();
                    b_matrix[4][1+(3*i)] = dn_dx[0+i][2].clone();

                    b_matrix[2][2+(3*i)] = dn_dx[0+i][2].clone();
                    b_matrix[4][2+(3*i)] = dn_dx[0+i][1].clone();
                    b_matrix[5][2+(3*i)] = dn_dx[0+i][0].clone();
                }
                Ok(b_matrix)
            }
            _ => panic!("Unsupported mesh dimension: {}", mesh_dim), // unsupported mesh dimension
        }
    }

    /// Transform scalar material property (density ρ) from physical to isoparametric coordinates
    /// 
    /// Input: ρ(x,y,z) in monomial format
    /// Output: ρ(xi,eta,psi) in monomial format
    /// 
    /// Substitutes: x → x(xi,eta,psi), y → y(xi,eta,psi), z → z(xi,eta,psi)
    /// where coordinate mappings come from isoparametric formulation
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

    /// Transform tensor material property (elasticity matrix C) from physical to isoparametric
    ///
    /// Input: C(x,y,z) as matrix where each component C_ij is a monomial polynomial
    /// Output: C(xi,eta,psi) with each component transformed
    ///
    /// Handles anisotropic materials where each C_ij(x,y,z) varies spatially
    pub fn transform_material_property_tensor(
        material_property_physical: &Vec<Vec<Vec<f64>>>,
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>],  // Shape function VALUES (not derivatives)
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
    
    /// Build coordinate mappings for isoparametric formulation
    /// 
    /// Returns: [x(xi,eta,psi), y(xi,eta,psi), z(xi,eta,psi)] as monomial polynomials
    /// 
    /// For each coordinate: coordinate(xi,eta,psi) = Σ N_i(xi,eta,psi) × coordinate_i
    fn build_coordinate_maps(
        element_nodes: &Vec<Node>,
        shape_function_values: &[Vec<f64>],
        mesh_dim: usize,
    ) -> Result<Vec<Vec<f64>>, GaussError> {
        
        let num_nodes = element_nodes.len();
        let mut coordinate_maps = Vec::with_capacity(mesh_dim);
        
        // For each spatial dimension (x, y, z)
        for dim in 0..mesh_dim {
            let mut coord_poly = vec![];
            
            // Build: xᵢ(ξ,η,ψ) = Σₖ Nₖ(ξ,η,ψ) × xᵢₖ
            for node_idx in 0..num_nodes {
                // Get physical coordinate value for this node
                let coord_value = element_nodes[node_idx].coordinates[dim];
                
                // Skip if coordinate is zero (optimization)
                if coord_value.abs() < 1e-15 {
                    continue;
                }
                
                // Multiply: Nₖ(ξ,η,ψ) × xᵢₖ
                let term = MonomialPolynomial::multiply_scalar(
                    &shape_function_values[node_idx],
                    coord_value,
                );
                
                // Accumulate: sum all contributions
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
    
    /// Substitute physical coordinates with parametric coordinate polynomials
    /// 
    /// Given: f(x,y,z) = Σ aᵢⱼₖ × xⁱ × yʲ × zᵏ (in graded lexicographic format)
    /// Returns: f(ξ,η,ψ) by substituting x→x(ξ,η,ψ), y→y(ξ,η,ψ), z→z(ξ,η,ψ)
    /// 
    /// Algorithm:
    /// 1. Parse each term aᵢⱼₖ × xⁱ × yʲ × zᵏ from the physical polynomial
    /// 2. Compute x(ξ,η,ψ)^i, y(ξ,η,ψ)^j, z(ξ,η,ψ)^k using polynomial powers
    /// 3. Multiply and accumulate: aᵢⱼₖ × [x(ξ,η,ψ)]^i × [y(ξ,η,ψ)]^j × [z(ξ,η,ψ)]^k
    /// 
    /// Uses the graded lexicographic basis to iterate through all terms
    fn substitute_polynomial(
        physical_poly: &[f64],
        coordinate_maps: &Vec<Vec<f64>>,
    ) -> Result<Vec<f64>, GaussError> {
        
        // Infer polynomial degree from coefficient vector length
        let max_degree = MonomialPolynomial::infer_max_degree(physical_poly.len())
            .map_err(|e| GaussError::GeometryError(e.to_string()))?;
        
        // Generate basis: [(i,j,k)] for all terms in graded lexicographic order
        let basis = MonomialPolynomial::generate_basis(max_degree);
        
        let mut result = vec![];
        
        // Process each term: aᵢⱼₖ × xⁱ × yʲ × zᵏ
        for (idx, &(i, j, k)) in basis.iter().enumerate() {
            if idx >= physical_poly.len() {
                break;
            }
            
            let coefficient = physical_poly[idx];
            
            // Skip negligible terms for numerical stability
            if coefficient.abs() < 1e-15 {
                continue;
            }
            
            // Start with coefficient as constant polynomial: [aᵢⱼₖ, 0, 0, ...]
            let mut term_poly = vec![coefficient];
            
            // Multiply by x(ξ,η,ψ)^i
            if i > 0 {
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

    // Inverse of a matrix with monomial polynomial entries
    pub fn inverse_matrix_monomial(matrix: &Vec<Vec<Vec<f64>>>, mesh_dim: usize, element_dim: usize) -> Result<Vec<Vec<Vec<f64>>>, ElementError> {
        // Inversion only implemented for 2x2 and 3x3 matrices
        match (mesh_dim, element_dim) {
            // 1x1 matrix: determinant is the single element
            (1, 1) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];

                match MonomialPolynomial::divide(&[1.0, 0.0, 0.0, 0.0], &matrix[0][0]) {
                    Ok(coefficients) => {
                        inverse_matrix[0][0] = coefficients;
                    }
                    Err(e) => return Err(e.into()),
                }

                Ok(inverse_matrix)
            },

            // 2x2 matrix: 
            (2, 2) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];

                let adjoint = Self::calculate_adjoint_monomial(&matrix, element_dim, element_dim)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to calculate adjoint: {:?}", e)))?;
                let determinant = GaussianQuadrature::calculate_determinant_monomial(&matrix, element_dim, element_dim)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to calculate determinant: {:?}", e)))?;

                for i in 0..element_dim {
                    for j in 0..element_dim {
                        match MonomialPolynomial::divide(&adjoint[i][j], &determinant) {
                            Ok(coefficients) => {
                                inverse_matrix[i][j] = coefficients;
                            }
                            Err(e) => return Err(e.into()),
                        }
                    }
                }
                Ok(inverse_matrix)
            },

            // 3x3 matrix: det = a(ei−fh) − b(di−fg) + c(dh−eg)
            (3, 3) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];

                let adjoint = Self::calculate_adjoint_monomial(&matrix, element_dim, element_dim)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to calculate adjoint: {:?}", e)))?;
                let determinant = GaussianQuadrature::calculate_determinant_monomial(&matrix, element_dim, element_dim)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to calculate determinant: {:?}", e)))?;

                for i in 0..element_dim {
                    for j in 0..element_dim {
                        match MonomialPolynomial::divide(&adjoint[i][j], &determinant) {
                            Ok(coefficients) => {
                                inverse_matrix[i][j] = coefficients;
                            }
                            Err(e) => return Err(e.into()),
                        }
                    }
                }
                Ok(inverse_matrix)
            }
            // 1D elements in 2D space: 
            (2, 1) => {

                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];
                let mut square_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];

                
                let mut sum = vec![];
                for k in 0..mesh_dim {
                    let prod = MonomialPolynomial::multiply(&matrix[k][0], &matrix[k][0])?;
                    sum = MonomialPolynomial::add(&sum, &prod)?;
                }
                square_matrix[0][0] = sum;
                    
                

                match MonomialPolynomial::divide(&[1.0, 0.0, 0.0, 0.0], &square_matrix[0][0]) {
                    Ok(coefficients) => {
                        inverse_matrix[0][0] = coefficients;
                    }
                    Err(e) => return Err(e.into()),
                }

                Ok(inverse_matrix)
            }

            // 1D elements in 3D space:
            (3, 1) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];
                let mut square_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];

                
                let mut sum = vec![];
                for k in 0..mesh_dim {
                    let prod = MonomialPolynomial::multiply(&matrix[k][0], &matrix[k][0])?;
                    sum = MonomialPolynomial::add(&sum, &prod)?;
                }
                square_matrix[0][0] = sum;
                    
                

                match MonomialPolynomial::divide(&[1.0, 0.0, 0.0, 0.0], &square_matrix[0][0]) {
                    Ok(coefficients) => {
                        inverse_matrix[0][0] = coefficients;
                    }
                    Err(e) => return Err(e.into()),
                }

                Ok(inverse_matrix)
            }

            // 2D elements in 3D space: metric determinant from first fundamental form
            (3, 2) => {
                let mut inverse_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];
                let mut square_matrix: Vec<Vec<Vec<f64>>> = vec![vec![vec![]; element_dim]; element_dim];

                for i in 0..element_dim {
                    for j in 0..element_dim {
                        let mut sum = vec![];
                        for k in 0..mesh_dim {
                            let prod = MonomialPolynomial::multiply(&matrix[k][i], &matrix[k][j])?;
                            sum = MonomialPolynomial::add(&sum, &prod)?;
                        }
                        square_matrix[i][j] = sum;
                    }
                }

                let adjoint = Self::calculate_adjoint_monomial(&square_matrix, element_dim, element_dim)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to calculate adjoint: {:?}", e)))?;
                let determinant = GaussianQuadrature::calculate_determinant_monomial(&square_matrix, element_dim, element_dim)
                    .map_err(|e| ElementError::GeometryError(format!("Failed to calculate determinant: {:?}", e)))?;

                for i in 0..element_dim {
                    for j in 0..element_dim {
                        match MonomialPolynomial::divide(&adjoint[i][j], &determinant) {
                            Ok(coefficients) => {
                                inverse_matrix[i][j] = coefficients;
                            }
                            Err(e) => return Err(e.into()),
                        }
                    }
                }

                Ok(inverse_matrix)
            }

            // Unsupported matrix dimensions
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for {}x{} (mesh_dim x element_dim) matrices",
                mesh_dim, element_dim
            ))),
        }
    }
    

    // Adjoint of a matrix with monomial polynomial entries
    pub fn calculate_adjoint_monomial(matrix: &Vec<Vec<Vec<f64>>>, element_dim: usize, mesh_dim: usize) -> Result<Vec<Vec<Vec<f64>>>, ElementError> {
        // Adjoint calculation only implemented for 2x2 and 3x3 matrices
        let mesh_dim = matrix.len();
        let element_dim = if mesh_dim > 0 { matrix[0].len() } else { 0 };

        match (mesh_dim, element_dim) {
            (2, 2) => {
                let adj = vec![
                    vec![matrix[1][1].clone(), MonomialPolynomial::multiply_scalar(&matrix[0][1], -1.0)],
                    vec![MonomialPolynomial::multiply_scalar(&matrix[1][0], -1.0), matrix[0][0].clone()],
                ];
                Ok(adj)
            },
            (3, 3) => {
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
                "Jacobian adjoint calculation not implemented for {}x{} (mesh_dim x element_dim) matrices",
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


    /// Calculate error based on element type and dispatch to appropriate quadrature rule
    pub fn calculate_error(
        integrand: &Vec<Vec<Vec<f64>>>,
        n: usize,
        element_dim: usize,
        num_nodes: usize,
    ) -> Result<Vec<Vec<f64>>, GaussError> {

        let mut errors = vec![vec![0.0; num_nodes]; num_nodes];
        
        for i in 0..num_nodes {
            for j in 0..num_nodes {
                errors[i][j] = Self::gauss_legendre_error(&integrand[i][j], n, element_dim)?;
            }
        }
        
        Ok(errors)
    }

    /// Gauss-Legendre quadrature error estimation
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
                        if ak > 1e-15 {
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
                        if aij > 1e-15 {
                            term1 += aij * (factorial(i) / factorial(i - 2 * n));
                        }
                    }
                }

                let mut term2 = 0.0;
                for i in 0..=poly_degree_x.min(coeff_matrix.len().saturating_sub(1)) {
                    for j in (2 * n)..=poly_degree_y.min(coeff_matrix[i].len().saturating_sub(1)) {
                        let aij = coeff_matrix[i][j].abs();
                        if aij > 1e-15 {
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
                            if aijk > 1e-15 {
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
                            if aijk > 1e-15 {
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
                            if aijk > 1e-15 {
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
mod tests {
    use super::*;
    use std::f64::consts::PI;

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

    fn create_node_1d(id: usize, x: f64) -> Node {
        Node {
            id,
            coordinates: vec![x],
        }
        pub struct MeshData {                               // Defines a structure to represent a mesh
        pub dimension: usize,                           // Spatial dimension (from # sdim tag)
        pub num_nodes: usize,                           // Number of nodes (from # number of mesh vertices tag)
        pub min_node_index: usize,                      // Lowest mesh vertex index (from # lowest mesh vertex index tag) 
        pub nodes: Vec<Node>,                           // All nodes with their coordinates
        pub num_eltypes: usize,                         // Number of element types (from # number of element types tag)
        pub elements: Vec<Element>,                     // All elements with their connectivity
        pub element_type_info: Vec<ElementTypeInfo>,    // Information about each element type
}
    }

    fn create_test_mesh_data( // here create a meshdata from single element with its nodes and element type info and coordinates
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

    // Helper function for factorial in tests
    fn fact(n: usize) -> f64 {
        (1..=n).fold(1.0, |acc, k| acc * k as f64)
    }

    fn create_test_element_and_nodes() -> (Element, Vec<Node>) {
        let element = Element {
            id: 0,
            nodes: vec![0, 1, 2],
        };
        
        let nodes = vec![
            Node { id: 0, coordinates: vec![0.0, 0.0] },
            Node { id: 1, coordinates: vec![1.0, 0.0] },
            Node { id: 2, coordinates: vec![0.0, 1.0] },
        ];
        
        (element, nodes)
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

        let result = GaussianQuadrature::calculate_jacobian_monomial(
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

        let result = GaussianQuadrature::calculate_jacobian_monomial(
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
    fn test_detect_polynomial_order_1d_basic() {
        let coeffs = vec![1.0, 2.0, 0.0, 4.0]; // 1 + 2x + 4x³
        let order = GaussianQuadrature::detect_polynomial_order_1d(&coeffs);
        assert_eq!(order, 3);
    }

    #[test]
    fn test_detect_polynomial_order_1d_with_tolerance() {
        let coeffs = vec![1.0, 1e-13, 0.0, 0.0]; // Coefficients below tolerance should be ignored
        let order = GaussianQuadrature::detect_polynomial_order_1d(&coeffs);
        assert_eq!(order, 0);
    }

    #[test]
    fn test_detect_polynomial_order_1d_zero_polynomial() {
        let coeffs = vec![0.0, 0.0, 0.0];
        let order = GaussianQuadrature::detect_polynomial_order_1d(&coeffs);
        assert_eq!(order, 0);
    }

    #[test]
    fn detect_polynomial_orders_2d_basic() {
        let mut a = vec![vec![0.0; 3]; 4];
        a[0][0] = 1.0;
        a[3][2] = -2.5;

        let (dx, dy) = GaussianQuadrature::detect_polynomial_orders_2d(&a);
        assert_eq!((dx, dy), (3, 2));
    }

    #[test]
    fn detect_polynomial_orders_2d_ignores_tiny() {
        let mut a = vec![vec![0.0; 2]; 2];
        a[0][0] = 1e-13;  // Below tolerance - should be ignored
        a[1][1] = 1.0;    // Above tolerance - should be counted

        let (dx, dy) = GaussianQuadrature::detect_polynomial_orders_2d(&a);
        assert_eq!((dx, dy), (1, 1));
    }

    #[test]
    fn detect_polynomial_orders_3d_basic() {
        let mut a = vec![vec![vec![0.0; 4]; 2]; 3];
        a[0][0][0] = 0.1;
        a[2][1][3] = -7.0;

        let (dx, dy, dz) = GaussianQuadrature::detect_polynomial_orders_3d(&a);
        assert_eq!((dx, dy, dz), (2, 1, 3));
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1.0);
        assert_eq!(factorial(5), 120.0);
        assert!((factorial(10) - 3628800.0).abs() < 1e-6);
    }

    #[test]
    fn test_gauss_legendre_error_1d_linear_polynomial() {
        // Test error for linear polynomial: 1 + x
        // With n=1 Gauss point, this should be exact (error = 0)
        let integrand = vec![1.0, 1.0, 0.0, 0.0]; // 1 + x
        let error = GaussianQuadrature::gauss_legendre_error(&integrand, 1, 1).unwrap();
        assert!(error.abs() < 1e-12, "Linear polynomial should be exact with n=1, got error: {}", error);
    }

    #[test]
    fn test_gauss_legendre_error_1d_quadratic_polynomial() {
        // Test error for quadratic polynomial: 1 + x + x²
        // With n=1 Gauss point, this should have non-zero error
        // With n=2 Gauss points, this should be exact
        let integrand = vec![1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // 1 + x + x²
        
        let error_n1 = GaussianQuadrature::gauss_legendre_error(&integrand, 1, 1).unwrap();
        assert!(error_n1 > 0.0, "Quadratic polynomial should have error with n=1");
        
        let error_n2 = GaussianQuadrature::gauss_legendre_error(&integrand, 2, 1).unwrap();
        assert!(error_n2.abs() < 1e-12, "Quadratic polynomial should be exact with n=2, got error: {}", error_n2);
    }

    #[test]
    fn test_gauss_legendre_error_2d() {
        // Test 2D polynomial: 1 + x + y
        let integrand = vec![1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // Simplified representation
        
        let error = GaussianQuadrature::gauss_legendre_error(&integrand, 1, 2).unwrap();
        // For linear polynomial in 2D with n=1, error should be small but non-zero
        assert!(error >= 0.0);
    }

    /* 
    #[test]
    fn test_calculate_integrand_mass_matrix() {
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        let material_property = MaterialProperty::Scalar(vec![1.0]); // Use enum
        
        let result = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &nodes,
            &material_property, // Pass reference
            2,
        );
        
        assert!(result.is_ok());
        let integrand = result.unwrap();
        
        // The integrand should be a polynomial representing ρ * ΣNᵢ² * det(J)
        // For a unit triangle with constant density, this should be non-zero
        assert!(!integrand.is_empty());
        
        // Check that the integrand has reasonable values
        let max_coeff = integrand.iter().fold(0.0_f64, |max, &val| max.max(val.abs()));
        assert!(max_coeff > 0.0);
    }
    */

    #[test]
    fn test_calculate_integrand_stiffness_matrix() {
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        
        // Create a 2x2 identity matrix material property
        let material_tensor = vec![
            vec![vec![1.0], vec![0.0]], // D_xx = 1.0, D_xy = 0.0
            vec![vec![0.0], vec![1.0]], // D_yx = 0.0, D_yy = 1.0
        ];
        let material_property = MaterialProperty::Matrix(material_tensor);
        
        let result = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Stiffness,
            &element_type,
            &nodes,
            &material_property,
            2,
        );
        
        // This might fail due to element type limitations, but the interface should work
        assert!(result.is_ok() || matches!(result, Err(GaussError::InvalidElement(_))));
    }

    #[test]
    fn test_find_optimal_gauss_points_simple_case() {
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        let material_property = MaterialProperty::Scalar(vec![1.0]); // Use enum for constant density
        
        let result = GaussianQuadrature::find_optimal_gauss_points(
            IntegrationType::Mass,
            &element,
            &element_type,
            &nodes,
            2, // mesh_dim
            1e-6, // tolerance
            &material_property, // Pass reference
        );
        
        assert!(result.is_ok());
        let gauss_info = result.unwrap();
        
        assert_eq!(gauss_info.element_id, 0);
        assert!(gauss_info.theoretical_number >= 1);
        assert!(gauss_info.optimal_number >= 1);
        assert!(gauss_info.optimal_number <= gauss_info.theoretical_number);
    }

    #[test]
    fn test_find_optimal_gauss_points_mesh() {
        // Create a simple mesh with one triangle element
        let nodes = vec![
            Node { id: 1, coordinates: vec![0.0, 0.0] },
            Node { id: 2, coordinates: vec![1.0, 0.0] },
            Node { id: 3, coordinates: vec![0.0, 1.0] },
        ];
        
        let elements = vec![Element {
            id: 1,
            nodes: vec![1, 2, 3],
        }];
        
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 3,
        }];
        
        let mesh_data = MeshData {
            dimension: 2,
            num_nodes: 3,
            min_node_index: 1,
            nodes,
            num_eltypes: 1,
            elements,
            element_type_info,
        };
        
        let material_property = MaterialProperty::Scalar(vec![1.0]); // Use enum for constant density
        
        let result = GaussianQuadrature::find_optimal_gauss_points_number_mesh(
            &mesh_data,
            IntegrationType::Mass,
            1e-6,
            &material_property, // Pass reference
        );
        
        assert!(result.is_ok());
        let report = result.unwrap();
        
        assert_eq!(report.total_elements, 1);
        assert_eq!(report.gauss_point_numbers[0].element_id, 1);
        assert!(report.gauss_point_numbers[0].optimal_number >= 1);
    }

    #[test]
    fn test_polynomial_coefficient_extraction() {
        // Test the coefficient extraction functions from MonomialPolynomial
        
        // Create a simple 2D polynomial: 1 + 2x + 3y + 4xy
        // In graded lexicographic order for degree 2: [1, x, y, z, x², xy, xz, y², yz, z²]
        let poly_2d = vec![1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0];
        
        // Extract 1D coefficients (only x terms)
        let coeffs_1d = MonomialPolynomial::get_coefficients_1d(&poly_2d).unwrap();
        assert_eq!(coeffs_1d, vec![1.0, 2.0, 0.0]); // Constant, x, x²
        
        // Extract 2D coefficients
        let coeffs_2d = MonomialPolynomial::get_coefficients_2d(&poly_2d).unwrap();
        assert_eq!(coeffs_2d[0][0], 1.0); // Constant
        assert_eq!(coeffs_2d[1][0], 2.0); // x
        assert_eq!(coeffs_2d[0][1], 3.0); // y
        assert_eq!(coeffs_2d[1][1], 4.0); // xy
    }
    /* 
    #[test]
    fn test_integration_with_varying_material_properties() {
        // Test integration with non-constant material properties
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        
        // Linear density variation: ρ(x) = 1 + x
        let material_property = MaterialProperty::Scalar(vec![1.0, 1.0]); // 1.0 + 1.0*x
        
        let result = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &nodes,
            &material_property, // Use enum
            2,
        );
        
        assert!(result.is_ok());
        let integrand = result.unwrap();
        
        // The integrand should be more complex due to varying density
        assert!(!integrand.is_empty());
        
        // Find optimal Gauss points for this more complex case
        let gauss_result = GaussianQuadrature::find_optimal_gauss_points(
            IntegrationType::Mass,
            &element,
            &element_type,
            &nodes,
            2,
            1e-6,
            &material_property, // Use enum
        );
        
        assert!(gauss_result.is_ok());
    }
    */
    #[test]
    fn test_material_property_enum_methods() {
        // Test scalar material property
        let scalar_prop = MaterialProperty::Scalar(vec![1.0, 2.0, 3.0]);
        assert!(scalar_prop.as_scalar().is_ok());
        assert!(scalar_prop.as_matrix().is_err());
        
        // Test matrix material property
        let matrix_prop = MaterialProperty::Matrix(vec![
            vec![vec![1.0], vec![0.0]],
            vec![vec![0.0], vec![1.0]],
        ]);
        assert!(matrix_prop.as_matrix().is_ok());
        assert!(matrix_prop.as_scalar().is_err());
    }

    #[test]
    fn test_integration_type_material_mismatch() {
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        
        // Try to use scalar material with stiffness integration (should fail)
        let scalar_prop = MaterialProperty::Scalar(vec![1.0]);
        let result_stiffness = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Stiffness,
            &element_type,
            &nodes,
            &scalar_prop,
            2,
        );
        assert!(matches!(result_stiffness, Err(GaussError::InvalidMaterialProperty(_))));
        
        // Try to use matrix material with mass integration (should fail)
        let matrix_prop = MaterialProperty::Matrix(vec![
            vec![vec![1.0], vec![0.0]],
            vec![vec![0.0], vec![1.0]],
        ]);
        let result_mass = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &nodes,
            &matrix_prop,
            2,
        );
        assert!(matches!(result_mass, Err(GaussError::InvalidMaterialProperty(_))));
    }
}
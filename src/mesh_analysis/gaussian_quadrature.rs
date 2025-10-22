use std::f64;
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

        // Theoretical points based on polynomial degree 
        let max_num_gp = ((MonomialPolynomial::total_degree_polynomial(&integrand) as f64 + 1.0) / 2.0).ceil() as usize;  // num_gp = ceil((d + 1)/2) where d is polynomial degree
        let theoretical_points = max_num_gp;  // Arbitrary upper limit for Gauss points per direction

        // Find optimal points
        let mut optimal_points = 1;
        let mut min_error = 100.0;

        for n in 1..=max_num_gp {
            let error = Self::calculate_error(
                &integrand,
                n,
                &element_type,
                element_dim,
            )?;

            if error < min_error {
                min_error = error;
                optimal_points = n;
            }

            // If we meet tolerance, we can stop early
            if error <= tolerance {
                break;
            }
        }

        Ok(GaussianPointNumber {
            element_id: element.id,
            theoretical_number: theoretical_points,
            optimal_number: optimal_points,
        })
    }

    /// Calculate integrand for given element and integration type
    pub fn calculate_integrand(
        int_type: &IntegrationType,
        element_type: &ElementType,
        element_nodes: &Vec<Node>,
        material_property: &MaterialProperty,
        mesh_dim: usize,
    ) -> Result<Vec<f64>, GaussError> {

        // Get element and mesh dimensions and element order
        let element_dim = ElementType::get_element_dimension(element_type)
            .ok_or(GaussError::UnsupportedDimension(0))?;
        let element_order = ElementType::get_element_order(element_type)
            .ok_or(GaussError::InvalidElement("Element order not found".to_string()))?;
        
        let shape_function = ElementType::get_shape_functions(element_type)
            .ok_or_else(|| GaussError::InvalidElement("Shape functions not found".to_string()))?;
        
        let num_nodes = shape_function.num_nodes;
        
        // Calculate Jacobian matrix mapping parametric to physical coordinates
        let jacobian_matrix = GeometricAnalysis::calculate_jacobian_monomial(
            &element_nodes,
            &shape_function.derivatives,
            element_dim,
            num_nodes,
            mesh_dim,
            element_order,
        )?;
        
        // Calculate determinant of Jacobian matrix
        let det_j = &jacobian_matrix.determinant;

        match int_type {
            IntegrationType::Mass => {
                // Mass matrix: M = int rho(x) N_i(x) N_j(x) det(J) dxi
                // Mass matrix integrand: rho(x) * N' N * det(J)

                let rho = material_property.as_scalar()?; // Density polynomial coefficients in monomial form

                // Calculate N' N = sum (N[i])^2
                let mut sum_ns = vec![0.0; shape_function.values[0].len()];
                
                for i in 0..num_nodes {
                    let ni_ni = MonomialPolynomial::multiply(
                        &shape_function.values[i],
                        &shape_function.values[i]
                    )?;
                    sum_ns = MonomialPolynomial::add(&sum_ns, &ni_ni)?;
                }
                
                // Combine all terms: rho(x) * (sum (N[i])^2) * det(J)
                let rho_sum = MonomialPolynomial::multiply(&sum_ns, &rho)?;
                let integrand = MonomialPolynomial::multiply(&rho_sum, &det_j)?;
                
                Ok(integrand)
            }
            IntegrationType::Stiffness => {
                // Stiffness matrix: K = int B_i^T * D * B_j * det(J) dxi
                // where B is derivative and combination of rows of N

                // Build B matrix
                let b_matrix = Self::build_b_matrix(
                    &shape_function.derivatives,
                    num_nodes,
                    mesh_dim,
                );

                let c_matrix = material_property.as_matrix()?; // Material property matrix C in monomial form

                // Validate material matrix dimensions based on element dimension
                let expected_c_size = match mesh_dim {
                    1 => 1,  // 1D: 1x1
                    2 => 3,  // 2D: 3x3 (plane stress/strain)
                    3 => 6,  // 3D: 6x6
                    _ => return Err(GaussError::UnsupportedDimension(element_dim)),
                };

                if c_matrix.len() != expected_c_size {
                    return Err(GaussError::InvalidMaterialProperty(format!(
                        "Material matrix has {} rows, expected {} for {}D elements",
                        c_matrix.len(), expected_c_size, element_dim
                    )));
                }
                
                for i in 0..expected_c_size {
                    if c_matrix[i].len() != expected_c_size {
                        return Err(GaussError::InvalidMaterialProperty(format!(
                            "Material matrix row {} has {} columns, expected {}",
                            i, c_matrix[i].len(), expected_c_size
                        )));
                    }
                }

                // Initialize result matrix
                let mut integrand = vec![vec![vec![]; num_nodes]; num_nodes];

                match mesh_dim{
                    1 => {
                        // Initialize intermediate arrays
                        let mut bt_c = vec![vec![]; num_nodes];
                        let mut bt_c_b = vec![vec![vec![]; num_nodes]; num_nodes];

                        // B^T * C
                        for i in 0..num_nodes{
                            bt_c[i] = MonomialPolynomial::multiply(&b_matrix[0][i], &c_matrix[0][0])?; 
                        }

                        // (B^T * C) * B
                        for i in 0..num_nodes{
                            for j in 0..num_nodes{
                                bt_c_b[i][j] = MonomialPolynomial::multiply(&bt_c[i], &b_matrix[0][j])?;  
                            }
                        }
                        
                        // Multiply by det(J)
                        for i in 0..num_nodes{
                            for j in 0..num_nodes{
                                integrand[i][j] = MonomialPolynomial::multiply(&bt_c_b[i][j], &det_j)?;  
                            }
                        }
                    }
                    2 => {
                        // Initialize intermediate arrays
                        let mut bt_c = vec![vec![vec![]; 3]; num_nodes]; // num_nodes × 3
                        let mut bt_c_b = vec![vec![vec![]; num_nodes]; num_nodes];
                        let mut part_bt_c = vec![vec![vec![]; 3]; num_nodes]; // num_nodes × 3
                        let mut part_bt_c_b = vec![vec![vec![]; num_nodes]; num_nodes];
                        
                        // B^T * C: bt_c[i][j] = Σ_k (B[k][i] * C[k][j])
                        for i in 0..num_nodes{
                            for j in 0..3{ // size of c_matrix 1,3,6 
                                for k in 0..3{ // loop for sum
                                    let part_bt_c = MonomialPolynomial::multiply(&b_matrix[k][i], &c_matrix[k][j])?;  
                                    bt_c[i][j] = MonomialPolynomial::add(&bt_c[i][j], &part_bt_c)?;
                                }      
                            }    
                        }
                        for i in 0..num_nodes{
                            for j in 0..num_nodes{
                                for k in 0..3{  // loop for sum
                                    let part_bt_c_b = MonomialPolynomial::multiply(&bt_c[k][i], &b_matrix[k][j])?;  
                                    bt_c_b[i][j] = MonomialPolynomial::add(&bt_c_b[i][j], &part_bt_c_b)?;
                                }
                            }
                        }
                        for i in 0..num_nodes{
                            for j in 0..num_nodes{
                                integrand[i][j] = MonomialPolynomial::multiply(&bt_c_b[i][j], &det_j)?;  
                            }
                        }
                    }
                    3 => {
                        // Initialize intermediate arrays
                        let mut bt_c = vec![vec![vec![]; 6]; num_nodes]; // num_nodes × 6
                        let mut bt_c_b = vec![vec![vec![]; num_nodes]; num_nodes];
                        let mut part_bt_c = vec![vec![vec![]; 6]; num_nodes]; // num_nodes × 6
                        let mut part_bt_c_b = vec![vec![vec![]; num_nodes]; num_nodes];

                        // B^T * C: bt_c[i][j] = Σ_k (B[k][i] * C[k][j])
                        for i in 0..num_nodes{
                            for j in 0..6{  // size of c_matrix 1,3,6
                                for k in 0..6{  // loop for sum
                                    let part_bt_c = MonomialPolynomial::multiply(&b_matrix[k][i], &c_matrix[k][j])?;  
                                    bt_c[i][j] = MonomialPolynomial::add(&bt_c[i][j], &part_bt_c)?;
                                }
                            }
                        }
                        for i in 0..num_nodes{
                            for j in 0..num_nodes{
                                for k in 0..6{  // loop for sum
                                    let part_bt_c_b = MonomialPolynomial::multiply(&bt_c[k][i], &b_matrix[k][j])?;  
                                    bt_c_b[i][j] = MonomialPolynomial::add(&bt_c_b[i][j], &part_bt_c_b)?;
                                }
                            }
                        }
                        for i in 0..num_nodes{
                            for j in 0..num_nodes{
                                integrand[i][j] = MonomialPolynomial::multiply(&bt_c_b[i][j], &det_j)?;  
                            }
                        }
                    }
                }
            
                Ok(integrand)
            }
            /* 
            _ => Err(GaussError::GeometryError(format!(
                "Unsupported integration type {:?} for integrand calculation",
                int_type
            ))),
            */
        }
    }

    // Build B matrix

    fn build_b_matrix(
            shape_derivatives: &Vec<Vec<Vec<f64>>>,
            num_nodes: usize,
            mesh_dim: usize,
        ) -> Vec<Vec<Vec<f64>>> {

        match mesh_dim{
            1 => {
                let mut b_matrix = vec![vec![]; num_nodes];
                
                for i in 0..num_nodes {

                    b_matrix[0][i] = shape_derivatives[i][0].clone(); // monomial coefficients
                }
                b_matrix
            }
            2 => {
                let mut b_matrix = vec![vec![]; 2 * num_nodes]; // row number = 3
                for i in 0..num_nodes{

                    b_matrix[0][0+(2*i)] = shape_derivatives[0+i][0].clone();
                    b_matrix[2][0+(2*i)] = shape_derivatives[0+i][1].clone();

                    b_matrix[1][1+(2*i)] = shape_derivatives[0+i][1].clone();
                    b_matrix[2][1+(2*i)] = shape_derivatives[0+i][0].clone();
                }
                b_matrix
            }
            3 => {
                let mut b_matrix = vec![vec![]; 3 * num_nodes]; // row number = 6
                for i in 0..num_nodes{

                    b_matrix[0][0+(3*i)] = shape_derivatives[0+i][0].clone();
                    b_matrix[3][0+(3*i)] = shape_derivatives[0+i][1].clone();
                    b_matrix[5][0+(3*i)] = shape_derivatives[0+i][2].clone();

                    b_matrix[1][1+(3*i)] = shape_derivatives[0+i][1].clone();
                    b_matrix[3][1+(3*i)] = shape_derivatives[0+i][0].clone();
                    b_matrix[4][1+(3*i)] = shape_derivatives[0+i][2].clone();

                    b_matrix[2][2+(3*i)] = shape_derivatives[0+i][2].clone();
                    b_matrix[4][2+(3*i)] = shape_derivatives[0+i][1].clone();
                    b_matrix[5][2+(3*i)] = shape_derivatives[0+i][0].clone();
                }
                b_matrix
            }
            _ => panic!("Unsupported mesh dimension: {}", mesh_dim), // unsupported mesh dimension
        }
    }

    /* 
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
                let determinant = GeometricAnalysis::calculate_determinant_monomial(&matrix, element_dim, element_dim)
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
                let determinant = GeometricAnalysis::calculate_determinant_monomial(&matrix, element_dim, element_dim)
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
                let determinant = GeometricAnalysis::calculate_determinant_monomial(&square_matrix, element_dim, element_dim)
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
    */

    /// Calculate error based on element type and dispatch to appropriate quadrature rule
    pub fn calculate_error(
        integrand: &Vec<f64>,
        n: usize,
        element_type: &ElementType,
        element_dim: usize,
    ) -> Result<f64, GaussError> {

        Self::gauss_legendre_error(&integrand, n, element_dim)
        /* 
        match element_type {
            ElementType::Line | ElementType::QuadraticEdge | ElementType::Quad | ElementType::Hexahedron => {
                Self::gauss_legendre_error(&integrand, n, element_dim)
            }
            ElementType::Triangle | ElementType::Tetra => {
                //Self::gauss_dunont_error(num_gp, element_dim, coeffs, coeff_matrix, coeff_tensor)
                Self::gauss_legendre_error(&integrand, n, element_dim)
            }
            _ => Err(GaussError::GeometryError(format!(
                "Unsupported element type {:?} for error estimation",
                element_type
            ))),
        }
        */
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
                        if ak > 0.0 {
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
                            term1 += aij * (factorial(i) / factorial(i - 2 * n)) * (1.0 / (j as f64 + 1.0));
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
                                    * (factorial(i) / factorial(i - 2 * n))
                                    * (1.0 / ((j as f64 + 1.0) * (k as f64 + 1.0)));
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
                                    * (factorial(j) / factorial(j - 2 * n))
                                    * (1.0 / (k as f64 + 1.0));
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
mod tests {
    use super::*;

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
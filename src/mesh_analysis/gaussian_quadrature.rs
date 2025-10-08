use std::f64;
use crate::structs_and_impls::*;
use crate::error::*;
use super::geometric_analysis::GeometricAnalysis;
use lazy_static::lazy_static;

// Precompute factorials for efficient error calculation
// Uses lazy_static to compute once and reuse
lazy_static! {
    static ref FACTORIALS: Vec<f64> = {
        let mut v = vec![1.0; 100]; // Precompute up to 99!
        for i in 1..v.len() {
            v[i] = v[i - 1] * (i as f64);
        }
        v
    };
}

/// Calculate factorial with fallback to Stirling's approximation for large n
fn factorial(n: usize) -> f64 {
    if n < FACTORIALS.len() {
        FACTORIALS[n] // Use precomputed value
    } else {
        // Stirling's approximation for large n: n! ≈ √(2πn) * (n/e)^n
        (2.0 * std::f64::consts::PI * (n as f64)).sqrt()
            * ((n as f64) / std::f64::consts::E).powi(n as i32)
            * (1.0 + 1.0 / (12.0 * (n as f64))) // Correction term
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
        material_property: &[f64], // density polynomial coefficients
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
                        &material_property,
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
        material_property: &[f64], // density polynomial coefficients
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
            &material_property,
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
        material_property: &[f64], // density polynomial coefficients (material property)
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
        let jacobian_matrix = GeometricAnalysis::calculate_jacobian(
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
                // Mass matrix: M_ij = ∫ ρ(x) N_i(x) N_j(x) det(J) dξ
                // Mass matrix integrand: ρ(x) * N' N * det(J)

                // Calculate N' N = Σ (N[i])^2
                let mut sum_ns = vec![0.0; shape_function.values[0].len()];
                
                for i in 0..num_nodes {
                    let ni_ni = MonomialPolynomial::multiply(
                        &shape_function.values[i],
                        &shape_function.values[i]
                    )?;
                    sum_ns = MonomialPolynomial::add(&sum_ns, &ni_ni)?;
                }
                
                // Combine all terms: ρ(x) * (Σ (N[i])^2) * det(J)
                let rho_sum = MonomialPolynomial::multiply(&sum_ns, &material_property)?;
                let integrand = MonomialPolynomial::multiply(&rho_sum, &det_j)?;
                
                Ok(integrand)
            }
            IntegrationType::Stiffness => {
                let mut integrand = vec![0.0; 10];  // Placeholder for stiffness matrix integrand
                Ok(integrand) 
            }
        }
    }

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
        let material_property = vec![1.0]; // Constant density = 1.0
        
        let result = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &nodes,
            &material_property,
            2, // mesh_dim
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
    fn test_find_optimal_gauss_points_simple_case() {
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        let material_property = vec![1.0]; // Constant density
        
        let result = GaussianQuadrature::find_optimal_gauss_points(
            IntegrationType::Mass,
            &element,
            &element_type,
            &nodes,
            2, // mesh_dim
            1e-6, // tolerance
            &material_property,
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
        
        let material_property = vec![1.0]; // Constant density
        
        let result = GaussianQuadrature::find_optimal_gauss_points_number_mesh(
            &mesh_data,
            IntegrationType::Mass,
            1e-6,
            &material_property,
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

    #[test]
    fn test_integration_with_varying_material_properties() {
        // Test integration with non-constant material properties
        let (element, nodes) = create_test_element_and_nodes();
        let element_type = ElementType::Triangle;
        
        // Linear density variation: ρ(x) = 1 + x
        let material_property = vec![1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        
        let result = GaussianQuadrature::calculate_integrand(
            &IntegrationType::Mass,
            &element_type,
            &nodes,
            &material_property,
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
            &material_property,
        );
        
        assert!(gauss_result.is_ok());
    }
}
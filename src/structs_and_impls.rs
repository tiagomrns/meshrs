use vtkio::model::CellType;  
use crate::error::*;                     // Import mesh data structures and error types from error module


#[derive(Debug, Clone)]             // Auto-implement Debug for printing and Clone for copying 
                                    // Clone is used to copy in the case of deformed mesh. (for undeformed mesh we dont need to clone because we maintain the originaal values)
pub struct Node {                   // Defines a structure to represent a mesh node/vertex
    pub id: usize,                  // Unique identifier for the node
    pub coordinates: Vec<f64>,      // Spatial coordinates (x, y for 2D, x,y,z for 3D)
}


#[derive(Debug, Clone)]
pub struct Element {                // Defines a structure to represent a mesh element
    pub id: usize,                  // Unique identifier for the element
    pub nodes: Vec<usize>,          // List of node IDs forming this element
}

/// Represents polynomials in 3D using graded lexicographic ordering
/// This is a compact representation for multivariate polynomials
#[derive(Debug, Clone)]
pub struct MonomialPolynomial {
    pub coefficients: Vec<f64>, // Coefficients in graded lexicographic order for 3D
    // Order: [1, x, y, z, x², xy, xz, y², yz, z², x³, x²y, x²z, xy², xyz, xz², y³, y²z, yz², z³, ...]
}


impl MonomialPolynomial {
    /// Create a new polynomial from coefficient vector
    /// Automatically trims to a valid length matching the 3D polynomial structure
    pub fn new(coefficients: Vec<f64>) -> Self {
        // Always create a valid polynomial, even if trimming fails
        let trimmed = Self::trim_to_valid_length_fallback(&coefficients);
        Self { coefficients: trimmed }
    }

    /// Infer maximum polynomial degree from coefficient vector length
    /// For 3D polynomials, uses formula: length = (max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6
    /// This is the number of monomials in a complete 3D polynomial of the given degree
    pub fn infer_max_degree(len: usize) -> Result<u32, PolynomialError> {
        for max_degree in 0..=100 {
            let expected = ((max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6) as usize;
            if expected == len {
                return Ok(max_degree);
            }
            if expected > len {
                break;
            }
        }
        Err(PolynomialError::InvalidCoefficientLength(len))
    }

    /// Calculate expected coefficient vector length for a given maximum degree in 3D
    /// Returns the number of monomials: (d+1)(d+2)(d+3)/6 for degree d
    pub fn expected_length(max_degree: u32) -> usize {
        ((max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6) as usize
    }

    /// Generate graded lexicographic basis for 3D with given maximum degree
    /// Returns a vector of (i, j, k) tuples representing x^i * y^j * z^k terms
    /// Order: [1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, ...]
    /// Sorted first by total degree (i+j+k), then by x (descending), then by y (descending)
    pub fn generate_basis(max_degree: u32) -> Vec<(u32, u32, u32)> {
        let mut basis = Vec::new();
        
        // Iterate through all possible total degrees (0 to max_degree)
        for total_degree in 0..=max_degree {
            // Generate all combinations (i, j, k) such that i + j + k = total_degree
            // Using reverse iteration to get lexicographic ordering
            for i in (0..=total_degree).rev() {
                for j in (0..=(total_degree - i)).rev() {
                    let k = total_degree - i - j;  // k is determined by constraint i+j+k=total_degree
                    basis.push((i, j, k));
                }
            }
        }   
        
        basis
    }

    // Map 3D exponent tuple (i,j,k) to its index in graded lexicographic ordering
    // This uses a combinatorial formula based on:
    // 1. Count all monomials with total degree < n (where n = i+j+k)
    // 2. Within degree n, count monomials that come before (i,j,k) in reverse lex order
    fn map_index(exponent: (u32, u32, u32)) -> usize {
        let (i, j, k) = exponent;
        let n = i + j + k;  // Total degree of the monomial
        
        // Count all monomials with total degree < n
        // This uses the formula for number of non-negative integer solutions to i+j+k<n
        // which equals (n)(n+1)(n+2)/6
        let index_before_degree_n = if n > 0 {
            (n as usize) * (n as usize + 1) * (n as usize + 2) / 6
        } else {
            0
        };
        
        // Within degree n, monomials are ordered in reverse lexicographic order:
        // x from n down to 0, then for each x, y from (n-x) down to 0
        
        // Count monomials with x > i (these come before our monomial)
        let mut count = 0;
        for x in (i + 1)..=n {
            // For this x value, number of possible y values = n - x + 1
            count += (n - x + 1) as usize;
        }
        
        // For x = i, count monomials with y > j (these also come before our monomial)
        count += (n - i - j) as usize;
        
        index_before_degree_n + count
    }

    /// Find the index of an exponent tuple in the basis
    /// Returns the position of (i,j,k) in the graded lexicographic ordering
    pub fn find_basis_index(basis: &[(u32, u32, u32)], exponent: (u32, u32, u32)) -> Result<usize, &'static str> {
        // Try fast O(1) lookup using the combinatorial formula
        let index = MonomialPolynomial::map_index(exponent);
        if index < basis.len() && basis[index] == exponent {
            Ok(index)
        } else {
            // Fall back to linear search if formula fails (should be rare)
            basis.iter()
                .position(|&exp| exp == exponent)
                .ok_or("Exponent not found in basis")
        }
    }

    /// Fallback version of trim_to_valid_length that never fails - used in constructor
    /// Ensures polynomial always has a valid length for 3D graded lexicographic ordering
    pub fn trim_to_valid_length_fallback(coeffs: &[f64]) -> Vec<f64> {
        // Handle empty input immediately
        if coeffs.is_empty() {
            return vec![0.0];
        }
        
        // Try the normal trimming first
        match Self::trim_to_valid_length(coeffs) {
            Ok(trimmed) => trimmed,
            Err(_) => {
                // If trimming fails, use a simple fallback
                // Find the last non-zero coefficient and return up to that point
                let mut last_non_zero = 0;
                for (i, &coeff) in coeffs.iter().enumerate() {
                    if coeff.abs() > 1e-12 {
                        last_non_zero = i;
                    }
                }
                
                // Return a constant polynomial if everything is zero or we couldn't find anything
                if last_non_zero == 0 {
                    if coeffs[0].abs() < 1e-12 {
                        vec![0.0]  // All zeros case
                    } else {
                        vec![coeffs[0]]  // Just constant term
                    }
                } else {
                    coeffs[..=last_non_zero].to_vec()
                }
            }
        }
    }

    /// Trim coefficient vector to a valid length for 3D polynomial representation
    /// Removes trailing zeros and adjusts length to match the expected size for some degree
    pub fn trim_to_valid_length(coeffs: &[f64]) -> Result<Vec<f64>, PolynomialError> {
        if coeffs.is_empty() {
            return Ok(vec![0.0]);
        }
        
        match Self::infer_max_degree(coeffs.len()) {
            Ok(max_degree) => {
                let basis = Self::generate_basis(max_degree);
                let mut actual_degree = 0;
                
                for (idx, &(i, j, k)) in basis.iter().enumerate() {
                    if idx < coeffs.len() && coeffs[idx].abs() > 1e-12 {
                        actual_degree = actual_degree.max(i + j + k);
                    }
                }
                
                let valid_length = Self::expected_length(actual_degree);
                let mut result = vec![0.0; valid_length];
                let copy_len = valid_length.min(coeffs.len());
                result[..copy_len].copy_from_slice(&coeffs[..copy_len]);
                Ok(result)
            },
            Err(_) => {
                // Fallback: find last non-zero coefficient and create valid polynomial
                let mut last_non_zero = 0;
                for (i, &coeff) in coeffs.iter().enumerate() {
                    if coeff.abs() > 1e-12 {
                        last_non_zero = i;
                    }
                }
                
                // Find the minimum degree that can accommodate the last non-zero coefficient index
                // We search for the smallest degree where expected_length(degree) > last_non_zero
                let mut degree = 0;
                while Self::expected_length(degree) <= last_non_zero + 1 {
                    degree += 1;
                    // Safety limit: polynomials above degree 20 are rarely needed in FEM applications
                    if degree > 20 {
                        return Err(PolynomialError::InvalidCoefficientLength(coeffs.len()));
                    }
                }
                
                let valid_length = Self::expected_length(degree);
                let mut result = vec![0.0; valid_length];
                let copy_len = valid_length.min(coeffs.len());
                result[..copy_len].copy_from_slice(&coeffs[..copy_len]);
                Ok(result)
            }
        }
    }

    /// Optimized multiplication with intelligent constant/zero detection
    /// This optimization significantly speeds up polynomial operations in FEM applications
    /// by handling common special cases before resorting to full polynomial multiplication
    pub fn multiply(first: &[f64], second: &[f64]) -> Result<Vec<f64>, PolynomialError> {
        // Early zero detection - if either polynomial is zero, return zero
        if Self::is_zero_polynomial(first) || Self::is_zero_polynomial(second) {
            return Ok(vec![0.0]);
        }

        // Constant detection and optimization - handles three cases:
        // Case 1: Both polynomials are constants (e.g., [5.0] * [3.0] = [15.0])
        // Case 2: First is constant, second is polynomial (e.g., [2.0] * [1, x, y] = [2, 2x, 2y])
        // Case 3: Second is constant, first is polynomial (symmetric to case 2)
        if let Some(first_const) = Self::get_constant_value(first) {
            if let Some(second_const) = Self::get_constant_value(second) {
                // Both are constants - return their product as constant
                return Ok(vec![first_const * second_const]);
            } else {
                // First is constant, second is not - use scalar multiplication
                return Ok(Self::multiply_scalar(second, first_const));
            }
        } else if let Some(second_const) = Self::get_constant_value(second) {
            // Second is constant, first is not - use scalar multiplication
            return Ok(Self::multiply_scalar(first, second_const));
        }

        // Identity optimization - multiplying by polynomial 1 returns the other polynomial
        if first.len() == 1 && first[0] == 1.0 {
            return Ok(second.to_vec());
        }
        if second.len() == 1 && second[0] == 1.0 {
            return Ok(first.to_vec());
        }

        // Fall back to standard multiplication for non-constant, non-identity polynomials
        // This performs the full O(n*m) multiplication using basis expansion
        let max_deg_first = Self::infer_max_degree(first.len())?;
        let max_deg_second = Self::infer_max_degree(second.len())?;

        let basis_first = Self::generate_basis(max_deg_first);
        let basis_second = Self::generate_basis(max_deg_second);

        // Result degree is sum of input degrees (e.g., degree 2 * degree 3 = degree 5)
        let result_max_degree = max_deg_first + max_deg_second;
        let result_len = Self::expected_length(result_max_degree);
        
        let mut result = vec![0.0; result_len];
        
        // Multiply each term in first by each term in second
        // (c1 * x^i * y^j * z^k) * (c2 * x^i' * y^j' * z^k') = (c1*c2) * x^(i+i') * y^(j+j') * z^(k+k')
        for (i, &(exp_i_x, exp_i_y, exp_i_z)) in basis_first.iter().enumerate() {
            for (j, &(exp_j_x, exp_j_y, exp_j_z)) in basis_second.iter().enumerate() {
                let coeff_product = first[i] * second[j];
                
                // Skip negligible terms to improve performance
                if coeff_product.abs() < 1e-12 {
                    continue;
                }
                
                // Add exponents to get result monomial: x^i * x^j = x^(i+j)
                let result_exponent = (
                    exp_i_x + exp_j_x,
                    exp_i_y + exp_j_y,
                    exp_i_z + exp_j_z,
                );
                
                let result_index = Self::map_index(result_exponent);
                if result_index < result_len {
                    result[result_index] += coeff_product;
                }
            }
        }
        
        Self::trim_to_valid_length(&result)
    }

    /// Optimized addition with intelligent constant/zero detection
    /// Handles polynomial addition: p1(x,y,z) + p2(x,y,z)
    pub fn add(first: &[f64], second: &[f64]) -> Result<Vec<f64>, PolynomialError> {
        // Early zero detection - adding zero doesn't change the polynomial
        if Self::is_zero_polynomial(first) {
            return Ok(second.to_vec());
        }
        if Self::is_zero_polynomial(second) {
            return Ok(first.to_vec());
        }

        // Constant detection and optimization
        if let Some(first_const) = Self::get_constant_value(first) {
            if let Some(second_const) = Self::get_constant_value(second) {
                // Both are constants - return their sum as constant
                return Ok(vec![first_const + second_const]);
            } else {
                // First is constant, second is polynomial - add constant to polynomial's constant term
                let mut result = second.to_vec();
                if !result.is_empty() {
                    result[0] += first_const;
                }
                return Self::trim_to_valid_length(&result);
            }
        } else if let Some(second_const) = Self::get_constant_value(second) {
            // Second is constant, first is polynomial - add constant to polynomial's constant term
            let mut result = first.to_vec();
            if !result.is_empty() {
                result[0] += second_const;
            }
            return Self::trim_to_valid_length(&result);
        }

        // Standard polynomial addition for non-constant polynomials
        // Simply add coefficients term-by-term
        let max_len = first.len().max(second.len());
        let result: Vec<f64> = (0..max_len)
            .map(|i| {
                first.get(i).unwrap_or(&0.0) + second.get(i).unwrap_or(&0.0)
            })
            .collect();
        
        Self::trim_to_valid_length(&result)
    }

    /// Multiply polynomial by scalar with optimization
    /// Computes: scalar * p(x,y,z) = scalar * Σ(coeff[i] * monomial[i])
    pub fn multiply_scalar(coeffs: &[f64], scalar: f64) -> Vec<f64> {
        // Early detection for zero scalar or zero polynomial
        if scalar.abs() < 1e-12 || Self::is_zero_polynomial(coeffs) {
            return vec![0.0];
        }
        
        // Identity optimization - multiplying by 1 returns the original polynomial
        if scalar == 1.0 {
            return coeffs.to_vec();
        }

        // For constant polynomial, just multiply the constant
        if let Some(const_val) = Self::get_constant_value(coeffs) {
            return vec![const_val * scalar];
        }

        // General case: multiply each coefficient by the scalar
        coeffs.iter()
            .map(|c| c * scalar)
            .collect()
    }

    // Evaluate polynomial at a point (x, y, z)
    // Computes: sum of coeff[idx] * x^i * y^j * z^k for each basis element (i,j,k)
    pub fn evaluate(coeffs: &[f64], point: (f64, f64, f64)) -> Result<f64, PolynomialError> {
        let max_degree = Self::infer_max_degree(coeffs.len())?;
        let (x, y, z) = point;
        let basis = Self::generate_basis(max_degree);

        // Handle case where coeffs is shorter than basis
        let evaluation_length = coeffs.len().min(basis.len());

        // Evaluate polynomial: Σ coeff[idx] * x^i * y^j * z^k
        Ok(coeffs[..evaluation_length].iter()
            .zip(&basis[..evaluation_length])
            .map(|(coeff, &(i, j, k))| {
                coeff * x.powi(i as i32) * y.powi(j as i32) * z.powi(k as i32)
            })
            .sum())
    }

    // Safe evaluation that never fails
    pub fn evaluate_safe(coeffs: &[f64], point: (f64, f64, f64)) -> f64 {
        Self::evaluate(coeffs, point).unwrap_or(0.0)
    }

    // Safe multiplication that never fails
    pub fn multiply_optimized(first: &[f64], second: &[f64]) -> Vec<f64> {
        Self::multiply(first, second).unwrap_or_else(|_| {
            // Fallback: return a constant zero polynomial
            vec![0.0]
        })
    }

    // Safe addition that never fails  
    pub fn add_optimized(first: &[f64], second: &[f64]) -> Vec<f64> {
        Self::add(first, second).unwrap_or_else(|_| {
            // Fallback: return the longer polynomial
            if first.len() >= second.len() {
                first.to_vec()
            } else {
                second.to_vec()
            }
        })
    }

    // ========== OPTIMIZED OPERATIONS WITH CONSTANT/ZERO DETECTION ==========

    /// Check if polynomial is effectively zero (all coefficients near zero)
    /// Uses tolerance of 1e-12 for floating-point comparison
    pub fn is_zero_polynomial(poly: &[f64]) -> bool {
        poly.iter().all(|&c| c.abs() < 1e-12)
    }
    
    /// Check if polynomial is constant (only first coefficient non-zero)
    /// A constant polynomial has the form: f(x,y,z) = c where c is a scalar
    pub fn is_constant_polynomial(poly: &[f64]) -> bool {
        if poly.is_empty() { 
            return true; 
        }
        
        // For trimmed polynomials, constant means length == 1
        if poly.len() == 1 {
            return true;
        }
        
        // For untrimmed polynomials, check if all coefficients beyond first are zero
        poly.iter().enumerate().all(|(i, &c)| i == 0 || c.abs() < 1e-12)
    }
    
    /// Get constant value from polynomial if it's constant
    /// Returns Some(value) if polynomial is constant, None otherwise
    pub fn get_constant_value(poly: &[f64]) -> Option<f64> {
        if Self::is_constant_polynomial(poly) && !poly.is_empty() {
            Some(poly[0])
        } else {
            None
        }
    }

    /// Get 1D coefficients (only x powers): [a₀, a₁, a₂, ...] for a₀ + a₁x + a₂x² + ...
    /// Extracts coefficients where y and z exponents are zero
    pub fn get_coefficients_1d(coeffs: &[f64]) -> Result<Vec<f64>, &'static str> {
        // Handle empty input
        if coeffs.is_empty() {
            return Ok(vec![0.0]);
        }
        
        let degree = match Self::infer_max_degree(coeffs.len()) {
            Ok(d) => d,
            Err(_) => return Ok(vec![0.0]), // Return zero polynomial on error
        };
        let mut result = vec![0.0; (degree + 1) as usize];
        let basis = Self::generate_basis(degree);
        
        for (idx, &(i, j, k)) in basis.iter().enumerate() {
            if j == 0 && k == 0 && idx < coeffs.len() {
                result[i as usize] = coeffs[idx];
            }
        }
        Ok(result)
    }

    /// Get 2D coefficients as matrix: result[i][j] = coefficient of x^i * y^j
    /// Extracts coefficients where z exponent is zero
    pub fn get_coefficients_2d(coeffs: &[f64]) -> Result<Vec<Vec<f64>>, &'static str> {
        // Handle empty input
        if coeffs.is_empty() {
            return Ok(vec![vec![0.0]]);
        }
        
        let degree = match Self::infer_max_degree(coeffs.len()) {
            Ok(d) => d,
            Err(_) => return Ok(vec![vec![0.0]]), // Return zero polynomial on error
        };
        let size = (degree + 1) as usize;
        let mut result = vec![vec![0.0; size]; size];
        let basis = Self::generate_basis(degree);
        
        for (idx, &(i, j, k)) in basis.iter().enumerate() {
            if k == 0 && idx < coeffs.len() {
                result[i as usize][j as usize] = coeffs[idx];
            }
        }
        Ok(result)
    }

    /// Get 3D coefficients as tensor: result[i][j][k] = coefficient of x^i * y^j * z^k
    /// Full 3D coefficient extraction
    pub fn get_coefficients_3d(coeffs: &[f64]) -> Result<Vec<Vec<Vec<f64>>>, &'static str> {
        // Handle empty input
        if coeffs.is_empty() {
            return Ok(vec![vec![vec![0.0]]]);
        }
        
        let degree = match Self::infer_max_degree(coeffs.len()) {
            Ok(d) => d,
            Err(_) => return Ok(vec![vec![vec![0.0]]]), // Return zero polynomial on error
        };
        let size = (degree + 1) as usize;
        let mut result = vec![vec![vec![0.0; size]; size]; size];
        let basis = Self::generate_basis(degree);
        
        for (idx, &(i, j, k)) in basis.iter().enumerate() {
            if idx < coeffs.len() {
                result[i as usize][j as usize][k as usize] = coeffs[idx];
            }
        }
        Ok(result)
    }

    /// Calculate the total degree with tolerance for floating point comparison
    pub fn total_degree_polynomial(polynomial: &[f64]) -> u32 {
        // Handle empty input
        if polynomial.is_empty() {
            return 0;
        }
        
        let max_degree = match Self::infer_max_degree(polynomial.len()) {
            Ok(deg) => deg,
            Err(_) => return 0, // Return 0 degree on error
        };
        
        let basis = Self::generate_basis(max_degree);
        let mut max_total_degree = 0;
        
        for (idx, &(i, j, k)) in basis.iter().enumerate() {
            if idx < polynomial.len() && polynomial[idx].abs() > 1e-12 {
                let total_degree = i + j + k;
                if total_degree > max_total_degree {
                    max_total_degree = total_degree;
                }
            }
        }
        
        max_total_degree
    }
}

// Convenience: Allow creating from Vec<f64> directly
impl From<Vec<f64>> for MonomialPolynomial {
    fn from(coefficients: Vec<f64>) -> Self {
        Self::new(coefficients)
    }
}

/// Shape functions and their derivatives for finite elements
/// Shape functions N(ξ,η,ζ) map natural coordinates to physical space
/// Used in finite element analysis for element interpolation
#[derive(Debug, Clone)]
pub struct ShapeFunction {
    pub values: Vec<Vec<f64>>,              // Shape function coefficients: values[node_i] = coefficients for N_i
    pub derivatives: Vec<Vec<Vec<f64>>>,    // Derivatives: derivatives[direction][node][coefficients]
                                            // direction: 0=∂/∂ξ, 1=∂/∂η, 2=∂/∂ζ (natural coordinates)
                                            // node: index of the shape function (0 to num_nodes-1)
    pub num_nodes: usize,                   // Number of nodes in the element
}

/// Element types available in VTK and supported FEM software
/// Covers 1D (lines), 2D (triangles, quads), and 3D (tetra, hex, pyramid, wedge) elements
/// Each type supports both linear (first-order) and higher-order (quadratic, biquadratic, triquadratic) variants
#[derive(Debug, Clone)]
pub enum ElementType {              // Enumeration of supported finite element types
    Vertex,                         // Point element (0D)   Comsol
    Line,                           // First order edge element (2 nodes)   Comsol
    QuadraticEdge,                  // Second order edge element (3 nodes)  Comsol
    Triangle,                       // First order triangular element (3 nodes)   Comsol
    QuadraticTriangle,              // Second order triangular element (6 nodes)   Comsol
    Quad,                           // First order quadrilateral element (4 nodes)   Comsol
    QuadraticQuad,                  // Second order quadrilateral element (8 nodes)  Abaqus
    BiquadraticQuad,                // Second order quadrilateral element (9 nodes with center)  Comsol
    Tetra,                          // First order tetrahedral element (4 nodes)   Comsol
    QuadraticTetra,                 // Second order tetrahedral element (10 nodes)   Comsol
    Pyramid,                        // First order pyramid element (5 nodes)   Comsol
    QuadraticPyramid,               // Second order pyramid element (13 nodes used from 14-node formulation)   Comsol
                                    // Note: Comsol uses 13 nodes, omitting the middle face node from full 14-node pyramid
    Wedge,                          // First order prism/wedge element (6 nodes)   Comsol
    QuadraticWedge,                 // Second order prism element (15 nodes)   Abaqus
    BiquadraticQuadraticWedge,      // Second order prism element   Comsol
    Hexahedron,                     // First order hexahedral element   Comsol
    QuadraticHexahedron,            // Second order hexahedral element  20 nodes   Abaqus
    BiquadraticQuadraticHexahedron, // Second order hexahedral element  24 nodes   Abaqus
    TriquadraticHexahedron,         // Second order hexahedral element  27 nodes   Comsol                
}

impl ElementType {
    
    /* 
    pub fn from_str_ansys(s: &str) -> Option<ElementType> { // Converts string from input data to its element type correspond - Ansys to xml vtk (probably no need to reorder nodes for xml vtk) 

        // No match found
        None
    }
    */

    pub fn eltype_vtk(&self) -> CellType { //converts element type from ElementType to VTK element type ID
        
        match self {
            ElementType::Vertex => CellType::Vertex,
            ElementType::Line => CellType::Line,  
            ElementType::QuadraticEdge => CellType::QuadraticEdge,  
            ElementType::Triangle => CellType::Triangle, 
            ElementType::QuadraticTriangle => CellType::QuadraticTriangle,  
            ElementType::Quad => CellType::Quad,
            ElementType::QuadraticQuad => CellType::QuadraticQuad, 
            ElementType::BiquadraticQuad => CellType::BiquadraticQuad, 
            ElementType::Tetra => CellType::Tetra,  
            ElementType::QuadraticTetra => CellType::QuadraticTetra,  
            ElementType::Pyramid => CellType::Pyramid,
            ElementType::QuadraticPyramid => CellType::QuadraticPyramid,   
            ElementType::Wedge => CellType::Wedge,   
            ElementType::QuadraticWedge => CellType::QuadraticWedge,
            ElementType::BiquadraticQuadraticWedge => CellType::BiquadraticQuadraticWedge,  
            ElementType::Hexahedron => CellType::Hexahedron,   
            ElementType::QuadraticHexahedron => CellType::QuadraticHexahedron,  
            ElementType::BiquadraticQuadraticHexahedron => CellType::BiquadraticQuadraticHexahedron, 
            ElementType::TriquadraticHexahedron => CellType::TriquadraticHexahedron,  
            
        }
    }

    /// Monomial shape functions and their derivatives for standard elements N and dN/dxi
    /// Monomial rule for Graded Lexicographic Order:
    pub fn get_shape_functions(
        element_type: &ElementType
    ) -> Option<ShapeFunction> {                  

        match element_type {
            ElementType::Vertex => {
                // Not implemented for Vertex element
                None
            },

            ElementType::Line => {
                let num_nodes = 2;
                // 1D: N0 = 1 - xi, N1 = xi
                // Padding to 3D degree 1: [1, xi, eta, psi]
                
                let values = vec![
                    vec![1.0, -1.0, 0.0, 0.0],  // N0 = 1 - xi
                    vec![0.0, 1.0, 0.0, 0.0],   // N1 = xi
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0],   // dN0/dxi = -1 , dN1/dxi = 1
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 0.0],   // dN0/deta = 0 , dN1/deta = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 0.0],  // dN0/dpsi = 0 , dN1/dpsi = 0
                    ],
                ]; 

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticEdge => {
                let num_nodes = 3;
                // 1D degree 2: N0 = 1 - 3xi + 2xi², N1 = -xi + 2xi², N2 = 4xi - 4xi²
                // Padding to 3D degree 2: [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]
                
                let values = vec![
                    vec![1.0, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N0 = 1 - 3xi + 2xi²
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N1 = -xi + 2xi²  
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N2 = 4xi - 4xi²
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dxi = -3 + 4xi 
                        vec![-1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dxi = -1 + 4xi
                        vec![4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dxi = 4 - 8xi 
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN0/deta = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/deta = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/deta = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dpsi = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Triangle => {
                let num_nodes = 3;
                // 2D degree 1: N0 = 1 - xi - eta, N1 = xi, N2 = eta
                // Padding to 3D degree 1: [1, xi, eta, psi]
                
                let values = vec![
                    vec![1.0, -1.0, -1.0, 0.0],  // N0 = 1 - xi - eta
                    vec![0.0, 1.0, 0.0, 0.0],    // N1 = xi
                    vec![0.0, 0.0, 1.0, 0.0],    // N2 = eta
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN0/dxi = -1
                        vec![1.0, 0.0, 0.0, 0.0],   // dN1/dxi = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dxi = 0
                    ],
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN0/deta = -1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/deta = 0
                        vec![1.0, 0.0, 0.0, 0.0],   // dN2/deta = 1
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0],   // dN0/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dpsi = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticTriangle => {
                let num_nodes = 6;
                // 2D degree 2, padding to 3D degree 2: [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]
                
                let values = vec![
                    vec![1.0, -3.0, -3.0, 0.0, 2.0, 4.0, 0.0, 2.0, 0.0, 0.0],  // N0
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N1
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0],   // N2
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0],  // N3
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0],    // N4
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0],  // N5
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dxi
                        vec![-1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dxi
                        vec![4.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dxi
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dxi
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN5/dxi
                    ],
                    vec![
                        vec![-3.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/deta
                        vec![-1.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/deta
                        vec![0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/deta
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/deta
                        vec![4.0, -4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN5/deta
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dpsi
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Quad => {
                let num_nodes = 4;
                // 2D degree 1 with bilinear term, padding to 3D degree 2
                // [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]

                let values = vec![
                    vec![1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  // N0 = (1-xi)(1-eta) = 1 - xi - eta + xi*eta
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],   // N1 = xi(1-eta) = xi - xi*eta
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],    // N2 = xi*eta
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],   // N3 = (1-xi)eta = eta - xi*eta
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dxi = -1 + eta
                        vec![1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dxi = 1 - eta
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dxi = eta
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dxi = -eta
                    ],
                    vec![
                        vec![-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/deta = -1 + xi
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/deta = -xi
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/deta = xi
                        vec![1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/deta = 1 - xi
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dpsi = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticQuad => {
                let num_nodes = 8;
                // 2D degree 2 serendipity element
                // Basis: [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]
                
                let values = vec![
                    vec![1.0, -3.0, -3.0, 0.0, 2.0, 5.0, 0.0, 2.0, 0.0, 0.0, 0.0, -2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // N0 = -2*eta^2*xi + 2*eta^2 - 2*eta*xi^2 + 5*eta*xi - 3*eta + 2*xi^2 - 3*xi + 1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // N1 = 2*eta^2*xi - 2*eta*xi^2 - eta*xi + 2*xi^2 - xi  
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          // N2 = 2*eta^2*xi + 2*eta*xi^2 - 3*eta*xi
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // N3 = -2*eta^2*xi + 2*eta^2 + 2*eta*xi^2 - eta*xi - eta
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // N4 = 4*xi - 4*eta*xi + 4*eta*xi^2 - 4*xi^2
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          // N5 = -4*xi*eta^2 + 4*xi*eta
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          // N6 = -4*eta*xi^2 + 4*eta*xi
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // N7 = 4*eta - 4*eta*xi + 4*eta^2*xi - 4*eta^2
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 5.0, 0.0, 0.0, -4.0, 0.0, -2.0, 0.0, 0.0],      // dN0/dxi = 5*eta + 4*xi - 4*eta*xi - 2*eta^2 - 3
                        vec![-1.0, 4.0, -1.0, 0.0, 0.0, -4.0, 0.0, 2.0, 0.0, 0.0],      // dN1/dxi = 4*xi - eta - 4*eta*xi + 2*eta^2 - 1
                        vec![0.0, 0.0, -3.0, 0.0, 0.0, 4.0, 0.0, 2.0, 0.0, 0.0],        // dN2/dxi = 4*eta*xi - 3*eta + 2*eta^2
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, -2.0, 0.0, 0.0],       // dN3/dxi = 4*eta*xi - eta - 2*eta^2
                        vec![4.0, -8.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0],       // dN4/dxi = 8*eta*xi - 8*xi - 4*eta + 4
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0],        // dN5/dxi = -4*eta^2 + 4*eta
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0],        // dN6/dxi = 4*eta - 8*eta*xi
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0],        // dN7/dxi = 4*eta^2 - 4*eta
                    ],
                    vec![
                        vec![-3.0, 5.0, 4.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0],      // dN0/deta = 4*eta + 5*xi - 4*eta*xi - 2*xi^2 - 3
                        vec![0.0, -1.0, 0.0, 0.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0],       // dN1/deta = 4*eta*xi - xi - 2*xi^2
                        vec![0.0, -3.0, 0.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0],        // dN2/deta = 4*eta*xi - 3*xi + 2*xi^2
                        vec![-1.0, -1.0, 4.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0],      // dN3/deta = 4*eta - xi - 4*eta*xi + 2*xi^2 - 1
                        vec![0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // dN4/deta = 4*xi^2 - 4*xi
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0],        // dN5/deta = 4*xi - 8*eta*xi
                        vec![0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // dN6/deta = -4*xi^2 + 4*xi
                        vec![4.0, -4.0, -8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0],       // dN7/deta = 8*eta*xi - 4*xi - 8*eta + 4
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN0/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN2/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN3/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/dpsi = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::BiquadraticQuad => {
                let num_nodes = 9;
                // 2D biquadratic element with center node
                // Basis: [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]
                
                let values = vec![
                    vec![1.0, -3.0, -3.0, 0.0, 2.0, 9.0, 0.0, 2.0, 0.0, 0.0, 0.0, -6.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // N0 = 4*eta^2*xi^2 - 6*eta^2*xi + 2*eta^2 - 6*eta*xi^2 + 9*eta*xi - 3*eta + 2*xi^2 - 3*xi + 1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // N1 = 4*eta^2*xi^2 - 2*eta^2*xi - 6*eta*xi^2 + 3*eta*xi + 2*xi^2 - xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          // N2 = 4*eta^2*xi^2 - 2*eta^2*xi - 2*eta*xi^2 + eta*xi
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 3.0, 0.0, 2.0, 0.0, 0.0, 0.0, -2.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // N3 = 4*eta^2*xi^2 - 6*eta^2*xi + 2*eta^2 - 2*eta*xi^2 + 3*eta*xi - eta
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // N4 = -8*eta^2*xi^2 + 8*eta^2*xi + 12*eta*xi^2 - 12*eta*xi - 4*xi^2 + 4*xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          // N5 = -8*eta^2*xi^2 + 4*eta^2*xi + 8*eta*xi^2 - 4*eta*xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          // N6 = -8*eta^2*xi^2 + 8*eta^2*xi + 4*eta*xi^2 - 4*eta*xi
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -12.0, 0.0, -4.0, 0.0, 0.0, 0.0, 8.0, 0.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // N7 = -8*eta^2*xi^2 + 12*eta^2*xi - 4*eta^2 + 8*eta*xi^2 - 12*eta*xi + 4*eta
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],      // N8 = 16*eta^2*xi^2 - 16*eta^2*xi - 16*eta*xi^2 + 16*eta*xi
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 9.0, 0.0, 0.0, -12.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // dN0/dxi = 9*eta + 4*xi - 12*eta*xi + 8*eta^2*xi - 6*eta^2 - 3
                        vec![-1.0, 4.0, 3.0, 0.0, 0.0, -12.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // dN1/dxi = 3*eta + 4*xi - 12*eta*xi + 8*eta^2*xi - 2*eta^2 - 1
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, -4.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // dN2/dxi = eta - 4*eta*xi + 8*eta^2*xi - 2*eta^2
                        vec![0.0, 0.0, 3.0, 0.0, 0.0, -4.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // dN3/dxi = 3*eta - 4*eta*xi + 8*eta^2*xi - 6*eta^2
                        vec![4.0, -8.0, -12.0, 0.0, 0.0, 24.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dxi = 24*eta*xi - 8*xi - 12*eta - 16*eta^2*xi + 8*eta^2 + 4
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 16.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // dN5/dxi = 16*eta*xi - 4*eta - 16*eta^2*xi + 4*eta^2
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // dN6/dxi = 8*eta*xi - 4*eta - 16*eta^2*xi + 8*eta^2
                        vec![0.0, 0.0, -12.0, 0.0, 0.0, 16.0, 0.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/dxi = 16*eta*xi - 12*eta - 16*eta^2*xi + 12*eta^2
                        vec![0.0, 0.0, 16.0, 0.0, 0.0, -32.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN8/dxi = 16*eta - 32*eta*xi + 32*eta^2*xi - 16*eta^2
                    ],
                    vec![
                        vec![-3.0, 9.0, 4.0, 0.0, -6.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // dN0/deta = 4*eta + 9*xi - 12*eta*xi + 8*eta*xi^2 - 6*xi^2 - 3
                        vec![0.0, 3.0, 0.0, 0.0, -6.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // dN1/deta = 3*xi - 4*eta*xi + 8*eta*xi^2 - 6*xi^2
                        vec![0.0, 1.0, 0.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // dN2/deta = xi - 4*eta*xi + 8*eta*xi^2 - 2*xi^2
                        vec![-1.0, 3.0, 4.0, 0.0, -2.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // dN3/deta = 4*eta + 3*xi - 12*eta*xi + 8*eta*xi^2 - 2*xi^2 - 1
                        vec![0.0, -12.0, 0.0, 0.0, 12.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/deta = 16*eta*xi - 12*xi - 16*eta*xi^2 + 12*xi^2
                        vec![0.0, -4.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],        // dN5/deta = 8*eta*xi - 4*xi - 16*eta*xi^2 + 8*xi^2
                        vec![0.0, -4.0, 0.0, 0.0, 4.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       // dN6/deta = 16*eta*xi - 4*xi - 16*eta*xi^2 + 4*xi^2
                        vec![4.0, -12.0, -8.0, 0.0, 8.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/deta = 24*eta*xi - 12*xi - 8*eta - 16*eta*xi^2 + 8*xi^2 + 4
                        vec![0.0, 16.0, 0.0, 0.0, -16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN8/deta = 16*xi - 32*eta*xi + 32*eta*xi^2 - 16*xi^2
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN0/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN2/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN3/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN8/dpsi = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Tetra => {
                let num_nodes = 4;
                // 3D degree 1: [1, xi, eta, psi]

                let values = vec![
                    vec![1.0, -1.0, -1.0, -1.0],  // N0 = 1 - xi - eta - psi
                    vec![0.0, 1.0, 0.0, 0.0],     // N1 = xi
                    vec![0.0, 0.0, 1.0, 0.0],     // N2 = eta
                    vec![0.0, 0.0, 0.0, 1.0],     // N3 = psi
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN0/dxi = -1
                        vec![1.0, 0.0, 0.0, 0.0],   // dN1/dxi = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dxi = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN3/dxi = 0
                    ],
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN0/deta = -1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/deta = 0
                        vec![1.0, 0.0, 0.0, 0.0],   // dN2/deta = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN3/deta = 0
                    ],
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN0/dpsi = -1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dpsi = 0
                        vec![1.0, 0.0, 0.0, 0.0],   // dN3/dpsi = 1
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticTetra => {
                let num_nodes = 10;
                // 3D degree 2: [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]

                let values = vec![
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 4.0, 4.0, 2.0, 4.0, 2.0],  // N0 = 2*eta^2 + 4*eta*psi + 4*eta*xi - 3*eta + 2*psi^2 + 4*psi*xi - 3*psi + 2*xi^2 - 3*xi + 1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N1 = 2*xi^2 - xi
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0],    // N2 = 2*eta^2 - eta
                    vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0],    // N3 = 2*psi^2 - psi
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, -4.0, 0.0, 0.0, 0.0],  // N4 = 4*xi - 4*eta*xi - 4*psi*xi - 4*xi^2
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0],     // N5 = 4*eta*xi
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0, 0.0],  // N6 = 4*eta - 4*eta*psi - 4*eta*xi - 4*eta^2
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0],   // N7 = 4*psi - 4*eta*psi - 4*psi*xi - 4*psi^2
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0],     // N8 = 4*psi - 4*eta*psi - 4*psi*xi - 4*psi^2
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0],     // N9 = 4*psi - 4*eta*psi - 4*psi*xi - 4*psi^2
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dxi = 4*eta + 4*psi + 4*xi - 3
                        vec![-1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dxi = 4*xi - 1
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dxi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dxi = 0
                        vec![4.0, -8.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dxi = 4 - 4*psi - 8*xi - 4*eta
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dxi = 4*eta
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dxi = -4*eta
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN7/dxi = -4*psi
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8/dxi = 4*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dxi = 0
                    ],
                    vec![
                        vec![-3.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/deta = 4*eta + 4*psi + 4*xi - 3
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/deta = 0
                        vec![-1.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/deta = 4*eta - 1
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/deta = 0
                        vec![0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/deta = -4*xi
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/deta = 4*xi
                        vec![4.0, -4.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN6/deta = 4 - 4*psi - 4*xi - 8*eta
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN7/deta = -4*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8/deta = 0
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/deta = 4*psi
                    ],
                    vec![
                        vec![-3.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dpsi = 4*eta + 4*psi + 4*xi - 3
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dpsi = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dpsi = 0
                        vec![-1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dpsi = 4*psi - 1
                        vec![0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dpsi = -4*xi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dpsi = 0
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dpsi = -4*eta
                        vec![4.0, -4.0, -4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN7/dpsi = 4 - 8*psi - 4*xi - 4*eta
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8/dpsi = 4*xi
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dpsi = 4*eta
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Pyramid => {
                let num_nodes = 5;
                // [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi², xi³, xi²eta, xi²psi, xi*eta², xi*eta*psi, xi*psi², eta³, eta²psi, eta*psi², psi³]

                let values = vec![
                    vec![1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // N0 = (1-xi)(1-eta)(1-psi)
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N1 = xi(1-eta)(1-psi)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // N2 = xi*eta(1-psi)
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N3 = (1-xi)eta(1-psi)
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // N4 = psi
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],   // dN0/dxi = -(1-eta)(1-psi) = -1 + eta + psi - eta*psi
                        vec![1.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],   // dN1/dxi = (1-eta)(1-psi) = 1 - eta - psi + eta*psi
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],    // dN2/dxi = eta(1-psi) = eta - eta*psi
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],    // dN3/dxi = -eta(1-psi) = -eta + eta*psi
                        vec![0.0],                                                  // dN4/dxi = 0
                    ],
                    vec![
                        vec![-1.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],   // dN0/deta = -(1-xi)(1-psi) = -1 + xi + psi - xi*psi
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],    // dN1/deta = -xi(1-psi) = -xi + xi*psi
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],    // dN2/deta = xi(1-psi) = xi - xi*psi
                        vec![1.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],   // dN3/deta = (1-xi)(1-psi) = 1 - xi - psi + xi*psi
                        vec![0.0],                                                  // dN4/deta = 0
                    ],
                    vec![
                        vec![-1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dpsi = -(1-xi)(1-eta) = -1 + xi + eta - xi*eta
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],    // dN1/dpsi = -xi(1-eta) = -xi + xi*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],    // dN2/dpsi = -xi*eta
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],    // dN3/dpsi = -(1-xi)eta = -eta + xi*eta
                        vec![1.0],                                                  // dN4/dpsi = 1
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Wedge => {
                let num_nodes = 6;
                // [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]
                
                let values = vec![
                    vec![1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0],  // N0 = (1-xi-eta)(1-psi)
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],   // N1 = xi(1-psi)
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],   // N2 = eta(1-psi)
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0], // N3 = (1-xi-eta)psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],     // N4 = xi*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],     // N5 = eta*psi
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dxi = -(1-psi) = -1 + psi
                        vec![1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dxi = 1 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dxi = 0
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dxi = -psi
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dxi = psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dxi = 0
                    ],
                    vec![
                        vec![-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/deta = -(1-psi) = -1 + psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/deta = 0
                        vec![1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/deta = 1 - psi
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/deta = -psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/deta = 0
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/deta = psi
                    ],
                    vec![
                        vec![-1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN0/dpsi = -(1-xi-eta) = -1 + xi + eta
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dpsi = -xi
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dpsi = -eta
                        vec![1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dpsi = 1 - xi - eta
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dpsi = xi
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dpsi = eta
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticWedge => {
                let num_nodes = 15;

                let values = vec![
                    // N0 = 2*(1-xi-eta)*(1-psi)*(0.5-xi-eta-psi)
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 4.0, 5.0, 2.0, 5.0, 2.0, 0.0, 0.0, -2.0, 0.0, -4.0, -2.0, 0.0, -2.0, -2.0, 0.0],
                    // N1 = 2*xi*(1-psi)*(xi-psi-0.5)
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
                    // N2 = 2*eta*(1-psi)*(eta-psi-0.5)
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 2.0, 0.0],
                    // N3 = 2*(1-xi-eta)*psi*(psi-xi-eta-0.5)
                    vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 3.0, 0.0, 3.0, 2.0, 0.0, 0.0, 2.0, 0.0, 4.0, -2.0, 0.0, 2.0, -2.0, 0.0],
                    // N4 = 2*xi*psi*(xi+psi-1.5)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
                    // N5 = 2*eta*psi*(eta+psi-1.5)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0],
                    // N6 = 4*xi*(1-xi-eta)*(1-psi)
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    // N7 = 4*xi*eta*(1-psi)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    // N8 = 4*(1-xi-eta)*eta*(1-psi)
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0],
                    // N9 = 4*xi*(1-xi-eta)*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    // N10 = 4*xi*eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    // N11 = 4*(1-xi-eta)*eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -4.0, 0.0, 0.0],
                    // N12 = 4*psi*(1-xi-eta)*(1-psi)
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0],
                    // N13 = 4*psi*xi*(1-psi)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0],
                    // N14 = 4*psi*eta*(1-psi)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0],
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 4.0, 5.0, 0.0, 0.0, -4.0, 0.0, -4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dxi = 4*xi + 4*eta + 5*psi - 4*xi*psi - 4*eta*psi - 2*psi^2 - 3
                        vec![-1.0, 4.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dxi = 4*xi - 3*psi - 1
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN2/dxi = 0
                        vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN3/dxi = 3*psi
                        vec![0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN4/dxi = -3*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dxi = 0
                        vec![4.0, -8.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN6/dxi = -8*xi - 4*eta - 4*psi + 8*xi*psi + 4*eta*psi + 4
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN7/dxi = 4*eta - 4*eta*psi
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN8/dxi = -4*eta + 4*eta*psi
                        vec![0.0, 0.0, 0.0, 4.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dxi = 4*psi - 8*xi*psi - 4*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN10/dxi = 4*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN11/dxi = -4*eta*psi
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN12/dxi = 8*xi*psi + 8*eta*psi + 8*psi^2 - 12*psi + 4
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN13/dxi = 4*psi - 8*xi*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN14/dxi = -8*eta*psi
                    ],
                    vec![
                        vec![-3.0, 4.0, 4.0, 5.0, 0.0, 0.0, -4.0, 0.0, -4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/deta = 4*xi + 4*eta + 5*psi - 4*xi*psi - 4*eta*psi - 2*psi^2 - 3
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN1/deta = 0
                        vec![-1.0, 0.0, 4.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/deta = 4*eta - 3*psi - 1
                        vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN3/deta = 3*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/deta = 0
                        vec![0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN5/deta = -3*psi
                        vec![0.0, -4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN6/deta = -4*xi + 4*xi*psi
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN7/deta = 4*xi - 4*xi*psi
                        vec![4.0, -4.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN8/deta = -4*xi - 8*eta - 4*psi + 4*xi*psi + 8*eta*psi + 4
                        vec![0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN9/deta = -4*xi*psi
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN10/deta = 4*xi*psi
                        vec![0.0, 0.0, 0.0, 4.0, -4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN11/deta = 4*psi - 4*xi*psi - 8*eta*psi
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN12/deta = 8*xi*psi + 8*eta*psi + 8*psi^2 - 12*psi + 4
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN13/deta = -8*xi*psi
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN14/deta = 4*psi - 8*eta*psi
                    ],
                    vec![
                        vec![-3.0, 5.0, 5.0, 4.0, -2.0, -2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dpsi = 5*xi + 5*eta + 4*psi - 2*xi^2 - 2*xi*eta - 2*eta^2 - 3
                        vec![0.0, -3.0, 0.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN1/dpsi = -3*xi + 4*psi + 2*xi^2
                        vec![0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN2/dpsi = -3*eta + 4*psi + 2*eta^2
                        vec![0.0, 3.0, 3.0, 4.0, 2.0, 4.0, -4.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dpsi = 3*xi + 3*eta + 4*psi + 2*xi^2 + 4*xi*eta - 4*xi*psi + 2*eta^2 - 4*eta*psi
                        vec![0.0, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN4/dpsi = -3*xi + 2*xi^2
                        vec![0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN5/dpsi = -3*eta + 2*eta^2
                        vec![0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN6/dpsi = -4*xi + 4*xi^2 + 4*xi*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN7/dpsi = -4*xi*eta
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN8/dpsi = -4*eta + 4*xi*eta + 4*eta^2
                        vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dpsi = 4*xi - 4*xi^2 - 4*xi*eta
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN10/dpsi = 4*xi*eta
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN11/dpsi = 4*eta - 4*xi*eta - 4*eta^2
                        vec![4.0, -4.0, -4.0, -8.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN12/dpsi = 4 - 4*xi - 4*eta - 8*psi + 8*xi*psi + 8*eta*psi
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN13/dpsi = 4*xi - 8*xi*psi
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN14/dpsi = 4*eta - 8*eta*psi
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })

            },

            ElementType::BiquadraticQuadraticWedge => {
                let num_nodes = 18;

                let values = vec![
                    // Node 0
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 4.0, 9.0, 2.0, 9.0, 2.0, 0.0, 0.0, -6.0, 0.0, -12.0, -6.0, 0.0, -6.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // N0 = 4.0*eta^2*psi^2 + 8.0*eta*psi^2*xi - 6.0*eta^2*psi - 6.0*eta*psi^2 + 4.0*psi^2*xi^2 - 12.0*eta*psi*xi - 6.0*psi^2*xi + 2.0*eta^2 + 9.0*eta*psi - 6.0*psi*xi^2 + 4.0*eta*xi + 2.0*psi^2 + 9.0*psi*xi + 2.0*xi^2 - 3.0*eta - 3.0*psi - 3.0*xi + 1.0
                    // Node 1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N1 = -2.0*psi^2*xi - 2.0*psi*xi^2 + 3.0*psi*xi + 2.0*xi^2 - xi
                    // Node 2
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N2 = -2.0*eta^2*psi - 2.0*eta*psi^2 + 2.0*eta^2 + 3.0*eta*psi - eta
                    // Node 3
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 2.0, 0.0, 0.0, 2.0, 0.0, 4.0, -3.0, 0.0, 2.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // N3 = 4.0*eta^2*psi^2 + 8.0*eta*psi^2*xi + 2.0*eta^2*psi - 3.0*eta*psi^2 + 4.0*psi^2*xi^2 + 4.0*eta*psi*xi - 3.0*psi^2*xi + 2.0*psi*xi^2 - eta*psi + 2.0*psi^2 - psi*xi
                    // Node 4
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N4 = 2.0*psi^2*xi^2 + psi^2*xi
                    // Node 5
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0], // N5 = 2.0*eta^2*psi^2 + eta*psi^2
                    // Node 6
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 12.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N6 = -8.0*eta*psi^2*xi + 12.0*eta*psi*xi - 4.0*eta*xi - 8.0*psi^2*xi^2 + 8.0*psi*xi^2 - 12.0*psi*xi - 4.0*xi^2 + 4.0*xi
                    // Node 7
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N7 = 8.0*eta*psi^2*xi - 12.0*eta*psi*xi + 4.0*eta*xi - 4.0*eta*psi
                    // Node 8
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0], // N8 = -8.0*eta^2*psi^2 + 8.0*eta^2*psi - 4.0*eta^2 - 8.0*eta*psi^2*xi + 12.0*eta*psi*xi - 4.0*eta*xi + 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta
                    // Node 9
                    vec![0.0, 0.0, 0.0, 0.0, 2.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, -8.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N9 = 8.0*eta*psi^2*xi - 8.0*eta*psi*xi + 8.0*psi^2*xi^2 - 8.0*psi*xi^2 + 2.0*psi*xi^2 + 2.0*eta*xi + 2.0*xi^2 - 4.0*psi*xi
                    // Node 10
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N10 = 8.0*eta*psi^2*xi - 4.0*eta*psi*xi
                    // Node 11
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, -8.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0], // N11 = 8.0*eta^2*psi^2 - 8.0*eta^2*psi + 8.0*eta*psi^2*xi - 8.0*eta*psi*xi + 2.0*eta*psi^2 - 4.0*eta*psi + 2.0*eta^2 + 2.0*eta*xi
                    // Node 12
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -12.0, 0.0, -12.0, -4.0, 0.0, 0.0, 8.0, 0.0, 16.0, 12.0, 0.0, 8.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0], // N12 = -8.0*eta^2*psi^2 + 8.0*eta^2*psi - 16.0*eta*psi^2*xi + 16.0*eta*psi*xi - 8.0*psi^2*xi^2 + 8.0*psi*xi^2 + 12.0*eta*psi^2 - 12.0*eta*psi + 12.0*psi^2*xi - 12.0*psi*xi - 4.0*psi^2 + 4.0*psi
                    // Node 13
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N13 = -8.0*psi^2*xi^2 + 8.0*psi*xi^2 + 4.0*psi^2*xi - 4.0*psi*xi
                    // Node 14
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0], // N14 = -8.0*eta^2*psi^2 + 8.0*eta^2*psi + 4.0*eta*psi^2 - 4.0*eta*psi
                    // Node 15
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, -16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N15 = 16.0*eta*psi^2*xi - 16.0*eta*psi*xi + 16.0*psi^2*xi^2 - 16.0*psi*xi^2 - 16.0*psi^2*xi + 16.0*psi*xi
                    // Node 16
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N16 = -16.0*eta*psi^2*xi + 16.0*eta*psi*xi
                    // Node 17
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, -16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0], // N17 = 16.0*eta^2*psi^2 - 16.0*eta^2*psi + 16.0*eta*psi^2*xi - 16.0*eta*psi*xi - 16.0*eta*psi^2 + 16.0*eta*psi
                ];

                let derivatives = vec![
                    // d/dxi
                    vec![
                        vec![-3.0, 4.0, 4.0, 9.0, 0.0, 0.0, -12.0, 0.0, 0.0, -6.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dxi = 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta + 8.0*psi^2*xi - 12.0*psi*xi + 4.0*xi - 6.0*psi^2 + 9.0*psi - 3.0
                        vec![-1.0, 4.0, 0.0, 3.0, 0.0, 0.0, -4.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dxi = 4.0*xi - 4.0*psi^2*xi + 3.0*psi - 2.0*psi^2 - 1.0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dxi = 0
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, 4.0, -3.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dxi = 8.0*eta*psi^2 - 4.0*eta*psi + 8.0*psi^2*xi - 4.0*psi*xi - 3.0*psi^2 + 4.0*psi*xi + 2.0*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dxi = 4.0*psi^2*xi - 2.0*psi*xi + 2.0*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dxi = 0
                        vec![4.0, -8.0, -4.0, -12.0, 0.0, 0.0, 16.0, 0.0, -8.0, 8.0, 0.0, 0.0, -16.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dxi = -8.0*eta*psi^2 + 12.0*eta*psi - 4.0*eta - 16.0*psi^2*xi + 16.0*psi*xi - 8.0*xi + 8.0*psi^2 - 12.0*psi + 4.0
                        vec![0.0, 0.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dxi = 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dxi = -8.0*eta*psi^2 + 12.0*eta*psi - 4.0*eta
                        vec![2.0, 4.0, 2.0, -4.0, 0.0, 0.0, 16.0, 0.0, 8.0, -8.0, 0.0, 0.0, 16.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dxi = 8.0*eta*psi^2 - 8.0*eta*psi + 2.0*eta + 16.0*psi^2*xi - 16.0*psi*xi + 4.0*xi - 8.0*psi^2 + 2.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dxi = 8.0*eta*psi^2 - 4.0*eta*psi
                        vec![0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dxi = 8.0*eta*psi^2 - 8.0*eta*psi + 2.0*eta
                        vec![0.0, 0.0, 0.0, -12.0, 0.0, 0.0, -16.0, 0.0, -16.0, 12.0, 0.0, 0.0, 32.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dxi = -16.0*eta*psi^2 + 16.0*eta*psi - 16.0*psi^2*xi + 16.0*psi*xi + 12.0*psi^2 - 12.0*psi
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -16.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dxi = -16.0*psi^2*xi + 16.0*psi*xi + 4.0*psi^2 - 4.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dxi = 0
                        vec![0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 32.0, 0.0, 16.0, -16.0, 0.0, 0.0, -32.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/dxi = 16.0*eta*psi^2 - 16.0*eta*psi + 32.0*psi^2*xi - 32.0*psi*xi - 16.0*psi^2 + 16.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/dxi = -16.0*eta*psi^2 + 16.0*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/dxi = 16.0*eta*psi^2 - 16.0*eta*psi
                    ],
                    vec![
                        // d/deta
                        vec![-3.0, 4.0, 4.0, 9.0, 0.0, 0.0, 0.0, -12.0, 0.0, -6.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/deta = 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta + 8.0*psi^2*xi - 12.0*psi*xi + 4.0*xi - 6.0*psi^2 + 9.0*psi - 3.0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/deta = 0
                        vec![-1.0, 0.0, 4.0, 3.0, 0.0, 0.0, 0.0, -4.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/deta = 4.0*eta - 4.0*eta*psi^2 + 3.0*psi - 2.0*psi^2 - 1.0
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, 4.0, -3.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/deta = 8.0*eta*psi^2 - 4.0*eta*psi + 8.0*psi^2*xi - 4.0*psi*xi - 3.0*psi^2 + 4.0*psi*xi + 2.0*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/deta = 0
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/deta = 4.0*eta*psi^2 - 2.0*eta*psi + 2.0*psi^2 - psi
                        vec![0.0, 4.0, -4.0, -12.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/deta = -8.0*psi^2*xi + 12.0*psi*xi - 4.0*xi
                        vec![0.0, 4.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/deta = 8.0*psi^2*xi - 12.0*psi*xi + 4.0*xi - 4.0*psi^2 + 6.0*psi - 2.0
                        vec![4.0, -8.0, -4.0, -12.0, 0.0, 0.0, -8.0, 16.0, 8.0, 8.0, 0.0, 0.0, 0.0, -16.0, -8.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/deta = -16.0*eta*psi^2 + 16.0*eta*psi - 8.0*eta - 8.0*psi^2*xi + 12.0*psi*xi - 4.0*xi + 8.0*psi^2 - 12.0*psi + 4.0
                        vec![2.0, 2.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/deta = 8.0*psi^2*xi - 8.0*psi*xi + 2.0*xi
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/deta = 8.0*psi^2*xi - 4.0*psi*xi
                        vec![2.0, 2.0, 4.0, -4.0, 0.0, 0.0, 8.0, 16.0, 0.0, -8.0, 0.0, 0.0, 0.0, 16.0, 8.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/deta = 16.0*eta*psi^2 - 16.0*eta*psi + 4.0*eta + 8.0*psi^2*xi - 8.0*psi*xi + 2.0*xi - 8.0*psi^2 + 2.0*psi
                        vec![0.0, 0.0, 0.0, -12.0, 0.0, 0.0, -16.0, 0.0, -16.0, 12.0, 0.0, 0.0, 0.0, 32.0, 32.0, 0.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/deta = -16.0*eta*psi^2 + 16.0*eta*psi - 16.0*psi^2*xi + 16.0*psi*xi + 12.0*psi^2 - 12.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/deta = 0
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, -16.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/deta = -16.0*eta*psi^2 + 16.0*eta*psi + 4.0*psi^2 - 4.0*psi
                        vec![0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, -16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/deta = 16.0*psi^2*xi - 16.0*psi*xi - 16.0*psi^2 + 16.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/deta = -16.0*psi^2*xi + 16.0*psi*xi
                        vec![0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 32.0, 0.0, -16.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/deta = 32.0*eta*psi^2 - 32.0*eta*psi + 16.0*psi^2*xi - 16.0*psi*xi - 16.0*psi^2 + 16.0*psi
                    ],
                    vec![
                        // d/dpsi
                        vec![-3.0, 9.0, 9.0, -3.0, -6.0, -12.0, 8.0, -6.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dpsi = 8.0*eta^2*psi + 8.0*eta*xi*psi - 12.0*eta*xi - 6.0*eta^2 - 6.0*eta*psi + 9.0*eta + 8.0*xi^2*psi - 6.0*xi^2 - 6.0*xi*psi + 9.0*xi + 4.0*psi - 3.0
                        vec![3.0, 0.0, 0.0, -1.0, -2.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dpsi = -4.0*xi^2*psi + 3.0*xi - 2.0*xi^2 - 4.0*xi*psi - psi
                        vec![3.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dpsi = -4.0*eta^2*psi + 3.0*eta - 2.0*eta^2 - 4.0*eta*psi - psi
                        vec![-1.0, 4.0, 4.0, -1.0, 2.0, 4.0, 8.0, 2.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dpsi = 8.0*eta^2*psi + 8.0*eta*xi*psi + 2.0*eta^2 - 6.0*eta*psi + 4.0*eta + 8.0*xi^2*psi + 2.0*xi^2 - 6.0*xi*psi + 4.0*xi + 4.0*psi - 1.0
                        vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dpsi = 4.0*xi^2*psi + xi + 2.0*xi^2 + 4.0*xi*psi
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dpsi = 4.0*eta^2*psi + eta + 2.0*eta^2 + 4.0*eta*psi
                        vec![0.0, -12.0, -4.0, 0.0, 8.0, 12.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dpsi = -16.0*eta*xi*psi + 12.0*eta*xi - 16.0*xi^2*psi + 8.0*xi^2 + 12.0*xi*psi - 12.0*xi - 4.0*eta
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, -12.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dpsi = 16.0*eta*xi*psi - 12.0*eta*xi + 16.0*eta*psi - 4.0*eta
                        vec![0.0, -4.0, -12.0, 0.0, 0.0, 12.0, 0.0, 8.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dpsi = -16.0*eta^2*psi + 8.0*eta^2 - 16.0*eta*xi*psi + 12.0*eta*xi + 12.0*eta*psi - 12.0*eta - 4.0*xi
                        vec![0.0, -4.0, 2.0, 0.0, -8.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dpsi = 16.0*eta*xi*psi - 8.0*eta*xi + 16.0*xi^2*psi - 8.0*xi^2 + 2.0*xi*psi - 4.0*xi + 2.0*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dpsi = 16.0*eta*xi*psi - 4.0*eta*xi
                        vec![0.0, 2.0, -4.0, 0.0, 0.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dpsi = 16.0*eta^2*psi - 8.0*eta^2 + 16.0*eta*xi*psi - 8.0*eta*xi + 2.0*eta*psi - 4.0*eta + 2.0*xi
                        vec![4.0, -12.0, -12.0, 4.0, 8.0, 16.0, -16.0, 8.0, -16.0, -8.0, 0.0, 0.0, -16.0, 0.0, -32.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dpsi = -16.0*eta^2*psi + 8.0*eta^2 - 32.0*eta*xi*psi + 16.0*eta*xi - 16.0*eta*psi + 12.0*eta - 16.0*xi^2*psi + 8.0*xi^2 - 16.0*xi*psi + 12.0*xi - 8.0*psi + 4.0
                        vec![0.0, -4.0, 0.0, 4.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dpsi = -16.0*xi^2*psi + 8.0*xi^2 + 4.0*xi*psi - 4.0*xi + 4.0*psi
                        vec![0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dpsi = -16.0*eta^2*psi + 8.0*eta^2 + 4.0*eta*psi - 4.0*eta + 4.0*psi
                        vec![0.0, 16.0, -16.0, 0.0, -16.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/dpsi = 32.0*eta*xi*psi - 16.0*eta*xi + 32.0*xi^2*psi - 16.0*xi^2 - 16.0*xi*psi + 16.0*xi
                        vec![0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/dpsi = -32.0*eta*xi*psi + 16.0*eta*xi + 16.0*eta
                        vec![0.0, 16.0, -16.0, 0.0, 0.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32.0, 0.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/dpsi = 32.0*eta^2*psi - 16.0*eta^2 + 32.0*eta*xi*psi - 16.0*eta*xi - 16.0*eta*psi + 16.0*eta
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Hexahedron => {
                let num_nodes = 8;

                // Polynomial basis: [1, xi, eta, psi, xi², xi*eta, xi*psi, eta², eta*psi, psi²]
                let values= vec![
                    vec![1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // N0 = eta*psi - psi - xi - eta + eta*xi + psi*xi - eta*psi*xi + 1
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N1 = xi - eta*xi - psi*xi + eta*psi*xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N2 = eta*xi - eta*psi*xi
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N3 = eta - eta*psi - eta*xi + eta*psi*xi
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N4 = psi - eta*psi - psi*xi + eta*psi*xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N5 = psi*xi - eta*psi*xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N6 = eta*psi*xi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N7 = eta*psi - eta*psi*xi
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],   // dN0/dxi = eta + psi - eta*psi - 1
                        vec![1.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],   // dN1/dxi = eta*psi - psi - eta + 1
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],    // dN2/dxi = eta - eta*psi
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],    // dN3/dxi = eta*psi - eta
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],    // dN4/dxi = eta*psi - psi
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],    // dN5/dxi = psi - eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],     // dN6/dxi = eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],    // dN7/dxi = -eta*psi
                    ],
                    vec![
                        vec![-1.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],   // dN0/deta = psi + xi - psi*xi - 1
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],    // dN1/deta = psi*xi - xi
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],    // dN2/deta = xi - psi*xi
                        vec![1.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],   // dN3/deta = psi*xi - xi - psi + 1
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],    // dN4/deta = psi*xi - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],    // dN5/deta = -psi*xi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],     // dN6/deta = psi*xi
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],    // dN7/deta = psi - psi*xi
                    ],
                    vec![
                        vec![-1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dpsi = eta + xi - eta*xi - 1
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],    // dN1/dpsi = eta*xi - xi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],    // dN2/dpsi = -eta*xi
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],    // dN3/dpsi = eta*xi - eta
                        vec![1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dpsi = eta*xi - xi - eta + 1
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],    // dN5/dpsi = xi - eta*xi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],     // dN6/dpsi = eta*xi
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],    // dN7/dpsi = eta - eta*xi
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            }


            ElementType::QuadraticHexahedron => { // needs to be checked
                let num_nodes = 20;

                let values = vec![
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 5.0, 5.0, 2.0, 5.0, 2.0, 0.0, -2.0, -2.0, -2.0, -7.0, -2.0, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, -2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N0 = 2*xi*eta^2*psi + 2*xi*eta*psi^2 + 2*xi^2*eta*psi - 2*xi*eta^2 - 2*eta*psi^2 - 7*xi*eta*psi + 2*eta^2 + 2*psi^2 + 5*xi*eta + 5*eta*psi + 5*xi*psi - 2*xi^2*eta - 2*xi^2*psi - 2*eta^2*psi + 2*xi^2 - 3*xi - 3*eta - 3*psi + 1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N1 = -2*xi*eta^2*psi - 2*xi^2*eta*psi - 2*xi^2*eta + 2*xi^2 - xi*eta*psi - 2*xi*eta^2 - 2*xi^2*psi - xi
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N2 = -2*xi*eta*psi^2 - 2*xi*eta^2*psi - 2*eta^2*psi - xi*eta*psi + 2*eta^2 - 2*eta*psi^2 - eta
                    vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N3 = -2*xi*eta*psi^2 - 2*xi^2*psi - xi*eta*psi + 2*psi^2 - 2*eta*psi^2 - psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N4 = 2*xi^2*eta - 2*xi*psi^2 - xi*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N5 = 2*eta^2*psi - 2*eta*psi^2 - eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // N6 = 2*xi*eta^2*psi + 2*xi*eta*psi^2 + xi*eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // N7 = xi*eta*psi
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N8 = -4*xi*eta*psi^2 + 4*xi*eta*psi - 4*xi^2 - 4*xi*eta + 4*xi
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N9 = -4*xi*eta^2*psi + 4*xi*eta*psi - 4*eta^2 + 4*eta - 4*eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N10 = 4*xi*eta
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N11 = 4*xi*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N12 = 4*eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N13 = 0
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N14 = 0
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N15 = 0
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N16 = -4*xi*eta*psi^2 + 4*xi*eta*psi - 4*psi^2 + 4*psi - 4*eta*psi^2
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N17 = -4*xi*eta*psi
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N18 = 0
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N19 = 0
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 5.0, 5.0, 4.0, -7.0, -3.0, -4.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0 = 2*eta^2*psi + 2*eta*psi^2 + 2*xi*eta*psi + 2*eta^2 + 2*psi^2 + 2*xi*eta + 2*xi*psi - 2*eta*psi^2 - 7*eta*psi - 4*xi*eta - 2*xi*psi - 4*eta^2 + 5*eta - 2*psi^2 + 5*psi + 4*xi - 3
                        vec![-1.0, 0.0, 0.0, 4.0, -1.0, -3.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1 = -2*eta^2*psi - 4*xi*eta*psi - 2*xi*psi^2 - 4*xi*eta - eta*psi - 2*eta^2 - 2*psi^2 + 4*xi - 1
                        vec![0.0, 0.0, 0.0, 0.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2 = -2*eta*psi^2 - 2*eta^2*psi - eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3 = -2*eta*psi^2 - 4*xi*psi - eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, -1.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4 = 4*xi*eta - 2*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN5 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6 = 2*eta^2*psi + 2*eta*psi^2 + eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7 = eta*psi
                        vec![4.0, -4.0, 0.0, -8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8 = -4*eta*psi^2 + 4*eta*psi - 4*eta + 4
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN9 = -4*eta^2*psi + 4*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN10 = 4*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN11 = 4*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN12 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN13 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN14 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN15 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN16 = -4*eta*psi^2 + 4*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN17 = -4*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN18 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN19 = 0
                    ],
                    vec![
                        vec![-3.0, 5.0, 5.0, -4.0, 5.0, -7.0, -2.0, 4.0, -3.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0 = 2*xi*eta^2 + 2*xi*psi^2 + 2*eta*psi^2 + 2*xi^2 + 2*psi^2 + 2*xi*eta + 2*eta*psi - 2*xi^2*psi - 7*xi*psi - 2*xi^2 - 4*eta*psi - 4*xi*eta + 5*psi + 5*xi + 4*eta - 3
                        vec![0.0, 0.0, 0.0, -1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1 = -2*xi*psi^2 - 2*xi^2*psi - xi*psi - 4*xi*eta
                        vec![-1.0, 0.0, 0.0, 0.0, -1.0, -3.0, 0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2 = -2*xi*psi^2 - 4*xi*eta*psi - 2*eta*psi^2 - xi*psi - 2*xi*eta - 2*psi^2 + 4*eta - 1
                        vec![0.0, 0.0, 0.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3 = -2*xi*psi^2 - 2*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN4 = 2*xi^2
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN5 = 4*eta*psi - 2*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6 = 2*xi*eta*psi + 2*xi*psi^2 + xi*psi
                        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7 = xi*psi
                        vec![0.0, -4.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN8 = -4*xi*psi^2 + 4*xi*psi - 4*xi
                        vec![4.0, 0.0, -8.0, 0.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN9 = -4*xi*eta*psi + 4*xi*psi - 8*eta + 4
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN10 = 4*xi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN11 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN12 = 4*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN13 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN14 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN15 = 0
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN16 = -4*xi*psi^2 + 4*xi*psi - 4*psi^2
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN17 = -4*xi*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN18 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN19 = 0
                    ],
                    vec![
                        vec![-3.0, 5.0, 5.0, -3.0, -7.0, 4.0, 2.0, -2.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0 = 2*xi*eta^2 + 2*xi*eta*psi + 2*xi^2*eta + 2*xi^2 + 2*eta^2 + 2*xi*eta + 2*xi*psi - 2*xi*eta^2 - 7*xi*eta - 2*eta^2*psi - 4*xi*psi - 4*eta*psi + 5*xi + 5*eta + 4*psi - 3
                        vec![0.0, 0.0, 0.0, -1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1 = -2*xi*eta^2 - 2*xi^2*eta - xi*eta - 4*xi*psi - 2*xi^2
                        vec![0.0, 0.0, 0.0, 0.0, -1.0, -3.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2 = -4*xi*eta*psi - 2*xi*eta^2 - 2*eta^2*psi - xi*eta - 2*eta^2 - 4*eta*psi + 3*eta
                        vec![-1.0, 0.0, 0.0, 0.0, -1.0, -4.0, -2.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3 = -4*xi*eta*psi - 2*xi^2*psi - xi*eta - 4*eta*psi + 4*psi - 1
                        vec![0.0, 0.0, 0.0, -1.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN4 = 2*xi^2 - 4*xi*psi - xi
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN5 = 2*eta^2 - 4*eta*psi - eta
                        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6 = 2*xi*eta^2 + 2*xi*eta*psi + xi*eta
                        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7 = xi*eta
                        vec![0.0, 4.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8 = -8*xi*eta*psi + 4*xi*eta + 4*xi
                        vec![0.0, 0.0, 0.0, 4.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9 = -4*xi*eta^2 + 4*xi*eta - 4*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN10 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN11 = 4*xi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN12 = 4*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN13 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN14 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN15 = 0
                        vec![4.0, -4.0, -4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16 = -8*eta*psi - 4*xi*eta + 4*xi + 4*eta - 4
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN17 = -4*xi*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN18 = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN19 = 0
                    ],
                ];


                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::BiquadraticQuadraticHexahedron => { // needs to be checked
                let num_nodes = 24;

                let values = vec![
                    // Node 0
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 4.0, 9.0, 2.0, 9.0, 2.0, 0.0, 0.0, -6.0, 0.0, -12.0, -6.0, 0.0, -6.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // N0 = 4.0*eta^2*psi^2 + 8.0*eta*psi^2*xi - 6.0*eta^2*psi - 6.0*eta*psi^2 + 4.0*psi^2*xi^2 - 12.0*eta*psi*xi - 6.0*psi^2*xi + 2.0*eta^2 + 9.0*eta*psi - 6.0*psi*xi^2 + 4.0*eta*xi + 2.0*psi^2 + 9.0*psi*xi + 2.0*xi^2 - 3.0*eta - 3.0*psi - 3.0*xi + 1.0
                    // Node 1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N1 = -2.0*psi^2*xi - 2.0*psi*xi^2 + 3.0*psi*xi + 2.0*xi^2 - xi
                    // Node 2
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N2 = -2.0*eta^2*psi - 2.0*eta*psi^2 + 2.0*eta^2 + 3.0*eta*psi - eta
                    // Node 3
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 2.0, 0.0, 0.0, 2.0, 0.0, 4.0, -3.0, 0.0, 2.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // N3 = 4.0*eta^2*psi^2 + 8.0*eta*psi^2*xi + 2.0*eta^2*psi - 3.0*eta*psi^2 + 4.0*psi^2*xi^2 + 4.0*eta*psi*xi - 3.0*psi^2*xi + 2.0*psi*xi^2 - eta*psi + 2.0*psi^2 - psi*xi
                    // Node 4
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N4 = 2.0*psi^2*xi^2 + psi^2*xi
                    // Node 5
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0], // N5 = 2.0*eta^2*psi^2 + eta*psi^2
                    // Node 6
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 12.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N6 = -8.0*eta*psi^2*xi + 12.0*eta*psi*xi - 4.0*eta*xi - 8.0*psi^2*xi^2 + 8.0*psi*xi^2 - 12.0*psi*xi - 4.0*xi^2 + 4.0*xi
                    // Node 7
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N7 = 8.0*eta*psi^2*xi - 12.0*eta*psi*xi + 4.0*eta*xi - 4.0*eta*psi
                    // Node 8
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0], // N8 = -8.0*eta^2*psi^2 + 8.0*eta^2*psi - 4.0*eta^2 - 8.0*eta*psi^2*xi + 12.0*eta*psi*xi - 4.0*eta*xi + 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta
                    // Node 9
                    vec![0.0, 0.0, 0.0, 0.0, 2.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, -8.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N9 = 8.0*eta*psi^2*xi - 8.0*eta*psi*xi + 8.0*psi^2*xi^2 - 8.0*psi*xi^2 + 2.0*psi*xi^2 + 2.0*eta*xi + 2.0*xi^2 - 4.0*psi*xi
                    // Node 10
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N10 = 8.0*eta*psi^2*xi - 4.0*eta*psi*xi
                    // Node 11
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, -8.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0], // N11 = 8.0*eta^2*psi^2 - 8.0*eta^2*psi + 8.0*eta*psi^2*xi - 8.0*eta*psi*xi + 2.0*eta*psi^2 - 4.0*eta*psi + 2.0*eta^2 + 2.0*eta*xi
                    // Node 12
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -12.0, 0.0, -12.0, -4.0, 0.0, 0.0, 8.0, 0.0, 16.0, 12.0, 0.0, 8.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0], // N12 = -8.0*eta^2*psi^2 + 8.0*eta^2*psi - 16.0*eta*psi^2*xi + 16.0*eta*psi*xi - 8.0*psi^2*xi^2 + 8.0*psi*xi^2 + 12.0*eta*psi^2 - 12.0*eta*psi + 12.0*psi^2*xi - 12.0*psi*xi - 4.0*psi^2 + 4.0*psi
                    // Node 13
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N13 = -8.0*psi^2*xi^2 + 8.0*psi*xi^2 + 4.0*psi^2*xi - 4.0*psi*xi
                    // Node 14
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0], // N14 = -8.0*eta^2*psi^2 + 8.0*eta^2*psi + 4.0*eta*psi^2 - 4.0*eta*psi
                    // Node 15
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, -16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N15 = 16.0*eta*psi^2*xi - 16.0*eta*psi*xi + 16.0*psi^2*xi^2 - 16.0*psi*xi^2 - 16.0*psi^2*xi + 16.0*psi*xi
                    // Node 16
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N16 = -16.0*eta*psi^2*xi + 16.0*eta*psi*xi
                    // Node 17
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, -16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0], // N17 = 16.0*eta^2*psi^2 - 16.0*eta^2*psi + 16.0*eta*psi^2*xi - 16.0*eta*psi*xi - 16.0*eta*psi^2 + 16.0*eta*psi
                ];

                let derivatives = vec![
                    // d/dxi
                    vec![
                        vec![-3.0, 4.0, 4.0, 9.0, 0.0, 0.0, -12.0, 0.0, 0.0, -6.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dxi = 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta + 8.0*psi^2*xi - 12.0*psi*xi + 4.0*xi - 6.0*psi^2 + 9.0*psi - 3.0
                        vec![-1.0, 4.0, 0.0, 3.0, 0.0, 0.0, -4.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dxi = 4.0*xi - 4.0*psi^2*xi + 3.0*psi - 2.0*psi^2 - 1.0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dxi = 0
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, 4.0, -3.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dxi = 8.0*eta*psi^2 - 4.0*eta*psi + 8.0*psi^2*xi - 4.0*psi*xi - 3.0*psi^2 + 4.0*psi*xi + 2.0*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dxi = 4.0*psi^2*xi - 2.0*psi*xi + 2.0*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dxi = 0
                        vec![4.0, -8.0, -4.0, -12.0, 0.0, 0.0, 16.0, 0.0, -8.0, 8.0, 0.0, 0.0, -16.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dxi = -8.0*eta*psi^2 + 12.0*eta*psi - 4.0*eta - 16.0*psi^2*xi + 16.0*psi*xi - 8.0*xi + 8.0*psi^2 - 12.0*psi + 4.0
                        vec![0.0, 0.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dxi = 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dxi = -8.0*eta*psi^2 + 12.0*eta*psi - 4.0*eta
                        vec![2.0, 4.0, 2.0, -4.0, 0.0, 0.0, 16.0, 0.0, 8.0, -8.0, 0.0, 0.0, 16.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dxi = 8.0*eta*psi^2 - 8.0*eta*psi + 2.0*eta + 16.0*psi^2*xi - 16.0*psi*xi + 4.0*xi - 8.0*psi^2 + 2.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dxi = 8.0*eta*psi^2 - 4.0*eta*psi
                        vec![0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dxi = 8.0*eta*psi^2 - 8.0*eta*psi + 2.0*eta
                        vec![0.0, 0.0, 0.0, -12.0, 0.0, 0.0, -16.0, 0.0, -16.0, 12.0, 0.0, 0.0, 32.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dxi = -16.0*eta*psi^2 + 16.0*eta*psi - 16.0*psi^2*xi + 16.0*psi*xi + 12.0*psi^2 - 12.0*psi
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -16.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dxi = -16.0*psi^2*xi + 16.0*psi*xi + 4.0*psi^2 - 4.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dxi = 0
                        vec![0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 32.0, 0.0, 16.0, -16.0, 0.0, 0.0, -32.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/dxi = 16.0*eta*psi^2 - 16.0*eta*psi + 32.0*psi^2*xi - 32.0*psi*xi - 16.0*psi^2 + 16.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/dxi = -16.0*eta*psi^2 + 16.0*eta*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/dxi = 16.0*eta*psi^2 - 16.0*eta*psi
                    ],
                    vec![
                        // d/deta
                        vec![-3.0, 4.0, 4.0, 9.0, 0.0, 0.0, 0.0, -12.0, 0.0, -6.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/deta = 8.0*eta*psi^2 - 12.0*eta*psi + 4.0*eta + 8.0*psi^2*xi - 12.0*psi*xi + 4.0*xi - 6.0*psi^2 + 9.0*psi - 3.0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/deta = 0
                        vec![-1.0, 0.0, 4.0, 3.0, 0.0, 0.0, 0.0, -4.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/deta = 4.0*eta - 4.0*eta*psi^2 + 3.0*psi - 2.0*psi^2 - 1.0
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, 4.0, -3.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/deta = 8.0*eta*psi^2 - 4.0*eta*psi + 8.0*psi^2*xi - 4.0*psi*xi - 3.0*psi^2 + 4.0*psi*xi + 2.0*psi^2 - psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/deta = 0
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/deta = 4.0*eta*psi^2 - 2.0*eta*psi + 2.0*psi^2 - psi
                        vec![0.0, 4.0, -4.0, -12.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/deta = -8.0*psi^2*xi + 12.0*psi*xi - 4.0*xi
                        vec![0.0, 4.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/deta = 8.0*psi^2*xi - 12.0*psi*xi + 4.0*xi - 4.0*psi^2 + 6.0*psi - 2.0
                        vec![4.0, -8.0, -4.0, -12.0, 0.0, 0.0, -8.0, 16.0, 8.0, 8.0, 0.0, 0.0, 0.0, -16.0, -8.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/deta = -16.0*eta*psi^2 + 16.0*eta*psi - 8.0*eta - 8.0*psi^2*xi + 12.0*psi*xi - 4.0*xi + 8.0*psi^2 - 12.0*psi + 4.0
                        vec![2.0, 2.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/deta = 8.0*psi^2*xi - 8.0*psi*xi + 2.0*xi
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/deta = 8.0*psi^2*xi - 4.0*psi*xi
                        vec![2.0, 2.0, 4.0, -4.0, 0.0, 0.0, 8.0, 16.0, 0.0, -8.0, 0.0, 0.0, 0.0, 16.0, 8.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/deta = 16.0*eta*psi^2 - 16.0*eta*psi + 4.0*eta + 8.0*psi^2*xi - 8.0*psi*xi + 2.0*xi - 8.0*psi^2 + 2.0*psi
                        vec![0.0, 0.0, 0.0, -12.0, 0.0, 0.0, -16.0, 0.0, -16.0, 12.0, 0.0, 0.0, 0.0, 32.0, 32.0, 0.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/deta = -16.0*eta*psi^2 + 16.0*eta*psi - 16.0*psi^2*xi + 16.0*psi*xi + 12.0*psi^2 - 12.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/deta = 0
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, -16.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/deta = -16.0*eta*psi^2 + 16.0*eta*psi + 4.0*psi^2 - 4.0*psi
                        vec![0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, -16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/deta = 16.0*psi^2*xi - 16.0*psi*xi - 16.0*psi^2 + 16.0*psi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/deta = -16.0*psi^2*xi + 16.0*psi*xi
                        vec![0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 32.0, 0.0, -16.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/deta = 32.0*eta*psi^2 - 32.0*eta*psi + 16.0*psi^2*xi - 16.0*psi*xi - 16.0*psi^2 + 16.0*psi
                    ],
                    vec![
                        // d/dpsi
                        vec![-3.0, 9.0, 9.0, -3.0, -6.0, -12.0, 8.0, -6.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dpsi = 8.0*eta^2*psi + 8.0*eta*xi*psi - 12.0*eta*xi - 6.0*eta^2 - 6.0*eta*psi + 9.0*eta + 8.0*xi^2*psi - 6.0*xi^2 - 6.0*xi*psi + 9.0*xi + 4.0*psi - 3.0
                        vec![3.0, 0.0, 0.0, -1.0, -2.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dpsi = -4.0*xi^2*psi + 3.0*xi - 2.0*xi^2 - 4.0*xi*psi - psi
                        vec![3.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dpsi = -4.0*eta^2*psi + 3.0*eta - 2.0*eta^2 - 4.0*eta*psi - psi
                        vec![-1.0, 4.0, 4.0, -1.0, 2.0, 4.0, 8.0, 2.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dpsi = 8.0*eta^2*psi + 8.0*eta*xi*psi + 2.0*eta^2 - 6.0*eta*psi + 4.0*eta + 8.0*xi^2*psi + 2.0*xi^2 - 6.0*xi*psi + 4.0*xi + 4.0*psi - 1.0
                        vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dpsi = 4.0*xi^2*psi + xi + 2.0*xi^2 + 4.0*xi*psi
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dpsi = 4.0*eta^2*psi + eta + 2.0*eta^2 + 4.0*eta*psi
                        vec![0.0, -12.0, -4.0, 0.0, 8.0, 12.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dpsi = -16.0*eta*xi*psi + 12.0*eta*xi - 16.0*xi^2*psi + 8.0*xi^2 + 12.0*xi*psi - 12.0*xi - 4.0*eta
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, -12.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dpsi = 16.0*eta*xi*psi - 12.0*eta*xi + 16.0*eta*psi - 4.0*eta
                        vec![0.0, -4.0, -12.0, 0.0, 0.0, 12.0, 0.0, 8.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dpsi = -16.0*eta^2*psi + 8.0*eta^2 - 16.0*eta*xi*psi + 12.0*eta*xi + 12.0*eta*psi - 12.0*eta - 4.0*xi
                        vec![0.0, -4.0, 2.0, 0.0, -8.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dpsi = 16.0*eta*xi*psi - 8.0*eta*xi + 16.0*xi^2*psi - 8.0*xi^2 + 2.0*xi*psi - 4.0*xi + 2.0*eta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dpsi = 16.0*eta*xi*psi - 4.0*eta*xi
                        vec![0.0, 2.0, -4.0, 0.0, 0.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dpsi = 16.0*eta^2*psi - 8.0*eta^2 + 16.0*eta*xi*psi - 8.0*eta*xi + 2.0*eta*psi - 4.0*eta + 2.0*xi
                        vec![4.0, -12.0, -12.0, 4.0, 8.0, 16.0, -16.0, 8.0, -16.0, -8.0, 0.0, 0.0, -16.0, 0.0, -32.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dpsi = -16.0*eta^2*psi + 8.0*eta^2 - 32.0*eta*xi*psi + 16.0*eta*xi - 16.0*eta*psi + 12.0*eta - 16.0*xi^2*psi + 8.0*xi^2 - 16.0*xi*psi + 12.0*xi - 8.0*psi + 4.0
                        vec![0.0, -4.0, 0.0, 4.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dpsi = -16.0*xi^2*psi + 8.0*xi^2 + 4.0*xi*psi - 4.0*xi + 4.0*psi
                        vec![0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dpsi = -16.0*eta^2*psi + 8.0*eta^2 + 4.0*eta*psi - 4.0*eta + 4.0*psi
                        vec![0.0, 16.0, -16.0, 0.0, -16.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/dpsi = 32.0*eta*xi*psi - 16.0*eta*xi + 32.0*xi^2*psi - 16.0*xi^2 - 16.0*xi*psi + 16.0*xi
                        vec![0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/dpsi = -32.0*eta*xi*psi + 16.0*eta*xi + 16.0*eta
                        vec![0.0, 16.0, -16.0, 0.0, 0.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32.0, 0.0, 0.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/dpsi = 32.0*eta^2*psi - 16.0*eta^2 + 32.0*eta*xi*psi - 16.0*eta*xi - 16.0*eta*psi + 16.0*eta
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::TriquadraticHexahedron => {  // needs to be checked
                let num_nodes = 27;

                let values = vec![
                    vec![1.0, -3.0, 2.0, -3.0, 9.0, -6.0, 2.0, -6.0, 4.0, -3.0, 9.0, -6.0, 9.0, -27.0, 18.0, -6.0, 18.0, -12.0, 2.0, -6.0, 4.0, -6.0, 18.0, -12.0, 4.0, -12.0, 8.0], // N0
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 3.0, -2.0, 3.0, -9.0, 6.0, -2.0, 6.0, -4.0, 2.0, -6.0, 4.0, -6.0, 18.0, -12.0, 4.0, -12.0, 8.0], // N1
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -3.0, 2.0, -2.0, 6.0, -4.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 4.0, -12.0, 8.0], // N2
                    vec![0.0, 0.0, 0.0, -1.0, 3.0, -2.0, 2.0, -6.0, 4.0, 0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -6.0, 18.0, -12.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 4.0, -12.0, 8.0], // N3
                    vec![0.0, -1.0, 2.0, 0.0, 3.0, -6.0, 0.0, -2.0, 4.0, 0.0, 3.0, -6.0, 0.0, -9.0, 18.0, 0.0, 6.0, -12.0, 0.0, -2.0, 4.0, 0.0, 6.0, -12.0, 0.0, -4.0, 8.0], // N4
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 0.0, -3.0, 6.0, 0.0, 2.0, -4.0, 0.0, -2.0, 4.0, 0.0, 6.0, -12.0, 0.0, -4.0, 8.0], // N5
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -4.0, 8.0], // N6
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 0.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 6.0, -12.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -4.0, 8.0], // N7
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, -4.0, 12.0, -8.0, 12.0, -36.0, 24.0, -8.0, 24.0, -16.0], // N8
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -8.0, 24.0, -16.0], // N9
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -8.0, 24.0, -16.0], // N10
                    vec![0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -4.0, 12.0, -8.0, 0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 12.0, -36.0, 24.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -8.0, 24.0, -16.0], // N11
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 4.0, -8.0, 0.0, -12.0, 24.0, 0.0, 8.0, -16.0], // N12
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0], // N13
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 8.0, -16.0], // N14
                    vec![0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -12.0, 24.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0], // N15
                    vec![0.0, 4.0, -4.0, 0.0, -12.0, 12.0, 0.0, 8.0, -8.0, 0.0, -12.0, 12.0, 0.0, 36.0, -36.0, 0.0, -24.0, 24.0, 0.0, 8.0, -8.0, 0.0, -24.0, 24.0, 0.0, 16.0, -16.0], // N16
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 12.0, -12.0, 0.0, -8.0, 8.0, 0.0, 8.0, -8.0, 0.0, -24.0, 24.0, 0.0, 16.0, -16.0], // N17
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -8.0, 8.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 16.0, -16.0], // N18
                    vec![0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 8.0, -8.0, 0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -24.0, 24.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 16.0, -16.0], // N19
                    vec![0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -16.0, 16.0, 0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 48.0, -48.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -32.0, 32.0], // N20
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -32.0, 32.0], // N21
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, -16.0, 16.0, 0.0, 48.0, -48.0, 0.0, -32.0, 32.0], // N22
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -32.0, 32.0], // N23
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, -16.0, 48.0, -32.0, 16.0, -48.0, 32.0], // N24
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -16.0, 32.0], // N25
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, -64.0, 64.0, 0.0, 64.0, -64.0], // N26
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 9.0, -6.0, 9.0, -27.0, 18.0, -6.0, 18.0, -12.0, 4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dxi
                        vec![-1.0, 3.0, -2.0, 3.0, -9.0, 6.0, -2.0, 6.0, -4.0, 4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dxi
                        vec![0.0, 0.0, 0.0, 1.0, -3.0, 2.0, -2.0, 6.0, -4.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dxi
                        vec![0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -6.0, 18.0, -12.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dxi
                        vec![0.0, 3.0, -6.0, 0.0, -9.0, 18.0, 0.0, 6.0, -12.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dxi
                        vec![0.0, 1.0, -2.0, 0.0, -3.0, 6.0, 0.0, 2.0, -4.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dxi
                        vec![0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dxi
                        vec![0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 6.0, -12.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dxi
                        vec![4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, -8.0, 24.0, -16.0, 24.0, -72.0, 48.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dxi
                        vec![0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dxi
                        vec![0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dxi
                        vec![0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 12.0, -36.0, 24.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN11/dxi
                        vec![0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0, 0.0, -24.0, 48.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN12/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN13/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN14/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -12.0, 24.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN15/dxi
                        vec![0.0, -12.0, 12.0, 0.0, 36.0, -36.0, 0.0, -24.0, 24.0, 0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN16/dxi
                        vec![0.0, -4.0, 4.0, 0.0, 12.0, -12.0, 0.0, -8.0, 8.0, 0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN17/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -8.0, 8.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN18/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -24.0, 24.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN19/dxi
                        vec![0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 48.0, -48.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN20/dxi
                        vec![0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN21/dxi
                        vec![0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, -32.0, 32.0, 0.0, 96.0, -96.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN22/dxi
                        vec![0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN23/dxi
                        vec![0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, -32.0, 96.0, -64.0, 32.0, -96.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN24/dxi
                        vec![0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -64.0, 0.0, -32.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN25/dxi
                        vec![0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, -128.0, 128.0, 0.0, 128.0, -128.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN26/dxi
                    ],
                    vec![
                        vec![-3.0, 9.0, -6.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 9.0, -27.0, 18.0, -12.0, 36.0, -24.0, 0.0, 0.0, 0.0, -6.0, 18.0, -12.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0], // dN0/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -4.0, 12.0, -8.0, 0.0, 0.0, 0.0, -6.0, 18.0, -12.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0], // dN1/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -3.0, 2.0, -4.0, 12.0, -8.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0], // dN2/deta
                        vec![-1.0, 3.0, -2.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -12.0, 36.0, -24.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0], // dN3/deta
                        vec![0.0, 3.0, -6.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -9.0, 18.0, 0.0, 12.0, -24.0, 0.0, 0.0, 0.0, 0.0, 6.0, -12.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0], // dN4/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 6.0, -12.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0], // dN5/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0], // dN6/deta
                        vec![0.0, 1.0, -2.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 12.0, -24.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0], // dN7/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 16.0, -48.0, 32.0, 0.0, 0.0, 0.0, 12.0, -36.0, 24.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0], // dN8/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0], // dN9/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 16.0, -48.0, 32.0, 0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0], // dN10/deta
                        vec![4.0, -12.0, 8.0, -8.0, 24.0, -16.0, 0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 24.0, -72.0, 48.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0], // dN11/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, -12.0, 24.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],    // dN12/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],    // dN13/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],    // dN14/deta
                        vec![0.0, -4.0, 8.0, 0.0, 8.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -24.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],  // dN15/deta
                        vec![0.0, -12.0, 12.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 36.0, -36.0, 0.0, -48.0, 48.0, 0.0, 0.0, 0.0, 0.0, -24.0, 24.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],  // dN16/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -16.0, 16.0, 0.0, 0.0, 0.0, 0.0, -24.0, 24.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],    // dN17/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -16.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],    // dN18/deta
                        vec![0.0, -4.0, 4.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -48.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],  // dN19/deta
                        vec![0.0, 16.0, -16.0, 0.0, -32.0, 32.0, 0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 96.0, -96.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],  // dN20/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],    // dN21/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 64.0, -64.0, 0.0, 0.0, 0.0, 0.0, 48.0, -48.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],    // dN22/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 64.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],    // dN23/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -32.0, 96.0, -64.0, 0.0, 0.0, 0.0, -16.0, 48.0, -32.0, 32.0, -96.0, 64.0, 0.0, 0.0, 0.0],  // dN24/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 32.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -32.0, 64.0, 0.0, 0.0, 0.0],    // dN25/deta
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -128.0, 128.0, 0.0, 0.0, 0.0, 0.0, -64.0, 64.0, 0.0, 128.0, -128.0, 0.0, 0.0, 0.0],    // dN26/deta
                    ],
                    vec![
                        vec![-3.0, 4.0, 0.0, 9.0, -12.0, 0.0, -6.0, 8.0, 0.0, 9.0, -12.0, 0.0, -27.0, 36.0, 0.0, 18.0, -24.0, 0.0, -6.0, 8.0, 0.0, 18.0, -24.0, 0.0, -12.0, 16.0, 0.0], // dN0/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -4.0, 0.0, -9.0, 12.0, 0.0, 6.0, -8.0, 0.0, -6.0, 8.0, 0.0, 18.0, -24.0, 0.0, -12.0, 16.0, 0.0], // dN1/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 6.0, -8.0, 0.0, 0.0, 0.0, 0.0, 6.0, -8.0, 0.0, -12.0, 16.0, 0.0], // dN2/dpsi
                        vec![0.0, 0.0, 0.0, 3.0, -4.0, 0.0, -6.0, 8.0, 0.0, 0.0, 0.0, 0.0, -9.0, 12.0, 0.0, 18.0, -24.0, 0.0, 0.0, 0.0, 0.0, 6.0, -8.0, 0.0, -12.0, 16.0, 0.0], // dN3/dpsi
                        vec![-1.0, 4.0, 0.0, 3.0, -12.0, 0.0, -2.0, 8.0, 0.0, 3.0, -12.0, 0.0, -9.0, 36.0, 0.0, 6.0, -24.0, 0.0, -2.0, 8.0, 0.0, 6.0, -24.0, 0.0, -4.0, 16.0, 0.0], // dN4/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -4.0, 0.0, -3.0, 12.0, 0.0, 2.0, -8.0, 0.0, -2.0, 8.0, 0.0, 6.0, -24.0, 0.0, -4.0, 16.0, 0.0], // dN5/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 4.0, 0.0, 2.0, -8.0, 0.0, 0.0, 0.0, 0.0, 2.0, -8.0, 0.0, -4.0, 16.0, 0.0], // dN6/dpsi
                        vec![0.0, 0.0, 0.0, 1.0, -4.0, 0.0, -2.0, 8.0, 0.0, 0.0, 0.0, 0.0, -3.0, 12.0, 0.0, 6.0, -24.0, 0.0, 0.0, 0.0, 0.0, 2.0, -8.0, 0.0, -4.0, 16.0, 0.0], // dN7/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 16.0, 0.0, 36.0, -48.0, 0.0, -24.0, 32.0, 0.0, 12.0, -16.0, 0.0, -36.0, 48.0, 0.0, 24.0, -32.0, 0.0], // dN8/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -16.0, 0.0, -12.0, 16.0, 0.0, 0.0, 0.0, 0.0, -24.0, 32.0, 0.0, 24.0, -32.0, 0.0], // dN9/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -16.0, 0.0, -24.0, 32.0, 0.0, 0.0, 0.0, 0.0, -12.0, 16.0, 0.0, 24.0, -32.0, 0.0], // dN10/dpsi
                        vec![0.0, 0.0, 0.0, -12.0, 16.0, 0.0, 12.0, -16.0, 0.0, 0.0, 0.0, 0.0, 36.0, -48.0, 0.0, -36.0, 48.0, 0.0, 0.0, 0.0, 0.0, -24.0, 32.0, 0.0, 24.0, -32.0, 0.0],  // dN11/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 16.0, 0.0, 12.0, -48.0, 0.0, -8.0, 32.0, 0.0, 4.0, -16.0, 0.0, -12.0, 48.0, 0.0, 8.0, -32.0, 0.0],  // dN12/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -16.0, 0.0, -4.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 32.0, 0.0, 8.0, -32.0, 0.0],    // dN13/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -16.0, 0.0, -8.0, 32.0, 0.0, 0.0, 0.0, 0.0, -4.0, 16.0, 0.0, 8.0, -32.0, 0.0],    // dN14/dpsi
                        vec![0.0, 0.0, 0.0, -4.0, 16.0, 0.0, 4.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -48.0, 0.0, -12.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 32.0, 0.0, 8.0, -32.0, 0.0],  // dN15/dpsi
                        vec![4.0, -8.0, 0.0, -12.0, 24.0, 0.0, 8.0, -16.0, 0.0, -12.0, 24.0, 0.0, 36.0, -72.0, 0.0, -24.0, 48.0, 0.0, 8.0, -16.0, 0.0, -24.0, 48.0, 0.0, 16.0, -32.0, 0.0], // dN16/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0, 0.0, -24.0, 48.0, 0.0, 16.0, -32.0, 0.0],  // dN17/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0],    // dN18/dpsi
                        vec![0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 8.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -24.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0],  // dN19/dpsi
                        vec![0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, -48.0, 96.0, 0.0, 48.0, -96.0, 0.0, 0.0, 0.0, 0.0, 32.0, -64.0, 0.0, -32.0, 64.0, 0.0],  // dN20/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -64.0, 0.0, -32.0, 64.0, 0.0],    // dN21/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -48.0, 96.0, 0.0, 32.0, -64.0, 0.0, -16.0, 32.0, 0.0, 48.0, -96.0, 0.0, -32.0, 64.0, 0.0],  // dN22/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 32.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -32.0, 64.0, 0.0],    // dN23/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -48.0, 64.0, 0.0, 48.0, -64.0, 0.0, 0.0, 0.0, 0.0, 48.0, -64.0, 0.0, -48.0, 64.0, 0.0],    // dN24/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 64.0, 0.0, 16.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -64.0, 0.0, -16.0, 64.0, 0.0],    // dN25/dpsi
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.0, -128.0, 0.0, -64.0, 128.0, 0.0, 0.0, 0.0, 0.0, -64.0, 128.0, 0.0, 64.0, -128.0, 0.0],    // dN26/dpsi
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            // Add other element types as needed...
            _ => None,
        }
    }

    // Get element dimension
    pub fn get_element_dimension(element_type: &ElementType) -> Option<usize> {
        match element_type {
            ElementType::Line | ElementType::QuadraticEdge => Some(1),
            ElementType::Triangle | ElementType::QuadraticTriangle | 
            ElementType::Quad | ElementType::QuadraticQuad | ElementType::BiquadraticQuad => Some(2),
            ElementType::Tetra | ElementType::QuadraticTetra | ElementType::Pyramid | 
            ElementType::Wedge | ElementType::QuadraticWedge | ElementType::BiquadraticQuadraticWedge | 
            ElementType::Hexahedron | ElementType::QuadraticHexahedron | 
            ElementType::BiquadraticQuadraticHexahedron | ElementType::TriquadraticHexahedron => Some(3),
            ElementType::Vertex => Some(0),
            _ => None,
        }
    }
    // Get element order
    pub fn get_element_order(element_type: &ElementType) -> Option<usize> {
        match element_type {
            ElementType::Line | ElementType::Triangle | ElementType::Quad | ElementType::Tetra |
            ElementType::Pyramid | ElementType::Wedge | ElementType::Hexahedron => Some(1),

            ElementType::QuadraticEdge | ElementType::QuadraticTriangle | ElementType::QuadraticQuad |
            ElementType::BiquadraticQuad | ElementType::QuadraticTetra | ElementType::QuadraticWedge |
            ElementType::BiquadraticQuadraticWedge | ElementType::QuadraticHexahedron |
            ElementType::BiquadraticQuadraticHexahedron | ElementType::TriquadraticHexahedron => Some(2),
            
            // Add higher order elements as needed
            
            ElementType::Vertex => Some(0),
            _ => None,
        }
    }

}


#[derive(Debug)]                                    // Auto-implement Debug for printing
pub struct MeshData {                               // Defines a structure to represent a mesh
    pub dimension: usize,                           // Spatial dimension (from # sdim tag)
    pub num_nodes: usize,                           // Number of nodes (from # number of mesh vertices tag)
    pub min_node_index: usize,                      // Lowest mesh vertex index (from # lowest mesh vertex index tag) 
    pub nodes: Vec<Node>,                           // All nodes with their coordinates
    pub num_eltypes: usize,                         // Number of element types (from # number of element types tag)
    pub elements: Vec<Element>,                     // All elements with their connectivity
    pub element_type_info: Vec<ElementTypeInfo>,    // Information about each element type
}

#[derive(Debug)]
pub struct ElementTypeInfo {
    pub element_type: ElementType,                  // The actual element type enum
    pub num_elements: usize,                        // Number of elements of this type
    pub start_index: usize,                         // Starting index in the main elements vector
    pub nodes_per_element: usize,                   // Nodes per element for this type
}

/// Enumeration of physical value types that can be stored at nodes/elements
/// These represent different physical quantities in simulation results
#[derive(Debug, Clone, PartialEq)]
pub enum ValueType {
    Displacement,   // Structural displacement vector [ux, uy, uz]
    Pressure,       // Scalar pressure field
    Stress,         // Stress tensor (6 components in Voigt notation)
    Strain,         // Strain tensor (6 components in Voigt notation)
    Velocity,       // Velocity vector [vx, vy, vz]
    Acceleration,   // Acceleration vector [ax, ay, az]
    Temperature,    // Scalar temperature field

    // Add more physical quantity types as needed
}

impl ValueType {
    /// Parse value type from text column headers in result files
    /// Converts common abbreviations and full names to ValueType enum
    pub fn from_str_txt(s: &str) -> Option<ValueType> {
        match s.to_lowercase().as_str() {
            // Displacement variants
            "u" | "v" | "w" | "displacement" => Some(ValueType::Displacement),
            // Pressure variants  
            // "p" | "pressure" => Some(ValueType::Pressure),
            // Stress variants
            // "s" | "stress" => Some(ValueType::Stress),
            // Strain variants
            // "strain" => Some(ValueType::Strain),
            // Velocity variants
            "vel" | "velocity" => Some(ValueType::Velocity),
            // Acceleration variants
            "accel" | "acceleration" => Some(ValueType::Acceleration),
            // Temperature variants
            // "t" | "temp" | "temperature" => Some(ValueType::Temperature),

            // Add more variants as needed

            // Unknown type
            _ => None,
        }
    }

    // Helper function to convert ValueType to attribute name
    pub fn get_attribute_name(value_type: &ValueType) -> &'static str {
        match value_type {
            ValueType::Displacement => "Displacement",
            ValueType::Velocity => "Velocity",
            ValueType::Acceleration => "Acceleration",
            // ValueType::Stress => "Stress",
            // ValueType::Strain => "Strain",
            // ValueType::Temperature => "Temperature",
            // ValueType::Pressure => "Pressure",

            // Add more mappings as needed for your ValueType enum

            _ => "Unknown",
        }
    }
}

#[derive(Debug)]
pub struct NodeValue {                  // Defines a structure to represent node value
    pub id: usize,                      // Unique identifier for the node
    pub values: Vec<f64>,               // node values (x, y for 2D, x,y,z for 3D)
}
 
#[derive(Debug)]
pub struct NodeValueTypeInfo {
    pub dimension: usize,               // Physical dimension (1=scalar, 3=vector, 6=tensor)
    pub num_nodes: usize,               // Total number of nodes
    pub nodevalue_type: ValueType,      // Type of physical quantity
    pub num_nodevalue_type: usize,      // Total number of node values type
    pub start_index: usize,             // Starting index in the node_values array
}

#[derive(Debug)]
pub struct NodeValueData {
    pub node_values: Vec<NodeValue>,                          // All node values: [type1_allnodes_alltime, type2_allnodes_alltime, ...]
    pub node_value_type_info: Vec<NodeValueTypeInfo>,   // Metadata describing node value organization
}


/// Geometric analysis module

#[derive(Debug, Clone)]
pub struct ElementQuality {
    pub element_id: usize,          // ID of the element being analyzed
    pub shape_metric: f64,          // shape metric value
    pub skewness_metric: f64,       // skewness metric value
    pub length_ratio: f64,          // length ratio metric value
    pub orientation_metric: f64,    // orientation metric value
    pub volume_metric: f64,         // volume metric
    pub jacobian_ratio: f64,        // max(detJ) / min(detJ)

    // more quality metrics can be added here
}

// Structure to hold the complete mesh quality analysis results
#[derive(Debug, Clone)]
pub struct MeshQualityReport {
    pub total_elements: usize,                    // Total number of elements that were successfully analyzed
    pub element_qualities: Vec<ElementQuality>,   // Quality metrics for each individual element
    // pub statistics: QualityStatistics,           // Overall statistical summary of mesh quality
}

/// Type of finite element integration to perform
#[derive(Debug, Clone, Copy)]
pub enum IntegrationType {
    Mass,      // Mass matrix integration: ∫ ρ N_i N_j dΩ
    Stiffness, // Stiffness matrix integration
}

#[derive(Debug, Clone)]
pub enum MaterialProperty {
    Scalar(Vec<f64>),           // For mass matrix (density)
    Matrix(Vec<Vec<Vec<f64>>>), // For stiffness matrix (anisotropic material)
}

impl MaterialProperty {
    pub fn as_scalar(&self) -> Result<&[f64], GaussError> {
        match self {
            MaterialProperty::Scalar(coeffs) => Ok(coeffs),
            _ => Err(GaussError::InvalidMaterialProperty(
                "Expected scalar material property".to_string(),
            )),
        }
    }

    pub fn as_matrix(&self) -> Result<&Vec<Vec<Vec<f64>>>, GaussError> {
        match self {
            MaterialProperty::Matrix(matrix) => Ok(matrix),
            _ => Err(GaussError::InvalidMaterialProperty(
                "Expected matrix material property".to_string(),
            )),
        }
    }
}


#[derive(Debug, Clone)]
pub struct GaussianPointNumber {
    pub element_id: usize,              // ID of the element being analyzed
    pub theoretical_number: usize,      // Theoretical number of Gauss points for the element
    pub optimal_number: usize,          // Optimal number of Gauss points used in the element
}

#[derive(Debug, Clone)]
pub struct GaussianPointNumberReport {
    pub total_elements: usize,                    // Total number of elements that were successfully analyzed
    pub gauss_point_numbers: Vec<GaussianPointNumber>,   // Gaussian point numbers for each individual element
}


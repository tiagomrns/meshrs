use std::collections::HashMap;

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
    pub fn new(coefficients: Vec<f64>) -> Self {
        Self { coefficients }
    }

    // Infer max_degree from coefficient vector length (always 3D)
    // Uses formula: length = (max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6
    fn infer_max_degree(len: usize) -> Result<u32, &'static str> {
        // For 3D: length = (max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6
        for max_degree in 0..=100 {
            let expected = ((max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6) as usize;
            if expected == len {
                return Ok(max_degree);
            }
            if expected > len {
                break;
            }
        }
        
        Err("Coefficient vector length does not correspond to a valid 3D polynomial degree")
    }

    // Calculate expected length for given max_degree (always 3D)
    pub fn expected_length(max_degree: u32) -> usize {
        ((max_degree + 1) * (max_degree + 2) * (max_degree + 3) / 6) as usize
    }

    // Generate graded lexicographic basis for 3D with given maximum degree
    // Order: [1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, ...]
    // Sorted by total degree, then by x, then by y
    fn generate_basis(max_degree: u32) -> Vec<(u32, u32, u32)> {
        let mut basis = Vec::new();
        
        // Iterate through all possible total degrees
        for total_degree in 0..=max_degree {
            // Generate all combinations (i, j, k) such that i + j + k = total_degree
            for i in (0..=total_degree).rev() {
                for j in (0..=(total_degree - i)).rev() {
                    let k = total_degree - i - j;
                    basis.push((i, j, k));
                }
            }
        }   
        
        basis
    }

    // Map exponent to index: (i+j+k)*(i+j+k+1)*(i+j+k+2)/6 + (j+k)*(j+k+1)/2 + k
    // This is the combinatorial formula for graded lexicographic ordering in 3D
    fn map_index(exponent: (u32, u32, u32)) -> usize {
        let (i, j, k) = exponent;
        let n = i + j + k;
        ((n * (n + 1) * (n + 2)) / 6 + ((j + k) * (j + k + 1)) / 2 + k) as usize
    }

    // Add two polynomials
    pub fn add(first: &[f64], second: &[f64]) -> Result<Vec<f64>, &'static str> {
    let max_len = first.len().max(second.len());
    let result = (0..max_len)
        .map(|i| {
            first.get(i).unwrap_or(&0.0) + second.get(i).unwrap_or(&0.0)
        })
        .collect();
    
    Ok(result)
}

    // Multiply polynomial by scalar
    pub fn multiply_scalar(coeffs: &[f64], scalar: f64) -> Vec<f64> {
        coeffs.iter()
            .map(|c| c * scalar)
            .collect()
    }

    // Multiply two polynomials using distributive property
    pub fn multiply(first: &[f64], second: &[f64]) -> Result<Vec<f64>, &'static str> {
        let max_deg_first = Self::infer_max_degree(first.len())?;
        let max_deg_second = Self::infer_max_degree(second.len())?;

        // Generate basis for both polynomials
        let basis_first = Self::generate_basis(max_deg_first);
        let basis_second = Self::generate_basis(max_deg_second);

        // The resulting max degree is the sum of the two max degrees
        let result_max_degree = max_deg_first + max_deg_second;
        let result_len = Self::expected_length(result_max_degree);
        
        let mut result = vec![0.0; result_len];
        
        // Multiply each term from first with each term from second
        for (i, &(exp_i_x, exp_i_y, exp_i_z)) in basis_first.iter().enumerate() {
            for (j, &(exp_j_x, exp_j_y, exp_j_z)) in basis_second.iter().enumerate() {
                let coeff_product = first[i] * second[j];
                
                if coeff_product.abs() < 1e-15 {
                    continue; // Skip negligible terms for numerical stability
                }
                
                // Add exponents to get the resulting term
                let result_exponent = (
                    exp_i_x + exp_j_x,
                    exp_i_y + exp_j_y,
                    exp_i_z + exp_j_z,
                );
                
                // Find the index in the result basis
                let result_index = Self::map_index(result_exponent);
                result[result_index] += coeff_product;
            }
        }
        
        Ok(result)
    }

    // Evaluate polynomial at a point (x, y, z)
    pub fn evaluate(coeffs: &[f64], point: (f64, f64, f64)) -> Result<f64, &'static str> {
        let max_degree = Self::infer_max_degree(coeffs.len())?;
        let (x, y, z) = point;
        let basis = Self::generate_basis(max_degree);

        // Sum over all terms: coefficient * x^i * y^j * z^k
        Ok(coeffs.iter()
            .zip(&basis)
            .map(|(coeff, &(i, j, k))| {
                coeff * x.powi(i as i32) * y.powi(j as i32) * z.powi(k as i32)
            })
            .sum())
    }

    // Helper to pad polynomial to a higher degree (fills with zeros)
    pub fn pad_to_degree(coeffs: &[f64], target_degree: u32) -> Result<Vec<f64>, &'static str> {
        let current_degree = Self::infer_max_degree(coeffs.len())?;
        
        if target_degree < current_degree {
            return Err("Target degree must be >= current degree");
        }
        
        if target_degree == current_degree {
            return Ok(coeffs.to_vec());
        }
        
        let target_len = Self::expected_length(target_degree);
        let mut result = coeffs.to_vec();
        result.resize(target_len, 0.0);
        
        Ok(result)
    }

    /// Get 1D coefficients (only x powers): [a₀, a₁, a₂, ...] for a₀ + a₁x + a₂x² + ...
    /// Extracts coefficients where y and z exponents are zero
    pub fn get_coefficients_1d(coeffs: &[f64]) -> Result<Vec<f64>, &'static str> {
        let degree = Self::infer_max_degree(coeffs.len())?;
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
        let degree = Self::infer_max_degree(coeffs.len())?;
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
        let degree = Self::infer_max_degree(coeffs.len())?;
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
        let max_degree = match Self::infer_max_degree(polynomial.len()) {
            Ok(deg) => deg,
            Err(_) => return 0,
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
/// Represented as polynomials in the natural coordinate system
#[derive(Debug, Clone)]
pub struct ShapeFunction {
    pub values: Vec<Vec<f64>>, // Each inner vec is coefficients for one shape function
    pub derivatives: Vec<Vec<Vec<f64>>>, // [node][direction][coefficients]
    pub num_nodes: usize,
}

/// Jacobian matrix mapping from natural to physical coordinates
/// J[i][j] = ∂x_i/∂ξ_j where x are physical coords, ξ are natural coords
#[derive(Debug, Clone)]
pub struct Jacobian {
    pub matrix: Vec<Vec<Vec<f64>>>, // Jacobian matrix J[i][j] = dx_i/dxi_j
    pub determinant: Vec<f64>, // Determinant as polynomial (for curved elements)
}

impl Jacobian {
    pub fn evaluate_determinant_at_point(&self, point: (f64, f64, f64)) -> Result<f64, &'static str> {
    MonomialPolynomial::evaluate(&self.determinant, point)
}
}

/// Element types available in VTK 
#[derive(Debug, Clone)]
pub enum ElementType {              // Enumeration of supported finite element types
    Vertex,                         // Vertex element   Comsol
    Line,                           // First order edge element   Comsol
    QuadraticEdge,                  // Second order edge element  Comsol
    Triangle,                       // First order triangular element   Comsol
    QuadraticTriangle,              // Second order triangular element   Comsol
    Quad,                           // First order quadrilateral element   Comsol
    QuadraticQuad,                  // Second order quadrilateral element 9 nodes  Abaqus
    BiquadraticQuad,                // Second order quadrilateral element 9 nodes  Comsol
    Tetra,                          // First order tetrahedral element   Comsol
    QuadraticTetra,                 // Second order tetrahedral element   Comsol
    Pyramid,                        // First order pyramid element   Comsol
    //QuadraticPyramid,            // Second order pyramid element   13 nodes   Comsol  ???
    // TriquadraticPyramid,            // Second order pyramid element   with padded nodes  19 nodes   Comsol
    Wedge,                          // First order prism element   Comsol
    QuadraticWedge,                 // Second order prism element   Abaqus
    BiquadraticQuadraticWedge,      // Second order prism element   Comsol
    Hexahedron,                     // First order hexahedral element   Comsol
    QuadraticHexahedron,            // Second order hexahedral element  20 nodes   Abaqus
    BiquadraticQuadraticHexahedron, // Second order hexahedral element  24 nodes   Abaqus
    TriquadraticHexahedron,         // Second order hexahedral element  27 nodes   Comsol

    //QuadraticHexahedron27,        // Second order hexahedral element  e.g. only available in comsol and vtk                     
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
            //ElementType::QuadraticPyramid => CellType::QuadraticPyramid,
            // ElementType::TriquadraticPyramid => CellType::TriquadraticPyramid,   
            ElementType::Wedge => CellType::Wedge,   
            ElementType::QuadraticWedge => CellType::QuadraticWedge,
            ElementType::BiquadraticQuadraticWedge => CellType::BiquadraticQuadraticWedge,  
            ElementType::Hexahedron => CellType::Hexahedron,   
            ElementType::QuadraticHexahedron => CellType::QuadraticHexahedron,  
            ElementType::BiquadraticQuadraticHexahedron => CellType::BiquadraticQuadraticHexahedron, 
            ElementType::TriquadraticHexahedron => CellType::TriquadraticHexahedron,  
            
        }
    }

    /// Monomial shape functions and their derivatives for standard elements
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
                // 1D: N1 = 1 - x, N2 = x
                // Padding to 3D degree 1: [1, x, y, z]
                
                let values = vec![
                    vec![1.0, -1.0, 0.0, 0.0],  // N1 = 1 - x
                    vec![0.0, 1.0, 0.0, 0.0],   // N2 = x
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN1/dx = -1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/dz = 0
                    ],
                    vec![
                        vec![1.0, 0.0, 0.0, 0.0],   // dN2/dx = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dz = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticEdge => {
                let num_nodes = 3;
                // 1D degree 2: N1 = 1 - 3x + 2x², N2 = -x + 2x², N3 = 4x - 4x²
                // Padding to 3D degree 2: [1, x, y, z, x², xy, xz, y², yz, z²]
                
                let values = vec![
                    vec![1.0, -3.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N1 = 1 - 3x + 2x²
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N2 = -x + 2x²  
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N3 = 4x - 4x²
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dx = -3 + 4x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dz = 0
                    ],
                    vec![
                        vec![-1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dx = -1 + 4x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dz = 0
                    ],
                    vec![
                        vec![4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dx = 4 - 8x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dz = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Triangle => {
                let num_nodes = 3;
                // 2D degree 1: N1 = 1 - x - y, N2 = x, N3 = y
                // Padding to 3D degree 1: [1, x, y, z]
                
                let values = vec![
                    vec![1.0, -1.0, -1.0, 0.0],  // N1 = 1 - x - y
                    vec![0.0, 1.0, 0.0, 0.0],    // N2 = x
                    vec![0.0, 0.0, 1.0, 0.0],    // N3 = y
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN1/dx = -1
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN1/dy = -1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN1/dz = 0
                    ],
                    vec![
                        vec![1.0, 0.0, 0.0, 0.0],   // dN2/dx = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0],   // dN3/dx = 0
                        vec![1.0, 0.0, 0.0, 0.0],   // dN3/dy = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN3/dz = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticTriangle => {
                let num_nodes = 6;
                // 2D degree 2, padding to 3D degree 2: [1, x, y, z, x², xy, xz, y², yz, z²]
                
                let values = vec![
                    vec![1.0, -3.0, -3.0, 0.0, 2.0, 4.0, 0.0, 2.0, 0.0, 0.0],  // N1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // N2
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0],   // N3
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0],  // N4
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0],    // N5
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0],  // N6
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dx
                        vec![0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dz
                    ],
                    vec![
                        vec![-1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dx
                        vec![0.0, -1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dz
                    ],
                    vec![
                        vec![4.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dx
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dx
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dz
                    ],
                    vec![
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN6/dx
                        vec![4.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN6/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dz
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Quad => {
                let num_nodes = 4;
                // 2D degree 1 with bilinear term, padding to 3D degree 2
                // [1, x, y, z, x², xy, xz, y², yz, z²]

                let values = vec![
                    vec![1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  // N1 = (1-x)(1-y) = 1 - x - y + xy
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],   // N2 = x(1-y) = x - xy
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],    // N3 = xy
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],   // N4 = (1-x)y = y - xy
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dx = -1 + y
                        vec![-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dy = -1 + x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dz = 0
                    ],
                    vec![
                        vec![1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dx = 1 - y
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dy = -x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dx = y
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dy = x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dx = -y
                        vec![1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dy = 1 - x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dz = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },


            ElementType::QuadraticQuad => {
                let num_nodes = 8;
                // 2D degree 2 serendipity element
                // Basis: [1, x, y, z, x², xy, xz, y², yz, z²]
                
                let values = vec![
                    vec![1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],        // N0 = 1 - x - y + xy
                    vec![0.0, 1.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0],        // N1 = x - x² - xy  
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, -1.0, 0.0, 0.0],         // N2 = x² + xy - y²
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0],        // N3 = y - xy - y²
                    vec![0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0],         // N4 = 4x² - 4xy
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, -4.0, 0.0, 0.0],         // N5 = 4xy - 4y²
                    vec![0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0],         // N6 = -4x² + 4xy
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 4.0, 0.0, 0.0],         // N7 = -4xy + 4y²
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN0/dx = -1 + y
                        vec![-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN0/dy = -1 + x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN0/dz = 0
                    ],
                    vec![
                        vec![1.0, -2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dx = 1 - 2x - y
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN1/dy = -x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN1/dz = 0
                    ],
                    vec![
                        vec![0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN2/dx = 2x + y
                        vec![0.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN2/dy = x - 2y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN2/dz = 0
                    ],
                    vec![
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN3/dx = -y
                        vec![1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dy = 1 - x - 2y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN3/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 8.0, -4.0, 0.0, 0.0, 0.0, 0.0],    // dN4/dx = 8x - 4y
                        vec![0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN4/dy = -4x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dx = 4y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, -8.0, 0.0, 0.0],    // dN5/dy = 4x - 8y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -8.0, 4.0, 0.0, 0.0, 0.0, 0.0],    // dN6/dx = -8x + 4y
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6/dy = 4x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0],    // dN7/dx = -4y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 8.0, 0.0, 0.0],    // dN7/dy = -4x + 8y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/dz = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::BiquadraticQuad => {
                let num_nodes = 9;
                // 2D biquadratic element with center node
                // Basis: [1, x, y, z, x², xy, xz, y², yz, z²]
                
                let values = vec![
                    vec![1.0, -1.5, -1.5, 0.0, 0.5, 2.25, 0.0, 0.5, 0.0, 0.0],      // N0 = 1 - 1.5x - 1.5y + 0.5x² + 2.25xy + 0.5y²
                    vec![0.0, 2.0, 0.0, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0],       // N1 = 2x - 2x² - 2xy
                    vec![0.0, -0.5, 0.0, 0.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0],       // N2 = -0.5x + 0.5x² + 0.25xy
                    vec![0.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0, -2.0, 0.0, 0.0],       // N3 = 2y - 2xy - 2y²
                    vec![0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],         // N4 = 4x²
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0],         // N5 = 4xy
                    vec![0.0, 0.0, -0.5, 0.0, 0.0, 0.25, 0.0, 0.5, 0.0, 0.0],       // N6 = -0.5y + 0.25xy + 0.5y²
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0],         // N7 = 4y²
                    vec![0.0, 0.0, 0.0, 0.0, -4.0, -4.0, 0.0, -4.0, 0.0, 0.0],      // N8 = -4x² - 4xy - 4y²
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.5, 1.0, 2.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dx = -1.5 + x + 2.25y
                        vec![-1.5, 2.25, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN0/dy = -1.5 + 2.25x + y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN0/dz = 0
                    ],
                    vec![
                        vec![2.0, -4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN1/dx = 2 - 4x - 2y
                        vec![0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN1/dy = -2x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN1/dz = 0
                    ],
                    vec![
                        vec![-0.5, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dx = -0.5 + x + 0.25y
                        vec![0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN2/dy = 0.25x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN2/dz = 0
                    ],
                    vec![
                        vec![0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN3/dx = -2y
                        vec![2.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dy = 2 - 2x - 4y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN3/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dx = 8x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN4/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dx = 4y
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dy = 4x
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN5/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // dN6/dx = 0.25y
                        vec![-0.5, 0.25, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dy = -0.5 + 0.25x + y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN6/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/dx = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0],     // dN7/dy = 8y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN7/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0],   // dN8/dx = -8x - 4y
                        vec![0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -8.0, 0.0, 0.0],   // dN8/dy = -4x - 8y
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // dN8/dz = 0
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Tetra => {
                let num_nodes = 4;
                // 3D degree 1: [1, x, y, z]

                let values = vec![
                    vec![1.0, -1.0, -1.0, -1.0],  // N1 = 1 - x - y - z
                    vec![0.0, 1.0, 0.0, 0.0],     // N2 = x
                    vec![0.0, 0.0, 1.0, 0.0],     // N3 = y
                    vec![0.0, 0.0, 0.0, 1.0],     // N4 = z
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN1/dx = -1
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN1/dy = -1
                        vec![-1.0, 0.0, 0.0, 0.0],  // dN1/dz = -1
                    ],
                    vec![
                        vec![1.0, 0.0, 0.0, 0.0],   // dN2/dx = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dy = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN2/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0],   // dN3/dx = 0
                        vec![1.0, 0.0, 0.0, 0.0],   // dN3/dy = 1
                        vec![0.0, 0.0, 0.0, 0.0],   // dN3/dz = 0
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0],   // dN4/dx = 0
                        vec![0.0, 0.0, 0.0, 0.0],   // dN4/dy = 0
                        vec![1.0, 0.0, 0.0, 0.0],   // dN4/dz = 1
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticTetra => {
                let num_nodes = 10;
                // 3D degree 2: [1, x, y, z, x², xy, xz, y², yz, z²]

                let values = vec![
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 4.0, 4.0, 2.0, 4.0, 2.0],  // N1
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // N2
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0],    // N3
                    vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0],    // N4
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, -4.0, 0.0, 0.0, 0.0],  // N5
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0],     // N6
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0, 0.0],  // N7
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, -4.0],   // N8
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0],     // N9
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0],     // N10
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dx
                        vec![0.0, -3.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dy
                        vec![0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dz
                    ],
                    vec![
                        vec![-1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dx
                        vec![0.0, -1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN4/dy
                        vec![0.0, 0.0, -1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dz
                    ],
                    vec![
                        vec![4.0, -8.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dx
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dy
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dx
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dz
                    ],
                    vec![
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN7/dx
                        vec![4.0, 0.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN7/dy
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN7/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8/dx
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8/dy
                        vec![4.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN8/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dy
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN9/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN10/dx
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN10/dy
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN10/dz
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Pyramid => {
                let num_nodes = 5;
                // 3D degree 1 with trilinear term, padding to 3D degree 2
                // [1, x, y, z, x², xy, xz, y², yz, z²]

                let values = vec![
                    vec![1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  // N1 = (1-x)(1-y)(1-z)
                    vec![0.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0],   // N2 = x(1-y)(1-z)
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],     // N3 = xy(1-z)
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0],   // N4 = (1-x)y(1-z)
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],     // N5 = z
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dx = -(1-y)(1-z) = -1 + y + z - yz
                        vec![0.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dy = -(1-x)(1-z) = -1 + x + z - xz
                        vec![0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dz = -(1-x)(1-y) = -1 + x + y - xy
                    ],
                    vec![
                        vec![1.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dx = (1-y)(1-z) = 1 - y - z + yz
                        vec![0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dy = -x(1-z) = -x + xz
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dz = -x(1-y) = -x + xy
                    ],
                    vec![
                        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dx = y(1-z) = y - yz
                        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dy = x(1-z) = x - xz
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dz = -xy
                    ],
                    vec![
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dx = -y(1-z) = -y + yz
                        vec![1.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dy = (1-x)(1-z) = 1 - x - z + xz
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dz = -(1-x)y = -y + xy
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dx = 0
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dy = 0
                        vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dz = 1
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Wedge => {
                let num_nodes = 6;
                // 3D degree 1 with bilinear terms, padding to 3D degree 2
                // [1, x, y, z, x², xy, xz, y², yz, z²]
                
                let values = vec![
                    vec![1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  // N1 = (1-x-y)(1-z)
                    vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],   // N2 = x(1-z)
                    vec![0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],   // N3 = y(1-z)
                    vec![0.0, -1.0, -1.0, 1.0, 0.0, 1.0, -1.0, 0.0, -1.0, 0.0], // N4 = (1-x-y)z
                    vec![0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],     // N5 = xz
                    vec![0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],     // N6 = yz
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dx = -(1-z) = -1 + z
                        vec![0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dy = -(1-z) = -1 + z
                        vec![0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN1/dz = -(1-x-y) = -1 + x + y
                    ],
                    vec![
                        vec![1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dx = 1 - z
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN2/dy = 0
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN2/dz = -x
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN3/dx = 0
                        vec![1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dy = 1 - z
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN3/dz = -y
                    ],
                    vec![
                        vec![-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dx = -z
                        vec![0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dy = -z
                        vec![0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  // dN4/dz = 1 - x - y
                    ],
                    vec![
                        vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dx = z
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dy = 0
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN5/dz = x
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dx = 0
                        vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dy = z
                        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   // dN6/dz = y
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::QuadraticWedge => {
                
                let num_nodes = 15;

                let values: Vec<Vec<f64>> = vec![
                    vec![1.0,-3.0,-3.0,-3.0, 2.0, 4.0, 5.0, 2.0, 5.0, 2.0, 0.0, 0.0,-2.0, 0.0,-4.0,-2.0, 0.0,-2.0,-2.0, 0.0], // N0
                    vec![0.0,-1.0, 0.0, 0.0, 2.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0], // N1
                    vec![0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 2.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-2.0, 2.0, 0.0], // N2
                    vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0,-1.0, 0.0,-1.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0,-2.0, 0.0, 2.0,-2.0, 0.0], // N3
                    vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0,-2.0, 0.0, 0.0, 2.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0], // N4
                    vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0,-2.0, 0.0], // N5
                    vec![0.0, 2.0, 0.0,-4.0,-4.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0], // N6
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N7
                    vec![0.0,-2.0, 2.0, 0.0, 0.0,-4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N8
                    vec![0.0, 2.0, 0.0, 4.0,-4.0, 0.0,-2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0], // N9
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N10
                    vec![0.0,-2.0,-2.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N11
                    vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0,-4.0, 0.0,-4.0,-4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 4.0, 4.0, 0.0], // N12
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0], // N13
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0,-4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 0.0], // N14
                ];

                let derivatives: Vec<Vec<Vec<f64>>> = vec![
                    vec![
                        vec![-3.0, 4.0, 4.0, 5.0, 0.0, 0.0,-4.0, 0.0,-4.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dx
                        vec![-3.0, 4.0, 4.0, 5.0, 0.0, 0.0, 0.0, 0.0,-4.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dy
                        vec![-3.0, 5.0, 5.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dz
                    ],
                    vec![
                        vec![-1.0, 4.0, 0.0, 0.0, 4.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0,-6.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0], // dN1/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dy
                        vec![0.0,-2.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dx
                        vec![-1.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-6.0, 6.0, 0.0], // dN2/dy
                        vec![0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0], // dN2/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0,-1.0, 0.0, 0.0,-2.0, 0.0,-2.0, 4.0, 0.0, 0.0, 6.0, 0.0, 2.0,-6.0, 0.0, 2.0,-6.0, 0.0], // dN3/dx
                        vec![0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0,-2.0, 4.0, 0.0, 0.0, 2.0, 0.0, 2.0,-2.0, 0.0, 6.0,-6.0, 0.0], // dN3/dy
                        vec![1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 0.0, 2.0,-6.0, 0.0, 0.0,-2.0, 0.0,-2.0, 6.0, 0.0,-6.0, 6.0, 0.0], // dN3/dz
                    ],
                    vec![
                        vec![1.0, 4.0, 0.0, 0.0, 4.0, 0.0, 2.0, 0.0, 0.0,-6.0, 0.0, 0.0, 6.0, 0.0, 0.0,-6.0, 0.0, 0.0, 0.0], // dN4/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0], // dN4/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dx
                        vec![1.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0,-6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0,-6.0, 0.0], // dN5/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0,-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0], // dN5/dz
                    ],
                    vec![
                        vec![2.0,-8.0, 0.0, 4.0,-8.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 12.0, 0.0, 0.0,-12.0, 0.0, 0.0, 0.0], // dN6/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dy
                        vec![0.0, 0.0, 0.0,-4.0, 0.0, 0.0,-4.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0], // dN6/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dx
                        vec![0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dz
                    ],
                    vec![
                        vec![-2.0, 0.0,-4.0, 0.0, 0.0,-8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dx
                        vec![2.0,-4.0, 0.0, 0.0, 0.0,-8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dz
                    ],
                    vec![
                        vec![2.0,-8.0, 0.0,-4.0,-8.0, 0.0,-4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 12.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0], // dN9/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dy
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0, 0.0, 0.0,-8.0, 0.0, 0.0, 0.0], // dN9/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dx
                        vec![0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dz
                    ],
                    vec![
                        vec![-2.0, 0.0,-4.0, 0.0, 0.0,-8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dx
                        vec![-2.0, 0.0,-4.0, 0.0, 0.0,-8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0,-4.0, 0.0, 0.0,-8.0, 0.0,-8.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 8.0, 8.0, 0.0], // dN12/dx
                        vec![0.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0, 0.0,-8.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 0.0], // dN12/dy
                        vec![2.0, 0.0, 0.0, 4.0, 0.0, 0.0, 8.0, 0.0, 8.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0,-8.0, 0.0,-8.0,-8.0, 0.0], // dN12/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0,-8.0, 0.0, 0.0, 0.0], // dN13/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 0.0, 0.0, 0.0], // dN13/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-8.0, 0.0], // dN14/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-8.0, 0.0], // dN14/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 0.0], // dN14/dz
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })

            },

            ElementType::BiquadraticQuadraticWedge => {
                let num_nodes = 18;

                // Polynomial basis: [1, x, y, z, x², xy, xz, y², yz, z²]
                let values = vec![
                    vec![0.487654320988, -1.0, -1.0, -0.611111111111, 0.444444444444, 0.888888888889, 0.666666666667, 0.444444444444, 0.666666666667, 0.111111111111],
                    vec![-0.0679012345679, 0.111111111111, -5.27713734955e-16, 0.0555555555556, 0.444444444444, 1.15719499557e-15, -0.666666666667, 4.38088614426e-16, -4.01691975728e-16, 0.111111111111],
                    vec![-0.0679012345679, -4.4408920985e-16, 0.111111111111, 0.0555555555556, 4.94473784955e-16, 9.51211520478e-16, 4.61619375175e-16, 0.444444444444, -0.666666666667, 0.111111111111],
                    vec![-0.0123456790123, -0.333333333333, -0.333333333333, 0.388888888889, 0.444444444444, 0.888888888889, -0.666666666667, 0.444444444444, -0.666666666667, 0.111111111111],
                    vec![0.0987654320988, -0.555555555556, 1.16029566851e-15, -0.277777777778, 0.444444444444, -2.50469749212e-15, 0.666666666667, -7.71778206233e-16, 5.43726392727e-16, 0.111111111111],
                    vec![0.0987654320988, 1.12410081243e-15, -0.555555555556, -0.277777777778, -7.08887847936e-16, -2.1437535206e-15, -3.89301814306e-16, 0.444444444444, 0.666666666667, 0.111111111111],
                    vec![0.327160493827, 0.888888888889, -0.333333333333, -1.05555555556, -0.888888888889, -0.888888888889, -4.46313478191e-16, -1.13353847745e-15, 0.666666666667, 0.555555555556],
                    vec![-0.00617283950617, 0.333333333333, 0.333333333333, -0.388888888889, 6.65968160308e-16, 0.888888888889, -0.666666666667, 4.06435840368e-16, -0.666666666667, 0.555555555556],
                    vec![0.327160493827, -0.333333333333, 0.888888888889, -1.05555555556, -2.12841598697e-15, -0.888888888889, 0.666666666667, -0.888888888889, 7.16424909104e-16, 0.555555555556],
                    vec![-0.172839506173, 0.888888888889, 0.333333333333, -0.0555555555556, -0.888888888889, -0.888888888889, 6.49577903159e-16, 1.33180141384e-15, -0.666666666667, 0.555555555556],
                    vec![0.16049382716, -0.333333333333, -0.333333333333, -0.722222222222, -2.22508824011e-15, 0.888888888889, 0.666666666667, -1.7197428204e-15, 0.666666666667, 0.555555555556],
                    vec![-0.172839506173, 0.333333333333, 0.888888888889, -0.0555555555556, 1.70374972981e-15, -0.888888888889, -0.666666666667, -0.888888888889, 2.11480270434e-16, 0.555555555556],
                    vec![0.524691358025, -1.66666666667, -1.66666666667, 0.222222222222, 1.11111111111, 2.22222222222, -7.39552089594e-16, 1.11111111111, -8.50574392057e-16, -0.222222222222],
                    vec![-0.0308641975309, -0.555555555556, 1.09004825172e-16, 0.222222222222, 1.11111111111, 5.64281736051e-16, -5.561676196e-17, 3.92197763271e-16, -2.08272427846e-16, -0.222222222222],
                    vec![-0.0308641975309, 8.881784197e-16, -0.555555555556, 0.222222222222, 6.7733316273e-17, -1.87266424387e-16, -9.10477155061e-16, 1.11111111111, -4.31693475692e-16, -0.222222222222],
                    vec![-0.154320987654, 2.22222222222, -4.98050535729e-15, 1.11111111111, -2.22222222222, -2.22222222222, -7.37135775276e-16, 2.99361821225e-15, -9.0366922897e-16, -1.11111111111],
                    vec![-0.154320987654, 0.0, -8.52024511066e-16, 1.11111111111, 8.77960459904e-16, 2.22222222222, -9.54249470128e-16, 6.81908874798e-16, -1.20404965067e-15, -1.11111111111],
                    vec![-0.154320987654, -3.60822483003e-15, 2.22222222222, 1.11111111111, 2.14156515294e-15, -2.22222222222, 4.46046845223e-16, -2.22222222222, -9.83365298982e-16, -1.11111111111],
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0, 0.888888888889, 0.888888888889, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-1.0, 0.888888888889, 0.888888888889, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.611111111111, 0.666666666667, 0.666666666667, 0.222222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.111111111111, 0.888888888889, 1.15719499557e-15, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-5.27713734955e-16, 1.15719499557e-15, 8.76177228852e-16, -4.01691975728e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0555555555556, -0.666666666667, -4.01691975728e-16, 0.222222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-4.4408920985e-16, 9.88947569911e-16, 9.51211520478e-16, 4.61619375175e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.111111111111, 9.51211520478e-16, 0.888888888889, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0555555555556, 4.61619375175e-16, -0.666666666667, 0.222222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-0.333333333333, 0.888888888889, 0.888888888889, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.333333333333, 0.888888888889, 0.888888888889, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.388888888889, -0.666666666667, -0.666666666667, 0.222222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-0.555555555556, 0.888888888889, -2.50469749212e-15, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![1.16029566851e-15, -2.50469749212e-15, -1.54355641247e-15, 5.43726392727e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.277777777778, 0.666666666667, 5.43726392727e-16, 0.222222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![1.12410081243e-15, -1.41777569587e-15, -2.1437535206e-15, -3.89301814306e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.555555555556, -2.1437535206e-15, 0.888888888889, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.277777777778, -3.89301814306e-16, 0.666666666667, 0.222222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.888888888889, -1.77777777778, -0.888888888889, -4.46313478191e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.333333333333, -0.888888888889, -2.2670769549e-15, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-1.05555555556, -4.46313478191e-16, 0.666666666667, 1.11111111111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.333333333333, 1.33193632062e-15, 0.888888888889, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.333333333333, 0.888888888889, 8.12871680736e-16, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.388888888889, -0.666666666667, -0.666666666667, 1.11111111111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-0.333333333333, -4.25683197395e-15, -0.888888888889, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.888888888889, -0.888888888889, -1.77777777778, 7.16424909104e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-1.05555555556, 0.666666666667, 7.16424909104e-16, 1.11111111111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.888888888889, -1.77777777778, -0.888888888889, 6.49577903159e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.333333333333, -0.888888888889, 2.66360282769e-15, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.0555555555556, 6.49577903159e-16, -0.666666666667, 1.11111111111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-0.333333333333, -4.45017648023e-15, 0.888888888889, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.333333333333, 0.888888888889, -3.43948564079e-15, 0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.722222222222, 0.666666666667, 0.666666666667, 1.11111111111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.333333333333, 3.40749945962e-15, -0.888888888889, -0.666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.888888888889, -0.888888888889, -1.77777777778, 2.11480270434e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.0555555555556, -0.666666666667, 2.11480270434e-16, 1.11111111111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-1.66666666667, 2.22222222222, 2.22222222222, -7.39552089594e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-1.66666666667, 2.22222222222, 2.22222222222, -8.50574392057e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.222222222222, -7.39552089594e-16, -8.50574392057e-16, -0.444444444444, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-0.555555555556, 2.22222222222, 5.64281736051e-16, -5.561676196e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![1.09004825172e-16, 5.64281736051e-16, 7.84395526542e-16, -2.08272427846e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.222222222222, -5.561676196e-17, -2.08272427846e-16, -0.444444444444, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![8.881784197e-16, 1.35466632546e-16, -1.87266424387e-16, -9.10477155061e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.555555555556, -1.87266424387e-16, 2.22222222222, -4.31693475692e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.222222222222, -9.10477155061e-16, -4.31693475692e-16, -0.444444444444, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![2.22222222222, -4.44444444444, -2.22222222222, -7.37135775276e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-4.98050535729e-15, -2.22222222222, 5.9872364245e-15, -9.0366922897e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![1.11111111111, -7.37135775276e-16, -9.0366922897e-16, -2.22222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 1.75592091981e-15, 2.22222222222, -9.54249470128e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-8.52024511066e-16, 2.22222222222, 1.3638177496e-15, -1.20404965067e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![1.11111111111, -9.54249470128e-16, -1.20404965067e-15, -2.22222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![-3.60822483003e-15, 4.28313030588e-15, -2.22222222222, 4.46046845223e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![2.22222222222, -2.22222222222, -4.44444444444, -9.83365298982e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![1.11111111111, 4.46046845223e-16, -9.83365298982e-16, -2.22222222222, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::Hexahedron => {
                let num_nodes = 8;

                // Polynomial basis: [1, x, y, z, x², xy, xz, y², yz, z²]
                let values= vec![
                    vec![ 1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0],  // N0 = (1-x)(1-y)(1-z)
                    vec![ 0.0,  1.0,  0.0, -1.0, 0.0, 0.0, -1.0, 0.0,  0.0, 0.0], // N1 = x(1-y)(1-z)
                    vec![ 0.0,  0.0,  1.0, -1.0, 0.0, 0.0,  0.0, 0.0, -1.0, 0.0], // N2 = xy(1-z)
                    vec![ 0.0, -1.0,  1.0, -1.0, 0.0, 1.0,  1.0, 0.0, -1.0, 0.0], // N3 = (1-x)y(1-z)
                    vec![ 0.0, -1.0, -1.0,  1.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0], // N4 = (1-x)(1-y)z
                    vec![ 0.0,  1.0,  0.0,  1.0, 0.0, 0.0,  1.0, 0.0,  0.0, 0.0], // N5 = x(1-y)z
                    vec![ 0.0,  0.0,  1.0,  1.0, 0.0, 0.0,  0.0, 0.0,  1.0, 0.0], // N6 = xyz
                    vec![ 0.0, -1.0,  1.0,  1.0, 0.0, 1.0, -1.0, 0.0,  1.0, 0.0], // N7 = (1-x)yz
                ];

                let derivatives = vec![
                    vec![
                        vec![-1.0,  0.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN0/dx = -(1-y)(1-z)
                        vec![ 0.0, -1.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN0/dy = -(1-x)(1-z)
                        vec![ 0.0,  0.0, -1.0,  1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN0/dz = -(1-x)(1-y)
                    ],
                    vec![
                        vec![ 1.0,  0.0,  0.0, -1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN1/dx = (1-y)(1-z)
                        vec![ 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN1/dy = -x(1-z)
                        vec![ 0.0,  0.0,  0.0, -1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN1/dz = -x(1-y)
                    ],
                    vec![
                        vec![ 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN2/dx = y(1-z)
                        vec![ 1.0,  0.0,  0.0, -1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN2/dy = x(1-z)
                        vec![ 0.0,  0.0,  0.0, -1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN2/dz = -xy
                    ],
                    vec![
                        vec![-1.0,  0.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN3/dx = -y(1-z)
                        vec![ 1.0,  0.0,  0.0, -1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN3/dy = (1-x)(1-z)
                        vec![ 0.0,  0.0,  0.0, -1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN3/dz = -(1-x)y
                    ],
                    vec![
                        vec![-1.0,  0.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN4/dx = -(1-y)z
                        vec![ 0.0, -1.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN4/dy = -(1-x)z
                        vec![ 1.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN4/dz = (1-x)(1-y)
                    ],
                    vec![
                        vec![ 1.0,  0.0,  0.0,  1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN5/dx = (1-y)z
                        vec![ 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN5/dy = -xz
                        vec![ 0.0,  0.0,  0.0,  1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN5/dz = x(1-y)
                    ],
                    vec![
                        vec![ 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN6/dx = yz
                        vec![ 1.0,  0.0,  0.0,  1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN6/dy = xz
                        vec![ 0.0,  0.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN6/dz = xy
                    ],
                    vec![
                        vec![-1.0,  0.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN7/dx = -yz
                        vec![ 1.0,  0.0,  0.0,  1.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN7/dy = (1-x)z
                        vec![ 0.0,  0.0,  1.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0], // dN7/dz = (1-x)y
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            }


            ElementType::QuadraticHexahedron => {
                let num_nodes = 20;

                let values = vec![
                    vec![1.0, -3.0, -3.0, -3.0, 2.0, 5.0, 5.0, 2.0, 5.0, 2.0, 0.0, -2.0, -2.0, -2.0, -7.0, -2.0, 0.0, -2.0, -2.0, 0.0], // N0
                    vec![0.0, -1.0, 0.0, 0.0, 2.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -2.0, -2.0, 2.0, 3.0, 2.0, 0.0, 0.0, 0.0, 0.0], // N1
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N2
                    vec![0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 2.0, 0.0, -2.0, 3.0, 0.0, 0.0, -2.0, 2.0, 0.0], // N3
                    vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0, 2.0, 0.0, 0.0, 2.0, 0.0, 3.0, -2.0, 0.0, 2.0, -2.0, 0.0], // N4
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0], // N5
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N6
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0], // N7
                    vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N8
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N9
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N10
                    vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0], // N11
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N12
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N13
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N14
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -4.0, 0.0, 0.0], // N15
                    vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 4.0, 0.0], // N16
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0], // N17
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // N18
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -4.0, 0.0, 0.0], // N19
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 4.0, 5.0, 5.0, 0.0, -4.0, -4.0, -2.0, -7.0, -2.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 2.0, 0.0], // dN0/dx
                        vec![-3.0, 5.0, 4.0, 5.0, -2.0, -4.0, -7.0, 0.0, -4.0, -2.0, 0.0, 0.0, 2.0, 0.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0], // dN0/dy
                        vec![-3.0, 5.0, 5.0, 4.0, -2.0, -7.0, -4.0, -2.0, -4.0, 0.0, 0.0, 2.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN0/dz
                    ],
                    vec![
                        vec![-1.0, 4.0, -1.0, -1.0, 0.0, -4.0, -4.0, 2.0, 3.0, 2.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -2.0, -2.0, 0.0], // dN1/dx
                        vec![0.0, -1.0, 0.0, 0.0, -2.0, 4.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, -4.0, -2.0, 0.0, 0.0, 0.0, 0.0], // dN1/dy
                        vec![0.0, -1.0, 0.0, 0.0, -2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN1/dz
                    ],
                    vec![
                        vec![0.0, 0.0, -3.0, 0.0, 0.0, 4.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -2.0, 2.0, 0.0], // dN2/dx
                        vec![0.0, -3.0, 0.0, 0.0, 2.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -4.0, 2.0, 0.0, 0.0, 0.0, 0.0], // dN2/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN2/dz
                    ],
                    vec![
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, -2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 2.0, -2.0, 0.0], // dN3/dx
                        vec![-1.0, -1.0, 4.0, -1.0, 2.0, -4.0, 3.0, 0.0, -4.0, 2.0, 0.0, 0.0, -2.0, 0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0], // dN3/dy
                        vec![0.0, 0.0, -1.0, 0.0, 0.0, 3.0, 0.0, -2.0, 4.0, 0.0, 0.0, -2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN3/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 4.0, 0.0, 3.0, -2.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, -2.0, 2.0, 0.0], // dN4/dx
                        vec![0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 3.0, 0.0, 4.0, -2.0, 0.0, 0.0, -2.0, 0.0, -4.0, 2.0, 0.0, 0.0, 0.0, 0.0], // dN4/dy
                        vec![-1.0, -1.0, -1.0, 4.0, 2.0, 3.0, -4.0, 2.0, -4.0, 0.0, 0.0, -2.0, 0.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN4/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 4.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 2.0, -2.0, 0.0], // dN5/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0], // dN5/dy
                        vec![0.0, -3.0, 0.0, 0.0, 2.0, 1.0, 4.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN5/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 2.0, 0.0], // dN6/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0], // dN6/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN6/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -2.0, -2.0, 0.0], // dN7/dx
                        vec![0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 1.0, 0.0, 4.0, 2.0, 0.0, 0.0, 2.0, 0.0, -4.0, -2.0, 0.0, 0.0, 0.0, 0.0], // dN7/dy
                        vec![0.0, 0.0, -3.0, 0.0, 0.0, 1.0, 0.0, 2.0, 4.0, 0.0, 0.0, 2.0, 0.0, -2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN7/dz
                    ],
                    vec![
                        vec![4.0, -8.0, -4.0, -4.0, 0.0, 8.0, 8.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dx
                        vec![0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dy
                        vec![0.0, -4.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN8/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // dN9/dx
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, -8.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN9/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, -8.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dx
                        vec![0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN10/dz
                    ],
                    vec![
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0], // dN11/dx
                        vec![4.0, -4.0, -8.0, -4.0, 0.0, 8.0, 4.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dy
                        vec![0.0, 0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN11/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -8.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dy
                        vec![0.0, 4.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN12/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0], // dN13/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN13/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN14/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // dN15/dx
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/dy
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN15/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0], // dN16/dx
                        vec![0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0], // dN16/dy
                        vec![4.0, -4.0, -4.0, -8.0, 0.0, 4.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN16/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0], // dN17/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0], // dN17/dy
                        vec![0.0, 4.0, 0.0, 0.0, 0.0, -4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN17/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0], // dN18/dx
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0], // dN18/dy
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN18/dz
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0], // dN19/dx
                        vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0], // dN19/dy
                        vec![0.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0], // dN19/dz
                    ],
                ];


                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::BiquadraticQuadraticHexahedron => {
                let num_nodes = 24;

                let values = vec![
                    vec![0.0625, -0.125, 0.0625, -0.125, 0.25, -0.125, 0.0625, -0.125, 0.0625, -0.125, 0.25, -0.125, 0.25, -0.5, 0.25, -0.125, 0.25, -0.125, 0.0625, -0.125, 0.0625, 0.25, -0.5, 0.25, 0.0625, -0.125, 0.0625],
                    vec![-0.0625, 0.125, -0.0625, -0.125, 0.25, -0.125, -0.0625, 0.125, -0.0625, 0.125, -0.25, 0.125, 0.25, -0.5, 0.25, -0.125, 0.25, -0.125, -0.0625, 0.125, -0.0625, -0.25, 0.5, -0.25, 0.0625, -0.125, 0.0625],
                    vec![0.0625, -0.125, 0.0625, 0.125, -0.25, 0.125, -0.0625, 0.125, -0.0625, -0.125, 0.25, -0.125, 0.25, -0.5, 0.25, 0.125, -0.25, 0.125, -0.0625, 0.125, -0.0625, 0.25, -0.5, 0.25, 0.0625, -0.125, 0.0625],
                    vec![-0.0625, 0.125, -0.0625, 0.125, -0.25, 0.125, 0.0625, -0.125, 0.0625, 0.125, -0.25, 0.125, 0.25, -0.5, 0.25, 0.125, -0.25, 0.125, 0.0625, -0.125, 0.0625, -0.25, 0.5, -0.25, 0.0625, -0.125, 0.0625],
                    vec![-0.0625, -0.125, -0.0625, -0.125, -0.25, -0.125, -0.0625, -0.125, -0.0625, -0.125, -0.25, -0.125, -0.25, -0.5, -0.25, -0.125, -0.25, -0.125, -0.0625, -0.125, -0.0625, -0.25, -0.5, -0.25, -0.0625, -0.125, -0.0625],
                    vec![0.0625, 0.125, 0.0625, -0.125, -0.25, -0.125, -0.0625, -0.125, -0.0625, 0.125, 0.25, 0.125, -0.25, -0.5, -0.25, -0.125, -0.25, -0.125, -0.0625, -0.125, -0.0625, 0.25, 0.5, 0.25, -0.0625, -0.125, -0.0625],
                    vec![-0.0625, -0.125, -0.0625, 0.125, 0.25, 0.125, -0.0625, -0.125, -0.0625, -0.125, -0.25, -0.125, -0.25, -0.5, -0.25, 0.125, 0.25, 0.125, -0.0625, -0.125, -0.0625, -0.25, -0.5, -0.25, -0.0625, -0.125, -0.0625],
                    vec![0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, -0.25, -0.5, -0.25, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625, 0.25, 0.5, 0.25, -0.0625, -0.125, -0.0625],
                    vec![-0.25, 0.5, -0.25, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, -0.5, 1.0, -0.5, -0.25, 0.5, -0.25, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25],
                    vec![0.0, 0.0, 0.0, -0.25, 0.5, -0.25, -0.25, 0.5, -0.25, -0.5, 1.0, -0.5, 0.5, -1.0, 0.5, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.25, -0.5, 0.25],
                    vec![-0.25, 0.5, -0.25, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, -0.5, 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, -0.25, 0.5, -0.25, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25],
                    vec![0.0, 0.0, 0.0, -0.25, 0.5, -0.25, -0.25, 0.5, -0.25, 0.5, -1.0, 0.5, -0.5, 1.0, -0.5, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.25, -0.5, 0.25],
                    vec![0.25, 0.5, 0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.25, 0.5, 0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25],
                    vec![0.0, 0.0, 0.0, 0.25, 0.5, 0.25, 0.25, 0.5, 0.25, -0.5, -1.0, -0.5, -0.5, -1.0, -0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.25, -0.5, -0.25],
                    vec![0.25, 0.5, 0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.25, 0.5, 0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25],
                    vec![0.0, 0.0, 0.0, 0.25, 0.5, 0.25, 0.25, 0.5, 0.25, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.25, -0.5, -0.25],
                    vec![0.0625, 0.0, -0.0625, 0.0, 0.0, 0.0, -0.0625, 0.0, 0.0625, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0625, 0.0, -0.0625, 0.0, 0.0, 0.0, -0.0625, 0.0, 0.0625],
                    vec![-0.0625, 0.0, 0.0625, 0.0, 0.0, 0.0, 0.0625, 0.0, -0.0625, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, -0.0625, 0.0, 0.0625, 0.0, 0.0, 0.0, 0.0625, 0.0, -0.0625],
                    vec![0.0625, 0.0, -0.0625, 0.0, 0.0, 0.0, -0.0625, 0.0, 0.0625, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, -0.0625, 0.0, 0.0625, 0.0, 0.0, 0.0, 0.0625, 0.0, -0.0625],
                    vec![-0.0625, 0.0, 0.0625, 0.0, 0.0, 0.0, 0.0625, 0.0, -0.0625, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0625, 0.0, -0.0625, 0.0, 0.0, 0.0, -0.0625, 0.0, 0.0625],
                    vec![0.0, 0.0, 0.0, 0.125, 0.0, -0.125, -0.125, 0.0, 0.125, -0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.125, 0.0, -0.125],
                    vec![0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.125, 0.0, -0.125, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, -0.125, 0.0, 0.125],
                    vec![0.0, -0.5, 0.0, 0.0, 1.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0, 0.5, 0.0, 0.0, -1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                ];

                let derivatives = vec![
                    vec![
                        vec![-0.125, 0.25, -0.125, 0.25, -0.5, 0.25, -0.125, 0.25, -0.125, 0.125, -0.25, 0.125, 0.5, -1.0, 0.5, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.125, 0.25, -0.125, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, -0.25, 0.5, -0.25, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0],
                        vec![-0.125, 0.125, 0.0, 0.25, -0.25, 0.0, -0.125, 0.125, 0.0, 0.25, -0.25, 0.0, -0.5, 0.5, 0.0, 0.25, -0.25, 0.0, -0.125, 0.125, 0.0, -0.5, 0.5, 0.0, -0.125, 0.125, 0.0],
                    ],
                    vec![
                        vec![0.125, -0.25, 0.125, 0.25, -0.5, 0.25, -0.125, 0.25, -0.125, -0.125, 0.25, -0.125, -0.5, 1.0, -0.5, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.125, 0.25, -0.125, -0.125, 0.25, -0.125, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, -0.25, 0.5, -0.25, 0.0, 0.0, 0.0, -0.25, 0.5, -0.25, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0],
                        vec![0.125, -0.125, 0.0, 0.25, -0.25, 0.0, 0.125, -0.125, 0.0, -0.25, 0.25, 0.0, -0.5, 0.5, 0.0, 0.25, -0.25, 0.0, 0.125, -0.125, 0.0, 0.5, -0.5, 0.0, -0.125, 0.125, 0.0],
                    ],
                    vec![
                        vec![-0.125, 0.25, -0.125, 0.25, -0.5, 0.25, 0.125, -0.25, 0.125, -0.125, 0.25, -0.125, 0.5, -1.0, 0.5, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.125, -0.25, 0.125, -0.125, 0.25, -0.125, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.25, -0.5, 0.25, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0],
                        vec![-0.125, 0.125, 0.0, -0.25, 0.25, 0.0, 0.125, -0.125, 0.0, 0.25, -0.25, 0.0, -0.5, 0.5, 0.0, -0.25, 0.25, 0.0, 0.125, -0.125, 0.0, -0.5, 0.5, 0.0, -0.125, 0.125, 0.0],
                    ],
                    vec![
                        vec![0.125, -0.25, 0.125, 0.25, -0.5, 0.25, 0.125, -0.25, 0.125, 0.125, -0.25, 0.125, -0.5, 1.0, -0.5, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.125, -0.25, 0.125, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.25, -0.5, 0.25, 0.0, 0.0, 0.0, -0.25, 0.5, -0.25, 0.125, -0.25, 0.125, 0.0, 0.0, 0.0],
                        vec![0.125, -0.125, 0.0, -0.25, 0.25, 0.0, -0.125, 0.125, 0.0, -0.25, 0.25, 0.0, -0.5, 0.5, 0.0, -0.25, 0.25, 0.0, -0.125, 0.125, 0.0, 0.5, -0.5, 0.0, -0.125, 0.125, 0.0],
                    ],
                    vec![
                        vec![-0.125, -0.25, -0.125, -0.25, -0.5, -0.25, -0.125, -0.25, -0.125, -0.125, -0.25, -0.125, -0.5, -1.0, -0.5, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.125, -0.25, -0.125, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.25, -0.5, -0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0],
                        vec![-0.125, -0.125, 0.0, -0.25, -0.25, 0.0, -0.125, -0.125, 0.0, -0.25, -0.25, 0.0, -0.5, -0.5, 0.0, -0.25, -0.25, 0.0, -0.125, -0.125, 0.0, -0.5, -0.5, 0.0, -0.125, -0.125, 0.0],
                    ],
                    vec![
                        vec![0.125, 0.25, 0.125, -0.25, -0.5, -0.25, -0.125, -0.25, -0.125, -0.125, -0.25, -0.125, 0.5, 1.0, 0.5, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.125, -0.25, -0.125, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.25, -0.5, -0.25, 0.0, 0.0, 0.0, 0.25, 0.5, 0.25, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0],
                        vec![0.125, 0.125, 0.0, -0.25, -0.25, 0.0, -0.125, -0.125, 0.0, 0.25, 0.25, 0.0, -0.5, -0.5, 0.0, -0.25, -0.25, 0.0, -0.125, -0.125, 0.0, 0.5, 0.5, 0.0, -0.125, -0.125, 0.0],
                    ],
                    vec![
                        vec![-0.125, -0.25, -0.125, -0.25, -0.5, -0.25, 0.125, 0.25, 0.125, -0.125, -0.25, -0.125, -0.5, -1.0, -0.5, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.125, 0.25, 0.125, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, 0.25, 0.5, 0.25, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0],
                        vec![-0.125, -0.125, 0.0, 0.25, 0.25, 0.0, -0.125, -0.125, 0.0, -0.25, -0.25, 0.0, -0.5, -0.5, 0.0, 0.25, 0.25, 0.0, -0.125, -0.125, 0.0, -0.5, -0.5, 0.0, -0.125, -0.125, 0.0],
                    ],
                    vec![
                        vec![0.125, 0.25, 0.125, -0.25, -0.5, -0.25, 0.125, 0.25, 0.125, 0.125, 0.25, 0.125, 0.5, 1.0, 0.5, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.125, 0.25, 0.125, 0.125, 0.25, 0.125, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, 0.25, 0.5, 0.25, 0.0, 0.0, 0.0, 0.25, 0.5, 0.25, -0.125, -0.25, -0.125, 0.0, 0.0, 0.0],
                        vec![0.125, 0.125, 0.0, 0.25, 0.25, 0.0, 0.125, 0.125, 0.0, 0.25, 0.25, 0.0, -0.5, -0.5, 0.0, 0.25, 0.25, 0.0, 0.125, 0.125, 0.0, 0.5, 0.5, 0.0, -0.125, -0.125, 0.0],
                    ],
                    vec![
                        vec![0.5, -1.0, 0.5, 0.0, 0.0, 0.0, -0.5, 1.0, -0.5, -0.5, 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0],
                        vec![0.5, -0.5, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.0],
                    ],
                    vec![
                        vec![-0.5, 1.0, -0.5, 0.5, -1.0, 0.5, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.25, 0.5, -0.25, -0.5, 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.0, -0.5, 0.5, 0.0],
                    ],
                    vec![
                        vec![-0.5, 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, -0.5, 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0],
                        vec![0.5, -0.5, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.0],
                    ],
                    vec![
                        vec![0.5, -1.0, 0.5, -0.5, 1.0, -0.5, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.25, 0.5, -0.25, -0.5, 1.0, -0.5, 0.0, 0.0, 0.0, -0.5, 1.0, -0.5, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.25, -0.5, 0.25, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.5, -0.5, 0.0, 0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.0, -0.5, 0.5, 0.0],
                    ],
                    vec![
                        vec![-0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0],
                        vec![0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.0],
                    ],
                    vec![
                        vec![-0.5, -1.0, -0.5, -0.5, -1.0, -0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.25, 0.5, 0.25, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, -1.0, -1.0, 0.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.0, -0.5, -0.5, 0.0],
                    ],
                    vec![
                        vec![0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0],
                        vec![0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.0],
                    ],
                    vec![
                        vec![0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.25, 0.5, 0.25, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, -0.25, -0.5, -0.25, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.0, -0.5, -0.5, 0.0],
                    ],
                    vec![
                        vec![-0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0],
                        vec![0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0],
                    ],
                    vec![
                        vec![0.125, 0.0, -0.125, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0],
                        vec![0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0],
                    ],
                    vec![
                        vec![0.125, 0.0, -0.125, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0],
                        vec![0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0],
                    ],
                    vec![
                        vec![-0.125, 0.0, 0.125, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.125, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.0, 0.0, 0.0],
                        vec![0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0],
                    ],
                    vec![
                        vec![-0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, -0.25, 0.0, 0.25, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.125, 0.0, -0.125, -0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, -0.5, 0.0, 0.0, 0.0, -0.125, 0.0, 0.125, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, -0.25, 0.0],
                    ],
                    vec![
                        vec![0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.25, 0.0, -0.25, -0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.125, 0.0, 0.125, 0.25, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, -0.5, 0.0, 0.0, 0.0, 0.125, 0.0, -0.125, -0.25, 0.0, 0.25, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, -0.25, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-0.5, 0.0, 0.0, 1.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.5, 0.0, 0.0, -1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                ];

                Some(ShapeFunction { values, derivatives, num_nodes })
            },

            ElementType::TriquadraticHexahedron => {
                let num_nodes = 27;

                let values = vec![
                    vec![1.0, -3.0, 2.0, -3.0, 9.0, -6.0, 2.0, -6.0, 4.0, -3.0, 9.0, -6.0, 9.0, -27.0, 18.0, -6.0, 18.0, -12.0, 2.0, -6.0, 4.0, -6.0, 18.0, -12.0, 4.0, -12.0, 8.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 3.0, -2.0, 3.0, -9.0, 6.0, -2.0, 6.0, -4.0, 2.0, -6.0, 4.0, -6.0, 18.0, -12.0, 4.0, -12.0, 8.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -3.0, 2.0, -2.0, 6.0, -4.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 4.0, -12.0, 8.0],
                    vec![0.0, 0.0, 0.0, -1.0, 3.0, -2.0, 2.0, -6.0, 4.0, 0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -6.0, 18.0, -12.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 4.0, -12.0, 8.0],
                    vec![0.0, -1.0, 2.0, 0.0, 3.0, -6.0, 0.0, -2.0, 4.0, 0.0, 3.0, -6.0, 0.0, -9.0, 18.0, 0.0, 6.0, -12.0, 0.0, -2.0, 4.0, 0.0, 6.0, -12.0, 0.0, -4.0, 8.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 0.0, -3.0, 6.0, 0.0, 2.0, -4.0, 0.0, -2.0, 4.0, 0.0, 6.0, -12.0, 0.0, -4.0, 8.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -4.0, 8.0],
                    vec![0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 0.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 6.0, -12.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -4.0, 8.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, -4.0, 12.0, -8.0, 12.0, -36.0, 24.0, -8.0, 24.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -8.0, 24.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -8.0, 24.0, -16.0],
                    vec![0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -4.0, 12.0, -8.0, 0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 12.0, -36.0, 24.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -8.0, 24.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 4.0, -8.0, 0.0, -12.0, 24.0, 0.0, 8.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 8.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -12.0, 24.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0],
                    vec![0.0, 4.0, -4.0, 0.0, -12.0, 12.0, 0.0, 8.0, -8.0, 0.0, -12.0, 12.0, 0.0, 36.0, -36.0, 0.0, -24.0, 24.0, 0.0, 8.0, -8.0, 0.0, -24.0, 24.0, 0.0, 16.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 12.0, -12.0, 0.0, -8.0, 8.0, 0.0, 8.0, -8.0, 0.0, -24.0, 24.0, 0.0, 16.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -8.0, 8.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 16.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 8.0, -8.0, 0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -24.0, 24.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 16.0, -16.0],
                    vec![0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -16.0, 16.0, 0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 48.0, -48.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -32.0, 32.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -32.0, 32.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, -16.0, 16.0, 0.0, 48.0, -48.0, 0.0, -32.0, 32.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -32.0, 32.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, -16.0, 48.0, -32.0, 16.0, -48.0, 32.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -16.0, 32.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, -64.0, 64.0, 0.0, 64.0, -64.0],
                ];

                let derivatives = vec![
                    vec![
                        vec![-3.0, 9.0, -6.0, 9.0, -27.0, 18.0, -6.0, 18.0, -12.0, 4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-3.0, 9.0, -6.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 9.0, -27.0, 18.0, -12.0, 36.0, -24.0, 0.0, 0.0, 0.0, -6.0, 18.0, -12.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0],
                        vec![-3.0, 4.0, 0.0, 9.0, -12.0, 0.0, -6.0, 8.0, 0.0, 9.0, -12.0, 0.0, -27.0, 36.0, 0.0, 18.0, -24.0, 0.0, -6.0, 8.0, 0.0, 18.0, -24.0, 0.0, -12.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![-1.0, 3.0, -2.0, 3.0, -9.0, 6.0, -2.0, 6.0, -4.0, 4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -4.0, 12.0, -8.0, 0.0, 0.0, 0.0, -6.0, 18.0, -12.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -4.0, 0.0, -9.0, 12.0, 0.0, 6.0, -8.0, 0.0, -6.0, 8.0, 0.0, 18.0, -24.0, 0.0, -12.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 1.0, -3.0, 2.0, -2.0, 6.0, -4.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -3.0, 2.0, -4.0, 12.0, -8.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 6.0, -8.0, 0.0, 0.0, 0.0, 0.0, 6.0, -8.0, 0.0, -12.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -6.0, 18.0, -12.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![-1.0, 3.0, -2.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 3.0, -9.0, 6.0, -12.0, 36.0, -24.0, 0.0, 0.0, 0.0, -2.0, 6.0, -4.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 3.0, -4.0, 0.0, -6.0, 8.0, 0.0, 0.0, 0.0, 0.0, -9.0, 12.0, 0.0, 18.0, -24.0, 0.0, 0.0, 0.0, 0.0, 6.0, -8.0, 0.0, -12.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 3.0, -6.0, 0.0, -9.0, 18.0, 0.0, 6.0, -12.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 3.0, -6.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -9.0, 18.0, 0.0, 12.0, -24.0, 0.0, 0.0, 0.0, 0.0, 6.0, -12.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0],
                        vec![-1.0, 4.0, 0.0, 3.0, -12.0, 0.0, -2.0, 8.0, 0.0, 3.0, -12.0, 0.0, -9.0, 36.0, 0.0, 6.0, -24.0, 0.0, -2.0, 8.0, 0.0, 6.0, -24.0, 0.0, -4.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 1.0, -2.0, 0.0, -3.0, 6.0, 0.0, 2.0, -4.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 6.0, -12.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -4.0, 0.0, -3.0, 12.0, 0.0, 2.0, -8.0, 0.0, -2.0, 8.0, 0.0, 6.0, -24.0, 0.0, -4.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0, 4.0, -8.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 4.0, 0.0, 2.0, -8.0, 0.0, 0.0, 0.0, 0.0, 2.0, -8.0, 0.0, -4.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 6.0, -12.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 1.0, -2.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -3.0, 6.0, 0.0, 12.0, -24.0, 0.0, 0.0, 0.0, 0.0, 2.0, -4.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 1.0, -4.0, 0.0, -2.0, 8.0, 0.0, 0.0, 0.0, 0.0, -3.0, 12.0, 0.0, 6.0, -24.0, 0.0, 0.0, 0.0, 0.0, 2.0, -8.0, 0.0, -4.0, 16.0, 0.0],
                    ],
                    vec![
                        vec![4.0, -12.0, 8.0, -12.0, 36.0, -24.0, 8.0, -24.0, 16.0, -8.0, 24.0, -16.0, 24.0, -72.0, 48.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 16.0, -48.0, 32.0, 0.0, 0.0, 0.0, 12.0, -36.0, 24.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.0, 16.0, 0.0, 36.0, -48.0, 0.0, -24.0, 32.0, 0.0, 12.0, -16.0, 0.0, -36.0, 48.0, 0.0, 24.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 4.0, -12.0, 8.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -16.0, 0.0, -12.0, 16.0, 0.0, 0.0, 0.0, 0.0, -24.0, 32.0, 0.0, 24.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 8.0, -24.0, 16.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 12.0, -8.0, 16.0, -48.0, 32.0, 0.0, 0.0, 0.0, 4.0, -12.0, 8.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -16.0, 0.0, -24.0, 32.0, 0.0, 0.0, 0.0, 0.0, -12.0, 16.0, 0.0, 24.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 12.0, -36.0, 24.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![4.0, -12.0, 8.0, -8.0, 24.0, -16.0, 0.0, 0.0, 0.0, -12.0, 36.0, -24.0, 24.0, -72.0, 48.0, 0.0, 0.0, 0.0, 8.0, -24.0, 16.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -12.0, 16.0, 0.0, 12.0, -16.0, 0.0, 0.0, 0.0, 0.0, 36.0, -48.0, 0.0, -36.0, 48.0, 0.0, 0.0, 0.0, 0.0, -24.0, 32.0, 0.0, 24.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0, 0.0, -24.0, 48.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, -12.0, 24.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 16.0, 0.0, 12.0, -48.0, 0.0, -8.0, 32.0, 0.0, 4.0, -16.0, 0.0, -12.0, 48.0, 0.0, 8.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -4.0, 8.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -16.0, 0.0, -4.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 32.0, 0.0, 8.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -16.0, 0.0, -8.0, 32.0, 0.0, 0.0, 0.0, 0.0, -4.0, 16.0, 0.0, 8.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -12.0, 24.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, -4.0, 8.0, 0.0, 8.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -24.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -4.0, 16.0, 0.0, 4.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -48.0, 0.0, -12.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 32.0, 0.0, 8.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, -12.0, 12.0, 0.0, 36.0, -36.0, 0.0, -24.0, 24.0, 0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, -12.0, 12.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 36.0, -36.0, 0.0, -48.0, 48.0, 0.0, 0.0, 0.0, 0.0, -24.0, 24.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],
                        vec![4.0, -8.0, 0.0, -12.0, 24.0, 0.0, 8.0, -16.0, 0.0, -12.0, 24.0, 0.0, 36.0, -72.0, 0.0, -24.0, 48.0, 0.0, 8.0, -16.0, 0.0, -24.0, 48.0, 0.0, 16.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, -4.0, 4.0, 0.0, 12.0, -12.0, 0.0, -8.0, 8.0, 0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -16.0, 16.0, 0.0, 0.0, 0.0, 0.0, -24.0, 24.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 12.0, -24.0, 0.0, -8.0, 16.0, 0.0, 8.0, -16.0, 0.0, -24.0, 48.0, 0.0, 16.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -8.0, 8.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -4.0, 0.0, -16.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, -8.0, 0.0, -8.0, 16.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -24.0, 24.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, -4.0, 4.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -12.0, 0.0, -48.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 8.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, -4.0, 8.0, 0.0, 8.0, -16.0, 0.0, 0.0, 0.0, 0.0, 12.0, -24.0, 0.0, -24.0, 48.0, 0.0, 0.0, 0.0, 0.0, -8.0, 16.0, 0.0, 16.0, -32.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 48.0, -48.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 16.0, -16.0, 0.0, -32.0, 32.0, 0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 96.0, -96.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -16.0, 32.0, 0.0, 0.0, 0.0, 0.0, -48.0, 96.0, 0.0, 48.0, -96.0, 0.0, 0.0, 0.0, 0.0, 32.0, -64.0, 0.0, -32.0, 64.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 16.0, -16.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -64.0, 0.0, -32.0, 64.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 16.0, -16.0, 0.0, -48.0, 48.0, 0.0, 32.0, -32.0, 0.0, -32.0, 32.0, 0.0, 96.0, -96.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -48.0, 48.0, 0.0, 64.0, -64.0, 0.0, 0.0, 0.0, 0.0, 48.0, -48.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -48.0, 96.0, 0.0, 32.0, -64.0, 0.0, -16.0, 32.0, 0.0, 48.0, -96.0, 0.0, -32.0, 64.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 32.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -32.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 16.0, 0.0, 64.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -16.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 32.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -32.0, 64.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -16.0, 48.0, -32.0, 0.0, 0.0, 0.0, -32.0, 96.0, -64.0, 32.0, -96.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, -48.0, 32.0, -32.0, 96.0, -64.0, 0.0, 0.0, 0.0, -16.0, 48.0, -32.0, 32.0, -96.0, 64.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -48.0, 64.0, 0.0, 48.0, -64.0, 0.0, 0.0, 0.0, 0.0, 48.0, -64.0, 0.0, -48.0, 64.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 16.0, -32.0, 0.0, 0.0, 0.0, 0.0, 32.0, -64.0, 0.0, -32.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 32.0, 0.0, 32.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -32.0, 0.0, -32.0, 64.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 64.0, 0.0, 16.0, -64.0, 0.0, 0.0, 0.0, 0.0, 16.0, -64.0, 0.0, -16.0, 64.0, 0.0],
                    ],
                    vec![
                        vec![0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -64.0, 64.0, 0.0, 0.0, 0.0, 0.0, -128.0, 128.0, 0.0, 128.0, -128.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.0, -64.0, 0.0, -128.0, 128.0, 0.0, 0.0, 0.0, 0.0, -64.0, 64.0, 0.0, 128.0, -128.0, 0.0, 0.0, 0.0],
                        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64.0, -128.0, 0.0, -64.0, 128.0, 0.0, 0.0, 0.0, 0.0, -64.0, 128.0, 0.0, 64.0, -128.0, 0.0],
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
            "p" | "pressure" => Some(ValueType::Pressure),
            // Stress variants
            "s" | "stress" => Some(ValueType::Stress),
            // Strain variants
            "strain" => Some(ValueType::Strain),
            // Velocity variants
            "vel" | "velocity" => Some(ValueType::Velocity),
            // Acceleration variants
            "accel" | "acceleration" => Some(ValueType::Acceleration),
            // Temperature variants
            "t" | "temp" | "temperature" => Some(ValueType::Temperature),

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
            ValueType::Stress => "Stress",
            ValueType::Strain => "Strain",
            ValueType::Temperature => "Temperature",
            ValueType::Pressure => "Pressure",

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

#[derive(Debug)]
pub struct ElementValue {               // Defines a structure to represent element value
    pub id: usize,                      // Unique identifier for the element
    pub values: Vec<f64>,               // element values (x, y for 2D, x,y,z for 3D)
}

#[derive(Debug)]
pub struct ElementValueTypeInfo {
    pub dimension: usize,           // Physical dimension (1=scalar, 3=vector, 6=tensor)
    pub num_elements: usize,      // Total number of values for this type across all time steps
    pub elementvalue_type: ValueType,  // Type of physical quantity
    pub num_elementvalue_type: usize,      // Total number of element values type
    pub start_index: usize,         // Starting index in the element_values array
}

#[derive(Debug)]
pub struct ElementValueData {
    pub element_values: Vec<ElementValue>,                          // All element values: [type1_allnodes_alltime, type2_allnodes_alltime, ...]
    pub element_value_type_info: Vec<ElementValueTypeInfo>,   // Metadata describing element value organization
}


/// Geometric analysis module

#[derive(Debug, Clone)]
pub struct ElementQuality {
    pub element_id: usize,      // ID of the element being analyzed
    pub det_jacobian_value: f64,      // determinant of the Jacobian matrix

    // more quality metrics can be added here
}

// Structure to hold the complete mesh quality analysis results
#[derive(Debug, Clone)]
pub struct MeshQualityReport {
    pub total_elements: usize,                    // Total number of elements that were successfully analyzed
    pub element_qualities: Vec<ElementQuality>,   // Quality metrics for each individual element
    // pub statistics: QualityStatistics,           // Overall statistical summary of mesh quality
}
 
/*  
// Structure to hold statistical summary of mesh quality metrics
#[derive(Debug, Clone)]
pub struct QualityStatistics {
    pub min_jacobian: f64,              // Minimum Jacobian determinant in the mesh
    pub max_jacobian: f64,              // Maximum Jacobian determinant in the mesh
    pub avg_jacobian: f64,              // Average Jacobian determinant across all elements
    pub negative_jacobian_count: usize, // Number of elements with negative Jacobian (invalid elements)
}
*/

/// Type of finite element integration to perform
#[derive(Debug, Clone, Copy)]
pub enum IntegrationType {
    Mass,      // Mass matrix integration: ∫ ρ N_i N_j dΩ
    Stiffness, // Stiffness matrix integration
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



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monomial_polynomial_basic_operations() {
        // Test with valid 3D polynomials
        // For degree 1 in 3D: [1, x, y, z] - length should be 4
        let poly1 = vec![1.0, 2.0, 0.0, 0.0]; // Represents 1 + 2x in 3D basis
        let poly2 = vec![2.0, 1.0, 0.0, 0.0]; // Represents 2 + x in 3D basis
        
        // Test addition
        let sum = MonomialPolynomial::add(&poly1, &poly2).unwrap();
        assert_eq!(sum, vec![3.0, 3.0, 0.0, 0.0]);
        
        // Test scalar multiplication
        let scaled = MonomialPolynomial::multiply_scalar(&poly1, 2.0);
        assert_eq!(scaled, vec![2.0, 4.0, 0.0, 0.0]);
        
        // Test polynomial multiplication - this will create degree 2 polynomial
        let product = MonomialPolynomial::multiply(&poly1, &poly2).unwrap();
        // (1 + 2x)(2 + x) = 2 + x + 4x + 2x² = 2 + 5x + 2x²
        // In 3D degree 2 basis: [1, x, y, z, x², xy, xz, y², yz, z²]
        let expected_length = MonomialPolynomial::expected_length(2); // Should be 10
        assert_eq!(product.len(), expected_length);
        
        // Check specific coefficients
        let coeffs_1d = MonomialPolynomial::get_coefficients_1d(&product).unwrap();
        assert!((coeffs_1d[0] - 2.0).abs() < 1e-12); // constant term
        assert!((coeffs_1d[1] - 5.0).abs() < 1e-12); // x term
        assert!((coeffs_1d[2] - 2.0).abs() < 1e-12); // x² term
    }

    #[test]
    fn test_monomial_polynomial_evaluation() {
        // Test with a simple 3D polynomial: 1 + 2x + 3y
        let poly = vec![1.0, 2.0, 3.0, 0.0]; // [1, x, y, z] basis
        
        // Test evaluation at (1.0, 0.0, 0.0)
        let result = MonomialPolynomial::evaluate(&poly, (1.0, 0.0, 0.0)).unwrap();
        assert!((result - 3.0).abs() < 1e-12); // 1 + 2*1 + 3*0 = 3
        
        // Test evaluation at (2.0, 1.0, 0.0)
        let result = MonomialPolynomial::evaluate(&poly, (2.0, 1.0, 0.0)).unwrap();
        assert!((result - 8.0).abs() < 1e-12); // 1 + 2*2 + 3*1 = 8
    }

    #[test]
    fn test_monomial_polynomial_degree_inference() {
        // Test degree inference for various polynomial degrees
        assert_eq!(MonomialPolynomial::infer_max_degree(1).unwrap(), 0); // constant [1]
        assert_eq!(MonomialPolynomial::infer_max_degree(4).unwrap(), 1); // linear [1, x, y, z]
        assert_eq!(MonomialPolynomial::infer_max_degree(10).unwrap(), 2); // quadratic [1, x, y, z, x², xy, xz, y², yz, z²]
        assert_eq!(MonomialPolynomial::infer_max_degree(20).unwrap(), 3); // cubic [1, x, y, z, x², xy, xz, y², yz, z², x³, x²y, x²z, xy², xyz, xz², y³, y²z, yz², z³]

        // Test invalid length
        assert!(MonomialPolynomial::infer_max_degree(5).is_err());
    }

    #[test]
    fn test_monomial_polynomial_coefficient_extraction() {
        // Create a polynomial: 1 + 2x + 3y + 5x² + 6xy
        // For degree 2 in 3D: [1, x, y, z, x², xy, xz, y², yz, z²]
        let poly = vec![1.0, 2.0, 3.0, 0.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0];

        // Test 2D coefficient extraction
        let coeffs_2d = MonomialPolynomial::get_coefficients_2d(&poly).unwrap();
        assert_eq!(coeffs_2d[0][0], 1.0); // constant
        assert_eq!(coeffs_2d[1][0], 2.0); // x
        assert_eq!(coeffs_2d[0][1], 3.0); // y
        assert_eq!(coeffs_2d[1][1], 6.0); // xy
        assert_eq!(coeffs_2d[2][0], 5.0); // x²
    }
}
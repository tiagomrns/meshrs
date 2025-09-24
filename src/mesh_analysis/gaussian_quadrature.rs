use std::collections::HashMap;

use crate::structs_and_impls::*;                  // Import mesh data structures and error types from lib module
use crate::error::*;
use super::geometric_analysis::GeometricAnalysis;

use std::ops::{Add, Sub, Mul, Div, Neg};
use num_traits::{Zero, One, FromPrimitive};
use lazy_static::lazy_static;


// Symbolic type for coefficient extraction
#[derive(Debug, Clone, PartialEq)]
pub enum Symbolic {
    Constant(f64),
    Variable(String),
    Add(Box<Symbolic>, Box<Symbolic>),
    Sub(Box<Symbolic>, Box<Symbolic>),
    Mul(Box<Symbolic>, Box<Symbolic>),
    Div(Box<Symbolic>, Box<Symbolic>),
    Neg(Box<Symbolic>),
    Sqrt(Box<Symbolic>),
}

impl Symbolic {
    pub fn var(name: &str) -> Self {
        Symbolic::Variable(name.to_string())
    }
    
    pub fn constant(value: f64) -> Self {
        Symbolic::Constant(value)
    }
    
    // Extract polynomial coefficients for a given variable
    pub fn extract_coefficients(&self, variable: &str) -> HashMap<i32, f64> {
        let mut coefficients = HashMap::new();
        self._extract_coefficients(variable, 1.0, &mut coefficients);
        coefficients
    }
    
    fn _extract_coefficients(&self, variable: &str, multiplier: f64, coefficients: &mut HashMap<i32, f64>) {
        match self {
            Symbolic::Constant(c) => {
                *coefficients.entry(0).or_insert(0.0) += multiplier * c;
            }
            Symbolic::Variable(name) => {
                if name == variable {
                    *coefficients.entry(1).or_insert(0.0) += multiplier;
                } else {
                    *coefficients.entry(0).or_insert(0.0) += multiplier;
                }
            }
            Symbolic::Add(a, b) => {
                a._extract_coefficients(variable, multiplier, coefficients);
                b._extract_coefficients(variable, multiplier, coefficients);
            }
            Symbolic::Sub(a, b) => {
                a._extract_coefficients(variable, multiplier, coefficients);
                b._extract_coefficients(variable, -multiplier, coefficients);
            }
            Symbolic::Mul(a, b) => {
                // For multiplication, we need to handle polynomial expansion
                // This is simplified - for full implementation you'd need proper polynomial multiplication
                if let Symbolic::Constant(c) = **a {
                    b._extract_coefficients(variable, multiplier * c, coefficients);
                } else if let Symbolic::Constant(c) = **b {
                    a._extract_coefficients(variable, multiplier * c, coefficients);
                } else {
                    // Placeholder for proper polynomial multiplication
                    *coefficients.entry(0).or_insert(0.0) += multiplier;
                }
            }
            Symbolic::Div(a, b) => {
                // Division is complex - for now, treat as constant
                if let Symbolic::Constant(c) = **b {
                    a._extract_coefficients(variable, multiplier / c, coefficients);
                } else {
                    *coefficients.entry(0).or_insert(0.0) += multiplier;
                }
            }
            Symbolic::Neg(a) => {
                a._extract_coefficients(variable, -multiplier, coefficients);
            }
            Symbolic::Sqrt(a) => {

                if let Symbolic::Constant(c) = **a {
                    *coefficients.entry(0).or_insert(0.0) += multiplier * c.sqrt();
                }
                
            }
        }
    }
}

// Implement FloatLike for Symbolic
impl FloatLike for Symbolic {}

impl Add for Symbolic {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        // Simplify constant addition
        match (&self, &other) {
            (Symbolic::Constant(a), Symbolic::Constant(b)) => Symbolic::Constant(a + b),
            (Symbolic::Constant(a), _) if *a == 0.0 => other,
            (_, Symbolic::Constant(b)) if *b == 0.0 => self,
            _ => Symbolic::Add(Box::new(self), Box::new(other)),
        }
    }
}

impl Sub for Symbolic {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        // Simplify constant subtraction
        match (&self, &other) {
            (Symbolic::Constant(a), Symbolic::Constant(b)) => Symbolic::Constant(a - b),
            (_, Symbolic::Constant(b)) if *b == 0.0 => self,
            _ => Symbolic::Sub(Box::new(self), Box::new(other)),
        }
    }
}

impl Mul for Symbolic {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        // Simplify constant multiplication
        match (&self, &other) {
            (Symbolic::Constant(a), Symbolic::Constant(b)) => Symbolic::Constant(a * b),
            (Symbolic::Constant(a), _) if *a == 1.0 => other,
            (_, Symbolic::Constant(b)) if *b == 1.0 => self,
            (Symbolic::Constant(a), _) if *a == 0.0 => Symbolic::Constant(0.0),
            (_, Symbolic::Constant(b)) if *b == 0.0 => Symbolic::Constant(0.0),
            _ => Symbolic::Mul(Box::new(self), Box::new(other)),
        }
    }
}

impl Div for Symbolic {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        // Simplify constant division
        match (&self, &other) {
            (Symbolic::Constant(a), Symbolic::Constant(b)) => {
                if *b != 0.0 {
                    Symbolic::Constant(a / b)
                } else {
                    panic!("Division by zero")
                }
            }
            (_, Symbolic::Constant(b)) if *b == 1.0 => self,
            (Symbolic::Constant(a), _) if *a == 0.0 => Symbolic::Constant(0.0),
            _ => Symbolic::Div(Box::new(self), Box::new(other)),
        }
    }
}

impl Neg for Symbolic {
    type Output = Self;
    fn neg(self) -> Self {
        // Simplify negation
        match self {
            Symbolic::Constant(a) => Symbolic::Constant(-a),
            Symbolic::Neg(a) => *a, // Double negation
            _ => Symbolic::Neg(Box::new(self)),
        }
    }
}

impl Zero for Symbolic {
    fn zero() -> Self {
        Symbolic::Constant(0.0)
    }
    
    fn is_zero(&self) -> bool {
        match self {
            Symbolic::Constant(c) => *c == 0.0,
            _ => false,
        }
    }
}

impl One for Symbolic {
    fn one() -> Self {
        Symbolic::Constant(1.0)
    }
}

impl FromPrimitive for Symbolic {
    fn from_i64(n: i64) -> Option<Self> {
        Some(Symbolic::Constant(n as f64))
    }
    
    fn from_u64(n: u64) -> Option<Self> {
        Some(Symbolic::Constant(n as f64))
    }
    
    fn from_f64(n: f64) -> Option<Self> {
        Some(Symbolic::Constant(n))
    }
}

impl Sqrt for Symbolic {
    type Output = Symbolic;
    fn sqrt(self) -> Symbolic {
        match self {
            Symbolic::Constant(c) => Symbolic::Constant(c.sqrt()),
            other => Symbolic::Sqrt(Box::new(other)),
        }
    }
}


#[derive(Debug, Clone)]
pub enum IntegrationType {
    Mass,
    Stiffness,
}

impl IntegrationType {
    pub fn get_polynomial_order(int_type: IntegrationType, element_type: &ElementType, point: &[f64]) -> Result<usize, GaussError> {
        match int_type {
            IntegrationType::Mass => {

                let shape_function = ElementType::get_shape_functions(element_type, point).unwrap();
                let spatial_dim_element_type = ElementType::get_element_dimension(element_type).unwrap(); //1, 2, or 3

                let polynomial_order = 2 * shape_function.p_order + spatial_dim_element_type * shape_function.p_order_deriv; // or spatialdimension * det_jacobian.grad
                
                
                Ok(polynomial_order)
            },
             
            IntegrationType::Stiffness => {

                Ok(0)   // Placeholder, implement actual logic as needed

            },

            _ => return Err(GaussError::UnsupportedIntegrationType),
        }
    }

    pub fn get_coefficient_a(
        int_type: IntegrationType, 
        element_type: &ElementType,
        element_nodes: &[Node],
        material_density: f64,
    ) -> Result<Vec<Vec<f64>>, GaussError> {
        match int_type {
            IntegrationType::Mass => {
                let dim = ElementType::get_element_dimension(element_type).unwrap();

                match dim {
                    1 => {
                        let symbolic_coords: Vec<Symbolic> = vec![Symbolic::var("xi")]; 
                        let shape_function = ElementType::get_shape_functions(element_type, &symbolic_coords)
                            .ok_or_else(|| GaussError::InvalidElement("Element type not supported".to_string()))?;
                        
                        // Calculate Jacobian symbolically
                        let jacobian = GeometricAnalysis::calculate_jacobian(element_nodes, &shape_function.derivatives)?;
                        
                        let num_nodes = shape_function.values.len();
                        let mut coeff_matrix = vec![vec![0.0; num_nodes]; num_nodes];

                        for i in 0..num_nodes {
                            for j in 0..num_nodes {
                                // Calculate the integrand symbolically: N_i * N_j * detJ * density
                                let integrand = shape_function.values[i].clone() * 
                                            shape_function.values[j].clone() *  
                                            jacobian.determinant.clone();  // add material density
                                
                                // Extract the coefficients from the polynomial
                                let coefficient = extract_constant_coefficient(&integrand);
                                coeff_matrix[i][j] = coefficient;
                            }
                        }
                        
                        Ok(coeff_matrix)
                    },
                    2 => {
                        let symbolic_coords: Vec<Symbolic> = vec![Symbolic::var("xi"), Symbolic::var("eta")]; 
                        let shape_function = ElementType::get_shape_functions(element_type, &symbolic_coords)
                            .ok_or_else(|| GaussError::InvalidElement("Element type not supported".to_string()))?;
                        
                        // Calculate Jacobian symbolically
                        let jacobian = GeometricAnalysis::calculate_jacobian(element_nodes, &shape_function.derivatives)?;
                        
                        let num_nodes = shape_function.values.len();
                        let num_dof = 2; // 2 DOF per node in 2D (u, v)
                        let total_dof = num_nodes * num_dof;
                        let mut coeff_matrix = vec![vec![0.0; total_dof]; total_dof];

                        for i in 0..num_nodes {
                            for j in 0..num_nodes {
                                // Calculate the integrand symbolically: N_i * N_j * detJ * density
                                let integrand = shape_function.values[i].clone() * 
                                            shape_function.values[j].clone() *  
                                            jacobian.determinant.clone();  // add material density
                                
                                // Extract the constant coefficient (a_00) from the polynomial
                                let coefficient = extract_constant_coefficient(&integrand);
                                
                                // Fill the block matrix structure
                                // For node i, DOF α and node j, DOF β: M[i*num_dof + α][j*num_dof + β]
                                for dof_i in 0..num_dof {
                                    for dof_j in 0..num_dof {
                                        let global_i = i * num_dof + dof_i;
                                        let global_j = j * num_dof + dof_j;
                                        
                                        // Only diagonal blocks are non-zero for mass matrix (u-u and v-v coupling)
                                        if dof_i == dof_j {
                                            coeff_matrix[global_i][global_j] = coefficient;
                                        }
                                        // Off-diagonal blocks (u-v coupling) are typically zero for mass matrices
                                        // coeff_matrix[global_i][global_j] = 0.0; (already initialized to 0)
                                    }
                                }
                            }
                        }
                        
                        Ok(coeff_matrix)
                    },
                    3 => {
                        let symbolic_coords: Vec<Symbolic> = vec![Symbolic::var("xi"), Symbolic::var("eta"), Symbolic::var("psi")];
                        let shape_function = ElementType::get_shape_functions(element_type, &symbolic_coords)
                            .ok_or_else(|| GaussError::InvalidElement("Element type not supported".to_string()))?;
                        
                        // Calculate Jacobian symbolically
                        let jacobian = GeometricAnalysis::calculate_jacobian(element_nodes, &shape_function.derivatives)?;
                        
                        let num_nodes = shape_function.values.len();
                        let num_dof = 3; // 3 DOF per node in 3D (u, v, w)
                        let total_dof = num_nodes * num_dof;
                        let mut coeff_matrix = vec![vec![0.0; total_dof]; total_dof];

                        for i in 0..num_nodes {
                            for j in 0..num_nodes {
                                // Calculate the integrand symbolically: N_i * N_j * detJ * density
                                let integrand = shape_function.values[i].clone() * 
                                            shape_function.values[j].clone() *  
                                            jacobian.determinant.clone();  // add material density
                                
                                // Extract the constant coefficient (a_00) from the polynomial
                                let coefficient = extract_constant_coefficient(&integrand);
                                
                                // Fill the block matrix structure
                                // For node i, DOF α and node j, DOF β: M[i*num_dof + α][j*num_dof + β]
                                for dof_i in 0..num_dof {
                                    for dof_j in 0..num_dof {
                                        let global_i = i * num_dof + dof_i;
                                        let global_j = j * num_dof + dof_j;
                                        
                                        // Only diagonal blocks are non-zero for mass matrix (u-u, v-v, w-w coupling)
                                        if dof_i == dof_j {
                                            coeff_matrix[global_i][global_j] = coefficient;
                                        }
                                        // Off-diagonal blocks (u-v, u-w, v-w coupling) are typically zero for mass matrices
                                        // coeff_matrix[global_i][global_j] = 0.0; (already initialized to 0)
                                    }
                                }
                            }
                        }
                        
                        Ok(coeff_matrix)
                    },
                    _ => return Err(GaussError::UnsupportedDimension(dim)),
                }
            },
             
            IntegrationType::Stiffness => {
                let coeff_a = 0.0; // Placeholder, implement actual logic as needed
                Ok(vec![vec![coeff_a]])
            }, 
        }
    }
}

// Helper function to extract the constant coefficient from a symbolic expression
fn extract_constant_coefficient(expr: &Symbolic) -> f64 {
    match expr {
        Symbolic::Constant(c) => *c,
        Symbolic::Variable(_) => 0.0,
        Symbolic::Add(a, b) => extract_constant_coefficient(a) + extract_constant_coefficient(b),
        Symbolic::Sub(a, b) => extract_constant_coefficient(a) - extract_constant_coefficient(b),
        Symbolic::Mul(a, b) => {
            let coeff_a = extract_constant_coefficient(a);
            let coeff_b = extract_constant_coefficient(b);
            coeff_a * coeff_b
        }
        Symbolic::Div(a, b) => {
            let coeff_a = extract_constant_coefficient(a);
            let coeff_b = extract_constant_coefficient(b);
            if coeff_b != 0.0 { coeff_a / coeff_b } else { 0.0 }
        }
        Symbolic::Neg(a) => -extract_constant_coefficient(a),
        Symbolic::Sqrt(a) => {
            let inner = extract_constant_coefficient(a);
            if inner >= 0.0 {
                inner.sqrt()
            } else {
                // Non-real constant sqrt → treat as 0 in purely real polynomial context
                0.0
            }
        }
    }
}

// More advanced coefficient extraction for polynomial terms
fn extract_polynomial_coefficients(expr: &Symbolic, variable: &str) -> HashMap<i32, f64> {
    let mut coefficients = HashMap::new();
    extract_polynomial_coefficients_recursive(expr, variable, 1.0, 0, &mut coefficients);
    coefficients
}

fn extract_polynomial_coefficients_recursive(
    expr: &Symbolic, 
    variable: &str, 
    multiplier: f64, 
    current_power: i32, 
    coefficients: &mut HashMap<i32, f64>
) {
    match expr {
        Symbolic::Constant(c) => {
            *coefficients.entry(current_power).or_insert(0.0) += multiplier * c;
        }
        Symbolic::Variable(name) => {
            if name == variable {
                *coefficients.entry(current_power + 1).or_insert(0.0) += multiplier;
            } else {
                *coefficients.entry(current_power).or_insert(0.0) += multiplier;
            }
        }
        Symbolic::Add(a, b) => {
            extract_polynomial_coefficients_recursive(a, variable, multiplier, current_power, coefficients);
            extract_polynomial_coefficients_recursive(b, variable, multiplier, current_power, coefficients);
        }
        Symbolic::Sub(a, b) => {
            extract_polynomial_coefficients_recursive(a, variable, multiplier, current_power, coefficients);
            extract_polynomial_coefficients_recursive(b, variable, -multiplier, current_power, coefficients);
        }
        Symbolic::Mul(a, b) => {
            if let Symbolic::Constant(c) = **a {
                extract_polynomial_coefficients_recursive(b, variable, multiplier * c, current_power, coefficients);
            } else if let Symbolic::Constant(c) = **b {
                extract_polynomial_coefficients_recursive(a, variable, multiplier * c, current_power, coefficients);
            } else {
                // Still a simplified product – exact multiplication requires convolving coefficient maps.
                extract_polynomial_coefficients_recursive(a, variable, multiplier, current_power, coefficients);
                extract_polynomial_coefficients_recursive(b, variable, multiplier, current_power, coefficients);
            }
        }
        Symbolic::Div(a, b) => {
            if let Symbolic::Constant(c) = **b {
                extract_polynomial_coefficients_recursive(a, variable, multiplier / c, current_power, coefficients);
            } else {
                // Non-constant denominator → not a polynomial; best effort passthrough
                extract_polynomial_coefficients_recursive(a, variable, multiplier, current_power, coefficients);
            }
        }
        Symbolic::Neg(a) => {
            extract_polynomial_coefficients_recursive(a, variable, -multiplier, current_power, coefficients);
        }
        Symbolic::Sqrt(a) => {
            // sqrt of a polynomial is generally not a polynomial.
            // If inner is constant, propagate; else ignore as non-polynomial.
            if let Symbolic::Constant(c) = **a {
                *coefficients.entry(current_power).or_insert(0.0) += multiplier * c.sqrt();
            }
            // Otherwise do nothing (no polynomial contribution).
        }
    }
}

// Factorial calculation with memoization
lazy_static! {
    static ref FACTORIALS: Vec<f64> = {
        let mut facts = vec![1.0; 100]; // Precompute up to 100!
        for i in 1..facts.len() {
            facts[i] = facts[i-1] * (i as f64);
        }
        facts
    };
}

fn factorial(n: usize) -> f64 {
    if n < FACTORIALS.len() {
        FACTORIALS[n]
    } else {
        // Fallback for large n (using Gamma function approximation)
        (2.0 * std::f64::consts::PI * (n as f64)).sqrt() * 
        ((n as f64) / std::f64::consts::E).powi(n as i32) * 
        (1.0 + 1.0 / (12.0 * (n as f64)))
    }
}

impl GaussianQuadrature {

    /// Finds minimal Gauss points n such that error(n) <= tolerance.
    /// Stops at max_n if tolerance is not reached.
    /// Finds minimal Gauss points n such that error(n) <= tolerance.
    /// Supports dim = 1, 2, 3. Returns None if tolerance not reached.
    pub fn optimize_gauss_points(
        tolerance: f64,
        polynomial_order: usize,
        dim: usize,
        int_type: IntegrationType,
        element_type: &ElementType,
        element_nodes: &[Node],
    ) -> Result<Option<usize>, GaussError> {
        if !tolerance.is_finite() || tolerance <= 0.0 || polynomial_order == 0 {
            return Err(GaussError::InvalidTolerance);
        }

        let max_n = 2 * polynomial_order - 1;

        // Get the coefficient matrix for this element
        let material_density = 1.0; // Default value, you might want to get this from element
        let coeff_matrix = IntegrationType::get_coefficient_a(
            int_type, 
            element_type, 
            element_nodes,
            material_density,
        )?;

        // Select error function depending on dimension
        match dim {
            1 => {
                // For 1D, flatten the coefficient matrix
                let coeffs: Vec<f64> = coeff_matrix.into_iter()
                    .flat_map(|row| row.into_iter())
                    .collect();
                
                for n in 1..=max_n {
                    let err = Self::calculate_1d_error(n, &coeffs);
                    if err <= tolerance {
                        return Ok(Some(n));
                    }
                }
            }
            2 => {
                for n in 1..=max_n {
                    let err = Self::calculate_2d_error(n, &coeff_matrix);
                    if err <= tolerance {
                        return Ok(Some(n));
                    }
                }
            }
            3 => {
                // This is a simplified version
                for n in 1..=max_n {
                    let err = Self::calculate_3d_error(n, &coeff_matrix);
                    if err <= tolerance {
                        return Ok(Some(n));
                    }
                }
            }
            _ => return Err(GaussError::UnsupportedDimension(dim)),
        }

        Ok(None)
    }

    // Helper functions for error calculation
    fn calculate_1d_error(n: usize, coeffs: &[f64]) -> f64 {
        let numerator = factorial(n).powi(4);
        let denominator = (2 * n + 1) as f64 * factorial(2 * n).powi(3);
        let constant = numerator / denominator;

        let mut sum = 0.0;
        for k in (2 * n)..coeffs.len() {
            if k >= 2 * n {
                let term = coeffs[k].abs() * factorial(k) as f64 / factorial(k - 2 * n) as f64;
                sum += term;
            }
        }
        constant * sum
    }

    fn calculate_2d_error(n: usize, coeff_matrix: &[Vec<f64>]) -> f64 {
        let numerator = factorial(n).powi(4);
        let denominator = (2 * n + 1) as f64 * factorial(2 * n).powi(3);
        let constant = numerator / denominator;

        let mut sum1 = 0.0;
        let mut sum2 = 0.0;

        for i in (2 * n)..coeff_matrix.len() {
            for j in 0..coeff_matrix[i].len() {
                let term = coeff_matrix[i][j].abs() * 
                    factorial(i) as f64 / factorial(i - 2 * n) as f64 * 
                    (1.0 / (j as f64 + 1.0));
                sum1 += term;
            }
        }

        for i in 0..coeff_matrix.len() {
            for j in (2 * n)..coeff_matrix[i].len() {
                let term = coeff_matrix[i][j].abs() * 
                    factorial(j) as f64 / factorial(j - 2 * n) as f64;
                sum2 += term;
            }
        }

        constant * (sum1 + sum2)
    }

    fn calculate_3d_error(n: usize, coeff_matrix: &[Vec<f64>]) -> f64 {
        // Simplified 3D error calculation
        let numerator = factorial(n).powi(4);
        let denominator = (2 * n + 1) as f64 * factorial(2 * n).powi(3);
        numerator / denominator * 1.0 // Placeholder
    }    

    // Simplified material property getters
    fn get_material_density(_element: &Element) -> Result<f64, GaussError> {
        Ok(1.0) // Default density
    }

    fn get_elastic_modulus(_element: &Element) -> Result<f64, GaussError> {
        Ok(1.0) // Default elastic modulus
    }


}

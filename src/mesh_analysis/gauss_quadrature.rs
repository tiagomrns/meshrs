

//analyse calculate number of optimal gaussian quadrature points

use std::fs::File;                          
use std::io::{self};    // For input/output operations 

use crate::database::*;                // Import mesh data structures and error types from database module



#[derive(Debug, Clone)]
pub struct GaussianQuadrature {
    pub points: Vec<f64>, 
    pub weights: Vec<f64>,
}


impl GaussianQuadrature { // Implementation for Gauss-Legendre quadrature
    fn gauss_legendre(num_points: usize) -> Result<GaussianQuadrature, GaussError> { //gives back the needed points and weights
        match num_points {
            1 => Ok(GaussianQuadrature {
                points: vec![0.0],
                weights: vec![2.0],
            }),
            2 => Ok(GaussianQuadrature {
                points: vec![-1.0 / 3f64.sqrt(), 1.0 / 3f64.sqrt()],
                weights: vec![1.0, 1.0],
            }),
            3 => Ok(GaussianQuadrature {
                points: vec![-(3.0 / 5.0).sqrt(), 0.0, (3.0 / 5.0).sqrt()],
                weights: vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0],
            }),
            _ => Err(GaussError::UnsupportedOrder(num_points)),
        }
    }
    fn gaussian_quadrature(num_points: usize) -> Result<GaussianQuadrature, GaussError> {
        // Implementation for Gaussian quadrature (currently only supports 1D)
        Self::gauss_legendre(num_points)
    }
}
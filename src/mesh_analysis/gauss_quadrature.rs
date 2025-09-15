use std::io::{self};
use std::collections::HashMap;
use std::f64;

use crate::lib::*;
use crate::error::*;

// Import the geometric analysis module to reuse existing functionality
use super::geometric_analysis::GeometricAnalysis;

// Custom error type for Gaussian quadrature operations
#[derive(Debug, Clone)]
pub enum GaussError {
    UnsupportedOrder(usize),
    UnsupportedDimension(usize),
    IntegrationError(String),
    ElementError(ElementError),
    GeometryError(String),
    InvalidElement(String),
}

impl From<ElementError> for GaussError {
    fn from(err: ElementError) -> Self {
        GaussError::ElementError(err)
    }
}

impl From<ParseError> for GaussError {
    fn from(err: ParseError) -> Self {
        GaussError::GeometryError(format!("Parse error: {:?}", err))
    }
}

#[derive(Debug, Clone)]
pub struct GaussianQuadrature {
    pub points: Vec<f64>,
    pub weights: Vec<f64>,
}

#[derive(Debug, Clone)]
pub enum IntegrationType {
    Mass,
    Stiffness,
}

impl GaussianQuadrature {
    // Keep the existing Gauss-Legendre implementation as it's core functionality
    fn gauss_legendre(num_points: usize) -> Result<GaussianQuadrature, GaussError> {
        match num_points {
            1 => Ok(GaussianQuadrature {
                points: vec![0.0],
                weights: vec![2.0],
            }),
            2 => Ok(GaussianQuadrature {
                points: vec![-1.0_f64 / 3.0_f64.sqrt(), 1.0_f64 / 3.0_f64.sqrt()],
                weights: vec![1.0, 1.0],
            }),
            3 => Ok(GaussianQuadrature {
                points: vec![-(3.0_f64 / 5.0_f64).sqrt(), 0.0, (3.0_f64 / 5.0_f64).sqrt()],
                weights: vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0],
            }),
            4 => Ok(GaussianQuadrature {
                points: vec![
                    -((3.0_f64 + 2.0_f64 * (6.0_f64 / 5.0_f64).sqrt()) / 7.0_f64).sqrt(),
                    -((3.0_f64 - 2.0_f64 * (6.0_f64 / 5.0_f64).sqrt()) / 7.0_f64).sqrt(),
                    ((3.0_f64 - 2.0_f64 * (6.0_f64 / 5.0_f64).sqrt()) / 7.0_f64).sqrt(),
                    ((3.0_f64 + 2.0_f64 * (6.0_f64 / 5.0_f64).sqrt()) / 7.0_f64).sqrt(),
                ],
                weights: vec![
                    (18.0_f64 - 30.0_f64.sqrt()) / 36.0,
                    (18.0_f64 + 30.0_f64.sqrt()) / 36.0,
                    (18.0_f64 + 30.0_f64.sqrt()) / 36.0,
                    (18.0_f64 - 30.0_f64.sqrt()) / 36.0,
                ],
            }),
            5 => Ok(GaussianQuadrature {
                points: vec![
                    -((5.0_f64 + 2.0_f64 * (10.0_f64 / 7.0_f64).sqrt()) / 9.0_f64).sqrt(),
                    -((5.0_f64 - 2.0_f64 * (10.0_f64 / 7.0_f64).sqrt()) / 9.0_f64).sqrt(),
                    0.0,
                    ((5.0_f64 - 2.0_f64 * (10.0_f64 / 7.0_f64).sqrt()) / 9.0_f64).sqrt(),
                    ((5.0_f64 + 2.0_f64 * (10.0_f64 / 7.0_f64).sqrt()) / 9.0_f64).sqrt(),
                ],
                weights: vec![
                    (322.0_f64 - 13.0_f64 * 70.0_f64.sqrt()) / 900.0,
                    (322.0_f64 + 13.0_f64 * 70.0_f64.sqrt()) / 900.0,
                    128.0 / 225.0,
                    (322.0_f64 + 13.0_f64 * 70.0_f64.sqrt()) / 900.0,
                    (322.0_f64 - 13.0_f64 * 70.0_f64.sqrt()) / 900.0,
                ],
            }),
            _ => Err(GaussError::UnsupportedOrder(num_points)),
        }
    }

    fn gaussian_quadrature(num_points: usize) -> Result<GaussianQuadrature, GaussError> {
        Self::gauss_legendre(num_points)
    }

    // Use the existing dimension detection from geometric_analysis
    fn get_element_dimension(element_type: &ElementType) -> Result<usize, GaussError> {
        match element_type {
            ElementType::Line | ElementType::QuadraticEdge => Ok(1),
            ElementType::Triangle | ElementType::QuadraticTriangle | 
            ElementType::Quad | ElementType::QuadraticQuad | ElementType::BiquadraticQuad => Ok(2),
            ElementType::Tetra | ElementType::QuadraticTetra | ElementType::Pyramid | 
            ElementType::Wedge | ElementType::QuadraticWedge | ElementType::BiquadraticQuadraticWedge | 
            ElementType::Hexahedron | ElementType::QuadraticHexahedron | 
            ElementType::BiquadraticQuadraticHexahedron | ElementType::TriquadraticHexahedron => Ok(3),
            ElementType::Vertex => Ok(0),
            _ => Err(GaussError::UnsupportedDimension(0)),
        }
    }

    // Simplified material property getters (you can extend these based on your material system)
    fn get_material_density(_element: &Element) -> Result<f64, GaussError> {
        Ok(1.0) // Default density - replace with actual material property lookup
    }

    fn get_elastic_modulus(_element: &Element) -> Result<f64, GaussError> {
        Ok(1.0) // Default elastic modulus - replace with actual material property lookup
    }

    // Main integration function - now uses existing Jacobian from geometric_analysis
    pub fn integrate(
        &self,
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        integration_type: &IntegrationType,
    ) -> Result<f64, GaussError> {
        let dimension = Self::get_element_dimension(element_type)?;
        
        match dimension {
            1 => self.integrate_1d(element, element_type, nodes, integration_type),
            2 => self.integrate_2d(element, element_type, nodes, integration_type),
            3 => self.integrate_3d(element, element_type, nodes, integration_type),
            _ => Err(GaussError::UnsupportedDimension(dimension)),
        }
    }

    // Integration for 1D elements - now leverages GeometricAnalysis
    fn integrate_1d(
        &self,
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        integration_type: &IntegrationType,
    ) -> Result<f64, GaussError> {
        let mut integral = 0.0;
        
        // Get element nodes using existing functionality
        let element_nodes = GeometricAnalysis::get_element_nodes(element, nodes)
            .map_err(|e| GaussError::ElementError(e))?;
        
        // Get shape functions using existing functionality
        let shape_function = GeometricAnalysis::get_shape_functions(element_type)?;
        
        for (i, &xi) in self.points.iter().enumerate() {
            let weight = self.weights[i];
            
            // Calculate Jacobian using existing functionality from geometric_analysis
            let jacobian = GeometricAnalysis::calculate_jacobian(&element_nodes, &shape_function.derivatives)
                .map_err(|e| GaussError::ElementError(e))?;
            
            let f_value = match integration_type {
                IntegrationType::Mass => {
                    let mut mass_contrib = 0.0;
                    for j in 0..shape_function.values.len() {
                        for k in 0..shape_function.values.len() {
                            mass_contrib += shape_function.values[j] * shape_function.values[k];
                        }
                    }
                    let density = Self::get_material_density(element)?;
                    mass_contrib * density * jacobian.determinant.abs()
                },
                IntegrationType::Stiffness => {
                    let mut stiff_contrib = 0.0;
                    for j in 0..shape_function.derivatives.len() {
                        for k in 0..shape_function.derivatives.len() {
                            if !shape_function.derivatives[j].is_empty() && !shape_function.derivatives[k].is_empty() {
                                stiff_contrib += shape_function.derivatives[j][0] * shape_function.derivatives[k][0];
                            }
                        }
                    }
                    let elastic_modulus = Self::get_elastic_modulus(element)?;
                    stiff_contrib * elastic_modulus * jacobian.determinant.abs()
                },
            };
            
            integral += weight * f_value;
        }
        
        Ok(integral)
    }

    // Integration for 2D elements - now leverages GeometricAnalysis
    fn integrate_2d(
        &self,
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        integration_type: &IntegrationType,
    ) -> Result<f64, GaussError> {
        let mut integral = 0.0;
        
        // Get element nodes using existing functionality
        let element_nodes = GeometricAnalysis::get_element_nodes(element, nodes)
            .map_err(|e| GaussError::ElementError(e))?;
        
        // Get shape functions using existing functionality
        let shape_function = GeometricAnalysis::get_shape_functions(element_type)?;
        
        for (i, &_xi) in self.points.iter().enumerate() {
            for (j, &_eta) in self.points.iter().enumerate() {
                let weight = self.weights[i] * self.weights[j];
                
                // Calculate Jacobian using existing functionality from geometric_analysis
                let jacobian = GeometricAnalysis::calculate_jacobian(&element_nodes, &shape_function.derivatives)
                    .map_err(|e| GaussError::ElementError(e))?;
                
                let f_value = match integration_type {
                    IntegrationType::Mass => {
                        let mut mass_contrib = 0.0;
                        for k in 0..shape_function.values.len() {
                            for l in 0..shape_function.values.len() {
                                mass_contrib += shape_function.values[k] * shape_function.values[l];
                            }
                        }
                        let density = Self::get_material_density(element)?;
                        mass_contrib * density * jacobian.determinant.abs()
                    },
                    IntegrationType::Stiffness => {
                        let mut stiff_contrib = 0.0;
                        for dim in 0..2 {
                            for k in 0..shape_function.derivatives.len() {
                                for l in 0..shape_function.derivatives.len() {
                                    if dim < shape_function.derivatives[k].len() && dim < shape_function.derivatives[l].len() {
                                        stiff_contrib += shape_function.derivatives[k][dim] * shape_function.derivatives[l][dim];
                                    }
                                }
                            }
                        }
                        let elastic_modulus = Self::get_elastic_modulus(element)?;
                        stiff_contrib * elastic_modulus * jacobian.determinant.abs()
                    },
                };
                
                integral += weight * f_value;
            }
        }
        
        Ok(integral)
    }

    // Integration for 3D elements - now leverages GeometricAnalysis
    fn integrate_3d(
        &self,
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        integration_type: &IntegrationType,
    ) -> Result<f64, GaussError> {
        let mut integral = 0.0;
        
        // Get element nodes using existing functionality
        let element_nodes = GeometricAnalysis::get_element_nodes(element, nodes)
            .map_err(|e| GaussError::ElementError(e))?;
        
        // Get shape functions using existing functionality
        let shape_function = GeometricAnalysis::get_shape_functions(element_type)?;
        
        for (i, &_xi) in self.points.iter().enumerate() {
            for (j, &_eta) in self.points.iter().enumerate() {
                for (k, &_zeta) in self.points.iter().enumerate() {
                    let weight = self.weights[i] * self.weights[j] * self.weights[k];
                    
                    // Calculate Jacobian using existing functionality from geometric_analysis
                    let jacobian = GeometricAnalysis::calculate_jacobian(&element_nodes, &shape_function.derivatives)
                        .map_err(|e| GaussError::ElementError(e))?;
                    
                    let f_value = match integration_type {
                        IntegrationType::Mass => {
                            let mut mass_contrib = 0.0;
                            for l in 0..shape_function.values.len() {
                                for m in 0..shape_function.values.len() {
                                    mass_contrib += shape_function.values[l] * shape_function.values[m];
                                }
                            }
                            let density = Self::get_material_density(element)?;
                            mass_contrib * density * jacobian.determinant.abs()
                        },
                        IntegrationType::Stiffness => {
                            let mut stiff_contrib = 0.0;
                            for dim in 0..3 {
                                for l in 0..shape_function.derivatives.len() {
                                    for m in 0..shape_function.derivatives.len() {
                                        if dim < shape_function.derivatives[l].len() && dim < shape_function.derivatives[m].len() {
                                            stiff_contrib += shape_function.derivatives[l][dim] * shape_function.derivatives[m][dim];
                                        }
                                    }
                                }
                            }
                            let elastic_modulus = Self::get_elastic_modulus(element)?;
                            stiff_contrib * elastic_modulus * jacobian.determinant.abs()
                        },
                    };
                    
                    integral += weight * f_value;
                }
            }
        }
        
        Ok(integral)
    }

    // Optimization functions remain the same but now use the improved integrate method
    fn find_optimal_points_for_element(
        element: &Element,
        element_type: &ElementType,
        nodes: &[Node],
        max_points: usize,
        integration_type: &IntegrationType,
        tolerance: f64,
    ) -> Result<usize, GaussError> {
        
        let reference_quad = Self::gaussian_quadrature(5)?; 
        let reference_value = reference_quad.integrate(element, element_type, nodes, integration_type)?;
        
        for num_points in 1..=max_points {
            if let Ok(quad) = Self::gaussian_quadrature(num_points) {
                if let Ok(computed_value) = quad.integrate(element, element_type, nodes, integration_type) {
                    let error = (computed_value - reference_value).abs();
                    let relative_error = if reference_value.abs() > 1e-14 {
                        error / reference_value.abs()
                    } else {
                        error
                    };
                    
                    if relative_error < tolerance {
                        return Ok(num_points);
                    }
                }
            }
        }
        
        Ok(max_points)
    }

    pub fn optimize_mesh_integration(
        mesh_data: &MeshData,
        integration_type: &IntegrationType,
        tolerance: f64,
        max_points: usize,
    ) -> Result<HashMap<usize, usize>, GaussError> {
        
        let mut element_gauss_points = HashMap::new();
        
        for type_info in &mesh_data.element_type_info {
            if matches!(type_info.element_type, ElementType::Vertex) {
                continue;
            }
            
            let start_idx = type_info.start_index;
            let end_idx = start_idx + type_info.num_elements;
            
            for element_idx in start_idx..end_idx {
                if element_idx < mesh_data.elements.len() {
                    let element = &mesh_data.elements[element_idx];
                    
                    match Self::find_optimal_points_for_element(
                        element,
                        &type_info.element_type,
                        &mesh_data.nodes,
                        max_points,
                        integration_type,
                        tolerance,
                    ) {
                        Ok(optimal_points) => {
                            element_gauss_points.insert(element.id, optimal_points);
                        },
                        Err(e) => {
                            println!("Warning: Failed to optimize element {}: {:?}", element.id, e);
                        }
                    }
                }
            }
        }
        
        Ok(element_gauss_points)
    }

    pub fn print_optimization_results(
        optimization_results: &HashMap<usize, usize>,
        integration_type: &IntegrationType,
    ) {
        println!("Optimal Gaussian Points for {:?} Integration:", integration_type);
        println!("Element ID -> Gaussian Points");
        println!("==============================");
        
        let mut sorted_results: Vec<_> = optimization_results.iter().collect();
        sorted_results.sort_by_key(|&(id, _)| id);
        
        for (element_id, gauss_points) in sorted_results {
            println!("{:8} -> {}", element_id, gauss_points);
        }
        
        let total_elements = optimization_results.len();
        let total_points: usize = optimization_results.values().sum();
        let avg_points = total_points as f64 / total_elements as f64;
        
        println!("\nStatistics:");
        println!("Total elements: {}", total_elements);
        println!("Total Gaussian points: {}", total_points);
        println!("Average points per element: {:.2}", avg_points);
        
        let mut distribution = HashMap::new();
        for &points in optimization_results.values() {
            *distribution.entry(points).or_insert(0) += 1;
        }
        
        println!("\nDistribution:");
        for points in 1..=5 {
            if let Some(&count) = distribution.get(&points) {
                let percentage = 100.0 * count as f64 / total_elements as f64;
                println!("{} points: {} elements ({:.1}%)", points, count, percentage);
            }
        }
    }
}
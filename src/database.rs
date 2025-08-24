use std::io;

use vtkio::model::CellType;                        // Import I/O module for error handling

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
    QuadraticPyramid,            // Second order pyramid element   13 nodes   Comsol  ???
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
            ElementType::QuadraticPyramid => CellType::QuadraticPyramid,
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

#[derive(Debug)]
pub enum ParseError {                               // Defines custom error types for better error handling
    IoError(io::Error),                             // File I/O errors (e.g., file not found)
    FormatError(String),                            // File format errors (malformed data, unexpected structure)        
    NumberParseError(String),                       // Failed number conversions (invalid float/int strings)
}

// Implement automatic conversion from std::io::Error to our ParseError
// This allows us to use the ? operator with file operations
impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> Self {
        ParseError::IoError(err)  // Wrap the I/O error in our custom error type
    }
}

// Implement automatic conversion from float parsing errors
// This allows us to use ? when parsing floating point numbers
impl From<std::num::ParseFloatError> for ParseError {
    fn from(err: std::num::ParseFloatError) -> Self {
        ParseError::NumberParseError(format!("Float parse error: {}", err))
    }
}

// Implement automatic conversion from integer parsing errors
// This allows us to use ? when parsing integers
impl From<std::num::ParseIntError> for ParseError {
    fn from(err: std::num::ParseIntError) -> Self {
        ParseError::NumberParseError(format!("Int parse error: {}", err))
    }
}
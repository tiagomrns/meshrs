use std::io;    // Import I/O module for error handling

// use vtkio::model::CellType;                        


#[derive(Debug)]
pub enum ParseError {                               // Defines custom error types for better error handling
    IoError(io::Error),                             // File I/O errors (e.g., file not found)
    FormatError(String),                            // File format errors (malformed data, unexpected structure)        
    NumberParseError(String),                       // Failed number conversions (invalid float/int strings)
    ElementError(String),                           // Element-specific errors (unsupported types, bad Jacobians, etc.)
    JacobianError(String),                          // Jacobian calculation errors (singular matrices, invalid mappings)
    ShapeFunctionError(String),                     // Shape function evaluation errors
    PolynomialError(String),                        
    MeshError(String),                              
    FileFormatError(String),  
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

impl From<std::num::TryFromIntError> for ParseError {
    fn from(err: std::num::TryFromIntError) -> Self {
        ParseError::NumberParseError(format!("Integer conversion error: {}", err))
    }
}

// Define a custom error type for element analysis operations
#[derive(Debug, Clone)]
pub enum ElementError {
    InvalidElement(String),   // Error for invalid or problematic elements
    JacobianError(String),    // Error for Jacobian calculation issues
    GeometryError(String),    // Error for geometric calculation problems
    PolynomialError(String),  
    DegenerateElement(String), 
    UnsupportedElementType(String), 
}

// Implement automatic conversion from ParseError to ElementError
impl From<ParseError> for ElementError {
    fn from(err: ParseError) -> Self {
        match err {
            // Convert specific ParseError types to appropriate ElementError types
            ParseError::JacobianError(msg) => ElementError::JacobianError(msg),
            ParseError::ElementError(msg) => ElementError::InvalidElement(msg),
            // For other ParseError types, wrap them as GeometryError
            _ => ElementError::GeometryError(format!("Parse error: {:?}", err)),
        }
    }
}

// Convert &str errors to ElementError
impl From<&str> for ElementError {
    fn from(err: &str) -> Self {
        ElementError::GeometryError(err.to_string())
    }
}

// Convert String errors to ElementError
impl From<String> for ElementError {
    fn from(err: String) -> Self {
        ElementError::GeometryError(err)
    }
}

// Add this conversion for polynomial errors
impl From<PolynomialError> for ElementError {
    fn from(err: PolynomialError) -> Self {
        ElementError::PolynomialError(format!("Polynomial error: {:?}", err))
    }
}

// Polynomial error type
#[derive(Debug, Clone)]
pub enum PolynomialError {
    InvalidCoefficientLength(usize),
    InvalidDegree(u32),
    EvaluationError(String),
    ArithmeticError(String),
    BasisError(String),
}

impl std::fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PolynomialError::InvalidCoefficientLength(len) => 
                write!(f, "Invalid coefficient length: {} (not a valid 3D polynomial length)", len),
            PolynomialError::InvalidDegree(deg) => 
                write!(f, "Invalid polynomial degree: {}", deg),
            PolynomialError::EvaluationError(msg) => 
                write!(f, "Polynomial evaluation error: {}", msg),
            PolynomialError::ArithmeticError(msg) => 
                write!(f, "Polynomial arithmetic error: {}", msg),
            PolynomialError::BasisError(msg) => 
                write!(f, "Polynomial basis error: {}", msg),
        }
    }
}

// Mesh analysis specific errors
#[derive(Debug, Clone)]
pub enum AnalysisError {
    MeshQualityError(String),
    IntegrationError(String),
    ConvergenceError(String),
    InvalidTolerance(String),
}

#[derive(Debug, Clone)]
pub enum GaussError {
    UnsupportedOrder(usize),
    UnsupportedDimension(usize),
    IntegrationError(String),
    UnsupportedIntegrationType,
    ElementError(ElementError),
    GeometryError(String),
    InvalidElement(String),
    InvalidTolerance,
    InvalidMaterialProperty(String),
    PolynomialError(String),
    InvalidIntegrationType(String),
    AnalysisError(AnalysisError), 
    ConvergenceError(String),     
}

impl std::fmt::Display for GaussError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GaussError::PolynomialError(msg) => write!(f, "Polynomial error in Gauss quadrature: {}", msg),
            GaussError::ElementError(err) => write!(f, "Element error: {:?}", err),
            GaussError::GeometryError(msg) => write!(f, "Geometry error: {}", msg),
            GaussError::AnalysisError(err) => write!(f, "Analysis error: {:?}", err),
            _ => write!(f, "{:?}", self),
        }
    }
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

impl From<PolynomialError> for GaussError {
    fn from(err: PolynomialError) -> Self {
        GaussError::PolynomialError(format!("Polynomial error: {:?}", err))
    }
}

impl From<AnalysisError> for GaussError {
    fn from(err: AnalysisError) -> Self {
        GaussError::AnalysisError(err)
    }
}

impl From<&str> for GaussError {
    fn from(err: &str) -> Self {
        GaussError::GeometryError(err.to_string())
    }
}

// Implement conversion from PolynomialError to ParseError
impl From<PolynomialError> for ParseError {
    fn from(err: PolynomialError) -> Self {
        ParseError::PolynomialError(format!("{:?}", err))
    }
}

// Writer errors for output operations
#[derive(Debug)]
pub enum WriterError {
    IoError(io::Error),
    FormatError(String),
    InvalidData(String),
    VtkError(String),
}

impl From<io::Error> for WriterError {
    fn from(err: io::Error) -> Self {
        WriterError::IoError(err)
    }
}

// Parser-specific errors
#[derive(Debug)]
pub enum ParserError {
    IoError(io::Error),
    FormatError(String),
    InvalidFormat(String),
    UnsupportedElement(String),
    MissingData(String),
}

impl From<io::Error> for ParserError {
    fn from(err: io::Error) -> Self {
        ParserError::IoError(err)
    }
}

impl From<ParseError> for ParserError {
    fn from(err: ParseError) -> Self {
        match err {
            ParseError::IoError(e) => ParserError::IoError(e),
            ParseError::FormatError(s) => ParserError::FormatError(s),
            _ => ParserError::InvalidFormat(format!("{:?}", err)),
        }
    }
}
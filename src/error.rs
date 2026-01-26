// Import standard I/O module for file operation error handling
use std::io;

// Define custom error type for parsing operations
// This enum covers various parsing errors that can occur during file reading
#[derive(Debug)]
pub enum ParseError {
    IoError(io::Error),                // File I/O errors (e.g., file not found, permission denied)
    FormatError(String),               // File format errors (malformed data, unexpected structure)        
    NumberParseError(String),          // Failed number conversions (invalid float/int strings)
    ElementError(String),              // Element-specific errors (unsupported types, bad Jacobians, etc.)
    JacobianError(String),             // Jacobian calculation errors (singular matrices, invalid mappings)
    ShapeFunctionError(String),        // Shape function evaluation errors
    PolynomialError(String),           // Polynomial operation errors (invalid degree, coefficient issues)
    MeshError(String),                 // General mesh structure errors
    FileFormatError(String),           // File format specific errors
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

// Implement automatic conversion from integer overflow/underflow errors
// This allows use of ? operator when converting between integer types
impl From<std::num::TryFromIntError> for ParseError {
    fn from(err: std::num::TryFromIntError) -> Self {
        // Wrap the conversion error with context
        ParseError::NumberParseError(format!("Integer conversion error: {}", err))
    }
}

// Define a custom error type for element analysis operations
// These errors occur during finite element calculations and mesh processing
#[derive(Debug, Clone)]
pub enum ElementError {
    InvalidElement(String),         // Error for invalid or problematic elements
    JacobianError(String),         // Error for Jacobian matrix calculation issues
    GeometryError(String),         // Error for geometric calculation problems
    PolynomialError(String),       // Error for polynomial evaluation issues
    DegenerateElement(String),     // Error for degenerate elements (zero volume, etc.)
    UnsupportedElementType(String), // Error for unsupported element types
}

// Implement automatic conversion from ParseError to ElementError
// This allows seamless error propagation between parsing and element analysis
impl From<ParseError> for ElementError {
    fn from(err: ParseError) -> Self {
        match err {
            // Map specific ParseError types to appropriate ElementError variants
            ParseError::JacobianError(msg) => ElementError::JacobianError(msg),
            ParseError::ElementError(msg) => ElementError::InvalidElement(msg),
            // For other ParseError types, wrap them as GeometryError with debug info
            _ => ElementError::GeometryError(format!("Parse error: {:?}", err)),
        }
    }
}

// Convert string slice errors to ElementError for convenience
impl From<&str> for ElementError {
    fn from(err: &str) -> Self {
        // Wrap string error messages as GeometryError
        ElementError::GeometryError(err.to_string())
    }
}

// Convert owned String errors to ElementError for convenience
impl From<String> for ElementError {
    fn from(err: String) -> Self {
        // Wrap string error messages as GeometryError
        ElementError::GeometryError(err)
    }
}

// Implement conversion from PolynomialError to ElementError
// This allows polynomial errors to be handled in element analysis context
impl From<PolynomialError> for ElementError {
    fn from(err: PolynomialError) -> Self {
        // Wrap polynomial error with context information
        ElementError::PolynomialError(format!("Polynomial error: {:?}", err))
    }
}

// Define polynomial-specific error type
// Used for operations on monomial polynomials and polynomial arithmetic
#[derive(Debug, Clone)]
pub enum PolynomialError {
    InvalidCoefficientLength(usize), // Coefficient vector has invalid length for 3D polynomial
    InvalidDegree(u32),              // Polynomial degree is invalid or out of bounds
    EvaluationError(String),         // Error during polynomial evaluation
    ArithmeticError(String),         // Error during polynomial arithmetic operations
    BasisError(String),              // Error with polynomial basis functions
}

// Implement Display trait for user-friendly error messages
// This allows PolynomialError to be printed in a human-readable format
impl std::fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // Format coefficient length errors with context
            PolynomialError::InvalidCoefficientLength(len) => 
                write!(f, "Invalid coefficient length: {} (not a valid 3D polynomial length)", len),
            // Format degree validation errors
            PolynomialError::InvalidDegree(deg) => 
                write!(f, "Invalid polynomial degree: {}", deg),
            // Format evaluation errors with message
            PolynomialError::EvaluationError(msg) => 
                write!(f, "Polynomial evaluation error: {}", msg),
            // Format arithmetic operation errors
            PolynomialError::ArithmeticError(msg) => 
                write!(f, "Polynomial arithmetic error: {}", msg),
            // Format basis function errors
            PolynomialError::BasisError(msg) => 
                write!(f, "Polynomial basis error: {}", msg),
        }
    }
}

// Define mesh analysis specific error types
// Used for mesh quality analysis and integration operations
#[derive(Debug, Clone)]
pub enum AnalysisError {
    MeshQualityError(String),   // Error during mesh quality metric calculation
    IntegrationError(String),   // Error during numerical integration
    ConvergenceError(String),   // Error when numerical method fails to converge
    InvalidTolerance(String),   // Error for invalid tolerance values
}

// Define Gaussian quadrature specific error types
// Used for Gaussian integration point optimization and analysis
#[derive(Debug, Clone)]
pub enum GaussError {
    UnsupportedOrder(usize),              // Gaussian quadrature order not supported
    UnsupportedDimension(usize),          // Element dimension not supported
    IntegrationError(String),             // Error during integration
    UnsupportedIntegrationType,           // Integration type (mass/stiffness) not supported
    ElementError(ElementError),           // Wrapped element error
    GeometryError(String),                // Geometry-related error
    InvalidElement(String),               // Invalid element for integration
    InvalidTolerance,                     // Invalid tolerance value
    InvalidMaterialProperty(String),      // Invalid material property specification
    PolynomialError(String),              // Polynomial-related error
    InvalidIntegrationType(String),       // Invalid integration type specified
    AnalysisError(AnalysisError),         // Wrapped analysis error
    ConvergenceError(String),             // Convergence failure during optimization
}

// Implement Display trait for GaussError to provide readable error messages
impl std::fmt::Display for GaussError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // Format polynomial errors with context
            GaussError::PolynomialError(msg) => write!(f, "Polynomial error in Gauss quadrature: {}", msg),
            // Format element errors with debug information
            GaussError::ElementError(err) => write!(f, "Element error: {:?}", err),
            // Format geometry errors
            GaussError::GeometryError(msg) => write!(f, "Geometry error: {}", msg),
            // Format analysis errors
            GaussError::AnalysisError(err) => write!(f, "Analysis error: {:?}", err),
            // Use default debug formatting for other variants
            _ => write!(f, "{:?}", self),
        }
    }
}

// Implement conversion from ElementError to GaussError
// Allows element errors to be propagated through Gaussian quadrature operations
impl From<ElementError> for GaussError {
    fn from(err: ElementError) -> Self {
        GaussError::ElementError(err)
    }
}

// Implement conversion from ParseError to GaussError
// Allows parsing errors to be handled in Gaussian quadrature context
impl From<ParseError> for GaussError {
    fn from(err: ParseError) -> Self {
        // Wrap parse error as geometry error with debug information
        GaussError::GeometryError(format!("Parse error: {:?}", err))
    }
}

// Implement conversion from PolynomialError to GaussError
// Allows polynomial errors to be handled in Gaussian quadrature context
impl From<PolynomialError> for GaussError {
    fn from(err: PolynomialError) -> Self {
        // Wrap polynomial error with context
        GaussError::PolynomialError(format!("Polynomial error: {:?}", err))
    }
}

// Implement conversion from AnalysisError to GaussError
// Allows analysis errors to be propagated through Gaussian quadrature operations
impl From<AnalysisError> for GaussError {
    fn from(err: AnalysisError) -> Self {
        GaussError::AnalysisError(err)
    }
}

// Implement conversion from string slice to GaussError for convenience
impl From<&str> for GaussError {
    fn from(err: &str) -> Self {
        // Convert string error to GeometryError variant
        GaussError::GeometryError(err.to_string())
    }
}

// Implement conversion from PolynomialError to ParseError
// Allows polynomial errors to be handled in parsing context
impl From<PolynomialError> for ParseError {
    fn from(err: PolynomialError) -> Self {
        // Wrap polynomial error with debug formatting
        ParseError::PolynomialError(format!("{:?}", err))
    }
}

// Define writer-specific error types for output operations
#[derive(Debug)]
pub enum WriterError {
    IoError(io::Error),        // I/O error during file writing
    FormatError(String),       // Format error in output data
    InvalidData(String),       // Invalid data for output format
    VtkError(String),          // VTK-specific output errors
}

// Implement conversion from I/O errors to WriterError
impl From<io::Error> for WriterError {
    fn from(err: io::Error) -> Self {
        WriterError::IoError(err)
    }
}

// Define parser-specific error types for input operations
#[derive(Debug)]
pub enum ParserError {
    IoError(io::Error),           // I/O error during file reading
    FormatError(String),          // File format error
    InvalidFormat(String),        // Invalid file format specification
    UnsupportedElement(String),   // Unsupported element type in input file
    MissingData(String),          // Required data missing from input file
}

// Implement conversion from I/O errors to ParserError
impl From<io::Error> for ParserError {
    fn from(err: io::Error) -> Self {
        ParserError::IoError(err)
    }
}

// Implement conversion from ParseError to ParserError
// Allows parse errors to be handled in parser context
impl From<ParseError> for ParserError {
    fn from(err: ParseError) -> Self {
        match err {
            // Map specific ParseError variants to ParserError equivalents
            ParseError::IoError(e) => ParserError::IoError(e),
            ParseError::FormatError(s) => ParserError::FormatError(s),
            // Wrap other parse errors as InvalidFormat with debug info
            _ => ParserError::InvalidFormat(format!("{:?}", err)),
        }
    }
}
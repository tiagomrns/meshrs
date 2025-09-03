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

// Define a custom error type for element analysis operations
#[derive(Debug, Clone)]
pub enum ElementError {
    InvalidElement(String),   // Error for invalid or problematic elements
    JacobianError(String),    // Error for Jacobian calculation issues
    GeometryError(String),    // Error for geometric calculation problems
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

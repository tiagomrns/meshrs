// Standard library imports for file I/O operations
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

// Import custom data structures and error types
use crate::structs_and_impls::*;  // NodeValueData, NodeValue, ValueType, etc.
use crate::error::*;               // ParseError and related error types

/// Parser for COMSOL text export files containing node field values.
/// 
/// This parser handles static analysis results exported from COMSOL as text files.
/// It can parse field values such as displacement, velocity, stress, strain, etc.
/// The text files are expected to have a specific format:
/// - Header lines starting with '%' containing column names
/// - Data lines with whitespace-separated numeric values
/// - Optional unit specifications in parentheses (e.g., "(Pa)")
/// 
/// # Example File Format
/// ```text
/// % Node   displacement_x (m)   displacement_y (m)   stress_xx (Pa)
/// 0.0      1.2e-5               3.4e-6               1.5e6
/// 0.1      2.3e-5               4.5e-6               1.7e6
/// ```
pub struct ComsolTxtNodeParser;

impl ComsolTxtNodeParser {
    /// Parse a COMSOL text result file from disk containing static analysis data.
    /// 
    /// This is the main entry point for parsing COMSOL text files. It opens the file,
    /// creates a buffered reader, and delegates to the main parse function.
    /// 
    /// # Arguments
    /// * `filename` - Path to the COMSOL text file on disk
    /// 
    /// # Returns
    /// * `Ok(NodeValueData)` - Successfully parsed node values with metadata
    /// * `Err(ParseError)` - File I/O error or parsing failure
    /// 
    /// # Example
    /// ```no_run
    /// let data = ComsolTxtNodeParser::parse_file("results.txt")?;
    /// ```
    pub fn parse_file(filename: &str) -> Result<NodeValueData, ParseError> {
        // Open the file from disk (may fail if file doesn't exist or no permissions)
        let file = File::open(filename)?;
        
        // Wrap in buffered reader for efficient line-by-line reading
        let reader = BufReader::new(file);
        
        // Delegate to main parse function that works with any reader
        Self::parse(reader)
    }

    /// Main parsing function that processes COMSOL text data from any buffered reader.
    /// 
    /// This function orchestrates the three-phase parsing process:
    /// 1. Parse headers to identify column names and their physical meanings
    /// 2. Parse all numeric data rows
    /// 3. Build structured NodeValueData with proper organization and metadata
    /// 
    /// # Type Parameters
    /// * `R` - Any type that implements BufRead (file, string buffer, network stream, etc.)
    /// 
    /// # Arguments
    /// * `reader` - Buffered reader providing access to the COMSOL text data
    /// 
    /// # Returns
    /// * `Ok(NodeValueData)` - Structured node values organized by type
    /// * `Err(ParseError)` - Header parsing error, invalid data, or format mismatch
    pub fn parse<R: BufRead>(reader: R) -> Result<NodeValueData, ParseError> {
        // Make reader mutable so we can consume lines during parsing
        let mut reader = reader;
        
        // PHASE 1: Parse headers to find column layout
        // Returns: column indices mapping, original column names, and optionally the first data line
        let (column_indices, column_names, first_data_line) = Self::parse_headers(&mut reader)?;
        
        // PHASE 2: Parse all data rows into raw numeric vectors
        // The first_data_line may already contain the first row if read during header parsing
        let data = Self::parse_data(&mut reader, &column_indices, first_data_line)?;
        
        // PHASE 3: Build the final structured NodeValueData
        // Converts raw numeric data into properly typed and organized NodeValue structs
        Self::build_node_value_data(data, &column_indices, &column_names)
    }

    /// Parse the header section of a COMSOL text file to extract column information.
    /// 
    /// This function reads lines until it finds the first data line, collecting header
    /// information along the way. COMSOL text files have headers starting with '%'.
    /// The last header line before data contains the column names (and optional units).
    /// 
    /// # Arguments
    /// * `reader` - Mutable reference to buffered reader (will consume header lines)
    /// 
    /// # Returns
    /// A tuple containing:
    /// * `ColumnIndices` - Mapping of columns to value types and dimensions
    /// * `Vec<String>` - List of original column names in order
    /// * `Option<String>` - First data line if encountered (to avoid re-reading)
    /// 
    /// # Errors
    /// * Returns `ParseError::FormatError` if no headers found or invalid format
    fn parse_headers<R: BufRead>(reader: &mut R) -> Result<(ColumnIndices, Vec<String>, Option<String>), ParseError> {
        // Buffer for reading each line from the file
        let mut line = String::new();
        
        // Store the last header line we encounter (contains column names)
        let mut last_header_line = String::new();
        
        // Read file line by line until we find the first data line
        while reader.read_line(&mut line)? > 0 {
            // Remove leading/trailing whitespace for easier processing
            let trimmed = line.trim();
            
            // Check if this is a header line (COMSOL headers start with '%')
            if trimmed.starts_with('%') {
                // Remove the '%' prefix and any extra whitespace, then save as last header
                last_header_line = trimmed.trim_start_matches('%').trim().to_string();
                
                // Clear the line buffer for reuse in next iteration
                line.clear();
            } 
            // Check if we found a non-empty, non-header line (this is data)
            else if !trimmed.is_empty() {
                // Save the first data line so we don't lose it
                let first_data_line = Some(trimmed.to_string());
                
                // Parse column names from the last header line we saved
                let column_names = Self::parse_header_line(&last_header_line)?;
                
                // Map column names to their types (displacement, stress, etc.) and indices
                let column_indices = Self::map_columns(&column_names)?;
                
                // Return all header info plus the first data line
                return Ok((column_indices, column_names, first_data_line));
            }
            
            // Clear line buffer for next iteration
            line.clear();
        }
        
        // If we reach here, file ended without any data lines
        // Check if we at least found a header
        if last_header_line.is_empty() {
            return Err(ParseError::FormatError("No header found".to_string()));
        }
        
        // Parse the header we found even though there's no data
        let column_names = Self::parse_header_line(&last_header_line)?;
        let column_indices = Self::map_columns(&column_names)?;
        
        // Return header info without a first data line
        Ok((column_indices, column_names, None))
    }

    /// Parse a single header line to extract column names, filtering out unit specifications.
    /// 
    /// COMSOL text files may include unit specifications in parentheses, such as:
    /// - Standalone: "displacement_x (m) displacement_y (m)"
    /// - Inline: "stress_xx(Pa) stress_yy(Pa)"
    /// 
    /// This function extracts just the column names, converting them to lowercase
    /// for case-insensitive matching.
    /// 
    /// # Arguments
    /// * `header_line` - The header line with '%' already removed
    /// 
    /// # Returns
    /// * `Ok(Vec<String>)` - List of column names in lowercase
    /// 
    /// # Example
    /// Input: "Node displacement_x (m) displacement_y (m) stress_xx(Pa)"
    /// Output: ["node", "displacement_x", "displacement_y", "stress_xx"]
    fn parse_header_line(header_line: &str) -> Result<Vec<String>, ParseError> {
        // Collection to store the extracted column names
        let mut cols = Vec::new();
        
        // Create an iterator over whitespace-separated tokens with peek capability
        let mut it = header_line.split_whitespace().peekable();

        // Process each token in the header line
        while let Some(tok) = it.next() {
            // Skip standalone unit tokens like "(Pa)" or "(m)"
            // These are separate tokens that we don't want as column names
            if tok.starts_with('(') && tok.ends_with(')') {
                continue;
            }

            // Handle inline units like "stress_xx(Pa)" by splitting on '('
            // Takes only the part before '(' as the column name
            let name = tok.split('(').next().unwrap_or("").trim().to_lowercase();
            
            // Skip empty names (can happen with malformed input)
            if name.is_empty() { 
                continue; 
            }

            // Check if the next token is a standalone unit specification
            // If so, consume it (skip it) so we don't process it as a column name
            if let Some(peek) = it.peek() {
                if peek.starts_with('(') && peek.ends_with(')') {
                    it.next(); // Consume the unit token
                }
            }

            // Add the extracted column name to our collection
            cols.push(name);
        }
        
        // Return the extracted column names
        Ok(cols)
    }


    /// Map column names to their indices and determine physical value types and dimensions.
    /// 
    /// This function performs a two-pass analysis:
    /// 1. Group columns by their base name (e.g., "displacement_x" → base: "displacement", component: "x")
    /// 2. Determine the dimension based on number of components:
    ///    - 1 component: Scalar value
    ///    - 3 components: Vector (x, y, z)
    ///    - 6 components: Tensor (xx, yy, zz, xy, xz, yz)
    /// 
    /// # Arguments
    /// * `column_names` - Slice of column names from the header
    /// 
    /// # Returns
    /// * `Ok(ColumnIndices)` - Mapping of columns to value types with dimensions
    /// * `Err(ParseError)` - No recognized value columns found
    /// 
    /// # Example
    /// Input: ["node", "displacement_x", "displacement_y", "displacement_z"]
    /// Output: ColumnIndices with one entry: ([1,2,3], ValueType::Displacement, 3)
    fn map_columns(column_names: &[String]) -> Result<ColumnIndices, ParseError> {
        // Initialize the structure that will hold our column mapping
        let mut indices = ColumnIndices {
            // Will store tuples of (column_indices, value_type, dimension)
            value_columns: Vec::new(),
        };

        // PASS 1: Identify value columns and their components
        // HashMap groups components by their base name
        // Key: base name (e.g., "displacement"), Value: list of (column_index, component) pairs
        let mut value_groups: HashMap<String, Vec<(usize, String)>> = HashMap::new();
        
        // Iterate through each column name with its position index
        for (idx, name) in column_names.iter().enumerate() {
            // Try to parse as a component name (e.g., "displacement_x" → ("displacement", "x"))
            if let Some((base_name, component)) = Self::parse_component_name(name) {
                // Add this component to the group for its base name
                value_groups.entry(base_name).or_insert_with(Vec::new).push((idx, component));
            }
            // Also check for simple scalar names without components (e.g., just "temperature")
            else if ValueType::from_str_txt(name).is_some() {
                // Treat as a scalar value with a pseudo-component named "scalar"
                value_groups.entry(name.clone()).or_insert_with(Vec::new).push((idx, "scalar".to_string()));
            }
            // Note: Column names that don't match any value type are silently ignored
            // (e.g., "node", "element", coordinate columns, etc.)
        }

        // PASS 2: Organize components and determine dimensions for each value type
        for (base_name, mut components) in value_groups {
            // Verify that the base name corresponds to a recognized value type
            if let Some(vtype) = ValueType::from_str_txt(&base_name) {
                // Sort components in standard order (x, y, z for vectors; xx, yy, zz, xy, xz, yz for tensors)
                // This ensures consistent data layout regardless of column order in file
                components.sort_by(|a, b| Self::component_order(&a.1).cmp(&Self::component_order(&b.1)));
                
                // Dimension is determined by the number of components:
                // 1 = scalar, 3 = vector, 6 = tensor (symmetric)
                let dim = components.len();
                
                // Extract just the column indices, discarding component names
                let column_indices: Vec<usize> = components.into_iter().map(|(idx, _)| idx).collect();
                
                // Store the complete mapping for this value type
                indices.value_columns.push((column_indices, vtype, dim));
            }
        }

        // Validation: Ensure we found at least one valid value column
        if indices.value_columns.is_empty() {
            return Err(ParseError::FormatError(
                "No recognized value columns found".to_string()
            ));
        }

        // Return the complete column mapping
        Ok(indices)
    }

    /// Parse component names with underscores into base name and component suffix.
    /// 
    /// COMSOL exports field components with naming conventions like:
    /// - Vector components: "displacement_x", "velocity_y", "displacement_z"
    /// - Tensor components: "stress_xx", "stress_yy", "strain_xy"
    /// 
    /// This function splits on the last underscore and validates that the base
    /// name is a recognized value type.
    /// 
    /// # Arguments
    /// * `name` - Column name to parse
    /// 
    /// # Returns
    /// * `Some((base, component))` - If successfully parsed and base is valid
    /// * `None` - If no underscore found or base name not recognized
    /// 
    /// # Examples
    /// - "displacement_x" → Some(("displacement", "x"))
    /// - "stress_xy" → Some(("stress", "xy"))
    /// - "temperature" → None (no underscore)
    /// - "unknown_x" → None (base not a recognized value type)
    fn parse_component_name(name: &str) -> Option<(String, String)> {
        // Find the last underscore in the name (rightmost occurrence)
        // This handles cases like "some_field_name_x" correctly
        if let Some(pos) = name.rfind('_') {
            // Split the name: everything before underscore is the base name
            let base = &name[..pos];
            
            // Everything after underscore is the component identifier (x, y, z, xx, xy, etc.)
            let component = &name[pos + 1..];
            
            // Verify that the base name is a recognized physical value type
            // (e.g., "displacement", "stress", "velocity", etc.)
            if ValueType::from_str_txt(base).is_some() {
                // Return the parsed components as owned strings
                return Some((base.to_string(), component.to_string()));
            }
        }
        
        // Return None if no underscore found or base name not recognized
        None
    }

    /// Define the standard ordering for field components.
    /// 
    /// This ensures components are always stored in a consistent order regardless
    /// of how they appear in the COMSOL export file. This is important for:
    /// - Consistent data layout across different files
    /// - Proper interpretation of vector/tensor data
    /// - Matching expected component order in downstream processing
    /// 
    /// # Component Order
    /// - Vector: x=0, y=1, z=2
    /// - Tensor (symmetric, Voigt notation): xx=0, yy=1, zz=2, xy=3, xz=4, yz=5
    /// - Scalar: scalar=0
    /// - Unknown: 999 (sorted last)
    /// 
    /// # Arguments
    /// * `component` - Component identifier string (e.g., "x", "yy", "xy")
    /// 
    /// # Returns
    /// * `usize` - Numeric order value for sorting
    fn component_order(component: &str) -> usize {
        match component {
            // Vector components (standard 3D Cartesian)
            "x" => 0,   // X-component
            "y" => 1,   // Y-component
            "z" => 2,   // Z-component
            
            // Symmetric tensor components in Voigt notation
            // Normal stresses/strains
            "xx" => 0,  // X-normal component
            "yy" => 1,  // Y-normal component
            "zz" => 2,  // Z-normal component
            // Shear stresses/strains
            "xy" => 3,  // XY-shear component
            "xz" => 4,  // XZ-shear component
            "yz" => 5,  // YZ-shear component
            
            // Scalar quantities (single value)
            "scalar" => 0,
            
            // Unrecognized components get sorted to the end
            _ => 999,
        }
    }

    /// Parse all numeric data rows from the file.
    /// 
    /// This function reads all remaining lines after the headers and converts them
    /// to vectors of f64 values. Each line becomes one row in the output matrix.
    /// The function handles:
    /// - Scientific notation (e.g., 1.23e-5)
    /// - Standard decimal notation
    /// - Whitespace-separated values
    /// - Empty lines and comment lines (skipped)
    /// 
    /// # Arguments
    /// * `reader` - Mutable reference to buffered reader (positioned after headers)
    /// * `_indices` - Column indices structure (currently unused but kept for future extensions)
    /// * `first_data_line` - Optional first data line already read during header parsing
    /// 
    /// # Returns
    /// * `Ok(Vec<Vec<f64>>)` - Matrix of numeric data (rows × columns)
    /// * `Err(ParseError)` - Number parsing error or no data found
    /// 
    /// # Note
    /// Each inner Vec<f64> represents one node's data across all columns
    fn parse_data<R: BufRead>(
        reader: &mut R,
        _indices: &ColumnIndices,
        first_data_line: Option<String>,
    ) -> Result<Vec<Vec<f64>>, ParseError> {
        // Storage for all data rows (each row is a Vec<f64>)
        let mut data = Vec::new();
        
        // Buffer for reading lines from the file
        let mut line = String::new();

        // Process the first data line if we already read it during header parsing
        // This avoids losing the first row of data
        if let Some(first_line) = first_data_line {
            // Parse the first line: split by whitespace and convert each token to f64
            let values: Result<Vec<f64>, _> = first_line
                .split_whitespace()                     // Split on any whitespace
                .map(|s| s.parse().map_err(|e|         // Try to parse each token as f64
                    ParseError::NumberParseError(       // Convert parse error to our error type
                        format!("Invalid number '{}': {}", s, e)
                    )
                ))
                .collect();                             // Collect all Results into one Result<Vec<f64>>
            
            // Add the parsed first row to our data collection
            // The ? operator will early-return if parsing failed
            data.push(values?);
        }

        // Continue reading remaining data lines from the file
        while reader.read_line(&mut line)? > 0 {
            // Remove leading/trailing whitespace
            let trimmed = line.trim();
            
            // Skip empty lines and comment lines (starting with '%')
            if trimmed.is_empty() || trimmed.starts_with('%') {
                line.clear(); // Clear buffer for next iteration
                continue;
            }

            // Parse this line: split by whitespace and convert to f64
            let values: Result<Vec<f64>, _> = trimmed
                .split_whitespace()                     // Split on any whitespace
                .map(|s| s.parse().map_err(|e|         // Try to parse each token as f64
                    ParseError::NumberParseError(       // Convert parse error to our error type
                        format!("Invalid number '{}': {}", s, e)
                    )
                ))
                .collect();                             // Collect all Results into one Result<Vec<f64>>
            
            // Add the parsed row to our data collection
            data.push(values?);
            
            // Clear line buffer for next iteration
            line.clear();
        }

        // Validation: Ensure we found at least one data row
        if data.is_empty() {
            return Err(ParseError::FormatError("No data rows found".to_string()));
        }

        // Return the complete data matrix
        Ok(data)
    }

    /// Build the final NodeValueData structure from raw numeric data.
    /// 
    /// This function transforms the raw data matrix into a structured format with:
    /// - Individual NodeValue structs for each node and value type
    /// - Metadata describing the organization and dimensions
    /// - Proper indexing for efficient access
    /// 
    /// The output structure organizes data by value type (all displacement values
    /// together, all stress values together, etc.) rather than by node. This
    /// improves cache locality when processing all values of one type.
    /// 
    /// # Arguments
    /// * `data` - Raw numeric data as a matrix (rows = nodes, columns = all fields)
    /// * `indices` - Column mapping that identifies which columns belong to which value types
    /// * `_column_names` - Original column names (currently unused)
    /// 
    /// # Returns
    /// * `Ok(NodeValueData)` - Structured data with NodeValue objects and metadata
    /// * `Err(ParseError)` - If no data found
    /// 
    /// # Data Layout
    /// The resulting NodeValueData has all values organized as:
    /// [displacement_node0, displacement_node1, ..., stress_node0, stress_node1, ...]
    /// With metadata indicating where each value type's section begins.
    fn build_node_value_data(
        data: Vec<Vec<f64>>,
        indices: &ColumnIndices,
        _column_names: &[String],
    ) -> Result<NodeValueData, ParseError> {
        // Number of nodes equals the number of data rows
        let num_nodes = data.len();
        
        // Validation: Ensure we have data to process
        if num_nodes == 0 {
            return Err(ParseError::FormatError("No data found".to_string()));
        }

        // Master collection of all NodeValue structs (organized by value type)
        let mut all_node_values = Vec::new();
        
        // Metadata describing each value type's location and properties
        let mut node_value_type_info = Vec::new();
        
        // Track where each value type's section starts in all_node_values
        let mut start_index = 0;

        // Process each value type (displacement, stress, velocity, etc.)
        for (column_indices, vtype, dim) in &indices.value_columns {
            // Temporary collection for this value type's NodeValue structs
            let mut type_node_values = Vec::new();
            
            // Create one NodeValue struct for each node
            for node_idx in 0..num_nodes {
                // Get the data row for this node
                let row = &data[node_idx];
                
                // Collect the component values for this value type
                let mut values = Vec::new();
                
                // Extract each component value from the appropriate columns
                for &col_idx in column_indices {
                    // Check if this column exists in the current row
                    if col_idx < row.len() {
                        // Use the actual value from the data
                        values.push(row[col_idx]);
                    } else {
                        // Use 0.0 as default if column is missing (shouldn't happen normally)
                        values.push(0.0);
                    }
                }
                
                // Create a NodeValue struct for this node and value type
                type_node_values.push(NodeValue {
                    id: node_idx,  // Node identifier (0-indexed row number)
                    values,        // All component values for this value type
                });
            }
            
            // Only add this value type if we successfully created NodeValues for it
            if !type_node_values.is_empty() {
                // Create metadata describing this value type's properties and location
                node_value_type_info.push(NodeValueTypeInfo {
                    dimension: *dim,                // Number of components (1/3/6)
                    num_nodes: num_nodes,           // Total number of nodes
                    nodevalue_type: vtype.clone(),  // Type (Displacement, Stress, etc.)
                    num_nodevalue_type: 1,          // Always 1 for static analysis
                    start_index,                     // Where this type's data starts
                });
                
                // Add all NodeValues for this type to the master collection
                all_node_values.extend(type_node_values);
                
                // Update start_index for the next value type
                start_index = all_node_values.len();
            }
        }

        // Create and return the final structured data
        Ok(NodeValueData {
            node_values: all_node_values,           // All NodeValue structs
            node_value_type_info,                   // Metadata for each value type
        })
    }
}

/// Internal structure for mapping column indices to physical value types.
/// 
/// This structure is created during header parsing and used during data parsing
/// to correctly interpret which columns contain which physical quantities.
/// 
/// # Fields
/// * `value_columns` - Vector of tuples, each containing:
///   - `Vec<usize>` - Column indices for all components of this value type
///   - `ValueType` - The physical quantity (Displacement, Stress, etc.)
///   - `usize` - Dimension (number of components: 1=scalar, 3=vector, 6=tensor)
/// 
/// # Example
/// For a file with columns: [Node, displacement_x, displacement_y, displacement_z, stress_xx]
/// This would contain:
/// - ([1, 2, 3], ValueType::Displacement, 3)
/// - ([4], ValueType::Stress, 1)
struct ColumnIndices {
    value_columns: Vec<(Vec<usize>, ValueType, usize)>,
}


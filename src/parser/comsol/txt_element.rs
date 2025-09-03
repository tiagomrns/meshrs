use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

//use crate::database::*;                     // Import mesh data structures and error types from database module

use crate::lib::*;                     // Import mesh data structures and error types from lib module
use crate::error::*;                     // Import mesh data structures and error types from error module

pub struct ComsolTxtElementParser;

impl ComsolTxtElementParser {
    /// Parse a COMSOL txt result file from disk for element data (static data only)
    pub fn parse_file(filename: &str) -> Result<ElementValueData, ParseError> {
        let file = File::open(filename)?; // Open the file from disk
        let reader = BufReader::new(file); // Wrap in buffered reader for efficient line reading
        Self::parse(reader) // Call main parse function with the reader
    }

    /// Main parsing function for static element data
    pub fn parse<R: BufRead>(reader: R) -> Result<ElementValueData, ParseError> {
        let mut reader = reader; // Make reader mutable for line-by-line reading
        
        // Step 1: Parse headers to find column layout
        let (column_indices, column_names, first_data_line) = Self::parse_headers(&mut reader)?; // Extract column structure from file headers
        
        // Step 2: Parse all data rows
        let data = Self::parse_data(&mut reader, &column_indices, first_data_line)?; // Read all numeric data rows
        
        // Step 3: Build the final ElementValueData structure
        Self::build_element_value_data(data, &column_indices, &column_names) // Convert raw data into structured ElementValue objects
    }

    /// Parse header section to find column names and their indices
    fn parse_headers<R: BufRead>(reader: &mut R) -> Result<(ColumnIndices, Vec<String>, Option<String>), ParseError> {
        let mut line = String::new(); // Buffer to read each line
        let mut last_header_line = String::new(); // Store the last header line containing column names
        
        // Read file line by line until we find data
        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim(); // Remove whitespace from line
            
            if trimmed.starts_with('%') { // COMSOL header lines start with %
                last_header_line = trimmed.trim_start_matches('%').trim().to_string(); // Remove % and store header content
                line.clear(); // Clear buffer for next line
            } else if !trimmed.is_empty() { // Found first non-empty, non-header line (data line)
                let first_data_line = Some(trimmed.to_string()); // Save first data line to process later
                let column_names = Self::parse_header_line(&last_header_line)?; // Extract column names from last header
                let column_indices = Self::map_columns(&column_names)?; // Map column names to their types and indices
                return Ok((column_indices, column_names, first_data_line)); // Return parsed header info
            }
            line.clear(); // Clear buffer for next iteration
        }
        
        // Handle case where file ends without data lines
        if last_header_line.is_empty() {
            return Err(ParseError::FormatError("No header found".to_string())); // No header found at all
        }
        
        let column_names = Self::parse_header_line(&last_header_line)?; // Parse the header we found
        let column_indices = Self::map_columns(&column_names)?; // Map columns even without data
        
        Ok((column_indices, column_names, None)) // Return header info without first data line
    }

    /// Parse a single header line to extract column names
    fn parse_header_line(header_line: &str) -> Result<Vec<String>, ParseError> {
        let mut cols = Vec::new();
        let mut it = header_line.split_whitespace().peekable();

        while let Some(tok) = it.next() {
            // Skip standalone unit tokens like "(Pa)"
            if tok.starts_with('(') && tok.ends_with(')') {
                continue;
            }

            // Handle inline units like "stress_xx(Pa)" just in case
            let name = tok.split('(').next().unwrap_or("").trim().to_lowercase();
            if name.is_empty() { continue; }

            // If the next token is a standalone units token, skip it
            if let Some(peek) = it.peek() {
                if peek.starts_with('(') && peek.ends_with(')') {
                    it.next();
                }
            }

            cols.push(name);
        }
        Ok(cols)
    }


    /// Map column names to their indices and determine column types
    fn map_columns(column_names: &[String]) -> Result<ColumnIndices, ParseError> {
        let mut indices = ColumnIndices { // Initialize structure to hold column mapping
            value_columns: Vec::new(), // Will store (column_indices, value_type, dimension) tuples
        };

        // First pass: identify value columns and their components
        let mut value_groups: HashMap<String, Vec<(usize, String)>> = HashMap::new(); // Group components by base name
        
        for (idx, name) in column_names.iter().enumerate() { // Iterate through each column name with its index
            // Check for vector components like "displacement_x", "stress_xx", etc.
            if let Some((base_name, component)) = Self::parse_component_name(name) { // Try to split name into base + component
                value_groups.entry(base_name).or_insert_with(Vec::new).push((idx, component)); // Group component with its base name
            }
            // Also check for simple names without components
            else if ValueType::from_str_txt(name).is_some() { // Check if name is a recognized scalar value type
                value_groups.entry(name.clone()).or_insert_with(Vec::new).push((idx, "scalar".to_string())); // Treat as scalar component
            }
        }

        // Second pass: organize components and determine dimensions
        for (base_name, mut components) in value_groups { // Process each group of components
            if let Some(vtype) = ValueType::from_str_txt(&base_name) { // Verify base name is a valid value type
                // Sort components for consistent ordering (x, y, z for vectors; xx, yy, zz, xy, xz, yz for tensors)
                components.sort_by(|a, b| Self::component_order(&a.1).cmp(&Self::component_order(&b.1))); // Sort by predefined component order
                
                let dim = components.len(); // Number of components determines dimension (1=scalar, 3=vector, 6=tensor)
                let column_indices: Vec<usize> = components.into_iter().map(|(idx, _)| idx).collect(); // Extract just the column indices
                
                indices.value_columns.push((column_indices, vtype, dim)); // Store mapping for this value type
            }
        }

        // Require at least one value column
        if indices.value_columns.is_empty() {
            return Err(ParseError::FormatError( // Error if no valid columns found
                "No recognized value columns found".to_string()
            ));
        }

        Ok(indices) // Return the column mapping
    }

    /// Parse component names like "displacement_x" -> ("displacement", "x")
    fn parse_component_name(name: &str) -> Option<(String, String)> {
        if let Some(pos) = name.rfind('_') { // Find last underscore in name
            let base = &name[..pos]; // Everything before underscore is base name
            let component = &name[pos + 1..]; // Everything after underscore is component
            
            // Check if the base name is a recognized value type
            if ValueType::from_str_txt(base).is_some() { // Verify base is valid value type
                return Some((base.to_string(), component.to_string())); // Return parsed components
            }
        }
        None // Return None if parsing failed
    }

    /// Define component ordering for consistent arrangement
    fn component_order(component: &str) -> usize {
        match component { // Define standard ordering for components
            "x" => 0, "y" => 1, "z" => 2, // Vector components
            "xx" => 0, "yy" => 1, "zz" => 2, "xy" => 3, "xz" => 4, "yz" => 5, // Tensor components
            "scalar" => 0, // Scalar component
            _ => 999, // Unknown components go last
        }
    }

    /// Parse all data rows (static data only)
    fn parse_data<R: BufRead>(
        reader: &mut R, // Buffered reader for file access
        _indices: &ColumnIndices, // Column indices (not used in current implementation)
        first_data_line: Option<String>, // First data line if already read during header parsing
    ) -> Result<Vec<Vec<f64>>, ParseError> {
        let mut data = Vec::new(); // Store all data rows as Vec<Vec<f64>>
        let mut line = String::new(); // Buffer for reading lines

        // Process first data line if we have one
        if let Some(first_line) = first_data_line { // If we already read first data line during header parsing
            let values: Result<Vec<f64>, _> = first_line // Parse first line into numbers
                .split_whitespace() // Split by whitespace
                .map(|s| s.parse().map_err(|e| // Convert each token to f64
                    ParseError::NumberParseError(format!("Invalid number '{}': {}", s, e)) // Custom error for parse failures
                ))
                .collect(); // Collect all parsed numbers
            data.push(values?); // Add first row to data collection
        }

        // Continue reading data lines
        while reader.read_line(&mut line)? > 0 { // Read each remaining line
            let trimmed = line.trim(); // Remove whitespace
            
            if trimmed.is_empty() || trimmed.starts_with('%') { // Skip empty lines and comments
                line.clear(); // Clear buffer
                continue; // Skip to next line
            }

            let values: Result<Vec<f64>, _> = trimmed // Parse current line into numbers
                .split_whitespace() // Split by whitespace
                .map(|s| s.parse().map_err(|e| // Convert each token to f64
                    ParseError::NumberParseError(format!("Invalid number '{}': {}", s, e)) // Custom error for parse failures
                ))
                .collect(); // Collect all parsed numbers
            
            data.push(values?); // Add row to data collection
            line.clear(); // Clear buffer for next iteration
        }

        if data.is_empty() { // Check if any data was found
            return Err(ParseError::FormatError("No data rows found".to_string())); // Error if no data
        }

        Ok(data) // Return all parsed data rows
    }

    /// Build the final ElementValueData structure with ElementValue structs
    fn build_element_value_data(
        data: Vec<Vec<f64>>, // Raw numeric data from file
        indices: &ColumnIndices, // Column mapping information
        _column_names: &[String], // Original column names (not used)
    ) -> Result<ElementValueData, ParseError> {
        let num_elements = data.len(); // Number of elements equals number of data rows
        
        if num_elements == 0 { // Verify we have data to process
            return Err(ParseError::FormatError("No data found".to_string())); // Error if no data
        }

        let mut all_element_values = Vec::new(); // Collection of all ElementValue structs
        let mut element_value_type_info = Vec::new(); // Metadata for each value type
        let mut start_index = 0; // Track starting position for each value type in the collection

        // Process each value type group
        for (column_indices, vtype, dim) in &indices.value_columns { // Iterate through each identified value type
            let mut type_element_values = Vec::new(); // ElementValues for current value type
            
            // For each element (row in data), create an ElementValue struct
            for element_idx in 0..num_elements { // Process each element/row
                let row = &data[element_idx]; // Get data for current element
                let mut values = Vec::new(); // Values for current element and value type
                
                // Collect all components for this value type
                for &col_idx in column_indices { // Iterate through columns for this value type
                    if col_idx < row.len() { // Check if column exists in this row
                        values.push(row[col_idx]); // Add actual value from data
                    } else {
                        values.push(0.0); // Default value if column missing
                    }
                }
                
                // Create ElementValue struct
                type_element_values.push(ElementValue { // Create new ElementValue for this element and value type
                    id: element_idx, // Element identifier (row index)
                    values, // Component values for this value type
                });
            }
            
            if !type_element_values.is_empty() { // If we have valid ElementValues for this type
                // Add metadata for this value type
                element_value_type_info.push(ElementValueTypeInfo { // Create metadata for this value type
                    dimension: *dim, // Physical dimension (1=scalar, 3=vector, 6=tensor)
                    num_elements: num_elements, // Total number of elements
                    elementvalue_type: vtype.clone(), // Type of physical quantity
                    num_elementvalue_type: 1, // Static analysis = 1 type per value
                    start_index, // Starting position in all_element_values collection
                });
                
                // Add all ElementValue structs for this type to the main collection
                all_element_values.extend(type_element_values); // Add all ElementValues for this type
                start_index = all_element_values.len(); // Update start index for next value type
            }
        }

        Ok(ElementValueData { // Create final data structure
            element_values: all_element_values, // All ElementValue structs organized by type
            element_value_type_info, // Metadata describing the organization
        })
    }
}

/// Helper structure for column indices (simplified for static data)
struct ColumnIndices {
    value_columns: Vec<(Vec<usize>, ValueType, usize)>, // (column_indices, value_type, dimension) for each value type
}


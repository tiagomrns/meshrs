// Import standard library modules for file I/O operations
use std::fs::File;              // For opening files
use std::io::{self, BufRead, BufReader}; // For buffered reading of files line by line

//use crate::database::*;                     // Import mesh data structures and error types from database module

//use crate::lib::*;                  // Import mesh data structures and error types from lib module
use crate::structs_and_impls::*;                  // Import mesh data structures and error types from lib module
use crate::error::*;                     // Import mesh data structures and error types from error module

// Implementation block for ElementType enum - adds methods to the ElementType enum
impl ElementType {
    
    /// Convert Abaqus element type string to our internal ElementType enum
    /// Takes an Abaqus element name (like "C3D8", "S4R") and returns corresponding ElementType
    fn from_str_abaqus(s: &str) -> Option<ElementType> {
        // Convert input to uppercase for case-insensitive matching
        
        //1D elements (line elements like beams, trusses)
        if s.contains("C1D2") { // C1D2 = 2-node link element in Abaqus
            return Some(ElementType::Line);         // Maps to simple line element in VTK
        } else if s.contains("C1D3") { // C1D3 = 3-node link element (quadratic)
            return Some(ElementType::QuadraticEdge); // Maps to quadratic edge in VTK
        }

        // Check 3D elements first to avoid substring conflicts with 2D elements
        // (e.g., "C3D4" contains "D4" but should be treated as 3D, not 2D)
        else if s.contains("3D4") { // 4-node linear tetrahedron (most common 3D element)
            return Some(ElementType::Tetra);
        } else if s.contains("3D5") { // 5-node linear pyramid (has square base, triangular sides)
            return Some(ElementType::Pyramid);
        } else if s.contains("3D6") { // 6-node linear triangular prism (wedge-shaped)
            return Some(ElementType::Wedge);
        } else if s.contains("3D8") { // 8-node linear brick/hexahedron (cube-like)
            return Some(ElementType::Hexahedron);
        } else if s.contains("3D10") { // 10-node quadratic tetrahedron (curved edges)
            return Some(ElementType::QuadraticTetra);
        } else if s.contains("3D15") { // 15-node quadratic triangular prism
            return Some(ElementType::QuadraticWedge);
        } else if s.contains("3D20") { // 20-node quadratic brick (curved edges)
            return Some(ElementType::QuadraticHexahedron);
        } 

        // Check 2D elements (shells, membranes)
        // These can have various prefixes: E=general, EG=enhanced, S=shell, D=membrane
        else if s.contains("E3") || s.contains("EG3") || s.contains("S3") || s.contains("D3") { 
            // 3-node linear triangular elements
            return Some(ElementType::Triangle);
        } else if s.contains("E4") || s.contains("EG4") || s.contains("S4") || s.contains("D4") { 
            // 4-node bilinear quadrilateral elements
            return Some(ElementType::Quad);
        } else if s.contains("E6") || s.contains("EG6") || s.contains("S6") || s.contains("D6") { 
            // 6-node quadratic triangular elements
            return Some(ElementType::QuadraticTriangle);
        } else if s.contains("E8") || s.contains("EG8") || s.contains("S8") || s.contains("D8") { 
            // 8-node biquadratic quadrilateral elements
            return Some(ElementType::QuadraticQuad);
        }

        //Cylindrical Element types (specialized for cylindrical coordinate systems)
        else if s.contains("CL9") { // 9-node cylindrical prism
            return Some(ElementType::QuadraticWedge);    // Map to closest VTK equivalent
        } else if s.contains("CL12") { // 12-node cylindrical brick
            return Some(ElementType::QuadraticHexahedron);
        } else if s.contains("CL18") { // 18-node cylindrical prism (high-order)
            return Some(ElementType::BiquadraticQuadraticWedge);
        } else if s.contains("CL24") { // 24-node cylindrical brick (high-order)
            return Some(ElementType::BiquadraticQuadraticHexahedron);
        }

        //Axisymmetric Element types (for problems with rotational symmetry)
        else if s.contains("AX2") {   // 2-node axisymmetric line element
            return Some(ElementType::Line);
        } else if s.contains("AX3") { // 3-node axisymmetric triangle
            return Some(ElementType::Triangle);
        } else if s.contains("AX4") || s.contains("AXA4") { // 4-node axisymmetric quad
            return Some(ElementType::Quad);
        } else if s.contains("AX6") { // 6-node axisymmetric quadratic triangle
            return Some(ElementType::QuadraticTriangle);
        } else if s.contains("AX8") || s.contains("AXA8") { // 8-node axisymmetric quadratic quad
            return Some(ElementType::QuadraticQuad);
        }

        // If no match found, return None (will cause parsing error)
        None
    }
}

// Public structure for the Abaqus INP file parser
pub struct AbaqusInpParser;

// Implementation of parsing methods for AbaqusInpParser
impl AbaqusInpParser {
    /// Main entry point for parsing an Abaqus .inp file
    /// Opens the file and coordinates the entire parsing process
    pub fn parse_file(filename: &str) -> Result<MeshData, ParseError> {
        // Open the file for reading
        let abaqus_inp_file = File::open(filename)?;           // ? operator propagates any file opening errors
        // Wrap in BufReader for efficient line-by-line reading
        let abaqus_inp_reader = BufReader::new(abaqus_inp_file);
        // Parse the file content and get initial mesh data
        let mut mesh_data = Self::parse_abaqus_inp(abaqus_inp_reader)?;
        
        // Post-processing step 1: Convert from Abaqus 1-based indexing to VTK 0-based indexing
        Self::convert_to_zero_based_indexing(&mut mesh_data);
        
        // Post-processing step 2: Reorder nodes within elements to match VTK conventions
        // (Abaqus and VTK may have different node ordering for the same element type)
        Self::reorder_nodes_abaqus(&mut mesh_data)?;

        Ok(mesh_data)
    }

    /// Parse the actual content of an Abaqus INP file
    /// Processes the file line by line, extracting nodes and elements
    pub fn parse_abaqus_inp(abaqus_inp_reader: BufReader<File>) -> Result<MeshData, ParseError> {
        // Initialize empty mesh data structure
        let mut mesh_data = MeshData {
            dimension: 3,                // Default to 3D (will be updated based on actual node coordinates)
            num_nodes: 0,               // Will be set after parsing nodes
            min_node_index: 1,          // ABAQUS uses 1-based indexing initially
            nodes: Vec::new(),          // Empty vector to store parsed nodes
            num_eltypes: 0,             // Will be incremented as we find different element types
            elements: Vec::new(),       // Empty vector to store parsed elements
            element_type_info: Vec::new(), // Information about each element type found
        };

        // Read all lines into memory for easier processing
        let lines: Vec<String> = abaqus_inp_reader.lines().collect::<Result<Vec<_>, _>>()?;
        // Create iterator that tracks both line number and content
        let mut line_iter = lines.iter().enumerate();
        let mut element_id_counter = 0;  // Track element IDs across different element types
        let mut parsed_eltypes = 0;      // Count how many element types we've parsed

        // Main parsing loop - process each line
        while let Some((line_num, current_line)) = line_iter.next() {
            let trimmed = current_line.trim();  // Remove leading/trailing whitespace

            // Skip empty lines
            if trimmed.is_empty() {
                continue;
            }

            // Pattern match on line content to identify sections
            match trimmed {
                // Check if this line starts a NODE section (case-insensitive)
                tag if tag.to_uppercase().starts_with("*NODE") => {
                    // Collect remaining lines for node parsing
                    let remaining_lines: Vec<&String> = lines.iter().skip(line_num + 1).collect();
                    // Parse all nodes starting from the next line
                    let (nodes_parsed, lines_consumed) = Self::parse_nodes(
                        &remaining_lines,     // Lines to parse
                        &mut mesh_data.nodes, // Where to store parsed nodes
                    )?;
                    
                    mesh_data.num_nodes = nodes_parsed;  // Update node count
                    // Set mesh dimension based on first node's coordinate count
                    if !mesh_data.nodes.is_empty() {
                        mesh_data.dimension = mesh_data.nodes[0].coordinates.len();
                    }
                    
                    // Skip the lines we just processed
                    for _ in 0..lines_consumed {
                        line_iter.next();
                    }
                }

                // Check if this line starts an ELEMENT section
                tag if tag.to_uppercase().starts_with("*ELEMENT") => {
                    // Collect remaining lines for element parsing
                    let remaining_lines: Vec<&String> = lines.iter().skip(line_num + 1).collect();
                    let start_index = mesh_data.elements.len();  // Where this element type starts in elements array
                    // Parse elements of this specific type
                    let (elements, type_info, lines_consumed) = Self::parse_single_element_type(
                        &remaining_lines,     // Lines to parse
                        trimmed,             // The *ELEMENT header line (contains type info)
                        element_id_counter,  // Starting element ID
                        start_index,         // Starting index in elements array
                    )?;
                    
                    // Update counters and add parsed data to mesh
                    element_id_counter += elements.len();        // Update element ID counter
                    mesh_data.elements.extend(elements);         // Add elements to main vector
                    mesh_data.element_type_info.push(type_info); // Add type info
                    parsed_eltypes += 1;                        // Increment element type counter
                    mesh_data.num_eltypes = parsed_eltypes;     // Update total count
                    
                    // Skip the lines we just processed
                    for _ in 0..lines_consumed {
                        line_iter.next();
                    }
                }
                
                // Ignore all other lines (comments, other keywords, etc.)
                _ => {}
            }
        }

        Ok(mesh_data)
    }

    /// Parse node data from NODE section
    /// Extracts node IDs and coordinates until next keyword or end of file
    fn parse_nodes(
        lines: &[&String],          // Lines to parse
        nodes: &mut Vec<Node>,      // Where to store parsed nodes
    ) -> Result<(usize, usize), ParseError> {  // Returns (nodes_parsed, lines_consumed)
        let mut nodes_parsed = 0;   // Count successful node parses
        let mut lines_consumed = 0; // Count lines processed
        
        // Process lines until we hit a keyword (starting with *) or end of data
        for line in lines {
            let trimmed = line.trim();
            
            // Stop if we hit empty line or new keyword section
            if trimmed.is_empty() || trimmed.starts_with('*') {
                break;
            }
            
            lines_consumed += 1;  // Count this line as processed
            
            // Split line by commas to get node ID and coordinates
            let parts: Vec<&str> = trimmed.split(',').map(|s| s.trim()).collect();
            
            // Need at least node ID and one coordinate
            if parts.len() < 2 {
                continue;  // Skip malformed lines
            }
            
            // Parse node ID (first field)
            let node_id = parts[0].parse::<usize>()
                .map_err(|e| ParseError::NumberParseError(format!(
                    "Invalid node ID '{}': {}", parts[0], e
                )))?;
            
            // Parse coordinates (remaining fields)
            let coordinates: Result<Vec<f64>, _> = parts[1..]  // Skip first element (node ID)
                .iter()
                .map(|s| s.parse::<f64>())  // Convert each string to f64
                .collect();
            
            // Handle coordinate parsing errors
            let coords = coordinates.map_err(|e| ParseError::NumberParseError(format!(
                "Invalid coordinate in node {}: {}", node_id, e
            )))?;
            
            // Create and store the node
            nodes.push(Node {
                id: node_id,
                coordinates: coords,
            });
            
            nodes_parsed += 1;  // Count successful parse
        }
        
        Ok((nodes_parsed, lines_consumed))
    }

    /// Parse a single element type section (one *ELEMENT keyword)
    /// Extracts element type from header and parses all elements of that type
    fn parse_single_element_type(
        lines: &[&String],          // Lines to parse
        header_line: &str,          // The *ELEMENT header line
        start_element_id: usize,    // Starting element ID for this type
        start_index: usize,         // Starting index in elements array
    ) -> Result<(Vec<Element>, ElementTypeInfo, usize), ParseError> {
        
        // Extract Abaqus element type name from header (e.g., "C3D8", "S4R")
        let full_element_type = Self::extract_full_element_type_from_header(header_line)?;
        // Convert Abaqus element type to our internal ElementType enum
        let element_type = ElementType::from_str_abaqus((&full_element_type))
            .ok_or_else(|| ParseError::FormatError(format!(
                "Unknown element type: {}", &full_element_type
            )))?;

        let mut elements = Vec::new();      // Store parsed elements
        let mut lines_consumed = 0;        // Track lines processed
        let mut nodes_per_element = 0;     // Validate consistent node count per element
        
        // Process element data lines
        for line in lines {
            let trimmed = line.trim();
            
            // Stop at empty line or new keyword section
            if trimmed.is_empty() || trimmed.starts_with('*') {
                break;
            }
            
            lines_consumed += 1;
            
            // Split line by commas to get element ID and node connectivity
            let parts: Vec<&str> = trimmed.split(',').map(|s| s.trim()).collect();
            
            // Need at least element ID and one node
            if parts.len() < 2 {
                continue;  // Skip malformed lines
            }
            
            // Parse element ID (first field)
            let element_id = parts[0].parse::<usize>()
                .map_err(|e| ParseError::NumberParseError(format!(
                    "Invalid element ID '{}': {}", parts[0], e
                )))?;
            
            // Parse node connectivity (remaining fields)
            let node_ids: Result<Vec<usize>, _> = parts[1..]  // Skip element ID
                .iter()
                .map(|s| s.parse::<usize>())  // Convert each node ID string to number
                .collect();
            
            // Handle node ID parsing errors
            let nodes = node_ids.map_err(|e| ParseError::NumberParseError(format!(
                "Invalid node ID in element {}: {}", element_id, e
            )))?;
            
            // Validate that all elements of this type have same number of nodes
            if nodes_per_element == 0 {
                nodes_per_element = nodes.len();  // Set expected count from first element
            } else if nodes.len() != nodes_per_element {
                return Err(ParseError::FormatError(format!(
                    "Element {} has {} nodes, expected {}",
                    element_id,
                    nodes.len(),
                    nodes_per_element
                )));
            }
            
            // Create and store the element
            elements.push(Element {
                id: element_id,
                nodes,
            });
        }

        let num_elements = elements.len();  // Count elements of this type

        // Create type information structure
        let type_info = ElementTypeInfo { 
            element_type,      // The ElementType enum value
            num_elements,      // How many elements of this type
            nodes_per_element, // Nodes per element for this type
            start_index,       // Where this type starts in the main elements array
        };

        Ok((elements, type_info, lines_consumed))
    }

    /// Extract the Abaqus element type name from the *ELEMENT header line
    /// Looks for "TYPE=" and extracts the element type that follows
    fn extract_full_element_type_from_header(header_line: &str) -> Result<String, ParseError> {
        let upper_line = header_line.to_uppercase();  // Convert to uppercase for searching
        
        // Look for "TYPE=" in the header line
        if let Some(type_start) = upper_line.find("TYPE=") {
            let type_part = &header_line[type_start + 5..]; // Skip past "TYPE="
            // Extract element type name (up to next comma or end of line)
            let full_type_name = type_part.split(',').next().unwrap_or("").trim();
            
            // Validate that we found a non-empty type name
            if full_type_name.is_empty() {
                return Err(ParseError::FormatError(format!(
                    "Empty element type in header: {}", header_line
                )));
            }
            
            Ok(full_type_name.to_string())
        } else {
            // No TYPE= found in header
            Err(ParseError::FormatError(format!(
                "Could not extract element type from header: {}", header_line
            )))
        }
    }

    /// Convert all node and element IDs from Abaqus 1-based to VTK 0-based indexing
    /// Abaqus numbers nodes/elements starting from 1, VTK starts from 0
    fn convert_to_zero_based_indexing(mesh_data: &mut MeshData) {
        // Update the minimum node index to reflect 0-based indexing
        mesh_data.min_node_index = 0;
        
        // Convert all node IDs from 1-based to 0-based
        for node in &mut mesh_data.nodes {
            if node.id > 0 {      // Safety check to avoid underflow
                node.id -= 1;     // Subtract 1 to convert to 0-based
            }
        }
        
        // Convert element IDs and their node references to 0-based
        for element in &mut mesh_data.elements {
            if element.id > 0 {   // Convert element ID
                element.id -= 1;
            }
            
            // Convert all node references within this element
            for node_id in &mut element.nodes {
                if *node_id > 0 { // Convert each node reference
                    *node_id -= 1;
                }
            }
        }
    }
    
    /// Apply node reordering for all element types in the mesh
    /// Different software packages may order nodes differently for the same element type
    fn reorder_nodes_abaqus(mesh_data: &mut MeshData) -> Result<(), ParseError> {
        // Process each element type separately
        for type_info in &mesh_data.element_type_info {
            let start_idx = type_info.start_index;           // First element of this type
            let end_idx = start_idx + type_info.num_elements; // One past last element of this type
            
            // Get mutable slice of elements for this type only
            let elements_slice = &mut mesh_data.elements[start_idx..end_idx];
            
            // Apply type-specific reordering to this group of elements
            Self::reorder_nodes_abaqus_for_eltype(elements_slice, &type_info.element_type)?;
        }
        
        Ok(())
    }

    /// Reorder nodes within elements of a specific type to match VTK conventions
    /// Each element type may have different node ordering between Abaqus and VTK
    fn reorder_nodes_abaqus_for_eltype(elements: &mut [Element], element_type: &ElementType) -> Result<(), ParseError> {
        // Process each element of this type
        for element in elements.iter_mut() {
            let abaqus_nodes = element.nodes.clone(); // Save original Abaqus node ordering
            
            // Apply type-specific reordering (most Abaqus types match VTK, but some don't)
            match element_type {
                ElementType::Line => {
                    // Line elements: Abaqus and VTK use identical node ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::QuadraticEdge => {
                    // Quadratic edge: 3-node line element with middle node
                    // Abaqus: [end1, end2, middle] -> VTK: [end1, middle, end2]
                    element.nodes = vec![abaqus_nodes[0], abaqus_nodes[2], abaqus_nodes[1]];
                }
                ElementType::Triangle => {
                    // Triangle elements: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::QuadraticTriangle => {
                    // 6-node quadratic triangle: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::Quad => {
                    // 4-node quadrilateral: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::QuadraticQuad => {
                    // 8-node quadratic quad: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::Tetra => {
                    // 4-node tetrahedron: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::QuadraticTetra => {
                    // 10-node quadratic tetrahedron: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::Pyramid => {
                    // 5-node pyramid: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::Wedge => {
                    // 6-node wedge/prism: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::QuadraticWedge => {
                    if abaqus_nodes.len() == 15 {           // Standard 15-node quadratic wedge
                        element.nodes = abaqus_nodes;       // Same ordering as VTK
                    } else if abaqus_nodes.len() == 9 {     // Special case: 9-node cylindrical prism
                        // Need to pad with additional nodes or remap - this might need adjustment
                        // based on your specific VTK writer expectations
                        element.nodes = vec![
                            abaqus_nodes[0], abaqus_nodes[1], abaqus_nodes[2], abaqus_nodes[3],
                            abaqus_nodes[4], abaqus_nodes[5], abaqus_nodes[6], abaqus_nodes[7],
                            abaqus_nodes[8], abaqus_nodes[9], abaqus_nodes[10], abaqus_nodes[11],
                            abaqus_nodes[12], abaqus_nodes[13], abaqus_nodes[14]
                        ];
                    }
                }
                ElementType::BiquadraticQuadraticWedge => {
                    // 18-node biquadratic wedge: needs specific reordering
                    element.nodes = vec![
                            abaqus_nodes[0], abaqus_nodes[1], abaqus_nodes[2], abaqus_nodes[3],
                            abaqus_nodes[4], abaqus_nodes[5], abaqus_nodes[9], abaqus_nodes[10],
                            abaqus_nodes[11], abaqus_nodes[12], abaqus_nodes[13], abaqus_nodes[14],
                            abaqus_nodes[6], abaqus_nodes[7], abaqus_nodes[8], abaqus_nodes[15],
                            abaqus_nodes[16], abaqus_nodes[17]
                        ];
                }
                ElementType::Hexahedron => {
                    // 8-node hexahedron/brick: Abaqus and VTK use same ordering
                    element.nodes = abaqus_nodes;
                }
                ElementType::QuadraticHexahedron => {
                    if abaqus_nodes.len() == 20 {           // Standard 20-node quadratic hexahedron
                        element.nodes = abaqus_nodes;       // Same ordering as VTK
                    } else if abaqus_nodes.len() == 12 {    // Special case: 12-node cylindrical brick
                        // Pad to 20 nodes for VTK compatibility - might need adjustment
                        element.nodes = vec![
                            abaqus_nodes[0], abaqus_nodes[1], abaqus_nodes[2], abaqus_nodes[3],
                            abaqus_nodes[4], abaqus_nodes[5], abaqus_nodes[6], abaqus_nodes[7],
                            abaqus_nodes[8], abaqus_nodes[9], abaqus_nodes[10], abaqus_nodes[11],
                            abaqus_nodes[12], abaqus_nodes[13], abaqus_nodes[14], abaqus_nodes[15],
                            abaqus_nodes[16], abaqus_nodes[17], abaqus_nodes[18], abaqus_nodes[19]
                        ];
                    }
                }
                ElementType::BiquadraticQuadraticHexahedron => {
                    // 24-node biquadratic hexahedron: complex reordering needed
                    element.nodes = vec![
                            abaqus_nodes[0], abaqus_nodes[1], abaqus_nodes[2], abaqus_nodes[3],
                            abaqus_nodes[4], abaqus_nodes[5], abaqus_nodes[6], abaqus_nodes[7],
                            abaqus_nodes[12], abaqus_nodes[13], abaqus_nodes[14], abaqus_nodes[15],
                            abaqus_nodes[16], abaqus_nodes[17], abaqus_nodes[18], abaqus_nodes[19],
                            abaqus_nodes[8], abaqus_nodes[9], abaqus_nodes[10], abaqus_nodes[11],
                            abaqus_nodes[23], abaqus_nodes[21], abaqus_nodes[20], abaqus_nodes[22]
                        ];
                }

                // Handle unsupported element types
                _ => {
                    return Err(ParseError::FormatError(format!(
                        "Node reordering not implemented for Abaqus element type: {:?}", 
                        element_type
                    )));
                }
            }
        }
        
        Ok(())
    }
}
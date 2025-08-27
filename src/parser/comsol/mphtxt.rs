use std::fs::File;                          
//use std::intrinsics::typed_swap_nonoverlapping;
// For file handling (opening files)
use std::io::{self, BufRead, BufReader};    // For input/output operations and buffered reading (buffered reading is more efficient for line-by-line processing)


use crate::database::*;                     // Import mesh data structures and error types from database module

//use crate::database::ElementType;

impl ElementType {
    fn from_str_comsol(s: &str) -> Option<ElementType> {  // Converts string from input data to its element type correspond - Comsol mphtxt to xml vtk
        match s {
            
            "vtx" => Some(ElementType::Vertex), 
            "edg" => Some(ElementType::Line),  
            "edg2" => Some(ElementType::QuadraticEdge),  
            "tri" => Some(ElementType::Triangle),  
            "tri2" => Some(ElementType::QuadraticTriangle),   
            "quad" => Some(ElementType::Quad),    
            "quad2" => Some(ElementType::BiquadraticQuad),   
            "tet" => Some(ElementType::Tetra),  
            "tet2" => Some(ElementType::QuadraticTetra),   
            "pyr" => Some(ElementType::Pyramid),   
            "pyr2" => Some(ElementType::QuadraticPyramid),   //?
            "prism" => Some(ElementType::Wedge),  
            "prism2" => Some(ElementType::BiquadraticQuadraticWedge),   
            "hex" => Some(ElementType::Hexahedron),  
            "hex2" => Some(ElementType::TriquadraticHexahedron),  

            _ => None,   
        }
    }
}

pub struct MphtxtParser;  // Defines a structure to represent a mesh parser for MPHTXT files

impl MphtxtParser {
    // Public method to parse an MPHTXT file from disk
    // Takes a filename and returns MeshData or ParseError
    pub fn parse_file(filename: &str) -> Result<MeshData, ParseError> {
        let mphtxt_file = File::open(filename)?;                        // Opens file and ? propagates errors
        let mphtxt_reader = BufReader::new(mphtxt_file);    // Creates a buffered reader for efficient line reading
        Self::parse_mphtxt(mphtxt_reader)                                           // Delegates to mphtxt parsing function
    }

    // Public method to parse MPHTXT file 
    pub fn parse_mphtxt(mphtxt_reader: BufReader<File>) -> Result<MeshData, ParseError> {
        // Initialize empty mesh data structure
        let mut mesh_data = MeshData {
            dimension: 0,                       // Will be set when we find # sdim
            num_nodes: 0,                       // Will be set when we find # number of mesh vertices
            min_node_index: 0,                  // Will be set when we find # lowest mesh vertex index 
            nodes: Vec::new(),                  // Will be populated when we find mesh points
            num_eltypes: 0,                     // Will be set when we find # number of element types 
            elements: Vec::new(),               // Will be populated when we find elements
            element_type_info: Vec::new(),      // Will be populated with element type information 
        };

        // Split mphtxt file into lines and store them in vector for line-by-line processing
        let mut lines = mphtxt_reader.lines();         // .lines() splits the string into lines
        let mut element_id_counter = 0;                                 // Global element ID counter
        let mut parsed_eltypes = 0;                                     // Counter for parsed element types

        // Main parsing loop to go through each line
        while let Some(line) = lines.next() {
            let current_line = line?;             // Reads next line, ? propagates errors
            let trimmed = current_line.trim();      // Trims leading and trailing whitespaces from the line

            if trimmed.is_empty() {                       // Skips empty lines
                continue;
            }

            // Pattern match on different MPHTXT section headers. 
            // Here we identify what type of data we are dealing with to easily parse it with relevent functions and store it in the correct place in the mesh_data structure. 
            // Each match branch handles a different type of data section
            match trimmed {
                // Extract spatial dimension from lines like "3 # sdim"
                tag if tag.contains("# sdim") => {
                    mesh_data.dimension = Self::parse_first_number(tag, "spatial dimension")?;  // ? propagates parsing errors
                }

                // Extract number of nodes from lines like "2607 # number of mesh vertices"
                tag if tag.contains("# number of mesh vertices") => {
                    mesh_data.num_nodes = Self::parse_first_number(tag, "number of nodes")?;  // ? propagates parsing errors
                }

                // Extract lowest mesh vertex index from lines like "0 # lowest mesh vertex index"
                tag if tag.contains("# lowest mesh vertex index") => {
                    mesh_data.min_node_index = Self::parse_first_number(tag, "lowest mesh vertex index")?;  // ? propagates parsing errors
                }
                
                // Parse mesh points (vertex/node coordinates) section
                // Search for line "# Mesh vertex coordinates"  
                tag if tag.contains("# Mesh vertex coordinates") => {
                    
                    // Parse nodes 
                    Self::parse_nodes(
                        &mut lines,                              // Remaining lines (Lines iterator for line-by-line processing after the header # Mesh vertex coordinates)
                        mesh_data.dimension,                     // Spatial dimension
                        mesh_data.min_node_index,     // Starting node index (0 for 0-based indexing)
                        mesh_data.num_nodes,                     // Expected number of nodes
                        &mut mesh_data.nodes,     
                        )?;
                    }
                    
                
                // Extract number of element types from lines like "5 # number of element types"
                tag if tag.contains("# number of element types") => {
                    mesh_data.num_eltypes = Self::parse_first_number(tag, "number of element types")?;      // ? propagates parsing errors
                }

                tag if tag.contains("# Type #") => {
                    // Parse type ID from header
                    //let type_id = Self::parse_type_id(tag)?;                                                // # Type #"type_id"

                    let start_index = mesh_data.elements.len();                                             // Starting index for this element type in the main elements vector. It is equal to the current length of elements vector beause we will append new elements to it.
                    let (elements, type_info) = Self::parse_single_element_type(
                        &mut lines,                                                                                // Remaining lines (Lines iterator for line-by-line processing after the header # Type #)
                        //type_id,
                        element_id_counter,                                                       // Number of elements already parsed (element ID counter) to assign correct IDs to new elements 
                        start_index,
                    )?;
                    
                    // Updates mesh_data structure with parsed elements and type information
                    element_id_counter += elements.len();                                                           // Increments element ID counter by the number of elements parsed in this type
                    mesh_data.elements.extend(elements);                                                      // Appends parsed elements to the main elements vector
                    mesh_data.element_type_info.push(type_info);                                                    // Adds element type information to the mesh data structure element_type_info vector
                    parsed_eltypes += 1;                                                                            // Increments parsed element types counter
                }
                
                
                // If line doesn't match any known pattern, ignore it
                _ => {}
            }
            
        }

        // Verify we parsed all expected element types
        if parsed_eltypes != mesh_data.num_eltypes {
            return Err(ParseError::FormatError(format!(
                "Expected {} element types, parsed {}",
                mesh_data.num_eltypes, parsed_eltypes
            )));
        }


        Ok(mesh_data) // Ok() returns successfully parsed mesh data. Returns ParseError if any error occurs during parsing
    }


    // Parse node coordinates section
    // Returns nodes with coordinates in vector format and updates the mesh_data structure
    fn parse_nodes(
        lines: &mut impl Iterator<Item = io::Result<String>>,
        dimension: usize,
        min_index: usize,
        num_nodes: usize,
        nodes: &mut Vec<Node>,
    ) -> Result<(), ParseError> {
        for node_id in min_index..(min_index + num_nodes) {
            let line = lines.next().ok_or(ParseError::FormatError(     // Checks if there is a next line after the header
                "Unexpected end of node data".to_string(),
            ))??;

            let coords: Result<Vec<f64>, _> = line                           // Processes the line to extract coordinates, returns an error (< , _>) if parsing fails
                .split_whitespace()                             // Splits line into parts separated by whitespace
                .take(dimension)                         // Takes only the first number of parts equal to the dimension (e.g., first 3 coordinates for 3D)
                .map(|s| s.parse::<f64>()) // Parses each part as f64
                .collect();                                                         // Collects parsed coordinates into a vector

            let coordinates = coords?;                    // Collects parsed coordinates (coords) into a vector (coordinates) or returns an error if parsing fails
            if coordinates.len() != dimension {                     // Checks if the number of coordinates matches the expected dimension
                return Err(ParseError::FormatError(format!(
                    "Node {} has {} coordinates, expected {}",
                    node_id,
                    coordinates.len(),
                    dimension
                )));
            }

            nodes.push(Node {                                       // Appends a new Node to the nodes vector
                id: node_id,
                coordinates,
            });
        }
        Ok(())
    }


    // Parse a single element type section
    // This function reads the element type header, number of vertices per element, number of elements, and connectivity data.
    // It returns a vector of Element and ElementTypeInfo. It returns an error if the format is incorrect or if any parsing fails.
    fn parse_single_element_type(
        lines: &mut impl Iterator<Item = io::Result<String>>,
        // type_id: usize,
        start_element_id: usize,
        start_index: usize,
    ) -> Result<(Vec<Element>, ElementTypeInfo), ParseError> {
        // Read element type tag (e.g., "4 edg2 # type name")
        let type_tag_line = Self::next_non_empty_line(lines)?;      // Gets the next non-empty line which should be the element type tag
        let type_tag_line = type_tag_line.trim();                     // Trims leading and trailing whitespaces from the line
        
        if type_tag_line.is_empty() {                                       // Checks if the line is empty. Returns an error if it is.
            return Err(ParseError::FormatError(
                "Element type tag line is empty".to_string()
            ));
        }
        
        let parts: Vec<&str> = type_tag_line.split_whitespace().collect();  // Splits the line into parts based on whitespace and collects them into a vector
        
        // At least 2 parts are expected (number and element type name: i.e. "4 edg2")
        if parts.len() < 2 {                                                // Checks if the line has at least 2 parts. Returns an error if it doesn't.
            return Err(ParseError::FormatError(format!(
                "Invalid element type tag format: '{}'", type_tag_line
            )));
        }

        let element_type = ElementType::from_str_comsol(parts[1])  // Converts the second part of the line (element type name) to ElementType enum
            .ok_or_else(|| ParseError::FormatError(format!(                             // Returns an error if the conversion fails
                "Unknown element type: {}", parts[1]
            )))?;

        // Find "# number of vertices per element"
        let vertices_line = Self::find_line_containing(                                                                       // Finds the line containing "# number of vertices per element"
            lines,
            "# number of vertices per element",
        )?;
        let nodes_per_element = Self::parse_first_number(&vertices_line, "number of vertices per element")?;    // Parses the number of vertices per element from the line

        // Get number of elements
        let count_line = Self::next_non_empty_line(lines)?;                                                 // Gets the next non-empty line which should contain the number of elements
        let num_elements = Self::parse_first_number(&count_line, "number of elements")?;        // Parses the number of elements from the line

        // Find "# Elements" header
        Self::find_line_containing(lines, "# Elements")?;       // Finds the line containing "# Elements" header, which indicates the start of connectivity data

        // Parse connectivity data
        let mut elements = Vec::new();
        for elem_idx in 0..num_elements {
            let line = Self::next_non_empty_line(lines)?;                       // Gets the next non-empty line which should contain connectivity data for the element

            let node_ids: Result<Vec<usize>, _> = line
                .split_whitespace()                                  // Splits the line into parts based on whitespace
                .map(|s| s.parse::<usize>())    // Parses each part as usize
                .collect();                                                             // Collects parsed node IDs into a vector

            let nodes = node_ids?;
            if nodes.len() != nodes_per_element {  // Checks if the number of nodes matches the expected number of vertices per element
                return Err(ParseError::FormatError(format!(
                    "Element {} has {} nodes, expected {}",
                    elem_idx,
                    nodes.len(),
                    nodes_per_element
                )));
            }

            elements.push(Element {                     // Appends a new Element to the elements vector
                id: start_element_id + elem_idx,        // Assigns a unique ID to the element based on the starting element ID and the current index
                nodes,                                  // Stores the node IDs for this element
                // element_type: element_type.clone(),     // Clones the element type to avoid ownership issues
            });
        }

        // Reorder nodes for the entire block of elements based on element type
        Self::reorder_nodes_comsol(&mut elements, &element_type)?;

        // Create element type information
        let type_info = ElementTypeInfo { 
            element_type,
            num_elements,
            nodes_per_element,
            start_index,
        };

        Ok((elements, type_info))
    }

    // Helper to get next non-empty line
    fn next_non_empty_line(
        lines: &mut impl Iterator<Item = io::Result<String>>
    ) -> Result<String, ParseError> {
        loop {
            let line = lines.next()
                .ok_or(ParseError::FormatError("Unexpected end of file".to_string()))??;
            
            if !line.trim().is_empty() {
                return Ok(line);
            }
        }
    }

    // Helper to find a line containing a specific pattern
    fn find_line_containing(
        lines: &mut impl Iterator<Item = io::Result<String>>,
        pattern: &str,
    ) -> Result<String, ParseError> {
        loop {
            let line = lines.next().ok_or(ParseError::FormatError(format!(
                "Could not find line: {}",
                pattern
            )))??;
            if line.contains(pattern) {
                return Ok(line);
            }
        }
    }

    // Helper to parse first number from a line
    fn parse_first_number(tag: &str, context: &str) -> Result<usize, ParseError> {
        // Split the line and find the first token that can be parsed as a number
        for part in tag.split_whitespace() {
            if let Ok(num) = part.parse::<usize>() {
                return Ok(num);
            }
        }
        
        Err(ParseError::FormatError(format!(
        "Failed to parse number for '{}' in line: '{}'", 
        context, tag
    )))
    }

    /// Helper to adapt node order in element connectivity
    fn reorder_nodes_comsol(elements: &mut Vec<Element>, element_type: &ElementType) -> Result<(), ParseError> {
    for element in elements.iter_mut() {
        let comsol_nodes = element.nodes.clone(); // Save original COMSOL ordering
        
        match element_type {
            ElementType::Vertex => {
                // Vertex elements don't need reordering
                element.nodes = comsol_nodes;
            }
            ElementType::Line => {
                // Line elements: COMSOL and VTK have same ordering
                element.nodes = comsol_nodes; // Keep original ordering
            }
            ElementType::QuadraticEdge => {
                // Quadratic edge elements: COMSOL and VTK have same ordering
                element.nodes = comsol_nodes; // Keep original ordering
            }
            ElementType::Triangle => {
                // Triangle elements: COMSOL and VTK have same ordering
                element.nodes = comsol_nodes; // Keep original ordering
            }
            ElementType::QuadraticTriangle => {
                // Quadratic triangle elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[2],
                        comsol_nodes[3], comsol_nodes[5], comsol_nodes[4]
                    ];
            }
            ElementType::Quad => {
                // Quad elements
                element.nodes = vec![comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2]];
            }
            ElementType::BiquadraticQuad => {
                // Biquadratic quad elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4], comsol_nodes[7], comsol_nodes[8], comsol_nodes[5],
                        comsol_nodes[6]
                    ];
            }
            ElementType::Tetra => {
                // Tetrahedron elements: COMSOL and VTK have same ordering
                element.nodes = vec![comsol_nodes[0], comsol_nodes[1], comsol_nodes[2], comsol_nodes[3]];
            }
            ElementType::QuadraticTetra => {
                // Quadratic tetrahedron elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[2], comsol_nodes[3],
                        comsol_nodes[4], comsol_nodes[6], comsol_nodes[5], comsol_nodes[7],
                        comsol_nodes[8], comsol_nodes[9]
                    ];
            }
            ElementType::Pyramid => {
                // Pyramid elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4]
                    ];
            }
            ElementType::QuadraticPyramid => {   // ?? node padding or ignoring etc.
                // Quadratic pyramid: 13 nodes
                element.nodes = comsol_nodes; // Keep original for now, adjust as needed
            }
            ElementType::Wedge => {
                // Wedge (prism) elements: COMSOL and VTK have same ordering
                element.nodes = comsol_nodes; // Keep original for now, adjust as needed
            }
            ElementType::BiquadraticQuadraticWedge => {
                // Biquadratic quadratic wedge elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[2], comsol_nodes[3],
                        comsol_nodes[4], comsol_nodes[5], comsol_nodes[6], comsol_nodes[8],
                        comsol_nodes[7], comsol_nodes[15], comsol_nodes[17], comsol_nodes[16],
                        comsol_nodes[9], comsol_nodes[11], comsol_nodes[14], comsol_nodes[10],
                        comsol_nodes[13], comsol_nodes[12]
                    ];
            }
            ElementType::Hexahedron => {
                // Hexahedron elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4], comsol_nodes[5], comsol_nodes[7], comsol_nodes[6]
                    ];
            }
            ElementType::TriquadraticHexahedron => {
                // Triquadratic hexahedron elements
                element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4], comsol_nodes[5], comsol_nodes[7], comsol_nodes[6],
                        comsol_nodes[8], comsol_nodes[11], comsol_nodes[12], comsol_nodes[9],
                        comsol_nodes[22], comsol_nodes[25], comsol_nodes[26], comsol_nodes[23],
                        comsol_nodes[13], comsol_nodes[15], comsol_nodes[21], comsol_nodes[19],
                        comsol_nodes[16], comsol_nodes[18], comsol_nodes[14], comsol_nodes[20],
                        comsol_nodes[10], comsol_nodes[24], comsol_nodes[17]
                    ];
            }
            _ => {
                return Err(ParseError::FormatError(format!(
                    "Node reordering not implemented for element type: {:?}", 
                    element_type
                )));
            }
        }
    }
    
    Ok(())
    }
}
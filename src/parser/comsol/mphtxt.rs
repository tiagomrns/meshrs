// File handling for opening and reading mesh files
use std::fs::File;
// I/O operations: buffered reading is efficient for line-by-line processing of large mesh files
use std::io::{self, BufRead, BufReader};

// Import mesh data structures (MeshData, Node, Element, ElementType, ElementTypeInfo)
use crate::structs_and_impls::*;
// Import error types for mesh parsing (ParseError)
use crate::error::*;

impl ElementType {
    /// Converts COMSOL element type strings to internal ElementType enum
    /// 
    /// COMSOL .mphtxt files use specific strings to represent element types.
    /// This function maps those strings to our internal representation that's compatible with VTK format.
    /// 
    /// # Arguments
    /// * `s` - Element type string from COMSOL file (e.g., "tri", "tet", "hex2")
    /// 
    /// # Returns
    /// * `Some(ElementType)` if the string is recognized
    /// * `None` if the element type is unknown
    fn from_str_comsol(s: &str) -> Option<ElementType> {
        match s {
            // 0D elements: point/vertex
            "vtx" => Some(ElementType::Vertex),
            
            // 1D elements: lines (linear and quadratic)
            "edg" => Some(ElementType::Line),
            "edg2" => Some(ElementType::QuadraticEdge),
            
            // 2D elements: triangles (linear and quadratic)
            "tri" => Some(ElementType::Triangle),
            "tri2" => Some(ElementType::QuadraticTriangle),
            
            // 2D elements: quadrilaterals (linear and biquadratic)
            "quad" => Some(ElementType::Quad),
            "quad2" => Some(ElementType::BiquadraticQuad),
            
            // 3D elements: tetrahedra (linear and quadratic)
            "tet" => Some(ElementType::Tetra),
            "tet2" => Some(ElementType::QuadraticTetra),
            
            // 3D elements: pyramids (linear and quadratic)
            "pyr" => Some(ElementType::Pyramid),
            "pyr2" => Some(ElementType::QuadraticPyramid),
            
            // 3D elements: prisms/wedges (linear and biquadratic-quadratic)
            "prism" => Some(ElementType::Wedge),
            "prism2" => Some(ElementType::BiquadraticQuadraticWedge),
            
            // 3D elements: hexahedra (linear and triquadratic)
            "hex" => Some(ElementType::Hexahedron),
            "hex2" => Some(ElementType::TriquadraticHexahedron),

            // Unknown element type
            _ => None,
        }
    }
}

/// Parser for COMSOL .mphtxt mesh files
/// 
/// The .mphtxt format is a text-based representation of finite element meshes exported from COMSOL.
/// It contains:
/// - Spatial dimension (2D or 3D)
/// - Node/vertex coordinates
/// - Element connectivity information (which nodes form each element)
/// - Multiple element types in the same mesh
pub struct MphtxtParser;

impl MphtxtParser {
    /// Parses a COMSOL .mphtxt mesh file from disk
    /// 
    /// This is the main entry point for parsing. It opens the file and delegates
    /// to the internal parsing function.
    /// 
    /// # Arguments
    /// * `filename` - Path to the .mphtxt file
    /// 
    /// # Returns
    /// * `Ok(MeshData)` - Successfully parsed mesh data structure
    /// * `Err(ParseError)` - File I/O error or format parsing error
    /// 
    /// # Example
    /// ```
    /// let mesh = MphtxtParser::parse_file("mesh.mphtxt")?;
    /// ```
    pub fn parse_file(filename: &str) -> Result<MeshData, ParseError> {
        // Open the file (returns I/O error if file doesn't exist or can't be read)
        let mphtxt_file = File::open(filename)?;
        
        // Create a buffered reader for efficient line-by-line reading
        let mphtxt_reader = BufReader::new(mphtxt_file);
        
        // Parse the mphtxt format
        Self::parse_mphtxt(mphtxt_reader)
    }

    /// Parses COMSOL .mphtxt format from a buffered reader
    /// 
    /// This is the core parsing function that processes the .mphtxt file line by line.
    /// 
    /// # COMSOL .mphtxt File Format Structure:
    /// 1. Header with spatial dimension: "3 # sdim" (for 3D meshes)
    /// 2. Number of nodes: "2607 # number of mesh vertices"
    /// 3. Lowest node index: "0 # lowest mesh vertex index"
    /// 4. Node coordinates section: "# Mesh vertex coordinates" followed by coordinate lines
    /// 5. Number of element types: "5 # number of element types"
    /// 6. For each element type:
    ///    - Type header: "# Type #1", "# Type #2", etc.
    ///    - Element type tag: "4 edg2 # type name"
    ///    - Vertices per element: "3 # number of vertices per element"
    ///    - Number of elements: "450 # number of elements"
    ///    - Elements header: "# Elements"
    ///    - Connectivity data: one line per element with node indices
    /// 
    /// # Arguments
    /// * `mphtxt_reader` - Buffered reader positioned at start of file
    /// 
    /// # Returns
    /// * `Ok(MeshData)` - Complete parsed mesh with nodes and elements
    /// * `Err(ParseError)` - Format error or unexpected end of file
    pub fn parse_mphtxt(mphtxt_reader: BufReader<File>) -> Result<MeshData, ParseError> {
        // Initialize empty mesh data structure that will be populated during parsing
        let mut mesh_data = MeshData {
            dimension: 0,                       // Spatial dimension (2D or 3D) - set from "# sdim" line
            num_nodes: 0,                       // Total node count - set from "# number of mesh vertices" line
            min_node_index: 0,                  // Starting index for nodes (usually 0) - set from "# lowest mesh vertex index" line
            nodes: Vec::new(),                  // Node coordinates - populated from "# Mesh vertex coordinates" section
            num_eltypes: 0,                     // Number of different element types - set from "# number of element types" line
            elements: Vec::new(),               // All elements from all types - populated from "# Type" sections
            element_type_info: Vec::new(),      // Metadata for each element type - populated from "# Type" sections
        };

        // Convert buffered reader into an iterator over lines
        let mut lines = mphtxt_reader.lines();
        
        // Global element ID counter - ensures unique IDs across all element types
        let mut element_id_counter = 0;
        
        // Track how many element types we've successfully parsed
        let mut parsed_eltypes = 0;

        // Main parsing loop: process file line by line
        while let Some(line) = lines.next() {
            // Extract line content, propagate I/O errors with ?
            let current_line = line?;
            
            // Remove leading/trailing whitespace
            let trimmed = current_line.trim();

            // Skip empty lines (common in .mphtxt format for readability)
            if trimmed.is_empty() {
                continue;
            }

            // Pattern matching on COMSOL section headers
            // Each branch handles a different section of the .mphtxt format
            match trimmed {
                // Parse spatial dimension (2D or 3D)
                // Example line: "3 # sdim" means 3D mesh
                tag if tag.contains("# sdim") => {
                    mesh_data.dimension = Self::parse_first_number(tag, "spatial dimension")?;
                }

                // Parse total number of nodes/vertices
                // Example line: "2607 # number of mesh vertices" means 2607 nodes
                tag if tag.contains("# number of mesh vertices") => {
                    mesh_data.num_nodes = Self::parse_first_number(tag, "number of nodes")?;
                }

                // Parse the starting index for node numbering (usually 0 for 0-based indexing)
                // Example line: "0 # lowest mesh vertex index"
                tag if tag.contains("# lowest mesh vertex index") => {
                    mesh_data.min_node_index = Self::parse_first_number(tag, "lowest mesh vertex index")?;
                }
                
                // Parse the node coordinates section
                // This section starts with "# Mesh vertex coordinates" followed by one line per node
                // Each line contains space-separated coordinates (x, y for 2D or x, y, z for 3D)
                tag if tag.contains("# Mesh vertex coordinates") => {
                    Self::parse_nodes(
                        &mut lines,                    // Iterator over remaining lines
                        mesh_data.dimension,           // How many coordinates per node (2 or 3)
                        mesh_data.min_node_index,      // Starting node ID (usually 0)
                        mesh_data.num_nodes,           // How many nodes to read
                        &mut mesh_data.nodes,          // Where to store parsed nodes
                    )?;
                }
                
                // Parse total number of different element types in the mesh
                // Example line: "5 # number of element types" means 5 different types (e.g., triangles, quads, tets, etc.)
                tag if tag.contains("# number of element types") => {
                    mesh_data.num_eltypes = Self::parse_first_number(tag, "number of element types")?;
                }

                // Parse a single element type section
                // Each section starts with "# Type #N" where N is the type number (1, 2, 3, ...)
                tag if tag.contains("# Type #") => {
                    // Calculate where this element type's elements will start in the global elements vector
                    let start_index = mesh_data.elements.len();
                    
                    // Parse this element type's data (type, nodes per element, connectivity)
                    let (elements, type_info) = Self::parse_single_element_type(
                        &mut lines,                 // Iterator over remaining lines
                        element_id_counter,         // Current global element ID for unique numbering
                        start_index,                // Where these elements start in global vector
                    )?;
                    
                    // Update global element ID counter for next element type
                    element_id_counter += elements.len();
                    
                    // Append this type's elements to the global elements vector
                    mesh_data.elements.extend(elements);
                    
                    // Store metadata about this element type (count, nodes per element, etc.)
                    mesh_data.element_type_info.push(type_info);
                    
                    // Track that we've successfully parsed one more element type
                    parsed_eltypes += 1;
                }
                
                // Ignore any other lines (comments, empty sections, etc.)
                _ => {}
            }
        }

        // Validation: ensure we parsed the expected number of element types
        if parsed_eltypes != mesh_data.num_eltypes {
            return Err(ParseError::FormatError(format!(
                "Expected {} element types, parsed {}",
                mesh_data.num_eltypes, parsed_eltypes
            )));
        }

        // Return successfully parsed mesh data
        Ok(mesh_data)
    }


    /// Parses the node coordinates section of a .mphtxt file
    /// 
    /// After the "# Mesh vertex coordinates" header, the file contains one line per node
    /// with space-separated coordinates. For 2D meshes: "x y", for 3D meshes: "x y z"
    /// 
    /// # Arguments
    /// * `lines` - Iterator over remaining file lines
    /// * `dimension` - Number of coordinates per node (2 for 2D, 3 for 3D)
    /// * `min_index` - Starting node ID (usually 0)
    /// * `num_nodes` - Expected number of nodes to parse
    /// * `nodes` - Vector to store parsed nodes
    /// 
    /// # Returns
    /// * `Ok(())` - All nodes parsed successfully
    /// * `Err(ParseError)` - Unexpected end of file or coordinate parsing error
    fn parse_nodes(
        lines: &mut impl Iterator<Item = io::Result<String>>,
        dimension: usize,
        min_index: usize,
        num_nodes: usize,
        nodes: &mut Vec<Node>,
    ) -> Result<(), ParseError> {
        // Read coordinates for each node
        for node_id in min_index..(min_index + num_nodes) {
            // Get next line, return error if file ends prematurely
            let line = lines.next().ok_or(ParseError::FormatError(
                "Unexpected end of node data".to_string(),
            ))??;  // Double ? handles both Option and Result

            // Parse coordinates from the line
            let coords: Result<Vec<f64>, _> = line
                .split_whitespace()              // Split by spaces/tabs
                .take(dimension)                 // Take only dimension number of values (2 for 2D, 3 for 3D)
                .map(|s| s.parse::<f64>())       // Convert each string to f64
                .collect();                      // Collect into Vec<f64> or return parse error

            // Extract coordinates or return error if parsing failed
            let coordinates = coords?;
            
            // Validate: ensure we got the expected number of coordinates
            if coordinates.len() != dimension {
                return Err(ParseError::FormatError(format!(
                    "Node {} has {} coordinates, expected {}",
                    node_id,
                    coordinates.len(),
                    dimension
                )));
            }

            // Add the node to our collection
            nodes.push(Node {
                id: node_id,
                coordinates,
            });
        }
        Ok(())
    }


    /// Parses a single element type section from the .mphtxt file
    /// 
    /// Each element type section has this structure:
    /// ```text
    /// # Type #N
    /// 4 edg2 # type name
    /// 3 # number of vertices per element
    /// 450 # number of elements
    /// # Elements
    /// 0 1 2    <- connectivity line (node indices for element 0)
    /// 1 3 4    <- connectivity line (node indices for element 1)
    /// ...
    /// ```
    /// 
    /// # Arguments
    /// * `lines` - Iterator over remaining file lines
    /// * `start_element_id` - Global element ID to assign to first element of this type
    /// * `start_index` - Index in global elements vector where this type's elements begin
    /// 
    /// # Returns
    /// * `Ok((elements, type_info))` - Parsed elements and metadata
    /// * `Err(ParseError)` - Format error or unexpected end of file
    fn parse_single_element_type(
        lines: &mut impl Iterator<Item = io::Result<String>>,
        start_element_id: usize,
        start_index: usize,
    ) -> Result<(Vec<Element>, ElementTypeInfo), ParseError> {
        // Read the element type tag line (e.g., "4 edg2 # type name")
        let type_tag_line = Self::next_non_empty_line(lines)?;
        let type_tag_line = type_tag_line.trim();
        
        // Validate that we got a non-empty line
        if type_tag_line.is_empty() {
            return Err(ParseError::FormatError(
                "Element type tag line is empty".to_string()
            ));
        }
        
        // Split line into parts: ["4", "edg2", "#", "type", "name"]
        let parts: Vec<&str> = type_tag_line.split_whitespace().collect();
        
        // Need at least number and type name (e.g., "4 edg2")
        if parts.len() < 2 {
            return Err(ParseError::FormatError(format!(
                "Invalid element type tag format: '{}'", type_tag_line
            )));
        }

        // Convert COMSOL element type string to our internal enum (e.g., "edg2" -> QuadraticEdge)
        let element_type = ElementType::from_str_comsol(parts[1])
            .ok_or_else(|| ParseError::FormatError(format!(
                "Unknown element type: {}", parts[1]
            )))?;

        // Find and parse the "# number of vertices per element" line
        let vertices_line = Self::find_line_containing(
            lines,
            "# number of vertices per element",
        )?;
        let nodes_per_element = Self::parse_first_number(&vertices_line, "number of vertices per element")?;

        // Get the number of elements of this type (next non-empty line)
        let count_line = Self::next_non_empty_line(lines)?;
        let num_elements = Self::parse_first_number(&count_line, "number of elements")?;

        // Find the "# Elements" header that precedes the connectivity data
        Self::find_line_containing(lines, "# Elements")?;

        // Parse connectivity data: one line per element with node indices
        let mut elements = Vec::new();
        for elem_idx in 0..num_elements {
            // Get connectivity line
            let line = Self::next_non_empty_line(lines)?;

            // Parse node indices from the line
            let node_ids: Result<Vec<usize>, _> = line
                .split_whitespace()              // Split by spaces
                .map(|s| s.parse::<usize>())     // Parse each as integer
                .collect();                      // Collect into Vec<usize>

            let nodes = node_ids?;
            
            // Validate: ensure we got the expected number of nodes
            if nodes.len() != nodes_per_element {
                return Err(ParseError::FormatError(format!(
                    "Element {} has {} nodes, expected {}",
                    elem_idx,
                    nodes.len(),
                    nodes_per_element
                )));
            }

            // Create element with unique global ID
            elements.push(Element {
                id: start_element_id + elem_idx,
                nodes,
            });
        }

        // Reorder node indices to match VTK format (COMSOL and VTK use different orderings)
        Self::reorder_nodes_comsol(&mut elements, &element_type)?;

        // Create metadata for this element type
        let type_info = ElementTypeInfo {
            element_type,
            num_elements,
            nodes_per_element,
            start_index,
        };

        Ok((elements, type_info))
    }

    /// Gets the next non-empty line from the file
    /// 
    /// Skips over empty lines and whitespace-only lines.
    /// 
    /// # Arguments
    /// * `lines` - Iterator over file lines
    /// 
    /// # Returns
    /// * `Ok(String)` - Next non-empty line
    /// * `Err(ParseError)` - Unexpected end of file
    fn next_non_empty_line(
        lines: &mut impl Iterator<Item = io::Result<String>>
    ) -> Result<String, ParseError> {
        loop {
            // Get next line, return error if file ends
            let line = lines.next()
                .ok_or(ParseError::FormatError("Unexpected end of file".to_string()))??;
            
            // Return line if it contains non-whitespace characters
            if !line.trim().is_empty() {
                return Ok(line);
            }
            // Otherwise continue looping to find next non-empty line
        }
    }

    /// Finds and returns the first line containing a specific pattern
    /// 
    /// Continues reading lines until one contains the search pattern.
    /// 
    /// # Arguments
    /// * `lines` - Iterator over file lines
    /// * `pattern` - Text to search for in lines
    /// 
    /// # Returns
    /// * `Ok(String)` - First line containing the pattern
    /// * `Err(ParseError)` - Pattern not found before end of file
    fn find_line_containing(
        lines: &mut impl Iterator<Item = io::Result<String>>,
        pattern: &str,
    ) -> Result<String, ParseError> {
        loop {
            // Get next line, return error if file ends before finding pattern
            let line = lines.next().ok_or(ParseError::FormatError(format!(
                "Could not find line: {}",
                pattern
            )))??;
            
            // Return line if it contains the search pattern
            if line.contains(pattern) {
                return Ok(line);
            }
            // Otherwise continue searching
        }
    }

    /// Parses the first number found in a line
    /// 
    /// Splits the line by whitespace and returns the first token that can be parsed as a number.
    /// 
    /// # Arguments
    /// * `tag` - Line to parse
    /// * `context` - Description of what we're parsing (for error messages)
    /// 
    /// # Returns
    /// * `Ok(usize)` - First number found in the line
    /// * `Err(ParseError)` - No valid number found
    /// 
    /// # Example
    /// ```
    /// parse_first_number("3 # sdim", "spatial dimension") -> Ok(3)
    /// parse_first_number("2607 # number of mesh vertices", "nodes") -> Ok(2607)
    /// ```
    fn parse_first_number(tag: &str, context: &str) -> Result<usize, ParseError> {
        // Iterate through space-separated tokens
        for part in tag.split_whitespace() {
            // Try to parse token as number
            if let Ok(num) = part.parse::<usize>() {
                return Ok(num);
            }
        }
        
        // No valid number found
        Err(ParseError::FormatError(format!(
            "Failed to parse number for '{}' in line: '{}'", 
            context, tag
        )))
    }

    /// Reorders element node indices from COMSOL format to VTK format
    /// 
    /// COMSOL and VTK use different conventions for ordering the nodes in an element.
    /// This function converts COMSOL's node ordering to VTK's expected ordering.
    /// 
    /// For example, a quad in COMSOL might be ordered: [0,1,2,3]
    /// But VTK expects: [0,1,3,2]
    /// 
    /// This ensures compatibility when exporting to VTK formats.
    /// 
    /// # Arguments
    /// * `elements` - Mutable reference to elements to reorder
    /// * `element_type` - Type of elements (determines reordering pattern)
    /// 
    /// # Returns
    /// * `Ok(())` - Reordering successful
    /// * `Err(ParseError)` - Reordering not implemented for this element type
    fn reorder_nodes_comsol(elements: &mut Vec<Element>, element_type: &ElementType) -> Result<(), ParseError> {
        // Process each element of this type
        for element in elements.iter_mut() {
            // Save original COMSOL node ordering
            let comsol_nodes = element.nodes.clone();
            
            // Apply element-type-specific reordering
            match element_type {
                // Vertex elements: single node, no reordering needed
                ElementType::Vertex => {
                    element.nodes = comsol_nodes;
                }
                
                // Line elements: 2 nodes, COMSOL and VTK use same ordering
                ElementType::Line => {
                    element.nodes = comsol_nodes;
                }
                
                // Quadratic edge elements: 3 nodes, same ordering
                ElementType::QuadraticEdge => {
                    element.nodes = comsol_nodes;
                }
                
                // Triangle elements: 3 nodes, same ordering
                ElementType::Triangle => {
                    element.nodes = comsol_nodes;
                }
                
                // Quadratic triangle elements: 6 nodes, middle edge nodes reordered
                // COMSOL: [corner0, corner1, corner2, edge0-1, edge1-2, edge2-0]
                // VTK:    [corner0, corner1, corner2, edge0-1, edge2-0, edge1-2]
                ElementType::QuadraticTriangle => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[2],
                        comsol_nodes[3], comsol_nodes[5], comsol_nodes[4]
                    ];
                }
                
                // Quad elements: 4 nodes, diagonal pair swapped
                // COMSOL: [0, 1, 2, 3] (counter-clockwise from one corner)
                // VTK:    [0, 1, 3, 2] (different winding order)
                ElementType::Quad => {
                    element.nodes = vec![comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2]];
                }
                
                // Biquadratic quad elements: 9 nodes (4 corners, 4 edge midpoints, 1 center)
                // Reorder edge midpoints and center node
                ElementType::BiquadraticQuad => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4], comsol_nodes[7], comsol_nodes[8], comsol_nodes[5],
                        comsol_nodes[6]
                    ];
                }
                
                // Tetrahedron elements: 4 nodes, same ordering
                ElementType::Tetra => {
                    element.nodes = vec![comsol_nodes[0], comsol_nodes[1], comsol_nodes[2], comsol_nodes[3]];
                }
                
                // Quadratic tetrahedron elements: 10 nodes (4 corners + 6 edge midpoints)
                // Reorder edge midpoints to match VTK convention
                ElementType::QuadraticTetra => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[2], comsol_nodes[3],
                        comsol_nodes[4], comsol_nodes[6], comsol_nodes[5], comsol_nodes[7],
                        comsol_nodes[8], comsol_nodes[9]
                    ];
                }
                
                // Pyramid elements: 5 nodes (4-sided base + 1 apex)
                // Reorder base nodes to match VTK quad ordering
                ElementType::Pyramid => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4]
                    ];
                }
                
                // Quadratic pyramid elements: 14 nodes in COMSOL -> 13 nodes in VTK
                // Note: One node is dropped (middle node) to match VTK format
                ElementType::QuadraticPyramid => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4], comsol_nodes[5], comsol_nodes[8], comsol_nodes[9],
                        comsol_nodes[6], comsol_nodes[10], comsol_nodes[11], comsol_nodes[13],
                        comsol_nodes[12]
                    ];
                }
                
                // Wedge (prism) elements: 6 nodes (2 triangular faces)
                // COMSOL and VTK use same ordering
                ElementType::Wedge => {
                    element.nodes = comsol_nodes;
                }
                
                // Biquadratic quadratic wedge elements: 18 nodes
                // Reorder edge and face midpoints
                ElementType::BiquadraticQuadraticWedge => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[2], comsol_nodes[3],
                        comsol_nodes[4], comsol_nodes[5], comsol_nodes[6], comsol_nodes[8],
                        comsol_nodes[7], comsol_nodes[15], comsol_nodes[17], comsol_nodes[16],
                        comsol_nodes[9], comsol_nodes[11], comsol_nodes[14], comsol_nodes[10],
                        comsol_nodes[13], comsol_nodes[12]
                    ];
                }
                
                // Hexahedron elements: 8 nodes (brick/cube)
                // Reorder to match VTK's expected node sequence
                ElementType::Hexahedron => {
                    element.nodes = vec![
                        comsol_nodes[0], comsol_nodes[1], comsol_nodes[3], comsol_nodes[2],
                        comsol_nodes[4], comsol_nodes[5], comsol_nodes[7], comsol_nodes[6]
                    ];
                }
                
                // Triquadratic hexahedron elements: 27 nodes (8 corners + 12 edge midpoints + 6 face centers + 1 volume center)
                // Complex reordering to match VTK convention
                ElementType::TriquadraticHexahedron => {
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
                
                // Unknown element type - return error
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
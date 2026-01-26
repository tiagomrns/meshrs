// VTK I/O library imports for reading and writing VTK files
use vtkio::model::*;
use vtkio::Error as WriterError;

// Standard library imports
use std::fs;

// Internal module imports
use crate::structs_and_impls::*;
use crate::error::*;

/// VTUWriter is responsible for converting mesh data structures to VTK Unstructured Grid (VTU) format.
/// VTU files are XML-based files that can be visualized in ParaView and other VTK-compatible tools.
/// 
/// The writer supports:
/// - Basic mesh geometry (nodes and elements)
/// - Field data (displacement, velocity, stress, etc.)
/// - Quality metrics and Gaussian point analysis
pub struct VTUWriter;

impl VTUWriter {

    /// Writes basic mesh geometry to a VTU file without any field data.
    /// 
    /// # Arguments
    /// * `mesh_data` - The mesh containing nodes and elements
    /// * `output_path` - File path where the VTU file will be written
    /// 
    /// # Returns
    /// * `Ok(())` on success
    /// * `Err(WriterError)` if file writing or VTK serialization fails
    pub fn write_vtu(
        mesh_data: &MeshData,
        output_path: &str,
    ) -> Result<(), WriterError> {

        // Buffer to hold the serialized VTU XML data before writing to file
        let mut vtu = Vec::new();

        // Step 1: Convert node coordinates to VTK points format
        // VTK requires 3D coordinates (x, y, z), so we pad 1D/2D meshes with zeros
        let points_data: Vec<f64> = mesh_data
            .nodes
            .iter()
            .flat_map(|node| {
                node.coordinates
                    .iter()
                    .copied()
                    .chain(std::iter::repeat(0.0)) // Pad with zeros for missing dimensions
                    .take(3)                       // Ensure exactly 3 components (x, y, z)
            })
            .collect();

        // Step 2: Pre-allocate vectors for cell topology data
        // Calculate total size to avoid reallocations during iteration
        let total_connectivity: usize = mesh_data.elements.iter().map(|e| e.nodes.len()).sum();

        // Connectivity: flat list of node IDs for all elements
        let mut connectivity = Vec::with_capacity(total_connectivity);
        // Offsets: cumulative count of nodes at the end of each element
        let mut offsets = Vec::with_capacity(mesh_data.elements.len());
        // Cell types: VTK cell type for each element (initialized as Vertex, updated later)
        let mut cell_types = vec![CellType::Vertex; mesh_data.elements.len()];
        let mut current_offset = 0;

        // Step 3: Build connectivity and offset arrays
        // VTK uses a flat connectivity array with offsets to delineate elements
        for element in &mesh_data.elements {
            // Append this element's node IDs to the connectivity array
            connectivity.extend(element.nodes.iter().map(|&id| id as u64));
            // Update the cumulative offset after adding this element's nodes
            current_offset += element.nodes.len() as u64;
            offsets.push(current_offset);
        }

        // Step 4: Assign VTK cell types based on element type groups
        // Each group defines a range of elements with the same element type
        for group in &mesh_data.element_type_info {
            // Convert our element type to VTK's cell type enumeration
            let vtk_type = group.element_type.eltype_vtk();
            
            // Assign the VTK cell type to all elements in this group's range
            for i in group.start_index..(group.start_index + group.num_elements) {
                if i < cell_types.len() {
                    cell_types[i] = vtk_type;
                } else {
                    // Return error instead of panicking for bounds violations
                    return Err(WriterError::IO(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Element index {} out of bounds (total elements: {})", i, cell_types.len())
                    )));
                }
            }
        }
        
        // Validation: Ensure consistency of element-related data structures
        println!("Total elements: {}", mesh_data.elements.len());
        println!("Cell types count: {}", cell_types.len());
        println!("Offsets count: {}", offsets.len());
        assert_eq!(cell_types.len(), mesh_data.elements.len());
        assert_eq!(offsets.len(), mesh_data.elements.len());

        // Step 5: Construct the VTK data structure
        Vtk {
            version: Version { major: 2, minor: 2 },   // VTK XML format version (compatible with ParaView)
            title: String::new(),                      // Optional title for the dataset
            byte_order: ByteOrder::LittleEndian,       // Byte order for binary data sections
            file_path: None,                           // File path is set during write operation
            data: DataSet::inline(UnstructuredGridPiece {
                points: IOBuffer::F64(points_data),    // Node coordinates as 64-bit floats
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity: connectivity,     // Flat array of node IDs
                        offsets: offsets,               // Cumulative node counts per element
                    },
                    types: cell_types,                  // VTK cell type for each element
                },
                data: Attributes {
                    ..Default::default()                // No point or cell data attributes
                },
            }),
        }.write_xml(&mut vtu)?;                         // Serialize to XML format

        // Step 6: Write the serialized XML data to the output file
        fs::write(output_path, &vtu)?;

        Ok(())
    }

    /// Writes mesh geometry along with nodal field values to a VTU file.
    /// Field values are attached as point data attributes.
    /// 
    /// # Arguments
    /// * `mesh_data` - The mesh containing nodes and elements
    /// * `node_value_data` - Field values at nodes (displacement, velocity, stress, etc.)
    /// * `output_path` - File path where the VTU file will be written
    /// 
    /// # Returns
    /// * `Ok(())` on success
    /// * `Err(WriterError)` if file writing or VTK serialization fails
    pub fn write_vtu_with_field_values(
        mesh_data: &MeshData,
        node_value_data: &NodeValueData,
        output_path: &str,
    ) -> Result<(), WriterError> {

        // Buffer to hold the serialized VTU XML data
        let mut vtu = Vec::new();

        /// Helper function to determine if a value type represents a vector quantity
        /// Vector quantities (displacement, velocity, acceleration) need padding to 3D
        fn is_vector_quantity(vtype: &ValueType) -> bool {
            match vtype {
                ValueType::Displacement | 
                ValueType::Velocity | 
                ValueType::Acceleration => true,
                _ => false, // Scalars and tensors handled differently
            }
        }

        // Step 1: Convert node coordinates to VTK points format
        let points_data: Vec<f64> = mesh_data
            .nodes
            .iter()
            .flat_map(|node| {
                node.coordinates
                    .iter()
                    .copied()
                    .chain(std::iter::repeat(0.0)) // Pad with zeros for missing dimensions
                    .take(3)                       // Ensure exactly 3 components (x, y, z)
            })
            .collect();

        // Step 2: Pre-allocate vectors for cell topology data
        let total_connectivity: usize = mesh_data.elements.iter().map(|e| e.nodes.len()).sum();

        let mut connectivity = Vec::with_capacity(total_connectivity);
        let mut offsets = Vec::with_capacity(mesh_data.elements.len());
        let mut cell_types = vec![CellType::Vertex; mesh_data.elements.len()];
        let mut current_offset = 0;

        // Step 3: Build connectivity and offset arrays
        for element in &mesh_data.elements {
            connectivity.extend(element.nodes.iter().map(|&id| id as u64));
            current_offset += element.nodes.len() as u64;
            offsets.push(current_offset);
        }

        // Step 4: Assign VTK cell types based on element type groups
        for group in &mesh_data.element_type_info {
            let vtk_type = group.element_type.eltype_vtk();
            
            // Assign the VTK cell type to all elements in this group's range
            for i in group.start_index..(group.start_index + group.num_elements) {
                if i < cell_types.len() {
                    cell_types[i] = vtk_type;
                } else {
                    return Err(WriterError::IO(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Element index {} out of bounds (total elements: {})", i, cell_types.len())
                    )));
                }
            }
        }

    // Step 5: Process node value data and convert to VTK point attributes
    // Each value type (displacement, stress, etc.) becomes a separate attribute

    let mut point_attributes = Vec::new();

    for type_info in &node_value_data.node_value_type_info {
        
        // Extract the slice of node values corresponding to this value type
        let start_idx = type_info.start_index;
        let end_idx = start_idx + type_info.num_nodes;
        
        // Flatten the values from all nodes of this type into a single array
        let mut raw_values = Vec::with_capacity(type_info.num_nodes * type_info.dimension);
        if end_idx <= node_value_data.node_values.len() {
            for node_value in &node_value_data.node_values[start_idx..end_idx] {
                raw_values.extend_from_slice(&node_value.values);
            }
        }

        // Get the human-readable attribute name for ParaView display
        let attr_name = ValueType::get_attribute_name(&type_info.nodevalue_type);
        // Determine if this value type needs vector padding
        let is_vector = is_vector_quantity(&type_info.nodevalue_type);

        // Create VTK attribute based on dimension and value type
        // VTK expects specific formats: scalars (1D), vectors (3D), tensors (9D)
        let attribute = match type_info.dimension {
            1 => {
                if is_vector {
                    // 1D vector (e.g., 1D displacement) -> pad to 3D: u -> (u, 0, 0)
                    let padded: Vec<f64> = raw_values.iter()
                        .flat_map(|&val| [val, 0.0, 0.0])
                        .collect();
                    Attribute::vectors(attr_name).with_data(IOBuffer::F64(padded))
                } else {
                    // True scalar value (e.g., temperature, pressure)
                    Attribute::scalars(attr_name, 1).with_data(IOBuffer::F64(raw_values))
                }
            },
            2 => {
                if is_vector {
                    // 2D vector -> pad to 3D: (u, v) -> (u, v, 0)
                    let padded: Vec<f64> = raw_values.chunks(2)
                        .flat_map(|c| [c[0], c[1], 0.0])
                        .collect();
                    Attribute::vectors(attr_name).with_data(IOBuffer::F64(padded))
                } else {
                    // Generic 2D field data (not a vector)
                    Attribute::field(attr_name).with_data(IOBuffer::F64(raw_values))
                }
            },
            3 => {
                // 3D vector - no padding needed
                Attribute::vectors(attr_name).with_data(IOBuffer::F64(raw_values))
            },
            4 => {
                // 2D stress/strain tensor (xx, yy, xy, yx) -> pad to 3D 9-component tensor
                // VTK tensor format: [xx, xy, xz, yx, yy, yz, zx, zy, zz]
                let padded: Vec<f64> = raw_values.chunks(4)
                    .flat_map(|c| [c[0], c[2], 0.0, c[3], c[1], 0.0, 0.0, 0.0, 0.0])
                    .collect();
                Attribute::tensors(attr_name).with_data(IOBuffer::F64(padded))
            },
            _ => {
                // Handle higher-dimensional data
                if type_info.dimension == 6 {
                    // 3D stress/strain tensor in Voigt notation (xx, yy, zz, xy, yz, xz)
                    Attribute::generic(attr_name, 6).with_data(IOBuffer::F64(raw_values))
                } else {
                    // Generic field with arbitrary dimension
                    Attribute::field(attr_name).with_data(IOBuffer::F64(raw_values))
                }
            }
        };
        point_attributes.push(attribute);
    }
    
        
        // Validation: Ensure consistency of element-related data structures
        println!("Total elements: {}", mesh_data.elements.len());
        println!("Cell types count: {}", cell_types.len());
        println!("Offsets count: {}", offsets.len());
        assert_eq!(cell_types.len(), mesh_data.elements.len());
        assert_eq!(offsets.len(), mesh_data.elements.len());

        // Step 6: Construct the VTK data structure with point attributes
        Vtk {
            version: Version { major: 2, minor: 2 },   // VTK XML format version
            title: String::new(),                      // Optional title for the dataset
            byte_order: ByteOrder::LittleEndian,       // Byte order for binary data
            file_path: None,                           // File path set during write
            data: DataSet::inline(UnstructuredGridPiece {
                points: IOBuffer::F64(points_data),    // Node coordinates
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity: connectivity,     // Flat array of node IDs
                        offsets: offsets,               // Cumulative node counts
                    },
                    types: cell_types,                  // VTK cell types
                },
                data: Attributes {
                    point: point_attributes,            // Field values at nodes
                    cell: Default::default(),           // No cell-centered data
                },
            }),
        }.write_xml(&mut vtu)?;

        // Step 7: Write the serialized XML data to the output file
        fs::write(output_path, &vtu)?;

        Ok(())
    }

    /// Writes mesh geometry along with quality metrics and Gaussian point analysis to a VTU file.
    /// Quality metrics and integration point data are attached as cell data attributes.
    /// 
    /// # Arguments
    /// * `mesh_data` - The mesh containing nodes and elements
    /// * `element_quality` - Mesh quality metrics (volume, shape, skewness, etc.)
    /// * `gaussian_points_report` - Analysis of optimal integration points per element
    /// * `output_path` - File path where the VTU file will be written
    /// 
    /// # Returns
    /// * `Ok(())` on success
    /// * `Err(WriterError)` if file writing or VTK serialization fails
    pub fn write_vtu_with_geometric_analysis_and_gaussian_points_number(
        mesh_data: &MeshData,
        element_quality: &MeshQualityReport,
        gaussian_points_report: &GaussianPointNumberReport,
        output_path: &str,
    ) -> Result<(), WriterError> {

        // Buffer to hold the serialized VTU XML data
        let mut vtu = Vec::new();

        // Step 1: Convert node coordinates to VTK points format
        let points_data: Vec<f64> = mesh_data
            .nodes
            .iter()
            .flat_map(|node| {
                node.coordinates
                    .iter()
                    .copied()
                    .chain(std::iter::repeat(0.0)) // Pad with zeros for missing dimensions
                    .take(3)                       // Ensure exactly 3 components (x, y, z)
            })
            .collect();

        // Step 2: Pre-allocate vectors for cell topology data
        let total_connectivity: usize = mesh_data.elements.iter().map(|e| e.nodes.len()).sum();

        let mut connectivity = Vec::with_capacity(total_connectivity);
        let mut offsets = Vec::with_capacity(mesh_data.elements.len());
        let mut cell_types = vec![CellType::Vertex; mesh_data.elements.len()];
        let mut current_offset = 0;

        // Step 3: Build connectivity and offset arrays
        for element in &mesh_data.elements {
            connectivity.extend(element.nodes.iter().map(|&id| id as u64));
            current_offset += element.nodes.len() as u64;
            offsets.push(current_offset);
        }

        // Step 4: Assign VTK cell types based on element type groups
        for group in &mesh_data.element_type_info {
            let vtk_type = group.element_type.eltype_vtk();
            
            // Assign the VTK cell type to all elements in this group's range
            for i in group.start_index..(group.start_index + group.num_elements) {
                if i < cell_types.len() {
                    cell_types[i] = vtk_type;
                } else {
                    return Err(WriterError::IO(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Element index {} out of bounds (total elements: {})", i, cell_types.len())
                    )));
                }
            }
        }

        // Step 5: Prepare cell attributes from quality metrics and Gaussian point analysis
        let mut cell_attributes = Vec::new();

        // Volume Metric: Measures element volume/area relative to ideal shape
        let mut volume_metric_values = vec![0.0; mesh_data.elements.len()];
        for quality in &element_quality.element_qualities {
            if quality.element_id < volume_metric_values.len() {
                volume_metric_values[quality.element_id] = quality.volume_metric;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Volume_Metric", 1)
                .with_data(IOBuffer::F64(volume_metric_values))
        );

        // Shape Metric: Measures deviation from ideal element shape
        let mut shape_metric_values = vec![0.0; mesh_data.elements.len()];
        for quality in &element_quality.element_qualities {
            if quality.element_id < shape_metric_values.len() {
                shape_metric_values[quality.element_id] = quality.shape_metric;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Shape_Metric", 1)
                .with_data(IOBuffer::F64(shape_metric_values))
        );

        // Skewness Metric: Measures element skewness/distortion
        let mut skewness_metric_values = vec![0.0; mesh_data.elements.len()];
        for quality in &element_quality.element_qualities {
            if quality.element_id < skewness_metric_values.len() {
                skewness_metric_values[quality.element_id] = quality.skewness_metric;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Skewness_Metric", 1)
                .with_data(IOBuffer::F64(skewness_metric_values))
        );

        // Length Ratio: Ratio of longest to shortest edge
        let mut length_ratio_values = vec![0.0; mesh_data.elements.len()];
        for quality in &element_quality.element_qualities {
            if quality.element_id < length_ratio_values.len() {
                length_ratio_values[quality.element_id] = quality.length_ratio;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Length_Ratio", 1)
                .with_data(IOBuffer::F64(length_ratio_values))
        );

        // Orientation Metric: Measures element alignment with coordinate axes
        let mut orientation_metric_values = vec![0.0; mesh_data.elements.len()];
        for quality in &element_quality.element_qualities {
            if quality.element_id < orientation_metric_values.len() {
                orientation_metric_values[quality.element_id] = quality.orientation_metric;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Orientation_Metric", 1)
                .with_data(IOBuffer::F64(orientation_metric_values))
        );

        // Jacobian Ratio: Ratio of min to max Jacobian determinant at integration points
        let mut jacobian_ratio_values = vec![0.0; mesh_data.elements.len()];
        for quality in &element_quality.element_qualities {
            if quality.element_id < jacobian_ratio_values.len() {
                jacobian_ratio_values[quality.element_id] = quality.jacobian_ratio;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Jacobian_Ratio", 1)
                .with_data(IOBuffer::F64(jacobian_ratio_values))
        );

        // Optimal Gaussian Points: Recommended number of integration points based on quality
        let mut gauss_points_values = vec![0.0; mesh_data.elements.len()];
        for gauss_point in &gaussian_points_report.gauss_point_numbers {
            if gauss_point.element_id < gauss_points_values.len() {
                gauss_points_values[gauss_point.element_id] = gauss_point.optimal_number as f64;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Optimal_Gaussian_Points", 1)
                .with_data(IOBuffer::F64(gauss_points_values))
        );

        // Theoretical Gaussian Points: Standard number of integration points for element type
        let mut theoretical_gauss_values = vec![0.0; mesh_data.elements.len()];
        for gauss_point in &gaussian_points_report.gauss_point_numbers {
            if gauss_point.element_id < theoretical_gauss_values.len() {
                theoretical_gauss_values[gauss_point.element_id] = gauss_point.theoretical_number as f64;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Theoretical_Gaussian_Points", 1)
                .with_data(IOBuffer::F64(theoretical_gauss_values))
        );

        // Gaussian Point Difference: Shows potential over/under-integration
        // Positive values indicate potential over-integration, negative indicate under-integration
        let mut gauss_point_diff = vec![0.0; mesh_data.elements.len()];
        for gauss_point in &gaussian_points_report.gauss_point_numbers {
            if gauss_point.element_id < gauss_point_diff.len() {
                gauss_point_diff[gauss_point.element_id] = 
                    gauss_point.theoretical_number as f64 - gauss_point.optimal_number as f64;
            }
        }
        cell_attributes.push(
            Attribute::scalars("Gaussian_Point_Difference", 1)
                .with_data(IOBuffer::F64(gauss_point_diff))
        );

        // Validation: Ensure consistency of element-related data structures
        println!("Total elements: {}", mesh_data.elements.len());
        println!("Cell types count: {}", cell_types.len());
        println!("Offsets count: {}", offsets.len());
        assert_eq!(cell_types.len(), mesh_data.elements.len());
        assert_eq!(offsets.len(), mesh_data.elements.len());

        // Step 6: Construct the VTK data structure with cell attributes
        Vtk {
            version: Version { major: 2, minor: 2 },   // VTK XML format version
            title: String::new(),                      // Optional title for the dataset
            byte_order: ByteOrder::LittleEndian,       // Byte order for binary data
            file_path: None,                           // File path set during write
            data: DataSet::inline(UnstructuredGridPiece {
                points: IOBuffer::F64(points_data),    // Node coordinates
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity: connectivity,     // Flat array of node IDs
                        offsets: offsets,               // Cumulative node counts
                    },
                    types: cell_types,                  // VTK cell types
                },
                data: Attributes {
                    point: Default::default(),          // No point-centered data
                    cell: cell_attributes,              // Quality and analysis metrics per element
                },
            }),
        }.write_xml(&mut vtu)?;

        // Step 7: Write the serialized XML data to the output file
        fs::write(output_path, &vtu)?;

        Ok(())
    }
}

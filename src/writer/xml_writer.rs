use vtkio::model::*; // import model definition of a VTK file
use vtkio::Error as WriterError;

use std::fs;
use std::vec;
//use std::io::Write;

use crate::database::*;  // Import model definitions (Node, Element, ElementType, etc.)


pub struct VTUWriter;  // Defines a structure to represent a mesh parser for MPHTXT files

impl VTUWriter {

    pub fn write_vtu(
        mesh_data: &MeshData,
        output_path: &str,
    ) -> Result<(), WriterError> {  //ParseError gibi bir sey lazim derive debug seyi olabilsin diye

        let mut vtu = Vec::new();  // ?? why vector?? Create a new vector to hold the VTU XML data ??

        // 1. Prepare points data
        let points_data: Vec<f64> = mesh_data  //u32
            .nodes
            .iter()
            .flat_map(|node| node.coordinates.iter().copied())
            .collect();

        // 2. Pre-calculate sizes
        let total_connectivity: usize = mesh_data.elements.iter().map(|e| e.nodes.len()).sum();

        let mut connectivity = Vec::with_capacity(total_connectivity);
        let mut offsets = Vec::with_capacity(mesh_data.elements.len());
        let mut cell_types = vec![CellType::Vertex; mesh_data.elements.len()];
        let mut current_offset = 0;

        // 3. Process elements individually (for connectivity and offsets)
        for element in &mesh_data.elements {
            connectivity.extend(element.nodes.iter().map(|&id| id as u64));
            current_offset += element.nodes.len() as u64;
            offsets.push(current_offset);
        }

        // 4. Process groups ONLY for cell types
        for group in &mesh_data.element_type_info {
            let vtk_type = group.element_type.eltype_vtk();
            
            // Fill cell types for this group
            for i in group.start_index..(group.start_index + group.num_elements) {
                if i < cell_types.len() {
                    cell_types[i] = vtk_type;
                } else {
                    panic!("Element index out of bounds: {} >= {}", i, cell_types.len());
                }
            }
        }
        
        println!("Total elements: {}", mesh_data.elements.len());
        println!("Cell types count: {}", cell_types.len());
        println!("Offsets count: {}", offsets.len());
        assert_eq!(cell_types.len(), mesh_data.elements.len());
        assert_eq!(offsets.len(), mesh_data.elements.len());

        Vtk {
            version: Version { major: 2, minor: 2 },   //Check the for paraview neeeded version
            title: String::new(),  //String::from("VTU example"),
            byte_order: ByteOrder::LittleEndian,
            file_path: None,
            data: DataSet::inline(UnstructuredGridPiece {
                points: IOBuffer::F64(points_data),
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity: connectivity,
                        offsets: offsets,
                    },
                    types: cell_types,  
                },
                data: Attributes {
                    ..Default::default()
                },
            }),
        }.write_xml(&mut vtu)?;

        // Write the vector to file
        fs::write(output_path, &vtu)?;

        Ok(())
    }

    pub fn write_vtu_with_field_values(
        mesh_data: &MeshData,
        node_value_data: &NodeValueData,
        element_value_data: &ElementValueData,
        output_path: &str,
    ) -> Result<(), WriterError> {  //ParseError gibi bir sey lazim derive debug seyi olabilsin diye

        let mut vtu = Vec::new();  // ?? why vector?? Create a new vector to hold the VTU XML data ??

        // 1. Prepare points data
        let points_data: Vec<f64> = mesh_data  //u32
            .nodes
            .iter()
            .flat_map(|node| node.coordinates.iter().copied())
            .collect();

        // 2. Pre-calculate sizes
        let total_connectivity: usize = mesh_data.elements.iter().map(|e| e.nodes.len()).sum();

        let mut connectivity = Vec::with_capacity(total_connectivity);
        let mut offsets = Vec::with_capacity(mesh_data.elements.len());
        let mut cell_types = vec![CellType::Vertex; mesh_data.elements.len()];
        let mut current_offset = 0;

        // 3. Process elements individually (for connectivity and offsets)
        for element in &mesh_data.elements {
            connectivity.extend(element.nodes.iter().map(|&id| id as u64));
            current_offset += element.nodes.len() as u64;
            offsets.push(current_offset);
        }

        // 4. Process groups ONLY for cell types
        for group in &mesh_data.element_type_info {
            let vtk_type = group.element_type.eltype_vtk();
            
            // Fill cell types for this group
            for i in group.start_index..(group.start_index + group.num_elements) {
                if i < cell_types.len() {
                    cell_types[i] = vtk_type;
                } else {
                    panic!("Element index out of bounds: {} >= {}", i, cell_types.len());
                }
            }
        }

    // Node value data preparation following the same pattern as element data

    // Build VTK Attributes from NodeValueData
    let mut point_attributes = Vec::new();

    for type_info in &node_value_data.node_value_type_info {
        // Create a vector to hold all values for this type
        let mut values_vec = Vec::with_capacity(type_info.num_nodes * type_info.dimension);
        
        // Extract values for this specific type from NodeValue structs
        let start_idx = type_info.start_index;
        let end_idx = start_idx + type_info.num_nodes;
        
        if end_idx <= node_value_data.node_values.len() {
            // Get the slice of NodeValue structs for this type
            let node_values_slice = &node_value_data.node_values[start_idx..end_idx];
            
            // Extract the actual values from each NodeValue struct
            for node_value in node_values_slice {
                values_vec.extend_from_slice(&node_value.values);
            }
            
            // Get the attribute name from ValueType
            let attr_name = ValueType::get_attribute_name(&type_info.nodevalue_type);
            
            // Create appropriate attribute based on dimension
            let attribute = match type_info.dimension {
                1 => {
                    // Scalar data
                    Attribute::scalars(attr_name, 1)
                        .with_data(IOBuffer::F64(values_vec))
                },
                3 => {
                    // Vector data (3D)
                    Attribute::vectors(attr_name)
                        .with_data(IOBuffer::F64(values_vec))
                },
                6 => {
                    // Tensor data (6 components for symmetric tensor)
                    Attribute::generic(attr_name, 6)
                        .with_data(IOBuffer::F64(values_vec))
                },
                9 => {
                    // Full tensor data (9 components)
                    Attribute::tensors(attr_name)
                        .with_data(IOBuffer::F64(values_vec))
                },
                _ => {
                    // Generic field data for other dimensions
                    Attribute::field(attr_name)
                        .with_data(IOBuffer::F64(values_vec))
                }
            };
            
            point_attributes.push(attribute);
        } else {
            panic!("Node value index out of bounds: {} > {}", end_idx, node_value_data.node_values.len());
        }
    }

    // Element value data preparation
    let mut cell_attributes = Vec::new();

    for type_info in &element_value_data.element_value_type_info {
        // Create a vector to hold all values for this type
        let mut values_vec = Vec::with_capacity(type_info.num_elements * type_info.dimension);
        
        // Extract values for this specific type from ElementValue structs
        let start_idx = type_info.start_index;
        let end_idx = start_idx + type_info.num_elements;
        
        if end_idx <= element_value_data.element_values.len() {
            // Get the slice of ElementValue structs for this type
            let element_values_slice = &element_value_data.element_values[start_idx..end_idx];
            
            // Extract the actual values from each ElementValue struct
            for element_value in element_values_slice {
                values_vec.extend_from_slice(&element_value.values);
            }
            
            let attr_name = ValueType::get_attribute_name(&type_info.elementvalue_type);
            
            let attribute = match type_info.dimension {
                1 => Attribute::scalars(attr_name, 1).with_data(IOBuffer::F64(values_vec)),
                3 => Attribute::vectors(attr_name).with_data(IOBuffer::F64(values_vec)),
                6 => Attribute::generic(attr_name, 6).with_data(IOBuffer::F64(values_vec)),
                9 => Attribute::tensors(attr_name).with_data(IOBuffer::F64(values_vec)),
                _ => Attribute::field(attr_name).with_data(IOBuffer::F64(values_vec)),
            };
            
            cell_attributes.push(attribute);
        } else {
            panic!("Element value index out of bounds: {} > {}", end_idx, element_value_data.element_values.len());
        }
    }
    
        
        println!("Total elements: {}", mesh_data.elements.len());
        println!("Cell types count: {}", cell_types.len());
        println!("Offsets count: {}", offsets.len());
        assert_eq!(cell_types.len(), mesh_data.elements.len());
        assert_eq!(offsets.len(), mesh_data.elements.len());

        Vtk {
            version: Version { major: 2, minor: 2 },   //Check the for paraview neeeded version
            title: String::new(),  //String::from("VTU example"),
            byte_order: ByteOrder::LittleEndian,
            file_path: None,
            data: DataSet::inline(UnstructuredGridPiece {
                points: IOBuffer::F64(points_data),
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity: connectivity,
                        offsets: offsets, 
                    },
                    types: cell_types,  // Cell types for each element 
                },
                data: Attributes {
                    point: point_attributes,
                    cell: cell_attributes,
                },
            }),
        }.write_xml(&mut vtu)?;

        // Write the vector to file
        fs::write(output_path, &vtu)?;

        Ok(())
    }

}

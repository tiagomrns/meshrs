use std::f64; // bring f64 constants (e.g., std::f64::consts)

use crate::structs_and_impls::*; // import Node, Element, ElementType, ShapeFunction, Jacobian, MeshData, etc.
use crate::error::*;             // import your ElementError and related error types

/// Geometry routines (Jacobian, quality metrics, etc.)
pub struct GeometricAnalysis;

impl GeometricAnalysis {
    /// Analyze mesh quality over all elements (skips Vertex elements).
    pub fn analyse_mesh_quality(mesh_data: &MeshData) -> Result<MeshQualityReport, ElementError> {
        let mut element_qualities = Vec::new(); // accumulate per-element quality results
        let mut processed_elements = 0;          // count of successfully processed elements

        // iterate over each element-type block in the mesh
        for type_info in &mesh_data.element_type_info {
            // skip vertex elements (0D) since Jacobian is not meaningful
            if matches!(type_info.element_type, ElementType::Vertex) {
                continue; // move to next type
            }

            let start_idx = type_info.start_index;                 // first element index in this type block
            let end_idx = start_idx + type_info.num_elements;      // one past the last element index

            // iterate elements of this type by index into mesh_data.elements
            for element_idx in start_idx..end_idx {
                // guard against out-of-bounds (robustness)
                if element_idx < mesh_data.elements.len() {
                    let element = &mesh_data.elements[element_idx]; // borrow the element

                    // compute element quality (determinant of Jacobian at center)
                    match Self::calculate_element_quality(element, &type_info.element_type, &mesh_data.nodes) {
                        Ok(quality) => {
                            element_qualities.push(quality); // store result
                            processed_elements += 1;         // increment success counter
                        }
                        Err(e) => {
                            // non-fatal: report but continue with other elements
                            println!("Warning: Failed to analyze element {}: {:?}", element.id, e);
                        }
                    }
                }
            }
        }

        // if nothing was processed, return an error
        if element_qualities.is_empty() {
            return Err(ElementError::GeometryError("No elements could be analyzed".to_string()));
        }

        // build and return the report
        Ok(MeshQualityReport {
            total_elements: processed_elements, // how many we processed
            element_qualities,                  // the list of per-element qualities
        })
    }

    /// Calculate per-element quality (currently determinant of J at element center).
    fn calculate_element_quality(
        element: &Element,            // which element
        element_type: &ElementType,   // its type (Line, Triangle, etc.)
        nodes: &[Node],               // all nodes in the mesh
    ) -> Result<ElementQuality, ElementError> {
        let element_nodes = Self::get_element_nodes(element, nodes)?; // gather Node structs for this element
        let center_point = Self::get_center_point(element_type);      // reference (natural) center point

        // evaluate shape functions (and derivatives) at center_point, using f64 arithmetic
        let shape_function = ElementType::get_shape_functions::<f64>(element_type, &center_point)
            .ok_or_else(|| ElementError::GeometryError("Failed to get shape functions".to_string()))?;

        // compute the Jacobian (and its determinant) using node coords and shape derivatives
        let jacobian = Self::calculate_jacobian(&element_nodes, &shape_function.derivatives)?;

        // return the element quality (currently just det(J))
        Ok(ElementQuality {
            element_id: element.id,            // pass through element id
            det_jacobian: jacobian.determinant // determinant at the center
        })
    }

    /// Return the natural/reference-space "center" point for each element type.
    fn get_center_point(element_type: &ElementType) -> Vec<f64> {
        match element_type {
            // 1D elements: midpoint in ξ
            ElementType::Line | ElementType::QuadraticEdge => vec![0.5],
            // 2D triangles: centroid (ξ=η=1/3)
            ElementType::Triangle | ElementType::QuadraticTriangle => vec![1.0 / 3.0, 1.0 / 3.0],
            // 2D quads: center (ξ=η=0.5)
            ElementType::Quad
            | ElementType::QuadraticQuad
            | ElementType::BiquadraticQuad => vec![0.5, 0.5],
            // 3D tetra: centroid (ξ=η=ζ=1/4)
            ElementType::Tetra | ElementType::QuadraticTetra => vec![0.25, 0.25, 0.25],
            // 3D pyramid: a typical interior point (ξ=η=ζ=0.5)
            ElementType::Pyramid => vec![0.5, 0.5, 0.5],
            // 3D wedge/prism: centroid in triangle (1/3,1/3) and mid through thickness (0.5)
            ElementType::Wedge
            | ElementType::QuadraticWedge
            | ElementType::BiquadraticQuadraticWedge => vec![1.0 / 3.0, 1.0 / 3.0, 0.5],
            // 3D hex family: center (0.5,0.5,0.5)
            ElementType::Hexahedron
            | ElementType::QuadraticHexahedron
            | ElementType::BiquadraticQuadraticHexahedron
            | ElementType::TriquadraticHexahedron => vec![0.5, 0.5, 0.5],
            // Vertex (0D) and fall-through default
            ElementType::Vertex => vec![0.0],
            _ => vec![0.0],
        }
    }

    /// Convert element’s node IDs to concrete Node structs (with basic validation).
    pub fn get_element_nodes(element: &Element, nodes: &[Node]) -> Result<Vec<Node>, ElementError> {
        let mut element_nodes = Vec::with_capacity(element.nodes.len()); // pre-allocate for speed

        // walk each node id in the element connectivity
        for &node_id in &element.nodes {
            // find the node with that id in the global nodes array
            if let Some(node) = nodes.iter().find(|n| n.id == node_id) {
                element_nodes.push(node.clone()); // push a clone (Node is small)
            } else {
                // if any id is missing, return a descriptive error
                return Err(ElementError::InvalidElement(format!(
                    "Node {} not found for element {}",
                    node_id, element.id
                )));
            }
        }

        // safety check: an element must have at least one node
        if element_nodes.is_empty() {
            return Err(ElementError::InvalidElement(format!(
                "No valid nodes found for element {}",
                element.id
            )));
        }

        Ok(element_nodes) // success
    }

    /// Compute Jacobian J and its determinant for an element mapping.
    ///
    /// J[i][j] = Σ_k ( dN_k/dξ_j * x_k_i )
    ///  - i indexes physical coordinates (x,y,z) → rows = mesh_dim
    ///  - j indexes natural coordinates (ξ,η,ψ)  → cols = element_dim
    ///  - k iterates element nodes
    pub fn calculate_jacobian<T: FloatLike>(
        element_nodes: &[Node],       // node coordinates in physical space
        shape_derivatives: &[Vec<T>], // dN_k/dξ_j for each node k and natural direction j
    ) -> Result<Jacobian<T>, ElementError> {
        // must have nodes
        if element_nodes.is_empty() {
            return Err(ElementError::GeometryError("No nodes provided".to_string()));
        }
        // node count must match number of derivative rows (one row per node)
        if element_nodes.len() != shape_derivatives.len() {
            return Err(ElementError::GeometryError(format!(
                "Node count ({}) doesn't match shape derivative count ({})",
                element_nodes.len(),
                shape_derivatives.len()
            )));
        }

        let mesh_dim = element_nodes[0].coordinates.len(); // x/y/z count
        let element_dim = shape_derivatives[0].len();       // ξ/η/ψ count

        // allocate J as mesh_dim × element_dim (initialized to zero)
        let mut jacobian_matrix = vec![vec![T::zero(); element_dim]; mesh_dim];

        // triple loop implementing J[i][j] = Σ_k dN_k/dξ_j * x_k_i
        for i in 0..mesh_dim {
            for j in 0..element_dim {
                for (k, node) in element_nodes.iter().enumerate() {
                    // guard against any inconsistent lengths
                    if j < shape_derivatives[k].len() && i < node.coordinates.len() {
                        let coord = T::from_f64(node.coordinates[i]).unwrap(); // cast f64 → T
                        jacobian_matrix[i][j] =
                            jacobian_matrix[i][j].clone() + shape_derivatives[k][j].clone() * coord;
                    }
                }
            }
        }

        // compute determinant (or geometric measure) of J
        let determinant = Self::calculate_jacobian_determinant(&jacobian_matrix)?;

        // package both matrix and determinant
        Ok(Jacobian {
            matrix: jacobian_matrix,
            determinant,
        })
    }

    /// Small-matrix determinant / geometric measure for Jacobians.
    ///
    /// Supported shapes:
    /// - 1×1, 2×2, 3×3 → standard determinant.
    /// - 2×1 → curve length in 2D: sqrt(dx^2 + dy^2).
    /// - 3×1 → curve length in 3D: sqrt(dx^2 + dy^2 + dz^2).
    /// - 3×2 → surface area factor in 3D: sqrt(det(JᵀJ)).
    fn calculate_jacobian_determinant<T: FloatLike>(matrix: &[Vec<T>]) -> Result<T, ElementError> {
        let rows = matrix.len(); // number of physical coords
        if rows == 0 {
            return Err(ElementError::GeometryError("Empty matrix".to_string())); // sanity check
        }

        let cols = matrix[0].len(); // number of natural coords
        // check rectangular shape (all rows same length)
        for row in matrix {
            if row.len() != cols {
                return Err(ElementError::GeometryError(
                    "Inconsistent matrix dimensions".to_string(),
                ));
            }
        }

        // branch by (rows, cols) to compute measure
        match (rows, cols) {
            // 1×1 determinant
            (1, 1) => Ok(matrix[0][0].clone()),

            // 2×2 determinant
            (2, 2) => Ok(matrix[0][0].clone() * matrix[1][1].clone()
                - matrix[0][1].clone() * matrix[1][0].clone()),

            // 3×3 determinant (Laplace expansion)
            (3, 3) => Ok(
                matrix[0][0].clone()
                    * (matrix[1][1].clone() * matrix[2][2].clone()
                        - matrix[1][2].clone() * matrix[2][1].clone())
                    - matrix[0][1].clone()
                        * (matrix[1][0].clone() * matrix[2][2].clone()
                            - matrix[1][2].clone() * matrix[2][0].clone())
                    + matrix[0][2].clone()
                        * (matrix[1][0].clone() * matrix[2][1].clone()
                            - matrix[1][1].clone() * matrix[2][0].clone()),
            ),

            // 2×1: length factor for a 1D element in 2D
            (2, 1) => {
                let g = matrix[0][0].clone() * matrix[0][0].clone()
                    + matrix[1][0].clone() * matrix[1][0].clone(); // dx^2 + dy^2
                Ok(g.sqrt()) // length = sqrt(dx^2 + dy^2)
            }

            // 3×1: length factor for a 1D element in 3D
            (3, 1) => {
                let g = matrix[0][0].clone() * matrix[0][0].clone()
                    + matrix[1][0].clone() * matrix[1][0].clone()
                    + matrix[2][0].clone() * matrix[2][0].clone(); // dx^2 + dy^2 + dz^2
                Ok(g.sqrt()) // length
            }

            // 3×2: area factor for a 2D element in 3D → sqrt(det(JᵀJ)) where G = JᵀJ is 2×2
            (3, 2) => {
                // compute G = JᵀJ (2×2)
                let mut g = vec![vec![T::zero(); 2]; 2]; // allocate 2×2 zero
                for i in 0..2 {                // i: column of J, row of G
                    for j in 0..2 {            // j: column of J, col of G
                        for k in 0..3 {        // k: row of J
                            g[i][j] = g[i][j].clone() + matrix[k][i].clone() * matrix[k][j].clone();
                        }
                    }
                }
                // determinant of 2×2 G
                let det_g = g[0][0].clone() * g[1][1].clone() - g[0][1].clone() * g[1][0].clone();
                Ok(det_g.sqrt()) // area factor = sqrt(det(G))
            }

            // else: not implemented for other shapes (you can extend here)
            _ => Err(ElementError::GeometryError(format!(
                "Jacobian determinant calculation not implemented for {}x{} matrices",
                rows, cols
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*; // bring GeometricAnalysis & helpers into scope

    // helper: 3D node
    fn create_node(id: usize, x: f64, y: f64, z: f64) -> Node {
        Node { id, coordinates: vec![x, y, z] }
    }
    // helper: 2D node
    fn create_node_2d(id: usize, x: f64, y: f64) -> Node {
        Node { id, coordinates: vec![x, y] }
    }
    // helper: 1D node
    fn create_node_1d(id: usize, x: f64) -> Node {
        Node { id, coordinates: vec![x] }
    }

    // assemble MeshData for tests (minimal fields)
    fn create_test_mesh_data(
        nodes: Vec<Node>,
        elements: Vec<Element>,
        element_type_info: Vec<ElementTypeInfo>,
    ) -> MeshData {
        MeshData {
            dimension: if !nodes.is_empty() { nodes[0].coordinates.len() } else { 3 }, // deduce dim
            num_nodes: nodes.len(),                 // count nodes
            min_node_index: nodes.iter().map(|n| n.id).min().unwrap_or(1), // min id
            nodes,                                  // all nodes
            num_eltypes: element_type_info.len(),   // how many element type blocks
            elements,                               // all elements
            element_type_info,                      // per-type metadata
        }
    }

    // ---- tests below are unchanged, but annotated by names ----

    #[test]
    pub fn test_1d_line_element_1d_space() {
        let nodes = vec![create_node_1d(1, 0.0), create_node_1d(2, 2.0)];
        let element = Element { id: 1, nodes: vec![1, 2] };
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Line,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 2,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        assert_eq!(report.total_elements, 1);
        assert_eq!(report.element_qualities[0].element_id, 1);
        assert!((report.element_qualities[0].det_jacobian - 2.0).abs() < 1e-10); // length = 2
    }

    #[test]
    pub fn test_1d_line_element_2d_space() {
        let nodes = vec![create_node_2d(1, 0.0, 0.0), create_node_2d(2, 1.0, 1.0)];
        let element = Element { id: 1, nodes: vec![1, 2] };
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Line,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 2,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        let expected = (2.0f64).sqrt(); // length of diagonal from (0,0) to (1,1)
        assert!((report.element_qualities[0].det_jacobian - expected).abs() < 1e-10);
    }

    #[test]
    pub fn test_2d_triangle_element_2d_space() {
        // right triangle with area 0.5 → standard mapping gives |J| = 1 at center for these nodes
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
        ];
        let element = Element { id: 1, nodes: vec![1, 2, 3] };
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 3,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
    }

    #[test]
    pub fn test_2d_quad_element_2d_space() {
        // unit square → |J| = 1 for bilinear mapping at center
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 1.0, 1.0),
            create_node_2d(4, 0.0, 1.0),
        ];
        let element = Element { id: 1, nodes: vec![1, 2, 3, 4] };
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Quad,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 4,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
    }

    #[test]
    pub fn test_3d_tetrahedron_element_3d_space() {
        // reference tetra → |J| = 1 at center for these nodes
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
            create_node(3, 0.0, 1.0, 0.0),
            create_node(4, 0.0, 0.0, 1.0),
        ];
        let element = Element { id: 1, nodes: vec![1, 2, 3, 4] };
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Tetra,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 4,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
    }

    #[test]
    pub fn test_3d_hexahedron_element_3d_space() {
        // unit cube → |J| = 1 at center for trilinear mapping
        let nodes = vec![
            create_node(1, 0.0, 0.0, 0.0),
            create_node(2, 1.0, 0.0, 0.0),
            create_node(3, 1.0, 1.0, 0.0),
            create_node(4, 0.0, 1.0, 0.0),
            create_node(5, 0.0, 0.0, 1.0),
            create_node(6, 1.0, 0.0, 1.0),
            create_node(7, 1.0, 1.0, 1.0),
            create_node(8, 0.0, 1.0, 1.0),
        ];
        let element = Element { id: 1, nodes: vec![1, 2, 3, 4, 5, 6, 7, 8] };
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Hexahedron,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 8,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element.clone()], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        assert!((report.element_qualities[0].det_jacobian - 1.0).abs() < 1e-10);
    }

    #[test]
    pub fn test_complete_mesh_analysis() {
        // triangle + 2x2 quad (area scaling 4)
        let nodes = vec![
            create_node_2d(1, 0.0, 0.0),
            create_node_2d(2, 1.0, 0.0),
            create_node_2d(3, 0.0, 1.0),
            create_node_2d(4, 0.0, 0.0),
            create_node_2d(5, 2.0, 0.0),
            create_node_2d(6, 2.0, 2.0),
            create_node_2d(7, 0.0, 2.0),
        ];
        let elements = vec![
            Element { id: 1, nodes: vec![1, 2, 3] },      // triangle
            Element { id: 2, nodes: vec![4, 5, 6, 7] },   // quad
        ];
        let element_type_info = vec![
            ElementTypeInfo { element_type: ElementType::Triangle, start_index: 0, num_elements: 1, nodes_per_element: 3 },
            ElementTypeInfo { element_type: ElementType::Quad,     start_index: 1, num_elements: 1, nodes_per_element: 4 },
        ];
        let mesh_data = MeshData {
            dimension: 2,
            num_nodes: nodes.len(),
            min_node_index: 1,
            nodes,
            num_eltypes: element_type_info.len(),
            elements,
            element_type_info,
        };

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        assert!(result.is_ok());
        let report = result.unwrap();

        assert_eq!(report.total_elements, 2);

        let triangle_quality = report.element_qualities.iter().find(|q| q.element_id == 1).unwrap();
        assert!((triangle_quality.det_jacobian - 1.0).abs() < 1e-10);

        let quad_quality = report.element_qualities.iter().find(|q| q.element_id == 2).unwrap();
        assert!((quad_quality.det_jacobian - 4.0).abs() < 1e-10); // 2×2 scaling → area factor 4
    }

    #[test]
    pub fn test_invalid_node_reference() {
        // node id 3 is missing, should not panic
        let nodes = vec![create_node_2d(1, 0.0, 0.0), create_node_2d(2, 1.0, 0.0)];
        let element = Element { id: 1, nodes: vec![1, 2, 3] }; // invalid reference
        let element_type_info = vec![ElementTypeInfo {
            element_type: ElementType::Triangle,
            num_elements: 1,
            start_index: 0,
            nodes_per_element: 3,
        }];
        let mesh_data = create_test_mesh_data(nodes, vec![element], element_type_info);

        let result = GeometricAnalysis::analyse_mesh_quality(&mesh_data);
        // we accept either Ok (skip) or Err (fail this element) — but must not panic
        assert!(result.is_err() || result.is_ok());
    }
}

#[cfg(test)]
mod test_runner {
    use super::tests::*; // import all tests

    #[test]
    fn run_all_geometric_analysis_tests() {
        println!("Running Geometric Analysis Tests...");
        println!("=====================================");

        // table of test names and functions to execute
        let tests = vec![
            ("1D Line in 1D space", test_1d_line_element_1d_space as fn()),
            ("1D Line in 2D space", test_1d_line_element_2d_space as fn()),
            ("2D Triangle in 2D space", test_2d_triangle_element_2d_space as fn()),
            ("2D Quad in 2D space", test_2d_quad_element_2d_space as fn()),
            ("3D Tetrahedron in 3D space", test_3d_tetrahedron_element_3d_space as fn()),
            ("3D Hexahedron in 3D space", test_3d_hexahedron_element_3d_space as fn()),
            ("Complete mesh analysis", test_complete_mesh_analysis as fn()),
            ("Invalid node reference", test_invalid_node_reference as fn()),
        ];

        let mut passed = 0usize; // counter for successes
        let mut failed = 0usize; // counter for failures

        for (name, f) in tests {
            print!("{:.<40}", name);                  // pretty print fixed width
            let result = std::panic::catch_unwind(|| f()); // catch panics to keep runner going
            if result.is_ok() {
                println!("PASSED");
                passed += 1;
            } else {
                println!("FAILED");
                failed += 1;
            }
        }

        println!("=====================================");
        println!("Total: {} tests, {} passed, {} failed", passed + failed, passed, failed);
        assert_eq!(failed, 0, "Some tests failed!"); // fail CI if any test failed
    }
}

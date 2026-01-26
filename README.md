# meshrs

![Rust](https://img.shields.io/badge/Language-Rust-blue) ![License](https://img.shields.io/github/license/tiagomrns/meshrs)

---

## Executive Summary

The meshrs addresses the interoperability gap between proprietary simulation software, academic codes, and in-house solvers by developing a modular framework in the Rust programming language. The framework enables the conversion of proprietary mesh formats (such as COMSOL `.mphtxt` and Abaqus `.inp`) to the open-standard VTU format for visualization in ParaView, while providing extensible modules for user-defined pre-processing analysis. Key contributions include mesh quality assessment using metrics like Jacobian ratio and skewness, integration point optimization based on material properties and a user-defined error tolerance, and deformation visualization of raw result data from `.txt` files. By leveraging Rust's memory safety, performance characteristics, and growing industry adoption, this work establishes a foundation for reliable, high-performance engineering software tools that bridge legacy systems and modern development practices.

---

## Key Features

- **Parsing**
  - Supports `.inp` (Abaqus) and `.mphtxt` (COMSOL) file formats for mesh parsing and `.txt` file format for nodal value parsing.
  
- **Analysis**
  - Calculates geometric quality metrics, including:
    - Knupp's algebraic quality metrics (Shape Metric, Skewness Metric, Length Ratio, Orientation Metric, Volume Metric)
    - Determinant check
    - Jacobian Ratio

- **Optimization**
  - Performs mass and stiffness matrix integration optimization of Gaussian quadrature points.
  
- **Output**
  - Generates VTU files compatible with ParaView for intuitive visualization.
  - Generates TXT files for detailed analysis.

---

## Installation

To get started with the meshrs, clone this repository and use **Cargo** (Rust's package manager) to build the project. Run the following commands:

```bash
git clone https://github.com/tiagomrns/meshrs.git
cd meshrs
cargo build --release
```

---

## Usage Examples

Below are various use cases for the tool:

1. **Basic Mesh Parsing**  
   Generate a VTU representation of a COMSOL mesh file:
   ```bash
   cargo run --release -- mesh.mphtxt
   ```

2. **Geometric Quality Analysis**
   Perform mesh quality analysis for an Abaqus mesh file:
   ```bash
   cargo run --release -- model.inp --analyze-mesh
   ```

3. **Gaussian Quadrature Optimization (Mass Matrix Integration)**  
   Optimize Gaussian integration points for the mass matrix:
   ```bash
   cargo run --release -- input_mesh.inp --optimize-gauss --gauss-tolerance 1e-4 --mass-integration --material-property 10.0,5.0,3.0,0.5
   ```

4. **Gaussian Quadrature Optimization (Stiffness Matrix Integration)**  
   Optimize Gaussian integration points for the stiffness matrix using provided material property for 3D mesh (6 at the beginning shows that the material matrix is 6x6):
   ```bash
   cargo run --release -- input_mesh.inp --optimize-gauss --stiffness-integration --material-property '6;10.0;1.3;4.0;0.6;1.0;0.3;1.3;0.8;0.3;0.1;0.5;0.7;4.0;0.3;0.7;0.2;0.8;0.4;0.6;0.1;0.2;0.9;0.0;0.1;1.0;0.5;0.8;0.0;0.6;0.1;0.3;0.7;0.4;0.1;0.1;0.3'
   ```

---

## Output Files

The tool generates the following output files for further analysis:

| File Extension                        | Description                                                    |
|---------------------------------------|----------------------------------------------------------------|
| `<mesh_name>_parsed_mesh_data.txt`    | Parsed mesh data, including nodes and elements, in a human-readable format. |
| `<mesh_name>_parsed_data.txt          | Raw parsed nodal values (if applicable).
| `<mesh_name>_mesh_quality_analysis.txt` | Detailed quality metrics of the parsed mesh, including shape, skewness, and Jacobian ratios. |
| `<mesh_name>_mesh_gauss_optimization.txt` | Report on the optimal number of Gaussian points found for each element. |
| `<mesh_name>_with_analysis.vtu` | VTU file containing the mesh geometry and computed quality metrics. |
| `<mesh_name>_quadrature_optimized_stiffness.vtu` | VTU file with optimized Gaussian integration points for Stiffness Matrix. |

---

## Developer Guide

### Extending the meshrs

### Adding a New ElementType

1. **Define a new `ElementType`**:
   - Open [structs_and_impls.rs](https://github.com/tiagomrns/meshrs/blob/main/src/structs_and_impls.rs).
   - Add a new variant to the `ElementType` enum.
   
2. **Shape Functions**:  
   - Implement the `get_shape_functions` method for your new `ElementType`.

3. **Integration Rules**:  
   - Extend the implementation in [gaussian_quadrature.rs](https://github.com/tiagomrns/meshrs/blob/main/src/mesh_analysis/gaussian_quadrature.rs) to define symbolic or numerical integration rules.

---

### Adding Support for a New Mesh Parser

1. **Create a New Parser Module**:  
   - Add your parser under the `src/parser` directory.  
   - Implement the necessary parsing logic (see [abaqus_inp.rs](https://github.com/tiagomrns/meshrs/blob/main/src/parser/abaqus/abaqus_inp.rs)).

2. **Connect Parser**:  
   - Register your parser in [parser/mod.rs](https://github.com/tiagomrns/meshrs/blob/main/src/parser/mod.rs).

---

## Known Limitations & Future Work

- **Hardcoded Output Path**:
  - The output directory is currently hardcoded in [main.rs](https://github.com/tiagomrns/meshrs/blob/main/src/main.rs).  
  - **Action for Users**: You must manually modify the `output_dir` variable in `main.rs` to match your local file system before running the tool. Failure to do so may result in the program panicking due to an invalid path.
 
- **High Order Shape Functions and Derivatives**:
  - Shape functions for Quadratic Hexahedron, Biquadratic Hexahedron, and Triquadratic Hexahedron require verification before use. 
  - **Action for Users**: Please review the implementation in structs_and_impls.rs if utilizing these specific element types for high-precision simulations.

- **Symbolic Integration Approximations**:
  - The symbolic math engine used in the tool has limitations when handling non-polynomial functions such as rational functions and square roots.
  
  - **Issue**: The following critical calculations rely on a Taylor Series expansion:
    - Calculation of the Inverse Jacobian (required for Stiffness Matrix in both square and non-square cases).
    - Computation of the Metric Determinant (required for Mass Matrices with non-square Jacobians).
    
  - **Consequence**:  
    - Using Taylor Series introduces approximation errors, which can accumulate during stiffness matrix calculations, leading to less precise results.
    
  - **Future Solution**:  
    - Replace the Taylor Series expansion with Chebyshev Interpolation or alternative polynomial interpolation methods to better approximate non-polynomial terms. This enhancement would significantly reduce numerical errors in the calculation process.

---

This README provides an overview of the **meshrs**, summarizing its capabilities, installation, usage, known limitations, and guidelines for extending its functionality.

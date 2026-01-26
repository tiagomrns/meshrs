# FEA Mesh Analysis & Optimization Tool

![Rust](https://img.shields.io/badge/Language-Rust-blue) ![License](https://img.shields.io/github/license/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik)

---

## Executive Summary

The FEA Mesh Analysis & Optimization Tool is a cutting-edge Rust-based solution designed to streamline the parsing, analysis, and optimization processes of finite element analysis meshes used in engineering simulations. This tool supports Abaqus (.inp) and COMSOL (.mphtxt) mesh formats for parsing, evaluates geometric quality using metrics like Knupp's quality metrics (Shape, Skewness, and Jacobian Ratio), and optimizes Gaussian quadrature points using symbolic mathematics. With output in the VTU format, the tool ensures seamless integration with ParaView for visualization and further processing.

---

## Key Features

- **Parsing**
  - Supports `.inp` (Abaqus) and `.mphtxt` (COMSOL) file formats.
  
- **Analysis**
  - Calculates geometric quality metrics, including:
    - Shape metric
    - Skewness
    - Jacobian Ratio

- **Optimization**
  - Performs mass and stiffness matrix integration optimization of Gaussian quadrature points.
  
- **Output**
  - Generates VTU files compatible with ParaView for intuitive visualization.

---

## Installation

To get started with the FEA Tool, clone this repository and use **Cargo** (Rust's package manager) to build the project. Run the following commands:

```bash
git clone https://github.com/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik.git
cd Master-Thesis-Caglar-Ayberk-Saglik
cargo build --release
```

---

## Usage Examples

Below are various use cases for the tool:

1. **Basic Mesh Parsing**  
   Generate a VTU representation of a COMSOL mesh file:
   ```bash
   cargo run --release -- mesh.mphtxt --print-vtu
   ```

2. **Geometric Quality Analysis**
   Perform mesh quality analysis for an Abaqus mesh file:
   ```bash
   cargo run --release -- model.inp --analyze-mesh
   ```

3. **Gaussian Quadrature Optimization (Mass Matrix Integration)**  
   Optimize Gaussian integration points for the mass matrix:
   ```bash
   cargo run --release -- input_mesh.inp --optimize-gauss --mass-integration
   ```

4. **Gaussian Quadrature Optimization (Stiffness Matrix Integration)**  
   Optimize Gaussian integration points for the stiffness matrix using provided material property:
   ```bash
   cargo run --release -- input_mesh.inp --optimize-gauss --stiffness-integration --material-property "[1.0, 0.0, 0.5]"
   ```

---

## Output Files

The tool generates the following output files for further analysis:

| File Extension                        | Description                                                    |
|---------------------------------------|----------------------------------------------------------------|
| `<mesh_name>_parsed_mesh_data.txt`    | Parsed mesh data, including nodes and elements, in a human-readable format. |
| `<mesh_name>_mesh_quality_analysis.txt` | Detailed quality metrics of the parsed mesh, including shape, skewness, and Jacobian ratios. |
| `<mesh_name>_quadrature_optimized_mass.vtu` | VTU file with optimized Gaussian integration points for Mass Matrix. |
| `<mesh_name>_quadrature_optimized_stiffness.vtu` | VTU file with optimized Gaussian integration points for Stiffness Matrix. |

---

## Developer Guide

### Extending the FEA Tool

### Adding a New ElementType

1. **Define a new `ElementType`**:
   - Open [structs_and_impls.rs](https://github.com/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik/blob/main/src/structs_and_impls.rs).
   - Add a new variant to the `ElementType` enum.
   
2. **Shape Functions**:  
   - Implement the `get_shape_functions` method for your new `ElementType`.

3. **Integration Rules**:  
   - Extend the implementation in [gaussian_quadrature.rs](https://github.com/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik/blob/main/src/mesh_analysis/gaussian_quadrature.rs) to define symbolic or numerical integration rules.

---

### Adding Support for a New Mesh Parser

1. **Create a New Parser Module**:  
   - Add your parser under the `src/parser` directory.  
   - Implement the necessary parsing logic (see [abaqus_inp.rs](https://github.com/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik/blob/main/src/parser/abaqus/abaqus_inp.rs)).

2. **Connect Parser**:  
   - Register your parser in [parser/mod.rs](https://github.com/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik/blob/main/src/parser/mod.rs).

---

## Known Limitations & Future Work

- **Hardcoded Output Path**:
  - The output directory is currently hardcoded in [main.rs](https://github.com/ayberksaglik/Master-Thesis-Caglar-Ayberk-Saglik/blob/main/src/main.rs).  
  - **Action for Users**: You must manually modify the `output_dir` variable in `main.rs` to match your local file system before running the tool. Failure to do so may result in the program panicking due to an invalid path.

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

This README provides an overview of the **FEA Mesh Analysis & Optimization Tool**, summarizing its capabilities, installation, usage, known limitations, and guidelines for extending its functionality.
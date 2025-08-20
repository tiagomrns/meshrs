# meshrs

[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg?maxAge=3600)](https://github.com/rust-lang/rust)

**A high-performance, safe, and modular toolkit for Finite Element mesh I/O, manipulation, and validation, written in Rust.**

`mesh-rs` aims to provide a robust foundation for scientific computing and CAE applications, leveraging Rust's strengths to avoid common pitfalls in mesh processing like memory unsafety, data races, and undefined behavior.

## Features

*   **I/O Support:** Read and write mesh files in various formats.
    *   [x] Gmsh (.msh, format 4.1)
    *   [ ] VTK Legacy (.vtk)
    *   [ ] Abaqus (.inp)
    *   [ ] XDMF (via `hdf5` crate)
    *   *... more to come! (Contributions welcome)*

*   **Core Data Structures:** Efficient, generic representations of mesh topology and geometry.
    *   `Mesh`: The primary structure holding nodes, elements, and physical groups.
    *   `ElementType`: Enumeration of standard 1D, 2D, and 3D element types (Beam, Tri, Quad, Tet, Hex, etc.).
    *   `PhysicalGroup`: Tag groups of entities (nodes, elements) with physical identifiers and names.

*   **Manipulation & Queries:**
    *   Iterate over nodes, elements, and boundaries.
    *   Calculate element quality metrics (aspect ratio, skewness, Jacobian).
    *   Find connected elements (node-to-element connectivity).
    *   *Planned: Mesh refinement, coarsening, and transformations.*

*   **Validation & Quality:** Tools to ensure mesh integrity before simulation.
    *   Check for duplicate nodes and elements.
    *   Validate element connectivity (e.g., no negative volumes).
    *   Detect and fix orphaned nodes.

*   **Performance:** Built with performance in mind from the ground up.
    *   Zero-cost abstractions.
    *   Parallel processing capabilities using `rayon`.

## Installation

Add `mesh-rs` to your `Cargo.toml`:

```toml
[dependencies]
mesh-rs = "0.1"

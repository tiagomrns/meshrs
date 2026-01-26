// Import standard library modules for system interaction
use std::env;                    // Provides access to command line arguments via env::args()
use std::process;               // Enables process control, particularly process::exit() for terminating with error codes
use std::path::{Path, PathBuf}; // Path manipulation utilities - Path for immutable path references, PathBuf for owned paths
use std::fs::File;              // For file operations
use std::io::Write;             // For writing to files

// Import custom modules that contain the core functionality
mod structs_and_impls;                   // Module containing data structures (MeshData, Node, Element, etc.)
use structs_and_impls::*;                // Wildcard import brings all public items from lib module into scope

// Import custom modules that contain the core functionality
mod error;                   // Module containing data structures (MeshData, Node, Element, etc.)
use error::*;                // Wildcard import brings all public items from error module into scope

mod parser;                     // Module containing file format parsers
use parser::comsol::mphtxt::MphtxtParser;       // Parser specifically for COMSOL .mphtxt mesh files
use parser::abaqus::abaqus_inp::AbaqusInpParser; // Parser for Abaqus .inp input files
use parser::comsol::txt_node::ComsolTxtNodeParser;       // Parser for nodal field values for COMSOL .txt result files

mod writer;                     // Module containing output file writers
use writer::xml_writer::VTUWriter; // Writer for VTK XML unstructured grid (.vtu) format

mod mesh_analysis;                   // Module containing mesh quality analysis  
use mesh_analysis::geometric_analysis::GeometricAnalysis; // Geometric analysis for mesh quality 
use mesh_analysis::gaussian_quadrature::*; // Gauss quadrature utilities for numerical integration 

/// Struct to hold parsed command line arguments in a structured format
/// This makes it easier to pass around and validate the user's input
struct CommandLineArgs {
    mesh_file: Option<String>,        // Path to mesh file (.mphtxt or .inp) - required
    node_txt_file: Option<String>,    // Path to txt file containing node field values - optional
    mesh_entity: MeshEntity,          // Specifies which type of values to parse from txt files
    print_vtu: bool,                  // Flag to control whether to print VTU file content to console
    analyze_mesh: bool,               // Flag to enable mesh quality analysis  
    optimize_gauss: bool,             // Flag to enable Gaussian quadrature optimization
    gauss_tolerance: f64,             // Tolerance for Gaussian quadrature optimization
    integration_type: IntegrationType, // Type of integration (Mass or Stiffness)
    material_property: MaterialProperty,      // Material property polynomial coefficients in monomial form
}

/// Enum to specify which type of values to parse from txt file(s)
/// This determines the parsing strategy and validation requirements
#[derive(Debug, Clone, PartialEq)] // Derive traits for debugging, cloning, and comparison
enum MeshEntity {
    Node,     // Parse txt file as containing node-based field values
    //Element,  // Parse txt file as containing element-based field values
    //Both,     // Parse two txt files (first=nodes, second=elements)
}

/// Parse command line arguments and validate them according to program requirements
/// Returns a structured CommandLineArgs object or an error message string
fn parse_arguments(args: Vec<String>) -> Result<CommandLineArgs, String> {
    // Check if user provided minimum number of arguments (program name + at least one file)
    if args.len() < 2 {
        return Err("Insufficient arguments".to_string()); // Return error if no arguments provided
    }

    // Initialize variables to store parsed arguments with default values
    let mut mesh_file: Option<String> = None;        // No mesh file found yet
    let mut node_txt_file: Option<String> = None;    // No node txt file found yet
    let mut mesh_entity = MeshEntity::Node;          // Default to node values (will be overridden if flags provided)
    let mut print_vtu = false;                        // Default to not printing VTU content
    let mut txt_files: Vec<String> = Vec::new();      // Collect all txt files to process based on mesh_entity later
    let mut analyze_mesh = false;                     // Default to not analyzing mesh quality
    let mut optimize_gauss = false;                   // Default to not optimizing Gaussian quadrature
    let mut gauss_tolerance = 1e-6;                   // Default tolerance for Gaussian quadrature
    let mut integration_type = IntegrationType::Mass; // Default integration type
    let mut material_property = MaterialProperty::Scalar(vec![1.0]); // Default material property polynomial coefficients in monomial form

    // Start parsing from index 1 (skip program name at index 0)
    let mut i = 1;
    while i < args.len() {                            // Loop through all arguments
        match args[i].as_str() {                      // Convert String to &str for pattern matching
            "--node-values" => {                      // User specified node values flag
                mesh_entity = MeshEntity::Node;         // Set parsing mode to node values
            },
            "--print-vtu" => {                        // User wants VTU content printed
                print_vtu = true;                     // Enable VTU content printing
            },
            "--analyze-mesh" => {                     // User wants mesh quality analysis  
                analyze_mesh = true;                  // Enable mesh analysis 
            },
            "--optimize-gauss" => {                   // User wants Gaussian quadrature optimization
                optimize_gauss = true;                // Enable Gauss optimization
            },
            "--gauss-tolerance" => {                  // User specified tolerance for Gauss optimization
                i += 1; // Move to next argument to read the tolerance value
                if i >= args.len() {                  // Check if value is missing
                    return Err("Missing value for --gauss-tolerance".to_string());
                }
                // Parse the tolerance value from string to f64, return error if invalid
                gauss_tolerance = args[i].parse::<f64>()
                    .map_err(|_| format!("Invalid tolerance value: {}", args[i]))?;
                // Validate that tolerance is positive (required for numerical methods)
                if gauss_tolerance <= 0.0 {
                    return Err("Tolerance must be positive".to_string());
                }
            },
            "--mass-integration" => {                 // User wants mass matrix integration
                integration_type = IntegrationType::Mass;
            },
            "--stiffness-integration" => {            // User wants stiffness matrix integration
                integration_type = IntegrationType::Stiffness;
            },
            "--material-property" => {                // User specified material property coefficients
                i += 1; // Move to next argument to read the coefficient values
                if i >= args.len() {                  // Check if values are missing
                    return Err("Missing value for --material-property".to_string());
                }
                
                // Parse based on integration type - material properties differ for mass vs stiffness
                material_property = match integration_type {
                    IntegrationType::Mass => {
                        // For mass integration: scalar polynomial coefficients (density, etc.)
                        // Format: comma-separated coefficients representing polynomial terms
                        let coefficients: Vec<f64> = args[i].split(',')
                            .map(|s| s.trim().parse::<f64>())
                            .collect::<Result<Vec<f64>, _>>()
                            .map_err(|_| format!("Invalid material property coefficients: {}", args[i]))?;
                        if coefficients.is_empty() {
                            return Err("Material property coefficients cannot be empty".to_string());
                        }
                        MaterialProperty::Scalar(coefficients)
                    }
                    IntegrationType::Stiffness => {
                        // For stiffness integration: tensor/matrix where each component is a polynomial
                        // New Format: "size;c00_coeffs;c01_coeffs;...;cNN_coeffs"
                        // Each component's coefficients are comma-separated polynomial terms
                        // Example: "3;1.0,0.5;0.2;0.0;0.2;1.0,0.3;0.0;0.0;0.0;0.5"
                        
                        let parts: Vec<&str> = args[i].split(';').collect();
                        if parts.is_empty() {
                            return Err("For stiffness matrix, provide: size;c00;c01;...".to_string());
                        }
                        
                        // First value is the matrix size (dimension-dependent: 1 for 1D, 3 for 2D, 6 for 3D)
                        let matrix_size: usize = parts[0].trim().parse()
                            .map_err(|_| format!("Invalid matrix size: {}", parts[0]))?;
                        
                        // Validate matrix size based on dimensionality
                        // 1D: 1x1, 2D: 3x3 (xx, xy, yy for plane stress/strain), 3D: 6x6 (full 3D elasticity)
                        if matrix_size != 1 && matrix_size != 3 && matrix_size != 6 {
                            return Err(format!(
                                "Matrix size must be 1 (1D), 3 (2D), or 6 (3D), got {}",
                                matrix_size
                            ));
                        }
                        
                        // Total number of matrix components (size x size)
                        let num_components = matrix_size * matrix_size;
                        
                        // Remaining parts are polynomial coefficients for each matrix component
                        let component_parts = &parts[1..];
                        
                        if component_parts.is_empty() {
                            return Err("No coefficients provided for stiffness matrix".to_string());
                        }
                        
                        // Verify we have exactly the right number of components for the matrix
                        if component_parts.len() != num_components {
                            return Err(format!(
                                "Expected {} components ({}x{} matrix), got {}. Format: size;c00;c01;...;cNN",
                                num_components, matrix_size, matrix_size, component_parts.len()
                            ));
                        }
                        
                        // Parse each component's polynomial coefficients into 2D structure
                        // material_tensor[i][j] = vector of polynomial coefficients for component (i,j)
                        let mut material_tensor = vec![vec![vec![]; matrix_size]; matrix_size];
                        
                        for (idx, component_str) in component_parts.iter().enumerate() {
                            // Parse comma-separated coefficients for this component
                            let coefficients: Vec<f64> = component_str
                                .split(',')
                                .map(|s| s.trim().parse::<f64>())
                                .collect::<Result<Vec<f64>, _>>()
                                .map_err(|_| format!("Invalid coefficient in component {}: '{}'", idx, component_str))?;
                            
                            if coefficients.is_empty() {
                                return Err(format!("Component {} has no coefficients", idx));
                            }
                            
                            // Convert flat index to (i, j) matrix indices using row-major order
                            // Row-major: idx = i * matrix_size + j
                            let i = idx / matrix_size;
                            let j = idx % matrix_size;
                            
                            material_tensor[i][j] = coefficients;
                        }
                        
                        MaterialProperty::Matrix(material_tensor)
                    }
                };
            },
            file if file.ends_with(".mphtxt") || file.ends_with(".inp") => { // Check if file is a mesh file
                if mesh_file.is_some() {              // Check if we already found a mesh file
                    return Err("Multiple mesh files specified".to_string()); // Error: only one mesh file allowed
                }
                mesh_file = Some(file.to_string());   // Store the mesh file path
            },
            file if file.ends_with(".txt") => {      // Check if file is a txt result file
                txt_files.push(file.to_string());     // Add to list of txt files for later processing
            },
            _ => {                                    // Unknown argument or unsupported file type
                return Err(format!("Unknown argument or unsupported file: {}", args[i])); // Return descriptive error
            }
        }
        i += 1; // Move to next argument
    }

    // Validate that we always have a mesh file (this is mandatory for all operations)
    if mesh_file.is_none() {
        return Err("Mesh file (.mphtxt or .inp) is required".to_string()); // Error if no mesh file provided
    }

    // Process txt files based on the specified mesh entity and validate file count requirements
    match mesh_entity {
        MeshEntity::Node => {                          // User wants to parse node values
            if txt_files.len() > 1 {                  // Check if more than one txt file provided
                return Err("Only one txt file allowed for --node-values".to_string()); // Error: too many files
            }
            if !txt_files.is_empty() {                // If at least one txt file was provided
                node_txt_file = Some(txt_files[0].clone()); // Assign first (and only) txt file as node file
            }
        }
    }

    // Validate that if txt files are provided, the user must specify what type they contain
    // This prevents ambiguity about how to parse the txt files
    if !txt_files.is_empty() && !args.iter().any(|arg| 
        arg.starts_with("--") && (arg.contains("values") || arg == "--print-vtu" || arg == "--analyze-mesh" || arg == "--optimize-gauss")
    ) { 
        return Err("When parsing txt files, you must specify --node-values".to_string());
    }

    // Return successfully parsed and validated arguments
    Ok(CommandLineArgs {
        mesh_file,        // Required mesh file path
        node_txt_file,    // Optional node values file path
        mesh_entity,       // Specified mesh entity parsing mode
        print_vtu,        // VTU printing preference
        analyze_mesh,     // Mesh analysis preference
        optimize_gauss,   // Gaussian quadrature optimization preference
        gauss_tolerance,  // Tolerance for Gauss optimization
        integration_type, // Integration type preference
        material_property,// Material property polynomial coefficients
    })
}

/// Print usage information to help users understand how to use the program
/// Shows required format, optional parameters, and example commands
fn print_usage(program_name: &str) {
    eprintln!("Usage: {} <mesh_file> [value_files...] [OPTIONS]", program_name);
    eprintln!("Required: mesh.mphtxt or model.inp");
    eprintln!("Optional: nodes.txt");
    eprintln!("Options:");
    eprintln!("  --node-values         Parse txt file as node values");
    eprintln!("  --print-vtu           Print VTU file content to console");
    eprintln!("  --analyze-mesh        Perform mesh quality analysis");
    eprintln!("  --optimize-gauss      Optimize Gaussian quadrature points");
    eprintln!("  --gauss-tolerance <val> Set tolerance for Gauss optimization (default: 1e-6)");
    eprintln!("  --mass-integration    Use mass matrix integration (default)");
    eprintln!("  --stiffness-integration Use stiffness matrix integration");
    eprintln!("  --material-property <spec> Material property specification:");
    eprintln!("     For mass: comma-separated polynomial coefficients (default: 1.0)");
    eprintln!("        Example: --material-property 1.0,0.5,0.1");
    eprintln!("     For stiffness: size;c00;c01;...;cNN where each component is polynomial coefficients");
    eprintln!("        Example: --material-property 2;1.0;0.0;0.0;1.0");
    eprintln!(); // Empty line for separation
    eprintln!("Note: Mesh file is always required. Value type must be specified when using txt files."); // Important usage note
    eprintln!(); // Empty line for separation
    eprintln!("Examples:"); // Section header for usage examples
    eprintln!("  {} mesh.mphtxt                              # Parse mesh only", program_name); // Basic usage example
    eprintln!("  {} mesh.mphtxt nodes.txt --node-values      # Parse mesh + node values", program_name); // Mesh + nodes example
    eprintln!("  {} mesh.mphtxt nodes.txt --node-values --print-vtu # Parse all with output", program_name); // Full featured example
    eprintln!("  {} mesh.mphtxt --analyze-mesh               # Parse mesh + quality analysis", program_name); // Analysis example  
    eprintln!("  {} mesh.mphtxt nodes.txt --node-values --analyze-mesh # Parse + analyze both meshes", program_name); // Full analysis example 
    eprintln!("  {} mesh.mphtxt --optimize-gauss --mass-integration --material-property 2.5", program_name); // Mass integration with scalar
    eprintln!("  {} mesh.mphtxt --optimize-gauss --mass-integration --material-property 1.0,0.5", program_name); // Mass integration with polynomial
    eprintln!("  {} mesh.mphtxt --optimize-gauss --gauss-tolerance 1e-8 # Custom tolerance", program_name); // Custom tolerance example
    eprintln!("  {} mesh.mphtxt --analyze-mesh --optimize-gauss # Both analysis and optimization", program_name); // Combined operations
    eprintln!("  {} mesh.mphtxt --optimize-gauss --stiffness-integration --material-property '3;2.0,2.0;0.5,0.0;0.5;1.0;0.0;0.0,1.2;0.7;1.0,1.2,2.0;1.0'", program_name); // Stiffness example
}

/// Perform mesh quality analysis and write results to text file 
/// Calculates geometric quality metrics for each element and generates a report
/// 
/// # Arguments
/// * `mesh_data` - The mesh to analyze
/// * `mesh_name` - Base name for output file
/// * `output_dir` - Directory where analysis file will be written
/// 
/// # Returns
/// * `Ok(())` if analysis succeeds
/// * `Err(String)` with error message if analysis or file writing fails
fn analyze_mesh_quality(mesh_data: &MeshData, mesh_name: &str, output_dir: &Path) -> Result<(), String> {
    println!("Performing mesh quality analysis...");
    
    // Construct output file path with descriptive name
    let output_path = output_dir.join(format!("{}_mesh_quality_analysis.txt", mesh_name));
    
    // Perform the geometric analysis using the specialized analysis module
    match GeometricAnalysis::analyse_mesh_quality(mesh_data) {
        Ok(quality_report) => {
            // Create output file for writing analysis results
            let mut file = File::create(&output_path)
                .map_err(|e| format!("Failed to create analysis file: {}", e))?;
            
            // Write the quality report in debug format for human readability
            writeln!(file, "{:#?}", quality_report)
                .map_err(|e| format!("Failed to write analysis results: {}", e))?;
            
            println!("Mesh quality analysis written to: {}", output_path.display());
            Ok(())
        },
        Err(e) => {
            // Return descriptive error message if analysis failed
            Err(format!("Mesh analysis failed: {:?}", e))
        }
    }
}

/// Perform Gaussian quadrature optimization and write results to text file
/// Finds optimal number of Gauss points for numerical integration on each element
/// 
/// # Arguments
/// * `mesh_data` - The mesh containing elements to optimize
/// * `integration_type` - Type of integration (Mass or Stiffness)
/// * `tolerance` - Convergence tolerance for optimization
/// * `material_property` - Material property (scalar or matrix with polynomial coefficients)
/// * `output_dir` - Directory where optimization results will be written
/// * `mesh_name` - Base name for output file
/// 
/// # Returns
/// * `Ok(())` if optimization succeeds
/// * `Err(String)` with error message if optimization or file writing fails
fn optimize_gaussian_quadrature(
    mesh_data: &MeshData,
    integration_type: IntegrationType,
    tolerance: f64,
    material_property: &MaterialProperty,
    output_dir: &Path,
    mesh_name: &str,
) -> Result<(), String> {
    println!("Performing Gaussian quadrature optimization...");
    
    // Construct output file path with descriptive name
    let output_path = output_dir.join(format!("{}_gauss_optimization.txt", mesh_name));
    
    // Find optimal Gauss point numbers for all elements in the mesh
    match GaussianQuadrature::find_optimal_gauss_points_number_mesh(
        mesh_data, integration_type, tolerance, material_property
    ) {
        Ok(optimization_results) => {
            // Create output file for writing optimization results
            let mut file = File::create(&output_path)
                .map_err(|e| format!("Failed to create optimization file: {}", e))?;
            
            // Write the optimization results in debug format for human readability
            writeln!(file, "{:#?}", optimization_results)
                .map_err(|e| format!("Failed to write optimization results: {}", e))?;
            
            println!("Gaussian quadrature optimization completed successfully");
            println!("Successfully analyzed {} elements", optimization_results.total_elements);
            Ok(())
        },
        Err(e) => {
            // Return descriptive error message if optimization failed
            let error_msg = format!("Gaussian quadrature optimization failed: {}", e);
            println!("{}", error_msg);
            Err(error_msg)
        }
    }
}

/// Write parsed mesh data to text file in debug format
/// Useful for verifying mesh parsing and inspecting mesh structure
/// 
/// # Arguments
/// * `mesh_data` - The mesh data to write
/// * `output_dir` - Directory where file will be written
/// * `mesh_name` - Base name for output file
/// 
/// # Returns
/// * `Ok(())` if writing succeeds
/// * `Err(String)` with error message if file creation or writing fails
fn write_mesh_data_to_file(mesh_data: &MeshData, output_dir: &Path, mesh_name: &str) -> Result<(), String> {
    println!("Writing parsed mesh data to file...");
    
    // Construct output file path with descriptive name
    let output_path = output_dir.join(format!("{}_parsed_mesh_data.txt", mesh_name));
    
    // Create output file
    let mut file = File::create(&output_path)
        .map_err(|e| format!("Failed to create mesh data file: {}", e))?;
    
    // Write mesh data in pretty-printed debug format for readability
    writeln!(file, "{:#?}", mesh_data)
        .map_err(|e| format!("Failed to write mesh data: {}", e))?;
    
    println!("Parsed mesh data written to: {}", output_path.display());
    Ok(())
}

/// Write parsed field data to text files in debug format
/// Useful for verifying field value parsing and inspecting node/element data
/// 
/// # Arguments
/// * `node_data` - Optional node value data to write
/// * `output_dir` - Directory where file will be written
/// * `mesh_name` - Base name for output file
/// 
/// # Returns
/// * `Ok(())` if writing succeeds (or if no data to write)
/// * `Err(String)` with error message if file creation or writing fails
fn write_field_data_to_file(
    node_data: &Option<NodeValueData>,
    output_dir: &Path, 
    mesh_name: &str
) -> Result<(), String> {
    // Only write if node data was provided
    if let Some(node_data) = node_data {
        // Construct output file path with descriptive name
        let output_path = output_dir.join(format!("{}_parsed_node_data.txt", mesh_name));
        
        // Create output file
        let mut file = File::create(&output_path)
            .map_err(|e| format!("Failed to create node data file: {}", e))?;
        
        // Write node data in pretty-printed debug format for readability
        writeln!(file, "{:#?}", node_data)
            .map_err(|e| format!("Failed to write node data: {}", e))?;
        
        println!("Parsed node data written to: {}", output_path.display());
    }
    
    Ok(())
}

/// Main entry point - the function that runs when the program starts
/// Orchestrates the entire workflow:
/// 1. Parse and validate command line arguments
/// 2. Load and parse mesh file (.mphtxt or .inp)
/// 3. Optionally load node field values from .txt file
/// 4. Optionally perform mesh quality analysis
/// 5. Optionally optimize Gaussian quadrature points
/// 6. Generate VTU output files for visualization
fn main() {
    // ========== STEP 1: PARSE COMMAND LINE ARGUMENTS ==========
    // Collect all command line arguments into a vector of strings
    // args[0] is always the program name, args[1+] are user-provided arguments
    let args: Vec<String> = env::args().collect();
    
    // Parse and validate command line arguments using our custom parser
    let cmd_args = match parse_arguments(args.clone()) {
        Ok(args) => args,           // If parsing succeeded, use the parsed arguments
        Err(e) => {                 // If parsing failed, show error and exit
            eprintln!("Error: {}", e);    // Print the specific error message
            print_usage(&args[0]);        // Show usage information to help user
            process::exit(1);             // Exit with error code 1 to indicate failure
        }
    };

    // ========== STEP 2: SET UP OUTPUT DIRECTORY ==========
    // Set up output directory where generated VTU files will be saved
    // This is a hardcoded path but could be made configurable in the future
    let output_dir = PathBuf::from(r"C:\Users\caybe\OneDrive\Vorlesungen\Master\Master Thesis\Master-Thesis-Caglar-Ayberk-Saglik\test_output");
    
    // Ensure the output directory exists, create it if it doesn't
    if !output_dir.exists() {                                   // Check if directory already exists
        if let Err(e) = std::fs::create_dir_all(&output_dir) {  // Try to create directory and all parent directories
            eprintln!("Error creating output directory: {:?}", e); // Print error if creation fails
            process::exit(1); // Exit with error code since we can't proceed without output directory
        }
    }

    // ========== STEP 3: PARSE MESH FILE ==========
    // Parse mesh file (always required) - unwrap is safe because we validated it exists
    let mesh_data = parse_mesh_file(cmd_args.mesh_file.as_ref().unwrap());

    // ========== STEP 4: PARSE OPTIONAL FIELD VALUE FILES ==========
    // Parse txt files for field values if provided - returns node_data
    let node_value_data = parse_txt_files(&cmd_args.node_txt_file);

    // ========== STEP 5: PRINT PARSING SUMMARY ==========
    // Print simple summaries of parsed data to inform user of successful parsing
    println!("{}", "=".repeat(50));                             // Print visual separator line
    println!("Mesh parsed successfully:");
    println!("  - Nodes: {}", mesh_data.num_nodes);
    println!("  - Elements: {}", mesh_data.elements.len());
    println!("  - Dimension: {}", mesh_data.dimension);
    println!("  - Material property coefficients: {:?}", cmd_args.material_property);

    // If node field values were parsed, show simple summary
    if node_value_data.is_some() {
        println!("Field values parsed: Node values: Yes");
    }

    // ========== STEP 6: GENERATE OUTPUT FILE NAMES ==========
    // Generate mesh name for output files by extracting filename without extension
    let mesh_name = cmd_args.mesh_file
        .as_ref()                                                // Get reference to Option<String>
        .map(|f| Path::new(f).file_stem().unwrap_or_default().to_string_lossy().to_string()) // Extract filename without extension
        .unwrap_or_else(|| "output".to_string());               // Fallback name if somehow no mesh file (shouldn't happen)

    // ========== STEP 7: WRITE DEBUG OUTPUT FILES ==========
    // Write parsed data to text files for inspection and debugging
    if let Err(e) = write_mesh_data_to_file(&mesh_data, &output_dir, &mesh_name) {
        eprintln!("Warning: Failed to write mesh data to file: {}", e);
    }

    if let Err(e) = write_field_data_to_file(&node_value_data, &output_dir, &mesh_name) {
        eprintln!("Warning: Failed to write field data to file: {}", e);
    }

    // ========== STEP 8: INITIALIZE ANALYSIS REPORT STRUCTURES ==========
    // Initialize empty report structures that will be filled if analysis is performed
    let mut mesh_quality_report = MeshQualityReport {
        total_elements: 0,
        element_qualities: Vec::new(),
    };
    
    let mut gauss_report = GaussianPointNumberReport {
        total_elements: 0,
        gauss_point_numbers: Vec::new(),
    };

    // ========== STEP 9: OPTIONAL MESH QUALITY ANALYSIS ==========
    // MESH QUALITY ANALYSIS SECTION
    if cmd_args.analyze_mesh {
        // Perform analysis and write results to file
        if let Err(e) = analyze_mesh_quality(&mesh_data, &mesh_name, &output_dir) {
            eprintln!("Warning: Mesh analysis failed: {}", e);
        } else {
            // Get the report for inclusion in VTU file
            match GeometricAnalysis::analyse_mesh_quality(&mesh_data) {
                Ok(report) => mesh_quality_report = report,
                Err(e) => eprintln!("Warning: Could not get mesh quality report for VTU: {:?}", e),
            }
        }
    }

    // ========== STEP 10: OPTIONAL GAUSSIAN QUADRATURE OPTIMIZATION ==========
    // GAUSSIAN QUADRATURE OPTIMIZATION SECTION
    if cmd_args.optimize_gauss {
        // Perform optimization and write results to file
        if let Err(e) = optimize_gaussian_quadrature(&mesh_data, cmd_args.integration_type, cmd_args.gauss_tolerance, &cmd_args.material_property, &output_dir, &mesh_name) {
            eprintln!("Warning: Gaussian quadrature optimization failed: {}", e);
        } else {
            // Get the report for inclusion in VTU file
            match GaussianQuadrature::find_optimal_gauss_points_number_mesh(
                &mesh_data, 
                cmd_args.integration_type, 
                cmd_args.gauss_tolerance, 
                &cmd_args.material_property
            ) {
                Ok(report) => gauss_report = report,
                Err(e) => eprintln!("Warning: Could not get Gauss report for VTU: {:?}", e),
            }
        }
    }

    // ========== STEP 11: GENERATE PRIMARY VTU FILE ==========
    // Generate output file name based on input mesh file name
    let output_filename = mesh_name;
    
    // Create full output path by combining directory and filename with .vtu extension
    let output_path = output_dir.join(format!("{}.vtu", output_filename));
    let output_str = output_path.to_str().expect("Invalid UTF-8 in output path"); // Convert PathBuf to &str for writer

    // Inform user about VTU file generation
    println!("\nWriting VTU file: {}", output_path.display());
    
    // Write VTU file with or without field values based on what data is available
    let write_result = match &node_value_data {
        Some(node_data) => {
            // Write VTU with both mesh and node field values
            VTUWriter::write_vtu_with_field_values(&mesh_data, node_data, output_str)
        },
        None => {
            // Write VTU with mesh data only
            VTUWriter::write_vtu(&mesh_data, output_str)
        }
    };

    // Check if VTU file writing was successful and inform user
    match write_result {
        Ok(_) => println!("Successfully wrote VTU file"),      // Success message
        Err(e) => eprintln!("Error writing VTU file: {:?}", e), // Error message if writing fails
    }

    // ========== STEP 12: GENERATE ANALYSIS VTU FILE IF NEEDED ==========
    // If analysis or optimization was performed, create additional VTU with those results
    if cmd_args.analyze_mesh || cmd_args.optimize_gauss {
        let analysis_output_path = output_dir.join(format!("{}_with_analysis.vtu", output_filename));
        let analysis_output_str = analysis_output_path.to_str().expect("Invalid UTF-8 in output path");

        println!("\nWriting analysis VTU file: {}", analysis_output_path.display());
        
        // Write VTU file that includes mesh quality and/or Gauss point number data
        match VTUWriter::write_vtu_with_geometric_analysis_and_gaussian_points_number(
            &mesh_data,
            &mesh_quality_report,
            &gauss_report,
            analysis_output_str,
        ) {
            Ok(_) => println!("Successfully wrote VTU file with analysis data: {}", analysis_output_path.display()),
            Err(e) => eprintln!("Error writing VTU file with analysis: {:?}", e),
        }
    }

    // ========== STEP 13: OPTIONAL CONSOLE OUTPUT ==========
    // Optionally print the VTU file content to terminal (useful for debugging or verification)
    if cmd_args.print_vtu {                                    // Only if user specified --print-vtu flag
        println!("\nVTU FILE CONTENT");                      // Header for file content section
        println!("{}", "=".repeat(50));                      // Visual separator
        
        // Read the generated VTU file and display its contents
        match std::fs::read_to_string(&output_path) {        // Read entire file into a string
            Ok(content) => println!("{}", content),          // Print file content to console
            Err(e) => eprintln!("Error reading VTU file: {:?}", e), // Print error if file read fails
        }

        // If analysis file was generated, print that too
        if cmd_args.analyze_mesh || cmd_args.optimize_gauss {
            let analysis_output_path = output_dir.join(format!("{}_with_analysis.vtu", output_filename));
            println!("\nANALYSIS VTU FILE CONTENT");
            println!("{}", "=".repeat(50));
            match std::fs::read_to_string(&analysis_output_path) {
                Ok(content) => println!("{}", content),
                Err(e) => eprintln!("Error reading analysis VTU file: {:?}", e),
            }
        }
    }

    // ========== STEP 14: COMPLETION MESSAGE ==========
    // Print final completion message to indicate successful program execution
    println!("\nPROCESSING COMPLETED!");
    println!("{}", "=".repeat(50));                          // Visual separator for clean output
}

/// Parse mesh file based on its extension and return structured mesh data
/// Supports both COMSOL (.mphtxt) and Abaqus (.inp) mesh formats
/// 
/// # Arguments
/// * `filename` - Path to the mesh file to parse
/// 
/// # Returns
/// * `MeshData` - Parsed mesh data structure containing nodes and elements
/// 
/// # Exits
/// * Exits with code 1 if file format is unsupported or parsing fails
fn parse_mesh_file(filename: &str) -> MeshData {
    let path = Path::new(filename);                           // Create Path object for easy extension checking
    
    // Inform user about which file is being processed
    println!("Parsing mesh file: {}", filename);
    
    // Select appropriate parser based on file extension
    let mesh_data = match path.extension().and_then(|ext| ext.to_str()) { // Extract file extension as string
        Some("mphtxt") => {                                   // COMSOL mesh file format
            println!("Using COMSOL MPHTXT parser");          // Inform user which parser is being used
            MphtxtParser::parse_file(filename)               // Use COMSOL-specific parser
        }
        Some("inp") => {                                      // Abaqus input file format
            println!("Using Abaqus INP parser");             // Inform user which parser is being used
            AbaqusInpParser::parse_file(filename)            // Use Abaqus-specific parser
        }
        _ => {                                                // Unsupported file extension
            eprintln!("Unsupported mesh file format: {}", filename); // Print error message
            process::exit(1);                                // Exit since we can't parse this format
        }
    };
    
    // Handle the result of the parsing operation
    match mesh_data {
        Ok(data) => data,   // If parsing succeeded, return the mesh data
        Err(e) => {         // If parsing failed
            eprintln!("Error parsing mesh file: {:?}", e);   // Print error details
            process::exit(1); // Exit with error code to indicate failure
        }
    }
}

/// Parse txt files for field values 
/// Currently only supports node value files from COMSOL
/// 
/// # Arguments
/// * `node_txt_file` - Optional path to txt file containing node field values
/// 
/// # Returns
/// * `Some(NodeValueData)` if file was provided and parsed successfully
/// * `None` if no file was provided
/// 
/// # Exits
/// * Exits with code 1 if file parsing fails
fn parse_txt_files(node_txt_file: &Option<String>) -> Option<NodeValueData> {
    if let Some(node_file) = node_txt_file {
        // Inform user about parsing operation
        println!("Parsing node values from: {}", node_file);
        
        // Attempt to parse the node values file
        match ComsolTxtNodeParser::parse_file(node_file) {
            Ok(data) => {
                println!("Successfully parsed node values");
                Some(data)  // Return parsed data
            },
            Err(e) => {
                // Print error and exit if parsing fails
                eprintln!("Error parsing node values from {}: {:?}", node_file, e);
                process::exit(1);
            }
        }
    } else {
        None  // No file provided, return None
    }
}
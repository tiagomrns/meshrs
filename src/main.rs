// Import standard library modules for system interaction
use std::env;                    // Provides access to command line arguments via env::args()
use std::process;               // Enables process control, particularly process::exit() for terminating with error codes
use std::path::{Path, PathBuf}; // Path manipulation utilities - Path for immutable path references, PathBuf for owned paths
use std::collections::HashMap;  // For storing Gaussian quadrature optimization results
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
use parser::comsol::txt_element::ComsolTxtElementParser; // Parser for element field values for COMSOL .txt result files

mod writer;                     // Module containing output file writers
use writer::xml_writer::VTUWriter; // Writer for VTK XML unstructured grid (.vtu) format

mod mesh_analysis;                   // Module containing mesh quality analysis  
use mesh_analysis::geometric_analysis::GeometricAnalysis; // Geometric analysis for mesh quality 
use mesh_analysis::gaussian_quadrature::*; // Gauss quadrature utilities for numerical integration 

/// Struct to hold parsed command line arguments in a structured format
/// This makes it easier to pass around and validate the user's input
struct CommandLineArgs {
    mesh_file: Option<String>,        // Path to mesh file (.mphtxt or .inp) - always required but Option for parsing stage
    node_txt_file: Option<String>,    // Path to txt file containing node field values - optional
    element_txt_file: Option<String>, // Path to txt file containing element field values - optional
    mesh_entity: MeshEntity,            // Specifies which type of values to parse from txt files
    print_vtu: bool,                  // Flag to control whether to print VTU file content to console
    analyze_mesh: bool,               // Flag to enable mesh quality analysis  
    optimize_gauss: bool,             // Flag to enable Gaussian quadrature optimization
    gauss_tolerance: f64,             // Tolerance for Gaussian quadrature optimization
    integration_type: IntegrationType, // Type of integration (Mass or Stiffness)
    material_property: Vec<f64>,      // Material property polynomial coefficients in monomial form
}

/// Enum to specify which type of values to parse from txt file(s)
/// This determines the parsing strategy and validation requirements
#[derive(Debug, Clone, PartialEq)] // Derive traits for debugging, cloning, and comparison
enum MeshEntity {
    Node,     // Parse txt file as containing node-based field values
    Element,  // Parse txt file as containing element-based field values
    Both,     // Parse two txt files (first=nodes, second=elements)
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
    let mut element_txt_file: Option<String> = None; // No element txt file found yet
    let mut mesh_entity = MeshEntity::Node;          // Default to node values (will be overridden if flags provided)
    let mut print_vtu = false;                        // Default to not printing VTU content
    let mut txt_files: Vec<String> = Vec::new();      // Collect all txt files to process based on mesh_entity later
    let mut analyze_mesh = false;                     // Default to not analyzing mesh quality
    let mut optimize_gauss = false;                   // Default to not optimizing Gaussian quadrature
    let mut gauss_tolerance = 1e-6;                   // Default tolerance for Gaussian quadrature
    let mut integration_type = IntegrationType::Mass; // Default integration type
    let mut material_property = vec![1.0]; // Default material property polynomial coefficients in monomial form

    // Start parsing from index 1 (skip program name at index 0)
    let mut i = 1;
    while i < args.len() {                            // Loop through all arguments
        match args[i].as_str() {                      // Convert String to &str for pattern matching
            "--node-values" => {                      // User specified node values flag
                mesh_entity = MeshEntity::Node;         // Set parsing mode to node values
            },
            "--element-values" => {                   // User specified element values flag
                mesh_entity = MeshEntity::Element;      // Set parsing mode to element values
            },
            "--both-values" => {                      // User specified both values flag
                mesh_entity = MeshEntity::Both;         // Set parsing mode to both types
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
                i += 1; // Move to next argument for the value
                if i >= args.len() {
                    return Err("Missing value for --gauss-tolerance".to_string());
                }
                gauss_tolerance = args[i].parse::<f64>()
                    .map_err(|_| format!("Invalid tolerance value: {}", args[i]))?;
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
                i += 1; // Move to next argument for the values
                if i >= args.len() {
                    return Err("Missing value for --material-property".to_string());
                }
                material_property = args[i].split(',')
                    .map(|s| s.trim().parse::<f64>())
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| format!("Invalid material property coefficients: {}", args[i]))?;
                if material_property.is_empty() {
                    return Err("Material property coefficients cannot be empty".to_string());
                }
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
        },
        MeshEntity::Element => {                       // User wants to parse element values
            if txt_files.len() > 1 {                  // Check if more than one txt file provided
                return Err("Only one txt file allowed for --element-values".to_string()); // Error: too many files
            }
            if !txt_files.is_empty() {                // If at least one txt file was provided
                element_txt_file = Some(txt_files[0].clone()); // Assign first (and only) txt file as element file
            }
        },
        MeshEntity::Both => {                          // User wants to parse both node and element values
            if txt_files.len() != 2 {                 // Check if exactly 2 txt files provided
                return Err("Exactly two txt files required for --both-values (first for nodes, second for elements)".to_string());
            }
            node_txt_file = Some(txt_files[0].clone());    // First file is for node values
            element_txt_file = Some(txt_files[1].clone()); // Second file is for element values
        }
    }

    // Validate that if txt files are provided, the user must specify what type they contain
    // This prevents ambiguity about how to parse the txt files
    if !txt_files.is_empty() && !args.iter().any(|arg| arg.starts_with("--") && (arg.contains("values") || arg == "--print-vtu" || arg == "--analyze-mesh" || arg == "--optimize-gauss")) { 
        return Err("When parsing txt files, you must specify --node-values, --element-values, or --both-values".to_string());
    }

    // Return successfully parsed and validated arguments
    Ok(CommandLineArgs {
        mesh_file,        // Required mesh file path
        node_txt_file,    // Optional node values file path
        element_txt_file, // Optional element values file path
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
    eprintln!("Usage: {} <mesh_file> [value_files...] [OPTIONS]", program_name); // Show general usage format
    eprintln!(); // Empty line for readability
    eprintln!("Required:"); // Section header for required files
    eprintln!("  mesh.mphtxt           COMSOL mesh file (required)");     // COMSOL mesh file description
    eprintln!("  model.inp             Abaqus mesh file (required)");     // Abaqus mesh file description
    eprintln!(); // Empty line for separation
    eprintln!("Optional Value Files:"); // Section header for optional files
    eprintln!("  nodes.txt             COMSOL node result file");         // Node values file description
    eprintln!("  elements.txt          COMSOL element result file");      // Element values file description
    eprintln!(); // Empty line for separation
    eprintln!("Options:"); // Section header for command line options
    eprintln!("  --node-values         Parse txt file as node values");          // Node values option description
    eprintln!("  --element-values      Parse txt file as element values");       // Element values option description
    eprintln!("  --both-values         Parse two txt files (first=nodes, second=elements)"); // Both values option description
    eprintln!("  --print-vtu           Print VTU file content to console");      // Print option description
    eprintln!("  --analyze-mesh        Perform mesh quality analysis");          // Analysis option description 
    eprintln!("  --optimize-gauss      Optimize Gaussian quadrature points");    // Gauss optimization description
    eprintln!("  --gauss-tolerance <val> Set tolerance for Gauss optimization (default: 1e-6)"); // Tolerance description
    eprintln!("  --mass-integration    Use mass matrix integration (default)");   // Mass integration description
    eprintln!("  --stiffness-integration Use stiffness matrix integration");     // Stiffness integration description
    eprintln!("  --material-property <coeffs> Material property polynomial coefficients as comma-separated values (default: 1.0)"); // Material property description
    eprintln!(); // Empty line for separation
    eprintln!("Note: Mesh file is always required. Value type must be specified when using txt files."); // Important usage note
    eprintln!(); // Empty line for separation
    eprintln!("Examples:"); // Section header for usage examples
    eprintln!("  {} mesh.mphtxt                              # Parse mesh only", program_name); // Basic usage example
    eprintln!("  {} mesh.mphtxt nodes.txt --node-values      # Parse mesh + node values", program_name); // Mesh + nodes example
    eprintln!("  {} mesh.mphtxt elements.txt --element-values # Parse mesh + element values", program_name); // Mesh + elements example
    eprintln!("  {} mesh.mphtxt nodes.txt elements.txt --both-values # Parse mesh + both values", program_name); // Complete example
    eprintln!("  {} model.inp nodes.txt elements.txt --both-values --print-vtu # Parse all with output", program_name); // Full featured example
    eprintln!("  {} mesh.mphtxt --analyze-mesh               # Parse mesh + quality analysis", program_name); // Analysis example  
    eprintln!("  {} mesh.mphtxt nodes.txt --node-values --analyze-mesh # Parse + analyze both meshes", program_name); // Full analysis example 
    eprintln!("  {} mesh.mphtxt --optimize-gauss --mass-integration # Optimize Gauss points for mass matrix", program_name);
    eprintln!("  {} mesh.mphtxt --optimize-gauss --gauss-tolerance 1e-8 # Custom tolerance", program_name);
    eprintln!("  {} mesh.mphtxt --analyze-mesh --optimize-gauss # Both analysis and optimization", program_name);
    eprintln!("  {} mesh.mphtxt --optimize-gauss --material-property 1.0 # Constant density variation", program_name);
}

/// Perform mesh quality analysis and write results to text file 
fn analyze_mesh_quality(mesh_data: &MeshData, mesh_name: &str, output_dir: &Path) -> Result<(), String> {
    println!("Performing mesh quality analysis...");
    
    let output_path = output_dir.join(format!("{}_mesh_quality_analysis.txt", mesh_name));
    
    match GeometricAnalysis::analyse_mesh_quality(mesh_data) {
        Ok(quality_report) => {
            let mut file = File::create(&output_path)
                .map_err(|e| format!("Failed to create analysis file: {}", e))?;
            
            writeln!(file, "{:#?}", quality_report)
                .map_err(|e| format!("Failed to write analysis results: {}", e))?;
            
            println!("Mesh quality analysis written to: {}", output_path.display());
            Ok(())
        },
        Err(e) => {
            Err(format!("Mesh analysis failed: {:?}", e))
        }
    }
}

/// Perform Gaussian quadrature optimization and write results to text file
fn optimize_gaussian_quadrature(
    mesh_data: &MeshData,
    integration_type: IntegrationType,
    tolerance: f64,
    material_property: &[f64],
    output_dir: &Path,
    mesh_name: &str,
) -> Result<(), String> {
    println!("Performing Gaussian quadrature optimization...");
    
    let output_path = output_dir.join(format!("{}_gauss_optimization.txt", mesh_name));
    
    match GaussianQuadrature::find_optimal_gauss_points_number_mesh(mesh_data, integration_type, tolerance, material_property) {
        Ok(optimization_results) => {
            let mut file = File::create(&output_path)
                .map_err(|e| format!("Failed to create optimization file: {}", e))?;
            
            writeln!(file, "{:#?}", optimization_results)
                .map_err(|e| format!("Failed to write optimization results: {}", e))?;
            
            println!("Gaussian quadrature optimization written to: {}", output_path.display());
            Ok(())
        },
        Err(e) => {
            Err(format!("Gaussian quadrature optimization failed: {:?}", e))
        }
    }
}

/// Write parsed mesh data to text file in debug format
fn write_mesh_data_to_file(mesh_data: &MeshData, output_dir: &Path, mesh_name: &str) -> Result<(), String> {
    println!("Writing parsed mesh data to file...");
    
    let output_path = output_dir.join(format!("{}_parsed_mesh_data.txt", mesh_name));
    
    let mut file = File::create(&output_path)
        .map_err(|e| format!("Failed to create mesh data file: {}", e))?;
    
    writeln!(file, "{:#?}", mesh_data)
        .map_err(|e| format!("Failed to write mesh data: {}", e))?;
    
    println!("Parsed mesh data written to: {}", output_path.display());
    Ok(())
}

/// Write parsed field data to text files in debug format
fn write_field_data_to_file(
    node_data: &Option<NodeValueData>, 
    element_data: &Option<ElementValueData>, 
    output_dir: &Path, 
    mesh_name: &str
) -> Result<(), String> {
    if let Some(node_data) = node_data {
        let output_path = output_dir.join(format!("{}_parsed_node_data.txt", mesh_name));
        let mut file = File::create(&output_path)
            .map_err(|e| format!("Failed to create node data file: {}", e))?;
        
        writeln!(file, "{:#?}", node_data)
            .map_err(|e| format!("Failed to write node data: {}", e))?;
        
        println!("Parsed node data written to: {}", output_path.display());
    }
    
    if let Some(element_data) = element_data {
        let output_path = output_dir.join(format!("{}_parsed_element_data.txt", mesh_name));
        let mut file = File::create(&output_path)
            .map_err(|e| format!("Failed to create element data file: {}", e))?;
        
        writeln!(file, "{:#?}", element_data)
            .map_err(|e| format!("Failed to write element data: {}", e))?;
        
        println!("Parsed element data written to: {}", output_path.display());
    }
    
    Ok(())
}

/// Main entry point - the function that runs when the program starts
/// Coordinates the entire parsing and conversion process from command line to VTU output
fn main() {
    // Collect all command line arguments into a vector of strings
    // args[0] is always the program name, args[1+] are user-provided arguments
    let args: Vec<String> = env::args().collect();
    
    // Parse and validate command line arguments using our custom parser
    let cmd_args = match parse_arguments(args.clone()) {
        Ok(args) => args,           // If parsing succeeded, use the parsed arguments
        Err(e) => {                 // If parsing failed, show error and exit
            eprintln!("Error: {}", e);    // Print the specific error message
            print_usage(&args[0]);        // Show usage information to help user
            process::exit(1);             // Exit with error code 1
        }
    };

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

    // Parse mesh file (always required) - unwrap is safe because we validated it exists
    let mesh_data = parse_mesh_file(cmd_args.mesh_file.as_ref().unwrap());

    // Parse txt files for field values if provided - returns (node_data, element_data) tuple
    let (node_value_data, element_value_data) = parse_txt_files(&cmd_args.node_txt_file, &cmd_args.element_txt_file);

    // Print simple summaries of parsed data
    println!("{}", "=".repeat(50));                             // Print visual separator line
    println!("Mesh parsed successfully:");
    println!("  - Nodes: {}", mesh_data.num_nodes);
    println!("  - Elements: {}", mesh_data.elements.len());
    println!("  - Dimension: {}", mesh_data.dimension);
    println!("  - Material property coefficients: {:?}", cmd_args.material_property);

    // If any field values were parsed, show simple summary
    if node_value_data.is_some() || element_value_data.is_some() {
        println!("Field values parsed:");
        if node_value_data.is_some() {
            println!("  - Node values: Yes");
        }
        if element_value_data.is_some() {
            println!("  - Element values: Yes");
        }
    }

    // Generate mesh name for output files
    let mesh_name = cmd_args.mesh_file
        .as_ref()                                                // Get reference to Option<String>
        .map(|f| Path::new(f).file_stem().unwrap_or_default().to_string_lossy().to_string()) // Extract filename without extension
        .unwrap_or_else(|| "output".to_string());               // Fallback name if somehow no mesh file (shouldn't happen)

    // Write parsed data to text files for inspection
    if let Err(e) = write_mesh_data_to_file(&mesh_data, &output_dir, &mesh_name) {
        eprintln!("Warning: Failed to write mesh data to file: {}", e);
    }

    if let Err(e) = write_field_data_to_file(&node_value_data, &element_value_data, &output_dir, &mesh_name) {
        eprintln!("Warning: Failed to write field data to file: {}", e);
    }

    let mut mesh_quality_report = MeshQualityReport {
        total_elements: 0,
        element_qualities: Vec::new(),
    };
    
    let mut gauss_report = GaussianPointNumberReport {
        total_elements: 0,
        gauss_point_numbers: Vec::new(),
    };

    // MESH QUALITY ANALYSIS SECTION
    if cmd_args.analyze_mesh {
        if let Err(e) = analyze_mesh_quality(&mesh_data, &mesh_name, &output_dir) {
            eprintln!("Warning: Mesh analysis failed: {}", e);
        } else {
            match GeometricAnalysis::analyse_mesh_quality(&mesh_data) {
                Ok(report) => mesh_quality_report = report,
                Err(e) => eprintln!("Warning: Could not get mesh quality report for VTU: {:?}", e),
            }
        }
    }

    // GAUSSIAN QUADRATURE OPTIMIZATION SECTION
    if cmd_args.optimize_gauss {
        if let Err(e) = optimize_gaussian_quadrature(&mesh_data, cmd_args.integration_type, cmd_args.gauss_tolerance, &cmd_args.material_property, &output_dir, &mesh_name) {
            eprintln!("Warning: Gaussian quadrature optimization failed: {}", e);
        } else {
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

    // Generate output file name based on input mesh file name
    let output_filename = mesh_name;
    
    // Create full output path by combining directory and filename with .vtu extension
    let output_path = output_dir.join(format!("{}.vtu", output_filename));
    let output_str = output_path.to_str().expect("Invalid UTF-8 in output path"); // Convert PathBuf to &str for writer

    // Inform user about VTU file generation
    println!("\nWriting VTU file: {}", output_path.display());
    
    // Write VTU file with or without field values based on what data is available
    let write_result = match (&node_value_data, &element_value_data) {
        (Some(node_data), Some(element_data)) => {              // Both node and element values available
            VTUWriter::write_vtu_with_field_values(&mesh_data, node_data, element_data, output_str)
        },
        (Some(node_data), None) => {                            // Only node values available
            // Create empty element data structure to satisfy writer interface
            let empty_element_data = ElementValueData {
                element_values: Vec::new(),                      // Empty values array
                element_value_type_info: Vec::new(),             // Empty type info array
            };
            VTUWriter::write_vtu_with_field_values(&mesh_data, node_data, &empty_element_data, output_str)
        },
        (None, Some(element_data)) => {                         // Only element values available
            // Create empty node data structure to satisfy writer interface
            let empty_node_data = NodeValueData {
                node_values: Vec::new(),                         // Empty values array
                node_value_type_info: Vec::new(),                // Empty type info array
            };
            VTUWriter::write_vtu_with_field_values(&mesh_data, &empty_node_data, element_data, output_str)
        },
        (None, None) => {                                       // No field values, just mesh geometry
            VTUWriter::write_vtu(&mesh_data, output_str)        // Use simpler writer method for mesh only
        }
    };

    // Check if VTU file writing was successful and inform user
    match write_result {
        Ok(_) => println!("Successfully wrote VTU file"),      // Success message
        Err(e) => eprintln!("Error writing VTU file: {:?}", e), // Error message if writing fails
    }

    if cmd_args.analyze_mesh || cmd_args.optimize_gauss {
        let analysis_output_path = output_dir.join(format!("{}_with_analysis.vtu", output_filename));
        let analysis_output_str = analysis_output_path.to_str().expect("Invalid UTF-8 in output path");

        println!("\nWriting analysis VTU file: {}", analysis_output_path.display());
        
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

    // Optionally print the VTU file content to terminal (useful for debugging or verification)
    if cmd_args.print_vtu {                                    // Only if user specified --print-vtu flag
        println!("\nVTU FILE CONTENT");                      // Header for file content section
        println!("{}", "=".repeat(50));                      // Visual separator
        
        // Read the generated VTU file and display its contents
        match std::fs::read_to_string(&output_path) {        // Read entire file into a string
            Ok(content) => println!("{}", content),          // Print file content to console
            Err(e) => eprintln!("Error reading VTU file: {:?}", e), // Print error if file read fails
        }

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

    // Print final completion message to indicate successful program execution
    println!("\nPROCESSING COMPLETED!");
    println!("{}", "=".repeat(50));                          // Visual separator for clean output
}

/// Parse mesh file based on its extension and return structured mesh data
/// Supports both COMSOL (.mphtxt) and Abaqus (.inp) mesh formats
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
            process::exit(1); // Exit with error code
        }
    }
}

/// Parse txt files for field values using appropriate parsers
/// Takes optional file paths and returns tuple of (node_data, element_data)
fn parse_txt_files(node_txt_file: &Option<String>, element_txt_file: &Option<String>) -> (Option<NodeValueData>, Option<ElementValueData>) {
    let mut node_value_data = None;       // Initialize as None - will be populated if file provided
    let mut element_value_data = None;    // Initialize as None - will be populated if file provided
    
    // Parse node values if file path provided
    if let Some(node_file) = node_txt_file {                 // Check if node file path exists
        println!("Parsing node values from: {}", node_file); // Inform user about node parsing
        match ComsolTxtNodeParser::parse_file(node_file) {   // Attempt to parse node values
            Ok(data) => {                                     // If parsing succeeded
                println!("Successfully parsed node values"); // Success message
                node_value_data = Some(data);                 // Store parsed node data
            },
            Err(e) => {                                       // If parsing failed
                eprintln!("Error parsing node values from {}: {:?}", node_file, e); // Error message with filename
                process::exit(1);                             // Exit since we can't continue with invalid data
            }
        }
    }
    
    // Parse element values if file path provided
    if let Some(element_file) = element_txt_file {           // Check if element file path exists
        println!("Parsing element values from: {}", element_file); // Inform user about element parsing
        match ComsolTxtElementParser::parse_file(element_file) {   // Attempt to parse element values
            Ok(data) => {                                     // If parsing succeeded
                println!("Successfully parsed element values"); // Success message
                element_value_data = Some(data);              // Store parsed element data
            },
            Err(e) => {                                       // If parsing failed
                eprintln!("Error parsing element values from {}: {:?}", element_file, e); // Error message with filename
                process::exit(1);                             // Exit since we can't continue with invalid data
            }
        }
    }
    
    // Return tuple of parsed data (both may be None if no files provided)
    (node_value_data, element_value_data)
}
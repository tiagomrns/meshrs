// Import standard library modules for system interaction
use std::env;                    // Provides access to command line arguments via env::args()
use std::process;               // Enables process control, particularly process::exit() for terminating with error codes
use std::path::{Path, PathBuf}; // Path manipulation utilities - Path for immutable path references, PathBuf for owned paths
use std::collections::HashMap;  // For storing Gaussian quadrature optimization results

// Import custom modules that contain the core functionality
//mod lib;                   // Module containing data structures (MeshData, Node, Element, etc.)
//use lib::*;                // Wildcard import brings all public items from lib module into scope

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
}

/* 
/// Create a deformed mesh by applying nodal displacements to the original mesh 
/// Returns a new MeshData structure with updated node coordinates
fn create_deformed_mesh(original_mesh: &MeshData, node_value_data: &NodeValueData) -> Result<MeshData, String> {
    // Find displacement data in the node values
    let displacement_info = node_value_data.node_value_type_info
        .iter()
        .find(|info| matches!(info.nodevalue_type, ValueType::Displacement))
        .ok_or("No displacement data found in node values")?;

    // Validate that we have displacement data for all mesh nodes
    if displacement_info.num_nodes != original_mesh.num_nodes {
        return Err(format!(
            "Displacement data node count ({}) doesn't match mesh node count ({})",
            displacement_info.num_nodes, original_mesh.num_nodes
        ));
    }

    // Create new mesh data with deformed coordinates
    let mut deformed_mesh = original_mesh.clone();

    // Apply displacements to each node
    for (i, original_node) in original_mesh.nodes.iter().enumerate() {
        // Find corresponding displacement values
        if let Some(displacement_node) = node_value_data.node_values
            .iter()
            .skip(displacement_info.start_index)
            .take(displacement_info.num_nodes)
            .find(|nv| nv.id == original_node.id) {
            
            // Apply displacement to original coordinates
            let mut new_coordinates = original_node.coordinates.clone();
            
            // Add displacement values (ensure we don't exceed coordinate dimensions)
            let disp_components = displacement_node.values.len().min(new_coordinates.len());
            for j in 0..disp_components {
                new_coordinates[j] += displacement_node.values[j];
            }
            
            // Update the deformed mesh node coordinates
            deformed_mesh.nodes[i].coordinates = new_coordinates;
        } else {
            return Err(format!("No displacement data found for node {}", original_node.id));
        }
    }

    Ok(deformed_mesh)
}
*/
 
/// Perform mesh quality analysis and print detailed results 
fn analyze_mesh_quality(mesh_data: &MeshData, mesh_name: &str) -> Result<(), String> {
    println!("\n{}", "=".repeat(60));
    println!("MESH QUALITY ANALYSIS: {}", mesh_name);
    println!("{}", "=".repeat(60));

    match GeometricAnalysis::analyse_mesh_quality(mesh_data) {
        Ok(quality_report) => {
            println!("Analysis completed successfully!");
            println!("Total elements analyzed: {}", quality_report.total_elements);
            
            // Print detailed statistics from the quality report
            println!("\nJacobian Determinant Statistics:");
            println!("  Minimum:     {:.6e}", quality_report.statistics.min_jacobian);
            println!("  Maximum:     {:.6e}", quality_report.statistics.max_jacobian);
            println!("  Average:     {:.6e}", quality_report.statistics.avg_jacobian);
            println!("  Negative:    {} elements ({:.2}%)", 
                quality_report.statistics.negative_jacobian_count,
                (quality_report.statistics.negative_jacobian_count as f64 / quality_report.total_elements as f64) * 100.0
            );

            // Print determinant for every element
            println!("\nDeterminant Results for All Elements:");
            println!("{}", "-".repeat(50));
            for (i, elem_quality) in quality_report.element_qualities.iter().enumerate() {
                println!("Element {}: {:.6e}", elem_quality.element_id, elem_quality.det_jacobian);
                
                // Optional: Add a separator every 10 elements for readability
                if (i + 1) % 10 == 0 && i + 1 < quality_report.element_qualities.len() {
                    println!("{}", "-".repeat(20));
                }
            }
            println!("{}", "-".repeat(50));

            // Additional quality analysis
            if quality_report.statistics.negative_jacobian_count > 0 {
                println!("  WARNING: {} elements have negative Jacobian determinants (inverted elements)!", 
                    quality_report.statistics.negative_jacobian_count);
                
                // Show details of the worst elements
                let mut worst_elements: Vec<_> = quality_report.element_qualities
                    .iter()
                    .filter(|q| q.det_jacobian < 0.0)
                    .collect();
                worst_elements.sort_by(|a, b| a.det_jacobian.partial_cmp(&b.det_jacobian).unwrap());
                
                println!("\n  Worst inverted elements (up to 5):");
                for (i, elem) in worst_elements.iter().take(5).enumerate() {
                    println!("    {}: Element {} - Jacobian det: {:.6e}", 
                        i + 1, elem.element_id, elem.det_jacobian);
                }
            }

            // Quality classification (for console output)
            let very_poor_count = quality_report.element_qualities.iter()
                .filter(|q| q.det_jacobian > 0.0 && q.det_jacobian < 0.1)
                .count();
            let poor_count = quality_report.element_qualities.iter()
                .filter(|q| q.det_jacobian >= 0.1 && q.det_jacobian < 0.5)
                .count();
            let good_count = quality_report.element_qualities.iter()
                .filter(|q| q.det_jacobian >= 0.5)
                .count();

            println!("\nElement Quality Distribution:");
            println!("  Excellent (det ≥ 0.5):     {} elements ({:.1}%)", 
                good_count, (good_count as f64 / quality_report.total_elements as f64) * 100.0);
            println!("  Poor (0.1 ≤ det < 0.5):    {} elements ({:.1}%)", 
                poor_count, (poor_count as f64 / quality_report.total_elements as f64) * 100.0);
            println!("  Very Poor (0 < det < 0.1): {} elements ({:.1}%)", 
                very_poor_count, (very_poor_count as f64 / quality_report.total_elements as f64) * 100.0);
            println!("  Inverted (det < 0):        {} elements ({:.1}%)", 
                quality_report.statistics.negative_jacobian_count,
                (quality_report.statistics.negative_jacobian_count as f64 / quality_report.total_elements as f64) * 100.0);

            // Overall mesh quality assessment
            let quality_percentage = ((quality_report.total_elements - quality_report.statistics.negative_jacobian_count - very_poor_count) as f64 
                / quality_report.total_elements as f64) * 100.0;
            
            println!("\nOverall Mesh Quality Assessment:");
            if quality_percentage >= 90.0 {
                println!("  ✓ EXCELLENT: {:.1}% of elements have acceptable quality", quality_percentage);
            } else if quality_percentage >= 75.0 {
                println!("  ⚠ GOOD: {:.1}% of elements have acceptable quality", quality_percentage);
            } else if quality_percentage >= 50.0 {
                println!("  ⚠ FAIR: {:.1}% of elements have acceptable quality", quality_percentage);
            } else {
                println!("  ✗ POOR: Only {:.1}% of elements have acceptable quality", quality_percentage);
            }

            println!("{}", "=".repeat(60));
            Ok(())
        },
        Err(e) => {
            println!("Analysis failed: {:?}", e);
            println!("{}", "=".repeat(60));
            Err(format!("Mesh analysis failed: {:?}", e))
        }
    }
}

/// Perform Gaussian quadrature optimization and print results
fn optimize_gaussian_quadrature(
    mesh_data: &MeshData,
    integration_type: &IntegrationType,
    tolerance: f64,
) -> Result<(), String> {
    println!("\n{}", "=".repeat(60));
    println!("GAUSSIAN QUADRATURE OPTIMIZATION");
    println!("{}", "=".repeat(60));
    
    println!("Integration type: {:?}", integration_type);
    println!("Error tolerance: {:.2e}", tolerance);
    println!("Analyzing {} elements...", mesh_data.elements.len());
    
    // First, calculate errorless integration points for comparison
    let mut errorless_points_map = HashMap::new();
    let mut errorless_total = 0usize;
    
    for type_info in &mesh_data.element_type_info {
        if matches!(type_info.element_type, ElementType::Vertex) {
            continue;  // Skip vertex elements
        }
        
        match GaussianQuadrature::get_required_polynomial_order(&type_info.element_type, integration_type) {
            Ok(required_poly_order) => {
                let errorless_points = match GaussianQuadrature::for_element(&type_info.element_type, required_poly_order) {
                    Ok(quad_rule) => quad_rule.points.len(),
                    Err(_) => 0,
                };
                
                let start_idx = type_info.start_index;
                let end_idx = start_idx + type_info.num_elements;
                for element_idx in start_idx..end_idx {
                    if element_idx < mesh_data.elements.len() {
                        let element = &mesh_data.elements[element_idx];
                        errorless_points_map.insert(element.id, errorless_points);
                        errorless_total += errorless_points;
                    }
                }
            },
            Err(_) => {
                let start_idx = type_info.start_index;
                let end_idx = start_idx + type_info.num_elements;
                for element_idx in start_idx..end_idx {
                    if element_idx < mesh_data.elements.len() {
                        let element = &mesh_data.elements[element_idx];
                        errorless_points_map.insert(element.id, 0);
                    }
                }
            }
        }
    }
    
    match GaussianQuadrature::optimize_mesh_gauss_points(mesh_data, integration_type, tolerance) {
        Ok(optimization_results) => {
            println!("Optimization completed successfully!");
            
            // Print ALL elements (no truncation, no ellipsis)
            println!("\nERRORLESS vs OPTIMIZED COMPARISON:");
            println!("{}", "=".repeat(60));
            println!("{:<12} | {:<15} | {:<15} | {:<10}", "Element ID", "Errorless Pts", "Optimized Pts", "Savings");
            println!("{}", "-".repeat(60));
            
            let mut total_savings: i32 = 0;
            let mut sorted_elements: Vec<_> = optimization_results.keys().copied().collect();
            sorted_elements.sort_unstable();

            for element_id in sorted_elements {
                let errorless_pts = *errorless_points_map.get(&element_id).unwrap_or(&0);
                let optimized_pts = *optimization_results.get(&element_id).unwrap_or(&0);
                let savings = errorless_pts as i32 - optimized_pts as i32;
                total_savings += savings;
                
                let savings_str = if savings > 0 {
                    format!("-{}", savings)
                } else if savings < 0 {
                    format!("+{}", -savings)
                } else {
                    "0".to_string()
                };
                
                println!("{:<12} | {:<15} | {:<15} | {:<10}", 
                    element_id, errorless_pts, optimized_pts, savings_str);
            }
            
            println!("{}", "-".repeat(60));
            
            // Summary statistics (kept)
            let total_elements = optimization_results.len();
            let optimized_total: usize = optimization_results.values().sum();
            
            println!("\nSUMMARY STATISTICS:");
            println!("{}", "=".repeat(30));
            println!("Total elements:           {}", total_elements);
            println!("Errorless total points:   {}", errorless_total);
            println!("Optimized total points:   {}", optimized_total);
            println!("Total point savings:      {}", total_savings);
            
            let efficiency = if errorless_total > 0 {
                (total_savings as f64 / errorless_total as f64) * 100.0
            } else {
                0.0
            };
            
            if efficiency > 0.0 {
                println!("Computational savings:    {:.1}%", efficiency);
            } else if efficiency < 0.0 {
                println!("Additional cost:          {:.1}%", -efficiency);
            } else {
                println!("No change in points");
            }

            // NOTE: “QUALITY METRICS” section intentionally removed.
            println!("{}", "=".repeat(60));
            Ok(())
        },
        Err(e) => {
            println!("Optimization failed: {:?}", e);
            println!("{}", "=".repeat(60));
            Err(format!("Gaussian quadrature optimization failed: {:?}", e))
        }
    }
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

    // Print summaries of parsed data to inform user about what was processed
    println!("{}", "=".repeat(50));                             // Print visual separator line
    print_mesh_summary(&mesh_data);                             // Display detailed mesh information

    // If any field values were parsed, show their summary as well
    if node_value_data.is_some() || element_value_data.is_some() {
        println!("{}", "=".repeat(50));                         // Print visual separator line
        print_result_summary(&node_value_data, &element_value_data); // Display field values information
    }

    // MESH QUALITY ANALYSIS SECTION
    // Always perform undeformed mesh analysis if requested
    if cmd_args.analyze_mesh {
        // Analyze original (undeformed) mesh
        if let Err(e) = analyze_mesh_quality(&mesh_data, "ORIGINAL MESH") {
            eprintln!("Warning: Original mesh analysis failed: {}", e);
        }

        // FUTURE: Deformed mesh analysis (commented out for now)
        // Uncomment the section below when you want to add deformed mesh analysis
        /*
        // If we have displacement data, also analyze the deformed mesh
        if let Some(ref node_data) = node_value_data {
            // Check if displacement data is available
            let has_displacement = node_data.node_value_type_info
                .iter()
                .any(|info| matches!(info.nodevalue_type, ValueType::Displacement));

            if has_displacement {
                println!("\nDisplacement data detected - analyzing deformed mesh...");
                match create_deformed_mesh(&mesh_data, node_data) {
                    Ok(deformed_mesh) => {
                        if let Err(e) = analyze_mesh_quality(&deformed_mesh, "DEFORMED MESH") {
                            eprintln!("Warning: Deformed mesh analysis failed: {}", e);
                        }
                    },
                    Err(e) => {
                        eprintln!("Warning: Could not create deformed mesh: {}", e);
                    }
                }
            } else {
                println!("\nNote: No displacement data found, skipping deformed mesh analysis.");
            }
        } else {
            println!("\nNote: No node data provided, analyzing undeformed mesh only.");
        }
        */
    }

    // GAUSSIAN QUADRATURE OPTIMIZATION SECTION
    if cmd_args.optimize_gauss {
        // Perform mesh-wide Gaussian quadrature optimization
        if let Err(e) = optimize_gaussian_quadrature(&mesh_data, &cmd_args.integration_type, cmd_args.gauss_tolerance) {
            eprintln!("Warning: Gaussian quadrature optimization failed: {}", e);
        }
    }

    // Generate output file name based on input mesh file name
    let output_filename = cmd_args.mesh_file
        .as_ref()                                                // Get reference to Option<String>
        .map(|f| Path::new(f).file_stem().unwrap_or_default().to_string_lossy().to_string()) // Extract filename without extension
        .unwrap_or_else(|| "output".to_string());               // Fallback name if somehow no mesh file (shouldn't happen)
    
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

    // Optionally print the VTU file content to terminal (useful for debugging or verification)
    if cmd_args.print_vtu {                                    // Only if user specified --print-vtu flag
        println!("\nVTU FILE CONTENT");                      // Header for file content section
        println!("{}", "=".repeat(50));                      // Visual separator
        
        // Read the generated VTU file and display its contents
        match std::fs::read_to_string(&output_path) {        // Read entire file into a string
            Ok(content) => println!("{}", content),          // Print file content to console
            Err(e) => eprintln!("Error reading VTU file: {:?}", e), // Print error if file read fails
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

/// Print a comprehensive summary of parsed result data (from COMSOL .txt files)
/// Shows value types, dimensions, and sample data for verification purposes
fn print_result_summary(node_data: &Option<NodeValueData>, element_data: &Option<ElementValueData>) {
    println!("RESULT DATA SUMMARY");                         // Section header
    println!("{}", "=".repeat(50));                          // Visual separator
    
    // Display node values information if available
    if let Some(node_result) = node_data {                   // Check if node data exists
        println!("NODE VALUES:");                            // Subsection header for node data
        println!("Number of value types: {}", node_result.node_value_type_info.len()); // Different types of values (temp, displacement, etc.)
        
        // Display detailed information about each value type
        println!("\nNode value types:");                     // Subsection for value type details
        for type_info in &node_result.node_value_type_info { // Loop through each value type
            println!("  {:?}:", type_info.nodevalue_type);   // Print the value type name (Temperature, Displacement, etc.)
            println!("    - Dimension: {}", type_info.dimension); // Scalar (1), vector (3), tensor (6), etc.
            println!("    - Number of nodes: {}", type_info.num_nodes); // Total number of nodes for this type
            println!("    - Number of value types: {}", type_info.num_nodevalue_type); // Number of value types
            println!("    - Start index: {}", type_info.start_index); // Where this value type starts in the flattened array
        }
        
        // Show sample node data for verification (helps confirm parsing worked correctly)
        if !node_result.node_values.is_empty() {             // Check if we have any nodes to display
            let show_count = 3;                               // Limit to 3 nodes to avoid overwhelming output
            let total_nodes = node_result.node_values.len();  // Get total number of node values
            
            println!("\nSample node values (first {} entries):", show_count.min(total_nodes)); // Header for sample data
            
            // Loop through first few node values to show sample data
            for i in 0..show_count.min(total_nodes) {         // Show up to show_count nodes or all nodes if fewer
                let node_value = &node_result.node_values[i]; // Get reference to current node value
                println!("  Node {}: values {:?}", node_value.id, node_value.values); // Print node ID and values
            }
        }
    } else {
        println!("NODE VALUES: None");
    }
    
    // Display element values information if available
    if let Some(element_result) = element_data {             // Check if element data exists
        println!("\nELEMENT VALUES:");                      // Subsection header for element data
        println!("Number of value types: {}", element_result.element_value_type_info.len()); // Different types of element values
        
        // Display detailed information about each element value type
        println!("\nElement value types:");                 // Subsection for element value type details
        for type_info in &element_result.element_value_type_info { // Loop through each element value type
            println!("  {:?}:", type_info.elementvalue_type); // Print the element value type name
            println!("    - Dimension: {}", type_info.dimension); // Scalar (1), vector (3), tensor (6), etc.
            println!("    - Number of elements: {}", type_info.num_elements); // Total number of elements for this type
            println!("    - Number of value types: {}", type_info.num_elementvalue_type); // Number of value types
            println!("    - Start index: {}", type_info.start_index); // Where this value type starts in the array
        }
        
        // Show sample element data for verification
        if !element_result.element_values.is_empty() {      // Check if we have any elements to display
            let show_count = 3;                             // Limit to 3 elements to avoid overwhelming output
            let total_elements = element_result.element_values.len(); // Get total number of element values
            
            println!("\nSample element values (first {} entries):", show_count.min(total_elements)); // Header for sample data
            
            // Loop through first few element values to show sample data
            for i in 0..show_count.min(total_elements) {     // Show up to show_count elements or all if fewer
                let element_value = &element_result.element_values[i]; // Get reference to current element value
                println!("  Element {}: values {:?}", element_value.id, element_value.values); // Print element ID and values
            }
        }
    } else {
        println!("\nELEMENT VALUES: None");
    }
}

/// Print detailed summary of parsed mesh data (from .mphtxt or .inp files)
/// Shows mesh statistics, sample node coordinates, and element connectivity information
fn print_mesh_summary(mesh: &MeshData) {
    println!("MESH PARSING RESULTS");                        // Main section header
    println!("{}", "=".repeat(50));                          // Visual separator line
    
    // Print basic mesh statistics to give user overview of parsed data
    println!("Spatial Dimension: {}", mesh.dimension);       // 2D or 3D mesh
    println!("Number of Nodes: {}", mesh.num_nodes);         // Total node count in mesh
    println!("Minimum Node Index: {}", mesh.min_node_index); // Starting index (usually 0 or 1)
    println!("Number of Element Types: {}", mesh.num_eltypes); // How many different element types (triangles, quads, etc.)
    println!("Total Elements: {}", mesh.elements.len());     // Total element count across all types
    
    // Node information section - show sample nodes to verify parsing
    println!("\nNODE INFORMATION");                          // Subsection header
    println!("{}", "=".repeat(30));                          // Subsection separator
    println!("Parsed {} nodes", mesh.nodes.len());           // Confirm total nodes parsed
    
    // Display sample nodes to verify parsing worked correctly
    if !mesh.nodes.is_empty() {                              // Check if we have nodes to display
        let show_count = 5;                                   // Show 5 nodes from beginning and end
        let total_nodes = mesh.nodes.len();                   // Get total number of nodes
    
        // Show first few nodes (or all if we have fewer than show_count)
        let first_end = show_count.min(total_nodes);          // Don't exceed total node count
        println!("\nFirst {} nodes:", first_end);            // Header for first nodes section
        for i in 0..first_end {                              // Loop through first nodes
            let node = &mesh.nodes[i];                       // Get reference to current node
            println!("  Node {}: {:?}", node.id, node.coordinates); // Print node ID and spatial coordinates
        }

        // Show last few nodes if we have many nodes (gives sense of node ID range)
        if total_nodes > show_count * 2 {                    // Only if we have more than 10 total nodes
            let last_start = total_nodes - show_count;       // Calculate starting index for last nodes
            println!("\nLast {} nodes:", show_count);        // Header for last nodes section
            for i in last_start..total_nodes {               // Loop through last nodes
                let node = &mesh.nodes[i];                   // Get reference to current node
                println!("  Node {}: {:?}", node.id, node.coordinates); // Print node ID and coordinates
            }
        } 
    }

    // Element connectivity information section - show how elements connect nodes
    println!("\nELEMENT CONNECTIVITY");                      // Subsection header
    println!("{}", "=".repeat(30));                          // Subsection separator

    // Display elements organized by type (triangles, quads, tetrahedra, etc.)
    for type_info in &mesh.element_type_info {               // Process each element type separately
        let total_elements = type_info.num_elements;          // Get count of elements of this type
        let show_count = 5;                                   // Show 5 elements from beginning and end for each type
    
        // Header showing element type and count
        println!("\n{:?} elements ({} total):", type_info.element_type, total_elements);
    
        // Show first few elements of this type
        let first_end = show_count.min(total_elements);       // Don't exceed total elements of this type
        for i in 0..first_end {                              // Loop through first elements of this type
            // Calculate the actual index in the global elements array
            let elem = &mesh.elements[type_info.start_index + i]; // Get element from global array using offset
            println!("  Element {}: nodes {:?}", elem.id, elem.nodes); // Print element ID and connected node IDs
        }
    
        // Show last few elements if we have many of this type
        if total_elements > show_count * 2 {                 // Only if we have more than 10 elements of this type
            let last_start = total_elements - show_count;    // Calculate starting index for last elements
            println!("  ...");                               // Indicate that we're skipping elements in between
            for i in last_start..total_elements {            // Loop through last elements of this type
                let elem = &mesh.elements[type_info.start_index + i]; // Get element from global array
                println!("  Element {}: nodes {:?}", elem.id, elem.nodes); // Print element ID and connectivity
            }
        }
    }
    
    // Final success message for mesh parsing
    println!("\nPARSING COMPLETED SUCCESSFULLY!");
    println!("{}", "=".repeat(50));                          // Visual separator for clean output
}
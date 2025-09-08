import argparse
import sys
import json
from pathlib import Path
from typing import Optional

from config.optimization_config import ProteinOptimizationConfig
from core.optimizer import ProteinOptimizer
from utils.logging_utils import setup_logging
from analysis.scoring import ScoreAnalyzer
from analysis.energy_calculator import EnergyAnalyzer


def create_argument_parser() -> argparse.ArgumentParser:
   
    # Input files
    parser.add_argument('--input', '-i', type=str,
                      help='Input PDB structure file')
    parser.add_argument('--protocol', '-p', type=str,
                      help='Rosetta protocol XML file')
    parser.add_argument('--config', '-c', type=str,
                      help='Configuration JSON file')
    
    # Optimization parameters
    parser.add_argument('--max-iterations', type=int, default=1000,
                      help='Maximum number of optimization iterations (default: 1000)')
    parser.add_argument('--n-structures', type=int, default=5,
                      help='Number of structures per iteration (default: 5)')
    parser.add_argument('--ph-value', type=float, default=8.0,
                      help='pH value for electrostatic calculations (default: 8.0)')
    
    # Energy goals
    parser.add_argument('--n-terminal-goal', type=float, default=26678.18,
                      help='N-terminal energy goal (default: 26678.18)')
    parser.add_argument('--c-terminal-goal', type=float, default=79814.70,
                      help='C-terminal energy goal (default: 79814.70)')
    parser.add_argument('--initial-score', type=float, default=-658.0,
                      help='Initial Rosetta score (default: -658.0)')
    
    # Tool paths
    parser.add_argument('--rosetta-path', type=str,
                      help='Path to Rosetta executable')
    parser.add_argument('--pdb2pqr-path', type=str,
                      help='Path to PDB2PQR executable')
    parser.add_argument('--apbs-path', type=str,
                      help='Path to APBS executable')
    
    # Logging and output
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                      default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path (default: auto-generated)')
    parser.add_argument('--no-console', action='store_true',
                      help='Disable console output')
    parser.add_argument('--output-dir', type=str,
                      help='Output directory for results')
    
    # Special modes
    parser.add_argument('--resume', type=str,
                      help='Resume optimization from results file')
    parser.add_argument('--additional-iterations', type=int, default=500,
                      help='Additional iterations when resuming (default: 500)')
    parser.add_argument('--analyze', type=str,
                      help='Analyze existing results file')
    parser.add_argument('--validate-only', action='store_true',
                      help='Only validate configuration and tools, do not run')
    
    # Advanced options
    parser.add_argument('--targeted-positions', type=str,
                      help='Comma-separated list of positions for targeted optimization')
    parser.add_argument('--targeted-residues', type=str,
                      help='Comma-separated list of residues for targeted optimization')
    
    return parser


def load_config_from_file(config_file: str) -> ProteinOptimizationConfig:

    with open(config_file, 'r') as f:
        config_data = json.load(f)
    
    return ProteinOptimizationConfig.from_dict(config_data)


def create_config_from_args(args) -> ProteinOptimizationConfig:

    config_kwargs = {}
    
    # Basic parameters
    if args.input:
        config_kwargs['input_pdb'] = args.input
    if args.protocol:
        config_kwargs['protocol_file'] = args.protocol
    
    # Optimization parameters
    config_kwargs['max_iterations'] = args.max_iterations
    config_kwargs['n_structures'] = args.n_structures
    config_kwargs['ph_value'] = args.ph_value
    
    # Energy goals
    config_kwargs['initial_score'] = args.initial_score
    config_kwargs['n_terminal_goal'] = args.n_terminal_goal
    config_kwargs['c_terminal_goal'] = args.c_terminal_goal
    
    # Tool paths
    if args.rosetta_path:
        config_kwargs['rosetta_path'] = args.rosetta_path
    if args.pdb2pqr_path:
        config_kwargs['pdb2pqr_path'] = args.pdb2pqr_path
    if args.apbs_path:
        config_kwargs['apbs_path'] = args.apbs_path
    
    return ProteinOptimizationConfig(**config_kwargs)


def run_optimization(config: ProteinOptimizationConfig, 
                   main_logger, 
                   args) -> None:

    optimizer = ProteinOptimizer(config, main_logger)
    
    if args.targeted_positions and args.targeted_residues:
        # Targeted optimization
        positions = [int(x.strip()) for x in args.targeted_positions.split(',')]
        residues = [x.strip() for x in args.targeted_residues.split(',')]
        
        results = optimizer.run_targeted_optimization(
            positions, residues, config.max_iterations
        )
    else:
        # Standard optimization
        results = optimizer.run_optimization()
    
    # Print summary
    print_optimization_summary(results, optimizer.get_current_state())


def resume_optimization(resume_file: str, 
                      additional_iterations: int,
                      main_logger) -> None:

    # Load configuration from the same directory
    results_dir = Path(resume_file).parent
    config_file = results_dir / "config.json"
    
    if not config_file.exists():
        print(f"Error: Configuration file not found: {config_file}")
        sys.exit(1)
    
    config = load_config_from_file(str(config_file))
    optimizer = ProteinOptimizer(config, main_logger)
    
    results = optimizer.resume_optimization(resume_file, additional_iterations)
    
    print_optimization_summary(results, optimizer.get_current_state())


def analyze_results(results_file: str) -> None:

    from utils.file_utils import FileManager
    
    file_manager = FileManager()
    
    try:
        results = file_manager.load_optimization_results(results_file)
    except Exception as e:
        print(f"Error loading results: {e}")
        sys.exit(1)
    
    if not results:
        print("No results found in file.")
        return
    
    # Load configuration if available
    results_dir = Path(results_file).parent
    config_file = results_dir / "config.json"
    
    if config_file.exists():
        config = load_config_from_file(str(config_file))
    else:
        print("Warning: Configuration file not found, using defaults for analysis")
        config = ProteinOptimizationConfig()
    
    # Perform analysis
    score_analyzer = ScoreAnalyzer()
    energy_analyzer = EnergyAnalyzer(config)
    
    print("\n" + "="*60)
    print("OPTIMIZATION RESULTS ANALYSIS")
    print("="*60)
    
    # Score analysis
    print(score_analyzer.generate_score_report(results))
    print("\n")
    
    # Energy analysis
    print(energy_analyzer.generate_energy_report(results))
    
    # Find Pareto optimal solutions
    pareto_solutions = energy_analyzer.find_pareto_optimal_solutions(results)
    if pareto_solutions:
        print(f"\nPARETTO OPTIMAL SOLUTIONS ({len(pareto_solutions)} found):")
        print("-" * 40)
        for i, sol in enumerate(pareto_solutions[:10], 1):
            print(f"{i:2d}. Score: {sol.score:7.2f}, N-diff: {sol.n_diff:7.1f}, C-diff: {sol.c_diff:7.1f}")


def print_optimization_summary(results, current_state) -> None:

    if not results:
        print("No optimization results to summarize.")
        return
    
    accepted_count = sum(1 for r in results if r.accepted)
    acceptance_rate = accepted_count / len(results) * 100
    
    best_score = min(r.score for r in results)
    final_score = current_state['current_score']
    score_improvement = results[0].score - final_score
    
    print("\n" + "="*50)
    print("OPTIMIZATION SUMMARY")
    print("="*50)
    print(f"Total iterations:     {len(results)}")
    print(f"Accepted mutations:   {accepted_count}")
    print(f"Acceptance rate:      {acceptance_rate:.1f}%")
    print(f"Initial score:        {results[0].score:.2f}")
    print(f"Final score:          {final_score:.2f}")
    print(f"Best score:           {best_score:.2f}")
    print(f"Score improvement:    {score_improvement:.2f}")
    print("="*50)


def validate_configuration(config: ProteinOptimizationConfig) -> bool:

    print("Validating configuration and environment...")
    
    # Check configuration
    errors = config.validate_paths()
    if errors:
        print("Configuration errors:")
        for error in errors:
            print(f"  - {error}")
        return False
    
    # Check tools
    from external_tools.rosetta import RosettaTool
    from external_tools.apbs import ElectrostaticCalculator
    
    rosetta_tool = RosettaTool(config)
    electrostatic_calc = ElectrostaticCalculator(config)
    
    tools_status = {
        'Rosetta': rosetta_tool.check_availability(),
        'PDB2PQR': electrostatic_calc.pdb2pqr.check_availability(),
        'APBS': electrostatic_calc.apbs.check_availability()
    }
    
    print("Tool availability:")
    all_available = True
    for tool, available in tools_status.items():
        status = "✓" if available else "✗"
        print(f"  {status} {tool}")
        if not available:
            all_available = False
    
    if all_available:
        print("✓ All tools are available")
    else:
        print("✗ Some tools are missing")
    
    return all_available


def main():
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Setup logging
    main_logger, performance_logger, structured_logger = setup_logging(
        log_level=args.log_level,
        log_file=args.log_file,
        console_output=not args.no_console
    )
    
    try:
        # Handle special modes
        if args.analyze:
            analyze_results(args.analyze)
            return
        
        if args.resume:
            resume_optimization(args.resume, args.additional_iterations, main_logger)
            return
        
        # Load or create configuration
        if args.config:
            config = load_config_from_file(args.config)
            print(f"Loaded configuration from: {args.config}")
        else:
            config = create_config_from_args(args)
        
        # Validate configuration
        if not validate_configuration(config):
            print("Configuration validation failed. Please check your setup.")
            sys.exit(1)
        
        if args.validate_only:
            print("Configuration validation completed successfully.")
            return
        
        # Ensure required arguments for optimization
        if not args.config and (not args.input or not args.protocol):
            print("Error: --input and --protocol are required (or use --config)")
            sys.exit(1)
        
        # Run optimization
        run_optimization(config, main_logger, args)
        
    except KeyboardInterrupt:
        print("\nOptimization interrupted by user.")
        sys.exit(1)
    except Exception as e:
        main_logger.log_error(e, "main execution")
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

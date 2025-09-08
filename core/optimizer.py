import os
import time
from typing import List, Optional
from pathlib import Path

from config.optimization_config import ProteinOptimizationConfig
from core.mutation import MutationResult, MutationGenerator, MutationHistory
from core.acceptance import AcceptanceCriteria
from external_tools.rosetta import RosettaTool
from external_tools.apbs import ElectrostaticCalculator
from external_tools.base import ToolChain
from utils.file_utils import FileManager, DataExporter
from utils.logging_utils import OptimizationLogger, PerformanceLogger, StructuredLogger


class ProteinOptimizer:
   
    def __init__(self, 
                 config: ProteinOptimizationConfig,
                 logger: Optional[OptimizationLogger] = None):

        self.config = config
        
        # Setup logging
        if logger is None:
            self.logger = OptimizationLogger(log_file=f"optimization_{int(time.time())}.log")
        else:
            self.logger = logger
        
        self.performance_logger = PerformanceLogger(self.logger.get_logger())
        
        # Initialize current state
        self.current_score = config.initial_score
        self.current_n_energy = 28487  # Initial N-terminal energy
        self.current_c_energy = 73823  # Initial C-terminal energy
        
        # Initialize components
        self.mutation_generator = MutationGenerator(config)
        self.acceptance_criteria = AcceptanceCriteria(config)
        self.mutation_history = MutationHistory()
        
        # Initialize tools
        self.rosetta_tool = RosettaTool(config)
        self.electrostatic_calculator = ElectrostaticCalculator(config)
        
        # Initialize file management
        self.file_manager = FileManager()
        self.data_exporter = DataExporter()
        
        # Validate environment
        self._validate_environment()
    
    def _validate_environment(self) -> None:
        """Validate that all required tools and files are available."""
        self.logger.log_optimization_start(self.config.to_dict())
        
        # Validate configuration
        config_errors = self.config.validate_paths()
        if config_errors:
            for error in config_errors:
                self.logger.log_warning(error)
        
        # Validate tools
        tool_status = self.electrostatic_calculator.validate_tools()
        self.logger.log_tool_status("Rosetta", self.rosetta_tool.check_availability(), self.config.rosetta_path)
        self.logger.log_tool_status("PDB2PQR", tool_status['pdb2pqr_available'], self.config.pdb2pqr_path)
        self.logger.log_tool_status("APBS", tool_status['apbs_available'], self.config.apbs_path)
        
        # Check if all tools are available
        all_tools_available = (
            self.rosetta_tool.check_availability() and
            tool_status['pdb2pqr_available'] and
            tool_status['apbs_available']
        )
        
        if not all_tools_available:
            self.logger.log_warning("Some tools are not available. Optimization may fail.")
    
    def run_optimization(self) -> List[MutationResult]:

        self.logger.get_logger().info(f"Starting optimization for {self.config.max_iterations} iterations")
        self.performance_logger.start_timer("total_optimization")
        
        for iteration in range(self.config.max_iterations):
            try:
                self.logger.log_iteration_start(iteration, self.config.max_iterations)
                self.performance_logger.start_timer(f"iteration_{iteration}")
                
                # Run single iteration
                result = self._run_single_iteration(iteration)
                
                # Add to history
                self.mutation_history.add_result(result)
                
                # Update current state if accepted
                if result.accepted:
                    self._accept_mutation(result, iteration)
                else:
                    self._reject_mutation(result, iteration)
                
                # Log iteration completion
                iteration_time = self.performance_logger.stop_timer(f"iteration_{iteration}")
                self.performance_logger.log_iteration_timing(iteration, iteration_time)
                
                # Cleanup temporary files
                self.file_manager.cleanup_temporary_files()
                
            except Exception as e:
                self.logger.log_error(e, f"iteration {iteration}")
                self.file_manager.cleanup_temporary_files()
                continue
        
        # Complete optimization
        total_time = self.performance_logger.stop_timer("total_optimization")
        self._finalize_optimization(total_time)
        
        return self.mutation_history.results
    
    def _run_single_iteration(self, iteration: int) -> MutationResult:

        # Create mutation
        self.performance_logger.start_timer("mutation_generation")
        position, residue = self.mutation_generator.create_mutation_resfile(iteration)
        self.performance_logger.stop_timer("mutation_generation")
        
        self.logger.log_mutation_details(position, residue, "random")
        
        # Run Rosetta relaxation
        self.performance_logger.start_timer("rosetta_execution")
        self.rosetta_tool.run_relaxation()
        rosetta_time = self.performance_logger.stop_timer("rosetta_execution")
        self.performance_logger.log_tool_timing("Rosetta", rosetta_time)
        
        # Analyze Rosetta results
        self.performance_logger.start_timer("score_analysis")
        score, best_structure = self.rosetta_tool.analyze_scores()
        self.performance_logger.stop_timer("score_analysis")
        
        # Calculate electrostatic energies
        self.performance_logger.start_timer("electrostatic_calculation")
        n_energy, c_energy = self.electrostatic_calculator.calculate_terminal_energies(best_structure)
        electrostatic_time = self.performance_logger.stop_timer("electrostatic_calculation")
        self.performance_logger.log_tool_timing("APBS/PDB2PQR", electrostatic_time)
        
        self.logger.log_scores(score, n_energy, c_energy)
        
        # Calculate energy differences from goals
        n_diff = abs(n_energy - self.config.n_terminal_goal)
        c_diff = abs(c_energy - self.config.c_terminal_goal)
        
        # Determine acceptance
        self.performance_logger.start_timer("acceptance_calculation")
        probability, mutation_type = self.acceptance_criteria.calculate_acceptance_probability(
            score, n_energy, c_energy,
            self.current_score, self.current_n_energy, self.current_c_energy
        )
        accepted = self.acceptance_criteria.should_accept_mutation(probability)
        self.performance_logger.stop_timer("acceptance_calculation")
        
        self.logger.log_acceptance_decision(accepted, probability, mutation_type)
        
        # Create result object
        result = MutationResult(
            iteration=iteration,
            score=score,
            n_terminal_energy=n_energy,
            c_terminal_energy=c_energy,
            n_diff=n_diff,
            c_diff=c_diff,
            accepted=accepted,
            acceptance_probability=probability,
            mutation_type=mutation_type,
            best_structure=best_structure,
            position=position,
            new_residue=residue
        )
        
        return result
    
    def _accept_mutation(self, result: MutationResult, iteration: int) -> None:

        # Update current state
        self.current_score = result.score
        self.current_n_energy = result.n_terminal_energy
        self.current_c_energy = result.c_terminal_energy
        
        # Save accepted structure
        output_name = self.file_manager.save_accepted_structure(
            result.best_structure, iteration, result.mutation_type
        )
        
        # Update input structure for next iteration
        os.rename("output_chain.pdb", self.config.input_pdb)
        
        # Remove the temporary best structure
        if Path(result.best_structure).exists():
            os.remove(result.best_structure)
        
        self.logger.get_logger().info(f"Mutation accepted: {result.mutation_type}")
    
    def _reject_mutation(self, result: MutationResult, iteration: int) -> None:

        # Save rejected structure for analysis (optional)
        self.file_manager.save_rejected_structure(
            result.best_structure, iteration, result.mutation_type
        )
        
        # Clean up temporary files
        temp_files = ["output_chain.pdb", result.best_structure]
        for temp_file in temp_files:
            if Path(temp_file).exists():
                os.remove(temp_file)
        
        self.logger.get_logger().debug(f"Mutation rejected: {result.mutation_type}")
    
    def _finalize_optimization(self, total_time: float) -> None:

        results = self.mutation_history.results
        accepted_count = len(self.mutation_history.get_accepted_mutations())
        
        self.logger.log_optimization_complete(
            self.current_score, accepted_count, len(results)
        )
        
        self.logger.get_logger().info(f"Total optimization time: {total_time:.1f} seconds")
        
        # Save results
        self._save_results()
    
    def _save_results(self) -> None:

        # Create results directory
        results_dir = self.file_manager.create_results_directory()
        
        # Archive the optimization run
        self.file_manager.archive_optimization_run(
            self.mutation_history.results, self.config, results_dir
        )
        
        # Export data in various formats
        csv_file = results_dir / "results.csv"
        self.data_exporter.export_to_csv(self.mutation_history.results, str(csv_file))
        
        summary_file = results_dir / "summary_stats.json"
        self.data_exporter.export_summary_stats(self.mutation_history, str(summary_file))
        
        trajectory_file = results_dir / "trajectory_data.json"
        self.data_exporter.export_trajectory_data(self.mutation_history.results, str(trajectory_file))
        
        self.logger.get_logger().info(f"Results saved to: {results_dir}")
    
    def get_current_state(self) -> dict:

        return {
            'current_score': self.current_score,
            'current_n_energy': self.current_n_energy,
            'current_c_energy': self.current_c_energy,
            'total_iterations': len(self.mutation_history.results),
            'accepted_mutations': len(self.mutation_history.get_accepted_mutations()),
            'acceptance_rate': self.mutation_history.get_acceptance_rate(),
            'best_score': self.mutation_history.get_best_score()
        }
    
    def resume_optimization(self, 
                          checkpoint_file: str, 
                          additional_iterations: int) -> List[MutationResult]:

        # Load previous results
        previous_results = self.file_manager.load_optimization_results(checkpoint_file)
        
        # Restore state from last accepted mutation
        accepted_results = [r for r in previous_results if r.accepted]
        if accepted_results:
            last_accepted = accepted_results[-1]
            self.current_score = last_accepted.score
            self.current_n_energy = last_accepted.n_terminal_energy
            self.current_c_energy = last_accepted.c_terminal_energy
        
        # Add previous results to history
        for result in previous_results:
            self.mutation_history.add_result(result)
        
        self.logger.get_logger().info(f"Resuming optimization from iteration {len(previous_results)}")
        
        # Update configuration for additional iterations
        original_max = self.config.max_iterations
        self.config.max_iterations = additional_iterations
        
        # Run additional iterations
        new_results = self.run_optimization()
        
        # Restore original max iterations
        self.config.max_iterations = original_max
        
        return self.mutation_history.results
    
    def run_targeted_optimization(self, 
                                target_positions: List[int],
                                target_residues: List[str],
                                max_iterations: int) -> List[MutationResult]:
        self.logger.get_logger().info(f"Starting targeted optimization: {len(target_positions)} positions")
        
        original_max = self.config.max_iterations
        self.config.max_iterations = max_iterations
        
        # Override mutation generator for targeted mutations
        original_generator = self.mutation_generator
        
        try:
            for iteration in range(max_iterations):
                # Select position and residue from targets
                position = target_positions[iteration % len(target_positions)]
                residue = target_residues[iteration % len(target_residues)]
                
                # Create targeted mutation
                self.mutation_generator.create_targeted_mutation(position, residue)
                
                # Run the rest of the iteration normally
                result = self._run_single_iteration(iteration)
                self.mutation_history.add_result(result)
                
                if result.accepted:
                    self._accept_mutation(result, iteration)
                else:
                    self._reject_mutation(result, iteration)
                
                self.file_manager.cleanup_temporary_files()
        
        finally:
            # Restore original settings
            self.config.max_iterations = original_max
            self.mutation_generator = original_generator
        
        return self.mutation_history.results
    
    def analyze_optimization_progress(self) -> dict:

        from analysis.scoring import ScoreAnalyzer
        from analysis.energy_calculator import EnergyAnalyzer
        
        if not self.mutation_history.results:
            return {"status": "no_data"}
        
        score_analyzer = ScoreAnalyzer()
        energy_analyzer = EnergyAnalyzer(self.config)
        
        analysis = {
            "current_state": self.get_current_state(),
            "trajectory_analysis": score_analyzer.analyze_score_trajectory(self.mutation_history.results),
            "energy_analysis": energy_analyzer.calculate_energy_deviations(self.mutation_history.results),
            "acceptance_patterns": score_analyzer.analyze_acceptance_patterns(self.mutation_history.results),
            "optimization_phases": score_analyzer.detect_optimization_phases(self.mutation_history.results)
        }
        
        return analysis

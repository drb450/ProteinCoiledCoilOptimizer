"""
Coiled-coil design. Trimodal MC search of coiled-coil proteins using criteria for rosetta score, N-terminal EEbcf, and C-terminal EEbcf

Citation: "Computational prediction of coiled-coil protein gelation dynamics and structure" 
by Dustin Britton.

Version: 3.0
"""

import os
import sys
import glob
import random
import math
import time
import logging
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
from shutil import copyfile
import json
from datetime import datetime


class ProteinOptimizationConfig:
    # Configuration parameters for protein optimization
    # Initialize configuration with default values
    # Important: change paths for rosetta, pdb2pqr, and apbs
    def __init__(self, 
                 initial_score: float = -658.0,
                 n_terminal_goal: float = 26678.18,
                 c_terminal_goal: float = 79814.70,
                 initial_n_diff: float = None,
                 initial_c_diff: float = None,
                 a_d_residues: List[int] = None,
                 e_g_residues: List[int] = None,
                 favorable_residues: List[str] = None,
                 secondary_residues: List[str] = None,
                 zeta: float = 3.93e-5,
                 alpha: float = 1.96e-4,
                 beta: float = 1.31e-4,
                 max_iterations: int = 1000,
                 n_structures: int = 5,
                 ph_value: float = 8.0,
                 rosetta_path: str = "/share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease",
                 pdb2pqr_path: str = "/share/apps/pdb2pqr/3.1.0/bin/pdb2pqr30",
                 apbs_apth: str = "/share/apps/apbs",
                 protocol_file: str = "symandrelax.xml",
                 input_pdb: str = "input.pdb"):
        
        # Energy targets
        self.initial_score = initial_score
        self.n_terminal_goal = n_terminal_goal
        self.c_terminal_goal = c_terminal_goal
        
        # Initial energy differences
        self.initial_n_diff = initial_n_diff if initial_n_diff is not None else abs(28487 - n_terminal_goal)
        self.initial_c_diff = initial_c_diff if initial_c_diff is not None else abs(73823 - c_terminal_goal)
        
        # Residue definitions
        self.a_d_residues = a_d_residues if a_d_residues is not None else [20, 23, 27, 30, 34, 37, 41, 44, 51]
        self.e_g_residues = e_g_residues if e_g_residues is not None else [24, 26, 31, 33, 38, 40, 45, 47, 52, 54]
        self.favorable_residues = favorable_residues if favorable_residues is not None else ['V', 'I', 'M', 'T', 'Q', 'L']
        self.secondary_residues = secondary_residues if secondary_residues is not None else ['A', 'E', 'K', 'Q', 'N', 'T', 'D', 'S']
        
        # Scaling factors for acceptance probability
        self.zeta = zeta
        self.alpha = alpha
        self.beta = beta
        
        # Simulation parameters
        self.max_iterations = max_iterations
        self.n_structures = n_structures
        self.ph_value = ph_value
        
        # File paths
        self.rosetta_path = rosetta_path
        self.pdb2pqr_path = pdb2pqr_path
        self.protocol_file = protocol_file
        self.apbs_path = apbs_path
        self.input_pdb = input_pdb
    
    def __repr__(self):
        # String representation for debugging
        return (f"ProteinOptimizationConfig(initial_score={self.initial_score}, "
                f"max_iterations={self.max_iterations}, "
                f"n_terminal_goal={self.n_terminal_goal})")
    
    def __eq__(self, other):
        # Compare two configuration objects
        if not isinstance(other, ProteinOptimizationConfig):
            return False
        return (self.initial_score == other.initial_score and
                self.n_terminal_goal == other.n_terminal_goal and
                self.c_terminal_goal == other.c_terminal_goal and
                self.max_iterations == other.max_iterations)


class MutationResult:
    # Results from a single mutation iteration
    
    def __init__(self, iteration: int, score: float, n_terminal_energy: float, 
                 c_terminal_energy: float, n_diff: float, c_diff: float,
                 accepted: bool, acceptance_probability: float, 
                 mutation_type: str, best_structure: str):
        # Initialize mutation result
        self.iteration = iteration
        self.score = score
        self.n_terminal_energy = n_terminal_energy
        self.c_terminal_energy = c_terminal_energy
        self.n_diff = n_diff
        self.c_diff = c_diff
        self.accepted = accepted
        self.acceptance_probability = acceptance_probability
        self.mutation_type = mutation_type
        self.best_structure = best_structure
    
    def __repr__(self):
        # String representation for debugging
        return (f"MutationResult(iteration={self.iteration}, "
                f"score={self.score:.2f}, accepted={self.accepted})")
    
    def to_dict(self) -> Dict[str, Any]:
        # Convert to dictionary for JSON serialization
        return {
            'iteration': self.iteration,
            'score': self.score,
            'n_terminal_energy': self.n_terminal_energy,
            'c_terminal_energy': self.c_terminal_energy,
            'n_diff': self.n_diff,
            'c_diff': self.c_diff,
            'accepted': self.accepted,
            'acceptance_probability': self.acceptance_probability,
            'mutation_type': self.mutation_type,
            'best_structure': self.best_structure
        }


class ProteinOptimizer:
    # Monte Carlo optimizer for coiled-coil proteins.
    
    def __init__(self, config: ProteinOptimizationConfig):
        # Initialize the optimizer with configuration.
        self.config = config
        self.current_score = config.initial_score
        self.current_n_energy = 28487  # Initial N-terminal energy
        self.current_c_energy = 73823  # Initial C-terminal energy
        
        self.results: List[MutationResult] = []
        self.index_range = list(range(100))
        
        # Setup logging
        self._setup_logging()
        
        # Validate environment
        self._validate_environment()
    
    def _setup_logging(self) -> None:
        # Configure logging for the optimization process
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler(f'protein_optimization_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def _validate_environment(self) -> None:
        # Validate tools
        required_files = [
            self.config.protocol_file,
            self.config.input_pdb
        ]
        
        required_executables = [
            self.config.rosetta_path,
            self.config.pdb2pqr_path,
            self.config.apbs_path
        ]
        
        for file_path in required_files:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"Required file not found: {file_path}")
        
        for executable in required_executables:
            if not Path(executable).exists() and not self._is_executable_in_path(executable):
                self.logger.warning(f"Executable not found in expected location: {executable}")
    
    def _is_executable_in_path(self, executable: str) -> bool:
        # Check exe path
        try:
            subprocess.run(['which', executable], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError:
            return False
    
    def create_mutation_resfile(self, iteration: int) -> Tuple[int, str]:
        # Create a resfile for the current mutation.
        list_choice = random.choice([0, 1])
        residue_positions = (self.config.a_d_residues if list_choice == 0 
                           else self.config.e_g_residues)
        residue_options = (self.config.favorable_residues if list_choice == 0 
                         else self.config.secondary_residues)
        
        target_residue = random.choice(residue_positions)
        new_residue = random.choice(residue_options)
        
        resfile_content = f"NATAA \nstart \n\n{target_residue} A PIKAA {new_residue}\n"
        
        with open("mut.resfile", 'w') as resfile:
            resfile.write(resfile_content)
        
        self.logger.debug(f"Mutation {iteration}: Position {target_residue} -> {new_residue}")
        return target_residue, new_residue
    
    def run_rosetta_relaxation(self) -> None:
        # Execute Rosetta protein relaxation.
        cmd = [
            self.config.rosetta_path,
            '-parser:protocol', self.config.protocol_file,
            '-s', self.config.input_pdb,
            '-nstruct', str(self.config.n_structures),
            '-out:prefix', 'temp_'
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.debug("Rosetta relaxation completed successfully")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Rosetta execution failed: {e}")
            raise
        
        time.sleep(2)
    
    def analyze_scores(self) -> Tuple[float, str]:
        # Analyze scores from Rosetta output and return best score and structure
        try:
            with open('temp_score.sc', 'r') as score_file:
                lines = score_file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError("Score file not found. Rosetta may have failed.")
        
        if len(lines) < 3:
            raise ValueError("Insufficient score data in output file")
        
        scores = []
        for line in lines[2:]:  # Skip header lines
            try:
                score = float(line.split()[1])
                scores.append(score)
            except (IndexError, ValueError) as e:
                self.logger.warning(f"Could not parse score from line: {line.strip()}")
                continue
        
        if not scores:
            raise ValueError("No valid scores found in output file")
        
        best_score = min(scores)
        
        # Find the corresponding structure file
        best_structure = ""
        for line in lines[2:]:
            try:
                if float(line.split()[1]) == best_score:
                    best_structure = line.split()[21] + '.pdb'
                    break
            except (IndexError, ValueError):
                continue
        
        if not best_structure:
            raise ValueError("Could not identify best scoring structure")
        
        self.logger.debug(f"Best score: {best_score}, Structure: {best_structure}")
        return best_score, best_structure
    
    def prepare_structure_for_apbs(self, structure_file: str) -> None:
        # Prepare structure for APBS calculation
        # Extract single chain
        try:
            with open(structure_file, 'r') as full_protein, \
                 open('output_chain.pdb', 'w') as chain_protein:
                for line in full_protein:
                    if 'TER' not in line:
                        chain_protein.write(line)
                    else:
                        break
        except FileNotFoundError:
            raise FileNotFoundError(f"Structure file not found: {structure_file}")
        
        # Run PDB2PQR
        pdb2pqr_cmd = [
            self.config.pdb2pqr_path,
            f'--apbs-input=output_pqr.in',
            f'--with-ph={self.config.ph_value}',
            '--titration-state-method=propka',
            '--ff=AMBER',
            structure_file,
            'output_pqr'
        ]
        
        try:
            subprocess.run(pdb2pqr_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError:
            # Try with --noopts flag if first attempt fails
            self.logger.warning("PDB2PQR failed, retrying with --noopts")
            pdb2pqr_cmd.insert(-2, '--noopts')
            subprocess.run(pdb2pqr_cmd, check=True, capture_output=True, text=True)
        
        # Handle potential PDB2PQR errors
        self._handle_pdb2pqr_errors(structure_file)
    
    def _handle_pdb2pqr_errors(self, structure_file: str) -> None:
        # Handle common PDB2PQR errors.
        for err_file in glob.glob("*.err"):
            try:
                with open(err_file, 'r') as check_err:
                    if "TypeError: '<' not supported between instances of" in check_err.read():
                        self.logger.warning("PDB2PQR TypeError detected, rerunning with --noopts")
                        subprocess.run([
                            self.config.pdb2pqr_path,
                            f'--apbs-input=output_pqr.in',
                            f'--with-ph={self.config.ph_value}',
                            '--titration-state-method=propka',
                            '--ff=AMBER',
                            '--noopts',
                            structure_file,
                            'output_pqr'
                        ], check=True)
                        break
            except Exception as e:
                self.logger.warning(f"Error handling PDB2PQR error file {err_file}: {e}")
    
    def calculate_electrostatic_energies(self) -> Tuple[float, float]:
        # Calculate electrostatic energies using APBS
        # Modify APBS input file
        try:
            with open('output_pqr.in', 'r') as input_file, \
                 open('output.in', 'w') as output_file:
                for line in input_file:
                    output_file.write(line.replace('total', 'comps'))
        except FileNotFoundError:
            raise FileNotFoundError("APBS input file not found")
        
        # Run APBS
        try:
            subprocess.run([self.config.apbs_path, '--output-file=output_apbs', 'output.in'], 
                         check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"APBS execution failed: {e}")
            raise
        
        # Parse APBS output
        try:
            with open('output_apbs', 'r') as apbs_file:
                lines = apbs_file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError("APBS output file not found")
        
        elec_energies = []
        for line in lines:
            if line.startswith('atom'):
                try:
                    energy = float(line.split()[2])
                    elec_energies.append(energy)
                except (IndexError, ValueError):
                    continue
        
        if not elec_energies:
            raise ValueError("No electrostatic energies found in APBS output")
        
        # Split energies
        mid_point = len(elec_energies) // 2
        polar_energies = elec_energies[mid_point:]
        
        n_terminal_energy = sum(polar_energies[:5])
        c_terminal_energy = sum(polar_energies[-5:])
        
        return n_terminal_energy, c_terminal_energy
    
    def calculate_acceptance_probability(self, score: float, n_diff: float, c_diff: float) -> Tuple[float, str]:
        # Calculate acceptance probability based on Monte Carlo criteria
        mutation_type = ""
        probability = 0.0
        
        score_improved = score <= self.current_score
        n_improved = n_diff <= abs(self.current_n_energy - self.config.n_terminal_goal)
        c_improved = c_diff <= abs(self.current_c_energy - self.config.c_terminal_goal)
        
        current_n_diff = abs(self.current_n_energy - self.config.n_terminal_goal)
        current_c_diff = abs(self.current_c_energy - self.config.c_terminal_goal)
        
        if score_improved and n_improved and c_improved:
            # Accept unconditionally
            probability = 1.0
            mutation_type = "accept"
        
        elif not score_improved and n_improved and c_improved:
            # R condition - score degraded
            delta_e = score - self.current_score
            exp_factor = delta_e / (self.current_score * self.config.zeta)
            probability = self._safe_exp(exp_factor)
            mutation_type = "score_penalty"
        
        elif score_improved and not n_improved and c_improved:
            # P condition - N-terminal degraded
            delta_e = n_diff - current_n_diff
            exp_factor = -delta_e / (current_n_diff * self.config.beta)
            probability = self._safe_exp(exp_factor)
            mutation_type = "n_terminal_penalty"
        
        elif score_improved and n_improved and not c_improved:
            # N condition - C-terminal degraded
            delta_e = c_diff - current_c_diff
            exp_factor = -delta_e / (current_c_diff * self.config.alpha)
            probability = self._safe_exp(exp_factor)
            mutation_type = "c_terminal_penalty"
        
        else:
            # Multiple penalties - combine probabilities
            prob_components = []
            
            if not score_improved:
                delta_e = score - self.current_score
                exp_factor = delta_e / (self.current_score * self.config.zeta)
                prob_components.append(self._safe_exp(exp_factor))
                mutation_type += "R"
            
            if not n_improved:
                delta_e = n_diff - current_n_diff
                exp_factor = -delta_e / (current_n_diff * self.config.beta)
                prob_components.append(self._safe_exp(exp_factor))
                mutation_type += "P"
            
            if not c_improved:
                delta_e = c_diff - current_c_diff
                exp_factor = -delta_e / (current_c_diff * self.config.alpha)
                prob_components.append(self._safe_exp(exp_factor))
                mutation_type += "N"
            
            probability = 1.0 if not prob_components else math.prod(prob_components)
            mutation_type = f"multi_penalty_{mutation_type}"
        
        return probability, mutation_type
    
    def _safe_exp(self, x: float) -> float:
        # Safely calculate exponential, handling overflow
        try:
            return math.exp(x)
        except OverflowError:
            return 0.0
    
    def should_accept_mutation(self, probability: float) -> bool:
        # Determine whether to accept mutation based on probability
        if probability >= 1.0:
            return True
        
        prob_percentage = min(int(round(probability * 100, 0)), 100)
        accept_indices = list(range(prob_percentage))
        random_index = random.choice(self.index_range)
        
        return random_index in accept_indices
    
    def cleanup_temporary_files(self) -> None:
        # Clean up temporary files from the current iteration
        patterns = ["output*", "temp_*", "*.err"]
        
        for pattern in patterns:
            for file_path in glob.glob(pattern):
                try:
                    os.remove(file_path)
                except OSError as e:
                    self.logger.warning(f"Could not remove {file_path}: {e}")
    
    def save_results(self, output_dir: str = "results") -> None:
        # Save optimization results to files
        Path(output_dir).mkdir(exist_ok=True)
        
        # Save results as JSON
        results_data = {
            "config": {
                "initial_score": self.config.initial_score,
                "n_terminal_goal": self.config.n_terminal_goal,
                "c_terminal_goal": self.config.c_terminal_goal,
                "max_iterations": self.config.max_iterations
            },
            "results": [r.to_dict() for r in self.results]
        }
        
        with open(f"{output_dir}/optimization_results.json", 'w') as f:
            json.dump(results_data, f, indent=2)
        
        # Save summary statistics
        accepted_results = [r for r in self.results if r.accepted]
        
        summary = {
            "total_iterations": len(self.results),
            "accepted_mutations": len(accepted_results),
            "acceptance_rate": len(accepted_results) / len(self.results) if self.results else 0,
            "final_score": self.current_score,
            "best_score": min((r.score for r in self.results), default=self.config.initial_score),
            "final_n_energy": self.current_n_energy,
            "final_c_energy": self.current_c_energy
        }
        
        with open(f"{output_dir}/optimization_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.logger.info(f"Results saved to {output_dir}/")
    
    def run_optimization(self) -> List[MutationResult]:
        # Run the complete Monte Carlo optimization
        self.logger.info(f"Starting protein optimization for {self.config.max_iterations} iterations")
        
        for iteration in range(self.config.max_iterations):
            try:
                self.logger.info(f"Iteration {iteration + 1}/{self.config.max_iterations}")
                
                # Create mutation
                position, residue = self.create_mutation_resfile(iteration)
                
                # Run Rosetta
                self.run_rosetta_relaxation()
                
                # Analyze results
                score, best_structure = self.analyze_scores()
                
                # Prepare for APBS
                self.prepare_structure_for_apbs(best_structure)
                
                # Calculate energies
                n_energy, c_energy = self.calculate_electrostatic_energies()
                
                # Calculate differences from goals
                n_diff = abs(n_energy - self.config.n_terminal_goal)
                c_diff = abs(c_energy - self.config.c_terminal_goal)
                
                # Determine acceptance
                probability, mutation_type = self.calculate_acceptance_probability(score, n_diff, c_diff)
                accepted = self.should_accept_mutation(probability)
                
                # Create result record
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
                    best_structure=best_structure
                )
                
                self.results.append(result)
                
                # Update current state if accepted
                if accepted:
                    self.current_score = score
                    self.current_n_energy = n_energy
                    self.current_c_energy = c_energy
                    
                    # Save accepted structure
                    output_name = f"{iteration}_{mutation_type}_output.pdb"
                    copyfile(best_structure, output_name)
                    
                    # Update input for next iteration
                    os.rename("output_chain.pdb", self.config.input_pdb)
                    os.remove(best_structure)
                    
                    self.logger.info(f"Mutation accepted: {mutation_type}, probability: {probability:.4f}")
                else:
                    # Save rejected structure for analysis
                    rejected_name = f"{iteration}_rejected_output.pdb"
                    copyfile(best_structure, rejected_name)
                    os.remove("output_chain.pdb")
                    
                    self.logger.info(f"Mutation rejected: {mutation_type}, probability: {probability:.4f}")
                
                # Cleanup
                self.cleanup_temporary_files()
                
            except Exception as e:
                self.logger.error(f"Error in iteration {iteration}: {e}")
                self.cleanup_temporary_files()
                continue
        
        self.logger.info("Optimization completed")
        self.save_results()
        
        return self.results


def main():
    # Main execution function
    # Load configuration (could be from file in production)
    config = ProteinOptimizationConfig()
    
    # Create optimizer
    optimizer = ProteinOptimizer(config)
    
    # Run optimization
    results = optimizer.run_optimization()
    
    # Print summary
    accepted_count = sum(1 for r in results if r.accepted)
    print(f"\nOptimization Summary:")
    print(f"Total iterations: {len(results)}")
    print(f"Accepted mutations: {accepted_count}")
    print(f"Acceptance rate: {accepted_count/len(results)*100:.1f}%")
    print(f"Final score: {optimizer.current_score:.2f}")
    print(f"Best score achieved: {min(r.score for r in results):.2f}")


if __name__ == "__main__":
    main()


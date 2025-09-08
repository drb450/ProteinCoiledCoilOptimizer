import os
import glob
import json
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
import logging

from core.mutation import MutationResult, MutationHistory
from config.optimization_config import ProteinOptimizationConfig


class FileManager:
    
    def __init__(self, working_directory: str = "."):
        self.logger.info(f"Archived optimization run to: {results_dir}")
        return results_dir
    
    def load_optimization_results(self, results_file: str) -> List[MutationResult]:
        with open(results_file, 'r') as f:
            results_data = json.load(f)
        
        results = []
        for data in results_data:
            result = MutationResult(
                iteration=data['iteration'],
                score=data['score'],
                n_terminal_energy=data['n_terminal_energy'],
                c_terminal_energy=data['c_terminal_energy'],
                n_diff=data['n_diff'],
                c_diff=data['c_diff'],
                accepted=data['accepted'],
                acceptance_probability=data['acceptance_probability'],
                mutation_type=data['mutation_type'],
                best_structure=data['best_structure'],
                position=data.get('position', 0),
                new_residue=data.get('new_residue', '')
            )
            results.append(result)
        
        return results
    
    def load_configuration(self, config_file: str) -> ProteinOptimizationConfig:
        with open(config_file, 'r') as f:
            config_data = json.load(f)
        
        return ProteinOptimizationConfig.from_dict(config_data)


class DataExporter:
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def export_to_csv(self, 
                     results: List[MutationResult], 
                     output_file: str) -> None:
        import csv
        
        if not results:
            self.logger.warning("No results to export")
            return
        
        fieldnames = [
            'iteration', 'score', 'n_terminal_energy', 'c_terminal_energy',
            'n_diff', 'c_diff', 'accepted', 'acceptance_probability',
            'mutation_type', 'best_structure', 'position', 'new_residue'
        ]
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in results:
                writer.writerow(result.to_dict())
        
        self.logger.info(f"Exported {len(results)} results to {output_file}")
    
    def export_summary_stats(self, 
                           history: MutationHistory, 
                           output_file: str) -> None:
        summary = history.get_summary_stats()
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.logger.info(f"Exported summary statistics to {output_file}")
    
    def export_trajectory_data(self, 
                              results: List[MutationResult], 
                              output_file: str) -> None:
        trajectory_data = {
            'iterations': [r.iteration for r in results],
            'scores': [r.score for r in results],
            'n_terminal_energies': [r.n_terminal_energy for r in results],
            'c_terminal_energies': [r.c_terminal_energy for r in results],
            'n_diffs': [r.n_diff for r in results],
            'c_diffs': [r.c_diff for r in results],
            'accepted': [r.accepted for r in results],
            'acceptance_probabilities': [r.acceptance_probability for r in results]
        }
        
        with open(output_file, 'w') as f:
            json.dump(trajectory_data, f, indent=2)
        
        self.logger.info(f"Exported trajectory data to {output_file}")


class StructureManager:
    
    def __init__(self, working_directory: str = "."):
        self.working_dir = Path(working_directory)
        self.logger = logging.getLogger(__name__)
    
    def extract_chain(self, 
                     input_pdb: str, 
                     chain_id: str = "A", 
                     output_pdb: str = "output_chain.pdb") -> str:
        output_path = self.working_dir / output_pdb
        
        with open(input_pdb, 'r') as infile, \
             open(output_path, 'w') as outfile:
            
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    if len(line) > 21 and line[21] == chain_id:
                        outfile.write(line)
                elif line.startswith('TER'):
                    break
                elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'REMARK')):
                    outfile.write(line)
        
        self.logger.debug(f"Extracted chain {chain_id} to {output_path}")
        return str(output_path)
    
    def validate_pdb_structure(self, pdb_file: str) -> Dict[str, Any]:
        validation = {
            'valid': False,
            'num_atoms': 0,
            'num_residues': 0,
            'chains': set(),
            'errors': []
        }
        
        if not Path(pdb_file).exists():
            validation['errors'].append(f"File does not exist: {pdb_file}")
            return validation
        
        try:
            residue_numbers = set()
            
            with open(pdb_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith(('ATOM', 'HETATM')):
                        if len(line) < 54:
                            validation['errors'].append(f"Line {line_num}: Invalid ATOM record length")
                            continue
                        
                        validation['num_atoms'] += 1
                        
                        # Extract chain and residue number
                        if len(line) > 21:
                            validation['chains'].add(line[21])
                        
                        if len(line) > 26:
                            try:
                                res_num = int(line[22:26].strip())
                                residue_numbers.add(res_num)
                            except ValueError:
                                validation['errors'].append(f"Line {line_num}: Invalid residue number")
            
            validation['num_residues'] = len(residue_numbers)
            validation['chains'] = list(validation['chains'])
            
            # Basic validation checks
            if validation['num_atoms'] == 0:
                validation['errors'].append("No atoms found in PDB file")
            elif validation['num_atoms'] < 10:
                validation['errors'].append("Very few atoms found (< 10)")
            
            if validation['num_residues'] == 0:
                validation['errors'].append("No residues found")
            
            validation['valid'] = len(validation['errors']) == 0
            
        except Exception as e:
            validation['errors'].append(f"Error reading PDB file: {e}")
        
        return validation
    
    def renumber_residues(self, 
                         input_pdb: str, 
                         output_pdb: str,
                         start_number: int = 1) -> str:
        output_path = self.working_dir / output_pdb
        
        current_res_num = None
        new_res_num = start_number - 1
        
        with open(input_pdb, 'r') as infile, \
             open(output_path, 'w') as outfile:
            
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    # Get current residue number
                    old_res_num = line[22:26].strip()
                    
                    # Check if we've moved to a new residue
                    if old_res_num != current_res_num:
                        current_res_num = old_res_num
                        new_res_num += 1
                    
                    # Replace residue number in the line
                    new_line = line[:22] + f"{new_res_num:4d}" + line[26:]
                    outfile.write(new_line)
                else:
                    outfile.write(line)
        
        self.logger.debug(f"Renumbered residues in {output_path}")
        return str(output_path)


class LogFileManager:
    
    def __init__(self, log_directory: str = "logs"):
        self.log_dir = Path(log_directory)
        self.log_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger(__name__)
    
    def create_optimization_log(self, run_id: str) -> str:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = self.log_dir / f"optimization_{run_id}_{timestamp}.log"
        
        return str(log_file)
    
    def archive_logs(self, run_id: str, results_dir: Path) -> None:
        log_archive_dir = results_dir / "logs"
        log_archive_dir.mkdir(exist_ok=True)
        
        # Copy all log files for this run
        for log_file in self.log_dir.glob(f"*{run_id}*"):
            shutil.copy2(log_file, log_archive_dir / log_file.name)
        
        self.logger.info(f"Archived logs to {log_archive_dir}")
    
    def cleanup_old_logs(self, days_to_keep: int = 30) -> int:
        import time
        
        cutoff_time = time.time() - (days_to_keep * 24 * 3600)
        removed_count = 0
        
        for log_file in self.log_dir.glob("*.log"):
            if log_file.stat().st_mtime < cutoff_time:
                try:
                    log_file.unlink()
                    removed_count += 1
                    self.logger.debug(f"Removed old log file: {log_file}")
                except OSError as e:
                    self.logger.warning(f"Could not remove {log_file}: {e}")
        
        return removed_count.working_dir = Path(working_directory)
        self.logger = logging.getLogger(__name__)
    
    def cleanup_temporary_files(self, patterns: Optional[List[str]] = None) -> List[str]:
        if patterns is None:
            patterns = [
                "output*",
                "temp_*",
                "*.err",
                "*.fasc",
                "mut.resfile"
            ]
        
        removed_files = []
        
        for pattern in patterns:
            file_pattern = self.working_dir / pattern
            for file_path in glob.glob(str(file_pattern)):
                try:
                    os.remove(file_path)
                    removed_files.append(file_path)
                    self.logger.debug(f"Removed temporary file: {file_path}")
                except OSError as e:
                    self.logger.warning(f"Could not remove {file_path}: {e}")
        
        return removed_files
    
    def backup_input_structure(self, input_pdb: str, iteration: int) -> str:
        backup_name = f"backup_iteration_{iteration}_{Path(input_pdb).name}"
        backup_path = self.working_dir / backup_name
        
        try:
            shutil.copy2(input_pdb, backup_path)
            self.logger.debug(f"Backed up {input_pdb} to {backup_path}")
            return str(backup_path)
        except Exception as e:
            self.logger.error(f"Failed to backup {input_pdb}: {e}")
            raise
    
    def save_accepted_structure(self, 
                              structure_file: str, 
                              iteration: int, 
                              mutation_type: str) -> str:
        output_name = f"accepted_{iteration:04d}_{mutation_type}.pdb"
        output_path = self.working_dir / output_name
        
        try:
            shutil.copy2(structure_file, output_path)
            self.logger.info(f"Saved accepted structure: {output_path}")
            return str(output_path)
        except Exception as e:
            self.logger.error(f"Failed to save accepted structure: {e}")
            raise
    
    def save_rejected_structure(self, 
                              structure_file: str, 
                              iteration: int, 
                              mutation_type: str) -> str:
        output_name = f"rejected_{iteration:04d}_{mutation_type}.pdb"
        output_path = self.working_dir / output_name
        
        try:
            shutil.copy2(structure_file, output_path)
            self.logger.debug(f"Saved rejected structure: {output_path}")
            return str(output_path)
        except Exception as e:
            self.logger.warning(f"Failed to save rejected structure: {e}")
            return ""
    
    def create_results_directory(self, base_name: str = "results") -> Path:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = self.working_dir / f"{base_name}_{timestamp}"
        
        results_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Created results directory: {results_dir}")
        
        return results_dir
    
    def archive_optimization_run(self, 
                                results: List[MutationResult], 
                                config: ProteinOptimizationConfig,
                                results_dir: Optional[Path] = None) -> Path:
        if results_dir is None:
            results_dir = self.create_results_directory()
        
        # Save configuration
        config_path = results_dir / "config.json"
        with open(config_path, 'w') as f:
            json.dump(config.to_dict(), f, indent=2)
        
        # Save results
        results_data = [r.to_dict() for r in results]
        results_path = results_dir / "optimization_results.json"
        with open(results_path, 'w') as f:
            json.dump(results_data, f, indent=2)
        
        # Copy accepted structures
        accepted_dir = results_dir / "accepted_structures"
        accepted_dir.mkdir(exist_ok=True)
        
        # Copy rejected structures
        rejected_dir = results_dir / "rejected_structures"
        rejected_dir.mkdir(exist_ok=True)
        
        # Move structure files to appropriate directories
        for pattern in ["accepted_*.pdb", "rejected_*.pdb"]:
            for file_path in glob.glob(pattern):
                if file_path.startswith("accepted_"):
                    shutil.move(file_path, accepted_dir / Path(file_path).name)
                elif file_path.startswith("rejected_"):
                    shutil.move(file_path, rejected_dir / Path(file_path).name)
        
        self

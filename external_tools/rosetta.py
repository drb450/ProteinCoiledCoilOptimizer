import subprocess
import time
from typing import Tuple, List
from pathlib import Path
import logging

from config.optimization_config import ProteinOptimizationConfig
from .base import ExternalTool


class RosettaTool(ExternalTool):
    
    def __init__(self, config: ProteinOptimizationConfig):
        super().__init__(config.rosetta_path)
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run_relaxation(self) -> None:
        # Validate required files exist
        self._validate_input_files()
        
        # Build command
        cmd = self._build_relaxation_command()
        
        # Execute Rosetta
        try:
            self.logger.debug(f"Running Rosetta command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, 
                check=True, 
                capture_output=True, 
                text=True,
                timeout=3600  # 1 hour timeout
            )
            self.logger.debug("Rosetta relaxation completed successfully")
            
            # Give Rosetta time to write output files
            time.sleep(2)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Rosetta execution failed: {e}")
            self.logger.error(f"Stderr: {e.stderr}")
            raise
        except subprocess.TimeoutExpired:
            self.logger.error("Rosetta execution timed out")
            raise
    
    def _validate_input_files(self) -> None:
        required_files = [
            self.config.protocol_file,
            self.config.input_pdb,
            "mut.resfile"  # Should be created by mutation generator
        ]
        
        for file_path in required_files:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"Required file not found: {file_path}")
    
    def _build_relaxation_command(self) -> List[str]:
        cmd = [
            self.executable_path,
            '-parser:protocol', self.config.protocol_file,
            '-s', self.config.input_pdb,
            '-nstruct', str(self.config.n_structures),
            '-out:prefix', 'temp_',
            '-packing:resfile', 'mut.resfile'
        ]
        
        # Add additional common Rosetta flags
        cmd.extend([
            '-out:no_nstruct_label',
            '-out:overwrite',
            '-mute', 'all'  # Reduce output verbosity
        ])
        
        return cmd
  
    def analyze_scores(self) -> Tuple[float, str]:

        score_file = 'temp_score.sc'
        
        try:
            with open(score_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Score file not found: {score_file}. Rosetta may have failed.")
        
        if len(lines) < 3:
            raise ValueError("Insufficient score data in output file")
        
        # Parse scores from score file
        scores_data = []
        for line in lines[2:]:  # Skip header lines
            try:
                parts = line.split()
                if len(parts) < 22:  # Ensure we have enough columns
                    continue
                    
                score = float(parts[1])
                structure_name = parts[21] + '.pdb'
                scores_data.append((score, structure_name))
                
            except (IndexError, ValueError) as e:
                self.logger.warning(f"Could not parse score from line: {line.strip()}")
                continue
        
        if not scores_data:
            raise ValueError("No valid scores found in output file")
        
        # Find best score
        best_score, best_structure = min(scores_data, key=lambda x: x[0])
        
        self.logger.debug(f"Best score: {best_score}, Structure: {best_structure}")
        return best_score, best_structure
    
    def get_score_summary(self) -> dict:
        score_file = 'temp_score.sc'
        
        try:
            with open(score_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            return {}
        
        scores = []
        for line in lines[2:]:  # Skip headers
            try:
                score = float(line.split()[1])
                scores.append(score)
            except (IndexError, ValueError):
                continue
        
        if not scores:
            return {}
        
        return {
            'num_structures': len(scores),
            'best_score': min(scores),
            'worst_score': max(scores),
            'mean_score': sum(scores) / len(scores),
            'score_range': max(scores) - min(scores)
        }
    
    def cleanup_output_files(self) -> None:
        import glob
        import os
        
        patterns = ["temp_*.pdb", "temp_score.sc", "*.fasc"]
        
        for pattern in patterns:
            for file_path in glob.glob(pattern):
                try:
                    os.remove(file_path)
                    self.logger.debug(f"Removed file: {file_path}")
                except OSError as e:
                    self.logger.warning(f"Could not remove {file_path}: {e}")
    
    def validate_output(self) -> bool:

        score_file = Path('temp_score.sc')
        if not score_file.exists():
            return False
        
        # Check if we have any structure files
        import glob
        structure_files = glob.glob("temp_*.pdb")
        
        return len(structure_files) > 0

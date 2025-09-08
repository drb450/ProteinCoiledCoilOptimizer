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
                lines = f.

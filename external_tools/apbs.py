import subprocess
import glob
import os
from typing import Tuple, List
from pathlib import Path
import logging

from config.optimization_config import ProteinOptimizationConfig
from .base import ExternalTool


class PDB2PQRTool(ExternalTool):
    
    def __init__(self, config: ProteinOptimizationConfig):
        super().__init__(config.pdb2pqr_path)
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def prepare_structure(self, structure_file: str) -> None:
        if not Path(structure_file).exists():
            raise FileNotFoundError(f"Structure file not found: {structure_file}")
        
        # Extract single chain first
        self._extract_single_chain(structure_file)
        
        # Run PDB2PQR
        self._run_pdb2pqr(structure_file)
        
        # Handle potential errors
        self._handle_pdb2pqr_errors(structure_file)
    
    def _extract_single_chain(self, structure_file: str) -> None:
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
    
    def _run_pdb2pqr(self, structure_file: str) -> None:
        """Run PDB2PQR with standard parameters."""
        cmd = [
            self.executable_path,
            f'--apbs-input=output_pqr.in',
            f'--with-ph={self.config.ph_value}',
            '--titration-state-method=propka',
            '--ff=AMBER',
            structure_file,
            'output_pqr'
        ]
        
        try:
            self.logger.debug(f"Running PDB2PQR: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.debug("PDB2PQR completed successfully")
            
        except subprocess.CalledProcessError as e:
            self.logger.warning("PDB2PQR failed, retrying with --noopts")
            # Try with --noopts flag if first attempt fails
            cmd.insert(-2, '--noopts')
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                self.logger.debug("PDB2PQR completed with --noopts")
            except subprocess.CalledProcessError as e2:
                self.logger.error(f"PDB2PQR failed even with --noopts: {e2}")
                raise
    
    def _handle_pdb2pqr_errors(self, structure_file: str) -> None:
        for err_file in glob.glob("*.err"):
            try:
                with open(err_file, 'r') as check_err:
                    content = check_err.read()
                    if "TypeError: '<' not supported between instances of" in content:
                        self.logger.warning("PDB2PQR TypeError detected, rerunning with --noopts")
                        self._run_pdb2pqr_with_noopts(structure_file)
                        break
            except Exception as e:
                self.logger.warning(f"Error handling PDB2PQR error file {err_file}: {e}")
    
    def _run_pdb2pqr_with_noopts(self, structure_file: str) -> None:
          cmd = [
            self.executable_path,
            f'--apbs-input=output_pqr.in',
            f'--with-ph={self.config.ph_value}',
            '--titration-state-method=propka',
            '--ff=AMBER',
            '--noopts',
            structure_file,
            'output_pqr'
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    
    def validate_output(self) -> bool:
        required_files = ['output_pqr.pqr', 'output_pqr.in']
        return all(Path(f).exists() for f in required_files)


class APBSTool(ExternalTool):
   
    def __init__(self, config: ProteinOptimizationConfig):
        super().__init__(config.apbs_path)
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def calculate_electrostatic_energies(self) -> Tuple[float, float]:
        # Prepare APBS input file
        self._prepare_apbs_input()
        
        # Run APBS
        self._run_apbs()
        
        # Parse results
        return self._parse_apbs_output()
    
    def _prepare_apbs_input(self) -> None:
        input_file = 'output_pqr.in'
        output_file = 'output.in'
        
        try:
            with open(input_file, 'r') as f_in, \
                 open(output_file, 'w') as f_out:
                for line in f_in:
                    # Replace 'total' with 'comps' to get component energies
                    modified_line = line.replace('total', 'comps')
                    f_out.write(modified_line)
        except FileNotFoundError:
            raise FileNotFoundError(f"APBS input file not found: {input_file}")
    
    def _run_apbs(self) -> None:
        cmd = [self.executable_path, '--output-file=output_apbs', 'output.in']
        
        try:
            self.logger.debug(f"Running APBS: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.debug("APBS calculation completed successfully")
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"APBS execution failed: {e}")
            self.logger.error(f"Stderr: {e.stderr}")
            raise
    
    def _parse_apbs_output(self) -> Tuple[float, float]:
        output_file = 'output_apbs'
        
        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"APBS output file not found: {output_file}")
        
        # Extract electrostatic energies
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
        
        # Split energies and calculate terminal energies
        mid_point = len(elec_energies) // 2
        polar_energies = elec_energies[mid_point:]
        
        # Take first 5 and last 5 residues for terminal calculations
        n_terminal_energy = sum(polar_energies[:5])
        c_terminal_energy = sum(polar_energies[-5:])
        
        self.logger.debug(f"N-terminal energy: {n_terminal_energy}")
        self.logger.debug(f"C-terminal energy: {c_terminal_energy}")
        
        return n_terminal_energy, c_terminal_energy
    
    def get_all_energies(self) -> List[float]:
        output_file = 'output_apbs'
        
        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            return []
        
        energies = []
        for line in lines:
            if line.startswith('atom'):
                try:
                    energy = float(line.split()[2])
                    energies.append(energy)
                except (IndexError, ValueError):
                    continue
        
        return energies
    
    def validate_output(self) -> bool:
        return Path('output_apbs').exists()


class ElectrostaticCalculator:
   
    def __init__(self, config: ProteinOptimizationConfig):
        self.config = config
        self.pdb2pqr = PDB2PQRTool(config)
        self.apbs = APBSTool(config)
        self.logger = logging.getLogger(__name__)
    
    def calculate_terminal_energies(self, structure_file: str) -> Tuple[float, float]:
        self.logger.debug(f"Calculating electrostatic energies for {structure_file}")
        
        # Prepare structure
        self.pdb2pqr.prepare_structure(structure_file)
        
        # Calculate energies
        n_energy, c_energy = self.apbs.calculate_electrostatic_energies()
        
        return n_energy, c_energy
    
    def validate_tools(self) -> dict:
        return {
            'pdb2pqr_available': self.pdb2pqr.check_availability(),
            'apbs_available': self.apbs.check_availability()
        }
    
    def cleanup_files(self) -> None:
        patterns = ["output*", "*.err"]
        
        for pattern in patterns:
            for file_path in glob.glob(pattern):
                try:
                    os.remove(file_path)
                    self.logger.debug(f"Removed file: {file_path}")
                except OSError as e:
                    self.logger.warning(f"Could not remove {file_path}: {e}")

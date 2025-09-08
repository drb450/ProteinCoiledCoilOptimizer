from typing import List, Optional, Dict, Any
from pathlib import Path


class ProteinOptimizationConfig:
   
    def __init__(self, 
                 initial_score: float = -658.0,
                 n_terminal_goal: float = 26678.18,
                 c_terminal_goal: float = 79814.70,
                 initial_n_diff: Optional[float] = None,
                 initial_c_diff: Optional[float] = None,
                 a_d_residues: Optional[List[int]] = None,
                 e_g_residues: Optional[List[int]] = None,
                 favorable_residues: Optional[List[str]] = None,
                 secondary_residues: Optional[List[str]] = None,
                 zeta: float = 3.93e-5,
                 alpha: float = 1.96e-4,
                 beta: float = 1.31e-4,
                 max_iterations: int = 1000,
                 n_structures: int = 5,
                 ph_value: float = 8.0,
                 rosetta_path: str = "/share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease",
                 pdb2pqr_path: str = "/share/apps/pdb2pqr/3.1.0/bin/pdb2pqr30",
                 apbs_path: str = "/share/apps/apbs",
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
        self.apbs_path = apbs_path
        self.protocol_file = protocol_file
        self.input_pdb = input_pdb
    
    def __repr__(self) -> str:

        return (f"ProteinOptimizationConfig(initial_score={self.initial_score}, "
                f"max_iterations={self.max_iterations}, "
                f"n_terminal_goal={self.n_terminal_goal})")
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, ProteinOptimizationConfig):
            return False
        return (self.initial_score == other.initial_score and
                self.n_terminal_goal == other.n_terminal_goal and
                self.c_terminal_goal == other.c_terminal_goal and
                self.max_iterations == other.max_iterations)
    
    def validate_paths(self) -> List[str]:
        errors = []
        
        # Check required files
        required_files = [self.protocol_file, self.input_pdb]
        for file_path in required_files:
            if not Path(file_path).exists():
                errors.append(f"Required file not found: {file_path}")
        
        # Check executables (basic existence check)
        required_executables = [self.rosetta_path, self.pdb2pqr_path, self.apbs_path]
        for executable in required_executables:
            if not Path(executable).exists():
                errors.append(f"Executable not found: {executable}")
        
        return errors
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'initial_score': self.initial_score,
            'n_terminal_goal': self.n_terminal_goal,
            'c_terminal_goal': self.c_terminal_goal,
            'initial_n_diff': self.initial_n_diff,
            'initial_c_diff': self.initial_c_diff,
            'a_d_residues': self.a_d_residues,
            'e_g_residues': self.e_g_residues,
            'favorable_residues': self.favorable_residues,
            'secondary_residues': self.secondary_residues,
            'zeta': self.zeta,
            'alpha': self.alpha,
            'beta': self.beta,
            'max_iterations': self.max_iterations,
            'n_structures': self.n_structures,
            'ph_value': self.ph_value,
            'rosetta_path': self.rosetta_path,
            'pdb2pqr_path': self.pdb2pqr_path,
            'apbs_path': self.apbs_path,
            'protocol_file': self.protocol_file,
            'input_pdb': self.input_pdb
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'ProteinOptimizationConfig':
        """Create configuration from dictionary."""
        return cls(**config_dict)

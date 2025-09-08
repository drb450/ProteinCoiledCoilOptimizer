import random
from typing import Dict, Any, Tuple, List
from config.optimization_config import ProteinOptimizationConfig


class MutationResult:
    """Results from a single mutation iteration."""
    
    def __init__(self, 
                 iteration: int, 
                 score: float, 
                 n_terminal_energy: float, 
                 c_terminal_energy: float, 
                 n_diff: float, 
                 c_diff: float,
                 accepted: bool, 
                 acceptance_probability: float, 
                 mutation_type: str, 
                 best_structure: str,
                 position: int,
                 new_residue: str):

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
        self.position = position
        self.new_residue = new_residue
    
    def __repr__(self) -> str:
        return (f"MutationResult(iteration={self.iteration}, "
                f"score={self.score:.2f}, accepted={self.accepted}, "
                f"mutation={self.position}{self.new_residue})")
    
    def to_dict(self) -> Dict[str, Any]:
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
            'best_structure': self.best_structure,
            'position': self.position,
            'new_residue': self.new_residue
        }


class MutationGenerator:

    
    def __init__(self, config: ProteinOptimizationConfig):

        self.config = config
    
    def create_mutation_resfile(self, iteration: int) -> Tuple[int, str]:

        # Choose between a/d and e/g positions
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
        
        return target_residue, new_residue
    
    def create_targeted_mutation(self, 
                               position: int, 
                               residue: str) -> Tuple[int, str]:

        resfile_content = f"NATAA \nstart \n\n{position} A PIKAA {residue}\n"
        
        with open("mut.resfile", 'w') as resfile:
            resfile.write(resfile_content)
        
        return position, residue
    
    def get_mutation_context(self, position: int) -> Dict[str, Any]:

        is_ad_position = position in self.config.a_d_residues
        is_eg_position = position in self.config.e_g_residues
        
        if is_ad_position:
            position_type = "a/d"
            favorable_residues = self.config.favorable_residues
        elif is_eg_position:
            position_type = "e/g"
            favorable_residues = self.config.secondary_residues
        else:
            position_type = "other"
            favorable_residues = []
        
        return {
            'position': position,
            'position_type': position_type,
            'is_ad_position': is_ad_position,
            'is_eg_position': is_eg_position,
            'favorable_residues': favorable_residues
        }


class MutationHistory:
    
    def __init__(self):
        self.results: List[MutationResult] = []
    
    def add_result(self, result: MutationResult) -> None:
        self.results.append(result)
    
    def get_accepted_mutations(self) -> List[MutationResult]:
        return [r for r in self.results if r.accepted]
    
    def get_rejected_mutations(self) -> List[MutationResult]:
        return [r for r in self.results if not r.accepted]
    
    def get_acceptance_rate(self) -> float:
        if not self.results:
            return 0.0
        return len(self.get_accepted_mutations()) / len(self.results)
    
    def get_best_score(self) -> float:
        if not self.results:
            return float('inf')
        return min(r.score for r in self.results)
    
    def get_mutation_frequency(self) -> Dict[int, int]:
        frequency = {}
        for result in self.results:
            position = result.position
            frequency[position] = frequency.get(position, 0) + 1
        return frequency
    
    def get_residue_frequency(self) -> Dict[str, int]:
        frequency = {}
        for result in self.results:
            residue = result.new_residue
            frequency[residue] = frequency.get(residue, 0) + 1
        return frequency
    
    def get_acceptance_by_mutation_type(self) -> Dict[str, float]:
        type_counts = {}
        type_accepted = {}
        
        for result in self.results:
            mut_type = result.mutation_type
            type_counts[mut_type] = type_counts.get(mut_type, 0) + 1
            if result.accepted:
                type_accepted[mut_type] = type_accepted.get(mut_type, 0) + 1
        
        return {mut_type: type_accepted.get(mut_type, 0) / count 
                for mut_type, count in type_counts.items()}
    
    def get_summary_stats(self) -> Dict[str, Any]:
        if not self.results:
            return {"total_iterations": 0}
        
        accepted = self.get_accepted_mutations()
        
        return {
            "total_iterations": len(self.results),
            "accepted_mutations": len(accepted),
            "acceptance_rate": self.get_acceptance_rate(),
            "best_score": self.get_best_score(),
            "final_score": self.results[-1].score,
            "score_improvement": self.results[0].score - self.results[-1].score,
            "mutation_frequency": self.get_mutation_frequency(),
            "residue_frequency": self.get_residue_frequency(),
            "acceptance_by_type": self.get_acceptance_by_mutation_type()
        }

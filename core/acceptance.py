import math
import random
from typing import Tuple
from config.optimization_config import ProteinOptimizationConfig


class AcceptanceCriteria:
    
    def __init__(self, config: ProteinOptimizationConfig):

        self.config = config
        self.index_range = list(range(100))
    
    def calculate_acceptance_probability(self, 
                                       score: float, 
                                       n_energy: float, 
                                       c_energy: float,
                                       current_score: float,
                                       current_n_energy: float,
                                       current_c_energy: float) -> Tuple[float, str]:

        # Calculate energy differences from goals
        n_diff = abs(n_energy - self.config.n_terminal_goal)
        c_diff = abs(c_energy - self.config.c_terminal_goal)
        current_n_diff = abs(current_n_energy - self.config.n_terminal_goal)
        current_c_diff = abs(current_c_energy - self.config.c_terminal_goal)
        
        # Determine improvements
        score_improved = score <= current_score
        n_improved = n_diff <= current_n_diff
        c_improved = c_diff <= current_c_diff
        
        return self._calculate_probability_by_conditions(
            score, n_diff, c_diff, current_score, current_n_diff, current_c_diff,
            score_improved, n_improved, c_improved
        )
    
    def _calculate_probability_by_conditions(self,
                                           score: float,
                                           n_diff: float,
                                           c_diff: float,
                                           current_score: float,
                                           current_n_diff: float,
                                           current_c_diff: float,
                                           score_improved: bool,
                                           n_improved: bool,
                                           c_improved: bool) -> Tuple[float, str]:
        
        if score_improved and n_improved and c_improved:
            # Accept unconditionally - all criteria improved
            return 1.0, "accept"
        
        elif not score_improved and n_improved and c_improved:
            # R condition - score degraded but energies improved
            delta_e = score - current_score
            exp_factor = delta_e / (current_score * self.config.zeta)
            probability = self._safe_exp(exp_factor)
            return probability, "score_penalty"
        
        elif score_improved and not n_improved and c_improved:
            # P condition - N-terminal degraded but others improved
            delta_e = n_diff - current_n_diff
            exp_factor = -delta_e / (current_n_diff * self.config.beta)
            probability = self._safe_exp(exp_factor)
            return probability, "n_terminal_penalty"
        
        elif score_improved and n_improved and not c_improved:
            # N condition - C-terminal degraded but others improved
            delta_e = c_diff - current_c_diff
            exp_factor = -delta_e / (current_c_diff * self.config.alpha)
            probability = self._safe_exp(exp_factor)
            return probability, "c_terminal_penalty"
        
        else:
            # Multiple penalties - combine probabilities
            return self._calculate_multi_penalty_probability(
                score, n_diff, c_diff, current_score, current_n_diff, current_c_diff,
                score_improved, n_improved, c_improved
            )
    
    def _calculate_multi_penalty_probability(self,
                                           score: float,
                                           n_diff: float,
                                           c_diff: float,
                                           current_score: float,
                                           current_n_diff: float,
                                           current_c_diff: float,
                                           score_improved: bool,
                                           n_improved: bool,
                                           c_improved: bool) -> Tuple[float, str]:
        prob_components = []
        mutation_type = "multi_penalty"
        
        if not score_improved:
            delta_e = score - current_score
            exp_factor = delta_e / (current_score * self.config.zeta)
            prob_components.append(self._safe_exp(exp_factor))
            mutation_type += "_R"
        
        if not n_improved:
            delta_e = n_diff - current_n_diff
            exp_factor = -delta_e / (current_n_diff * self.config.beta)
            prob_components.append(self._safe_exp(exp_factor))
            mutation_type += "_P"
        
        if not c_improved:
            delta_e = c_diff - current_c_diff
            exp_factor = -delta_e / (current_c_diff * self.config.alpha)
            prob_components.append(self._safe_exp(exp_factor))
            mutation_type += "_N"
        
        probability = 1.0 if not prob_components else math.prod(prob_components)
        return probability, mutation_type
    
    def _safe_exp(self, x: float) -> float:

        try:
            return math.exp(x)
        except OverflowError:
            return 0.0
    
    def should_accept_mutation(self, probability: float) -> bool:

        if probability >= 1.0:
            return True
        
        # Convert probability to percentage and use random sampling
        prob_percentage = min(int(round(probability * 100, 0)), 100)
        accept_indices = list(range(prob_percentage))
        random_index = random.choice(self.index_range)
        
        return random_index in accept_indices
    
    def get_acceptance_details(self,
                             score: float,
                             n_energy: float,
                             c_energy: float,
                             current_score: float,
                             current_n_energy: float,
                             current_c_energy: float) -> dict:

        n_diff = abs(n_energy - self.config.n_terminal_goal)
        c_diff = abs(c_energy - self.config.c_terminal_goal)
        current_n_diff = abs(current_n_energy - self.config.n_terminal_goal)
        current_c_diff = abs(current_c_energy - self.config.c_terminal_goal)
        
        score_improved = score <= current_score
        n_improved = n_diff <= current_n_diff
        c_improved = c_diff <= current_c_diff
        
        probability, mutation_type = self.calculate_acceptance_probability(
            score, n_energy, c_energy, current_score, current_n_energy, current_c_energy
        )
        
        return {
            'score_change': score - current_score,
            'n_energy_change': n_energy - current_n_energy,
            'c_energy_change': c_energy - current_c_energy,
            'n_diff_change': n_diff - current_n_diff,
            'c_diff_change': c_diff - current_c_diff,
            'score_improved': score_improved,
            'n_improved': n_improved,
            'c_improved': c_improved,
            'acceptance_probability': probability,
            'mutation_type': mutation_type,
            'will_accept': self.should_accept_mutation(probability)
        }

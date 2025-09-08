import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from pathlib import Path
import logging

from core.mutation import MutationResult


class ScoreAnalyzer:
   
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def analyze_score_trajectory(self, results: List[MutationResult]) -> Dict[str, Any]:
        if not results:
            return {}
        
        scores = [r.score for r in results]
        accepted_results = [r for r in results if r.accepted]
        accepted_scores = [r.score for r in accepted_results]
        
        analysis = {
            'total_iterations': len(results),
            'accepted_iterations': len(accepted_results),
            'acceptance_rate': len(accepted_results) / len(results) if results else 0,
            
            # Score statistics
            'initial_score': scores[0],
            'final_score': scores[-1],
            'best_score': min(scores),
            'worst_score': max(scores),
            'score_improvement': scores[0] - scores[-1],
            'score_range': max(scores) - min(scores),
            'mean_score': np.mean(scores),
            'std_score': np.std(scores),
            
            # Accepted score statistics
            'accepted_score_improvement': (accepted_scores[0] - accepted_scores[-1] 
                                         if accepted_scores else 0),
            'accepted_mean_score': np.mean(accepted_scores) if accepted_scores else 0,
            'accepted_std_score': np.std(accepted_scores) if accepted_scores else 0,
        }
        
        # Add convergence analysis
        analysis.update(self._analyze_convergence(scores, accepted_scores))
        
        return analysis
    
    def _analyze_convergence(self, all_scores: List[float], 
                           accepted_scores: List[float]) -> Dict[str, Any]:
        convergence = {}
        
        if len(all_scores) > 10:
            # Calculate moving averages
            window_size = min(10, len(all_scores) // 4)
            moving_avg = self._calculate_moving_average(all_scores, window_size)
            convergence['moving_average_final'] = moving_avg[-1] if moving_avg else 0

            plateau_threshold = 0.1  # REU threshold for plateau detection
            plateau_length = self._detect_plateau(all_scores, plateau_threshold)
            convergence['plateau_length'] = plateau_length

            if len(all_scores) > 1:
                improvement_rate = (all_scores[0] - all_scores[-1]) / len(all_scores)
                convergence['improvement_rate_per_iteration'] = improvement_rate
        
        return convergence
    
    def _calculate_moving_average(self, scores: List[float], window_size: int) -> List[float]:
        if len(scores) < window_size:
            return []
        
        moving_avg = []
        for i in range(window_size - 1, len(scores)):
            window = scores[i - window_size + 1:i + 1]
            moving_avg.append(np.mean(window))
        
        return moving_avg
    
    def _detect_plateau(self, scores: List[float], threshold: float) -> int:
        if len(scores) < 5:
            return 0
        
        plateau_length = 0
        last_score = scores[-1]
        
        for i in range(len(scores) - 1, -1, -1):
            if abs(scores[i] - last_score) <= threshold:
                plateau_length += 1
            else:
                break
        
        return plateau_length
    
    def analyze_energy_trajectory(self, results: List[MutationResult]) -> Dict[str, Any]:
        if not results:
            return {}
        
        n_energies = [r.n_terminal_energy for r in results]
        c_energies = [r.c_terminal_energy for r in results]
        n_diffs = [r.n_diff for r in results]
        c_diffs = [r.c_diff for r in results]
        
        return {
            # N-terminal analysis
            'n_terminal': {
                'initial': n_energies[0],
                'final': n_energies[-1],
                'mean': np.mean(n_energies),
                'std': np.std(n_energies),
                'min': min(n_energies),
                'max': max(n_energies),
                'improvement': n_energies[0] - n_energies[-1]
            },
            
            # C-terminal analysis
            'c_terminal': {
                'initial': c_energies[0],
                'final': c_energies[-1],
                'mean': np.mean(c_energies),
                'std': np.std(c_energies),
                'min': min(c_energies),
                'max': max(c_energies),
                'improvement': c_energies[0] - c_energies[-1]
            },
            
            # Goal differences
            'n_diff_improvement': n_diffs[0] - n_diffs[-1],
            'c_diff_improvement': c_diffs[0] - c_diffs[-1],
            'final_n_diff': n_diffs[-1],
            'final_c_diff': c_diffs[-1]
        }
    
    def analyze_acceptance_patterns(self, results: List[MutationResult]) -> Dict[str, Any]:
        if not results:
            return {}
        
        # Group by mutation type
        type_analysis = {}
        for result in results:
            mut_type = result.mutation_type
            if mut_type not in type_analysis:
                type_analysis[mut_type] = {
                    'count': 0,
                    'accepted': 0,
                    'avg_probability': 0,
                    'scores': []
                }
            
            type_analysis[mut_type]['count'] += 1
            type_analysis[mut_type]['scores'].append(result.score)
            type_analysis[mut_type]['avg_probability'] += result.acceptance_probability
            
            if result.accepted:
                type_analysis[mut_type]['accepted'] += 1
        
        # Calculate final statistics
        for mut_type in type_analysis:
            data = type_analysis[mut_type]
            data['acceptance_rate'] = data['accepted'] / data['count']
            data['avg_probability'] /= data['count']
            data['avg_score'] = np.mean(data['scores'])
            data['std_score'] = np.std(data['scores'])
        
        return type_analysis
    
    def detect_optimization_phases(self, results: List[MutationResult]) -> List[Dict[str, Any]]:
        if len(results) < 20:
            return [{'phase': 'insufficient_data', 'length': len(results)}]
        
        phases = []
        scores = [r.score for r in results]
        
        window_size = 20
        improvement_threshold = 0.5  # REU per iteration
        
        i = 0
        while i < len(scores) - window_size:
            window_scores = scores[i:i + window_size]
            improvement_rate = (window_scores[0] - window_scores[-1]) / window_size
            
            if improvement_rate > improvement_threshold:
                phase_type = 'rapid_improvement'
            elif improvement_rate > 0.1:
                phase_type = 'gradual_improvement'
            elif improvement_rate > -0.1:
                phase_type = 'plateau'
            else:
                phase_type = 'degradation'
            
            phase_end = i + window_size
            while (phase_end < len(scores) - window_size and 
                   self._same_phase_type(scores[phase_end:phase_end + window_size], phase_type)):
                phase_end += window_size
            
            phases.append({
                'phase': phase_type,
                'start': i,
                'end': min(phase_end, len(scores)),
                'length': min(phase_end,

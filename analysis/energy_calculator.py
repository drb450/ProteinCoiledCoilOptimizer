import numpy as np
from typing import List, Dict, Any, Tuple, Optional
import logging

from core.mutation import MutationResult
from config.optimization_config import ProteinOptimizationConfig


class EnergyAnalyzer:
    
    def __init__(self, config: ProteinOptimizationConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def calculate_energy_deviations(self, results: List[MutationResult]) -> Dict[str, Any]:
        if not results:
            return {}
        
        n_deviations = []
        c_deviations = []
        n_energies = []
        c_energies = []
        
        for result in results:
            n_dev = abs(result.n_terminal_energy - self.config.n_terminal_goal)
            c_dev = abs(result.c_terminal_energy - self.config.c_terminal_goal)
            
            n_deviations.append(n_dev)
            c_deviations.append(c_dev)
            n_energies.append(result.n_terminal_energy)
            c_energies.append(result.c_terminal_energy)
        
        return {
            'n_terminal_deviations': {
                'initial': n_deviations[0],
                'final': n_deviations[-1],
                'best': min(n_deviations),
                'worst': max(n_deviations),
                'mean': np.mean(n_deviations),
                'std': np.std(n_deviations),
                'improvement': n_deviations[0] - n_deviations[-1],
                'goal': self.config.n_terminal_goal
            },
            'c_terminal_deviations': {
                'initial': c_deviations[0],
                'final': c_deviations[-1],
                'best': min(c_deviations),
                'worst': max(c_deviations),
                'mean': np.mean(c_deviations),
                'std': np.std(c_deviations),
                'improvement': c_deviations[0] - c_deviations[-1],
                'goal': self.config.c_terminal_goal
            },
            'combined_deviation': {
                'initial': n_deviations[0] + c_deviations[0],
                'final': n_deviations[-1] + c_deviations[-1],
                'improvement': (n_deviations[0] + c_deviations[0]) - (n_deviations[-1] + c_deviations[-1])
            }
        }
    
    def analyze_energy_correlations(self, results: List[MutationResult]) -> Dict[str, float]:
        if len(results) < 3:
            return {}
        
        scores = [r.score for r in results]
        n_energies = [r.n_terminal_energy for r in results]
        c_energies = [r.c_terminal_energy for r in results]
        n_diffs = [r.n_diff for r in results]
        c_diffs = [r.c_diff for r in results]
        
        correlations = {}
        
        try:
            correlations['score_vs_n_energy'] = np.corrcoef(scores, n_energies)[0, 1]
            correlations['score_vs_c_energy'] = np.corrcoef(scores, c_energies)[0, 1]
            correlations['n_energy_vs_c_energy'] = np.corrcoef(n_energies, c_energies)[0, 1]
            correlations['score_vs_n_diff'] = np.corrcoef(scores, n_diffs)[0, 1]
            correlations['score_vs_c_diff'] = np.corrcoef(scores, c_diffs)[0, 1]
            correlations['n_diff_vs_c_diff'] = np.corrcoef(n_diffs, c_diffs)[0, 1]
        except Exception as e:
            self.logger.warning(f"Error calculating correlations: {e}")
        
        return correlations
    
    def find_pareto_optimal_solutions(self, results: List[MutationResult]) -> List[MutationResult]:
        if not results:
            return []
        
        pareto_solutions = []
        
        for i, result_i in enumerate(results):
            is_dominated = False
            
            for j, result_j in enumerate(results):
                if i == j:
                    continue
                
                # Check if result_j dominates result_i
                if (result_j.score <= result_i.score and
                    result_j.n_diff <= result_i.n_diff and
                    result_j.c_diff <= result_i.c_diff and
                    (result_j.score < result_i.score or
                     result_j.n_diff < result_i.n_diff or
                     result_j.c_diff < result_i.c_diff)):
                    is_dominated = True
                    break
            
            if not is_dominated:
                pareto_solutions.append(result_i)
        
        # Sort by total score for easier interpretation
        pareto_solutions.sort(key=lambda x: x.score)
        
        return pareto_solutions
    
    def calculate_energy_landscape_metrics(self, results: List[MutationResult]) -> Dict[str, Any]:
        if len(results) < 10:
            return {}
        
        scores = [r.score for r in results]
        n_energies = [r.n_terminal_energy for r in results]
        c_energies = [r.c_terminal_energy for r in results]
        
        metrics = {}
        
        # Ruggedness metrics
        metrics['score_ruggedness'] = self._calculate_ruggedness(scores)
        metrics['n_energy_ruggedness'] = self._calculate_ruggedness(n_energies)
        metrics['c_energy_ruggedness'] = self._calculate_ruggedness(c_energies)
        
        # Smoothness metrics (inverse of ruggedness)
        metrics['score_smoothness'] = 1.0 / (1.0 + metrics['score_ruggedness'])
        
        # Gradient analysis
        score_gradients = np.diff(scores)
        metrics['score_gradient_mean'] = np.mean(score_gradients)
        metrics['score_gradient_std'] = np.std(score_gradients)
        
        # Local minima detection
        metrics['local_minima_count'] = self._count_local_minima(scores)
        
        return metrics
    
    def _calculate_ruggedness(self, values: List[float]) -> float:
        if len(values) < 2:
            return 0.0
        
        differences = [abs(values[i+1] - values[i]) for i in range(len(values) - 1)]
        return np.mean(differences)
    
    def _count_local_minima(self, values: List[float], window_size: int = 3) -> int:
        if len(values) < 2 * window_size + 1:
            return 0
        
        minima_count = 0
        for i in range(window_size, len(values) - window_size):
            is_minimum = True
            center_value = values[i]
            
            # Check if center is smaller than all neighbors in window
            for j in range(i - window_size, i + window_size + 1):
                if j != i and values[j] <= center_value:
                    is_minimum = False
                    break
            
            if is_minimum:
                minima_count += 1
        
        return minima_count
    
    def generate_energy_report(self, results: List[MutationResult]) -> str:
        if not results:
            return "No results to analyze."
        
        deviations = self.calculate_energy_deviations(results)
        correlations = self.analyze_energy_correlations(results)
        pareto_solutions = self.find_pareto_optimal_solutions(results)
        landscape = self.calculate_energy_landscape_metrics(results)
        
        report = []
        report.append("ENERGY ANALYSIS REPORT")
        report.append("=" * 40)
        report.append("")
        
        # Energy deviations
        report.append("ENERGY GOAL DEVIATIONS:")
        n_dev = deviations['n_terminal_deviations']
        c_dev = deviations['c_terminal_deviations']
        
        report.append(f"  N-terminal (goal: {n_dev['goal']:.1f}):")
        report.append(f"    Initial deviation: {n_dev['initial']:.1f}")
        report.append(f"    Final deviation: {n_dev['final']:.1f}")
        report.append(f"    Best deviation: {n_dev['best']:.1f}")
        report.append(f"    Improvement: {n_dev['improvement']:.1f}")
        report.append("")
        
        report.append(f"  C-terminal (goal: {c_dev['goal']:.1f}):")
        report.append(f"    Initial deviation: {c_dev['initial']:.1f}")
        report.append(f"    Final deviation: {c_dev['final']:.1f}")
        report.append(f"    Best deviation: {c_dev['best']:.1f}")
        report.append(f"    Improvement: {c_dev['improvement']:.1f}")
        report.append("")
        
        combined = deviations['combined_deviation']
        report.append(f"  Combined deviation improvement: {combined['improvement']:.1f}")
        report.append("")
        
        # Correlations
        if correlations:
            report.append("ENERGY CORRELATIONS:")
            for corr_name, corr_value in correlations.items():
                if not np.isnan(corr_value):
                    report.append(f"  {corr_name}: {corr_value:.3f}")
            report.append("")
        
        # Pareto solutions
        report.append(f"PARETO OPTIMAL SOLUTIONS: {len(pareto_solutions)} found")
        if pareto_solutions:
            report.append("  Best solutions (score, n_diff, c_diff):")
            for i, sol in enumerate(pareto_solutions[:5], 1):
                report.append(f"    {i}. Score: {sol.score:.2f}, N-diff: {sol.n_diff:.1f}, C-diff: {sol.c_diff:.1f}")
            report.append("")
        
        # Landscape metrics
        if landscape:
            report.append("ENERGY LANDSCAPE:")
            report.append(f"  Score ruggedness: {landscape['score_ruggedness']:.3f}")
            report.append(f"  Score smoothness: {landscape['score_smoothness']:.3f}")
            report.append(f"  Local minima count: {landscape['local_minima_count']}")
            if 'score_gradient_mean' in landscape:
                report.append(f"  Average score gradient: {landscape['score_gradient_mean']:.3f}")
        
        return "\n".join(report)


class EnergyOptimizationMetrics:
   
    def __init__(self, config: ProteinOptimizationConfig):
        """Initialize with configuration goals."""
        self.config = config
    
    def calculate_goal_achievement(self, results: List[MutationResult]) -> Dict[str, Any]:
        if not results:
            return {}
        
        final_result = results[-1]
        best_n_diff = min(r.n_diff for r in results)
        best_c_diff = min(r.c_diff for r in results)
        
        # Define achievement thresholds (could be made configurable)
        n_threshold = self.config.initial_n_diff * 0.1  # 10% of initial difference
        c_threshold = self.config.initial_c_diff * 0.1
        
        return {
            'n_terminal_achieved': final_result.n_diff < n_threshold,
            'c_terminal_achieved': final_result.c_diff < c_threshold,
            'both_achieved': (final_result.n_diff < n_threshold and 
                            final_result.c_diff < c_threshold),
            'n_improvement_percent': (self.config.initial_n_diff - final_result.n_diff) / self.config.initial_n_diff * 100,
            'c_improvement_percent': (self.config.initial_c_diff - final_result.c_diff) / self.config.initial_c_diff * 100,
            'best_n_achieved': best_n_diff < n_threshold,
            'best_c_achieved': best_c_diff < c_threshold,
            'thresholds': {
                'n_threshold': n_threshold,
                'c_threshold': c_threshold
            }
        }

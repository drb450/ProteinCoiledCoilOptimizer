import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any
import json


class OptimizationLogger:
  
    def __init__(self, 
                 name: str = "protein_optimization",
                 log_level: int = logging.INFO,
                 log_file: Optional[str] = None,
                 console_output: bool = True):
        self.name = name
        self.log_level = log_level
        self.log_file = log_file
        self.console_output = console_output
        
        self.logger = self._setup_logger()
    
    def _setup_logger(self) -> logging.Logger:

        logger = logging.getLogger(self.name)
        logger.setLevel(self.log_level)
        
        # Remove any existing handlers
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Console handler
        if self.console_output:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(self.log_level)
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
        
        # File handler
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(self.log_level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        
        return logger
    
    def get_logger(self) -> logging.Logger:
        return self.logger
    
    def log_optimization_start(self, config_dict: Dict[str, Any]) -> None:

        self.logger.info("=" * 60)
        self.logger.info("PROTEIN OPTIMIZATION STARTED")
        self.logger.info("=" * 60)
        self.logger.info(f"Configuration:")
        
        for key, value in config_dict.items():
            if isinstance(value, (list, dict)):
                self.logger.info(f"  {key}: {json.dumps(value, indent=2)}")
            else:
                self.logger.info(f"  {key}: {value}")
    
    def log_iteration_start(self, iteration: int, total_iterations: int) -> None:
        self.logger.info(f"--- Iteration {iteration + 1}/{total_iterations} ---")
    
    def log_mutation_details(self, position: int, residue: str, mutation_type: str) -> None:
        self.logger.debug(f"Mutation: Position {position} -> {residue} (type: {mutation_type})")
    
    def log_scores(self, score: float, n_energy: float, c_energy: float) -> None:
        self.logger.debug(f"Scores - Rosetta: {score:.2f}, N-term: {n_energy:.1f}, C-term: {c_energy:.1f}")
    
    def log_acceptance_decision(self, 
                              accepted: bool, 
                              probability: float, 
                              mutation_type: str) -> None:
        status = "ACCEPTED" if accepted else "REJECTED"
        self.logger.info(f"{status}: {mutation_type}, probability: {probability:.4f}")
    
    def log_optimization_complete(self, 
                                final_score: float, 
                                total_accepted: int, 
                                total_iterations: int) -> None:
        acceptance_rate = total_accepted / total_iterations * 100
        
        self.logger.info("=" * 60)
        self.logger.info("PROTEIN OPTIMIZATION COMPLETED")
        self.logger.info("=" * 60)
        self.logger.info(f"Final score: {final_score:.2f}")
        self.logger.info(f"Total iterations: {total_iterations}")
        self.logger.info(f"Accepted mutations: {total_accepted}")
        self.logger.info(f"Acceptance rate: {acceptance_rate:.1f}%")
    
    def log_error(self, error: Exception, context: str = "") -> None:
        if context:
            self.logger.error(f"Error in {context}: {str(error)}")
        else:
            self.logger.error(f"Error: {str(error)}")
        
        # Log full traceback at debug level
        self.logger.debug("Full traceback:", exc_info=True)
    
    def log_warning(self, message: str) -> None:
        self.logger.warning(message)
    
    def log_tool_status(self, tool_name: str, available: bool, path: str = "") -> None:
        status = "AVAILABLE" if available else "NOT AVAILABLE"
        if path:
            self.logger.info(f"Tool {tool_name}: {status} at {path}")
        else:
            self.logger.info(f"Tool {tool_name}: {status}")


class PerformanceLogger:
   
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.timers = {}
    
    def start_timer(self, name: str) -> None:
        import time
        self.timers[name] = time.time()
        self.logger.debug(f"Started timer: {name}")
    
    def stop_timer(self, name: str) -> float:
        import time
        
        if name not in self.timers:
            self.logger.warning(f"Timer '{name}' was not started")
            return 0.0
        
        duration = time.time() - self.timers[name]
        del self.timers[name]
        
        self.logger.debug(f"Timer {name}: {duration:.2f} seconds")
        return duration
    
    def log_iteration_timing(self, iteration: int, duration: float) -> None:
        self.logger.info(f"Iteration {iteration} completed in {duration:.1f} seconds")
    
    def log_tool_timing(self, tool_name: str, duration: float) -> None:
        self.logger.debug(f"{tool_name} execution time: {duration:.2f} seconds")
    
    def log_memory_usage(self) -> None:
        try:
            import psutil
            process = psutil.Process()
            memory_mb = process.memory_info().rss / 1024 / 1024
            self.logger.debug(f"Memory usage: {memory_mb:.1f} MB")
        except ImportError:
            self.logger.debug("Memory monitoring not available (psutil not installed)")


class StructuredLogger:
    
    def __init__(self, output_file: str):
        self.output_file = Path(output_file)
        self.events = []
        
        # Ensure directory exists
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    def log_event(self, event_type: str, data: Dict[str, Any]) -> None:
        event = {
            'timestamp': datetime.now().isoformat(),
            'event_type': event_type,
            'data': data
        }
        
        self.events.append(event)
    
    def log_iteration_event(self, 
                          iteration: int, 
                          score: float, 
                          n_energy: float,
                          c_energy: float,
                          accepted: bool,
                          probability: float) -> None:
        self.log_event('iteration_complete', {
            'iteration': iteration,
            'score': score,
            'n_terminal_energy': n_energy,
            'c_terminal_energy': c_energy,
            'accepted': accepted,
            'acceptance_probability': probability
        })
    
    def log_tool_event(self, tool_name: str, duration: float, success: bool) -> None:
        self.log_event('tool_execution', {
            'tool': tool_name,
            'duration_seconds': duration,
            'success': success
        })
    
    def flush_to_file(self) -> None:
        with open(self.output_file, 'w') as f:
            json.dump(self.events, f, indent=2)
    
    def get_events_by_type(self, event_type: str) -> list:
        return [event for event in self.events if event['event_type'] == event_type]


def setup_logging(log_level: str = "INFO",
                 log_file: Optional[str] = None,
                 console_output: bool = True,
                 structured_log_file: Optional[str] = None) -> tuple:
    # Convert string level to logging constant
    level_map = {
        'DEBUG': logging.DEBUG,
        'INFO': logging.INFO,
        'WARNING': logging.WARNING,
        'ERROR': logging.ERROR
    }
    
    numeric_level = level_map.get(log_level.upper(), logging.INFO)
    
    # Create default log file if none provided
    if log_file is None and console_output is False:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = f"protein_optimization_{timestamp}.log"
    
    # Setup main logger
    main_logger = OptimizationLogger(
        log_level=numeric_level,
        log_file=log_file,
        console_output=console_output
    )
    
    # Setup performance logger
    performance_logger = PerformanceLogger(main_logger.get_logger())
    
    # Setup structured logger if requested
    structured_logger = None
    if structured_log_file:
        structured_logger = StructuredLogger(structured_log_file)
    
    return main_logger, performance_logger, structured_logger

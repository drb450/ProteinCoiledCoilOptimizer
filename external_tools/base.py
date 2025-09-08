import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional
import logging


class ExternalTool(ABC):
    
    def __init__(self, executable_path: str):

        self.executable_path = executable_path
        self.logger = logging.getLogger(self.__class__.__name__)
        self._validate_executable()
    
    def _validate_executable(self) -> None:

        if not Path(self.executable_path).exists():
            # Try to find in PATH
            if not self._is_executable_in_path():
                self.logger.warning(f"Executable not found: {self.executable_path}")
    
    def _is_executable_in_path(self) -> bool:

        try:
            subprocess.run(['which', self.executable_path], 
                         check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError:
            return False
    
    def check_availability(self) -> bool:

        return (Path(self.executable_path).exists() or 
                self._is_executable_in_path())
    
    @abstractmethod
    def validate_output(self) -> bool:

        pass


class ToolChain:
    def __init__(self):

        self.tools = []
        self.logger = logging.getLogger(__name__)
    
    def add_tool(self, tool: ExternalTool, name: str) -> None:

        self.tools.append((name, tool))
    
    def validate_all_tools(self) -> dict:

        status = {}
        for name, tool in self.tools:
            available = tool.check_availability()
            status[name] = available
            
            if not available:
                self.logger.warning(f"Tool '{name}' is not available")
            else:
                self.logger.debug(f"Tool '{name}' is available")
        
        return status
    
    def get_missing_tools(self) -> list:

        missing = []
        for name, tool in self.tools:
            if not tool.check_availability():
                missing.append(name)
        return missing
    
    def is_ready(self) -> bool:

        return len(self.get_missing_tools()) == 0


class CommandBuilder:
    def __init__(self, executable: str):

        self.executable = executable
        self.args = []
    
    def add_flag(self, flag: str) -> 'CommandBuilder':

        self.args.append(flag)
        return self
    
    def add_option(self, option: str, value: str) -> 'CommandBuilder':

        self.args.extend([option, str(value)])
        return self
    
    def add_argument(self, arg: str) -> 'CommandBuilder':

        self.args.append(str(arg))
        return self
    
    def build(self) -> list:

        return [self.executable] + self.args
    
    def reset(self) -> 'CommandBuilder':

        self.args = []
        return self

"""
CAT-Surface Python Wrapper

A Python subprocess wrapper for CAT-Surface command-line tools.
This module provides a Pythonic interface to the CAT-Surface neuroimaging
tools for cortical surface analysis.

Author: GitHub Copilot
License: Same as CAT-Surface (GNU General Public License)
"""

import os
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
from dataclasses import dataclass
from abc import ABC, abstractmethod


class CATSurfaceError(Exception):
    """Base exception for CAT-Surface wrapper errors."""
    pass


class CATExecutableNotFoundError(CATSurfaceError):
    """Raised when CAT-Surface executable is not found."""
    pass


class CATProcessError(CATSurfaceError):
    """Raised when CAT-Surface process fails."""
    
    def __init__(self, message: str, returncode: int, cmd: List[str], 
                 stdout: str = "", stderr: str = ""):
        super().__init__(message)
        self.returncode = returncode
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr


@dataclass
class CATResult:
    """Result from a CAT-Surface operation."""
    returncode: int
    stdout: str
    stderr: str
    cmd: List[str]
    output_files: List[str]
    
    @property
    def success(self) -> bool:
        """True if the operation was successful."""
        return self.returncode == 0


class CATConfig:
    """Configuration for CAT-Surface wrapper."""
    
    def __init__(self, executable_path: Optional[str] = None):
        """
        Initialize CAT configuration.
        
        Args:
            executable_path: Path to CAT-Surface executables directory.
                           If None, will search standard locations.
        """
        self.executable_path = self._find_executable_path(executable_path)
        
    def _find_executable_path(self, provided_path: Optional[str] = None) -> str:
        """Find CAT-Surface executable directory."""
        if provided_path:
            if os.path.isdir(provided_path):
                return provided_path
            else:
                raise CATExecutableNotFoundError(
                    f"Provided executable path does not exist: {provided_path}"
                )
        
        # Search standard locations
        search_paths = [
            # Current project build directories
            "/Users/gaser/Library/CloudStorage/Dropbox/GitHub/CAT-Surface/build-native-arm64/Progs",
            "/Users/gaser/Library/CloudStorage/Dropbox/GitHub/CAT-Surface/build-native/Progs",
            "/Users/gaser/Library/CloudStorage/Dropbox/GitHub/CAT-Surface/build-x86_64-pc-linux/Progs",
            # Standard installation locations
            "/usr/local/bin",
            "/usr/bin",
            "/opt/CAT-Surface/bin",
        ]
        
        for path in search_paths:
            if os.path.isdir(path):
                # Check if at least one CAT executable exists
                test_executables = ["CAT_SurfConvert", "CAT_VolAmap", "CAT_SurfDistance"]
                for exe in test_executables:
                    if os.path.isfile(os.path.join(path, exe)):
                        return path
        
        raise CATExecutableNotFoundError(
            "CAT-Surface executables not found. Please specify executable_path."
        )
    
    def get_executable(self, name: str) -> str:
        """Get full path to a CAT executable."""
        exe_path = os.path.join(self.executable_path, name)
        if not os.path.isfile(exe_path):
            raise CATExecutableNotFoundError(f"Executable not found: {exe_path}")
        return exe_path


class CATToolBase(ABC):
    """Base class for CAT-Surface tool wrappers."""
    
    def __init__(self, config: Optional[CATConfig] = None):
        """
        Initialize CAT tool wrapper.
        
        Args:
            config: CAT configuration object. If None, uses default config.
        """
        self.config = config or CATConfig()
        
    @property
    @abstractmethod
    def executable_name(self) -> str:
        """Name of the CAT executable."""
        pass
        
    def _validate_file_exists(self, filepath: str, description: str = "File"):
        """Validate that a file exists."""
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"{description} not found: {filepath}")
            
    def _validate_output_dir(self, filepath: str):
        """Validate that output directory exists or can be created."""
        output_dir = os.path.dirname(filepath)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            
    def _run_command(self, args: List[str], timeout: Optional[int] = None) -> CATResult:
        """
        Execute CAT command with given arguments.
        
        Args:
            args: Command line arguments
            timeout: Command timeout in seconds
            
        Returns:
            CATResult object with execution results
        """
        executable = self.config.get_executable(self.executable_name)
        cmd = [executable] + args
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False
            )
            
            # Determine output files from arguments (basic heuristic)
            output_files = self._extract_output_files(args)
            
            cat_result = CATResult(
                returncode=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
                cmd=cmd,
                output_files=output_files
            )
            
            if result.returncode != 0:
                error_msg = f"CAT command failed with return code {result.returncode}"
                if result.stderr:
                    error_msg += f"\nStderr: {result.stderr}"
                raise CATProcessError(
                    error_msg, result.returncode, cmd, result.stdout, result.stderr
                )
                
            return cat_result
            
        except subprocess.TimeoutExpired:
            raise CATProcessError(f"Command timed out after {timeout} seconds", -1, cmd)
        except FileNotFoundError:
            raise CATExecutableNotFoundError(f"Executable not found: {executable}")
            
    def _extract_output_files(self, args: List[str]) -> List[str]:
        """Extract output file paths from command arguments (basic implementation)."""
        # This is a simple heuristic - could be improved per tool
        output_files = []
        for i, arg in enumerate(args):
            if not arg.startswith('-') and i > 0:
                # Skip first non-option argument (usually input)
                if '.' in arg and not args[i-1].startswith('-'):
                    output_files.append(arg)
        return output_files


# Global configuration instance
_default_config = None


def get_default_config() -> CATConfig:
    """Get the default CAT configuration."""
    global _default_config
    if _default_config is None:
        _default_config = CATConfig()
    return _default_config


def set_executable_path(path: str):
    """Set the global executable path for CAT-Surface tools."""
    global _default_config
    _default_config = CATConfig(path)


# Import convenience functions from submodules
from .volume import (
    segment_with_amap, denoise_sanlm, compute_thickness_pbt,
    extract_surface_marching_cubes, project_volume_to_surface
)
from .surface import (
    average_surfaces, create_pial_white_surfaces, calculate_surface_distance,
    warp_surface, create_spherical_surface, resample_surface
)

# Define what gets imported with "from cat_surface import *"
__all__ = [
    # Core classes and exceptions
    'CATSurfaceError', 'CATExecutableNotFoundError', 'CATProcessError',
    'CATResult', 'CATConfig', 'CATToolBase',
    
    # Configuration functions
    'get_default_config', 'set_executable_path',
    
    # Volume processing convenience functions
    'segment_with_amap', 'denoise_sanlm', 'compute_thickness_pbt',
    'extract_surface_marching_cubes', 'project_volume_to_surface',
    
    # Surface processing convenience functions
    'average_surfaces', 'create_pial_white_surfaces', 'calculate_surface_distance',
    'warp_surface', 'create_spherical_surface', 'resample_surface',
]
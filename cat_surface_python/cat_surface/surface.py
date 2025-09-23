"""
Surface processing tools for CAT-Surface Python wrapper.

This module contains wrappers for CAT-Surface surface processing and manipulation tools.
"""

from typing import Optional, List, Union
from . import CATToolBase, CATResult, get_default_config, CATConfig


class SurfAverage(CATToolBase):
    """Wrapper for CAT_SurfAverage - Calculate average and RMS of surfaces."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfAverage"
    
    def average_surfaces(self,
                        surface_files: List[str],
                        output_file: str) -> CATResult:
        """
        Calculate average of multiple surfaces.
        
        Args:
            surface_files: List of input surface files
            output_file: Output average surface file
            
        Returns:
            CATResult with execution details
        """
        for f in surface_files:
            self._validate_file_exists(f, "Surface file")
        self._validate_output_dir(output_file)
        
        args = ["-avg", output_file] + surface_files
        
        return self._run_command(args)
    
    def rms_surfaces(self,
                    surface_files: List[str],
                    rms_output: str) -> CATResult:
        """
        Calculate root mean square error of surfaces.
        
        Args:
            surface_files: List of input surface files
            rms_output: Output RMS values file
            
        Returns:
            CATResult with execution details
        """
        for f in surface_files:
            self._validate_file_exists(f, "Surface file")
        self._validate_output_dir(rms_output)
        
        args = ["-rms", rms_output] + surface_files
        
        return self._run_command(args)


class Surf2PialWhite(CATToolBase):
    """Wrapper for CAT_Surf2PialWhite - Create pial and white matter surfaces."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_Surf2PialWhite"
    
    def create_surfaces(self,
                       central_surface: str,
                       thickness_file: str,
                       pial_output: str,
                       white_output: str,
                       weight: float = 1.0,
                       equivolume: bool = False,
                       check_intersect: bool = False,
                       verbose: bool = False) -> CATResult:
        """
        Create pial and white matter surfaces from central surface and thickness.
        
        Args:
            central_surface: Input central surface file
            thickness_file: Thickness values file
            pial_output: Output pial surface file
            white_output: Output white matter surface file
            weight: Weight factor for thickness displacement
            equivolume: Use equi-volume model correction
            check_intersect: Check for self-intersections
            verbose: Enable verbose mode
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(central_surface, "Central surface")
        self._validate_file_exists(thickness_file, "Thickness file")
        self._validate_output_dir(pial_output)
        self._validate_output_dir(white_output)
        
        args = []
        
        if equivolume:
            args.append("-equivolume")
        if check_intersect:
            args.append("-check-intersect")
        if weight != 1.0:
            args.extend(["-weight", str(weight)])
        if verbose:
            args.append("-verbose")
        
        args.extend([central_surface, thickness_file, pial_output, white_output])
        
        return self._run_command(args)


class SurfDistance(CATToolBase):
    """Wrapper for CAT_SurfDistance - Calculate distance between surfaces."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfDistance"
    
    def calculate_distance(self,
                          surface1: str,
                          surface2: str,
                          output_values: str,
                          thickness_file: Optional[str] = None) -> CATResult:
        """
        Calculate distance between two surfaces or using thickness file.
        
        Args:
            surface1: First surface file (or central surface if using thickness)
            surface2: Second surface file (ignored if using thickness)
            output_values: Output distance values file
            thickness_file: Optional thickness file for internal surface generation
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(surface1, "First surface")
        if thickness_file:
            self._validate_file_exists(thickness_file, "Thickness file")
        else:
            self._validate_file_exists(surface2, "Second surface")
        self._validate_output_dir(output_values)
        
        args = []
        
        if thickness_file:
            args.extend(["-thickness", thickness_file, surface1, output_values])
        else:
            args.extend([surface1, surface2, output_values])
        
        return self._run_command(args)


class SurfWarp(CATToolBase):
    """Wrapper for CAT_SurfWarp - Warp surfaces using deformation fields."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfWarp"
    
    def warp_surface(self,
                    input_surface: str,
                    output_surface: str,
                    deformation_file: str) -> CATResult:
        """
        Warp surface using deformation field.
        
        Args:
            input_surface: Input surface file
            output_surface: Output warped surface file
            deformation_file: Deformation field file
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_surface, "Input surface")
        self._validate_file_exists(deformation_file, "Deformation file")
        self._validate_output_dir(output_surface)
        
        args = [input_surface, output_surface, deformation_file]
        
        return self._run_command(args)


class Surf2Sphere(CATToolBase):
    """Wrapper for CAT_Surf2Sphere - Create spherical parameterization."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_Surf2Sphere"
    
    def create_sphere(self,
                     input_surface: str,
                     output_sphere: str,
                     stop_at: Optional[float] = None,
                     verbose: bool = False) -> CATResult:
        """
        Create spherical parameterization of surface.
        
        Args:
            input_surface: Input surface file
            output_sphere: Output spherical surface file
            stop_at: Optional stopping criterion
            verbose: Enable verbose mode
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_surface, "Input surface")
        self._validate_output_dir(output_sphere)
        
        args = [input_surface, output_sphere]
        
        if stop_at is not None:
            args.append(str(stop_at))
        if verbose:
            args.append("verbose")
        
        return self._run_command(args)


class SurfDeform(CATToolBase):
    """Wrapper for CAT_SurfDeform - Deform surfaces."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfDeform"
    
    def deform_surface(self,
                      input_surface: str,
                      output_surface: str,
                      deformation_params: dict) -> CATResult:
        """
        Deform surface with specified parameters.
        
        Args:
            input_surface: Input surface file
            output_surface: Output deformed surface file
            deformation_params: Dictionary of deformation parameters
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_surface, "Input surface")
        self._validate_output_dir(output_surface)
        
        args = [input_surface, output_surface]
        
        # Add deformation parameters (implementation depends on tool specifics)
        for key, value in deformation_params.items():
            args.extend([f"-{key}", str(value)])
        
        return self._run_command(args)


class SurfCorrectThicknessFolding(CATToolBase):
    """Wrapper for CAT_SurfCorrectThicknessFolding - Correct thickness in folded regions."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfCorrectThicknessFolding"
    
    def correct_thickness(self,
                         input_surface: str,
                         thickness_file: str,
                         output_thickness: str,
                         correction_factor: float = 1.0) -> CATResult:
        """
        Correct thickness values in folded regions.
        
        Args:
            input_surface: Input surface file
            thickness_file: Input thickness values file
            output_thickness: Output corrected thickness file
            correction_factor: Correction factor for folded regions
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_surface, "Input surface")
        self._validate_file_exists(thickness_file, "Thickness file")
        self._validate_output_dir(output_thickness)
        
        args = [input_surface, thickness_file, output_thickness]
        
        if correction_factor != 1.0:
            args.extend(["-factor", str(correction_factor)])
        
        return self._run_command(args)


class SurfReduce(CATToolBase):
    """Wrapper for CAT_SurfReduce - Reduce surface mesh complexity."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfReduce"
    
    def reduce_surface(self,
                      input_surface: str,
                      output_surface: str,
                      reduction_factor: float = 0.5,
                      preserve_topology: bool = True) -> CATResult:
        """
        Reduce surface mesh complexity.
        
        Args:
            input_surface: Input surface file
            output_surface: Output reduced surface file
            reduction_factor: Factor for mesh reduction (0-1)
            preserve_topology: Preserve surface topology
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_surface, "Input surface")
        self._validate_output_dir(output_surface)
        
        args = [input_surface, output_surface, str(reduction_factor)]
        
        if not preserve_topology:
            args.append("-no-topology")
        
        return self._run_command(args)


class SurfResample(CATToolBase):
    """Wrapper for CAT_SurfResample - Resample surface data."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_SurfResample"
    
    def resample_surface(self,
                        input_surface: str,
                        template_surface: str,
                        output_surface: str,
                        values_file: Optional[str] = None,
                        output_values: Optional[str] = None,
                        interpolation: str = "linear") -> CATResult:
        """
        Resample surface to match template resolution.
        
        Args:
            input_surface: Input surface file
            template_surface: Template surface file for resampling
            output_surface: Output resampled surface file
            values_file: Optional input values file to resample
            output_values: Optional output values file
            interpolation: Interpolation method
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_surface, "Input surface")
        self._validate_file_exists(template_surface, "Template surface")
        self._validate_output_dir(output_surface)
        
        args = [input_surface, template_surface, output_surface]
        
        if values_file:
            self._validate_file_exists(values_file, "Values file")
            args.append(values_file)
            
        if output_values:
            self._validate_output_dir(output_values)
            args.append(output_values)
        
        return self._run_command(args)


# Convenience functions for surface processing
def average_surfaces(surface_files: List[str], output_file: str,
                    config: Optional[CATConfig] = None) -> CATResult:
    """Convenience function for surface averaging."""
    tool = SurfAverage(config or get_default_config())
    return tool.average_surfaces(surface_files, output_file)


def create_pial_white_surfaces(central_surface: str, thickness_file: str,
                              pial_output: str, white_output: str,
                              config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for creating pial/white surfaces."""
    tool = Surf2PialWhite(config or get_default_config())
    return tool.create_surfaces(central_surface, thickness_file, pial_output, white_output, **kwargs)


def calculate_surface_distance(surface1: str, surface2: str, output_values: str,
                              config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for surface distance calculation."""
    tool = SurfDistance(config or get_default_config())
    return tool.calculate_distance(surface1, surface2, output_values, **kwargs)


def warp_surface(input_surface: str, output_surface: str, deformation_file: str,
                config: Optional[CATConfig] = None) -> CATResult:
    """Convenience function for surface warping."""
    tool = SurfWarp(config or get_default_config())
    return tool.warp_surface(input_surface, output_surface, deformation_file)


def create_spherical_surface(input_surface: str, output_sphere: str,
                           config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for spherical parameterization."""
    tool = Surf2Sphere(config or get_default_config())
    return tool.create_sphere(input_surface, output_sphere, **kwargs)


def resample_surface(input_surface: str, template_surface: str, output_surface: str,
                    config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for surface resampling."""
    tool = SurfResample(config or get_default_config())
    return tool.resample_surface(input_surface, template_surface, output_surface, **kwargs)
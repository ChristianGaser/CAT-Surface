"""
Volume processing tools for CAT-Surface Python wrapper.

This module contains wrappers for CAT-Surface volume processing tools.
"""

from typing import Optional, List, Tuple, Union
from . import CATToolBase, CATResult, get_default_config, CATConfig


class VolAmap(CATToolBase):
    """Wrapper for CAT_VolAmap - Segmentation with adaptive MAP."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_VolAmap"
    
    def segment(self,
                input_file: str,
                label_file: str,
                output_file: Optional[str] = None,
                iters: int = 50,
                subsample: int = 96,
                iters_icm: int = 50,
                mrf_weight: float = 0.0,
                las_weight: float = 0.5,
                bias_fwhm: float = 0.0,
                h_ornlm: float = 0.05,
                sigma_ornlm: float = -1.0,
                pve: bool = True,
                cleanup: int = 2,
                write_seg: Tuple[bool, bool, bool] = (False, True, False),
                write_label: bool = True,
                write_corr: bool = True,
                write_bias: bool = False,
                use_median: bool = False,
                use_bmap: bool = False,
                verbose: bool = False) -> CATResult:
        """
        Perform segmentation with adaptive MAP.
        
        Args:
            input_file: Input NIfTI volume file
            label_file: File containing segmentation labels for initialization
            output_file: Output basename (if None, uses input filename)
            iters: Number of iterations for Amap approach
            subsample: Subsampling factor for Amap approach
            iters_icm: Number of iterations for ICM algorithm
            mrf_weight: Weight of Markov Random Field prior (0-1)
            las_weight: Weight of local adaptive segmentation (0-1)
            bias_fwhm: FWHM for bias correction smoothing kernel
            h_ornlm: Smoothing parameter for exponential function
            sigma_ornlm: Sigma for Rician noise correction (default: -1.0 = use h_ornlm)
            pve: Use Partial Volume Estimation with 5 classes
            cleanup: Clean-up level (0=none, 1=medium, 2=strong)
            write_seg: Tuple of (CSF, GM, WM) flags for writing separate images
            write_label: Write label image
            write_corr: Write nu-corrected image
            write_bias: Write bias field
            use_median: Use local mean instead of median for peaks
            use_bmap: Use BMAP instead of AMAP (experimental)
            verbose: Enable verbose mode
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_file, "Input file")
        self._validate_file_exists(label_file, "Label file")
        
        args = ["-label", label_file]
        
        # Add optional parameters
        if iters != 50:
            args.extend(["-iters", str(iters)])
        if subsample != 96:
            args.extend(["-sub", str(subsample)])
        if iters_icm != 50:
            args.extend(["-iters-icm", str(iters_icm)])
        if mrf_weight != 0.0:
            args.extend(["-mrf", str(mrf_weight)])
        if las_weight != 0.5:
            args.extend(["-las", str(las_weight)])
        if bias_fwhm != 0.0:
            args.extend(["-bias-fwhm", str(bias_fwhm)])
        if h_ornlm != 0.05:
            args.extend(["-h-ornlm", str(h_ornlm)])
        if sigma_ornlm >= 0.0:
            args.extend(["-sigma-ornlm", str(sigma_ornlm)])
        if not pve:
            args.extend(["-pve", "0"])
        if cleanup != 2:
            args.extend(["-cleanup", str(cleanup)])
        
        # Write segmentation flags
        if write_seg != (False, True, False):
            seg_values = [int(x) for x in write_seg]
            args.extend(["-write-seg"] + [str(x) for x in seg_values])
        
        # Boolean flags
        if not write_label:
            args.append("-nowrite-label")
        if not write_corr:
            args.append("-nowrite-corr")
        if write_bias:
            args.append("-write-bias")
        if use_median:
            args.append("-use-median")
        if use_bmap:
            args.append("-use-bmap")
        if verbose:
            args.append("-verbose")
        
        # Add input and output files
        args.append(input_file)
        if output_file:
            self._validate_output_dir(output_file)
            args.append(output_file)
        
        return self._run_command(args)


class VolSanlm(CATToolBase):
    """Wrapper for CAT_VolSanlm - Spatial Adaptive Non-Local Means denoising."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_VolSanlm"
    
    def denoise(self,
                input_file: str,
                output_file: str,
                rician: bool = False) -> CATResult:
        """
        Perform SANLM denoising.
        
        Args:
            input_file: Input NIfTI volume file
            output_file: Output denoised volume file
            rician: Use Rician noise model instead of Gaussian
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_file, "Input file")
        self._validate_output_dir(output_file)
        
        args = []
        if rician:
            args.append("-rician")
        
        args.extend([input_file, output_file])
        
        return self._run_command(args)


class VolThicknessPbt(CATToolBase):
    """Wrapper for CAT_VolThicknessPbt - Projection-based cortical thickness estimation."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_VolThicknessPbt"
    
    def compute_thickness(self,
                         input_file: str,
                         gmt_output: str,
                         ppm_output: str,
                         wmd_output: Optional[str] = None,
                         csd_output: Optional[str] = None,
                         gmt1_output: Optional[str] = None,
                         gmt2_output: Optional[str] = None,
                         n_avgs: int = 8,
                         fill_holes: float = 0.5,
                         downsample: float = 0.0,
                         sharpen: float = 0.0,
                         correct_voxelsize: float = 0.0,
                         verbose: bool = False) -> CATResult:
        """
        Perform projection-based cortical thickness estimation.
        
        Args:
            input_file: Input PVE label image
            gmt_output: Output GMT (thickness) file
            ppm_output: Output PPM (percentage position map) file
            wmd_output: Optional WMD output file
            csd_output: Optional CSD output file
            gmt1_output: Optional GMT1 output file
            gmt2_output: Optional GMT2 output file
            n_avgs: Number of averages for distance estimation
            fill_holes: Threshold to fill holes in PPM image
            downsample: Downsample resolution for PPM and GMT
            sharpen: Amount of sharpening for PPM map
            correct_voxelsize: Amount of thickness correction for voxel-size
            verbose: Enable verbose mode
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_file, "Input file")
        self._validate_output_dir(gmt_output)
        self._validate_output_dir(ppm_output)
        
        args = []
        
        # Add optional parameters
        if verbose:
            args.append("-verbose")
        if n_avgs != 8:
            args.extend(["-n-avgs", str(n_avgs)])
        if fill_holes != 0.5:
            args.extend(["-fill-holes", str(fill_holes)])
        if downsample != 0.0:
            args.extend(["-downsample", str(downsample)])
        if sharpen != 0.0:
            args.extend(["-sharpen", str(sharpen)])
        if correct_voxelsize != 0.0:
            args.extend(["-correct-voxelsize", str(correct_voxelsize)])
        
        # Add required files
        args.extend([input_file, gmt_output, ppm_output])
        
        # Add optional output files
        if wmd_output:
            self._validate_output_dir(wmd_output)
            args.append(wmd_output)
        if csd_output:
            self._validate_output_dir(csd_output)
            args.append(csd_output)
        if gmt1_output:
            self._validate_output_dir(gmt1_output)
            args.append(gmt1_output)
        if gmt2_output:
            self._validate_output_dir(gmt2_output)
            args.append(gmt2_output)
        
        return self._run_command(args)


class VolMarchingCubes(CATToolBase):
    """Wrapper for CAT_VolMarchingCubes - Surface extraction using marching cubes."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_VolMarchingCubes"
    
    def extract_surface(self,
                       input_file: str,
                       output_file: str,
                       threshold: float = 0.5) -> CATResult:
        """
        Extract surface using marching cubes algorithm.
        
        Args:
            input_file: Input volume file
            output_file: Output surface file
            threshold: Isosurface threshold value
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(input_file, "Input file")
        self._validate_output_dir(output_file)
        
        args = [input_file, output_file, str(threshold)]
        
        return self._run_command(args)


class Vol2Surf(CATToolBase):
    """Wrapper for CAT_Vol2Surf - Project volume data onto surface."""
    
    @property
    def executable_name(self) -> str:
        return "CAT_Vol2Surf"
    
    def project_to_surface(self,
                          volume_file: str,
                          surface_file: str,
                          output_file: str,
                          interpolation: str = "linear") -> CATResult:
        """
        Project volume data onto surface.
        
        Args:
            volume_file: Input volume file
            surface_file: Input surface file
            output_file: Output surface values file
            interpolation: Interpolation method
            
        Returns:
            CATResult with execution details
        """
        self._validate_file_exists(volume_file, "Volume file")
        self._validate_file_exists(surface_file, "Surface file")
        self._validate_output_dir(output_file)
        
        args = [volume_file, surface_file, output_file]
        
        return self._run_command(args)


# Convenience functions for volume processing
def segment_with_amap(input_file: str, label_file: str, output_file: Optional[str] = None,
                     config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for CAT_VolAmap segmentation."""
    tool = VolAmap(config or get_default_config())
    return tool.segment(input_file, label_file, output_file, **kwargs)


def denoise_sanlm(input_file: str, output_file: str, rician: bool = False,
                 config: Optional[CATConfig] = None) -> CATResult:
    """Convenience function for CAT_VolSanlm denoising."""
    tool = VolSanlm(config or get_default_config())
    return tool.denoise(input_file, output_file, rician)


def compute_thickness_pbt(input_file: str, gmt_output: str, ppm_output: str,
                         config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for CAT_VolThicknessPbt."""
    tool = VolThicknessPbt(config or get_default_config())
    return tool.compute_thickness(input_file, gmt_output, ppm_output, **kwargs)


def extract_surface_marching_cubes(input_file: str, output_file: str, threshold: float = 0.5,
                                  config: Optional[CATConfig] = None) -> CATResult:
    """Convenience function for CAT_VolMarchingCubes."""
    tool = VolMarchingCubes(config or get_default_config())
    return tool.extract_surface(input_file, output_file, threshold)


def project_volume_to_surface(volume_file: str, surface_file: str, output_file: str,
                             config: Optional[CATConfig] = None, **kwargs) -> CATResult:
    """Convenience function for CAT_Vol2Surf."""
    tool = Vol2Surf(config or get_default_config())
    return tool.project_to_surface(volume_file, surface_file, output_file, **kwargs)
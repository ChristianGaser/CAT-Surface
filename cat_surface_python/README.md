# CAT-Surface Python Wrapper

A comprehensive Python wrapper for CAT-Surface neuroimaging tools, providing a Pythonic interface to cortical surface analysis functionality.

## Overview

CAT-Surface is a collection of high-performance command-line tools for surface-based neuroimaging analysis. This Python wrapper provides:

- **Subprocess-based interface** to all major CAT-Surface tools
- **Pythonic API** with proper error handling and validation
- **Type hints** for better development experience
- **Comprehensive examples** and documentation
- **Batch processing capabilities** for multiple subjects

## Features

### Volume Processing Tools
- **CAT_VolAmap**: Segmentation with adaptive MAP
- **CAT_VolSanlm**: Spatial Adaptive Non-Local Means denoising
- **CAT_VolThicknessPbt**: Projection-based cortical thickness estimation
- **CAT_VolMarchingCubes**: Surface extraction using marching cubes
- **CAT_Vol2Surf**: Project volume data onto surfaces

### Surface Processing Tools
- **CAT_SurfAverage**: Calculate average and RMS of surfaces
- **CAT_Surf2PialWhite**: Create pial and white matter surfaces
- **CAT_SurfDistance**: Calculate distances between surfaces
- **CAT_SurfWarp**: Warp surfaces using deformation fields
- **CAT_Surf2Sphere**: Create spherical parameterizations
- **CAT_SurfDeform**: Deform surfaces
- **CAT_SurfCorrectThicknessFolding**: Correct thickness in folded regions
- **CAT_SurfReduce**: Reduce surface mesh complexity
- **CAT_SurfResample**: Resample surface data

## Installation

### Prerequisites

1. **CAT-Surface**: Ensure CAT-Surface tools are compiled and available
2. **Python 3.7+**: Required for type hints and modern Python features

### Setup

1. **Clone or download** this wrapper to your local machine
2. **Add to Python path** or install as a package:

```bash
# Option 1: Add to PYTHONPATH
export PYTHONPATH="/path/to/cat_surface_python:$PYTHONPATH"

# Option 2: Install as editable package
cd cat_surface_python
pip install -e .
```

## Quick Start

### Basic Usage

```python
import cat_surface

# Configure executable path
cat_surface.set_executable_path("/path/to/CAT-Surface/build/Progs")

# Segment brain volume
result = cat_surface.segment_with_amap(
    input_file="brain_t1.nii",
    label_file="initial_labels.nii",
    output_file="segmented",
    verbose=True
)

# Compute cortical thickness
thickness_result = cat_surface.compute_thickness_pbt(
    input_file="segmented_label.nii",
    gmt_output="thickness.nii",
    ppm_output="position_map.nii"
)

# Extract surface
surface_result = cat_surface.extract_surface_marching_cubes(
    input_file="segmented_label.nii",
    output_file="central_surface.gii",
    threshold=1.5
)
```

### Class-based Interface

```python
from cat_surface.volume import VolAmap, VolThicknessPbt
from cat_surface.surface import SurfDistance

# Use class interface for more control
segmenter = VolAmap()
result = segmenter.segment(
    input_file="brain.nii",
    label_file="labels.nii",
    iters=100,
    pve=True,
    cleanup=2
)

# Calculate surface distances
distance_tool = SurfDistance()
dist_result = distance_tool.calculate_distance(
    surface1="white.gii",
    surface2="pial.gii",
    output_values="distances.txt"
)
```

## Configuration

### Executable Path Setup

The wrapper needs to know where CAT-Surface executables are located:

```python
# Method 1: Set global path
cat_surface.set_executable_path("/path/to/executables")

# Method 2: Create custom configuration
from cat_surface import CATConfig
config = CATConfig("/path/to/executables")

# Method 3: Auto-detection (searches standard locations)
config = CATConfig()  # Will search common installation paths
```

### Standard Installation Paths

The wrapper automatically searches these locations:
- `/usr/local/bin`
- `/usr/bin`
- `/opt/CAT-Surface/bin`
- Project build directories (for development)

## API Reference

### Volume Processing

#### segment_with_amap(input_file, label_file, ...)
Segment brain volume using adaptive MAP algorithm.

**Parameters:**
- `input_file` (str): Input NIfTI volume
- `label_file` (str): Initial segmentation labels
- `output_file` (str, optional): Output basename
- `iters` (int): Number of iterations (default: 50)
- `pve` (bool): Use partial volume estimation (default: True)
- `cleanup` (int): Cleanup level 0-2 (default: 2)
- `verbose` (bool): Enable verbose output

#### compute_thickness_pbt(input_file, gmt_output, ppm_output, ...)
Compute cortical thickness using projection-based method.

**Parameters:**
- `input_file` (str): Input PVE label image
- `gmt_output` (str): Output thickness file
- `ppm_output` (str): Output position map file
- `n_avgs` (int): Number of averages (default: 8)
- `fill_holes` (float): Hole filling threshold (default: 0.5)

### Surface Processing

#### create_pial_white_surfaces(central_surface, thickness_file, ...)
Create pial and white matter surfaces from central surface.

**Parameters:**
- `central_surface` (str): Input central surface
- `thickness_file` (str): Thickness values
- `pial_output` (str): Output pial surface
- `white_output` (str): Output white surface
- `weight` (float): Displacement weight (default: 1.0)
- `equivolume` (bool): Apply equi-volume correction

#### calculate_surface_distance(surface1, surface2, output_values, ...)
Calculate distances between two surfaces.

**Parameters:**
- `surface1` (str): First surface file
- `surface2` (str): Second surface file  
- `output_values` (str): Output distance values
- `thickness_file` (str, optional): Use thickness instead of second surface

## Error Handling

The wrapper provides comprehensive error handling:

```python
from cat_surface import CATSurfaceError, CATProcessError

try:
    result = cat_surface.segment_with_amap(
        input_file="brain.nii",
        label_file="labels.nii"
    )
except FileNotFoundError as e:
    print(f"Input file not found: {e}")
except CATProcessError as e:
    print(f"CAT process failed: {e}")
    print(f"Return code: {e.returncode}")
    print(f"Error output: {e.stderr}")
except CATSurfaceError as e:
    print(f"General CAT error: {e}")
```

## Batch Processing

Process multiple subjects efficiently:

```python
subjects = ["subj001", "subj002", "subj003"]

for subject in subjects:
    try:
        # Segmentation
        seg_result = cat_surface.segment_with_amap(
            input_file=f"data/{subject}/t1.nii",
            label_file=f"data/{subject}/labels.nii",
            output_file=f"output/{subject}/segmented"
        )
        
        # Thickness estimation
        thickness_result = cat_surface.compute_thickness_pbt(
            input_file=f"output/{subject}/segmented_label.nii",
            gmt_output=f"output/{subject}/thickness.nii",
            ppm_output=f"output/{subject}/position.nii"
        )
        
        print(f"{subject} processed successfully!")
        
    except CATSurfaceError as e:
        print(f"Error processing {subject}: {e}")
        continue
```

## Examples

See `examples.py` for comprehensive usage examples including:
- Volume segmentation and thickness estimation
- Surface extraction and processing
- Batch processing workflows
- Error handling patterns
- Advanced surface analysis

## File Format Support

The wrapper supports all file formats handled by CAT-Surface:
- **NIfTI** (.nii, .nii.gz): Volume data
- **GIFTI** (.gii): Surface geometry and data
- **FreeSurfer** formats: Surface files
- **BIC objects** (.obj): Surface meshes
- **Text files**: Surface values and coordinates

## Development

### Adding New Tools

To add support for additional CAT-Surface tools:

1. **Examine the tool's CLI interface** in the source code
2. **Create a wrapper class** inheriting from `CATToolBase`
3. **Implement the `executable_name` property** and main methods
4. **Add convenience functions** for easy access
5. **Write tests and examples**

Example:
```python
class NewTool(CATToolBase):
    @property
    def executable_name(self) -> str:
        return "CAT_NewTool"
    
    def process(self, input_file: str, output_file: str) -> CATResult:
        args = [input_file, output_file]
        return self._run_command(args)
```

### Testing

Run the examples to test functionality:

```python
python examples.py
```

## Contributing

1. **Fork the repository**
2. **Create a feature branch**
3. **Add your improvements**
4. **Write tests and documentation**
5. **Submit a pull request**

## License

This wrapper follows the same license as CAT-Surface: dual-licensed under the
GNU General Public License v3 (or later) for open-source use, with a commercial
license option available. See the main [LICENSE](../LICENSE) file for details.

## Support

For issues specific to this Python wrapper:
- Check the examples and documentation
- Verify CAT-Surface executables are properly installed
- Ensure file paths and permissions are correct

For CAT-Surface tool issues:
- Refer to the original CAT-Surface documentation
- Contact the CAT-Surface developers

## Acknowledgments

- **Christian Gaser** and the CAT-Surface development team
- **Original CAT-Surface project**: https://github.com/ChristianGaser/CAT-Surface
- **CAT12 toolbox**: https://github.com/ChristianGaser/cat12

---

**Note**: This is a subprocess wrapper that calls the original CAT-Surface command-line tools. Ensure CAT-Surface is properly compiled and accessible before using this wrapper.
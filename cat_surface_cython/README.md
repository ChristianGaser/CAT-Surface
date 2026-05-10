# cat-surf — Python bindings for CAT-Surface

[![PyPI](https://img.shields.io/pypi/v/cat-surf)](https://pypi.org/project/cat-surf/)
[![Python](https://img.shields.io/pypi/pyversions/cat-surf)](https://pypi.org/project/cat-surf/)
[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](https://github.com/ChristianGaser/CAT-Surface/blob/master/LICENSE)

`cat-surf` provides Python access to [CAT-Surface](https://github.com/ChristianGaser/CAT-Surface),
a mature C/C++ toolkit for surface-based neuroimaging analysis, focusing on the processing and analysis of cortical surface meshes.

CAT-Surface has been used internally for more than 10 years as part of
[CAT12](https://github.com/ChristianGaser/cat12), the SPM toolbox for computational anatomy.
With this package, it is now distributed as an independent Python package for direct integration into external workflows to make CAT-Surface easier to use:

- in **T1Prep** workflows (currently relying on CAT-Surface binaries)
- in **Python-based** pipelines without subprocess-heavy wrappers
- in reproducible environments via published wheels

---

## Installation

```bash
pip install cat-surf
```

Pre-built wheels are available for:

- **macOS** — arm64 (Apple Silicon), x86_64 (Intel)
- **Linux** — x86_64, aarch64 (manylinux)
- **Python** 3.9 – 3.13

---

## Basic usage

```python
import cat_surf

# Check version
print(cat_surf.__version__)

# Load a surface file (GIFTI, FreeSurfer, BIC/MNI formats)
vertices, faces = cat_surf.read_surface("lh.central.gii")

# Per-vertex area
area, total_area = cat_surf.get_area(vertices, faces)

# Euler characteristic
chi = cat_surf.euler_characteristic(vertices, faces)

# Smooth per-vertex data (heat kernel, FWHM in mm)
smoothed = cat_surf.smooth_heatkernel(vertices, faces, area, fwhm=20.0)
```

### Surface operations

| Function | Description |
|---|---|
| `read_surface` / `write_surface` | Multi-format surface I/O (GIFTI, FreeSurfer, BIC) |
| `read_values` / `write_values` | Per-vertex scalar I/O |
| `get_area` | Per-vertex and total surface area |
| `get_area_normalized` | Area normalized by a reference sphere |
| `euler_characteristic` | Topological integrity check |
| `smooth_heatkernel` | Heat-kernel smoothing of per-vertex data |
| `smooth_mesh` | Laplacian/Taubin mesh vertex smoothing |
| `smoothed_curvatures` | Mean curvature estimation |
| `sulcus_depth` | Sulcal depth via depth potential |
| `reduce_mesh` | Quadric mesh decimation |
| `remove_intersections` | Self-intersection repair |
| `count_intersections` | Count self-intersecting faces |
| `surf_average` | Vertex-wise averaging across surfaces |
| `surf_to_sphere` | Inflate surface to sphere |
| `sphere_radius` | Mean radius of a spherical surface |
| `correct_thickness_folding` | Folding-based thickness correction |
| `point_distance` | Vertex-to-surface point distances (symmetric) |
| `point_distance_mean` | Mean vertex-to-surface distance |
| `hausdorff_distance` | Hausdorff distance between two surfaces |

### Volume operations

| Function | Description |
|---|---|
| `vol_sanlm` | Structure-adaptive non-local means denoising |
| `vol_blood_vessel_correction` | Blood vessel intensity correction |
| `vol_thickness_pbt` | Cortical thickness via projection-based method |
| `vol_amap` | Adaptive maximum a posteriori tissue segmentation |
| `vol_marching_cubes` | Isosurface extraction from a volume file |

### Registration

| Function | Description |
|---|---|
| `bbreg` | Full BBR pipeline: optional NMI init → boundary-based surface registration |
| `bbreg_detect_contrast` | Auto-detect T1/FLAIR vs T2/BOLD from WM/GM intensity ratio |
| `volume_register_nmi` | Cross-modal rigid registration via Normalised Mutual Information (≈ `mri_coreg`) |
| `volume_register_robust` | Same-modality rigid registration via Tukey biweight M-estimation (≈ `mri_robust_register`) |

#### Basic BBR usage

```python
import cat_surf
import nibabel as nib

# Load surfaces (GIFTI or any format supported by cat_surf.read_surface)
lh_verts, lh_faces = cat_surf.read_surface("lh.white.surf.gii")
rh_verts, rh_faces = cat_surf.read_surface("rh.white.surf.gii")

# Full pipeline: NMI init from T1w reference, then BBR
matrix, cost = cat_surf.bbreg(
    "bold_mean.nii.gz",
    lh_surface=(lh_verts, lh_faces),
    rh_surface=(rh_verts, rh_faces),
    ref_file="T1w.nii.gz",     # NMI initialisation
    verbose=True,
)
print(f"BBR cost: {cost:.4f}")
print("EPI → T1 matrix:\n", matrix)

# Save the 4×4 transform for use with FSL / ANTs
import numpy as np
np.savetxt("epi_to_t1.txt", matrix)
```

#### Standalone volume registration

```python
# Cross-modal (EPI ↔ T1w) — NMI
matrix, nmi = cat_surf.volume_register_nmi("T1w.nii.gz", "bold_mean.nii.gz")

# Same-modality (T1w ↔ T1w) — robust IRLS
matrix, res = cat_surf.volume_register_robust("t1_ref.nii.gz", "t1_moving.nii.gz")
```

#### Contrast auto-detection

```python
# 0 = T1/FLAIR, 1 = T2/BOLD, -1 = undetermined
contrast = cat_surf.bbreg_detect_contrast(
    "bold_mean.nii.gz",
    lh_surface=(lh_verts, lh_faces),
    rh_surface=(rh_verts, rh_faces),
)
```

### Conversion utilities

| Function | Description |
|---|---|
| `arrays_to_polygons` | Convert NumPy vertex/face arrays to internal polygon mesh |
| `polygons_to_arrays` | Convert internal polygon mesh back to NumPy arrays |

---

## Typical use cases

- Cortical mesh processing (resampling, smoothing, metrics)
- Thickness and folding related computations
- Volume-to-surface projection
- Denoising and volume preprocessing for structural MRI

---

## Citation / provenance

If you use this package in research, please cite [Dahnke et al., 2013](https://doi.org/10.1016/j.neuroimage.2012.09.050) and mention the `cat-surf` package version for reproducibility.

---

## Source

- Source: <https://github.com/ChristianGaser/CAT-Surface>

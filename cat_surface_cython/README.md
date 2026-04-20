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
| `get_area` | Per-vertex and total surface area |
| `euler_characteristic` | Topological integrity check |
| `smooth_heatkernel` | Heat-kernel smoothing of per-vertex data |
| `smoothed_curvatures` | Mean curvature estimation |
| `sulcus_depth` | Sulcal depth via depth potential |
| `reduce_mesh` | Quadric mesh decimation |
| `remove_intersections` | Self-intersection repair |
| `surf_average` | Vertex-wise averaging across surfaces |
| `correct_thickness_folding` | Folding-based thickness correction |

### Volume operations

| Function | Description |
|---|---|
| `read_values` / `write_values` | Per-vertex scalar I/O |
| `sanlm_filter` | Spatial-adaptive non-local means denoising |

---

## Typical use cases

- Cortical mesh processing (resampling, smoothing, metrics)
- Thickness and folding related computations
- Volume-to-surface projection
- Denoising and volume preprocessing for structural MRI

---

## Citation / provenance

If you use this package in research, please cite the CAT12/CAT-Surface methodology
used in your workflow and mention the `cat-surf` package version for reproducibility.

---

## Source

- Source: <https://github.com/ChristianGaser/CAT-Surface>

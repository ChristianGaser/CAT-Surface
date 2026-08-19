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

### Building from source

The extension modules are generated from the `.pyx` sources at build time;
the intermediate `.c` files are build artifacts and are not tracked in the
repository. Cython is declared in `build-system.requires`, so any PEP 517
build installs it for you:

```bash
# build libCAT first (see the CAT-Surface README), then:
pip install ./cat_surface_cython
```

For an in-place development build, Cython has to be in the ambient
environment and libCAT has to be findable:

```bash
pip install cython
CAT_BUILD_DIR=/path/to/CAT-Surface python setup.py build_ext --inplace
```

Note that this regenerates every `.c` file, so a Cython version that differs
from the one used elsewhere changes all of them — harmless now that they are
untracked, but it is why they are not kept in git.

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

| Function | Description | Mirrors |
| --- | --- | --- |
| `read_surface` / `write_surface` | Multi-format surface I/O (GIFTI, FreeSurfer, BIC) | — |
| `read_values` / `write_values` | Per-vertex scalar I/O | — |
| `get_area` | Per-vertex and total surface area | `CAT_SurfArea` |
| `get_area_normalized` | Area normalized by a reference sphere | `CAT_SurfArea -sphere` |
| `euler_characteristic` | Topological integrity check | — |
| `smooth_heatkernel` | Heat-kernel smoothing of per-vertex data | — |
| `smooth_mesh` | Laplacian/Taubin mesh vertex smoothing | — |
| `smoothed_curvatures` | Mean curvature estimation | — |
| `sulcus_depth` | Sulcal depth via depth potential | — |
| `reduce_mesh` | Quadric (QEM) mesh decimation (`aggressiveness`, `preserve_sharp`) | `CAT_SurfReduce` |
| `fix_self_intersect` | Self-intersection repair (topology preserving) | `CAT_SurfFixSelfIntersect` |
| `count_intersections` | Count intersecting triangle pairs | `CAT_SurfSelfIntersect` |
| `surf_average` | Vertex-wise averaging across surfaces (`return_rms` for std dev) | `CAT_SurfAverage` |
| `surf_to_sphere` | Inflate surface to sphere | `CAT_Surf2Sphere` |
| `sphere_radius` | Mean radius of a spherical surface | — |
| `correct_thickness_folding` | Folding-based thickness correction (`slope`, `max_dist` clipping) | `CAT_SurfCorrectThicknessFolding` |
| `point_distance` | Linked vertex distance between two meshes | `CAT_SurfDistance -link` |
| `point_distance_mean` | Mean (Tfs) vertex distance | `CAT_SurfDistance -mean` |
| `hausdorff_distance` | Hausdorff distance between two surfaces | — |
| `resample_to_sphere` | Resample surface/values onto a target sphere | `CAT_SurfResample` |
| `surf_deform` | Deform a surface toward a volume isovalue | `CAT_SurfDeform` |
| `surf_to_pial_white` | Estimate pial + white surfaces from a central surface (`remove_intersect` repairs both) | `CAT_Surf2PialWhite` |
| `central_to_pial` | Generate a pial surface from central + thickness | — |
| `surf_warp` | DARTEL-based spherical registration | `CAT_SurfWarp` |
| `spherical_demon` | Spherical Demons spherical registration | `CAT_SurfSphericalDemon` |

### Volume operations

| Function | Description | Mirrors |
| --- | --- | --- |
| `vol_sanlm` | Structure-adaptive non-local means denoising | `CAT_VolSanlm` |
| `vol_blood_vessel_correction` | Blood vessel intensity correction | *(no CLI; the correction `CAT_VolThicknessPbt` applies unless `-no-bvc`)* |
| `vol_thickness_pbt` | Cortical thickness via projection-based method | `CAT_VolThicknessPbt` |
| `vol_thickness_qc` | Shape triage of implausibly thick cortex | `CAT_VolThicknessQC` |
| `vol_amap` | Adaptive maximum a posteriori tissue segmentation | `CAT_VolAmap` (core only) |
| `vol_marching_cubes` | Isosurface extraction with genus-0 topology correction | `CAT_VolMarchingCubes` |
| `vol_smooth` | Isotropic Gaussian volume smoothing | `CAT_VolSmooth` |
| `vol_sheetness` | Multi-scale Hessian sheetness (plate) filter | `CAT_VolSheetness` |
| `vol_oriented_median` | Median over a sheetness-oriented neighbourhood | `CAT_VolLocalStat -oriented` |
| `vol_sulcus_repair` | Pre-PBT repair of glued sulci and lost WM blades | `CAT_VolSulcusRepair` |
| `vol2surf` | Map a volume to a surface along inward normals | `CAT_Vol2Surf` |

#### Thin structures and the shrinking bias

Every isotropic regularizer — a local median, a Potts MRF, total variation —
penalizes boundary area, and a thin structure has an extreme area-to-volume
ratio, so deleting it is always the cheaper labelling. That is why one and the
same median filter opens glued sulci in one place and closes cerebellar fissures
in another, and why tuning its strength only trades one failure for the other.

`vol_sheetness` replaces the smoothness prior with a *shape* prior read off the
Hessian eigenvalues, so it keeps thin sheets, ignores blobs, and shrinks nothing.
Everything below is built on it, and all of it is a **no-op where the sheetness is
zero** — the oriented operators are numerically identical to the isotropic ones
they replace away from thin structures, so they are safe to switch on:

```python
import cat_surf, nibabel as nib

t1  = nib.load("t1_corr.nii").get_fdata().astype("float32")
lab = nib.load("label.nii").get_fdata().astype("float32")
vx  = nib.load("label.nii").header.get_zooms()[:3]

# Inspect the evidence first: dark sheets are sulcal CSF, bright ones WM blades
sheet, normal = cat_surf.vol_sheetness(t1, voxelsize=vx, polarity=-1,
                                       return_normal=True)

# A median that cannot close a sulcus (orientation taken from the intensity)
clean = cat_surf.vol_oriented_median(lab, guide=t1, voxelsize=vx)

# Recover sulcal CSF the classifier discarded, then run PBT on the repaired map
lab = cat_surf.vol_sulcus_repair(t1, lab, voxelsize=vx, verbose=True)
gmt, ppm, _, _ = cat_surf.vol_thickness_pbt(lab, voxelsize=vx,
                                            oriented_filter=True)

```

### Volume input convention

Every volume-consuming function accepts **three input forms** interchangeably:

```python
import cat_surf, nibabel as nib
import numpy as np

# 1. file path
v, f = cat_surf.vol_marching_cubes("brain.nii.gz", threshold=0.5)

# 2. (ndarray, affine) tuple — supply a 4×4 RAS+ affine
v, f = cat_surf.vol_marching_cubes((array, affine), threshold=0.5)

# 3. any nibabel-image-like object (.affine + .get_fdata())
v, f = cat_surf.vol_marching_cubes(nib.load("brain.nii.gz"), threshold=0.5)
```

The data is auto-converted to float32 (matching the C library's expectation),
and 4-D series default to the middle frame. No copy is made if the array is
already float32/Fortran-order. Internally a minimal `nifti_image` is built
from the affine — only the fields libCAT actually reads (`dims`, `voxel
sizes`, `sto_xyz`).

### Registration

| Function | Description |
| --- | --- |
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
| --- | --- |
| `arrays_to_polygons` | Convert NumPy vertex/face arrays to internal polygon mesh |
| `polygons_to_arrays` | Convert internal polygon mesh back to NumPy arrays |

---

## `cat_surf.cli` — drop-in replacement for the CAT binaries

The `cat_surf.cli` subpackage mirrors the `CAT_*` command-line binaries
one-to-one: same names (snake_case, `CAT_` prefix dropped), same
positional argument order, same option semantics and defaults. Each
function reads its inputs from disk, calls the in-memory numpy wrapper,
and writes the outputs — ideal for porting shell scripts to Python.

```python
from cat_surf import cli

# CAT_VolMarchingCubes brain.nii.gz brain.gii -thresh 0.5
cli.vol_marching_cubes("brain.nii.gz", "brain.gii", threshold=0.5)

# CAT_Surf2PialWhite -remove_intersect central.gii thickness.txt labels.nii \
#     pial.gii white.gii
cli.surf2pial_white("central.gii", "thickness.txt", "labels.nii",
                    "pial.gii", "white.gii", remove_intersect=True)

# CAT_SurfDistance -mean surf1.gii surf2.gii out.txt
cli.surf_distance("surf1.gii", "surf2.gii", "out.txt", mode="mean")

# CAT_VolSanlm in.nii out.nii -strength 1.0
cli.vol_sanlm("in.nii", "out.nii", strength=1.0)
```

The full mapping:

| CAT binary | `cat_surf.cli.<name>` |
| --- | --- |
| `CAT_Surf2PialWhite` | `surf2pial_white` |
| `CAT_Surf2Sphere` | `surf2sphere` |
| `CAT_SurfArea` | `surf_area` |
| `CAT_SurfAverage` | `surf_average` |
| `CAT_SurfCorrectThicknessFolding` | `surf_correct_thickness_folding` |
| `CAT_SurfDeform` | `surf_deform` |
| `CAT_SurfDistance` | `surf_distance` |
| `CAT_SurfReduce` | `surf_reduce` |
| `CAT_SurfFixSelfIntersect` | `surf_fix_self_intersect` |
| `CAT_SurfResample` | `surf_resample` |
| `CAT_SurfSphericalDemon` | `surf_spherical_demon` |
| `CAT_SurfWarp` | `surf_warp` |
| `CAT_Vol2Surf` | `vol2surf` |
| `CAT_VolAmap` | `vol_amap` |
| `CAT_VolLocalStat` | `vol_local_stat` (`-oriented` only) |
| `CAT_VolMarchingCubes` | `vol_marching_cubes` |
| `CAT_VolThicknessQC` | `vol_thickness_qc` |
| `CAT_VolSanlm` | `vol_sanlm` |
| `CAT_VolSheetness` | `vol_sheetness` |
| `CAT_VolSmooth` | `vol_smooth` |
| `CAT_VolSulcusRepair` | `vol_sulcus_repair` |
| `CAT_VolThicknessPbt` | `vol_thickness_pbt` |

For composable in-memory pipelines, prefer the lower-level `cat_surf`
API directly — the CLI shims are just thin convenience wrappers.

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

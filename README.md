# CAT: Cortex Analysis Tools
Christian Gaser christian.gaser@uni-jena.de Jena University Hospital, Germany.

CAT-Surface provides high-performance command-line tools for surface-based neuroimaging analysis, focusing on the processing and analysis of cortical surface meshes.
These tools are integral to the [CAT toolbox](https://github.com/ChristianGaser/cat12) for SPM, and serve as the backend for most of CAT’s surface-based processing, including extraction, smoothing, morphometric analysis, and format conversion. The tools are written in portable C and optimized for efficient batch processing of large cohorts in neuroimaging studies.

The repository also contains additional third-party libraries (see `3rdparty/`) for mesh handling and image I/O, bundled to simplify compilation and avoid dependency issues. These packages are not maintained within this project but are bundled for convenience so that no extra downloads are needed during compilation.


## Build Instructions

**1. Build the tools:**
```bash
./autogen.sh       # Generate the configure script
./configure        # Configure the build (add --prefix if needed)
make               # Build all tools
make install       # Optional: install system-wide
```

## API Documentation (Doxygen)

CAT-Surface includes a Doxygen configuration to generate API documentation from
the C/C++ sources and public headers.

Install Doxygen:

- macOS: `brew install doxygen`
- Linux: `sudo apt-get install doxygen`

Generate docs:

```bash
./scripts/generate_doxygen_docs.sh
```

The generated HTML entry point is:

- `docs/doxygen/html/index.html`

Warnings are written to:

- `docs/doxygen/warnings.log`

---
## Tools and Functions

Below is a summary of the available command-line programs in CAT-Surface, each designed for a distinct step in the cortical surface processing pipeline.

| Tool                        | Description |
|-----------------------------|-------------|
| **CAT_Vol2Surf**                | Projects values from a 3D image (volume) onto the cortical surface mesh vertices. |
| **CAT_VolAmap**                 | Performs adaptive maximum a posteriori tissue classification/segmentation on volumetric MRI data. `-mrf-aniso` relaxes the MRF prior along thin sheets so it stops erasing them (see below). |
| **CAT_VolCalc**                 | Voxel-wise image calculator in the spirit of SPM's `spm_imcalc`: evaluates a formula over one or more co-registered volumes. In matrix mode the variable `X` stands for the vector of all inputs at a voxel, so `mean(X)`, `median(X)` and `std(X)` reduce across any number of images. |
| **CAT_VolBloodVesselCorrection**| Reduces blood-vessel-related intensity artifacts in volumetric scans before downstream analysis. |
| **CAT_VolLocalStat**            | Applies a local statistic (mean, min, max, std, median, grey open/close) over a voxel neighbourhood. `-oriented` runs the median over a sheetness-oriented neighbourhood so it cannot close a thin structure (see below). |
| **CAT_VolMarchingCubes**        | Extracts a surface mesh from volumetric data using a marching cubes isosurface algorithm. |
| **CAT_VolSanlm**                | Applies spatially adaptive non-local means denoising to volumetric MRI images. |
| **CAT_VolSheetness**            | Multi-scale Hessian sheetness (plate) filter: detects thin sheet-like structures — sulcal CSF, gyral white-matter blades — and ignores blobs (see below). |
| **CAT_VolSmooth**               | Smooths a volume with an isotropic Gaussian kernel. `-oriented` instead diffuses *along* thin sheets and not across them (see below). |
| **CAT_VolSulcusRepair**         | Anatomy-aware repair of a PVE label map before thickness estimation: opens glued sulci and reconnects gyri broken by small missegmentations, using the intensity image as evidence (see below). |
| **CAT_VolThicknessPbt**         | Estimates cortical thickness from volumetric tissue maps using a projection-based thickness method. `-oriented-filter` replaces its internal isotropic medians with sheetness-oriented ones (see below). |
| **CAT_SurfApplyWarp**           | Applies deformation fields (from CAT_ApplySurf) to transform surface meshes. |
| **CAT_SurfApplyWarpValues**     | Applies surface deformations to vertex-wise data arrays (e.g., morphometric parameters). |
| **CAT_SurfSmooth**              | Performs heat kernel smoothing on surface meshes or vertex-wise data, using an exact spectral method. |
| **CAT_SurfDeform**              | Legacy algorithm (after David MacDonald) for cortical surface extraction from 3D anatomical images. Mainly for archival/historical purposes. |
| **CAT_SurfCurvature**           | Extracts folding-related surface parameters (e.g., mean curvature, Gaussian curvature, sulcal depth) and optionally smooths results using the diffusion heat kernel. |
| **CAT_SurfCorrectThicknessFolding** | Corrects cortical thickness values for folding-related variation (optionally using linear thickness-dependent weighting via `-slope`). |
| **CAT_DumpCurv_ui**             | GUI batch interface for CAT_DumpCurv to process large collections of surface files. |
| **CAT_SurfArea**                | Calculates total and/or local surface area metrics from cortical meshes. |
| **CAT_SurfAverage**             | Computes an average surface from multiple aligned input meshes. |
| **CAT_SurfConvert**             | Converts between surface mesh file formats: BIC (obj), Freesurfer, and OOGL. |
| **CAT_SurfDistance**            | Computes pointwise and/or aggregate geometric distances between surface meshes. |
| **CAT_SurfMeasure2Txt**         | Converts Freesurfer curvature (`.curv`) files to plain text format. |
| **CAT_SurfPlotValuesAtMaximum** | For a given reference, extracts or plots the values at the vertex of maximum value for each surface file. |
| **CAT_SurfPlotValuesAtPoint**   | Extracts or plots the values at specified (x, y, z) coordinates for each input surface. |
| **CAT_SurfReduce**              | Simplifies a surface mesh by reducing vertex/triangle count while preserving geometry. |
| **CAT_SurfRemoveIntersections** | Detects and repairs self-intersections in surface meshes. |
| **CAT_SurfResample**            | Resamples surface geometry or vertex-wise values onto a target mesh/grid resolution. |
| **CAT_SurfResampleSpherical**   | Resamples a spherical surface mesh onto a standard sphere (for surface-based morphometry or group analysis). |
| **CAT_SurfSheet2Surf**          | Maps a 2D image (PGM format) onto a surface mesh. |
| **CAT_Surf2PialWhite**          | Derives pial and white matter surfaces from a central cortical surface representation. |
| **CAT_Surf2Sheet**              | Flattens surface data (e.g., curvature, morphometry) onto a 2D sheet (PGM image), for visualization or further analysis. |
| **CAT_Surf2Sphere**             | Inflates a cortical surface mesh onto a sphere using the Caret/Van Essen inflation approach. |
| **CAT_SurfWarp**                | Warps one surface to another using non-linear surface-based registration. |
| **CAT_SurfBBReg**               | Boundary-Based Registration (BBR): rigid co-registration of a functional volume to cortical surfaces. Includes NMI-based volume initialisation, automatic T1/T2 contrast detection, and optional pre-smoothing. |

### Thin structures and the sheetness family

Several of the tools above share one idea, so they are described together.

Every isotropic regularizer — a local median, a Markov random field, total
variation — penalizes boundary area. A thin structure has an extreme
area-to-volume ratio, so removing it is always the cheaper labelling. This
*shrinking bias* is why one and the same median filter opens a glued sulcus in
one place and closes a cerebellar fissure in another, and why turning its
strength up or down only trades one failure for the other.

**CAT_VolSheetness** replaces the smoothness prior with a *shape* prior taken
from the Hessian eigenvalues |l1| <= |l2| <= |l3|:

    sheet (plate):  |l1| ~ |l2| ~ 0,  |l3| large
    tube  (line):   |l1| ~ 0,         |l2| ~ |l3| large
    blob:           |l1| ~ |l2| ~ |l3| large

It therefore keeps thin sheets, ignores blobs, and shrinks nothing, because it
makes no statement about boundary length. Sulcal CSF is a dark sheet in a T1
(`-polarity -1`), a gyral white-matter blade is a bright sheet (`-polarity 1`);
the same operator finds both. Run it on its own to check the scale range and
polarity on a new protocol before switching on anything that uses it.

The field is consumed by:

| Tool | Option |
|------|--------|
| **CAT_VolLocalStat** | `-oriented` (with `-stat 7`) — median over a sheet-oriented neighbourhood |
| **CAT_VolSmooth** | `-oriented` — coherence-enhancing smoothing along the sheet |
| **CAT_VolAmap** | `-mrf-aniso 1\|2` — locally varying beta, or direction-weighted Potts |
| **CAT_VolThicknessPbt** | `-oriented-filter` — replaces the three isotropic medians inside PBT |
| **CAT_VolSulcusRepair** | always — the evidence term for the pre-PBT repair |

Every one of these is **a no-op where no sheet is detected**: the oriented
operator is then numerically identical to the isotropic one it replaces, which
is what makes them safe to enable. All are off by default.

**CAT_VolSulcusRepair** addresses three failure modes that break central-surface
extraction, all of which are failures of *evidence* rather than of smoothness —
which is why no local filter fixes them:

1. **Glued sulci.** Two banks of a tight sulcus end up as one thick grey-matter
   band because no CSF was detected between them. Typical in the occipital
   midline, where cortex is thin and contrast poorest.
2. **Broken gyri.** A thin gyral white-matter blade is interrupted by a small
   missegmentation, which corrupts the distance map and with it the thickness
   and the central surface along the whole blade.
3. **Residual partial-volume error** where (1) happens: the label map reports no
   CSF nearby while the intensity image still shows a dip across the sulcus.

Regularization cannot repair any of these, because it cannot create evidence —
it only redistributes what the classifier already committed to. Each step goes
back to the intensity image instead and recovers evidence the classifier
discarded. Every operation is one-sided and gated on several independent pieces
of evidence, so none of them can cause the failure mode it is meant to fix.

A typical pre-PBT sequence:

```bash
# Inspect the evidence first (dark sheets = sulcal CSF)
CAT_VolSheetness -polarity -1 -v t1_corr.nii sheetness.nii

# Repair the label map, then estimate thickness from the repaired map
CAT_VolSulcusRepair -verbose t1_corr.nii label.nii label_repaired.nii
CAT_VolThicknessPbt -oriented-filter label_repaired.nii gmt.nii ppm.nii
```

References: Han et al., *Proc SPIE Med Imag* 4322:194–203, 2001 (ACE); Han et
al., *NeuroImage* 23(3):997–1012, 2004 (CRUISE); Kim et al., *NeuroImage*
27(1):210–221, 2005 (CLASP); Descoteaux et al., *Med Image Anal*
10(4):638–651, 2006 (sheetness); Sato et al., *Med Image Anal* 2(2):143–168,
1998.

### External binary GIFTI files

Specifying an output name with the `.dat` extension causes CAT-Surface to write
a `.gii` header alongside a `.dat` binary. The header uses the GIFTI
`ExternalFileBinary` encoding so the resulting pair is compatible with SPM12.

## Continuous Integration

A GitHub Actions workflow located in `.github/workflows/ci.yml` automatically
runs `./autogen.sh`, `./configure`, `make` and `make check` on every push or
pull request. It then builds the Python bindings against that tree and runs
`cat_surface_cython/tests/smoke_test.py`, which imports every extension module
and re-checks the numeric contracts — the generated `cat_surf/_*.c` sources are
untracked build artifacts, so only an actual build proves the Cython sources
still compile and link.

## Python bindings for CAT-Surface

[`cat-surf`](cat_surface_cython/README.md) provides Python access to CAT-Surface.

## License

CAT-Surface is dual-licensed:

1. **GNU General Public License v3** (or later) — for open-source use
2. **Commercial License** — for proprietary applications

See the [LICENSE](LICENSE) file for full details. For commercial licensing
inquiries, contact christian.gaser@uni-jena.de.


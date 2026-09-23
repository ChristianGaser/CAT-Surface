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
| **CAT_VolAmap**                 | Performs adaptive maximum a posteriori tissue classification/segmentation on volumetric MRI data. |
| **CAT_VolCalc**                 | Voxel-wise image calculator in the spirit of SPM's `spm_imcalc`: evaluates a formula over one or more co-registered volumes. In matrix mode the variable `X` stands for the vector of all inputs at a voxel, so `mean(X)`, `median(X)` and `std(X)` reduce across any number of images. |
| **CAT_VolLocalStat**            | Applies a local statistic (mean, min, max, std, median, grey open/close) over a voxel neighbourhood. `-oriented` runs the median over a sheetness-oriented neighbourhood so it cannot close a thin structure (see below). |
| **CAT_VolMarchingCubes**        | Extracts a surface mesh from volumetric data using a marching cubes isosurface algorithm. |
| **CAT_VolSanlm**                | Applies spatially adaptive non-local means denoising to volumetric MRI images. |
| **CAT_VolSheetness**            | Multi-scale Hessian sheetness (plate) filter: detects thin sheet-like structures — sulcal CSF, gyral white-matter blades — and ignores blobs (see below). |
| **CAT_VolSmooth**               | Smooths a volume with an isotropic Gaussian kernel. |
| **CAT_VolThicknessPbt**         | Estimates cortical thickness from volumetric tissue maps using a projection-based thickness method. Its internal median filters are sheetness-oriented, so they cannot close a thin structure (see below). |
| **CAT_SurfSmooth**              | Performs heat kernel smoothing on surface meshes or vertex-wise data, using an exact spectral method. |
| **CAT_SurfDeform**              | Deforms a surface mesh towards an isovalue of a volume, e.g. the central surface onto the 0.5 level of the PPM (used by T1Prep to refine the marching-cubes surface). |
| **CAT_SurfCurvature**           | Extracts folding-related surface parameters (e.g., mean curvature, Gaussian curvature, sulcal depth) and optionally smooths results using the diffusion heat kernel. |
| **CAT_SurfCorrectThicknessFolding** | Corrects cortical thickness values for folding-related variation (optionally using linear thickness-dependent weighting via `-slope`). |
| **CAT_SurfArea**                | Calculates total and/or local surface area metrics from cortical meshes. |
| **CAT_SurfAverage**             | Computes an average surface from multiple aligned input meshes. |
| **CAT_SurfConvert**             | Converts between surface mesh file formats: BIC (obj), Freesurfer, and OOGL. |
| **CAT_SurfDistance**            | Computes pointwise and/or aggregate geometric distances between surface meshes. |
| **CAT_SurfInfo**                | Prints a summary of a surface: vertices, faces and edges, surface area, enclosed volume, bounding box, Euler number, genus, connected components and self-intersections. For a GIFTI file it also lists the DataArrays embedded next to the mesh (thickness, curvature, labels) with their value range and NaN count. `-tab` gives the same numbers as key/value pairs. |
| **CAT_SurfMeasure2Txt**         | Converts Freesurfer curvature (`.curv`) files to plain text format. |
| **CAT_SurfPlotValuesAtMaximum** | For a given reference, extracts or plots the values at the vertex of maximum value for each surface file. |
| **CAT_SurfPlotValuesAtPoint**   | Extracts or plots the values at specified (x, y, z) coordinates for each input surface. |
| **CAT_SurfReduce**              | Simplifies a surface mesh by reducing vertex/triangle count while preserving geometry. |
| **CAT_SurfResample**            | Resamples surface geometry or vertex-wise values onto a target mesh/grid resolution. |
| **CAT_SurfResampleSpherical**   | Resamples a spherical surface mesh onto a standard sphere (for surface-based morphometry or group analysis). |
| **CAT_SurfSheet2Surf**          | Maps a 2D image (PGM format) onto a surface mesh. |
| **CAT_Surf2PialWhite**          | Derives pial and white matter surfaces from a central cortical surface representation. |
| **CAT_Surf2Sheet**              | Flattens surface data (e.g., curvature, morphometry) onto a 2D sheet (PGM image), for visualization or further analysis. |
| **CAT_Surf2Sphere**             | Inflates a cortical surface mesh onto a sphere using the Caret/Van Essen inflation approach. |
| **CAT_SurfBBReg**               | Boundary-Based Registration (BBR): rigid co-registration of a functional volume to cortical surfaces. Includes NMI-based volume initialisation, automatic T1/T2 contrast detection, and optional pre-smoothing. |
| **CAT_SurfSphericalDemon**      | Spherical Demons registration of a surface to a template. |
| **CAT_SurfFixSelfIntersect**    | Removes self-intersections by locally smoothing the intersecting regions; `-reference` retreats towards the surface a deformation started from where smoothing cannot separate crossed sheets. |
| **CAT_SurfSelfIntersect**       | Counts self-intersections and marks the intersecting triangles. |
| **CAT_SurfFixTopology**         | Corrects the topology of a brain surface. |
| **CAT_SurfMarkDefects**         | Locates and marks topological defects using a spherical mapping (text output). |
| **CAT_SurfArtifacts**           | Locates and marks artifacts by comparing a surface with a smoothed version of it (text output). |
| **CAT_SurfCentral2Pial**        | Estimates the pial or white surface from the central surface and thickness, optionally with an equi-volume correction. |
| **CAT_SurfSulcusDepth**         | Sulcal depth as the Euclidean distance between the central surface and its convex hull. |
| **CAT_SurfDepthPotential**      | Depth potential of a surface. |
| **CAT_SurfConvexity**           | Convexity values of a surface. |
| **CAT_SurfSharpness**           | Sharpness values of a surface. |
| **CAT_SurfFractalDimension**    | Local fractal dimension of a surface. |
| **CAT_SurfRatio**               | Surface ratio after Toro et al. (2008), normalized for surface area so that it is scaling invariant. |
| **CAT_Surf2ConvexHull**         | Extracts the convex hull of a surface. |
| **CAT_SurfHausdorff**           | Point-by-point Hausdorff distance between two surfaces. |
| **CAT_SurfAreaDistortion**      | Area distortion between two surfaces. |
| **CAT_SurfAngularDistortion**   | Angular distortion between two surfaces. |
| **CAT_SurfMetricDistortion**    | Metric distortion between a surface and its spherical map. |
| **CAT_SurfIsometize**           | Adjusts a spherical map to be more isometric (area preserving). |
| **CAT_SurfRefine**              | Refines a mesh to a maximum edge length, optionally shorter where the absolute mean curvature is large. |
| **CAT_SurfSeparatePolygons**    | Splits a mesh into its disjoint parts; index -1 writes the largest one. |
| **CAT_SurfSeparateClusters**    | Finds the connected clusters of vertex values on a surface. |
| **CAT_SurfAddValues**           | Stores vertex-wise values together with the mesh in a GIFTI file. |
| **CAT_SurfValuesAverage**       | Averages values (and meshes) of resampled GIFTI files. |
| **CAT_SurfResampleMulti**       | Resamples, and optionally smooths, values of several surfaces in one call. |
| **CAT_Surf2ROIMulti**           | Resamples annotations and computes per-ROI means for several surfaces in one call (JSON output). |
| **CAT_Surf2SPH**                | Spherical-harmonic coefficients of a surface. |
| **CAT_SurfSPH2Surf**            | Reconstructs a surface from spherical-harmonic coefficients. |
| **CAT_SurfResampleSphericalSPH**| Creates an equally sampled surface from spherical-harmonic coefficients. |
| **CAT_SurfSmoothAreal**         | Areal smoothing of surface points. |
| **CAT_SurfSmoothConvexity**     | Diffusion smoothing of surface points weighted by negative convexity. |
| **CAT_SurfSmoothSharpness**     | Diffusion smoothing of surface points weighted by sharpness. |
| **CAT_SurfSmoothDiffusion**     | Heat-kernel diffusion smoothing of values or surface points. |
| **CAT_SurfSmoothLaplacian**     | HC Laplacian smoothing of surface meshes (Vollmer et al., 1999). |
| **CAT_GlmEstimate**             | Estimates a general linear model on vertex- or voxel-wise data, from groups and covariate files or an R-style model formula (`-formula`). |
| **CAT_VolAverage**              | Averages volumes, optionally writing the standard deviation and z-scores. |
| **CAT_VolLayerSmooth**          | Smooths within cortical layers, along iso-depth contours, without crossing the GM/WM or GM/CSF boundary. |
| **CAT_VolOrnlm**                | Optimized blockwise non-local means (ORNLM) denoising. |

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
| **CAT_VolThicknessPbt** | always — the medians inside PBT are oriented (`-oriented-cutoff` sets their cutoff) |
| **CAT_VolMarchingCubes** | `-strength-sulci` — opens buried sulci in the PPM itself |

Every one of these is **a no-op where no sheet is detected**: the oriented
operator is then numerically identical to the isotropic one it replaces, which
is what makes them safe to enable. The PBT medians are oriented by default; the
other two are off until asked for.

**Buried sulci at the surface stage.** A glued sulcus — two banks of a tight
sulcus that end up as one thick grey-matter band because no CSF was detected
between them — is typical in the occipital midline, where cortex is thin and
contrast poorest. Marching cubes sees only the PPM, but no intensity image is
required to find one there, because the PPM
carries the geometry itself. Crossing a sulcus it runs 1 (WM) → 0.5 → ~0 (pial)
→ 0.5 → 1, and crossing a gyral blade it runs 0 → 0.5 → ~1 → 0.5 → 0, so a
sulcus is a valley and a blade a ridge. A buried sulcus is simply a valley whose
floor never reaches the isovalue: the shape is intact, only the amplitude is
missing, which is why the isosurface fuses the banks while the valley stays
plainly visible to a shape filter.

`CAT_VolMarchingCubes -strength-sulci` uses that field three times: valley
floors sitting just above the isovalue are pushed below it; the gyral boost of
`-strength-gyrimask` is damped there, because strengthening a thin white-matter
finger otherwise lifts the neighbouring sulcal floor back over the isovalue; and
the median filter is oriented along the sheet so it cannot re-close what was
opened.

**Getting a response at all.** Every sheetness-gated step compares the response
against a threshold, and the filter's automatic noise scale is half the *largest*
Hessian norm in the volume. On a T1 that scale is set by the strongest edges in
the image, so a thin sulcal sheet lands an order of magnitude below it and the
defaults can do nothing. The response is therefore anchored (its p99.9 is scaled
to 1) before any threshold is applied; measure it before tuning anything else:

```bash
CAT_VolSheetness -polarity -1 -v t1_corr.nii sheet.nii   # look at p99 and max
```

Pick the gain that puts the p99 of the response near twice the threshold. The
gain does **not** carry over between tools, because each measures a different
image: `CAT_VolLocalStat -sheet-strength` is measured on the image it filters,
`CAT_VolMarchingCubes -sulci-sheet-strength` on the PPM. The PPM is the
better-conditioned of the two — its dynamic range is bounded and its structures
are all of comparable scale — so it generally needs far less gain than a T1.
Run `CAT_VolMarchingCubes -strength-sulci 1 -verbose` and it prints the p99 and
maximum of the response next to the threshold, and warns outright when the gain
is too low to have any effect.

A typical sequence:

```bash
# Inspect the evidence first (dark sheets = sulcal CSF)
CAT_VolSheetness -polarity -1 -v t1_corr.nii sheetness.nii

CAT_VolThicknessPbt label.nii gmt.nii ppm.nii

# ... and again at the surface stage, on the PPM
CAT_VolMarchingCubes -strength-sulci 1 -thresh 0.5 -verbose ppm.nii central.gii
```

**Inspect that map before enabling anything downstream**, because two numbers
have to match for any of this to be visible: the response your data produces,
and the cutoff the consumer gates on.

The oriented median admits a neighbour when `s*(dhat.n)^2 < cutoff`. The 9
offsets lying in the sheet plane are always admitted; the 6 face neighbours drop
out at `s = cutoff`, the 12 edge neighbours at `2*cutoff`, the 8 corners at
`3*cutoff`. A one-voxel-thick sheet carries the sheet value on the in-plane
offsets only, so it survives the median from `s = 2*cutoff` upwards — **the
cutoff is half the sheetness at which a thin structure starts being protected**.
It defaults to 0.10 and is set with `-sheet-cutoff` (`-oriented-cutoff` in
`CAT_VolThicknessPbt`).

The response is usually much lower than expected, because the automatic noise
scale is half the largest Hessian norm in the volume — on a whole head the
scalp/air step, nothing a sulcal dip approaches. On a 0.5 mm MPRAGE the
dark-sheet map has `p99 = 0.20` and a maximum of 0.56. Raise the gain
(`-strength` in `CAT_VolSheetness`, `-sheet-strength` in `CAT_VolLocalStat`,
`-sulci-sheet-strength` in `CAT_VolMarchingCubes`), lower `-c`, or widen the
scale range so it brackets your voxel size — the scale defaults assume 0.5 mm
data:

```bash
# 1 mm data, and a gain chosen by looking at the map above
CAT_VolSheetness -polarity -1 -sigma-min 0.5 -sigma-max 1.5 -strength 8 \
    t1_corr.nii sheetness.nii
```

The gain multiplies the response and clamps it to `[0,1]`. Because it is linear
and fixes zero, no value of it can break the no-op guarantee above. It does
amplify the noise floor along with the sheets, and it lifts the strongest
responses first, so raise it while watching the map rather than blind.

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


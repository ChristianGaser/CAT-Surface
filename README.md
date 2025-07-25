# CAT: Cortex Analysis Tools
Christian Gaser christian.gaser@uni-jena.de Jena University Hospital, Germany.

CAT-Surface provides high-performance command-line tools for surface-based neuroimaging analysis, focusing on the processing, conversion, and analysis of cortical surface meshes.
These tools are integral to the [CAT12 toolbox](https://github.com/ChristianGaser/cat12) for SPM, and serve as the backend for most of CAT12’s surface-based processing, including extraction, smoothing, morphometric analysis, and format conversion. The tools are written in portable C and optimized for efficient batch processing of large cohorts in neuroimaging studies.

The repository also contains additional third-party libraries (see `3rdparty/`) for mesh handling and image I/O, bundled to simplify compilation and avoid dependency issues. These packages are not maintained within this project but are bundled for convenience so that no extra downloads are needed during compilation.

## External Dependencies

The only external library required to build CAT-Surface is:

- [FFTW3](https://www.fftw.org) (Fastest Fourier Transform in the West)  
  FFTW3 is used for efficient spectral filtering and smoothing operations on surface data.

If you do not have root access (e.g., in cloud or Codespaces environments), you can [build FFTW3 locally](#installing-fftw3-without-root-access).

---

## Build Instructions

**1. Install FFTW3:**

- **Linux:**  
  With root: `sudo apt-get install libfftw3-dev`  
  Without root: see below.
- **macOS:**  
  `brew install fftw`
- **Windows:**  
  Use MSYS2 or Windows Subsystem for Linux, then install as on Linux.

If you lack sudo rights, build FFTW3 in your home directory:
```
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar xzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
./configure --prefix=$HOME/.local
make && make install
```
Then compile CAT-Surface with:
```
export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH
./configure CPPFLAGS="-I$HOME/.local/include" LDFLAGS="-L$HOME/.local/lib"
make
```

**2. Build the tools:**
```
./autogen.sh       # Generate the configure script
./configure        # Configure the build (add --prefix if needed)
make               # Build all tools
make install       # Optional: install system-wide
```
---
## Tools and Functions

Below is a summary of the available command-line programs in CAT-Surface, each designed for a distinct step in the cortical surface processing pipeline.

| Tool                        | Description |
|-----------------------------|-------------|
| **CAT_Vol2Surf**                | Projects values from a 3D image (volume) onto the cortical surface mesh vertices. |
| **CAT_SurfApplyWarp**           | Applies deformation fields (from CAT_ApplySurf) to transform surface meshes. |
| **CAT_SurfApplyWarpValues**     | Applies surface deformations to vertex-wise data arrays (e.g., morphometric parameters). |
| **CAT_SurfSmooth**              | Performs heat kernel smoothing on surface meshes or vertex-wise data, using an exact spectral method. |
| **CAT_SurfDeform**              | Legacy algorithm (after David MacDonald) for cortical surface extraction from 3D anatomical images. Mainly for archival/historical purposes. |
| **CAT_SurfCurvature**           | Extracts folding-related surface parameters (e.g., mean curvature, Gaussian curvature, sulcal depth) and optionally smooths results using the diffusion heat kernel. |
| **CAT_DumpCurv_ui**             | GUI batch interface for CAT_DumpCurv to process large collections of surface files. |
| **CAT_SurfArea**                | Calculates total and/or local surface area metrics from cortical meshes. |
| **CAT_SurfConvert**             | Converts between surface mesh file formats: BIC (obj), Freesurfer, and OOGL. |
| **CAT_SurfMeasure2Txt**         | Converts Freesurfer curvature (`.curv`) files to plain text format. |
| **CAT_SurfPlotValuesAtMaximum** | For a given reference, extracts or plots the values at the vertex of maximum value for each surface file. |
| **CAT_SurfPlotValuesAtPoint**   | Extracts or plots the values at specified (x, y, z) coordinates for each input surface. |
| **CAT_SurfResampleSpherical**   | Resamples a spherical surface mesh onto a standard sphere (for surface-based morphometry or group analysis). |
| **CAT_SurfSheet2Surf**          | Maps a 2D image (PGM format) onto a surface mesh. |
| **CAT_Surf2Sheet**              | Flattens surface data (e.g., curvature, morphometry) onto a 2D sheet (PGM image), for visualization or further analysis. |
| **CAT_Surf2Sphere**             | Inflates a cortical surface mesh onto a sphere using the Caret/Van Essen inflation approach. |

### External binary GIFTI files

Specifying an output name with the `.dat` extension causes CAT-Surface to write
a `.gii` header alongside a `.dat` binary. The header uses the GIFTI
`ExternalFileBinary` encoding so the resulting pair is compatible with SPM12.

## Continuous Integration

A GitHub Actions workflow located in `.github/workflows/ci.yml` automatically
runs `./autogen.sh`, `./configure`, `make` and a small test on every push or
pull request.


# Cython-based BICPL Surface IO

This module provides a minimal Cython extension to read BICPL `.obj` surfaces
from the 3rdparty bicpl-surface library bundled in this repository and expose
them as NumPy arrays (vertices, faces) with a small nibabel-like wrapper.

Status: experimental, triangles-only (non-tri polygons are skipped).

## Build

Requirements:
- Python 3.8+
- Cython, NumPy
- A built CAT-Surface tree (so that `build-native-arm64/.libs/libCAT.a` exists),
  or set `CAT_BUILD_DIR` to a different build folder with `.libs/libCAT.a`.

Build and develop:

```
python -m pip install -e ./cat_surface_cython[build]
python ./cat_surface_cython/setup.py build_ext --inplace
```

If your build dir is different, set:

```
export CAT_BUILD_DIR=/path/to/CAT-Surface/build-native
```

## Usage

```
from cat_surface import load_bic_surface
surf = load_bic_surface("/path/to/surface.obj")
print(surf)
print(surf.vertices.shape, surf.triangles.shape)
```

## Notes
- Only the first POLYGONS object in the file is returned.
- Faces must be triangles; others are skipped.
- This extension links against `libCAT.a`, which already links the required
  3rdparty code.

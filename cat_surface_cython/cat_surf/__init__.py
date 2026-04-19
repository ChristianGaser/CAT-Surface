"""
cat_surf — Cython bindings to the CAT-Surface library (libCAT).

Provides fast, in-memory access to cortical surface mesh algorithms
without subprocess overhead.  All functions accept and return plain
numpy arrays so they work seamlessly with nibabel, surfa, trimesh,
or any other mesh library.

Quick start
-----------
>>> import numpy as np
>>> import cat_surf
>>>
>>> # Load a surface (uses libCAT multi-format I/O)
>>> vertices, faces = cat_surf.read_surface("lh.central.gii")
>>>
>>> # Or bring your own arrays via nibabel:
>>> # import nibabel as nib
>>> # gii = nib.load("lh.central.gii")
>>> # vertices = gii.agg_data("pointset")
>>> # faces = gii.agg_data("triangle").astype(np.int32)
>>>
>>> # Compute per-vertex area
>>> area, total = cat_surf.get_area(vertices, faces)
>>>
>>> # Euler characteristic
>>> chi = cat_surf.euler_characteristic(vertices, faces)
>>>
>>> # Smooth per-vertex data
>>> import cat_surf
>>> smoothed = cat_surf.smooth_heatkernel(vertices, faces, area, fwhm=3.0)
"""

__version__ = "0.1.0"

# --- I/O functions ---
from cat_surf._io import (
    read_surface,
    write_surface,
    read_values,
    write_values,
)

# --- Core surface operations ---
from cat_surf._surf import (
    get_area,
    euler_characteristic,
    point_distance,
    point_distance_mean,
    hausdorff_distance,
    smooth_heatkernel,
    smoothed_curvatures,
    sulcus_depth,
    reduce_mesh,
    sphere_radius,
    smooth_mesh,
)

__all__ = [
    # I/O
    "read_surface",
    "write_surface",
    "read_values",
    "write_values",
    # Surface operations
    "get_area",
    "euler_characteristic",
    "point_distance",
    "point_distance_mean",
    "hausdorff_distance",
    "smooth_heatkernel",
    "smoothed_curvatures",
    "sulcus_depth",
    "reduce_mesh",
    "sphere_radius",
    "smooth_mesh",
]

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

try:
    from importlib.metadata import version
    __version__ = version("cat-surf")
except Exception:  # pragma: no cover
    __version__ = "0.0.dev0"

import os as _os
import sys as _sys


def _redirect_libomp_to_torch_if_present():
    """
    macOS only: when PyTorch is installed in the same env, rewrite each
    cat_surf .so file's libomp.dylib dependency from our bundled
    ``@loader_path/.dylibs/libomp.dylib`` to PyTorch's
    ``<torch>/lib/libomp.dylib``.

    Why this exists.  Two libomp instances in the same process corrupt
    each other's thread-pool TLS in ways that survive
    ``KMP_DUPLICATE_LIB_OK=TRUE`` (silent crash inside the first
    parallel region we run).  The robust fix is to ensure exactly one
    libomp.dylib per process by pointing cat_surf at the same dylib
    PyTorch already loads.

    The rewrite is done once per install via ``install_name_tool`` and
    recorded in a marker file inside ``cat_surf/.dylibs/`` so subsequent
    imports skip the work.  If the rewrite fails (read-only install,
    torch unavailable, ...), we fall back to the bundled libomp and let
    ``KMP_DUPLICATE_LIB_OK`` paper over the conflict.
    """
    if _sys.platform != "darwin":
        return

    pkg_dir = _os.path.dirname(__file__)
    dylibs_dir = _os.path.join(pkg_dir, ".dylibs")
    marker = _os.path.join(dylibs_dir, ".libomp_redirected_to_torch")
    if _os.path.isfile(marker):
        return

    # Locate torch's libomp.dylib without importing torch (importing
    # torch initialises its libomp; that's fine but slow if the user
    # doesn't actually need it yet).  We rely on torch's installed
    # location being discoverable from sys.path.
    torch_libomp = None
    for sp in _sys.path:
        candidate = _os.path.join(sp, "torch", "lib", "libomp.dylib")
        if _os.path.isfile(candidate):
            torch_libomp = candidate
            break
    if torch_libomp is None:
        return

    import glob as _glob
    import subprocess as _sp
    bundled_dep = "@loader_path/.dylibs/libomp.dylib"
    so_files = _glob.glob(_os.path.join(pkg_dir, "*.so"))
    rewrote_any = False
    for so in so_files:
        try:
            _sp.check_call([
                "install_name_tool", "-change",
                bundled_dep, torch_libomp, so,
            ], stderr=_sp.DEVNULL)
        except (_sp.CalledProcessError, FileNotFoundError, OSError):
            continue  # install_name_tool missing or .so read-only
        # install_name_tool invalidates the code signature.  macOS's
        # hardened runtime then SIGKILLs the process at load time.
        # Re-sign ad-hoc to restore loadability.
        try:
            _sp.check_call(
                ["codesign", "--remove-signature", so],
                stderr=_sp.DEVNULL,
            )
        except (_sp.CalledProcessError, FileNotFoundError, OSError):
            pass
        try:
            _sp.check_call(
                ["codesign", "--sign", "-", so],
                stderr=_sp.DEVNULL,
            )
            rewrote_any = True
        except (_sp.CalledProcessError, FileNotFoundError, OSError):
            pass

    if rewrote_any:
        try:
            _os.makedirs(dylibs_dir, exist_ok=True)
            with open(marker, "w") as _f:
                _f.write(torch_libomp + "\n")
        except OSError:
            pass


_redirect_libomp_to_torch_if_present()

# Even after the redirect above, keep the duplicate-init workaround as a
# safety net for installs where the rewrite couldn't run (read-only
# install) or for users who haven't installed torch.
_os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

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
    surf_to_sphere,
    get_area_normalized,
    surf_average,
    correct_thickness_folding,
    remove_intersections,
    count_intersections,
    resample_to_sphere,
    surf_deform,
    surf_to_pial_white,
    central_to_pial,
)

# --- Volume operations ---
from cat_surf._vol import (
    vol_sanlm,
    vol_blood_vessel_correction,
    vol_thickness_pbt,
    vol_amap,
    vol_marching_cubes,
)

# --- Registration ---
from cat_surf._bbreg import (
    bbreg,
    bbreg_detect_contrast,
    volume_register_nmi,
    volume_register_robust,
)

# --- Volume-to-surface mapping ---
from cat_surf._vol2surf import vol2surf

# --- DARTEL spherical surface registration ---
from cat_surf._surf_warp import surf_warp

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
    "surf_to_sphere",
    "get_area_normalized",
    "surf_average",
    "correct_thickness_folding",
    "remove_intersections",
    "count_intersections",
    "resample_to_sphere",
    "surf_deform",
    "surf_to_pial_white",
    "central_to_pial",
    # Volume operations
    "vol_sanlm",
    "vol_blood_vessel_correction",
    "vol_thickness_pbt",
    "vol_amap",
    "vol_marching_cubes",
    # Registration
    "bbreg",
    "bbreg_detect_contrast",
    "volume_register_nmi",
    "volume_register_robust",
    # Volume-to-surface
    "vol2surf",
    # DARTEL surface warp
    "surf_warp",
]

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

# Belt-and-suspenders for macOS users who somehow end up with two
# libomp instances in the same process despite the wheel's @rpath
# resolution preferring PyTorch's libomp (see setup.py rpath list).
# Without this env var, Intel's OpenMP runtime aborts with
# "OMP: Error #15: ... libomp.dylib already initialized."
import os as _os
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
    resample_annot,
    surf_curvature,
    surf_deform,
    surf_to_pial_white,
    central_to_pial,
    surf2roi_unit,
    surf_ratio,
    surf_fractal_dimension,
)

# --- Volume operations ---
from cat_surf._vol import (
    vol_sanlm,
    vol_blood_vessel_correction,
    vol_thickness_pbt,
    vol_amap,
    vol_marching_cubes,
    vol_smooth,
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

# --- Spherical Demons spherical surface registration ---
from cat_surf._spherical_demon import spherical_demon


def surf2roi_multi(units):
    """Resample annotations and compute per-ROI means for multiple units.

    Mirrors ``CAT_Surf2ROIMulti``.  Each unit dict specifies a source sphere,
    target sphere, annotation file, and values file.  Units sharing the same
    hemisphere are merged into a single per-ROI table.

    Parameters
    ----------
    units : list of dict
        Each dict has:

        ``src_sphere`` : str
            Source sphere matching the annotation vertex count.
        ``trg_sphere`` : str
            Target sphere to resample onto.
        ``annot`` : str
            FreeSurfer ``.annot`` annotation file.
        ``vals`` : str
            Per-vertex values file on the *target* sphere.
        ``hemi`` : str, optional
            ``"lh"`` or ``"rh"``; guessed from file-name patterns when
            omitted.

    Returns
    -------
    roi_stats : dict
        Keys ``"lh"``, ``"rh"``, ``"unknown"``.  Each value is a list of
        dicts ``{"id": int, "name": str, "mean": float, "n": int}``,
        sorted by ROI id.
    """
    def _guess_hemi(*paths):
        for p in paths:
            if p is None:
                continue
            for m in ("lh.", "_lh", "/lh"):
                if m in p:
                    return "lh"
            for m in ("rh.", "_rh", "/rh"):
                if m in p:
                    return "rh"
        return None

    buckets = {"lh": {}, "rh": {}, "unknown": {}}

    for unit in units:
        hemi = unit.get("hemi") or _guess_hemi(
            unit.get("annot"), unit.get("vals"),
            unit.get("src_sphere"), unit.get("trg_sphere"),
        )
        hemi_key = (
            "lh" if hemi and hemi.lower() in ("lh", "left") else
            "rh" if hemi and hemi.lower() in ("rh", "right") else
            "unknown"
        )

        stats = surf2roi_unit(
            unit["src_sphere"], unit["trg_sphere"],
            unit["annot"], unit["vals"]
        )
        bucket = buckets[hemi_key]
        for s in stats:
            roi_id = s["id"]
            if roi_id not in bucket:
                bucket[roi_id] = {"id": roi_id, "name": s["name"],
                                   "sum": 0.0, "n": 0}
            bucket[roi_id]["sum"] += s["sum"]
            bucket[roi_id]["n"]   += s["n"]

    out = {}
    for hk in ("lh", "rh", "unknown"):
        rows = []
        for d in buckets[hk].values():
            n = d["n"]
            rows.append({
                "id":   d["id"],
                "name": d["name"],
                "mean": d["sum"] / n if n > 0 else float("nan"),
                "n":    n,
            })
        rows.sort(key=lambda x: x["id"])
        out[hk] = rows
    return out


def resample_multi(units, fwhm=0.0, areal=False):
    """Resample and optionally smooth multiple surface units, then concatenate.

    Mirrors ``CAT_SurfResampleMulti``.  Each unit resamples one surface
    (and optional per-vertex values) from its own source sphere to its own
    target sphere, then all units are concatenated into a single mesh.

    Parameters
    ----------
    units : list of dict
        Each dict has:

        ``surf`` : (vertices, faces)
            Source mesh arrays.
        ``src_sphere`` : (vertices, faces)
            Source sphere arrays (same topology as ``surf``).
        ``trg_sphere`` : (vertices, faces)
            Target sphere arrays.
        ``vals`` : array_like or None, optional
            Per-vertex values on the source sphere.
        ``mask`` : array_like or None, optional
            Binary mask on the *target* sphere; values at masked-out
            vertices are set to ``NaN`` before optional smoothing.
        ``fwhm`` : float, optional
            Per-unit smoothing FWHM in mm.  Overrides the global *fwhm*
            argument when provided and >= 0.

    fwhm : float
        Global default smoothing FWHM in mm (used for units without
        ``fwhm`` key).  0 = no smoothing.
    areal : bool
        Use areal (sum-preserving) interpolation for values.

    Returns
    -------
    vertices : ndarray, shape (V_total, 3), float64
        Concatenated vertex positions.
    faces : ndarray, shape (F_total, 3), int32
        Concatenated faces (indices shifted per unit).
    values : ndarray or None, shape (V_total,), float64
        Concatenated values if any unit had ``vals``; otherwise ``None``.
    """
    import numpy as np
    all_verts = []
    all_faces = []
    all_vals  = []
    has_vals  = False
    v_offset  = 0

    for unit in units:
        sv, sf = unit["surf"]
        ssv, ssf = unit["src_sphere"]
        tsv, tsf = unit["trg_sphere"]
        u_vals = unit.get("vals")
        u_mask = unit.get("mask")
        u_fwhm = unit.get("fwhm", None)
        if u_fwhm is None or float(u_fwhm) < 0.0:
            u_fwhm = fwhm

        nv, nf, nvals = resample_to_sphere(
            sv, sf, ssv, ssf, tsv, tsf,
            values=u_vals,
            label_interpolation=False,
            areal_interpolation=areal,
        )

        if nvals is not None:
            has_vals = True
            if u_mask is not None:
                import numpy as _np
                mask_arr = _np.asarray(u_mask, dtype=bool)
                nvals = nvals.copy()
                nvals[~mask_arr] = float("nan")
            if float(u_fwhm) > 0.0:
                nvals = smooth_heatkernel(nv, nf, nvals, fwhm=float(u_fwhm))
        else:
            nvals = None

        nf_shifted = np.asarray(nf, dtype=np.int32) + v_offset
        all_verts.append(nv)
        all_faces.append(nf_shifted)
        all_vals.append(nvals)
        v_offset += len(nv)

    cat_verts = np.concatenate(all_verts, axis=0)
    cat_faces = np.concatenate(all_faces, axis=0).astype(np.int32)
    if has_vals:
        cat_vals = np.concatenate(
            [v if v is not None else np.full(len(all_verts[i]), float("nan"))
             for i, v in enumerate(all_vals)],
            axis=0,
        )
    else:
        cat_vals = None

    return cat_verts, cat_faces, cat_vals


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
    "resample_annot",
    "surf_curvature",
    "surf_deform",
    "surf_to_pial_white",
    "central_to_pial",
    "surf2roi_unit",
    "surf2roi_multi",
    "resample_multi",
    "surf_ratio",
    "surf_fractal_dimension",
    # Volume operations
    "vol_smooth",
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
    # Spherical Demons surface registration
    "spherical_demon",
]

"""
cat_surf.cli — file-path mirror of the CAT_* command-line binaries.

Each function in this module mirrors exactly one ``CAT_<Tool>`` binary
in name, positional-argument order, option semantics, and defaults.
The implementation is a thin shim:

1. read inputs from disk via :mod:`cat_surf.read_surface` / ``read_values``,
2. call the in-memory numpy wrapper in :mod:`cat_surf`,
3. write the outputs back to disk.

This is the right entry point if you are porting shell scripts that
call the CAT binaries, or if you want a drop-in Python replacement
for the binaries.  For composable in-memory pipelines, use the
lower-level :mod:`cat_surf` API directly.

Naming convention
-----------------
The binary ``CAT_<X>`` maps to ``cat_surf.cli.<x>`` where ``<x>`` is
``<X>`` in snake_case with the ``CAT_`` prefix dropped::

    CAT_Surf2PialWhite                -> surf2pial_white
    CAT_Surf2Sphere                   -> surf2sphere
    CAT_SurfArea                      -> surf_area
    CAT_SurfAverage                   -> surf_average
    CAT_SurfCorrectThicknessFolding   -> surf_correct_thickness_folding
    CAT_SurfCurvature                 -> surf_curvature
    CAT_SurfDeform                    -> surf_deform
    CAT_SurfDistance                  -> surf_distance
    CAT_SurfReduce                    -> surf_reduce
    CAT_SurfRemoveIntersections       -> surf_remove_intersections
    CAT_SurfResample                  -> surf_resample
    CAT_SurfResample -label <annot>   -> surf_resample_annot
    CAT_SurfResampleMulti             -> surf_resample_multi
    CAT_Surf2ROIMulti                 -> surf2roi_multi
    CAT_SurfWarp                      -> surf_warp  (use avg=True for -avg)
    CAT_Vol2Surf                      -> vol2surf
    CAT_VolAmap                       -> vol_amap
    CAT_VolBloodVesselCorrection      -> vol_blood_vessel_correction
    CAT_VolMarchingCubes              -> vol_marching_cubes
    CAT_VolSanlm                      -> vol_sanlm
    CAT_VolSmooth                     -> vol_smooth
    CAT_VolThicknessPbt               -> vol_thickness_pbt
"""
from __future__ import annotations

import os
import numpy as np

from cat_surf._io import (
    read_surface,
    write_surface,
    read_values,
    write_values,
)
from cat_surf import (
    surf_to_pial_white as _surf_to_pial_white,
    surf_to_sphere as _surf_to_sphere,
    get_area as _get_area,
    get_area_normalized as _get_area_normalized,
    surf_average as _surf_average,
    correct_thickness_folding as _correct_thickness_folding,
    surf_deform as _surf_deform,
    point_distance as _point_distance,
    point_distance_mean as _point_distance_mean,
    reduce_mesh as _reduce_mesh,
    remove_intersections as _remove_intersections,
    count_intersections as _count_intersections,
    resample_to_sphere as _resample_to_sphere,
    resample_annot as _resample_annot,
    surf_curvature as _surf_curvature,
    surf_warp as _surf_warp,
    vol2surf as _vol2surf,
    vol_amap as _vol_amap,
    vol_blood_vessel_correction as _vol_blood_vessel_correction,
    vol_marching_cubes as _vol_marching_cubes,
    vol_sanlm as _vol_sanlm,
    vol_thickness_pbt as _vol_thickness_pbt,
    vol_smooth as _vol_smooth,
    surf2roi_multi as _surf2roi_multi,
    resample_multi as _resample_multi,
)

# ---------------------------------------------------------------------------
# Volume I/O helpers (read float32 array + affine, write back)
# ---------------------------------------------------------------------------

def _load_volume_meta(path):
    """Load a NIfTI volume as (float32_array, affine, voxelsize)."""
    try:
        import nibabel as nib
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "cat_surf.cli volume operations require nibabel for writing back "
            "results.  Install nibabel or use the lower-level cat_surf API."
        ) from exc
    img = nib.load(path)
    return img


def _save_volume_like(out_path, data, like_img, dtype=None):
    """Write *data* to *out_path* using *like_img*'s header and affine."""
    import nibabel as nib
    if dtype is not None:
        data = data.astype(dtype, copy=False)
    out = nib.Nifti1Image(data, like_img.affine, like_img.header)
    out.set_data_dtype(data.dtype)
    nib.save(out, out_path)


# ===========================================================================
# Surface tools
# ===========================================================================

def surf2pial_white(surface_file, thickness_file, label_file,
                    pial_file, white_file, **kwargs):
    """Mirror of ``CAT_Surf2PialWhite``.

    Estimate pial and white-matter surfaces from a central surface and
    cortical thickness, writing both surfaces to disk.

    Positional arguments
    --------------------
    surface_file : str
        Input central surface.
    thickness_file : str
        Per-vertex thickness values.
    label_file : str
        NIfTI tissue label volume.
    pial_file : str
        Output pial surface.
    white_file : str
        Output white-matter surface.

    Keyword arguments are forwarded to :func:`cat_surf.surf_to_pial_white`.
    """
    v, f = read_surface(surface_file)
    t = read_values(thickness_file)
    pv, pf, wv, wf = _surf_to_pial_white(v, f, t, label_file, **kwargs)
    write_surface(pial_file, pv, pf)
    write_surface(white_file, wv, wf)


def surf2sphere(surface_file, output_file, stop_at=5, verbose=False):
    """Mirror of ``CAT_Surf2Sphere``.

    Maps a surface to a sphere using iterative inflation.
    """
    v, f = read_surface(surface_file)
    sv, sf = _surf_to_sphere(v, f, stop_at=stop_at, verbose=verbose)
    write_surface(output_file, sv, sf)


def surf_area(surface_file, output_values_file=None, sphere_file=None,
              use_log=False):
    """Mirror of ``CAT_SurfArea``.

    Returns the total surface area.  If ``output_values_file`` is given,
    also writes per-vertex area values.
    """
    v, f = read_surface(surface_file)
    if sphere_file is not None:
        sv, sf = read_surface(sphere_file)
        area, total = _get_area_normalized(v, f, sv, sf)
    else:
        area, total = _get_area(v, f)
    if use_log:
        area = np.log10(area)
        total = float(np.log10(total))
    if output_values_file is not None:
        write_values(output_values_file, area)
    return total


def surf_average(output_avg_file, *surface_files, rms_file=None,
                 verbose=False):
    """Mirror of ``CAT_SurfAverage``.

    Average vertex-wise across surfaces and write the result.  All
    surfaces must share topology.  ``rms_file`` (optional) gets per-vertex
    root-mean-square deviation.
    """
    if not surface_files:
        raise ValueError("at least one surface_file is required")
    surfaces = [read_surface(p) for p in surface_files]
    if verbose:
        for i, p in enumerate(surface_files):
            print(f"{i}:  {p}")
    if rms_file is not None:
        avg_v, faces, rms = _surf_average(surfaces, return_rms=True)
        write_values(rms_file, rms)
    else:
        avg_v, faces = _surf_average(surfaces)
    write_surface(output_avg_file, avg_v, faces)


def surf_correct_thickness_folding(surface_file, thickness_file,
                                   output_thickness_file,
                                   slope=0.0, max_dist=6.0):
    """Mirror of ``CAT_SurfCorrectThicknessFolding``."""
    v, f = read_surface(surface_file)
    t = read_values(thickness_file)
    out = _correct_thickness_folding(v, f, t, slope=slope, max_dist=max_dist)
    write_values(output_thickness_file, out)


def surf_deform(volume_file, surface_file, output_surface_file, **kwargs):
    """Mirror of ``CAT_SurfDeform``.

    Note: argument order matches the binary (volume first, then surface).
    """
    v, f = read_surface(surface_file)
    nv, nf = _surf_deform(v, f, volume_file, **kwargs)
    write_surface(output_surface_file, nv, nf)


def surf_distance(surface_file, surface_file2_or_None,
                  output_values_file, thickness_file=None,
                  mode="mean", max_dist=6.0,
                  check_intersect=False, verbose=False):
    """Mirror of ``CAT_SurfDistance``.

    Two modes:
    - Two surfaces: pass both surface files.
    - Central + thickness: pass ``surface_file`` and set
      ``thickness_file``; ``surface_file2_or_None`` is ignored.

    ``mode`` is ``"mean"`` (Tfs, default) or ``"link"``.
    """
    if mode not in ("mean", "link"):
        raise ValueError("mode must be 'mean' or 'link'")
    v1, f1 = read_surface(surface_file)
    if thickness_file is not None:
        from cat_surf import central_to_pial
        t = read_values(thickness_file)
        pv, pf = central_to_pial(v1, f1, t, sigma=0.0, iterations=0,
                                  check_intersect=check_intersect,
                                  verbose=verbose)
        wv, wf = central_to_pial(v1, f1, -t, sigma=0.0, iterations=0,
                                  check_intersect=check_intersect,
                                  verbose=verbose)
        v1, f1, v2, f2 = pv, pf, wv, wf
    else:
        if surface_file2_or_None is None:
            raise ValueError("second surface required when thickness_file is None")
        v2, f2 = read_surface(surface_file2_or_None)
    if mode == "mean":
        d, _ = _point_distance_mean(v1, f1, v2, f2, symmetric=False,
                                     max_dist=max_dist)
    else:
        d, _ = _point_distance(v1, f1, v2, f2, symmetric=False,
                                max_dist=max_dist)
    write_values(output_values_file, d)


def surf_reduce(input_file, output_file, ratio=0.5, aggressiveness=7.0,
                preserve_sharp=True, verbose=False):
    """Mirror of ``CAT_SurfReduce``.

    ``ratio`` < 1 keeps that fraction of faces; ``ratio`` >= 1 is the
    absolute target face count.
    """
    v, f = read_surface(input_file)
    if ratio <= 0.0:
        raise ValueError("ratio must be > 0")
    target = int(round(ratio * f.shape[0])) if ratio < 1.0 else int(ratio)
    nv, nf = _reduce_mesh(v, f, target_faces=target,
                          aggressiveness=aggressiveness,
                          preserve_sharp=preserve_sharp, verbose=verbose)
    write_surface(output_file, nv, nf)


def surf_remove_intersections(input_file, output_file=None, *,
                              count_only=False, max_iters=10,
                              inner_loops=3, fill_holes=True, verbose=False):
    """Mirror of ``CAT_SurfRemoveIntersections``.

    When ``count_only=True``, returns the number of self-intersections
    without writing any output (matches ``-count``).
    """
    v, f = read_surface(input_file)
    if count_only:
        return _count_intersections(v, f)
    if output_file is None:
        raise ValueError("output_file required unless count_only=True")
    nv, nf = _remove_intersections(v, f, max_iters=max_iters,
                                    inner_loops=inner_loops,
                                    fill_holes=fill_holes, verbose=verbose)
    write_surface(output_file, nv, nf)


def surf_resample(surface_file_or_None, sphere_file_or_None,
                  target_sphere_file, output_surface_file_or_None,
                  input_values_file=None, output_values_file=None,
                  label=False, areal=False):
    """Mirror of ``CAT_SurfResample``.

    Resamples a surface (and optional per-vertex values) from one
    spherical parameterization to another.  Use ``"NULL"`` or ``None``
    for unused file slots, matching the CLI.
    """
    def _is_null(p):
        return p is None or p == "NULL"

    src_v = src_f = None
    if not _is_null(surface_file_or_None):
        src_v, src_f = read_surface(surface_file_or_None)
    sph_v = sph_f = None
    if not _is_null(sphere_file_or_None):
        sph_v, sph_f = read_surface(sphere_file_or_None)
    tsph_v, tsph_f = read_surface(target_sphere_file)

    vals = None
    if input_values_file is not None:
        vals = read_values(input_values_file)

    if src_v is None or sph_v is None:
        raise ValueError(
            "resample_to_sphere requires both surface and sphere; pass them "
            "explicitly in this Python API")
    nv, nf, nvals = _resample_to_sphere(
        src_v, src_f, sph_v, sph_f, tsph_v, tsph_f,
        values=vals, label_interpolation=label, areal_interpolation=areal)

    if not _is_null(output_surface_file_or_None):
        write_surface(output_surface_file_or_None, nv, nf)
    if output_values_file is not None and nvals is not None:
        write_values(output_values_file, nvals)


def surf_warp(source_file, source_sphere_file,
              target_file, target_sphere_file,
              output_sphere_file, **kwargs):
    """Mirror of ``CAT_SurfWarp`` (DARTEL spherical registration).

    Writes the warped source sphere.  The Jacobian-determinant / PGM
    outputs of the binary are not surfaced.  Pass ``avg=True`` to
    enable the CLI's ``-avg`` flag (pole-rotated double run averaged).
    """
    sv, sf = read_surface(source_file)
    ssv, ssf = read_surface(source_sphere_file)
    tv, tf = read_surface(target_file)
    tsv, tsf = read_surface(target_sphere_file)
    wv, wf = _surf_warp((sv, sf), (ssv, ssf), (tv, tf), (tsv, tsf), **kwargs)
    write_surface(output_sphere_file, wv, wf)


def surf_curvature(surface_file, output_values_file, curvtype=0,
                   fwhm=0.0, use_abs_values=False, invert_values=False):
    """Mirror of ``CAT_SurfCurvature``.

    Compute a per-vertex curvature map and write it to disk.  See
    :func:`cat_surf.surf_curvature` for the full ``curvtype`` table.
    """
    v, f = read_surface(surface_file)
    curv = _surf_curvature(v, f, curvtype=curvtype, fwhm=fwhm,
                           use_abs_values=use_abs_values,
                           invert_values=invert_values)
    write_values(output_values_file, curv)


def surf_resample_annot(source_surface_file, source_sphere_file,
                        target_sphere_file,
                        annot_in_file, annot_out_file):
    """Mirror of ``CAT_SurfResample -label … <annot_in> <annot_out>``.

    Resample a FreeSurfer ``.annot`` label file from one cortical
    spherical parameterisation to another.  Used by T1Prep's atlas
    label-resampling step.
    """
    return _resample_annot(source_surface_file, source_sphere_file,
                           target_sphere_file,
                           annot_in_file, annot_out_file)


# ===========================================================================
# Volume tools
# ===========================================================================

def vol2surf(surface_file, volume_file, output_values_file,
             thickness_file=None, offset_file=None, **kwargs):
    """Mirror of ``CAT_Vol2Surf``.

    Map a volume to a surface and write per-vertex values.
    """
    v, f = read_surface(surface_file)
    t = read_values(thickness_file) if thickness_file else None
    o = read_values(offset_file) if offset_file else None
    values, _ = _vol2surf(volume_file, v, f, thickness=t, offset=o, **kwargs)
    if values.ndim == 1:
        write_values(output_values_file, values)
    else:
        # multi-grid output: write per-step file with grid value in name
        base, ext = os.path.splitext(output_values_file)
        if ext.lower() == ".txt":
            ext = ".txt"
        else:
            ext = ""
        for j in range(values.shape[1]):
            write_values(f"{base}_g{j:02d}{ext}", values[:, j])


def vol_amap(input_file, label_file, output_file=None, **kwargs):
    """Mirror of ``CAT_VolAmap``.

    Note: only the segmentation label volume is written (matching the
    binary's ``_seg`` output).  For probability maps / bias field use
    :func:`cat_surf.vol_amap` directly.
    """
    import nibabel as nib
    img = nib.load(input_file)
    label_img = nib.load(label_file)
    vol = img.get_fdata().astype(np.float32)
    lab = label_img.get_fdata().astype(np.uint8)
    vx = img.header.get_zooms()[:3]
    prob, lab_out, mean = _vol_amap(vol, lab, voxelsize=vx, **kwargs)
    if output_file is None:
        base, _ext = os.path.splitext(os.path.splitext(input_file)[0])
        output_file = base + "_seg.nii"
    _save_volume_like(output_file, lab_out, img, dtype=np.uint8)


def vol_blood_vessel_correction(input_file, output_file=None):
    """Mirror of ``CAT_VolBloodVesselCorrection``."""
    import nibabel as nib
    img = nib.load(input_file)
    vol = img.get_fdata().astype(np.float32)
    vx = img.header.get_zooms()[:3]
    out = _vol_blood_vessel_correction(vol, voxelsize=vx)
    if output_file is None:
        d = os.path.dirname(input_file) or "."
        b = os.path.basename(input_file)
        output_file = os.path.join(d, "bvc_" + b)
    _save_volume_like(output_file, out, img)


def vol_marching_cubes(input_file, output_surface_file, label_file=None,
                       **kwargs):
    """Mirror of ``CAT_VolMarchingCubes``."""
    v, f = _vol_marching_cubes(input_file, label=label_file, **kwargs)
    write_surface(output_surface_file, v, f)


def vol_sanlm(input_file, output_file=None, is_rician=False, strength=1.0):
    """Mirror of ``CAT_VolSanlm``."""
    import nibabel as nib
    img = nib.load(input_file)
    vol = img.get_fdata().astype(np.float32)
    out = _vol_sanlm(vol, is_rician=is_rician, strength=strength)
    if output_file is None:
        d = os.path.dirname(input_file) or "."
        b = os.path.basename(input_file)
        output_file = os.path.join(d, "n" + b)
    _save_volume_like(output_file, out, img, dtype=np.float32)


def vol_thickness_pbt(input_file, gmt_file=None, ppm_file=None,
                      dist_csf_file=None, dist_wm_file=None, **kwargs):
    """Mirror of ``CAT_VolThicknessPbt``.

    Writes any of the four outputs whose file path is provided; pass
    ``None`` to skip an output.
    """
    import nibabel as nib
    img = nib.load(input_file)
    vol = img.get_fdata().astype(np.float32)
    vx = img.header.get_zooms()[:3]
    gmt, ppm, dcsf, dwm = _vol_thickness_pbt(vol, voxelsize=vx, **kwargs)
    if gmt_file:
        _save_volume_like(gmt_file, gmt, img, dtype=np.float32)
    if ppm_file:
        _save_volume_like(ppm_file, ppm, img, dtype=np.float32)
    if dist_csf_file:
        _save_volume_like(dist_csf_file, dcsf, img, dtype=np.float32)
    if dist_wm_file:
        _save_volume_like(dist_wm_file, dwm, img, dtype=np.float32)


def vol_smooth(input_file, output_file=None, fwhm=8.0, use_mask=False):
    """Mirror of ``CAT_VolSmooth``.

    Smooth a NIfTI volume with an isotropic Gaussian kernel and write
    the result to disk.

    Parameters
    ----------
    input_file : str
        Input NIfTI file.
    output_file : str, optional
        Output path.  Defaults to the same directory as ``input_file``
        with a ``s<fwhm>`` prefix prepended to the filename, matching
        the ``CAT_VolSmooth`` default naming convention.
    fwhm : float
        Full-width at half-maximum in mm (default 8.0).
    use_mask : bool
        Use masked smoothing (default False).
    """
    import nibabel as nib
    img = nib.load(input_file)
    vol = img.get_fdata().astype(np.float32)
    vx = img.header.get_zooms()[:3]
    out = _vol_smooth(vol, voxelsize=vx, fwhm=fwhm, use_mask=use_mask)
    if output_file is None:
        d = os.path.dirname(input_file) or "."
        b = os.path.basename(input_file)
        output_file = os.path.join(d, f"s{fwhm}{b}")
    _save_volume_like(output_file, out, img, dtype=np.float32)


def surf_resample_multi(units, output_surface_file, output_values_file=None,
                        fwhm=0.0, areal=False):
    """Mirror of ``CAT_SurfResampleMulti``.

    Resample and optionally smooth multiple surface units, concatenate
    them, and write the combined mesh (and values) to disk.

    Parameters
    ----------
    units : list of dict
        Each dict has:

        ``surf`` : str
            Source surface file path.
        ``src_sphere`` : str
            Source sphere file path (matches ``surf`` topology).
        ``trg_sphere`` : str
            Target sphere file path.
        ``vals`` : str or None, optional
            Per-vertex values file on the source sphere.
        ``mask`` : str or None, optional
            Binary mask file on the *target* sphere.
        ``fwhm`` : float, optional
            Per-unit smoothing FWHM in mm.  Overrides global *fwhm*.

    output_surface_file : str
        Output file for the concatenated surface.
    output_values_file : str, optional
        Output file for the concatenated values.  Written only when at
        least one unit provided ``vals``.
    fwhm : float
        Global default smoothing FWHM for units without ``fwhm``.
    areal : bool
        Use areal (sum-preserving) interpolation.
    """
    loaded_units = []
    for u in units:
        sv, sf = read_surface(u["surf"])
        ssv, ssf = read_surface(u["src_sphere"])
        tsv, tsf = read_surface(u["trg_sphere"])
        vals = read_values(u["vals"]) if u.get("vals") else None
        mask = read_values(u["mask"]) if u.get("mask") else None
        loaded_units.append({
            "surf":       (sv, sf),
            "src_sphere": (ssv, ssf),
            "trg_sphere": (tsv, tsf),
            "vals":       vals,
            "mask":       mask,
            "fwhm":       u.get("fwhm", None),
        })

    cat_v, cat_f, cat_vals = _resample_multi(loaded_units, fwhm=fwhm, areal=areal)
    write_surface(output_surface_file, cat_v, cat_f)
    if output_values_file is not None and cat_vals is not None:
        write_values(output_values_file, cat_vals)


def surf2roi_multi(units, output_json_file):
    """Mirror of ``CAT_Surf2ROIMulti``.

    Resample annotation labels onto target spheres and compute per-ROI
    means, writing the results as a JSON file.

    Parameters
    ----------
    units : list of dict
        Each dict has:

        ``src_sphere`` : str
            Source sphere file path (matches annotation vertex count).
        ``trg_sphere`` : str
            Target sphere file path.
        ``annot`` : str
            FreeSurfer ``.annot`` annotation file.
        ``vals`` : str
            Per-vertex values file on the *target* sphere.
        ``hemi`` : str, optional
            ``"lh"`` or ``"rh"``; guessed from file-name patterns when
            omitted.

    output_json_file : str
        Output JSON file path.  The JSON structure mirrors the binary::

            {
              "lh": [{"id": <int>, "name": "<str>", "mean": <float>, "n": <int>}, ...],
              "rh": [...],
              "unknown": [...]
            }
    """
    import json
    roi_stats = _surf2roi_multi(units)
    with open(output_json_file, "w") as fp:
        json.dump(roi_stats, fp, indent=2)


__all__ = [
    # Surface tools
    "surf2pial_white",
    "surf2sphere",
    "surf_area",
    "surf_average",
    "surf_correct_thickness_folding",
    "surf_deform",
    "surf_distance",
    "surf_reduce",
    "surf_remove_intersections",
    "surf_resample",
    "surf_resample_annot",
    "surf_resample_multi",
    "surf_curvature",
    "surf_warp",
    "surf2roi_multi",
    # Volume tools
    "vol2surf",
    "vol_amap",
    "vol_blood_vessel_correction",
    "vol_marching_cubes",
    "vol_sanlm",
    "vol_smooth",
    "vol_thickness_pbt",
]

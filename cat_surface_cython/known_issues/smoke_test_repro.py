"""
Smoke-test every public cat_surf function against a PPM/seg pair.

This is the reproducer for the get_area_normalized SIGABRT documented
in this directory's README.md.  All earlier functions pass; the abort
hits at the get_area_normalized step.

Usage
-----
    PYTHONPATH=/path/to/cat_surface_cython \\
        OMP_NUM_THREADS=1 \\
        python -u -X faulthandler smoke_test_repro.py \\
            --ppm /path/to/lh.ppm.<subject>.nii \\
            --seg /path/to/lh.seg.<subject>.nii
"""
import argparse
import os
import sys
import time
import tempfile
import traceback

import numpy as np
import nibabel as nib

os.environ.setdefault("OMP_NUM_THREADS", "1")

import cat_surf as cs  # noqa: E402

_ap = argparse.ArgumentParser()
_ap.add_argument("--ppm", required=True,
                 help="PPM_volume NIfTI (e.g. lh.ppm.<subj>.nii)")
_ap.add_argument("--seg", required=True,
                 help="Hemi_volume NIfTI (e.g. lh.seg.<subj>.nii)")
_args = _ap.parse_args()

PPM = _args.ppm
SEG = _args.seg

results = []


def run(name, fn):
    print(f"-> {name} ...", flush=True)
    t0 = time.perf_counter()
    try:
        rv = fn()
        dt = time.perf_counter() - t0
        line = f"PASS {name:<30} {dt:6.2f}s  {rv if rv else ''}"
        print(line, flush=True)
        results.append((name, "PASS", f"{dt:6.2f}s  {rv if rv else ''}"))
        return rv
    except Exception as exc:
        dt = time.perf_counter() - t0
        print(f"FAIL {name}: {type(exc).__name__}: {exc}", flush=True)
        results.append((name, "FAIL", f"{dt:6.2f}s  {type(exc).__name__}: {exc}"))
        traceback.print_exc()
        return None


def skip(name, reason):
    results.append((name, "SKIP", reason))


ppm_img = nib.load(PPM)
seg_img = nib.load(SEG)
ppm = ppm_img.get_fdata().astype(np.float32)
seg = seg_img.get_fdata().astype(np.float32)
vx = ppm_img.header.get_zooms()[:3]
print(f"PPM shape {ppm.shape}, dtype on disk {ppm_img.get_data_dtype()}, "
      f"voxelsize {vx}, min/max {ppm.min():.3f}/{ppm.max():.3f}")

tmp = tempfile.mkdtemp(prefix="cat_surf_test_")
print("tmp dir:", tmp)


# ---------------------------------------------------------------- volume
def t_sanlm():
    out = cs.vol_sanlm(seg, is_rician=False, strength=1.0)
    assert out.shape == seg.shape and out.dtype == np.float32
    return f"shape {out.shape}"


def t_bvc():
    out = cs.vol_blood_vessel_correction(seg, voxelsize=vx)
    return f"diff_voxels={int(((out!=seg).sum()))}"


def t_pbt():
    gmt, ppm2, dcsf, dwm = cs.vol_thickness_pbt(
        seg, voxelsize=vx, n_avgs=2, n_median_filter=2,
        median_subsample=2, range_val=0.45,
        correct_voxelsize=-0.75, sulcal_width=5.0,
    )
    return f"gmt mean {gmt.mean():.3f}, ppm>=0.5 frac {(ppm2>=0.5).mean():.3f}"


def t_amap():
    lab = (seg.round().clip(0, 3)).astype(np.uint8)
    prob, lab_out, mean = cs.vol_amap(
        seg, lab, voxelsize=vx,
        n_pure_classes=3, n_iters=5, sub=96, pve=True,
        weight_mrf=0.0, n_iters_icm=5,
    )
    return f"prob {prob.shape}, mean {mean}"


def t_mc():
    v, f = cs.vol_marching_cubes(
        PPM, label=SEG, threshold=0.5,
        pre_fwhm=1.0, iter_laplacian=50,
        n_median_filter=2, strength_gyri_mask=0.1,
    )
    return f"V={len(v)} F={len(f)}"


def t_mc_fast():
    v, f = cs.vol_marching_cubes(PPM, threshold=0.5, fast=True)
    return f"V={len(v)} F={len(f)}"


run("vol_sanlm",                    t_sanlm)
run("vol_blood_vessel_correction",  t_bvc)
run("vol_thickness_pbt",            t_pbt)
run("vol_amap",                     t_amap)
run("vol_marching_cubes",           t_mc)
run("vol_marching_cubes(fast)",     t_mc_fast)

# Single canonical MC mesh used for all surface tests below.
v_arr, f_arr = cs.vol_marching_cubes(
    PPM, label=SEG, threshold=0.5,
    pre_fwhm=1.0, iter_laplacian=50,
    n_median_filter=2, strength_gyri_mask=0.1,
)
print(f"\nMC surface: V={len(v_arr)} F={len(f_arr)}\n")


# ---------------------------------------------------------------- surface
def t_get_area():
    area, total = cs.get_area(v_arr, f_arr)
    return f"total_area={total:.1f} mm^2"


def t_euler():
    chi = cs.euler_characteristic(v_arr, f_arr)
    return f"chi={chi}"


def t_point_distance():
    d, maxd = cs.point_distance(v_arr, f_arr, v_arr, f_arr)
    return f"max={maxd:.4f}"


def t_point_distance_mean():
    d, meand = cs.point_distance_mean(v_arr, f_arr, v_arr, f_arr)
    return f"mean={meand:.4f}"


def t_hausdorff():
    d, hd = cs.hausdorff_distance(v_arr, f_arr, v_arr, f_arr)
    return f"hausdorff={hd:.4f}"


def t_smooth_heatkernel():
    vals = np.linalg.norm(v_arr, axis=1)
    out = cs.smooth_heatkernel(v_arr, f_arr, vals, fwhm=3.0)
    return f"shape={out.shape}"


def t_smoothed_curvatures():
    out = cs.smoothed_curvatures(v_arr, f_arr, fwhm=3.0, n_iter=10)
    return f"shape={out.shape}, mean={out.mean():.4f}"


def t_sulcus_depth():
    out = cs.sulcus_depth(v_arr, f_arr)
    return f"shape={out.shape}"


def t_reduce_mesh():
    target = max(100, int(0.25 * len(f_arr)))
    v2, f2 = cs.reduce_mesh(v_arr, f_arr, target_faces=target,
                            aggressiveness=7.0, preserve_sharp=True)
    return f"F {len(f_arr)}->{len(f2)}"


def t_sphere_radius():
    r = cs.sphere_radius(v_arr, f_arr)
    return f"r={r:.4f}"


def t_smooth_mesh():
    out_v, out_f = cs.smooth_mesh(v_arr, f_arr, iterations=5,
                                  alpha=0.5, beta=-0.53)
    return f"V shape={out_v.shape}"


def t_count_intersections():
    n = cs.count_intersections(v_arr, f_arr)
    return f"intersections={n}"


def t_remove_intersections():
    v2, f2 = cs.remove_intersections(v_arr, f_arr)
    return f"V {len(v_arr)}->{len(v2)}"


def t_correct_thickness_folding():
    thick = np.full(len(v_arr), 2.5, dtype=np.float32)
    out = cs.correct_thickness_folding(v_arr, f_arr, thick,
                                       slope=1.0, max_dist=6.0)
    return f"thick shape {out.shape}, mean {out.mean():.3f}"


run("get_area",                  t_get_area)
run("euler_characteristic",      t_euler)
run("point_distance",            t_point_distance)
run("point_distance_mean",       t_point_distance_mean)
run("hausdorff_distance",        t_hausdorff)
run("smooth_heatkernel",         t_smooth_heatkernel)
run("smoothed_curvatures",       t_smoothed_curvatures)
run("sulcus_depth",              t_sulcus_depth)
run("reduce_mesh",               t_reduce_mesh)
run("sphere_radius",             t_sphere_radius)
run("smooth_mesh",               t_smooth_mesh)
run("count_intersections",       t_count_intersections)
run("remove_intersections",      t_remove_intersections)
run("correct_thickness_folding", t_correct_thickness_folding)


# ---- spherical / multi-surface — heavy; run on smaller mesh
def t_surf_to_sphere():
    target = max(2000, int(0.10 * len(f_arr)))
    v2, f2 = cs.reduce_mesh(v_arr, f_arr, target_faces=target,
                            aggressiveness=7.0, preserve_sharp=True)
    sv, sf = cs.surf_to_sphere(v2, f2, stop_at=4)
    return f"sphere V={len(sv)} F={len(sf)}"


sphere = run("surf_to_sphere",  t_surf_to_sphere)


def t_get_area_normalized():
    # Reduce mesh + sphere of matching topology required.
    target = max(2000, int(0.10 * len(f_arr)))
    v2, f2 = cs.reduce_mesh(v_arr, f_arr, target_faces=target,
                            aggressiveness=7.0, preserve_sharp=True)
    sv, sf = cs.surf_to_sphere(v2, f2, stop_at=4)
    area, total = cs.get_area_normalized(v2, f2, sv, sf)
    return f"total={total:.2f}"


run("get_area_normalized",  t_get_area_normalized)
skip("resample_to_sphere", "needs source+target sphere pairs of matching topology")
skip("surf_average",       "needs >=2 input surfaces of matching topology")
skip("surf_to_pial_white", "needs thickness + matching PPM; exercised in T1Prep main loop")
skip("central_to_pial",    "needs thickness; exercised in T1Prep main loop")
skip("surf_deform",        "exercised separately in T1Prep via PPM_volume")
skip("surf_warp",          "needs source+target sphere/values pairs")


# ---------------------------------------------------------------- I/O
def t_write_read_surface():
    p = os.path.join(tmp, "test.gii")
    cs.write_surface(p, v_arr, f_arr)
    v2, f2 = cs.read_surface(p)
    assert v2.shape == v_arr.shape and f2.shape == f_arr.shape
    return f"roundtrip V={len(v2)} F={len(f2)}"


def t_write_read_values():
    p = os.path.join(tmp, "test.values")
    vals = np.arange(len(v_arr), dtype=np.float32) * 0.01
    cs.write_values(p, vals)
    v2 = cs.read_values(p)
    return f"len={len(v2)} dtype={v2.dtype}"


run("write/read_surface",  t_write_read_surface)
run("write/read_values",   t_write_read_values)


# ---------------------------------------------------------------- vol2surf / bbreg
def t_vol2surf():
    out = cs.vol2surf(PPM, v_arr, f_arr,
                      grid_start=-0.4, grid_end=0.4, grid_steps=5,
                      map_func="mean")
    # API returns (values, grid_positions)
    if isinstance(out, tuple):
        vals = out[0]
    else:
        vals = out
    return f"values shape={vals.shape}, mean={vals.mean():.4f}"


def t_bbreg_detect():
    contrast = cs.bbreg_detect_contrast(PPM, lh_surface=(v_arr, f_arr))
    return f"contrast={contrast}"


run("vol2surf",                 t_vol2surf)
run("bbreg_detect_contrast",    t_bbreg_detect)
skip("bbreg",                   "needs source volume + target volume + surface pair")
skip("volume_register_nmi",     "needs fixed+moving volume pair")
skip("volume_register_robust",  "needs fixed+moving volume pair")


# ---------------------------------------------------------------- report
w_name = max(len(n) for n, *_ in results)
print("\n" + "=" * (w_name + 50))
print(f"{'Function':<{w_name}}  Status  Detail")
print("=" * (w_name + 50))
for n, s, d in results:
    print(f"{n:<{w_name}}  {s:6}  {d}")

n_pass = sum(1 for _, s, _ in results if s == "PASS")
n_fail = sum(1 for _, s, _ in results if s == "FAIL")
n_skip = sum(1 for _, s, _ in results if s == "SKIP")
print(f"\n{n_pass} PASS  {n_fail} FAIL  {n_skip} SKIP")
sys.exit(1 if n_fail else 0)

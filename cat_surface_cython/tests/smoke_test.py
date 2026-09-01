"""
Build gate for the cat_surf Cython bindings.

The generated ``cat_surf/_*.c`` sources are build artifacts and are not
tracked (see ``.gitignore``), so nothing in the repository proves that the
``.pyx`` sources still cythonize, compile, link against libCAT, and import.
This script is that proof, and it runs in CI right after the extensions are
built.

It deliberately depends on numpy only -- no pytest, no nibabel -- so it can
run in the same environment the wheel is built in.

Usage
-----
    CAT_BUILD_DIR=/path/to/CAT-Surface pip install ./cat_surface_cython
    python cat_surface_cython/tests/smoke_test.py
"""
import os
import sys
import tempfile
import traceback

import numpy as np

import cat_surf
from cat_surf import cli

N = 24
VX = (1.0, 1.0, 1.0)

_failures = []


def check(name, condition, detail=""):
    """Record one assertion; keep going so a run reports every failure."""
    if condition:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name}{(' -- ' + detail) if detail else ''}")
        _failures.append(name)


def section(title):
    print(f"\n{title}")


# ---------------------------------------------------------------------------
# 1. Every extension module imports, and the advertised API resolves.
#    This is what catches a .pyx that no longer cythonizes or links.
# ---------------------------------------------------------------------------
def test_api_surface():
    section("API surface")

    missing = [n for n in cat_surf.__all__ if not hasattr(cat_surf, n)]
    check("every name in cat_surf.__all__ resolves", not missing, str(missing))

    missing = [n for n in cli.__all__ if not hasattr(cli, n)]
    check("every name in cat_surf.cli.__all__ resolves", not missing, str(missing))

    # Importing the package pulls in every extension module; name them
    # explicitly so a failure says which one broke.
    for mod in ("_io", "_surf", "_vol", "_convert", "_volume", "_bbreg",
                "_vol2surf", "_surf_warp", "_spherical_demon"):
        try:
            __import__(f"cat_surf.{mod}")
            ok = True
        except Exception as exc:  # pragma: no cover
            ok = False
            print(f"       {exc}")
        check(f"cat_surf.{mod} imports", ok)


# ---------------------------------------------------------------------------
# 2. Sheetness separates a thin sheet from a blob and orients across it.
# ---------------------------------------------------------------------------
def test_sheetness():
    section("Hessian sheetness")

    vol = np.full((N, N, N), 2.0, np.float32)
    vol[8, :, :] = 1.0                      # thin dark sheet
    vol[15:18, 15:18, 15:18] = 1.0          # dark blob

    S, nrm = cat_surf.vol_sheetness(vol, voxelsize=VX, polarity=-1,
                                    sigma_min=0.5, sigma_max=1.0, n_scales=2,
                                    return_normal=True)

    check("response has the input shape", S.shape == (N, N, N))
    check("normals are (X, Y, Z, 3)", nrm.shape == (N, N, N, 3))
    check("thin dark sheet responds", S[8, 12, 12] > 0.5, f"{S[8, 12, 12]:.3f}")
    check("dark blob is rejected", S[16, 16, 16] < 0.2, f"{S[16, 16, 16]:.3f}")
    check("normal points across the sheet",
          abs(float(nrm[8, 12, 12, 0])) > 0.95)

    S_bright = cat_surf.vol_sheetness(vol, voxelsize=VX, polarity=1,
                                      sigma_min=0.5, sigma_max=1.0, n_scales=2)
    check("opposite polarity ignores a dark sheet", S_bright[8, 12, 12] < 0.05)


# ---------------------------------------------------------------------------
# 3. The oriented filters keep what the isotropic ones destroy -- and are
#    exact no-ops where no sheet was found.  That invariant is what makes
#    every -oriented / -oriented-filter option safe to enable,
#    so it is guarded here as well as in tests/test_sheetness.c.
# ---------------------------------------------------------------------------
def test_oriented_filters():
    section("Oriented filters")

    vol = np.full((N, N, N), 2.0, np.float32)
    vol[12, :, :] = 1.0

    ones = np.ones((N, N, N), np.float32)
    zeros = np.zeros((N, N, N), np.float32)
    nrm = np.zeros((N, N, N, 3), np.float32)
    nrm[..., 0] = 1.0

    oriented = cat_surf.vol_oriented_median(vol, sheetness=ones, normal=nrm)
    isotropic = cat_surf.vol_oriented_median(vol, sheetness=zeros, normal=nrm)

    check("isotropic median erases the sheet", isotropic[12, 12, 12] > 1.9)
    check("oriented median keeps the sheet", oriented[12, 12, 12] < 1.1)

    plain = cat_surf.vol_oriented_median(vol, sheetness=zeros, normal=nrm)
    check("median at zero sheetness is bit-identical to isotropic",
          np.array_equal(plain, isotropic))

    # The C side reads normals as normal[3*i + c] with i a Fortran index;
    # a round trip through a real estimated field checks that packing.
    _, est = cat_surf.vol_sheetness(vol, voxelsize=VX, sigma_min=0.5,
                                    sigma_max=1.0, n_scales=2,
                                    return_normal=True)
    out = cat_surf.vol_oriented_median(vol, sheetness=ones, normal=est)
    check("estimated normals round-trip through the filter",
          out.shape == vol.shape and np.isfinite(out).all())


# ---------------------------------------------------------------------------
# 4. Pre-PBT repair: opens a glued sulcus, bridges a broken blade, and does
#    neither where the evidence is absent.
# ---------------------------------------------------------------------------
def test_open_ppm_sulci():
    section("Buried sulci in the PPM")

    # Two one-voxel structures with the same amplitudes and opposite curvature:
    # a sulcal valley whose floor stops short of the isovalue, and a gyral
    # ridge.  A filter keying on "thin sheet" alone would cut both.
    ppm = np.full((N, N, N), 0.95, np.float32)
    ppm[8, :, :] = 0.60
    ppm[15, :, :] = 0.60
    ppm[16, :, :] = 0.95
    ppm[17, :, :] = 0.60

    out = cat_surf.vol_open_ppm_sulci(ppm, voxelsize=VX, isovalue=0.5,
                                      sigma_min=0.5, sigma_max=1.0, n_scales=2)

    check("buried valley is pushed below the isovalue", out[8, 12, 12] < 0.5,
          f"{out[8, 12, 12]:.3f}")
    check("gyral crest is left alone", out[16, 12, 12] == ppm[16, 12, 12])
    check("solid values are left alone", out[4, 12, 12] == ppm[4, 12, 12])
    check("no value is ever raised", bool((out <= ppm + 1e-7).all()))

    noop = cat_surf.vol_open_ppm_sulci(ppm, voxelsize=VX, isovalue=0.5,
                                       strength=0.0)
    check("strength 0 is an exact no-op", np.array_equal(noop, ppm))


def test_marching_cubes_sulci_kwargs():
    section("Marching-cubes buried-sulcus kwargs")

    import inspect
    params = inspect.signature(cat_surf.vol_marching_cubes).parameters
    for k in ("strength_sulci", "sulci_cutoff", "sulci_sheet_strength",
              "sulci_thresh", "sulci_band", "sulci_normalize", "sulci_skeleton"):
        check(f"vol_marching_cubes accepts {k}", k in params)

    # A PPM-like blob with a buried valley: the correction must run end to end
    # with the skeleton both off and on, and leave the mesh usable either way.
    zz, yy, xx = np.mgrid[0:N, 0:N, 0:N]
    r = np.sqrt((xx - N / 2.) ** 2 + (yy - N / 2.) ** 2 + (zz - N / 2.) ** 2)
    ppm = np.zeros((N, N, N), np.float32)
    ppm[r < N / 3.] = 0.95
    ppm[(np.abs(xx - N // 2) <= 0) & (r < N / 4.)] = 0.62
    aff = np.eye(4)

    for flag in (False, True):
        try:
            v, f = cat_surf.vol_marching_cubes((ppm, aff), threshold=0.5,
                                               strength_sulci=1.0,
                                               sulci_skeleton=flag)
            ok = len(v) > 0 and len(f) > 0
        except Exception as exc:  # pragma: no cover
            ok = False
            print(f"       {exc}")
        check(f"marching cubes runs with sulci_skeleton={flag}", ok)


def test_sulcal_barrier():
    section("Sulcal barrier in PBT")

    M = 64
    def build(gap_is_csf):
        v = np.ones((M, M, M), np.float32)
        v[24:40, 10:54, 10:54] = 2.0
        v[16:24, 10:54, 10:54] = 3.0
        v[40:48, 10:54, 10:54] = 3.0
        v[10:16, 10:54, 10:54] = 1.0
        v[48:54, 10:54, 10:54] = 1.0
        if gap_is_csf:
            v[30:34, 10:54, 10:54] = 1.0
        return v

    vx = (0.5, 0.5, 0.5)

    def run(vol, barrier, gmtmax=0.0):
        g, ppm, _, _ = cat_surf.vol_thickness_pbt(vol, voxelsize=vx,
                                                  sulcal_barrier=barrier,
                                                  barrier_gmtmax=gmtmax)
        prof = ppm[24:41, 32, 32]
        xs = [24 + i for i in range(len(prof) - 1)
              if (prof[i] - 0.5) * (prof[i + 1] - 0.5) < 0]
        return g, ppm, xs

    # A properly segmented sulcus splits the cortical band, so no two fronts
    # ever meet and there is nothing for the barrier to do.  That has to be an
    # exact no-op, not merely a small one -- it is what makes the option safe
    # to leave on for every subject.
    g_off, p_off, x_off = run(build(True), False)
    g_on, p_on, x_on = run(build(True), True)
    check("a segmented sulcus is left bit-identical",
          np.array_equal(g_off, g_on) and np.array_equal(p_off, p_on))

    # With the CSF lost, the distance runs across the fused grey matter and the
    # thickness follows it; the barrier has to put the missing boundary back.
    # The gate is pinned here rather than derived.  This phantom is fused
    # everywhere, so the median thickness it would derive the gate from is
    # itself measuring the fused band -- the estimator assumes glued sulci are
    # a minority, which holds on a brain (1.6% of the cortex on the subject it
    # was checked against) but not on a phantom built entirely out of them.
    gf_off, pf_off, xf_off = run(build(False), False)
    gf_on, pf_on, xf_on = run(build(False), True, gmtmax=5.0)

    check("fusing the banks overestimates the thickness",
          gf_off[27, 32, 32] > g_off[27, 32, 32] + 1.0,
          f"{gf_off[27, 32, 32]:.2f} vs {g_off[27, 32, 32]:.2f}")
    check("the barrier recovers it",
          abs(gf_on[27, 32, 32] - g_off[27, 32, 32]) < 0.3,
          f"{gf_on[27, 32, 32]:.2f} vs {g_off[27, 32, 32]:.2f}")
    check("and never inflates a thickness", bool((gf_on <= gf_off + 1e-4).all()))
    check("the central surfaces move back apart",
          len(xf_on) == 2 and len(xf_off) == 2 and
          (xf_on[1] - xf_on[0]) > (xf_off[1] - xf_off[0]),
          f"{xf_on} vs {xf_off}")


def test_barrier_gate_scales_with_thickness():
    section("Barrier gate follows the cortex")

    M = 64

    def banks(bank, fused):
        """Two banks of `bank` voxels each, fused or properly separated."""
        v = np.ones((M, M, M), np.float32)
        mid = 32
        a0, b1 = mid - 2 * bank, mid + 2 * bank
        v[a0:b1, 10:54, 10:54] = 2.0
        v[a0 - 8:a0, 10:54, 10:54] = 3.0
        v[b1:b1 + 8, 10:54, 10:54] = 3.0
        v[a0 - 14:a0 - 8, 10:54, 10:54] = 1.0
        v[b1 + 8:b1 + 14, 10:54, 10:54] = 1.0
        if not fused:
            v[mid - 1:mid + 1, 10:54, 10:54] = 1.0
        return v

    vx = (0.5, 0.5, 0.5)

    # The same factor has to work on both a thin and a thick cortex, because
    # the threshold it produces is derived from each one separately.  A fixed
    # millimetre gate can only be right for the thickness it was tuned on.
    for bank, label in ((3, "thin"), (6, "thick")):
        g_ok, _, _, _ = cat_surf.vol_thickness_pbt(banks(bank, False), voxelsize=vx)
        g_off, _, _, _ = cat_surf.vol_thickness_pbt(banks(bank, True), voxelsize=vx)
        g_on, _, _, _ = cat_surf.vol_thickness_pbt(banks(bank, True), voxelsize=vx,
                                                   sulcal_barrier=True)
        ref = g_ok[g_ok > 0.3].mean()
        bad = g_off[g_off > 0.3].mean()
        fixed = g_on[g_on > 0.3].mean()
        check(f"{label} cortex: fusing overestimates", bad > ref + 0.2,
              f"{bad:.2f} vs {ref:.2f}")
        check(f"{label} cortex: the derived gate still catches it",
              abs(fixed - ref) < abs(bad - ref),
              f"{fixed:.2f} -> ref {ref:.2f}, was {bad:.2f}")

    # and a properly segmented cortex is untouched whatever its thickness
    for bank, label in ((3, "thin"), (6, "thick")):
        v = banks(bank, False)
        g0, p0, _, _ = cat_surf.vol_thickness_pbt(v, voxelsize=vx)
        g1, p1, _, _ = cat_surf.vol_thickness_pbt(v, voxelsize=vx, sulcal_barrier=True)
        check(f"{label} cortex: segmented case is bit-identical",
              np.array_equal(g0, g1) and np.array_equal(p0, p1))


def test_option_no_ops():
    section("Backward compatibility of the new options")

    x = np.arange(N)[:, None, None]
    lab = np.where((x < 6) | (x > 17), 3.0, 2.0).astype(np.float32) \
        * np.ones((N, N, N), np.float32)
    lab[0, :, :] = 1.0

    g1, p1, _, _ = cat_surf.vol_thickness_pbt(lab, voxelsize=VX, fast=True,
                                              oriented_filter=False)
    g2, p2, _, _ = cat_surf.vol_thickness_pbt(lab, voxelsize=VX, fast=True,
                                              oriented_filter=True,
                                              oriented_strength=0.0)
    check("PBT at oriented_strength=0 is bit-identical to isotropic",
          np.array_equal(g1, g2) and np.array_equal(p1, p2))

    # The gain has to be able to amplify, not just attenuate: every consumer
    # gates on the response through a threshold, and the automatic noise scale
    # routinely leaves a real sulcus an order of magnitude below the one in
    # vol_oriented_median, which is hard at 0.5.
    sheet_vol = np.full((N, N, N), 2.0, np.float32)
    sheet_vol[12, :, :] = 1.0
    s_low = cat_surf.vol_sheetness(sheet_vol, voxelsize=VX, polarity=-1,
                                   sigma_min=0.5, sigma_max=1.0, n_scales=2,
                                   gain=0.25)
    s_high = cat_surf.vol_sheetness(sheet_vol, voxelsize=VX, polarity=-1,
                                    sigma_min=0.5, sigma_max=1.0, n_scales=2,
                                    gain=0.5)
    unclamped = s_high < 1.0
    check("sheetness gain scales the response linearly",
          np.allclose(s_high[unclamped], 2.0 * s_low[unclamped], atol=1e-5))
    check("sheetness gain never leaves [0, 1]",
          s_high.min() >= 0.0 and s_high.max() <= 1.0)
    check("a zero response stays zero at any gain",
          np.array_equal(s_high[s_low == 0.0],
                         np.zeros_like(s_high[s_low == 0.0])))

    s_up = cat_surf.vol_sheetness(sheet_vol, voxelsize=VX, polarity=-1,
                                  sigma_min=0.5, sigma_max=1.0, n_scales=2,
                                  gain=4.0)
    check("a gain above 1 lifts the response over the median threshold",
          s_low.max() <= 0.5 and s_up.max() > 0.5)


# ---------------------------------------------------------------------------
# Surface info: every reported quantity has a closed-form value on a regular
# octahedron, and the GIFTI listing has to survive a round-trip through disk.
# ---------------------------------------------------------------------------
def test_surf_info():
    section("Surface info")

    pts = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                    [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.float64)
    tri = np.array([[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
                    [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]],
                   dtype=np.int32)

    info = cat_surf.surf_info(pts, tri)

    check("V, F and E of an octahedron",
          info["n_points"] == 6 and info["n_polygons"] == 8
          and info["n_edges"] == 12)
    check("V - E + F = 2", info["euler"] == 2)
    check("genus 0, one component, closed",
          info["genus"] == 0 and info["n_components"] == 1 and info["closed"])
    check("no boundary or non-manifold edges",
          info["n_boundary_edges"] == 0 and info["n_nonmanifold_edges"] == 0)
    check("surface area is 4*sqrt(3)",
          abs(info["surface_area"] - 4.0 * np.sqrt(3.0)) < 1e-6,
          f"{info['surface_area']:.6f}")
    check("volume is 4/3", abs(info["volume"] - 4.0 / 3.0) < 1e-6,
          f"{info['volume']:.6f}")
    check("all edges have length sqrt(2)",
          abs(info["edge_mean"] - np.sqrt(2.0)) < 1e-6
          and abs(info["edge_max"] - info["edge_min"]) < 1e-9)
    check("centroid at the origin",
          np.allclose(info["centroid"], 0.0, atol=1e-9))
    check("a convex body has no self-intersections",
          info["n_self_intersections"] == 0)

    skipped = cat_surf.surf_info(pts, tri, check_intersections=False)
    check("a skipped test reports None, not zero",
          skipped["n_self_intersections"] is None
          and skipped["n_triangle_hits"] is None
          and skipped["n_intersecting_polygons"] is None)

    # an open patch: two triangles sharing one edge
    patch_pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                         dtype=np.float64)
    patch_tri = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
    patch = cat_surf.surf_info(patch_pts, patch_tri,
                               check_intersections=False)
    check("chi of a disc is 1", patch["euler"] == 1)
    check("the patch is open, with four boundary edges",
          not patch["closed"] and patch["n_boundary_edges"] == 4)
    check("genus and volume are None when the surface is open",
          patch["genus"] is None and patch["volume"] is None)
    check("area of the unit square",
          abs(patch["surface_area"] - 1.0) < 1e-9)

    # a .gii carries the mesh in two DataArrays and may embed data next to it
    tmp = os.path.join(tempfile.gettempdir(), "cat_surf_smoke_info.gii")
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0, np.nan])
    try:
        cat_surf.write_surface(tmp, pts, tri, values=values)
        cli_info = cat_surf.cli.surf_info(tmp, check_intersections=False)

        check("the CLI mirror reproduces the array API",
              cli_info["euler"] == 2 and cli_info["n_edges"] == 12
              and abs(cli_info["surface_area"] - info["surface_area"]) < 1e-6)

        darrays = cli_info["gifti_darrays"]
        check("the embedded array is listed next to the mesh",
              darrays is not None and len(darrays) == 3, str(darrays and len(darrays)))

        shape = [d for d in (darrays or [])
                 if d["intent_name"] == "NIFTI_INTENT_SHAPE"]
        check("the shape array is found", len(shape) == 1)
        if shape:
            da = shape[0]
            check("one value per vertex", da["n_values"] == 6)
            check("the NaN is counted, not silently dropped",
                  da["n_nonfinite"] == 1)
            check("statistics skip the NaN instead of becoming nan",
                  da["mean"] is not None and abs(da["mean"] - 3.0) < 1e-5
                  and abs(da["min"] - 1.0) < 1e-5
                  and abs(da["max"] - 5.0) < 1e-5,
                  str(da["mean"]))
        mesh_arrays = [d for d in (darrays or [])
                       if d["intent_name"] in ("NIFTI_INTENT_POINTSET",
                                               "NIFTI_INTENT_TRIANGLE")]
        check("the mesh arrays carry no statistics",
              len(mesh_arrays) == 2
              and all(d["mean"] is None for d in mesh_arrays))
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


# ---------------------------------------------------------------------------
# 6. Thickness QC tells a plate from a solid mass, which is the whole point:
#    a glued sulcus is recoverable, a subcortical mass is not, and thickness
#    alone cannot separate them.
# ---------------------------------------------------------------------------
def main():
    print(f"cat_surf {cat_surf.__version__} — binding smoke test")
    for test in (test_api_surface, test_sheetness, test_oriented_filters,
                 test_open_ppm_sulci, test_marching_cubes_sulci_kwargs, test_sulcal_barrier,
                 test_barrier_gate_scales_with_thickness, test_option_no_ops,
                 test_surf_info):
        try:
            test()
        except Exception:  # pragma: no cover
            print(f"  FAIL {test.__name__} raised:")
            traceback.print_exc()
            _failures.append(test.__name__)

    print()
    if _failures:
        print(f"FAILED ({len(_failures)}): {', '.join(_failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

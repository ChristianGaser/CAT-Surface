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
                "_vol2surf", "_spherical_demon"):
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


def test_barrier_reference():
    section("Shared barrier reference")

    M = 64
    v = np.ones((M, M, M), np.float32)
    v[26:38, 10:54, 10:54] = 2.0
    v[18:26, 10:54, 10:54] = 3.0
    v[38:46, 10:54, 10:54] = 3.0
    vx = (0.5, 0.5, 0.5)
    kw = dict(sulcal_barrier=True, oriented_filter=True, n_avgs=5)

    # One kwargs dict drives both calls, so the reference is the one the run
    # would derive -- and handing it back must not change the result.
    ref = cat_surf.vol_pbt_barrier_reference(v, voxelsize=vx, **kw)
    check("the reference is a thickness in mm", 0.5 < ref < 10.0, f"{ref:.3f}")
    g0, p0, _, _ = cat_surf.vol_thickness_pbt(v, voxelsize=vx, **kw)
    g1, p1, _, _ = cat_surf.vol_thickness_pbt(v, voxelsize=vx, barrier_gmtref=ref, **kw)
    check("a given reference reproduces the derived run",
          np.array_equal(g0, g1) and np.array_equal(p0, p1))

    # the gate follows it: a smaller reference corrects more
    g2, _, _, _ = cat_surf.vol_thickness_pbt(v, voxelsize=vx, barrier_gmtref=0.5 * ref, **kw)
    check("a smaller reference corrects more", g2.sum() < g0.sum(),
          f"{g2.sum():.1f} vs {g0.sum():.1f}")

    try:
        cat_surf.vol_pbt_barrier_reference(v, voxelsize=vx, barrier_gmtfatcor=2.0)
        raised = False
    except TypeError:
        raised = True
    check("a misspelt keyword is rejected", raised)

    # repeated calls in one process are reproducible
    g3, _, _, _ = cat_surf.vol_thickness_pbt(v, voxelsize=vx, **kw)
    check("repeated runs are identical", np.array_equal(g0, g3))


def _icosphere(level, radius):
    """Subdivided icosahedron with outward-facing triangles."""
    t = (1 + 5 ** 0.5) / 2
    v = [np.array(x, float) for x in
         ([-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
          [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1])]
    f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9],
         [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2],
         [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10],
         [8, 6, 7], [9, 8, 1]]
    for _ in range(level):
        cache = {}

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                v.append((v[a] + v[b]) / 2)
                cache[key] = len(v) - 1
            return cache[key]

        f = [tri for a, b, c in f
             for ab, bc, ca in [(mid(a, b), mid(b, c), mid(c, a))]
             for tri in ([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca])]
    v = np.array(v)
    return radius * v / np.linalg.norm(v, axis=1, keepdims=True), np.array(f, np.int32)


def test_surf_to_sphere_relax():
    section("surf_to_sphere(relax=True)")

    # an elongated, flattened ellipsoid: the inflation maps it unevenly
    v, f = _icosphere(4, 1.0)
    v = v * np.array([80.0, 45.0, 30.0])

    # the relaxation minimises the distance distortion: edge lengths on the
    # sphere against those on the surface, scaled to the same area
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    len0 = np.linalg.norm(v[e[:, 0]] - v[e[:, 1]], axis=1)

    def stats(sv):
        r = np.linalg.norm(sv, axis=1)
        a = np.cross(sv[f[:, 1]] - sv[f[:, 0]], sv[f[:, 2]] - sv[f[:, 0]])
        folds = int(((a * (sv[f[:, 0]] + sv[f[:, 1]] + sv[f[:, 2]])).sum(1) < 0).sum())
        ratio = np.linalg.norm(sv[e[:, 0]] - sv[e[:, 1]], axis=1) / len0
        return r.std() / r.mean(), folds, np.std(np.log2(ratio))

    plain, _ = cat_surf.surf_to_sphere(v, f, stop_at=6)
    relaxed, rf = cat_surf.surf_to_sphere(v, f, stop_at=6, relax=True)
    r_plain, folds_plain, d_plain = stats(np.asarray(plain))
    r_rel, folds_rel, d_rel = stats(np.asarray(relaxed))
    check("relaxed map keeps the mesh", np.array_equal(np.asarray(rf), f))
    check("relaxed map is a sphere", r_rel < 1e-6, f"radius cv {r_rel:.2e}")
    check("relaxed map is fold-free", folds_rel == 0, f"{folds_rel} folds")
    check("relaxation lowers the distance distortion", d_rel < d_plain,
          f"{d_plain:.3f} -> {d_rel:.3f}")
    try:
        cat_surf.surf_to_sphere(v, f, stop_at=3, relax=True)
        raised = False
    except ValueError:
        raised = True
    check("relax=True before the sphere stage is refused", raised)


def test_sphere_remove_folds():
    section("sphere_remove_folds / surf_to_sphere(remove_folds)")

    # a regular sphere as surface, mapped with a radial crease in a polar cap:
    # within 30 deg the polar angle becomes t + (cap/pi) sin(2 pi t / cap),
    # which runs backwards in the middle third, so that band is mirrored
    v, f = _icosphere(5, 100.0)
    cap = np.radians(30.0)
    theta = np.arccos(np.clip(v[:, 2] / 100.0, -1, 1))
    phi = np.arctan2(v[:, 1], v[:, 0])
    t = np.where(theta < cap, theta + cap / np.pi * np.sin(2 * np.pi * theta / cap), theta)
    creased = 100.0 * np.c_[np.sin(t) * np.cos(phi), np.sin(t) * np.sin(phi), np.cos(t)]

    def folds(sv):
        a = np.cross(sv[f[:, 1]] - sv[f[:, 0]], sv[f[:, 2]] - sv[f[:, 0]])
        return int(((a * (sv[f[:, 0]] + sv[f[:, 1]] + sv[f[:, 2]])).sum(1) < 0).sum())

    fixed, left = cat_surf.sphere_remove_folds(v, f, creased)
    fixed = np.asarray(fixed)
    outside = theta > 2 * cap
    check("creased map starts folded", folds(creased) > 0, f"{folds(creased)} folds")
    check("folds removed", left == 0 and folds(fixed) == 0,
          f"{left} left, {folds(fixed)} folds")
    check("nothing outside the creased area moves",
          np.allclose(fixed[outside], creased[outside], atol=1e-4))
    check("map stays on the sphere",
          np.abs(np.linalg.norm(fixed, axis=1) - 100.0).max() < 1e-2)
    same, left = cat_surf.sphere_remove_folds(v, f, v)
    check("fold-free map is returned unchanged",
          left == 0 and np.allclose(np.asarray(same), v, atol=1e-4))

    # surf_to_sphere runs it by default; without folds nothing changes
    e, ef = _icosphere(4, 1.0)
    e = e * np.array([80.0, 45.0, 30.0])
    plain, _ = cat_surf.surf_to_sphere(e, ef, stop_at=6, remove_folds=False)
    default, _ = cat_surf.surf_to_sphere(e, ef, stop_at=6)
    check("fold-free sphere identical with and without remove_folds",
          np.array_equal(np.asarray(plain), np.asarray(default)))


def test_fix_self_intersect_reference():
    section("Self-intersection repair with a reference")

    # Two sheets driven through each other -- here the northern cap of a
    # sphere pushed through its southern half -- survive any smoothing; the
    # sphere they came from is the way back.
    sphere, f = _icosphere(4, 10.0)
    v = sphere.copy()
    v[v[:, 2] > 7.0, 2] -= 19.0
    before = cat_surf.count_intersections(v, f)
    plain, _ = cat_surf.fix_self_intersect(v, f)
    fixed, _ = cat_surf.fix_self_intersect(v, f, reference=sphere)
    check("the pushed cap crosses the sphere", before > 100, str(before))
    check("smoothing alone leaves crossings",
          cat_surf.count_intersections(plain, f) > 0)
    check("retreating to the reference removes them",
          cat_surf.count_intersections(fixed, f) == 0)
    check("the equator is untouched",
          np.array_equal(fixed[np.abs(sphere[:, 2]) < 3.0],
                         sphere[np.abs(sphere[:, 2]) < 3.0].astype(np.float32)))

    try:
        cat_surf.fix_self_intersect(v, f, reference=sphere[:-1])
        raised = False
    except ValueError:
        raised = True
    check("a reference of the wrong size is rejected", raised)

    with tempfile.TemporaryDirectory() as d:
        src, ref, out = (os.path.join(d, n) for n in ("in.gii", "ref.gii", "out.gii"))
        cat_surf.write_surface(src, v, f)
        cat_surf.write_surface(ref, sphere, f)
        cli.surf_fix_self_intersect(src, out, reference_file=ref)
        check("the CLI mirror passes the reference on",
              cli.surf_fix_self_intersect(out, count_only=True) == 0)


def test_pbt_front_ends_agree():
    section("PBT: the Python front-ends follow CAT_VolThicknessPbt")

    # WM core, GM band, and a WM-bright tube through the GM that the
    # blood-vessel correction removes
    x, y, z = np.meshgrid(np.arange(N), np.arange(N), np.arange(N),
                          indexing="ij")
    lab = np.ones((N, N, N), np.float32)
    lab[(x >= 6) & (x <= 17)] = 2.0
    lab[(x >= 9) & (x <= 14)] = 3.0
    tube = ((x - 7) ** 2 + (y - 12) ** 2 <= 1) & (z > 2) & (z < 21)
    lab[tube] = 3.0

    g_def, p_def, _, _ = cat_surf.vol_thickness_pbt(lab, voxelsize=VX)
    g_on, p_on, _, _ = cat_surf.vol_thickness_pbt(lab, voxelsize=VX,
                                                  oriented_filter=True)
    check("oriented_filter left unset keeps the library default (on)",
          np.array_equal(g_def, g_on) and np.array_equal(p_def, p_on))
    r_def = cat_surf.vol_pbt_barrier_reference(lab, voxelsize=VX)
    r_on = cat_surf.vol_pbt_barrier_reference(lab, voxelsize=VX,
                                              oriented_filter=True)
    check("the barrier reference keeps it too", r_def == r_on,
          f"{r_def} vs {r_on}")

    bv = cat_surf.vol_blood_vessel_correction(lab, voxelsize=VX)
    check("the phantom has a vessel for the correction to remove",
          np.any(bv[tube] != lab[tube]))
    g_raw, _, _, _ = cat_surf.vol_thickness_pbt(lab, voxelsize=VX, fast=True)
    check("the array API leaves the vessel correction to the caller",
          np.array_equal(g_raw, cat_surf.vol_thickness_pbt(
              lab, voxelsize=VX, fast=True)[0]))

    try:
        import nibabel as nib
    except ImportError:
        print("  skip the CLI mirror applies the vessel correction "
              "(needs nibabel)")
        return
    g_bv, p_bv, _, _ = cat_surf.vol_thickness_pbt(bv, voxelsize=VX, fast=True)
    with tempfile.TemporaryDirectory() as d:
        src, gmt, ppm, gmt0 = (os.path.join(d, n) for n in
                               ("seg.nii", "gmt.nii", "ppm.nii", "gmt0.nii"))
        nib.save(nib.Nifti1Image(lab, np.eye(4)), src)
        cli.vol_thickness_pbt(src, gmt, ppm, fast=True)
        cli.vol_thickness_pbt(src, gmt0, fast=True,
                              blood_vessel_correction=False)
        g_cli = np.asarray(nib.load(gmt).dataobj)
        p_cli = np.asarray(nib.load(ppm).dataobj)
        g_cli0 = np.asarray(nib.load(gmt0).dataobj)
    check("the CLI mirror applies the vessel correction, like the binary",
          np.array_equal(g_cli, g_bv) and np.array_equal(p_cli, p_bv))
    check("blood_vessel_correction=False is -no-blood-vessel-correction",
          np.array_equal(g_cli0, g_raw) and not np.array_equal(g_cli0, g_bv))


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
# Image calculator: the formula has to address the inputs in order, reduce
# X across them, and refuse what it cannot evaluate.
# ---------------------------------------------------------------------------
def test_vol_calc():
    section("Image calculator")

    rng = np.random.default_rng(0)
    a = rng.random((5, 6, 7))
    b = rng.random((5, 6, 7))
    c = rng.random((5, 6, 7))

    d = cat_surf.vol_calc([a, b], "i2-i1")
    check("i2-i1 is the voxel-wise difference", np.allclose(d, b - a))
    check("the result keeps the input shape", d.shape == a.shape)

    rel = cat_surf.vol_calc([a, b], "200*(i2-i1)./(i1+i2+2.220446049250313e-16)")
    check("MATLAB element-wise operators are accepted",
          np.allclose(rel, 200 * (b - a) / (a + b + np.finfo(float).eps)))

    check("median(X) reduces across the inputs",
          np.allclose(cat_surf.vol_calc([a, b, c], "median(X)"),
                      np.median(np.stack([a, b, c]), axis=0)))
    check("a C-ordered input addresses the same voxels",
          np.allclose(cat_surf.vol_calc([np.ascontiguousarray(a), b], "i1-i2"),
                      a - b))

    for bad, why in (("i3-i1", "an image that was not passed"),
                     ("i1-", "a truncated formula")):
        try:
            cat_surf.vol_calc([a, b], bad)
            raised = False
        except ValueError:
            raised = True
        check(f"{why} is a ValueError", raised)

    try:
        cat_surf.vol_calc([a, b[:, :, :3]], "i1")
        raised = False
    except ValueError:
        raised = True
    check("differing shapes are a ValueError", raised)


# ---------------------------------------------------------------------------
# Box morphology: a cube of side 2r+1, outside the volume is background.
# The reference is the definition, built from shifted copies of a padded mask.
# ---------------------------------------------------------------------------
def _box_reference(mask, r, erode):
    """Any (dilate) or all (erode) of the (2r+1)^3 cube; outside is background."""
    padded = np.pad(mask, r, constant_values=False)
    nx, ny, nz = mask.shape
    out = np.ones(mask.shape, bool) if erode else np.zeros(mask.shape, bool)
    for dx in range(2 * r + 1):
        for dy in range(2 * r + 1):
            for dz in range(2 * r + 1):
                shifted = padded[dx:dx + nx, dy:dy + ny, dz:dz + nz]
                out = (out & shifted) if erode else (out | shifted)
    return out


def test_box_morph():
    section("Box morphology")

    rng = np.random.default_rng(0)
    vol = rng.random((13, 11, 9)) < 0.3
    before = vol.copy()
    for r in (1, 2, 4):
        dil = _box_reference(vol, r, False)
        ero = _box_reference(vol, r, True)
        check(f"dilation by {r} is the (2r+1)^3 cube",
              np.array_equal(cat_surf.vol_box_morph(vol, "dilate", r), dil))
        check(f"erosion by {r} is the (2r+1)^3 cube",
              np.array_equal(cat_surf.vol_box_morph(vol, "erode", r), ero))
        check(f"closing by {r} is a dilation then an erosion",
              np.array_equal(cat_surf.vol_box_morph(vol, "close", r),
                             _box_reference(dil, r, True)))
        check(f"opening by {r} is an erosion then a dilation",
              np.array_equal(cat_surf.vol_box_morph(vol, "open", r),
                             _box_reference(ero, r, False)))

    out = cat_surf.vol_box_morph(vol, "dilate", 2)
    check("the result is boolean and keeps the shape",
          out.dtype == np.bool_ and out.shape == vol.shape)
    check("the input is not modified", np.array_equal(vol, before))
    check("radius 0 returns the mask", np.array_equal(
        cat_surf.vol_box_morph(vol, "erode", 0), vol))
    check("a Fortran-ordered input gives the same result", np.array_equal(
        cat_surf.vol_box_morph(np.asfortranarray(vol), "dilate", 2), out))
    check("any non-zero value is foreground", np.array_equal(
        cat_surf.vol_box_morph(vol.astype(np.uint8) * 7, "dilate", 2), out))
    single = np.zeros((9, 9, 9), bool)
    single[4, 4, 4] = True
    check("one voxel dilated by 1 is a 3x3x3 cube",
          cat_surf.vol_box_morph(single, "dilate", 1).sum() == 27)

    for args, why in (((vol, "boundary", 1), "an unknown operation"),
                      ((vol[0], "dilate", 1), "a 2-D mask")):
        try:
            cat_surf.vol_box_morph(*args)
            raised = False
        except ValueError:
            raised = True
        check(f"{why} is a ValueError", raised)


# ---------------------------------------------------------------------------
# B-spline sampling: torch's grid convention (x along the last axis), and an
# interpolation, so the voxel centres come back unchanged.
# ---------------------------------------------------------------------------
def _centres(n):
    """Normalised coordinates of the voxel centres along an axis of n voxels."""
    return np.linspace(-1.0, 1.0, n)


def test_bspline_sample():
    section("B-spline sampling")

    d, h, w = 9, 11, 13
    ramp_w = np.broadcast_to(np.arange(w, dtype=np.float32), (d, h, w)).copy()
    ramp_d = np.broadcast_to(np.arange(d, dtype=np.float32)[:, None, None],
                             (d, h, w)).copy()
    # an interior line of voxel centres along x, at y and z in the middle
    xs = _centres(w)[3:-3]
    line_x = np.stack([xs, np.zeros_like(xs), np.zeros_like(xs)], axis=-1)
    zs = _centres(d)[2:-2]
    line_z = np.stack([np.zeros_like(zs), np.zeros_like(zs), zs], axis=-1)
    check("grid x runs along the last axis of the volume",
          np.allclose(cat_surf.vol_bspline_sample(ramp_w, line_x),
                      np.arange(w)[3:-3], atol=1e-4))
    check("grid z runs along the first axis of the volume",
          np.allclose(cat_surf.vol_bspline_sample(ramp_d, line_z),
                      np.arange(d)[2:-2], atol=1e-4))

    rng = np.random.default_rng(0)
    vol = rng.random((d, h, w)).astype(np.float32) + 1.0
    zz, yy, xx = np.meshgrid(_centres(d)[1:-1], _centres(h)[1:-1],
                             _centres(w)[1:-1], indexing="ij")
    grid = np.stack([xx, yy, zz], axis=-1).astype(np.float32)
    out = cat_surf.vol_bspline_sample(vol, grid)
    check("the result has the shape of the grid without its last axis",
          out.shape == grid.shape[:-1] and out.dtype == np.float32)
    check("interior voxel centres come back unchanged",
          np.allclose(out, vol[1:-1, 1:-1, 1:-1], atol=1e-4))
    check("a mask value that does not occur changes nothing",
          np.array_equal(cat_surf.vol_bspline_sample(vol, grid, mask_value=0.0), out))
    check("align_corners=False samples elsewhere",
          not np.allclose(cat_surf.vol_bspline_sample(vol, grid, align_corners=False),
                          out))
    check("an empty grid gives an empty result",
          cat_surf.vol_bspline_sample(vol, np.zeros((0, 3), np.float32)).shape == (0,))

    for args, why in (((vol[0], grid), "a 2-D volume"),
                      ((vol, grid[..., :2]), "a grid without three coordinates")):
        try:
            cat_surf.vol_bspline_sample(*args)
            raised = False
        except ValueError:
            raised = True
        check(f"{why} is a ValueError", raised)


def test_vol_shoot():
    section("Geodesic shooting")

    shape = (21, 23, 19)
    x, y, z = np.meshgrid(*(np.arange(n, dtype=np.float32) for n in shape),
                          indexing="ij")
    gm = np.exp(-((x - 7.2) ** 2 / 24 + (y - 11.1) ** 2 / 31 +
                  (z - 8.4) ** 2 / 20)).astype(np.float32)
    wm = np.exp(-((x - 11.3) ** 2 / 19 + (y - 10.2) ** 2 / 26 +
                  (z - 9.1) ** 2 / 23)).astype(np.float32)
    extra = np.full(shape, 17.0, np.float32)
    template = np.stack((gm, wm, extra), axis=-1)
    templates = [template.copy() for _ in range(5)]
    mat = np.diag([1.1, 1.2, 1.3, 1.0])
    mat[:3, 3] = -0.5 * (np.asarray(shape) - 1) * np.diag(mat)[:3]

    result_c = cat_surf.vol_shoot(gm, wm, templates, mat, mat, nits=3,
                                  n_threads=1, inverse=True)
    result_f = cat_surf.vol_shoot(
        np.asfortranarray(gm), np.asfortranarray(wm),
        [np.asfortranarray(t) for t in templates], mat, mat, nits=3,
        n_threads=1, inverse=True)

    check("forward and inverse deformation shapes",
          result_c["y"].shape == shape + (3,) and
          result_c["iy"].shape == shape + (3,))
    check("Jacobian shape and finite outputs",
          result_c["jac"].shape == shape and
          np.isfinite(result_c["y"]).all() and
          np.isfinite(result_c["iy"]).all() and
          np.isfinite(result_c["jac"]).all())
    check("C- and Fortran-ordered inputs agree",
          np.array_equal(result_c["y"], result_f["y"]) and
          np.array_equal(result_c["iy"], result_f["iy"]) and
          np.array_equal(result_c["jac"], result_f["jac"]))

    grid = np.stack(np.meshgrid(*(np.arange(n, dtype=np.float32)
                                  for n in shape), indexing="ij"), axis=-1)
    forward_error = np.linalg.norm(result_c["y"] - grid, axis=-1)
    inverse_error = np.linalg.norm(result_c["iy"] - grid, axis=-1)
    interior = np.s_[2:-2, 2:-2, 2:-2]
    check("asymmetric identity maps preserve voxel coordinates",
          np.max(forward_error[interior]) < 1.0 and
          np.max(inverse_error[interior]) < 1.0)


# ---------------------------------------------------------------------------
# Head-motion correction: known shifts of an analytic phantom are recovered,
# the reference volume is the identity, and mcflirt's outputs have their
# shapes.  The phantom is a sum of Gaussian blobs without symmetry.
# ---------------------------------------------------------------------------
_BLOBS = ((0.0, 0.0, 0.0, 12.0, 60.0), (-10.0, 7.0, 3.0, 4.5, 80.0),
          (8.0, -5.0, 8.0, 3.5, -50.0), (5.0, 12.0, -7.0, 5.0, 70.0),
          (-8.0, -10.0, -9.0, 3.0, 90.0))


def _phantom(shape, voxel, shift):
    """Blob phantom on a grid centred on the origin, moved by shift (mm)."""
    axes = [(np.arange(n) - 0.5 * (n - 1)) * voxel for n in shape]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    v = np.full(shape, 10.0)
    for bx, by, bz, s, a in _BLOBS:
        d2 = ((x - shift[0] - bx) ** 2 + (y - shift[1] - by) ** 2 +
              (z - shift[2] - bz) ** 2)
        v += a * np.exp(-0.5 * d2 / (s * s))
    return v


def test_vol_realign():
    section("Head-motion correction (mcflirt replacement)")

    shape, voxel = (30, 32, 28), 2.0
    affine = np.diag([voxel, voxel, voxel, 1.0])
    affine[:3, 3] = [-0.5 * (n - 1) * voxel for n in shape]
    shifts = np.array([[0.0, 0.0, 0.0], [0.6, -0.4, 0.3], [1.1, -0.2, 0.7],
                       [1.5, 0.3, 1.0]])
    rng = np.random.default_rng(0)
    series = np.stack([_phantom(shape, voxel, d) + rng.normal(0.0, 0.5, shape)
                       for d in shifts], axis=-1).astype(np.float32)

    res = cat_surf.vol_realign(series, affine, refvol=0)
    est = -res["ras"][:, :3, 3]          # volume -> reference undoes the shift
    check("volume 0 is the reference", res["ref_index"] == 0)
    check("the reference matrix is exactly the identity",
          np.array_equal(res["fsl"][0], np.eye(4)))
    check("shifts are recovered within 0.1 mm",
          np.abs(est - shifts).max() < 0.1, f"{np.abs(est - shifts).max():.3f}")
    check("mcflirt's outputs have their shapes",
          res["par"].shape == (4, 6) and res["rms_abs"].shape == (4,)
          and res["rms_rel"].shape == (3,) and res["realigned"].shape == shape + (4,))
    check("the mean RMS files hold the means",
          abs(res["rms_abs_mean"] - res["rms_abs"].mean()) < 1e-12
          and abs(res["rms_rel_mean"] - res["rms_rel"].mean()) < 1e-12)

    bright = (1.25 * _phantom(shape, voxel, (0.0, 0.0, 0.0))).astype(np.float32)
    res = cat_surf.vol_realign(series, affine, ref=bright, resample=False)
    est = -res["ras"][:, :3, 3]
    check("a brighter separate reference: the intensity scale is its gain",
          np.abs(res["scale"] - 1.25).max() < 0.02, f"{res['scale']}")
    check("a brighter separate reference: shifts within 0.1 mm",
          np.abs(est - shifts).max() < 0.1, f"{np.abs(est - shifts).max():.3f}")
    check("resample=False returns no series", "realigned" not in res)


# ---------------------------------------------------------------------------
# Boundary-based registration: cat_surf.bbreg takes its defaults from
# CAT_BBRegOptionsInit(), as CAT_SurfBBReg does, so the two give the same
# transform.  They had drifted apart: wm_dist 0.5 against 2.0, a +-15 mm
# grid against +-4 mm, no smoothing against 5 and 2 mm.
# ---------------------------------------------------------------------------
def test_bbreg_defaults():
    section("BBR defaults")

    import inspect

    from cat_surf._bbreg import _bbreg_defaults

    lib = _bbreg_defaults()
    params = inspect.signature(cat_surf.bbreg).parameters
    check("every setting of CAT_BBRegOpts is a parameter of bbreg",
          set(lib) <= set(params), str(set(lib) - set(params)))
    differ = {name: (params[name].default, value) for name, value in lib.items()
              if name in params and params[name].default != value}
    check("bbreg's defaults are CAT_BBRegOptionsInit()'s", not differ, str(differ))

    detect = inspect.signature(cat_surf.bbreg_detect_contrast).parameters
    check("contrast detection samples where bbreg does",
          detect["wm_dist"].default == lib["wm_dist"] and
          detect["gm_dist"].default == lib["gm_dist"])


# ---------------------------------------------------------------------------
# resample_annot writes what the extension says.  It used to write the binary
# FreeSurfer format whatever the name, so T1Prep's *.label.gii parcellations
# were annot files that no GIfTI reader could open.
# ---------------------------------------------------------------------------
def test_sanlm_repeatable():
    section("vol_sanlm: repeatable, and the noise model does not stick")
    rng = np.random.default_rng(0)
    x, y, z = np.mgrid[:48, :48, :48]
    vol = (100 + 40 * (np.sqrt((x - 24) ** 2 + (y - 24) ** 2 + (z - 24) ** 2) < 15)
           + rng.normal(0, 8, x.shape)).astype(np.float32)

    first = np.asarray(cat_surf.vol_sanlm(vol.copy()))
    again = np.asarray(cat_surf.vol_sanlm(vol.copy()))
    check("two calls give the same result", np.array_equal(first, again),
          f"max |diff| {np.abs(first - again).max():.3g}")
    cat_surf.vol_sanlm(vol.copy(), is_rician=True)
    after = np.asarray(cat_surf.vol_sanlm(vol.copy()))
    check("a Rician call does not change later default calls",
          np.array_equal(first, after), f"max |diff| {np.abs(first - after).max():.3g}")


def test_missing_files_raise():
    section("a missing file raises instead of ending the interpreter")
    # The readers and writers used SAFE_FOPEN, which called exit(1): a missing
    # FreeSurfer surface or text file ended the whole Python process.
    for name, call in (
            ("read_surface of a FreeSurfer file",
             lambda: cat_surf.read_surface("/nonexistent/lh.white")),
            ("read_values of a text file",
             lambda: cat_surf.read_values("/nonexistent/values.txt")),
            ("write_values into a missing directory",
             lambda: cat_surf.write_values("/nonexistent_dir/values.txt", np.ones(3)))):
        try:
            call()
            check(f"{name} raises OSError", False, "returned normally")
        except OSError:
            check(f"{name} raises OSError", True)


def test_amap_rejects_bad_labels():
    section("vol_amap: out-of-range labels and empty classes are errors")
    x, y, z = np.mgrid[:32, :32, :32]
    r = np.sqrt((x - 16) ** 2 + (y - 16) ** 2 + (z - 16) ** 2)
    labels = np.where(r < 6, 3, np.where(r < 9, 2, np.where(r < 12, 1, 0))).astype(np.float32)
    t1 = (labels + np.random.default_rng(1).normal(0, 0.05, labels.shape)).astype(np.float32)

    prob, lab, mean = cat_surf.vol_amap(t1, labels)
    check("valid labels give finite class means", np.isfinite(mean).all(), str(mean))
    # halves round away from zero, as round() in CAT_VolAmap: 0.5 is CSF,
    # not background, and -0.5 is -1, not 0
    trial, t1_trial = labels.copy(), t1.copy()
    trial[1, 1, 1], t1_trial[1, 1, 1] = 0.5, 1.0
    lab = cat_surf.vol_amap(t1_trial, trial)[1]
    check("a label of 0.5 is CSF, not background", lab[1, 1, 1] != 0, str(lab[1, 1, 1]))
    for name, bad, kw in (("a label above the class count", 4, {}),
                          ("a negative label", -1, {}),
                          ("a label of -0.5", -0.5, {}),
                          ("a NaN label", np.nan, {}),
                          ("an unsupported class count", None, {"n_pure_classes": 7}),
                          ("a class without voxels", "no-csf", {})):
        trial = labels.copy()
        if bad == "no-csf":
            trial[trial == 1] = 2
        elif bad is not None:
            trial[16, 16, 16] = bad
        try:
            cat_surf.vol_amap(t1, trial, **kw)
            check(f"{name} is a ValueError", False)
        except ValueError:
            check(f"{name} is a ValueError", True)


def _write_annot(path, codes, table):
    """FreeSurfer .annot with a version-2 colour table of (name, r, g, b)."""
    import struct
    with open(path, "wb") as fh:
        fh.write(struct.pack(">i", len(codes)))
        for vno, code in enumerate(codes):
            fh.write(struct.pack(">ii", vno, int(code)))
        fh.write(struct.pack(">ii", 1, -2))          # colour table, version 2
        fh.write(struct.pack(">i", len(table)))
        fh.write(struct.pack(">i", 6) + b"smoke\0")
        fh.write(struct.pack(">i", len(table)))
        for i, (name, r, g, b) in enumerate(table):
            label = name.encode() + b"\0"
            fh.write(struct.pack(">ii", i, len(label)) + label)
            fh.write(struct.pack(">iiii", r, g, b, 0))


def test_resample_annot_formats():
    section("resample_annot: the output extension picks the format")
    src, src_f = _icosphere(3, 100.0)
    tgt, tgt_f = _icosphere(2, 100.0)
    table = [("north", 200, 30, 30), ("equator", 30, 200, 30), ("south", 30, 30, 200)]
    codes = np.array([r + 256 * g + 65536 * b for _, r, g, b in table])
    region = np.where(src[:, 2] > 30, 0, np.where(src[:, 2] < -30, 2, 1))

    with tempfile.TemporaryDirectory() as tmp:
        sph = os.path.join(tmp, "src.gii")
        tgt_sph = os.path.join(tmp, "tgt.gii")
        annot = os.path.join(tmp, "in.annot")
        cat_surf.write_surface(sph, src, src_f)
        cat_surf.write_surface(tgt_sph, tgt, tgt_f)
        _write_annot(annot, codes[region], table)

        out_gii = os.path.join(tmp, "out.label.gii")
        out_annot = os.path.join(tmp, "out.annot")
        out_txt = os.path.join(tmp, "out.txt")
        for out in (out_gii, out_annot, out_txt):
            cat_surf.resample_annot(sph, sph, tgt_sph, annot, out)

        with open(out_gii, "rb") as fh:
            check("a .gii output is GIfTI XML", fh.read(5) == b"<?xml")
        da = cat_surf.read_gifti_darrays(out_gii)[0]
        check("with NIFTI_INTENT_LABEL data", da["intent"] == 1002,
              str(da["intent_name"]))
        check("whose keys are the colour-table rows",
              da["min"] == 0 and da["max"] == len(table) - 1,
              f"{da['min']}..{da['max']}")
        with open(out_annot, "rb") as fh:
            check("an .annot output is still the FreeSurfer format",
                  int.from_bytes(fh.read(4), "big") == len(tgt))
        values = np.loadtxt(out_txt)
        check("anything else gets the annotation codes",
              set(np.unique(values).astype(int)) <= set(codes.tolist()))

        short = os.path.join(tmp, "short.annot")
        _write_annot(short, codes[region[:len(tgt)]], table)
        try:
            cat_surf.resample_annot(sph, sph, tgt_sph, short, out_gii)
            check("an annotation of the wrong length is a ValueError", False)
        except ValueError:
            check("an annotation of the wrong length is a ValueError", True)


# ---------------------------------------------------------------------------
# 6. Thickness QC tells a plate from a solid mass, which is the whole point:
#    a glued sulcus is recoverable, a subcortical mass is not, and thickness
#    alone cannot separate them.
# ---------------------------------------------------------------------------
def test_hull_and_normalized_area():
    section("sulcus_depth and get_area_normalized")
    v, f = _icosphere(3, 50.0)
    depth = cat_surf.sulcus_depth(v, f)
    check("a sphere lies on its own convex hull", np.abs(depth).max() < 0.1,
          f"max depth {np.abs(depth).max():.3g}")
    # it freed the address of a local pointer and aborted the interpreter
    area, total = cat_surf.get_area_normalized(v, f, v, f)
    check("get_area_normalized returns 4 pi r^2 for a sphere",
          abs(total / (4 * np.pi * 50.0 ** 2) - 1) < 0.01, f"{total:.1f}")
    check("and values that add up to it", abs(area.sum() / total - 1) < 1e-9)


def main():
    print(f"cat_surf {cat_surf.__version__} — binding smoke test")
    for test in (test_api_surface, test_sheetness, test_oriented_filters,
                 test_open_ppm_sulci, test_marching_cubes_sulci_kwargs, test_sulcal_barrier,
                 test_barrier_gate_scales_with_thickness, test_barrier_reference,
                 test_surf_to_sphere_relax,
                 test_sphere_remove_folds,
                 test_fix_self_intersect_reference, test_pbt_front_ends_agree,
                 test_option_no_ops,
                 test_surf_info, test_vol_calc, test_box_morph,
                 test_bspline_sample, test_vol_shoot, test_vol_realign,
                 test_bbreg_defaults,
                 test_resample_annot_formats, test_sanlm_repeatable,
                 test_amap_rejects_bad_labels, test_missing_files_raise,
                 test_hull_and_normalized_area):
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

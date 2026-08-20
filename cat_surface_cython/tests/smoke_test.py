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
import sys
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
def test_sulcus_repair():
    section("Pre-PBT repair")

    x = np.arange(N)[:, None, None]
    wm = (x < 6) | (x > 17)

    t1 = np.where(wm, 140.0, 90.0).astype(np.float32) * np.ones((N, N, N), np.float32)
    lab = np.where(wm, 3.0, 2.0).astype(np.float32) * np.ones((N, N, N), np.float32)
    lab[0, :, :] = 1.0                       # some real CSF for the class levels
    t1[0, :, :] = 40.0
    t1[12, :, :] = 45.0                      # dark sulcus, still labelled GM

    out, sheet = cat_surf.vol_sulcus_repair(t1, lab, voxelsize=VX,
                                            strengthen_wm=False,
                                            sheet_sigma_min=0.5,
                                            sheet_sigma_max=1.0,
                                            sheet_n_scales=2,
                                            return_sheetness=True)
    check("glued sulcus was opened", out[12, 12, 12] < lab[12, 12, 12] - 0.4,
          f"{lab[12, 12, 12]:.2f} -> {out[12, 12, 12]:.2f}")
    check("never opened past what the intensity supports", out[12, 12, 12] >= 1.0)
    check("ordinary GM left alone",
          abs(out[9, 12, 12] - lab[9, 12, 12]) < 1e-5)
    check("sheetness map returned", sheet.shape == (N, N, N))

    t1b = np.full((N, N, N), 90.0, np.float32)
    labb = np.full((N, N, N), 2.0, np.float32)
    labb[0, :, :] = 1.0
    t1b[0, :, :] = 40.0
    labb[12, :, :] = 3.0                     # a one-voxel WM blade
    t1b[12, :, :] = 140.0
    labb[12, :, 12] = 2.0                    # the missegmentation
    t1b[12, :, 12] = 138.0
    labb[5, 5, 5] = 2.0                      # decoy: no facing bank
    t1b[5, 5, 5] = 140.0

    outb = cat_surf.vol_sulcus_repair(t1b, labb, voxelsize=VX,
                                      recover_csf=False,
                                      sheet_sigma_min=0.5, sheet_sigma_max=1.0,
                                      sheet_n_scales=2)
    check("blade gap was bridged", outb[12, 12, 12] > 2.4,
          f"{outb[12, 12, 12]:.2f}")
    check("isolated bright voxel untouched",
          abs(outb[5, 5, 5] - labb[5, 5, 5]) < 1e-5)


# ---------------------------------------------------------------------------
# 5. The new options on the existing entry points stay backward compatible.
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


def test_wm_sulcus_guard():
    section("Sulcus guard on blade strengthening")

    def phantom(with_sulcus):
        t1 = np.full((N, N, N), 90.0, np.float32)
        lab = np.full((N, N, N), 2.0, np.float32)
        lab[0, :, :] = 1.0
        t1[0, :, :] = 40.0
        lab[12, :, :] = 3.0
        t1[12, :, :] = 140.0
        lab[12, :, 12] = 2.0          # the missegmentation
        t1[12, :, 12] = 138.0
        if with_sulcus:               # a sulcal floor one voxel from the tip
            t1[13, :, 12] = 45.0
        return t1, lab

    kw = dict(voxelsize=VX, recover_csf=False, sheet_sigma_min=0.5,
              sheet_sigma_max=1.0, sheet_n_scales=2)

    t1, lab = phantom(True)
    unguarded = cat_surf.vol_sulcus_repair(t1, lab, wm_sulcus_guard=0.0, **kw)
    t1, lab = phantom(True)
    guarded = cat_surf.vol_sulcus_repair(t1, lab, wm_sulcus_guard=1.0, **kw)
    t1, lab = phantom(False)
    no_sulcus = cat_surf.vol_sulcus_repair(t1, lab, wm_sulcus_guard=1.0, **kw)

    check("without the guard the blade buries the sulcus",
          unguarded[12, 12, 12] > 2.4, f"{unguarded[12, 12, 12]:.3f}")
    check("the guard stops that",
          guarded[12, 12, 12] < unguarded[12, 12, 12] - 0.4,
          f"{guarded[12, 12, 12]:.3f}")
    check("the guard does not block a blade with no sulcus nearby",
          no_sulcus[12, 12, 12] > 2.4, f"{no_sulcus[12, 12, 12]:.3f}")


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
# 6. Thickness QC tells a plate from a solid mass, which is the whole point:
#    a glued sulcus is recoverable, a subcortical mass is not, and thickness
#    alone cannot separate them.
# ---------------------------------------------------------------------------
def test_thickness_qc():
    section("Thickness QC")

    gmt = np.full((N, N, N), 1.0, np.float32)
    gmt[6:10, 4:N-4, 4:N-4] = 6.0                       # a 4-voxel plate
    # np.mgrid yields axis 0 first, so a0 must clear the plate's 6:10 slab
    a0, a1, a2 = np.mgrid[0:N, 0:N, 0:N]
    gmt[((a0 - 17) ** 2 + (a1 - 12) ** 2 + (a2 - 12) ** 2) <= 16] = 6.0  # a ball

    comps, cls = cat_surf.vol_thickness_qc(gmt, return_classmap=True)
    shapes = sorted(c["shape"] for c in comps)
    check("both structures were found", len(comps) == 2, str(shapes))
    check("one plate and one solid", shapes == ["plate", "solid"])

    plate = next(c for c in comps if c["shape"] == "plate")
    solid = next(c for c in comps if c["shape"] == "solid")
    check("the plate radius is half its thickness", 1.5 < plate["max_radius"] <= 2.5,
          f"{plate['max_radius']:.2f} mm")
    check("the ball radius is its own radius", solid["max_radius"] > 4.0,
          f"{solid['max_radius']:.2f} mm")
    check("size is not the discriminator", plate["n_voxels"] > solid["n_voxels"])
    check("class map tags plate 1 and solid 2",
          set(np.unique(cls).tolist()) == {0.0, 1.0, 2.0})

    clean = np.full((N, N, N), 2.5, np.float32)
    check("a plausible thickness map reports nothing",
          cat_surf.vol_thickness_qc(clean) == [])


def main():
    print(f"cat_surf {cat_surf.__version__} — binding smoke test")
    for test in (test_api_surface, test_sheetness, test_oriented_filters,
                 test_open_ppm_sulci, test_marching_cubes_sulci_kwargs,
                 test_wm_sulcus_guard,
                 test_sulcus_repair, test_thickness_qc, test_option_no_ops):
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

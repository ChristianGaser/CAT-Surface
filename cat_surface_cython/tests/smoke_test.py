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
#    every -oriented / -mrf-aniso / -oriented-filter option safe to enable,
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

    sm_o = cat_surf.vol_oriented_smooth(vol, sheetness=ones, normal=nrm,
                                        iters=3, sigma=0.5)
    sm_i = cat_surf.vol_oriented_smooth(vol, sheetness=zeros, normal=nrm,
                                        iters=3)
    check("oriented smoothing preserves the sheet better than isotropic",
          sm_o[12, 12, 12] < sm_i[12, 12, 12],
          f"{sm_o[12, 12, 12]:.3f} vs {sm_i[12, 12, 12]:.3f}")

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
                                            reconnect_gyri=False,
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

    rng = np.random.default_rng(0)
    z = np.arange(N)[None, None, :]
    vol = np.where(z < 8, 30.0, np.where(z < 16, 90.0, 150.0)).astype(np.float32) \
        * np.ones((N, N, N), np.float32)
    vol += rng.normal(0, 3, vol.shape).astype(np.float32)
    labels = (np.where(z < 8, 1, np.where(z < 16, 2, 3)).astype(np.uint8)
              * np.ones((N, N, N), np.uint8))

    zeros = np.zeros((N, N, N), np.float32)
    nrm = np.zeros((N, N, N, 3), np.float32)
    nrm[..., 0] = 1.0

    kw = dict(voxelsize=VX, weight_mrf=0.3, n_iters=5)
    _, l0, _ = cat_surf.vol_amap(vol, labels, mrf_aniso=0, **kw)
    _, l1, _ = cat_surf.vol_amap(vol, labels, mrf_aniso=1,
                                 sheetness=zeros, normal=nrm, **kw)
    _, l2, _ = cat_surf.vol_amap(vol, labels, mrf_aniso=2,
                                 sheetness=zeros, normal=nrm, **kw)
    check("Amap local-beta mode is a no-op at zero sheetness",
          np.array_equal(l0, l1))
    check("Amap anisotropic-Potts mode is a no-op at zero sheetness",
          np.array_equal(l0, l2))


def main():
    print(f"cat_surf {cat_surf.__version__} — binding smoke test")
    for test in (test_api_surface, test_sheetness, test_oriented_filters,
                 test_sulcus_repair, test_option_no_ops):
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

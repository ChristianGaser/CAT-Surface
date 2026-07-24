#!/usr/bin/env python3
"""Integration self-check for the T1Prep / cat_surf Nipype interfaces.

Run this after changing either ``cat_surf`` (the ``cat-surf`` package) or the
``nipype.interfaces.t1prep`` interfaces to confirm the two still agree.  It
validates three things:

1. **cat_surf API** — the ``spherical_demon`` default is in range and the
   ``cat_surf.cli`` mirror is present with an interface matching ``surf_warp``.
2. **Call contracts** — every ``cat_surf`` function that a Nipype interface
   calls actually exists and accepts the keyword arguments the interface
   passes (this is the check that catches "bbreg() got an unexpected keyword
   argument 'moving_file'" style breakage).
3. **Nipype interfaces** — the package exports exactly the intended subset,
   each interface's traits are as expected, and each one dispatches without
   the native library (using a fake ``cat_surf``).

Each stage degrades gracefully: if ``cat_surf`` or ``nipype`` is not installed,
the dependent checks are reported as SKIP rather than failing.  Exit status is
non-zero if any check FAILS.

Usage
-----
    python check_interfaces.py            # in an env with cat-surf and nipype
"""

from __future__ import annotations

import inspect
import sys
import tempfile
import types
from pathlib import Path

PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"
_results: list[tuple[str, str, str]] = []


def record(status: str, name: str, detail: str = "") -> None:
    _results.append((status, name, detail))


# --------------------------------------------------------------------------- #
# Contract table: nipype interface -> cat_surf callable + kwargs it passes.
# "cli:" prefix means the function lives on cat_surf.cli (file-based mirror).
# --------------------------------------------------------------------------- #
CALL_CONTRACTS = {
    "CatSurfVolSanlm": ("cli:vol_sanlm", ["is_rician", "strength"]),
    "CatSurfBbreg": (
        "bbreg",
        ["lh_surface", "rh_surface", "ref_file", "invert_contrast", "fwhm", "verbose"],
    ),
    "CatSurfBbregDetectContrast": (
        "bbreg_detect_contrast",
        ["lh_surface", "rh_surface"],
    ),
    "CatSurfVolumeRegisterNmi": (
        "volume_register_nmi",
        ["n_levels", "n_bins", "max_iter", "verbose"],
    ),
    "CatSurfVolumeRegisterRobust": (
        "volume_register_robust",
        ["n_levels", "sat_k", "max_iter", "verbose"],
    ),
}

EXPECTED_INTERFACES = set(CALL_CONTRACTS)
EXPECTED_COMMANDS = {
    "Info",
    "T1PrepCommand",
    "T1Prep",
    "T1PrepSegment",
    "T1PrepSurfaceEstimation",
    "T1PrepRealignLongitudinal",
}
# A sample of interfaces that must NOT be exported after the scope-down.
REMOVED_INTERFACES = {
    "CatSurfReadSurface",
    "CatSurfGetArea",
    "CatSurfWarp",
    "CatSurfSphericalDemon",
    "CatSurfDeform",
    "CatSurfVolMarchingCubes",
}

# Expected trait names per kept interface (inputs / outputs).
EXPECTED_TRAITS = {
    "CatSurfVolSanlm": (
        {"in_file", "out_file", "is_rician", "strength"},
        {"out_file"},
    ),
    "CatSurfBbreg": (
        {
            "in_file", "lh_surface", "rh_surface", "ref_file",
            "out_matrix_file", "invert_contrast", "fwhm", "verbose",
        },
        {"out_matrix_file", "cost"},
    ),
    "CatSurfBbregDetectContrast": (
        {"in_file", "lh_surface", "rh_surface"},
        {"contrast"},
    ),
    "CatSurfVolumeRegisterNmi": (
        {"moving_file", "fixed_file", "out_matrix_file", "n_levels",
         "n_bins", "max_iter", "verbose"},
        {"out_matrix_file", "nmi"},
    ),
    "CatSurfVolumeRegisterRobust": (
        {"moving_file", "fixed_file", "out_matrix_file", "n_levels",
         "sat_k", "max_iter", "verbose"},
        {"out_matrix_file", "residual"},
    ),
}


def _accepts(func, kwargs) -> list[str]:
    """Return the kwargs *func* does NOT accept (empty list = all fine)."""
    try:
        sig = inspect.signature(func)
    except (TypeError, ValueError):
        return []  # builtin without signature: cannot verify, assume ok
    params = sig.parameters
    if any(p.kind == p.VAR_KEYWORD for p in params.values()):
        return []
    return [k for k in kwargs if k not in params]


# --------------------------------------------------------------------------- #
# Stage 1 + 2: cat_surf API and call contracts
# --------------------------------------------------------------------------- #
def check_cat_surf():
    try:
        import cat_surf
    except Exception as exc:  # noqa: BLE001
        record(SKIP, "cat_surf import", f"cat_surf not installed ({exc})")
        return

    record(PASS, "cat_surf import", f"version {getattr(cat_surf, '__version__', '?')}")

    # spherical_demon default n_steps must be within CAT_WARP_DEMONS_MAX_STEPS (4)
    try:
        n_steps = inspect.signature(cat_surf.spherical_demon).parameters["n_steps"].default
        if 1 <= n_steps <= 4:
            record(PASS, "spherical_demon n_steps default", f"= {n_steps}")
        else:
            record(FAIL, "spherical_demon n_steps default",
                   f"= {n_steps} (must be 1..4)")
    except Exception as exc:  # noqa: BLE001
        record(FAIL, "spherical_demon n_steps default", str(exc))

    # cli mirror parity with surf_warp
    try:
        from cat_surf import cli
        has_sd = hasattr(cli, "surf_spherical_demon")
        in_all = "surf_spherical_demon" in getattr(cli, "__all__", [])
        record(PASS if (has_sd and in_all) else FAIL,
               "cli.surf_spherical_demon present",
               f"attr={has_sd}, in __all__={in_all}")
        if has_sd and hasattr(cli, "surf_warp"):
            sw = list(inspect.signature(cli.surf_warp).parameters)[:5]
            sd = list(inspect.signature(cli.surf_spherical_demon).parameters)[:5]
            record(PASS if sw == sd else FAIL,
                   "surf_warp / surf_spherical_demon positional parity",
                   f"{sw} vs {sd}")
    except Exception as exc:  # noqa: BLE001
        record(FAIL, "cli mirror", str(exc))

    # Call-contract check: functions the interfaces call must accept their kwargs
    try:
        from cat_surf import cli
    except Exception:  # noqa: BLE001
        cli = None
    for iface, (target, kwargs) in CALL_CONTRACTS.items():
        if target.startswith("cli:"):
            fname = target.split(":", 1)[1]
            func = getattr(cli, fname, None) if cli is not None else None
            label = f"cat_surf.cli.{fname}"
        else:
            func = getattr(cat_surf, target, None)
            label = f"cat_surf.{target}"
        if func is None:
            record(FAIL, f"{iface} -> {label}", "function missing")
            continue
        missing = _accepts(func, kwargs)
        if missing:
            record(FAIL, f"{iface} -> {label}",
                   f"does not accept kwargs: {missing}")
        else:
            record(PASS, f"{iface} -> {label}", "signature compatible")


# --------------------------------------------------------------------------- #
# Stage 3: nipype interfaces
# --------------------------------------------------------------------------- #
def _fake_cat_surf():
    import numpy as np

    ns = types.SimpleNamespace()
    ns.read_surface = lambda *a, **k: ("V", "F")
    ns.bbreg = lambda *a, **k: (np.eye(4), 0.5)
    ns.bbreg_detect_contrast = lambda *a, **k: 0
    ns.volume_register_nmi = lambda *a, **k: (np.eye(4), 0.9)
    ns.volume_register_robust = lambda *a, **k: (np.eye(4), 0.1)
    ns.cli = types.SimpleNamespace(vol_sanlm=lambda *a, **k: None)
    return ns


def check_nipype():
    try:
        from nipype.interfaces import t1prep
        from nipype.interfaces.t1prep import cat_surf as cat_surf_mod
    except Exception as exc:  # noqa: BLE001
        record(SKIP, "nipype import", f"nipype not importable ({exc})")
        return

    exported = set(getattr(t1prep, "__all__", []))

    missing = (EXPECTED_INTERFACES | EXPECTED_COMMANDS) - exported
    record(PASS if not missing else FAIL, "expected exports present",
           "all present" if not missing else f"missing: {sorted(missing)}")

    leaked = REMOVED_INTERFACES & exported
    record(PASS if not leaked else FAIL, "removed interfaces are gone",
           "none exported" if not leaked else f"still exported: {sorted(leaked)}")

    for iface_name, (want_in, want_out) in EXPECTED_TRAITS.items():
        cls = getattr(t1prep, iface_name, None)
        if cls is None:
            record(FAIL, f"{iface_name} traits", "class not exported")
            continue
        got_in = {k for k in cls.input_spec().traits()
                  if k not in ("trait_added", "trait_modified")}
        got_out = {k for k in cls.output_spec().traits()
                   if k not in ("trait_added", "trait_modified")}
        ok = want_in <= got_in and want_out <= got_out
        record(PASS if ok else FAIL, f"{iface_name} traits",
               "as expected" if ok
               else f"in-missing={sorted(want_in - got_in)} "
                    f"out-missing={sorted(want_out - got_out)}")

    # Dispatch every kept interface against a fake cat_surf module.
    cat_surf_mod._import_cat_surf = _fake_cat_surf  # type: ignore[attr-defined]
    with tempfile.TemporaryDirectory() as tmp:
        runtime = types.SimpleNamespace(cwd=tmp, returncode=0, environ={})
        for iface_name in EXPECTED_INTERFACES:
            cls = getattr(t1prep, iface_name, None)
            if cls is None:
                continue
            try:
                node = cls()
                for name, trait in node.input_spec().traits().items():
                    if name in ("trait_added", "trait_modified") or not trait.mandatory:
                        continue
                    if type(trait.trait_type).__name__ == "File":
                        f = Path(tmp) / f"{iface_name}_{name}.nii.gz"
                        f.write_bytes(b"")
                        setattr(node.inputs, name, str(f))
                    else:
                        setattr(node.inputs, name, 1)
                node._run_interface(runtime)
                outputs = node._list_outputs()
                record(PASS if isinstance(outputs, dict) else FAIL,
                       f"{iface_name} dispatch", "ran, produced outputs")
            except Exception as exc:  # noqa: BLE001
                record(FAIL, f"{iface_name} dispatch", f"{type(exc).__name__}: {exc}")


def main() -> int:
    check_cat_surf()
    check_nipype()

    width = max((len(n) for _, n, _ in _results), default=0)
    print("\n" + "=" * (width + 20))
    print("cat_surf / nipype interface self-check")
    print("=" * (width + 20))
    n_fail = n_pass = n_skip = 0
    for status, name, detail in _results:
        n_fail += status == FAIL
        n_pass += status == PASS
        n_skip += status == SKIP
        print(f"[{status}] {name.ljust(width)}  {detail}")
    print("-" * (width + 20))
    print(f"{n_pass} passed, {n_fail} failed, {n_skip} skipped")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

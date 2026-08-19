"""
Build script for cat_surf — Cython bindings to libCAT.

Usage (local development):
    python setup.py build_ext --inplace

Usage (CI / pip wheel):
    export CAT_BUILD_DIR=/path/to/build    # autotools build tree
    pip wheel .

Cython is required: the generated .c sources are build artifacts and are
not tracked in the repository.  pyproject.toml declares Cython in
build-system.requires, so any PEP 517 build installs it automatically; a
direct "python setup.py build_ext" needs it in the ambient environment.

CI wheel builds via cibuildwheel use a vendored layout: headers and
static libraries are staged into ``_vendor/`` inside this directory
before running cibuildwheel.  See python-wheels.yml for details.

Environment variables:
    CAT_SURFACE_ROOT  Source tree root   (default: parent of this file)
    CAT_BUILD_DIR     Autotools build tree containing .libs/libCAT.a
"""
import os
import re
import subprocess
import sys
import sysconfig

import numpy as np
from setuptools import setup, Extension

_HERE = os.path.abspath(os.path.dirname(__file__))

# Linking a libCAT.a that no longer matches the sources — or that does not
# cover every architecture the extension is built for — still produces a
# wheel that installs cleanly: unresolved symbols in a Python extension are
# left for dynamic lookup.  The failure only surfaces much later, at import,
# as "symbol not found in flat namespace '_CAT_...'".  The checks below turn
# those silent mismatches into visible warnings at build time.
_STRICT = os.environ.get("CAT_STRICT_BUILD", "") not in ("", "0")


def _warn(msg):
    """Report a build-time mismatch; abort instead if CAT_STRICT_BUILD is set."""
    if _STRICT:
        sys.exit(f"ERROR: {msg}")
    print(f"WARNING: {msg}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Paths — CI puts everything into _vendor/; local dev uses the source tree
# ---------------------------------------------------------------------------
_VENDOR = os.path.join(_HERE, "_vendor")

if os.path.isdir(_VENDOR):
    # ---- CI / vendored mode ----
    CAT_ROOT = _VENDOR
    BUILD_DIR = os.path.join(_VENDOR, "build")
else:
    # ---- Local development mode ----
    CAT_ROOT = os.environ.get(
        "CAT_SURFACE_ROOT",
        os.path.abspath(os.path.join(_HERE, os.pardir)),
    )
    # "." covers an in-source autotools build, which lives at the source root
    # rather than in a build-*/ subdirectory.  Omitting it used to make an
    # abandoned out-of-tree tree shadow the build the developer is actually
    # working in.
    _platform_dirs = [
        "build-native-arm64",
        "build-native",
        "build-x86_64-pc-linux",
        "build-x86_64-w64-mingw32",
        "build",
        ".",
    ]
    BUILD_DIR = os.environ.get("CAT_BUILD_DIR", "")
    if not BUILD_DIR:
        # Rank the candidates by archive mtime and take the newest.  Taking
        # the *first* match instead silently preferred whichever tree came
        # earliest in the list, however stale it was.
        _found = []
        for d in _platform_dirs:
            candidate = os.path.normpath(os.path.join(CAT_ROOT, d))
            for _rel in (os.path.join(".libs", "libCAT.a"), "libCAT.a", "libCAT.la"):
                _archive = os.path.join(candidate, _rel)
                if os.path.isfile(_archive):
                    _found.append((os.path.getmtime(_archive), candidate))
                    break
        if _found:
            _found.sort(reverse=True)
            BUILD_DIR = _found[0][1]
            if len(_found) > 1:
                print(
                    f"Using newest libCAT build tree: {BUILD_DIR}\n"
                    "  (ignoring older: "
                    + ", ".join(c for _, c in _found[1:])
                    + ")",
                    file=sys.stderr,
                )
    if not BUILD_DIR:
        sys.exit(
            "ERROR: Could not locate a build directory with libCAT.  "
            "Set CAT_BUILD_DIR or build CAT-Surface first."
        )

LIBS_DIR = os.path.join(BUILD_DIR, ".libs")
if not os.path.isdir(LIBS_DIR):
    LIBS_DIR = BUILD_DIR

# ---------------------------------------------------------------------------
# Static archives
# ---------------------------------------------------------------------------
LIBCAT_A = os.path.join(LIBS_DIR, "libCAT.a")
if not os.path.isfile(LIBCAT_A):
    sys.exit(f"ERROR: Cannot find {LIBCAT_A}")

LIBFFTW3_A = os.path.join(BUILD_DIR, "3rdparty", "fftw-build", ".libs", "libfftw3.a")
if not os.path.isfile(LIBFFTW3_A):
    sys.exit(f"ERROR: Cannot find {LIBFFTW3_A}")


# ---------------------------------------------------------------------------
# Sanity checks on the archive we are about to link
# ---------------------------------------------------------------------------
def _newest_source_mtime(root):
    """Return the mtime of the most recently edited libCAT C source/header."""
    newest = 0.0
    for sub in ("Lib", "Include"):
        for dirpath, _dirnames, filenames in os.walk(os.path.join(root, sub)):
            for name in filenames:
                if name.endswith((".c", ".h")):
                    newest = max(newest, os.path.getmtime(os.path.join(dirpath, name)))
    return newest


_src_mtime = _newest_source_mtime(CAT_ROOT)
if _src_mtime and _src_mtime > os.path.getmtime(LIBCAT_A):
    _warn(
        f"{LIBCAT_A}\n  is older than the newest C source in {CAT_ROOT}.\n"
        "  The bindings will link against a stale library and any symbol added\n"
        "  since that build will fail at import, not at build time.\n"
        "  Run 'make' in the build tree first, or point CAT_BUILD_DIR at the\n"
        "  tree you actually build in."
    )


def _archive_archs(path):
    """Architectures present in a macOS archive, or an empty set elsewhere."""
    if sys.platform != "darwin":
        return set()
    try:
        out = subprocess.run(
            ["lipo", "-archs", path], capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return set()
    return set(out.split())


def _target_archs():
    """Architectures this extension is being compiled for on macOS."""
    if sys.platform != "darwin":
        return set()
    flags = os.environ.get("ARCHFLAGS") or sysconfig.get_config_var("CFLAGS") or ""
    return set(re.findall(r"-arch\s+(\S+)", flags))


_have = _archive_archs(LIBCAT_A)
_want = _target_archs()
_missing = _want - _have if (_have and _want) else set()
if _missing:
    _warn(
        f"{LIBCAT_A}\n  provides {'/'.join(sorted(_have))} but the extension is "
        f"being built for {'/'.join(sorted(_want))}.\n"
        f"  The {'/'.join(sorted(_missing))} slice(s) would link with every libCAT\n"
        "  symbol unresolved and fail at import.  Build a universal libCAT, or\n"
        f"  restrict this build:  ARCHFLAGS=\"{' '.join('-arch ' + a for a in sorted(_have))}\" "
        "pip install ."
    )

# ---------------------------------------------------------------------------
# Include directories
# ---------------------------------------------------------------------------
include_dirs = [
    np.get_include(),
    os.path.join(CAT_ROOT, "Include"),
    os.path.join(CAT_ROOT, "3rdparty", "bicpl-surface", "Include"),
    os.path.join(CAT_ROOT, "3rdparty", "volume_io", "Include"),
    os.path.join(CAT_ROOT, "3rdparty", "nifti"),
    os.path.join(CAT_ROOT, "3rdparty", "gifticlib"),
    os.path.join(CAT_ROOT, "3rdparty", "nii2mesh"),
    os.path.join(CAT_ROOT, "3rdparty", "zlib"),
    os.path.join(CAT_ROOT, "3rdparty", "expat"),
    os.path.join(CAT_ROOT, "3rdparty", "dartel"),
    BUILD_DIR,                       # for config.h
]

# ---------------------------------------------------------------------------
# Link flags (platform-specific)
# ---------------------------------------------------------------------------
extra_link_args = [LIBCAT_A, LIBFFTW3_A, "-lm", "-lz"]

if sys.platform == "darwin":
    # macOS: link the (system) expat the bundled headers expect.
    #
    # No OpenMP anywhere: libCAT parallelises with pthreads (resolved via
    # libSystem on macOS), so cat_surf's .so files pull in no libomp.  This
    # is what keeps cat_surf from becoming a second OpenMP runtime in a
    # process that already hosts one (e.g. PyTorch's libomp.dylib), which
    # used to cause "OMP: Error #15" / thread-pool TLS corruption.
    extra_link_args.append("-lexpat")
elif sys.platform == "linux":
    # Linux: pthread for libCAT's threads.  libCAT is pure C, so no C++
    # runtime has to be linked in.
    extra_link_args += ["-lpthread"]

# ---------------------------------------------------------------------------
# Cython is required
# ---------------------------------------------------------------------------
# The generated .c files are build artifacts and are not tracked (see
# .gitignore); every module is cythonized from its .pyx here.  Cython is
# declared in pyproject.toml's build-system.requires, so PEP 517 builds --
# pip install, pip wheel, cibuildwheel -- always have it in the isolated build
# environment.  Only a direct "python setup.py build_ext" in an environment
# without Cython can reach the error below.
try:
    from Cython.Build import cythonize
except ImportError:  # pragma: no cover
    sys.exit(
        "ERROR: Cython is required to build cat-surf.\n"
        "       Install it with 'pip install cython', or build through pip\n"
        "       ('pip install .'), which installs it automatically."
    )

ext_suffix = ".pyx"

# ---------------------------------------------------------------------------
# Extension modules
# ---------------------------------------------------------------------------
common_kwargs = dict(
    include_dirs=include_dirs,
    extra_link_args=extra_link_args,
    language="c",
    define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
                   ("BICAPI", "")],
)

extensions = [
    Extension(
        "cat_surf._convert",
        [os.path.join("cat_surf", "_convert" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._volume",
        [os.path.join("cat_surf", "_volume" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._io",
        [os.path.join("cat_surf", "_io" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._surf",
        [os.path.join("cat_surf", "_surf" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._vol",
        [os.path.join("cat_surf", "_vol" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._bbreg",
        [os.path.join("cat_surf", "_bbreg" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._vol2surf",
        [os.path.join("cat_surf", "_vol2surf" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._surf_warp",
        [os.path.join("cat_surf", "_surf_warp" + ext_suffix)],
        **common_kwargs,
    ),
    Extension(
        "cat_surf._spherical_demon",
        [os.path.join("cat_surf", "_spherical_demon" + ext_suffix)],
        **common_kwargs,
    ),
]

extensions = cythonize(
    extensions,
    compiler_directives={
        "language_level": "3",
        "boundscheck": False,
        "wraparound": False,
        "cdivision": True,
    },
)

# ---------------------------------------------------------------------------
# Package setup (metadata lives in pyproject.toml)
# ---------------------------------------------------------------------------
setup(
    ext_modules=extensions,
    packages=["cat_surf", "cat_surf.cli"],
)

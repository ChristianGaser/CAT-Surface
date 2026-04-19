"""
Build script for cat_surf — Cython bindings to libCAT.

Usage (local development):
    python setup.py build_ext --inplace

Usage (CI / pip wheel):
    export CAT_BUILD_DIR=/path/to/build    # autotools build tree
    pip wheel .

CI wheel builds via cibuildwheel use a vendored layout: headers and
static libraries are staged into ``_vendor/`` inside this directory
before running cibuildwheel.  See python-wheels.yml for details.

Environment variables:
    CAT_SURFACE_ROOT  Source tree root   (default: parent of this file)
    CAT_BUILD_DIR     Autotools build tree containing .libs/libCAT.a
"""
import os
import platform
import sys

import numpy as np
from setuptools import setup, Extension

_HERE = os.path.abspath(os.path.dirname(__file__))

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
    _platform_dirs = [
        "build-native-arm64",
        "build-native",
        "build-x86_64-pc-linux",
        "build-x86_64-w64-mingw32",
        "build",
    ]
    BUILD_DIR = os.environ.get("CAT_BUILD_DIR", "")
    if not BUILD_DIR:
        for d in _platform_dirs:
            candidate = os.path.join(CAT_ROOT, d)
            if os.path.isfile(os.path.join(candidate, ".libs", "libCAT.a")) or \
               os.path.isfile(os.path.join(candidate, "libCAT.la")):
                BUILD_DIR = candidate
                break
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

_machine = platform.machine().lower()

if sys.platform == "darwin":
    # macOS: system expat + bundled static libomp
    extra_link_args.append("-lexpat")
    _omp_a = os.path.join(
        CAT_ROOT, "3rdparty", "libomp", "lib",
        f"libomp-{_machine}.a",
    )
    if os.path.isfile(_omp_a):
        extra_link_args.append(_omp_a)
    else:
        extra_link_args.append("-lomp")
elif sys.platform == "linux":
    # Linux: system GOMP + pthread + stdc++/gcc for MeshFix C++ objects
    extra_link_args += ["-lgomp", "-lstdc++", "-lpthread"]
elif sys.platform == "win32":
    # Windows (MinGW): GOMP + stdc++
    extra_link_args += ["-lgomp", "-lstdc++"]

# ---------------------------------------------------------------------------
# Use Cython if available, else fall back to pre-generated .c files
# ---------------------------------------------------------------------------
try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

ext_suffix = ".pyx" if USE_CYTHON else ".c"

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
]

if USE_CYTHON:
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
    packages=["cat_surf"],
)

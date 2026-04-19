"""
Build script for cat_surf — Cython bindings to libCAT.

Usage:
    python setup.py build_ext --inplace
    pip install -e .
"""
import os
import sys
import numpy as np
from setuptools import setup, Extension

# ---------------------------------------------------------------------------
# Paths — adjust CAT_SURFACE_ROOT if building out-of-tree
# ---------------------------------------------------------------------------
CAT_ROOT = os.environ.get(
    "CAT_SURFACE_ROOT",
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)),
)

# Detect platform-specific build directory for libCAT.a
_platform_dirs = [
    "build-native-arm64",
    "build-native",
    "build-x86_64-pc-linux",
    "build-x86_64-w64-mingw32",
]
BUILD_DIR = os.environ.get("CAT_BUILD_DIR", "")
if not BUILD_DIR:
    for d in _platform_dirs:
        candidate = os.path.join(CAT_ROOT, d)
        # Look for the libtool archive or the static lib
        if os.path.isfile(os.path.join(candidate, "libCAT.la")) or \
           os.path.isfile(os.path.join(candidate, ".libs", "libCAT.a")):
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
    BUILD_DIR,                       # for config.h
]

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
# Link against the *static* archive directly so that all symbols from
# bundled 3rdparty code (bicpl, volume_io, gifticlib, …) are pulled in.
# Using -lCAT with -undefined dynamic_lookup would defer resolution and
# fail at import time.
LIBCAT_A = os.path.join(LIBS_DIR, "libCAT.a")
if not os.path.isfile(LIBCAT_A):
    sys.exit(f"ERROR: Cannot find {LIBCAT_A}")

library_dirs = []
libraries = []

# Extra link flags — we need math, zlib, and possibly fftw3
extra_link_args = [LIBCAT_A, "-lm", "-lz"]

# FFTW3: try pkg-config first, fall back to -lfftw3
try:
    import subprocess
    fftw_libs = subprocess.check_output(
        ["pkg-config", "--libs", "fftw3"], text=True,
    ).strip().split()
    extra_link_args.extend(fftw_libs)
except Exception:
    extra_link_args.append("-lfftw3")

# macOS needs expat from system or bundled
if sys.platform == "darwin":
    extra_link_args.append("-lexpat")

# OpenMP: libCAT uses OpenMP — link the bundled static libomp on macOS
import platform
_machine = platform.machine().lower()
if sys.platform == "darwin":
    _omp_a = os.path.join(
        CAT_ROOT, "3rdparty", "libomp", "lib",
        f"libomp-{_machine}.a",
    )
    if os.path.isfile(_omp_a):
        extra_link_args.append(_omp_a)
    else:
        extra_link_args.append("-lomp")

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
    library_dirs=library_dirs,
    libraries=libraries,
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
# Package metadata
# ---------------------------------------------------------------------------
setup(
    name="cat_surf",
    version="0.1.0",
    description="Cython bindings to the CAT-Surface library (libCAT)",
    author="Christian Gaser",
    author_email="christian.gaser@uni-jena.de",
    license="GPL-2.0",
    packages=["cat_surf"],
    ext_modules=extensions,
    python_requires=">=3.8",
    install_requires=["numpy>=1.20"],
    extras_require={
        "io": ["nibabel>=3.0"],
    },
)

#!/usr/bin/env bash
# conda-recipe/build.sh — build script for the cat-surf conda package.
#
# This script is run by conda-build inside the activated conda build
# environment.  It performs two steps:
#
#   1. Build libCAT.a (+ bundled libfftw3.a) from the C source tree using
#      the autotools build system.  The conda compilers and flags are
#      inherited automatically via CC, CFLAGS, etc.
#
#   2. Build and install the Cython extension by running
#      `pip install cat_surface_cython/` pointing setup.py at the freshly
#      built libCAT.a.
#
# Environment variables set by conda-build that are used here:
#   SRC_DIR   — unpacked source tree (= repo root from the GitHub archive)
#   PREFIX    — conda install prefix (libraries, headers go here)
#   CPU_COUNT — number of CPUs available for parallel make
#   CC, CXX, CFLAGS, CXXFLAGS, LDFLAGS — conda toolchain flags

set -euxo pipefail

# ---------------------------------------------------------------------------
# Step 1 — build libCAT.a inside a subdirectory to keep SRC_DIR clean
# ---------------------------------------------------------------------------
CAT_BUILD="${SRC_DIR}/build-conda"
mkdir -p "${CAT_BUILD}"
cd "${CAT_BUILD}"

autoupdate || true
autoreconf -fi "${SRC_DIR}"

"${SRC_DIR}/configure" \
    --srcdir="${SRC_DIR}" \
    --enable-static \
    --disable-shared \
    --disable-doc \
    CFLAGS="${CFLAGS:--O2} -fPIC" \
    CXXFLAGS="${CXXFLAGS:--O2} -fPIC" \
    CPPFLAGS="${CPPFLAGS:-}"

make -j"${CPU_COUNT:-1}" libCAT.la

# ---------------------------------------------------------------------------
# Step 2 — build the Cython extension against the freshly built libCAT.a
#
# On macOS, CAT_LIBOMP_PREFIX points setup.py to conda-forge's llvm-openmp
# (installed in $PREFIX by the host dependency).  conda-build then patches
# the baked-in @rpath entries to $PREFIX/lib via install_name_tool so
# dyld resolves libomp.dylib from the active conda environment at runtime.
# ---------------------------------------------------------------------------
export CAT_SURFACE_ROOT="${SRC_DIR}"
export CAT_BUILD_DIR="${CAT_BUILD}"

if [[ "$(uname)" == "Darwin" ]]; then
    export CAT_LIBOMP_PREFIX="${PREFIX}"
fi

cd "${SRC_DIR}/cat_surface_cython"
"${PYTHON}" -m pip install . \
    --no-deps \
    --no-build-isolation \
    --prefix="${PREFIX}" \
    -vv

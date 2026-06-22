AC_INIT([Cortex Analysis Tools],[1.0],[Christian Gaser <christian.gaser@uni-jena.de>])
AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
LT_INIT

AC_ENABLE_STATIC
AC_DISABLE_SHARED

AC_CANONICAL_HOST

build_linux=no
build_windows=no
build_mac=no

# Detect the target system
build_emscripten=no
case "${host_os}" in
    linux*)
        build_linux=yes
        ;;
    cygwin*|mingw*)
        build_windows=yes
        ;;
    darwin*)
        build_mac=yes
        ;;
    *emscripten*)
        build_emscripten=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

AM_CONDITIONAL([EMSCRIPTEN], [test "$build_emscripten" = "yes"])

# Pass the conditionals to automake
AM_CONDITIONAL([LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([OSX], [test "$build_mac" = "yes"])

# Set the linker driver used in Progs/Makefile.am.
# On Windows, MeshFix (C++) is excluded from libCAT, so we can use the plain
# C linker (gcc) — this avoids pulling in libstdc++ and the winpthread chain.
# On all other platforms, libCAT includes C++ objects and must be linked with g++.
if test "$build_windows" = "yes"; then
    PROG_CCLD='$(CC)'
    PROG_LINK_TAG='CC'
else
    PROG_CCLD='$(CXX)'
    PROG_LINK_TAG='CXX'
fi
AC_SUBST([PROG_CCLD])
AC_SUBST([PROG_LINK_TAG])

# --- FFTW3 (bundled source in 3rdparty/fftw/) ---
# Configure the bundled FFTW3 directly from this configure script so that
# autoreconf does not recurse into the FFTW3 source tree (its Makefile.am
# uses GNU-make extensions that cause automake warnings/errors).
fftw_abs_srcdir=$(cd "$srcdir" && pwd)/3rdparty/fftw
FFTW3_CFLAGS="-I${fftw_abs_srcdir}/api"
FFTW3_LIBS='$(top_builddir)/3rdparty/fftw-build/.libs/libfftw3.a'
AC_SUBST([FFTW3_CFLAGS])
AC_SUBST([FFTW3_LIBS])

# Run FFTW3's own pre-generated configure in a build sub-directory.
# Pass the same CC, CFLAGS, host/build triplets for cross-compilation.
AC_MSG_NOTICE([configuring bundled FFTW3...])
fftw_builddir="$ac_pwd/3rdparty/fftw-build"
AS_MKDIR_P(["$fftw_builddir"])
fftw_configure_args="--enable-static --disable-shared --disable-fortran --disable-doc"
test -n "$host_alias" && fftw_configure_args="$fftw_configure_args --host=$host_alias"
test -n "$build_alias" && fftw_configure_args="$fftw_configure_args --build=$build_alias"
(cd "$fftw_builddir" && \
 "$fftw_abs_srcdir/configure" \
   CC="$CC" CFLAGS="$CFLAGS" \
   $fftw_configure_args --quiet) || \
  AC_MSG_ERROR([bundled FFTW3 configure failed])

AC_CHECK_HEADERS(values.h malloc.h)

# Parallelism is provided by pthreads (see Lib/CAT_Vol.c); OpenMP is no
# longer used, so no OpenMP runtime detection is needed here.

# On older glibc (< 2.17), clock_gettime lives in librt; pull in -lrt when
# the symbol is not already available in libc.
AC_SEARCH_LIBS([clock_gettime], [rt])

AC_CONFIG_FILES([Makefile Progs/Makefile tests/Makefile])
AC_OUTPUT

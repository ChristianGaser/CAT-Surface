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

# --- OpenMP detection (optional) ---
# Try the standard AC_OPENMP first.  On macOS with Apple Clang this will
# fail because Apple ships Clang without -fopenmp support.  In that case
# we look for the bundled static libomp in 3rdparty/libomp/ (arch-specific),
# then fall back to /usr/local/lib.
AC_OPENMP

AS_IF([test "x$ac_cv_prog_c_openmp" = "xunsupported" && test "$build_mac" = "yes"],
  [AC_MSG_CHECKING([for static libomp on macOS])
   omp_found=no
   # Derive arch from host_cpu for the arch-specific library name
   case "${host_cpu}" in
     x86_64|amd64) omp_arch=x86_64 ;;
     aarch64|arm64|arm) omp_arch=arm64 ;;
     *) omp_arch="" ;;
   esac
   omp_srcdir=$(cd "$srcdir" && pwd)/3rdparty/libomp
   # 1) Bundled arch-specific library in 3rdparty/libomp/lib/
   if test -f "$omp_srcdir/lib/libomp-${omp_arch}.a"; then
     OPENMP_CFLAGS="-Xpreprocessor -fopenmp -I${omp_srcdir}/include"
     OPENMP_LIBS="${omp_srcdir}/lib/libomp-${omp_arch}.a"
     omp_found=yes
     AC_MSG_RESULT([yes (bundled $omp_arch)])
   fi
   # 2) System-wide static library
   AS_IF([test "$omp_found" = "no"],
     [for d in /usr/local/lib /usr/lib; do
        if test -f "$d/libomp.a"; then
          OPENMP_CFLAGS="-Xpreprocessor -fopenmp"
          OPENMP_LIBS="$d/libomp.a"
          omp_found=yes
          AC_MSG_RESULT([yes ($d/libomp.a)])
          break
        fi
      done])
   AS_IF([test "$omp_found" = "no"],
     [OPENMP_CFLAGS=""
      OPENMP_LIBS=""
      AC_MSG_RESULT([no – OpenMP disabled])])])

# When AC_OPENMP succeeds on Linux/GCC the flags are in $OPENMP_CFLAGS and
# the library is linked implicitly.  Provide an empty OPENMP_LIBS fallback.
AS_IF([test "x$OPENMP_LIBS" = "x"], [OPENMP_LIBS=""])

# On GCC (Linux, Windows/Cygwin/MinGW), the compiler driver links libgomp
# automatically when -fopenmp is passed at link time.  Propagate the flag to
# OPENMP_LIBS so it reaches LDADD.
if test -n "$OPENMP_CFLAGS" && test -z "$OPENMP_LIBS"; then
  OPENMP_LIBS="$OPENMP_CFLAGS"
fi

AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_LIBS])

# On older glibc (< 2.17), clock_gettime lives in librt.  libgomp (OpenMP)
# needs it, so pull in -lrt when required.
AC_SEARCH_LIBS([clock_gettime], [rt])

AC_CONFIG_FILES([Makefile Progs/Makefile tests/Makefile])
AC_OUTPUT

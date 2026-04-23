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

# --- OpenMP detection (bundled-static only) ---
# Policy: enable OpenMP only when a matching static runtime archive is present
# in 3rdparty and links successfully for the current target.
AC_OPENMP

OPENMP_BUNDLED_FOUND=no
OPENMP_BUNDLED_CFLAGS=""
OPENMP_BUNDLED_LIBS=""

case "${host_cpu}" in
  x86_64|amd64) omp_arch=x86_64 ;;
  aarch64|arm64|arm) omp_arch=arm64 ;;
  *) omp_arch="" ;;
esac

# Windows: libstdc++ (posix-threads model in MSYS2 MinGW64) and libgomp both
# reference winpthread symbols.  The binaries.yml CI configure step passes
# -static in LDFLAGS so the linker automatically prefers .a over .dll.a for
# every -l reference, including ones injected by the C++ driver after LDADD.
# No extra detection is needed here.

if test "$build_mac" = "yes"; then
  omp_srcdir=$(cd "$srcdir" && pwd)/3rdparty/libomp
  if test -n "$omp_arch" && test -f "$omp_srcdir/lib/libomp-${omp_arch}.a"; then
    OPENMP_BUNDLED_CFLAGS="-Xpreprocessor -fopenmp -I${omp_srcdir}/include"
    OPENMP_BUNDLED_LIBS="${omp_srcdir}/lib/libomp-${omp_arch}.a"
    OPENMP_BUNDLED_FOUND=yes
  fi
elif test "$build_windows" = "yes"; then
  # On Windows/MinGW, libgomp.a is used directly from the active compiler.
  # Static linking of libgomp (and its libwinpthread dependency) is achieved
  # by passing -static in LDFLAGS at configure time (see binaries.yml), which
  # makes ld prefer .a over .dll.a for every library reference globally.
  gomp_lib=`$CC -print-file-name=libgomp.a 2>/dev/null`
  if test -n "$omp_arch" && test "$gomp_lib" != "libgomp.a" && test -f "$gomp_lib"; then
    OPENMP_BUNDLED_CFLAGS="-fopenmp"
    OPENMP_BUNDLED_LIBS="$gomp_lib"
    OPENMP_BUNDLED_FOUND=yes
  else
    AC_MSG_WARN([libgomp.a not found; OpenMP disabled on Windows])
  fi
elif test "$build_linux" = "yes"; then
  gomp_srcdir=$(cd "$srcdir" && pwd)/3rdparty/libgomp
  if test -n "$omp_arch" && test -f "$gomp_srcdir/lib/libgomp-${omp_arch}.a"; then
    OPENMP_BUNDLED_CFLAGS="-fopenmp"
    OPENMP_BUNDLED_LIBS="${gomp_srcdir}/lib/libgomp-${omp_arch}.a"
    OPENMP_BUNDLED_FOUND=yes
  fi
fi

if test "$OPENMP_BUNDLED_FOUND" = "yes"; then
  AC_MSG_NOTICE([found bundled OpenMP runtime: $OPENMP_BUNDLED_LIBS])
  AC_MSG_CHECKING([whether bundled OpenMP runtime links with real OpenMP code])
  # Use an actual OpenMP translation unit so the linker must resolve gomp
  # symbols from the archive.  A trivial int-main would not pull in any
  # archive objects and would silently succeed even with a non-PIC archive.
  cat > conftest_openmp_runtime.c <<'EOF'
#ifdef _MSC_VER
# include <omp.h>
#else
# include <omp.h>
#endif
int main(void)
{
    int n = 0;
#pragma omp parallel reduction(+:n)
    n += omp_get_num_threads();
    return n == 0;
}
EOF
  # Do not pass -pthread on Windows: MinGW maps it to -lwinpthread which
  # would re-introduce the dynamic libwinpthread-1.dll dependency we removed.
  openmp_test_pthread_flag=""
  test "$build_windows" != "yes" && openmp_test_pthread_flag="-pthread"
  if $CC $CPPFLAGS $CFLAGS $OPENMP_BUNDLED_CFLAGS $LDFLAGS \
        conftest_openmp_runtime.c $OPENMP_BUNDLED_LIBS $openmp_test_pthread_flag \
        -o conftest_openmp_runtime$ac_exeext >/dev/null 2>&1; then
    OPENMP_CFLAGS="$OPENMP_BUNDLED_CFLAGS"
    OPENMP_LIBS="$OPENMP_BUNDLED_LIBS"
    AC_MSG_RESULT([yes])
  else
    OPENMP_CFLAGS=""
    OPENMP_LIBS=""
    AC_MSG_RESULT([no])
    AC_MSG_WARN([bundled OpenMP runtime failed to link; OpenMP disabled for this target])
  fi
  rm -f conftest_openmp_runtime.c conftest_openmp_runtime$ac_exeext
else
  OPENMP_CFLAGS=""
  OPENMP_LIBS=""
  AC_MSG_NOTICE([no bundled static OpenMP runtime for this target; OpenMP disabled])
fi

AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_LIBS])

# On older glibc (< 2.17), clock_gettime lives in librt.  libgomp (OpenMP)
# needs it, so pull in -lrt when required.
AC_SEARCH_LIBS([clock_gettime], [rt])

AC_CONFIG_FILES([Makefile Progs/Makefile tests/Makefile])
AC_OUTPUT

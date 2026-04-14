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

AC_CHECK_LIB(fftw3,fftw_execute,,AC_MSG_ERROR(The fftw3 library is required.))

AC_CHECK_HEADERS(values.h malloc.h)

# --- OpenMP detection (optional) ---
# Try the standard AC_OPENMP first.  On macOS with Apple Clang this will
# fail because Apple ships Clang without -fopenmp support.  In that case
# we fall back to the Homebrew-installed libomp (brew install libomp).
AC_OPENMP

AS_IF([test "x$ac_cv_prog_c_openmp" = "xunsupported" && test "$build_mac" = "yes"],
  [AC_MSG_CHECKING([for Homebrew libomp on macOS])
   omp_prefix=""
   for d in /opt/homebrew /usr/local; do
     if test -f "$d/opt/libomp/lib/libomp.dylib"; then
       omp_prefix="$d/opt/libomp"
       break
     fi
   done
   AS_IF([test -n "$omp_prefix"],
     [OPENMP_CFLAGS="-Xpreprocessor -fopenmp -I${omp_prefix}/include"
      OPENMP_LIBS="-L${omp_prefix}/lib -lomp"
      AC_MSG_RESULT([yes ($omp_prefix)])],
     [OPENMP_CFLAGS=""
      OPENMP_LIBS=""
      AC_MSG_RESULT([no – install with: brew install libomp])])])

# When AC_OPENMP succeeds on Linux/GCC the flags are in $OPENMP_CFLAGS and
# the library is linked implicitly.  Provide an empty OPENMP_LIBS fallback.
AS_IF([test "x$OPENMP_LIBS" = "x"], [OPENMP_LIBS=""])

AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_LIBS])

AC_CONFIG_FILES([Makefile Progs/Makefile tests/Makefile])
AC_OUTPUT

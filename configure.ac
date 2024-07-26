dnl $Id$
dnl Process this file with autoconf to produce a configure script.

AC_INIT([Cortex Analysis Tools],[1.0],[Christian Gaser <christian.gaser@uni-jena.de>])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
LT_INIT

AC_ENABLE_STATIC
AC_DISABLE_SHARED

smr_WITH_BUILD_PATH

AC_CANONICAL_HOST

build_linux=no
build_windows=no
build_mac=no

# Detect the target system
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
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

# Pass the conditionals to automake
AM_CONDITIONAL([LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([OSX], [test "$build_mac" = "yes"])

AC_CHECK_LIB(fftw3,fftw_execute,,AC_MSG_ERROR(The fftw3 library is required.))
AC_CHECK_LIB(expat, XML_ParserCreate,, AC_MSG_ERROR(The expat library is required.))

# AC_CHECK_HEADERS(expat.h)
# AC_CHECK_HEADERS(fftw3.h)
AC_CHECK_HEADERS(values.h malloc.h)

AC_CHECK_HEADERS()

AC_CONFIG_FILES( Makefile Progs/Makefile)
AC_OUTPUT

dnl $Id$
dnl Process this file with autoconf to produce a configure script.

AC_INIT([Cortex Analysis Tools], 0.9, Christian Gaser <christian.gaser@uni-jena.de>)
AM_INIT_AUTOMAKE

AC_CONFIG_HEADER(config.h)

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
LT_INIT

AC_PROG_LIBTOOL
AC_ENABLE_STATIC
AC_DISABLE_SHARED

smr_WITH_BUILD_PATH
mni_REQUIRE_BICPL

AC_CHECK_LIB(fftw3,fftw_execute,,AC_MSG_ERROR(The fftw3 library is required.))
AC_CHECK_LIB(expat, XML_ParserCreate,, AC_MSG_ERROR(The expat library is required.))

AC_CHECK_HEADERS(expat.h)
AC_CHECK_HEADERS(fftw3.h)
AC_CHECK_HEADERS(values.h malloc.h)

AC_CHECK_HEADERS()

AC_CONFIG_FILES( Makefile Progs/Makefile)
AC_CONFIG_MACRO_DIR([m4])
AC_OUTPUT

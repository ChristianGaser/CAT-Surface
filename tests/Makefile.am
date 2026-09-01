AM_CPPFLAGS = -I$(top_srcdir)/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/volume_io/Include \
    -I$(top_srcdir)/3rdparty/gifticlib \
    -I$(top_srcdir)/3rdparty/nii2mesh \
    -I$(top_srcdir)/3rdparty/nifti

check_PROGRAMS = test_math test_glmformula test_sheetness test_calc test_surfinfo
TESTS = test_math test_glmformula test_sheetness test_calc test_surfinfo

test_math_SOURCES = test_math.c

test_math_LDADD = -lm

test_glmformula_SOURCES = test_glmformula.c

test_glmformula_LDADD = -lm

test_sheetness_SOURCES = test_sheetness.c

test_sheetness_LDADD = $(top_builddir)/libCAT.la -lm

test_calc_SOURCES = test_calc.c

test_calc_LDADD = $(top_builddir)/libCAT.la -lm

test_surfinfo_SOURCES = test_surfinfo.c

test_surfinfo_LDADD = $(top_builddir)/libCAT.la -lm

AM_CPPFLAGS = -I$(top_srcdir)/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/volume_io/Include \
    -I$(top_srcdir)/3rdparty/gifticlib \
    -I$(top_srcdir)/3rdparty/nii2mesh \
    -I$(top_srcdir)/3rdparty/nifti

check_PROGRAMS = test_math test_glmformula test_sheetness test_calc test_surfinfo test_pialprofile test_deform test_folding test_intersect
TESTS = test_math test_glmformula test_sheetness test_calc test_surfinfo test_pialprofile test_deform test_folding test_intersect

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

test_pialprofile_SOURCES = test_pialprofile.c test_phantoms.h

test_pialprofile_LDADD = $(top_builddir)/libCAT.la -lm

test_deform_SOURCES = test_deform.c test_phantoms.h

test_deform_LDADD = $(top_builddir)/libCAT.la -lm

test_folding_SOURCES = test_folding.c

test_folding_LDADD = $(top_builddir)/libCAT.la -lm

test_intersect_SOURCES = test_intersect.c

test_intersect_LDADD = $(top_builddir)/libCAT.la -lm

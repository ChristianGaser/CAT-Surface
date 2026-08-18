AM_CPPFLAGS = -I$(top_srcdir)/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/nifti

check_PROGRAMS = test_math test_glmformula test_sheetness
TESTS = test_math test_glmformula test_sheetness

test_math_SOURCES = test_math.c

test_math_LDADD = -lm

test_glmformula_SOURCES = test_glmformula.c

test_glmformula_LDADD = -lm

test_sheetness_SOURCES = test_sheetness.c

test_sheetness_LDADD = $(top_builddir)/libCAT.la -lm

AM_CPPFLAGS = -I$(top_srcdir)/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/nifti

check_PROGRAMS = test_math test_glmformula test_sheetness test_calc test_myelin
TESTS = test_math test_glmformula test_sheetness test_calc test_myelin

test_math_SOURCES = test_math.c

test_math_LDADD = -lm

test_glmformula_SOURCES = test_glmformula.c

test_glmformula_LDADD = -lm

test_sheetness_SOURCES = test_sheetness.c

test_sheetness_LDADD = $(top_builddir)/libCAT.la -lm

test_calc_SOURCES = test_calc.c

test_calc_LDADD = $(top_builddir)/libCAT.la -lm

test_myelin_SOURCES = test_myelin.c

test_myelin_LDADD = $(top_builddir)/libCAT.la -lm

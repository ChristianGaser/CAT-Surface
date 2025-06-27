AM_CPPFLAGS = -I$(top_srcdir)/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/nifti

check_PROGRAMS = test_math
TESTS = test_math

test_math_SOURCES = test_math.c

test_math_LDADD = -lm

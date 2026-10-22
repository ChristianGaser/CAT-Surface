AM_CPPFLAGS = -I$(top_srcdir)/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/volume_io/Include \
    -I$(top_srcdir)/3rdparty/gifticlib \
    -I$(top_srcdir)/3rdparty/expat/lib \
    -I$(top_srcdir)/3rdparty/nii2mesh \
    -I$(top_srcdir)/3rdparty/nifti

# Included by the tests rather than compiled on their own.
EXTRA_DIST = minunit.h math_subset.c

check_PROGRAMS = test_math test_glmformula test_sheetness test_calc test_surfinfo test_pialprofile test_deform test_folding test_intersect test_annot test_morph test_resample test_nifti test_volaverage test_bspline test_volumereg test_volrealign test_sulc test_separate
TESTS = test_math test_glmformula test_sheetness test_calc test_surfinfo test_pialprofile test_deform test_folding test_intersect test_annot test_morph test_resample test_nifti test_volaverage test_bspline test_volumereg test_volrealign test_sulc test_separate

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

test_sulc_SOURCES = test_sulc.c

test_sulc_LDADD = $(top_builddir)/libCAT.la -lm

test_separate_SOURCES = test_separate.c

test_separate_LDADD = $(top_builddir)/libCAT.la -lm

test_intersect_SOURCES = test_intersect.c

test_intersect_LDADD = $(top_builddir)/libCAT.la -lm

test_annot_SOURCES = test_annot.c

test_annot_LDADD = $(top_builddir)/libCAT.la -lm

test_morph_SOURCES = test_morph.c

test_morph_LDADD = $(top_builddir)/libCAT.la -lm

test_resample_SOURCES = test_resample.c

test_resample_LDADD = $(top_builddir)/libCAT.la -lm

test_nifti_SOURCES = test_nifti.c

test_nifti_LDADD = $(top_builddir)/libCAT.la -lm

test_volaverage_SOURCES = test_volaverage.c

test_volaverage_LDADD = $(top_builddir)/libCAT.la -lm

test_bspline_SOURCES = test_bspline.c

test_bspline_LDADD = $(top_builddir)/libCAT.la -lm

test_volumereg_SOURCES = test_volumereg.c

test_volumereg_LDADD = $(top_builddir)/libCAT.la -lm

test_volrealign_SOURCES = test_volrealign.c

test_volrealign_LDADD = $(top_builddir)/libCAT.la -lm

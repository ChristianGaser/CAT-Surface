SUBDIRS = . Progs

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/MarchingCubes -I$(top_srcdir)/niftilib -I$(top_srcdir)/zlib

lib_LTLIBRARIES = libCAT.la

EXTRA_DIST = \
    autogen.sh \
    MarchingCubes/README \
    gifticlib/LICENSE.gifti \
    gifticlib/README.gifti \
    $(m4_files)

s2kit10_SRC = \
    s2kit10/FST_semi_fly.c \
    s2kit10/cospmls.c \
    s2kit10/makeweights.c \
    s2kit10/pmls.c \
    s2kit10/seminaive.c \
    s2kit10/naive_synthesis.c \
    s2kit10/primitive.c
    
nl_SRC = \
    nl/nl_api.c \
    nl/nl_blas.c \
    nl/nl_context.c \
    nl/nl_iterative_solvers.c \
    nl/nl_matrix.c \
    nl/nl_os.c \
    nl/nl_preconditioners.c \
    nl/nl_superlu.c
    
zlib_SRC = \
    zlib/adler32.c \
    zlib/crc32.c \
    zlib/gzio.c  \
    zlib/inffast.c \
    zlib/inftrees.c \
    zlib/trees.c \
    zlib/zutil.c \
    zlib/compress.c \
    zlib/deflate.c \
    zlib/infback.c \
    zlib/inflate.c \
    zlib/uncompr.c

nifti_SRC = \
    niftilib/nifti1_io.c \
    niftilib/znzlib.c

gifti_SRC = \
    gifticlib/gifti_io.c \
    gifticlib/gifti_tool.c \
    gifticlib/gifti_xml.c

dartel_SRC = \
    dartel/diffeosphere.c \
    dartel/diffeopoly.c \
    dartel/optimizersphere.c
    
noinst_HEADERS = \
    Include/CAT_Map.h \
    Include/CAT_Smooth.h \
    Include/CAT_Curvature.h \
    Include/CAT_Surf.h \
    Include/CAT_Pinv.h \
    Include/CAT_SurfaceIO.h \
    Include/CAT_SPH.h \
    Include/CAT_Interpolate.h \
    Include/CAT_Isometrics.h \
    Include/CAT_Patch.h \
    Include/CAT_Resample.h \
    Include/CAT_Complexity.h \
    Include/CAT_Defect.h \
    Include/CAT_Intersect.h \
    Include/CAT_FixTopology.h \
    Include/CAT_Refine.h \
    Include/CAT_DeformPolygons.h \
    Include/CAT_Separate.h \
    Include/CAT_Octree.h \
    Include/CAT_NiftiIO.h \
    Include/CAT_ConvexHull.h \
    Include/CAT_Gyrification.h \
    Include/CAT_DepthPotential.h \
    MarchingCubes/MarchingCubes.h \
    MarchingCubes/LookUpTable.h \
    MarchingCubes/genus0.h \
    MarchingCubes/tricases.h \
    dartel/dartel.h \
    niftilib/analyze75.h \
    niftilib/nifti1.h \
    niftilib/nifti1_io.h \
    niftilib/nifti1_local.h \
    niftilib/znzlib.h \
    gifticlib/gifti_io.h \
    gifticlib/gifti_test.h \
    gifticlib/gifti_tool.h \
    gifticlib/gifti_xml.h \
    zlib/crc32.h \
    zlib/deflate.h \
    zlib/inffast.h \
    zlib/inffixed.h \
    zlib/inflate.h \
    zlib/inftrees.h \
    zlib/trees.h \
    zlib/zconf.h \
    zlib/zlib.h \
    zlib/zutil.h \
    s2kit10/FST_semi_fly.h \
    s2kit10/naive_synthesis.h \
    s2kit10/primitive.h \
    s2kit10/cospmls.h \
    s2kit10/makeweights.h \
    s2kit10/pmls.h \
    s2kit10/seminaive.h \
    nl/nl.h \
    nl/nl_context.h \
    nl/nl_matrix.h \
    nl/nl_private.h \
    nl/nl_blas.h \
    nl/nl_iterative_solvers.h \
    nl/nl_preconditioners.h \
    nl/nl_superlu.h

libCAT_la_SOURCES = \
    Lib/CAT_Surf.c \
    Lib/CAT_Curvature.c \
    Lib/CAT_SurfaceIO.c \
    Lib/CAT_Map.c \
    Lib/CAT_Smooth.c \
    Lib/CAT_SPH.c \
    Lib/CAT_Pinv.c \
    Lib/CAT_SelfIntersect.c \
    Lib/CAT_Intersect.c \
    Lib/CAT_Octree.c \
    Lib/CAT_Resample.c \
    Lib/CAT_DeformPolygons.c \
    Lib/CAT_Topology.c \
    Lib/CAT_Defect.c \
    Lib/CAT_Patch.c \
    Lib/CAT_Complexity.c \
    Lib/CAT_Refine.c \
    Lib/CAT_Interpolate.c \
    Lib/CAT_LaplaceBeltrami.c \
    Lib/CAT_Isometrics.c \
    Lib/CAT_Separate.c \
    Lib/CAT_NiftiIO.c \
    Lib/CAT_ConvexHull.c \
    Lib/CAT_Gyrification.c \
    Lib/CAT_DepthPotential.c \
    MarchingCubes/MarchingCubes.c \
    MarchingCubes/genus0.c \
    $(zlib_SRC) \
    $(nl_SRC) \
    $(nifti_SRC) \
    $(gifti_SRC) \
    $(dartel_SRC) \
    $(s2kit10_SRC)

m4_files = m4/mni_REQUIRE_LIB.m4 \
    m4/mni_REQUIRE_MNILIBS.m4 \
    m4/mni_REQUIRE_OPENINVENTOR.m4 \
    m4/mni_cxx_have_koenig_lookup.m4 \
    m4/smr_CGAL_MAKEFILE.m4 \
    m4/smr_OPTIONAL_LIB.m4 \
    m4/smr_REQUIRED_LIB.m4 \
    m4/smr_WITH_BUILD_PATH.m4

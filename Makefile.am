SUBDIRS = . Progs

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS = -pthread -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/3rdparty/MarchingCubes -I$(top_srcdir)/3rdparty/nifti -I$(top_srcdir)/3rdparty/zlib
AM_CPPFLAGS = -pthread -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/3rdparty/MarchingCubes -I$(top_srcdir)/3rdparty/nifti -I$(top_srcdir)/3rdparty/zlib

lib_LTLIBRARIES = libCAT.la

EXTRA_DIST = \
    autogen.sh \
    3rdparty/MarchingCubes/README \
    3rdparty/gifticlib/LICENSE.gifti \
    3rdparty/gifticlib/README.gifti \
    $(m4_files)

s2kit10_SRC = \
    3rdparty/s2kit10/FST_semi_fly.c \
    3rdparty/s2kit10/cospmls.c \
    3rdparty/s2kit10/makeweights.c \
    3rdparty/s2kit10/pmls.c \
    3rdparty/s2kit10/seminaive.c \
    3rdparty/s2kit10/naive_synthesis.c \
    3rdparty/s2kit10/primitive.c
    
zlib_SRC = \
    3rdparty/zlib/adler32.c \
    3rdparty/zlib/crc32.c \
    3rdparty/zlib/gzio.c  \
    3rdparty/zlib/inffast.c \
    3rdparty/zlib/inftrees.c \
    3rdparty/zlib/trees.c \
    3rdparty/zlib/zutil.c \
    3rdparty/zlib/compress.c \
    3rdparty/zlib/deflate.c \
    3rdparty/zlib/infback.c \
    3rdparty/zlib/inflate.c \
    3rdparty/zlib/uncompr.c

nifti_SRC = \
    3rdparty/nifti/nifti1_io.c \
    3rdparty/nifti/znzlib.c

gifti_SRC = \
    3rdparty/gifticlib/gifti_io.c \
    3rdparty/gifticlib/gifti_tool.c \
    3rdparty/gifticlib/gifti_xml.c

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
    Include/CAT_Separate.h \
    Include/CAT_Octree.h \
    Include/CAT_ConvexHull.h \
    Include/CAT_Gyrification.h \
    Include/CAT_DepthPotential.h \
    Include/CAT_Vol.h \
    Include/CAT_Kmeans.h \
    Include/CAT_Amap.h \
    Include/CAT_NiftiLib.h \
    3rdparty/MarchingCubes/MarchingCubes.h \
    3rdparty/MarchingCubes/LookUpTable.h \
    3rdparty/MarchingCubes/genus0.h \
    3rdparty/MarchingCubes/tricases.h \
    dartel/dartel.h \
    3rdparty/nifti/analyze75.h \
    3rdparty/nifti/nifti1.h \
    3rdparty/nifti/nifti1_io.h \
    3rdparty/nifti/nifti1_local.h \
    3rdparty/nifti/znzlib.h \
    3rdparty/gifticlib/gifti_io.h \
    3rdparty/gifticlib/gifti_test.h \
    3rdparty/gifticlib/gifti_tool.h \
    3rdparty/gifticlib/gifti_xml.h \
    3rdparty/zlib/crc32.h \
    3rdparty/zlib/deflate.h \
    3rdparty/zlib/inffast.h \
    3rdparty/zlib/inffixed.h \
    3rdparty/zlib/inflate.h \
    3rdparty/zlib/inftrees.h \
    3rdparty/zlib/trees.h \
    3rdparty/zlib/zconf.h \
    3rdparty/zlib/zlib.h \
    3rdparty/zlib/zutil.h \
    3rdparty/s2kit10/FST_semi_fly.h \
    3rdparty/s2kit10/naive_synthesis.h \
    3rdparty/s2kit10/primitive.h \
    3rdparty/s2kit10/cospmls.h \
    3rdparty/s2kit10/makeweights.h \
    3rdparty/s2kit10/pmls.h \
    3rdparty/s2kit10/seminaive.h

libCAT_la_SOURCES = \
    Lib/CAT_Surf.c \
    Lib/CAT_Curvature.c \
    Lib/CAT_SurfaceIO.c \
    Lib/CAT_Map.c \
    Lib/CAT_Smooth.c \
    Lib/CAT_SPH.c \
    Lib/CAT_Pinv.c \
    Lib/CAT_Intersect.c \
    Lib/CAT_Octree.c \
    Lib/CAT_Resample.c \
    Lib/CAT_Topology.c \
    Lib/CAT_Defect.c \
    Lib/CAT_Patch.c \
    Lib/CAT_Complexity.c \
    Lib/CAT_Refine.c \
    Lib/CAT_Interpolate.c \
    Lib/CAT_Isometrics.c \
    Lib/CAT_Separate.c \
    Lib/CAT_ConvexHull.c \
    Lib/CAT_Gyrification.c \
    Lib/CAT_DepthPotential.c \
    Lib/CAT_Vol.c \
    Lib/CAT_Kmeans.c \
    Lib/CAT_Amap.c \
    Lib/CAT_NiftiLib.c \
    Lib/CAT_Sanlm.c \
    3rdparty/MarchingCubes/MarchingCubes.c \
    3rdparty/MarchingCubes/genus0.c \
    $(zlib_SRC) \
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

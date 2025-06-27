SUBDIRS = . Progs tests

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

CFLAGS0 = -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/dartel \
    -I$(top_srcdir)/3rdparty/genus0 -I$(top_srcdir)/3rdparty/nifti \
    -I$(top_srcdir)/3rdparty/nifti -I$(top_srcdir)/3rdparty/s2kit10 \
    -I$(top_srcdir)/3rdparty/gifticlib -I$(top_srcdir)/3rdparty/zlib \
    -I$(top_srcdir)/3rdparty/volume_io/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/MeshFix-V2.1/include/Kernel \
    -I$(top_srcdir)/3rdparty/MeshFix-V2.1/include/TMesh

if EMSCRIPTEN
CFLAGS0 += -s USE_ZLIB=1 -s MODULARIZE=1 -s FETCH=1 -s WASM=1 -s EXPORT_ES6=1 \
    -s FORCE_FILESYSTEM=1 -s USE_PTHREADS=1 -s ENVIRONMENT=web
endif

if WINDOWS
  AM_CFLAGS = ${CFLAGS0}
else
  AM_CFLAGS = -pthread ${CFLAGS0}
endif
AM_CPPFLAGS=${AM_CFLAGS}
AM_CXXFLAGS=${AM_CFLAGS} 

lib_LTLIBRARIES = libCAT.la

EXTRA_DIST = \
    autogen.sh \
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

MeshFix_SRC = \
    3rdparty/MeshFix-V2.1/src/TMesh/edge.cpp \
    3rdparty/MeshFix-V2.1/src/TMesh/vertex.cpp \
    3rdparty/MeshFix-V2.1/src/TMesh/triangle.cpp \
    3rdparty/MeshFix-V2.1/src/TMesh/tin.cpp \
    3rdparty/MeshFix-V2.1/src/Algorithms/detectIntersections.cpp \
    3rdparty/MeshFix-V2.1/src/Algorithms/holeFilling.cpp \
    3rdparty/MeshFix-V2.1/src/Algorithms/marchIntersections.cpp \
    3rdparty/MeshFix-V2.1/src/Algorithms/checkAndRepair.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/heap.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/matrix.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/orientation.c \
    3rdparty/MeshFix-V2.1/src/Kernel/list.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/coordinates.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/tmesh.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/graph.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/point.cpp \
    3rdparty/MeshFix-V2.1/src/Kernel/jqsort.cpp

expat_SRC = \
    3rdparty/expat/lib/xmlparse.c \
    3rdparty/expat/lib/xmlrole.c \
    3rdparty/expat/lib/xmltok.c  \
    3rdparty/expat/lib/xmltok_impl.c \
    3rdparty/expat/lib/xmltok_ns.c

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
    3rdparty/gifticlib/gifti_xml.c

dartel_SRC = \
    dartel/diffeosphere.c \
    dartel/optimizersphere.c \
    dartel/diffeopoly.c 

bic_SRC = \
    3rdparty/volume_io/Geometry/colour.c  \
    3rdparty/volume_io/Geometry/colour_def.c  \
    3rdparty/volume_io/Geometry/gaussian.c  \
    3rdparty/volume_io/Geometry/inverse.c  \
    3rdparty/volume_io/Geometry/newton.c  \
    3rdparty/volume_io/Geometry/points.c  \
    3rdparty/volume_io/Geometry/splines.c  \
    3rdparty/volume_io/Geometry/tensors.c  \
    3rdparty/volume_io/Geometry/transforms.c  \
    3rdparty/volume_io/Prog_utils/alloc.c  \
    3rdparty/volume_io/Prog_utils/alloc_check.c  \
    3rdparty/volume_io/Prog_utils/arrays.c  \
    3rdparty/volume_io/Prog_utils/files.c  \
    3rdparty/volume_io/Prog_utils/print.c  \
    3rdparty/volume_io/Prog_utils/progress.c  \
    3rdparty/volume_io/Prog_utils/string.c  \
    3rdparty/volume_io/Prog_utils/time.c  \
    3rdparty/bicpl-surface/ParseArgv.c  \
    3rdparty/bicpl-surface/Data_structures/bintree.c  \
    3rdparty/bicpl-surface/Data_structures/bitlist.c  \
    3rdparty/bicpl-surface/Data_structures/build_bintree.c  \
    3rdparty/bicpl-surface/Data_structures/hash2_table.c  \
    3rdparty/bicpl-surface/Data_structures/hash_table.c  \
    3rdparty/bicpl-surface/Data_structures/object_bintrees.c  \
    3rdparty/bicpl-surface/Data_structures/point_bintree.c  \
    3rdparty/bicpl-surface/Data_structures/ray_bintree.c  \
    3rdparty/bicpl-surface/Data_structures/skiplist.c  \
    3rdparty/bicpl-surface/Geometry/approx_sqrt.c  \
    3rdparty/bicpl-surface/Geometry/clip_3d.c  \
    3rdparty/bicpl-surface/Geometry/closest_point.c  \
    3rdparty/bicpl-surface/Geometry/curvature.c  \
    3rdparty/bicpl-surface/Geometry/flatten.c  \
    3rdparty/bicpl-surface/Geometry/geodesic_distance.c  \
    3rdparty/bicpl-surface/Geometry/geometry.c  \
    3rdparty/bicpl-surface/Geometry/intersect.c  \
    3rdparty/bicpl-surface/Geometry/line_circle.c  \
    3rdparty/bicpl-surface/Geometry/map_polygons.c  \
    3rdparty/bicpl-surface/Geometry/path_surface.c  \
    3rdparty/bicpl-surface/Geometry/plane_polygon_intersect.c  \
    3rdparty/bicpl-surface/Geometry/platonics.c  \
    3rdparty/bicpl-surface/Geometry/points.c  \
    3rdparty/bicpl-surface/Geometry/poly_dist.c  \
    3rdparty/bicpl-surface/Geometry/polygon_sphere.c  \
    3rdparty/bicpl-surface/Geometry/ray_intersect.c  \
    3rdparty/bicpl-surface/Geometry/segment_polygons.c  \
    3rdparty/bicpl-surface/Geometry/smooth_curvature.c  \
    3rdparty/bicpl-surface/Geometry/smooth_lines.c  \
    3rdparty/bicpl-surface/Geometry/smooth_polygons.c  \
    3rdparty/bicpl-surface/Geometry/solve_plane.c  \
    3rdparty/bicpl-surface/Geometry/subdivide_lines.c  \
    3rdparty/bicpl-surface/Geometry/subdivide_polygons.c  \
    3rdparty/bicpl-surface/Geometry/surface_area.c  \
    3rdparty/bicpl-surface/Geometry/tetrahedrons.c  \
    3rdparty/bicpl-surface/Geometry/tubes.c  \
    3rdparty/bicpl-surface/Marching_cubes/isosurfaces.c  \
    3rdparty/bicpl-surface/Marching_cubes/marching_cubes.c  \
    3rdparty/bicpl-surface/Marching_cubes/marching_no_holes.c  \
    3rdparty/bicpl-surface/Marching_cubes/marching_tetra.c  \
    3rdparty/bicpl-surface/Numerical/amoeba.c  \
    3rdparty/bicpl-surface/Numerical/gaussian.c  \
    3rdparty/bicpl-surface/Numerical/gradient_minimize.c  \
    3rdparty/bicpl-surface/Numerical/histogram.c  \
    3rdparty/bicpl-surface/Numerical/least_squares.c  \
    3rdparty/bicpl-surface/Numerical/matrix_svd.c  \
    3rdparty/bicpl-surface/Numerical/minimize_lsq.c  \
    3rdparty/bicpl-surface/Numerical/minimize_lsq_include.c  \
    3rdparty/bicpl-surface/Numerical/numerical.c  \
    3rdparty/bicpl-surface/Numerical/quadratic.c  \
    3rdparty/bicpl-surface/Numerical/real_quadratic.c  \
    3rdparty/bicpl-surface/Numerical/statistics.c  \
    3rdparty/bicpl-surface/Numerical/t_stat.c  \
    3rdparty/bicpl-surface/Objects/coalesce.c  \
    3rdparty/bicpl-surface/Objects/colours.c  \
    3rdparty/bicpl-surface/Objects/graphics_io.c  \
    3rdparty/bicpl-surface/Objects/lines.c  \
    3rdparty/bicpl-surface/Objects/markers.c  \
    3rdparty/bicpl-surface/Objects/models.c  \
    3rdparty/bicpl-surface/Objects/object_io.c  \
    3rdparty/bicpl-surface/Objects/objects.c  \
    3rdparty/bicpl-surface/Objects/pixels.c  \
    3rdparty/bicpl-surface/Objects/poly_neighs.c  \
    3rdparty/bicpl-surface/Objects/polygons.c  \
    3rdparty/bicpl-surface/Objects/quadmesh.c  \
    3rdparty/bicpl-surface/Objects/rgb_lookup.c  \
    3rdparty/bicpl-surface/Objects/tag_objects.c  \
    3rdparty/bicpl-surface/Objects/text.c  \
    3rdparty/bicpl-surface/Prog_utils/arguments.c  \
    3rdparty/bicpl-surface/Prog_utils/globals.c  \
    3rdparty/bicpl-surface/Prog_utils/time.c  \
    3rdparty/bicpl-surface/Transforms/compute_tps.c  \
    3rdparty/bicpl-surface/Transforms/compute_xfm.c  \
    3rdparty/bicpl-surface/Transforms/matrix_basics.c  \
    3rdparty/bicpl-surface/Transforms/optimize.c  \
    3rdparty/bicpl-surface/Transforms/procrustes.c  \
    3rdparty/bicpl-surface/Transforms/rotmat_to_ang.c  \
    3rdparty/bicpl-surface/Transforms/safe_compute_xfm.c  \
    3rdparty/bicpl-surface/Transforms/transform_io.c  \
    3rdparty/bicpl-surface/Transforms/transforms.c  \
    3rdparty/bicpl-surface/bicpl_clapack/F77_aloc.c  \
    3rdparty/bicpl-surface/bicpl_clapack/d_sign.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dbdsqr.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dcopy.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgebd2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgebrd.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgelq2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgelqf.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgemm.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgemv.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgeqr2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgeqrf.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dger.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dgesvd.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlabrd.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlacpy.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlamch.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlange.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlapy2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlarf.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlarfb.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlarfg.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlarft.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlartg.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlas2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlascl.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlaset.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasq1.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasq2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasq3.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasq4.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasq5.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasq6.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasr.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasrt.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlassq.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dlasv2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dnrm2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorg2r.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorgbr.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorgl2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorglq.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorgqr.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorm2r.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dormbr.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dorml2.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dormlq.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dormqr.c  \
    3rdparty/bicpl-surface/bicpl_clapack/drot.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dscal.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dswap.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dtrmm.c  \
    3rdparty/bicpl-surface/bicpl_clapack/dtrmv.c  \
    3rdparty/bicpl-surface/bicpl_clapack/exit_.c  \
    3rdparty/bicpl-surface/bicpl_clapack/ieeeck.c  \
    3rdparty/bicpl-surface/bicpl_clapack/ilaenv.c  \
    3rdparty/bicpl-surface/bicpl_clapack/lsame.c  \
    3rdparty/bicpl-surface/bicpl_clapack/pow_dd.c  \
    3rdparty/bicpl-surface/bicpl_clapack/pow_di.c  \
    3rdparty/bicpl-surface/bicpl_clapack/s_cat.c  \
    3rdparty/bicpl-surface/bicpl_clapack/s_cmp.c  \
    3rdparty/bicpl-surface/bicpl_clapack/s_copy.c  \
    3rdparty/bicpl-surface/bicpl_clapack/xerbla.c
    
noinst_HEADERS = \
    Include/CAT_Map.h \
    Include/CAT_Smooth.h \
    Include/CAT_Curvature.h \
    Include/CAT_Surf.h \
    Include/CAT_Warp.h \
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
    Include/CAT_Math.h \
    Include/CAT_Kmeans.h \
    Include/CAT_Amap.h \
    Include/CAT_Bmap.h \
    Include/CAT_NiftiLib.h \
    Include/CAT_MarchingCubes.h \
    3rdparty/bicpl-surface/Include/bicpl.h \
    3rdparty/genus0/LookUpTable.h \
    3rdparty/genus0/genus0.h \
    3rdparty/genus0/tricases.h \
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
    Lib/CAT_Math.c \
    Lib/CAT_Surf.c \
    Lib/CAT_Warp.c \
    Lib/CAT_Curvature.c \
    Lib/CAT_SurfaceIO.c \
    Lib/CAT_Map.c \
    Lib/CAT_Smooth.c \
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
    Lib/CAT_Bmap.c \
    Lib/CAT_NiftiLib.c \
    Lib/CAT_Sanlm.c \
    Lib/CAT_Deform.c \
    Lib//CAT_MarchingCubes.c \
    Lib/CAT_SPH.c \
    3rdparty/genus0/genus0.c \
    $(expat_SRC) \
    $(zlib_SRC) \
    $(nifti_SRC) \
    $(gifti_SRC) \
    $(dartel_SRC) \
    $(bic_SRC) \
    $(s2kit10_SRC)



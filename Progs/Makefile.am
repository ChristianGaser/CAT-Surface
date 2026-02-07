# Compiler flags
AM_CFLAGS = -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/3rdparty/dartel \
    -I$(top_srcdir)/3rdparty/genus0 -I$(top_srcdir)/3rdparty/nifti \
    -I$(top_srcdir)/3rdparty/nifti -I$(top_srcdir)/3rdparty/s2kit10 \
    -I$(top_srcdir)/3rdparty/gifticlib -I$(top_srcdir)/3rdparty/zlib \
    -I$(top_srcdir)/3rdparty/volume_io/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/nii2mesh \
    -I$(top_srcdir)/3rdparty/MeshFix-V2.1/include/Kernel \
    -I$(top_srcdir)/3rdparty/MeshFix-V2.1/include/TMesh

AM_CPPFLAGS = ${AM_CFLAGS}

# When libCAT includes C++ objects (MeshFix), programs must link with the C++
# runtime. In some cross toolchains the shared C++ runtime depends on
# libgcc_s.so.1 which may not be discoverable at link time. If a static
# libstdc++ is available, prefer static libstdc++/libgcc to avoid rpath-link
# issues.
CAT_LIBSTDCXX_A = $(shell $(CXX) -print-file-name=libstdc++.a)
ifneq ($(CAT_LIBSTDCXX_A),libstdc++.a)
AM_LDFLAGS += -static-libstdc++ -static-libgcc
endif

# Help the linker discover libgcc_s when a shared libstdc++ is pulled in via
# libtool's handling of libstdc++.la (common in cross toolchains).
CAT_LIBGCC_S_SO = $(shell $(CXX) -print-file-name=libgcc_s.so.1)
ifneq ($(CAT_LIBGCC_S_SO),libgcc_s.so.1)
CAT_LIBGCC_S_DIR = $(shell dirname $(CAT_LIBGCC_S_SO) 2>/dev/null)
AM_LDFLAGS += -Wl,-rpath-link,$(CAT_LIBGCC_S_DIR) -L$(CAT_LIBGCC_S_DIR)
endif

# libCAT contains C++ objects (MeshFix), so we must link programs with the C++
# linker to ensure libstdc++/unwind symbols are resolved correctly.
CCLD = $(CXX)
LINK = $(LIBTOOL) $(AM_V_lt) --tag=CXX $(AM_LIBTOOLFLAGS) \
    $(LIBTOOLFLAGS) --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
    $(AM_LDFLAGS) $(LDFLAGS) -o $@

# Linker flags
LDADD = ../libCAT.la -lm -lfftw3

# Distribute scripts
dist_bin_SCRIPTS = \
    CAT_ResampleFS_ui \
    CAT_FractalDimensionFS_ui \
    CAT_Resample_ui \
    CAT_FractalDimension_ui \
    CAT_DumpCurvL4FS_ui \
    CAT_DumpSurfArea_ui \
    CAT_SulcusDepth_ui \
    CAT_DumpCurv_ui

# Programs to build
bin_PROGRAMS = \
    CAT_Surf2Sheet \
    CAT_SurfSheet2Surf \
    CAT_Vol2Surf \
    CAT_SurfCurvature \
    CAT_SurfApplyWarp \
    CAT_SurfApplyWarpValues \
    CAT_SurfArea \
    CAT_SurfRatio \
    CAT_SurfSmooth \
    CAT_SurfPlotValuesAtPoint \
    CAT_SurfPlotValuesAtMaximum \
    CAT_SurfMeasure2Txt \
    CAT_SurfConvert \
    CAT_Surf2Sphere \
    CAT_SurfSeparateClusters \
    CAT_SurfSeparatePolygons \
    CAT_SurfResampleSpherical \
    CAT_SurfResample \
    CAT_SurfResampleMulti \
    CAT_SurfAreaDistortion \
    CAT_SurfAngularDistortion \
    CAT_SurfMetricDistortion \
    CAT_SurfAverage \
    CAT_SurfFixTopology \
    CAT_SurfIsometize \
    CAT_SurfWarp \
    CAT_SurfWarpDemon \
    CAT_SurfSmoothAreal \
    CAT_SurfSmoothLaplacian \
    CAT_SurfSmoothConvexity \
    CAT_SurfSmoothSharpness \
    CAT_SurfHausdorff \
    CAT_SurfSharpness \
    CAT_SurfSelfIntersect \
    CAT_SurfMarkDefects \
    CAT_SurfArtifacts \
    CAT_SurfFractalDimension \
    CAT_SurfRefine \
    CAT_SurfConvexity \
    CAT_VolMarchingCubes \
    CAT_SurfAddValues \
    CAT_SurfSulcusDepth \
    CAT_SurfDepthPotential \
    CAT_SurfCentral2Pial \
    CAT_Surf2ConvexHull \
    CAT_SurfDistance \
    CAT_Test \
    CAT_SurfSmoothDiffusion \
    CAT_VolAmap \
    CAT_VolAverage \
    CAT_VolSmooth \
    CAT_VolSanlm \
    CAT_VolOrnlm \
    CAT_VolThicknessPbt \
    CAT_VolLayerSmooth \
    CAT_SurfResampleSphericalSPH \
    CAT_Surf2SPH \
    CAT_SurfDeform \
    CAT_Surf2PialWhite \
    CAT_SurfCorrectThicknessFolding \
    CAT_SurfReduce \
    CAT_SurfRemoveIntersections \
    CAT_Surf2ROIMulti \
    CAT_SurfSPH2Surf

# Programs requiring C++ standard library
CAT_SurfRemoveIntersections_LDADD = $(LDADD) -lstdc++

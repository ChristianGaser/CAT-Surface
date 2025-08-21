# Compiler flags
AM_CFLAGS = -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/3rdparty/dartel \
    -I$(top_srcdir)/3rdparty/genus0 -I$(top_srcdir)/3rdparty/nifti \
    -I$(top_srcdir)/3rdparty/nifti -I$(top_srcdir)/3rdparty/s2kit10 \
    -I$(top_srcdir)/3rdparty/gifticlib -I$(top_srcdir)/3rdparty/zlib \
    -I$(top_srcdir)/3rdparty/volume_io/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include \
    -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl \
    -I$(top_srcdir)/3rdparty/nii2mesh

AM_CPPFLAGS = ${AM_CFLAGS}

# Linker flags
LDADD = ../libCAT.la -lm -lfftw3 -pthread

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
    CAT_VolThicknessPbt \
    CAT_SurfResampleSphericalSPH \
    CAT_Surf2SPH \
    CAT_SurfDeform \
    CAT_Surf2PialWhite \
    CAT_SurfCorrectThicknessFolding \
    CAT_SurfReduce \
    CAT_SurfSPH2Surf

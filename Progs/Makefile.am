AM_CFLAGS = -DHAVE_ZLIB -I$(top_srcdir)/Include -I$(top_srcdir)/3rdparty/genus0 -I$(top_srcdir)/3rdparty/nifti -I$(top_srcdir)/3rdparty/zlib -I$(top_srcdir)/3rdparty/volume_io/Include -I$(top_srcdir)/3rdparty/bicpl-surface/Include -I$(top_srcdir)/3rdparty/bicpl-surface/Include/bicpl
AM_CPPFLAGS=${AM_CFLAGS}

LDADD = ../libCAT.la -lm

dist_bin_SCRIPTS = \
    CAT_ResampleFS_ui \
    CAT_FractalDimensionFS_ui \
    CAT_Resample_ui \
    CAT_FractalDimension_ui \
    CAT_DumpCurvL4FS_ui \
    CAT_DumpSurfArea_ui \
    CAT_SulcusDepth_ui \
    CAT_DumpCurv_ui
    
bin_PROGRAMS = \
    CAT_Surf2Sheet \
    CAT_Sheet2Surf \
    CAT_3dVol2Surf \
    CAT_DumpCurv \
    CAT_ApplyWarpSurf \
    CAT_ApplyWarpValues \
    CAT_DumpSurfArea \
    CAT_DumpSurfaceRatio \
    CAT_BlurSurfHK \
    CAT_PlotValuesAtPoint \
    CAT_PlotValuesAtMaximum \
    CAT_FreesurferCurv2Txt \
    CAT_SurfConvert \
    CAT_Surf2Sphere \
    CAT_SeparateCluster \
    CAT_SeparatePolygon \
    CAT_ResampleSphericalSurf \
    CAT_ResampleSurf \
    CAT_AreaDistortion \
    CAT_AngularDistortion \
    CAT_MetricDistortion \
    CAT_AverageSurfaces \
    CAT_FixTopology \
    CAT_Isometize \
    CAT_WarpSurf \
    CAT_WarpSurfDemon \
    CAT_ResampleSphericalSurfSPH \
    CAT_ArealSmoothing \
    CAT_SmoothConvexity \
    CAT_SmoothSharpness \
    CAT_Surf2SPH \
    CAT_SPH2Surf \
    CAT_Hausdorff \
    CAT_Sharpness \
    CAT_SelfIntersect \
    CAT_MarkDefects \
    CAT_Artifacts \
    CAT_FractalDimension \
    CAT_RefineMesh \
    CAT_DumpConvexity \
    CAT_VolMarchingCubes \
    CAT_AddValuesToSurf \
    CAT_SulcusDepth \
    CAT_DumpDepthPotential \
    CAT_Central2Pial \
    CAT_Surf2ConvexHull \
    CAT_SurfDistance \
    CAT_Test \
    CAT_SmoothDiffusion \
    CAT_VolAmap \
    CAT_VolAverage \
    CAT_VolSmooth \
    CAT_VolSanlm \
    CAT_VolThicknessPbt

# experimental code, not yet working or for testing only!
##  CAT_Test \
##  CAT_DumpGI \
##  CAT_DumpGyrificationIndex \ 
##  CAT_DeformMesh

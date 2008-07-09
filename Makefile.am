AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

dist_bin_SCRIPTS = \
	CAT_GLM_tmp \
	CAT_DumpGI_ui \
	CAT_DeformSurf_ui \
	CAT_Obj2Neurolens \
	CAT_DumpCurv_ui

bin_PROGRAMS = \
	CAT_DeformSurf \
	CAT_Surf2Sheet \
	CAT_Sheet2Surf \
	CAT_3dVol2Surf \
	CAT_DumpCurv \
	CAT_DumpGI \
	CAT_ApplyWarpSurf \
	CAT_ApplyWarpValues \
	CAT_Color2Surf \
	CAT_DumpSurfArea \
	CAT_DumpSurfaceRatio \
	CAT_BlurSurfHK \
	CAT_GLM_Estimate \
	CAT_PlotValuesAtPoint \
	CAT_PlotValuesAtMaximum \
	CAT_FreesurferCurv2Txt \
	CAT_SurfConvert \
	CAT_Surf2Sphere \
	CAT_SeparateCluster \
	CAT_ResampleSphericalSurf \
	CAT_AreaDistortion \
	CAT_AngularDistortion \
	CAT_AverageSurfaces \
	CAT_GLM_Math \
	CAT_FindFlats \
	CAT_Isometize \
	CAT_LBConformalMap \
	CAT_WarpSurf \
	CAT_ResampleSphericalSurfSPH \
	CAT_Surf2SPH \
	CAT_SPH2Surf

s2kit10_SRC = \
	s2kit10/FST_semi_fly.c \
	s2kit10/cospmls.c \
	s2kit10/makeweights.c \
	s2kit10/pmls.c \
	s2kit10/seminaive.c \
	s2kit10/csecond.c \
	s2kit10/naive_synthesis.c \
	s2kit10/primitive.c

extra_headers = \
	CAT_Deform.h \
	CAT_SheetIO.h \
	CAT_Map2d.h \
	CAT_Blur2d.h \
	CAT_Curvature.h \
	CAT_Surf.h \
	CAT_Pinv.h \
	CAT_SurfaceIO.h \
	CAT_SPH.h \
	CAT_DeformPrototypes.h	

EXTRA_DIST = $(extra_headers)

CAT_Surf2Sheet_SOURCES = CAT_Surf2Sheet.c CAT_Map2d.c CAT_Blur2d.c CAT_SheetIO.c CAT_SurfaceIO.c
CAT_Sheet2Surf_SOURCES = CAT_Sheet2Surf.c CAT_Map2d.c CAT_SheetIO.c CAT_SurfaceIO.c
CAT_3dVol2Surf_SOURCES = CAT_3dVol2Surf.c CAT_SurfaceIO.c
CAT_DumpCurv_SOURCES = CAT_DumpCurv.c CAT_Blur2d.c CAT_Curvature.c CAT_SurfaceIO.c
CAT_DumpGI_SOURCES = CAT_DumpGI.c CAT_Blur2d.c CAT_Curvature.c CAT_SurfaceIO.c
CAT_DeformSurf_SOURCES = CAT_DeformSurf.c CAT_DeformPolygons.c CAT_SurfaceIO.c
CAT_ApplyWarpSurf_SOURCES = CAT_ApplyWarpSurf.c CAT_Map2d.c CAT_SurfaceIO.c CAT_Surf.c diffeosphere.c optimizersphere.c
CAT_ApplyWarpValues_SOURCES = CAT_ApplyWarpValues.c CAT_Map2d.c CAT_SurfaceIO.c CAT_Surf.c diffeosphere.c optimizersphere.c
CAT_DumpSurfArea_SOURCES = CAT_DumpSurfArea.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_AverageSurfaces_SOURCES = CAT_AverageSurfaces.c CAT_SurfaceIO.c
CAT_DumpSurfaceRatio_SOURCES = CAT_DumpSurfaceRatio.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_Color2Surf_SOURCES = CAT_Color2Surf.c CAT_SurfaceIO.c
CAT_BlurSurfHK_SOURCES = CAT_BlurSurfHK.c CAT_Blur2d.c CAT_SurfaceIO.c
CAT_GLM_Estimate_SOURCES = CAT_GLM_Estimate.c CAT_Pinv.c CAT_SliceIO.c
CAT_GLM_Math_SOURCES = CAT_GLM_Math.c
CAT_SurfConvert_SOURCES = CAT_SurfConvert.c CAT_SurfaceIO.c 
CAT_FreesurferCurv2Txt_SOURCES = CAT_FreesurferCurv2Txt.c CAT_SurfaceIO.c 
CAT_PlotValuesAtPoint_SOURCES = CAT_PlotValuesAtPoint.c CAT_SurfaceIO.c
CAT_PlotValuesAtMaximum_SOURCES = CAT_PlotValuesAtMaximum.c CAT_SurfaceIO.c
CAT_SeparateCluster_SOURCES = CAT_SeparateCluster.c CAT_Blur2d.c CAT_SurfaceIO.c
CAT_Surf2Sphere_SOURCES = CAT_Surf2Sphere.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_ResampleSphericalSurf_SOURCES = CAT_ResampleSphericalSurf.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_AreaDistortion_SOURCES = CAT_AreaDistortion.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_AngularDistortion_SOURCES = CAT_AngularDistortion.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_FindFlats_SOURCES = CAT_FindFlats.c CAT_Curvature.c CAT_Blur2d.c
CAT_Isometize_SOURCES = CAT_Isometize.c CAT_Surf.c CAT_Blur2d.c CAT_Curvature.c CAT_SurfaceIO.c CAT_Map2d.c
CAT_WarpSurf_SOURCES = CAT_WarpSurf.c CAT_Map2d.c CAT_Blur2d.c CAT_SheetIO.c CAT_Curvature.c CAT_SurfaceIO.c CAT_Surf.c diffeosphere.c optimizersphere.c
CAT_Surf2SPH_SOURCES = CAT_Surf2SPH.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c CAT_SPH.c $(s2kit10_SRC)
CAT_SPH2Surf_SOURCES = CAT_SPH2Surf.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c CAT_SPH.c $(s2kit10_SRC)
CAT_ResampleSphericalSurfSPH_SOURCES = CAT_ResampleSphericalSurfSPH.c CAT_Surf.c CAT_SurfaceIO.c CAT_Map2d.c CAT_SPH.c $(s2kit10_SRC)
CAT_LBConformalMap_SOURCES = CAT_LBConformalMap.cpp CAT_Surf.c CAT_Blur2d.c CAT_Curvature.c CAT_SurfaceIO.c CAT_Map2d.c \
	nl/nl_api.c nl/nl_blas.c nl/nl_context.c nl/nl_iterative_solvers.c nl/nl_matrix.c nl/nl_os.c nl/nl_preconditioners.c nl/nl_superlu.c

m4_files = m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4

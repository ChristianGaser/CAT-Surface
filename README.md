# CAT: Cortex Analysis Tools
Christian Gaser (christian.gaser@uni-jena.de), Jena University Hospital, Germany.

These tools are used within [CAT12](https://github.com/ChristianGaser/cat12) for almost all of the surfaced-based functions.

## The following libraries are necessary to compile CAT:
- [netcdf-4.2](https://github.com/Unidata/netcdf-c)
- [minc-1.5.1](https://github.com/BIC-MNI/minc)
- [bicp-1.4.6](https://github.com/BIC-MNI/minc)
- [fftw-3.3.2](https://github.com/FFTW/fftw3)
- [expat-2.0](https://github.com/libexpat/libexpat)

## CAT_3dVol2Surf
Map values form 3D volume to surface

## CAT_ApplyWarpSurf
Apply deformations obtained with CAT_ApplySurf to surface

## CAT_ApplyWarpValues
Apply deformations obtained with CAT_ApplySurf to values form surface

## CAT_BlurSurfHK
Heat kernel smoothing of either surface or surface values

## CAT_Color2Surf
Colorize surface using surface value files

## CAT_DeformSurf
Extract cortical surface from 3D volume (very old algorithm from David MacDonald)

## CAT_DumpCurv
Extract different folding parameters and optionally post-process and/or smooth
values using diffusion heat kernel

## CAT_DumpCurv_ui
User interface for CAT_DumpCurv for dealing with many files

## CAT_DumpGI
Calculate gyrification index using the ratio of local surface area and local
inflated surface area. Local surface area can be approximated by use of different
curve types (default: mean curvature averaged over 3mm).

## CAT_DumpGI_ui
User interface for CAT_DumpGI for dealing with many files.

## CAT_DumpSurfArea
Extract surface area

## CAT_SurfConvert
Convert BIC obj to Freesurfer or OOGL and vice versa.

## CAT_FreesurferCurv2Txt
Convert freesurfer curv to txt-files

## CAT_GLM_Estimate
WARNING: error in estimating designs with more than one sample.
Estimate GLM parameters to analyze surface values or minc 3D volumes.

## CAT_GLM_Math
Perform simple math operations to value files (from CAT_GLM_Estimate)

## CAT_LBConformalMap
Laplace-Beltrami Conformal Mapping using the approach of Haker et al. "Conformal 
Surface Parameterization for Texture Mapping" 2000

## CAT_Obj2Neurolens
Display BIC obj-files with Neurolens

## CAT_PlotValuesAtMaximum
Plot values at maximum of a reference for each file

## CAT_PlotValuesAtPoint
Plot values at coordinate x y z for each file

## CAT_ResampleSphericalSurf
Resamples a spherical inflated surface to a sphere.

## CAT_ResampleSphericalSurf
Resamples a spherical inflated surface to an external given sphere.

## CAT_Sheet2Surf
Maps 2d image in pgm format to surface

## CAT_Surf2Sheet
Maps surface values or mean curvature to 2d sheet in pgm format

## CAT_Surf2Sphere
Maps a surface to a sphere by using the caret inflating approach. 

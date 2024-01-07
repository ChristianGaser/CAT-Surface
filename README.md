# CAT: Cortex Analysis Tools
Christian Gaser christian.gaser@uni-jena.de Jena University Hospital, Germany.

These tools are used within [CAT12](https://github.com/ChristianGaser/cat12) for almost all of the surfaced-based functions.

## The following libraries are necessary to compile CAT:
- [expat-2.0](https://github.com/libexpat/libexpat)
- [fftw-3](https://www.fftw.org)

## The following library was largely modified:
- [bicp-1.4.6](https://github.com/BIC-MNI/minc)
All dependencies in bicpl on minc and volume_io were removed and replaced by 
nifti-functions. The following functions from minc were added:
- ParseArgv.c
- alloc.c
- alloc_check.c
- progress.c
- arrays.c

## CAT_3dVol2Surf
Map values form 3D volume to surface

## CAT_ApplyWarpSurf
Apply deformations obtained with CAT_ApplySurf to surface

## CAT_ApplyWarpValues
Apply deformations obtained with CAT_ApplySurf to values form surface

## CAT_BlurSurfHK
Heat kernel smoothing of either surface or surface values

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

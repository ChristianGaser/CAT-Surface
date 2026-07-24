# Nipype interfaces for `cat_surf`

Reference Nipype interfaces wrapping the `cat_surf` / `t1prep.cat_surf` Python API.

## Scope

The upstream Nipype package `nipype.interfaces.t1prep` ships only the **subset of
interfaces needed to replace FreeSurfer / FSL (bbreg) / ANTs** in an
fMRIPrep-style pipeline:

| Interface | Replaces | Wraps |
| --- | --- | --- |
| `CatSurfVolSanlm` | ad-hoc denoising | `cat_surf.cli.vol_sanlm` |
| `CatSurfBbreg` | FSL `bbregister` / `flirt -bbr` | `cat_surf.bbreg` |
| `CatSurfBbregDetectContrast` | contrast heuristic | `cat_surf.bbreg_detect_contrast` |
| `CatSurfVolumeRegisterNmi` | ANTs affine (cross-modal, NMI) | `cat_surf.volume_register_nmi` |
| `CatSurfVolumeRegisterRobust` | ANTs affine (same-modal, robust) | `cat_surf.volume_register_robust` |

Plus the T1Prep command interfaces (`T1Prep`, `T1PrepSegment`,
`T1PrepSurfaceEstimation`, `T1PrepRealignLongitudinal`) that run
`python -m t1prep.<module>`.

## This folder

`cat_surf_interfaces_full.py` is an **archived reference** holding the *complete*
set of per-function interfaces (surface I/O, geometry, smoothing, resampling,
DARTEL / Spherical Demons registration, all volume ops). These were moved out of
the wired Nipype module to keep its PR surface and test matrix small, but are
kept here so any of them can be promoted back when a workflow needs it.

- Imports use absolute `nipype.*` paths and `_import_cat_surf` is inlined, so a
  class can be copied straight back into a live package.
- The registration/geometry interfaces demonstrating the recommended pattern
  (`TraitedSpec` + `SimpleInterface` + file-based I/O, e.g. `CatSurfWarp`,
  `CatSurfSphericalDemon`, `CatSurfGetArea`) should be the template when
  promoting an older `traits.Any` array-passing interface.

This file is **not** imported by the package and is not part of the build or the
Nipype test suite.

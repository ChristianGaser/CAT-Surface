# Third-Party Libraries

This directory contains libraries that are not maintained by the CAT-Surface
maintainers. These libraries are bundled for convenience so that no extra
downloads are needed during compilation.

## Components and Licenses

| Component | Description | License | Source |
|-----------|-------------|---------|--------|
| **bicpl-surface** | BIC Programming Library for surface mesh handling | MIT-like (MNI) | [BIC-MNI/bicpl](https://github.com/BIC-MNI/bicpl) |
| **volume_io** | Volume I/O and geometry utilities from MINC | MIT-like (MNI) | [BIC-MNI/minc](https://github.com/BIC-MNI/minc) |
| **nifti** | NIfTI-1 I/O library | Public Domain | [NIFTI](http://nifti.nimh.nih.gov/) |
| **gifticlib** | GIFTI I/O library | Public Domain | [NITRC](http://www.nitrc.org/projects/gifti) |
| **zlib** | Compression library | zlib License | [zlib.net](https://zlib.net/) |
| **expat** | XML parsing library | MIT License | [libexpat](https://libexpat.github.io/) |
| **nii2mesh** | Mesh simplification (quadric decimation) | BSD-2-Clause | [neurolabusc/nii2mesh](https://github.com/neurolabusc/nii2mesh) |
| **s2kit10** | Spherical harmonics transform | GPL-2+ | Kostelec & Rockmore |
| **dartel** | Diffeomorphic registration | GPL-2+ (SPM) | John Ashburner |
| **genus0** | Topology correction | Research use | Steve Haker (BWH/Harvard) |
| **MeshFix-V2.1** | Mesh repair (self-intersections) | **GPL-3** or Commercial | [MarcoAttene/MeshFix-V2.1](https://github.com/MarcoAttene/MeshFix-V2.1) |

## License Notes

### Permissive Licenses (bicpl-surface, volume_io, nifti, gifticlib, zlib, expat, nii2mesh)
These components have permissive licenses (MIT, BSD, public domain, zlib) that
allow use in both open-source and commercial applications without restrictions.

### GPL Components (s2kit10, dartel, MeshFix)
- **s2kit10** and **dartel**: GPL-2+, compatible with CAT-Surface GPL-3
- **MeshFix**: GPL-3 or requires commercial license from IMATI-GE/CNR

### Research-Use Components (genus0)
The genus0 library is distributed for research purposes. Contact the original
authors (haker@bwh.harvard.edu) for commercial licensing.

## Modifications

The following libraries have been modified from their original versions:

- **bicpl-surface**: Based on bicpl-1.4.6, with all volume function dependencies
  removed (using NIfTI functions instead). ParseArgv.c added from minc-1.5.1.
- **volume_io**: From minc-1.5.1, geometry and utility functions only.

## References

- MeshFix: Attene, M. "A lightweight approach to repairing digitized polygon
  meshes." The Visual Computer, 2010. DOI: 10.1007/s00371-010-0416-3
- S2kit: Kostelec & Rockmore, "FFTs on the rotation group", J Fourier Anal Appl 
  14, 145–179, 2008. DOI: 10.1007/s00041-008-9013-5

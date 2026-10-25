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
| **gifticlib** | GIFTI I/O library (version 1.17, 2021) | Public Domain | [gifti_clib](https://github.com/NIFTI-Imaging/gifti_clib) |
| **fftw** | Library for computing the discrete Fourier transform (DFT)  | GPL/FFTW License | [FFTW](https://fftw.org/) |
| **zlib** | Compression library (version 1.3.2, 2026) | zlib License | [zlib.net](https://zlib.net/) |
| **expat** | XML parsing library (version 2.9.0, 2026) | MIT License | [libexpat](https://libexpat.github.io/) |
| **nii2mesh** | Mesh simplification (quadric decimation) | BSD-2-Clause | [neurolabusc/nii2mesh](https://github.com/neurolabusc/nii2mesh) |
| **s2kit10** | Spherical harmonics transform | GPL-2+ | Kostelec & Rockmore |
| **dartel** | Diffeomorphic registration (no longer bundled, see below) | GPL-2+ (SPM) | John Ashburner |
| **genus0** | Topology correction | No license grant in the source (see below) | Steve Haker (BWH/Harvard) |

## License Notes

### Permissive Licenses (bicpl-surface, volume_io, nifti, gifticlib, zlib, expat, nii2mesh)
These components have permissive licenses (MIT, BSD, public domain, zlib) that
allow use in both open-source and commercial applications. Most of them still
require that their notice is kept with every copy: see `zlib/LICENSE`,
`expat/COPYING`, `volume_io/LICENSE` (which also covers the bundled netCDF and
ParseArgv headers), `bicpl-surface/COPYING`, `nifti/LICENSE`,
`gifticlib/LICENSE.gifti` and `nii2mesh/LICENSE`.

### GPL Components (fftw, s2kit10, dartel)
- **fftw**: GPL-2+ (see `fftw/COPYING`). It is built as a static library and
  linked into `libCAT` and into every program in `Progs/`.
- **s2kit10**: GPL-2+, compiled into `libCAT`.
- **dartel**: no longer bundled; the DARTEL code was moved to `../deprecated/`.
  A small part of `diffeopoly.c` ((c) John Ashburner 2007, from SPM) remains in
  `../Lib/CAT_WarpDemons.c` (sphere sample points of Spherical Demons).

Because FFTW and s2kit10 are part of every CAT-Surface binary, the binaries,
`libCAT` and the `cat-surf` Python package can only be distributed under the
GPL. The commercial option in `../LICENSE` covers the CAT-Surface original
source code; a commercial binary additionally needs a commercial FFTW license
(available from MIT), a replacement for s2kit10 and a check of the code
derived from DARTEL.

### genus0 (no license grant)
The genus0 sources (`genus0.c`, `genus0.h`) name the authors (Steve Haker and
Andre Robatino, Surgical Planning Lab, Brigham and Women's Hospital) but contain
no license statement. Without a license grant from the authors, genus0 cannot
be relicensed under the GPL, and its inclusion in the GPL distribution and in
the commercial option is unclear. A written permission or license from the
authors (haker@bwh.harvard.edu) should be obtained and stored here as
`genus0/LICENSE`.

## Security

expat parses every GIfTI file, and zlib inflates and deflates the
gzip-base64 GIfTI arrays and `.nii.gz` files, so both read data that can come
from anywhere. They were updated from 2.1.0 (2012) and 1.2.3 (2005), which
carried known memory-corruption bugs on these paths (e.g. expat
CVE-2016-0718, zlib CVE-2018-25032), to expat 2.9.0 and zlib 1.3.2; keep them
current. Two classes do not apply here: gifticlib installs an expat default
handler, which turns internal entity expansion ("billion laughs") off, and
nothing calls zlib's `inflateGetHeader()` (CVE-2022-37434). nifti1_io (1.43,
2010) has nothing security-relevant in its upstream history since; note that
nifti1_io looks for an image file by trying `.nii` before `.nii.gz`, which
CAT's readers work around (`read_nifti_float/double` read the voxels from the
file that was asked for).

## Modifications

The following libraries have been modified from their original versions:

- **bicpl-surface**: Based on bicpl-1.4.6, with all volume function dependencies
  removed (using NIfTI functions instead). ParseArgv.c added from minc-1.5.1.
- **volume_io**: From minc-1.5.1, geometry and utility functions only.
- **dartel**: Added more OpenMP support (now in `../deprecated/`).
- **zlib**: 1.3.2 as released, except that `zconf.h` includes `unistd.h` on
  every non-Windows system, as zlib's own `./configure` would set it.
- **expat**: the `lib/` sources of 2.9.0 as released, with a hand-written
  `expat_config.h` in place of the one expat's configure generates (hash salt
  from `/dev/urandom` on POSIX systems, `rand_s()` on Windows); the unused
  `xmlwf` tool is not vendored.
- **gifticlib**: 1.17 (NIFTI-Imaging/gifti_clib at d3e873d) as released,
  except that `gifti_io.h` includes `nifti1_io.h`; upstream defaults to
  `nifti2_io.h`, which is not vendored. The test and tool programs are not
  vendored.

## References

- S2kit: Kostelec & Rockmore, "FFTs on the rotation group", J Fourier Anal Appl 
  14, 145–179, 2008. DOI: 10.1007/s00041-008-9013-5

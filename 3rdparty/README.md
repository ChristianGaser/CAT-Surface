This directory contains libraries that are not maintained by the CAT-Surface
maintainers. These libraries are necessary to compile CAT_Surface and are 
provided here as courtesy to the end-users.

The following library was largely modified:
- [bicp-1.4.6](https://github.com/BIC-MNI/minc)

All dependencies in bicpl on minc and volume_io were removed and replaced by 
nifti-functions. The following functions from minc were added:
- ParseArgv.c
- alloc.c
- alloc_check.c
- progress.c
- arrays.c

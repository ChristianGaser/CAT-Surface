This directory contains libraries that are not maintained by the CAT-Surface
maintainers. These libraries are needed to compile CAT_Surface and are provided 
here as a courtesy to end users.

The following libraries have undergone major changes and all dependencies on 
volume functions have been removed, since we are constantly using nifti functions
only:
- [bicp-1.4.6](https://github.com/BIC-MNI/minc)
- [minc-1.5.1/volume_io](https://github.com/BIC-MNI/minc)

The following function was added from minc-1.5.1:
- ParseArgv.c

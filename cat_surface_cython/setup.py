from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import sys
import numpy as np

CAT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_BUILD = os.environ.get("CAT_BUILD_DIR", os.path.join(CAT_ROOT, "build-native-arm64"))

include_dirs = [
    os.path.join(CAT_ROOT, "3rdparty", "volume_io", "Include"),
    os.path.join(CAT_ROOT, "3rdparty", "bicpl-surface", "Include"),
    np.get_include(),
]

library_dirs = [
    os.path.join(DEFAULT_BUILD, ".libs"),
]

extra_objects = [
    os.path.join(DEFAULT_BUILD, ".libs", "libCAT.a"),
]

extra_compile_args = ["-std=c99"]
extra_link_args = []

extensions = [
    Extension(
        name="cat_surface.surface_io",
        sources=[os.path.join(CAT_ROOT, "cat_surface_cython", "cat_surface", "surface_io.pyx")],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        extra_objects=extra_objects,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language="c",
    )
]

class BuildExt(build_ext):
    def build_extensions(self):
        try:
            from Cython.Build import cythonize
        except Exception:
            print("Cython is required to build this extension. pip install cython", file=sys.stderr)
            raise
        self.extensions = cythonize(self.extensions, compiler_directives={"language_level": 3})
        super().build_extensions()

setup(
    name="cat-surface-cython",
    version="0.1.0",
    description="Cython-based BICPL surface IO for CAT-Surface",
    packages=["cat_surface"],
    package_dir={"cat_surface": os.path.join("cat_surface_cython", "cat_surface")},
    ext_modules=extensions,
    cmdclass={"build_ext": BuildExt},
    python_requires=">=3.8",
    install_requires=["numpy"],
)

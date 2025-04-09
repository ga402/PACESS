from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("distance3d", ["distance3d.pyx"], include_dirs=[np.get_include()]),
    Extension("bikernel", ["bikernel.pyx"], include_dirs=[np.get_include()])
]

setup(
    name="cython_geometry_kernels",
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    include_dirs=[np.get_include()]
)

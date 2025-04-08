from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name="get_distances_utils",
    ext_modules=cythonize("get_distances_utils.pyx", compiler_directives={'language_level': "3"}),
    include_dirs=[numpy.get_include()]
)

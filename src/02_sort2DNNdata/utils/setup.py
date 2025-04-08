from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy


extensions = [
    Extension("image_loader", ["image_loader.pyx"], include_dirs=[numpy.get_include()]),
    Extension("subsample_mean", ["subsample_mean.pyx"], include_dirs=[numpy.get_include()])
]



setup(
    ext_modules=cythonize(extensions, language_level="3"),
    #include_dirs=[numpy.get_include()],
)
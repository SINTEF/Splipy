import os
import sys

from setuptools import Extension
from Cython.Build import cythonize
import numpy as np

extensions = cythonize(
    Extension(
        "splipy.basis_eval",
        ["splipy/basis_eval.pyx"],
        include_dirs=[np.get_include()],
    )
)

def build(setup_kwargs):
    setup_kwargs.update({
        "ext_modules": extensions,
    })

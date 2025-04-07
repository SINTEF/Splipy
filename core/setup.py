from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

extensions = cythonize(
    Extension(
        "splipy_core._core",
        ["src/splipy_core/core.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    )
)

setup(
    ext_modules=extensions,
)

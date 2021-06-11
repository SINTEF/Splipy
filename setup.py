#!/usr/bin/env python3

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np

with open('PyPI_text.md') as f:
    long_description = f.read()

setup(
    name='Splipy',
    version='1.5.3',
    description='Spline modelling library for Python',
    long_description_content_type='text/markdown',
    long_description=long_description,
    keywords=['Bspline', 'Splines', 'NURBS', 'Curve', 'Surface', 'Volume', 'Interpolation', 'Approximation', 'Fit', 'Integration', 'Differentiation'],
    url='https://github.com/sintefmath/Splipy',
    maintainer='Kjetil Andre Johannessen',
    maintainer_email='kjetijo@gmail.com',
    license='GNU public license v3',
    packages=find_packages(),
    package_data={
        'splipy': ['templates/*.bpt'],
    },
    install_requires=[
        'numpy    >= 1.15',
        'scipy    >= 1.2',
    ],
    extras_require={
        'FiniteElement': ["nutils>=4.0"],
        'Images':        ["opencv-python>=4.0"],
        'Rhino':         ["rhino3dm>=0.14"],
    },
    # ext_modules=cythonize("splipy/basis_eval.pyx"),
    ext_modules=cythonize([
        Extension(
            'splipy.basis_eval',
            ['splipy/basis_eval.pyx'],
            include_dirs=[np.get_include()],
        )
    ]),
    # include_dirs=[np.get_include()],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Multimedia :: Graphics :: 3D Modeling',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
    ],
)

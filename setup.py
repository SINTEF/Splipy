#!/usr/bin/env python

from setuptools import setup
from splipy import __version__

setup(
    name='Splipy',
    version='1.2.0',
    description='Spline modelling library for Python',
    long_description='This repository contains the Splipy packages. Splipy is a pure python library for the creation, evaluation and manipulation of B-spline and NURBS geometries. It supports n-variate splines of any dimension, but emphasis is made on the use of curves, surfaces and volumes. The library is designed primarily for analysis use, and therefore allows fine-grained control over many aspects which is not possible to achieve with conventional CAD tools.',
    keywords=['B-spline', 'Splines', 'NURBS', 'Curve', 'Surface', 'Volume', 'Interpolation', 'Approximation', 'Fit', 'Integration', 'Differentiation'],
    url='https://github.com/sintefmath/Splipy',
    maintainer='Kjetil Andre Johannessen',
    maintainer_email='kjetijo@gmail.com',
    license='GNU public license v3',
    packages=['splipy', 'splipy.utils', 'splipy.io'],
    package_data={
        'splipy': ['templates/*.bpt'],
    },
    install_requires=[
        'numpy >= 1.9',
        'scipy >= 0.17',
    ],
    extra_requires={
        'FiniteElement': ["nutils>=2.0"],
        'Images':        ["opencv-python>=3.3"],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Multimedia :: Graphics :: 3D Modeling',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
)

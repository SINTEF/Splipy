#!/usr/bin/env python

from setuptools import setup
from splipy import __version__

setup(
    name='Splipy',
    version='1.1.0',
    description='Spline modelling library for Python',
    url='https://github.com/sintefmath/Splipy',
    maintainer='Kjetil Andre Johannessen',
    maintainer_email='kjetijo@gmail.com',
    license='GNU public license v3',
    packages=['splipy', 'splipy.utils', 'splipy.IO'],
    package_data={
        'splipy': ['templates/*.bpt'],
    },
    install_requires=[
        'numpy >= 1.9',
        'scipy >= 0.17',
    ],
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

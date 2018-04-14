#!/usr/bin/env python3

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
# import numpy as np

class CustomBuildExtCommand(build_ext):

    def run(self):
        import numpy as np
        self.include_dirs.append(np.get_include())
        build_ext.run(self)

with open('PyPI_text.md') as f:
    long_description = f.read()

setup(
    name='Splipy',
    version='1.3.0',
    description='Spline modelling library for Python',
    long_description_content_type='text/markdown',
    long_description=long_description,
    keywords=['Bspline', 'Splines', 'NURBS', 'Curve', 'Surface', 'Volume', 'Interpolation', 'Approximation', 'Fit', 'Integration', 'Differentiation'],
    url='https://github.com/sintefmath/Splipy',
    maintainer='Kjetil Andre Johannessen',
    maintainer_email='kjetijo@gmail.com',
    license='GNU public license v3',
    packages=['splipy', 'splipy.utils', 'splipy.io'],
    package_data={
        'splipy': ['templates/*.bpt'],
    },
    setup_requires=['numpy >= 1.9'],
    install_requires=[
        'numpy >= 1.9',
        'scipy >= 0.17',
    ],
    extras_require={
        'FiniteElement': ["nutils>=3.0"],
        'Images':        ["opencv-python>=3.3"],
    },
    cmdclass={'build_ext': CustomBuildExtCommand},
    ext_modules=[Extension('splipy.basis_eval', ['splipy/basis_eval.c'])],
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

#!/usr/bin/env python

from setuptools import setup

# Creates the __version__ name. We can't import it because it will try to load
# dependencies before they are installed.
with open('splipy/__version__.py', 'r') as f:
    exec(f.read())

setup(
    name='Splipy',
    version='1.0.0',
    description='Spline modelling library for Python',
    maintainer='Arne Morten Kvarving',
    maintainer_email='arne.morten.kvarving@sintef.no',
    packages=['splipy', 'splipy.utils', 'splipy.IO'],
    package_data={
        'splipy': ['templates/*.bpt'],
    },
    install_requires=[
        'numpy >= 1.9',
        'scipy >= 0.17',
    ],
)

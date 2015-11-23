# GeoModeller

This repository contains the GeoModeller module for the IFEM project.

1. To compile this, you will need to add the IFEM repository, at https://launchpad.net/~ifem/

Remember to run `sudo apt-get update` after adding the repository and continuing.

2. Install the cmakerules (plan is to discontinue these at some point) from https://github.com/akva2/IFEM-cmakerules

3. Followed by installing these packages

    sudo apt-get install cmake python-dev libnewmat10ldbl libgotools-compositemodel-dev libgotools-core-dev libgotools-igeslib-dev libgotools-implicitization-dev libgotools-intersections-dev libgotools-isogeometricmodel-dev libgotools-qualitymodule-dev libgotools-qualitymodule1 libgotools-parametrization-dev libgotools-topology-dev libgotools-trivariate-dev libgotools-trivariatemodel-dev libopennurbs1-dev libttl-dev libsisl-dev zlib1g-dev python-numpy python-scipy python-yaml python-opencv python-epydoc


### Compiling the code

To compile type the following in the commandline

    cmake -DCMAKE_BUILD_TYPE=Release .
    make


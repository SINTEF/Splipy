===========
GeoModeller
===========

This repository contains the GeoModeller module for the IFEM project.

1. To compile this, you will need to add the IFEM repository, at https://launchpad.net/~ifem/
   Remember to run `sudo apt-get update` after adding the repository and continuing.

2. Followed by installing these packages::

    sudo apt-get install ifem-geomodeler-builddeps

Compiling the code
------------------

To compile type the following in the command-line::

    cmake -DCMAKE_BUILD_TYPE=Release .
    make

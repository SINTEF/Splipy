# Examples on how to write geometries to file
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian Univeristy of Science and Technology (NTNU)
# Date:      March 2016
#

from sys import path
path.append('../')
from splipy.io import *
import splipy.surface_factory as surfaces

# create a NURBS torus
torus = surfaces.torus(minor_r=1, major_r=4)


# STL files are tesselated linear triangles. View with i.e. meshlab
with STL('torus.stl') as my_file:
    my_file.write(torus, n=(50, 150)) # specify resolution of 50x150 evaluation pts


# G2 files are native GoTools (http://www.sintef.no/projectweb/geometry-toolkits/gotools/)
with G2('torus.g2') as my_file:
    my_file.write(torus)

# Examples on how to write geometries to file
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian University of Science and Technology (NTNU)
# Date:      March 2016
#

from pathlib import Path
from splipy.io import G2, STL
import splipy.surface_factory as surfaces

# create a NURBS torus
torus = surfaces.torus(minor_r=1, major_r=4)

# Path to file directory
path = Path(__file__).parent

# STL files are tessellated linear triangles. View with i.e. meshlab
with STL(str(path / "torus.stl")) as my_file:
    my_file.write(torus, n=(50, 150)) # specify resolution of 50x150 evaluation pts


# G2 files are native GoTools (http://www.sintef.no/projectweb/geometry-toolkits/gotools/)
with G2(str(path / "torus.g2")) as my_file:
    my_file.write(torus)

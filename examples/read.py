# Examples on how to write geometries to file
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian University of Science and Technology (NTNU)
# Date:      October 2017
#

from sys import path
path.append('../')
from splipy.io import *

# G2 files are native GoTools (http://www.sintef.no/projectweb/geometry-toolkits/gotools/)

# Read a single NURBS patch from the file 'sphere.g2'
with G2('sphere.g2') as my_file:
    my_sphere = my_file.read()

# Read multiple NURBS patches from the file 'teapot.g2'
with G2('teapot.g2') as my_file:
    my_teapot = my_file.read()

print(type(my_teapot)) # <class 'list'>
print(len(my_teapot))  # 32
print(type(my_sphere)) # <class 'list'>
print(len(my_sphere))  # 1

# dump knot vectors and controlpoints of all 33 Surface patches to screen
print(my_sphere[0])
for surf in my_teapot:
    print(surf)

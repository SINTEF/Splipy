# Examples on how to create a sweeping geometry
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian University of Science and Technology (NTNU)
# Date:      June 2017
#

from sys import path
path.append('../')

from splipy import *
from splipy.io import *
from numpy import pi, cos, sin
import numpy as np
import splipy.curve_factory   as curves
import splipy.surface_factory as surfaces


### create a parametric 3D-function
def trefoil(t):
    t = np.array(t)
    return np.array([sin(t) + 2*sin(2*t), cos(t) - 2*cos(2*t), -sin(3*t)]).T

### do an adaptive best cubic spline-fit of this function
path  = curves.fit(trefoil, 0, 2*pi)

### since we know it is a closed curve, enforce this on the path
path  = path.make_periodic(0,0)

### create a sweeping curve (either a circle or square)
shape = curves.circle(r=0.2)
# shape = 16*curves.n_gon(4)

### sweep *shape along *path
srf = surfaces.sweep(path, shape)

### write results to file. Use meshlab (www.meshlab.net) to view stl-files
with STL('trefoil.stl') as f:
    f.write(srf, n=(150,30))


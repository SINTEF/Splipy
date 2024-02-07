# Reuleaux triangles, or rollers, are objects of constant width that are not
# circles. More info: https://en.wikipedia.org/wiki/Reuleaux_triangle
#
# Author:    Kjetil Andre Johannessen
# Institute: Norwegian University of Science and Technology (NTNU)
# Date:      March 2016
#

from sys import path
path.append('../')
from splipy import *
import splipy.curve_factory as curves
import splipy.surface_factory as surfaces
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pi, cos, sin

# create the three sides of the triangle, each consisting of a circle segment
c1 = curves.circle_segment(pi/3)
c2 = c1.clone().rotate(2*pi/3) + [1,0]
c3 = c1.clone().rotate(4*pi/3) + [cos(pi/3), sin(pi/3)]

# merge the three circles into one, and center it at the origin
c = c1.append(c2).append(c3)
c -= c.center()

# plot the reuleaux triangle
t = np.linspace(c.start(0), c.end(0), 151) # 151 parametric evaluation points
x = c(t)                                   # evaluate (x,y)-coordinates
plt.plot(x[:,0], x[:,1])
plt.axis('equal')
plt.show()

# split the triangle in two, and align this with the y-axis
two_parts = c.split((c.start(0) + c.end(0)) / 2.0)
half_curve = two_parts[0].rotate(2*pi/3)

# create a surface by revolving around its symmetry axis (y-axis)
surf = surfaces.revolve(half_curve, axis=[0,1,0])

# plot the resulting surface on a uniform grid of 71x71 evaluation points
u = np.linspace(surf.start(0), surf.end(0), 71)
v = np.linspace(surf.start(1), surf.end(1), 71)
x = surf(u,v)
ax = plt.axes(projection='3d')
ax.plot_wireframe(x[:,:,0], x[:,:,1], x[:,:,2])
plt.show()

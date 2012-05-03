from math import *
from GoTools import * 
from GoTools.CurveFactory import * 
from GoTools.SurfaceFactory import * 
from GoTools.VolumeFactory import * 
from GoTools.SurfaceModelFactory import * 

SetFinalOutput(file="trimmed.g2")

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

rect     = Rectangle([-2,-2,0], # corner
                     x_axis,    # u-axis
                     y_axis,    # v-axis
                     4, 4)      # length

cylinder = CylinderSurface(-z_axis, # center
                           z_axis,  # axis
                           1,       # radius
                           2)       # height


trimmed_surface = TrimSurface(rect, cylinder)
model = RegularizeSurface(trimmed_surface)

# print model
s1 =  model[0]
s2 =  model[1]
s3 =  model[2]
s4 =  model[3]

FinalOutput([s1,s2,s3,s4], True)

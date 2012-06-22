from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from math import *

# global setup parameters
SetDimension(3)
SetFinalOutput("coons.g2")

# model parameters
N    = 10
xMax = 2*pi-0.9
yMax = 5.0

# First we create the point clouds to interpolate
p_bottom = list()
p_top = list()
p_left = list()
p_right = list()

for i in range(0,N):
	x = i*xMax/(N-1)
	p_bottom.append(Point( x, sin(x) - x/xMax*sin(xMax), 0))

for i in range(0,N):
	x = i*xMax/(N-1)
	p_top.append(Point( 2*x, yMax, 0))

for i in range(0,N):
	y = i*yMax/(N-1)
	p_left.append(Point( y*(y-yMax)/yMax/yMax*4, y, 0))

for i in range(0,N):
	y = i*yMax/(N-1)
	p_right.append(Point( xMax + y*xMax/yMax, y, 0))

# interpolate the four boundary curves
bottom = ApproximateCurve(p_bottom, range(0,N))
top    = ApproximateCurve(p_top,    range(0,N))
right  = ApproximateCurve(p_right,  range(0,N))
left   = ApproximateCurve(p_left,   range(0,N))

# make coons patch for the geometry interior
surface = CoonsSurfacePatch([bottom, top, left, right])

# write output
FinalOutput([surface], True)

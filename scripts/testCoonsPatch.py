from GoTools import *
from GoTools.SurfaceFactory import *
from GeoUtils.Interpolate import *
from math import *

# global setup parameters
SetDimension(3)
SetFinalOutput("coons.g2")

# model parameters
N    = 10
xMax = 2*pi-0.9
yMax = 5.0

# First we create the point clouds to interpolate
p_bottom = [[],[]]
p_top    = [[],[]]
p_left   = [[],[]]
p_right  = [[],[]]

for i in range(0,N):
	x = i*xMax/(N-1)
	p_bottom[0].append( x )
	p_bottom[1].append(sin(x) - x/xMax*sin(xMax))

for i in range(0,N):
	x = i*xMax/(N-1)
	p_top[0].append( 2*x )
	p_top[1].append(yMax)

for i in range(0,N):
	y = i*yMax/(N-1)
	p_left[0].append( y*(y-yMax)/yMax/yMax*4 )
	p_left[1].append(y)

for i in range(0,N):
	y = i*yMax/(N-1)
	p_right[0].append( xMax + y*xMax/yMax )
	p_right[1].append(y)

# interpolate the four boundary curves
bottom = UniformCubicCurve(p_bottom[0], p_bottom[1])
top    = UniformCubicCurve(p_top[0]   , p_top[1]   )
right  = UniformCubicCurve(p_right[0] , p_right[1] )
left   = UniformCubicCurve(p_left[0]  , p_left[1]  )
top.FlipParametrization()
left.FlipParametrization()

# make coons patch for the geometry interior
surface = CoonsSurfacePatch([bottom, right, top, left])

# write output
FinalOutput([surface], True)

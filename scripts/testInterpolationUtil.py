from math import *
from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GeoUtils.Interpolate import *

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

# allocate arrays and number of interpolation points N
N = 10;
xi  = []
xp  = []
yp  = []
crv = []

# create some dataset
for i in range(N):
	x = i*2*pi/(N-1)
	y = sin(x)
	xi.append(x)
	xp.append(x)
	yp.append(y)


# test the Interpolate entry points and optional parameters
crv.append(Linear(range(10)))
crv.append(Linear(range(10), range(10)))
crv.append(Linear(range(10), range(10), range(10)))
crv.append(Cubic(yp))
crv.append(Cubic(range(10), yp))
crv.append(Cubic(range(10), range(10), yp))

# test the Interpolate different entry points
crv.append(Linear(xp,yp))
crv.append(UniformCubic(xp, yp))

FinalOutput(crv)


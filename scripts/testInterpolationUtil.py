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
pts = []
xi  = []
xp  = []
yp  = []
crv = [None]*12

# create some dataset
for i in range(N):
	x = i*2*pi/(N-1)
	y = sin(x)
	pts.append(Point(x,y,0))
	xi.append(x)
	xp.append(x)
	yp.append(y)


# test the GoTools entry point for linear interpolation
crv[0]  = ApproximateCurve(pts, xi, 2)

# test the Interpolate entry points and optional parameters
crv[1] = Linear(range(10))
crv[2] = Linear(range(10), range(10))
crv[3] = Linear(range(10), range(10), range(10))
crv[4] = Cubic(yp)
crv[5] = Cubic(range(10), yp)
crv[6] = Cubic(range(10), range(10), yp)

# test the Interpolate different entry points
crv[7]  = Linear(xp,yp)
crv[8]  = LinearP(pts)
crv[9]  = CubicP( pts)
crv[10] = UniformCubic(xp, yp)
crv[11] = UniformCubicP( pts)

FinalOutput(crv)


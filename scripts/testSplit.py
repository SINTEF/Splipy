from math import *
from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GeoUtils.Interpolate import *
from GeoUtils.Elementary import *
from GeoUtils.Factory import *
from GeoUtils.Refinement import *

# standard nice things
o = Point(0,0,0)
x = Point(1,0,0)
y = Point(0,1,0)
z = Point(0,0,1)

# create a slightly more interesting volume than the unit cube
pts = []
for i in range(10):
	t = i*2.0*pi/9
	pts.append( Point(t, sin(t), 0))
crv1 = CubicP(pts)
crv2 = LineSegment(2*y, 5*y + 7*x)
surf = LoftBetween(crv1, crv2)
vol  = ExtrudeSurface(surf, z+y, 3.0)

# make it a somewhat more uniformly refined
vol.RaiseOrder(0,1,1)
UniformVolume(vol, 2, 3)
UniformVolume(vol, 3, 3)

# split in u-direction
vols   = SplitVolume(vol, 0, 0.3)
result = []
# split in v-direction
for v in vols:
	result += SplitVolume(v, 1, 0.3)
vols   = result[:]
result = []
# split in w-direction
for v in vols:
	result += SplitVolume(v, 2, 0.3)


# test the splitting of surfaces as well
surf2 = Translate(surf, 5*z)
surf2.RaiseOrder(0,2)
UniformSurface(surf2, 2, 3)
s       = SplitSurface(surf2, 0, 0.22)
result += SplitSurface(s[0],  1, 0.22)
result += SplitSurface(s[1],  1, 0.22)


# output file name
SetFinalOutput("split.g2")
FinalOutput( result , True)

 

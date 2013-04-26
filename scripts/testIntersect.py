from math import *
from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GeoUtils.Elementary import *
from GeoUtils.Refinement import *

# standard nice things
o = Point(0,0,0)
x = Point(1,0,0)
y = Point(0,1,0)
z = Point(0,0,1)

c1 = Circle(o, 1, z)
c2 = Circle(0.5*z, 0.8, z)
c3 = Circle(z, 0.7, z)
s  = Rectangle(-z-x-y*0.8, x, z, 2, 3)
s2 = CylinderSurface(-y, y, 0.1, 2.5)

s.RaiseOrder(1,1)
UniformCurve(c1, 2)
UniformCurve(c2, 2)
UniformCurve(c3, 2)
UniformSurface(s, 0, 3)
UniformSurface(s2, 0, 3)

result = []

if c1.Intersect(s):
	result.append(LineSegment(o, x))
if c2.Intersect(s):
	result.append(LineSegment(1*z, 1*z + x))
if c3.Intersect(s):
	result.append(LineSegment(2*z, 2*z + x))
if s.Intersect(s2):
	result.append(LineSegment(3*z, 3*z + x))
if s2.Intersect(s):
	result.append(LineSegment(4*z, 4*z + x))
if s.Intersect(c2):
	result.append(LineSegment(5*z, 5*z + x))
if s.Intersect(c1):
	result.append(LineSegment(6*z, 6*z + x))
if s.Intersect(c3):
	result.append(LineSegment(7*z, 7*z + x))

# output file name
SetFinalOutput("test.g2")
FinalOutput( result, True);

 

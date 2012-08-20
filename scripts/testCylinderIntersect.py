from math import *
from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

# make the unit cylinder around the z-axis which we will do intersections with
cyl       = CylinderSurface([0,0,-2],  z_axis , 1.0, 4)

# make 3 different objects creating different intersections
orthoCyl  = CylinderSurface([-2,0,0],  x_axis , 0.5, 4)
sphere    = SphereSurface([2,0,0], 1.0)
tiltedCyl = CylinderSurface(origin,  z_axis , 0.5, 7)

# fix the tilted cylinder
tiltedCyl.Rotate(y_axis, pi/4)
tiltedCyl.Translate(-(1.0-sqrt(2)/4)*x_axis)

# do the actual intersections
twoCurves, noPoints  = IntersectCylinder(orthoCyl,  cyl)
noCurves,  onePointA = IntersectCylinder(sphere,    cyl)
oneCurve,  onePointB = IntersectCylinder(tiltedCyl, cyl)

# accumulate the result curves
result = twoCurves
result.extend(oneCurve)
result.append(LineSegment(onePointA[0], onePointB[0]))

FinalOutput(result, True)


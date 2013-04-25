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

# output file name
SetFinalOutput("test.g2")
results = []

# create a rational curve and manipulate
rc1 = CircleSegment(o,x,2*pi/3,z)
rc2 = Translate(rc1, z)
rc3 = Rotate(   rc1, x, pi/2)
rc4 = Mirror(   rc2, o, x+y)
results.append(rc1)
results.append(rc2)
results.append(rc3)
results.append(rc4)

# create two non-rational and one rational surfaces
e1 = rc1.Rebuild(7,3) # to get a non-rational circle representation
e2 = LineSegment([cos(2*pi/3), sin(2*pi/3), 0], 0.4*y-x)
e3 = LineSegment(0.4*y-x, -2*y)
e4 = LineSegment(-2*y, x)
s1 = CoonsSurfacePatch([e1,e2,e3,e4])
s2 = Rectangle(x, -x-2*y, 2*x-y, sqrt(5), 3*sqrt(2))
s3 = CylinderSurface(x-z, 2*x-y, 1, 3) # should really have been a CylinderSEction, but that function is buggy

# build a surface model and manipulate this
sm1 = SurfaceModel(surfaces=[s1,s2,s3])
sm2 = Translate(sm1, 3*z)
sm3 = Rotate(   sm2, y, 2*pi/3)
sm4 = Mirror(   sm2, -3*y, x+2*y)
results.append(sm1)
results.append(sm2)
results.append(sm3)
results.append(sm4)

# create a few volumes and add the to a volume model
# vol1 = Sphere(x+y, 1.0)
# vol2 = Torus(x+y+2*z, z, 1.0, 0.5)
# vol3 = Parallelepiped(o, -x, -y-x, -x-y+z, 1, 1, 1)
vol1 = ExtrudeSurface(s1, -z, 3)
vol2 = Parallelepiped(-2*y, x+2*y, -x+2.4*y, y+z, 1, 1, 1)
vol3 = Cylinder(x-2*y-2*z, 1.0/2/sqrt(2)*x+(1.0/sqrt(2)-2)*y-2*z, z, 1)

vm1 = VolumeModel([vol1, vol2, vol3])
vm2 = Translate(vm1, 3*z)
vm3 = Rotate(   vm2, y, 2*pi/3)
#####        the mirror case actually fails... don't know why
# vm4 = Mirror(   vm1, 5*y, y) 
#################################################
results.append(vm1)
results.append(vm2)
results.append(vm3)


FinalOutput( results, True);

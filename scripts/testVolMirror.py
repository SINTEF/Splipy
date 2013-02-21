from math import *
from GoTools import *
from GeoUtils.Elementary import *
from GoTools.VolumeFactory import *

o = Point(0,0,0)
x = Point(1,0,0)
y = Point(0,1,0)
z = Point(0,0,1)


box  = Box(o, x, y, 1,1,1)
cone = Cone( [4,1,4], -z, pi/6, 4 )
ball = Sphere([1,5,5], 1.0)

mirBox  = Mirror(box,  o, x+y+z)
mirCone = Mirror(cone, o, x+y+z)
mirBall = Mirror(ball, o, x+y+z)
print 'BOX'
print box
print ''
print '**** Mirrored ****'
print ''
print mirBox
print 'CONE'
print ''
print '**** Mirrored ****'
print ''
print mirCone

FinalOutput([box, cone, ball, mirBox, mirCone, mirBall])

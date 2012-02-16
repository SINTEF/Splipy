from GoTools import *
from GoTools.CurveFactory import *

SetDimension(3)

# Infinite line
p0        = Point(1.2,1.2,1.2)
direction = Point(4,2,1)
infline = Line(p0,direction)

# Line segment
p0 = Point(1,1,1.1)
p1 = Point(2,0,1)
line = LineSegment(p0,p1)

# Circle
c = Point(1,0,0)
normal = Point(1,-1,0.5)
circle = Circle(c,1.5,normal)

# Circle segment
c = Point(1,0,0.2)
p0 = Point(2,1,0.2)
normal = Point(1,-1,0.5)
circleseg = CircleSegment(c,p0,1.1,normal)

# Ellipse
c = Point(-1,2,0)
axis = Point(3,1,0.4)
normal = Point(0.8,-2,-1)
ellipse = Ellipse(c,axis,3.0,1.4,normal)

# Elliptic segment
c = Point(-1,1.5,0)
axis = Point(3,1,0.4)
normal = Point(0.8,-2,-1)
ellipseg = EllipticSegment(c,axis,3.0,1.4,1.2,2.9,normal)

# Helix
c = Point(-1,3,2)
p0 = Point(0,1.1,2.5)
axis = Point(2,-1,0.7)
helix = Helix(c,p0,axis,0.6,15)

print infline
print line
print circle
print circleseg
print ellipse
print ellipseg
print helix

from GoTools import *
from GoTools.CurveFactory import *

SetDimension(3)

# Line segment
p0 = Point(1,1,1.1)
p1 = Point(2,0,1)
line = LineSegment(p0,p1)

l = line.Split([0.3,0.7])
l += line.Split(0.5)

FinalOutput(l,True)

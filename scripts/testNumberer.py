from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GeoUtils.IO import Numberer

SetDimension(3)

numberer = Numberer()
tol = 1e-6
output = []

# Edge -> Vertex
line = LineSegment(Point(0,0,0), Point(1,1,1))
getter = numberer.GetterFunction('edge', 'vertex')
assert(abs(getter(line, 0) - Point(0,0,0)) < tol)
assert(abs(getter(line, 1) - Point(1,1,1)) < tol)

# Face -> Vertex
face = Rectangle(Point(0,0,0), Point(1,0,0), Point(0,1,0), 1, 1)
getter = numberer.GetterFunction('face', 'vertex')
assert(abs(getter(face, 0) - Point(0,0,0)) < tol)
assert(abs(getter(face, 1) - Point(1,0,0)) < tol)
assert(abs(getter(face, 2) - Point(0,1,0)) < tol)
assert(abs(getter(face, 3) - Point(1,1,0)) < tol)

# Face -> Edge
getter = numberer.GetterFunction('face', 'edge')
for i in xrange(4):
    output.append(getter(face, i))

# Volume -> Vertex
volume = Box(Point(0,0,0), Point(1,0,0), Point(0,1,0), 1, 1, 1)
getter = numberer.GetterFunction('volume', 'vertex')
assert(abs(getter(volume, 0) - Point(0,0,0)) < tol)
assert(abs(getter(volume, 1) - Point(1,0,0)) < tol)
assert(abs(getter(volume, 2) - Point(0,1,0)) < tol)
assert(abs(getter(volume, 3) - Point(1,1,0)) < tol)
assert(abs(getter(volume, 4) - Point(0,0,1)) < tol)
assert(abs(getter(volume, 5) - Point(1,0,1)) < tol)
assert(abs(getter(volume, 6) - Point(0,1,1)) < tol)
assert(abs(getter(volume, 7) - Point(1,1,1)) < tol)

# Volume -> Edge
getter = numberer.GetterFunction('volume', 'edge')
for i in xrange(12):
    output.append(getter(volume, i))

# Volume -> Face
getter = numberer.GetterFunction('volume', 'face')
for i in xrange(6):
    output.append(getter(volume, i))

FinalOutput(output)

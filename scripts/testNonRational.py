from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

SetFinalOutput(file="RoundThings.g2")

crv  = Circle(      origin,     1,    z_axis)
surf = CircularDisc(z_axis,  [1,1,0], z_axis)
vol  = Sphere(      [0,0,3],    1)

# since convertion from Disc to NURBS is not implemented, we need to create
# surf in a different way:

surf = ContractTo(crv, origin)
surf.Translate(z_axis)

# convert them all to non-rational representations
c = NonRationalCurve(crv)
s = NonRationalSurface(surf)
v = NonRationalVolume(vol)

FinalOutput([c, s, v], True)

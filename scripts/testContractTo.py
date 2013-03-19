from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *

# standard nice things
origin = Point(0,0,0)
x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

# create unit disc by contracting unit circle to origin
crv  = Circle( origin, 1, z_axis)
surf = ContractCurveTo(crv, origin)

# create a unit ball by contracting unit shell to the origin
ball  = Sphere(origin, 1)
edges = ball.GetFaces()
shell = edges[1]

shell.Translate([0,0,2]) # move up for displaying purposes
vol   = ContractSurfaceTo(shell, [0,0,2])

FinalOutput([surf, vol], True)

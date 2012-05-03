
from GoTools import * 
from GoTools.SurfaceFactory import * 
from GoTools.VolumeFactory import * 


#   Make a cross of eight boxes, with the 
#   four center ones being degenerated to
#   extruded triangles (6 corner points)
#
#                    _____
#                   |     |
#                   |     |
#              _____|_____|_____
#             |     | \ / |     |
#             |     |  X  |     |
#             |_____|_/_\_|_____|
#                   |     |
#                   |     |
#                   |_____|
#            

SetFinalOutput("degen.g2")

x_axis = Point(1,0,0)
y_axis = Point(0,1,0)
z_axis = Point(0,0,1)

origin = Point(0,0,0)


Center = Rectangle(origin, x_axis, z_axis, 2, 2)
Center.Project("Z")
R1 = Rectangle([ 1,-1,0],  y_axis, z_axis, 2, 2)
R2 = Rectangle([ 1, 1,0], -x_axis, z_axis, 2, 2)
R3 = Rectangle([-1, 1,0], -y_axis, z_axis, 2, 2)
R4 = Rectangle([-1,-1,0],  x_axis, z_axis, 2, 2)

B1 = LoftSurfaces([Center, R1])
B2 = LoftSurfaces([Center, R2])
B3 = LoftSurfaces([Center, R3])
B4 = LoftSurfaces([Center, R4])
B5 = Box([ 1,-1,0], x_axis, y_axis, 2,2,2)
B6 = Box([-1, 1,0], x_axis, y_axis, 2,2,2)
B7 = Box([-3,-1,0], x_axis, y_axis, 2,2,2)
B8 = Box([-1,-3,0], x_axis, y_axis, 2,2,2)

B4.FlipParametrization(2)
B4.SwapParametrization(1,2)

B6.FlipParametrization(2)
B6.SwapParametrization(1,2)

B7.FlipParametrization(2)
B7.SwapParametrization(0,2)

B2.FlipParametrization(2)
B2.SwapParametrization(0,1)

FinalOutput([B1, B2, B3, B4, B5, B6, B7, B8], True)

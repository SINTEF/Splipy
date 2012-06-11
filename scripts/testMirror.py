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

# output file name
SetFinalOutput("helmet.g2")

# create a rational curve and refine
rc1 = CircleSegment(origin,  # center
                    x_axis,  # start
                    pi,      # angle 
                    z_axis)  # normal
rc1.Translate([0,0,1])

# create a rational surface and refine it
rs1 = RotationalCurveSweep(rc1,    # rotating curve
                           origin, # center point
                           x_axis, # rotational axis
                           pi)     # angle

surf  = ConvertNonRational(rs1);
surf2 = MirrorSurface(surf,  [0, 0, -1] , [0,0,1]);
surf3 = MirrorSurface(surf,  [0, 1, 0]  , [0,1,0]);
surf4 = MirrorSurface(surf2, [0, 1, 0]  , [0,1,0]);

FinalOutput( [surf, surf2, surf3, surf4], True);

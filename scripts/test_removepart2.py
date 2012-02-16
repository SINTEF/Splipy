from GoTools import *
from GoTools.SurfaceFactory import *
from GoTools.SurfaceModelFactory import *

SetDimension(3)
SetFinalOutput("removepart2_3.g2")

# First surface
surf1 = Rectangle(Point(0,-1,0),Point(1,0,0),Point(0,1,0),4,2)

# Second surface
surf2 = TorusSurface(Point(2,0,0),Point(0,1,0),1.4,0.1)
WriteG2("removepart2_1.g2",[surf1,surf2],True)

# Remove part
surf3 = TrimSurface(surf1,surf2)
WriteG2("removepart2_2.g2",surf3)

# Regularize
surf_mod = RegularizeSurface(surf3)
FinalOutput(surf_mod,True)

from math import *
from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *

# standard nice things
o = Point(0,0,0)
x = Point(1,0,0)
y = Point(0,1,0)
z = Point(0,0,1)

box1 = Box(o,    x,  y, 1,1,1)
box2 = Box(x+y,  y, -x, 1,1,1)
box3 = Box(2*z,  x, -z, 1,1,1)
vols = [box1, box2, box3]

mod   = VolumeModel([box1, box2, box3])
surfs = mod.GetShell()

# output file name
SetFinalOutput("vm.g2")
FinalOutput( surfs , True)


from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *

SetDimension(3)
SetTolerance(gap=1.e-5)
SetDebugLevel(1)

ls1 = LineSegment((3,0,0),(4,0,0))
disc_1 = RotationalCurveSweep(ls1,(0,0,0),(0,0,1))

ls2 = LineSegment((3,0,5),(4.3,0,5))
disc_2 = RotationalCurveSweep(ls2,(0,0,5),(0,0,1))

ls3 = LineSegment((2.1,0,10),(3.2,0,10))
disc_3 = RotationalCurveSweep(ls3,(0,0,10),(0,0,1))

ls4 = LineSegment((3,0,15),(4,0,15))
disc_4 = RotationalCurveSweep(ls4,(1.2,0,15),(0,0,1))

ls5 = LineSegment((3,0,20),(4,0,20))
disc_5 = RotationalCurveSweep(ls5,(0,0,20),(0.25,0,1))

tube = LoftSurfaces([disc_1,disc_2,disc_3,disc_4,disc_5])
tube2 = tube + (7,7,-2)

WriteG2("tube_sec.g2",(disc_1,disc_2,disc_3,disc_4,disc_5),False,2)

SetFinalOutput("tube.g2")
FinalOutput((tube,tube2))

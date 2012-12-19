
from math import *
from GeoUtils.Refinement import *
from GoTools import *
from GoTools.VolumeFactory import *

########################################

#######     Initial parameters    ######

########################################

SetFinalOutput(file="fourRefinedPawns.g2")

o = Point(0,0,0)
x = Point(1,0,0)
y = Point(0,1,0)
z = Point(0,0,1)

# read the pawn and uniformly refine it in all directions
vol = ReadG2("pawn.g2")
vol1 = vol.Clone()
vol2 = vol.Clone()
vol3 = vol.Clone()
vol4 = vol.Clone()
UniformVolume(vol1, 1)
UniformVolume(vol2, 2)
UniformVolume(vol3, 3)
UniformVolume(vol4, 0)
UniformVolume(vol4, 0)

# move the pawns around to make it easier to see if it worked
vol1 += 30*x
vol2 += 30*y
vol3 += 30*y + 30*x
vol  += 15*y + 15*x  # unrefined guy at the center

FinalOutput([vol,vol1,vol2,vol3,vol4])


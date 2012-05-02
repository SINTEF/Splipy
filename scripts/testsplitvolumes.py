from GoTools import *
from GoTools.VolumeFactory import *
from math import *

# Unit cube
box = Box([0,0,0],[1,0,0],[0,1,0],1,1,1)

l = box.Split([0.3, 0.7],0)
l += box.Split([0.3, 0.7],1)
l += box.Split([0.3, 0.7],2)
l += box.Split(0.5,0)

FinalOutput(l,True)

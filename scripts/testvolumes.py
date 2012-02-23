from GoTools import *
from GoTools.VolumeFactory import *
from math import *

# Unit cube
box = Box([0,0,0],[1,0,0],[0,1,0],1,1,1)

# Cone
cone = Cone([0,0,1],[0,0,1],pi/4,1,1)

# Cylinder
cylinder = Cylinder([0,0,0],[0,1,0],[0,0,1],1)

# Parallelepiped
pip = Parallelepiped([0,0,0],[1,0,0],[sqrt(2)/2,sqrt(2)/2,0],[0,0,1],1,1,1)

# Sphere
sphere = Sphere([0,0,0],1)

# Torus
torus = Torus([0,0,0],[0,0,1],1,0.7)

FinalOutput([box,cone,cylinder,pip,sphere,torus],True)

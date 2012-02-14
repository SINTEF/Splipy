from GoTools import *

SetDimension(3)

# Create the surface
origo = Point(0,0,0)
surf = Rectangle(origo,Point(1,0,0),Point(0,1,0),1,1)

# Create the hole
circ = Circle(Point(0.5,0.5,0),0.3,Point(0,0,1))
surf_w_hole = AddLoop(surf,circ)

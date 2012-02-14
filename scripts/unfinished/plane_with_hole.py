from GeoModeller import *

SetDimension(3)

# Create the plane
origo = Point(0,0,0)
surf = Rectangle(origo,Point(1,0,0),Point(0,1,0),3,5)

# Create the hole
circ = Circle(Point(1.2,2,0),0.6,Point(0,0,1))
print circ
WriteG2("foo.g2",circ)
surf_w_hole = AddLoop(surf,circ)

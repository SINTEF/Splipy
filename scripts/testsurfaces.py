from GoTools import *

SetDimension(3)

# Infinite plane
p0 = Point(x=1,y=1,z=1.1)
normal = Point(-1,1,-2)
plane = Plane(p0,normal)

# Sphere
p0 = Point(1.4,6,3)
sphere = SphereSurface(p0,4.5)

# Infinite cylinder surface
p0 = Point(0,1,-2)
axis = Point(1,1.3,-0.3)
infcyl = CylinderSurface(p0,axis,2.2)

# Limited cylinder surface
p0 = Point(2,1,-4)
axis = Point(1,1.3,-0.3)
cyl = CylinderSurface(p0,axis,2.4,4.1)

# Infinite cone surface
apex = Point(-1.4,2.1,4.1)
axis = Point(-1.4,0.6,1.9)
infcone = ConeSurface(apex,axis,0.55)

# Limited cone surface
apex = Point(-2.7,7.8,2.1)
axis = Point(-1.4,0.6,1.9)
cone = ConeSurface(apex,axis,0.55,0.0,3.7)

# Torus
p0 = Point(1,13,2)
axis = Point(0,1.4,-0.5)
torus = TorusSurface(p0,axis,3.7,2.2)

# Circular disc
p0 = Point(-2,0,1.1)
p1 = Point(-2.1,0.9,0.6)
normal = Point(0.9,0.6,0.9)
disc = CircularDisc(p0,p1,normal)

print plane
print sphere 
print infcyl
print cyl
print infcone
print cone
print torus
print disc

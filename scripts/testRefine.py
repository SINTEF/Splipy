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

# create a nonrational sin-curve and refine it
pList = list()
for i in range(0,10):
	x = i*2.0*pi/9
	pList.append(Point(x,sin(x), 0))
c1 = InterpolateCurve(pList, range(0,10))
c2 = c1.Clone()
# c1 should still be equal to c2 after refinement
c2.InsertKnot(0.5)
c2.RaiseOrder(2)
c2.InsertKnot(4.5)

# create a rational curve and refine
rc1 = CircleSegment(origin,  # center
                    x_axis,  # start
                    pi,      # angle 
                    z_axis)  # normal
rc2 = rc1.Clone()
rc2.InsertKnot(1.0)
rc2.RaiseOrder(2)
rc2.InsertKnot(1.5)


# create a nonrational surface and refine it
c3 = LineSegment([0, 0, 0], [0, 0, 3])
s1 = LinearCurveSweep(c3,      # first curve
                      c1,      # second curve
                      [0,0,1]) # start point
s2 = s1.Clone()
s2.InsertKnot(0, 0.33)
s2.InsertKnot(1, 4.0)
s2.RaiseOrder(1,3)
s2.InsertKnot(0, 0.66)
s2.InsertKnot(1, 6.0)

# create a rational surface and refine it
rs1 = RotationalCurveSweep(rc1,    # rotating curve
                           origin, # center point
                           x_axis, # rotational axis
                           pi)     # angle
rs1.Translate(z_axis)
rs2 = rs1.Clone()
rs2.InsertKnot(0, 1)
rs2.InsertKnot(1, 1)
rs2.RaiseOrder(1,2)
rs2.InsertKnot(0, 2)
rs2.InsertKnot(1, 2)

# create a nonrational volume and refine it
surf = s1.Clone()
surf.Translate([0,0,4])
v1 = ExtrudeSurface(surf,   # surface to extrude
                    y_axis, # direction to extrude
                    4)      # amount (thickness)
v2 = v1.Clone()
v2.InsertKnot(0, 0.333)
v2.InsertKnot(1, 1)
v2.InsertKnot(2, 0.333)
v2.RaiseOrder(1,2,3)
v2.InsertKnot(0, 0.666)
v2.InsertKnot(1, 5)
v2.InsertKnot(2, 0.666)

# create a rational volume and refine it
surf = rs1.Clone()
surf.Translate([-3,0,4])
rv1 = ExtrudeSurface(surf,   # surface to extrude
                     z_axis, # direction to extrude
                     1)      # amount (thickness)
rv2 = rv1.Clone()
rv2.InsertKnot(0, 0.333)
rv2.InsertKnot(1, 1)
rv2.InsertKnot(2, 0.333)
rv2.RaiseOrder(1,2,3)
rv2.InsertKnot(0, 0.666)
rv2.InsertKnot(1, 2)
rv2.InsertKnot(2, 0.666)

allSplines = [ c1,  c2,
              rc1, rc2,
               s1,  s2,
              rs1, rs2,
               v1,  v2,
              rv1, rv2]
FinalOutput(allSplines, True)

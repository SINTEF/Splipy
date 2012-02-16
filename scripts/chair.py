from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from math import *

SetDimension(3)
SetFinalOutput("chair.g2")

# Some basic geometry
origo      = Point(0.0,0.0,0.0)
x_axis     = Point(1.0,0.0,0.0)
y_axis     = Point(0.0,1.0,0.0)
z_axis     = Point(0.0,0.0,1.0)
neg_x_axis = Point(-1.0,0.0,0.0)

# Chair parameters
xz_thick        = 0.15  # Frame thickness in x- and z-direction
y_thick         = 0.25  # Frame thickness in y-direction
seat_thick      = 0.05  # Thickness of seat and back plate
arm_dist        = 1.7   # Inner distance between the two arm rests
beam_height     = 1.25  # Height from floor to bottom of beams
arm_height      = 1.9   # Height from floor to bottom of intersection of arm rest and back leg
front_back_dist = 1.5   # Inner distance between back and front leg

# Thickness vectors

x_thick_vec = Point(xz_thick,0.0,0.0)
y_thick_vec = Point(0.0,y_thick,0.0)
z_thick_vec = Point(0.0,0.0,xz_thick)

y_gap_vec = Point(0.0,arm_dist,0.0)

# Right front leg and arm rest
fl_1 = LineSegment(Point(front_back_dist+xz_thick,0,0),x_thick_vec,True)
fl_2 = LineSegment(Point(front_back_dist+0.5*xz_thick,0,0.5*beam_height),x_thick_vec,True)
fl_3 = LineSegment(Point(front_back_dist+xz_thick,0,beam_height),x_thick_vec,True)
fl_4 = LineSegment(Point(front_back_dist+xz_thick,0,beam_height+xz_thick),x_thick_vec,True)
fl_5 = LineSegment(Point(front_back_dist,0,arm_height-0.5*xz_thick),x_thick_vec+z_thick_vec,True)
fl_6 = LineSegment(Point(front_back_dist*0.7,0,arm_height),z_thick_vec,True)
fl_7 = LineSegment(Point(front_back_dist*0.4,0,arm_height-0.5*xz_thick),z_thick_vec,True)
fl_8 = LineSegment(Point(xz_thick,0,arm_height),z_thick_vec,True)

leg_rest_surf = LoftCurves([fl_1,fl_2,fl_3,fl_4,fl_5,fl_6,fl_7,fl_8])

ls1 = LineSegment(origo,y_thick_vec)
leg_rest_1 = LinearSurfaceSweep(ls1,leg_rest_surf,origo)

# Other front leg and arm rest
leg_rest_2 = leg_rest_1 + y_thick_vec + y_gap_vec

# Back legs
back_in_rad = 0.5*arm_dist
back_center_height = beam_height+xz_thick+seat_thick+back_in_rad
back_center = Point(xz_thick,y_thick+back_in_rad,back_center_height)

leg_back_1 = Box(origo,x_axis,y_axis,xz_thick,y_thick,back_center_height)
leg_back_2 = leg_back_1 + y_thick_vec + y_gap_vec

# Top half-circle
r = Rectangle(Point(0,0,back_center_height),x_axis,y_axis,xz_thick,y_thick)
up_bend = RotationalSurfaceSweep(r,back_center,neg_x_axis,pi)

# Beams
back_beam  = Box(Point(0.0,y_thick,beam_height),x_axis,y_axis,xz_thick,arm_dist,xz_thick)
front_beam = back_beam + Point(front_back_dist+xz_thick,0,0)
right_beam = Box(Point(xz_thick,y_thick,beam_height),x_axis,y_axis,front_back_dist,y_thick,xz_thick)
left_beam  = right_beam + y_gap_vec-y_thick_vec

# Seat
seat = Box(Point(0,y_thick,xz_thick+beam_height),x_axis,y_axis,front_back_dist+2*xz_thick,arm_dist,seat_thick)

# Back
back = Cylinder(back_center,back_center+(0,back_in_rad,0),neg_x_axis,seat_thick)

FinalOutput([leg_rest_1,leg_rest_2,leg_back_1,leg_back_2,up_bend,back_beam,front_beam,left_beam,right_beam,back,seat],True)

from GoTools import *

SetDimension(3)

# x-coordinate for the vertical stiffener lines
stiff_x_1 = 8.517894104809960
stiff_x_2 = 8.660254037844390
stiff_x_3 = 8.775534171775530
stiff_x_4 = 8.917894104809960

# y-coordinate for the horizontal stiffener lines
stiff_y_1 = 4.6
stiff_y_2 = 5.0

# rotation angle for inner and outer side of sphere

ang_in    = 0.872664625997165
ang_out   = 0.878371594299277

# Make som corners in stiffener sections

pt_2_top = Point(stiff_x_2,0,stiff_y_2)
pt_3_top = Point(stiff_x_3,0,stiff_y_2)
pt_1_bot = Point(stiff_x_1,0,stiff_y_1)
pt_2_bot = Point(stiff_x_2,0,stiff_y_1)
pt_3_bot = Point(stiff_x_3,0,stiff_y_1)

origo = Point(0,0,0)
x_axis = Point(1,0,0)
z_axis = Point(0,0,1)
neg_y_axis = Point(0,-1,0)


# Create surface to be rotated. Start with sphere
circ_in = CircleSegment(origo,pt_2_top,ang_in,neg_y_axis)
circ_out = CircleSegment(origo,pt_3_top,ang_out,neg_y_axis)

sphere_section = LoftCurves([circ_in,circ_out])

# Create stiffener sections
stiff_1_section = Rectangle(pt_1_bot,x_axis,z_axis,stiff_x_2-stiff_x_1,stiff_y_2-stiff_y_1)
stiff_2_section = Rectangle(pt_2_bot,x_axis,z_axis,stiff_x_3-stiff_x_2,stiff_y_2-stiff_y_1)
stiff_3_section = Rectangle(pt_3_bot,x_axis,z_axis,stiff_x_4-stiff_x_3,stiff_y_2-stiff_y_1)

# Store all into surface model
surf_model = SurfaceModel([sphere_section,stiff_1_section,stiff_2_section,stiff_3_section])

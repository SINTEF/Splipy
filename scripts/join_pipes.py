from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from math import *

SetDimension(3)
SetFinalOutput("pipes.g2")

# Input parameters
center = Point(0,0,0)
ax_1 = Point(1,0,0)
ax_2 = Point(2,1,0)
len_1 = 8.0
len_2 = 6.0
r_in = 1.0
r_out = 1.2

# Intersection curves a center
# Define n_ax_i for i = 1,2,3 to be a basis of normalized (but not perpendicular)
# vectors, where n_ax_1 is direction of main pipe, n_ax_2 direction of side pipe,
# and n_ax_3 direction of their cross product

ax_3 = ax_1 % ax_2
n_ax_1 = ax_1 / sqrt(ax_1*ax_1)
n_ax_2 = ax_2 / sqrt(ax_2*ax_2)
n_ax_3 = ax_3 / sqrt(ax_3*ax_3)

# Calculate relative radius along the two semi axis for the ellipse given as
# the intersection of the two tubes
cos_ax = n_ax_1 * n_ax_2
rad_rel_1 = sqrt(2/(1-cos_ax))
rad_rel_2 = sqrt(2/(1+cos_ax))

# Bottom points at intersection of tubes
c_p_in = center + r_in*n_ax_3
c_p_out = center + r_out*n_ax_3

# Intersection curves at centre
c_circ_in_1  = CircleSegment(center,c_p_in,pi,ax_1)
c_circ_out_1 = CircleSegment(center,c_p_out,pi,ax_1)
c_circ_in_2  = EllipticSegment(center,ax_3,r_in,r_in*rad_rel_1,0,pi,n_ax_2-n_ax_1)
c_circ_out_2 = EllipticSegment(center,ax_3,r_out,r_out*rad_rel_1,0,pi,n_ax_2-n_ax_1)
c_circ_in_3  = EllipticSegment(center,ax_3,r_in,r_in*rad_rel_2,0,pi,-(n_ax_1+n_ax_2))
c_circ_out_3 = EllipticSegment(center,ax_3,r_out,r_out*rad_rel_2,0,pi,-(n_ax_1+n_ax_2))

# End curves
end_vec_1 = ax_1 * len_1
end_vec_2 = ax_2 * len_2

e_circ_in_1  = CircleSegment(center+end_vec_1,c_p_in+end_vec_1, pi,ax_1)
e_circ_out_1 = CircleSegment(center+end_vec_1,c_p_out+end_vec_1,pi,ax_1)
e_circ_in_2  = CircleSegment(center+end_vec_1,c_p_in+end_vec_1, pi,-ax_1)
e_circ_out_2 = CircleSegment(center+end_vec_1,c_p_out+end_vec_1,pi,-ax_1)

e_circ_in_3  = CircleSegment(center-end_vec_1,c_p_in-end_vec_1, pi,ax_1)
e_circ_out_3 = CircleSegment(center-end_vec_1,c_p_out-end_vec_1,pi,ax_1)
e_circ_in_4  = CircleSegment(center-end_vec_1,c_p_in-end_vec_1, pi,-ax_1)
e_circ_out_4 = CircleSegment(center-end_vec_1,c_p_out-end_vec_1,pi,-ax_1)

e_circ_in_5  = CircleSegment(center+end_vec_2,c_p_in+end_vec_2, pi,-ax_2)
e_circ_out_5 = CircleSegment(center+end_vec_2,c_p_out+end_vec_2,pi,-ax_2)
e_circ_in_6  = CircleSegment(center+end_vec_2,c_p_in+end_vec_2, pi,ax_2)
e_circ_out_6 = CircleSegment(center+end_vec_2,c_p_out+end_vec_2,pi,ax_2)

# Intersection surfaces at center
c_circ_in_1 = NonRationalCurve(c_circ_in_1)
c_circ_out_1 = NonRationalCurve(c_circ_out_1)
c_face_1 = LoftCurves([c_circ_in_1,c_circ_out_1], order=3)
c_circ_in_2 = NonRationalCurve(c_circ_in_2)
c_circ_out_2 = NonRationalCurve(c_circ_out_2)
c_face_2 = LoftCurves([c_circ_in_2,c_circ_out_2], order=3)
c_circ_in_3 = NonRationalCurve(c_circ_in_3)
c_circ_out_3 = NonRationalCurve(c_circ_out_3)
c_face_3 = LoftCurves([c_circ_in_3,c_circ_out_3], order=3)

# End faces
e_circ_in_1 = NonRationalCurve(e_circ_in_1)
e_circ_out_1 = NonRationalCurve(e_circ_out_1)
e_face_1 = LoftCurves([e_circ_in_1,e_circ_out_1], order=3)
e_circ_in_2 = NonRationalCurve(e_circ_in_2)
e_circ_out_2 = NonRationalCurve(e_circ_out_2)
e_face_2 = LoftCurves([e_circ_in_2,e_circ_out_2], order=3)
e_circ_in_3 = NonRationalCurve(e_circ_in_3)
e_circ_out_3 = NonRationalCurve(e_circ_out_3)
e_face_3 = LoftCurves([e_circ_in_3,e_circ_out_3], order=3)
e_circ_in_4 = NonRationalCurve(e_circ_in_4)
e_circ_out_4 = NonRationalCurve(e_circ_out_4)
e_face_4 = LoftCurves([e_circ_in_4,e_circ_out_4], order=3)
e_circ_in_5 = NonRationalCurve(e_circ_in_5)
e_circ_out_5 = NonRationalCurve(e_circ_out_5)
e_face_5 = LoftCurves([e_circ_in_5,e_circ_out_5], order=3)
e_circ_in_6 = NonRationalCurve(e_circ_in_6)
e_circ_out_6 = NonRationalCurve(e_circ_out_6)
e_face_6 = LoftCurves([e_circ_in_6,e_circ_out_6], order=3)

# Pipe solids
pipe_1 = LoftSurfaces([c_face_1,e_face_1], order=2)
pipe_2 = LoftSurfaces([c_face_2,e_face_2], order=2)
pipe_3 = LoftSurfaces([c_face_1,e_face_3], order=2)
pipe_4 = LoftSurfaces([c_face_3,e_face_4], order=2)
pipe_5 = LoftSurfaces([c_face_3,e_face_5], order=2)
pipe_6 = LoftSurfaces([c_face_2,e_face_6], order=2)

FinalOutput([pipe_1,pipe_2,pipe_3,pipe_4,pipe_5,pipe_6],True)

__doc__ = 'Implementation of various elementary operations on a per-controlpoint level.'

from GoTools import *
from math import *
import numpy as np

def getNonWeightCP(cp):
	if len(cp) == 4:
		w = cp[3]
		return Point(cp[0]/w, cp[1]/w, cp[2]/w)
	else:
		return cp
	

def Mirror(obj, point, normal):
	"""Mirror a curve, surface or volume around a plane
	@param obj: The obj to mirror
	@type obj: Curve, Surface or Volume
	@param point: The point to mirror about
	@type point: Point or list of floats
	@param normal: The plane normal to mirror about
	@type normal: Point or list of floats
	@return: Mirrored result
	@rtype: Curve, Surface or Volume
	"""

	if isinstance(point,list):
		point = Point(list=point)
	if isinstance(normal, list):
		normal = Point(list=normal)
	
	normal.Normalize()
	mirroredCP = []

	rational = (len(obj[0]) == 4)

	for cp in obj:

		# fix rational control points
		if rational:
			w = cp[3]

		# compute the actual end points
		oCP = getNonWeightCP(cp) - point
		movLen = abs(oCP*normal)
		endPos = getNonWeightCP(cp) - 2*movLen*normal

		# fix storage of rational control points
		if rational:
			w = cp[3]
			mirroredCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
		else:
			mirroredCP.append(endPos)
	
	p = obj.GetOrder()
	if type(obj) is Curve:
		knots = obj.GetKnots(with_multiplicities=True)
		return Curve(p, knots, mirroredCP, rational)
	elif type(obj) is Surface:
		knots1, knots2 = obj.GetKnots(with_multiplicities=True)
		return Surface(p[0], p[1], knots1, knots2, mirroredCP, rational)
	elif type(obj) is Volume:
		knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
		return Volume(p[0], p[1], p[2], knots1, knots2, knots3, mirroredCP, rational)

def Translate(obj, vector):
	"""Translate a curve, surface or volume 
	@param obj: The obj to translate
	@type obj: Curve, Surface or Volume
	@param vector: The direction to move the object
	@type vector: Point or list of floats
	@return: Translated result
	@rtype: Curve, Surface or Volume
	"""
	
	if isinstance(vector,list):
		vector = Point(list=vector)

	translatedCP = []

	rational = (len(obj[0]) == 4)
	for cp in obj:

		# fix rational control points
		if rational:
			w = cp[3]

		# compute the actual end points
		endPos = getNonWeightCP(cp) + vector

		# fix storage of rational control points
		if rational:
			w = cp[3]
			translatedCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
		else:
			translatedCP.append(endPos)
	
	p = obj.GetOrder()
	if type(obj) is Curve:
		knots = obj.GetKnots(with_multiplicities=True)
		return Curve(p, knots, translatedCP, rational)
	elif type(obj) is Surface:
		knots1, knots2 = obj.GetKnots(with_multiplicities=True)
		return Surface(p[0], p[1], knots1, knots2, translatedCP, rational)
	elif type(obj) is Volume:
		knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
		return Volume(p[0], p[1], p[2], knots1, knots2, knots3, translatedCP, rational)

def Scale(obj, amount):
	"""Scale a curve, surface or volume 
	@param obj: The obj to scaled
	@type obj: Curve, Surface or Volume
	@param amount: The amount to scale the object with
	@type amount: float
	@return: Scaled result
	@rtype: Curve, Surface or Volume
	"""

	scaledCP = []

	rational = (len(obj[0]) == 4)
	for cp in obj:

		# fix rational control points
		if rational:
			w = cp[3]

		# compute the actual end points
		endPos = getNonWeightCP(cp) * amount

		# fix storage of rational control points
		if rational:
			w = cp[3]
			scaledCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
		else:
			scaledCP.append(endPos)
	
	p = obj.GetOrder()
	if type(obj) is Curve:
		knots = obj.GetKnots(with_multiplicities=True)
		return Curve(p, knots, scaledCP, rational)
	elif type(obj) is Surface:
		knots1, knots2 = obj.GetKnots(with_multiplicities=True)
		return Surface(p[0], p[1], knots1, knots2, scaledCP, rational)
	elif type(obj) is Volume:
		knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
		return Volume(p[0], p[1], p[2], knots1, knots2, knots3, scaledCP, rational)


def rotation_matrix(axis,theta):
	axis = axis/np.sqrt(np.dot(axis,axis))
	a = np.cos(theta/2)
	b,c,d = -axis*np.sin(theta/2)
	return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
	                 [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
	                 [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
	

def Rotate(obj, normal, degrees):
	"""Rotate a curve, surface or volume 
	@param obj: The obj to rotate
	@type obj: Curve, Surface or Volume
	@param amount: The amount to rotate the object with
	@type amount: float
	@return: Scaled result
	@rtype: Curve, Surface or Volume
	"""

	rotatedCP = []
	theta = degrees / 360.0 * 2 * pi;
	print 'theta = ', theta

	if isinstance(normal, list):
		normal = Point(list=normal)

	R = rotation_matrix(normal, theta);

	rational = (len(obj[0]) == 4)
	for cp in obj:

		# fix rational control points
		if rational:
			w = cp[3]

		# do the actual rotation
		c      = getNonWeightCP(cp) ;
		endPos = Point( R[0][0]*c[0] + R[0][1]*c[1] + R[0][2]*c[2] ,
		                R[1][0]*c[0] + R[1][1]*c[1] + R[1][2]*c[2] ,
		                R[2][0]*c[0] + R[2][1]*c[1] + R[2][2]*c[2] )

		# fix storage of rational control points
		if rational:
			w = cp[3]
			rotatedCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
		else:
			rotatedCP.append(endPos)
	
	p = obj.GetOrder()
	if type(obj) is Curve:
		knots = obj.GetKnots(with_multiplicities=True)
		return Curve(p, knots, rotatedCP, rational)
	elif type(obj) is Surface:
		knots1, knots2 = obj.GetKnots(with_multiplicities=True)
		return Surface(p[0], p[1], knots1, knots2, rotatedCP, rational)
	elif type(obj) is Volume:
		knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
		return Volume(p[0], p[1], p[2], knots1, knots2, knots3, rotatedCP, rational)

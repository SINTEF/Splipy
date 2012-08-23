__doc__ = 'Implementation of various interpolation schemes.'

from GoTools import *

def MirrorVolume(volume, point, normal):
	"""Mirror a volume around a plane
	@param volume: The volume to mirror
	@type volume: Volume
	@param point: The point to mirror about
	@type point: Point
	@param normal: The plane normal to mirror about
	@type normal: Point
	@return: Mirrored volume
	@rtype: Volume
	"""
	
	normal.Normalize()
	mirroredCP = []

	rational = (len(volume[0]) == 4)
	for cp in volume:

		# fix rational control points
		if rational:
			w = cp[3]
			cpNonWeight = Point(cp[0]/w, cp[1]/w, cp[2]/w)
		else:
			cpNonWeight = cp

		# compute the actual end points
		oCP = cpNonWeight - point
		movLen = abs(oCP*normal)
		endPos = cpNonWeight - 2*movLen*normal

		# fix storage of rational control points
		if rational:
			w = cp[3]
			mirroredCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
		else:
			mirroredCP.append(endPos)
	
	p = volume.GetOrder()
	knots1, knots2, knots3 = volume.GetKnots(with_multiplicities=True)
	return Volume(p[0], p[1], p[2], knots1, knots2, knots3, mirroredCP, rational)



__doc__ = 'Implementation of various refinement schemes.'

from math import atan
from math import pi

def UniformCurve(curve):
  """Uniformly refine a curve by halfing each knot interval
  @param curve: The curve to refine
  @type curve: Curve
  @return: None
  """

  knots = curve.GetKnots()
  for i in range(0,len(knots)-2):
    curve.InsertKnot((knots[i]+knots[i+1])/2)

# Chop each knot span in half
def UniformSurface(surface, direction=0):
  """Uniformly refine a surface by halfing each knot interval
  @param surface: The surface to refine
  @type surface: Surface 
  @param direction: The direction to refine in (0 = both, 1, 2)
  @type direction: integer
  @return: None
  """
  knots_u, knots_v = surface.GetKnots()
  if direction == 0 or direction == 1:
    for i in range(0,len(knots_u)-1):
      surface.InsertKnot(0,(knots_u[i]+knots_u[i+1])/2)
  if direction == 0 or direction == 2:
    for i in range(0,len(knots_v)-1):
      surface.InsertKnot(1,(knots_v[i]+knots_v[i+1])/2)

# Chop each knot span in half
def UniformVolume(volume, direction=0):
  """Uniformly refine a volume by halfing each knot interval
  @param volume: The volume to refine
  @type volume: Volume 
  @param direction: The direction to refine in (0 = both, 1, 2, 3)
  @type direction: integer
  @return: None
  """
  knots_u, knots_v, knots_w = volume.GetKnots()
  if direction == 0 or direction == 1:
    for i in range(0,len(knots_u)-1):
      volume.InsertKnot(0,(knots_u[i]+knots_u[i+1])/2)
  if direction == 0 or direction == 2:
    for i in range(0,len(knots_v)-1):
      volume.InsertKnot(1,(knots_v[i]+knots_v[i+1])/2)
  if direction == 0 or direction == 3:
    for i in range(0,len(knots_w)-1):
      volume.InsertKnot(2,(knots_w[i]+knots_w[i+1])/2)


def BoundaryLayerCurve(curve, start, scale, n):
  """Refine a curve with geometric grading
  @param curve: The curve to refine
  @type curve: Curve
  @param start: The end of the curve to have grading towards (1 or 2)
  @type start: integer
  @param scale: The geometric scaling factor
  @type scale: Float
  @param n: The total number of control points
  @type n: integer
  @return: None
  """
  knots = curve.GetKnots()
  pwr  = 1
  smm = 0.0
  for i in range(0,n):
    smm += pwr
    pwr *= scale
  if start == 1:
    minor = knots[1]
    major = knots[0]
  if start == 2:
    minor = knots[-2]
    major = knots[-1]

  alpha = 1/smm
  smm = 0.0
  pwr = 1.0
  for i in range(0,n-1):
    smm += alpha*pwr
    pwr *= scale
    curve.InsertKnot(smm*major+(1.0-smm)*minor)

# Geometric distribution of knots
def GeometricRefineSurface(surface, direction, alpha, n):
	"""Refine a surface by making a geometric distribution of element sizes
	Consider Surface.FlipParametrization if you need refinement towards the other edge
	@param surface: The surface to refine
	@type surface: Surface 
	@param direction: The direction to refine in (u=1 or v=2) 
	@type direction: int
	@param alpha: The ratio between two sequential knot segments
	@type alpha: float
	@param n: The number of knots to insert
	@type n: int
	@return: None
	"""
	
	# some error tests on input
	if n<=0:
		print 'n should be greater than 0'
		return None

	flipBack = False
	if direction < 0:
		surface.FlipParametrization(-direction-1)
		direction = -direction
		flipBack = True

	# fetch knots
	knots_u, knots_v = surface.GetKnots()
	if direction == 1:
		knotStart = knots_u[0]
		knotEnd   = knots_u[-1]
	elif direction == 2:
		knotStart = knots_v[0]
		knotEnd   = knots_v[-1]
	else:
		print 'Direction should be 1 or 2'
		return None
	dk = knotEnd - knotStart

	# evaluate the factors
	n = n+1 # redefine n to be knot spans instead of new internal knots
	totProd = 1.0
	totSum  = 0.0
	for i in range(n):
		totSum  += totProd
		totProd *= alpha
	d1 = 1.0 / totSum
	knot = d1

	# do the actual knot insertion
	for i in range(n-1):
		surface.InsertKnot(direction-1, knotStart + knot*dk)
		knot += alpha*d1
		d1   *= alpha

	if flipBack:
		surface.FlipParametrization(direction-1)


# Edge refinement
def EdgeRefineSurface(surface, direction, S, n):
	"""Refine a surface by both edges, by sampling a atan-function
	@param surface: The surface to refine
	@type surface: Surface 
	@param direction: The direction to refine in (u=1 or v=2) 
	@type direction: int
	@param S: The slope of the atan-function
	@type S: float
	@param n: The number of knots to insert
	@type n: int
	@return: None
	"""
	
	# some error tests on input
	if n<=0:
		print 'n should be greater than 0'
		return None

	# fetch knots
	knots_u, knots_v = surface.GetKnots()
	if direction == 1:
		knotStart = knots_u[0]
		knotEnd   = knots_u[-1]
	elif direction == 2:
		knotStart = knots_v[0]
		knotEnd   = knots_v[-1]
	else:
		print 'Direction should be 1 or 2'
		return None
	dk = knotEnd - knotStart

	# evaluate the factors
	newKnots = []
	maxAtan  = atan(S)
	for i in range(1,n+1):
		xi  = -1.0 + 2.0*i/(n+1)
		xi *= S
		newKnots.append(knotStart + (atan(xi)+maxAtan)/2/maxAtan*dk)

	# do the actual knot insertion
	for x in newKnots:
		surface.InsertKnot(direction-1, x)



__doc__ = 'Implementation of various refinement schemes.'

def UniformCurve(curve):
  """Uniformly refine a curve by halfing each knot interval
  @param curve: The curve to refine
  @type curve: Curve
  @return: None
  """

  knots = curve.GetKnots()
  for i in range(0,len(knots)-2):
    curve.InsertKnot((knots[i]+knots[i+1])/2)

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
    for i in range(0,len(knots_u)-2):
      surface.InsertKnot(0,(knots_u[i]+knots_u[i+1])/2)
  if direction == 0 or direction == 2:
    end = len(knots_v)-2
    if end == 0:
      end = 1;
    for i in range(0,end):
      surface.InsertKnot(1,(knots_v[i]+knots_v[i+1])/2)

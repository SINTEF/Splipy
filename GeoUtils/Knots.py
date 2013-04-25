__doc__ = 'Implementation of knot vector helpers.'

from GoTools import *

def FullContinuity(params, order):
  """ Create a full-continuity knot vector
      @param params: The parameter values
      @type params: List of float
      @param order: The order of the spline
      @type order: integer
      @return: The knot vector
      @rtype: List of float
  """
  result = []
  for i in range(order-1):
    result.append(params[0])
  result += params 
  for i in range(order-1):
    result.append(params[-1])

  return result

def KnotExist(knotVec, knot, ignoreEnds=True):
  """" Check if knot already exist in the object based on the global refine tolerance
       @param knotVec: The knot vector to check against
       @type knotVec:  List of floats
       @param knot: The new knot to insert
       @type knot:  Float
       @param ignoreEnds: If the first and last knot should be ignored 
       @type ignoreEnds:  Bool
       @return: If the knot already exist
       @rtype: Bool
       """

  tol = GetTolerance('refine')
  start = knotVec[0]
  end   = knotVec[-1]
  for k in knotVec:
    # note that one may pass knot vector with multiplicities here
    # which means that the endpoints may appear multiple places in
    # the knot vector
    if (k == start or k == end) and ignoreEnds:
      continue
    if abs(k-knot) < tol:
      return True
  return False

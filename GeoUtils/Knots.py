__doc__ = 'Implementation of knot vector helpers.'

from GoTools import *

def FullContinuity(params, order):
  """ Create a full-continuity knot vector
      @param params: The parameter values
      @type params: List of float
      @param order: The order of the spline
      @type order: Int
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
  """ Check if knot already exist in the object based on the global refine tolerance
      @param knotVec: The knot vector to check against
      @type knotVec:  List of floats
      @param knot: The new knot to insert
      @type knot:  Float
      @param ignoreEnds: If the first and last knot should be ignored 
      @type ignoreEnds:  Bool
      @return: If the knot already exist
      @rtype: Bool
  """

  tol    = GetTolerance('refine')
  start  = knotVec[0]
  end    = knotVec[-1]
  result = 0;
  for k in knotVec:
    # note that one may pass knot vector with multiplicities here
    # which means that the endpoints may appear multiple places in
    # the knot vector
    if (k == start or k == end) and ignoreEnds:
      continue
    if abs(k-knot) < tol:
      result += 1
  return result

def GetGrevillePoints(knotVec):
  """" Generate the list of all greville points ascociated with a knot vector. Polynomial degree of 
       the basis is implicitly given by the first knot multiplicity
       @param knotVec: The knot vector to check against
       @type knotVec:  List of floats
       @return: The greville points
       @rtype: List of floats
       """

  p = 1
  while knotVec[p] == knotVec[0]:
    p += 1
  n = len(knotVec) - p
  results = []
  for i in range(n):
    pt = 0
    for j in range(i+1, i+p):
      pt += knotVec[j]
    results.append(pt/(p-1))
  return results


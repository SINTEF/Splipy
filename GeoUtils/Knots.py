__doc__ = 'Implementation of knot vector helpers.'

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

__doc__ = 'Implementation of knot vector helpers.'

def FullContinuity(params, order):
  result = []
  for i in range(order-1):
    result.append(params[0])
  result += params 
  for i in range(order-1):
    result.append(params[-1])

  return result

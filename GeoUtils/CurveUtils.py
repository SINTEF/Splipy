__doc__ = 'Implementation of various curve utilities'

def CurveLengthParametrization(pts,normalize=False):
    """Get knots corresponding to a curvelength parametrization of a spline'
    @param pts      : The nodes of the control polygon
    @type  pts      : List of float
    @param normalize: Whether or not to normalize parametrization to 1.0
    @type  normalize: Boolean
    @return         : The parametrization
    @rtype          : List of float
    """
    # Number of points
    np = len(pts)
    
    # Curve length parametrization
    s = 0.0
    knots = []
    knots.append(s)
    for i in range(1,np):
        ds = abs(pts[i]-pts[i-1])
        s = s + ds
        knots.append(s)
    
    # Curve length (approximation)
    s = knots[-1]
    
    # Normalized curve length parametrization
    if (normalize):
        for i in range(0,np):
            knots[i] = knots[i]/s
    
    return knots

def GetCurvePoints(curve):
    """Get value of a given curve in all its knots
    @param curve: The curve
    @type  curve: Curve
    @return     : Value of curve in parameters
    @rtype      : List of points
    """
    knots = curve.GetKnots()
    
    pts = []
    for xi in knots:
        pts.append(curve.Evaluate(xi))
    
    return pts

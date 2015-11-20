__doc__ = 'Implementation of various smoothing operations on a per-controlpoint level.'

from GoTools import *

def getNonWeightCP(cp):
    if len(cp) == 4:
        w = cp[3]
        return Point(cp[0]/w, cp[1]/w, cp[2]/w)
    else:
        return cp

def SmoothSurfInterior(obj, nIter):
    """Smooth a surface by setting the interior control points to the average of itself and 8 neighbours
    @param obj  : The surface to smooth
    @type  obj  : Surface
    @param nIter: Number of smoothing iterations
    @type  nIter: Int
    @return     : Smoothed result
    @rtype      : Surface
    """

    p              = obj.GetOrder()
    knots1, knots2 = obj.GetKnots(with_multiplicities=True)
    n = [len(knots1) - p[0], len(knots2) - p[1]]
    smoothedCP = []
    newCP      = []
    for cp in obj:
        smoothedCP.append(cp)
        newCP.append(cp)

    rational = (len(obj[0]) == 4)
    for iteration in range(nIter):
        for j in range(1,n[1]-1):
            for i in range(1,n[0]-1):
                cpI = j*n[0] + i

                # fix rational control points
                if rational:
                    w = smoothedCP[cpI][3]

                cp  = Point(0,0,0)
                for jj in range(-n[0],n[0]+1, n[0]):
                    for ii in range(-1,2):
                        cp += getNonWeightCP(smoothedCP[cpI+ii+jj])
                cp  /= 9


                # fix storage of rational control points
                if rational:
                    w = cp[3]
                    newCP[cpI] = Point(list=[cp[0]*w, cp[1]*w, cp[2]*w, w])
                else:
                    newCP[cpI] = cp
        smoothedCP = newCP

    return Surface(p[0], p[1], knots1, knots2, smoothedCP, rational)

def SmoothVolInterior(obj, nIter):
    """Smooth a volume by setting the interior control points to the average of itself and 26 neighbours
    @param obj  : The volume to smooth
    @type  obj  : Volume
    @param nIter: Number of smoothing iterations
    @type  nIter: Int
    @return     : Smoothed result
    @rtype      : Volume
    """

    p              = obj.GetOrder()
    knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
    n = [len(knots1) - p[0], len(knots2) - p[1], len(knots3) - p[2]]
    smoothedCP = []
    newCP      = []
    for cp in obj:
        smoothedCP.append(cp)
        newCP.append(cp)

    rational = (len(obj[0]) == 4)
    for iteration in range(nIter):
        for k in range(1,n[2]-1):
            for j in range(1,n[1]-1):
                for i in range(1,n[0]-1):
                    cpI = k * n[0]*n[1] + j*n[0] + i

                    # fix rational control points
                    if rational:
                        w = smoothedCP[cpI][3]

                    cp  = Point(0,0,0)
                    for kk in range(-n[0]*n[1],n[0]*n[1]+1, n[0]*n[1]):
                        for jj in range(-n[0],n[0]+1, n[0]):
                            for ii in range(-1,2):
                                cp += getNonWeightCP(smoothedCP[cpI+ii+jj+kk])
                    cp  /= 27


                    # fix storage of rational control points
                    if rational:
                        w = cp[3]
                        newCP[cpI] = Point(list=[cp[0]*w, cp[1]*w, cp[2]*w, w])
                    else:
                        newCP[cpI] = cp
        smoothedCP = newCP

    return Volume(p[0], p[1], p[2], knots1, knots2, knots3, smoothedCP, rational)

def LaplaceSmoothing(crv,r=1.0,nsweep=1):
    """Smooth a curve by applying a laplacian
    @param crv   : The curve to smooth
    @type  crv   : Curve
    @param r     : Smoothing ratio
    @type  r     : Float
    @param nsweep: Number of sweeps to apply
    @type  nsweep: int
    @return      : Smoothed curve
    @rtype       : Curve
    """
    pts1   = GetCurvePoints(crv)
    knots1 = CurveLengthParametrization(pts1,False)
    n1     = len(knots1)
    
    dxi = []
    dxi1 = []
    for i in range(1,n1):
        dxi.append(knots1[i]-knots1[i-1])
        dxi1.append(knots1[i]-knots1[i-1])
    
    for n in range(0,nsweep):
        dxi[0] = r*dxi1[0] + 0.5*(1.0-r)*(dxi1[-1]+dxi1[1])
        for i in range(1,n1-2):
            dxi[i] = r*dxi1[i] + 0.5*(1.0-r)*(dxi1[i-1]+dxi1[i+1])
        dxi[n1-2] = r*dxi1[n1-2] + 0.5*(1.0-r)*(dxi1[n1-3] + dxi1[0])
        
        for i in range(0,n1-1):
            dxi1[i] = dxi[i]
    
    s1 = knots1[-1]
    s2 = 0.0
    for ds in dxi:
        s2 = s2 + ds
    
    knots = []
    s = 0.0
    knots.append(s)
    for i in range(0,n1-1):
        s = s + dxi[i]*s1/s2
        knots.append(s)
    
    pts = []
    pts.append(pts1[0])
    for i in range(1,n1-1):
        pts.append(crv.Evaluate(knots[i]))
    pts.append(pts1[-1])
    
    t1 = crv.EvaluateTangent(knots1[0])
    t2 = crv.EvaluateTangent(knots1[-1])
    curve = InterpolateCurve(pts,knots,t1,t2)
    return curve

__doc__ = 'Implementation of various factory methods which should have been in GoTools, but are missing or bugged'

from GoTools import *
from GoTools.SurfaceFactory import *
from GeoUtils.Knots import *
from GeoUtils.Elementary import *
from GoTools.CurveFactory import *

def SplitSurface(original, param, value):
    """Split a SplineSurface in two by continued knot insertion
    @param original: The Surface to split
    @type  original: Surface
    @param    param: The parametric direction to split along (u=0, v=1)
    @type     param: Int
    @param    value: The parametric value to split along
    @type     value: Float
    @return:         The two new surfaces
    @rtype:          List of Surfaces
    """

    surf      = original.Clone()
    dim       = GetDimension()
    rational  = (len(surf[0]) == (dim+1))
    knot      = surf.GetKnots(True)
    p         = surf.GetOrder()

    # error check input (Note: that rational should really fail unless you know what you're doing)
#         if rational:
#                 return None
    if( value <= knot[param][0] or knot[param][-1] <= value):
      return None

    # Figure out if knot already exist
    knotCount = KnotExist(knot[param], value)

    # insert the actual knot
    for i in range(knotCount, p[param]-1):
        surf.InsertKnot(param, value)

    knot = surf.GetKnots(True)
    # figure out which knot interval (ki) the split lies
    for ki in range(len(knot[param])):
        if abs(knot[param][ki]-value) < GetTolerance('refine'):
            break
    ki = ki-1

    n = [len(knot[0]) - p[0] , len(knot[1]) - p[1]]

    coefs1 = []
    coefs2 = []
    knot1u = []
    knot1v = []
    knot2u = []
    knot2v = []

    k = 0
    if param == 0:
        for y in range(n[1]):
            for x in range(n[0]):
                if x <= ki:
                    coefs1.append(surf[k])
                if x >= ki:
                    coefs2.append(surf[k])
                k = k+1

        for i in range(ki+p[0]):
                knot1u.append(knot[0][i])
        knot1u.append(knot1u[-1])
        for i in range(ki+1, len(knot[0])):
                knot2u.append(knot[0][i])
        knot2u = [knot2u[0]] + knot2u

        knot1v = knot[1]
        knot2v = knot[1]
    else:
        for y in range(n[1]):
            for x in range(n[0]):
                if y <= ki:
                    coefs1.append(surf[k])
                if y >= ki:
                    coefs2.append(surf[k])
                k = k+1

        for i in range(ki+p[1]):
            knot1v.append(knot[1][i])
        knot1v.append(knot1v[-1])
        for i in range(ki+1, len(knot[1])):
            knot2v.append(knot[1][i])
        knot2v = [knot2v[0]] + knot2v

        knot1u = knot[0]
        knot2u = knot[0]


    result = []
    result.append(Surface(p[0], p[1], knot1u, knot1v, coefs1, rational))
    result.append(Surface(p[0], p[1], knot2u, knot2v, coefs2, rational))
    return  result


def SplitVolume(original, param, value):
    """Split a SplineVolume in two by continued knot insertion
    @param original: The Volume to split
    @type  original: Volume
    @param    param: The parametric direction to split along (u=0, v=1, w=2)
    @type     param: Int
    @param    value: The parametric value to split along
    @type     value: Float
    @return:         The two new volumes
    @rtype:          List of Volumes
    """

    vol       = original.Clone()
# this actually fails for 2D models, but no Surface.GetDim() or
# Surface.IsRational()  :(
    rational  = (len(vol[0]) == 4)
    knot      = vol.GetKnots(True)
    p         = vol.GetOrder()

    # error check input (Note: that rational should really fail unless you know what you're doing)
#         if rational:
#                 return None
    if( value <= knot[param][0] or knot[param][-1] <= value):
      return None

    # Figure out if knot already exist
    knotCount = KnotExist(knot[param], value)

    # insert the actual knot
    for i in range(knotCount, p[param]-1):
        vol.InsertKnot(param, value)

    knot = vol.GetKnots(True)
    # figure out which knot interval (ki) the split lies
    for ki in range(len(knot[param])):
        if abs(knot[param][ki]-value) < GetTolerance('refine'):
            break
    ki = ki-1

    n = [len(knot[0]) - p[0] , len(knot[1]) - p[1], len(knot[2]) - p[2]]

    coefs1 = []
    coefs2 = []
    knot1u = []
    knot1v = []
    knot1w = []
    knot2u = []
    knot2v = []
    knot2w = []

    k = 0
    if param == 0:
        for z in range(n[2]):
            for y in range(n[1]):
                for x in range(n[0]):
                    if x <= ki:
                        coefs1.append(vol[k])
                    if x >= ki:
                        coefs2.append(vol[k])
                    k = k+1

        for i in range(ki+p[0]):
                knot1u.append(knot[0][i])
        knot1u.append(knot1u[-1])
        for i in range(ki+1, len(knot[0])):
                knot2u.append(knot[0][i])
        knot2u = [knot2u[0]] + knot2u

        knot1v = knot[1]
        knot2v = knot[1]
        knot1w = knot[2]
        knot2w = knot[2]
    elif param == 1:
        for z in range(n[2]):
            for y in range(n[1]):
                for x in range(n[0]):
                    if y <= ki:
                        coefs1.append(vol[k])
                    if y >= ki:
                        coefs2.append(vol[k])
                    k = k+1

        for i in range(ki+p[1]):
            knot1v.append(knot[1][i])
        knot1v.append(knot1v[-1])
        for i in range(ki+1, len(knot[1])):
            knot2v.append(knot[1][i])
        knot2v = [knot2v[0]] + knot2v

        knot1u = knot[0]
        knot2u = knot[0]
        knot1w = knot[2]
        knot2w = knot[2]
    else:
        for z in range(n[2]):
            for y in range(n[1]):
                for x in range(n[0]):
                    if z <= ki:
                        coefs1.append(vol[k])
                    if z >= ki:
                        coefs2.append(vol[k])
                    k = k+1

        for i in range(ki+p[2]):
            knot1w.append(knot[2][i])
        knot1w.append(knot1w[-1])
        for i in range(ki+1, len(knot[2])):
            knot2w.append(knot[2][i])
        knot2w = [knot2w[0]] + knot2w

        knot1u = knot[0]
        knot2u = knot[0]
        knot1v = knot[1]
        knot2v = knot[1]


    result = []
    result.append(Volume(p[0], p[1], p[2], knot1u, knot1v, knot1w, coefs1, rational))
    result.append(Volume(p[0], p[1], p[2], knot2u, knot2v, knot2w, coefs2, rational))
    return  result


def RotationalCrv2CrvSweep(crv1, crv2, axis):
    """ Rotational sweep the surface in between to curves,
    where both curves live on the same cylinder.
    CONSTRAINT: Curves must share the height-parametrization, i.e.
    share the z-coordinates if living on the unit cylinder
    @param crv1: The first curve on the cylinder
    @type  crv1: Curve
    @param crv2: The second curve on the cylinder
    @type  crv2: Curve
    @param axis: The cylinder axis; 0=x-axis, 1=y-axis, 2=z-axis
    @type  axis: Int
    @return:     Swept surface
    @rtype:      Surface
    """
    controlPoints = []

    listAxis = [0,0,0]
    listAxis[axis] = 1.0
    nAxis = Point(list=listAxis)

    for i in range(len(crv1)):

#                   A
#          , - ~ ~ ~ -0,
#      , '           /   ' ,
#    ,              /        ;,
#   ,           r1 /          ,'.  ControlPt 1
#  ,              /            , '0
#  ,             * ) Theta     , /
#  ,              \            ,/
#   ,              \          ,0  ControlPt 2
#    ,           r2 \        ,/
#      ,             \    , '/
#        ' - , _ _ _ ,0_"___0  ControlPt 3
#                     B
#
#       going for two knot spans. Should work for all theta
#       strictly less than 2pi (though bad results for close
#       to 2pi
#

        a = crv1[i]
        b = crv2[i]
        origin1 = (a*nAxis)*nAxis
        origin2 = (b*nAxis)*nAxis
        dist    = abs(origin1-origin2)
        r1      = abs(origin1-a)
        r2      = abs(origin2-b)

        # !!! HARDCODE FIX !!!
        # if the first curve starts or stops tangent to the cylinder radial parametrization
        # we try and fix this by tweaking the surrounding control points. In this case, setting
        # the second CP to the average height between the first two.
        if i==1 and abs(crv1[i][axis]-crv1[i-1][axis]) < 1e-3:
            # print 'fixing height'
            newHeight = (crv1[i+1][axis]+crv1[i][axis]) / 2.0
            origin1 = Point( newHeight, origin1[1], origin1[2])
            origin2 = Point( newHeight, origin2[1], origin2[2])
            a = Point(newHeight, a[1], a[2])
            b = Point(newHeight, b[1], b[2])
        elif i==len(crv1)-2 and abs(crv1[i][axis]-crv1[i+1][axis]) < 1e-3:
            # print 'fixing height'
            newHeight = (crv1[i-1][axis]+crv1[i][axis]) / 2.0
            origin1 = Point( newHeight, origin1[1], origin1[2])
            origin2 = Point( newHeight, origin2[1], origin2[2])
            a = Point(newHeight, a[1], a[2])
            b = Point(newHeight, b[1], b[2])


        ###   this test isn't really true since the control points are floating outside
        ###   the actual curve in geometry space (but nice to have to illustrate whats
        ###   going on)
        # if abs(r1-r2) > GetTolerance('gap'):
        #         print 'Error. Curves not on the same cylinder radius'
        #         return

        if dist > GetTolerance('gap'):
            print 'Error. Curves not on the same cylinder height'
            return

        l1 = a-origin1
        l2 = b-origin2
        if abs(l1*l2/abs(l1)/abs(l2)-1.0) < 1e-14 or \
           abs(l1*l2/abs(l1)/abs(l2)+1.0) < 1e-14:
            theta = pi
        else:
            theta = acos((l1*l2)/abs(l1)/abs(l2))
#                 if((l1 % l2)[axis] < 0):
#                         theta = - theta
        dt    = theta / 4.0
        dr    = (r2-r1)/4.0

        rotPoint = l1.Clone()
        rotPoint.Normalize()

        # kjetil hack for nasty patches
        controlPoints.extend([crv1[i][0], crv1[i][1], crv1[i][2], 1.0])
        rotPoint.Rotate(nAxis, dt)
        rotPoint = rotPoint / cos(dt)

        for j in range(1,5):
            r = r1 + j*dr

            w = cos(dt)*(j%2) + (1.0-(j%2))
            controlPoints.append((origin1[0] + r*rotPoint[0]) * w)
            controlPoints.append((origin1[1] + r*rotPoint[1]) * w)
            controlPoints.append((origin1[2] + r*rotPoint[2]) * w)
            controlPoints.append(w)

            rotPoint.Rotate(nAxis, dt)
            if j%2==0:
                rotPoint = rotPoint / cos(dt)
            else:
                rotPoint = rotPoint * cos(dt)

    knotXi = [0,0,0,1,1,2,2,2]

    surf   = Surface(3, crv1.GetOrder(), knotXi, crv1.GetKnots(True), controlPoints, True)

    return surf

def LoftBetween(obj1, obj2):
    """Linear lofing the domain between two objects. If the arguments
    are curves, the extruded domain will be a surface, if the arguments
    are surfaces, the domain will be the contained volume.
    @param obj1: The first object
    @type obj1: Curve or Surface
    @param obj2: The second object
    @type obj2: Curve or Surface
    @return: The contained domain
    @rtype: Surface or Volume
    """

    cObj1 = obj1.Clone()
    cObj2 = obj2.Clone()

    MakeCommonSplineSpace(cObj1, cObj2)

    newCP = []
    for c in cObj1:
        newCP.append(c)
    for c in cObj2:
        newCP.append(c)

    p = cObj1.GetOrder()
    rational = len(cObj1[0]) == 4
    newKnot = [0, 0, 1, 1]  # linear in the extrusion direction
    newP    = 2

    if type(cObj1) is Curve:
        knots = cObj1.GetKnots(with_multiplicities=True)
        return Surface(p, newP, knots, newKnot, newCP, rational)
    elif type(cObj1) is Surface:
        knots1, knots2 = cObj1.GetKnots(with_multiplicities=True)
        return Volume(p[0], p[1], newP, knots1, knots2, newKnot, newCP, rational)

def Thicken(obj, amount):
    """Extrudes a surface in the normal direction to create a volume representation of it.
    May cause self intersections for large amount or high curvature surfaces
    @param    obj: The surface to thicken
    @type     obj: Surface
    @param amount: The distance to extrude
    @type  amount: Float
    @return:       Thickened surface
    @rtype:        Volume
    """

    # this is not a hard restriction. Just have to come up with another
    # sampling point strategy
    p1,p2 = obj.GetOrder()
    if p1>4 or p2>4:
      print 'Surfaces of higher order than 4 not supported for thicken command'
      return None

    # sampling point choice: use the greville point of the corresponding bicubic surface representation
    sampleObj = obj.Clone()
    sampleObj.RaiseOrder(4-p1, 4-p2)
    knot1, knot2 = sampleObj.GetKnots(True)
    g1 = GetGrevillePoints(knot1)
    g2 = GetGrevillePoints(knot2)

    outPts = []
    inPts = []
    for v in g2:
      for u in g1:
        p = sampleObj.Evaluate(u,v)
        n = sampleObj.EvaluateNormal(u,v)
        outPts.append( p + 0.5*amount * n)
        inPts.append(  p - 0.5*amount * n)

    outShell = InterpolateSurface(outPts, g1, g2)
    inShell  = InterpolateSurface(inPts,  g1, g2)

    return LoftBetween(inShell, outShell)

def HyperEllipse(radius, center, N, order=4, quadrant=0):
  """Generate a hyperellipse in the X-Y plane
  @param radius: Radius of the original circle
  @type radius: Float
  @param center: Origo of the circle
  @type center: Point
  @param N: Number of spline coefficients on curve
  @type N: Integer
  @param order: The order of the hyperellipse. 2 for a Squircle (default)
  @type order: Integer
  @param quadrant: Quadrants to include. 0 for all quadrants
  @type quadrant: Integer
  @return: New (cubic) curve describing the hyperellipse
  @rtype: Curve
  """
  if quadrant == 0:
    N = N/4
  vals = []
  if quadrant == 0 or quadrant == 1:
    vals = [[center[0]+radius[0]*pow(cos(knot),2.0/order),center[1]+radius[1]*pow(sin(knot),2.0/order),center[2]] for knot in [i*pi/(2*(N-1)) for i in range(N)]]
  if quadrant == 0 or quadrant == 2:
    val2 = [[center[0]-radius[0]*pow(cos(knot),2.0/order),center[1]+radius[1]*pow(sin(knot),2.0/order),center[2]] for knot in [i*pi/(2*(N-1)) for i in range(N)]]
    val2.reverse()
    if quadrant == 0:
      del val2[0]
    vals += val2
  if quadrant == 0 or quadrant == 3:
    val2 = [[center[0]-radius[0]*pow(cos(knot),2.0/order),center[0]-radius[1]*pow(sin(knot),2.0/order),center[2]] for knot in [i*pi/(2*(N-1)) for i in range(N)]]
    if quadrant == 0:
      del val2[0]
    vals += val2
  if quadrant == 0 or quadrant == 4:
    val2 = [[center[0]+radius[0]*pow(cos(knot),2.0/order),center[0]-radius[1]*pow(sin(knot),2.0/order),center[2]] for knot in [i*pi/(2*(N-1)) for i in range(N)]]
    if quadrant == 0:
      del val2[0]
    val2.reverse()
    vals += val2

  # Make loop
  if quadrant==0:
    vals.append(vals[0])

  return InterpolateCurve(vals, range(len(vals)))

def CombineSurfaces(surf1, surf2):
    """Combine two surfaces into a single patch with internal C^0 knot
    @param surf1: The first Surface
    @type  surf1: Surface
    @param surf2: The second Surface
    @type  surf2: Surface
    @return:      The combined surface
    @rtype:       Surface
    """

    # this actually fails for 2D models, but no Surface.GetDim() or Surface.IsRational()  :(
    rational  = (len(surf1[0]) == 4)
    knot1 = surf1.GetKnots(True)
    p1    = surf1.GetOrder()
    n1    = [len(knot1[0]) - p1[0] , len(knot1[1]) - p1[1]]

    knot2 = surf2.GetKnots(True)
    p2    = surf2.GetOrder()
    n2    = [len(knot2[0]) - p2[0] , len(knot2[1]) - p2[1]]

    # swap parametrizations since we want the sem at vmax
    if n1[1]==n2[1] and p1[1]==p2[1]:
        surf1.SwapParametrization()
        surf2.SwapParametrization()
    elif n1[1]==n2[0] and p1[1]==p2[0]:
        surf1.SwapParametrization()
    elif n1[0]==n2[1] and p1[0]==p2[1]:
        surf2.SwapParametrization()

    # update discretization values
    knot1 = surf1.GetKnots(True)
    p1    = surf1.GetOrder()
    n1    = [len(knot1[0]) - p1[0] , len(knot1[1]) - p1[1]]
    knot2 = surf2.GetKnots(True)
    p2    = surf2.GetOrder()
    n2    = [len(knot2[0]) - p2[0] , len(knot2[1]) - p2[1]]

    # make sure that the sem is on top (vmax) of surf1
    tol = GetTolerance('neighbour')
    if n1[0]==n2[0] and p1[0]==p2[0]:
        n = n1[0]
        reverse = 0
        while reverse<4:
            done = True
            for i in range(n):
                if abs(surf2[i] - surf1[n1[0]*(n1[1]-1)+i]) > tol:
                    done = False
                    break
            if done:
                reverse = 5
            else:
                if reverse==0:
                    surf1.FlipParametrization(1)
                elif reverse==1:
                    surf2.FlipParametrization(1);
                elif reverse==2:
                    surf1.FlipParametrization(1);
                reverse = reverse+1
    else:
        print('Error: CombineSurfaces failed due to surfaces not having matching parametrization')
        return None


    # surfaces now organized as surf1 below surf2 with sem at vmax
    coefs =  []
    knotU = knot1[0]
    start = knot2[1][0]
    for i in range(len(knot2[1])):
        knot2[1][i] = knot2[1][i] + knot1[1][-1] - start
    knotV =  knot1[1][:-1] + knot2[1][p2[1]:]

    k = 0
    for y in range(n1[1]-1):
        for x in range(n1[0]):
            coefs.append(surf1[k])
            k = k+1
    k=0
    for y in range(n2[1]):
        for x in range(n2[0]):
            coefs.append(surf2[k])
            k = k+1

    return Surface(p1[0], p1[1], knotU, knotV, coefs, rational)

def CombineCurves(crv1, crv2):
    """Combine two curves into a single representation with internal C^0 knot. Assumes that the end of crv1 matches the start of crv2
    @param crv1: The first Curve
    @type  crv1: Curve
    @param crv2: The second Curve
    @type  crv2: Curve
    @return:     The combined Curve
    @rtype:      Curve
    """

    c1        = crv1.Clone()
    c2        = crv2.Clone()
    dim       = GetDimension()
    rational  = (len(c1[0]) == dim+1)
    if c1.GetOrder() > c2.GetOrder():
        c2.RaiseOrder(c1.GetOrder() - c2.GetOrder())
    else:
        c1.RaiseOrder(c2.GetOrder() - c1.GetOrder())

    knot1 = c1.GetKnots(True)
    p1    = c1.GetOrder()
    n1    = len(knot1) - p1

    knot2 = c2.GetKnots(True)
    p2    = c2.GetOrder()
    n2    = len(knot2) - p2


    coefs =  []
    knot =  knot1[:-1] + [k+knot1[-1] for k in knot2[p2:]] # merge 

    for i in range(n1-1):
        coefs.append(c1[i])
    for i in range(n2):
        coefs.append(c2[i])

    return Curve(p1, knot, coefs, rational)

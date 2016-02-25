__doc__ = 'Implementation of various factory methods which should have been in GoTools, but are missing or bugged'

from GoTools import *
from GoTools.SurfaceFactory import *
from GeoUtils.Knots import *
from GeoUtils.Elementary import *
from GoTools.CurveFactory import *


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


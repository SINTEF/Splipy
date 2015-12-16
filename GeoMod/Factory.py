from Curve   import *
from Surface import *
from Volume  import *

def line(a, b):
    """ Create a line between the points a and b
    @param a: start point
    @type  a: Point_like
    @param b: end point
    @type  b: Point_like
    @return : Linear spline curve from a to b
    @rtype  : Curve
    """
    return Curve(controlpoints=[a,b])

def circle(r=1):
    """ Create a circle at the origin
    @param r: circle radius
    @type  r: Float
    @return : A periodic, quadratic NURBS curve
    @rtype  : Curve
    """
    if r <= 0:
        raise ValueError('radius needs to be positive')

    w = 1.0/np.sqrt(2)
    controlpoints = [[r,0,1], [r*w,r*w,w], [0,r,1], [-r*w,r*w,w], [-r,0,1],
                     [-r*w,-r*w,w], [0,-r,1], [r*w,-r*w,w]];
    knot = np.array([0,0,0, 1,1, 2,2, 3,3, 4,4,4])/4.0*2*pi
    return Curve(3, knot, controlpoints, True, True)

def circle_segment(theta, r=1):
    """ Create a circle segment at the origin with start at (1,0)
    @param theta: circle angle in radians
    @type  theta: Float
    @param r    : circle radius
    @type  r    : Float
    @return     : A quadratic NURBS curve
    @rtype      : Curve
    """
    # error test input
    if abs(theta) > 2*pi:
        raise ValueError('theta needs to be in range [-2pi,2pi]')
    if r <= 0:
        raise ValueError('radius needs to be positive')

    # build knot vector
    knot_spans = int(theta / (2*pi/3) )
    knot = [0]
    for i in range(knot_spans+1):
        knot += [i]*2
    knot += [knot_spans] # knot vector [0,0,0,1,1,2,2,..,n,n,n]
    knot = np.array(knot) / float(knot[-1]) * theta # set parametic space to [0,theta]

    n = (knot_spans-1)*2+3 # number of control points needed
    cp = []
    t  = 0                             # current angle
    dt = float(theta) / knot_spans / 2 # angle step

    # build control points
    for i in range(n):
        w = 1 - (i%2)*(1-cos(dt))      # weights = 1 and cos(dt) every other i
        x = r*cos(t)
        y = r*sin(t)
        cp += [[x,y,w]]
        t  += dt

    return Curve(3, knot, cp, True)

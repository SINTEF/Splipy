from Curve   import *
from Surface import *
import CurveFactory

def square(size=(1,1)):
    """ Create a 2D square with lower right corner at (0,0)
    @param size: size in all directions, or (width,height)
    @type  size: Float or List of Floats
    @return    : a square
    @rtype     : Surface
    """
    result = Surface() # unit square
    result.scale(size)
    return result
    
def disc(r=1, type='radial'):
    """ Create surface representation of a circular disc with center at (0,0)
    @param r   : radius
    @type  r   : Float
    @param type: 'radial' or 'square'
    @type  type: String
    @return    : a circular disc
    @rtype     : Surface
    """
    if type is 'radial':
        c = CurveFactory.circle(r)
        cp = np.zeros((16,3))
        cp[:,-1] = 1 
        cp[1::2,:] = c.controlpoints
        return Surface(BSplineBasis(), c.basis, cp, True)
    elif type is 'square':
        w = 1/np.sqrt(2)
        cp = [ [-r*w,-r*w,1], [0,-r,w], [r*w,-r*w,1],
               [-r,     0,w], [0, 0,1], [r,     0,w],
               [-r*w, r*w,1], [0, r,w], [r*w, r*w,1]]
        basis1 = BSplineBasis(3)
        basis2 = BSplineBasis(3)
        return Surface(basis1, basis2, cp, True)
    else:
        raise ValueError('invalid type argument')

def sphere(r=1):
    """ Create sphere shell
    @param r   : radius of sphere
    @type  r   : Float
    @return    : a sphere
    @rtype     : Surface
    """
    pi = np.pi
    circle = CurveFactory.circle_segment(pi, r)
    circle.rotate(-pi/2)
    circle.rotate(pi/2, (1,0,0))    # flip up into xz-plane
    return revolve(circle)
    

def extrude(curve, h):
    """ Extrude a curve by sweeping it straight up in the z-direction
    to a given height 
    @param curve : curve to extrude
    @type  curve : Curve
    @param h     : height in z-direction
    @type  h     : Float
    @return      : an extruded surface
    @rtype       : Surface
    """
    curve.set_dimension(3) # add z-components (if not already present)
    n  = len(curve)        # number of control points of the curve
    cp = np.zeros((2*n,4))
    cp[:n,:] = curve.controlpoints # the first control points form the bottom
    curve += (0,0,h)
    cp[n:,:] = curve.controlpoints # the last control points form the top
    return Surface(curve.basis, BSplineBasis(2), cp, curve.rational)

def revolve(curve, theta=2*np.pi):
    """ Revolve a surface by sweeping a curve in a rotational fashion around
    the z-axis
    @param curve : curve to revolve
    @type  curve : Curve
    @param theta : angle in radians
    @type  theta : Float
    @return      : a revolved surface
    @rtype       : Surface
    """
    pi = np.pi
    curve.set_dimension(3) # add z-components (if not already present)
    curve.force_rational() # add weight (if not already present)
    n  = len(curve)        # number of control points of the curve
    cp = np.zeros((8*n,4))
    basis = BSplineBasis(3, [0,0,0,1,1,2,2,3,3,4,4,4], periodic=0)
    basis *= 2*pi/4        # set parametric domain to (0,2pi) in v-direction

    # loop around the circle and set control points by the traditional 9-point
    # circle curve with weights 1/sqrt(2), only here C0-periodic, so 8 points
    for i in range(8):
        if i%2 == 0:
            weight = 1
        else:
            weight = 1/sqrt(2)
        cp[i*n:(i+1)*n,:]  = curve.controlpoints
        cp[i*n:(i+1)*n,2] *= weight
        cp[i*n:(i+1)*n,3] *= weight
        curve.rotate(pi/4)
    return Surface(curve.basis, basis, cp, True)
    

def cylinder(r=1, h=1):
    """ Create cylinder shell with no top or bottom starting at the xy-plane,
    and the height increases in the z-direction
    @param r   : radius of cylinder
    @type  r   : Float
    @param h   : height in z-direction
    @type  h   : Float
    @return    : a cylinder shell
    @rtype     : Surface
    """
    return extrude(CurveFactory.circle(r), h)

def torus(minor_r=1, major_r=3):
    """ Create a torus (doughnut) by revolving a circle of size minor_r around
    the z-axis with radius major_r
    @param minor_r: the thickness of the torus, or radius in xz-plane
    @type  minor_r: Float
    @param major_r: the size of the torus, or radius in xy-plane
    @type  major_r: Float
    @return       : a periodic torus
    @rtype        : Surface
    """
    circle = CurveFactory.circle(minor_r)
    circle.rotate(pi/2, (1,0,0))    # flip up into xz-plane
    circle.translate((major_r,0,0)) # move into position to spin around z-axis
    return revolve(circle)

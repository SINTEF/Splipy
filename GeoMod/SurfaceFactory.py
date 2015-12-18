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
    if r <= 0:
        raise ValueError('radius needs to be positive')
    weight = 1/np.sqrt(2)
    cp = []
    for i in range(5):
        for j in range(8):
            theta = j*2*pi/9    # latitude  east/west
            phi   = i*pi/4-pi/2 # longitude north/south
            x = r * np.cos(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.cos(phi)
            z = r *                 np.sin(phi)
            w = 1
            if i%2==1:
                w *= weight
            if j%2==1:
                w *= weight
            cp.append([x,y,z,w])

    knot1  = np.array([0,0,0,1,1,2,2,3,3,4,4,4]) / 4.0 * 2*pi
    knot2  = np.array([0,0,0,1,1,2,2,2]        ) / 2.0 *   pi - pi/2
    basis1 = BSplineBasis(3, knot1, 0) # periodic
    basis2 = BSplineBasis(3, knot2)
    return Surface(basis1, basis2, cp, True)

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
    n = len(curve)        # number of control points of the circle
    cp = np.zeros((2*n,4))
    cp[:n,:] = curve.controlpoints # the first control points form the bottom
    curve += (0,0,h)
    cp[n:,:] = curve.controlpoints # the last control points form the top
    return Surface(curve.basis, BSplineBasis(2), cp, curve.rational)
    

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

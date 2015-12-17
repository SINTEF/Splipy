from Curve   import *
from Surface import *
import CurveFactory

def square(width=1, height=1):
    """ Create a 2D square with lower right corner at (0,0)
    @param width : width in x-direction
    @type  height: Float
    @return      : a square
    @rtype       : Surface
    """
    result = Surface() # unit square
    result.scale((width,height))
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
    


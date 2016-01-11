from Surface import *
from Volume import *
import CurveFactory
import SurfaceFactory

def cube(size=(1,1,1)):
    """ Create a volumetric cube with lower right corner at (0,0,0)
    @param size: size in all directions, or (width,depth,height)
    @type  size: Float or List of Floats
    @return    : a linear parametrized box
    @rtype     : Volume
    """
    result = Volume()
    result.scale(size)
    return result

def revolve(surf, theta=2*pi):
    """ Revolve a volume by sweeping a surface in a rotational fashion around
    the z-axis
    @param surf  : Surface to revolve
    @type  surf  : Surface
    @param theta : angle in radians
    @type  theta : Float
    @return      : a revolved surface
    @rtype       : Surface
    """
    surf = surf.clone()   # clone input surface, throw away old reference
    surf.set_dimension(3) # add z-components (if not already present)
    surf.force_rational() # add weight (if not already present)
    n  = len(surf)        # number of control points of the surface
    cp = np.zeros((8*n,4))
    basis = BSplineBasis(3, [0,0,0,1,1,2,2,3,3,4,4,4], periodic=0)
    basis *= 2*pi/4        # set parametric domain to (0,2pi) in w-direction

    # loop around the circle and set control points by the traditional 9-point
    # circle curve with weights 1/sqrt(2), only here C0-periodic, so 8 points
    for i in range(8):
        if i%2 == 0:
            weight = 1.0
        else:
            weight = 1.0/sqrt(2)
        cp[i*n:(i+1)*n,:]  = np.reshape(surf.controlpoints.transpose(1,0,2), (n,4))
        cp[i*n:(i+1)*n,2] *= weight
        cp[i*n:(i+1)*n,3] *= weight
        surf.rotate(pi/4)
    return Volume(surf.basis1, surf.basis2, basis, cp, True)

def cylinder(r=1, h=1):
    """ Create a solid cylinder with starting at the xy-plane,
    and the height increases in the z-direction
    @param r   : radius of cylinder
    @type  r   : Float
    @param h   : height in z-direction
    @type  h   : Float
    @return    : a solid cylinder
    @rtype     : Volume
    """
    shell = SurfaceFactory.cylinder(r,h)
    cp = []
    for controlpoint in shell:
        cp.append([0,0,controlpoint[2],controlpoint[3]]) # project to z-axis
    for controlpoint in shell:
        cp.append(list(controlpoint))

    return Volume(shell.basis1, shell.basis2, BSplineBasis(), cp, True)

def extrude(surf, h):
    """ Extrude a surface by sweeping it straight up in the z-direction
    to a given height 
    @param surf : surf to extrude
    @type  surf : Surface
    @param h     : height in z-direction
    @type  h     : Float
    @return      : an extruded surface
    @rtype       : Surface
    """
    surf.set_dimension(3) # add z-components (if not already present)
    cp = []
    for controlpoint in surf:
        cp.append(list(controlpoint))
    surf += (0,0,h)
    for controlpoint in surf:
        cp.append(list(controlpoint))
    surf -= (0,0,h)
    return Volume(surf.basis1, surf.basis2, BSplineBasis(2), cp, surf.rational)


def edge_surfaces(surfaces):
    """ Create the surface defined by the area between 2 or 6 input surfaces.
    In case of 6 input surfaces, then these must defined in the following
    order: bottom, top, left, right, back, front with opposing sides
    parametrized in the same directions
    @param surfaces: Two or six edge surfaces
    @type  surfaces: List of surfaces
    @return        : enclosed volume
    @rtype         : Volume
    """
    if len(surfaces)==2:
        surf1 = surfaces[0].clone()
        surf2 = surfaces[1].clone()
        Surface.make_surfaces_identical(surf1, surf2)
        (n1,n2,d) = surf1.controlpoints.shape # d = dimension + rational

        controlpoints          = np.zeros((n1, n2, 2, d))
        controlpoints[:,:,0,:] = surf1.controlpoints
        controlpoints[:,:,1,:] = surf2.controlpoints
        
        # Volume constructor orders control points in a different way, so we
        # create it from scratch here
        result = Volume()
        result.basis1       = surf1.basis1
        result.basis2       = surf1.basis2
        result.basis3       = BSplineBasis(2)
        result.dimension    = surf1.dimension
        result.rational     = surf1.rational
        result.controlpoints= controlpoints

        return result
    elif len(surfaces)==6:
        raise NotImplementedError('Should have been coons patch algorithm here. Come back later')
    else:
        raise ValueError('Requires two or six input surfaces')

from Curve   import *
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
    n = len(shell)         # number of control points of the cylinder shell
    cp = []
    for controlpoint in shell:
        cp.append([0,0,controlpoint[2],controlpoint[3]]) # project to z-axis
    for controlpoint in shell:
        cp.append(list(controlpoint))

    return Volume(shell.basis1, shell.basis2, BSplineBasis(), cp, True)


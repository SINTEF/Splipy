from Curve   import *
from Surface import *
from math    import pi,sin,cos,sqrt
import CurveFactory
import inspect
import copy
        

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
        w = 1/sqrt(2)
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

def revolve(curve, theta=2*pi):
    """ Revolve a surface by sweeping a curve in a rotational fashion around
    the z-axis
    @param curve : curve to revolve
    @type  curve : Curve
    @param theta : angle in radians
    @type  theta : Float
    @return      : a revolved surface
    @rtype       : Surface
    """
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
            weight = 1.0
        else:
            weight = 1.0/sqrt(2)
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

def edge_curves(curves):
    """ Create the surface defined by the area between 2 or 4 input curves. In
    case of 4 input curves, then these must define an ordered directional
    closed loop around the resulting surface.
    @param curves: Two or four edge curves
    @type  curves: List of curves
    @return      : intermediate surface
    @rtype       : Surface
    """
    if len(curves)==2:
        crv1 = copy.deepcopy(curves[0])
        crv2 = copy.deepcopy(curves[1])
        Curve.make_curves_identical(crv1, crv2)
        (n,d) = crv1.controlpoints.shape # d = dimension + rational

        controlpoints       = np.zeros((2*n, d))
        controlpoints[:n,:] = crv1.controlpoints
        controlpoints[n:,:] = crv2.controlpoints
        linear              = BSplineBasis(2)

        return Surface(crv1.basis, linear, controlpoints, crv1.rational)
    elif len(curves)==4:
        # coons patch (https://en.wikipedia.org/wiki/Coons_patch)
        bottom = curves[0]
        right  = curves[1]
        top    = copy.deepcopy(curves[2])
        left   = copy.deepcopy(curves[3]) # gonna change these two, so make copies
        top.flip_parametrization()
        left.flip_parametrization()
        # create linear interpolation between opposing sides
        s1     = edge_curves([bottom,top])
        s2     = edge_curves([left,right])
        s2.swap_parametrization()
        # create (linear,linear) corner parametrization
        linear = BSplineBasis(2)
        rat    = s1.rational     # using control-points from top/bottom, so need to know if these are rational
        s3     = Surface(linear,linear, [bottom[0],bottom[-1],top[0],top[-1]], rat)

        # in order to add spline surfaces, they need identical parametrization
        Surface.make_surfaces_identical(s1,s2)
        Surface.make_surfaces_identical(s1,s3)
        Surface.make_surfaces_identical(s2,s3)

        result = s1
        result.controlpoints += s2.controlpoints
        result.controlpoints -= s3.controlpoints
        return result
    else:
        raise ValueError('Requires two or four input curves')

def thicken(curve, amount):
    """ Add a thickness to a curve to generate a surface. For 2D-curves this
    is going to generate a 2D planar surface with the curve through the center,
    while 3D curves will generate a surface tube around it and with open ends.
    The resulting surface is an approximation generated by interpolating at
    the greville points. It will use the same discretization as the input
    curve. Function does not check for self-intersection.
    @param curve : generating curve
    @type  curve : Curve
    @param amount: Either a constant or variable thickness. For function input
                   one needs to specify parameters x,y,z and/or t.
    @type  amount: Float or function
    @return      : Surrounding surface
    @rtype       : Surface
    """
    # NOTES: There are several pitfalls with this function
    #  * self intersection:
    #     could be handled more gracefully, but is here ignored
    #  * choice of discretization:
    #     the offset curve is computed from the velocity (tangent) which is of
    #     one continuity less than the original curve. In particular C1
    #     quadratic curves will get very noticable C0-kinks in them. Currently
    #     this is completely ignored and we keep the high continuity of the
    #     original curve.
    #  * width given by function input
    #     could produce wild behaviour. Original discretization might not 
    #     produce a satisfactory result
    #  * 3D tube geometries: 
    #     unimplemented as of now. Would like to see the three points above
    #     resolved before this is implemented. Rough idea is to compute the
    #     acceleration and binormal vectors to the curve and sketch out a
    #     circle in the plane defined by these two vectors

    t = curve.basis.greville()
    if curve.dimension==2:
        # linear parametrization across domain
        n             = len(curve)
        left_points   = np.matrix(np.zeros((n,2)))
        right_points  = np.matrix(np.zeros((n,2)))
        linear        = BSplineBasis(2)

        x = curve.evaluate(t)            # curve at interpolation points
        v = curve.evaluate_tangent(t)    # velocity at interpolation points
        v = np.array(v)
        l = sqrt(v[:,0]**2+v[:,1]**2) # normalizing factor for velocity
        v[:,0] = v[:,0] / l
        v[:,1] = v[:,1] / l
        v = np.matrix(v)
        if inspect.isfunction(amount):
            arg_names = inspect.getargspec(amount).args
            argc      = len(arg_names)
            argv      = [0]*argc
            for i in range(n):
                # build up the list of arguments (in case not all of (x,y,t) are specified)
                for j in range(argc):
                    if arg_names[j] == 'x':
                        argv[j] = x[i,0]
                    elif arg_names[j] == 'y':
                        argv[j] = x[i,1]
                    elif arg_names[j] == 't':
                        argv[j] = t[i]
                # figure out the distane at this particular point
                dist = amount(*argv)

                # store interpolation points
                right_points[i,0] = x[i,0] - v[i,1]*dist # x at bottom
                right_points[i,1] = x[i,1] + v[i,0]*dist # y at bottom
                left_points[ i,0] = x[i,0] + v[i,1]*dist # x at top
                left_points[ i,1] = x[i,1] - v[i,0]*dist # y at top
        else:
            right_points[:,0] = x[:,0] - v[:,1]*amount # x at bottom
            right_points[:,1] = x[:,1] + v[:,0]*amount # y at bottom
            left_points[ :,0] = x[:,0] + v[:,1]*amount # x at top
            left_points[ :,1] = x[:,1] - v[:,0]*amount # y at top
        # perform interpolation on each side
        right = CurveFactory.interpolate(right_points, curve.basis)
        left  = CurveFactory.interpolate(left_points,  curve.basis)

        # draw the linear surface in between
        controlpoints = np.zeros((2*n,2))
        controlpoints[:n,:] = right.controlpoints
        controlpoints[n:,:] = left.controlpoints
        return Surface(curve.basis, linear, controlpoints)

    else: # dimension=3, we will create a surrounding tube
        raise NotImplementedError('Currently only 2D supported. See comments in source code')



__doc__ = 'Implementation of various interpolation schemes.'

import GoTools
import numpy as np
from operator import itemgetter

FREE=1
NATURAL=2
HERMITE=3
PERIODIC=4
TANGENT=5
TANGENTNATURAL=6
UNIFORM=8

def getBasis(t, knot, d=0, fromRight=True):
    """Get all basis functions evaluated in a given set of points
    @param t:         The parametric coordinates to perform evaluation
    @type  t:         List of floats (in non-decreasing order)
    @param knot:      Open knot vector corresponding to the basis functions
    @type  knot:      List of floats (in non-decreasing order)
    @param d:         Number of derivatives
    @type  d:         Integer
    @param fromRight: True if evaluation should be done in the limit from above
    @type  fromRight: Boolean
    @return:          Matrix of all basis functions in all points
    @rtype:           Numpy.Matrix
    """
    TOL = 1e-10 # knot vector snap tolerance
    p = 1       # knot vector order
    while knot[p]==knot[0]:
        p += 1
    n = len(knot)-p                   # number of basis functions
    N = np.matrix(np.zeros((len(t),n)))
    for i in range(len(t)):           # for all evaluation points t
        if t[i]<knot[0] or t[i]>knot[-1]:
            continue
        M = np.zeros(p+1)             # temp storage to keep all the function evaluations
        mu = p                        # index of last non-zero basis function
        if abs(t[i]-knot[-1]) < TOL:  # special case the endpoint
            mu = n
        else:
            if fromRight:
                while t[i]-knot[mu] > -TOL:
                    mu += 1
            else:
                while t[i] > knot[mu]:
                    mu += 1
        M[-2] = 1 # the last entry is a dummy-zero which is never used
        for q in range(1,p):
            for j in range(p-q-1, p):
                k = mu-p+j   # 'i'-index in global knot vector (ref Hughes book pg.21)
                if abs(knot[k+q]-knot[k])>TOL:
                    if q<p-d:                                               # in case of normal evaluation
                        M[j] = M[j]*float(t[i]-knot[k])/(knot[k+q]-knot[k])
                    else:                                                   # in case of derivative evaluation
                        M[j] = M[j]*float(q)/(knot[k+q]-knot[k])
                else:                                                       # in case of multiplicative knot
                    M[j] = 0
                if abs(knot[k+q+1]-knot[k+1])>TOL:                          # and the same for the second term in the sum
                    if q<p-d:
                        M[j] = M[j] + M[j+1]*float(knot[k+q+1]-t[i])/(knot[k+q+1]-knot[k+1])
                    else:
                        M[j] = M[j] - M[j+1]*float(q)/(knot[k+q+1]-knot[k+1])
        N[i,(mu-p):mu] = M[0:-1]
    return N

def InterpolateCurve(x, t, knot):
    """Interpolate a curve at a given set of points. User is responsible that the input problem is well posed
    @param x:    The physical coordinates of the points to interpolate. Size nxd, where d is the dimension and n is the number of points
    @type  x:    Numpy.Matrix
    @param t:    The parametric coordinates of the points to interpolate. Length n
    @type  t:    List of floats
    @param knot: Open knot vector to use for interpolation. Size n+p, where p is the order to use for interpolation
    @type  knot: List of floats
    @return:     Interpolating curve
    @rtype:      Curve
    """
    p = 1 # knot vector order
    while knot[p]==knot[0]:
        p += 1
    N = getBasis(t,knot)
    C = np.linalg.solve(N,x)
    return GoTools.Curve(p, knot, C.tolist(), False);

def ApproximateCurve(x, t, knot):
    """Approximate a curve of m discrete points using a spline of n control points, where n<m
    @param x:    The physical coordinates of the points to approximate. Size mxd, where d is the dimension and m is the number of points
    @type  x:    Numpy.Matrix
    @param t:    The parametric coordinates of the points above. Length m
    @type  t:    List of floats
    @param knot: Open knot vector to use for approximation. Size n+p, where p is the spline order and n is the number of control points
    @type knot:  List of floats
    @return:     Approximating curve
    @rtype:      Curve
    """
    p = 1 # knot vector order
    while knot[p]==knot[0]:
        p += 1
    N = getBasis(t,knot)
    C = np.linalg.solve(N.transpose()*N, N.transpose()*x)
    return GoTools.Curve(p, knot, C.tolist(), False);


def InterpolateSurface(x, u, v, knotU, knotV):
    """Interpolate a surface at a given set of points. User is responsible that the input problem is well posed
    @param x:     The physical coordinates of the points to interpolate. Size (NxM)xD, where D is the dimension and NxM is the number of points
    @type  x:     List of list of floats
    @param u:     The parametric coordinates in the first direction of the points to interpolate. Length N
    @type  u:     List of floats
    @param v:     The parametric coordinates in the second direction of the points to interpolate. Length M
    @type  v:     List of floats
    @param knotU: Open knot vector to use for interpolation. Size N+p, where p is the order in the first direction
    @type  knotU: List of floats
    @param knotV: Open knot vector to use for interpolation. Size M+q, where q is the order in the second direction
    @type  knotV: List of floats
    @return:      Interpolating surface
    @rtype:       Surface
    """
    Nu = getBasis(u,knotU)
    Nv = getBasis(v,knotV)
    invNu = np.linalg.inv(Nu)
    invNv = np.linalg.inv(Nv)
    p = 1 # knot vector order
    q = 1
    while knotU[p]==knotU[0]:
        p += 1
    while knotV[q]==knotV[0]:
        q += 1
    n   = len(Nu)
    m   = len(Nv)
    dim = len(x[0])
    controlpts = np.zeros((n*m,dim))
    for d in range(dim):
        X = np.matrix(np.zeros((n,m))) # rewrap input x-values to matrix format
        k = 0
        for j in range(m):
            for i in range(n):
                X[i,j] = x[k][d]
                k+=1
        C = invNu*X*(invNv.transpose()) # compute control points
        k=0
        for j in range(m):
            for i in range(n):
                controlpts[k,d] = C[i,j] # rewrap output c-values to list format
                k+=1
    return GoTools.Surface(p,q, knotU, knotV, controlpts.tolist(), False)

def ApproximateSurface(x, u, v, knotU, knotV):
    """Approximate a surface of nxm discrete points using a spline of NxM control points, where n>N and m>M
    @param x:     The physical coordinates of the points to interpolate. Size (nxm)xD, where D is the dimension and NxM is the number of points
    @type  x:     List of list of floats
    @param u:     The parametric coordinates in the first direction of the points to interpolate. Length N
    @type  u:     List of floats
    @param v:     The parametric coordinates in the second direction of the points to interpolate. Length M
    @type  v:     List of floats
    @param knotU: Open knot vector to use for interpolation. Size N+p, where p is the order in the first direction
    @type  knotU: List of floats
    @param knotV: Open knot vector to use for interpolation. Size M+q, where q is the order in the second direction
    @type  knotV: List of floats
    @return:      Interpolating surface
    @rtype:       Surface
    """
    Nu = getBasis(u,knotU)
    Nv = getBasis(v,knotV)
    NuT = Nu.transpose()
    NvT = Nv.transpose()
    invNu = np.linalg.inv(NuT*Nu) # symmetric matrix
    invNv = np.linalg.inv(NvT*Nv)
    p = 1 # knot vector order
    q = 1
    while knotU[p]==knotU[0]:
        p += 1
    while knotV[q]==knotV[0]:
        q += 1
    N   = len(u)
    M   = len(v)
    n   = len(knotU) -p
    m   = len(knotV) -q
    dim = len(x[0])
    controlpts = np.zeros((n*m,dim))
    for d in range(dim):
        X = np.matrix(np.zeros((N,M))) # rewrap input x-values to matrix format
        k = 0
        for j in range(M):
            for i in range(N):
                X[i,j] = x[k][d]
                k+=1
        C = invNu * NuT * X * Nv * invNv # compute control points
        k=0
        for j in range(m):
            for i in range(n):
                controlpts[k,d] = C[i,j] # rewrap output c-values to list format
                k+=1
    return GoTools.Surface(p,q, knotU, knotV, controlpts.tolist(), False)

def LinearCurve(x=[], y=[], z=[], pts=[], t=[]):
    """Linearly interpolate a list of points.
    @param x:   The x-coordinate of the points to interpolate
    @type  x:   List of floats
    @param y:   y-coordinates
    @type  y:   List of floats
    @param z:   z-coordinates
    @type  z:   List of floats
    @param pts: Coordinates of the points to interpolate (supersedes x, y and z)
    @type  pts: List of Point
    @param t:   The parametric coordinates of the points above. Length n. If not given,
                arclength parametrization is used.
    @type  t:   List of floats
    @return:    Linear interpolated curve in 3 dimensions
    @rtype:     Curve
    """
    if pts:
        x, y, z = [map(itemgetter(i), pts) for i in xrange(3)]

    pts = np.matrix(np.zeros((3,len(x))))
    pts[0,:] = x
    if(len(y) != 0):
        pts[1,:] = y
    if(len(z) != 0):
        pts[2,:] = z

    if not t:
        n = len(x)
        t   = [0]
        for i in range(1,n):
            x0 = pts[:,i-1]
            x1 = pts[:,i]
            dist = np.linalg.norm(x1-x0)  # eucledian distance between points i and (i-1)
            t.append(t[i-1]+dist)

    knot = [t[0]] + t + [t[-1]]

    return InterpolateCurve(pts.transpose(), t, knot)

def CubicCurve(x=[], y=[], z=[], boundary=FREE, derX=[], derY=[], derZ=[], pts=[], t=[], der=[]):
    """Cubic spline interpolation of a list of points.
    @param x:        The x-coordinate of the points to interpolate
    @type  x:        List of floats (n points)
    @param y:        y-coordinates
    @type  y:        List of floats
    @param z:        z-coordinates
    @type  z:        List of floats
    @param boundary: Boundary conditions
    @type  boundary: FREE, NATURAL, HERMITE, PERIODIC, TANGENT or TANGENTNATURAL
    @param derX:     In case of Hermite or Tangent boundary conditions, one must supply tangent information
    @type  derX:     List of floats (two points for Tangent, n points for Hermite)
    @param derY:     Y-component of tangent information
    @type  derY:     List of floats
    @param derZ:     Z-component of tangent information
    @type  derZ:     List of floats
    @param pts:      Coordinates of the points to interpolate (supersedes x, y and z)
    @type  pts:      List of Point
    @param t:        The parametric coordinates of the points above. Length n. If not given,
                     arclength parametrization is used.
    @type  t:        List of floats
    @param der:      Derivatives (supersedes derX, derY and derZ)
    @type  der:      List of Point
    @return:         Cubic interpolated curve in 3 dimensions
    @rtype:          Curve
    """
    if pts:
        x, y, z = [map(itemgetter(i), pts) for i in xrange(3)]
    if der:
        derX, derY, derZ = [map(itemgetter(i), der) for i in xrange(3)]

    pts = np.matrix(np.zeros((3,len(x))))
    pts[0,:] = x
    if(len(y) != 0):
        pts[1,:] = y
    if(len(z) != 0):
        pts[2,:] = z

    n = len(x)

    # create a knot vector based on arclength discretization
    if not t:
        t = [0]
        for i in range(1,n):
            x0 = pts[:,i-1]
            x1 = pts[:,i]
            dist = np.linalg.norm(x1-x0)  # eucledian distance between points i and (i-1)
            t.append(t[i-1]+dist)

    # modify knot vector for chosen boundary conditions
    knot = [t[0]]*3 + t + [t[-1]]*3
    if boundary == FREE:
        del knot[-5]
        del knot[4]
    elif boundary == HERMITE:
        knot = sorted(knot + t[1:-1])

    # fetch the main system matrix, will be modified to take into account boundary conditions
    N = getBasis(t,knot)  # matrix
    pts = pts.transpose() # right-hand-side vector

    # add boundary conditions
    if boundary==TANGENT or boundary==HERMITE or boundary==TANGENTNATURAL:
        if boundary == TANGENT:
            dn  = getBasis([t[0], t[-1]], knot, 1)
            N   = np.resize(N,   (N.shape[0]+2,   N.shape[1]))
            pts = np.resize(pts, (pts.shape[0]+2, pts.shape[1]))
        elif boundary == TANGENTNATURAL:
            dn  = getBasis([t[0]], knot, 1)
            N   = np.resize(N,   (N.shape[0]+1,   N.shape[1]))
            pts = np.resize(pts, (pts.shape[0]+1, pts.shape[1]))
        elif boundary == HERMITE:
            dn  = getBasis(t, knot, 1)
            N   = np.resize(N,   (N.shape[0]+n,   N.shape[1]))
            pts = np.resize(pts, (pts.shape[0]+n, pts.shape[1]))
        N[n:,:]   = dn
        pts[n:,0] = derX            # derivative given at start/end (TANGENT)
        if len(derY) != 0:          # or at all internal knots (HERMITE)
            pts[n:,1]  = derY
        else:
            pts[n:,1]  = 0
        if len(derZ) != 0:
            pts[n:,2]  = derZ
        else:
            pts[n:,2]  = 0

    if boundary == PERIODIC:
        dn   = getBasis([t[0], t[-1]], knot, 1)
        ddn  = getBasis([t[0], t[-1]], knot, 2)
        N    = np.resize(N,   (N.shape[0]+2,   N.shape[1]))
        pts  = np.resize(pts, (pts.shape[0]+2, pts.shape[1]))
        N[-2,:]    =  dn[0,:] -  dn[1,:] # first derivative matching at start/end
        N[-1,:]    = ddn[0,:] - ddn[1,:] # second derivative matching
        pts[-2:,:] = 0

    if boundary == NATURAL or boundary == TANGENTNATURAL:
        if boundary == NATURAL:
            ddn  = getBasis([t[0], t[-1]], knot, 2)
            new  = 2
        elif boundary == TANGENTNATURAL:
            ddn  = getBasis([t[-1]], knot, 2)
            new  = 1
        N    = np.resize(N,   (N.shape[0]+new,   N.shape[1]))
        pts  = np.resize(pts, (pts.shape[0]+new, pts.shape[1]))
        N[-new:,:] = ddn
        pts[-new:,:] = 0

    # solve system to get controlpoints
    C = np.linalg.solve(N,pts)

    # wrap it all into a GoTools curve and return
    return GoTools.Curve(4, knot, C.tolist(), False);

def UniformCubicCurve(x=[], y=[], z=[], pts=[]):
    """Cubic spline interpolation a list of points by uniform parametrization
    @param x:   The x-coordinate of the points to interpolate
    @type  x:   List of floats
    @param y:   y-coordinates
    @type  y:   List of floats
    @param z:   z-coordinates
    @type  z:   List of floats
    @param pts: Coordinates of the points to interpolate (supersedes x, y and z)
    @type pts:  List of Point
    @return:    Cubic interpolated curve in 3 dimensions
    @rtype:     Curve
    """
    if pts:
        x, y, z = [map(itemgetter(i), pts) for i in xrange(3)]

    pts = np.matrix(np.zeros((3,len(x))))
    pts[0,:] = x
    if(len(y) != 0):
        pts[1,:] = y
    if(len(z) != 0):
        pts[2,:] = z

    n = len(x)
    t = range(n)
    knot = [t[0]]*4 + t[2:-2] + [t[-1]]*4
    return InterpolateCurve(pts.transpose(), t, knot)

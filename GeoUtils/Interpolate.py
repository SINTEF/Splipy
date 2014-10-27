__doc__ = 'Implementation of various interpolation schemes.'

import GoTools
import GoTools.CurveFactory
import numpy as np

def getBasis(t, knot):
    """Get all basis functions evaluated in a given set of points
    @param t: The parametric coordinates to perform evaluation
    @type t:  List of floats
    @param t: Open knot vector corresponding to the basis functions
    @type t:  List of floats
    @return:  Matrix of all basis functions in all points
    @rtype:   Numpy.Matrix
    """
    p = 1 # knot vector order
    while knot[p]==knot[0]:
        p += 1
    n = len(knot)-p # number of basis functions
    N = np.matrix(np.zeros((len(t),n)))
    for i in range(len(t)):  # for all evaluation points t
        if t[i]<knot[0] or t[i]>knot[-1]:
            continue
        M = np.zeros(p+1)    # temp storage to keep all the function evaluations
        mu = p           # index of last non-zero basis function
        if t[i]==knot[-1]:   # special case the endpoint
            mu = n
        else:
            while knot[mu]<t[i]:
                mu += 1
        M[-2] = 1 # the last entry is a dummy-zero which is never used
        for q in range(1,p):
            for j in range(p-q-1, p):
                k = mu-p+j # 'i'-index in global knot vector (ref Hughes book pg.21)
                if not knot[k+q]==knot[k]:
                    M[j] = M[j]*float(t[i]-knot[k])/(knot[k+q]-knot[k])
                else:
                    M[j] = 0
                if not knot[k+q+1]==knot[k+1]:
                    M[j] = M[j] + M[j+1]*float(knot[k+q+1]-t[i])/(knot[k+q+1]-knot[k+1])
        N[i,(mu-p):mu] = M[0:-1]
    return N

def InterpolateCurve(x, t, knot):
    """Interpolate a curve at a given set of points. User is responsible that the input problem is well posed
    @param x:    The physical coordinates of the points to interpolate. Size nxd, where d is the dimension and n is the number of points
    @type x:     Numpy.Matrix
    @param t:    The parametric coordinates of the points to interpolate. Length n
    @type t:     List of floats
    @param knot: Open knot vector to use for interpolation. Size n+p, where p is the order to use for interpolation
    @type knot:  List of floats
    @return:     Control points corresponding to the interpolating curve
    @rtype:      Numpy.Matrix
    """
    N = getBasis(t,knot)
    C = np.linalg.solve(N,x)
    return C

def ApproximateCurve(x, t, knot):
    """Approximate a curve of m discrete points using a spline of n control points, where n<m
    @param x:    The physical coordinates of the points to approximate. Size mxd, where d is the dimension and m is the number of points
    @type x:     Numpy.Matrix
    @param t:    The parametric coordinates of the points above. Length m
    @type t:     List of floats
    @param knot: Open knot vector to use for approximation. Size n+p, where p is the spline order and n is the number of control points
    @type knot:  List of floats
    @return:     Control points corresponding to the approximation curve
    @rtype:      Numpy.Matrix
    """
    N = getBasis(t,knot)
    C = np.linalg.solve(N.transpose()*N, N.transpose()*x)
    return C

def InterpolateSurface(x, u, v, knotU, knotV):
    """Interpolate a surface at a given set of points. User is responsible that the input problem is well posed
    @param x:     The physical coordinates of the points to interpolate. Size (NxM)xD, where D is the dimension and NxM is the number of points
    @type x:      Numpy.Matrix
    @param u:     The parametric coordinates in the first direction of the points to interpolate. Length N
    @type u:      List of floats
    @param v:     The parametric coordinates in the second direction of the points to interpolate. Length M
    @type v:      List of floats
    @param knotU: Open knot vector to use for interpolation. Size N+p, where p is the order in the first direction
    @type knotU:  List of floats
    @param knotV: Open knot vector to use for interpolation. Size M+q, where q is the order in the second direction
    @type knotV:  List of floats
    @return:      Control points corresponding to the interpolating curve
    @rtype:       Numpy.Matrix
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

def Linear(x,y=[],z=[]):
    """Linear interpolate a list of points (arclength parametrization)
    @param x: The x-coordinate of the points to interpolate
    @type x:  List of floats
    @param y: y-coordinates
    @type y:  List of floats
    @param z: z-coordinates
    @type z:  List of floats
    @return:  Linear interpolated curve in 3 dimensions
    @rtype:   Curve
    """

    pts = np.matrix(np.zeros((3,len(x))))
    pts[0,:] = x
    if(len(y) != 0):
        pts[1,:] = y
    if(len(z) != 0):
        pts[2,:] = z

    n = len(x)
    t   = [0]
    for i in range(1,n):
        x0 = pts[:,i-1]
        x1 = pts[:,i]
        dist = np.linalg.norm(x1-x0)  # eucledian distance between points i and (i-1)
        t.append(t[i-1]+dist)
    knot = [0]+t+[t[-1]]
    C = InterpolateCurve(pts.transpose(), t, knot)
    return GoTools.Curve(2, knot, C.tolist(), False);

def Cubic(x, y=[], z=[]):
    """Cubic spline interpolation a list of points by arclength parametrization
    @param x: The x-coordinate of the points to interpolate
    @type x: List of floats
    @param y: y-coordinates
    @type y: List of floats
    @param z: z-coordinates
    @type z: List of floats
    @return: Cubic interpolated curve in 3 dimensions
    @rtype: Curve
    """
    pts = np.matrix(np.zeros((3,len(x))))
    pts[0,:] = x
    if(len(y) != 0):
        pts[1,:] = y
    if(len(z) != 0):
        pts[2,:] = z

    n = len(x)
    t   = [0]
    for i in range(1,n):
        x0 = pts[:,i-1]
        x1 = pts[:,i]
        dist = np.linalg.norm(x1-x0)  # eucledian distance between points i and (i-1)
        t.append(t[i-1]+dist)
    knot = [t[0]]*4 + t[2:-2] + [t[-1]]*4
    C = InterpolateCurve(pts.transpose(), t, knot)
    return GoTools.Curve(4, knot, C.tolist(), False);

def UniformCubic(x, y=[], z=[]):
    """Cubic spline interpolation a list of points by uniform parametrization
    @param x: The x-coordinate of the points to interpolate
    @type x:  List of floats
    @param y: y-coordinates
    @type y:  List of floats
    @param z: z-coordinates
    @type z:  List of floats
    @return:  Cubic interpolated curve in 3 dimensions
    @rtype:   Curve
    """
    pts = np.matrix(np.zeros((3,len(x))))
    pts[0,:] = x
    if(len(y) != 0):
        pts[1,:] = y
    if(len(z) != 0):
        pts[2,:] = z

    n = len(x)
    t = range(n)
    knot = [t[0]]*4 + t[2:-2] + [t[-1]]*4
    C = InterpolateCurve(pts.transpose(), t, knot)
    return GoTools.Curve(4, knot, C.tolist(), False);


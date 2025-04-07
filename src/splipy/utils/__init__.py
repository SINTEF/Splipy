# -*- coding: utf-8 -*-

from itertools import combinations, product
from math import atan2, sqrt
import numpy as np

try:
    from collections.abc import Sized
except ImportError:
    from collections import Sized

def is_right_hand(patch, tol=1e-3):
    param = tuple((a+b)/2 for a,b in zip(patch.start(), patch.end()))

    if patch.dimension == patch.pardim == 3:
        du = patch.derivative(*param, d=(1,0,0))
        dv = patch.derivative(*param, d=(0,1,0))
        dw = patch.derivative(*param, d=(0,0,1))

        # Normalize
        du = du / np.linalg.norm(du)
        dv = dv / np.linalg.norm(dv)
        dw = dw / np.linalg.norm(dw)

        # Compare cross product
        return np.dot(dw, np.cross(du, dv)) >= tol

    if patch.dimension == patch.pardim == 2:
        du = patch.derivative(*param, d=(1,0))
        dv = patch.derivative(*param, d=(0,1))

        # Normalize
        du = du / np.linalg.norm(du)
        dv = dv / np.linalg.norm(dv)

        return np.cross(du, dv) >= tol

    raise ValueError("Right-handedness only defined for 2D or 3D patches in 2D or 3D space, respectively")

def rotation_matrix(theta, axis):
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2)
    b, c, d = -axis*np.sin(theta / 2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)],
                      [2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c]])

def sections(src_dim, tgt_dim):
    """Generate all boundary sections from a source dimension to a target
    dimension. For example, `sections(3,1)` generates all edges on a volume.

    The return values are lists of length `src_dim` with each element either 0,
    --1 or `None`, which are suitable for passing to
    :func:`splipy.SplineObject.section`.
    """
    # Enumerate all combinations of fixed directions
    nfixed = src_dim - tgt_dim
    for fixed in combinations(range(src_dim), r=nfixed):
        # Enumerate all {0,-1}^n over the fixed directions
        for indices in product([0, -1], repeat=nfixed):
            args = [None] * src_dim
            for f, i in zip(fixed, indices[::-1]):
                args[f] = i
            yield args

def section_from_index(src_dim, tgt_dim, i):
    """Return the i'th section from a source dimension to a target dimension.

    See :func:`splipy.Utils.sections` for more information.
    """
    for j, s in enumerate(sections(src_dim, tgt_dim)):
        if i == j:
            return s

def section_to_index(section):
    """Return the index corresponding to a section."""
    src_dim = len(section)
    tgt_dim = sum(1 for s in section if s is None)
    for i, t in enumerate(sections(src_dim, tgt_dim)):
        if tuple(section) == tuple(t):
            return i

def check_section(*args, **kwargs):
    """check_section(u, v, ...)

    Parse arguments and return a section spec.

    The keyword argument `pardim` *must* be provided. The return value is a
    section as described in :func:`splipy.Utils.sections`.
    """
    pardim = kwargs['pardim']
    args = list(args)
    while len(args) < pardim:
        args.append(None)
    for k in set(kwargs.keys()) & set('uvw'):
        index = 'uvw'.index(k)
        args[index] = kwargs[k]
    return args

def check_direction(direction, pardim):
    if   direction in {0, 'u', 'U'} and 0 < pardim:
        return 0
    elif direction in {1, 'v', 'V'} and 1 < pardim:
        return 1
    elif direction in {2, 'w', 'W'} and 2 < pardim:
        return 2
    raise ValueError('Invalid direction')

def ensure_flatlist(x):
    """Flattens a multi-list x to a single index list."""
    if isinstance(x[0], Sized):
        return x[0]
    return x

def is_singleton(x):
    """Checks if x is list-like."""
    return not isinstance(x, Sized)

def ensure_listlike(x, dups=1):
    """Wraps x in a list if it's not list-like."""
    try:
        while len(x) < dups:
            x = list(x)
            x.append(x[-1])
        return x
    except TypeError:
        return [x] * dups
    except IndexError:
        return []

def rotate_local_x_axis(xaxis=(1,0,0), normal=(0,0,1)):
    # rotate xaxis vector back to reference domain (r=1, around origin)
    theta = atan2(normal[1], normal[0])
    phi   = atan2(sqrt(normal[0]**2+normal[1]**2), normal[2])
    R1 = rotation_matrix(-theta, (0,0,1))
    R2 = rotation_matrix(-phi,   (0,1,0))
    if len(xaxis) != 3: # typically 2D geometries
        xaxis = [xaxis[0], xaxis[1], 0]
    xaxis = np.array([xaxis])
    xaxis = xaxis.dot(R1).dot(R2)
    # if xaxis is orthogonal to normal, then xaxis[2]==0 now. If not then
    # treating it as such is the closest projection, which makes perfect sense
    return atan2(xaxis[0,1], xaxis[0,0])

def flip_and_move_plane_geometry(obj, center=(0,0,0), normal=(0,0,1)):
    """re-orients a planar geometry by moving it to a different location and
    tilting it"""
    # don't touch it if not needed. translate or scale operations may force
    # object into 3D space
    if not np.allclose(normal, np.array([0,0,1])):
        theta = atan2(normal[1], normal[0])
        phi   = atan2(sqrt(normal[0]**2+normal[1]**2), normal[2])
        obj.rotate(phi,   (0,1,0))
        obj.rotate(theta, (0,0,1))
    if not np.allclose(center, 0):
        obj.translate(center)
    return obj

def reshape(cps, newshape, order='C', ncomps=None):
    """Like numpy's reshape, but preserves control points of several dimensions
    that are stored contiguously.

    The return value has shape (*newshape, ncomps), where ncomps is the number
    of components per control point, as inferred by the size of `cps` and the
    desired shape.

    The `order` argument ('C' or 'F') determines the order in which control
    points are read, but does *not* affect the order in which each component of
    a control point is read.
    """
    npts = np.prod(newshape)
    if ncomps is None:
        try:
            ncomps = cps.size // npts
        except AttributeError:
            ncomps = len(cps) // npts

    if order == 'C':
        shape = list(newshape) + [ncomps]
    elif order == 'F':
        shape = list(newshape[::-1]) + [ncomps]
    cps = np.reshape(cps, shape)
    if order == 'F':
        spec = list(range(len(newshape)))[::-1] + [len(newshape)]
        cps = cps.transpose(spec)
    return cps

def uniquify(iterator):
    """Iterates over all elements in `iterator`, removing duplicates."""
    seen = set()
    for i in iterator:
        if i in seen:
            continue
        seen.add(i)
        yield i

def raise_order_1D(n, k, T, P, m, periodic):
    """ Implementation of method in "Efficient Degree Elevation and Knot Insertion
        for B-spline Curves using Derivatives" by Qi-Xing Huang a Shi-Min Hu, Ralph R Martin. Only the case of open knot vector is fully implemented
    :param int n: (n+1) is the number of initial basis functions
    :param int k: spline order
    :param T: knot vector
    :param P: weighted NURBS coefficients
    :param int m: number of degree elevations
    :param int periodic: Number of continuous derivatives at start and end. -1 is not periodic, 0 is continuous, etc.
    :return Q: new control points
    """

    u = np.unique(T[k-1:-k+1])
    S = u.size - 1
    d = P.shape[0] # dimension of spline

    # Find multiplicity of the knot vector T
    b = BSplineBasis(k, T)
    z = [k-1-b.continuity(t0) for t0 in b.knot_spans()]
    
    # Step 1: Find Pt_i^j
    Pt = np.zeros((d,n+1,k))
    Pt[:,:,0] = P
    Pt = np.concatenate((Pt,Pt[:,0:periodic+1,:]),axis=1) 
    n += periodic+1

    beta = np.cumsum(z[1:-1],dtype=int)
    beta = np.insert(beta,0,0) # include the empty sum (=0)
    for l in range(1,k):
        for i in range(0,n+1-l):
            if T[i+l] < T[i+k]:
                Pt[:,i,l] = (Pt[:,i+1,l-1] - Pt[:,i,l-1])/(T[i+k]-T[i+l])

    # Step 2: Create new knot vector Tb 
    nb = n + S*m
    Tb = np.zeros(nb+m+k+1)
    Tb[:k-1] = T[:k-1]
    Tb[-k+1:] = T[-k+1:]
    j = k-1 
    for i in range(0,len(z)):
        Tb[j:j+z[i]+m] = u[i]
        j = j+z[i]+m
    
    # Step 3: Find boundary values of Qt_i^j
    Qt = np.zeros((d,nb+1,k))
    l_arr = np.array(range(1,k))
    alpha = np.cumprod((k-l_arr)/(k+m-l_arr))
    alpha = np.insert(alpha,0,1) # include the empty product (=1)
    indices = range(0,k)
    Qt[:,0,indices] = np.multiply(Pt[:,0,indices],np.reshape(alpha[indices],(1,1,k))) # (21)
    for p in range(0,S):
        indices = range(k-z[p],k)
        Qt[:,beta[p]+p*m,indices] = np.multiply(Pt[:,beta[p],indices],np.reshape(alpha[indices],(1,1,z[p]))) # (22)
    
    for p in range(0,S):
        idx = beta[p]+p*m
        Qt[:,idx+1:m+idx+1,k-1] = np.repeat(Qt[:,idx:idx+1,k-1],m,axis=1) # (23)
    
    # Step 4: Find remaining values of Qt_i^j
    for j in range(k-1,0,-1):
        for i in range(0,nb):
            if Tb[i+k+m] > Tb[i+j]:
                Qt[:,i+1,j-1] = Qt[:,i,j-1] + (Tb[i+k+m] - Tb[i+j])*Qt[:,i,j] #(20) with Qt replacing Pt

    return Qt[:,:,0]

__all__ = [
    'nutils', 'refinement', 'image', 'NACA', 'curve', 'smooth',
    'rotation_matrix', 'sections', 'section_from_index', 'section_to_index',
    'check_section', 'check_direction', 'ensure_flatlist', 'is_singleton',
    'ensure_listlike', 'rotate_local_x_axis', 'flip_and_move_plane_geometry',
    'reshape','raise_order_1D'
]

import numpy as np
import copy
from itertools import product
from GeoMod import BSplineBasis
from GeoMod.Utils import ensure_listlike

__all__ = ['ControlPointOperations']


def get_rotation_matrix(theta, axis):
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2)
    b, c, d = -axis*np.sin(theta / 2)
    return np.matrix([[a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)],
                      [2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c]])


class ControlPointOperations(object):

    def __init__(self, bases, controlpoints, rational):
        bases = [(b if b else BSplineBasis()) for b in bases]
        self.bases = bases
        if controlpoints is None:
            # We want to generate tuples in column-major format, hence the [::-1]
            controlpoints = [c[::-1] for c in product(*(b.greville() for b in bases[::-1]))]
            if len(controlpoints[0]) == 1:
                controlpoints = [tuple(list(c) + [0.0]) for c in controlpoints]
        controlpoints = np.array(controlpoints)
        self.dimension = controlpoints.shape[-1] - rational
        self.rational = rational

        cp_shape = tuple([b.num_functions() for b in bases[::-1]] + [self.dimension + rational])
        controlpoints = np.reshape(controlpoints, cp_shape)

        # Compensate for numpy's row-major ordering
        spec = tuple(list(range(len(bases)))[::-1] + [len(bases)])
        self.controlpoints = controlpoints.transpose(spec)

    def _validate_domain(self, *params):
        """Check whether the given evaluation parameters are valid."""
        for b, p in zip(self.bases, params):
            if b.periodic < 0:
                if min(p) < b.start() or b.end() < max(p):
                    raise ValueError('Evaluation outside parametric domain')

    def evaluate(self, *params):
        """Evaluate the object at given parametric values
        @param params: Parametric coordinate point(s)
        @type  params: Float or list of Floats
        @return : Geometry coordinates. Matrix X(i1,i2,...,j) of component x(j) evaluated at t(i1,i2,...)
        @rtype  : numpy.array
        """
        params = [ensure_listlike(p) for p in params]

        self._validate_domain(*params)

        Ns = [b.evaluate(p) for b, p in zip(self.bases, params)]
        idx = len(self.bases) - 1
        result = self.controlpoints
        for N in Ns[::-1]:
            result = np.tensordot(N, result, axes=(1, idx))

        if self.rational:
            for i in range(self.dimension):
                result[..., i] /= result[..., -1]
            result = np.delete(result, self.dimension, -1)

        if all(s == 1 for s in result.shape[:-1]):
            result = result.reshape(self.dimension)

        return result

    def evaluate_derivative(self, *params, **kwargs):
        """Evaluate the derivative of the object at given parametric values

        This function accepts n parametric coordinate point(s) (which are
        floats or list of floats), and an optional order of the derivative
        to compute.

        @param params: Parametric coordinate point(s)
        @type  params: Float or list of Floats
        @param d: Derivative order (default (1,1,...))
        @type  d: Int or Tuple
        @return: Derivatives. Matrix X(i1,i2,...,j) of component x(j)
            differentiationed d(j) times, evaluated at t(i1,i2,....)

        Tangent of curve at single point
        curve.evaluate_derivative(1.0)

        Double derivative of curve at single point
        curve.evaluate_derivative(1.0, d=2)

        Third derivative of curve at several points
        curve.evaluate_derivative([0.0, 1.0, 2.0], d=3)

        Tangents of surface
        surface.evaluate_derivative(0.5, 0.7, d=(1,0))
        surface.evaluate_derivative(0.5, 0.7, d=(0,1))

        Cross-derivative of surface
        surface.evaluate_derivative(0.5, 0.7, d=(1,1))
        """
        params = [ensure_listlike(p) for p in params]

        derivs = kwargs.get('d', [1] * len(self.bases))
        derivs = ensure_listlike(derivs)

        self._validate_domain(*params)

        dNs = [b.evaluate(p, d) for b, p, d in zip(self.bases, params, derivs)]
        idx = len(self.bases) - 1
        result = self.controlpoints
        for dN in dNs[::-1]:
            result = np.tensordot(dN, result, axes=(1, idx))

        if self.rational:
            if sum(derivs) > 1:
                raise RuntimeError('Rational derivative not implemented for order %i' % sum(derivs))
            Ns = [b.evaluate(p) for b, p in zip(self.bases, params)]
            non_derivative = self.controlpoints
            for N in Ns[::-1]:
                non_derivative = np.tensordot(N, non_derivative, axes=(1, idx))
            W = non_derivative[..., -1]  # W(u,v)
            Wd = result[..., -1]         # dW(u,v)/du or /dv
            for i in range(self.dimension):
                result[..., i] = result[..., i] / W - non_derivative[..., i] * Wd / W / W;
            result = np.delete(result, self.dimension, 1)

        if all(s == 1 for s in result.shape[:-1]):
            result = result.reshape(self.dimension)

        return result

    def start(self, direction=None):
        """Return the start of the parametric domain"""
        if direction is None:
            return tuple(b.start() for b in self.bases)
        return self.bases[direction].start()

    def end(self, direction=None):
        """Return the end of the parametric domain"""
        if direction is None:
            return tuple(b.end() for b in self.bases)
        return self.bases[direction].end()

    def order(self, direction=None):
        """Return polynomial order (degree + 1)"""
        if direction is None:
            return tuple(b.order for b in self.bases)
        return self.bases[direction].order

    def translate(self, x):
        """Translate, i.e. move a B-spline object a given distance
        @param x: The direction and amount to move
        @type  x: Point_like
        """
        # 3D rational example: create a 4x4 translation matrix
        #
        #  |xw|      |  1   0   0  x1 |   |xw|
        #  |yw|   =  |  0   1   0  x2 | * |yw|
        #  |zw|      |  0   0   1  x3 |   |zw|
        #  | w|_new  |  0   0   0   1 |   | w|_old
        #
        #  PS: we even need a rational representation for non-rational splines
        #      in order to formulate translation as a matrix-matrix product
        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points
        if len(x) > dim:  # typical case: requesting movement in z-direction for 2D geometries
            self.set_dimension(len(x))
            dim = self.dimension

        # set up the translation matrix
        translation_matrix = np.matrix(np.identity(dim + 1))
        for i in range(dim):
            translation_matrix[i, -1] = x[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        if not self.rational:
            cp = np.matrix(np.ones((n, dim + 1)))  # pad with weights=1
            cp[:, :-1] = np.reshape(self.controlpoints, (n, dim))
        else:
            cp = np.matrix(np.reshape(self.controlpoints, (n, dim + rat)))

        # do the actual scaling by matrix-matrix multiplication
        cp = cp * translation_matrix.T  # right-mult, so we need transpose

        # store results
        if self.rational:
            self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)
        else:
            self.controlpoints = np.reshape(np.array(cp[:, :-1]), self.controlpoints.shape)

        return self

    def scale(self, s):
        """Scale, or magnify a B-spline object a given amount
        @param ax: Scaling factors, possible different in each direction
        @type  ax: Point_like or Float
        """
        # 3D rational example: create a 4x4 scaling matrix
        #
        #  |xw|      |  sx    0    0   0 |   |xw|
        #  |yw|   =  |   0   sy    0   0 | * |yw|
        #  |zw|      |   0    0   sz   0 |   |zw|
        #  | w|_new  |   0    0    0   1 |   | w|_old
        #
        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points
        s = ensure_listlike(s, dups=3)

        # set up the scaling matrix
        scale_matrix = np.matrix(np.identity(dim + rat))
        for i in range(dim):
            scale_matrix[i, i] = s[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.matrix(np.reshape(self.controlpoints, (n, dim + rat)))

        # do the actual scaling by matrix-matrix multiplication
        cp = cp * scale_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def rotate(self, theta, normal=[0, 0, 1]):
        """Rotate a B-spline object around an axis
        @param theta:  Angle to rotate, as measured in radians
        @type  theta:  Float
        @param normal: The normal axis (if 3D) to rotate object around
        @type  normal: Point_like
        """
        # 2D rational example: create a 3x3 rotation matrix
        #
        #  |xw|     | cos(t)  -sin(t)   0 |   |xw|
        #  |yw|   = | sin(t)   cos(t)   0 | * |yw|
        #  | w|_new |   0        0      1 |   | w|_old
        #

        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points
        if not (normal[0] == 0 and normal[1] == 0):  # rotating a 2D geometry out of the xy-plane
            self.set_dimension(3)
            dim = self.dimension

        # set up the rotation matrix
        if dim == 2:
            R = np.matrix([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]
                           ]).T  # we do right-multiplication, so we need a transpose
        elif dim == 3:
            normal = np.array(normal)
            R = get_rotation_matrix(theta, normal)
        else:
            raise RuntimeError('rotation undefined for geometries other than 2D and 3D')

        rotation_matrix = np.matrix(np.identity(dim + rat))
        rotation_matrix[0:dim, 0:dim] = R

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.matrix(np.reshape(self.controlpoints, (n, dim + rat)))

        # do the actual rotation by matrix-matrix multiplication
        cp = cp * rotation_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def mirror(self, normal):
        """Mirror a B-spline object around a plane through the origin
        @param normal: The plane normal to mirror about
        @type  normal: Point_like
        """
        # 3D rational example: create a 4x4 reflection matrix
        #
        #  normal = [a,b,c]
        #
        #  |xw|     |  1-2a^2   -2ab    -2ac   0 |   |xw|
        #  |yw|   = |   -2ab   1-2b^2   -2bc   0 | * |yw|
        #  |zw|   = |   -2ac    -2bc   1-2c^2  0 | * |zw|
        #  | w|_new |    0       0         0   1 |   | w|_old
        #
        #  PS: A reflection about a line is not a linear transformation; it is
        #      an affine transformation.

        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points

        if dim != 3:
            raise RuntimeError('reflection undefined for geometries other than 3D')

        # fixup the input normal to right form
        normal = np.array(normal)
        normal = normal / np.sqrt(np.dot(normal, normal))  # normalize it

        # set up the reflection matrix
        reflection_matrix = np.matrix(np.identity(dim + rat))
        reflection_matrix[0:dim, 0:dim] -= 2 * np.outer(normal, normal)

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.matrix(np.reshape(self.controlpoints, (n, dim + rat)))

        # do the actual rotation by matrix-matrix multiplication
        cp = cp * reflection_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def project(self, plane):
        """Projects geometry onto a plane or axis. project('xy') will project it
        onto the xy-plane, by setting all z-components to 0. project('y') will
        project it onto the y-axis by setting all x and z coordinates to 0.
        @param plane: Any combination of 'x', 'y' and/or 'z'
        @type  plane: String
        """
        keep = [False] * 3
        for s in plane:
            if s == 'x' or s == 'X':
                keep[0] = True
            if s == 'y' or s == 'Y':
                keep[1] = True
            if s == 'z' or s == 'Z':
                keep[2] = True

        dim = self.dimension
        for i in range(dim):
            if not keep[i]:
                self.controlpoints[..., i] = 0

        return self

    def bounding_box(self):
        """Get the bounding box of a B-spline computed off the control-point values.
        Might be inaccurate for rational splines with wild weight values
        @return: List of minimum and maximum values, sorted as xmin,xmax,ymin,ymax,...
        @rtype : list
        """
        dim = self.dimension

        result = []
        for i in range(dim):
            result.append(np.min(self.controlpoints[..., i]))
            result.append(np.max(self.controlpoints[..., i]))
        return result

    def set_dimension(self, new_dim):
        dim = self.dimension
        pardim = len(self.controlpoints.shape) - 1  # 1=Curve, 2=Surface, 3=Volume
        shape = self.controlpoints.shape
        while new_dim > dim:
            self.controlpoints = np.insert(self.controlpoints, dim, np.zeros(shape[:-1]), pardim)
            dim += 1
        while new_dim < dim:
            np.insert(self.controlpoints, dim, pardim)
            dim -= 1
        self.dimension = new_dim

        return self

    def force_rational(self):
        """Force a rational representation by including weights of all value 1"""
        if not self.rational:
            dim = self.dimension
            shape = self.controlpoints.shape
            pardim = len(self.controlpoints.shape) - 1  # 1=Curve, 2=Surface, 3=Volume
            self.controlpoints = np.insert(self.controlpoints, dim, np.ones(shape[:-1]), pardim)
            self.rational = 1

        return self

    def clone(self):
        return copy.deepcopy(self)

    __call__ = evaluate

    def __iadd__(self, x):
        self.translate(x)
        return self

    def __isub__(self, x):
        self.translate(-np.array(x))  # can't do -x if x is a list, so we rewrap it here
        return self

    def __imul__(self, x):
        self.scale(x)
        return self

    def __itruediv__(self, x):
        self.scale(1.0 / x)
        return self

    __ifloordiv__ = __itruediv__  # integer division (should not distinguish)
    __idiv__ = __itruediv__  # python2 compatibility

    def __add__(self, x):
        new_obj = copy.deepcopy(self)
        new_obj += x
        return new_obj

    def __sub__(self, x):
        new_obj = copy.deepcopy(self)
        new_obj -= x
        return new_obj

    def __mul__(self, x):
        new_obj = copy.deepcopy(self)
        new_obj *= x
        return new_obj

    def __div__(self, x):
        new_obj = copy.deepcopy(self)
        new_obj /= x
        return new_obj

import numpy as np

def get_rotation_matrix(theta,axis):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    return np.matrix([[a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)],
                     [2*(b*c+a*d),      a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c),      2*(c*d+a*b),     a*a+d*d-b*b-c*c]])

class ControlPointOperations:

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
        n   = len(self)     # number of control points

        # set up the translation matrix
        translation_matrix = np.matrix(np.identity(dim+1))
        for i in range(dim):
            translation_matrix[i,-1] = x[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        if not self.rational:
            cp        = np.matrix(np.ones((n, dim+1))) # pad with weights=1
            cp[:,:-1] = np.reshape(self.controlpoints, (n, dim))
        else:
            cp = np.matrix(np.reshape(self.controlpoints, (n, dim+rat)))

        # do the actual scaling by matrix-matrix multiplication
        cp  = cp * translation_matrix.T # right-mult, so we need transpose

        # store results
        if self.rational:
            self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)
        else:
            self.controlpoints = np.reshape(np.array(cp[:,:-1]), self.controlpoints.shape)

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
        n   = len(self)     # number of control points
        try:
            len(s)
        except TypeError:
            s = [s,s,s]

        # set up the scaling matrix
        scale_matrix = np.matrix(np.identity(dim+rat))
        for i in range(dim):
            scale_matrix[i,i] = s[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp  = np.matrix(np.reshape(self.controlpoints, (n, dim+rat)))

        # do the actual scaling by matrix-matrix multiplication
        cp  = cp * scale_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

    def rotate(self, theta, normal=[0,0,1]):
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
        n   = len(self)     # number of control points

        # set up the rotation matrix
        if dim==2:
            R = np.matrix([[np.cos(theta),  -np.sin(theta)],
                           [np.sin(theta),   np.cos(theta)]]).T # we do right-multiplication, so we need a transpose
        elif dim==3:
            normal = np.array(normal)
            R = get_rotation_matrix(theta, normal)
        else:
            raise RuntimeError('rotation undefined for geometries other than 2D and 3D')

        rotation_matrix = np.matrix(np.identity(dim+rat))
        rotation_matrix[0:dim,0:dim] = R 

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp  = np.matrix(np.reshape(self.controlpoints, (n, dim+rat)))

        # do the actual rotation by matrix-matrix multiplication
        cp  = cp * rotation_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

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
        n   = len(self)     # number of control points

        if dim != 3:
            raise RuntimeError('reflection undefined for geometries other than 3D')

        # fixup the input normal to right form
        normal = np.array(normal)
        normal = normal/np.sqrt(np.dot(normal,normal)) # normalize it

        # set up the reflection matrix
        reflection_matrix = np.matrix(np.identity(dim+rat))
        reflection_matrix[0:dim,0:dim] -= 2*np.outer(normal,normal)

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp  = np.matrix(np.reshape(self.controlpoints, (n, dim+rat)))

        # do the actual rotation by matrix-matrix multiplication
        cp  = cp * reflection_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

    def bounding_box(self):
        """Get the bounding box of a B-spline computed off the control-point values.
        Might be inaccurate for rational splines with wild weight values
        @return: List of minimum and maximum values, sorted as xmin,xmax,ymin,ymax,...
        @rtype : list
        """
        dim = self.dimension
        rat = self.rational
        n   = len(self)     # number of control points

        result       = []
        for i in range(dim):
            result.append(np.min(self.controlpoints[...,i]))
            result.append(np.max(self.controlpoints[...,i]))
        return result



    def __iadd__(self, x):
        self.translate(x)
        return self

    def __isub__(self, x):
        self.translate(-np.array(x)) # can't do -x if x is a list, so we rewrap it here
        return self

    def __imul__(self, x):
        self.scale(x)
        return self


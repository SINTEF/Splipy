import numpy as np

class ControlPointOperations:   

    def rotation_matrix(theta,axis):
        axis = axis/np.sqrt(np.dot(axis,axis))
        a = np.cos(theta/2)
        b,c,d = -axis*np.sin(theta/2)
        return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)],
                        [2*(b*c+a*d),      a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                        [2*(b*d-a*c),      2*(c*d+a*b),     a*a+d*d-b*b-c*c]])

    def Translate(self, x):
        """Translate, i.e. move a B-spline object a given distance
        @param x: The direction and amount to move
        @type  x: Point_like
        """
        x = np.array(x)
        for i in range(len(self)):
            self[i] = self[i] + x

    def Scale(self, x):
        """Scale, or magnify a B-spline object a given amount
        @param x: Scaling factor
        @type  x: Float
        """
        for i in range(len(self)):
            self[i] = self[i] * x

    def Rotate(self, theta, normal):
        """Rotate a B-spline object around an axis
        @param theta:  Angle to rotate, as measured in radians
        @type  theta:  Float
        @param normal: The normal axis to rotate object around
        @type  normal: Point_like
        """
        normal = np.array(normal)
        R = rotation_matrix(theta, normal);

        for i in range(len(self)):
            self[i] = R * self[i]

    def Mirror(self, point, normal):
        """Mirror a B-spline object around a plane
        @param point:  The point to mirror about
        @type  point:  Point_like
        @param normal: The plane normal to mirror about
        @type  normal: Point_like 
        """
        normal = np.array(normal)
        point  = np.array(point )
        for i in range(len(self)):
            move_length = 2 * np.inner(self[i] - point,  normal)
            self[i] = self[i] - move_length * normal

    def BoundingBox(self):
        """Get the bounding box of a B-spline computed off the control-point values.
        Might be inaccurate for rational splines with wild weight values
        @return: List of minimum and maximum values, sorted as xmin,xmax,ymin,ymax,...
        @rtype : list
        """
        controlpoint = self[0]
        dimension    = len(controlpoint)
        result       = []
        for i in range(dimension):
            result.append( np.inf)
            result.append(-np.inf)
        for i in range(len(self)):
            controlpoint = self[i]
            for j in range(dimension):
                result[2*j  ] = min(controlpoint[j], result[2*j  ])
                result[2*j+1] = max(controlpoint[j], result[2*j+1])
        return result



    def __iadd__(self, x):
        self.Translate(x)
        return self

    def __isub__(self, x):
        self.Translate(-np.array(x)) # can't do -x if x is a list, so we rewrap it here
        return self

    def __imul__(self, x):
        self.Scale(x)
        return self

    def __idiv__(self, x):
        self.Scale(1.0/x)
        return self

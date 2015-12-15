from BSplineBasis import *
from ControlPointOperations import *
import numpy as np

class Curve(ControlPointOperations):

    def __init__(self, order=None, knot=None, controlpoints=None, rational=False):

        self.basis = BSplineBasis(order, knot)
        
        # if none provided, create the default geometry which is the linear mapping onto the unit line [0,0]->[1,0]
        if controlpoints is None:
            controlpoints = []
            greville_points = self.basis.greville()
            for point in greville_points:
                controlpoints.append([point, 0])
        
        self.controlpoints = np.array(controlpoints)
        self.rational      = rational
        self.dimension     = len(controlpoints[0]) - rational

    def evaluate(self, t):
        """Evaluate the curve at given parametric values
        @param t: Parametric coordinate point(s)
        @type  t: Float or list of Floats
        @return : Geometry coordinates. Matrix X(i,j) of component x(j) evaluated at t(i)
        @rtype  : numpy.array
        """
        # compute basis functions for all points t. N(i,j) is a matrix of all functions j for all points i
        N = self.basis.evaluate(t)

        # compute physical points [x,y,z] for all points t[i]. For rational curves, compute [X,Y,Z,W] (in projective space)
        result = N * self.controlpoints

        # Project rational curves down to geometry space: x = X/W, y=Y/W, z=Z/W
        if self.rational: 
            for i in range(self.dimension):
                result[:,i] /= result[:,-1] 
            result = np.delete(result, self.dimension, 1) # remove the weight column

        if result.shape[0] == 1: # in case of single value input t, return vector instead of matrix
            result = np.array(result[0,:]).reshape(self.dimension)

        return result

    def evaluate_tangent(self, t):
        """Evaluate the tangent of the curve at given parametric values
        @param t: Parametric coordinate point(s)
        @type  t: Float or list of Floats
        @return : Tangent evaluation. Matrix DX(i,j) of component x(j) evaluated at t(i)
        @rtype  : numpy.array
        """
        # compute basis functions for all points t. dN(i,j) is a matrix of the derivative of all functions j for all points i
        dN = self.basis.evaluate(t, 1)

        # compute physical points [dx/dt,dy/dt,dz/dt] for all points t[i]
        result = dN * self.controlpoints

        # Rational curves need the quotient rule to compute derivatives (controlpoints are stored as x_i*w_i)
        # x(t) = sum_i N_i(t) * w_i*x_i / W(t)
        # W(t) = sum_j N_j(t) * w_j
        # dx/dt =  sum_i N_i'(t)*w_i*x_i / W(t) - sum_i N_i(t)*w_i*x_i* W'(t)/W(t)^2
        if self.rational: 
            N = self.basis.evaluate(t)
            non_derivative = N*self.controlpoints
            W    = non_derivative[:,-1]  # W(t)
            Wder = result[:,-1]          # W'(t)
            for i in range(self.dimension):
                result[:,i] = result[:,i]/W - non_derivative[i,:]*Wder/W/W

            result = np.delete(result, self.dimension, 1) # remove the weight column

        return result

    def flip_parametrization(self):
        """Swap direction of the curve by making it go in the reverse direction. Parametric domain remain unchanged"""
        self.basis.reverse()
        self.controlpoints = self.controlpoints[::-1]

    def get_order(self):
        """Return polynomial order (degree + 1) of spline curve"""
        return self.basis.order

    def get_knots(self, with_multiplicities=False):
        """Get the knots of a spline curve
        @param with_multiplicities: Set to true to obtain the knot vector with multiplicities
        @type with_multiplicities : Boolean
        @return:                    List with the knot values
        @rtype :                    List of float
        """
        if with_multiplicities:
            return self.basis.knots
        else:
            return self.basis.get_knot_spans()

    def force_rational(self):
        """Force a rational representation by including weights of all value 1"""
        if not self.rational:
            n,d = self.controlpoints.shape # n controlpoints of dimension d
            self.controlpoints = np.insert(self.controlpoints, d, np.array([1]*n), 1)
            self.rational = 1

    def reparametrize(self, start=0, end=1):
        """Redefine the parametric domain to be (start,end)"""
        if end <= start:
            raise ValueError('end must be larger than start')
        self.basis.normalize()     # set domain to (0,1)
        self.basis *= (end-start)
        self.basis += start

    def raise_order(self, amount):
        # create the new basis
        newKnot  = self.basis.get_raise_order_knot(amount)
        newBasis = BSplineBasis(self.basis.order + amount, newKnot, self.basis.periodic)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_t = newBasis.greville()        # parametric interpolation points (t)
        N_old = self.basis.evaluate(interpolation_pts_t) 
        N_new = newBasis.evaluate(  interpolation_pts_t) 
        interpolation_pts_x = N_old * self.controlpoints # projective interpolation points (x,y,z,w)

        # solve the interpolation problem
        self.controlpoints = np.linalg.solve(N_new, interpolation_pts_x)
        self.basis         = newBasis
        




    def __len__(self):
        """return the number of control points (basis functions) for this curve"""
        return len(self.basis)

    def __getitem__(self, i):
        if self.rational:
            return self.controlpoints[i,:-1] / self.controlpoints[i,-1]
        else:
            return self.controlpoints[i,:]

    def __setitem__(self, i, newCP):
        if self.rational:
            self.controlpoints[i,:-1] = newCP * self.controlpoints[i,-1]
        else:
            self.controlpoints[i,:]   = newCP
        return self

    def __repr__(self):
        return str(self.basis) + '\n' + str(self.controlpoints)


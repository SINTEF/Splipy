from BSplineBasis import *
from ControlPointOperations import *
import numpy as np

class Surface(ControlPointOperations):

    def __init__(self, order1=2, order2=2, knot1=None, knot2=None, controlpoints=None, rational=False):

        self.basis1 = BSplineBasis(order1, knot1)
        self.basis2 = BSplineBasis(order2, knot2)
        
        # if none provided, create the default geometry which is the linear mapping onto the unit square (0,1)^2
        if controlpoints == None:
            controlpoints = []
            greville_points1 = self.basis1.Greville()
            greville_points2 = self.basis2.Greville()
            for p2 in greville_points2:
                for p1 in greville_points1:
                    controlpoints.append([p1,p2])
        
        self.dimension     = len(controlpoints[0]) - rational
        self.rational      = rational

        # controlpoints are given in at 2-index (ji,k) for u[i], v[j], x[k]
        # reshape theese into 3-index (j,i,k)
        self.controlpoints = np.reshape(controlpoints, (len(self.basis2), len(self.basis1), self.dimension + self.rational ))
        # swap axis 0 and 1, to make it (i,j,k) 
        self.controlpoints = self.controlpoints.transpose((1,0,2))

    def Evaluate(self, u, v):
        """Evaluate the surface at given parametric values
        @param u: Parametric coordinate point(s) in first direction
        @type  u: Float or list of Floats
        @param v: Parametric coordinate point(s) in second direction
        @type  v: Float or list of Floats
        @return : Geometry coordinates. 3D-array X(i,j,k) of component x(k) evaluated at (u[i],v[j])
        @rtype  : numpy.array
        """
        # compute basis functions for all points t. Nu(i,j) is a matrix of all functions j for all points u[i]
        Nu = self.basis1.Evaluate(u)
        Nv = self.basis2.Evaluate(v)

        # compute physical points [x,y,z] for all points (u[i],v[j]). For rational surfaces, compute [X,Y,Z,W] (in projective space)
        result = np.tensordot(Nu, self.controlpoints, axes=(1,0))
        result = np.tensordot(Nv, result,             axes=(1,1))

        # Project rational surfaces down to geometry space: x = X/W, y=Y/W, z=Z/W
        if self.rational: 
            for i in range(self.dimension):
                result[:,:,i] /= result[:,:,-1] 
            result = np.delete(result, self.dimension, 2) # remove the matrix of weight evaluation

        if result.shape[0] == 1 and result.shape[1] == 1: # in case of single value input (u,v), return vector instead of 3D-matrix
            result = np.array(result[0,0,:]).reshape(self.dimension)
        return result

    def FlipParametrization(self, direction):
        """Swap direction of the surface by making it go in the reverse direction. Parametric domain remain unchanged
           @param direction: The parametric direction to flip (0=u, 1=v)
           @type  direction: Int
        """
        if direction==0:
            self.basis1.Reverse()
            self.controlpoints = self.controlpoints[::-1, :,:]
        elif direction==1:
            self.basis2.Reverse()
            self.controlpoints = self.controlpoints[:, ::-1 ,:]
        else:
            raise ValueError('direction must be 0 or 1')

    def GetOrder(self):
        """Return spline surface order (polynomial degree + 1) in all parametric directions"""
        return (self.basis1.order, self.basis2.order)

    def GetKnots(self, with_multiplicities=False):
        """Get the knots of the spline surface
        @param with_multiplicities: Set to true to obtain the knot vector with multiplicities
        @type with_multiplicities : Boolean
        @return:                    List with the knot values
        @rtype :                    Tuple with List of float
        """
        if with_multiplicities:
            return (self.basis1.knots, self.basis2.knots)
        else:
            return (self.basis1.GetKnotSpans(), self.basis2.GetKnotSpans())

    def ForceRational(self):
        """Force a rational representation by including weights of all value 1"""
        if not self.rational:
            n1,n2,d = self.controlpoints.shape # n1 x n2 controlpoints of dimension d
            self.controlpoints = np.insert(self.controlpoints, d, np.ones((n1,n2)), 2)
            self.rational = 1

    def SwapParametrization(self):
        """Swaps the two surface parameter directions"""
        self.controlpoints = self.controlpoints.transpose((1,0,2))  # re-order controlpoints
        self.basis1, self.basis2 = self.basis2, self.basis1         # swap knot vectors

    def ReParametrize(self, umin=0, umax=1, vmin=0, vmax=1):
        """Redefine the parametric domain to be (umin,umax) x (vmin,vmax)"""
        if umax <= umin or vmax <= vmin:
            raise ValueError('end must be larger than start')
        self.basis1.Normalize()     # set domain to (0,1)
        self.basis1 *= (umax-umin)
        self.basis1 += umin
        self.basis2.Normalize()
        self.basis2 *= (vmax-vmin)
        self.basis2 += vmin

    def GetEdges(self):
        """Return the four edge curves in (parametric) order: bottom, right, top, left"""
        # ASSUMPTION: open knot vectors
        (p1,p2)     = self.GetOrder()
        (n1,n2,dim) = self.controlpoints.shape
        rat         = self.rational
        umin = Curve(p2, self.basis2, np.reshape(self.controlpoints[ 0,:,:], (n2,dim), rat))
        umax = Curve(p2, self.basis2, np.reshape(self.controlpoints[-1,:,:], (n2,dim), rat))
        vmin = Curve(p1, self.basis1, np.reshape(self.controlpoints[:, 0,:], (n1,dim), rat))
        vmax = Curve(p1, self.basis1, np.reshape(self.controlpoints[:,-1,:], (n1,dim), rat))
        # make the curves form a clockwise oriented closed loop around surface
        umax.FlipParametrization()
        vmax.FlipParametrization()
        return [vmin, umax, vmax, umin]




    def __len__(self):
        """return the number of control points (basis functions) for this surface"""
        return len(self.basis1)*len(self.basis2)

    def __getitem__(self, i):
        (n1,n2,dim) = self.controlpoints.shape
        i1 =     i % n1
        i2 = int(i / n1)
        if self.rational:
            return self.controlpoints[i1,i2,:-1] / self.controlpoints[i1,i2,-1]
        else:
            return self.controlpoints[i1,i2,:]

    def __setitem__(self, i, newCP):
        (n1,n2,dim) = self.controlpoints.shape
        i1 =     i % n1
        i2 = int(i / n1)
        if self.rational:
            self.controlpoints[i1,i2,:-1] = newCP * self.controlpoints[i1,i2,-1]
        else:
            self.controlpoints[i1,i2,:]   = newCP
        return self

    def __repr__(self):
        result = str(self.basis1) + '\n' + str(self.basis2) + '\n'
        # print legacy controlpoint enumeration
        n1,n2,n3 = self.controlpoints.shape
        for j in range(n2):
            for i in range(n1):
                result += str(self.controlpoints[i,j,:]) + '\n'
        return result


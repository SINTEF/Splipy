from BSplineBasis import *
from ControlPointOperations import *
import numpy as np

class Surface(ControlPointOperations):

    def __init__(self, basis1=None, basis2=None, controlpoints=None, rational=False):

        if basis1 is None:
            basis1 = BSplineBasis()
        if basis2 is None:
            basis2 = BSplineBasis()

        self.basis1 = basis1
        self.basis2 = basis2
        
        # if none provided, create the default geometry which is the linear mapping onto the unit square (0,1)^2
        if controlpoints is None:
            controlpoints = []
            greville_points1 = self.basis1.greville()
            greville_points2 = self.basis2.greville()
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

    def evaluate(self, u, v):
        """Evaluate the surface at given parametric values
        @param u: Parametric coordinate point(s) in first direction
        @type  u: Float or list of Floats
        @param v: Parametric coordinate point(s) in second direction
        @type  v: Float or list of Floats
        @return : Geometry coordinates. 3D-array X(i,j,k) of component x(k) evaluated at (u[i],v[j])
        @rtype  : numpy.array
        """
        # for single-value input, wrap it into a list
        try:
            len(u)
        except TypeError:
            u = [u]
        try:
            len(v)
        except TypeError:
            v = [v]

        # error test input
        if self.basis1.periodic < 0: # periodic functions can evaluate everywhere
            if min(u) < self.basis1.start() or self.basis1.end() < max(u):
                raise ValueError('evaluation outside parametric domain')
        if self.basis2.periodic < 0: # periodic functions can evaluate everywhere
            if min(v) < self.basis2.start() or self.basis2.end() < max(v):
                raise ValueError('evaluation outside parametric domain')

        # compute basis functions for all points t. Nu(i,j) is a matrix of all functions j for all points u[i]
        Nu = self.basis1.evaluate(u)
        Nv = self.basis2.evaluate(v)

        # compute physical points [x,y,z] for all points (u[i],v[j]). For rational surfaces, compute [X,Y,Z,W] (in projective space)
        result = np.tensordot(Nv, self.controlpoints, axes=(1,1))
        result = np.tensordot(Nu, result,             axes=(1,1))

        # Project rational surfaces down to geometry space: x = X/W, y=Y/W, z=Z/W
        if self.rational: 
            for i in range(self.dimension):
                result[:,:,i] /= result[:,:,-1] 
            result = np.delete(result, self.dimension, 2) # remove the matrix of weight evaluation

        if result.shape[0] == 1 and result.shape[1] == 1: # in case of single value input (u,v), return vector instead of 3D-matrix
            result = np.array(result[0,0,:]).reshape(self.dimension)
        return result

    def flip_parametrization(self, direction):
        """Swap direction of the surface by making it go in the reverse direction. Parametric domain remain unchanged
           @param direction: The parametric direction to flip (0=u, 1=v)
           @type  direction: Int
        """
        if direction==0:
            self.basis1.reverse()
            self.controlpoints = self.controlpoints[::-1, :,:]
        elif direction==1:
            self.basis2.reverse()
            self.controlpoints = self.controlpoints[:, ::-1 ,:]
        else:
            raise ValueError('direction must be 0 or 1')

    def get_order(self):
        """Return spline surface order (polynomial degree + 1) in all parametric directions"""
        return (self.basis1.order, self.basis2.order)

    def get_knots(self, with_multiplicities=False):
        """Get the knots of the spline surface
        @param with_multiplicities: Set to true to obtain the knot vector with multiplicities
        @type with_multiplicities : Boolean
        @return:                    List with the knot values
        @rtype :                    Tuple with List of float
        """
        if with_multiplicities:
            return (self.basis1.knots, self.basis2.knots)
        else:
            return (self.basis1.get_knot_spans(), self.basis2.get_knot_spans())

    def force_rational(self):
        """Force a rational representation by including weights of all value 1"""
        if not self.rational:
            n1,n2,d = self.controlpoints.shape # n1 x n2 controlpoints of dimension d
            self.controlpoints = np.insert(self.controlpoints, d, np.ones((n1,n2)), 2)
            self.rational = 1

    def swap_parametrization(self):
        """Swaps the two surface parameter directions"""
        self.controlpoints = self.controlpoints.transpose((1,0,2))  # re-order controlpoints
        self.basis1, self.basis2 = self.basis2, self.basis1         # swap knot vectors

    def reparametrize(self, umin=0, umax=1, vmin=0, vmax=1):
        """Redefine the parametric domain to be (umin,umax) x (vmin,vmax)"""
        if umax <= umin or vmax <= vmin:
            raise ValueError('end must be larger than start')
        self.basis1.normalize()     # set domain to (0,1)
        self.basis1 *= (umax-umin)
        self.basis1 += umin
        self.basis2.normalize()
        self.basis2 *= (vmax-vmin)
        self.basis2 += vmin

    def get_edges(self):
        """Return the four edge curves in (parametric) order: bottom, right, top, left"""
        # ASSUMPTION: open knot vectors
        (p1,p2)     = self.get_order()
        (n1,n2,dim) = self.controlpoints.shape
        rat         = self.rational
        umin = Curve(p2, self.basis2, np.reshape(self.controlpoints[ 0,:,:], (n2,dim), rat))
        umax = Curve(p2, self.basis2, np.reshape(self.controlpoints[-1,:,:], (n2,dim), rat))
        vmin = Curve(p1, self.basis1, np.reshape(self.controlpoints[:, 0,:], (n1,dim), rat))
        vmax = Curve(p1, self.basis1, np.reshape(self.controlpoints[:,-1,:], (n1,dim), rat))
        # make the curves form a clockwise oriented closed loop around surface
        umax.flip_parametrization()
        vmax.flip_parametrization()
        return [vmin, umax, vmax, umin]

    def raise_order(self, raise_u, raise_v):
        """Raise the order of a spline surface
        @param raise_u: Number of polynomial degrees to increase in u
        @type  raise_u: Int
        @param raise_v: Number of polynomial degrees to increase in v
        @type  raise_v: Int
        """
        # create the new basis
        newKnot1  = self.basis1.get_raise_order_knot(raise_u)
        newBasis1 = BSplineBasis(self.basis1.order + raise_u, newKnot1, self.basis1.periodic)
        newKnot2  = self.basis2.get_raise_order_knot(raise_v)
        newBasis2 = BSplineBasis(self.basis2.order + raise_v, newKnot2, self.basis2.periodic)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_u = newBasis1.greville()        # parametric interpolation points u
        interpolation_pts_v = newBasis2.greville()        # parametric interpolation points v
        N_u_old = self.basis1.evaluate( interpolation_pts_u )
        N_u_new =   newBasis1.evaluate( interpolation_pts_u )
        N_v_old = self.basis2.evaluate( interpolation_pts_v )
        N_v_new =   newBasis2.evaluate( interpolation_pts_v )
        tmp = np.tensordot(N_v_old, self.controlpoints, axes=(1,1))
        tmp = np.tensordot(N_u_old, tmp,                axes=(1,1)) # projective interpolation points (x,y,z,w)
        interpolation_pts_x = tmp # 3D-tensor with elements (i,j,k) of component x[k] evaluated at u[i] v[j]

        # solve the interpolation problem
        N_u_inv = np.linalg.inv(N_u_new)
        N_v_inv = np.linalg.inv(N_v_new) # these are inverses of the 1D problems, and small compared to the total number of unknowns
        tmp = np.tensordot(N_v_inv, interpolation_pts_x, axes=(1,1))
        tmp = np.tensordot(N_u_inv, tmp,                 axes=(1,1))

        # update the basis and controlpoints of the surface
        self.controlpoints = tmp
        self.basis1        = newBasis1
        self.basis2        = newBasis2
        




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


from BSplineBasis import *
from ControlPointOperations import *
import numpy as np

class Volume(ControlPointOperations):

    def __init__(self, order1=2, order2=2, order3=2, knot1=None, knot2=None, knot3=None, controlpoints=None, rational=False):

        self.basis1 = BSplineBasis(order1, knot1)
        self.basis2 = BSplineBasis(order2, knot2)
        self.basis3 = BSplineBasis(order3, knot3)
        
        # if none provided, create the default geometry which is the linear mapping onto the unit cube (0,1)^3
        if controlpoints == None:
            controlpoints = []
            greville_points1 = self.basis1.greville()
            greville_points2 = self.basis2.greville()
            greville_points3 = self.basis3.greville()
            for p3 in greville_points3:
              for p2 in greville_points2:
                  for p1 in greville_points1:
                      controlpoints.append([p1,p2,p3])
        
        self.dimension     = len(controlpoints[0]) - rational
        self.rational      = rational

        # controlpoints are given in as 2-index (kji,l) for u[i], v[j], w[k], x[l]
        # reshape these into 4-index (k,j,i,l)
        self.controlpoints = np.reshape(controlpoints, (len(self.basis3), len(self.basis2), len(self.basis1), self.dimension + self.rational ))
        # swap axis 0 and 2, to make it (i,j,k,l) 
        self.controlpoints = self.controlpoints.transpose((2,1,0,3))

    def evaluate(self, u, v, w):
        """Evaluate the volume at given parametric values
        @param u: Parametric coordinate point(s) in first direction
        @type  u: Float or list of Floats
        @param v: Parametric coordinate point(s) in second direction
        @type  v: Float or list of Floats
        @param w: Parametric coordinate point(s) in second direction
        @type  w: Float or list of Floats
        @return : Geometry coordinates. 4D-array X(i,j,k,l) of component x(l) evaluated at (u[i],v[j],w[k])
        @rtype  : numpy.array
        """
        # compute basis functions for all points t. Nu(i,j) is a matrix of all functions j for all points u[i]
        Nu = self.basis1.evaluate(u)
        Nv = self.basis2.evaluate(v)
        Nw = self.basis3.evaluate(w)

        # compute physical points [x,y,z] for all points (u[i],v[j],w[k]). For rational volumes, compute [X,Y,Z,W] (in projective space)
        result = np.tensordot(Nu, self.controlpoints, axes=(1,0))
        result = np.tensordot(Nv, result,             axes=(1,1))
        result = np.tensordot(Nw, result,             axes=(1,2))

        # Project rational volumes down to geometry space: x = X/W, y=Y/W, z=Z/W
        if self.rational: 
            for i in range(self.dimension):
                result[:,:,:,i] /= result[:,:,:,-1] 
            result = np.delete(result, self.dimension, 3) # remove all weight evaluations

        # in case of single value input (u,v,w), return vector instead of 4D-matrix
        if result.shape[0] == 1 and result.shape[1] == 1 and result.shape[2] == 1:
            result = np.array(result[0,0,0,:]).reshape(self.dimension)
        return result

    def flip_parametrization(self, direction):
        """Swap direction of the volume by making it go in the reverse direction. Parametric domain remain unchanged
           @param direction: The parametric direction to flip (0=u, 1=v, 2=w)
           @type  direction: Int
        """
        if direction==0:
            self.basis1.reverse()
            self.controlpoints = self.controlpoints[::-1, :,:,:]
        elif direction==1:
            self.basis2.reverse()
            self.controlpoints = self.controlpoints[:, ::-1 ,:,:]
        elif direction==2:
            self.basis3.reverse()
            self.controlpoints = self.controlpoints[:,:, ::-1 ,:]
        else:
            raise ValueError('direction must be 0,1 or 2')

    def get_order(self):
        """Return spline volume order (polynomial degree + 1) in all parametric directions"""
        return (self.basis1.order, self.basis2.order, self.basis3.order)

    def get_knots(self, with_multiplicities=False):
        """Get the knots of the spline volume
        @param with_multiplicities: Set to true to obtain the knot vector with multiplicities
        @type with_multiplicities : Boolean
        @return:                    List with the knot values
        @rtype :                    Tuple with List of float
        """
        if with_multiplicities:
            return (self.basis1.knots, self.basis2.knots, self.basis3.knots)
        else:
            return (self.basis1.get_knot_spans(), self.basis2.get_knot_spans(), self.basis3.get_knot_spans())

    def force_rational(self):
        """Force a rational representation by including weights of all value 1"""
        if not self.rational:
            n1,n2,n3,d = self.controlpoints.shape # n1 x n2 x n3 controlpoints of dimension d
            self.controlpoints = np.insert(self.controlpoints, d, np.ones((n1,n2,n3)), 3)
            self.rational = 1

    def swap_parametrization(self, pardir1, pardir2):
        """Swaps the two volume parameter directions"""
        if (pardir1==0 and pardir2==1) or (pardir1==1 and pardir2==0):
            self.controlpoints = self.controlpoints.transpose((1,0,2,3))  # re-order controlpoints
            self.basis1, self.basis2 = self.basis2, self.basis1           # swap knot vectors
        elif (pardir1==0 and pardir2==2) or (pardir1==2 and pardir2==0):
            self.controlpoints = self.controlpoints.transpose((2,1,0,3))  # re-order controlpoints
            self.basis1, self.basis3 = self.basis3, self.basis1           # swap knot vectors
        elif (pardir1==1 and pardir2==2) or (pardir1==2 and pardir2==1):
            self.controlpoints = self.controlpoints.transpose((0,2,1,3))  # re-order controlpoints
            self.basis2, self.basis3 = self.basis3, self.basis2           # swap knot vectors
        else:
            raise ValueError('pardir1 and pardir2 must be different from each other and either 0,1 or 2')


    def reparametrize(self, umin=0, umax=1, vmin=0, vmax=1, wmin=0, wmax=1):
        """Redefine the parametric domain to be (umin,umax) x (vmin,vmax) x (wmin,wmax)"""
        if umax <= umin or vmax <= vmin or wmax <= wmin:
            raise ValueError('end must be larger than start')
        self.basis1.normalize()     # set domain to (0,1)
        self.basis1 *= (umax-umin)
        self.basis1 +=  umin
        self.basis2.normalize()
        self.basis2 *= (vmax-vmin)
        self.basis2 +=  vmin
        self.basis3.normalize()
        self.basis3 *= (wmax-wmin)
        self.basis3 +=  wmin

    def get_faces(self):
        """Return a list of the 6 boundary faces of this volume (with outward normal vector). They are ordered as (umin,umax,vmin,vmax,wmin,wmax)"""
        # ASSUMPTION: open knot vectors
        (p1,p2,p3)     = self.get_order()
        (n1,n2,n3,dim) = self.controlpoints.shape
        rat            = self.rational
        umin = Surface(p3,p2, self.basis3, self.basis2, np.reshape(self.controlpoints[ 0,:,:,:], (n2*n3,dim), rat))
        umax = Surface(p3,p2, self.basis3, self.basis2, np.reshape(self.controlpoints[-1,:,:,:], (n2*n3,dim), rat))
        vmin = Surface(p3,p1, self.basis3, self.basis1, np.reshape(self.controlpoints[:, 0,:,:], (n1*n3,dim), rat))
        vmax = Surface(p3,p1, self.basis3, self.basis1, np.reshape(self.controlpoints[:,-1,:,:], (n1*n3,dim), rat))
        wmin = Surface(p2,p1, self.basis2, self.basis1, np.reshape(self.controlpoints[:,:, 0,:], (n1*n2,dim), rat))
        wmax = Surface(p2,p1, self.basis2, self.basis1, np.reshape(self.controlpoints[:,:,-1,:], (n1*n2,dim), rat))
        umax.swap_parametrization()
        vmax.swap_parametrization()
        wmax.swap_parametrization()
        return [umin, umax, vmin, vmax, wmin, wmax]

        




    def __len__(self):
        """return the number of control points (basis functions) for this volume"""
        return len(self.basis1)*len(self.basis2)*len(self.basis3)

    def __getitem__(self, i):
        (n1,n2,n3,dim) = self.controlpoints.shape
        i1 = int(i % n1     )
        i2 = int(i / n1)% n2
        i3 = int(i / n1 / n2)
        if self.rational:
            return self.controlpoints[i1,i2,i3,:-1] / self.controlpoints[i1,i2,i3,-1]
        else:
            return self.controlpoints[i1,i2,i3,:]

    def __setitem__(self, i, newCP):
        (n1,n2,n3,dim) = self.controlpoints.shape
        i1 = int(i % n1     )
        i2 = int(i / n1)% n2
        i3 = int(i / n1 / n2)
        if self.rational:
            self.controlpoints[i1,i2,i3,:-1] = newCP * self.controlpoints[i1,i2,i3,-1]
        else:
            self.controlpoints[i1,i2,i3,:]   = newCP
        return self

    def __repr__(self):
        result  = str(self.basis1) + '\n';
        result += str(self.basis2) + '\n';
        result += str(self.basis3) + '\n';
        # print legacy controlpoint enumeration
        n1,n2,n3,dim = self.controlpoints.shape
        for k in range(n3):
            for j in range(n2):
                for i in range(n1):
                    result += str(self.controlpoints[i,j,k,:]) + '\n'
        return result


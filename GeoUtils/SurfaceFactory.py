__doc__ = 'Implementation of various methods which should have been in GoTools.SurfaceFactory or GoTools.Surface.'

from GoTools import *

def SplitSurface(original, param, value):
        """Split a SplineSurface in two by continued knot insertion
        @param original: The Surface to split
        @type original: Surface
        @param param: The parametric direction to split along (u=0, v=1)
        @type param: Int
        @param value: The parametric value to split along
        @type value: Float
        @return: The two new surfaces
        @rtype: List of Surfaces
        """
        
        surf         = original.Clone()
 # this actually fails for 2D models, but no Surface.GetDim() or 
 # Surface.IsRational()  :(
        rational     = (len(surf[0]) == 4)
        knot1, knot2 = surf.GetKnots(True);
        p            = surf.GetOrder();

        # error check input (Note: that rational should really fail unless you know what you're doing)
#         if rational:
#                 return None
        if(param == 0):
                if(value < knot1[0] or knot1[-1] < value):
                        return None
                knot = knot1
        elif(param == 1):
                if(value < knot2[0] or knot2[-1] < value):
                        return None
                knot = knot2
        else:
                return None

        # Figure out if knot already exist
        knotCount = 0
        for i in range(len(knot)):
                if abs(knot[i]-value) < GetTolerance('neighbour'):
                        value = knot[i] # snap to existing knot value
                        knotCount = knotCount + 1

        # insert the actual knot
        for i in range(knotCount, p[param]-1):
                surf.InsertKnot(param, value)

        knot1, knot2 = surf.GetKnots(True)
        if param==0:
                knot = knot1
        else:
                knot = knot2
        # figure out which knot interval (ki) the split lies
        for ki in range(len(knot)):
                if abs(knot[ki]-value) < GetTolerance('neighbour'):
                        break
        ki = ki-1
        
        n = [len(knot1) - p[0] , len(knot2) - p[1]]

        coefs1 = []
        coefs2 = []
        knot1u = []
        knot1v = []
        knot2u = []
        knot2v = []
        
        k = 0
        if param == 0:
                for y in range(n[1]):
                        for x in range(n[0]):
                                if x <= ki:
                                        coefs1.append(surf[k])
                                if x >= ki:
                                        coefs2.append(surf[k])
                                k = k+1

                for i in range(ki+p[0]):
                        knot1u.append(knot1[i])
                knot1u.append(knot1u[-1])
                for i in range(ki+1, len(knot1)):
                        knot2u.append(knot1[i])
                knot2u = [knot2u[0]] + knot2u

                knot1v = knot2
                knot2v = knot2
        else:
                for y in range(n[1]):
                        for x in range(n[0]):
                                if y <= ki:
                                        coefs1.append(surf[k])
                                if y >= ki:
                                        coefs2.append(surf[k])
                                k = k+1

                for i in range(ki+p[1]):
                        knot1v.append(knot2[i])
                knot1v.append(knot1v[-1])
                for i in range(ki+1, len(knot2)):
                        knot2v.append(knot2[i])
                knot2v = [knot2v[0]] + knot2v

                knot1u = knot1
                knot2u = knot1

#         print knot1u
#         print knot1v
#         print knot2u
#         print knot2v
# 
#         for point in surf:
#                 print point
#         print '****     SURFACE 1    ****'
#         for point in coefs1:
#                 print point
#         print '****     SURFACE 2    ****'
#         for point in coefs2:
#                 print point
        
        result = []
        result.append(Surface(p[0], p[1], knot1u, knot1v, coefs1, rational))
        result.append(Surface(p[0], p[1], knot2u, knot2v, coefs2, rational))
        return  result
        
        # p = original.GetOrder()
        # knots1, knots2, knots3 = volume.GetKnotsWithMult()
        # return Volume(p[0], p[1], p[2], knots1, knots2, knots3, mirroredCP, rational)




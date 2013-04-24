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
        knot1, knot2 = surf.GetKnots(True)
        p            = surf.GetOrder()

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


def RotationalCrv2CrvSweep(crv1, crv2, axis):
        """ Rotational sweep the surface in between to curves,
        where both curves live on the same cylinder.
        CONSTRAINT: Curves must share the height-parametrization, i.e.
        lie share z-coordinates if living on the unit cylinder
        @param crv1: The first curve on the cylinder
        @type  crv1: Curve
        @param crv2: The second curve on the cylinder
        @type  crv2: Curve
        @param axis: The cylinder axis (0=x-axis, 1=y, 2=z)
        @type  axis: int
        @return : Swept surface
        @rtype  : Surface
        """
        controlPoints = []

        listAxis = [0,0,0]
        listAxis[axis] = 1.0
        nAxis = Point(list=listAxis)
        
        for i in range(len(crv1)):

#                   A   
#          , - ~ ~ ~ -0,
#      , '           /   ' ,
#    ,              /        ;,
#   ,           r1 /          ,'.  ControlPt 1
#  ,              /            , '0
#  ,             * ) Theta     , /
#  ,              \            ,/
#   ,              \          ,0  ControlPt 2
#    ,           r2 \        ,/
#      ,             \    , '/
#        ' - , _ _ _ ,0_"___0  ControlPt 3
#                     B
#
#       going for two knot spans. Should work for all theta 
#       strictly less than 2pi (though bad results for close
#       to 2pi
#   

                a = crv1[i]
                b = crv2[i]
                origin1 = (a*nAxis)*nAxis
                origin2 = (b*nAxis)*nAxis
                dist    = abs(origin1-origin2)
                r1      = abs(origin1-a)
                r2      = abs(origin2-b)

                # !!! HARDCODE FIX !!!
                # if the first curve starts or stops tangent to the cylinder radial parametrization
                # we try and fix this by tweaking the surrounding control points. In this case, setting
                # the second CP to the average height between the first two.
                if i==1 and abs(crv1[i][axis]-crv1[i-1][axis]) < 1e-3:
                        # print 'fixing height'
                        newHeight = (crv1[i+1][axis]+crv1[i][axis]) / 2.0
                        origin1 = Point( newHeight, origin1[1], origin1[2])
                        origin2 = Point( newHeight, origin2[1], origin2[2])
                        a = Point(newHeight, a[1], a[2])
                        b = Point(newHeight, b[1], b[2])
                elif i==len(crv1)-2 and abs(crv1[i][axis]-crv1[i+1][axis]) < 1e-3:
                        # print 'fixing height'
                        newHeight = (crv1[i-1][axis]+crv1[i][axis]) / 2.0
                        origin1 = Point( newHeight, origin1[1], origin1[2])
                        origin2 = Point( newHeight, origin2[1], origin2[2])
                        a = Point(newHeight, a[1], a[2])
                        b = Point(newHeight, b[1], b[2])
                        

                ###   this test isn't really true since the control points are floating outside
                ###   the actual curve in geometry space (but nice to have to illustrate whats 
                ###   going on)
                # if abs(r1-r2) > GetTolerance('gap'):
                #         print 'Error. Curves not on the same cylinder radius'
                #         return

                if dist > GetTolerance('gap'):
                        print 'Error. Curves not on the same cylinder height'
                        return
                
                l1 = a-origin1
                l2 = b-origin2
                if abs(l1*l2/abs(l1)/abs(l2)-1.0) < 1e-14 or \
                   abs(l1*l2/abs(l1)/abs(l2)+1.0) < 1e-14:
                        theta = pi
                else:
                        theta = acos((l1*l2)/abs(l1)/abs(l2))
#                 if((l1 % l2)[axis] < 0):
#                         theta = - theta
                dt    = theta / 4.0
                dr    = (r2-r1)/4.0

                rotPoint = l1.Clone()
                rotPoint.Normalize()

                # kjetil hack for nasty patches
                controlPoints.extend([crv1[i][0], crv1[i][1], crv1[i][2], 1.0])
                rotPoint.Rotate(nAxis, dt)
                rotPoint = rotPoint / cos(dt)
                
                for j in range(1,5):
                        r = r1 + j*dr
                        
                        w = cos(dt)*(j%2) + (1.0-(j%2))
                        controlPoints.append((origin1[0] + r*rotPoint[0]) * w)
                        controlPoints.append((origin1[1] + r*rotPoint[1]) * w)
                        controlPoints.append((origin1[2] + r*rotPoint[2]) * w)
                        controlPoints.append(w)

                        rotPoint.Rotate(nAxis, dt)
                        if j%2==0:
                                rotPoint = rotPoint / cos(dt)
                        else:
                                rotPoint = rotPoint * cos(dt)
        
        knotXi = [0,0,0,1,1,2,2,2]

        surf   = Surface(3, crv1.GetOrder(), knotXi, crv1.GetKnots(True), controlPoints, True)

        return surf

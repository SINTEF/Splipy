__doc__ = 'Function for connecting curves which live on cylinders'

from math import *
from GoTools import *
from GoTools.CurveFactory import *

gap_epsilon = 1e-4;

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
        listAxis[axis] = 1.0;
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
                # if abs(r1-r2) > gap_epsilon:
                #         print 'Error. Curves not on the same cylinder radius'
                #         return

                if dist > gap_epsilon:
                        print 'Error. Curves not on the same cylinder height'
                        return
                
                l1 = a-origin1;
                l2 = b-origin2;
                if abs(l1*l2/abs(l1)/abs(l2)-1.0) < 1e-14 or \
                   abs(l1*l2/abs(l1)/abs(l2)+1.0) < 1e-14:
                        theta = pi
                else:
                        theta = acos((l1*l2)/abs(l1)/abs(l2))
#                 if((l1 % l2)[axis] < 0):
#                         theta = - theta;
                dt    = theta / 4.0;
                dr    = (r2-r1)/4.0;

                rotPoint = l1.Clone()
                rotPoint.Normalize()

                # kjetil hack for nasty patches
                controlPoints.extend([crv1[i][0], crv1[i][1], crv1[i][2], 1.0])
                rotPoint.Rotate(nAxis, dt)
                rotPoint = rotPoint / cos(dt)
                
                for j in range(1,5):
                        r = r1 + j*dr
                        
                        w = cos(dt)*(j%2) + (1.0-(j%2));
                        controlPoints.append((origin1[0] + r*rotPoint[0]) * w)
                        controlPoints.append((origin1[1] + r*rotPoint[1]) * w)
                        controlPoints.append((origin1[2] + r*rotPoint[2]) * w)
                        controlPoints.append(w);

                        rotPoint.Rotate(nAxis, dt);
                        if j%2==0:
                                rotPoint = rotPoint / cos(dt);
                        else:
                                rotPoint = rotPoint * cos(dt);
        
        knotXi = [0,0,0,1,1,2,2,2]

        surf   = Surface(3, crv1.GetOrder(), knotXi, crv1.GetKnots(True), controlPoints, True)

        return surf;



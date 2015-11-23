__doc__ = 'Implementation of various elementary operations on a per-controlpoint level.'

from GoTools import *
from math import *
import numpy as np

def getNonWeightCP(cp):
    if len(cp) == 4:
        w = cp[3]
        return Point(cp[0]/w, cp[1]/w, cp[2]/w)
    else:
        return cp


def Mirror(obj, point, normal):
    """Mirror a curve, surface or volume around a plane
    @param obj: The obj to mirror
    @type obj: Curve, Surface, Volume, SurfaceModel or VolumeModel
    @param point: The point to mirror about
    @type point: Point or list of floats
    @param normal: The plane normal to mirror about
    @type normal: Point or list of floats
    @return: Mirrored result
    @rtype: Curve, Surface, Volume, SurfaceModel or VolumeModel
    """

    # convert to input to proper usable vectors
    if isinstance(point,list):
        point = Point(list=point)
    if isinstance(normal, list):
        normal = Point(list=normal)
    
    # for input collections, manipulate all objects independently
    if type(obj) is SurfaceModel or type(obj) is VolumeModel:
        ans = []
        for o in obj:
            ans.append(Mirror(o, point, normal))
        if type(obj) is SurfaceModel:
            return SurfaceModel(ans)
        if type(obj) is VolumeModel:
            return VolumeModel(ans)

    normal.Normalize()
    mirroredCP = []

    rational = (len(obj[0]) == 4)

    for cp in obj:

        # fix rational control points
        if rational:
            w = cp[3]

        # compute the actual end points
        oCP = getNonWeightCP(cp) - point
        movLen = abs(oCP*normal) 
        endPos = getNonWeightCP(cp) - 2*movLen*normal

        # fix storage of rational control points
        if rational:
            w = cp[3]
            mirroredCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
        else:
            mirroredCP.append(endPos)
    
    p = obj.GetOrder()
    if type(obj) is Curve:
        knots = obj.GetKnots(with_multiplicities=True)
        return Curve(p, knots, mirroredCP, rational)
    elif type(obj) is Surface:
        knots1, knots2 = obj.GetKnots(with_multiplicities=True)
        return Surface(p[0], p[1], knots1, knots2, mirroredCP, rational)
    elif type(obj) is Volume:
        knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
        return Volume(p[0], p[1], p[2], knots1, knots2, knots3, mirroredCP, rational)

def Translate(obj, vector):
    """Translate a curve, surface or volume 
    @param obj: The obj to translate
    @type obj: Curve, Surface, Volume, SurfaceModel or VolumeModel
    @param vector: The direction to move the object
    @type vector: Point or list of floats
    @return: Translated result
    @rtype: Curve, Surface, Volume, SurfaceModel or VolumeModel
    """
    
    # convert to input to proper usable vectors
    if isinstance(vector,list):
        vector = Point(list=vector)

    # for input collections, manipulate all objects independently
    if type(obj) is SurfaceModel or type(obj) is VolumeModel:
        ans = []
        for o in obj:
            ans.append(Translate(o, vector))
        if type(obj) is SurfaceModel:
            return SurfaceModel(ans)
        if type(obj) is VolumeModel:
            return VolumeModel(ans)

    translatedCP = []

    rational = (len(obj[0]) == 4)
    for cp in obj:

        # fix rational control points
        if rational:
            w = cp[3]

        # compute the actual end points
        endPos = getNonWeightCP(cp) + vector

        # fix storage of rational control points
        if rational:
            w = cp[3]
            translatedCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
        else:
            translatedCP.append(endPos)
    
    p = obj.GetOrder()
    if type(obj) is Curve:
        knots = obj.GetKnots(with_multiplicities=True)
        return Curve(p, knots, translatedCP, rational)
    elif type(obj) is Surface:
        knots1, knots2 = obj.GetKnots(with_multiplicities=True)
        return Surface(p[0], p[1], knots1, knots2, translatedCP, rational)
    elif type(obj) is Volume:
        knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
        return Volume(p[0], p[1], p[2], knots1, knots2, knots3, translatedCP, rational)

def Scale(obj, amount):
    """Scale a curve, surface or volume 
    @param obj: The obj to scaled
    @type obj: Curve, Surface, Volume, SurfaceModel or VolumeModel
    @param amount: The amount to scale the object with
    @type amount: float
    @return: Scaled result
    @rtype: Curve, Surface, Volume, SurfaceModel or VolumeModel
    """

    # for input collections, manipulate all objects independently
    if type(obj) is SurfaceModel or type(obj) is VolumeModel:
        ans = []
        for o in obj:
            ans.append(Scale(o, amount))
        if type(obj) is SurfaceModel:
            return SurfaceModel(ans)
        if type(obj) is VolumeModel:
            return VolumeModel(ans)

    scaledCP = []

    rational = (len(obj[0]) == 4)
    for cp in obj:

        # fix rational control points
        if rational:
            w = cp[3]

        # compute the actual end points
        endPos = getNonWeightCP(cp) * amount

        # fix storage of rational control points
        if rational:
            w = cp[3]
            scaledCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
        else:
            scaledCP.append(endPos)
    
    p = obj.GetOrder()
    if type(obj) is Curve:
        knots = obj.GetKnots(with_multiplicities=True)
        return Curve(p, knots, scaledCP, rational)
    elif type(obj) is Surface:
        knots1, knots2 = obj.GetKnots(with_multiplicities=True)
        return Surface(p[0], p[1], knots1, knots2, scaledCP, rational)
    elif type(obj) is Volume:
        knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
        return Volume(p[0], p[1], p[2], knots1, knots2, knots3, scaledCP, rational)


def rotation_matrix(axis,theta):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
    

def Rotate(obj, normal, theta):
    """Rotate a curve, surface or volume 
    @param obj: The obj to rotate
    @type obj: Curve, Surface, Volume, SurfaceModel or VolumeModel
    @param normal: The normal to rotate object around
    @type normal: float
    @param theta: Angle to rotate, as measured in radians
    @type theta: float
    @return: Scaled result
    @rtype: Curve, Surface, Volume, SurfaceModel or VolumeModel
    """

    # convert to input to proper usable vectors
    if isinstance(normal, list):
        normal = Point(list=normal)

    # for input collections, manipulate all objects independently
    if type(obj) is SurfaceModel or type(obj) is VolumeModel:
        ans = []
        for o in obj:
            ans.append(Rotate(o, normal, theta))
        if type(obj) is SurfaceModel:
            return SurfaceModel(ans)
        if type(obj) is VolumeModel:
            return VolumeModel(ans)

    rotatedCP = []

    R = rotation_matrix(normal, theta);

    rational = (len(obj[0]) == 4)
    for cp in obj:

        # fix rational control points
        if rational:
            w = cp[3]

        # do the actual rotation
        c      = getNonWeightCP(cp) ;
        endPos = Point( R[0][0]*c[0] + R[0][1]*c[1] + R[0][2]*c[2] ,
                        R[1][0]*c[0] + R[1][1]*c[1] + R[1][2]*c[2] ,
                        R[2][0]*c[0] + R[2][1]*c[1] + R[2][2]*c[2] )

        # fix storage of rational control points
        if rational:
            w = cp[3]
            rotatedCP.append(Point(list=[endPos[0]*w, endPos[1]*w, endPos[2]*w, w]))
        else:
            rotatedCP.append(endPos)
    
    p = obj.GetOrder()
    if type(obj) is Curve:
        knots = obj.GetKnots(with_multiplicities=True)
        return Curve(p, knots, rotatedCP, rational)
    elif type(obj) is Surface:
        knots1, knots2 = obj.GetKnots(with_multiplicities=True)
        return Surface(p[0], p[1], knots1, knots2, rotatedCP, rational)
    elif type(obj) is Volume:
        knots1, knots2, knots3 = obj.GetKnots(with_multiplicities=True)
        return Volume(p[0], p[1], p[2], knots1, knots2, knots3, rotatedCP, rational)

def MakeCommonSplineSpace(obj1, obj2):
    """ Make obj1 and obj2 have the same degree, and same knot vector and same 
    rational representation. Uses 'Knot' tolerance (see GoTools.SetTolerance).
    @param obj1: The first object
    @type obj1: Curve or Surface 
    @param obj2: The second object
    @type obj2: Curve or Surface 
    @return: none
    """
    
    tol  = GetTolerance('knot');
    p1   = obj1.GetOrder()
    p2   = obj2.GetOrder()
    rat1 = len(obj1[0]) == 4 # here assuming that obj1 & obj2 lies in 3D space
    rat2 = len(obj2[0]) == 4
    if rat1 or rat2:
        obj1.ForceRational()
        obj2.ForceRational()
    obj1.ReParametrize()  # set knot vectors to span [0,1]
    obj2.ReParametrize()

    if type(obj1) is Curve:
        if p1 > p2:
            obj2.RaiseOrder(p1-p2)
        else:
            obj1.RaiseOrder(p2-p1)
        k1 = obj1.GetKnots(True)
        k2 = obj2.GetKnots(True)
        done = False
        i = 0
        j = 0
        while not done:
            if abs(k1[i] - k2[j]) <= tol:
                i = i+1
                j = j+1
            elif k1[i] < k2[j]:
                obj2.InsertKnot(k1[i])
                i = i+1
            elif k1[i] > k2[j]:
                obj1.InsertKnot(k2[j])
                j = j+1
            if i == len(k1) or j == len(k2):
                done = True


    else:
        if type(obj1) is Surface:
            if p1[0] > p2[0]:
                obj2.RaiseOrder(p1[0]-p2[0], 0)
            else:
                obj1.RaiseOrder(p2[0]-p1[0], 0)
            if p1[1] > p2[1]:
                obj2.RaiseOrder(0, p1[1]-p2[1])
            else:
                obj1.RaiseOrder(0, p2[1]-p1[1])
        elif type(obj1) is Volume:
            if p1[0] > p2[0]:
                obj2.RaiseOrder(p1[0]-p2[0], 0, 0)
            else:
                obj1.RaiseOrder(p2[0]-p1[0], 0, 0)
            if p1[1] > p2[1]:
                obj2.RaiseOrder(0, p1[1]-p2[1], 0)
            else:
                obj1.RaiseOrder(0, p2[1]-p1[1], 0)
            if p1[2] > p2[2]:
                obj2.RaiseOrder(0, 0, p1[2]-p2[2])
            else:
                obj1.RaiseOrder(0, 0, p2[2]-p1[2])

        k1 = obj1.GetKnots(True)
        k2 = obj2.GetKnots(True)
        done = False
        for d in range(len(p1)):
            i = 0
            j = 0
            while not done:
                if abs(k1[d][i] - k2[d][j]) <= tol:
                    i = i+1
                    j = j+1
                elif k1[d][i] < k2[d][j]:
                    obj2.InsertKnot(d, k1[d][i])
                    i = i+1
                elif k1[d][i] > k2[d][j]:
                    obj1.InsertKnot(d, k2[d][i])
                    j = j+1
                if i == len(k1[d]) or j == len(k2[d]):
                    done = True

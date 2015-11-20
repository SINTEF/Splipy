__doc__ = 'Implementation of various refinement schemes.'

from math import atan
from math import pi
from GeoUtils.Knots import *
from GeoUtils.Factory import *

def UniformCurve(curve, n=1):
    """Uniformly refine a curve by inserting n new knots into each knot interval
    @param curve: The curve to refine
    @type  curve: Curve
    @param n    : number of new knots in each interval
    @type  n    : Int
    @return     : None
    """
    
    knots = curve.GetKnots()
    for i in range(0,len(knots)-1):
        for j in range(n):
            newKnot = knots[i] + 1.0*(knots[i+1]-knots[i])/(n+1)*(j+1)
            curve.InsertKnot(newKnot);

# Chop each knot span in half
def UniformSurface(surface, direction=0, n=1):
    """Uniformly refine a surface by inserting n new knots into each knot interval
    @param surface  : The surface to refine
    @type  surface  : Surface 
    @param direction: The direction to refine in (0 = both, 1, 2)
    @type  direction: Int
    @param n        : number of new knots in each interval
    @type  n        : Int
    @return         : None
    """
    knots_u, knots_v = surface.GetKnots()
    if direction == 0 or direction == 1:
        for i in range(0,len(knots_u)-1):
            for j in range(n):
                newKnot = knots_u[i] + 1.0*(knots_u[i+1]-knots_u[i])/(n+1)*(j+1)
                surface.InsertKnot(0, newKnot);
    if direction == 0 or direction == 2:
        for i in range(0,len(knots_v)-1):
            for j in range(n):
                newKnot = knots_v[i] + 1.0*(knots_v[i+1]-knots_v[i])/(n+1)*(j+1)
                surface.InsertKnot(1, newKnot);

# Chop each knot span in half
def UniformVolume(volume, direction=0, n=1):
    """Uniformly refine a volume by inserting n new knots into each knot interval
    @param volume   : The volume to refine
    @type  volume   : Volume 
    @param direction: The direction to refine in (0 = both, 1, 2, 3)
    @type  direction: Int
    @param n        : number of new knots in each interval
    @type  n        : Int
    @return         : None
    """
    knots_u, knots_v, knots_w = volume.GetKnots()
    if direction == 0 or direction == 1:
        for i in range(0,len(knots_u)-1):
            for j in range(n):
                newKnot = knots_u[i] + 1.0*(knots_u[i+1]-knots_u[i])/(n+1)*(j+1)
                volume.InsertKnot(0, newKnot);
    if direction == 0 or direction == 2:
        for i in range(0,len(knots_v)-1):
            for j in range(n):
                newKnot = knots_v[i] + 1.0*(knots_v[i+1]-knots_v[i])/(n+1)*(j+1)
                volume.InsertKnot(1, newKnot);
    if direction == 0 or direction == 3:
        for i in range(0,len(knots_w)-1):
            for j in range(n):
                newKnot = knots_w[i] + 1.0*(knots_w[i+1]-knots_w[i])/(n+1)*(j+1)
                volume.InsertKnot(2, newKnot);


# Geometric distribution of knots
def GeometricRefineCurve(curve, alpha, n):
    """Refine a curve by making a geometric distribution of element sizes
    Consider Curve.FlipParametrization if you need refinement towards the other edge
    @param curve: The curve to refine
    @type  curve: Curve 
    @param alpha: The ratio between two sequential knot segments
    @type  alpha: Float
    @param n    : The number of knots to insert
    @type  n    : Int
    @return     : None
    """
    
    # some error tests on input
    if n<=0:
        print 'n should be greater than 0'
        return None

    # fetch knots
    knots = curve.GetKnots()
    knotStart = knots[0]
    knotEnd   = knots[-1]
    dk = knotEnd - knotStart

    # evaluate the factors
    n = n+1 # redefine n to be knot spans instead of new internal knots
    totProd = 1.0
    totSum  = 0.0
    for i in range(n):
        totSum  += totProd
        totProd *= alpha
    d1 = 1.0 / totSum
    knot = d1

    # do the actual knot insertion
    oldKnots = curve.GetKnots()
    for i in range(n-1):
        newKnot = knotStart + knot*dk
        if not KnotExist(oldKnots, newKnot):
            curve.InsertKnot(newKnot)
        knot += alpha*d1
        d1   *= alpha


# Geometric distribution of knots
def GeometricRefineSurface(surface, direction, alpha, n):
    """Refine a surface by making a geometric distribution of element sizes
    Consider Surface.FlipParametrization if you need refinement towards the other edge
    @param surface  : The surface to refine
    @type  surface  : Surface 
    @param direction: The direction to refine in (u=1 or v=2) 
    @type  direction: Int
    @param alpha    : The ratio between two sequential knot segments
    @type  alpha    : Float
    @param n        : The number of knots to insert
    @type  n        : Int
    @return         : None
    """
    
    # some error tests on input
    if n<=0:
        print 'n should be greater than 0'
        return None

    flipBack = False
    if direction < 0:
        surface.FlipParametrization(-direction-1)
        direction = -direction
        flipBack = True

    # fetch knots
    knots_u, knots_v = surface.GetKnots()
    if direction == 1:
        knotStart = knots_u[0]
        knotEnd   = knots_u[-1]
    elif direction == 2:
        knotStart = knots_v[0]
        knotEnd   = knots_v[-1]
    else:
        print 'Direction should be 1 or 2'
        return None
    dk = knotEnd - knotStart

    # evaluate the factors
    n = n+1 # redefine n to be knot spans instead of new internal knots
    totProd = 1.0
    totSum  = 0.0
    for i in range(n):
        totSum  += totProd
        totProd *= alpha
    d1 = 1.0 / totSum
    knot = d1

    # do the actual knot insertion
    oldKnots = surface.GetKnots()[direction-1]
    for i in range(n-1):
        newKnot = knotStart + knot*dk
        if not KnotExist(oldKnots, newKnot):
            surface.InsertKnot(direction-1, newKnot)
        knot += alpha*d1
        d1   *= alpha

    if flipBack:
        surface.FlipParametrization(direction-1)

# Geometric distribution of knots
def GeometricRefineVolume(volume, direction, alpha, n):
    """Refine a volume by making a geometric distribution of element sizes
    Consider Volume.FlipParametrization if you need refinement towards the other edge
    @param volume   : The volume to refine
    @type  volume   : Volume 
    @param direction: The direction to refine in (u=1, v=2 or w=3) 
    @type  direction: Int
    @param alpha    : The ratio between two sequential knot segments
    @type  alpha    : Float
    @param n        : The number of knots to insert
    @type  n        : Int
    @return         : None
    """
    
    # some error tests on input
    if n<=0:
        print 'n should be greater than 0'
        return None

    flipBack = False
    if direction < 0:
        volume.FlipParametrization(-direction-1)
        direction = -direction
        flipBack = True

    # fetch knots
    knots_u, knots_v, knots_w = volume.GetKnots()
    if direction == 1:
        knotStart = knots_u[0]
        knotEnd   = knots_u[-1]
    elif direction == 2:
        knotStart = knots_v[0]
        knotEnd   = knots_v[-1]
    elif direction == 3:
        knotStart = knots_w[0]
        knotEnd   = knots_w[-1]
    else:
        print 'Direction should be 1, 2 or 3' 
        return None
    dk = knotEnd - knotStart

    # evaluate the factors
    n = n+1 # redefine n to be knot spans instead of new internal knots
    totProd = 1.0
    totSum  = 0.0
    for i in range(n):
        totSum  += totProd
        totProd *= alpha
    d1 = 1.0 / totSum
    knot = d1

    # do the actual knot insertion
    oldKnots = volume.GetKnots()[direction-1]
    for i in range(n-1):
        newKnot = knotStart + knot*dk
        if not KnotExist(oldKnots, newKnot):
            volume.InsertKnot(direction-1, newKnot)
        knot += alpha*d1
        d1   *= alpha

    if flipBack:
        volume.FlipParametrization(direction-1)


# Edge refinement
def EdgeRefineSurface(surface, direction, S, n):
    """Refine a surface by both edges, by sampling a atan-function
    @param surface  : The surface to refine
    @type  surface  : Surface 
    @param direction: The direction to refine in (u=1 or v=2) 
    @type  direction: Int
    @param S        : The slope of the atan-function
    @type  S        : Float
    @param n        : The number of knots to insert
    @type  n        : Int
    @return         : None
    """
    
    # some error tests on input
    if n<=0:
        print 'n should be greater than 0'
        return None
    
    # fetch knots
    knots_u, knots_v = surface.GetKnots()
    if direction == 1:
        knotStart = knots_u[0]
        knotEnd   = knots_u[-1]
    elif direction == 2:
        knotStart = knots_v[0]
        knotEnd   = knots_v[-1]
    else:
        print 'Direction should be 1 or 2'
        return None
    dk = knotEnd - knotStart
    
    # evaluate the factors
    newKnots = []
    maxAtan  = atan(S)
    for i in range(1,n+1):
        xi  = -1.0 + 2.0*i/(n+1)
        xi *= S
        newKnots.append(knotStart + (atan(xi)+maxAtan)/2/maxAtan*dk)
    
    # do the actual knot insertion
    oldKnots = surface.getKnots()[direction-1]
    for x in newKnots:
        if not KnotExist(oldKNots, x):
            surface.InsertKnot(direction-1, x)

def _splitvector(len, parts):
    delta = len/(parts)
    sizes = [delta for i in range(parts)]
    remainder = len-parts*delta
    for i in range(parts-remainder+1, parts):
        sizes[i] = sizes[i]+1
    result = [0]
    for i in range(1,parts):
        result.append(sizes[i]+result[i-1])
    return result

def SubdivideSurface(srfs, subdivisions):
    """Subdivide a list of surfaces
       The X-Y plane is subdivided in subdivisions+1 blocks.
       @param vols        : Volumes/Surfaces to split
       @type  vols        : List of Volume/Surfaces
       @param subdivisions: Number of subdivisions to perform in the X-Y plane
       @type  subdivisions: Int
       @return            : List of new surfaces
    """

    # Subdivisions
    if subdivisions > 0:
      volsx = []
      for v in srfs:
        v.ReParametrize()
        knotx, knoty = v.GetKnots()
        splits = _splitvector(len(knotx), subdivisions+1)
        splitval = []
        for i in range(1,len(splits)):
          splitval.append(knotx[splits[i]])
        volstmp = [v,v]
        for val in splitval:
          volstmp = SplitSurface(volstmp[1], 0, val)
          volsx.append(volstmp[0])
          if val == splitval[-1]:
            volsx.append(volstmp[1])

      volsy = []
      for v in volsx:
        v.ReParametrize()
        knotx, knoty = v.GetKnots()
        splits = _splitvector(len(knoty), subdivisions+1)
        splitval = []
        for i in range(1,len(splits)):
          splitval.append(knoty[splits[i]])
        volstmp = [v,v]
        for val in splitval:
          volstmp = SplitSurface(volstmp[1], 1, val)
          volsy.append(volstmp[0])
          if val == splitval[-1]:
            volsy.append(volstmp[1])
    else:
      volsy = srfs

    return volsy

def SubdivideVolume(vols, subdivisions, zsplit):
    """Subdivide a list of volumes
       First the X-Y plane is subdivided in subdivisions+1 blocks.
       Then the result is split in zsplit+1 parts in the z direction.
       @param vols        : Volumes to split
       @type  vols        : List of Volume
       @param subdivisions: Number of subdivisions to perform in the X-Y plane
       @type  subdivisions: Int
       @param zsplit      : Number of times to split in the z direction
       @type  zsplit      : Int
       @return            : List of new volumes
    """

    # Subdivisions
    if subdivisions > 0:
      volsx = []
      for v in vols:
        v.ReParametrize()
        knotx, knoty, knotz = v.GetKnots()
        splits = _splitvector(len(knotx), subdivisions+1)
        splitval = []
        for i in range(1,len(splits)):
          splitval.append(knotx[splits[i]])
        volsx += v.Split(splitval, 0)

      volsy = []
      for v in volsx:
        v.ReParametrize()
        knotx, knoty, knotz = v.GetKnots()
        splits = _splitvector(len(knoty), subdivisions+1)
        splitval = []
        for i in range(1,len(splits)):
          splitval.append(knoty[splits[i]])
        volsy += v.Split(splitval, 1)
    else:
      volsy = vols

    if zsplit > 0:
      volsz = []
      for v in volsy:
        v.ReParametrize()
        knotx, knoty, knotz = v.GetKnots()
        splits = _splitvector(len(knotz), zsplit+1)
        splitval = []
        for i in range(1,len(splits)):
          splitval.append(knotz[splits[i]])
        volsz += v.Split(splitval, 2)
    else:
      volsz = volsy

    return volsz

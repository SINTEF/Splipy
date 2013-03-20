from GoTools import * 
import os
import numpy as np
import math

def distanceFunction(surfaces, wallset, patches=[]):
    """
    Calculate shortest wall distance
        patches = list of all patches for which the wall distance is to be calculated
        wallsrfid = info about wall (see InputFile.GetTopologSet)
    """
    
    worksurfaces = []
    if len(patches):
      for patch in patches:
          worksurfaces.append(surfaces[patch-1])
    else:
      worksurfaces = surfaces
    
    wallsurfaces = []
    for idx in wallset:
        wallsurfaces.append(surfaces[idx-1])
    
    wallcurves = []
    
    for idx in wallset:
      for edge in wallset[idx].edge:
        wallcurves.append(getWallCurve(surfaces[idx-1], edge-1))

    D4 = calcDistQuadMinNewton(wallcurves, worksurfaces)
    
    return D4
    
def calcDistStraightLine(worksurfaces):
    """
    Calculate minimum wall distance based on straight line approximation (assumes perpendicular grid lines)
    """
    
    D = []
    
    for surface in worksurfaces:
        wallcurve = getWallCurve(surface)
        (knots_xi, knots_eta) = surface.GetKnots()
        
        wdist = np.zeros((len(knots_xi), len(knots_eta)))
        
        wcrvpts = getCurvePoints(surface, knots_xi, knots_eta[0])
        
        j = 0
        
        for knot_eta in knots_eta:
            crvpts = getCurvePoints(surface, knots_xi, knot_eta)
            i = 0
            for (wallpt, pt) in zip(wcrvpts, crvpts):
                wdist[i,j] = calcPtsDistance(wallpt, pt)
                i = i + 1
            j = j+1
        D.append(wdist)
    
    return D



def calcDistQuadMinNewton(wallcurves, worksurfaces):
    """
    Calculate minimum wall distance based on a combination of quadratic minimization and Newton's method
    """
    
    D = []
    srfID = 0
    for surface in worksurfaces:
        print 'Working on surface number ' + str(srfID+1)
        
        D_surface = []
        
        (knots_xi, knots_eta) = surface.GetKnots()
        
        wdist = np.zeros((len(knots_xi), len(knots_eta)))
        s = np.zeros(4)
        curveID = 0
        
        for curve in wallcurves:
            print 'Calculating shortest wall distance from wallcurve ' + str(curveID+1)
            crv_knots_xi = curve.GetKnots()
            
            i = 0
            
            length = len(crv_knots_xi)
            #s[0] = crv_knots_xi[int(length*0.1)]
            tmp_Newton = np.infty
            tmp = crv_knots_xi[0]
            
            for knot_xi in knots_xi:
                j = 0
                for knot_eta in knots_eta:
                    
                    pt = surface.Evaluate(knot_xi, knot_eta)
                    k = 0

                    s[0] = crv_knots_xi[int(length*0.1)]
                    s[1] = crv_knots_xi[int(length*0.5)]
                    s[2] = crv_knots_xi[int(length*0.9)]
                    
                    while True:
                        m = 0
                        while m <= 4:
                            s[3] = calcQuadMinParam(curve, pt, s[0], s[1], s[2], tmp)
                            tol = np.absolute(s[3] - tmp)
                            tmp = s[3]
                        
                            s = determineSs(curve, pt, s)
                            m = m+1

                        tmp = s[0]
                    
                        s_star = calcSingleNewton(curve, crv_knots_xi, tmp, pt) 

                        tol = np.absolute(s_star - tmp_Newton)
                        tmp_Newton = s_star
                        tmp = tmp_Newton
                        
                        tolerance = calcTolerance(crv_knots_xi, tmp_Newton, 1e-10)

                        if tol < tolerance:
                            wdist[i,j] = calcPtsDistance(curve.Evaluate(tmp_Newton), pt)
                            j = j + 1
                            break;
                            
                        if k > 100:  # Fix slow or no convergence
                            tmp = tmp*1.1
                            s = s*0.9
                            k = 0
                        k = k+1

                i =  i + 1
            curveID = curveID + 1
            D_surface.append(wdist)
        D.append(determineShortestDistance(D_surface))
        srfID = srfID + 1
    return D
    
def calcDistQuadMin(wallcurves, worksurfaces):
    """
    Calculate minimum wall distance based on quadratic minimization
    """
    
    D = []
    
    for surface in worksurfaces:
        D_surface = []
        
        s = np.zeros(4)
        
        (knots_xi, knots_eta) = surface.GetKnots()
        wdist = np.zeros((len(knots_xi), len(knots_eta)))
        
        for curve in wallcurves:
            crv_knots_xi = curve.GetKnots()
            length = len(crv_knots_xi)
            i = 0
            for knot_xi in knots_xi:
                j = 0
                s[0] = crv_knots_xi[int(length*0.1)]
                s[1] = crv_knots_xi[int(length*0.5)]
                s[2] = crv_knots_xi[int(length*0.9)]
                for knot_eta in knots_eta:
                    pt = surface.Evaluate(knot_xi, knot_eta)
                    
                    if j > 0:
                        idx = (np.abs(np.array(crv_knots_xi)-tmp)).argmin()
                        if idx == len(crv_knots_xi)-1:
                            s[0] = crv_knots_xi[idx-1]
                            s[1] = (crv_knots_xi[idx-1]  + crv_knots_xi[idx])/2
                            s[2] = crv_knots_xi[idx]
                        else:
                            s[0] = crv_knots_xi[idx]
                            s[1] = (crv_knots_xi[idx]  + crv_knots_xi[idx+1])/2
                            s[2] = crv_knots_xi[idx+1]
                    k = 1
                    tmp = np.infty
                    while True:
                        s[3] = calcQuadMinParam(curve, pt, s[0], s[1], s[2], tmp)
                        tol = np.absolute(s[3] - tmp)
                        tmp = s[3]
                        
                        tolerance = calcTolerance(crv_knots_xi, s[3], 1e-10)
                        
                        if tol < tolerance:
                            wdist[i,j] = calcPtsDistance(curve.Evaluate(s[3]), pt)
                            #print 'Converged ' + str(i) + ',  ' + str(j) 
                            j = j + 1
                            break;
                        
                        s = determineSs(curve, pt, s)
                        k = k+1
                        if k > 100:
                            # Modify in case of slow or no convergence
                            s[2] = s[1]
                            s[1] = s[0]
                            k = 0
                i =  i + 1
            D_surface.append(wdist)
        D.append(determineShortest(D_surface))
    return D


def calcDistMinNewton(wallcurves, worksurfaces):
    """
    Calculate minimum wall distance based on Newton's method
    """
    
    
    D = []
    
    for surface in worksurfaces:
        D_surface = []
        
        (knots_xi, knots_eta) = surface.GetKnots()
        
        wdist = np.zeros((len(knots_xi), len(knots_eta)))
        s = np.zeros(4)
        
        curveID = 0
        
        for curve in wallcurves:
            crv_knots_xi = curve.GetKnots()
            
            i = 0
            print 'Calculating shortest wall distance for wallcurve ' + str(curveID)
            
            length = len(crv_knots_xi)
            
            tmp = crv_knots_xi[0]
            
            for knot_xi in knots_xi:
                j = 0
                for knot_eta in knots_eta:
                    
                    pt = surface.Evaluate(knot_xi, knot_eta)
                    k = 0
    
                    while True:
                        s_star = calcSingleNewton(curve, crv_knots_xi, tmp, pt) 

                        tol = np.absolute(s_star - tmp)
                        tmp = s_star
                        
                        tolerance = calcTolerance(crv_knots_xi, tmp, 1e-10)
                        if tol < tolerance:
                            print str(i) + ', ' + str(j) + ' converged!'
                            wdist[i,j] = calcPtsDistance(curve.Evaluate(tmp), pt)
                            j = j + 1
                            break;
                            
                        if k > 50:  # Fix slow or no convergence
                            tmp = crv_knots_xi[int(length*np.random.rand())]
                        k = k+1


                i =  i + 1
            D_surface.append(wdist)
            curveID = curveID + 1
        D.append(determineShortestDistance(D_surface))
    return D


def calcSingleNewton(curve,knots,s, pt):
    """
    Do a single Newton iteration
    """
    
    (pt_x, pt_y) = parsePoint(pt)

    width = 0.0001
    
    D1 = calcPtsDistance(curve.Evaluate(s-width), pt)**2
    D2 = calcPtsDistance(curve.Evaluate(s), pt)**2
    D3 = calcPtsDistance(curve.Evaluate(s+width), pt)**2
    
    D_p = (D3-D1)/(2*width)
    D_pp = (D3-2*D2+D1)/(width**2)

    s_star = s - D_p/D_pp

    return limitS(s_star, knots[0], knots[-1])

    

def determineSs(wallcurve, pt, s,flag=False):
    P = np.zeros(4)

    for i in xrange(0,4):
        P[i] = calcDistPtParam(wallcurve, pt, s[i])

    idx = P.argmax()

    s[idx] = np.infty
    s.sort()
    
    return s

def limitS(s, smin, smax):

    if s < smin:
        s = smin
    elif s > smax:
        s = smax
    
    return s
    
def calcTolerance(knots, s_star, tolerance):
    
    idx = (np.abs(np.array(knots)-s_star)).argmin()
    
    tol = np.infty
    
    if knots[idx] <= s_star:
        if idx == len(knots)-1:
            tol = abs(knots[idx] - knots[idx-1])*tolerance
        else:
            tol = abs(knots[idx+1] - knots[idx])*tolerance
    elif knots[idx] >= s_star:
        if idx == 0:
            tol = abs(knots[idx+1] - knots[idx])*tolerance
        else:
            tol = abs(knots[idx] - knots[idx-1])*tolerance
    
    return tol
                    
def calcDistPtParam(wallcurve, pt, s):
    
    pt = parsePoint(pt)
    
    curvept = parsePoint(wallcurve.Evaluate(s))
    
    return (curvept[0] - pt[0])**2 + (curvept[1] - pt[1])**2
    

def calcQuadMinParam(wallcurve, pt, s1, s2, s3, s_star_last):
    y12 = s1**2 - s2**2; y23 = s2**2 - s3**2; y31 = s3**2 - s1**2
    s12 = s1 - s2; s23 = s2 - s3; s31 = s3 - s1
    
    D1 = calcDistPtParam(wallcurve, pt, s1)**2
    D2 = calcDistPtParam(wallcurve, pt, s2)**2
    D3 = calcDistPtParam(wallcurve, pt, s3)**2

    
    if abs(s23*D1 + s31*D2 + s12*D3) > 1e-20:
        s_star = 1.0/2* (y23*D1 + y31*D2 + y12*D3) / (s23*D1 + s31*D2 + s12*D3)
    else:
        s_star = np.mean([s1, s2, s3])
    
    knots = wallcurve.GetKnots()
    
    return limitS(s_star, knots[0], knots[-1])

def getWallCurve(surface, idx):
    edges = surface.GetEdges()
    return edges[idx]
    
def parsePoint(point):
    #3D point: 0.254421 -0.0743915 0
    strpoint = str(point).split()
    
    return (float(strpoint[2]), float(strpoint[3]))
    
def getCurvePoints(surface,knots_xi, knot_eta):
    crvpoints = []
    for knot_xi in knots_xi:
        crvpoints.append(parsePoint(surface.Evaluate(knot_xi, knot_eta)))
    
    return np.array(crvpoints)
        
def calcPtsDistance(pt1, pt2):
    """
    Calculate shortest distance between two points
    """
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2)

def determineShortestDistance(Ds):
    """
    Determine which array gives the shortes distance to the wall, based on sum of wall distances
    """
    
    dist = []

    for D in Ds:
        dist.append(np.sum(D))
    
    idx = np.array(dist).argmin()
    
    return Ds[idx]

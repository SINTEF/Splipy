from GoTools import *
import os
import numpy as np
import math
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize

def distanceFunction(surfaces, wallset, patches=[]):
    """Calculate shortest wall distance
       @param surfaces: Surfaces in model
       @type surfaces: List of Surface
       @param wallset: The edges defining the wall
       @type wallset: List of tuple of (patch, edge) numbers
       @param: Optional list of patches to process
       @type patches: List of integer
       @return: Coefficients of distance field
       @rtype: List of list of float
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
        wallcurves.append(getWallCurve(surfaces[idx-1], edgeNumber(edge-1)))

    D = calcDistScipy(wallcurves, worksurfaces)

    return D


def calcDistScipy(wallcurves, worksurfaces):
    """Calculate minimum wall distance using scipy and minimize scalar
       @param wallcurves: List of curves describing the wall
       @type wallcurves: List of curve
       @param worksurfaces: Surfaces to process
       @type worksurfaces: List of Surface
       @return Wall distance
       @rtype List of doubles
    """

    D = []
    srfID = 0
    for surface in worksurfaces:
        print 'Working on surface number ' + str(srfID+1)

        (knots_xi, knots_eta) = surface.GetKnots()
        Dp = [0]*len(knots_xi)*len(knots_eta)

        i = 0
        for knot_xi in knots_xi:
            j = 0
            for knot_eta in knots_eta:
                pt = surface.Evaluate(knot_xi, knot_eta)

                mindist = np.infty

                for curve in wallcurves:
                    curveID = 0
                    crv_knots_xi = curve.GetKnots()
                    s0 = (crv_knots_xi[-1] + crv_knots_xi[0])/2.0
                    #length = len(crv_knots_xi)
                    lb = crv_knots_xi[0]
                    ub = crv_knots_xi[-1]
                    res = minimize_scalar(calcPtsDistanceCurve, s0, args=(curve, pt), bounds=(lb,ub), method='bounded', options={'xtol': 1e-10, 'disp': False})
                    tmp = calcPtsDistance(curve.Evaluate(res.x), pt)
                    if tmp < mindist:
                        mindist = tmp
                Dp[j*len(knots_xi)+i] = mindist
                j = j+1
            i = i+1
        curveID = curveID + 1
        D.append(surface.Interpolate(Dp))
        srfID = srfID + 1
    return D


def getWallCurve(surface, idx):
    edges = surface.GetEdges()
    return edges[idx]

def getWallFace(volume, idx):
    faces = volume.GetFaces()
    return faces[idx]

def edgeNumber(edgeID):
    """
    Convert edge number to parametric order
    """

    if edgeID == 0:
        return 3
    elif edgeID == 1:
        return 1
    elif edgeID == 2:
        return 0
    elif edgeID == 3:
        return 2

def faceNumber(edgeID):
    """
    Convert face number to parametric order
    """

    remap = [1,2,5,6,3,4]
    return remap[edgeID]-1

def calcPtsDistance(pt1, pt2):
    """
    Calculate shortest distance between two points
    """
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2)

def calcPtsDistanceCurve(s, curve, pt2):
    """
    Calculate shortest distance between two points
    """
    pt1 = curve.Evaluate(s)
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2)

def calcPtsDistance3D(pt1, pt2):
    """
    Calculate shortest distance between two points
    """
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2)

def calcPtsDistanceSurface(s, surface, pt2):
    """
    Calculate shortest distance between two points
    """
    pt1 = surface.Evaluate(s[0], s[1])
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2)

def distanceFunction3D(volumes, wallset, patches=[]):
    """Calculate shortest wall distance
       @param volumes: Volumes in model
       @type volumes: List of Volume
       @param wallset: The faces defining the wall
       @type wallset: List of tuple of (patch, face) numbers
       @param: Optional list of patches to process
       @type patches: List of integer
       @return: Coefficients of distance field
       @rtype: List of list of float
    """

    workvolumes = []
    if len(patches):
      for patch in patches:
          workvolumes.append(volumes[patch-1])
    else:
      workvolumes = volumes

    wallvolumes = []
    for idx in wallset:
        wallvolumes.append(volumes[idx-1])

    wallfaces = []

    for idx in wallset:
      for face in wallset[idx].face:
        wallfaces.append(getWallFace(volumes[idx-1], faceNumber(face-1)))

    D = calcDistScipy3D(wallfaces, workvolumes)

    return D


def calcDistScipy3D(wallfaces, workvolumes):
    """Calculate minimum wall distance using scipy and minimize scalar
       @param wallcurves: List of faces describing the wall
       @type wallcurves: List of Surface
       @param workvolumes: Volumes to process
       @type workvolumes: List of Volume
       @return Wall distance
       @rtype List of doubles
    """

    D = []
    volID = 1
    for volume in workvolumes:
        print 'Working on volume number ' + str(volID)

        (knots_xi, knots_eta, knots_gamma) = volume.GetGreville()
        Dp = [0]*len(knots_gamma)*len(knots_eta)*len(knots_xi)

        i = 0
        for knot_xi in knots_xi:
            j = 0
            for knot_eta in knots_eta:
                k = 0
                for knot_gamma in knots_gamma:
                    pt = volume.Evaluate(knot_xi, knot_eta, knot_gamma)

                    mindist = np.infty

                    for face in wallfaces:
                        (crv_knots_xi, crv_knots_eta) = face.GetKnots()
                        s0 = ((crv_knots_xi[-1] + crv_knots_xi[0])/2.0,
                              (crv_knots_eta[-1] + crv_knots_eta[0])/2.0)
                        lb = (crv_knots_xi[0], crv_knots_xi[-1])
                        ub = (crv_knots_eta[0], crv_knots_eta[-1])
                        res = minimize(_calcPtsDistanceSurface, s0, args=(face, pt), bounds=(lb,ub), method='SLSQP', jac=False)
                        tmp = _calcPtsDistance3D(face.Evaluate(res.x[0], res.x[1]), pt)
                        if tmp < mindist:
                            mindist = tmp
                    Dp[k*len(knots_xi)*len(knots_eta)+j*len(knots_xi)+i] = mindist
                    k = k+1
                j = j+1
            i = i+1
        D.append(volume.Interpolate(Dp))
        volID = volID + 1
    return D

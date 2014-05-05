__doc__ = 'Implementation of transfinite interpolation algorithms.'

from CurveUtils import *
from GoTools import *
from GoTools.SurfaceFactory import *
from GoTools.CurveFactory import *

def TFIcurve(curve, j1=1, r=0.9):
  pts1 = GetCurvePoints(curve[0])
  pts2 = GetCurvePoints(curve[1])
  pts3 = GetCurvePoints(curve[2])
  pts4 = GetCurvePoints(curve[3])

  imax = len(pts1)
  jmax = len(pts3)

  ximin  = CurveLengthParametrization(pts1,True)
  ximax  = CurveLengthParametrization(pts2,True)

  etamin2 = CurveLengthParametrization(pts3,True)
  etamax2 = CurveLengthParametrization(pts4,True)

  knots3 = curve[2].GetKnots()
  knots4 = curve[3].GetKnots()
  etamin = []
  etamax = []
  smin = etamin2[jmax-1]-etamin2[j1-1]
  smax = etamax2[jmax-1]-etamax2[j1-1]
  for j in range(j1-1,jmax):
    etamin.append((etamin2[j]-etamin2[j1-1])/smin)
    etamax.append((etamax2[j]-etamax2[j1-1])/smax)

  pts  = []
  pts.append(pts3[j1])
  for i in range(1,imax-1):
    p_ij = 1.0 - (ximax[i]-ximin[i])*(etamax[1]-etamin[1])
    bxi  = (ximin[i] + etamin[1]*(ximax[i]-ximin[i]))/p_ij
    beta = (etamin[1] + ximin[i]*(etamax[1]-etamin[1]))/p_ij

    pt = (1.0-beta)*pts1[i] + beta*pts2[i] + (1-bxi)*pts3[j1] + bxi*pts4[j1]
    pt = pt - (1.0-beta)*(1.0-bxi)*pts1[0] - beta*(1.0-bxi)*pts3[-1]
    pt = pt - (1.0-beta)*bxi*pts1[-1] - bxi*beta*pts4[-1]
    pts.append(pt)

  pts.append(pts4[j1])

  knots = CurveLengthParametrization(pts,False)
  curve  = InterpolateCurve(pts,knots)
  return pts, curve

def TFISurface(curve):
  """Construct a surface from four edge curves by using transfinite interpolation
     @param curves: Surface edges
     @type curves: List of Curve
     @param perodic1: Whether or not surface is periodic in first direction
     @type periodic1: Boolean
     @param periodic2: Whether or not surface is periodic in second direction
     @type periodic2: Boolean
     @return: New surface
     @rtype: Surface
  """
  # Knots without multiplicity
  knotsu = curve[0].GetKnots()
  knotsv = curve[1].GetKnots()

  # Span of parameter spaces
  ulen = knotsu[-1]-knotsu[0]
  vlen = knotsv[-1]-knotsv[0]

  # Points for transfinite interpolation
  # Corner points
  P1 = curve[0].Evaluate(knotsu[0])
  P2 = curve[1].Evaluate(knotsu[-1])
  P3 = curve[2].Evaluate(knotsu[0])
  P4 = curve[3].Evaluate(knotsu[-1])

  pts = []
  for eta in knotsv:
    c2 = curve[1].Evaluate(eta)
    c4 = curve[3].Evaluate(eta)

    v = (eta-knotsv[0])/vlen
    for xi in knotsu:
      c1 = curve[0].Evaluate(xi)
      c3 = curve[2].Evaluate(xi)

      u = (xi-knotsu[0])/ulen
      c = (1.0-v)*c1 + v*c3 + (1.0-u)*c2 + u*c4
      c = c - (1.0-u)*(1.0-v)*P2 + u*v*P3 - u*(1.0-v)*P1 - (1.0-u)*v*P4
      pts.append(c)

  # Uniform knot vectors
  dxi = 1.0/(len(knotsu)-1)
  knots1 = [i*dxi for i in range(len(knotsu))]
  deta = 1.0/(len(knotsv)-1)
  knots2 = [i*deta for i in range(len(knotsu))]

  surface = InterpolateSurface(pts, knots1, knots2)
  return surface

def TFIGradedSurface(curve, periodic1 = False, periodic2 = False):
  """Construct a graded surface from four edge curves by using transfinite interpolation
     @param curve: Surface edges
     @type curve: List of Curve
     @param perodic1: Whether or not surface is periodic in first direction
     @type periodic1: Boolean
     @param periodic2: Whether or not surface is periodic in second direction
     @type periodic2: Boolean
     @return: New surface
     @rtype: Surface
  """
#     - 4 4 4 4 4 4 -
#     1             2
#     1             2
#     1             2
#     1             2
#     - 3 3 3 3 3 3 -

  # Smoothing factor for grid size
  factor = 0.8

  # Get knot vectots
  knots1 = curve[0].GetKnots()
  knots2 = curve[1].GetKnots()
  knots3 = curve[2].GetKnots()
  knots4 = curve[3].GetKnots()

  # Get points
  x1 = []
  y1 = []
  pts1 = []
  for xi in knots1:
    pt = curve[0].Evaluate(xi)
    pts1.append(pt)
    x1.append(pt[0])
    y1.append(pt[1])
  x2 = []
  y2 = []
  pts2 = []
  for xi in knots2:
    pt = curve[1].Evaluate(xi)
    pts2.append(pt)
    x2.append(pt[0])
    y2.append(pt[1])
  x3 = []
  y3 = []
  pts3 = []
  for eta in knots3:
    pt = curve[2].Evaluate(eta)
    pts3.append(pt)
    x3.append(pt[0])
    y3.append(pt[1])
  x4 = []
  y4 = []
  pts4 = []
  for eta in knots4:
    pt = curve[3].Evaluate(eta)
    pts4.append(pt)
    x4.append(pt[0])
    y4.append(pt[1])

  # Curve length parametrization
  ximin  = CurveLengthParametrization(pts1,True)
  ximax  = CurveLengthParametrization(pts2,True)
  etamin = CurveLengthParametrization(pts3,True)
  etamax = CurveLengthParametrization(pts4,True)

  x01 = x1
  y01 = y1
  x02 = x2
  y02 = y2
  x03 = x3
  y03 = y3
  x04 = x4
  y04 = y4

  imax = len(x1)
  jmax = len(x3)

  xvec = []
  yvec = []
  pts  = []
  for j in range(0,jmax):
    x01m = x01
    y01m = y01
    x02m = x02
    y02m = y02
    for i in range(0,imax):
      x03m = x03
      y03m = y03
      x04m = x04
      y04m = y04

      p_ij = 1.0 - (ximax[i]-ximin[i])*(etamax[j]-etamin[j])
      bxi  = (ximin[i] + etamin[j]*(ximax[i]-ximin[i]))/p_ij
      beta = (etamin[j] + ximin[i]*(etamax[j]-etamin[j]))/p_ij

      x = (1.0-beta)*x01m[i] + beta*x02m[i] + (1-bxi)*x03m[j] + bxi*x04m[j]
      x = x - (1.0-beta)*(1.0-bxi)*x01m[0] - beta*(1.0-bxi)*x03m[-1]
      x = x - (1.0-beta)*bxi*x01m[-1] - bxi*beta*x04m[-1]
      xvec.append(x)

      y = (1.0-beta)*y01m[i] + beta*y02m[i] + (1-bxi)*y03m[j] + bxi*y04m[j]
      y = y - (1.0-beta)*(1.0-bxi)*y01m[0] - beta*(1.0-bxi)*y03m[-1]
      y = y - (1.0-beta)*bxi*y01m[-1] - bxi*beta*y04m[-1]
      yvec.append(y)

      pt = Point(list=[x,y,0.0])
      pts.append(pt)

  # Curve length parametrization
  xi = []
  s = 0.0
  xi.append(s)
  for i in range(0,imax-1):
    ds = abs(pts1[i+1]-pts1[i])
    s = s + ds
    xi.append(s)

  eta = []
  s = 0.0
  eta.append(s)
  for j in range(0,jmax-1):
    ds = abs(pts3[j+1]-pts3[j])
    s = s + ds
    eta.append(s)

  # Interpolate points
  surface  = InterpolateSurface(pts,xi,eta)
  return surface

def _GetCurvePoints(curve):
  knots = curve.GetKnots()

  pts = []
  for xi in knots:
    pts.append(curve.Evaluate(xi))

  return pts

def _GetTangents(curve, normalize = False):
  tangents = []
  knots = curve.GetKnots()
  for xi in knots:
    tau = curve.EvaluateTangent(xi)
    if (normalize):
      tau = tau/abs(tau)
    tangents.append(tau)

  return tangents

def _GetNormals(curve,righthandsystem=True,normalize=False):
  tangents = _GetTangents(curve,normalize)
  normals = []
  for t in tangents:
    if (righthandsystem):
      normals.append(Point(list=[-t[1],t[0],0.0]))
    else:
      normals.append(Point(list=[t[1],-t[0],0.0]))

  return normals

def _OrthogonalCurve(curve1,curve2,curve3,curve4,j=1):
  pts1 = GetCurvePoints(curve1)
  pts2 = GetCurvePoints(curve2)
  pts3 = GetCurvePoints(curve3)
  pts4 = GetCurvePoints(curve4)

  n1 = len(pts1)

  dy1 = abs(pts3[j]-pts3[j-1])
  dy2 = abs(pts4[j]-pts4[j-1])
  s1  = abs(pts3[-1]-pts3[j-1])
  s2  = abs(pts4[-1]-pts4[j-1])

  sn1 = CurveLengthParametrization(pts1,True)
  normal1 = _GetNormals(curve1,True,True)

  pts = []
  pts.append(pts3[j])
  for i in range(1,n1-1):
    dy = sn1[i]*dy2 + (1.0-sn1[i])*dy1
    pts.append(pts1[i]+normal1[i]*dy)
  pts.append(pts4[j])

  knots = CurveLengthParametrization(pts)
  curve = InterpolateCurve(pts,knots)
  return pts, curve

def _UniformParametrization(pts):
  n = len(pts)
  knots = range(0,n)

  return knots

def TFIOrthogonalSurface(curve ,rorto=0.9, no=0, r=0.0, nr=0, re=0.05, nre=0, ne1=0, ne2=0):
  pts1 = GetCurvePoints(curve[0])
  pts2 = GetCurvePoints(curve[1])
  pts3 = GetCurvePoints(curve[2])
  pts4 = GetCurvePoints(curve[3])

  n1 = len(pts1)
  n2 = len(pts3)

  pts = pts1

  factoro = 1.0
  factor = 1.0
  factore = 1.0
  crv = curve[0].Clone()
  for j in range(1,n2-1):
    ptsi, crvi = TFIcurve([crv,curve[1],curve[2],curve[3]],j,1.0)
    WriteG2('crv.g2', crvi)
    ptso, crvo = _OrthogonalCurve(crv,curve[1],curve[2],curve[3],j)

    if (j > no):
      factoro = factoro*rorto*rorto
      factor = factor*r
      factore = factore*re*re

    ptsc = []
    for i in range(0,n1):
      pt = factoro*ptso[i] + (1.0-factoro)*ptsi[i]
      ptsc.append(pt)

    knotsc = CurveLengthParametrization(ptsc)
    ds = knotsc[-1]/(len(knotsc)-1)
    crv = InterpolateCurve(ptsc,knotsc)

    ptsm = GetCurvePoints(crv)
    kn = crv.GetKnots()

    if (j > no):
      if ((nre > 0) or (nr>0)):
        if (nre > 0):
          crv1 = crv.Clone()
          crv = LaplaceSmoothingEnds(crv1,ne1,ne2,factore,nre)

        params = crv.GetKnots()
        for xi in params:
          pts.append(crv.Evaluate(xi))
    else:
      for pt in ptsc:
        pts.append(pt)

  for pt in pts2:
    pts.append(pt)

  xi = _UniformParametrization(pts1)
  eta = _UniformParametrization(pts3)

  surface = InterpolateSurface(pts,xi,eta)
  return surface

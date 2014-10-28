from GoTools import *
# from GoTools.CurveFactory import *
# from GoTools.SurfaceFactory import *
from GeoUtils.CurveUtils import *
import GeoUtils.Interpolate as ip

from itertools import product, chain
from numpy import cumsum
from operator import itemgetter


def LinearSurface(curves, eval_xi=None, eval_eta=None):
  """Performs linear transfinite interpolation between four curves to form a surface.
  @param curves: The curves. Curves 1 and 2 run from curve 3 to curve 4, and curves
  3 and 4 run from curve 1 to curve 2.
  @type curves: List of Curve
  @param eval_xi: Evaluation indices in the u-direction (optional, default all)
  @type eval_xi: List of integer, or None
  @param eval_eta: Evaluation indices in the v-direction (optional, default all)
  @type eval_eta: List of integer, or None
  @return: The generated surface (if eval_xi and eval_eta are not given) or the point cloud
  @rtype: Surface or List of Point
  """

  pts_ximin, pts_ximax, pts_etamin, pts_etamax = map(_GetExtPoints, curves)
  nxi  = len(pts_ximin)  - 2
  neta = len(pts_etamin) - 2
  ximin, ximax, etamin, etamax = map(
    lambda p: CurveLengthParametrization(p, True),
    [pts_ximin, pts_ximax, pts_etamin, pts_etamax])

  def ConvertEval(kts, length):
    ret = []
    for k in kts:
      if k == 0:
        ret.append(k)
      elif k == length - 1:
        ret.append(k+2)
      else:
        ret.append(k+1)
    return ret

  # Evaluate at all points if not specified
  interpolate = not eval_xi and not eval_eta
  if eval_xi is None:
    eval_xi = range(0, nxi+2)
  else:
    eval_xi = ConvertEval(eval_xi, nxi)
  if eval_eta is None:
    eval_eta = range(0, neta+2)
  else:
    eval_eta = ConvertEval(eval_eta, neta)

  points = []
  for eta, xi in product(eval_eta, eval_xi):
    emin, emax, pemin, pemax = map(itemgetter(eta), [etamin, etamax, pts_etamin, pts_etamax])
    xmin, xmax, pxmin, pxmax = map(itemgetter(xi),  [ximin,  ximax,  pts_ximin,  pts_ximax])

    # Solve a small linear system to determine the optimal placement of the node.
    # Make a drawing, it will become clear.
    # [    1.0     xmin-xmax ] [ x ]     [ xmin ]
    # [                      ] [   ]  =  [      ]
    # [ emin-emax     1.0    ] [ y ]     [ xmax ]
    det = 1.0 - (xmax - xmin) * (emax - emin)
    x = (xmin + emin * (xmax - xmin)) / det
    y = (emin + xmin * (emax - emin)) / det

    pt = (1-y)*pxmin + y*pxmax + (1-x)*pemin + x*pemax
    pt -= (1-x)*(1-y)*pts_ximin[0] + x*(1-y)*pts_ximin[-1]
    pt -= (1-x)*y*pts_ximax[0] + x*y*pts_ximax[-1]
    points.append(pt)

  if interpolate:
    return ip.InterpolateSurface(points,
                                 [0, 0.5] + range(1, nxi-1) + [nxi-1.5, nxi-1],
                                 [0, 0.5] + range(1, neta-1) + [neta-1.5, neta-1],
                                 [0]*3 + range(nxi) + [nxi-1]*3,
                                 [0]*3 + range(neta) + [neta-1]*3)
  return points


def OrthogonalSurface(curves, init_normalf, normalf=None, ranges=[],
                      fac_blend=0.9, fac_smooth=0.98, nsweeps=0, bnd_layer=0):
  """Performs orthogonalized transfinite interpolation between four curves to
  form a surface.
  @param curves: The curves, Curves 1 and 2 run from curve 3 to curve 4, and curves
  3 and 4 run from curve 1 to curve 2. Orthogonalization is performed near curve 1.
  @type curves: List of Curve
  @param init_normalf: Function to compute the intial normal vectors.
  @type init_normalf: Function
  @param normalf: (Optional) Function to compute curve normals during the intermediate
  stages Use this for performance if possible. Must be a function that takes two curves,
  crv_prev and crv (the last two iterations) and returns the normals of crv, pointing
  away from crv_prev.
  @type normalf: Function or None
  @param ranges: (Optional) A list of point ranges within which Laplacian smoothing
  will be performed. Each range is a tuple of the form (start, end) following
  regular Python range convention (i.e. start is included, end is not).
  @type ranges: List of (Integer, Integer)
  @param fac_blend: Determines how quickly the surface blends into linear TFI. If close
  to 1, it will happen slowly. If close to 0, it will happen (very) quickly.
  @type fac_blend: Float
  @param fac_smooth: Determines how quickly the strength of the Laplacian smoothing will
  increase. If close to 1, it will happen slowly. If close to 0, it will happen (very)
  quickly.
  @param nsweeps: (Optional, default 0) Number of Laplacian smoothing sweeps to perform.
  @type nsweeps: Integer
  @param bnd_layer: (Optional, default 0) Number of elements in the boundary layer.
  Within this layer, there will be no blending with linear TFI and no smoothing.
  @type bnd_layer: Integer
  @return: The generated surface
  @rtype: Surface
  """

  if normalf is None:
    # The default normalf function computes (t % v) % t, where
    # - t is the tangent vector at a point of the current curve
    # - v is the difference between two corresponding points of
    #   the current and previous curves
    # - % is the cross product
    def normalf(crv_prev, crv, kts):
      pre = [crv_prev.Evaluate(k) for k in kts]
      post = [crv.Evaluate(k) for k in kts]
      outward = [(pb - pa).Normalize() for pa, pb in zip(pre, post)]
      tangents = [crv.EvaluateTangent(k) for k in crv.GetKnots()]
      return [((tng % out) % tng).Normalize() for tng, out in zip(tangents, outward)]

  crv, endcurve, left, right = curves
  ekts = _GetExtKnots(crv)
  startpts, endpts = map(_GetExtPoints, [crv, endcurve])

  left = left.Clone()
  lkts = left.GetKnots()
  left.InsertKnot((lkts[0] + lkts[1])/2)
  left.InsertKnot((lkts[-2] + lkts[-1])/2)

  right = right.Clone()
  rkts = right.GetKnots()
  right.InsertKnot((rkts[0] + rkts[1])/2)
  right.InsertKnot((rkts[-2] + rkts[-1])/2)

  pts = [startpts]
  normals = init_normalf(crv, _GetExtKnots(crv))

  r_blend, r_smooth = 1.0, 1.0
  for j in xrange(1, len(left.GetKnots()) - 1):

    # Increase strength of blending and smoothing
    if j > bnd_layer:
      r_blend *= fac_blend**2
      r_smooth *= fac_smooth**2

    cl = left.GetSubCurve(left.GetKnots()[j-1], left.GetKnots()[-1])
    cr = right.GetSubCurve(right.GetKnots()[j-1], right.GetKnots()[-1])

    # Simple orthogonal projection
    ptso = _OrthogonalCurve(crv, cl.Evaluate(cr.GetKnots()[1]),
                            cr.Evaluate(cr.GetKnots()[1]), normals)

    # Perform smoothing if needed
    if j > bnd_layer and nsweeps > 0:
      ptso = _CurveSmoothing(ptso, ranges, nsweeps=nsweeps, r=1.0-r_smooth)

    # Linear TFI
    ptsi = LinearSurface([crv, endcurve, cl, cr], eval_eta=[1])

    # Blend the two together
    ptsc = [r_blend*pto + (1.0-r_blend)*pti for pto, pti in zip(ptso, ptsi)]

    # Create the new and the past curves
    crv_prev = crv

    crv = ip.InterpolateCurve(ptsc, ekts, [ekts[0]]*4 + ekts[2:-2] + [ekts[-1]]*4)

    # Update the normals and save the new points
    normals = normalf(crv_prev, crv, ekts)
    pts.append(ptsc)

  pts.append(endpts)

  nxi = len(ekts) - 2
  neta = len(lkts)
  return ip.InterpolateSurface(list(chain.from_iterable(pts)),
                               [0, 0.5] + range(1, nxi-1) + [nxi-1.5, nxi-1],
                               [0, 0.5] + range(1, neta-1) + [neta-1.5, neta-1],
                               [0]*3 + range(nxi) + [nxi-1]*3,
                               [0]*3 + range(neta) + [neta-1]*3)


# def LinearVolume(surfaces, eval_xi=None, eval_eta=None, eval_zeta=None, interpolate=True, order=4):
#   """Performs linear transfinite interpolation between six surfaces to form a volume.
#   @param surfaces: The surfaces. Surfaces 1 and 2, 3 and 4, as well as 5 and 6 are
#   opposite from each other and parametrized in the same direction. The first parameter
#   direction of surfaces 1 through 4 runs from surface 5 to surface 6. The first parameter
#   direction of surfaces 5 and 6 runs from surface 3 to surface 4.
#   @type curves: List of Surface
#   @param eval_xi: Evaluation indices in the u-direction (optional, default all)
#   @type eval_xi: List of integer, or None
#   @param eval_eta: Evaluation indices in the v-direction (optional, default all)
#   @type eval_eta: List of integer, or None
#   @param eval_zeta: Evaluation indices in the w-direction (optional, default all)
#   @type eval_zeta: List of integer, or None
#   @param interpolate: If true, returns an interpolated surface (optional, default true)
#   @type interpolate: Bool
#   @param order: Interpolation order, if interpolate is true (optional, default 4)
#   @type order: Integer
#   @return: The generated volume (if interpolate is true) or the point cloud
#   @rtype: Volume or List of Point
#   """


def _UniformParametrization(pts):
  return range(0, len(pts))


def _OrthogonalCurve(curve, pa, pb, normals):
  pre = _GetExtPoints(curve)
  knots = CurveLengthParametrization(pre, True)
  da, db = abs(pa - pre[0]), abs(pb - pre[-1])

  post = [p + n * ((1.0-k)*da + k*db) for p, n, k in zip(pre, normals, knots)[1:-1]]
  return [pa] + post + [pb]


def _CurveSmoothing(pts, ranges, nsweeps=1, r=0.2):
  orig_knots = CurveLengthParametrization(pts, False)
  curve = ip.InterpolateCurve(pts, orig_knots, [orig_knots[0]]*4 + orig_knots[2:-2] + [orig_knots[-1]]*4)
  knots = list(orig_knots)

  N = len(orig_knots)

  for _ in xrange(0, nsweeps):
    ks = list(knots)

    for rng in ranges:
      indices = range(rng[0], rng[1])[1:-1]
      for i in indices:
        ks[i%N] = knots[i%N] + r/2 * (knots[(i+1)%N] + knots[(i-1)%N] - 2*knots[i%N])

    knots = ks

  return [pts[0]] + [curve.Evaluate(k) for k in knots[1:-1]] + [pts[-1]]


def _GetSurfacePoints(srf):
  kus, kvs = srf.GetKnots()
  return [[srf.Evaluate(ku, kv) for kv in kvs] for ku in kus]


def _PtsToLists(pts):
  return [map(itemgetter(i), pts) for i in xrange(3)]


def _GetExtPoints(curve):
  return [curve.Evaluate(k) for k in _GetExtKnots(curve)]


def _GetExtKnots(thing):
  def fix(kts):
    return [kts[0], (kts[0]+kts[1])/2] + kts[1:-1] + [(kts[-2]+kts[-1])/2, kts[-1]]
  kts = thing.GetKnots()
  if isinstance(kts, tuple):
    return tuple([fix(ks) for ks in kts])
  return fix(kts)

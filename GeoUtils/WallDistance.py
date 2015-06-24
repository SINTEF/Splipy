from GoTools import *
import numpy as np
import math, sys
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize

# Parallel code copied from the nutils.org project.
import os, multiprocessing, time

def progress(title, i, tot):
    args = (title, i, tot, 100 * float(i) / tot)
    s = '\r%s: %i/%i (%.2f%%)' % args
    if i == tot:
        s += '\n'
    sys.stdout.write(s)
    sys.stdout.flush()

class Fork( object ):
  'nested fork context, unwinds at exit'

  def __init__( self, nprocs ):
    'constructor'

    self.nprocs = nprocs

  def __enter__( self ):
    'fork and return iproc'

    for self.iproc in range( self.nprocs-1 ):
      self.child_pid = os.fork()
      if self.child_pid:
        break
    else:
      self.child_pid = None
      self.iproc = self.nprocs-1
    return self.iproc

  def __exit__( self, *exc_info ):
    'kill all processes but first one'

    exctype, excvalue, tb = exc_info
    status = 0
    try:
      if exctype == KeyboardInterrupt:
        status = 1
      elif exctype == GeneratorExit:
        if self.iproc:
          log.stack( 'generator failed with unknown exception' )
        status = 1
      elif exctype:
        if self.iproc:
          log.stack( exc_info )
        status = 1
      if self.child_pid:
        child_pid, child_status = os.waitpid( self.child_pid, 0 )
        if child_pid != self.child_pid:
          log.error( 'pid failure! got %s, was waiting for %s' % (child_pid,self.child_pid) )
          status = 1
        elif child_status:
          status = 1
    except: # should not happen.. but just to be sure
      status = 1
    if self.iproc:
      os._exit( status )
    if not exctype:
      assert status == 0, 'one or more subprocesses failed'

def pariter( iterable, nprocs ):
  'iterate parallel, helper generator'
  shared_iter = multiprocessing.RawValue( 'i', nprocs )
  lock = multiprocessing.Lock()
  with Fork( nprocs ) as iproc:
    iiter = iproc
    for n, it in enumerate( iterable ):
      if n < iiter:
        continue
      assert n == iiter
      yield it
      with lock:
        iiter = shared_iter.value
        shared_iter.value = iiter + 1

def shzeros( shape, dtype=float ):
  'create zero-initialized array in shared memory'
  if isinstance( shape, int ):
    shape = shape,
  else:
    assert all( isinstance(sh,int) for sh in shape )
  size = np.product( shape ) if shape else 1
  if dtype == float:
    typecode = 'd'
  elif dtype == int:
    typecode = 'i'
  else:
    raise Exception, 'invalid dtype: %r' % dtype
  buf = multiprocessing.RawArray( typecode, size )
  return np.frombuffer( buf, dtype ).reshape( shape )

def distanceFunction(patches, wallset, process=[], display=True):
    """Calculate shortest wall distance
       @param surfaces: Patches in model
       @type surfaces: List of Surface or Volume
       @param wallset: The edges/faces defining the wall
       @type wallset: List of tuple of (patch, edge/face) numbers
       @param process: Optional list of patches to process
       @type process: List of integer
       @param display: (Optional) Write progress to stdout
       @type display: Boolean
       @return: Coefficients of distance field
       @rtype: List of list of float
    """

    assert len(wallset), 'Empty wallset'
    if isinstance(patches[0], Surface):
      return _distanceFunction2D(patches, wallset, process, display=display)

    if isinstance(patches[0], Volume):
      return _distanceFunction3D(patches, wallset, process, display=display)

    raise ValueError('Can only process surfaces or volumes')

def _distanceFunction2D(surfaces, wallset, patches=[], display=True):
    """Calculate shortest wall distance in 2D
       @param surfaces: Surfaces in model
       @type surfaces: List of Surface
       @param wallset: The edges defining the wall
       @type wallset: List of tuple of (patch, edge) numbers
       @param: Optional list of patches to process
       @type patches: List of integer
       @param display: (Optional) Write progress to stdout
       @type display: Boolean
       @return: Coefficients of distance field
       @rtype: List of list of float
    """

    worksurfaces = []
    if len(patches):
      for patch in patches:
          worksurfaces.append(surfaces[patch-1])
    else:
      worksurfaces = surfaces

    walledges = []
    for idx, info in wallset.iteritems():
      for edge in info.edge:
        walledges.append( surfaces[idx-1].GetEdges()[edge-1] )

    D = _calcDistScipy(walledges, worksurfaces, display=display)
    return D

def _calcDistScipy(wallcurves, worksurfaces, display=True):
    """Calculate minimum wall distance using scipy and minimize scalar
       @param wallcurves: List of curves describing the wall
       @type wallcurves: List of curve
       @param worksurfaces: Surfaces to process
       @type worksurfaces: List of Surface
       @param display: (Optional) Write progress to stdout
       @type display: Boolean
       @return Wall distance
       @rtype List of doubles
    """

    patch_size = lambda patch: np.prod( map( len, patch.GetGreville() ) )
    max_patch_size = max(map(patch_size, worksurfaces))
    D = shzeros((len(worksurfaces), max_patch_size))
    lenD = shzeros( len(worksurfaces) )
    t0 = time.time()
    for srfID, surface in pariter(enumerate(worksurfaces), GetProcessorCount()):
        if display:
            progress('Computing wall distance', srfID+1, len(D))

        (knots_xi, knots_eta) = surface.GetGreville()
        wdist = np.zeros((len(knots_xi), len(knots_eta)))

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
                    lb = crv_knots_xi[0]
                    ub = crv_knots_xi[-1]
                    res = minimize_scalar(_calcPtsDistanceCurve, s0, args=(curve, pt), bounds=(lb,ub), method='bounded', options={'xtol': 1e-10, 'xatol': 1e-10, 'disp': False})
                    tmp = _calcPtsDistance(curve.Evaluate(res.x), pt)
                    if tmp < mindist:
                        mindist = tmp
                wdist[i,j] = mindist
                j = j+1
            i = i+1
        size = wdist.size
        data = surface.Interpolate(list(wdist.flatten('F')))
        D[srfID,:size] = data
        lenD[srfID] = size

    return [list(Di[:lenDi]) for (Di,lenDi) in zip(D,lenD)]


def _getWallCurve(surface, idx):
    edges = surface.GetEdges()
    return edges[idx]

def _getWallFace(volume, idx):
    faces = volume.GetFaces()
    return faces[idx]

def _calcPtsDistance(pt1, pt2):
    """
    Calculate shortest distance between two points
    """
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2)

def _calcPtsDistanceCurve(s, curve, pt2):
    """
    Calculate shortest distance between two points
    """
    pt1 = curve.Evaluate(s)
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2)

def _calcPtsDistance3D(pt1, pt2):
    """
    Calculate shortest distance between two points
    """
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2)

def _calcPtsDistanceSurface(s, surface, pt2):
    """
    Calculate shortest distance between two points
    """
    pt1 = surface.Evaluate(s[0], s[1])
    return np.sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2)

def _distanceFunction3D(volumes, wallset, patches=[], display=True):
    """Calculate shortest wall distance in 3D
       @param volumes: Volumes in model
       @type volumes: List of Volume
       @param wallset: The faces defining the wall
       @type wallset: List of tuple of (patch, face) numbers
       @param: Optional list of patches to process
       @type patches: List of integer
       @param display: (Optional) Write progress to stdout
       @type display: Boolean
       @return: Coefficients of distance field
       @rtype: List of list of float
    """

    workvolumes = []
    if len(patches):
      for patch in patches:
        workvolumes.append(volumes[patch-1])
    else:
      workvolumes = volumes

    wallfaces = []
    for idx, info in wallset.iteritems():
      for face in info.face:
        wallfaces.append( volumes[idx-1].GetFaces()[face-1] )

    D = _calcDistScipy3D(wallfaces, workvolumes, display=display)

    return D


def _calcDistScipy3D(wallfaces, workvolumes, display=True):
    """Calculate minimum wall distance using scipy and minimize scalar
       @param wallcurves: List of faces describing the wall
       @type wallcurves: List of Surface
       @param workvolumes: Volumes to process
       @type workvolumes: List of Volume
       @param display: (Optional) Write progress to stdout
       @type display: Boolean
       @return Wall distance
       @rtype List of doubles
    """

    patch_size = lambda patch: np.prod( map( len, patch.GetGreville() ) )
    max_patch_size = max( map( patch_size, workvolumes ) )
    D = shzeros( (len(workvolumes), max_patch_size) )
    lenD = shzeros( len(workvolumes) )
    t0 = time.time()
    for volID, volume in pariter( enumerate(workvolumes),GetProcessorCount()):
        if display:
            progress('Computing wall distance', volID+1, len(D))
        (knots_xi, knots_eta, knots_gamma) = volume.GetGreville()
        wdist = np.zeros((len(knots_xi), len(knots_eta), len(knots_gamma)))

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
                        b1 = (crv_knots_xi[0], crv_knots_xi[-1])
                        b2 = (crv_knots_eta[0], crv_knots_eta[-1])
                        res = minimize(_calcPtsDistanceSurface, s0, args=(face, pt), method='L-BFGS-B', bounds=(b1,b2), jac=False)
                        tmp = _calcPtsDistance3D(face.Evaluate(res.x[0], res.x[1]), pt)
                        if not np.isfinite(tmp): print 'min dist @ vol_%i[%i,%i,%i] = '%(volID,i,j,k), tmp
                        if tmp < mindist:
                            mindist = tmp
                    wdist[i,j,k] = mindist
                    k = k+1
                j = j+1
            i = i+1
        size = wdist.size
        data = volume.Interpolate(list(wdist.flatten('F')))
        D[volID,:size] = data
        lenD[volID] = size
    return [list(Di[:lenDi]) for (Di,lenDi) in zip(D,lenD)]

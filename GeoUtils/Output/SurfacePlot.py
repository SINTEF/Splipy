__doc__ = 'Class for plots from HDF5 data'

from GeoUtils.IO import InputFile, IFEMResultDatabase
from GoTools import *
from itertools import product
from matplotlib import pyplot
from numpy import array, inf, log10, sqrt
from operator import attrgetter
from xml.dom import minidom


def Plot2DPatches(patches, tess=3, bbox=None, *args, **kwargs):
  """ Plots a collection of 2D patches (in the xy-plane) using matplotlib.
      Remaining arguments are passed to the matplotlib plot() function.
      @param patches: The patches to plot.
      @type patches: Iterable collection of Surface.
      @param tess: Number of tesselation points to use per knotspan.
      @type tess: Integer > 0.
      @param bbox: Optional bounding box (xmin, xmax, ymin and ymax)
      @type bbox: List of float
  """
  def plot(xs, ys):
    for i, (x, y) in enumerate(zip(xs, ys)):
      lw = 0.2
      if i == 0 or i == len(xs) - 1:
        lw = 1.0
      pyplot.plot(x, y, 'k-', linewidth=lw, *args, **kwargs)

  for patch in patches:
    ku, kv = patch.GetTesselationParams(tess, 1)
    points = patch.EvaluateGrid(ku, kv)
    points = [points[i*len(ku) : (i+1)*len(ku)] for i in xrange(len(kv))]

    xs = [[p[0] for p in ps] for ps in points]
    ys = [[p[1] for p in ps] for ps in points]
    plot(xs, ys)

    ku, kv = patch.GetTesselationParams(1, tess)
    points = patch.EvaluateGrid(ku, kv)
    points = [points[i*len(ku) : (i+1)*len(ku)] for i in xrange(len(kv))]

    xs = zip(*[[p[0] for p in ps] for ps in points])
    ys = zip(*[[p[1] for p in ps] for ps in points])
    plot(xs, ys)

  if not bbox:
    model = SurfaceModel(patches)
    xmin, xmax, ymin, ymax, _, _ = model.GetBoundingBox()
  else:
    xmin, xmax, ymin, ymax = bbox
  dx = (xmax - xmin) * 0.02
  dy = (ymax - ymin) * 0.02

  pyplot.xlim([xmin - dx, xmax + dx])
  pyplot.ylim([ymin - dy, ymax + dy])

  pyplot.axes().set_aspect('equal', 'datalim')
  pyplot.show()


class SurfacePlot:
  """Plot (flow) data from one topologyset in an HDF5 file."""

  # Basic operations
  def __init__( self, hdf5, topo, surf, levels=-1, setup=None, opts={} ):
    """Constructor
       @param hdf5: HDF5-file with data
       @type hdf5: str
       @param topo: XINP-file with <topologysets> tag
       @type topo: str
       @param surf: surface to extract data at
       @type surf: list of tuple of (patch, face) numbers
       @param levels: time levels in hdf5, default to last
       @type levels: int, list of ints, None
       @param dim: 2D or 3D
       @type dim: int
       @param setup: XINP-file with <fluidproperties> tag
       @type setup: str
       @param opts: options
       @type opts: dict
    """
    # fix filenames if necessary
    hdf5 = hdf5[:-5] if hdf5.endswith('hdf5') else hdf5
    topo = topo if topo.endswith('.xinp') else topo+'.xinp'
    setup = topo if not setup else \
            setup if setup.endswith('.xinp') else \
            setup+'.xinp'

    # find fieldname
    dom = minidom.parse( hdf5+'.xml' )
    for elem in dom.firstChild.getElementsByTagName( 'entry' ):
      if elem.getAttribute('description')=='velocity' and not \
         elem.getAttribute('name')=='restart': break
    tag = elem.getAttribute('name')
    ischorin = not tag.endswith('+p')

    # prepare data
    output = IFEMResultDatabase( hdf5 )
    basis = output.GetBasisForField( tag )
    # if levels is iterable/None, average over specified/all time levels
    levels = levels if hasattr( levels, '__iter__' ) else \
            range(output.GetTimeLevels()+1) if levels is None else \
            [output.GetTimeLevels() if levels<0 else levels]
    self.setup = minidom.parse( setup )
    geom = output.GetGeometry( basis, 0 ) # domain geometry, assert level=0
    self.dim = len(geom[0].GetOrder()) # dimension of domain

    # average if multiple time levels (TODO: inefficiently for whole volume)
    def _deinterleave( level ):
      if ischorin:
        return output.GetFieldCoefs(tag,level), output.GetFieldCoefs('pressure',level)
      else:
        vdofs, pdofs = [], []
        for coef in output.GetFieldCoefs(tag,level):
          vdofs.append( array( coef ).reshape( -1, self.dim+1 )[:,:-1].flatten().tolist() )
          pdofs.append( array( coef ).reshape( -1, self.dim+1 )[:,-1:].flatten().tolist() )
        return vdofs, pdofs
    vdofs, pdofs = [0]*len(geom), [0]*len(geom)
    for level in levels:
      vdofs_l, pdofs_l = _deinterleave( level )
      for i, (vdofs_li, pdofs_li) in enumerate(zip(vdofs_l, pdofs_l)):
        vdofs[i] += array( vdofs_li )
        pdofs[i] += array( pdofs_li )
    _average = lambda arr: ( arr/float(len(levels)) ).tolist()
    velo, pres = [], []
    for i, (geom_i, vdofs_i, pdofs_i) in enumerate(zip(geom, vdofs, pdofs)):
      velo.append( geom_i.Clone(_average(vdofs_i)) )
      pres.append( geom_i.Clone(_average(pdofs_i)) )

    # translate to boundary data
    self.surf = InputFile(topo).GetTopologySet( surf, convention="gotools" )
    name = 'Face' if self.dim==3 else 'Edge'
    self.geom, self.velo, self.pres = [], [], []
    for i, patchinfo in self.surf.iteritems():
      for face in attrgetter(name.lower())( patchinfo ):
        faces = lambda field: attrgetter('Get%ss'%name)( field )()
        if self.dim==3 and opts.has_key('edge#'):
          for j in opts['edge#']:
            self.geom.append( faces(geom[i-1])[face-1].GetEdges()[j] )
            self.velo.append( faces(velo[i-1])[face-1].GetEdges()[j] )
            self.pres.append( faces(pres[i-1])[face-1].GetEdges()[j] )
        else:
          self.geom.append( faces(geom[i-1])[face-1] )
          self.velo.append( faces(velo[i-1])[face-1] )
          self.pres.append( faces(pres[i-1])[face-1] )

  def __iter__( self ):
    "Iterate over control points of 2D/3D surface"
    # Order reversal [::-1] in _astuple makes __iter__ run fastest along the
    # first direction as is the convention in geomodeler.py (but not itertools)
    _astuple = lambda arg: arg[::-1] if isinstance(arg,tuple) else (arg,)
    funcs = self.geom, self.velo, self.pres
    for i, patch in enumerate(self.geom):
      for knot in product( *_astuple(patch.GetGreville()) ):
        args = [func[i].Evaluate( *_astuple(knot) ) for func in funcs] \
             + [func[i].EvaluateTangent( *_astuple(knot) ) for func in funcs]
        yield i, args

  def min( self, f ):
    """The minimum value of f (considering only Greville points)"""
    min_f = inf
    for i, arg in self:
      min_f = min( f(*arg), min_f )
    return min_f

  def max( self, f ):
    """The maximum value of f (considering only Greville points)"""
    max_f = -inf
    for i, arg in self:
      max_f = max( f(*arg), max_f )
    return max_f

  def interpolate( self, f ):
    """Interpolate a function onto the basis of self.geom, not necessary
       in general.
       @param f: function of primary variables to be interpolated
       @type f: function
       @return: interpolated function, a clone of self.geom
       @rtype: list of Curve/Surface
    """
    i0, data = None, len(self.geom)*[[]]
    for i, args in self:
      if not i0==i: # for every new patch i
        if isinstance( i0, int ): # skip i0==None
          data[i0] = self.geom[i0].Clone( self.geom[i0].Interpolate(datai) )
        i0, datai = i, []
      datai.append( evaluable(*args) )
    data[i] = self.geom[i].Clone( self.geom[i].Interpolate(datai) )
    return data

  # Built-in plottable fields
  x = lambda self: lambda g, v, p, dg, dv, dp: g[0]
  y = lambda self: lambda g, v, p, dg, dv, dp: g[1]
  z = lambda self: lambda g, v, p, dg, dv, dp: g[2]

  def normal( self ):
    """Normal vector assuming 2D or extrusion in 3rd direction, as well as
       left-handed patches."""
    return (lambda g, v, p, dg, dv, dp:  dg[1]/sqrt(dg[0]**2+g[1]**2),
            lambda g, v, p, dg, dv, dp: -dg[0]/sqrt(dg[0]**2+g[1]**2))

  def pressureCoefficient( self, pinf=0, Uinf=1 ):
    params = self.setup.getElementsByTagName('fluidproperties')[0]
    rho = float( params.getAttribute('rho') )
    return lambda g, v, p, dg, dv, dp: (p[0]-pinf)/(.5*rho*Uinf**2)

  def logCp( self, pinf=0, Uinf=1 ):
    params = self.setup.getElementsByTagName('fluidproperties')[0]
    rho = float( params.getAttribute('rho') )
    denom = .5*rho*Uinf**2
    def logCp( g, v, p, dg, dv, dp ):
      num = p[0]-pinf
      return log10(1+num/denom) if num>0 else -log10(1-num/denom)
    return logCp

  # Built-in plot commands
  def plot( self, x, y, fmt='k.-', axis=None, **kwargs ):
    """Built-in plotting function similar but not identical in behavior to
       matplotlib.pyplot.plot, namely, requires both x and y arguments, and
       accepts an axis to plot inside. The x, y arguments here are moreover
       simply function objects that take the primary variables as arguments,
       see the built-in examples.
       @param data: x-coordinate
       @type data: function
       @param data: y-coordinate
       @type data: function
       @param fmt: line format string, optional (as pyplot.plot)
       @type fmt: str
       @param axis: axis to plot inside, optional (default creates new axis)
       @type axis: matplotlib.axes.AxesSubplot
       @param kwargs: any other options for matplotlib.pyplot.plot
       @type kwargs: dict
    """
    # Note, in function call axis cannot default to pyplot.figure().gca() as
    # this is created upon loading the class, not upon calling the function.
    # This way a new axis is made for each function call (without axis arg).
    axis = pyplot.figure().gca() if axis is None else axis
    label = kwargs.pop( 'label', '_nolabel_' )
    funcs = self.geom, self.velo, self.pres
    i0 = None
    for i, args in self:
      if not i0==i: # for every new patch i
        if isinstance( i0, int ): # skip i0==None
          axis.plot( xi, yi, fmt, label='_nolabel_', **kwargs )
        i0, xi, yi = i, [], []
      xi.append( x(*args) )
      yi.append( y(*args) )
    axis.plot( xi, yi, fmt, label=label, **kwargs )

  def quiver( self, x, y, u, v, axis=None, **kwargs ):
    """Built-in plotting function similar but not identical in behavior to
       matplotlib.pyplot.quiver, namely, requires x, y, u and v arguments,
       and accepts an axis to plot inside. The x, y, u, v arguments here
       are moreover simply function objects that take the primary variables
       as arguments, see the built-in examples.
       @param data: x-coordinate
       @type data: function
       @param data: y-coordinate
       @type data: function
       @param data: u-component of arrow
       @type data: function
       @param data: y-component of arrow
       @type data: function
       @param axis: axis to plot inside, optional (default creates new axis)
       @type axis: matplotlib.axes.AxesSubplot
       @param kwargs: any other options for matplotlib.pyplot.plot
       @type kwargs: dict
    """
    axis = pyplot.figure().gca() if axis is None else axis
    label = kwargs.pop( 'label', '_nolabel_' )
    funcs = self.geom, self.velo, self.pres
    i0 = None
    for i, args in self:
      if not i0==i: # for every new patch i
        if isinstance( i0, int ): # skip i0==None
          axis.quiver( xi, yi, ui, vi, label='_nolabel_', **kwargs )
        i0, xi, yi, ui, vi = i, [], [], [], []
      xi.append( x(*args) )
      yi.append( y(*args) )
      ui.append( u(*args) )
      vi.append( v(*args) )
    axis.quiver( xi, yi, ui, vi, label=label, **kwargs )

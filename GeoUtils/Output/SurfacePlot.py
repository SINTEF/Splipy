__doc__ = 'Class for plots from HDF5 data'

import xml.dom.minidom
from GoTools import *

from itertools import product
from operator import attrgetter
from InputFile import *
from matplotlib import pyplot

class SurfacePlot:
  """Plot (flow) data from one topologyset in an HDF5 file."""

  def __init__( self, hdf5, topo, surf, level=-1, setup=None, opts={} ):
    """Constructor
       @param hdf5: HDF5-file with data
       @type hdf5: str
       @param topo: XINP-file with <topologysets> tag
       @type topo: str
       @param surf: surface to extract data at
       @type surf: list of tuple of (patch, face) numbers
       @param level: time level in hdf5, default to last
       @type level: int
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

    # prepare data
    output = IFEMResultDatabase( hdf5 )
    basis = output.GetBasisForField( 'pressure' )
    level = output.GetTimeLevels() if level<0 else level
    self.setup = xml.dom.minidom.parse( setup )

    # translate to boundary data
    geom = output.GetGeometry( basis, 0 ) # domain geometry, assert level=0
    self.dim = len(geom[0].GetOrder()) # dimension of domain
    self.surf = InputFile(topo).GetTopologySet( surf, convention="gotools" )
    field = 'u_x+u_y+u_z' if self.dim==3 else 'u_x+u_y' # assert Chorin
    velo = output.GetField( field, level, geom )
    pres = output.GetField( 'pressure', level, geom )
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
    for i, patch in enumerate(self.geom):
      for knot in product( *_astuple(patch.GetGreville()) ):
        yield i, _astuple(knot) # reverse order back

  def getField( self, func ):
    """Get generic evaluable field for use in self.plot(),
       for examples of use, see the built-ins below.
       @param: a function of geometry, velocity and pressure
       @type: function
    """
    i0, data = None, len(self.geom)*[[]]
    for i, coord in self:
      if not i0==i: # for every new patch i
        if isinstance( i0, int ): # skip i0==None
          data[i0] = self.geom[i0].Clone( self.geom[i0].Interpolate(datai) )
        i0, coordi, datai = i, [], []
      g = self.geom[i].Evaluate(*coord)
      v = self.velo[i].Evaluate(*coord)
      p = self.pres[i].Evaluate(*coord)[0]
      coordi.append( map(float,coord) )
      datai.append( func(g,v,p) )
    data[i] = self.geom[i].Clone( self.geom[i].Interpolate(datai) )
    return data

  # Built-in plottable fields
  x = lambda self: self.getField( lambda g, v, p: g[0] )
  y = lambda self: self.getField( lambda g, v, p: g[1] )
  z = lambda self: self.getField( lambda g, v, p: g[2] )
  def pressureCoefficient( self, pinf=0, Uinf=1 ):
    params = self.setup.getElementsByTagName('fluidproperties')[0]
    rho = float( params.getAttribute('rho') )
    Cp = lambda g, v, p: (p-pinf)/(.5*rho*Uinf**2)
    return self.getField( Cp )

  def plot( self, x, y, fmt='k.-', axis=None, **kwargs ):
    """Built-in plotting function for SurfacePlot.getField output similar
       but not identical in behavior to matplotlib.pyplot.plot, namely,
       requires both x and y data, and accepts an axis to plot inside.
       @param data: x-coordinates, output of self.getField()
       @type data: list of Surface/Curve
       @param data: y-coordinates, output of self.getField()
       @type data: list of Surface/Curve
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
    i0 = None
    for i, coord in self:
      if not i0==i: # for every new patch i
        if isinstance( i0, int ): # skip i0==None
          axis.plot( xi, yi, fmt, label='_nolabel_', **kwargs )
        i0, xi, yi = i, [], []
        evalf = lambda func: func[i].Evaluate( *coord ) # can go up-scope as coord a pointer
      xi.append( evalf(x) )
      yi.append( evalf(y) )
    axis.plot( xi, yi, fmt, label=label, **kwargs )

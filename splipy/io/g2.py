import numpy as np
from itertools import chain, product
from splipy import BSplineBasis, Curve, Surface, Volume, SplineObject
from .master import MasterIO
import splipy.surface_factory as SurfaceFactory
import splipy.curve_factory   as CurveFactory
import splipy.state as state
from numpy import sqrt, pi


class G2(MasterIO):

    def circle(self):
        dim   = int(     next(self.fstream).strip())
        r     = float(   next(self.fstream).strip())
        center= np.array(next(self.fstream).split(' '), dtype=float)
        normal= np.array(next(self.fstream).split(' '), dtype=float)
        xaxis = np.array(next(self.fstream).split(' '), dtype=float)
        param = np.array(next(self.fstream).split(' '), dtype=float)
        reverse =        next(self.fstream).strip() != '0'

        result = CurveFactory.circle(r=r, center=center, normal=normal, xaxis=xaxis)
        result.reparam(param)
        if reverse:
            result.reverse()
        return result

    def ellipse(self):
        dim   = int(     next(self.fstream).strip())
        r1    = float(   next(self.fstream).strip())
        r2    = float(   next(self.fstream).strip())
        center= np.array(next(self.fstream).split(' '), dtype=float)
        normal= np.array(next(self.fstream).split(' '), dtype=float)
        xaxis = np.array(next(self.fstream).split(' '), dtype=float)
        param = np.array(next(self.fstream).split(' '), dtype=float)
        reverse =        next(self.fstream).strip() != '0'

        result = CurveFactory.ellipse(r1=r1, r2=r2, center=center, normal=normal, xaxis=xaxis)
        result.reparam(param)
        if reverse:
            result.reverse()
        return result

    def line(self):
        dim      = int(     next(self.fstream).strip())
        start    = np.array(next(self.fstream).split(' '), dtype=float)
        direction= np.array(next(self.fstream).split(' '), dtype=float)
        finite   =          next(self.fstream).strip() != '0'
        param    = np.array(next(self.fstream).split(' '), dtype=float)
        reverse  =          next(self.fstream).strip() != '0'
        d = np.array(direction)
        s = np.array(start)
        d /= np.linalg.norm(d)
        if not finite:
            param = [-state.unlimited, +state.unlimited]

        result = CurveFactory.line(s-d*param[0], s+d*param[1])
        if reverse:
            result.reverse()
        return result


#   def cone(self):
#       dim      = int(     next(self.fstream).strip())
#       r        = float(   next(self.fstream).strip())
#       center   = np.array(next(self.fstream).split(' '), dtype=float)
#       z_axis   = np.array(next(self.fstream).split(' '), dtype=float)
#       x_axis   = np.array(next(self.fstream).split(' '), dtype=float)
#       angle    = float(   next(self.fstream).strip())
#       finite   =          next(self.fstream).strip() != '0'
#       param_u  = np.array(next(self.fstream).split(' '), dtype=float)
#       if finite:
#           param_v=np.array(next(self.fstream).split(' '), dtype=float)


    def cylinder(self):
        dim      = int(     next(self.fstream).strip())
        r        = float(   next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(' '), dtype=float)
        z_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        x_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        finite   =          next(self.fstream).strip() != '0'
        param_u  = np.array(next(self.fstream).split(' '), dtype=float)
        if finite:
            param_v=np.array(next(self.fstream).split(' '), dtype=float)
        else:
            param_v=[-state.unlimited, 0]
        swap     =          next(self.fstream).strip() != '0'

        center = center-z_axis*param_v[0]
        h      = param_v[1] - param_v[0]
        result = SurfaceFactory.cylinder(r=r, center=center, xaxis=x_axis, axis=z_axis, h=h)
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def disc(self):
        dim      = int(     next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(' '), dtype=float)
        r        = float(   next(self.fstream).strip())
        z_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        x_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        degen    =          next(self.fstream).strip() != '0'
        angles   = [float(  next(self.fstream).strip()) for i in range(4)]
        param_u  = np.array(next(self.fstream).split(' '), dtype=float)
        param_v  = np.array(next(self.fstream).split(' '), dtype=float)
        swap     =          next(self.fstream).strip() != '0'

        if degen:
            result = SurfaceFactory.disc(r=r, center=center, xaxis=x_axis, normal=z_axis, type='radial')
        else:
            if not(np.allclose(np.diff(angles), pi/2, atol=1e-10)):
                raise RuntimeError('Unknown square parametrization of disc elementary surface')
            result = SurfaceFactory.disc(r=r, center=center, xaxis=x_axis, normal=z_axis, type='square')
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def torus(self):
        dim      = int(     next(self.fstream).strip())
        r2       = float(   next(self.fstream).strip())
        r1       = float(   next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(' '), dtype=float)
        z_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        x_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        select_out=         next(self.fstream).strip() != '0' # I have no idea what this does :(
        param_u  = np.array(next(self.fstream).split(' '), dtype=float)
        param_v  = np.array(next(self.fstream).split(' '), dtype=float)
        swap     =          next(self.fstream).strip() != '0'

        result = SurfaceFactory.torus(minor_r=r1, major_r=r2, center=center, normal=z_axis, xaxis=x_axis)
        result.reparam(param_u, param_v)
        if(swap):
            result.swap()
        return result

    def sphere(self):
        dim      = int(     next(self.fstream).strip())
        r        = float(   next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(' '), dtype=float)
        z_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        x_axis   = np.array(next(self.fstream).split(' '), dtype=float)
        param_u  = np.array(next(self.fstream).split(' '), dtype=float)
        param_v  = np.array(next(self.fstream).split(' '), dtype=float)
        swap     =          next(self.fstream).strip() != '0'

        result = SurfaceFactory.sphere(r=r, center=center, xaxis=x_axis, zaxis=z_axis).swap()
        if swap:
            result.swap()
        result.reparam(param_u, param_v)
        return result

    g2_type = [100, 200, 700] # curve, surface, volume identifiers
    classes = [Curve, Surface, Volume]

    g2_generators = {120:line, 130:circle, 140:ellipse,
                     260:cylinder, 292:disc, 270:sphere, 290:torus} #, 280:cone

    def __init__(self, filename):
        if filename[-3:] != '.g2':
            filename += '.g2'
        self.filename = filename

    def __enter__(self):
        return self

    def write(self, obj):
        if not hasattr(self, 'fstream'):
            self.onlywrite = True
            self.fstream   = open(self.filename, 'w')
        if not self.onlywrite:
            raise IOError('Could not write to file %s' % (self.filename))

        """Write the object in GoTools format. """
        if isinstance(obj[0], SplineObject): # input SplineModel or list
            for o in obj:
                self.write(o)
            return

        for i in range(obj.pardim):
            if obj.periodic(i):
                obj = obj.split(obj.start(i), i)

        self.fstream.write('{} 1 0 0\n'.format(G2.g2_type[obj.pardim-1]))
        self.fstream.write('{} {}\n'.format(obj.dimension, int(obj.rational)))
        for b in obj.bases:
            self.fstream.write('%i %i\n' % (len(b.knots) - b.order, b.order))
            self.fstream.write(' '.join('%f' % k for k in b.knots))
            self.fstream.write('\n')

        for cp in obj:
            self.fstream.write(' '.join('%f' % x for x in cp))
            self.fstream.write('\n')

    def read(self):
        if not hasattr(self, 'fstream'):
            self.onlywrite = False
            self.fstream   = open(self.filename, 'r')

        if self.onlywrite:
            raise IOError('Could not read from file %s' % (self.filename))

        result = []

        for line in self.fstream:
            line = line.strip()
            if not line:
                continue

            # read object type
            objtype, major, minor, patch = map(int, line.split(' '))
            if (major, minor, patch) != (1, 0, 0):
                raise IOError('Unknown G2 format')

            # if obj type is in factory methods (cicle, torus etc), create it now
            if objtype in G2.g2_generators:
                print(objtype)
                constructor = getattr(self, G2.g2_generators[objtype].__name__)
                result.append( constructor() )
                continue

            # for "normal" splines (Curves, Surfaces, Volumes) create it now
            pardim = [i for i in range(len(G2.g2_type)) if G2.g2_type[i] == objtype]
            if not pardim:
                raise IOError('Unknown G2 object type {}'.format(objtype))
            pardim = pardim[0] + 1
            cls    = G2.classes[pardim-1]

            _, rational = next(self.fstream).strip().split(' ')
            rational = bool(int(rational))

            bases = [self.read_basis() for _ in range(pardim)]
            ncps = 1
            for b in bases:
                ncps *= b.num_functions()

            cps = [tuple(map(float, next(self.fstream).split(' ')))
                   for _ in range(ncps)]

            args = bases + [cps, rational]
            result.append(cls(*args))

        return result

    def read_basis(self):
        ncps, order = map(int, next(self.fstream).split(' '))
        kts = list(map(float, next(self.fstream).split(' ')))
        return BSplineBasis(order, kts, -1)

    def __exit__(self, exc_type, exc_value, traceback):
        pass


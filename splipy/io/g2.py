from itertools import chain, product

import numpy as np
from numpy import sqrt, pi, savetxt

from ..curve import Curve
from ..surface import Surface
from ..volume import Volume
from ..splineobject import SplineObject
from ..basis import BSplineBasis
from ..trimmedsurface import TrimmedSurface
from ..utils import flip_and_move_plane_geometry, rotate_local_x_axis
from .. import surface_factory, curve_factory, state

from .master import MasterIO


class G2(MasterIO):

    def read_next_non_whitespace(self):
        line = next(self.fstream).strip()
        while not line:
            line = next(self.fstream).strip()
        return line

    def circle(self):
        dim   = int(     self.read_next_non_whitespace().strip())
        r     = float(   next(self.fstream).strip())
        center= np.array(next(self.fstream).split(), dtype=float)
        normal= np.array(next(self.fstream).split(), dtype=float)
        xaxis = np.array(next(self.fstream).split(), dtype=float)
        param = np.array(next(self.fstream).split(), dtype=float)
        reverse =        next(self.fstream).strip() != '0'

        result = curve_factory.circle(r=r, center=center, normal=normal, xaxis=xaxis)
        result.reparam(param)
        if reverse:
            result.reverse()
        return result

    def ellipse(self):
        dim   = int(     self.read_next_non_whitespace().strip())
        r1    = float(   next(self.fstream).strip())
        r2    = float(   next(self.fstream).strip())
        center= np.array(next(self.fstream).split(), dtype=float)
        normal= np.array(next(self.fstream).split(), dtype=float)
        xaxis = np.array(next(self.fstream).split(), dtype=float)
        param = np.array(next(self.fstream).split(), dtype=float)
        reverse =        next(self.fstream).strip() != '0'

        result = curve_factory.ellipse(r1=r1, r2=r2, center=center, normal=normal, xaxis=xaxis)
        result.reparam(param)
        if reverse:
            result.reverse()
        return result

    def line(self):
        dim      = int(     self.read_next_non_whitespace().strip())
        start    = np.array(next(self.fstream).split(), dtype=float)
        direction= np.array(next(self.fstream).split(), dtype=float)
        finite   =          next(self.fstream).strip() != '0'
        param    = np.array(next(self.fstream).split(), dtype=float)
        reverse  =          next(self.fstream).strip() != '0'
        d = np.array(direction)
        s = np.array(start)
        # d /= np.linalg.norm(d)
        if not finite:
            param = [-state.unlimited, +state.unlimited]

        result = curve_factory.line(s+d*param[0], s+d*param[1])
        if reverse:
            result.reverse()
        return result


#   def cone(self):
#       dim      = int(     self.read_next_non_whitespace().strip())
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
        dim      = int(     self.read_next_non_whitespace().strip())
        r        = float(   next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(), dtype=float)
        z_axis   = np.array(next(self.fstream).split(), dtype=float)
        x_axis   = np.array(next(self.fstream).split(), dtype=float)
        finite   =          next(self.fstream).strip() != '0'
        param_u  = np.array(next(self.fstream).split(), dtype=float)
        if finite:
            param_v=np.array(next(self.fstream).split(), dtype=float)
        else:
            param_v=[-state.unlimited, state.unlimited]
        swap     =          next(self.fstream).strip() != '0'

        center = center + z_axis*param_v[0]
        h      = param_v[1] - param_v[0]
        result = surface_factory.cylinder(r=r, center=center, xaxis=x_axis, axis=z_axis, h=h)
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def disc(self):
        dim      = int(     self.read_next_non_whitespace().strip())
        center   = np.array(next(self.fstream).split(), dtype=float)
        r        = float(   next(self.fstream).strip())
        z_axis   = np.array(next(self.fstream).split(), dtype=float)
        x_axis   = np.array(next(self.fstream).split(), dtype=float)
        degen    =          next(self.fstream).strip() != '0'
        angles   = [float(  next(self.fstream).strip()) for i in range(4)]
        param_u  = np.array(next(self.fstream).split(), dtype=float)
        param_v  = np.array(next(self.fstream).split(), dtype=float)
        swap     =          next(self.fstream).strip() != '0'

        if degen:
            result = surface_factory.disc(r=r, center=center, xaxis=x_axis, normal=z_axis, type='radial')
        else:
            if not(np.allclose(np.diff(angles), pi/2, atol=1e-10)):
                raise RuntimeError('Unknown square parametrization of disc elementary surface')
            result = surface_factory.disc(r=r, center=center, xaxis=x_axis, normal=z_axis, type='square')
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def plane(self):
        dim        = int(     self.read_next_non_whitespace().strip())
        center     = np.array(next(self.fstream).split(), dtype=float)
        normal     = np.array(next(self.fstream).split(), dtype=float)
        x_axis     = np.array(next(self.fstream).split(), dtype=float)
        finite     =          next(self.fstream).strip() != '0'
        if finite:
            param_u= np.array(next(self.fstream).split(), dtype=float)
            param_v= np.array(next(self.fstream).split(), dtype=float)
        else:
            param_u= [-state.unlimited, +state.unlimited]
            param_v= [-state.unlimited, +state.unlimited]
        swap       =          next(self.fstream).strip() != '0'

        result = Surface() * [param_u[1]-param_u[0], param_v[1]-param_v[0]] + [param_u[0],param_v[0]]
        result.rotate(rotate_local_x_axis(x_axis, normal))
        result = flip_and_move_plane_geometry(result,center,normal)
        result.reparam(param_u, param_v)
        if(swap):
            result.swap()
        return result

    def torus(self):
        dim      = int(     self.read_next_non_whitespace().strip())
        r2       = float(   next(self.fstream).strip())
        r1       = float(   next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(), dtype=float)
        z_axis   = np.array(next(self.fstream).split(), dtype=float)
        x_axis   = np.array(next(self.fstream).split(), dtype=float)
        select_out=         next(self.fstream).strip() != '0' # I have no idea what this does :(
        param_u  = np.array(next(self.fstream).split(), dtype=float)
        param_v  = np.array(next(self.fstream).split(), dtype=float)
        swap     =          next(self.fstream).strip() != '0'

        result = surface_factory.torus(minor_r=r1, major_r=r2, center=center, normal=z_axis, xaxis=x_axis)
        result.reparam(param_u, param_v)
        if(swap):
            result.swap()
        return result

    def sphere(self):
        dim      = int(     self.read_next_non_whitespace().strip())
        r        = float(   next(self.fstream).strip())
        center   = np.array(next(self.fstream).split(), dtype=float)
        z_axis   = np.array(next(self.fstream).split(), dtype=float)
        x_axis   = np.array(next(self.fstream).split(), dtype=float)
        param_u  = np.array(next(self.fstream).split(), dtype=float)
        param_v  = np.array(next(self.fstream).split(), dtype=float)
        swap     =          next(self.fstream).strip() != '0'

        result = surface_factory.sphere(r=r, center=center, xaxis=x_axis, zaxis=z_axis).swap()
        if swap:
            result.swap()
        result.reparam(param_u, param_v)
        return result

    def splines(self, pardim):
        cls    = G2.classes[pardim-1]

        _, rational = self.read_next_non_whitespace().strip().split()
        rational = bool(int(rational))

        bases = [self.read_basis() for _ in range(pardim)]
        ncps = 1
        for b in bases:
            ncps *= b.num_functions()

        cps = [tuple(map(float, next(self.fstream).split()))
               for _ in range(ncps)]

        args = bases + [cps, rational]
        return cls(*args)

    def surface_of_linear_extrusion(self):
        dim      = int(      self.read_next_non_whitespace().strip())
        crv      = self.splines(1)
        normal   = np.array(self.read_next_non_whitespace().split(), dtype=float)
        finite   =          next(self.fstream).strip() != '0'
        param_u  = np.array(next(self.fstream).split(), dtype=float)
        if finite:
            param_v=np.array(next(self.fstream).split(), dtype=float)
        else:
            param_v=[-state.unlimited, +state.unlimited]
        swap     =          next(self.fstream).strip() != '0'

        result = surface_factory.extrude(crv + normal*param_v[0], normal*(param_v[1]-param_v[0]))
        result.reparam(param_u, param_v)

        if swap:
            result.swap()
        return result


    def bounded_surface(self):
        objtype = int( next(self.fstream).strip() )

        # create the underlying surface which all trimming curves are to be applied
        if objtype in G2.g2_generators:
            constructor = getattr(self, G2.g2_generators[objtype].__name__)
            surface = constructor()
        elif objtype == 200:
            surface = self.splines(2)
        else:
            raise IOError('Unsopported trimmed surface or malformed input file')

        # for all trimming loops
        numb_loops = int( self.read_next_non_whitespace() )
        all_loops  = []
        for i in range(numb_loops):

            # for all cuve pieces of that loop
            numb_crvs, space_epsilon = next(self.fstream).split()
            state.parametric_absolute_tolerance = float(space_epsilon)
            one_loop = []
            for j in range(int(numb_crvs)):

                # read a physical and parametric representation of the same curve
                _, parameter_curve_type, space_curve_type = map(int, self.read_next_non_whitespace().split())
                two_curves = []
                for crv_type in [parameter_curve_type, space_curve_type]:
                    if crv_type in G2.g2_generators:
                        constructor = getattr(self, G2.g2_generators[crv_type].__name__)
                        crv = constructor()
                    elif crv_type == 100:
                        crv = self.splines(1)
                    else:
                        raise IOError('Unsopported trimming curve or malformed input file')
                    two_curves.append(crv)

                # only keep the parametric version (re-generate physical one if we need it)
                one_loop.append(two_curves[0])
                self.trimming_curves.append(two_curves[1])
            all_loops.append(one_loop)

        return TrimmedSurface(surface.bases[0], surface.bases[1], surface.controlpoints, surface.rational, all_loops, raw=True)

    g2_type = [100, 200, 700] # curve, surface, volume identifiers
    classes = [Curve, Surface, Volume]
    g2_generators = {120:line, 130:circle, 140:ellipse,
                     260:cylinder, 292:disc, 270:sphere, 290:torus, 250:plane,
                     210:bounded_surface, 261:surface_of_linear_extrusion} #, 280:cone

    def __init__(self, filename):
        if filename[-3:] != '.g2':
            filename += '.g2'
        self.filename = filename
        self.trimming_curves = []

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
            self.fstream.write(' '.join('%.16g' % k for k in b.knots))
            self.fstream.write('\n')

        savetxt(self.fstream, obj.controlpoints.reshape(-1, obj.dimension + obj.rational, order='F'),
                fmt='%.16g', delimiter=' ', newline='\n')

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
            objtype, major, minor, patch = map(int, line.split())
            if (major, minor, patch) != (1, 0, 0):
                raise IOError('Unknown G2 format')

            # if obj type is in factory methods (cicle, torus etc), create it now
            if objtype in G2.g2_generators:
                constructor = getattr(self, G2.g2_generators[objtype].__name__)
                result.append( constructor() )
                continue

            # for "normal" splines (Curves, Surfaces, Volumes) create it now
            pardim = [i for i in range(len(G2.g2_type)) if G2.g2_type[i] == objtype]
            if not pardim:
                raise IOError('Unknown G2 object type {}'.format(objtype))
            pardim = pardim[0] + 1
            result.append(self.splines(pardim))

        return result

    def read_basis(self):
        ncps, order = map(int, next(self.fstream).split())
        kts = list(map(float, next(self.fstream).split()))
        return BSplineBasis(order, kts, -1)

    def __exit__(self, exc_type, exc_value, traceback):
        if hasattr(self, 'fstream'):
            self.fstream.close()

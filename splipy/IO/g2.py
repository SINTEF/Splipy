import numpy as np
from itertools import chain, product
from splipy import BSplineBasis, Curve, Surface, Volume, SplineObject
from .master import MasterIO


class G2(MasterIO):

    g2_type = [100, 200, 700] # curve, surface, volume identifiers
    classes = [Curve, Surface, Volume]

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

            objtype, major, minor, patch = map(int, line.split(' '))
            if (major, minor, patch) != (1, 0, 0):
                raise IOError('Unknown G2 format')

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


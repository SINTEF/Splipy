import numpy as np
from itertools import chain, product
from splipy import Curve, Surface, BSplineBasis
from .master import MasterIO
import splipy.state as state
from rhino3dm import Brep, File3dm, NurbsCurve, NurbsSurface

class ThreeDM(MasterIO):

    def __init__(self, filename):
        if filename[-4:] != '.3dm':
            filename += '.3dm'
        self.filename = filename
        self.trimming_curves = []

    def __enter__(self):
        return self

    def write(self, obj):
        raise IOError('Writing to 3DM not supported')

    def read(self):
        if not hasattr(self, 'fstream'):
            self.onlywrite = False
            self.fstream   = File3dm.Read(self.filename)

        if self.onlywrite:
            raise IOError('Could not read from file %s' % (self.filename))

        result = []

        for obj in self.fstream.Objects:
            if type(obj.Geometry) is Brep:
                for idx in range(0,len(obj.Geometry.Faces)):
                    nsrf = obj.Geometry.Faces[idx].UnderlyingSurface().ToNurbsSurface()
                    result.append(self.read_surface(nsrf))

            if type(obj.Geometry) is NurbsCurve:
                result.append(self.read_curve(obj.Geometry))

            if type(obj.Geometry) is NurbsSurface:
                result.append(self.read_surface(obj.Geometry))

        return result

    def read_surface(self, nsrf):
        knotsu = [0]
        for i in nsrf.KnotsU:
            knotsu.append(i)
        knotsu.append(knotsu[len(knotsu)-1])
        knotsu[0] = knotsu[1]

        knotsv = [0]
        for i in nsrf.KnotsV:
            knotsv.append(i)
        knotsv.append(knotsv[len(knotsv)-1])
        knotsv[0] = knotsv[1]

        basisu = BSplineBasis(nsrf.OrderU, knotsu, -1)
        basisv = BSplineBasis(nsrf.OrderV, knotsv, -1)
        cpts = []

        cpts = np.ndarray((nsrf.Points.CountU*nsrf.Points.CountV, 3 + nsrf.IsRational))
        for v in range(0,nsrf.Points.CountV):
            for u in range(0,nsrf.Points.CountU):
                cpts[u+v*nsrf.Points.CountU,0] = nsrf.Points[u,v].X
                cpts[u+v*nsrf.Points.CountU,1] = nsrf.Points[u,v].Y
                cpts[u+v*nsrf.Points.CountU,2] = nsrf.Points[u,v].Z
                if nsrf.IsRational:
                    cpts[u+v*nsrf.Points.CountU,3] = 1
                    cpts[u+v*nsrf.Points.CountU,:] *= nsrf.Points[u,v].W

        return Surface(basisu, basisv, cpts, nsrf.IsRational)

    def read_curve(self, ncrv):
        knots = [0]
        for i in ncrv.Knots:
            knots.append(i)
        knots[0] = knots[1]
        knots.append(knots[len(knots)-1])
        basis = BSplineBasis(ncrv.Order, knots, -1)
        cpts = []

        cpts = np.ndarray((len(ncrv.Points), ncrv.Dimension + ncrv.IsRational))
        for u in range(0,len(ncrv.Points)):
            cpts[u,0] = ncrv.Points[u].X
            cpts[u,1] = ncrv.Points[u].Y
            if ncrv.Dimension > 2:
                cpts[u,2] = ncrv.Points[u].Z
            if ncrv.IsRational:
                cpts[u,-1] = 1
                cpts[u,:] *= ncrv.Points[u].W

        return Curve(basis, cpts, ncrv.IsRational)

    def __exit__(self, exc_type, exc_value, traceback):
        # Apperently File3DM objects don't need to dedicated cleanup/close code
        pass

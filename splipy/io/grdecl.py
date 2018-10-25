import numpy as np
from itertools import islice, product, chain
from splipy import Curve, Surface, Volume, SplineObject, BSplineBasis
from splipy import curve_factory, surface_factory, volume_factory
from .master import MasterIO
import re
import warnings
from scipy.spatial import Delaunay
from scipy.spatial.qhull import QhullError
from tqdm import tqdm


class Box(object):

    def __init__(self, x):
        self.x = x


class DiscontBoxMesh(object):

    def __init__(self, n, coord, zcorn):
        nx, ny, nz = n

        X  = np.empty(n + 1, dtype=object)
        Xz = np.zeros((nx + 1, ny + 1, 2 * nz, 3))
        cells = np.empty(n, dtype=object)

        for i, j, k in product(range(nx), range(ny), range(nz)):
            x = []
            for k0, j0, i0 in product(range(2), repeat=3):
                # Interpolate to find the x,y values of this point
                zmin, zmax = coord[i+i0, j+j0, :, 2]
                z  = zcorn[2*i+i0, 2*j+j0, 2*k+k0]
                t = (z - zmax) / (zmin - zmax)
                point = coord[i+i0, j+j0, 0] * t + coord[i+i0, j+j0, 1] * (1 - t)
                x.append(point)

                if X[i+i0,j+j0,k+k0] is None:
                    X[i+i0,j+j0,k+k0] = [point]
                else:
                    X[i+i0,j+j0,k+k0].append(point)
                Xz[i+i0,j+j0,2*k+k0,:] = point

            cells[i,j,k] = Box(x)

        self.X = X
        self.Xz = Xz
        self.n = n

        def hull_or_none(x):
            try:
                return Delaunay(x)
            except QhullError:
                return None

        self.plane_hull = np.array([
            [Delaunay(np.reshape(coord[i:i+2, j:j+2, :, :], (8,3))) for j in range(ny)]
            for i in range(nx)
        ], dtype=object)

        self.hull = np.array([
            [[hull_or_none(cell.x) for cell in cell_tower] for cell_tower in cells_tmp]
            for cells_tmp in cells
        ], dtype=object)

    def cell_at(self, x, guess=None):
        # First, find the 'tower' containing x
        check = -1
        if guess is not None:
            i, j, _ = guess
            check = self.plane_hull[i,j].find_simplex(x)
        if check >= -1:
            for (i, j), hull in np.ndenumerate(self.plane_hull):
                check = hull.find_simplex(x)
                if check >= 0: break

        assert check >= 0

        # Find the correct cell in the 'tower'
        check = -1
        if guess is not None:
            _, _, k = guess
            check = self.hull[i,j,k].find_simplex(x)
        if check >= -1:
            for k, hull in enumerate(self.hull[i,j]):
                if hull is None: continue
                check = hull.find_simplex(x)
                if check >= 0: break

        assert check >= -1
        return i, j, k

    def get_c0_avg(self):
        """Compute best-approximation vertices for a continuous mesh by averaging the location of all
        corners that 'should' coincide.
        """
        return np.array([[[np.mean(k,axis=0) for k in j] for j in i] for i in self.X])

    def get_discontinuous_all(self):
        """Return a list of vertices suitable for a fully discontinuous mesh."""
        return list(chain.from_iterable(xs[::-1] for xs in self.X.T.flat))

    def get_discontinuous_z(self):
        """Return a list of vertices suitable for a mixed continuity mesh."""
        return self.Xz


class GRDECL(MasterIO):

    def __init__(self, filename):
        if not filename.endswith('.grdecl'):
            filename += '.grdecl'
        self.filename = filename
        self.attribute = {}

    def __enter__(self):
        self.fstream    = open(self.filename, 'r')
        self.line_number = 0
        return self

    def read_specgrid(self):
        args = next(self.fstream).strip().split()
        return np.array(args[:3], dtype=np.int32)

    def read_coord(self):
        nx, ny = self.n[:2]
        ans = np.zeros((nx + 1, ny + 1, 2, 3))
        for j, i in product(range(ny+1), range(nx+1)):
            args = next(self.fstream).split()
            ans[i,j,0,:] = np.array(args[:3], dtype=np.float64)
            ans[i,j,1,:] = np.array(args[3:], dtype=np.float64)
        return ans

    def read_zcorn(self):
        ntot = np.prod(self.n)*8
        numbers = []
        while len(numbers) < ntot:
            numbers += next(self.fstream).split()
        return np.reshape(np.array(numbers, dtype=np.float64), self.n*2, order='F')

    def cell_property(self, dtype=np.float64):
        ntot = np.prod(self.n)
        numbers = []
        while len(numbers) < ntot:
            numbers += next(self.fstream).split()
        return np.array(numbers, dtype=dtype)

    def read(self):
        for line in self.fstream:
            line = line.strip().lower()
            if line == 'specgrid':
                self.n = self.read_specgrid()
            elif line == 'coord':
                self.coord = self.read_coord()
            elif line == 'zcorn':
                self.zcorn = self.read_zcorn()
            elif line in {'actnum', 'permx', 'permy', 'permz', 'poro', 'satnum', 'rho'}:
                dtype = np.int32 if line in {'actnum', 'satnum'} else np.float64
                self.attribute[line] = self.cell_property(dtype)
            elif line in {'grid', '/', ''} or line.startswith('--'):
                pass
            elif not re.match('[0-9]', line[0]):
                warnings.showwarning(
                    'Unkown keyword "{}" encountered in file'.format(line.split()[0]),
                    SyntaxWarning, self.filename, self.line_number, line=[],
                )
            else:
                pass # silently skip large number blocks

        self.raw = DiscontBoxMesh(self.n, self.coord, self.zcorn)

        nx, ny, nz = self.n

        # Create the C0-mesh
        X = self.raw.get_c0_avg()
        b1 = BSplineBasis(2, [0] + [i/nx for i in range(nx+1)] + [1])
        b2 = BSplineBasis(2, [0] + [i/ny for i in range(ny+1)] + [1])
        b3 = BSplineBasis(2, [0] + [i/nz for i in range(nz+1)] + [1])
        self.c0_vol = volume_factory.interpolate(X, [b1, b2, b3])

        # Create the C^{-1} mesh
        Xm1 = self.raw.get_discontinuous_all()
        b1 = BSplineBasis(2, sorted(list(range(self.n[0]+1))*2))
        b2 = BSplineBasis(2, sorted(list(range(self.n[1]+1))*2))
        b3 = BSplineBasis(2, sorted(list(range(self.n[2]+1))*2))
        discont_vol = Volume(b1, b2, b3, Xm1)

        # Create mixed discontinuity mesh: C^0, C^0, C^{-1}
        Xz = self.raw.get_discontinuous_z()
        b1 = BSplineBasis(2, sorted(list(range(self.n[0]+1))+[0,self.n[0]]))
        b2 = BSplineBasis(2, sorted(list(range(self.n[1]+1))+[0,self.n[1]]))
        b3 = BSplineBasis(2, sorted(list(range(self.n[2]+1))*2))
        self.true_vol = Volume(b1, b2, b3, Xz, raw=True)

        return self.true_vol, discont_vol, self.c0_vol

    def texture(self, p, n, method='full'):
        # Set the dimensions of geometry and texture map
        ngeom    = np.floor(self.n / (p-1))
        ntexture = np.floor(self.n * n)
        ngeom    = ngeom.astype(np.int32)
        ntexture = ntexture.astype(np.int32)

        # Create the geometry
        ngx, ngy, ngz = ngeom
        b1 = BSplineBasis(p, [0]*(p-1) + [i/ngx for i in range(ngx+1)] + [1]*(p-1))
        b2 = BSplineBasis(p, [0]*(p-1) + [i/ngy for i in range(ngy+1)] + [1]*(p-1))
        b3 = BSplineBasis(p, [0]*(p-1) + [i/ngz for i in range(ngz+1)] + [1]*(p-1))

        nx, ny, nz = self.n
        l2_fit = surface_factory.least_square_fit
        u, v, w = np.linspace(0, 1, nx+1), np.linspace(0, 1, ny+1), np.linspace(0, 1, nz+1)

        if method == 'full':
            bottom = l2_fit(self.c0_vol[:,:, 0,:].squeeze(), [b1, b2], [u, v])
            top    = l2_fit(self.c0_vol[:,:,-1,:].squeeze(), [b1, b2], [u, v])
            left   = l2_fit(self.c0_vol[ 0,:,:,:].squeeze(), [b2, b3], [v, w])
            right  = l2_fit(self.c0_vol[-1,:,:,:].squeeze(), [b2, b3], [v, w])
            front  = l2_fit(self.c0_vol[:, 0,:,:].squeeze(), [b1, b3], [u, w])
            back   = l2_fit(self.c0_vol[:,-1,:,:].squeeze(), [b1, b3], [u, w])
            volume = volume_factory.edge_surfaces([left, right, front, back, bottom, top])

        elif method == 'z':
            bottom = l2_fit(self.c0_vol[:,:, 0,:].squeeze(), [b1, b2], [u, v])
            top    = l2_fit(self.c0_vol[:,:,-1,:].squeeze(), [b1, b2], [u, v])
            volume = volume_factory.edge_surfaces([bottom, top])
            volume.set_order(p)
            volume.refine(ngx - 1, direction='w')

        # Point-to-cell mapping
        # TODO: Optimize more
        u = [np.linspace(0,1,n) for n in ntexture]
        points = volume(*u).reshape(-1, 3)
        cellids = np.zeros(points.shape[:-1], dtype=int)
        cell = None
        for ptid, point in enumerate(tqdm(points, desc='Inverse mapping')):
            i, j, k = cell = self.raw.cell_at(point, guess=cell)
            cellid = i*ny*nz + j*nz + k
            cellids[ptid] = cellid

        cellids = cellids.reshape(tuple(len(param) for param in u))

        all_textures = {}
        for name in self.attribute:
            data = self.attribute[name][cellids]

            # TODO: This flattens the image if it happens to be 3D (or higher...)
            # do we need a way to communicate the structure back to the caller?
            data = data.reshape(-1, data.shape[-1])

            # TODO: This normalizes the image,
            # but we need a way to communicate the ranges back to the caller
            a, b = min(data.flat), max(data.flat)
            data = ((data - a) / (b - a) * 255).astype(np.uint8)

            all_textures[name] = data

        return volume, all_textures

    def __exit__(self, exc_type, exc_value, traceback):
        self.fstream.close()

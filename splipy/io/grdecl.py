import numpy as np
from itertools import islice
from splipy import Curve, Surface, Volume, SplineObject, BSplineBasis
from splipy import curve_factory, surface_factory, volume_factory
from .master import MasterIO
import re
import warnings
from scipy.spatial import Delaunay
from scipy.spatial.qhull import QhullError

class box(object):

    def __init__(self, x):
        self.x = x

class disconnected_box_mesh(object):

    def __init__(self, n):
        self.n = n
        self.cells = []

    def populate(self, coord, zcorn):
        # self.X = [[[[None]]*(self.n[2]+1)]*(self.n[1]+1)]*(self.n[0]+1)
        self.X  = np.empty(self.n+1, dtype=object)
        self.Xz = np.zeros((self.n[0]+1, self.n[1]+1, 2*self.n[2], 3))
        # print(self.X.shape)
        for i in range(self.n[0]):
            self.cells.append([[]])
            for j in range(self.n[1]):
                self.cells[i].append([])
                for k in range(self.n[2]):
                    x = []
                    for k0 in range(2):
                        for j0 in range(2):
                            for i0 in range(2):
                                z0 = coord[i+i0, j+j0, 0, 2]
                                z1 = coord[i+i0, j+j0, 1, 2]
                                z  = zcorn[2*i+i0, 2*j+j0, 2*k+k0]
                                t  = (z-z1) / (z0-z1)
                                dx = coord[i+i0, j+j0, 1, :2] - coord[i+i0, j+j0, 0, :2]
                                xy = coord[i+i0, j+j0, 0, :2] + dx*t
                                x.append(np.array([xy[0], xy[1], z]))
                                if self.X[i+i0,j+j0,k+k0] == None:
                                    self.X[i+i0,j+j0,k+k0] = [x[-1]]
                                else:
                                    self.X[i+i0,j+j0,k+k0].append(x[-1])
                                self.Xz[i+i0,j+j0,2*k+k0,:] = x[-1]
                    self.cells[i][j].append(box(x))
        self.X = np.array(self.X)

        # print('Building convex hulls for searching')
        # self.hull    = [[[Delaunay(mycell.x) for mycell in heights] for heights in columns] for columns in self.cells]
        # self.ij_hull = [[Delaunay(np.reshape(coord[i:i+2,j:j+2,:,:], (8,3))) for j in range(self.n[1])] for i in range(self.n[0])]
        self.hull    = [[[]]]
        self.ij_hull =  [[]]
        for i in range(self.n[0]):
            self.hull.append([[]])
            self.ij_hull.append([])
            for j in range(self.n[1]):
                self.hull[i].append([])
                for k in range(self.n[2]):
                    try:
                        self.hull[i][j].append(Delaunay(self.cells[i][j][k].x))
                    except QhullError: # degenerate volume
                        self.hull[i][j].append(None)
                self.ij_hull[i].append(Delaunay(np.reshape(coord[i:i+2,j:j+2,:,:], (8,3))))

    def cell_at(self, x):
        i,j = np.where(np.array([[hull.find_simplex(x) for hull in columns] for columns in self.ij_hull]) >= 0)
        i,j = i[0], j[0]
        k   = np.where(np.array([ hull.find_simplex(x) if hull is not None else -1 for hull in self.hull[i][j]]) >= 0)
        return (i,j,k[0][0])

    def get_c0_avg(self):
        return np.array([[[np.mean(k,axis=0) for k in j] for j in i] for i in self.X])

    def get_discontinuous_all(self):
        ans = []
        for k in range(self.n[2]+1):
            for j in range(self.n[1]+1):
                for i in range(self.n[0]+1):
                    while len(self.X[i][j][k])>0:
                        ans.append(self.X[i][j][k].pop())
        return ans

    def get_discontinuous_z(self):
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

    def readline(self):
        self.line_number += 1
        return self.fstream.readline()

    def specgrid(self):
        args = self.readline()
        args = args.strip().split()
        return np.array(args[:3], dtype=np.int32)

    def coord(self):
        n = (self.n[0]+1, self.n[1]+1,2,3)
        ans = np.zeros(n)
        for j in range(self.n[1]+1):
            for i in range(self.n[0]+1):
                args = self.readline().split()
                ans[i,j,0,:] = np.array(args[:3], dtype=np.float64)
                ans[i,j,1,:] = np.array(args[3:], dtype=np.float64)
        return ans

    def zcorn(self):
        ntot = np.prod(self.n)*8
        ans  = np.zeros(self.n*2)
        numbers = []
        while len(numbers) < ntot:
            numbers += self.readline().split()
        # this tripple nested for loops can probably be exchanged by a numpy
        # reshape-call
        n = 0
        for k in range(2*self.n[2]):
            for j in range(2*self.n[1]):
                for i in range(2*self.n[0]):
                    ans[i,j,k] = float(numbers[n])
                    n += 1
        return ans

    def cell_property(self, dtype=np.float64):
        ntot = np.prod(self.n)
        numbers = []
        while len(numbers) < ntot:
            numbers += self.readline().split()
        return np.array(numbers, dtype=dtype)

    def read(self):
        line = self.readline()
        while not line == '':
            line = line.strip().lower()
            if   line == 'specgrid': self.n      = self.specgrid()
            elif line == 'coord':    self.coord  = self.coord()
            elif line == 'zcorn':    self.zcorn  = self.zcorn()
            elif line == 'actnum':   self.attribute['actnum'] = self.cell_property(np.int32)
            elif line == 'permx':    self.attribute['permx']  = self.cell_property()
            elif line == 'permy':    self.attribute['permy']  = self.cell_property()
            elif line == 'permz':    self.attribute['permz']  = self.cell_property()
            elif line == 'poro':     self.attribute['poro']   = self.cell_property()
            elif line == 'satnum':   self.attribute['satnum'] = self.cell_property(np.int32)
            elif line == 'rho':      self.attribute['rho']    = self.cell_property()
            elif line == 'grid':     pass # this signals the start of interesting stuff
            elif line == '/':        pass # end read dataset
            elif line == '' :        pass # endline
            elif line[:2] == '--' :  pass # comment line
            elif not re.match('[0-9]',line[0]):
                # k = self.line_number
                warnings.showwarning('Unkown keyword "{}" encountered in file'.format(line.split()[0]), SyntaxWarning, self.filename, self.line_number, line=[])
            else:
                pass # silently skip large number blocks
            line = self.readline()

        self.raw = disconnected_box_mesh(self.n)
        self.raw.populate(self.coord, self.zcorn)

        # create the C0-mesh
        X = self.raw.get_c0_avg()
        # print(X.shape)
        # print(X)
        b1 = BSplineBasis(2, [0]+[i/self.n[0] for i in range(self.n[0]+1)]+[1])
        b2 = BSplineBasis(2, [0]+[i/self.n[1] for i in range(self.n[1]+1)]+[1])
        b3 = BSplineBasis(2, [0]+[i/self.n[2] for i in range(self.n[2]+1)]+[1])
        c0_vol = volume_factory.interpolate(X, [b1,b2,b3])

        # create the C^{-1} mesh
        Xm1 = self.raw.get_discontinuous_all()
        b1 = BSplineBasis(2, sorted(list(range(self.n[0]+1))*2))
        b2 = BSplineBasis(2, sorted(list(range(self.n[1]+1))*2))
        b3 = BSplineBasis(2, sorted(list(range(self.n[2]+1))*2))
        discont_vol = Volume(b1,b2,b3, Xm1)

        # create mixed discontinuity mesh: C^0, C^0, C^{-1}
        Xz = self.raw.get_discontinuous_z()
        b1 = BSplineBasis(2, sorted(list(range(self.n[0]+1))+[0,self.n[0]]))
        b2 = BSplineBasis(2, sorted(list(range(self.n[1]+1))+[0,self.n[1]]))
        b3 = BSplineBasis(2, sorted(list(range(self.n[2]+1))*2))
        self.true_vol = Volume(b1,b2,b3, Xz, raw=True)

        return self.true_vol, discont_vol, c0_vol

    def texture(self, p, n):
        # set the dimensions of geometry and texture map
        ngeom    = np.floor(self.n / (p-1))
        ntexture = np.floor(self.n * n)
        ngeom    = ngeom.astype(   np.int32)
        ntexture = ntexture.astype(np.int32)

        # create the geometry
        b1 = BSplineBasis(p, [0]*(p-1)+[i/ngeom[0] for i in range(ngeom[0]+1)]+[1]*(p-1))
        b2 = BSplineBasis(p, [0]*(p-1)+[i/ngeom[1] for i in range(ngeom[1]+1)]+[1]*(p-1))
        b3 = BSplineBasis(p, [0]*(p-1)+[i/ngeom[2] for i in range(ngeom[2]+1)]+[1]*(p-1))
        l2_fit = surface_factory.least_square_fit
        v,u = np.linspace(0,1,self.n[1]+1), np.linspace(0,1,self.n[0]+1)
        w   = sorted(list(np.linspace(0,1,self.n[2]-1, endpoint=False))*2+[0,1])
        X = self.true_vol[:,:, 0,:].squeeze()
        bottom = l2_fit(self.true_vol[:,:, 0,:].squeeze(), [b1,b2], [u, v])
        top    = l2_fit(self.true_vol[:,:,-1,:].squeeze(), [b1,b2], [u, v])
        left   = l2_fit(self.true_vol[ 0,:,:,:].squeeze(), [b2,b3], [v, w])
        right  = l2_fit(self.true_vol[-1,:,:,:].squeeze(), [b2,b3], [v, w])
        front  = l2_fit(self.true_vol[:, 0,:,:].squeeze(), [b1,b3], [u, w])
        back   = l2_fit(self.true_vol[:,-1,:,:].squeeze(), [b1,b3], [u, w])
        volume = volume_factory.edge_surfaces([left, right, front, back, bottom, top])
        # volume = volume_factory.edge_surfaces([bottom, top])
        # volume.raise_order(0,0,p-2)
        # volume.refine(0,0,ngeom[2]-1)

        return [volume]
        # return [top, bottom, front, back, left, right]


        all_textures = {}
        for name in self.attribute:
            one_texture = np.zeros(ntexture)
            u = [np.linspace(0,1,n) for n in ntexture]
            x = volume(*u)
            x = np.reshape(x, [np.prod(ntexture),3])
            for I in range(np.prod(ntexture)):
                print(x[I,:])
                i,j,k = self.raw.cell_at(x[I,:])
                J = i*self.n[1]*self.n[2] + j*self.n[2] + k
                one_texture[I] = self.attribute[name][J]
            all_textures[name] = one_texture
        return volume, all_textures



    def __exit__(self, exc_type, exc_value, traceback):
        self.fstream.close()

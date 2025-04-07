# --- Automatic generated test file  ---
# Generator    : generate_derivative_spline.py
# Date         : 2017-01-23
# Git revision : b458827cf91658b8da98744b230af11b9c95dd94

import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis
from math import sqrt
import unittest


class TestDerivativeSpline(unittest.TestCase):
    def test_curve_2D_p2(self):
        controlpoints = np.array([[  1.,  -1.],
       [ 35.,   0.],
       [ 61.,  -4.],
       [ 95.,  -5.]])
        basis0 = BSplineBasis(2, np.array([ 0. ,  0. ,  0.7,  1.7,  3. ,  3. ]))
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 1)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_2D_p23(self):
        controlpoints = np.array([[   1.,   -3.],
       [  52.,   -2.],
       [  96.,   -5.],
       [   2.,   20.],
       [  50.,   21.],
       [  96.,   22.],
       [   2.,   47.],
       [  48.,   52.],
       [ 100.,   54.],
       [  -2.,   75.],
       [  49.,   78.],
       [  98.,   73.],
       [  -2.,  104.],
       [  45.,   98.],
       [  98.,   95.]])
        basis0 = BSplineBasis(2, np.array([ 0. ,  0. ,  1.2,  2. ,  2. ]))
        basis1 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  1.1,  1.8,  3. ,  3. ,  3. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 1)
        self.assertEqual(surf2.order(direction=1), 3)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 2)
        self.assertEqual(surf3.order(direction=1), 2)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_curve_2D_p4_C0_periodic(self):
        controlpoints = np.array([[  99.,    1.],
       [  73.,   71.],
       [  -2.,  102.],
       [ -73.,   66.],
       [-100.,    3.],
       [ -69.,  -68.],
       [   3., -104.],
       [  70.,  -73.]])
        basis0 = BSplineBasis(4, np.array([-0.7,  0. ,  0. ,  0. ,  0.6,  2.2,  2.8,  4.1,  5.3,  6. ,  6. ,  6. ,  6.6]),0)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 3)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_2D_p22_C0_periodic(self):
        controlpoints = np.array([[  64.,    3.],
       [   2.,   58.],
       [ -56.,   -5.],
       [   0.,  -65.],
       [  76.,   -2.],
       [   3.,   77.],
       [ -71.,   -2.],
       [   0.,  -74.],
       [  86.,    0.],
       [  -4.,   87.],
       [ -91.,   -4.],
       [   2.,  -92.],
       [ 104.,   -2.],
       [  -5.,   96.],
       [ -98.,   -4.],
       [  -2.,  -98.]])
        basis0 = BSplineBasis(2, np.array([-0.6,  0. ,  1.1,  2. ,  3.4,  4. ,  5.1]),0)
        basis1 = BSplineBasis(2, np.array([ 0. ,  0. ,  1.2,  1.7,  3. ,  3. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 1)
        self.assertEqual(surf2.order(direction=1), 2)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 2)
        self.assertEqual(surf3.order(direction=1), 1)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_curve_3D_p3(self):
        controlpoints = np.array([[  -1.,   -1.,    0.],
       [  29.,    1.,   -5.],
       [  68.,    0.,    1.],
       [ 104.,   -1.,    4.]])
        basis0 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  1.3,  2. ,  2. ,  2. ]))
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 2)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p43(self):
        controlpoints = np.array([[   4.,    4.,   -4.],       [  13.,    1.,    3.],       [  36.,   -5.,    3.],       [  48.,    0.,   -4.],       [  64.,    2.,    2.],       [  86.,    3.,   -4.],       [  95.,    4.,    1.],       [  -4.,   23.,   -1.],       [  17.,   20.,    0.],       [  34.,   18.,   -3.],       [  47.,   22.,    3.],       [  66.,   19.,   -5.],       [  80.,   18.,    0.],       [  95.,   17.,   -1.],       [   3.,   42.,    3.],       [  14.,   41.,    0.],       [  30.,   37.,   -2.],       [  52.,   41.,   -5.],       [  68.,   37.,    4.],       [  86.,   40.,    0.],       [ 104.,   36.,    2.],       [   3.,   57.,    4.],       [  12.,   57.,    1.],       [  37.,   60.,    0.],       [  53.,   55.,   -4.],       [  63.,   60.,    2.],       [  86.,   59.,    2.],       [  95.,   61.,    3.],       [  -4.,   81.,   -4.],       [  18.,   84.,   -1.],       [  28.,   79.,    0.],       [  48.,   75.,   -2.],       [  71.,   83.,    0.],       [  83.,   80.,   -4.],       [ 104.,   84.,   -4.],       [  -1.,  104.,   -3.],       [  13.,  100.,    2.],       [  37.,  102.,   -5.],       [  51.,  101.,   -4.],       [  66.,   97.,   -3.],       [  81.,   95.,   -4.],       [ 103.,   97.,    3.]])
        basis0 = BSplineBasis(4, np.array([ 0. ,  0. ,  0. ,  0. ,  0.8,  1.7,  2.7,  4. ,  4. ,  4. ,  4. ]))
        basis1 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  1.2,  1.9,  2.8,  4. ,  4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 3)
        self.assertEqual(surf2.order(direction=1), 3)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 4)
        self.assertEqual(surf3.order(direction=1), 2)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p434(self):
        controlpoints = np.array([[   3.,    1.,   -3.],       [  23.,   -1.,   -2.],       [  41.,    0.,    1.],       [  56.,    3.,   -3.],       [  75.,    1.,    4.],       [ 101.,    0.,   -1.],       [   2.,   33.,    3.],       [  22.,   35.,   -2.],       [  38.,   35.,   -1.],       [  63.,   37.,   -2.],       [  83.,   31.,    4.],       [ 103.,   36.,   -2.],       [  -4.,   64.,    1.],       [  20.,   66.,   -3.],       [  43.,   68.,   -5.],       [  63.,   63.,   -1.],       [  76.,   68.,    2.],       [ 104.,   70.,   -5.],       [   4.,   96.,   -4.],       [  19.,  101.,    0.],       [  42.,  101.,    2.],       [  64.,  100.,   -2.],       [  82.,   97.,    3.],       [ 102.,  101.,    3.],       [   0.,   -4.,   21.],       [  17.,   -5.,   17.],       [  42.,    1.,   20.],       [  56.,    0.,   14.],       [  78.,   -4.,   16.],       [ 103.,   -1.,   13.],       [   3.,   38.,   21.],       [  22.,   32.,   14.],       [  35.,   29.,   16.],       [  63.,   35.,   15.],       [  76.,   36.,   13.],       [ 102.,   33.,   18.],       [  -4.,   68.,   16.],       [  23.,   67.,   20.],       [  37.,   67.,   18.],       [  59.,   69.,   15.],       [  75.,   70.,   12.],       [  95.,   66.,   19.],       [   1.,   98.,   19.],       [  24.,  103.,   12.],       [  40.,   98.,   18.],       [  61.,  102.,   21.],       [  84.,   97.,   18.],       [ 102.,   97.,   12.],       [   3.,   -4.,   28.],       [  15.,   -4.,   29.],       [  38.,    4.,   35.],       [  62.,   -1.,   31.],       [  79.,   -1.,   30.],       [  96.,    3.,   28.],       [   3.,   35.,   31.],       [  19.,   32.,   37.],       [  42.,   28.,   30.],       [  60.,   29.,   37.],       [  83.,   38.,   37.],       [ 103.,   35.,   37.],       [  -5.,   69.,   33.],       [  21.,   70.,   35.],       [  39.,   62.,   36.],       [  63.,   68.,   28.],       [  77.,   67.,   34.],       [  97.,   66.,   28.],       [   4.,  102.,   31.],       [  23.,  104.,   37.],       [  42.,   95.,   31.],       [  58.,   97.,   35.],       [  79.,   95.,   32.],       [ 100.,  101.,   31.],       [  -3.,    1.,   54.],       [  18.,   -5.,   50.],       [  35.,   -5.,   49.],       [  55.,    4.,   52.],       [  82.,    1.,   51.],       [  95.,   -3.,   48.],       [  -2.,   30.,   54.],       [  23.,   33.,   53.],       [  41.,   37.,   50.],       [  63.,   29.,   51.],       [  83.,   32.,   54.],       [  98.,   33.,   49.],       [  -3.,   64.,   48.],       [  22.,   70.,   50.],       [  40.,   66.,   47.],       [  59.,   69.,   46.],       [  75.,   63.,   46.],       [  95.,   66.,   54.],       [  -4.,  101.,   47.],       [  19.,   97.,   48.],       [  44.,   99.,   46.],       [  55.,   98.,   46.],       [  75.,   96.,   46.],       [  95.,   98.,   54.],       [   3.,    1.,   64.],       [  16.,   -5.,   70.],       [  43.,   -1.,   62.],       [  55.,   -4.,   68.],       [  81.,    0.,   69.],       [ 100.,   -3.,   68.],       [   1.,   30.,   66.],       [  21.,   36.,   62.],       [  39.,   29.,   62.],       [  56.,   33.,   62.],       [  84.,   28.,   66.],       [  96.,   32.,   63.],       [   2.,   61.,   69.],       [  18.,   66.,   67.],       [  36.,   65.,   68.],       [  60.,   70.,   65.],       [  77.,   71.,   66.],       [  97.,   68.,   65.],       [  -1.,  101.,   63.],       [  22.,   96.,   71.],       [  44.,  101.,   66.],       [  55.,   97.,   62.],       [  84.,  103.,   69.],       [  97.,  101.,   64.],       [  -4.,    2.,   78.],       [  24.,   -1.,   85.],       [  35.,   -3.,   82.],       [  63.,   -2.,   80.],       [  82.,    4.,   82.],       [  98.,   -1.,   82.],       [  -2.,   32.,   85.],       [  18.,   36.,   84.],       [  35.,   34.,   87.],       [  60.,   31.,   86.],       [  81.,   32.,   84.],       [ 104.,   35.,   78.],       [   2.,   65.,   84.],       [  24.,   62.,   86.],       [  42.,   70.,   87.],       [  57.,   69.,   84.],       [  81.,   66.,   86.],       [ 104.,   65.,   80.],       [  -3.,   97.,   88.],       [  17.,  101.,   86.],       [  43.,   96.,   78.],       [  58.,   95.,   84.],       [  75.,   96.,   86.],       [  98.,   97.,   83.],       [   2.,    3.,   95.],       [  21.,   -3.,   97.],       [  38.,    1.,  102.],       [  56.,   -5.,   96.],       [  79.,   -1.,   97.],       [  95.,   -3.,   96.],       [  -3.,   29.,  100.],       [  18.,   32.,  101.],       [  42.,   36.,  101.],       [  60.,   38.,   96.],       [  80.,   33.,   95.],       [  96.,   30.,  103.],       [  -4.,   64.,  100.],       [  18.,   66.,   99.],       [  44.,   63.,  100.],       [  59.,   61.,  104.],       [  81.,   62.,   95.],       [  98.,   66.,  103.],       [  -2.,   96.,  101.],       [  22.,  100.,   98.],       [  36.,   99.,  103.],       [  60.,  100.,   98.],       [  76.,   98.,   96.],       [  95.,   98.,   98.]])
        basis0 = BSplineBasis(4, np.array([ 0. ,  0. ,  0. ,  0. ,  1. ,  1.7,  3. ,  3. ,  3. ,  3. ]))
        basis1 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  0.8,  2. ,  2. ,  2. ]))
        basis2 = BSplineBasis(4, np.array([ 0. ,  0. ,  0. ,  0. ,  1.2,  1.9,  2.8,  4. ,  4. ,  4. ,  4. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 3)
        self.assertEqual(vol2.order(direction=1), 3)
        self.assertEqual(vol2.order(direction=2), 4)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 4)
        self.assertEqual(vol3.order(direction=1), 2)
        self.assertEqual(vol3.order(direction=2), 4)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 4)
        self.assertEqual(vol4.order(direction=1), 3)
        self.assertEqual(vol4.order(direction=2), 3)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)
    def test_curve_3D_p4_C0_periodic(self):
        controlpoints = np.array([[ 103.,    0.,   -1.],
       [  67.,   74.,    0.],
       [  -2.,  101.,   -4.],
       [ -76.,   72.,    0.],
       [ -97.,   -3.,   -1.],
       [ -75.,  -69.,   -2.],
       [   0.,  -98.,    3.],
       [  65.,  -70.,    1.]])
        basis0 = BSplineBasis(4, np.array([-0.9,  0. ,  0. ,  0. ,  0.6,  2.3,  3.1,  4.1,  5.1,  6. ,  6. ,  6. ,  6.6]),0)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 3)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p43_C0_periodic(self):
        controlpoints = np.array([[  58.,   -2.,    1.],       [  47.,   43.,    3.],       [   1.,   57.,   -4.],       [ -38.,   40.,   -3.],       [ -63.,   -1.,   -3.],       [ -41.,  -41.,    4.],       [  -4.,  -63.,    1.],       [  40.,  -43.,    2.],       [  68.,    1.,   -5.],       [  52.,   56.,    3.],       [  -4.,   73.,    1.],       [ -57.,   52.,   -2.],       [ -73.,   -2.,    3.],       [ -49.,  -57.,    1.],       [  -5.,  -69.,   -4.],       [  48.,  -57.,    0.],       [  82.,   -4.,    2.],       [  59.,   57.,    1.],       [   4.,   88.,   -2.],       [ -63.,   58.,   -4.],       [ -86.,    3.,   -4.],       [ -57.,  -59.,    4.],       [   3.,  -88.,   -1.],       [  63.,  -59.,   -3.],       [ 101.,    2.,   -2.],       [  69.,   74.,   -5.],       [  -1.,  103.,   -3.],       [ -71.,   71.,   -2.],       [-105.,   -5.,    4.],       [ -68.,  -68.,    1.],       [   2.,  -99.,    3.],       [  68.,  -67.,   -2.]])
        basis0 = BSplineBasis(4, np.array([-1.3,  0. ,  0. ,  0. ,  1.2,  1.7,  3.3,  3.9,  4.7,  6. ,  6. ,  6. ,  7.2]),0)
        basis1 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  0.6,  2. ,  2. ,  2. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 3)
        self.assertEqual(surf2.order(direction=1), 3)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 4)
        self.assertEqual(surf3.order(direction=1), 2)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p334_C0_periodic(self):
        controlpoints = np.array([[  54.,    1.,   -1.],       [  25.,   38.,    4.],       [ -22.,   41.,    4.],       [ -46.,    4.,    4.],       [ -28.,  -41.,    1.],       [  21.,  -47.,    0.],       [  59.,   -3.,   -2.],       [  27.,   59.,   -3.],       [ -31.,   57.,    1.],       [ -63.,    3.,    2.],       [ -36.,  -55.,    0.],       [  35.,  -51.,    0.],       [  71.,    2.,    2.],       [  41.,   63.,   -3.],       [ -40.,   62.,    1.],       [ -72.,    3.,    4.],       [ -36.,  -65.,    0.],       [  34.,  -69.,   -5.],       [  83.,    1.,   -2.],       [  41.,   78.,    3.],       [ -45.,   79.,   -1.],       [ -86.,    2.,   -2.],       [ -42.,  -81.,   -1.],       [  43.,  -76.,    1.],       [ 102.,   -2.,   -1.],       [  48.,   90.,    4.],       [ -51.,   90.,   -5.],       [-104.,   -4.,   -2.],       [ -51.,  -85.,    3.],       [  45.,  -86.,   -5.],       [  47.,   -3.,   12.],       [  28.,   45.,   14.],       [ -30.,   41.,    7.],       [ -48.,    1.,   13.],       [ -27.,  -44.,    5.],       [  21.,  -43.,   11.],       [  60.,    1.,    6.],       [  29.,   58.,    8.],       [ -35.,   50.,    6.],       [ -63.,    4.,    8.],       [ -34.,  -53.,   13.],       [  31.,  -54.,   12.],       [  79.,    2.,    8.],       [  37.,   63.,    7.],       [ -40.,   69.,   14.],       [ -73.,    4.,    8.],       [ -35.,  -67.,   14.],       [  42.,  -62.,   14.],       [  92.,   -4.,    5.],       [  41.,   78.,   11.],       [ -48.,   76.,   10.],       [ -92.,   -3.,    6.],       [ -46.,  -73.,    7.],       [  41.,  -78.,   14.],       [  96.,   -4.,   11.],       [  49.,   82.,   12.],       [ -52.,   85.,    7.],       [-102.,    3.,    7.],       [ -51.,  -90.,    9.],       [  48.,  -91.,    6.],       [  54.,   -4.,   15.],       [  23.,   47.,   16.],       [ -22.,   45.,   22.],       [ -54.,    2.,   21.],       [ -22.,  -46.,   20.],       [  22.,  -44.,   17.],       [  66.,   -4.,   21.],       [  29.,   50.,   17.],       [ -36.,   51.,   21.],       [ -66.,   -3.,   15.],       [ -32.,  -58.,   23.],       [  27.,  -54.,   23.],       [  72.,    1.,   23.],       [  40.,   69.,   23.],       [ -39.,   66.,   22.],       [ -77.,    0.,   19.],       [ -37.,  -60.,   17.],       [  41.,  -69.,   21.],       [  85.,   -1.,   18.],       [  47.,   80.,   19.],       [ -49.,   74.,   19.],       [ -89.,    0.,   22.],       [ -41.,  -77.,   19.],       [  41.,  -78.,   15.],       [ 100.,    1.,   16.],       [  51.,   89.,   21.],       [ -51.,   88.,   19.],       [ -98.,    0.,   15.],       [ -53.,  -85.,   18.],       [  47.,  -87.,   24.],       [  54.,    3.,   33.],       [  29.,   48.,   29.],       [ -24.,   39.,   27.],       [ -47.,   -2.,   32.],       [ -29.,  -41.,   33.],       [  20.,  -47.,   33.],       [  63.,   -5.,   28.],       [  26.,   51.,   34.],       [ -32.,   52.,   25.],       [ -64.,    0.,   32.],       [ -33.,  -52.,   25.],       [  34.,  -60.,   30.],       [  78.,   -5.,   31.],       [  33.,   64.,   34.],       [ -39.,   61.,   33.],       [ -72.,    4.,   29.],       [ -41.,  -70.,   26.],       [  42.,  -63.,   31.],       [  87.,    4.,   26.],       [  45.,   78.,   32.],       [ -49.,   76.,   31.],       [ -89.,    3.,   27.],       [ -46.,  -79.,   25.],       [  47.,  -72.,   34.],       [  97.,   -5.,   32.],       [  47.,   84.,   30.],       [ -51.,   90.,   25.],       [ -98.,   -2.,   34.],       [ -46.,  -86.,   27.],       [  52.,  -85.,   31.],       [  54.,   -2.,   35.],       [  29.,   42.,   43.],       [ -27.,   47.,   42.],       [ -48.,   -5.,   37.],       [ -29.,  -44.,   41.],       [  27.,  -40.,   43.],       [  60.,    4.,   36.],       [  35.,   56.,   40.],       [ -28.,   56.,   41.],       [ -63.,    4.,   41.],       [ -32.,  -53.,   40.],       [  31.,  -58.,   35.],       [  78.,    0.,   37.],       [  37.,   67.,   35.],       [ -39.,   67.,   39.],       [ -77.,    2.,   41.],       [ -33.,  -65.,   37.],       [  39.,  -67.,   39.],       [  89.,   -4.,   37.],       [  48.,   77.,   41.],       [ -42.,   77.,   43.],       [ -91.,    2.,   42.],       [ -41.,  -74.,   39.],       [  47.,  -74.,   36.],       [  95.,    0.,   42.],       [  50.,   88.,   41.],       [ -47.,   90.,   38.],       [-102.,   -5.,   38.],       [ -50.,  -82.,   37.],       [  46.,  -87.,   38.]])
        basis0 = BSplineBasis(3, np.array([-0.9,  0. ,  0. ,  1.4,  2. ,  3.4,  4.1,  5. ,  5. ,  6.4]),0)
        basis1 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  1.1,  2.1,  3. ,  3. ,  3. ]))
        basis2 = BSplineBasis(4, np.array([ 0. ,  0. ,  0. ,  0. ,  0.9,  2. ,  2. ,  2. ,  2. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 2)
        self.assertEqual(vol2.order(direction=1), 3)
        self.assertEqual(vol2.order(direction=2), 4)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 3)
        self.assertEqual(vol3.order(direction=1), 2)
        self.assertEqual(vol3.order(direction=2), 4)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 3)
        self.assertEqual(vol4.order(direction=1), 3)
        self.assertEqual(vol4.order(direction=2), 3)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)
    def test_curve_3D_p4_C1_periodic(self):
        controlpoints = np.array([[  96.,    2.,    0.],
       [  80.,   61.,    1.],
       [  28.,   93.,    1.],
       [ -29.,   96.,    0.],
       [ -86.,   53.,   -3.],
       [-105.,    1.,   -2.],
       [ -83.,  -63.,    3.],
       [ -28.,  -95.,   -5.],
       [  32.,  -95.,   -4.],
       [  82.,  -62.,    0.]])
        basis0 = BSplineBasis(4, np.array([ -2.2,  -1. ,   0. ,   0. ,   0.9,   2.2,   3.1,   4.3,   4.7,   5.8,   6.8,   8. ,   9. ,
         9. ,   9.9,  11.2]),1)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 3)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p44_C1_periodic(self):
        controlpoints = np.array([[  61.,    4.,   -3.],       [  46.,   37.,    0.],       [  19.,   55.,    3.],       [ -17.,   58.,   -4.],       [ -48.,   36.,    4.],       [ -64.,   -3.,    1.],       [ -48.,  -33.,   -4.],       [ -18.,  -62.,    4.],       [  15.,  -56.,   -5.],       [  52.,  -32.,    1.],       [  66.,    3.,    3.],       [  50.,   43.,   -3.],       [  15.,   63.,    1.],       [ -22.,   63.,   -2.],       [ -50.,   37.,   -3.],       [ -68.,    1.,    0.],       [ -59.,  -42.,   -2.],       [ -17.,  -61.,    0.],       [  16.,  -67.,   -5.],       [  51.,  -39.,    1.],       [  76.,   -3.,   -4.],       [  60.,   38.,   -5.],       [  26.,   69.,   -5.],       [ -28.,   66.,   -3.],       [ -63.,   41.,   -2.],       [ -74.,    1.,    2.],       [ -56.,  -47.,   -5.],       [ -28.,  -67.,   -4.],       [  24.,  -66.,   -2.],       [  55.,  -42.,   -3.],       [  77.,   -2.,   -5.],       [  65.,   51.,   -3.],       [  27.,   71.,    0.],       [ -28.,   73.,    3.],       [ -65.,   44.,    0.],       [ -77.,    1.,   -1.],       [ -68.,  -45.,   -3.],       [ -27.,  -79.,    3.],       [  24.,  -78.,   -4.],       [  68.,  -47.,    1.],       [  90.,   -3.,   -3.],       [  71.,   50.,   -3.],       [  31.,   79.,    4.],       [ -29.,   79.,   -4.],       [ -72.,   47.,    0.],       [ -85.,    3.,    0.],       [ -74.,  -52.,    0.],       [ -29.,  -84.,   -2.],       [  23.,  -87.,   -5.],       [  66.,  -49.,    2.],       [  93.,   -2.,   -1.],       [  71.,   57.,    3.],       [  33.,   93.,    0.],       [ -26.,   85.,    0.],       [ -72.,   59.,    3.],       [ -96.,    0.,    4.],       [ -75.,  -53.,    2.],       [ -32.,  -88.,   -1.],       [  24.,  -93.,    3.],       [  79.,  -56.,   -4.],       [ 101.,    3.,   -3.],       [  84.,   60.,    0.],       [  30.,   92.,   -1.],       [ -27.,   99.,    2.],       [ -81.,   58.,    2.],       [ -97.,   -3.,    1.],       [ -85.,  -55.,   -3.],       [ -33.,  -91.,    1.],       [  26.,  -91.,   -1.],       [  77.,  -60.,    0.]])
        basis0 = BSplineBasis(4, np.array([ -1.7,  -1.2,   0. ,   0. ,   0.7,   1.9,   2.8,   4. ,   5.2,   5.9,   7.3,   7.8,   9. ,
         9. ,   9.7,  10.9]),1)
        basis1 = BSplineBasis(4, np.array([ 0. ,  0. ,  0. ,  0. ,  0.9,  2.2,  3.4,  4. ,  4. ,  4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 3)
        self.assertEqual(surf2.order(direction=1), 4)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 4)
        self.assertEqual(surf3.order(direction=1), 3)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p432_C1_periodic(self):
        controlpoints = np.array([[  50.,    3.,   -1.],       [  41.,   29.,    2.],       [  12.,   49.,    3.],       [ -15.,   44.,   -2.],       [ -42.,   34.,   -3.],       [ -46.,   -4.,   -2.],       [ -39.,  -32.,    4.],       [ -17.,  -53.,    3.],       [  11.,  -52.,    0.],       [  38.,  -34.,   -2.],       [  59.,   -2.,    1.],       [  48.,   30.,   -3.],       [  22.,   61.,   -1.],       [ -22.,   59.,   -4.],       [ -50.,   33.,    4.],       [ -61.,    0.,    2.],       [ -52.,  -35.,   -5.],       [ -18.,  -59.,   -4.],       [  18.,  -61.,   -4.],       [  51.,  -34.,    3.],       [  72.,   -2.,   -2.],       [  58.,   39.,   -1.],       [  26.,   64.,    3.],       [ -19.,   65.,   -1.],       [ -61.,   39.,    0.],       [ -72.,   -4.,    0.],       [ -59.,  -40.,    3.],       [ -27.,  -68.,   -4.],       [  22.,  -65.,    1.],       [  60.,  -46.,   -1.],       [  84.,    4.,    1.],       [  69.,   44.,    4.],       [  25.,   76.,   -3.],       [ -23.,   76.,   -5.],       [ -70.,   51.,   -1.],       [ -77.,   -3.,   -1.],       [ -66.,  -46.,    1.],       [ -27.,  -81.,   -2.],       [  22.,  -75.,    3.],       [  68.,  -44.,   -3.],       [  86.,   -1.,    1.],       [  76.,   50.,    4.],       [  30.,   86.,    2.],       [ -26.,   87.,   -3.],       [ -71.,   57.,   -1.],       [ -88.,    2.,   -1.],       [ -74.,  -57.,    2.],       [ -26.,  -90.,    3.],       [  24.,  -86.,    0.],       [  76.,  -58.,    2.],       [ 100.,   -2.,    2.],       [  83.,   62.,    0.],       [  26.,   93.,    4.],       [ -34.,   94.,   -3.],       [ -86.,   61.,   -5.],       [-104.,   -3.,   -1.],       [ -78.,  -55.,   -5.],       [ -32.,  -99.,   -3.],       [  27.,  -95.,   -1.],       [  81.,  -56.,    0.],       [  52.,    3.,    9.],       [  38.,   34.,    8.],       [  13.,   51.,    9.],       [ -14.,   49.,   12.],       [ -39.,   33.,   13.],       [ -51.,    1.,   14.],       [ -39.,  -32.,   14.],       [ -21.,  -52.,    8.],       [  12.,  -51.,    6.],       [  44.,  -25.,    9.],       [  61.,    4.,    9.],       [  52.,   35.,    6.],       [  14.,   56.,   11.],       [ -16.,   55.,   11.],       [ -46.,   35.,    7.],       [ -64.,    4.,   12.],       [ -47.,  -33.,    7.],       [ -22.,  -61.,    8.],       [  14.,  -59.,    9.],       [  48.,  -36.,    5.],       [  69.,    2.,    5.],       [  55.,   39.,   13.],       [  18.,   64.,    6.],       [ -23.,   62.,   12.],       [ -56.,   37.,   14.],       [ -70.,   -1.,   13.],       [ -60.,  -44.,   10.],       [ -24.,  -70.,   10.],       [  25.,  -71.,   10.],       [  54.,  -44.,   12.],       [  80.,   -1.,   10.],       [  61.,   43.,   14.],       [  26.,   78.,   12.],       [ -30.,   79.,    7.],       [ -66.,   49.,   14.],       [ -82.,    4.,   11.],       [ -63.,  -43.,    8.],       [ -27.,  -76.,   10.],       [  26.,  -74.,    6.],       [  66.,  -43.,    9.],       [  85.,    4.,   13.],       [  70.,   57.,    9.],       [  25.,   84.,   14.],       [ -25.,   86.,    6.],       [ -71.,   50.,    9.],       [ -91.,    1.,   13.],       [ -70.,  -55.,   12.],       [ -30.,  -89.,   11.],       [  28.,  -91.,   12.],       [  72.,  -56.,    6.],       [ 102.,   -3.,    6.],       [  83.,   61.,   10.],       [  29.,   94.,   10.],       [ -29.,   97.,    5.],       [ -84.,   54.,   12.],       [ -97.,    4.,   12.],       [ -79.,  -56.,   12.],       [ -30.,  -98.,    7.],       [  29.,  -98.,   12.],       [  78.,  -60.,    6.],       [  48.,    1.,   18.],       [  40.,   24.,   18.],       [  11.,   52.,   19.],       [ -12.,   50.,   19.],       [ -43.,   34.,   19.],       [ -53.,   -3.,   24.],       [ -45.,  -34.,   23.],       [ -17.,  -49.,   18.],       [  17.,  -47.,   21.],       [  42.,  -26.,   17.],       [  62.,   -2.,   16.],       [  53.,   32.,   22.],       [  19.,   61.,   24.],       [ -21.,   52.,   20.],       [ -53.,   39.,   18.],       [ -59.,   -3.,   16.],       [ -49.,  -35.,   15.],       [ -17.,  -61.,   24.],       [  20.,  -62.,   16.],       [  44.,  -32.,   18.],       [  65.,    4.,   18.],       [  58.,   41.,   20.],       [  19.,   71.,   16.],       [ -23.,   64.,   17.],       [ -56.,   37.,   16.],       [ -71.,    3.,   23.],       [ -58.,  -42.,   22.],       [ -26.,  -69.,   15.],       [  19.,  -67.,   15.],       [  54.,  -37.,   17.],       [  84.,    1.,   22.],       [  66.,   42.,   15.],       [  23.,   76.,   20.],       [ -20.,   76.,   21.],       [ -69.,   50.,   15.],       [ -85.,    0.,   23.],       [ -69.,  -47.,   18.],       [ -22.,  -77.,   21.],       [  20.,  -81.,   24.],       [  69.,  -46.,   22.],       [  87.,   -4.,   21.],       [  77.,   53.,   19.],       [  24.,   89.,   22.],       [ -26.,   80.,   24.],       [ -77.,   54.,   22.],       [ -92.,    3.,   22.],       [ -72.,  -57.,   21.],       [ -30.,  -85.,   23.],       [  30.,  -83.,   18.],       [  68.,  -54.,   15.],       [ 102.,   -4.,   19.],       [  85.,   60.,   17.],       [  29.,   93.,   20.],       [ -33.,   95.,   18.],       [ -79.,   62.,   19.],       [-104.,    1.,   19.],       [ -86.,  -58.,   19.],       [ -28.,  -91.,   22.],       [  31.,  -99.,   15.],       [  83.,  -55.,   23.],       [  50.,   -4.,   30.],       [  42.,   28.,   32.],       [  20.,   52.,   26.],       [ -18.,   52.,   31.],       [ -42.,   26.,   33.],       [ -46.,   -4.,   34.],       [ -40.,  -31.,   27.],       [ -14.,  -50.,   32.],       [  15.,  -53.,   31.],       [  45.,  -28.,   29.],       [  64.,    3.,   29.],       [  46.,   31.,   28.],       [  17.,   59.,   25.],       [ -14.,   56.,   31.],       [ -47.,   32.,   34.],       [ -57.,   -2.,   31.],       [ -45.,  -37.,   31.],       [ -20.,  -53.,   31.],       [  18.,  -53.,   32.],       [  53.,  -35.,   29.],       [  70.,    2.,   26.],       [  60.,   36.,   32.],       [  17.,   65.,   31.],       [ -25.,   69.,   28.],       [ -56.,   43.,   26.],       [ -74.,   -3.,   33.],       [ -60.,  -37.,   27.],       [ -20.,  -66.,   26.],       [  17.,  -69.,   27.],       [  55.,  -41.,   34.],       [  78.,    4.,   27.],       [  63.,   42.,   34.],       [  26.,   75.,   29.],       [ -23.,   71.,   31.],       [ -60.,   44.,   26.],       [ -84.,    2.,   29.],       [ -62.,  -50.,   32.],       [ -25.,  -73.,   34.],       [  21.,  -81.,   27.],       [  60.,  -48.,   33.],       [  92.,    1.,   29.],       [  69.,   55.,   26.],       [  25.,   85.,   28.],       [ -29.,   81.,   33.],       [ -73.,   50.,   28.],       [ -88.,    0.,   25.],       [ -76.,  -52.,   30.],       [ -25.,  -86.,   25.],       [  23.,  -83.,   27.],       [  70.,  -53.,   27.],       [ 100.,    3.,   32.],       [  78.,   56.,   28.],       [  26.,   92.,   29.],       [ -32.,   99.,   25.],       [ -81.,   59.,   27.],       [ -98.,    2.,   27.],       [ -80.,  -56.,   33.],       [ -33.,  -92.,   29.],       [  26.,  -98.,   31.],       [  78.,  -64.,   33.],       [  53.,   -2.,   42.],       [  39.,   24.,   37.],       [  16.,   50.,   40.],       [ -16.,   48.,   40.],       [ -37.,   25.,   40.],       [ -53.,    0.,   35.],       [ -37.,  -34.,   35.],       [ -17.,  -46.,   44.],       [  14.,  -50.,   35.],       [  37.,  -27.,   39.],       [  64.,    4.,   39.],       [  46.,   36.,   44.],       [  14.,   60.,   43.],       [ -17.,   58.,   43.],       [ -48.,   31.,   42.],       [ -62.,   -4.,   44.],       [ -46.,  -40.,   37.],       [ -20.,  -58.,   37.],       [  14.,  -61.,   44.],       [  49.,  -37.,   35.],       [  65.,    1.,   36.],       [  52.,   40.,   41.],       [  18.,   63.,   35.],       [ -20.,   67.,   39.],       [ -61.,   40.,   42.],       [ -73.,    2.,   41.],       [ -60.,  -39.,   41.],       [ -20.,  -67.,   35.],       [  23.,  -72.,   36.],       [  53.,  -45.,   43.],       [  81.,   -4.,   36.],       [  64.,   49.,   38.],       [  23.,   73.,   41.],       [ -22.,   76.,   36.],       [ -63.,   49.,   38.],       [ -82.,   -2.,   38.],       [ -61.,  -43.,   43.],       [ -30.,  -72.,   44.],       [  28.,  -74.,   36.],       [  62.,  -45.,   36.],       [  85.,    0.,   36.],       [  71.,   50.,   40.],       [  27.,   90.,   44.],       [ -27.,   88.,   36.],       [ -72.,   49.,   38.],       [ -87.,    2.,   42.],       [ -69.,  -53.,   41.],       [ -31.,  -91.,   36.],       [  28.,  -83.,   42.],       [  69.,  -50.,   36.],       [  99.,    0.,   36.],       [  77.,   63.,   42.],       [  32.,   93.,   44.],       [ -28.,   93.,   41.],       [ -86.,   61.,   36.],       [ -97.,    3.,   36.],       [ -78.,  -56.,   44.],       [ -34.,  -99.,   36.],       [  29.,  -97.,   37.],       [  82.,  -63.,   39.]])
        basis0 = BSplineBasis(4, np.array([ -2.4,  -0.8,   0. ,   0. ,   0.6,   1.9,   3.2,   3.8,   5.1,   6.3,   6.6,   8.2,   9. ,
         9. ,   9.6,  10.9]),1)
        basis1 = BSplineBasis(3, np.array([ 0. ,  0. ,  0. ,  1.4,  2.2,  2.6,  4. ,  4. ,  4. ]))
        basis2 = BSplineBasis(2, np.array([ 0. ,  0. ,  0.9,  2. ,  3.4,  4. ,  4. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 3)
        self.assertEqual(vol2.order(direction=1), 3)
        self.assertEqual(vol2.order(direction=2), 2)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 4)
        self.assertEqual(vol3.order(direction=1), 2)
        self.assertEqual(vol3.order(direction=2), 2)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 4)
        self.assertEqual(vol4.order(direction=1), 3)
        self.assertEqual(vol4.order(direction=2), 1)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)
    def test_curve_2D_p5(self):
        controlpoints = np.array([[  -3.,   -2.],
       [  10.,    2.],
       [  25.,    4.],
       [  39.,    1.],
       [  54.,   -1.],
       [  74.,   -2.],
       [  81.,   -2.],
       [ 102.,    3.]])
        basis0 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0.8,  2.3,  2.8,  4. ,  4. ,  4. ,  4. ,  4. ]))
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 4)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_2D_p56(self):
        controlpoints = np.array([[   1.,    3.],       [  17.,   -3.],       [  25.,   -5.],       [  40.,   -2.],       [  60.,    2.],       [  72.,    0.],       [  85.,    3.],       [  95.,    2.],       [  -5.,   13.],       [  11.,    8.],       [  25.,   14.],       [  47.,   14.],       [  58.,    7.],       [  66.,    7.],       [  86.,   13.],       [ 104.,   10.],       [  -1.,   25.],       [  11.,   26.],       [  32.,   21.],       [  40.,   21.],       [  56.,   24.],       [  67.,   29.],       [  85.,   23.],       [  98.,   24.],       [  -2.,   38.],       [  17.,   39.],       [  24.,   34.],       [  41.,   38.],       [  61.,   33.],       [  67.,   39.],       [  87.,   36.],       [ 104.,   37.],       [  -4.,   54.],       [  16.,   49.],       [  27.,   54.],       [  41.,   51.],       [  57.,   46.],       [  67.,   52.],       [  82.,   52.],       [  99.,   48.],       [  -2.,   61.],       [  14.,   67.],       [  32.,   61.],       [  40.,   60.],       [  59.,   63.],       [  72.,   63.],       [  84.,   66.],       [ 103.,   67.],       [   4.,   76.],       [   9.,   78.],       [  32.,   76.],       [  41.,   79.],       [  59.,   79.],       [  76.,   78.],       [  81.,   78.],       [ 101.,   74.],       [  -4.,   88.],       [  15.,   84.],       [  27.,   83.],       [  43.,   84.],       [  57.,   91.],       [  71.,   87.],       [  83.,   91.],       [ 104.,   85.],       [  -3.,  103.],       [  12.,   96.],       [  25.,   98.],       [  42.,  104.],       [  60.,   97.],       [  70.,   98.],       [  89.,  101.],       [ 100.,   99.]])
        basis0 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0.7,  2.3,  2.9,  4. ,  4. ,  4. ,  4. ,  4. ]))
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.8,  2.1,  3.2,  4. ,  4. ,  4. ,  4. ,  4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 4)
        self.assertEqual(surf2.order(direction=1), 6)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 5)
        self.assertEqual(surf3.order(direction=1), 5)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_curve_2D_p5_C0_periodic(self):
        controlpoints = np.array([[ 100.,    1.],
       [  76.,   58.],
       [  35.,   95.],
       [ -28.,   99.],
       [ -82.,   57.],
       [-105.,   -1.],
       [ -84.,  -61.],
       [ -36.,  -91.],
       [  34.,  -97.],
       [  80.,  -62.]])
        basis0 = BSplineBasis(5, np.array([-1.4,  0. ,  0. ,  0. ,  0. ,  0.8,  2. ,  3.3,  4.1,  4.6,  5.6,  7. ,  7. ,  7. ,  7. ,
        7.8]),0)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 4)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_2D_p57_C0_periodic(self):
        controlpoints = np.array([[  56.,   -1.],       [  47.,   39.],       [  21.,   59.],       [ -19.,   61.],       [ -52.,   33.],       [ -65.,   -3.],       [ -54.,  -38.],       [ -24.,  -60.],       [  16.,  -53.],       [  53.,  -33.],       [  67.,   -2.],       [  47.,   35.],       [  23.,   57.],       [ -22.,   65.],       [ -56.,   39.],       [ -69.,    0.],       [ -52.,  -36.],       [ -25.,  -60.],       [  22.,  -64.],       [  50.,  -42.],       [  66.,   -4.],       [  51.,   40.],       [  20.,   62.],       [ -26.,   62.],       [ -56.,   37.],       [ -67.,    4.],       [ -53.,  -39.],       [ -26.,  -71.],       [  22.,  -63.],       [  57.,  -46.],       [  76.,   -4.],       [  60.,   39.],       [  26.,   68.],       [ -21.,   72.],       [ -59.,   42.],       [ -77.,    4.],       [ -56.,  -48.],       [ -23.,  -74.],       [  18.,  -70.],       [  64.,  -48.],       [  79.,   -1.],       [  66.,   40.],       [  21.,   70.],       [ -22.,   78.],       [ -60.,   50.],       [ -76.,    1.],       [ -67.,  -41.],       [ -24.,  -72.],       [  23.,  -73.],       [  67.,  -48.],       [  83.,    1.],       [  70.,   52.],       [  29.,   82.],       [ -31.,   80.],       [ -70.,   48.],       [ -86.,   -3.],       [ -63.,  -50.],       [ -30.,  -84.],       [  27.,  -80.],       [  69.,  -50.],       [  87.,   -1.],       [  73.,   55.],       [  26.,   83.],       [ -32.,   86.],       [ -67.,   53.],       [ -83.,    3.],       [ -74.,  -49.],       [ -22.,  -81.],       [  31.,  -86.],       [  73.,  -52.],       [  95.,    1.],       [  75.,   54.],       [  31.,   86.],       [ -30.,   84.],       [ -75.,   55.],       [ -89.,   -3.],       [ -70.,  -52.],       [ -31.,  -85.],       [  30.,  -87.],       [  72.,  -49.],       [  90.,   -1.],       [  77.,   55.],       [  29.,   93.],       [ -35.,   85.],       [ -76.,   55.],       [ -97.,    3.],       [ -74.,  -56.],       [ -28.,  -93.],       [  30.,  -92.],       [  81.,  -58.],       [ 101.,    4.],       [  83.,   63.],       [  27.,   96.],       [ -34.,   99.],       [ -80.,   54.],       [ -98.,   -5.],       [ -84.,  -55.],       [ -36.,  -91.],       [  33.,  -93.],       [  78.,  -55.]])
        basis0 = BSplineBasis(5, np.array([-0.6,  0. ,  0. ,  0. ,  0. ,  1.2,  2.4,  3.2,  4.2,  4.9,  6.4,  7. ,  7. ,  7. ,  7. ,
        8.2]),0)
        basis1 = BSplineBasis(7, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.6,  1.9,  3.4,  4. ,  4. ,  4. ,  4. ,  4. ,
        4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 4)
        self.assertEqual(surf2.order(direction=1), 7)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 5)
        self.assertEqual(surf3.order(direction=1), 6)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_curve_2D_p5_C1_periodic(self):
        controlpoints = np.array([[  99.,    2.],
       [  91.,   54.],
       [  50.,   82.],
       [  -2.,   96.],
       [ -47.,   89.],
       [ -88.,   53.],
       [-105.,   -5.],
       [ -88.,  -48.],
       [ -55.,  -83.],
       [   3., -102.],
       [  50.,  -85.],
       [  89.,  -52.]])
        basis0 = BSplineBasis(5, np.array([ -2.3,  -0.8,   0. ,   0. ,   0. ,   0.8,   2.1,   3.4,   3.9,   4.7,   5.9,   6.9,   7.7,
         9.2,  10. ,  10. ,  10. ,  10.8,  12.1]),1)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 4)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_2D_p66_C1_periodic(self):
        controlpoints = np.array([[  59.,   -4.],       [  52.,   22.],       [  36.,   48.],       [   8.,   54.],       [ -10.,   59.],       [ -39.,   49.],       [ -54.,   27.],       [ -57.,    1.],       [ -54.,  -26.],       [ -34.,  -47.],       [ -11.,  -59.],       [  17.,  -64.],       [  38.,  -49.],       [  58.,  -22.],       [  64.,    3.],       [  61.,   33.],       [  45.,   54.],       [  19.,   69.],       [ -17.,   62.],       [ -43.,   54.],       [ -56.,   32.],       [ -62.,   -4.],       [ -64.,  -26.],       [ -46.,  -56.],       [ -14.,  -70.],       [  19.,  -62.],       [  40.,  -49.],       [  63.,  -27.],       [  69.,   -1.],       [  68.,   29.],       [  50.,   55.],       [  18.,   67.],       [ -18.,   68.],       [ -50.,   59.],       [ -67.,   27.],       [ -74.,   -2.],       [ -68.,  -36.],       [ -48.,  -58.],       [ -12.,  -76.],       [  13.,  -69.],       [  42.,  -57.],       [  70.,  -35.],       [  76.,    4.],       [  68.,   38.],       [  54.,   58.],       [  16.,   81.],       [ -20.,   78.],       [ -46.,   61.],       [ -77.,   32.],       [ -77.,    0.],       [ -70.,  -30.],       [ -52.,  -65.],       [ -18.,  -75.],       [  19.,  -78.],       [  46.,  -59.],       [  73.,  -32.],       [  83.,    4.],       [  77.,   40.],       [  56.,   69.],       [  24.,   81.],       [ -22.,   80.],       [ -50.,   71.],       [ -77.,   39.],       [ -84.,   -5.],       [ -77.,  -34.],       [ -51.,  -71.],       [ -23.,  -89.],       [  23.,  -89.],       [  52.,  -69.],       [  74.,  -35.],       [  97.,   -2.],       [  83.,   43.],       [  61.,   74.],       [  21.,   91.],       [ -24.,   90.],       [ -54.,   77.],       [ -88.,   42.],       [ -89.,    0.],       [ -80.,  -44.],       [ -55.,  -71.],       [ -21.,  -93.],       [  25.,  -87.],       [  57.,  -77.],       [  86.,  -42.],       [ 104.,    4.],       [  88.,   39.],       [  61.,   74.],       [  21.,   95.],       [ -27.,  100.],       [ -59.,   73.],       [ -88.,   44.],       [ -96.,    3.],       [ -94.,  -40.],       [ -63.,  -81.],       [ -20., -102.],       [  22.,  -94.],       [  66.,  -81.],       [  94.,  -45.]])
        basis0 = BSplineBasis(6, np.array([ -1.9,  -0.7,   0. ,   0. ,   0. ,   0. ,   1.1,   1.8,   3. ,   4.2,   4.9,   6. ,   7.4,
         8.1,   9.1,  10.3,  11. ,  11. ,  11. ,  11. ,  12.1,  12.8]),1)
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.3,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 5)
        self.assertEqual(surf2.order(direction=1), 6)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 6)
        self.assertEqual(surf3.order(direction=1), 5)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_curve_2D_p7_C2_periodic(self):
        controlpoints = np.array([[  97.,    2.],
       [  94.,   37.],
       [  79.,   60.],
       [  50.,   91.],
       [  22.,  101.],
       [ -18.,  103.],
       [ -53.,   84.],
       [ -81.,   68.],
       [ -93.,   33.],
       [-100.,    1.],
       [ -92.,  -33.],
       [ -76.,  -68.],
       [ -55.,  -89.],
       [ -23.,  -98.],
       [  19.,  -98.],
       [  48.,  -89.],
       [  81.,  -65.],
       [  93.,  -38.]])
        basis0 = BSplineBasis(7, np.array([ -3.3,  -1.9,  -0.8,   0. ,   0. ,   0. ,   0. ,   1. ,   2.3,   2.7,   4.1,   5.2,   6.4,
         6.7,   7.9,   9. ,  10. ,  11.4,  11.7,  13.1,  14.2,  15. ,  15. ,  15. ,  15. ,  16. ,
        17.3,  17.7]),2)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 6)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_2D_p56_C2_periodic(self):
        controlpoints = np.array([[  63.,   -5.],       [  53.,   26.],       [  38.,   42.],       [   9.,   58.],       [ -13.,   61.],       [ -36.,   47.],       [ -53.,   25.],       [ -62.,    1.],       [ -53.,  -29.],       [ -38.,  -51.],       [ -11.,  -62.],       [  16.,  -57.],       [  36.,  -43.],       [  56.,  -25.],       [  62.,    0.],       [  57.,   26.],       [  39.,   53.],       [  14.,   62.],       [ -13.,   65.],       [ -36.,   54.],       [ -57.,   26.],       [ -65.,    3.],       [ -62.,  -26.],       [ -43.,  -53.],       [ -15.,  -63.],       [  15.,  -66.],       [  38.,  -56.],       [  56.,  -30.],       [  70.,    3.],       [  66.,   30.],       [  47.,   50.],       [  11.,   64.],       [ -19.,   65.],       [ -49.,   52.],       [ -60.,   31.],       [ -68.,   -3.],       [ -66.,  -36.],       [ -49.,  -53.],       [ -14.,  -74.],       [  15.,  -73.],       [  41.,  -52.],       [  63.,  -26.],       [  71.,    4.],       [  70.,   28.],       [  43.,   62.],       [  15.,   71.],       [ -22.,   73.],       [ -47.,   57.],       [ -65.,   28.],       [ -74.,   -4.],       [ -69.,  -31.],       [ -47.,  -63.],       [ -19.,  -73.],       [  13.,  -69.],       [  49.,  -56.],       [  62.,  -31.],       [  83.,    2.],       [  67.,   35.],       [  50.,   64.],       [  18.,   76.],       [ -15.,   81.],       [ -48.,   62.],       [ -77.,   35.],       [ -85.,   -5.],       [ -75.,  -36.],       [ -49.,  -68.],       [ -21.,  -77.],       [  20.,  -78.],       [  48.,  -67.],       [  71.,  -39.],       [  86.,    2.],       [  80.,   37.],       [  53.,   63.],       [  15.,   80.],       [ -18.,   85.],       [ -50.,   62.],       [ -72.,   35.],       [ -89.,    1.],       [ -73.,  -36.],       [ -57.,  -62.],       [ -22.,  -78.],       [  20.,  -78.],       [  49.,  -67.],       [  80.,  -40.],       [  85.,   -4.],       [  77.,   43.],       [  55.,   65.],       [  17.,   85.],       [ -20.,   86.],       [ -54.,   69.],       [ -79.,   42.],       [ -89.,   -2.],       [ -77.,  -44.],       [ -57.,  -76.],       [ -17.,  -83.],       [  15.,  -87.],       [  60.,  -67.],       [  82.,  -41.],       [  90.,    3.],       [  87.,   39.],       [  56.,   74.],       [  18.,   89.],       [ -25.,   91.],       [ -57.,   76.],       [ -84.,   44.],       [ -91.,    2.],       [ -86.,  -39.],       [ -55.,  -72.],       [ -26.,  -95.],       [  19.,  -88.],       [  61.,  -77.],       [  83.,  -45.],       [  96.,   -5.],       [  90.,   44.],       [  65.,   75.],       [  17.,   93.],       [ -20.,   99.],       [ -60.,   79.],       [ -89.,   48.],       [ -98.,   -1.],       [ -95.,  -43.],       [ -63.,  -81.],       [ -21.,  -99.],       [  18., -100.],       [  63.,  -76.],       [  87.,  -44.]])
        basis0 = BSplineBasis(5, np.array([ -2.9,  -2. ,  -1.3,   0. ,   0. ,   1.1,   1.8,   3. ,   4.3,   4.7,   5.7,   6.6,   8.4,
         8.8,  10.1,  11. ,  11.7,  13. ,  13. ,  14.1,  14.8,  16. ]),2)
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1. ,  2. ,  3.1,  4. ,  4. ,  4. ,  4. ,  4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 4)
        self.assertEqual(surf2.order(direction=1), 6)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 5)
        self.assertEqual(surf3.order(direction=1), 5)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_curve_3D_p5(self):
        controlpoints = np.array([[ -4.,  -4.,   4.],
       [ 19.,  -5.,  -5.],
       [ 26.,   1.,  -2.],
       [ 45.,  -1.,  -4.],
       [ 60.,  -3.,  -1.],
       [ 75.,   4.,  -5.],
       [ 89.,  -1.,   1.],
       [ 97.,   0.,   1.]])
        basis0 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0.7,  1.9,  2.9,  4. ,  4. ,  4. ,  4. ,  4. ]))
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 4)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p77(self):
        controlpoints = np.array([[   2.,    0.,    4.],       [  10.,    4.,   -1.],       [  23.,   -3.,    0.],       [  30.,    3.,   -5.],       [  40.,   -2.,   -4.],       [  44.,    0.,    1.],       [  57.,   -3.,   -3.],       [  61.,    1.,   -1.],       [  71.,    3.,    2.],       [  80.,   -5.,   -2.],       [  86.,   -4.,   -5.],       [ 101.,    3.,   -2.],       [   4.,   15.,   -3.],       [   8.,    9.,   -5.],       [  17.,   18.,    2.],       [  25.,   13.,    0.],       [  32.,   12.,   -4.],       [  48.,   14.,    0.],       [  53.,   17.,   -5.],       [  67.,   13.,    0.],       [  68.,   14.,    1.],       [  82.,   15.,   -5.],       [  95.,   16.,   -3.],       [  98.,   13.,    1.],       [  -1.,   30.,   -4.],       [  10.,   27.,    2.],       [  16.,   27.,    0.],       [  27.,   25.,    4.],       [  34.,   25.,    2.],       [  41.,   30.,   -4.],       [  56.,   32.,    1.],       [  63.,   32.,    1.],       [  68.,   23.,   -4.],       [  83.,   27.,   -4.],       [  89.,   26.,   -1.],       [  99.,   30.,    2.],       [   4.,   41.,    3.],       [  12.,   43.,   -2.],       [  18.,   42.,   -3.],       [  24.,   45.,    3.],       [  41.,   46.,   -4.],       [  49.,   41.,   -5.],       [  50.,   44.,    1.],       [  63.,   44.,    4.],       [  67.,   40.,   -1.],       [  77.,   43.,   -3.],       [  91.,   38.,    3.],       [ 104.,   41.,    2.],       [   0.,   57.,    0.],       [   8.,   60.,   -5.],       [  16.,   54.,    4.],       [  23.,   53.,    4.],       [  39.,   61.,    0.],       [  49.,   55.,   -5.],       [  54.,   61.,    0.],       [  60.,   60.,   -5.],       [  75.,   52.,   -1.],       [  78.,   58.,   -1.],       [  93.,   58.,    3.],       [  99.,   57.,   -2.],       [   4.,   75.,   -3.],       [   8.,   74.,    2.],       [  20.,   71.,    0.],       [  30.,   73.,   -4.],       [  40.,   73.,    0.],       [  47.,   72.,    3.],       [  58.,   70.,   -1.],       [  68.,   74.,   -5.],       [  71.,   67.,    1.],       [  86.,   75.,   -4.],       [  87.,   72.,   -3.],       [ 102.,   72.,    2.],       [   2.,   82.,    3.],       [  10.,   83.,    1.],       [  14.,   82.,    2.],       [  29.,   86.,    3.],       [  36.,   85.,   -4.],       [  48.,   81.,    0.],       [  59.,   88.,    0.],       [  68.,   82.,    2.],       [  68.,   80.,    2.],       [  82.,   81.,   -5.],       [  94.,   83.,    2.],       [ 101.,   83.,    3.],       [   1.,  104.,   -2.],       [   8.,  104.,   -4.],       [  17.,  101.,   -4.],       [  22.,  100.,    2.],       [  37.,  103.,   -2.],       [  47.,  104.,    1.],       [  55.,   99.,    1.],       [  63.,   98.,    3.],       [  69.,   99.,    3.],       [  85.,   95.,   -5.],       [  87.,   98.,    1.],       [ 101.,  102.,   -3.]])
        basis0 = BSplineBasis(7, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.9,  2.2,  3.2,  3.7,  4.8,  6. ,  6. ,  6. ,
        6. ,  6. ,  6. ,  6. ]))
        basis1 = BSplineBasis(7, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.3,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 6)
        self.assertEqual(surf2.order(direction=1), 7)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 7)
        self.assertEqual(surf3.order(direction=1), 6)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p565(self):
        controlpoints = np.array([[   3.,   -1.,   -4.],       [  18.,    1.,   -5.],       [  28.,   -2.,    1.],       [  41.,    0.,   -5.],       [  54.,   -1.,    0.],       [  70.,    1.,   -5.],       [  82.,   -1.,    1.],       [  96.,    0.,   -2.],       [  -3.,   20.,   -4.],       [  13.,   17.,    1.],       [  25.,   20.,    3.],       [  46.,   15.,    0.],       [  59.,   17.,    0.],       [  74.,   13.,    4.],       [  84.,   20.,   -1.],       [ 101.,   17.,    1.],       [  -1.,   28.,   -2.],       [  12.,   30.,   -4.],       [  24.,   33.,    0.],       [  46.,   33.,    3.],       [  60.,   30.,    1.],       [  75.,   34.,    4.],       [  87.,   30.,   -3.],       [  96.,   37.,    1.],       [  -3.,   50.,   -3.],       [  16.,   51.,   -5.],       [  28.,   46.,    0.],       [  40.,   45.,    2.],       [  61.,   51.,    2.],       [  73.,   45.,   -5.],       [  86.,   45.,    2.],       [  96.,   53.,   -5.],       [   0.,   70.,    4.],       [  11.,   70.,    0.],       [  23.,   62.,   -2.],       [  41.,   63.,   -5.],       [  53.,   63.,    2.],       [  70.,   65.,    1.],       [  89.,   64.,    3.],       [ 103.,   70.,   -4.],       [  -5.,   85.,    0.],       [  11.,   87.,    2.],       [  29.,   81.,    1.],       [  40.,   82.,    2.],       [  60.,   81.,    4.],       [  68.,   81.,    2.],       [  88.,   86.,    3.],       [  97.,   79.,   -4.],       [  -3.,   98.,   -2.],       [   9.,   98.,   -5.],       [  25.,  104.,    2.],       [  38.,   97.,    2.],       [  60.,   97.,    1.],       [  68.,  101.,    3.],       [  81.,   99.,    4.],       [ 101.,  100.,    4.],       [   2.,   -5.,   17.],       [  12.,   -1.,   21.],       [  24.,   -2.,   24.],       [  47.,    0.,   22.],       [  56.,   -1.,   23.],       [  70.,    4.,   15.],       [  87.,   -2.,   18.],       [ 100.,   -1.,   23.],       [  -5.,   12.,   17.],       [  17.,   13.,   21.],       [  31.,   16.,   17.],       [  47.,   20.,   17.],       [  58.,   21.,   21.],       [  74.,   19.,   24.],       [  86.,   15.,   22.],       [ 104.,   15.,   19.],       [  -2.,   37.,   18.],       [  17.,   34.,   15.],       [  24.,   32.,   24.],       [  42.,   37.,   16.],       [  53.,   28.,   23.],       [  72.,   33.,   20.],       [  89.,   37.,   19.],       [ 100.,   31.,   17.],       [   3.,   47.,   20.],       [  18.,   46.,   24.],       [  26.,   50.,   22.],       [  43.,   51.,   21.],       [  58.,   47.,   21.],       [  66.,   47.,   23.],       [  87.,   53.,   16.],       [  98.,   50.,   24.],       [  -4.,   71.,   17.],       [  13.,   63.,   21.],       [  26.,   68.,   21.],       [  39.,   71.,   22.],       [  53.,   63.,   21.],       [  69.,   67.,   22.],       [  88.,   61.,   16.],       [  95.,   64.,   20.],       [  -3.,   83.,   15.],       [  10.,   79.,   23.],       [  26.,   85.,   19.],       [  44.,   83.,   17.],       [  55.,   85.,   17.],       [  70.,   85.,   21.],       [  90.,   87.,   20.],       [  97.,   84.,   15.],       [  -1.,   96.,   16.],       [  11.,   97.,   21.],       [  25.,   98.,   22.],       [  41.,   95.,   21.],       [  58.,  104.,   22.],       [  74.,   96.,   15.],       [  86.,  101.,   18.],       [  95.,  102.,   20.],       [   4.,    4.,   42.],       [  11.,   -4.,   44.],       [  31.,    0.,   42.],       [  45.,   -4.,   36.],       [  58.,   -2.,   44.],       [  69.,   -4.,   35.],       [  83.,   -2.,   39.],       [  96.,    1.,   42.],       [   4.,   13.,   37.],       [  16.,   16.,   42.],       [  28.,   20.,   38.],       [  37.,   18.,   40.],       [  57.,   18.,   41.],       [  71.,   20.,   43.],       [  87.,   20.,   40.],       [  97.,   20.,   40.],       [   4.,   30.,   41.],       [  14.,   35.,   44.],       [  28.,   38.,   39.],       [  42.,   37.,   43.],       [  57.,   29.,   39.],       [  71.,   38.,   41.],       [  86.,   38.,   36.],       [  99.,   32.,   41.],       [   0.,   53.,   39.],       [  13.,   53.,   44.],       [  29.,   49.,   42.],       [  46.,   46.,   40.],       [  58.,   50.,   43.],       [  69.,   47.,   43.],       [  81.,   54.,   41.],       [ 101.,   45.,   38.],       [   0.,   66.,   36.],       [  14.,   70.,   44.],       [  32.,   69.,   38.],       [  39.,   67.,   44.],       [  53.,   63.,   35.],       [  74.,   65.,   37.],       [  89.,   62.,   44.],       [ 103.,   62.,   35.],       [   3.,   86.,   44.],       [  10.,   81.,   36.],       [  32.,   79.,   44.],       [  43.,   84.,   36.],       [  56.,   80.,   44.],       [  74.,   85.,   41.],       [  83.,   85.,   40.],       [  97.,   87.,   41.],       [   1.,  100.,   43.],       [  17.,   99.,   39.],       [  26.,  100.,   41.],       [  41.,   97.,   41.],       [  55.,   96.,   44.],       [  70.,  103.,   35.],       [  89.,  102.,   42.],       [  98.,  104.,   44.],       [  -4.,    2.,   61.],       [  13.,   -3.,   60.],       [  27.,    0.,   57.],       [  40.,   -4.,   61.],       [  54.,   -5.,   57.],       [  72.,    4.,   58.],       [  83.,   -1.,   61.],       [  96.,    1.,   59.],       [  -4.,   17.,   57.],       [  18.,   20.,   63.],       [  24.,   17.,   58.],       [  45.,   12.,   60.],       [  60.,   14.,   57.],       [  70.,   21.,   60.],       [  83.,   19.,   63.],       [ 100.,   13.,   60.],       [  -4.,   34.,   64.],       [  10.,   30.,   62.],       [  26.,   29.,   60.],       [  40.,   33.,   61.],       [  60.,   32.,   62.],       [  69.,   28.,   55.],       [  83.,   29.,   55.],       [  96.,   32.,   63.],       [   0.,   46.,   59.],       [  15.,   51.,   58.],       [  30.,   53.,   58.],       [  41.,   46.,   64.],       [  58.,   46.,   55.],       [  74.,   48.,   64.],       [  88.,   47.,   64.],       [ 101.,   48.,   55.],       [   1.,   71.,   61.],       [  12.,   70.,   58.],       [  27.,   69.,   57.],       [  47.,   65.,   60.],       [  53.,   69.,   64.],       [  67.,   64.,   59.],       [  85.,   64.,   62.],       [  99.,   62.,   57.],       [  -4.,   82.,   59.],       [  17.,   86.,   64.],       [  29.,   85.,   64.],       [  38.,   80.,   63.],       [  60.,   82.,   56.],       [  74.,   87.,   60.],       [  84.,   82.,   56.],       [ 102.,   82.,   64.],       [   4.,  103.,   64.],       [  14.,  101.,   59.],       [  29.,   98.,   62.],       [  46.,   96.,   63.],       [  52.,   95.,   55.],       [  74.,   97.,   64.],       [  89.,  102.,   61.],       [  97.,   96.,   59.],       [   2.,   -3.,   84.],       [  14.,    3.,   80.],       [  30.,    3.,   75.],       [  46.,   -1.,   82.],       [  61.,   -4.,   81.],       [  69.,    3.,   82.],       [  82.,   -1.,   83.],       [  95.,    4.,   77.],       [   2.,   19.,   75.],       [  14.,   16.,   82.],       [  31.,   15.,   82.],       [  46.,   18.,   82.],       [  61.,   21.,   84.],       [  67.,   14.,   84.],       [  81.,   17.,   79.],       [  97.,   16.,   78.],       [  -3.,   34.,   83.],       [  11.,   35.,   84.],       [  29.,   35.,   77.],       [  46.,   38.,   76.],       [  61.,   34.,   75.],       [  70.,   32.,   80.],       [  90.,   28.,   75.],       [ 101.,   37.,   79.],       [  -5.,   52.,   82.],       [  19.,   50.,   79.],       [  32.,   49.,   77.],       [  43.,   45.,   80.],       [  58.,   50.,   82.],       [  74.,   46.,   77.],       [  87.,   48.,   84.],       [  96.,   46.,   79.],       [   2.,   63.,   76.],       [  14.,   64.,   84.],       [  33.,   67.,   77.],       [  38.,   65.,   84.],       [  60.,   70.,   80.],       [  74.,   63.,   75.],       [  86.,   66.,   76.],       [  95.,   67.,   79.],       [   4.,   80.,   79.],       [  18.,   87.,   75.],       [  29.,   79.,   82.],       [  45.,   86.,   75.],       [  58.,   88.,   83.],       [  69.,   82.,   83.],       [  86.,   78.,   76.],       [  99.,   83.,   78.],       [  -5.,  103.,   80.],       [  12.,   98.,   75.],       [  26.,   99.,   79.],       [  39.,   96.,   75.],       [  55.,  103.,   81.],       [  67.,   98.,   77.],       [  88.,   99.,   81.],       [ 103.,  101.,   79.],       [  -2.,   -5.,  103.],       [  18.,   -4.,  100.],       [  31.,    1.,   95.],       [  45.,    3.,   99.],       [  61.,   -5.,  104.],       [  67.,   -2.,  104.],       [  90.,   -1.,   98.],       [ 104.,    1.,   97.],       [   0.,   13.,  101.],       [  18.,   16.,   98.],       [  30.,   14.,   95.],       [  38.,   19.,   95.],       [  60.,   18.,  104.],       [  67.,   14.,  103.],       [  87.,   19.,  100.],       [  99.,   13.,  102.],       [  -5.,   37.,   98.],       [  15.,   30.,   97.],       [  33.,   36.,  101.],       [  40.,   29.,  104.],       [  61.,   37.,  104.],       [  69.,   28.,   96.],       [  87.,   34.,   99.],       [ 100.,   31.,  100.],       [  -5.,   52.,  102.],       [  15.,   46.,  100.],       [  24.,   45.,   98.],       [  47.,   50.,   98.],       [  58.,   48.,   99.],       [  68.,   47.,   97.],       [  87.,   51.,   95.],       [ 104.,   47.,  100.],       [   1.,   68.,   95.],       [  10.,   67.,   99.],       [  24.,   66.,   97.],       [  40.,   61.,   99.],       [  60.,   66.,   95.],       [  70.,   67.,  102.],       [  84.,   62.,  101.],       [  98.,   71.,  103.],       [  -1.,   80.,  100.],       [  11.,   80.,   96.],       [  33.,   78.,  100.],       [  38.,   88.,   95.],       [  54.,   82.,  100.],       [  67.,   82.,   98.],       [  84.,   80.,  101.],       [  97.,   79.,  100.],       [  -1.,  101.,  102.],       [  17.,  101.,  101.],       [  29.,  104.,  101.],       [  39.,  102.,   97.],       [  61.,   97.,  101.],       [  70.,   96.,   96.],       [  82.,  100.,   98.],       [ 103.,  102.,   95.]])
        basis0 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0.7,  1.6,  2.9,  4. ,  4. ,  4. ,  4. ,  4. ]))
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.7,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        basis2 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0.8,  2. ,  2. ,  2. ,  2. ,  2. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 4)
        self.assertEqual(vol2.order(direction=1), 6)
        self.assertEqual(vol2.order(direction=2), 5)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 5)
        self.assertEqual(vol3.order(direction=1), 5)
        self.assertEqual(vol3.order(direction=2), 5)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 5)
        self.assertEqual(vol4.order(direction=1), 6)
        self.assertEqual(vol4.order(direction=2), 4)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)
    def test_curve_3D_p6_C0_periodic(self):
        controlpoints = np.array([[ 103.,    1.,   -1.],
       [  90.,   46.,   -1.],
       [  52.,   85.,    3.],
       [   1.,   98.,   -5.],
       [ -51.,   84.,    4.],
       [ -90.,   51.,    2.],
       [-104.,   -5.,   -3.],
       [ -87.,  -54.,    0.],
       [ -46.,  -88.,    1.],
       [  -5.,  -99.,    4.],
       [  45.,  -86.,    0.],
       [  90.,  -52.,    0.]])
        basis0 = BSplineBasis(6, np.array([-1.2,  0. ,  0. ,  0. ,  0. ,  0. ,  1.1,  2.2,  3.3,  4.3,  4.8,  6.1,  6.8,  8. ,  8. ,
        8. ,  8. ,  8. ,  9.1]),0)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 5)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p77_C0_periodic(self):
        controlpoints = np.array([[  55.,   -2.,    0.],       [  57.,   27.,   -1.],       [  36.,   47.,    3.],       [  16.,   55.,   -1.],       [ -10.,   53.,   -5.],       [ -39.,   46.,    0.],       [ -59.,   24.,    2.],       [ -61.,    0.,    3.],       [ -55.,  -22.,    4.],       [ -38.,  -44.,   -4.],       [ -18.,  -57.,   -4.],       [   8.,  -57.,   -3.],       [  33.,  -44.,    1.],       [  56.,  -22.,   -4.],       [  63.,   -4.,   -1.],       [  58.,   29.,   -3.],       [  38.,   51.,    0.],       [  12.,   58.,    0.],       [ -19.,   66.,    1.],       [ -38.,   50.,    1.],       [ -62.,   24.,   -5.],       [ -61.,   -5.,    0.],       [ -60.,  -32.,   -2.],       [ -45.,  -49.,    2.],       [ -12.,  -68.,   -3.],       [  17.,  -68.,    1.],       [  43.,  -49.,   -1.],       [  56.,  -31.,    3.],       [  70.,    4.,   -5.],       [  63.,   31.,    4.],       [  45.,   58.,    0.],       [  11.,   72.,   -5.],       [ -16.,   66.,    0.],       [ -44.,   54.,    1.],       [ -64.,   28.,    0.],       [ -66.,    2.,    1.],       [ -63.,  -32.,   -1.],       [ -43.,  -50.,   -4.],       [ -21.,  -66.,   -3.],       [  10.,  -67.,    1.],       [  41.,  -52.,    0.],       [  57.,  -33.,   -2.],       [  72.,   -3.,   -3.],       [  65.,   31.,   -3.],       [  48.,   54.,    0.],       [  13.,   68.,    0.],       [ -19.,   67.,    3.],       [ -45.,   54.,   -3.],       [ -63.,   28.,   -1.],       [ -72.,   -4.,   -3.],       [ -66.,  -37.,    4.],       [ -49.,  -59.,   -4.],       [ -16.,  -70.,    0.],       [  15.,  -72.,   -3.],       [  46.,  -54.,   -2.],       [  65.,  -37.,    1.],       [  73.,    3.,   -4.],       [  65.,   36.,   -2.],       [  44.,   63.,    4.],       [  21.,   73.,    0.],       [ -19.,   74.,   -2.],       [ -44.,   56.,   -3.],       [ -68.,   35.,   -2.],       [ -82.,    1.,   -2.],       [ -73.,  -35.,    2.],       [ -47.,  -58.,   -3.],       [ -17.,  -77.,   -4.],       [  15.,  -76.,   -5.],       [  45.,  -57.,    1.],       [  74.,  -38.,    1.],       [  86.,    4.,   -5.],       [  70.,   34.,    1.],       [  56.,   68.,   -2.],       [  16.,   83.,   -2.],       [ -19.,   75.,   -1.],       [ -53.,   61.,    0.],       [ -77.,   34.,    3.],       [ -81.,    3.,   -5.],       [ -73.,  -36.,    0.],       [ -55.,  -67.,    3.],       [ -21.,  -85.,   -2.],       [  16.,  -77.,   -2.],       [  52.,  -66.,   -4.],       [  70.,  -35.,    4.],       [  91.,    2.,   -4.],       [  73.,   38.,    0.],       [  54.,   66.,   -1.],       [  19.,   79.,   -4.],       [ -15.,   81.,    2.],       [ -53.,   63.,    0.],       [ -83.,   36.,   -1.],       [ -83.,   -5.,    0.],       [ -80.,  -41.,   -5.],       [ -54.,  -72.,    2.],       [ -18.,  -85.,    2.],       [  23.,  -87.,   -2.],       [  55.,  -73.,   -5.],       [  75.,  -38.,    4.],       [  89.,   -3.,   -1.],       [  79.,   35.,    3.],       [  60.,   68.,   -4.],       [  18.,   85.,   -2.],       [ -24.,   89.,   -2.],       [ -61.,   73.,   -3.],       [ -81.,   41.,   -3.],       [ -92.,    1.,   -4.],       [ -82.,  -36.,    3.],       [ -61.,  -67.,   -5.],       [ -20.,  -94.,    3.],       [  22.,  -89.,    1.],       [  61.,  -72.,   -1.],       [  81.,  -37.,    1.],       [  91.,   -5.,    1.],       [  83.,   46.,    3.],       [  56.,   76.,   -2.],       [  23.,   97.,   -2.],       [ -20.,   94.,   -4.],       [ -60.,   72.,   -2.],       [ -86.,   38.,    3.],       [ -99.,   -5.,   -1.],       [ -84.,  -43.,   -4.],       [ -59.,  -74.,   -1.],       [ -24.,  -90.,    3.],       [  20.,  -90.,    3.],       [  59.,  -78.,    1.],       [  90.,  -39.,   -3.],       [  97.,   -2.,    4.],       [  89.,   48.,    4.],       [  66.,   77.,   -3.],       [  20.,   95.,   -2.],       [ -25.,  102.,    0.],       [ -67.,   73.,    3.],       [ -90.,   41.,   -2.],       [-103.,    0.,   -1.],       [ -88.,  -43.,   -1.],       [ -67.,  -83.,    4.],       [ -23., -102.,   -5.],       [  25.,  -95.,   -5.],       [  59.,  -83.,   -3.],       [  95.,  -40.,    3.]])
        basis0 = BSplineBasis(7, np.array([ -0.7,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   1.2,   1.6,   3.1,   3.8,   5.2,   6.3,
         6.7,   8.3,   9. ,   9. ,   9. ,   9. ,   9. ,   9. ,  10.2]),0)
        basis1 = BSplineBasis(7, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1. ,  2. ,  3.3,  4. ,  4. ,  4. ,  4. ,  4. ,
        4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 6)
        self.assertEqual(surf2.order(direction=1), 7)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 7)
        self.assertEqual(surf3.order(direction=1), 6)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p766_C0_periodic(self):
        controlpoints = np.array([[  53.,   -5.,   -2.],       [  45.,   24.,    1.],       [  30.,   38.,   -4.],       [   6.,   50.,   -2.],       [  -9.,   46.,   -3.],       [ -36.,   41.,    2.],       [ -45.,   22.,    1.],       [ -51.,   -2.,    1.],       [ -48.,  -20.,   -4.],       [ -32.,  -36.,    4.],       [  -7.,  -47.,   -5.],       [  11.,  -45.,   -2.],       [  33.,  -40.,   -1.],       [  47.,  -20.,   -1.],       [  58.,   -1.,   -4.],       [  51.,   29.,    4.],       [  34.,   40.,   -3.],       [  10.,   58.,    1.],       [ -14.,   52.,    2.],       [ -36.,   41.,    4.],       [ -55.,   25.,    1.],       [ -58.,    3.,   -3.],       [ -47.,  -23.,   -1.],       [ -39.,  -49.,    0.],       [ -10.,  -58.,    0.],       [  10.,  -60.,    0.],       [  31.,  -42.,    3.],       [  54.,  -24.,   -4.],       [  64.,   -4.,    4.],       [  56.,   23.,   -2.],       [  36.,   47.,    3.],       [   9.,   62.,   -2.],       [ -18.,   62.,    4.],       [ -42.,   53.,    4.],       [ -56.,   28.,    1.],       [ -65.,   -5.,   -5.],       [ -53.,  -33.,   -4.],       [ -36.,  -51.,    0.],       [ -15.,  -62.,    2.],       [   9.,  -61.,    1.],       [  43.,  -47.,   -5.],       [  59.,  -31.,    1.],       [  66.,   -4.,   -4.],       [  62.,   33.,   -4.],       [  43.,   53.,   -5.],       [  17.,   64.,   -5.],       [ -15.,   65.,   -1.],       [ -46.,   49.,   -3.],       [ -67.,   31.,    3.],       [ -66.,    0.,    3.],       [ -67.,  -33.,    1.],       [ -40.,  -58.,    4.],       [ -12.,  -72.,    2.],       [  11.,  -71.,   -2.],       [  46.,  -58.,    1.],       [  57.,  -26.,   -4.],       [  71.,    4.,   -4.],       [  64.,   37.,   -4.],       [  50.,   62.,    4.],       [  14.,   69.,   -3.],       [ -22.,   70.,   -3.],       [ -48.,   63.,    2.],       [ -69.,   32.,   -4.],       [ -73.,   -3.,   -5.],       [ -68.,  -34.,    0.],       [ -50.,  -57.,   -2.],       [ -14.,  -77.,    1.],       [  15.,  -70.,    0.],       [  48.,  -58.,   -1.],       [  66.,  -28.,    2.],       [  82.,    3.,    4.],       [  69.,   32.,   -2.],       [  54.,   59.,   -1.],       [  13.,   83.,   -1.],       [ -18.,   84.,    1.],       [ -51.,   66.,    1.],       [ -77.,   37.,    4.],       [ -79.,   -4.,   -1.],       [ -70.,  -33.,   -5.],       [ -55.,  -59.,   -2.],       [ -18.,  -79.,   -5.],       [  15.,  -76.,   -4.],       [  47.,  -60.,    1.],       [  74.,  -32.,    2.],       [  92.,    1.,    1.],       [  77.,   39.,   -4.],       [  52.,   72.,    3.],       [  21.,   81.,    3.],       [ -17.,   85.,   -4.],       [ -60.,   63.,    1.],       [ -80.,   40.,   -5.],       [ -88.,    4.,    3.],       [ -82.,  -35.,   -3.],       [ -57.,  -73.,    4.],       [ -24.,  -86.,   -1.],       [  14.,  -86.,    0.],       [  50.,  -72.,   -1.],       [  77.,  -37.,   -2.],       [  96.,   -2.,    0.],       [  83.,   42.,   -4.],       [  57.,   71.,    1.],       [  22.,   93.,   -1.],       [ -26.,   94.,    0.],       [ -56.,   71.,    1.],       [ -87.,   41.,   -4.],       [ -96.,   -2.,    2.],       [ -84.,  -42.,   -3.],       [ -62.,  -78.,   -3.],       [ -21.,  -92.,   -4.],       [  24.,  -89.,   -2.],       [  53.,  -77.,   -2.],       [  85.,  -38.,    1.],       [ 101.,    2.,    0.],       [  92.,   40.,    3.],       [  57.,   80.,    2.],       [  24.,   97.,    0.],       [ -27.,   99.,   -4.],       [ -62.,   80.,   -3.],       [ -94.,   39.,    3.],       [-101.,    4.,   -1.],       [ -93.,  -41.,   -2.],       [ -62.,  -75.,    1.],       [ -25.,  -98.,    4.],       [  18., -102.,   -5.],       [  65.,  -79.,    2.],       [  85.,  -46.,   -1.],       [  51.,   -3.,    7.],       [  44.,   23.,    2.],       [  30.,   37.,   11.],       [  13.,   50.,    1.],       [ -10.,   48.,    4.],       [ -27.,   43.,    4.],       [ -41.,   21.,    1.],       [ -52.,   -5.,    6.],       [ -42.,  -23.,    6.],       [ -30.,  -40.,    7.],       [ -12.,  -49.,    2.],       [  13.,  -46.,   10.],       [  32.,  -40.,   10.],       [  48.,  -17.,    4.],       [  54.,    2.,    9.],       [  54.,   25.,    1.],       [  30.,   42.,    7.],       [  12.,   56.,   11.],       [ -16.,   57.,   11.],       [ -39.,   40.,    5.],       [ -55.,   23.,    8.],       [ -53.,    4.,    8.],       [ -52.,  -25.,    8.],       [ -32.,  -48.,    6.],       [ -14.,  -53.,    7.],       [  11.,  -57.,    5.],       [  32.,  -45.,    9.],       [  51.,  -23.,    5.],       [  63.,   -3.,    6.],       [  57.,   24.,    8.],       [  38.,   51.,    7.],       [  11.,   59.,    3.],       [ -14.,   62.,   10.],       [ -39.,   43.,    2.],       [ -54.,   23.,    3.],       [ -59.,    0.,    2.],       [ -54.,  -25.,    3.],       [ -39.,  -52.,    3.],       [ -15.,  -57.,    4.],       [  15.,  -60.,   11.],       [  37.,  -53.,    3.],       [  57.,  -25.,    5.],       [  71.,   -3.,    3.],       [  57.,   27.,    5.],       [  43.,   54.,    4.],       [  13.,   65.,    6.],       [ -18.,   66.,    3.],       [ -45.,   56.,    9.],       [ -58.,   25.,    4.],       [ -74.,    3.,    2.],       [ -66.,  -33.,    6.],       [ -46.,  -58.,    2.],       [ -14.,  -65.,    7.],       [  16.,  -65.,    2.],       [  39.,  -56.,   10.],       [  58.,  -28.,    5.],       [  78.,   -1.,    3.],       [  63.,   33.,    8.],       [  51.,   56.,    4.],       [  15.,   74.,    6.],       [ -12.,   68.,    2.],       [ -46.,   56.,    5.],       [ -64.,   37.,    2.],       [ -76.,   -3.,    4.],       [ -68.,  -36.,    4.],       [ -43.,  -57.,    3.],       [ -14.,  -72.,    6.],       [  17.,  -71.,    1.],       [  49.,  -57.,    1.],       [  66.,  -38.,    6.],       [  85.,    0.,    2.],       [  69.,   37.,   10.],       [  55.,   66.,    8.],       [  21.,   77.,    7.],       [ -22.,   77.,   10.],       [ -55.,   59.,    3.],       [ -77.,   33.,    9.],       [ -86.,   -2.,    9.],       [ -73.,  -31.,    6.],       [ -53.,  -64.,    5.],       [ -17.,  -83.,    3.],       [  15.,  -84.,    5.],       [  52.,  -63.,    7.],       [  71.,  -35.,    6.],       [  89.,    3.,    9.],       [  78.,   36.,    7.],       [  53.,   69.,    2.],       [  21.,   80.,   11.],       [ -23.,   83.,    6.],       [ -51.,   73.,    9.],       [ -80.,   34.,    3.],       [ -90.,    0.,    7.],       [ -82.,  -35.,    3.],       [ -60.,  -73.,   11.],       [ -20.,  -84.,    3.],       [  19.,  -84.,    7.],       [  55.,  -67.,   10.],       [  74.,  -33.,    4.],       [  96.,    1.,   10.],       [  84.,   37.,    4.],       [  59.,   76.,   10.],       [  24.,   88.,    8.],       [ -21.,   88.,    4.],       [ -57.,   73.,    3.],       [ -90.,   35.,   10.],       [ -91.,    3.,    2.],       [ -83.,  -42.,    3.],       [ -56.,  -72.,    4.],       [ -21.,  -88.,    6.],       [  15.,  -95.,    2.],       [  53.,  -72.,    6.],       [  82.,  -36.,   10.],       [ 103.,   -3.,   11.],       [  85.,   40.,    6.],       [  67.,   73.,    3.],       [  21.,  100.,    4.],       [ -23.,   95.,    3.],       [ -63.,   76.,    1.],       [ -95.,   44.,    5.],       [-103.,   -5.,    2.],       [ -94.,  -41.,    6.],       [ -64.,  -74.,    6.],       [ -24.,  -98.,    3.],       [  24.,  -93.,    8.],       [  60.,  -77.,    5.],       [  86.,  -48.,    5.],       [  46.,   -1.,   13.],       [  45.,   23.,   14.],       [  32.,   34.,   16.],       [  12.,   46.,   10.],       [  -9.,   52.,   14.],       [ -36.,   38.,    9.],       [ -47.,   21.,   13.],       [ -55.,   -5.,   12.],       [ -46.,  -22.,   17.],       [ -29.,  -36.,   15.],       [ -12.,  -45.,    8.],       [  11.,  -46.,   11.],       [  28.,  -44.,   10.],       [  45.,  -26.,   14.],       [  53.,   -1.,    9.],       [  46.,   25.,   15.],       [  37.,   42.,   17.],       [  10.,   59.,   13.],       [ -10.,   58.,   16.],       [ -33.,   40.,    9.],       [ -54.,   26.,    8.],       [ -54.,   -4.,   16.],       [ -50.,  -25.,   16.],       [ -37.,  -45.,   12.],       [ -14.,  -52.,   14.],       [   9.,  -58.,   18.],       [  35.,  -48.,   15.],       [  51.,  -29.,   16.],       [  57.,    1.,   18.],       [  52.,   24.,    8.],       [  35.,   53.,   11.],       [  16.,   56.,   10.],       [ -16.,   63.,   13.],       [ -38.,   51.,   14.],       [ -55.,   31.,   15.],       [ -66.,   -5.,   10.],       [ -56.,  -24.,   11.],       [ -35.,  -45.,   10.],       [ -16.,  -58.,   11.],       [  16.,  -57.,   10.],       [  38.,  -53.,   14.],       [  55.,  -30.,   17.],       [  70.,   -2.,   10.],       [  60.,   26.,    9.],       [  38.,   53.,   18.],       [  15.,   65.,    8.],       [ -14.,   70.,   10.],       [ -42.,   56.,   15.],       [ -62.,   34.,   12.],       [ -67.,   -4.,   17.],       [ -63.,  -27.,   12.],       [ -48.,  -59.,   18.],       [ -15.,  -70.,    9.],       [  17.,  -67.,   10.],       [  43.,  -51.,   13.],       [  62.,  -32.,    9.],       [  72.,    0.,    9.],       [  65.,   37.,   10.],       [  48.,   56.,    9.],       [  14.,   71.,   14.],       [ -12.,   77.,   17.],       [ -49.,   54.,   12.],       [ -67.,   28.,   15.],       [ -74.,   -3.,   13.],       [ -73.,  -30.,   14.],       [ -45.,  -58.,   12.],       [ -20.,  -77.,   13.],       [  18.,  -76.,   17.],       [  46.,  -62.,   15.],       [  66.,  -32.,    8.],       [  85.,   -5.,   15.],       [  71.,   32.,    8.],       [  51.,   60.,   16.],       [  14.,   80.,   15.],       [ -17.,   75.,   18.],       [ -52.,   59.,   14.],       [ -73.,   33.,    9.],       [ -84.,   -1.,   16.],       [ -76.,  -40.,    8.],       [ -53.,  -65.,   10.],       [ -23.,  -80.,   14.],       [  22.,  -83.,   17.],       [  45.,  -63.,   17.],       [  69.,  -41.,   10.],       [  88.,   -5.,   10.],       [  82.,   39.,   18.],       [  57.,   71.,   10.],       [  23.,   82.,   13.],       [ -24.,   81.,    9.],       [ -51.,   71.,   10.],       [ -78.,   34.,    9.],       [ -87.,   -4.,   10.],       [ -76.,  -37.,   15.],       [ -56.,  -66.,    9.],       [ -19.,  -86.,    9.],       [  15.,  -87.,   10.],       [  50.,  -69.,    8.],       [  75.,  -39.,   15.],       [  97.,    1.,   13.],       [  80.,   37.,   10.],       [  54.,   74.,    9.],       [  23.,   91.,   13.],       [ -25.,   95.,   16.],       [ -59.,   78.,   14.],       [ -88.,   42.,   15.],       [ -92.,   -2.,   10.],       [ -89.,  -44.,   15.],       [ -59.,  -75.,    8.],       [ -19.,  -97.,    9.],       [  24.,  -87.,   15.],       [  55.,  -78.,    9.],       [  87.,  -43.,   16.],       [  96.,    3.,   10.],       [  91.,   43.,   15.],       [  64.,   82.,   12.],       [  23.,   98.,   17.],       [ -23.,  100.,   17.],       [ -61.,   79.,   14.],       [ -89.,   47.,   15.],       [-104.,   -2.,   12.],       [ -91.,  -48.,    8.],       [ -58.,  -79.,   13.],       [ -19., -102.,   12.],       [  25., -101.,   13.],       [  62.,  -80.,   16.],       [  91.,  -42.,    8.],       [  54.,    2.,   18.],       [  49.,   24.,   23.],       [  27.,   35.,   15.],       [  11.,   46.,   15.],       [  -9.,   52.,   22.],       [ -35.,   40.,   18.],       [ -45.,   20.,   17.],       [ -50.,    1.,   24.],       [ -49.,  -23.,   18.],       [ -29.,  -36.,   18.],       [  -8.,  -46.,   17.],       [   9.,  -45.,   21.],       [  27.,  -36.,   24.],       [  47.,  -20.,   20.],       [  58.,    1.,   17.],       [  55.,   19.,   16.],       [  30.,   43.,   18.],       [  14.,   52.,   20.],       [ -18.,   58.,   22.],       [ -32.,   42.,   15.],       [ -56.,   22.,   18.],       [ -55.,    2.,   15.],       [ -51.,  -27.,   21.],       [ -37.,  -46.,   21.],       [ -11.,  -51.,   24.],       [  13.,  -57.,   19.],       [  30.,  -46.,   15.],       [  54.,  -23.,   23.],       [  66.,   -2.,   18.],       [  55.,   29.,   20.],       [  43.,   52.,   18.],       [  15.,   60.,   19.],       [ -14.,   62.,   23.],       [ -39.,   44.,   19.],       [ -56.,   22.,   24.],       [ -61.,   -2.,   19.],       [ -55.,  -27.,   17.],       [ -37.,  -52.,   22.],       [ -17.,  -63.,   23.],       [   8.,  -65.,   16.],       [  36.,  -54.,   15.],       [  53.,  -27.,   23.],       [  70.,    1.,   23.],       [  63.,   32.,   15.],       [  39.,   51.,   21.],       [  12.,   67.,   19.],       [ -12.,   67.,   19.],       [ -47.,   55.,   21.],       [ -66.,   33.,   22.],       [ -68.,   -3.,   19.],       [ -66.,  -32.,   21.],       [ -47.,  -54.,   24.],       [ -16.,  -65.,   18.],       [  15.,  -65.,   24.],       [  42.,  -50.,   17.],       [  61.,  -26.,   22.],       [  74.,   -5.,   20.],       [  70.,   29.,   18.],       [  43.,   61.,   21.],       [  17.,   72.,   17.],       [ -13.,   76.,   20.],       [ -48.,   58.,   19.],       [ -72.,   29.,   16.],       [ -77.,   -2.,   24.],       [ -69.,  -36.,   21.],       [ -48.,  -63.,   23.],       [ -14.,  -77.,   20.],       [  19.,  -78.,   16.],       [  42.,  -57.,   17.],       [  72.,  -33.,   15.],       [  78.,    3.,   24.],       [  74.,   40.,   22.],       [  52.,   58.,   16.],       [  18.,   77.,   17.],       [ -21.,   80.,   17.],       [ -53.,   58.,   23.],       [ -72.,   36.,   18.],       [ -85.,   -1.,   21.],       [ -76.,  -31.,   15.],       [ -53.,  -64.,   23.],       [ -16.,  -76.,   18.],       [  22.,  -83.,   15.],       [  46.,  -64.,   20.],       [  73.,  -35.,   18.],       [  83.,   -2.,   24.],       [  81.,   40.,   24.],       [  54.,   65.,   18.],       [  17.,   81.,   19.],       [ -25.,   81.,   22.],       [ -60.,   71.,   19.],       [ -81.,   40.,   24.],       [ -90.,   -4.,   23.],       [ -80.,  -36.,   20.],       [ -53.,  -71.,   17.],       [ -23.,  -91.,   24.],       [  17.,  -84.,   24.],       [  53.,  -65.,   23.],       [  79.,  -38.,   24.],       [  94.,    4.,   24.],       [  82.,   41.,   17.],       [  59.,   71.,   18.],       [  20.,   93.,   16.],       [ -24.,   90.,   24.],       [ -62.,   76.,   20.],       [ -88.,   41.,   20.],       [ -97.,    3.,   20.],       [ -89.,  -38.,   21.],       [ -57.,  -75.,   15.],       [ -22.,  -87.,   24.],       [  24.,  -90.,   19.],       [  60.,  -72.,   23.],       [  87.,  -45.,   24.],       [  95.,   -1.,   24.],       [  90.,   41.,   16.],       [  62.,   74.,   23.],       [  18.,   93.,   18.],       [ -27.,   94.,   21.],       [ -68.,   76.,   16.],       [ -87.,   42.,   20.],       [-103.,    3.,   18.],       [ -95.,  -39.,   19.],       [ -65.,  -74.,   15.],       [ -26.,  -97.,   24.],       [  27.,  -94.,   20.],       [  62.,  -82.,   21.],       [  93.,  -42.,   21.],       [  45.,   -4.,   28.],       [  46.,   25.,   22.],       [  29.,   36.,   28.],       [   8.,   51.,   29.],       [  -9.,   45.,   26.],       [ -31.,   41.,   25.],       [ -47.,   22.,   26.],       [ -52.,   -3.,   23.],       [ -45.,  -18.,   28.],       [ -31.,  -37.,   30.],       [ -17.,  -49.,   29.],       [   7.,  -51.,   28.],       [  26.,  -43.,   24.],       [  48.,  -20.,   28.],       [  52.,   -3.,   26.],       [  49.,   22.,   29.],       [  33.,   39.,   28.],       [  14.,   50.,   21.],       [ -11.,   55.,   30.],       [ -31.,   43.,   22.],       [ -49.,   20.,   30.],       [ -56.,   -4.,   23.],       [ -54.,  -23.,   27.],       [ -37.,  -48.,   31.],       [ -12.,  -54.,   24.],       [  12.,  -55.,   29.],       [  33.,  -44.,   31.],       [  49.,  -25.,   24.],       [  67.,   -2.,   30.],       [  54.,   24.,   29.],       [  38.,   50.,   28.],       [  16.,   59.,   21.],       [ -10.,   58.,   24.],       [ -39.,   52.,   25.],       [ -54.,   28.,   31.],       [ -63.,    2.,   25.],       [ -57.,  -32.,   23.],       [ -37.,  -46.,   26.],       [ -11.,  -59.,   27.],       [   9.,  -59.,   28.],       [  35.,  -54.,   29.],       [  58.,  -32.,   31.],       [  71.,    1.,   22.],       [  57.,   28.,   23.],       [  43.,   50.,   29.],       [  20.,   70.,   29.],       [ -11.,   68.,   29.],       [ -39.,   50.,   22.],       [ -63.,   27.,   30.],       [ -69.,   -2.,   30.],       [ -58.,  -31.,   23.],       [ -43.,  -54.,   25.],       [ -19.,  -64.,   28.],       [  14.,  -69.,   23.],       [  42.,  -54.,   30.],       [  59.,  -28.,   26.],       [  76.,   -5.,   25.],       [  72.,   34.,   23.],       [  45.,   60.,   30.],       [  12.,   77.,   27.],       [ -18.,   71.,   27.],       [ -47.,   56.,   28.],       [ -72.,   33.,   24.],       [ -74.,   -2.,   24.],       [ -67.,  -29.,   23.],       [ -50.,  -58.,   29.],       [ -21.,  -75.,   22.],       [  21.,  -74.,   30.],       [  48.,  -63.,   30.],       [  72.,  -31.,   25.],       [  78.,    1.,   31.],       [  74.,   30.,   30.],       [  52.,   64.,   28.],       [  16.,   80.,   31.],       [ -18.,   76.,   26.],       [ -50.,   64.,   28.],       [ -71.,   34.,   29.],       [ -80.,   -3.,   30.],       [ -77.,  -33.,   25.],       [ -47.,  -63.,   31.],       [ -18.,  -80.,   30.],       [  16.,  -83.,   28.],       [  48.,  -67.,   24.],       [  69.,  -37.,   28.],       [  87.,    0.,   28.],       [  76.,   38.,   27.],       [  49.,   70.,   26.],       [  17.,   83.,   23.],       [ -20.,   88.,   31.],       [ -55.,   72.,   29.],       [ -80.,   36.,   28.],       [ -93.,   -2.,   22.],       [ -79.,  -39.,   25.],       [ -60.,  -69.,   21.],       [ -19.,  -89.,   28.],       [  18.,  -84.,   28.],       [  52.,  -71.,   27.],       [  76.,  -36.,   28.],       [  90.,    0.,   29.],       [  83.,   39.,   27.],       [  55.,   77.,   28.],       [  24.,   89.,   22.],       [ -21.,   90.,   28.],       [ -56.,   76.,   22.],       [ -85.,   45.,   24.],       [ -95.,   -1.,   28.],       [ -89.,  -46.,   25.],       [ -60.,  -78.,   22.],       [ -18.,  -88.,   28.],       [  17.,  -90.,   23.],       [  62.,  -70.,   31.],       [  85.,  -37.,   27.],       [  97.,    2.,   28.],       [  88.,   42.,   28.],       [  58.,   83.,   26.],       [  24.,   93.,   30.],       [ -28.,   93.,   30.],       [ -64.,   80.,   28.],       [ -92.,   46.,   24.],       [-101.,    3.,   30.],       [ -87.,  -47.,   22.],       [ -60.,  -81.,   26.],       [ -25., -102.,   30.],       [  22., -103.,   27.],       [  60.,  -77.,   28.],       [  86.,  -46.,   22.],       [  47.,   -5.,   29.],       [  40.,   24.,   35.],       [  35.,   38.,   33.],       [  11.,   44.,   32.],       [ -14.,   45.,   36.],       [ -35.,   34.,   28.],       [ -42.,   20.,   31.],       [ -49.,    3.,   31.],       [ -42.,  -25.,   37.],       [ -30.,  -42.,   32.],       [ -16.,  -53.,   35.],       [   8.,  -51.,   30.],       [  28.,  -35.,   32.],       [  41.,  -24.,   35.],       [  52.,   -3.,   29.],       [  54.,   19.,   32.],       [  33.,   45.,   37.],       [  16.,   54.,   37.],       [ -18.,   51.,   29.],       [ -37.,   48.,   35.],       [ -52.,   26.,   32.],       [ -59.,    4.,   31.],       [ -54.,  -30.,   35.],       [ -36.,  -48.,   35.],       [ -14.,  -53.,   28.],       [   9.,  -57.,   36.],       [  33.,  -45.,   36.],       [  49.,  -25.,   38.],       [  58.,    4.,   35.],       [  58.,   25.,   30.],       [  42.,   50.,   31.],       [   9.,   62.,   28.],       [ -18.,   59.,   33.],       [ -38.,   44.,   34.],       [ -58.,   22.,   31.],       [ -68.,    4.,   33.],       [ -61.,  -25.,   30.],       [ -35.,  -52.,   32.],       [ -14.,  -57.,   36.],       [  16.,  -57.,   35.],       [  35.,  -50.,   30.],       [  54.,  -24.,   36.],       [  69.,   -2.,   36.],       [  59.,   27.,   32.],       [  42.,   58.,   32.],       [  10.,   69.,   29.],       [ -12.,   68.,   29.],       [ -46.,   50.,   34.],       [ -60.,   30.,   30.],       [ -70.,    3.,   36.],       [ -62.,  -29.,   31.],       [ -42.,  -56.,   37.],       [ -18.,  -69.,   34.],       [  15.,  -63.,   33.],       [  43.,  -56.,   35.],       [  64.,  -31.,   33.],       [  71.,    3.,   32.],       [  66.,   37.,   29.],       [  48.,   59.,   34.],       [  12.,   73.,   34.],       [ -17.,   73.,   32.],       [ -45.,   63.,   37.],       [ -69.,   31.,   32.],       [ -76.,    0.,   32.],       [ -69.,  -32.,   34.],       [ -49.,  -59.,   31.],       [ -17.,  -70.,   34.],       [  12.,  -73.,   34.],       [  45.,  -54.,   35.],       [  63.,  -37.,   35.],       [  85.,    2.,   32.],       [  76.,   32.,   32.],       [  51.,   64.,   34.],       [  22.,   80.,   30.],       [ -20.,   76.,   34.],       [ -50.,   66.,   35.],       [ -73.,   33.,   30.],       [ -86.,   -1.,   28.],       [ -73.,  -40.,   28.],       [ -52.,  -62.,   36.],       [ -16.,  -81.,   32.],       [  16.,  -77.,   38.],       [  47.,  -67.,   30.],       [  74.,  -39.,   34.],       [  92.,    2.,   35.],       [  76.,   42.,   34.],       [  53.,   68.,   37.],       [  16.,   89.,   31.],       [ -24.,   80.,   28.],       [ -50.,   63.,   34.],       [ -78.,   38.,   29.],       [ -85.,   -1.,   35.],       [ -80.,  -36.,   29.],       [ -59.,  -69.,   32.],       [ -22.,  -90.,   38.],       [  23.,  -83.,   37.],       [  51.,  -68.,   35.],       [  82.,  -41.,   36.],       [  90.,   -4.,   28.],       [  88.,   40.,   33.],       [  59.,   69.,   33.],       [  20.,   92.,   37.],       [ -19.,   94.,   32.],       [ -63.,   71.,   37.],       [ -82.,   37.,   29.],       [ -94.,   -3.,   32.],       [ -81.,  -43.,   37.],       [ -56.,  -71.,   31.],       [ -22.,  -93.,   36.],       [  18.,  -94.,   34.],       [  60.,  -69.,   31.],       [  88.,  -45.,   35.],       [ 102.,   -3.,   33.],       [  94.,   44.,   37.],       [  66.,   81.,   28.],       [  19.,   98.,   38.],       [ -27.,   98.,   34.],       [ -67.,   78.,   30.],       [ -88.,   47.,   38.],       [-101.,    0.,   28.],       [ -89.,  -41.,   29.],       [ -60.,  -81.,   35.],       [ -24.,  -99.,   32.],       [  26.,  -98.,   32.],       [  62.,  -83.,   30.],       [  94.,  -48.,   30.],       [  49.,    0.,   35.],       [  40.,   23.,   40.],       [  29.,   39.,   37.],       [  14.,   52.,   39.],       [  -8.,   43.,   38.],       [ -30.,   41.,   42.],       [ -49.,   25.,   41.],       [ -47.,    1.,   44.],       [ -47.,  -24.,   40.],       [ -27.,  -36.,   37.],       [ -12.,  -52.,   44.],       [  16.,  -53.,   42.],       [  34.,  -40.,   44.],       [  49.,  -23.,   43.],       [  59.,    1.,   43.],       [  47.,   23.,   38.],       [  34.,   39.,   39.],       [  11.,   50.,   43.],       [ -13.,   59.,   36.],       [ -35.,   47.,   41.],       [ -54.,   28.,   35.],       [ -58.,   -4.,   40.],       [ -50.,  -24.,   36.],       [ -32.,  -44.,   41.],       [ -10.,  -59.,   41.],       [  17.,  -55.,   35.],       [  33.,  -46.,   37.],       [  53.,  -23.,   39.],       [  58.,   -1.,   40.],       [  52.,   27.,   35.],       [  41.,   49.,   41.],       [  18.,   58.,   35.],       [ -17.,   64.,   36.],       [ -35.,   49.,   35.],       [ -53.,   27.,   39.],       [ -66.,   -4.,   44.],       [ -60.,  -32.,   36.],       [ -41.,  -50.,   44.],       [ -15.,  -58.,   37.],       [  17.,  -57.,   37.],       [  35.,  -47.,   44.],       [  57.,  -25.,   40.],       [  65.,   -5.,   43.],       [  60.,   26.,   44.],       [  43.,   51.,   37.],       [  15.,   63.,   40.],       [ -19.,   67.,   36.],       [ -46.,   53.,   44.],       [ -59.,   28.,   38.],       [ -68.,    1.,   44.],       [ -67.,  -29.,   43.],       [ -41.,  -54.,   43.],       [ -13.,  -66.,   35.],       [  19.,  -71.,   35.],       [  44.,  -55.,   38.],       [  66.,  -28.,   43.],       [  72.,   -3.,   42.],       [  72.,   28.,   37.],       [  47.,   60.,   36.],       [  17.,   70.,   39.],       [ -21.,   72.,   37.],       [ -46.,   62.,   41.],       [ -71.,   31.,   40.],       [ -73.,    3.,   42.],       [ -66.,  -38.,   36.],       [ -48.,  -60.,   44.],       [ -18.,  -72.,   40.],       [  17.,  -77.,   38.],       [  43.,  -59.,   38.],       [  72.,  -30.,   36.],       [  77.,   -5.,   35.],       [  70.,   35.,   35.],       [  49.,   61.,   36.],       [  14.,   82.,   36.],       [ -18.,   75.,   37.],       [ -46.,   68.,   36.],       [ -79.,   36.,   40.],       [ -82.,   -3.,   41.],       [ -72.,  -34.,   44.],       [ -50.,  -68.,   38.],       [ -21.,  -79.,   42.],       [  19.,  -75.,   38.],       [  55.,  -64.,   38.],       [  71.,  -39.,   43.],       [  89.,    1.,   44.],       [  81.,   35.,   43.],       [  56.,   68.,   38.],       [  16.,   81.,   41.],       [ -19.,   81.,   37.],       [ -52.,   69.,   36.],       [ -76.,   39.,   43.],       [ -87.,    3.,   39.],       [ -83.,  -43.,   39.],       [ -53.,  -65.,   39.],       [ -19.,  -82.,   36.],       [  22.,  -88.,   43.],       [  53.,  -73.,   40.],       [  83.,  -38.,   35.],       [  96.,    3.,   39.],       [  86.,   40.,   39.],       [  58.,   69.,   43.],       [  22.,   95.,   41.],       [ -22.,   87.,   43.],       [ -57.,   77.,   38.],       [ -82.,   41.,   39.],       [ -92.,    3.,   39.],       [ -86.,  -43.,   39.],       [ -59.,  -70.,   44.],       [ -18.,  -89.,   42.],       [  23.,  -93.,   42.],       [  60.,  -75.,   35.],       [  80.,  -39.,   35.],       [ 104.,    3.,   43.],       [  85.,   42.,   37.],       [  61.,   78.,   35.],       [  25.,   95.,   42.],       [ -25.,   93.,   44.],       [ -59.,   82.,   41.],       [ -86.,   43.,   40.],       [-100.,   -2.,   36.],       [ -96.,  -43.,   38.],       [ -64.,  -81.,   43.],       [ -18.,  -95.,   39.],       [  24.,  -96.,   38.],       [  65.,  -83.,   41.],       [  89.,  -45.,   38.]])
        basis0 = BSplineBasis(7, np.array([-1.4,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.7,  2. ,  2.7,  4.3,  4.7,  6.2,  7. ,  7.6,
        9. ,  9. ,  9. ,  9. ,  9. ,  9. ,  9.7]),0)
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.8,  2. ,  3. ,  4. ,  4. ,  4. ,  4. ,  4. ,  4. ]))
        basis2 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.2,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 6)
        self.assertEqual(vol2.order(direction=1), 6)
        self.assertEqual(vol2.order(direction=2), 6)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 7)
        self.assertEqual(vol3.order(direction=1), 5)
        self.assertEqual(vol3.order(direction=2), 6)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 7)
        self.assertEqual(vol4.order(direction=1), 6)
        self.assertEqual(vol4.order(direction=2), 5)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)
    def test_curve_3D_p5_C1_periodic(self):
        controlpoints = np.array([[  96.,   -5.,    0.],
       [  83.,   46.,   -1.],
       [  47.,   86.,    3.],
       [  -5.,   98.,   -2.],
       [ -52.,   87.,   -4.],
       [ -87.,   45.,   -2.],
       [-105.,    1.,   -1.],
       [ -85.,  -48.,    1.],
       [ -51.,  -90.,   -5.],
       [   4.,  -97.,    3.],
       [  45.,  -91.,   -3.],
       [  83.,  -49.,    2.]])
        basis0 = BSplineBasis(5, np.array([ -2. ,  -0.7,   0. ,   0. ,   0. ,   0.7,   2.2,   2.6,   3.8,   5.2,   6. ,   6.9,   8. ,
         9.3,  10. ,  10. ,  10. ,  10.7,  12.2]),1)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 4)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p67_C1_periodic(self):
        controlpoints = np.array([[  63.,   -5.,   -4.],       [  58.,   26.,   -3.],       [  39.,   44.,   -3.],       [  13.,   58.,    2.],       [ -12.,   58.,   -3.],       [ -35.,   50.,    3.],       [ -57.,   23.,   -4.],       [ -61.,    1.,   -4.],       [ -59.,  -24.,   -3.],       [ -42.,  -46.,    0.],       [ -13.,  -62.,   -2.],       [  16.,  -63.,    3.],       [  35.,  -47.,    2.],       [  50.,  -30.,    3.],       [  60.,    4.,    2.],       [  61.,   32.,    0.],       [  38.,   48.,   -1.],       [  16.,   60.,    0.],       [ -16.,   64.,   -4.],       [ -37.,   48.,   -2.],       [ -55.,   31.,    3.],       [ -70.,    0.,   -5.],       [ -54.,  -33.,   -1.],       [ -37.,  -49.,   -5.],       [ -10.,  -66.,   -3.],       [   9.,  -64.,    4.],       [  41.,  -46.,    3.],       [  62.,  -24.,   -4.],       [  65.,    0.,    3.],       [  58.,   30.,   -5.],       [  37.,   56.,    4.],       [  14.,   62.,   -5.],       [ -14.,   71.,    1.],       [ -45.,   50.,    4.],       [ -66.,   32.,    0.],       [ -65.,   -5.,    3.],       [ -67.,  -35.,    2.],       [ -43.,  -55.,    0.],       [ -16.,  -64.,   -4.],       [  19.,  -65.,    1.],       [  41.,  -57.,    0.],       [  57.,  -31.,    4.],       [  76.,   -4.,    1.],       [  70.,   30.,   -4.],       [  43.,   54.,   -4.],       [  19.,   72.,    2.],       [ -21.,   71.,   -4.],       [ -47.,   59.,   -5.],       [ -66.,   27.,    0.],       [ -71.,    4.,    4.],       [ -62.,  -30.,    2.],       [ -43.,  -56.,   -2.],       [ -21.,  -72.,    3.],       [  18.,  -71.,   -3.],       [  43.,  -61.,   -4.],       [  68.,  -35.,    2.],       [  81.,   -1.,   -5.],       [  66.,   33.,    2.],       [  47.,   57.,    1.],       [  14.,   76.,    4.],       [ -15.,   71.,   -5.],       [ -46.,   59.,   -2.],       [ -69.,   31.,   -2.],       [ -79.,   -5.,    1.],       [ -75.,  -35.,   -2.],       [ -48.,  -63.,   -5.],       [ -20.,  -72.,    0.],       [  13.,  -75.,   -3.],       [  46.,  -62.,   -2.],       [  65.,  -31.,   -5.],       [  77.,   -4.,   -1.],       [  71.,   31.,    2.],       [  49.,   61.,   -3.],       [  20.,   81.,    2.],       [ -17.,   75.,   -2.],       [ -47.,   68.,    4.],       [ -70.,   34.,    4.],       [ -82.,   -3.,   -5.],       [ -76.,  -36.,   -2.],       [ -51.,  -69.,   -5.],       [ -23.,  -82.,   -3.],       [  15.,  -76.,   -3.],       [  51.,  -69.,   -2.],       [  72.,  -36.,   -1.],       [  83.,   -2.,    2.],       [  79.,   36.,    2.],       [  53.,   65.,   -3.],       [  19.,   88.,    0.],       [ -18.,   86.,   -3.],       [ -52.,   67.,   -4.],       [ -83.,   36.,   -2.],       [ -85.,   -5.,   -3.],       [ -74.,  -34.,    4.],       [ -51.,  -72.,    4.],       [ -19.,  -87.,    3.],       [  18.,  -90.,   -4.],       [  55.,  -68.,    0.],       [  76.,  -36.,    1.],       [  91.,    0.,   -5.],       [  82.,   35.,    3.],       [  52.,   73.,   -5.],       [  15.,   87.,    4.],       [ -16.,   92.,    4.],       [ -56.,   73.,   -5.],       [ -83.,   39.,    4.],       [ -92.,    3.,   -5.],       [ -87.,  -38.,    4.],       [ -61.,  -75.,   -4.],       [ -23.,  -88.,   -4.],       [  19.,  -86.,    4.],       [  56.,  -69.,   -4.],       [  78.,  -42.,    1.],       [  97.,    3.,   -3.],       [  88.,   42.,    4.],       [  64.,   74.,   -4.],       [  17.,   97.,    0.],       [ -17.,   94.,   -3.],       [ -56.,   77.,   -3.],       [ -85.,   45.,   -1.],       [ -99.,    4.,    2.],       [ -82.,  -42.,    4.],       [ -60.,  -72.,   -3.],       [ -24.,  -92.,   -3.],       [  21.,  -93.,   -1.],       [  57.,  -76.,    0.],       [  88.,  -40.,    1.],       [ 103.,    0.,    0.],       [  92.,   41.,    3.],       [  63.,   76.,   -3.],       [  19.,   93.,    0.],       [ -18.,  101.,    4.],       [ -67.,   73.,    0.],       [ -91.,   46.,    2.],       [-102.,   -3.,   -5.],       [ -88.,  -48.,    1.],       [ -68.,  -78.,    1.],       [ -20.,  -95.,    2.],       [  19.,  -96.,    3.],       [  57.,  -77.,   -3.],       [  90.,  -39.,    4.]])
        basis0 = BSplineBasis(6, np.array([ -1.9,  -0.8,   0. ,   0. ,   0. ,   0. ,   0.6,   2. ,   3.4,   4.1,   4.7,   5.8,   7.4,
         7.9,   9.1,  10.2,  11. ,  11. ,  11. ,  11. ,  11.6,  13. ]),1)
        basis1 = BSplineBasis(7, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.3,  2.1,  3.3,  4. ,  4. ,  4. ,  4. ,  4. ,
        4. ,  4. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 5)
        self.assertEqual(surf2.order(direction=1), 7)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 6)
        self.assertEqual(surf3.order(direction=1), 6)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p765_C1_periodic(self):
        controlpoints = np.array([[  52.,    2.,    2.],       [  44.,   15.,    3.],       [  39.,   37.,   -5.],       [  18.,   41.,   -5.],       [   3.,   46.,    0.],       [ -24.,   42.,    3.],       [ -40.,   37.,   -2.],       [ -44.,   19.,    4.],       [ -53.,   -1.,    4.],       [ -48.,  -23.,   -5.],       [ -32.,  -35.,   -2.],       [ -22.,  -46.,    4.],       [  -4.,  -52.,   -5.],       [  16.,  -51.,   -1.],       [  32.,  -31.,    1.],       [  49.,  -23.,    0.],       [  62.,   -5.,   -4.],       [  51.,   24.,   -5.],       [  44.,   42.,   -5.],       [  26.,   53.,    0.],       [   0.,   59.,    0.],       [ -21.,   56.,    0.],       [ -38.,   40.,   -2.],       [ -54.,   20.,    3.],       [ -59.,    0.,    0.],       [ -50.,  -25.,   -4.],       [ -41.,  -38.,   -2.],       [ -25.,  -51.,    4.],       [   0.,  -58.,    1.],       [  17.,  -52.,    2.],       [  44.,  -39.,   -5.],       [  57.,  -21.,   -4.],       [  62.,    1.,   -4.],       [  60.,   27.,   -5.],       [  47.,   45.,    0.],       [  23.,   60.,    1.],       [  -1.,   65.,   -5.],       [ -27.,   58.,    0.],       [ -43.,   50.,    1.],       [ -63.,   25.,   -2.],       [ -69.,   -2.,   -2.],       [ -58.,  -22.,   -2.],       [ -46.,  -46.,    0.],       [ -24.,  -60.,   -2.],       [  -3.,  -68.,   -2.],       [  22.,  -60.,    0.],       [  50.,  -47.,   -2.],       [  61.,  -27.,    1.],       [  71.,    1.,   -5.],       [  66.,   32.,   -5.],       [  55.,   56.,   -3.],       [  32.,   72.,    0.],       [  -4.,   76.,    2.],       [ -30.,   69.,   -1.],       [ -57.,   50.,    1.],       [ -70.,   27.,    1.],       [ -79.,   -3.,   -4.],       [ -72.,  -29.,    3.],       [ -49.,  -57.,   -4.],       [ -29.,  -67.,    1.],       [  -2.,  -75.,   -4.],       [  24.,  -71.,   -5.],       [  48.,  -51.,   -3.],       [  73.,  -25.,   -4.],       [  80.,    2.,   -2.],       [  74.,   28.,    1.],       [  55.,   62.,    1.],       [  36.,   78.,   -1.],       [  -4.,   85.,   -4.],       [ -31.,   81.,    0.],       [ -58.,   62.,   -3.],       [ -75.,   32.,    3.],       [ -85.,    3.,    4.],       [ -77.,  -31.,   -1.],       [ -55.,  -56.,   -4.],       [ -36.,  -74.,   -5.],       [   1.,  -85.,    4.],       [  31.,  -82.,   -5.],       [  63.,  -63.,   -5.],       [  75.,  -29.,    4.],       [  88.,    1.,   -5.],       [  86.,   33.,   -2.],       [  69.,   69.,    0.],       [  35.,   86.,    0.],       [   0.,   88.,   -3.],       [ -39.,   88.,    3.],       [ -64.,   61.,    4.],       [ -88.,   39.,    0.],       [ -93.,    2.,   -3.],       [ -85.,  -32.,    4.],       [ -65.,  -62.,    3.],       [ -37.,  -81.,   -4.],       [   2.,  -94.,    3.],       [  30.,  -84.,   -1.],       [  60.,  -64.,   -4.],       [  89.,  -40.,   -1.],       [ 103.,   -1.,    2.],       [  92.,   42.,   -1.],       [  71.,   71.,   -3.],       [  38.,   91.,    2.],       [   1.,  104.,   -4.],       [ -36.,   94.,    0.],       [ -70.,   70.,   -2.],       [ -94.,   40.,    1.],       [ -99.,    0.,    0.],       [ -98.,  -34.,   -2.],       [ -70.,  -70.,    3.],       [ -40.,  -89.,   -3.],       [   1., -102.,   -2.],       [  35.,  -94.,   -1.],       [  67.,  -73.,   -2.],       [  89.,  -38.,    4.],       [  54.,    4.,    4.],       [  50.,   17.,    4.],       [  38.,   34.,    3.],       [  14.,   43.,    9.],       [  -5.,   51.,    4.],       [ -22.,   46.,   11.],       [ -36.,   31.,    8.],       [ -48.,   19.,    4.],       [ -46.,   -2.,   10.],       [ -52.,  -21.,   12.],       [ -39.,  -31.,    3.],       [ -15.,  -47.,    7.],       [   3.,  -48.,    7.],       [  16.,  -51.,    7.],       [  38.,  -36.,    6.],       [  46.,  -23.,    3.],       [  56.,    1.,    8.],       [  50.,   18.,   11.],       [  40.,   43.,   12.],       [  25.,   53.,   12.],       [  -4.,   59.,   10.],       [ -26.,   50.,    5.],       [ -39.,   42.,   10.],       [ -59.,   24.,    3.],       [ -54.,   -5.,   10.],       [ -57.,  -21.,    8.],       [ -42.,  -41.,    7.],       [ -22.,  -50.,   11.],       [  -3.,  -54.,    3.],       [  20.,  -57.,    4.],       [  36.,  -43.,    9.],       [  51.,  -20.,    3.],       [  64.,    3.,    6.],       [  61.,   23.,    7.],       [  45.,   45.,   10.],       [  30.,   58.,    8.],       [  -4.,   65.,   11.],       [ -23.,   59.,    6.],       [ -46.,   45.,    7.],       [ -66.,   29.,    6.],       [ -66.,   -1.,   11.],       [ -59.,  -26.,    5.],       [ -43.,  -48.,    4.],       [ -23.,  -64.,   11.],       [  -2.,  -67.,    4.],       [  29.,  -62.,    6.],       [  48.,  -48.,    5.],       [  64.,  -24.,    6.],       [  75.,   -3.,    5.],       [  65.,   32.,    3.],       [  52.,   54.,    8.],       [  25.,   69.,   12.],       [  -5.,   77.,    5.],       [ -30.,   72.,   12.],       [ -57.,   48.,    4.],       [ -71.,   27.,    8.],       [ -76.,    1.,   10.],       [ -71.,  -24.,    4.],       [ -52.,  -49.,    7.],       [ -33.,  -66.,    7.],       [   0.,  -76.,   12.],       [  25.,  -70.,    7.],       [  48.,  -51.,    5.],       [  73.,  -28.,    4.],       [  88.,    3.,    7.],       [  80.,   28.,    5.],       [  54.,   59.,    9.],       [  29.,   73.,    8.],       [  -4.,   86.,   10.],       [ -36.,   75.,    7.],       [ -64.,   59.,    5.],       [ -79.,   30.,   12.],       [ -84.,    0.,    9.],       [ -74.,  -37.,    9.],       [ -57.,  -61.,    6.],       [ -30.,  -80.,    9.],       [  -3.,  -87.,   10.],       [  29.,  -74.,    4.],       [  63.,  -58.,   12.],       [  78.,  -31.,    4.],       [  88.,   -1.,    8.],       [  86.,   37.,    5.],       [  65.,   66.,    7.],       [  39.,   80.,   12.],       [  -4.,   95.,    3.],       [ -32.,   88.,    3.],       [ -61.,   68.,    5.],       [ -83.,   35.,    8.],       [ -93.,   -5.,    9.],       [ -86.,  -38.,    8.],       [ -65.,  -69.,    9.],       [ -38.,  -84.,    9.],       [   3.,  -88.,    3.],       [  30.,  -81.,    8.],       [  65.,  -69.,    9.],       [  82.,  -34.,    4.],       [ 103.,   -1.,   12.],       [  93.,   42.,    5.],       [  68.,   65.,    7.],       [  40.,   89.,   11.],       [   0.,   99.,    7.],       [ -43.,   94.,    4.],       [ -72.,   75.,   10.],       [ -96.,   35.,    9.],       [-105.,    4.,   11.],       [ -93.,  -35.,   11.],       [ -68.,  -67.,    5.],       [ -38.,  -93.,   10.],       [  -1.,  -99.,   12.],       [  42.,  -89.,    8.],       [  74.,  -75.,    5.],       [  89.,  -36.,    6.],       [  54.,    0.,   13.],       [  47.,   19.,   11.],       [  38.,   32.,   17.],       [  20.,   41.,   15.],       [   1.,   50.,   11.],       [ -18.,   48.,   17.],       [ -39.,   38.,   20.],       [ -43.,   18.,   14.],       [ -55.,    4.,   16.],       [ -47.,  -23.,   14.],       [ -35.,  -40.,   13.],       [ -23.,  -51.,   13.],       [   2.,  -52.,   19.],       [  17.,  -50.,   11.],       [  32.,  -35.,   18.],       [  44.,  -19.,   13.],       [  58.,   -4.,   18.],       [  57.,   21.,   11.],       [  44.,   43.,   14.],       [  22.,   49.,   11.],       [  -2.,   53.,   13.],       [ -24.,   54.,   14.],       [ -40.,   44.,   17.],       [ -49.,   18.,   12.],       [ -63.,    4.,   18.],       [ -58.,  -27.,   16.],       [ -40.,  -39.,   13.],       [ -27.,  -59.,   14.],       [   3.,  -59.,   19.],       [  26.,  -51.,   17.],       [  43.,  -38.,   20.],       [  57.,  -21.,   11.],       [  65.,   -2.,   11.],       [  56.,   27.,   18.],       [  45.,   50.,   18.],       [  27.,   61.,   20.],       [  -3.,   62.,   14.],       [ -27.,   57.,   15.],       [ -48.,   44.,   15.],       [ -60.,   22.,   18.],       [ -65.,   -3.,   18.],       [ -57.,  -25.,   17.],       [ -45.,  -49.,   11.],       [ -24.,  -59.,   16.],       [  -1.,  -67.,   16.],       [  20.,  -65.,   13.],       [  42.,  -47.,   11.],       [  64.,  -22.,   19.],       [  73.,   -4.,   11.],       [  71.,   28.,   13.],       [  52.,   48.,   11.],       [  24.,   66.,   13.],       [   2.,   71.,   17.],       [ -30.,   74.,   17.],       [ -57.,   57.,   16.],       [ -67.,   27.,   19.],       [ -76.,    2.,   13.],       [ -74.,  -33.,   11.],       [ -54.,  -49.,   13.],       [ -31.,  -73.,   13.],       [   1.,  -78.,   18.],       [  33.,  -71.,   12.],       [  54.,  -55.,   12.],       [  70.,  -33.,   12.],       [  84.,   -2.,   11.],       [  75.,   35.,   11.],       [  55.,   58.,   16.],       [  30.,   81.,   16.],       [   4.,   85.,   17.],       [ -35.,   81.,   13.],       [ -62.,   62.,   13.],       [ -76.,   36.,   16.],       [ -82.,    1.,   14.],       [ -80.,  -37.,   11.],       [ -63.,  -62.,   18.],       [ -32.,  -77.,   16.],       [  -5.,  -81.,   12.],       [  33.,  -77.,   19.],       [  58.,  -57.,   12.],       [  75.,  -33.,   18.],       [  90.,   -4.,   12.],       [  87.,   30.,   13.],       [  68.,   60.,   17.],       [  39.,   87.,   19.],       [   3.,   92.,   13.],       [ -35.,   85.,   13.],       [ -68.,   69.,   13.],       [ -80.,   35.,   16.],       [ -96.,    3.,   17.],       [ -89.,  -37.,   13.],       [ -70.,  -68.,   14.],       [ -32.,  -84.,   14.],       [   0.,  -91.,   15.],       [  32.,  -88.,   11.],       [  68.,  -69.,   12.],       [  89.,  -32.,   18.],       [ 100.,   -2.,   14.],       [  96.,   38.,   17.],       [  68.,   71.,   14.],       [  40.,   88.,   18.],       [   0.,   98.,   20.],       [ -35.,   94.,   20.],       [ -74.,   70.,   13.],       [ -97.,   40.,   16.],       [ -96.,   -3.,   14.],       [ -95.,  -39.,   14.],       [ -75.,  -67.,   18.],       [ -40.,  -90.,   19.],       [  -3., -100.,   18.],       [  40.,  -93.,   15.],       [  70.,  -71.,   14.],       [  95.,  -39.,   17.],       [  47.,   -3.,   26.],       [  43.,   16.,   28.],       [  36.,   35.,   19.],       [  16.,   43.,   28.],       [   2.,   51.,   19.],       [ -19.,   48.,   24.],       [ -32.,   38.,   20.],       [ -43.,   21.,   28.],       [ -50.,   -5.,   25.],       [ -49.,  -17.,   21.],       [ -36.,  -33.,   21.],       [ -15.,  -46.,   27.],       [  -2.,  -50.,   22.],       [  22.,  -45.,   27.],       [  37.,  -34.,   21.],       [  43.,  -15.,   21.],       [  57.,   -2.,   21.],       [  56.,   21.,   28.],       [  44.,   40.,   25.],       [  26.,   56.,   26.],       [   0.,   55.,   24.],       [ -21.,   53.,   27.],       [ -42.,   43.,   23.],       [ -59.,   27.,   24.],       [ -61.,    3.,   23.],       [ -56.,  -19.,   26.],       [ -45.,  -46.,   28.],       [ -25.,  -50.,   22.],       [  -4.,  -60.,   25.],       [  19.,  -58.,   21.],       [  45.,  -38.,   24.],       [  49.,  -27.,   25.],       [  62.,    2.,   26.],       [  60.,   29.,   25.],       [  42.,   51.,   19.],       [  24.,   59.,   19.],       [  -5.,   62.,   22.],       [ -28.,   61.,   23.],       [ -50.,   44.,   24.],       [ -57.,   29.,   20.],       [ -69.,   -2.,   19.],       [ -66.,  -29.,   21.],       [ -52.,  -48.,   27.],       [ -24.,  -62.,   27.],       [   3.,  -67.,   20.],       [  29.,  -60.,   25.],       [  44.,  -43.,   19.],       [  56.,  -29.,   20.],       [  77.,    4.,   23.],       [  66.,   24.,   26.],       [  52.,   52.,   22.],       [  30.,   64.,   19.],       [  -1.,   75.,   26.],       [ -31.,   72.,   28.],       [ -52.,   51.,   23.],       [ -66.,   30.,   24.],       [ -73.,   -2.,   22.],       [ -75.,  -31.,   26.],       [ -52.,  -49.,   23.],       [ -24.,  -73.,   24.],       [   1.,  -73.,   21.],       [  29.,  -71.,   28.],       [  52.,  -50.,   27.],       [  65.,  -28.,   26.],       [  84.,    2.,   25.],       [  72.,   34.,   24.],       [  59.,   61.,   27.],       [  28.,   74.,   20.],       [  -5.,   81.,   26.],       [ -31.,   80.,   28.],       [ -55.,   55.,   19.],       [ -75.,   31.,   22.],       [ -88.,   -3.,   23.],       [ -76.,  -28.,   19.],       [ -55.,  -61.,   25.],       [ -34.,  -74.,   25.],       [   2.,  -81.,   26.],       [  32.,  -77.,   20.],       [  56.,  -57.,   25.],       [  77.,  -36.,   28.],       [  93.,   -3.,   25.],       [  88.,   33.,   24.],       [  61.,   66.,   27.],       [  31.,   88.,   20.],       [   3.,   91.,   26.],       [ -34.,   82.,   23.],       [ -68.,   68.,   25.],       [ -85.,   35.,   26.],       [ -94.,    1.,   27.],       [ -83.,  -40.,   25.],       [ -65.,  -70.,   20.],       [ -40.,  -90.,   25.],       [   3.,  -94.,   28.],       [  30.,  -84.,   19.],       [  61.,  -66.,   25.],       [  81.,  -38.,   28.],       [ 100.,   -2.,   20.],       [  95.,   37.,   23.],       [  68.,   66.,   28.],       [  39.,   96.,   27.],       [   4.,   99.,   23.],       [ -40.,   97.,   28.],       [ -73.,   71.,   28.],       [ -90.,   33.,   25.],       [-101.,    0.,   24.],       [ -96.,  -43.,   25.],       [ -67.,  -69.,   20.],       [ -39.,  -92.,   28.],       [  -2.,  -98.,   28.],       [  38.,  -90.,   22.],       [  72.,  -67.,   23.],       [  96.,  -41.,   27.],       [  46.,    4.,   27.],       [  48.,   21.,   27.],       [  36.,   33.,   32.],       [  19.,   42.,   36.],       [  -5.,   47.,   31.],       [ -22.,   50.,   32.],       [ -38.,   33.,   29.],       [ -48.,   23.,   27.],       [ -49.,   -3.,   29.],       [ -47.,  -23.,   31.],       [ -37.,  -35.,   28.],       [ -20.,  -47.,   30.],       [  -5.,  -50.,   30.],       [  16.,  -46.,   34.],       [  38.,  -39.,   34.],       [  44.,  -15.,   36.],       [  54.,   -1.,   34.],       [  51.,   27.,   34.],       [  45.,   37.,   34.],       [  22.,   55.,   30.],       [   1.,   59.,   27.],       [ -27.,   49.,   35.],       [ -44.,   43.,   36.],       [ -49.,   19.,   27.],       [ -60.,    4.,   36.],       [ -50.,  -18.,   27.],       [ -40.,  -39.,   34.],       [ -19.,  -52.,   28.],       [  -4.,  -61.,   33.],       [  22.,  -58.,   36.],       [  36.,  -38.,   27.],       [  50.,  -22.,   27.],       [  66.,   -2.,   27.],       [  65.,   21.,   36.],       [  43.,   45.,   28.],       [  21.,   60.,   34.],       [   1.,   67.,   30.],       [ -21.,   60.,   32.],       [ -43.,   43.,   28.],       [ -62.,   28.,   31.],       [ -66.,   -5.,   28.],       [ -59.,  -24.,   28.],       [ -45.,  -51.,   28.],       [ -26.,  -57.,   30.],       [  -1.,  -64.,   35.],       [  25.,  -64.,   33.],       [  47.,  -50.,   35.],       [  56.,  -29.,   36.],       [  77.,   -1.,   31.],       [  66.,   31.,   36.],       [  50.,   51.,   32.],       [  32.,   65.,   33.],       [   0.,   77.,   29.],       [ -29.,   69.,   34.],       [ -51.,   57.,   27.],       [ -66.,   25.,   35.],       [ -79.,   -4.,   30.],       [ -74.,  -31.,   35.],       [ -52.,  -52.,   36.],       [ -32.,  -66.,   36.],       [  -2.,  -73.,   35.],       [  29.,  -67.,   36.],       [  56.,  -52.,   31.],       [  71.,  -26.,   34.],       [  78.,   -5.,   28.],       [  81.,   28.,   29.],       [  55.,   58.,   30.],       [  30.,   73.,   30.],       [  -1.,   88.,   34.],       [ -33.,   81.,   30.],       [ -63.,   55.,   36.],       [ -73.,   35.,   31.],       [ -83.,    3.,   34.],       [ -77.,  -37.,   35.],       [ -62.,  -62.,   30.],       [ -31.,  -80.,   29.],       [   1.,  -87.,   29.],       [  31.,  -78.,   33.],       [  61.,  -59.,   27.],       [  80.,  -35.,   29.],       [  96.,   -2.,   33.],       [  87.,   33.,   30.],       [  65.,   67.,   36.],       [  31.,   86.,   31.],       [   3.,   96.,   32.],       [ -35.,   89.,   28.],       [ -62.,   67.,   35.],       [ -83.,   32.,   31.],       [ -88.,    2.,   27.],       [ -80.,  -34.,   33.],       [ -62.,  -66.,   30.],       [ -31.,  -90.,   30.],       [   3.,  -87.,   30.],       [  38.,  -84.,   28.],       [  67.,  -64.,   32.],       [  82.,  -37.,   28.],       [  98.,   -5.,   30.],       [  93.,   34.,   29.],       [  73.,   71.,   33.],       [  36.,   95.,   30.],       [   0.,   95.,   28.],       [ -34.,   94.,   36.],       [ -75.,   74.,   35.],       [ -97.,   41.,   35.],       [ -97.,    1.,   29.],       [ -95.,  -34.,   32.],       [ -71.,  -70.,   35.],       [ -41.,  -93.,   28.],       [  -1.,  -97.,   36.],       [  33.,  -90.,   35.],       [  65.,  -68.,   32.],       [  90.,  -39.,   35.],       [  46.,    0.,   39.],       [  44.,   21.,   41.],       [  33.,   37.,   37.],       [  23.,   47.,   39.],       [   2.,   54.,   44.],       [ -17.,   42.,   40.],       [ -40.,   35.,   35.],       [ -51.,   14.,   40.],       [ -51.,    4.,   42.],       [ -42.,  -23.,   35.],       [ -33.,  -39.,   43.],       [ -21.,  -45.,   35.],       [  -1.,  -55.,   36.],       [  21.,  -43.,   37.],       [  33.,  -38.,   37.],       [  44.,  -20.,   41.],       [  59.,    0.,   36.],       [  55.,   26.,   40.],       [  45.,   43.,   42.],       [  23.,   50.,   40.],       [  -4.,   55.,   39.],       [ -24.,   56.,   41.],       [ -39.,   37.,   35.],       [ -58.,   27.,   44.],       [ -60.,   -5.,   42.],       [ -52.,  -22.,   37.],       [ -38.,  -42.,   41.],       [ -23.,  -50.,   40.],       [  -4.,  -59.,   41.],       [  22.,  -54.,   40.],       [  38.,  -46.,   43.],       [  50.,  -20.,   35.],       [  68.,   -4.,   40.],       [  60.,   20.,   39.],       [  50.,   47.,   41.],       [  28.,   65.,   35.],       [  -2.,   68.,   40.],       [ -29.,   57.,   37.],       [ -44.,   46.,   40.],       [ -66.,   24.,   35.],       [ -63.,   -1.,   41.],       [ -64.,  -31.,   44.],       [ -45.,  -50.,   36.],       [ -26.,  -59.,   43.],       [  -2.,  -69.,   36.],       [  28.,  -61.,   36.],       [  51.,  -44.,   35.],       [  63.,  -26.,   35.],       [  79.,    3.,   42.],       [  65.,   33.,   37.],       [  51.,   50.,   42.],       [  24.,   69.,   40.],       [  -3.,   78.,   36.],       [ -28.,   64.,   39.],       [ -51.,   49.,   41.],       [ -74.,   26.,   35.],       [ -71.,   -2.,   41.],       [ -68.,  -26.,   36.],       [ -56.,  -52.,   42.],       [ -27.,  -73.,   37.],       [  -4.,  -80.,   44.],       [  29.,  -70.,   44.],       [  50.,  -57.,   41.],       [  67.,  -28.,   40.],       [  87.,    3.,   38.],       [  74.,   32.,   37.],       [  63.,   58.,   39.],       [  32.,   75.,   43.],       [   2.,   79.,   37.],       [ -30.,   75.,   42.],       [ -61.,   63.,   42.],       [ -80.,   32.,   37.],       [ -86.,   -4.,   39.],       [ -74.,  -31.,   38.],       [ -57.,  -55.,   42.],       [ -35.,  -77.,   37.],       [  -5.,  -84.,   43.],       [  29.,  -73.,   36.],       [  58.,  -56.,   42.],       [  76.,  -29.,   38.],       [  90.,    2.,   40.],       [  81.,   31.,   43.],       [  60.,   60.,   44.],       [  33.,   84.,   39.],       [   3.,   93.,   39.],       [ -32.,   83.,   43.],       [ -65.,   61.,   43.],       [ -83.,   31.,   44.],       [ -95.,   -3.,   41.],       [ -81.,  -38.,   39.],       [ -70.,  -62.,   37.],       [ -35.,  -81.,   44.],       [  -5.,  -95.,   35.],       [  30.,  -81.,   39.],       [  64.,  -66.,   43.],       [  85.,  -39.,   35.],       [  96.,   -4.,   35.],       [  94.,   34.,   38.],       [  73.,   73.,   36.],       [  40.,   95.,   39.],       [   4.,  101.,   39.],       [ -42.,   92.,   41.],       [ -70.,   72.,   39.],       [ -91.,   42.,   40.],       [-101.,    1.,   44.],       [ -95.,  -39.,   38.],       [ -70.,  -69.,   42.],       [ -39.,  -96.,   38.],       [   2.,  -98.,   39.],       [  36.,  -92.,   37.],       [  69.,  -74.,   38.],       [  95.,  -37.,   41.]])
        basis0 = BSplineBasis(7, np.array([ -2.4,  -1.4,   0. ,   0. ,   0. ,   0. ,   0. ,   0.6,   2.2,   2.8,   4.1,   4.6,   6. ,
         7.3,   8.1,   8.8,   9.6,  10.6,  12. ,  12. ,  12. ,  12. ,  12. ,  12.6,  14.2]),1)
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.1,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        basis2 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0.8,  2. ,  2. ,  2. ,  2. ,  2. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 6)
        self.assertEqual(vol2.order(direction=1), 6)
        self.assertEqual(vol2.order(direction=2), 5)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 7)
        self.assertEqual(vol3.order(direction=1), 5)
        self.assertEqual(vol3.order(direction=2), 5)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 7)
        self.assertEqual(vol4.order(direction=1), 6)
        self.assertEqual(vol4.order(direction=2), 4)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)
    def test_curve_3D_p6_C2_periodic(self):
        controlpoints = np.array([[ 102.,    0.,   -5.],
       [  89.,   33.,   -5.],
       [  75.,   66.,    4.],
       [  40.,   91.,    3.],
       [  -4.,  103.,   -4.],
       [ -40.,   91.,   -4.],
       [ -68.,   73.,   -1.],
       [ -90.,   37.,    0.],
       [-105.,    1.,    4.],
       [ -93.,  -40.,    1.],
       [ -75.,  -69.,   -5.],
       [ -35.,  -94.,    1.],
       [   4., -105.,    0.],
       [  35.,  -95.,    3.],
       [  70.,  -73.,   -3.],
       [  91.,  -40.,   -1.]])
        basis0 = BSplineBasis(6, np.array([ -2.8,  -2.4,  -1.4,   0. ,   0. ,   0. ,   1. ,   2.1,   3.1,   4.3,   5.4,   6. ,   7.3,
         7.7,   9.1,  10.3,  11.2,  11.6,  12.6,  14. ,  14. ,  14. ,  15. ,  16.1,  17.1]),2)
        crv  = Curve(basis0, controlpoints,False)
        crv2 = crv.clone()
        crv2 = crv.get_derivative_spline(0)
        self.assertEqual(crv2.order(direction=0), 5)

        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
    def test_surface_3D_p56_C2_periodic(self):
        controlpoints = np.array([[  63.,   -2.,    0.],       [  56.,   26.,   -2.],       [  35.,   47.,   -3.],       [  16.,   59.,   -1.],       [ -13.,   62.,    2.],       [ -35.,   47.,    2.],       [ -50.,   21.,   -2.],       [ -58.,    0.,   -2.],       [ -54.,  -31.,   -3.],       [ -40.,  -46.,    4.],       [ -14.,  -61.,    4.],       [  13.,  -61.,   -3.],       [  34.,  -52.,   -2.],       [  56.,  -26.,    0.],       [  64.,    4.,   -5.],       [  55.,   26.,   -2.],       [  44.,   56.,   -3.],       [  11.,   69.,   -1.],       [ -15.,   69.,   -2.],       [ -43.,   48.,   -3.],       [ -57.,   29.,   -3.],       [ -72.,    3.,   -5.],       [ -62.,  -33.,   -5.],       [ -44.,  -51.,    4.],       [ -13.,  -63.,    4.],       [  14.,  -61.,    1.],       [  39.,  -52.,   -2.],       [  62.,  -32.,   -2.],       [  78.,   -1.,   -2.],       [  67.,   29.,   -2.],       [  45.,   56.,    4.],       [  18.,   74.,   -2.],       [ -12.,   72.,   -3.],       [ -49.,   59.,   -2.],       [ -62.,   35.,   -4.],       [ -76.,    3.,   -5.],       [ -65.,  -29.,    1.],       [ -42.,  -60.,   -2.],       [ -17.,  -69.,   -2.],       [  16.,  -77.,    2.],       [  41.,  -57.,    2.],       [  70.,  -33.,   -2.],       [  80.,   -2.,    4.],       [  71.,   34.,    1.],       [  49.,   67.,    0.],       [  19.,   76.,   -4.],       [ -21.,   79.,   -3.],       [ -53.,   61.,    3.],       [ -75.,   38.,    3.],       [ -82.,   -3.,   -5.],       [ -68.,  -38.,    0.],       [ -54.,  -59.,    0.],       [ -18.,  -80.,   -3.],       [  21.,  -76.,   -2.],       [  50.,  -61.,   -1.],       [  69.,  -33.,   -5.],       [  86.,   -5.,   -4.],       [  82.,   39.,    0.],       [  55.,   64.,    1.],       [  20.,   80.,   -4.],       [ -19.,   81.,   -2.],       [ -52.,   65.,   -2.],       [ -81.,   32.,   -2.],       [ -91.,   -3.,   -2.],       [ -82.,  -39.,    0.],       [ -59.,  -68.,    1.],       [ -21.,  -82.,    0.],       [  19.,  -90.,    0.],       [  55.,  -68.,    2.],       [  77.,  -34.,   -1.],       [  95.,   -4.,    3.],       [  87.,   40.,   -4.],       [  59.,   75.,   -5.],       [  19.,   94.,   -4.],       [ -24.,   95.,   -3.],       [ -61.,   74.,    2.],       [ -84.,   44.,    1.],       [ -95.,    1.,   -5.],       [ -82.,  -42.,   -1.],       [ -61.,  -76.,   -2.],       [ -18.,  -93.,    0.],       [  22.,  -91.,    0.],       [  54.,  -76.,    1.],       [  88.,  -42.,   -4.],       [ 103.,   -1.,   -1.],       [  87.,   47.,   -2.],       [  66.,   77.,    1.],       [  26.,   95.,   -1.],       [ -24.,   99.,    4.],       [ -64.,   80.,    2.],       [ -87.,   45.,    0.],       [-105.,   -2.,    3.],       [ -86.,  -40.,    0.],       [ -65.,  -81.,    2.],       [ -21.,  -97.,    4.],       [  24.,  -95.,    3.],       [  58.,  -82.,    4.],       [  89.,  -46.,   -3.]])
        basis0 = BSplineBasis(5, np.array([ -3.2,  -1.8,  -0.6,   0. ,   0. ,   1.4,   2.2,   3.2,   3.8,   4.8,   5.9,   7. ,   7.8,
         8.7,   9.8,  11.2,  12.4,  13. ,  13. ,  14.4,  15.2,  16.2]),2)
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.8,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        surf  = Surface(basis0, basis1, controlpoints,False)
        surf2 = surf.clone()
        surf2 = surf.get_derivative_spline(0)
        self.assertEqual(surf2.order(direction=0), 4)
        self.assertEqual(surf2.order(direction=1), 6)
        surf3 = surf.get_derivative_spline(1)
        self.assertEqual(surf3.order(direction=0), 5)
        self.assertEqual(surf3.order(direction=1), 5)

        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
    def test_volume_3D_p565_C2_periodic(self):
        controlpoints = np.array([[  54.,    2.,   -1.],       [  41.,   18.,    2.],       [  29.,   37.,   -5.],       [   6.,   49.,   -5.],       [ -12.,   46.,   -2.],       [ -32.,   40.,   -2.],       [ -44.,   20.,    4.],       [ -51.,   -1.,    1.],       [ -43.,  -26.,   -2.],       [ -36.,  -39.,   -4.],       [  -9.,  -47.,    4.],       [   7.,  -45.,   -2.],       [  32.,  -37.,    1.],       [  40.,  -22.,   -4.],       [  60.,   -1.,   -5.],       [  50.,   23.,   -2.],       [  35.,   48.,    2.],       [  13.,   52.,   -3.],       [  -9.,   56.,   -1.],       [ -35.,   45.,    2.],       [ -55.,   22.,   -5.],       [ -58.,    0.,    1.],       [ -49.,  -27.,   -5.],       [ -34.,  -43.,    1.],       [ -12.,  -58.,    1.],       [  13.,  -60.,   -1.],       [  36.,  -43.,    3.],       [  48.,  -21.,   -1.],       [  64.,   -5.,   -2.],       [  56.,   29.,    0.],       [  41.,   48.,   -3.],       [  19.,   66.,   -4.],       [ -16.,   62.,   -5.],       [ -40.,   53.,   -3.],       [ -64.,   33.,   -3.],       [ -67.,    4.,    4.],       [ -65.,  -29.,    4.],       [ -46.,  -54.,    1.],       [ -19.,  -65.,    4.],       [  15.,  -69.,    3.],       [  42.,  -51.,   -5.],       [  62.,  -27.,    4.],       [  72.,   -2.,   -2.],       [  64.,   32.,   -1.],       [  47.,   54.,   -2.],       [  12.,   77.,    3.],       [ -14.,   75.,    2.],       [ -43.,   57.,   -4.],       [ -72.,   33.,    1.],       [ -79.,   -2.,    1.],       [ -66.,  -36.,   -4.],       [ -50.,  -61.,    3.],       [ -21.,  -74.,   -3.],       [  12.,  -69.,    2.],       [  49.,  -59.,   -1.],       [  65.,  -33.,   -5.],       [  81.,   -1.,    1.],       [  73.,   37.,    2.],       [  52.,   63.,   -4.],       [  17.,   77.,   -1.],       [ -21.,   77.,   -3.],       [ -56.,   66.,   -1.],       [ -79.,   40.,    1.],       [ -79.,    2.,   -1.],       [ -75.,  -40.,   -2.],       [ -54.,  -69.,    4.],       [ -16.,  -81.,   -1.],       [  16.,  -81.,    0.],       [  50.,  -70.,    3.],       [  70.,  -35.,   -3.],       [  93.,    2.,   -1.],       [  79.,   42.,   -1.],       [  59.,   74.,    0.],       [  19.,   89.,   -2.],       [ -24.,   86.,   -3.],       [ -60.,   73.,   -2.],       [ -88.,   38.,    4.],       [ -90.,   -2.,   -2.],       [ -86.,  -40.,    3.],       [ -61.,  -76.,    2.],       [ -24.,  -90.,   -3.],       [  23.,  -91.,    1.],       [  58.,  -70.,    3.],       [  78.,  -44.,   -2.],       [ 101.,   -5.,   -5.],       [  88.,   39.,   -1.],       [  61.,   74.,   -3.],       [  25.,   93.,    1.],       [ -18.,   99.,   -4.],       [ -67.,   77.,   -2.],       [ -90.,   45.,    4.],       [ -96.,    4.,   -3.],       [ -95.,  -40.,   -3.],       [ -61.,  -78.,   -3.],       [ -20.,  -97.,    1.],       [  24., -101.,   -5.],       [  64.,  -76.,    1.],       [  93.,  -45.,   -4.],       [  53.,    1.,    5.],       [  47.,   24.,    9.],       [  30.,   40.,    3.],       [  12.,   44.,    8.],       [  -7.,   51.,    6.],       [ -29.,   35.,    7.],       [ -48.,   23.,    3.],       [ -49.,    1.,    8.],       [ -43.,  -27.,    3.],       [ -32.,  -37.,    3.],       [ -16.,  -51.,    4.],       [  13.,  -48.,   11.],       [  28.,  -44.,    4.],       [  47.,  -21.,    9.],       [  59.,    2.,    8.],       [  51.,   27.,    4.],       [  32.,   49.,    6.],       [  17.,   58.,    8.],       [ -13.,   54.,    8.],       [ -42.,   45.,    3.],       [ -58.,   24.,    5.],       [ -62.,    0.,   11.],       [ -54.,  -26.,    9.],       [ -41.,  -48.,    5.],       [ -11.,  -56.,    5.],       [  16.,  -60.,    8.],       [  36.,  -45.,    3.],       [  57.,  -25.,    5.],       [  63.,    1.,   11.],       [  64.,   24.,    2.],       [  45.,   55.,    6.],       [  16.,   66.,    2.],       [ -20.,   66.,    6.],       [ -42.,   47.,    3.],       [ -60.,   27.,    6.],       [ -70.,    1.,    9.],       [ -60.,  -32.,    2.],       [ -40.,  -52.,    4.],       [ -19.,  -66.,    4.],       [  16.,  -66.,    3.],       [  39.,  -57.,    9.],       [  59.,  -30.,    4.],       [  72.,    4.,    7.],       [  68.,   29.,    7.],       [  45.,   62.,    7.],       [  17.,   72.,    2.],       [ -12.,   73.,    8.],       [ -48.,   54.,    9.],       [ -63.,   32.,    2.],       [ -79.,    1.,    6.],       [ -69.,  -28.,    5.],       [ -46.,  -62.,    8.],       [ -14.,  -74.,    6.],       [  20.,  -77.,    2.],       [  51.,  -62.,    5.],       [  69.,  -34.,    3.],       [  85.,    0.,    6.],       [  76.,   35.,    2.],       [  55.,   66.,    3.],       [  18.,   81.,   10.],       [ -18.,   82.,    2.],       [ -51.,   67.,    4.],       [ -78.,   33.,   11.],       [ -83.,    2.,    4.],       [ -78.,  -39.,    9.],       [ -51.,  -66.,    6.],       [ -21.,  -84.,    8.],       [  21.,  -83.,   11.],       [  49.,  -62.,    7.],       [  70.,  -38.,   11.],       [  94.,    2.,    5.],       [  80.,   44.,    3.],       [  60.,   69.,   11.],       [  16.,   87.,    9.],       [ -20.,   85.,    5.],       [ -62.,   75.,    2.],       [ -82.,   41.,    4.],       [ -90.,   -1.,   10.],       [ -88.,  -38.,    6.],       [ -60.,  -74.,    5.],       [ -21.,  -88.,    6.],       [  21.,  -92.,   10.],       [  59.,  -75.,   11.],       [  79.,  -43.,    2.],       [  96.,    4.,    8.],       [  87.,   38.,   10.],       [  65.,   83.,    1.],       [  19.,   96.,    5.],       [ -23.,   99.,    6.],       [ -63.,   78.,    2.],       [ -86.,   39.,    5.],       [-104.,    2.,    7.],       [ -86.,  -46.,    3.],       [ -59.,  -74.,    5.],       [ -27.,  -96.,    4.],       [  19.,  -98.,    3.],       [  59.,  -75.,    6.],       [  88.,  -48.,   10.],       [  52.,   -5.,   12.],       [  44.,   16.,   13.],       [  33.,   38.,    8.],       [   6.,   47.,   15.],       [ -11.,   49.,   17.],       [ -29.,   39.,   13.],       [ -43.,   17.,    8.],       [ -50.,   -2.,   17.],       [ -45.,  -24.,   13.],       [ -31.,  -38.,   15.],       [ -15.,  -48.,   17.],       [  15.,  -50.,   12.],       [  28.,  -35.,   11.],       [  42.,  -20.,   17.],       [  58.,    1.,   15.],       [  56.,   23.,   18.],       [  38.,   45.,   17.],       [  16.,   54.,   10.],       [  -9.,   52.,   10.],       [ -33.,   42.,   13.],       [ -50.,   26.,   15.],       [ -62.,    2.,   13.],       [ -48.,  -25.,   17.],       [ -37.,  -47.,   15.],       [  -9.,  -61.,   17.],       [  16.,  -56.,   10.],       [  40.,  -43.,   16.],       [  50.,  -22.,   14.],       [  68.,    0.,   12.],       [  62.,   30.,   10.],       [  42.,   51.,   18.],       [  15.,   65.,   10.],       [ -16.,   66.,   17.],       [ -44.,   50.,   16.],       [ -64.,   30.,    9.],       [ -67.,    1.,   10.],       [ -64.,  -34.,   15.],       [ -41.,  -53.,   10.],       [ -12.,  -70.,   15.],       [   9.,  -69.,   14.],       [  36.,  -54.,   10.],       [  60.,  -27.,   16.],       [  79.,    3.,   15.],       [  65.,   33.,   10.],       [  43.,   61.,   13.],       [  15.,   73.,   13.],       [ -13.,   69.,   11.],       [ -45.,   59.,   16.],       [ -64.,   28.,   18.],       [ -77.,   -3.,   12.],       [ -66.,  -31.,   14.],       [ -45.,  -57.,   13.],       [ -15.,  -70.,   14.],       [  20.,  -70.,   14.],       [  43.,  -59.,   15.],       [  64.,  -35.,   15.],       [  79.,    3.,   13.],       [  71.,   37.,    9.],       [  52.,   68.,   14.],       [  13.,   80.,   11.],       [ -18.,   82.,   11.],       [ -52.,   65.,   17.],       [ -80.,   38.,   17.],       [ -79.,    2.,    9.],       [ -77.,  -39.,   11.],       [ -57.,  -67.,   12.],       [ -21.,  -82.,   11.],       [  14.,  -82.,   10.],       [  52.,  -65.,   12.],       [  70.,  -40.,   13.],       [  91.,    4.,   10.],       [  78.,   38.,   15.],       [  58.,   74.,    8.],       [  20.,   89.,   11.],       [ -16.,   84.,   15.],       [ -53.,   75.,   17.],       [ -85.,   39.,   10.],       [ -88.,    0.,    8.],       [ -81.,  -45.,   10.],       [ -58.,  -72.,   11.],       [ -17.,  -92.,    9.],       [  16.,  -87.,   15.],       [  56.,  -67.,   14.],       [  81.,  -44.,   11.],       [  97.,   -1.,   16.],       [  88.,   45.,    8.],       [  59.,   79.,   13.],       [  18.,   95.,   13.],       [ -26.,   94.,   10.],       [ -61.,   77.,   14.],       [ -89.,   42.,   16.],       [ -98.,   -5.,   12.],       [ -86.,  -48.,   15.],       [ -65.,  -82.,   15.],       [ -27.,  -98.,   16.],       [  23.,  -94.,   14.],       [  61.,  -79.,   13.],       [  88.,  -49.,   13.],       [  54.,   -1.,   23.],       [  49.,   16.,   17.],       [  33.,   42.,   17.],       [  12.,   45.,   24.],       [  -8.,   53.,   15.],       [ -34.,   35.,   16.],       [ -49.,   21.,   19.],       [ -54.,    0.,   22.],       [ -42.,  -19.,   20.],       [ -32.,  -44.,   16.],       [ -14.,  -48.,   20.],       [   7.,  -45.,   24.],       [  32.,  -38.,   21.],       [  43.,  -23.,   23.],       [  60.,    4.,   19.],       [  47.,   29.,   22.],       [  36.,   50.,   20.],       [   8.,   58.,   22.],       [ -12.,   55.,   15.],       [ -40.,   42.,   22.],       [ -56.,   22.,   15.],       [ -55.,   -3.,   23.],       [ -51.,  -28.,   18.],       [ -37.,  -43.,   21.],       [ -16.,  -61.,   24.],       [  11.,  -58.,   22.],       [  39.,  -50.,   21.],       [  51.,  -29.,   21.],       [  70.,    3.,   21.],       [  56.,   29.,   17.],       [  46.,   48.,   21.],       [  17.,   63.,   18.],       [ -18.,   64.,   24.],       [ -42.,   56.,   24.],       [ -64.,   26.,   15.],       [ -71.,    4.,   15.],       [ -62.,  -33.,   19.],       [ -39.,  -51.,   17.],       [ -11.,  -64.,   22.],       [  12.,  -67.,   18.],       [  41.,  -51.,   16.],       [  63.,  -31.,   20.],       [  79.,    1.,   19.],       [  67.,   37.,   22.],       [  43.,   63.,   16.],       [  19.,   77.,   16.],       [ -18.,   77.,   19.],       [ -44.,   59.,   16.],       [ -70.,   31.,   21.],       [ -75.,   -3.,   22.],       [ -71.,  -32.,   24.],       [ -43.,  -55.,   20.],       [ -20.,  -72.,   17.],       [  17.,  -74.,   23.],       [  42.,  -61.,   15.],       [  70.,  -28.,   20.],       [  83.,    3.,   15.],       [  79.,   37.,   22.],       [  56.,   62.,   15.],       [  19.,   80.,   16.],       [ -17.,   78.,   23.],       [ -52.,   67.,   21.],       [ -73.,   31.,   22.],       [ -89.,   -2.,   15.],       [ -79.,  -34.,   17.],       [ -54.,  -68.,   15.],       [ -15.,  -83.,   24.],       [  18.,  -79.,   24.],       [  55.,  -67.,   16.],       [  71.,  -35.,   15.],       [  93.,   -1.,   22.],       [  79.,   39.,   15.],       [  56.,   70.,   23.],       [  16.,   90.,   18.],       [ -23.,   88.,   16.],       [ -59.,   67.,   21.],       [ -79.,   43.,   18.],       [ -91.,   -2.,   16.],       [ -79.,  -42.,   21.],       [ -63.,  -67.,   20.],       [ -21.,  -88.,   23.],       [  18.,  -91.,   16.],       [  61.,  -76.,   18.],       [  86.,  -44.,   24.],       [  96.,    1.,   15.],       [  88.,   38.,   16.],       [  58.,   77.,   20.],       [  25.,   99.,   24.],       [ -26.,   95.,   24.],       [ -63.,   81.,   20.],       [ -90.,   45.,   19.],       [-105.,   -1.,   21.],       [ -89.,  -48.,   19.],       [ -62.,  -75.,   23.],       [ -25.,  -99.,   15.],       [  26., -100.,   23.],       [  61.,  -83.,   24.],       [  94.,  -40.,   15.],       [  50.,   -3.,   31.],       [  49.,   24.,   23.],       [  26.,   37.,   25.],       [   7.,   50.,   27.],       [  -7.,   51.,   26.],       [ -35.,   42.,   26.],       [ -42.,   23.,   24.],       [ -53.,    3.,   22.],       [ -41.,  -27.,   29.],       [ -34.,  -43.,   30.],       [ -12.,  -52.,   26.],       [  12.,  -47.,   25.],       [  29.,  -37.,   27.],       [  41.,  -17.,   30.],       [  60.,    1.,   28.],       [  52.,   24.,   29.],       [  31.,   42.,   23.],       [   8.,   57.,   23.],       [ -12.,   58.,   28.],       [ -34.,   44.,   29.],       [ -54.,   28.,   30.],       [ -62.,   -3.,   25.],       [ -58.,  -22.,   25.],       [ -35.,  -46.,   24.],       [ -17.,  -57.,   26.],       [  12.,  -58.,   29.],       [  39.,  -46.,   29.],       [  54.,  -24.,   30.],       [  62.,    4.,   23.],       [  62.,   24.,   26.],       [  38.,   56.,   27.],       [  14.,   67.,   24.],       [ -12.,   63.,   22.],       [ -44.,   50.,   31.],       [ -58.,   28.,   23.],       [ -71.,    2.,   26.],       [ -62.,  -29.,   26.],       [ -42.,  -48.,   25.],       [ -14.,  -64.,   25.],       [  11.,  -61.,   23.],       [  46.,  -55.,   27.],       [  64.,  -26.,   26.],       [  73.,   -5.,   28.],       [  70.,   35.,   23.],       [  51.,   55.,   31.],       [  20.,   73.,   24.],       [ -14.,   73.,   31.],       [ -47.,   61.,   30.],       [ -66.,   35.,   28.],       [ -76.,   -5.,   25.],       [ -72.,  -32.,   23.],       [ -51.,  -57.,   26.],       [ -14.,  -70.,   26.],       [  13.,  -70.,   25.],       [  46.,  -58.,   23.],       [  65.,  -31.,   25.],       [  82.,    1.,   30.],       [  76.,   37.,   22.],       [  54.,   63.,   27.],       [  14.,   78.,   26.],       [ -15.,   79.,   22.],       [ -56.,   69.,   25.],       [ -78.,   37.,   27.],       [ -86.,    0.,   30.],       [ -72.,  -42.,   29.],       [ -51.,  -61.,   24.],       [ -17.,  -82.,   25.],       [  14.,  -79.,   24.],       [  52.,  -67.,   31.],       [  73.,  -39.,   27.],       [  92.,    1.,   26.],       [  86.,   41.,   25.],       [  59.,   69.,   31.],       [  22.,   92.,   27.],       [ -21.,   85.,   29.],       [ -55.,   72.,   23.],       [ -81.,   43.,   23.],       [ -94.,   -2.,   25.],       [ -86.,  -39.,   24.],       [ -55.,  -72.,   25.],       [ -17.,  -90.,   22.],       [  15.,  -88.,   23.],       [  56.,  -75.,   29.],       [  87.,  -36.,   21.],       [ 100.,   -3.,   30.],       [  89.,   45.,   31.],       [  58.,   80.,   27.],       [  26.,   94.,   24.],       [ -22.,   93.,   29.],       [ -66.,   77.,   25.],       [ -93.,   41.,   23.],       [-102.,    4.,   27.],       [ -90.,  -39.,   30.],       [ -59.,  -77.,   23.],       [ -25., -102.,   27.],       [  26., -103.,   24.],       [  65.,  -76.,   28.],       [  93.,  -43.,   23.],       [  53.,    1.,   31.],       [  45.,   20.,   33.],       [  33.,   34.,   33.],       [  11.,   46.,   31.],       [ -15.,   53.,   28.],       [ -28.,   40.,   35.],       [ -45.,   23.,   30.],       [ -52.,   -3.,   33.],       [ -50.,  -22.,   29.],       [ -34.,  -42.,   31.],       [ -14.,  -52.,   30.],       [  14.,  -54.,   28.],       [  34.,  -42.,   37.],       [  47.,  -20.,   34.],       [  55.,    2.,   31.],       [  51.,   25.,   30.],       [  32.,   45.,   29.],       [   9.,   54.,   28.],       [  -9.,   54.,   33.],       [ -38.,   49.,   30.],       [ -53.,   25.,   36.],       [ -59.,    3.,   33.],       [ -54.,  -23.,   29.],       [ -37.,  -42.,   31.],       [ -12.,  -58.,   36.],       [  14.,  -60.,   35.],       [  38.,  -47.,   36.],       [  52.,  -22.,   32.],       [  66.,   -5.,   31.],       [  60.,   31.,   34.],       [  36.,   51.,   31.],       [  14.,   69.,   32.],       [ -14.,   67.,   34.],       [ -42.,   55.,   29.],       [ -60.,   28.,   35.],       [ -69.,   -2.,   37.],       [ -64.,  -28.,   29.],       [ -44.,  -56.,   29.],       [ -20.,  -70.,   33.],       [  13.,  -66.,   34.],       [  39.,  -50.,   33.],       [  64.,  -25.,   30.],       [  74.,    2.,   33.],       [  67.,   33.,   30.],       [  46.,   58.,   32.],       [  19.,   76.,   30.],       [ -14.,   69.,   30.],       [ -51.,   61.,   29.],       [ -68.,   28.,   30.],       [ -78.,   -4.,   31.],       [ -66.,  -36.,   30.],       [ -51.,  -62.,   29.],       [ -14.,  -72.,   29.],       [  13.,  -71.,   31.],       [  44.,  -59.,   32.],       [  70.,  -35.,   37.],       [  84.,    0.,   30.],       [  77.,   41.,   34.],       [  52.,   62.,   32.],       [  14.,   77.,   34.],       [ -17.,   80.,   31.],       [ -52.,   68.,   32.],       [ -76.,   33.,   36.],       [ -83.,   -5.,   32.],       [ -75.,  -38.,   30.],       [ -56.,  -65.,   34.],       [ -14.,  -80.,   32.],       [  17.,  -79.,   37.],       [  53.,  -70.,   37.],       [  73.,  -37.,   37.],       [  89.,   -5.,   35.],       [  79.,   35.,   30.],       [  57.,   72.,   29.],       [  23.,   94.,   34.],       [ -24.,   86.,   29.],       [ -59.,   75.,   35.],       [ -84.,   40.,   30.],       [ -92.,    3.,   33.],       [ -80.,  -45.,   35.],       [ -59.,  -74.,   30.],       [ -17.,  -87.,   30.],       [  18.,  -91.,   32.],       [  58.,  -76.,   37.],       [  84.,  -37.,   29.],       [ 101.,    0.,   31.],       [  89.,   47.,   33.],       [  58.,   82.,   33.],       [  21.,   93.,   35.],       [ -24.,   94.,   34.],       [ -62.,   77.,   36.],       [ -90.,   47.,   28.],       [ -98.,    4.,   36.],       [ -95.,  -43.,   36.],       [ -65.,  -83.,   33.],       [ -26.,  -93.,   32.],       [  24.,  -93.,   33.],       [  61.,  -74.,   37.],       [  94.,  -46.,   35.],       [  54.,    2.,   35.],       [  48.,   26.,   43.],       [  30.,   40.,   41.],       [  14.,   50.,   35.],       [ -13.,   45.,   37.],       [ -32.,   43.,   37.],       [ -45.,   22.,   41.],       [ -52.,   -1.,   35.],       [ -44.,  -24.,   36.],       [ -32.,  -39.,   40.],       [ -15.,  -47.,   37.],       [  14.,  -49.,   42.],       [  29.,  -42.,   44.],       [  42.,  -18.,   39.],       [  59.,    3.,   37.],       [  47.,   29.,   42.],       [  38.,   49.,   44.],       [   8.,   53.,   43.],       [ -12.,   61.,   39.],       [ -38.,   45.,   43.],       [ -51.,   28.,   44.],       [ -63.,    2.,   44.],       [ -53.,  -26.,   39.],       [ -36.,  -47.,   43.],       [ -17.,  -61.,   39.],       [  10.,  -53.,   40.],       [  38.,  -50.,   37.],       [  54.,  -22.,   44.],       [  61.,    3.,   41.],       [  63.,   25.,   36.],       [  36.,   48.,   39.],       [  16.,   65.,   40.],       [ -12.,   62.,   38.],       [ -40.,   52.,   41.],       [ -57.,   33.,   38.],       [ -68.,    2.,   38.],       [ -58.,  -28.,   43.],       [ -38.,  -51.,   36.],       [ -13.,  -64.,   41.],       [  11.,  -64.,   42.],       [  42.,  -50.,   36.],       [  55.,  -32.,   36.],       [  71.,   -5.,   37.],       [  64.,   31.,   35.],       [  46.,   54.,   38.],       [  20.,   68.,   42.],       [ -12.,   68.,   35.],       [ -52.,   55.,   43.],       [ -64.,   35.,   44.],       [ -77.,    4.,   44.],       [ -67.,  -35.,   41.],       [ -47.,  -64.,   37.],       [ -22.,  -76.,   39.],       [  20.,  -79.,   37.],       [  48.,  -61.,   44.],       [  68.,  -29.,   36.],       [  84.,    3.,   35.],       [  77.,   39.,   40.],       [  54.,   60.,   43.],       [  17.,   82.,   41.],       [ -23.,   83.,   41.],       [ -53.,   60.,   36.],       [ -77.,   36.,   35.],       [ -84.,   -1.,   38.],       [ -73.,  -40.,   42.],       [ -49.,  -66.,   42.],       [ -17.,  -82.,   36.],       [  21.,  -86.,   42.],       [  49.,  -64.,   40.],       [  78.,  -39.,   39.],       [  93.,    1.,   37.],       [  86.,   36.,   41.],       [  58.,   70.,   38.],       [  15.,   86.,   37.],       [ -24.,   91.,   43.],       [ -62.,   72.,   39.],       [ -81.,   42.,   38.],       [ -95.,    0.,   38.],       [ -80.,  -37.,   35.],       [ -59.,  -71.,   40.],       [ -18.,  -91.,   41.],       [  21.,  -95.,   42.],       [  60.,  -76.,   35.],       [  80.,  -37.,   41.],       [  98.,    2.,   36.],       [  94.,   47.,   42.],       [  64.,   82.,   41.],       [  18.,   92.,   44.],       [ -19.,  102.,   42.],       [ -67.,   76.,   39.],       [ -95.,   43.,   43.],       [-101.,    0.,   35.],       [ -86.,  -42.,   44.],       [ -58.,  -82.,   43.],       [ -26., -101.,   43.],       [  22., -102.,   35.],       [  66.,  -84.,   37.],       [  90.,  -49.,   36.]])
        basis0 = BSplineBasis(5, np.array([ -2.9,  -2.1,  -0.8,   0. ,   0. ,   1.2,   1.8,   2.6,   4.4,   5. ,   5.9,   7.4,   8. ,
         9.2,  10.1,  10.9,  12.2,  13. ,  13. ,  14.2,  14.8,  15.6]),2)
        basis1 = BSplineBasis(6, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.2,  2. ,  2. ,  2. ,  2. ,  2. ,  2. ]))
        basis2 = BSplineBasis(5, np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  1.2,  1.9,  3. ,  3. ,  3. ,  3. ,  3. ]))
        vol  = Volume(basis0, basis1, basis2, controlpoints,False)
        vol2 = vol.clone()
        vol2 = vol.get_derivative_spline(0)
        self.assertEqual(vol2.order(direction=0), 4)
        self.assertEqual(vol2.order(direction=1), 6)
        self.assertEqual(vol2.order(direction=2), 5)
        vol3 = vol.get_derivative_spline(1)
        self.assertEqual(vol3.order(direction=0), 5)
        self.assertEqual(vol3.order(direction=1), 5)
        self.assertEqual(vol3.order(direction=2), 5)
        vol4 = vol.get_derivative_spline(2)
        self.assertEqual(vol4.order(direction=0), 5)
        self.assertEqual(vol4.order(direction=1), 6)
        self.assertEqual(vol4.order(direction=2), 4)

        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)

if __name__ == '__main__':
    unittest.main()

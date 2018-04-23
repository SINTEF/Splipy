# -*- coding: utf-8 -*-

from splipy.io import G2
import splipy.surface_factory as SurfaceFactory
import splipy.volume_factory as VolumeFactory
import numpy as np
import unittest
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestG2(unittest.TestCase):

    def test_read_rational_surf(self):
        with G2(THIS_DIR + '/test_geometries/torus.g2') as myfile:
            surf = myfile.read()
        self.assertEqual(len(surf), 1)
        surf = surf[0]
        self.assertTrue(surf.rational)


    def test_write_and_read_multipatch_surface(self):
        # write teapot to file and test if its there
        teapot = SurfaceFactory.teapot()
        with G2('teapot.g2') as myfile:
            myfile.write(teapot)
        self.assertTrue(os.path.isfile('teapot.g2'))

        # read the written file and compare that it's the same thing
        with G2('teapot.g2') as myfile:
            read_teapot = myfile.read()
        self.assertEqual(len(read_teapot), 32)
        for i in range(32):
            self.assertTrue(np.all(teapot[i].shape == read_teapot[i].shape))

        # clean up after us
        os.remove('teapot.g2')

    def test_read_doublespaced(self):
        with G2(THIS_DIR + '/test_geometries/lshape.g2') as myfile: # controlpoints are separated by two spaces
            one_surf = myfile.read()
        self.assertEqual(len(one_surf), 1)
        self.assertEqual(one_surf[0].shape[0], 3)
        self.assertEqual(one_surf[0].shape[1], 2)

    def test_write_and_read_surface(self):
        # write disc to file and test if its there
        disc = SurfaceFactory.disc(type='square')
        with G2('disc.g2') as myfile:
            myfile.write(disc)
        self.assertTrue(os.path.isfile('disc.g2'))

        # read the written file and compare that it's the same thing
        with G2('disc.g2') as myfile:
            read_disc = myfile.read()
        self.assertEqual(len(read_disc), 1)
        read_disc = read_disc[0]
        self.assertTrue(np.all(disc.shape == read_disc.shape))

        # clean up after us
        os.remove('disc.g2')


    def test_write_and_read_volume(self):
        # write sphere to file and test if its there
        sphere = VolumeFactory.sphere(type='square')
        with G2('sphere.g2') as myfile:
            myfile.write(sphere)
        self.assertTrue(os.path.isfile('sphere.g2'))

        # read the written file and compare that it's the same thing
        with G2('sphere.g2') as myfile:
            read_sphere = myfile.read()
        self.assertEqual(len(read_sphere), 1)
        read_sphere = read_sphere[0]
        self.assertTrue(np.all(sphere.shape == read_sphere.shape))
        self.assertTrue(sphere.rational == read_sphere.rational)

        # clean up after us
        os.remove('sphere.g2')

    def test_read_elementary_curves(self):
        with G2(THIS_DIR + '/test_geometries/elementary_curves.g2') as myfile:
            my_curves = myfile.read()

        self.assertEqual(len(my_curves), 3)

        # check circle (r=3, center=(1,0,0), xaxis=(1,1,0)
        circle = my_curves[0]
        t = np.linspace(circle.start(0), circle.end(0), 25)
        x = circle(t)
        self.assertTrue(np.allclose((x[:,0]-1)**2 + x[:,1]**2 + x[:,2]**2, 3**2))
        self.assertTrue(np.allclose(circle[0], [3/np.sqrt(2)+1,3/np.sqrt(2),0,1]))

        # check ellipse (r1=3, r2=5, center=(1,0,0), xaxis(0,1,0)
        ellipse = my_curves[1]
        t = np.linspace(ellipse.start(0), ellipse.end(0), 25)
        x = ellipse(t)
        self.assertTrue(np.allclose(ellipse[0], [1,3,0,1]))
        self.assertTrue(np.allclose(ellipse[2], [-4,0,0,1]))

        # check line piece (p1=(1,0,0), direction=(0,6,0), length=4 (times direction))
        line = my_curves[2]
        self.assertAlmostEqual(line.length(), 24)
        self.assertFalse(line.rational)
        self.assertTrue(np.allclose(line[0], [1,0,0]))

    def test_read_elementary_surfaces(self):
        with G2(THIS_DIR + '/test_geometries/elementary_surfaces.g2') as myfile:
            my_surfaces = myfile.read()

        self.assertEqual(len(my_surfaces), 4)

        # check cylinder ( center=(1,0,0), x-axis=(0,1,0), radius=2)
        cylinder = my_surfaces[0]
        self.assertTrue(cylinder.rational)
        self.assertTrue(np.allclose(cylinder[0,0], [1,2,0,1]))

        # check sphere ( radius=1.5, center=(4,0,0), z-axis=(0,1,1))
        sphere = my_surfaces[1]
        self.assertTrue(sphere.rational)
        self.assertTrue(np.allclose(sphere[0,-1], [4,3/np.sqrt(8),3/np.sqrt(8), 1])) # north pole

        # check disc ( radius=2.5, center=(6,0,0), x-axis=(0,0,-1 ), normal=(1,0,0)
        disc = my_surfaces[2]
        self.assertTrue(disc.rational)
        self.assertTrue(np.allclose(disc[ 0,0], [6,0,0,1]))    # center
        self.assertTrue(np.allclose(disc[-1,0], [6,0,-2.5,1])) # center+x_axis

        # check torus ( ...)
        torus = my_surfaces[3]
        self.assertTrue(torus.rational)

    def test_from_step(self):
        # quite large nasty g2 file which contains cylinders, planes, trimming etc

        with G2(THIS_DIR + '/test_geometries/winglet_from_step.g2') as myfile:
            my_surfaces = myfile.read()
            trim_curves = myfile.trimming_curves

        # we only test that we are able to parse this file. No claims as to the
        # correctness of this parsing
        self.assertEqual(len(my_surfaces), 19)
        self.assertEqual(len(trim_curves), 122)

if __name__ == '__main__':
    unittest.main()

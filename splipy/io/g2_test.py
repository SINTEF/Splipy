# -*- coding: utf-8 -*-

from splipy import *
from splipy.io import *
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
        t = np.linspace(circle.start(), circle.end(), 25)
        x = circle(t)
        self.assertTrue(np.allclose((x[:,0]-1)**2 + x[:,1]**2 + x[:,2]**2, 3**2))
        self.assertTrue(np.allclose(circle[0], [3/sqrt(2)+1,3/sqrt(2),0,1]))

        # check ellipse (r1=3, r2=5, center=(1,0,0), xaxis(0,1,0)
        ellipse = my_curves[1]
        t = np.linspace(ellipse.start(), ellipse.end(), 25)
        x = ellipse(t)
        self.assertTrue(np.allclose(ellipse[0], [1,3,0,1]))
        self.assertTrue(np.allclose(ellipse[2], [-4,0,0,1]))

        # check line piece (p1=(1,0,0), direction=(6,1,0), length=4)
        line = my_curves[2]
        self.assertAlmostEqual(line.length(), 4)
        self.assertFalse(line.rational)
        self.assertTrue(np.allclose(line[0], [1,0,0]))



if __name__ == '__main__':
    unittest.main()

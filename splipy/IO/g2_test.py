# -*- coding: utf-8 -*-

from splipy import *
from splipy.IO import *
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


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-

from operator import itemgetter
from splipy import Volume
from splipy.SplineModel import SplineModel, Orientation
import unittest
import numpy as np


class TestOrientation(unittest.TestCase):

    def test_identical(self):
        v1 = Volume()
        v2 = Volume()
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (0, 1, 2))
        self.assertEqual(ori.flip, (False, False, False))

    def test_permute(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap('u', 'v')
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (1, 0, 2))
        self.assertEqual(ori.flip, (False, False, False))

        v2.swap('v', 'w')
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (2, 0, 1))
        self.assertEqual(ori.flip, (False, False, False))

    def test_flip(self):
        v1 = Volume()
        v2 = Volume()

        v2.reverse('u')
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (0, 1, 2))
        self.assertEqual(ori.flip, (True, False, False))

        v2.reverse('w')
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (0, 1, 2))
        self.assertEqual(ori.flip, (True, False, True))

    def test_permute_and_flip(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap('u', 'v')
        v2.reverse('v')
        v2.swap('u', 'w')
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (1, 2, 0))
        self.assertEqual(ori.flip, (True, False, False))

    def test_map_section(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap('u', 'v')
        v2.reverse('v')
        v2.swap('u', 'w')
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.map_section((-1, 0, -1)), (-1, -1, -1))
        self.assertEqual(ori.map_section((0, 0, None)), (-1, None, 0))

    def test_map_array(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap('u', 'v')
        v2.reverse('v')
        v2.swap('u', 'w')
        ori = Orientation.compute(v1, v2)
        self.assertTrue((
            ori.map_array(np.reshape(np.arange(8, dtype=int), (2,2,2)))
            == np.reshape([2, 6, 3, 7, 0, 4, 1, 5], (2,2,2))
        ).all())


class TestModel(unittest.TestCase):

    def test_node_counts(self):
        model = SplineModel(3,3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1,0,0))
        model.add(v + (0,0,1))

        # 3 volumes
        cat = model.catalogue
        self.assertEqual(len(cat.top_nodes()), 3)

        # 16 faces
        cat = cat.lower
        self.assertEqual(len(cat.top_nodes()), 16)

        # 28 edges
        cat = cat.lower
        self.assertEqual(len(cat.top_nodes()), 28)

        # 16 vertices
        cat = cat.lower
        self.assertEqual(len(cat.top_nodes()), 16)

    def test_lookup(self):
        model = SplineModel(3,3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1,0,0))
        model.add(v + (0,0,1))

        q = Volume().refine(1, 1, 1)
        self.assertIs(model[q].node, model[v].node)

        q = Volume().refine(1, 1, 1).reverse('u').swap('u', 'v')
        self.assertIs(model[q].node, model[v].node)

        self.assertIs(model[v.section(w=-1)].node, model[(v+(0,0,1)).section(w=0)].node)

    def test_cp_numbering(self):
        model = SplineModel(3,3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1,0,0))
        model.add(v + (0,0,1))

        # Predictable numbering of control points
        model.generate_cp_numbers()
        cps = model.cps()
        self.assertEqual(cps.shape, (63, 3))
        np.testing.assert_almost_equal(cps[:27], v.controlpoints.reshape(-1,3))
        np.testing.assert_almost_equal(cps[27:45], v.controlpoints[1:].reshape(-1,3) + (1,0,0))
        np.testing.assert_almost_equal(cps[45:], v.controlpoints[:,:,1:].reshape(-1,3) + (0,0,1))

    def test_cell_numbering(self):
        model = SplineModel(3,3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1,0,0))
        model.add(v + (0,0,1))

        # Predictable numbering of cells
        model.generate_cell_numbers()
        self.assertTrue((model[v].node.cell_numbers.flatten() == np.arange(8, dtype=int)).all())
        self.assertTrue((model[v+(1,0,0)].node.cell_numbers.flatten() == np.arange(8, 16, dtype=int)).all())
        self.assertTrue((model[v+(0,0,1)].node.cell_numbers.flatten() == np.arange(16, 24, dtype=int)).all())

    def test_faces(self):
        model = SplineModel(3,3)
        v = Volume().refine(1, 1, 1)
        vr = v + (1,0,0)
        vu = v + (0,0,1)
        model.add([v, vr, vu])

        model[v.section(u=0)].name = 'uneg'
        model[vr.section(u=-1)].name = 'upos'
        model[vu.section(u=0)].name = 'uneg'
        model[vu.section(u=-1)].name = 'upos'

        for vv in [v, vr, vu]:
            model[vv.section(v=0)].name = 'vneg'
            model[vv.section(v=-1)].name = 'vpos'

        model[v.section(w=0)].name = 'wneg'
        model[vr.section(w=0)].name = 'wneg'
        model[vr.section(w=-1)].name = 'wpos'
        model[vu.section(w=-1)].name = 'wpos'

        model.generate_cp_numbers()
        model.generate_cell_numbers()
        faces = model.faces()

        self.assertEqual(len(faces), 100)

        # The owner should always be the lower numbered cell
        self.assertTrue(((faces['owner'] < faces['neighbor']) | (faces['neighbor'] == -1)).all())

        # Where the boundary is named, the neighbor is always -1
        self.assertTrue((faces[faces['name'] != None]['neighbor'] == -1).all())

        # On internal faces, the neighbor is always > 0
        self.assertTrue((faces[faces['name'] == None]['neighbor'] > 0).all())

        # Boundaries occur the expected number of times
        self.assertEqual(np.sum(faces['name'] == None), 44)
        self.assertEqual(np.sum(faces['name'] == 'uneg'), 8)
        self.assertEqual(np.sum(faces['name'] == 'upos'), 8)
        self.assertEqual(np.sum(faces['name'] == 'vneg'), 12)
        self.assertEqual(np.sum(faces['name'] == 'vpos'), 12)
        self.assertEqual(np.sum(faces['name'] == 'wneg'), 8)
        self.assertEqual(np.sum(faces['name'] == 'wpos'), 8)

        # Every cell is mentioned by exactly six faces
        nmentions = [0] * 24
        for face in faces:
            nmentions[face['owner']] += 1
            if face['neighbor'] > 0:
                nmentions[face['neighbor']] += 1
        self.assertEqual(nmentions, [6] * 24)

        # Face normals always point toward the neighbor
        # In this case, that means in a positive axial direction,
        # except possibly for boundary faces
        cps = model.cps()
        for face in faces:
            for I in [face['nodes'][:-1], face['nodes'][1:]]:
                vertices = cps[I]
                z = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[1])

                # Should be zero, zero and +/- 0.25
                self.assertAlmostEqual(sum(np.abs(z)), 0.25)

                # j = index of the 0.25
                j = np.argmax(np.abs(z))
                others = {0: [1, 2], 1: [0, 2], 2: [0, 1]}[j]
                np.testing.assert_almost_equal(z[others], 0.0)

                # Check that the 0.25 has the right sign
                if face['name'] is None or face['name'][1:] == 'pos':
                    np.testing.assert_almost_equal(z[j], 0.25)
                else:
                    np.testing.assert_almost_equal(z[j], -0.25)

                # And that it occurs in the expected index (for boundary faces)
                if face['name'] is not None:
                    self.assertEqual(j, 'uvw'.index(face['name'][0]))

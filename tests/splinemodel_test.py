# -*- coding: utf-8 -*-

from operator import itemgetter
from splipy import Volume, Surface
from splipy.splinemodel import SplineModel, Orientation, IFEMWriter, IFEMConnection
from splipy.io import G2
from splipy import curve_factory, surface_factory, volume_factory
import unittest
import numpy as np

import os


THIS_DIR = os.path.dirname(os.path.abspath(__file__))


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

        v2.swap("u", "v")
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (1, 0, 2))
        self.assertEqual(ori.flip, (False, False, False))

        v2.swap("v", "w")
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (2, 0, 1))
        self.assertEqual(ori.flip, (False, False, False))

    def test_flip(self):
        v1 = Volume()
        v2 = Volume()

        v2.reverse("u")
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (0, 1, 2))
        self.assertEqual(ori.flip, (True, False, False))

        v2.reverse("w")
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (0, 1, 2))
        self.assertEqual(ori.flip, (True, False, True))

    def test_permute_and_flip(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap("u", "v")
        v2.reverse("v")
        v2.swap("u", "w")
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.perm, (1, 2, 0))
        self.assertEqual(ori.flip, (True, False, False))

    def test_map_section(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap("u", "v")
        v2.reverse("v")
        v2.swap("u", "w")
        ori = Orientation.compute(v1, v2)
        self.assertEqual(ori.map_section((-1, 0, -1)), (-1, -1, -1))
        self.assertEqual(ori.map_section((0, 0, None)), (-1, None, 0))

    def test_map_array(self):
        v1 = Volume()
        v2 = Volume()

        v2.swap("u", "v")
        v2.reverse("v")
        v2.swap("u", "w")
        ori = Orientation.compute(v1, v2)
        self.assertTrue(
            (
                ori.map_array(np.reshape(np.arange(8, dtype=int), (2, 2, 2)))
                == np.reshape([2, 6, 3, 7, 0, 4, 1, 5], (2, 2, 2))
            ).all()
        )


class TestModel(unittest.TestCase):
    def test_node_counts(self):
        model = SplineModel(3, 3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1, 0, 0))
        model.add(v + (0, 0, 1))

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
        model = SplineModel(3, 3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1, 0, 0))
        model.add(v + (0, 0, 1))

        q = Volume().refine(1, 1, 1)
        self.assertIs(model[q].node, model[v].node)

        q = Volume().refine(1, 1, 1).reverse("u").swap("u", "v")
        self.assertIs(model[q].node, model[v].node)

        self.assertIs(model[v.section(w=-1)].node, model[(v + (0, 0, 1)).section(w=0)].node)

    def test_cp_numbering(self):
        model = SplineModel(3, 3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1, 0, 0))
        model.add(v + (0, 0, 1))

        # Predictable numbering of control points
        model.generate_cp_numbers()
        cps = model.cps()
        self.assertEqual(cps.shape, (63, 3))
        np.testing.assert_almost_equal(cps[:27], v.controlpoints.reshape(-1, 3))
        np.testing.assert_almost_equal(cps[27:45], v.controlpoints[1:].reshape(-1, 3) + (1, 0, 0))
        np.testing.assert_almost_equal(cps[45:], v.controlpoints[:, :, 1:].reshape(-1, 3) + (0, 0, 1))

    def test_cell_numbering(self):
        model = SplineModel(3, 3)
        v = Volume().refine(1, 1, 1)
        model.add(v)
        model.add(v + (1, 0, 0))
        model.add(v + (0, 0, 1))

        # Predictable numbering of cells
        model.generate_cell_numbers()
        self.assertTrue((model[v].node.cell_numbers.flatten() == np.arange(8, dtype=int)).all())
        self.assertTrue(
            (model[v + (1, 0, 0)].node.cell_numbers.flatten() == np.arange(8, 16, dtype=int)).all()
        )
        self.assertTrue(
            (model[v + (0, 0, 1)].node.cell_numbers.flatten() == np.arange(16, 24, dtype=int)).all()
        )

    def test_faces(self):
        model = SplineModel(3, 3)
        v = Volume().refine(1, 1, 1)
        vr = v + (1, 0, 0)
        vu = v + (0, 0, 1)
        model.add([v, vr, vu])

        model[v.section(u=0)].name = "uneg"
        model[vr.section(u=-1)].name = "upos"
        model[vu.section(u=0)].name = "uneg"
        model[vu.section(u=-1)].name = "upos"

        for vv in [v, vr, vu]:
            model[vv.section(v=0)].name = "vneg"
            model[vv.section(v=-1)].name = "vpos"

        model[v.section(w=0)].name = "wneg"
        model[vr.section(w=0)].name = "wneg"
        model[vr.section(w=-1)].name = "wpos"
        model[vu.section(w=-1)].name = "wpos"

        model.generate_cp_numbers()
        model.generate_cell_numbers()
        faces = model.faces()

        self.assertEqual(len(faces), 100)

        # The owner should always be the lower numbered cell
        self.assertTrue(((faces["owner"] < faces["neighbor"]) | (faces["neighbor"] == -1)).all())

        # Where the boundary is named, the neighbor is always -1
        self.assertTrue((faces[faces["name"] != None]["neighbor"] == -1).all())

        # On internal faces, the neighbor is always > 0
        self.assertTrue((faces[faces["name"] == None]["neighbor"] > 0).all())

        # Boundaries occur the expected number of times
        self.assertEqual(np.sum(faces["name"] == None), 44)
        self.assertEqual(np.sum(faces["name"] == "uneg"), 8)
        self.assertEqual(np.sum(faces["name"] == "upos"), 8)
        self.assertEqual(np.sum(faces["name"] == "vneg"), 12)
        self.assertEqual(np.sum(faces["name"] == "vpos"), 12)
        self.assertEqual(np.sum(faces["name"] == "wneg"), 8)
        self.assertEqual(np.sum(faces["name"] == "wpos"), 8)

        # Every cell is mentioned by exactly six faces
        nmentions = [0] * 24
        for face in faces:
            nmentions[face["owner"]] += 1
            if face["neighbor"] > 0:
                nmentions[face["neighbor"]] += 1
        self.assertEqual(nmentions, [6] * 24)

        # Face normals always point toward the neighbor
        # In this case, that means in a positive axial direction,
        # except possibly for boundary faces
        cps = model.cps()
        for face in faces:
            for I in [face["nodes"][:-1], face["nodes"][1:]]:
                vertices = cps[I]
                z = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[1])

                # Should be zero, zero and +/- 0.25
                self.assertAlmostEqual(sum(np.abs(z)), 0.25)

                # j = index of the 0.25
                j = np.argmax(np.abs(z))
                others = {0: [1, 2], 1: [0, 2], 2: [0, 1]}[j]
                np.testing.assert_almost_equal(z[others], 0.0)

                # Check that the 0.25 has the right sign
                if face["name"] is None or face["name"][1:] == "pos":
                    np.testing.assert_almost_equal(z[j], 0.25)
                else:
                    np.testing.assert_almost_equal(z[j], -0.25)

                # And that it occurs in the expected index (for boundary faces)
                if face["name"] is not None:
                    self.assertEqual(j, "uvw".index(face["name"][0]))

    def test_orient(self):
        connections = [
            [
                IFEMConnection(1, 2, 2, 1, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 4, 3, 0),
                IFEMConnection(2, 6, 6, 5, 0),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 2, 1, 3),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 3, 3, 1),
                IFEMConnection(2, 6, 5, 5, 1),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 2, 2, 2),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 3, 3, 2),
                IFEMConnection(2, 6, 6, 5, 3),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 2, 2, 1),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 4, 3, 3),
                IFEMConnection(2, 6, 5, 5, 2),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 2, 3, 1),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 2, 3, 1),
                IFEMConnection(2, 6, 5, 5, 4),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 2, 5, 5),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 3, 3, 4),
                IFEMConnection(2, 6, 2, 5, 5),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 2, 4, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 6, 5, 0),
                IFEMConnection(2, 4, 2, 3, 2),
                IFEMConnection(2, 6, 6, 5, 6),
                IFEMConnection(3, 4, 2, 1, 0),
                IFEMConnection(3, 7, 6, 5, 0),
                IFEMConnection(4, 8, 6, 5, 0),
                IFEMConnection(5, 6, 2, 1, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 2, 1, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 1, 1),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 6, 4),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 4, 1, 1),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 5, 4),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 1, 0),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 4, 1, 4),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 6, 0),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 3, 1),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 2, 1, 6),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 3, 5),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 1, 5),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 5, 1, 4),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 1, 7),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 4, 5),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 5, 1, 1),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 1, 2),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 5, 5),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 3, 1, 0),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 5, 7),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 2, 2),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 3, 1, 5),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 3, 6),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 2, 4),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 6, 1, 5),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
            [
                IFEMConnection(1, 2, 6, 5, 0),
                IFEMConnection(1, 3, 4, 3, 0),
                IFEMConnection(1, 5, 2, 1, 0),
                IFEMConnection(2, 4, 4, 4, 7),
                IFEMConnection(2, 6, 2, 1, 0),
                IFEMConnection(3, 4, 6, 2, 7),
                IFEMConnection(3, 7, 2, 1, 0),
                IFEMConnection(4, 8, 5, 1, 7),
                IFEMConnection(5, 7, 4, 3, 0),
                IFEMConnection(5, 6, 6, 5, 0),
                IFEMConnection(6, 8, 4, 3, 0),
                IFEMConnection(7, 8, 6, 5, 0),
            ],
        ]

        for i, ref_topo in enumerate(connections):
            model = SplineModel(3, 3)
            with G2(THIS_DIR + "/geometries/cube-8-orient{}.g2".format(i)) as myfile:
                model.add(myfile.read())

            writer = IFEMWriter(model)
            my_topo = list(writer.connections())

            my_topo = sorted(my_topo, key=lambda c: c.slave)
            my_topo = sorted(my_topo, key=lambda c: c.master)
            ref_topo = sorted(ref_topo, key=lambda c: c.slave)
            ref_topo = sorted(ref_topo, key=lambda c: c.master)

            for my_con, ref_con in zip(my_topo, ref_topo):
                assert my_con == ref_con
            assert len(my_topo) == len(ref_topo)

    def test_2d_self_connection(self):
        c = curve_factory.line([1, 0], [2, 0])
        surf = surface_factory.revolve(c)
        surf = surf.split(surf.knots("v")[0], direction="v")  # break periodicity
        surf.set_dimension(2)
        model = SplineModel(2, 2)
        model.add(surf)

        writer = IFEMWriter(model)
        expected = [IFEMConnection(1, 1, 3, 4, 0)]
        for connection, want in zip(writer.connections(), expected):
            self.assertEqual(connection, want)

    def test_3d_self_connection(self):
        square = Surface() + [1, 0]
        square = square.rotate(np.pi / 2, (1, 0, 0))
        vol = volume_factory.revolve(square)
        vol = vol.split(vol.knots("w")[0], direction="w")  # break periodicity
        model = SplineModel(3, 3)
        model.add(vol, raise_on_twins=False)

        writer = IFEMWriter(model)
        expected = [IFEMConnection(1, 1, 5, 6, 0)]
        for connection, want in zip(writer.connections(), expected):
            self.assertEqual(connection, want)

    def test_3d_torus_self_connection(self):
        torus = volume_factory.torus()
        torus = torus.split(torus.knots(1)[0], direction="v")
        torus = torus.split(torus.knots(2)[0], direction="w")
        model = SplineModel(3, 3)
        model.add(torus, raise_on_twins=False)

        writer = IFEMWriter(model)
        expected = [IFEMConnection(1, 1, 3, 4, 0), IFEMConnection(1, 1, 5, 6, 0)]
        for connection, want in zip(writer.connections(), expected):
            self.assertEqual(connection, want)

    def test_3d_self_double_connection(self):
        c1 = curve_factory.circle(r=1, center=(3, 0, 0), normal=(0, 1, 0))
        c2 = curve_factory.circle(r=2, center=(3, 0, 0), normal=(0, 1, 0))
        ring = surface_factory.edge_curves(c1, c2)
        vol = volume_factory.revolve(ring)
        vol = vol.split(vol.knots("u")[0], direction="u")  # break periodicity
        vol = vol.split(vol.knots("w")[0], direction="w")
        model = SplineModel(3, 3)
        model.add(vol, raise_on_twins=False)

        writer = IFEMWriter(model)
        expected = [IFEMConnection(1, 1, 1, 2, 0), IFEMConnection(1, 1, 5, 6, 0)]
        for connection, want in zip(writer.connections(), expected):
            self.assertEqual(connection, want)

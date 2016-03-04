# -*- coding: utf-8 -*-

from splipy import *
from splipy.utils import *
import splipy.state as state
import numpy as np
from collections import Counter
from itertools import chain, product, permutations

try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping


class VertexDict(MutableMapping):
    """A dictionary where the keys are numpy arrays, and where equality
    computed in an a pproximate sense for floating point numbers.

    All keys must have the same dimensions.
    """

    def __init__(self):
        # List of (key, value) pairs
        self.internal = []

    def _eq(self, a, b):
        """Check whether two numpy arrays are almost equal, according to the given
        tolerances.
        """
        return np.allclose(a, b,
                           rtol=state.controlpoint_relative_tolerance,
                           atol=state.controlpoint_absolute_tolerance)

    def _candidate(self, key):
        """Return the internal index for the first stored mapping that matches the
        given key.

        :param numpy.array key: The key to look for
        :raises KeyError: If the key is not found
        """
        for i, (k, v) in enumerate(self.internal):
            if self._eq(key, k):
                return i
        raise KeyError(key)

    def __setitem__(self, key, value):
        """Assign a key to a value."""
        try:
            c = self._candidate(key)
            k, _ = self.internal[c]
            self.internal[c] = (k, value)
        except KeyError:
            self.internal.append((key, value))

    def __getitem__(self, key):
        """Gets the value assigned to a key.

        :raises KeyError: If the key is not found
        """
        c = self._candidate(key)
        return self.internal[c][1]

    def __delitem__(self, key):
        """Deletes an assignment."""
        self.internal = [(k, v) for k, v in self.internal
                         if not self._eq(key, k)]

    def __iter__(self):
        """Iterate over all keys.

        .. note:: This generates all the stored keys, not all matching keys.
        """
        for k, _ in self.internal:
            yield k

    def items(self):
        """Return a list of key, value pairs."""
        return list(self.internal)

    def __len__(self):
        """Returns the number of stored assignments."""
        return len(self.internal)


class OrientationError(RuntimeError):
    """An `OrientationError` is raised by certain methods in
    :class:`splipy.SplineModel.Orientation` indicating an inability to match
    two objects.
    """
    pass

class Orientation(object):
    """An `Orientation` represents a mapping between two coordinate systems: the
    *reference* system and the *actual* or *mapped* system.

    The orientation is represented in terms of two tuples:

    - `perm` is a permutation: direction `d` in the reference system is mapped
      to direction `perm[d]` in the actual system.
    - `flip` is a tuple of bools: `flip[d]` will be `True` to indicate that
      direction `d` *in the reference system* should be reversed.
    """

    def __init__(self, perm, flip):
        """Initialize an Orientation object.

        .. warning:: This constructor is for internal use. Use
            :func:`splipy.SplineModel.Orientation.compute` to generate new
            orientation objects.
        """
        self.perm = perm
        self.flip = flip
        self.perm_inv = tuple(perm.index(d) for d in range(len(perm)))

    @classmethod
    def compute(cls, cpa, cpb=None):
        """Compute and return a new orientation object representing the mapping between
        `cpa` (the reference system) and `cpb` (the mapped system).

        If `cpb` is not given, the identity orientation is returned.

        :param SplineObject cpa: The reference system
        :param SplineObject cpb: The mapped system
        :return: The orientation relating the two systems
        :rtype: Orientation
        :raises OrientationError: If the two objects do not match
        """
        shape_a = cpa.controlpoints.shape
        pardim = len(shape_a) - 1

        # Return the identity orientation if no cpb
        if cpb is None:
            return cls(tuple(range(pardim)),
                       tuple(False for _ in range(pardim)))

        shape_b = cpb.controlpoints.shape

        # Deal with the easy cases: dimension mismatch, and
        # comparing the shapes as multisets
        if len(shape_a) != len(shape_b):
            raise OrientationError("Mismatching parametric dimensions")
        if shape_a[-1] != shape_b[-1]:
            raise OrientationError("Mismatching physical dimensions")
        if Counter(shape_a) != Counter(shape_b):
            raise OrientationError("Non-matching objects")

        # Enumerate all permutations of directions
        for perm in permutations(range(pardim)):
            transposed = cpb.controlpoints.transpose(perm + (pardim,))
            if transposed.shape != shape_a:
                continue
            # Enumerate all possible direction reversals
            for flip in product([False, True], repeat=pardim):
                slices = tuple(slice(None, None, -1) if f else slice(None) for f in flip)
                test_b = transposed[slices + (slice(None),)]
                if np.allclose(cpa.controlpoints, test_b,
                               rtol=state.controlpoint_relative_tolerance,
                               atol=state.controlpoint_absolute_tolerance):
                    if all([cpa.bases[i].matches(cpb.bases[perm[i]], reverse=flip[i]) for i in range(pardim)]):
                        return cls(perm, flip)

        raise OrientationError("Non-matching objects")

    @property
    def pardim(self):
        return len(self.perm)

    def __mul__(self, other):
        """Compose two mappings.

        If `ort_left` maps system `A` (reference) to system `B`, and
        `ort_right` maps system `B` (reference) to system `C`, then `ort_left *
        ort_right` maps system `A` (reference) to system `C`.
        """
        assert self.pardim == other.pardim
        perm = tuple(other.perm[self.perm[d]] for d in range(self.pardim))

        # Flip a direction if it is flipped in exactly one orientation
        flip = tuple(self.flip[d] != other.flip[self.perm[d]] for d in range(self.pardim))

        return Orientation(perm, flip)

    def map_section(self, section):
        """Map a section in the mapped system to the reference system.

        The input is a section tuple as described in
        :func:`splipy.Utils.section`. It should be described relative to the
        mapped coordinate system. The return value is the corresponding section
        relative to the reference system.
        """
        permuted = tuple(section[d] for d in self.perm)

        flipped = ()
        for s, f, in zip(permuted, self.flip):
            # Flipping only applies to indexed directions, not variable ones
            if f and s is not None:
                flipped += (0 if s == -1 else -1,)
            else:
                flipped += (s,)

        return flipped

    def view_section(self, section):
        """Reduce a mapping to a lower dimension.

        The input is a section tuple as described in
        :func:`splipy.Utils.section`. It should be described relative to the
        mapped coordinate system. The return value is an `Orientation` object
        that describes the same orentation, relative to the section only.

        For example, for a 2D orentation that flips only the second direction,
        `view_section` will return the identity orientation for edges along the
        first direction, but the inverse orientation for edges along the second
        (flipped) direction.
        """
        # The directions (in the mapped system) that are variable
        variable_dirs = [i for i, s in enumerate(section) if s is None]

        # The directions (in the reference system) that are variable
        # Sort to get the order that they have in the reference system
        actual_dirs = sorted(self.perm_inv[d] for d in variable_dirs)

        # perm[d] maps a direction to the mapped system, and its index in
        # variable_dirs reveals the permutation
        new_perm = tuple(variable_dirs.index(self.perm[d]) for d in actual_dirs)

        new_flip = tuple(self.flip[d] for d in actual_dirs)
        return self.__class__(new_perm, new_flip)


class TopologicalNode(object):
    """A `TopologicalNode` object refers to a single, persistent point in the
    topological graph. It represents some object of dimension `d` (that is, a
    point, an edge, etc.) and it has references to all the other objects it
    touches (that is, a face has a reference to four vertices, four edges, and
    any volumes it may be part of).

    Relevant attributes are:

    - obj: The underlying `SplineObject`
    - lower_nodes: A nested list of lower-order connections to other
      `TopologicalNode` objects. `lower_nodes[d]` is a list of all connected
      node objects with dimension `d`.
    - higher_nodes: A dictionary of higher-order connections to other
      `TopologicalNode` objects. `higher_nodes[d]` is a set of all connected
      node objects with dimension `d`.

    .. note:: Connections to lower order nodes are `ordered` corresponding to
        the natural ordering of sections (see :func:`splipy.Utils.sections`).
        Connections to higher order nodes are not ordered.

    .. warning:: This class is mostly for internal use. It performs no checks
        of any kind.
    """

    def __init__(self, obj, lower_nodes):
        """Initialize a `TopologicalNode` object associated with the given
        `SplineObject` and lower order nodes.

        :param SplineObject obj: The underlying spline object
        :param lower_nodes: A nested list of lower order nodes
        """
        self.obj = obj
        self.lower_nodes = lower_nodes
        self.higher_nodes = {}

        for dim_nodes in self.lower_nodes:
            for node in dim_nodes:
                node.assign_higher(self)

    @property
    def pardim(self):
        return self.obj.pardim

    def assign_higher(self, node):
        """Add a link to a node of higher dimension."""
        self.higher_nodes.setdefault(node.pardim, set()).add(node)

    def view(self, other_obj=None):
        """Return a `NodeView` object of this node.

        The returned view has an orientation that matches that of the input
        argument. If the argument is omitted, the identity orientation is used.

        :param SplineObject other_obj: The "view" object
        :raises OrientationError: If the object does not match this node's
            underlying object
        """
        if other_obj:
            orientation = Orientation.compute(self.obj, other_obj)
        else:
            orientation = Orientation.compute(self.obj)
        return NodeView(self, orientation)


class NodeView(object):
    """A `NodeView` object refers to a *view* to a point in the topological graph.
    It is composed of a node (:class:`splipy.SplineModel.TopologicalNode`) and
    an orientation (:class:`splipy.SplineModel.Orienation`).

    .. note:: Unlike `TopologicalNode` objects, `NodeView` objects are not
        persistent.
    """

    def __init__(self, node, orientation=None):
        """Initialize a `NodeView` object with the given node and orientation.

        .. warning:: This constructor is for internal use.
        """
        self.node = node
        self.orientation = orientation

    @property
    def pardim(self):
        return self.node.pardim

    def section(self, *args, **kwargs):
        """Return a section. See :func:`splipy.SplineObject.section` for more details
        on the input arguments.

        The return value is another `NodeView` object with a different
        underlying node, but whose orientation corresponds to the original
        view, as expected.
        """
        section = check_section(*args, pardim=self.pardim, **kwargs)
        tgt_dim = sum(1 for s in section if s is None)

        # The index of the section in the reference system
        ref_idx = section_to_index(self.orientation.map_section(section))

        # The underlying node
        node = self.node.lower_nodes[tgt_dim][ref_idx]

        # The underlying lower-order node may not have an orientation that
        # matches the higher-order node, so we need to compose two orientations
        ref_ori = Orientation.compute(node.obj, self.node.obj.section(*section, unwrap_points=False))
        my_ori = self.orientation.view_section(section)

        return NodeView(node, ref_ori * my_ori)

    def corner(self, i):
        """Return the i'th corner."""
        return self.section(*section_from_index(self.pardim, 0, i))

    @property
    def corners(self):
        """A tuple of all corners."""
        return tuple(self.section(s) for s in sections(self.pardim, 0))

    def edge(self, i):
        """Return the i'th edge."""
        return self.section(*section_from_index(self.pardim, 1, i))

    @property
    def edges(self):
        """A tuple of all edges."""
        return tuple(self.section(s) for s in sections(self.pardim, 1))

    def face(self, i):
        """Return the i'th face."""
        return self.section(*section_from_index(self.pardim, 2, i))

    @property
    def faces(self):
        """A tuple of all faces."""
        return tuple(self.section(s) for s in sections(self.pardim, 2))


class ObjectCatalogue(object):
    """An `ObjectCatalogue` maintains a complete topological graph of objects with
    at most `pardim` parametric directions.
    """

    def __init__(self, pardim):
        """Initialize a catalogue for objects of parametric dimension
        `pardim`.
        """
        self.pardim = pardim

        # Internal mapping from tuples of lower-order nodes to lists of nodes
        self.internal = {}

        # Each catalogue has a catalogue of lower dimension
        # For points, we use a VertexDict
        if pardim > 0:
            self.lower = ObjectCatalogue(pardim - 1)
        else:
            self.lower = VertexDict()

    def lookup(self, obj, add=False):
        """Obtain the `NodeView` object corresponding to a given object.

        If the keyword argument `add` is true, this function may generate one
        or more new nodes to accommodate the object.

        :param SplineObject obj: The object to look up
        :param bool add: Whether to allow adding new objects
        :return: A corresponding view
        :rtype: NodeView

        .. warning:: The object *must* be a `SplineObject`, even for points.
        """
        # Pass lower-dimensional objects through to the lower levels
        if self.pardim > obj.pardim:
            return self.lower.lookup(obj, add=add)

        # Special case for points: self.lower is a mapping from array to node
        if self.pardim == 0:
            if add:
                node = TopologicalNode(obj, [])
                return self.lower.setdefault(obj.controlpoints, node).view()
            return self.lower[obj.controlpoints].view()

        # Get all nodes of lower dimension (points, vertices, etc.)
        # This involves a recursive call to self.lower.__call__
        lower_nodes = []
        for i in range(0, self.pardim):
            nodes = tuple(self.lower.lookup(obj.section(*args, unwrap_points=False), add=add).node
                          for args in sections(self.pardim, i))
            lower_nodes.append(nodes)

        # Try looking up the lower-order nodes in the internal dictionary,
        # which maps tuples of nodes to lists of nodes. E.g. for volumes we
        # look up faces, for faces we look up edges, etc. Return the first one
        # we find for which the .view() function succeeds. This can throw a
        # KeyError (if this particular combination of lower-order nodes is new)
        # or an OrientationError (if it is not new, but the objects don't
        # match). If that happens, we generate a new node and return the
        # identity view on it.
        try:
            for candidate_node in self.internal[lower_nodes[-1]]:
                return candidate_node.view(obj)
        # FIXME: It might be useful to optionally not silence OrientationError,
        # since that more often than not indicates a real error
        except (KeyError, OrientationError):
            if not add:
                raise KeyError("No such object found")
            node = TopologicalNode(obj, lower_nodes)
            # Assign the new node to each possible permutation of lower-order
            # nodes. This is slight overkill since some of these permutations
            # are invalid, but c'est la vie.
            for p in permutations(lower_nodes[-1]):
                self.internal.setdefault(p, []).append(node)
            return node.view()

    def add(self, obj):
        """Add new nodes to the graph to accommodate the given object, then return the
        corresponding `NodeView` object.

        This is equivalent to calling
        :func:`splipy.SplineModel.ObjectCatalogue.lookup` with `add` set to
        true.

        :param SplineObject obj: The object to add
        :return: A corresponding view
        :rtype: NodeView

        .. warning:: The object *must* be a `SplineObject`, even for points.
        """
        return self.lookup(obj, add=True)

    __call__ = add

    def top_nodes(self):
        """Return all nodes of the highest parametric dimension."""
        return self.nodes(self.pardim)

    def nodes(self, pardim):
        """Return all nodes of a given parametric dimension."""
        if self.pardim == pardim:
            if self.pardim > 0:
                return set(chain.from_iterable(self.internal.values()))
            return set(self.lower.values())
        return self.lower.nodes(pardim)


# FIXME: This class is unfinished, and right now it doesn't do much other than
# wrap ObjectCatalogue

class SplineModel(object):

    def __init__(self, pardim=3, dimension=3, objs=[]):
        self.pardim = pardim
        self.dimension = dimension

        self.catalogue = ObjectCatalogue(pardim)
        self.add_patches(objs)

    def add_patch(self, obj):
        self.add_patches([obj])

    def add_patches(self, objs):
        self._validate(objs)
        self._generate(objs)

    def boundary(self):
        return [node for node in self.catalogue.nodes(self.pardim-1) if len(node.higher_nodes[self.pardim])==1]

    def _validate(self, objs):
        if any(p.dimension != self.dimension for p in objs):
            raise ValueError("Patches with different dimension added")
        if any(p.pardim != self.pardim for p in objs):
            raise ValueError("Patches with different parametric dimension added")

    def _generate(self, objs):
        for p in objs:
            self.catalogue(p)

    def summary(self):
        c = self.catalogue
        while isinstance(c, ObjectCatalogue):
            print('Dim {}: {}'.format(c.pardim, len(c.top_nodes())))
            c = c.lower

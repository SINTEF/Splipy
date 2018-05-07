# -*- coding: utf-8 -*-

from splipy import SplineObject
from splipy.utils import check_section, sections, section_from_index, uniquify
import splipy.state as state
import numpy as np
from collections import Counter, OrderedDict
from itertools import chain, product, permutations, islice

try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping


def _section_to_index(section):
    """Replace all `None` in `section` with `slice(None)`, so that it
    works as a numpy array indexing tuple.
    """
    return tuple(slice(None) if s is None else s for s in section)


face_t = np.dtype([('nodes', int, (4,)), ('owner', int, ()), ('neighbor', int, ()), ('name', object, ())])


class VertexDict(MutableMapping):
    """A dictionary where the keys are numpy arrays, and where equality
    is computed in an approximate sense for floating point numbers.

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

    def map_array(self, array):
        """Map an array in the mapped system to the reference system."""
        array = array.transpose(*self.perm)
        flips = tuple(slice(None, None, -1) if f else slice(None) for f in self.flip)
        return array[flips]

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
      `TopologicalNode` objects. `higher_nodes[d]` is a list of all connected
      node objects with dimension `d`.
    - owner: A top-level `TopologicalNode` object that "owns" this node.

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
        self.owner = None

        self.name = None
        self.cell_numbers = None
        self.cp_numbers = None

        for dim_nodes in self.lower_nodes:
            for node in dim_nodes:
                node.assign_higher(self)

        # Take ownership of lower nodes that are unaccounted for
        if self.pardim > 0:
            for node in lower_nodes[-1]:
                if node.owner is None:
                    node._transfer_ownership(self)

    @property
    def pardim(self):
        return self.obj.pardim

    @property
    def nhigher(self):
        return len(self.higher_nodes[self.pardim + 1])

    def assign_higher(self, node):
        """Add a link to a node of higher dimension."""
        self.higher_nodes.setdefault(node.pardim, list()).append(node)

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

    def _transfer_ownership(self, new_owner):
        """Transfers ownership of this node to a new owner. This operation is
        transitive, so all child nodes owned by this node, or who are
        owner-less will also be transferred.

        :param TopologicalNode new_owner: The new owner
        """
        self.owner = new_owner

        if self.pardim > 0:
            for child in self.lower_nodes[-1]:
                if child.owner is self or child.owner is None:
                    child._transfer_ownership(new_owner)

    def generate_cp_numbers(self, start=0):
        """Generate a control point numbering starting at `start`. Return the next unused index."""
        assert self.owner is None

        # Initialize a number array
        shape = self.obj.shape
        numbers = np.empty(shape, dtype=int)
        numbers[:] = 0

        # Flag control points owned by other top-level objects with -1
        for node, section in zip(self.lower_nodes[-1], sections(self.pardim, self.pardim-1)):
            if node.owner is not self:
                numbers[_section_to_index(section)] = -1

        # Fill in control point numbers for the ones we do own
        mask = np.where(numbers != -1)
        nowned = len(mask[0])
        numbers[mask] = np.arange(start, start + nowned, dtype=int)

        # This method takes care of communicating results to children
        self.assign_cp_numbers(numbers)
        return start + nowned

    def assign_cp_numbers(self, numbers):
        """Directly assign control point numbers."""
        self.cp_numbers = numbers

        # Control point numbers for owned children must be communicated to them
        if self.pardim > 0:
            for node, section in zip(self.lower_nodes[-1], sections(self.pardim, self.pardim-1)):
                if node.owner is self or node.owner is self.owner:
                    # Since this runs in a direct line of ownership, we don't need to be concerned with
                    # orientations not matching up.
                    node.assign_cp_numbers(numbers[_section_to_index(section)])

    def read_cp_numbers(self):
        """Read control point numbers for unowned control points from child nodes."""
        for node, section in zip(self.lower_nodes[-1], sections(self.pardim, self.pardim-1)):
            if node.owner is not self:
                # The two sections may not agree on orientation, so we fix this here.
                ori = Orientation.compute(self.obj.section(*section), node.obj)
                self.cp_numbers[_section_to_index(section)] = ori.map_array(node.cp_numbers)

        assert (self.cp_numbers != -1).all()

    def generate_cell_numbers(self, start=0):
        """Generate a cell numbering starting at `start`. Return the next unused index."""
        assert self.owner is None

        # Cells are never shared between top-level objects, so no need to worry about ownership here
        shape = [len(kvec) - 1 for kvec in self.obj.knots()]
        nelems = np.prod(shape)
        self.cell_numbers = np.reshape(np.arange(start, start + nelems, dtype=int), shape)
        return start + nelems

    def faces(self):
        """Return all faces owned by this node, as a list of numpy arrays with dtype `face_t`."""
        assert self.pardim == 3
        assert self.obj.order() == (2,2,2)
        shape = [len(kvec) - 1 for kvec in self.obj.knots()]
        ncells = np.prod(shape)
        retval = []

        def mkindex(dim, z, a, b):
            rval = [a, b] if dim != 1 else [b, a]
            rval.insert(dim, z)
            return tuple(rval)

        lower = iter(self.lower_nodes[-1])
        for d in range(self.pardim):
            # Number of faces in one "slice"
            nperslice = ncells // shape[d]

            # First, get all internal faces in this direction
            # The owner (lowest cell index) is guaranteed to be toward the lower end
            # TODO: We assume a right-hand coordinate system here
            nfaces = ncells - nperslice
            faces = np.empty((nfaces,), dtype=face_t)
            faces['nodes'][:,0] = self.cp_numbers[mkindex(d, np.s_[1:-1], np.s_[:-1], np.s_[:-1])].flatten()
            faces['nodes'][:,1] = self.cp_numbers[mkindex(d, np.s_[1:-1], np.s_[1:],  np.s_[:-1])].flatten()
            faces['nodes'][:,2] = self.cp_numbers[mkindex(d, np.s_[1:-1], np.s_[1:],  np.s_[1:])].flatten()
            faces['nodes'][:,3] = self.cp_numbers[mkindex(d, np.s_[1:-1], np.s_[:-1], np.s_[1:])].flatten()
            faces['owner'] = self.cell_numbers[mkindex(d, np.s_[:-1], np.s_[:], np.s_[:])].flatten()
            faces['neighbor'] = self.cell_numbers[mkindex(d, np.s_[1:],  np.s_[:], np.s_[:])].flatten()
            retval.append(faces)

            # Go through the two boundaries
            for bdnode, bdindex in zip(islice(lower, 2), (0, -1)):
                assert bdnode.nhigher in {1, 2}

                # Faces on an interface are only returned from the owner
                if bdnode.owner is not self:
                    continue

                faces = np.empty((nperslice,), dtype=face_t)
                faces['nodes'][:,0] = self.cp_numbers[mkindex(d, bdindex, np.s_[:-1], np.s_[:-1])].flatten()
                faces['nodes'][:,1] = self.cp_numbers[mkindex(d, bdindex, np.s_[1:], np.s_[:-1])].flatten()
                faces['nodes'][:,2] = self.cp_numbers[mkindex(d, bdindex, np.s_[1:], np.s_[1:])].flatten()
                faces['nodes'][:,3] = self.cp_numbers[mkindex(d, bdindex, np.s_[:-1], np.s_[1:])].flatten()
                faces['owner'] = self.cell_numbers[mkindex(d, bdindex, np.s_[:], np.s_[:])].flatten()
                faces['name'] = bdnode.name

                # If we're on the left boundary, the face normal must point in the other direction
                # NOTE: We copy when swapping here, since we are swapping values which are views into
                # a mutable array!
                if bdindex == 0:
                    faces['nodes'][:,1], faces['nodes'][:,3] = (
                        faces['nodes'][:,3].copy(), faces['nodes'][:,1].copy()
                    )

                # If there's a neighbor on the interface we need neighbouring cell numbers
                if bdnode.nhigher == 1:
                    faces['neighbor'] = -1
                else:
                    neighbor = next(c for c in bdnode.higher_nodes[3] if c is not self)

                    # Find out which face the interface is as numbered from the neighbor's perspective
                    nb_index = neighbor.lower_nodes[2].index(bdnode)

                    # Get the spline object on that interface as oriented from the neighbor's perspective
                    nb_sec = section_from_index(3, 2, nb_index)
                    nb_obj = neighbor.obj.section(*nb_sec)

                    # Compute the relative orientation
                    ori = Orientation.compute(bdnode.obj, nb_obj)

                    # Get the neighbor cell numbers from the neighbor's perspective, and map them to our system
                    cellidxs = neighbor.cell_numbers[_section_to_index(nb_sec)]
                    faces['neighbor'] = ori.map_array(cellidxs).flatten()

                retval.append(faces)

        for faces in retval:
            assert ((faces['owner'] < faces['neighbor']) | (faces['neighbor'] == -1)).all()
        return retval


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

    @property
    def name(self):
        return self.node.name

    @name.setter
    def name(self, value):
        self.node.name = value

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
        self.internal = OrderedDict()

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
    __getitem__ = lookup

    def top_nodes(self):
        """Return all nodes of the highest parametric dimension."""
        return self.nodes(self.pardim)

    def nodes(self, pardim):
        """Return all nodes of a given parametric dimension."""
        if self.pardim == pardim:
            if self.pardim > 0:
                return list(uniquify(chain.from_iterable(self.internal.values())))
            return list(uniquify(self.lower.values()))
        return self.lower.nodes(pardim)


# FIXME: This class is unfinished, and right now it doesn't do much other than
# wrap ObjectCatalogue

class SplineModel(object):

    def __init__(self, pardim=3, dimension=3, objs=[]):
        self.pardim = pardim
        self.dimension = dimension

        self.catalogue = ObjectCatalogue(pardim)
        self.add(objs)

    def add(self, obj):
        if isinstance(obj, SplineObject):
            obj = [obj]
        self._validate(obj)
        self._generate(obj)

    def __getitem__(self, obj):
        return self.catalogue[obj]

    def boundary(self, name=None):
        for node in self.catalogue.nodes(self.pardim-1):
            if node.nhigher == 1 and (name is None or name == node.name):
                yield node

    def assign_boundary(self, name):
        """Give a name to all unnamed boundary nodes."""
        for node in self.boundary():
            if node.name is None:
                node.name = name

    def _validate(self, objs):
        if any(p.dimension != self.dimension for p in objs):
            raise ValueError("Patches with different dimension added")
        if any(p.pardim != self.pardim for p in objs):
            raise ValueError("Patches with different parametric dimension added")

    def _generate(self, objs):
        for p in objs:
            self.catalogue.add(p)

    def generate_cp_numbers(self):
        index = 0
        for node in self.catalogue.top_nodes():
            index = node.generate_cp_numbers(index)
        self.ncps = index
        for node in self.catalogue.top_nodes():
            node.read_cp_numbers()

    def generate_cell_numbers(self):
        index = 0
        for node in self.catalogue.top_nodes():
            index = node.generate_cell_numbers(index)
        self.ncells = index

    def cps(self):
        cps = np.zeros((self.ncps, self.dimension))
        for node in self.catalogue.top_nodes():
            indices = node.cp_numbers.reshape(-1)
            values = node.obj.controlpoints.reshape(-1, self.dimension)
            cps[indices] = values
        return cps

    def faces(self):
        assert self.pardim == 3
        faces = list(chain.from_iterable(node.faces() for node in self.catalogue.top_nodes()))
        return np.hstack(faces)

    def summary(self):
        c = self.catalogue
        while isinstance(c, ObjectCatalogue):
            print('Dim {}: {}'.format(c.pardim, len(c.top_nodes())))
            c = c.lower

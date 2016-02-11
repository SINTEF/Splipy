# -*- coding: utf-8 -*-

from GeoMod import *
from GeoMod.Utils import *
import numpy as np
from collections import Counter
from itertools import product, permutations

try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping


class VertexDict(MutableMapping):
    """A dictionary where the keys are numpy arrays, and where equality
    computed in an a pproximate sense for floating point numbers.

    All keys must have the same dimensions.
    """

    def __init__(self, rtol=0, atol=1e-8):
        """Initialize a vertex dictionary.

        :param float rtol: Relative tolerance to use for equality check
        :param float atol: Absolute tolerance to use for equality check
        """
        self.rtol = rtol
        self.atol = atol

        # List of (key, value) pairs
        self.internal = []

    def _eq(self, a, b):
        """Check whether two numpy arrays are almost equal, according to the given
        tolerances.
        """
        return np.allclose(a, b, rtol=self.rtol, atol=self.atol)

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
    :class:`GeoMod.SplineModel.Orientation` indicating an inability to match
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

        This constructor is for internal use. Use
        :func:`GeoMod.SplineModel.Orientation.compute` to generate new
        orientation objects.
        """
        self.perm = perm
        self.flip = flip
        self.perm_inv = tuple(perm.index(d) for d in range(len(perm)))

    @classmethod
    def compute(cls, cpa, cpb=None):
        """Compute and return a new orientation object representing the mapping between
        `cpa` (the reference system) and `cpb` (the mapped system).

        Each argument can be either a `SplineObject` or a numpy array
        representing the control points.

        If `cpb` is not given, the identity orientation is returned.

        :param cpa: The reference system
        :param cpb: The mapped system
        :return: The orientation relating the two systems
        :rtype: Orientation
        :raises OrientationError: If the two objects do not match
        """
        if not isinstance(cpa, np.ndarray):
            cpa = cpa.controlpoints

        shape_a = cpa.shape
        pardim = len(shape_a) - 1

        # Return the identity orientation if no cpb
        if cpb is None:
            return cls(tuple(range(pardim)),
                       tuple(False for _ in range(pardim)))

        if not isinstance(cpb, np.ndarray):
            cpb = cpb.controlpoints

        shape_b = cpb.shape

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
            transposed = cpb.transpose(perm + (pardim,))
            if transposed.shape != shape_a:
                continue
            # Enumerate all possible direction reversals
            for flip in product([False, True], repeat=pardim):
                slices = tuple(slice(None, None, -1) if f else slice(None) for f in flip)
                test_b = transposed[slices + (slice(None),)]
                # FIXME: Hardcoded tolerances
                if np.allclose(cpa, test_b, rtol=0, atol=1e-8):
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
        :func:`GeoMod.Utils.section`. It should be described relative to the
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
        :func:`GeoMod.Utils.section`. It should be described relative to the
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

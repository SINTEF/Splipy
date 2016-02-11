# -*- coding: utf-8 -*-

from GeoMod import *
import numpy as np

try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping


class VertexDict(MutableMapping):

    def __init__(self, rtol=0, atol=1e-8):
        self.rtol = rtol
        self.atol = atol
        self.internal = []

    def _eq(self, a, b):
        return np.allclose(a, b, rtol=self.rtol, atol=self.atol)

    def _candidate(self, key):
        for i, (k, v) in enumerate(self.internal):
            if self._eq(key, k):
                return i
        raise KeyError(key)

    def __setitem__(self, key, value):
        try:
            c = self._candidate(key)
            k, _ = self.internal[c]
            self.internal[c] = (k, value)
        except KeyError:
            self.internal.append((key, value))

    def __getitem__(self, key):
        c = self._candidate(key)
        return self.internal[c][1]

    def __delitem__(self, key):
        self.internal = [(k, v) for k, v in self.internal
                         if not self._eq(key, k)]

    def __iter__(self):
        for k, _ in self.internal:
            yield k

    def items(self):
        return list(self.internal)

    def __len__(self):
        return len(self.internal)

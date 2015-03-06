__doc__ = 'Class for handling specification of topology sets and automatic renumbering of patches'

import os

from itertools import groupby
from numpy import cumsum, ones, prod, zeros
from operator import itemgetter, methodcaller

from . import HDF5File, InputFile
from .. import WallDistance
from GoTools import *
from GoTools.Preprocess import NaturalNodeNumbers

class Numberer(object):
  """ Handles specification of topology sets and automatic renumbering of patches.
  """

  # A group is a named collection of patches of a single kind
  # ('volume', 'face', etc.)
  class Group:

    def __init__(self, name, kind):
      self.kind = kind
      self.name = name
      self.patches = []


  # A boundary is a named collection of boundary components.
  class Boundary:

    def __init__(self, name):
      self.name = name
      self.components = []

    # Returns a dict of type kind -> list of components of that kind
    def KindsAndComponents(self):
      kinds = {}
      for comp in self.components:
        if not comp.kind in kinds:
          kinds[comp.kind] = []
        kinds[comp.kind].append(comp)
      return kinds


  # A boundary component has a group (applies to all patches in that group),
  # a kind ('face', 'edge', etc.), and a collection of indices.
  #
  # It corresponds to number of <item> entries in a topology set, e.g.
  # for patch in group:
  #   <item patch="patchnumber">indexes</item>
  class BoundaryComponent:

    def __init__(self, group, indexes, kind):
      self.group = group
      self.indexes = indexes
      self.kind = kind


  # An output group is a topology set that is not a boundary, i.e. it only
  # consists of patches.  Its components are simply groups.
  class OutputGroup:

    def __init__(self, name, kind):
      self.kind = kind
      self.name = name
      self.components = []


  # A processor is a collection of patches.
  class Proc:

    def __init__(self):
      self.items = []

    def ndofs(self):
      return sum(map(itemgetter('ndofs'), self.items))


  def __init__(self):
    self._patches = []
    self._groups = {}
    self._boundaries = {}
    self._outputgroups = {}
    self._numbering = None


  def AddPatches(self, patches, allow_left=False):
    """ Adds patches to the geometry.
        @param patches: The patches to add.
        @type patches: List of patches
        @param allow_left: Whether to allow left-handed coordinate systems. (Default False.)
        @type allow_left: Bool
    """
    if not allow_left:
      for i, p in enumerate(patches):
        try:
          lhs = p.MakeRHS()
        except:
          lhs = False
        if lhs:
          raise Exception("Left-handed patch added (index %i)." % i)

    self._patches += patches


  def AddGroup(self, name, kind, objects):
    """ Adds patches to a group. The patches must also be added to the model
        using `AddPatches`.
        @param name: The name of the group (will be created if it doesn't exist).
        Names beginning with two underscores are reserved!
        @type name: String
        @param kind: The kind of the group.  Optional.
        @type kind: 'volume', 'face' or 'edge'
        @param objects: The patches to add.
        @type objects: List of patches
    """
    if kind is None:
      kind = {Curve: 'edge', Surface: 'face', Volume: 'volume'}[type(objects[0])]
    if not name in self._groups:
      self._groups[name] = Numberer.Group(name, kind)
    for o in objects:
      self._groups[name].patches.append(o)


  def Groups(self):
    """ @return: The already defined group names.
        @rtype: List of String
    """
    return self._groups.keys()


  def WriteGroup(self, name, filename):
    """ Writes the patches of a group to a file.
        @param name: The group name.
        @type name: String
        @param filename: The filename to write to.
        @type filename: String
    """
    WriteG2(filename, self._groups[name].patches)


  def WriteGroups(self, filename):
    """ Writes all the groups to files.
        @param filename: The filename base to write to.
        @param filename: String
    """
    for grp in self._groups.keys():
      self.WriteGroup(grp, '%s-grp-%s.g2' % (filename, grp))


  def WriteOutputGroup(self, name, filename):
    """ Writes the patches of an output group to a file.
        @param name: The group name.
        @type name: String
        @param filename: The filename to write to.
        @type filename: String
    """
    patches = []
    for component in self._outputgroups[name].components:
      patches.extend(component.patches)
    WriteG2(filename, patches)


  def WriteOutputGroups(self, filename):
    """ Writes all the output groups to files.
        @param filename: The filename base to write to.
        @param filename: String
    """
    for grp in self._outputgroups.keys():
      self.WriteOutputGroup(grp, '%s-outgrp-%s.g2' % (filename, grp))


  def AddBoundary(self, name, components):
    """ Adds boundary components to a boundary. Each component must be a tuple on the
        form (groupname, kind, indexes), which will add, for each patch in the given group,
        the sub-components of the given kind with the given indexes.  Indexes may be a list
        or a single element.

        E.g. ('mygroup', 'edge', [1,2,3]) will add edges 1-3 from each patch in mygroup.

        Alternatively, the first element of the comoponent may be a single patch or a list of
        patches, in which case a group will automatically be made with the appropriate kind.

        The indices are zero-indexed and must conform to the IFEM convention. (NOT the
        GoTools convention.)

        @param name: The name of the topology set (will be created if it doesn't exist).
        @type name: String
        @param components: Components to add.
        @type components: List of Tuple of (String | Patch | List of Patch, String, Int | List of
        Int), or a single tuple.
    """
    if not name in self._boundaries:
      self._boundaries[name] = Numberer.Boundary(name)
    bnd = self._boundaries[name]

    if type(components) is tuple:
      components = [components]

    for cname, ckind, cidxs in components:
      if type(cidxs) is not list:
        cidxs = [cidxs]
      if type(cname) is str:
        group = self._groups[cname]
      else:
        if type(cname) in [Curve, Surface, Volume]:
          cname = [cname]
        groupname = '__%i' % len(self._groups)
        self.AddGroup(groupname, None, cname)
        group = self._groups[groupname]
      bnd.components.append(Numberer.BoundaryComponent(group, cidxs, ckind))


  def AddOutputGroup(self, name, kind, components):
    """ Adds groups to an output group. An output group is like a boundary, in that
        it will produce a topology set, but it produces a topology set of whole patches,
        not of subcomponents.

        The component may be either the name of a group or a list of such, or a patch or
        list of patches.  In the latter case, the group will automatically be made.

        @param name: The name of the topology set (will be created if it doesn't exist).
        @type name: String
        @param kind: The kind of the topology set.
        @type kind: 'volume', 'face' or 'edge'
        @param components: The groups to add.
        @type components: String | Patch | List of String | List of Patch
    """
    if type(components) in [Curve, Volume, Surface, str]:
      components = [components]
    if type(components[0]) in [Curve, Volume, Surface]:
      groupname = '__%i' % len(self._groups)
      self.AddGroup(groupname, None, components)
      components = [groupname]

    if not name in self._outputgroups:
      self._outputgroups[name] = Numberer.OutputGroup(name, kind)
    outgroup = self._outputgroups[name]
    for cname in components:
      group = self._groups[cname]
      if group.kind != outgroup.kind:
        raise Exception('Numberer: Kind mismatch: %s and %s' % (group.kind, outgroup.kind))
      outgroup.components.append(group)


  def AddWallGroup(self, wallgroup, suffix='_walldist' ):
    """Looks up the namesake boundary group and prepares to write wall distance upon invoking WriteEveryting
       @param wallgroup: name of an existing boundary group
       @type wallgroup: str
       @param suffix: suffix to file basename for wall dist file
       @type suffix: str"""
    assert self._boundaries.has_key( wallgroup ), 'Wall group has to be added to boundary groups first'
    self.wallgroup = wallgroup
    self.wallsuffix = suffix


  def Boundaries(self):
    """ @return: The already defined boundary names.
        @rtype: List of String
    """
    return self._boundaries.keys()


  def WriteBoundary(self, name, filename):
    """ Writes the geometry of a boundary to a file.
        @param name: The boundary name.
        @type name: String
        @param filename: The filename to write to.
        @type filename: String
    """
    items = []
    bnd = self._boundaries[name]
    for bcmp in bnd.components:
      fn = self.GetterFunction(bcmp.group.kind, bcmp.kind)
      for patch in bcmp.group.patches:
        for cidx in bcmp.indexes:
          items.append(fn(patch, cidx))
    WriteG2(filename, items)


  def WriteBoundaries(self, filename):
    """ Writes all the boundaries to files.
        @param filename: The filename base to write to.
        @type filename: String
    """
    for bnd in self._boundaries.keys():
      self.WriteBoundary(bnd, '%s-bnd-%s.g2' % (filename, bnd))


  def Renumber(self, nprocs, outprocs=None):
    """ Generates a numbering of all the patches and assigns a range of patches
        to each processor. Doing this optimally is NP-complete, so this method
        uses a quasi-optimal greedy algorithm.

        All patches must be added before calling Renumber(), but it's not necessary
        that all the topology sets are specified. Renumber() can be called several
        times if necessary.

        @param nprocs: Number of processors to optimize for
        @type nprocs: Int
        @param outprocs: Number of processors to write output for (optional)
        @type outprocs: Int
    """

    if outprocs is None:
      outprocs = nprocs

    if nprocs == 1:
      # Do no special renumbering for nprocs = 1
      self._numbering = [{'patch': p,
                          'ndofs': prod(map(len, p.GetKnots())),
                          'index': i}
                         for i, p in enumerate(self._patches)]
    else:
      items = [{'patch': p,
                'ndofs': prod(map(len, p.GetKnots()))}
               for p in self._patches]
      tot_ndofs = sum(i['ndofs'] for i in items)

      # List of temporary processor objects.  We will attempt to optimize for
      # nprocs processors, but will output for more.
      temp_procs = [Numberer.Proc() for _ in xrange(nprocs)]

      # Add the patches in order of decreasing ndofs to the processor with
      # lowest number of ndofs so far.  This produces a numbering that depends
      # only on the patch list and nprocs. This is a greedy heuristic algorithm
      # for an NP-complete problem.  It is guaranteed that no processor will
      # get more than 4/3 of the optimal load.
      items.sort(key=itemgetter('ndofs'), reverse=True)
      for n in items:
        temp_procs[0].items.append(n)
        temp_procs.sort(key=methodcaller('ndofs'))

      # Establish the numbering.
      self._numbering = []
      for procnum, proc in enumerate(temp_procs):
        for i in proc.items:
          i['index'] = len(self._numbering)
          self._numbering.append(i)

    # Distribute the patches to the correct number of processors, without
    # changing the order.  Due to the order restriction, this problem can
    # be solved optimally quite quickly.
    self._procs = [Numberer.Proc() for _ in xrange(outprocs)]
    ranges = _linpar([i['ndofs'] for i in self._numbering], outprocs)
    for procnum, (rng, proc) in enumerate(zip(ranges, self._procs)):
      for r in rng:
        proc.items.append(self._numbering[r])
        self._numbering[r]['procnum'] = procnum

    self._patchmap = {i['patch']: i['index'] for i in self._numbering}


  # Throws an exception if Renumber() has not been called.
  def CheckNumbering(self):
    if self._numbering is None:
      raise Exception('Renumbering has not been performed yet')


  def GetNumber(self, patch):
    """ Gets the new number of a patch.
        @param patch: The patch.
        @type patch: Patch
        @return: The new patch number (zero-indexed)
        @rtype: Int
    """
    self.CheckNumbering()
    return self._patchmap[patch]


  def Partitioning(self, base_indent=0, indent=2, include_header=True):
    """ Returns the XML representation of the patch partitioning.
        @param base_indent: The base number of spaces to indent with.
        @type base_indent: Int
        @param indent: The number of spaces to add for each indentation leve.
        @type indent: 2
        @param include_header: Whether to include the root element <partitioning>
        @type include_header: Bool
        @return: The XML.
        @rtype: String with trailing newline
    """
    self.CheckNumbering()

    base = lambda k: ' '*(base_indent + k*indent)

    items = map(len, (list(v) for _, v in groupby(self._numbering, itemgetter('procnum'))))
    starts = [0] + list(cumsum(items)[:-1])
    ends = list(cumsum(items) - 1)

    lst = []

    if include_header:
      lst.append(base(0) + '<partitioning procs="%i">' % len(starts))

    for i, (s, e) in enumerate(zip(starts, ends)):
      lst.append(base(1) + '<part proc="%i" lower="%i" upper="%i" />' % (i, s+1, e+1))

    if include_header:
      lst.append(base(0) + '</partitioning>')

    ret = ''
    for t in lst:
      ret += t + '\n'

    return ret


  def TopologySets(self, sets, base_indent=0, indent=2, include_header=True):
    """ Returns the XML representation the topology sets.
        @param sets: The names of the topology sets to write.
        @type sets: List of String
        @param base_indent: The base number of spaces to indent with.
        @type base_indent: Int
        @param indent: The number of spaces to add for each indentation leve.
        @type indent: 2
        @param include_header: Whether to include the root element <partitioning>
        @type include_header: Bool
        @return: The XML.
        @rtype: String with trailing newline
    """
    self.CheckNumbering()

    base = lambda k: ' '*(base_indent + k*indent)

    lst = []
    idxtostr = lambda i: str(i+1)

    if include_header:
      lst.append(base(0) + '<topologysets>')

    for s in sets:
      if s in self._boundaries:
        bnd = self._boundaries[s]
        for kind, comps in bnd.KindsAndComponents().iteritems():
          lst.append(base(1) + '<set name="%s" type="%s">' % (bnd.name, kind))
          for comp in comps:
            for patch in comp.group.patches:
              num = self.GetNumber(patch)
              params = (num + 1, ' '.join(map(idxtostr, comp.indexes)))
              lst.append(base(2) + '<item patch="%i">%s</item>' % params)
          lst.append(base(1) + '</set>')

      if s in self._outputgroups:
        outgroup = self._outputgroups[s]
        lst.append(base(1) + '<set name="%s" type="%s">' % (outgroup.name, outgroup.kind))
        for group in outgroup.components:
          for patch in group.patches:
            num = self.GetNumber(patch)
            lst.append(base(2) + '<item patch="%i"></item>' % (num + 1))
        lst.append(base(1) + '</set>')

    if include_header:
      lst.append(base(0) + '</topologysets>')

    ret = ''
    for t in lst:
      ret += t + '\n'

    return ret


  def WriteEverything(self, filename, indent=2, periodic={}, display=True):
    """ All-in-one method for writing out everything you might need (g2 geometry file,
        natural node numbers and IFEM xml file.
        @param filename: The base filename to write to (no extension).
        @type filename: String
        @param indent: Number of spaces to add for each indentation level in the XML file.
        @type indent: Int
        @param display: (Optional) Write walldistance progress to stdout
        @type display: Boolean
    """
    self.CheckNumbering()

    basename = os.path.basename(filename)
    patchlist = map(itemgetter('patch'), self._numbering)

    WriteG2('%s.g2' % filename, patchlist)

    numbers = NaturalNodeNumbers(patchlist,**periodic)
    f = HDF5File('%s_nodenumbers' % filename)
    for i, (p, n) in enumerate(zip(patchlist, numbers)):
      f.AddGeometry('Common', i+1, 0, p)
      f.AddField('Common', 'node numbers', i+1, 0, 1, n)

    with open('%s.xinp' % filename, 'w') as f:
      f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n\n')
      f.write('<geometry>\n')

      f.write(self.Partitioning(indent, indent, True))
      f.write(' '*indent + '<patchfile>%s.g2</patchfile>\n' % basename)
      f.write(' '*indent + '<nodefile>%s_nodenumbers.hdf5</nodefile>\n' % basename)
      f.write(self.TopologySets(self._boundaries.keys() + self._outputgroups.keys(),
                                indent, indent, True))

      f.write('</geometry>\n')

    if hasattr( self, 'wallsuffix' ):
      wallinfo = InputFile( filename+'.xinp' ).GetTopologySet( self.wallgroup )
      dist = WallDistance.distanceFunction( patchlist, wallinfo, display=display )
      g = HDF5File( filename+self.wallsuffix )
      for i, (p, d) in enumerate(zip(patchlist, dist)):
        g.AddGeometry('Common', i+1, 0, p)
        g.AddField('Common', 'wall distance', i+1, 0, 1, d)


  def PrintLoadBalance(self):
    """ Writes a summary of the load balancing results to the console.
    """
    self.CheckNumbering()

    tot = 0

    print 'Load balancing results:'
    for ndofs, procs in groupby(self._procs, methodcaller('ndofs')):
      npatches = len(list(procs))
      print '~%i DoFs: %i processors' % (ndofs, npatches)
      tot += ndofs * npatches

    print 'Total: ~%i DoFs' % tot


  # Returns the appropriate getter function for getting type <target>
  # from type <source>, where source and target are 'volume', 'face',
  # 'edge', 'vertex'.
  def GetterFunction(self, source, target):
    _mapFE = {0: 3, 1: 1, 2: 0, 3: 2}

    _mapVE = {0: [ 4, 10,  6,  8],
              1: [ 5, 11,  7,  9],
              2: [ 0,  9,  2,  8],
              3: [ 1, 11,  3, 10],
              4: [ 0,  5,  1,  4],
              5: [ 2,  7,  3,  6]}

    # Gets vertex #i from the edge
    def GetEdgeVertex(edge, i):
      return edge.Evaluate(edge.GetKnots()[-1 if i == 1 else 0])

    # Gets vertex #i from the face (IFEM convention, zero-indexed)
    def GetFaceVertex(face, i):
      edge = face.GetEdges()[0 if i < 2 else 2]
      return GetEdgeVertex(edge, i % 2)

    # Gets edge #i from the face (IFEM convention, zero-indexed)
    def GetFaceEdge(face, i):
      return face.GetEdges()[_mapFE[i]]

    # Gets vertex #i from the volume (IFEM convention, zero-indexed)
    def GetVolVertex(vol, i):
      face = vol.GetFaces()[4 if i < 4 else 5]
      return GetFaceVertex(face, i % 4)

    # Gets edge #i from the volume (IFEM convention, zero-indexed)
    def GetVolEdge(vol, i):
      for f_idx, e_idxs in _mapVE.iteritems():
        if i in e_idxs:
          return vol.GetFaces()[f_idx].GetEdges()[e_idxs.index(i)]

    # Gets face #i from the volume (IFEM convention, zero-indexed)
    def GetVolFace(vol, i):
      return vol.GetFaces()[i]


    # Pick the right getter
    if source == 'volume':
      if target == 'face':
        return GetVolFace
      elif target == 'edge':
        return GetVolEdge
      elif target == 'vertex':
        return GetVolVertex

    elif source == 'face':
      if target == 'edge':
        return GetFaceEdge
      elif target == 'vertex':
        return GetFaceVertex

    elif source == 'edge':
      if target == 'vertex':
        return GetEdgeVertex

    raise Exception("No getter function from %s to %s" % (source, target))



# This utility function solves the linear partitioning problem.
# Given an array of positive numbers [s_1, ..., s_n], find k
# subintervals, so that the maximal sum over them is minimized.
# This is a load balancing problem with restricted ordering.
# It is solved using a dynamic programming algorithm as explained
# here:
# http://www8.cs.umu.se/kurser/TDBAfl/VT06/algorithms/BOOK/BOOK2/NODE45.HTM
def _linpar(s, k):
  prefixes = cumsum(s)
  n = len(s)

  M = zeros((n, k), dtype=int)  # Cost of each partition, we are interested in M[-1,-1]
  D = -ones((n, k), dtype=int)  # Splits, used to reconstruct the solution afterwards
  M[:,0] = prefixes
  M[0,:] = s[0]

  for i in xrange(1, n):
    for j in xrange(1, k):
      M[i,j] = -1
      for x in xrange(0, i):
        cost = max(M[x,j-1], prefixes[i] - prefixes[x])
        if M[i,j] > cost or M[i,j] < 0:
          M[i,j] = cost
          D[i,j] = x

  ranges = []

  def reconstruct(n, k):
    if k == 0:
      ranges.append(xrange(0, n+1))
    else:
      split = D[n,k]
      assert(split >= 0)
      reconstruct(split, k-1)
      ranges.append(xrange(split+1, n+1))

  reconstruct(n-1, k-1)

  return ranges

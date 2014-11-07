__doc__ = 'Implementation of file I/O classes.'

import xml.dom.minidom
import os
import sys
import json

from collections import namedtuple
from itertools import groupby
from numpy import cumsum, prod
from operator import itemgetter, methodcaller
from . import WallDistance

from GoTools import *
# import scitools.filetable
from GoTools.Preprocess import NaturalNodeNumbers

def topologystring(geotype, name, patches, entries):
  xml = '  <set name="%s" type="%s">\n' % (name, geotype)
  for patch in patches:
    xml += '    <item patch="%i">%s</item>\n' % (patch, entries)
  xml += '  </set>\n'
  return xml

def topologystringfrompairs(geotype, name, entries):
  xml = '  <set name="%s" type="%s">\n' % (name, geotype)
  for (patch, face) in entries:
    xml += '    <item patch="%i">%s</item>\n' % (patch, face)
  xml += '  </set>\n'
  return xml

def writeAsc(X, U, fname):
  f = open(fname,'w')
  for i in range(0,len(U)):
    f.write('%f %f\n' % (X[i], U[i]))
  f.close()

def ParseArgs(args, defaults):
  """ Parses command line arguments in a standardized way. Takes a list of command line
      arguments and a dictionary mapping arguments to default values. The defaults dict
      will be modified. The arguments must be of the form: 'argname=newvalue'

      This works for string, integer and float values (as determined by the default).
      For boolean values, an argument 'argname' will set argname to True, while an argument
      'noargname' will set it to False.

      Arguments are parsed in order, so later values can override earlier ones.

      A special argument of form 'paramfile=somefile.json' will load 'somefile.json' and apply
      the arguments in that file, in order. To specify boolean arguments in JSON, give true
      or false as the value.

      @param args: Command-line arguments
      @type args: List of String
      @param defaults: Defaults
      @type defaults: Dict
      @returns: None
      @rtype: NoneType
  """

  def set_def(key, value):
    for k, v in defaults.iteritems():
      if value is None and key in [k, 'no' + k] and isinstance(v, bool):
        defaults[k] = k == key
      elif k == key:
        if isinstance(v, bool):
          defaults[k] = bool(value)
        elif isinstance(v, float):
          defaults[k] = float(value)
        elif isinstance(v, int):
          defaults[k] = int(value)
        else:
          defaults[k] = str(value)

  for arg in args:
    if arg.startswith('paramfile=') :
      with open(arg[10:], 'r') as f:
        dc = json.load(f)
      for k, v in dc.iteritems():
        set_def(k, v)

    split = arg.split('=', 1)
    set_def(split[0], split[1] if len(split) > 1 else None)

def readMCfile(filename,remlast=True, silent=False):
    """
    Read text file with multiple columns
    If remlast = True, the last row will be deleted
    Return as tuple (x,y)
    """

    # Open file for reading
    try:
        readfile = open(filename,'r')
    except IOError:
        if silent != True:
            print 'Utils.FileInputOutput.py::readMCfile: Could not open file ' + filename + ' for reading.'
        #sys.exit(1)
        return (None, None)

    try:
        data = scitools.filetable.read(readfile)
    except Exception, e:
        print 'Utils.FileInputOutput.py::readMCfile: Could not read from file ' + filename + '.'
        print str(e)
        data = None

        #sys.exit(1)

    readfile.close()
    if data == None:
        print 'Utils.FileInputOutput.py::readMCfile: Trying alterative file reading method for' + filename + '.'
        data =  readLinesWithSeparation(filename,remlast=True, silent=silent)
        #pdb.set_trace()
    try:
      if remlast == True:
          x = data[0:-1,0]        # Remove last line which often is problematic
          y = data[0:-1,1:]       # Remove last line which often is problematic
      else:
          x = data[:,0]
          y = data[:,1:]
      x.reshape(x.shape[0],1)
    except TypeError:
      return (None, None)

    return (x,y)

class InputFile:
  """Class for working with IFEM input (.xinp) files.
     @param path: The path of the xinp to open
     @type path: String
  """
  PatchInfo = namedtuple("PatchInfo","vertex edge face")
  def __init__(self, path):
    self.path = path
    self.abspath = os.path.abspath(path)
    self.dom = xml.dom.minidom.parse(path)

  def GetGeometryFile(self):
    """ Extract the geometry definition (.g2 file)
        @return: The file name
        @rtype: String
    """
    geometry = self.dom.getElementsByTagName('geometry')[0]
    result = geometry.getElementsByTagName('patchfile')[0]
    return result.childNodes[0].nodeValue

  def GetTopologySets(self, context="", nocontext=False):
    """ Return a list of all topology sets.
        @param context: App context to read from (optional)
        @type context: String
        @return: Topology sets.
        @rtype: List of String
    """

    if nocontext:
      topologyset = self.dom.getElementsByTagName('topologysets')[0]
    else:
      if len(context):
        ctx = self.dom.getElementsByTagName(context)[0]
        geometry = ctx.getElementsByTagName('geometry')[0]
      else:
        geometry = self.dom.getElementsByTagName('geometry')[0]
      topologyset = geometry.getElementsByTagName('topologysets')[0]

    names = [n.getAttributeNode('name').nodeValue
             for n in topologyset.getElementsByTagName('set')]
    return list(set(names))


  def GetTopologySet(self, name, convention="gotools", context="", nocontext=False, toptype=""):
    """ Extract a topology set from the input file.
        @param name: Name of topology set
        @type name: String
        @param convention: Which numbering convention to use (optional)
        @type convention: 'gotools' (default) or 'ifem'
        @param context: App context to read from (optional)
        @type context: String
        @param toptype: Get only items of a given type (optional)
        @type toptype: String
        @return: Requested topologyset
        @rtype: Dict with PatchInfo
        """
    result = {}

    if type(name) is list:
      names = name
    else:
      names = [name]

    if nocontext:
      topologyset = self.dom.getElementsByTagName('topologysets')[0]
    else:
      if len(context):
        ctx = self.dom.getElementsByTagName(context)[0]
        geometry = ctx.getElementsByTagName('geometry')[0]
      else:
        geometry = self.dom.getElementsByTagName('geometry')[0]
      topologyset = geometry.getElementsByTagName('topologysets')[0]

    for name in names:
      for topset in topologyset.getElementsByTagName('set'):
        if topset.getAttributeNode('name').nodeValue == name:
          typ = topset.getAttributeNode('type').nodeValue
          for item in topset.getElementsByTagName('item'):
            patch = int(item.getAttributeNode('patch').nodeValue)
            values = map(int, item.childNodes[0].nodeValue.split())
            if not result.has_key(patch):
              result[patch] = InputFile.PatchInfo([], [], [])

            if convention == 'gotools':
              if typ == 'edge':
                if toptype == 'edge' or len(toptype) == 0:
                  remap = [4,2,1,3]
                  result[patch].edge.extend([remap[v-1] for v in values])
              elif typ == 'face':
                if toptype == 'face' or len(toptype) == 0:
                  remap = [1,2,5,6,3,4]
                  result[patch].face.extend([remap[v-1] for v in values])
              else:
                result[patch].vertex.extend(values)

            elif convention == 'ifem':
              if typ == 'edge':
                target = result[patch].edge
              elif typ == 'face':
                target = result[patch].face
              else:
                target = result[patch].vertex
              target.extend(values)

    return result

  def write(self, name):
    """ Write the input file to disk.
        @param name: The filename
        @type name: String
        @return: None
    """
    f = open(name,'w')
    f.write(self.dom.toxml('utf-8'))
    f.close()

  def writeOpenFOAM(self, name):
    """ Writes an equivalent mesh in OpenFOAM format. The topology sets
        will be written as disjointed parts, e.g. the set 'myset' will
        yield OpenFOAM patches named 'myset-001', 'myset-002', etc. You
        can use regular expressions in OpenFOAM to group them together,
        or use the createPatch utility to create a new patch from a group
        of existing patches.

        @param name: The path to save to. This must be a directory.
                     Will be created if it does not exist.
        @type name: String
        @return: None
    """

    # Ensure that the output location is sane
    try:
      os.makedirs(name)
    except OSError:
      if not os.path.isdir(name):
        raise Exception("The path '%s' already exists and is not a directory." % name)

    # Header function
    def header(cls, obj, note=None):
      s = 'FoamFile\n{\n'
      s += '    version     2.0;\n'
      s += '    format      ascii;\n'
      s += '    class       %s;\n' % cls
      if note:
        s += '    note        "%s";\n' % note
      s += '    object      %s;\n' % obj
      s += '}\n'
      return s

    patchlist = ReadG2(os.path.join(os.path.dirname(self.abspath), self.GetGeometryFile()))
    for p in patchlist:
      if p.MakeRHS():
        raise Exception("Found left-handed patch. Aborting.")


    # Write points
    nodenumbers = NaturalNodeNumbers(patchlist)
    maxnum = max(map(max, nodenumbers))
    points = [None] * (maxnum + 1)
    for patchid, patch in enumerate(patchlist):
      kus, kvs, kws = patch.GetKnots()
      patchpoints = [patch.Evaluate(ku, kv, kw)
                     for kw in kws for kv in kvs for ku in kus]
      for localid, pt in enumerate(patchpoints):
        points[nodenumbers[patchid][localid]] = pt

    with open(os.path.join(name, 'points'), 'w') as out:
      out.write(header('vectorField', 'points'))
      out.write('%i\n' % len(points))
      out.write('(\n')
      for pt in points:
        out.write('(%f %f %f)\n' % tuple(pt))
      out.write(')\n')

    # Face-cell connectivity structure
    sort_tuple = lambda tup: tuple(sorted(tup))
    rev_tuple = lambda tup: tuple(reversed(tup))
    localid = lambda iu, iv, iw, nu, nv: iw*nu*nv + iv*nu + iu

    facedict = {}
    cellnum = 0
    face_idxs = [(0,2,6,4), (1,5,7,3), (0,1,3,2), (4,6,7,5), (0,4,5,1), (2,3,7,6)]
    for patchid, patch in enumerate(patchlist):
      nu, nv, nw = map(len, patch.GetKnots())
      localid_patch = lambda (iu, iv, iw): localid(iu, iv, iw, nu, nv)

      for iw in range(nw-1):
        for iv in range(nv-1):
          for iu in range(nu-1):
            globalids = [nodenumbers[patchid][localid_patch(tp)]
                         for tp in [(iu,iv,iw), (iu+1,iv,iw), (iu,iv,iw+1), (iu+1,iv,iw+1),
                                    (iu,iv+1,iw), (iu+1,iv+1,iw), (iu,iv+1,iw+1), (iu+1,iv+1,iw+1)]]
            tuples = [(iu,iv,iw), (iu+1,iv,iw), (iu,iv,iw+1), (iu+1,iv,iw+1),
                      (iu,iv+1,iw), (iu+1,iv+1,iw), (iu,iv+1,iw+1), (iu+1,iv+1,iw+1)]
            for fidxs in face_idxs:
              face = [globalids[fid] for fid in fidxs]
              if sort_tuple(face) in facedict:
                assert(facedict[sort_tuple(face)]['neighbour'] == -1)
                facedict[sort_tuple(face)]['neighbour'] = cellnum
              else:
                facedict[sort_tuple(face)] = {'tuple': tuple(face),
                                              'owner': cellnum,
                                              'neighbour': -1}
            cellnum += 1

    internal_faces = [f for f in facedict.values() if f['neighbour'] >= 0]
    boundary_faces = [f for f in facedict.values() if f['neighbour'] == -1]

    # Identify boundaries
    def get_faces(patchid, faceid):
      nu, nv, nw = map(len, patchlist[patchid].GetKnots())
      localid_patch = lambda iu, iv, iw: localid(iu, iv, iw, nu, nv)
      ret = []
      if faceid <= 2:
        iu = 0 if faceid == 1 else nu - 1
        return [sort_tuple((nodenumbers[patchid][localid_patch(iu,iv,iw)],
                            nodenumbers[patchid][localid_patch(iu,iv+1,iw)],
                            nodenumbers[patchid][localid_patch(iu,iv+1,iw+1)],
                            nodenumbers[patchid][localid_patch(iu,iv,iw+1)]))
                for iv in range(0, nv-1) for iw in range(0, nw-1)]
      if faceid <= 4:
        iv = 0 if faceid == 3 else nv - 1
        return [sort_tuple((nodenumbers[patchid][localid_patch(iu,iv,iw)],
                            nodenumbers[patchid][localid_patch(iu+1,iv,iw)],
                            nodenumbers[patchid][localid_patch(iu+1,iv,iw+1)],
                            nodenumbers[patchid][localid_patch(iu,iv,iw+1)]))
                for iu in range(0, nu-1) for iw in range(0, nw-1)]
      iw = 0 if faceid == 5 else nw - 1
      return [sort_tuple((nodenumbers[patchid][localid_patch(iu,iv,iw)],
                          nodenumbers[patchid][localid_patch(iu+1,iv,iw)],
                          nodenumbers[patchid][localid_patch(iu+1,iv+1,iw)],
                          nodenumbers[patchid][localid_patch(iu,iv+1,iw)]))
              for iu in range(0, nu-1) for iv in range(0, nv-1)]

    for f in boundary_faces:
      f['boundary'] = None

    for ts_name in self.GetTopologySets():
      ts = self.GetTopologySet(ts_name, convention='ifem', toptype='face')
      for patchid, patchinfo in ts.iteritems():
        for faceid in patchinfo.face:
          faces = get_faces(patchid-1, faceid)
          for f in faces:
            facedict[f]['boundary'] = ts_name

    # Order faces
    internal_faces.sort(key=itemgetter('neighbour'))
    internal_faces.sort(key=itemgetter('owner'))
    boundary_faces.sort(key=itemgetter('owner'))
    boundary_faces.sort(key=itemgetter('boundary'))

    nfaces = len(facedict)

    # Write faces
    facenum = 0
    with open(os.path.join(name, 'faces'), 'w') as out:
      out.write(header('faceList', 'faces'))
      out.write('%i\n' % nfaces)
      out.write('(\n')
      for f in internal_faces + boundary_faces:
        out.write('4(%i %i %i %i)\n' % f['tuple'])
        f['label'] = facenum
        facenum += 1
      out.write(')\n')

    # Write owners and neighbours
    note = 'nPoints: %i nCells: %i nFaces: %i nInternalFaces: %i' % (len(points),
                                                                     cellnum,
                                                                     len(facedict),
                                                                     len(internal_faces))

    for thing in ['owner', 'neighbour']:
      with open(os.path.join(name, thing), 'w') as out:
        out.write(header('labelList', thing, note=note))
        out.write('%i\n' % nfaces)
        out.write('(\n')
        for f in internal_faces + boundary_faces:
          out.write('%i\n' % f[thing])
        out.write(')\n')

    # Write boundaries
    bdpatches = [(n,list(f)) for n,f in groupby(boundary_faces, itemgetter('boundary'))]
    with open(os.path.join(name, 'boundary'), 'w') as out:
      out.write(header('polyBoundaryMesh', 'boundary'))
      out.write('%i\n' % len(bdpatches))
      out.write('(\n')
      for name, faces in bdpatches:
        if name is None:
          print >> sys.stderr, "Warning: Found boundary faces with no topology set."
          continue
        out.write('%s\n' % name)
        out.write('{\n')
        out.write('    type patch;\n')
        out.write('    nFaces %i;\n' % len(faces))
        out.write('    startFace %i;\n' % faces[0]['label'])
        out.write('}\n')
      out.write(')\n')



class HDF5File:
  """ Handle output of fields and geometries to HDF5+XML
      @param prefix: Prefix for filename (no extension)
      @type prefix: String
  """
  FieldInfo = namedtuple("FieldInfo","basis patches components fieldtype")
  def __init__(self, prefix):
    self.prefix = prefix
    self.create = True
    self.basismap = {}

  def __del__(self):
    """ The destructor writes the final .xml file to disk."""
    if (len(self.basismap)):
      f = open(self.prefix+'.xml','w')
      xml ='<info>\n'
      for entry in self.basismap:
        xml += '  <entry name="%s" description="%s" type="%s" basis="%s" patches="%i" components="%i"/>\n' \
                %(entry, entry, self.basismap[entry].fieldtype, self.basismap[entry].basis, self.basismap[entry].patches, self.basismap[entry].components)
      xml += '  <levels>0</levels>\n</info>\n'
      f.write(xml)
      f.close()

  def AddGeometry(self, name, patch, level, data):
    """ Add a geometry basis to the HDF5 file.
        @param name: Name of basis
        @type name: String
        @param patch: Patch number
        @type patch: Integer
        @param level: Time level to write at
        @type level: Integer
        @param data: The entity to write
        @type data: Curve, Surface or Volume
        @return: None
     """
    WriteHDF5Geometry(self.prefix+'.hdf5', name, patch, level, data, self.create)
    self.create = False

  def AddField(self, basis, name, patch, level, components, data):
    """ Add a field to the HDF5 file.
        @param basis: Name of basis of the field
        @type name: String
        @param name: Name of field
        @type name: String
        @param patch: Patch number
        @type patch: Integer
        @param level: Time level to write at
        @type level: Integer
        @param data: The entity to write
        @param components: Number of components in field
        @type components: Integer
        @param data: The coefficients of the field
        @type data: List of float
        @return: None
     """
    WriteHDF5Field(self.prefix+'.hdf5', name, patch, level, data, self.create)
    if not self.basismap.has_key(name):
      self.basismap[name] = HDF5File.FieldInfo('', 0, 1,'')
    typ = 'field'
    if type(data[0]) is int:
      typ = 'intfield'
    self.basismap[name] = HDF5File.FieldInfo(basis, max(self.basismap[name].patches, patch), components, typ)

class IFEMResultDatabase:
  """ Class for working with IFEM data output files
      @param path: The path to the xml/hdf5 pair (no extension)
      @type path: String
  """
  def __init__(self, path):
    self.dom = xml.dom.minidom.parse(path+'.xml')
    self.path = path

  def GetBasisForField(self, field):
    """ Extract the basis associated with a given field
          @param field: The field name
          @type  geom: String
          @return: Basis for field
          @rtype: String
    """
    # Find number of patches in field
    fields = self.dom.getElementsByTagName('entry')
    patches = 1
    for f in fields:
      if f.getAttributeNode('name').nodeValue == field:
        return f.getAttributeNode('basis').nodeValue

    return ''

  def _GetPatchesForBasis(self, basis):
    """ Get number of patches in a basis
          @param basis: The basis name
          @type basis: String
          @return: Number of patches
          @rtype: Integer
    """
    fields = self.dom.getElementsByTagName('entry')
    patches = 1
    for field in fields:
      if field.getAttributeNode('basis').nodeValue == basis:
        patches = int(field.getAttributeNode('patches').nodeValue)

    return patches

  def GetGeometry(self, basis, level):
    """ Extract the geometry definition for a given field
        @param geom: The basis name
        @type  geom: String
        @param level: Level to read geometry at
        @type level: Integer
        @return The geometry basis
        @rtype: List of curves, surfaces or volumes)
    """
    patches = self._GetPatchesForBasis(basis)

    res = []
    for i in range(1, patches+1):
      res.append(ReadHDF5Geometry(self.path+'.hdf5', basis, i, level))

    return res

  def GetFieldCoefs(self, field, level):
    """ Extract the coefficients for a given field
        @param field: The field name
        @type  field: String
        @param level: Level to read field at
        @type level: Integer
        @return The field coefficients
        @rtype: List of float
    """
    basis = self.GetBasisForField(field)
    patches = self._GetPatchesForBasis(basis)

    res = []
    for i in range(1, patches+1):
      res.append(ReadHDF5Field(self.path+'.hdf5', field, i, level))

    return res

  def GetField(self, field, level, geom=None):
    """ Extract a given field
        @param field: The field name
        @type  field: String
        @param level: Level to read field at
        @type level: Integer
        @param geom: Basis/geometry
        @type geom: Curve, Volume, Surface or NoneType
        @return The field
        @rtype: List of (curves, surfaces or volumes)
    """
    if geom is None:
      basis = self.GetBasisForField(field)
      geom = self.GetGeometry(basis, 0)

    res = []
    for i in range(len(geom)):
      coefs = ReadHDF5Field(self.path+'.hdf5', field, i+1, level)
      if isinstance(geom[0], Volume):
        k1, k2, k3 = geom[i].GetKnots(True)
        order = geom[i].GetOrder()
        res.append(Volume(order[0], order[1], order[2], k1, k2, k3, coefs))
      elif isinstance(geom[0], Surface):
        k1, k2 = geom[i].GetKnots(True)
        order = geom[i].GetOrder()
        res.append(Surface(order[0], order[1], k1, k2, coefs))
      elif isinstance(geom[0], Curve):
        k1 = geom[i].GetKnots(True)
        order = geom[i].GetOrder()
        res.append(Surface(order, k1, coefs))
      del coefs
    del geom

    return res

  def GetTimeLevels(self):
    """ Extract number of timelevels in result data set
        @return: Number of time levels
        @rtype: Integer
    """
    return int(self.dom.getElementsByTagName('levels')[0].firstChild.nodeValue)

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
      for p in patches:
        try:
          lhs = p.MakeRHS()
        except:
          lhs = False
        if lhs:
          raise Exception("Left-handed patch added.")

    self._patches += patches


  def AddGroup(self, name, kind, objects):
    """ Adds patches to a group. The patches must also be added to the model
        using `AddPatches`.
        @param name: The name of the group (will be created if it doesn't exist).
        @type name: String
        @param kind: The kind of the group.
        @type kind: 'volume', 'face' or 'edge'
        @param objects: The patches to add.
        @type objects: List of patches
    """
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


  def AddBoundary(self, name, components):
    """ Adds boundary components to a boundary. Each component must be a tuple on the
        form (groupname, kind, indexes), which will add, for each patch in the given group,
        the sub-components of the given kind with the given indexes.

        E.g. ('mygroup', 'edge', [1,2,3]) will add edges 1-3 from each patch in mygroup.

        The indices are zero-indexed and must conform to the IFEM convention. (NOT the
        GoTools convention.)

        @param name: The name of the topology set (will be created if it doesn't exist).
        @type name: String
        @param components: Components to add.
        @type components: List of Tuple of (String, String, List of Int)
    """
    if not name in self._boundaries:
      self._boundaries[name] = Numberer.Boundary(name)
    bnd = self._boundaries[name]
    for cname, ckind, cidxs in components:
      group = self._groups[cname]
      bnd.components.append(Numberer.BoundaryComponent(group, cidxs, ckind))


  def AddOutputGroup(self, name, kind, components):
    """ Adds groups to an output group. An output group is like a boundary, in that
        it will produce a topology set, but it produces a topology set of whole patches,
        not of subcomponents.

        @param name: The name of the topology set (will be created if it doesn't exist).
        @type name: String
        @param kind: The kind of the topology set.
        @type kind: 'volume', 'face' or 'edge'
        @param components: The groups to add.
        @type components: List of String
    """
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


  def Renumber(self, nprocs):
    """ Generates a numbering of all the patches and assigns a range of patches
        to each processor. Doing this optimally is NP-complete, so this method
        uses a quasi-optimal greedy algorithm.

        All patches must be added before calling Renumber(), but it's not necessary
        that all the topology sets are specified. Renumber() can be called several
        times if necessary.
        @param nprocs: Number of processors to optimize for
        @type nprocs: Int
    """
    items = [{'patch': p,
              'ndofs': prod(map(len, p.GetKnots()))}
             for p in self._patches]

    items.sort(key=itemgetter('ndofs'), reverse=True)

    self._procs = [Numberer.Proc() for _ in xrange(nprocs)]

    for n in items:
      self._procs[0].items.append(n)
      self._procs.sort(key=methodcaller('ndofs'))

    self._numbering = []
    for procnum, proc in enumerate(self._procs):
      for i in proc.items:
        i['procnum'] = procnum
        i['index'] = len(self._numbering)
        self._numbering.append(i)

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


  def WriteEverything(self, filename, indent=2):
    """ All-in-one method for writing out everything you might need (g2 geometry file,
        natural node numbers and IFEM xml file.
        @param filename: The base filename to write to (no extension).
        @type filename: String
        @param indent: Number of spaces to add for each indentation level in the XML file.
        @type indent: Int
    """
    self.CheckNumbering()

    basename = os.path.basename(filename)
    patchlist = map(itemgetter('patch'), self._numbering)

    WriteG2('%s.g2' % filename, patchlist)

    numbers = NaturalNodeNumbers(patchlist)
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
      dist = WallDistance.distanceFunction( patchlist, wallinfo )
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

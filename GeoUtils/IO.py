__doc__ = 'Implementation of file I/O classes.'

import xml.dom.minidom
from collections import namedtuple
from GoTools import *

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

class InputFile:
  """ Class for working with IFEM input (.xinp) files.
      @param path: The path of the xinp to open
      @type path: String
  """
  PatchInfo = namedtuple("PatchInfo","vertex edge face")
  def __init__(self, path):
    self.dom = xml.dom.minidom.parse(path)

  def GetGeometryFile(self):
    """ Extract the geometry definition (.g2 file)
        @return: The file name
        @rtype: String
    """
    geometry = self.dom.getElementsByTagName('geometry')[0]
    result = geometry.getElementsByTagName('patchfile')[0]
    return result.childNodes[0].nodeValue


  def GetTopologySet(self, name, context="", nocontext=False, toptype=""):
    """ Extract a topology set from the input file.
        @param name: Name of topology set
        @type name: String
        @param context: App context to read from (optional)
        @type context: String
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
            if not result.has_key(patch):
              result[patch] = InputFile.PatchInfo([], [], [])
            if typ == 'edge':
              if  toptype == 'edge' or len(toptype) == 0:
                remap = [4,2,1,3]
                ed = int(item.childNodes[0].nodeValue)
                result[patch].edge.append(remap[ed-1])
            elif typ == 'face':
              if toptype == 'face' or len(toptype) == 0:
                remap = [1,2,5,6,3,4]
                fa = int(item.childNodes[0].nodeValue)
                result[patch].face.append(remap[fa-1])
            else:
              result[int(item.getAttributeNode('patch').nodeValue)].vertex.append(int(item.childNodes[0].nodeValue))
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

  """ Extract the basis associated with a given field
      @param field: The field name
      @type  geom: String
      @return: Basis for field
      @rtype: String
  """
  def GetBasisForField(self, field):
    # Find number of patches in field
    fields = self.dom.getElementsByTagName('entry')
    patches = 1
    for f in fields:
      if f.getAttributeNode('name').nodeValue == field:
        return f.getAttributeNode('basis').nodeValue

    return ''

  """ Get number of patches in a basis
      @param basis: The basis name
      @type basis: String
      @return: Number of patches
      @rtype: Integer
  """
  def _GetPatchesForBasis(self, basis):
    fields = self.dom.getElementsByTagName('entry')
    patches = 1
    for field in fields:
      if field.getAttributeNode('basis').nodeValue == basis:
        patches = int(field.getAttributeNode('patches').nodeValue)

    return patches

  """ Extract the geometry definition for a given field
      @param geom: The basis name
      @type  geom: String
      @param level: Level to read geometry at
      @type level: Integer
      @return The geometry basis
      @rtype: List of curves, surfaces or volumes)
  """
  def GetGeometry(self, basis, level):
    patches = self._GetPatchesForBasis(basis)

    res = []
    for i in range(1, patches+1):
      res.append(ReadHDF5Geometry(self.path+'.hdf5', basis, i, level))

    return res

  """ Extract the coefficients for a given field
      @param field: The field name
      @type  field: String
      @param level: Level to read field at
      @type level: Integer
      @return The field coefficients
      @rtype: List of float
  """
  def GetFieldCoefs(self, field, level):
    basis = self.GetBasisForField(field)
    patches = self._GetPatchesForBasis(basis)

    res = []
    for i in range(1, patches+1):
      res.append(ReadHDF5Field(self.path+'.hdf5', field, i, level))

    return res

  """ Extract a given field
      @param field: The field name
      @type  field: String
      @param level: Level to read field at
      @type level: Integer
      @return The field
      @rtype: List of (curves, surfaces or volumes)
  """
  def GetField(self, field, level):
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
    return res

  """ Extract number of timelevels in result data set
      @return: Number of time levels
      @rtype: Integer
  """
  def GetTimeLevels(self):
    return int(self.dom.getElementsByTagName('levels').nodeValue)

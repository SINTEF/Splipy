__doc__ = 'Class for working with IFEM data output files'

from collections import namedtuple
import xml.dom.minidom
from GoTools import *

class IFEMResultDatabase:
  """ Class for working with IFEM data output files
      @param path: The path to the xml/hdf5 pair (no extension)
      @type path: String
  """
  FieldInfo = namedtuple("FieldInfo","basis patches components fieldtype once description")
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
      try:
        if field.getAttributeNode('basis').nodeValue == basis:
          patches = int(field.getAttributeNode('patches').nodeValue)
      except:
        pass

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
      group = "%i/%i" %(level,i)
      res.append(ReadHDF5Field(self.path+'.hdf5', field, group))

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
      group = "%i/%i" %(level,i+1)
      coefs = ReadHDF5Field(self.path+'.hdf5', field, group)
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

  def GetFieldInfo(self, field):
    """ Extract info for a given field
        @param field: The field name
        @type field: String
        @return: Field infos for given name
        @rtype: List of FieldInfo
    """
    result = []
    nodes = self.dom.getElementsByTagName('entry')
    for node in nodes:
      if node.getAttributeNode('name').nodeValue == field:
        name    = field
        desc    = node.getAttributeNode('description').nodeValue
        ftype   = node.getAttributeNode('type').nodeValue
        try:
          basis   = node.getAttributeNode('basis').nodeValue
        except:
           basis = ''
        try:
          patches = int(node.getAttributeNode('patches').nodeValue)
        except:
          patches = 1
        try:
          comp    = int(node.getAttributeNode('components').nodeValue)
        except:
          comp = 1
        once    = False
        if node.getAttributeNode('once'):
          once = node.getAttributeNode('once').nodeValue == 'true'
        result.append(IFEMResultDatabase.FieldInfo(basis,patches,comp,ftype,once,desc))

    return result

  def GetFields(self):
    result = []
    nodes = self.dom.getElementsByTagName('entry')
    for node in nodes:
      result.append(node.getAttributeNode('name').nodeValue)

    return list(set(result)) # to get unique values

  def GetTimeLevels(self):
    """ Extract number of timelevels in result data set
        @return: Number of time levels
        @rtype: Integer
    """
    return int(self.dom.getElementsByTagName('levels')[0].firstChild.nodeValue)

  def GetTime(self, level):
    """ Read the time value for a given level
        @return: The time value
        @rtype: Float
    """
    return ReadHDF5Field(self.path+'.hdf5', 'timeinfo/SIMbase-1', '%i'%(level))[0]

  def GetIntField(self, level, name):
    """ Read an integer field for a given level
        @return: The field
        @rtype: List of int
    """
    return ReadHDF5Field(self.path+'.hdf5', name, '%i' %(level), True)

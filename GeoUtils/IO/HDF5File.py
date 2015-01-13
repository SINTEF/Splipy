__doc__ = 'Handle output of fields and geometries to HDF5+XML'

from collections import namedtuple
from GoTools import *

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

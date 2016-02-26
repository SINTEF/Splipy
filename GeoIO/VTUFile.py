__doc__ = 'Handle output of fields and geometries to VTU'

from collections import namedtuple

class VTUFile:
    """ Constructor
        @param prefix: Prefix for files
        @type prefix: String
        @param binary: Dummy parameter to be compatible with VTF interface
        @type binary: (Ignored)
    """
    BlockInfo = namedtuple("BlockInfo", "dimension nodes elements fields")
    def __init__(self, prefix, binary=True):
        self._prefix = prefix
        self._data = {}
        self._step = 0

    def AddGeometryBlock(self, nodes, elements, blockid, dim):
      """ Add a geometry block to the VTU file
            @param nodes: Nodal coordinates (interleaved format)
            @type nodes: List of float
            @param elements: Element nodal connectivities
            @type elements: List of Integer
            @param blockid: ID of block being added
            @type blockid: Integer
            @param dim: Dimensionality of geometry block
            @type dim: Integer
            @return: Self
            @rtype: VTUFile
      """
      if not blockid in self._data:
        self._data[blockid] = VTUFile.BlockInfo(dim, nodes, elements, [])

      return self

    def AddGeometryDescriptor(self, size):
      """ Dummy method to be compatible with VTF interface
      """

    def AddGeometrySet(self, elems, block, name):
      """ Dummy method to be compatible with VTF interface
      """

    def AddField(self, coefs, block):
      """ Add a field to VTU file
          @param coefs: Field data
          @type coefs: List of float
          @param block: Block to map field to
          @type block: Integer
          @return: Block IDs
          @rtype: Tuple of (blockid, fieldid)
      """
      self._data[block].fields.append(("",coefs,-1))
      return (block, len(self._data[block].fields)-1)

    def AddFieldBlocks(self, blocks, name, comp=1, displacement=False):
      """ Add a description block for a field
          @param blocks: List with blocks corresponding to field
          @type blocks: List of tuple with (blockid, fieldid)
          @param name: Name of field
          @type name: String
          @param comp: Number of components in field
          @type comp: Integer
          @param displacement: Dummy parameter to be compatible with VTF interface
          @type displacement: (Ignored)
          @return: None
      """
      for block in blocks:
        self._data[block[0]].fields[block[1]] = (name, self._data[block[0]].fields[block[1]][1], comp)

    def AddState(self, time):
      """ Add a timestep to VTU file
          @param time: Time level
          @type time: Float
          @return: None
      """
      f = open('%s-%05i.vtu' %(self._prefix,self._step), 'w')
      f.write('<?xml version="1.0"?>\n')
      f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
      f.write('\t<UnstructuredGrid>\n')
      for idx in self._data:
        nperelm = pow(2, self._data[idx].dimension)
        f.write('\t\t<Piece NumberOfCells="%i" NumberOfPoints="%i">\n' % (len(self._data[idx].elements)/nperelm,len(self._data[idx].nodes)/3))
        self._WriteGeometry(f, self._data[idx])
        f.write('\t\t\t<PointData>\n')
        for field in self._data[idx].fields:
          # Disable writing of vector fields for now
          if field[2] == 1:
            self._WriteField(f, field)
        f.write('\t\t\t</PointData>')
        f.write('\t\t</Piece>\n')
        while self._data[idx].fields:
          del self._data[idx].fields[0]
      f.write('\t</UnstructuredGrid>\n')
      f.write('</VTKFile>')
      f.close()

      self._step = self._step + 1

    def _WriteGeometry(self, f, data):
      """ Write geometry piece to file
          @param f: File to write to
          @type f: File
          @param data: Block descriptor
          @type data: VTUFile.FieldInfo
      """
      nperelm = pow(2, data.dimension)
      f.write('\t\t\t<Points>\n')
      f.write('\t\t\t\t<DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="ascii">\n')
      for i in range(0, len(data.nodes)/3):
        f.write('\t\t\t\t\t%f %f %f\n' %(data.nodes[i*3],
                                         data.nodes[i*3+1],
                                         data.nodes[i*3+2]))
      f.write('\t\t\t\t</DataArray>\n')
      f.write('\t\t\t</Points>\n')
      f.write('\t\t\t<Cells>\n')
      f.write('\t\t\t\t<DataArray type="Int32" Name="connectivity" NumberOfComponents="1" format="ascii">\n')
      f.write('\t\t\t\t\t')
      for e in data.elements:
        f.write('%i ' %(e))
      f.write('\n')
      f.write('\t\t\t\t</DataArray>\n')
      if nperelm == 8:
        elmtype = 12
      if nperelm == 4:
        elmtype = 9
      if nperelm == 2:
        elmtype = 3
      f.write('\t\t\t\t<DataArray type="UInt8" Name="types" NumberOfComponents="1" format="ascii">\n')
      f.write('\t\t\t\t\t')
      for i in range(0, len(data.elements)/nperelm):
        f.write('%i ' %(elmtype))
      f.write('\n')
      f.write('\t\t\t\t</DataArray>\n')
      f.write('\t\t\t\t<DataArray type="Int32" Name="offsets" NumberOfComponents="1" format="ascii">\n')
      f.write('\t\t\t\t\t')
      for i in range(0, len(data.elements)/nperelm):
        f.write('%i ' %((i+1)*nperelm))
      f.write('\n')
      f.write('\t\t\t\t</DataArray>\n')
      f.write('\t\t\t</Cells>\n')

    def _WriteField(self, f, data):
      """ Write a field to file
          @param f: File to write to
          @type f: File
          @param data: Field descriptor
          @type data: List of tuple with (name, coefs, components)
          @return: None
      """
      f.write('\t\t\t\t<DataArray type="Float32" Name="%s" NumberOfComponents="%i" format="ascii">\n' %(data[0], data[2]))
      for coef in data[1]:
        if data[2] == 1:
          f.write('%f '%(coef))
        else:
          for d in coef:
            f.write('%f '%(d))
      f.write('\n')
      f.write('\t\t\t\t</DataArray>\n')

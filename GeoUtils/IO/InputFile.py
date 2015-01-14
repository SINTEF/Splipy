__doc__ == 'Class for working with IFEM input (.xinp) files'

from collections import namedtuple

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
                  remap = [1,2,3,4,5,6]
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

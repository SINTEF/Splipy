__doc__ == 'Class for working with IFEM input (.xinp) files'

import numpy as np
import os, subprocess, time, xml.dom.minidom, sys
from collections import namedtuple
from itertools import groupby, product
from operator import itemgetter
from GoTools import ReadG2
from GoTools.Preprocess import NaturalNodeNumbers

def progress(title, i, tot):
    """Progress output to stdout."""
    args = (title, i, tot, 100 * float(i) / tot)
    s = '\r%s: %i/%i (%.2f%%)' % args
    if i == tot:
        s += '\n'
    sys.stdout.write(s)
    sys.stdout.flush()

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

    # helper function
    def split_header( xml ):
      j = xml.find( '?>' )
      return (xml[:j+2],xml[j+2:]) if j>0 else ('',xml)

    # handle <include> tags
    for i in range(10): # max recursion depth
      elems = self.dom.getElementsByTagName('include')
      if not len(elems): break
      for elem in elems:
        try:
          include = xml.dom.minidom.parse( elem.firstChild.nodeValue ).firstChild
          elem.parentNode.replaceChild( include, elem )
        except xml.parsers.expat.ExpatError: # multiple root elements
          name = 'InputFile.%s.tmp' % time.strftime('%s')
          parent = elem.parentNode
          # Create copy and heal it
          f = open( elem.firstChild.nodeValue, 'r' )
          g = open( name, 'w' )
          g.write( '%s\n<xml>%s</xml>\n'%split_header(f.read()) )
          f.close()
          g.close()
          # Inject included xml
          parent.removeChild( elem )
          for include in xml.dom.minidom.parse( name ).firstChild.childNodes:
            if isinstance( include, xml.dom.minidom.Element ):
              parent.appendChild( include )
          # Remove copy
          subprocess.call( 'rm %s'%name, shell=True )

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
        typ = topset.getAttributeNode('type').nodeValue
        nam = topset.getAttributeNode('name').nodeValue
        if toptype and not typ == toptype or not nam == name: continue

        for item in topset.getElementsByTagName('item'):
          patch = int(item.getAttributeNode('patch').nodeValue)
          if not result.has_key(patch):
            result[patch] = InputFile.PatchInfo([], [], [])

          if not len(item.childNodes): continue # empty <item> tag
          values = map(int, item.childNodes[0].nodeValue.split())

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

  def writeOpenFOAM(self, name, display=False):
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

    # Find out how many faces in total
    def patch_num_faces(patch, faceid=0):
      nu, nv, nw = map(len, patch.GetKnots())
      if faceid == 0:
        return 3*nu*nv*nw - 2*(nu*nv + nu*nw + nv*nw) + nu + nv + nw
      elif faceid <= 2:
        return (nv-1)*(nw-1)
      elif faceid <= 4:
        return (nu-1)*(nw-1)
      return (nu-1)*(nv-1)

    internal_bnds = set((patchid, faceid)
                        for patchid in xrange(0, len(patchlist))
                        for faceid in xrange(1, 7))
    num_bnd_faces = 0
    for ts_name in self.GetTopologySets():
      ts = self.GetTopologySet(ts_name, convention='ifem', toptype='face')
      for patchid, patchinfo in ts.iteritems():
        internal_bnds.difference_update((patchid-1, faceid)
                                        for faceid in patchinfo.face)
        num_bnd_faces += sum(patch_num_faces(patchlist[patchid-1], faceid)
                             for faceid in patchinfo.face)
    internal_bnd_faces = sum(patch_num_faces(patchlist[p], faceid) for p, faceid in internal_bnds) / 2
    num_faces = sum(patch_num_faces(p) for p in patchlist) - internal_bnd_faces

    # Face-cell connectivity structure
    localid = lambda iu, iv, iw, nu, nv: iw*nu*nv + iv*nu + iu

    face_table = np.zeros((num_faces,),
                          dtype=[('owner', np.int32),
                                 ('neighbour', np.int32),
                                 ('boundary', np.int32),
                                 ('tup1', np.int32),
                                 ('tup2', np.int32),
                                 ('tup3', np.int32),
                                 ('tup4', np.int32)])
    face_table[:]['neighbour'] = -1
    face_table[:]['boundary'] = -1

    next_cell, next_face = 0, 0
    previous_faces = {}
    face_idxs = [(0,2,6,4), (1,5,7,3), (0,1,3,2), (4,6,7,5), (0,4,5,1), (2,3,7,6)]
    if display:
      progress('Processing faces', next_face, num_faces)
    for patchid, patch in enumerate(patchlist):
      nu, nv, nw = map(len, patch.GetKnots())
      localid_patch = lambda (iu, iv, iw): localid(iu, iv, iw, nu, nv)

      for iw, iv, iu in product(xrange(nw-1), xrange(nv-1), xrange(nu-1)):
        globalids = [nodenumbers[patchid][localid_patch(tp)]
                     for tp in [(iu,iv,iw), (iu+1,iv,iw), (iu,iv,iw+1), (iu+1,iv,iw+1),
                                (iu,iv+1,iw), (iu+1,iv+1,iw), (iu,iv+1,iw+1), (iu+1,iv+1,iw+1)]]
        for fidxs in face_idxs:
          face_ns = tuple(globalids[fid] for fid in fidxs)
          face_hash = frozenset(face_ns)

          if face_hash in previous_faces:
            face_id = previous_faces[face_hash]
            assert(face_table[face_id]['neighbour'] == -1)
            face_table[face_id]['neighbour'] = next_cell
            del previous_faces[face_hash]
          else:
            face_table[next_face]['owner'] = next_cell
            face_table[next_face]['tup1'] = face_ns[0]
            face_table[next_face]['tup2'] = face_ns[1]
            face_table[next_face]['tup3'] = face_ns[2]
            face_table[next_face]['tup4'] = face_ns[3]
            previous_faces[face_hash] = next_face
            next_face += 1

          if next_face % 100000 == 0 and display:
            progress('Processing faces', next_face, num_faces)

        next_cell += 1

    if display:
      progress('Processing faces', num_faces, num_faces)

    # Identify boundaries
    def get_faces(patchid, faceid):
      nu, nv, nw = map(len, patchlist[patchid].GetKnots())
      localid_patch = lambda iu, iv, iw: localid(iu, iv, iw, nu, nv)
      ret = []
      if faceid <= 2:
        iu = 0 if faceid == 1 else nu - 1
        return [frozenset((nodenumbers[patchid][localid_patch(iu,iv,iw)],
                           nodenumbers[patchid][localid_patch(iu,iv+1,iw)],
                           nodenumbers[patchid][localid_patch(iu,iv+1,iw+1)],
                           nodenumbers[patchid][localid_patch(iu,iv,iw+1)]))
                for iv in range(0, nv-1) for iw in range(0, nw-1)]
      if faceid <= 4:
        iv = 0 if faceid == 3 else nv - 1
        return [frozenset((nodenumbers[patchid][localid_patch(iu,iv,iw)],
                           nodenumbers[patchid][localid_patch(iu+1,iv,iw)],
                           nodenumbers[patchid][localid_patch(iu+1,iv,iw+1)],
                           nodenumbers[patchid][localid_patch(iu,iv,iw+1)]))
                for iu in range(0, nu-1) for iw in range(0, nw-1)]
      iw = 0 if faceid == 5 else nw - 1
      return [frozenset((nodenumbers[patchid][localid_patch(iu,iv,iw)],
                         nodenumbers[patchid][localid_patch(iu+1,iv,iw)],
                         nodenumbers[patchid][localid_patch(iu+1,iv+1,iw)],
                         nodenumbers[patchid][localid_patch(iu,iv+1,iw)]))
              for iu in range(0, nu-1) for iv in range(0, nv-1)]

    topologysets = self.GetTopologySets()

    for ts_id, ts_name in enumerate(topologysets):
      ts = self.GetTopologySet(ts_name, convention='ifem', toptype='face')
      for patchid, patchinfo in ts.iteritems():
        for faceid in patchinfo.face:
          faces = get_faces(patchid-1, faceid)
          for f in faces:
            face_id = previous_faces[f]
            face_table[face_id]['boundary'] = ts_id

    # Order faces
    face_table.sort(order=['boundary', 'owner', 'neighbour'])

    # Write faces, simultaneously compute boundary boundaries (heh)
    prev_boundary = -1
    bnd_bnds = []
    with open(os.path.join(name, 'faces'), 'w') as out:
      out.write(header('faceList', 'faces'))
      out.write('%i\n' % num_faces)
      out.write('(\n')
      for i in xrange(num_faces):
        out.write('4(%i %i %i %i)\n' % (face_table[i]['tup1'],
                                        face_table[i]['tup2'],
                                        face_table[i]['tup3'],
                                        face_table[i]['tup4']))
        if face_table[i]['boundary'] != prev_boundary:
          bnd_bnds.append(i)
          prev_boundary = face_table[i]['boundary']
      if face_table[-1]['boundary'] != -1:
        bnd_bnds.append(num_faces)
      out.write(')\n')

    # Write owners and neighbours
    note = 'nPoints: %i nCells: %i nFaces: %i nInternalFaces: %i' % (len(points),
                                                                     next_cell,
                                                                     num_faces,
                                                                     num_faces - num_bnd_faces)
    for thing in ['owner', 'neighbour']:
      with open(os.path.join(name, thing), 'w') as out:
        out.write(header('labelList', thing, note=note))
        out.write('%i\n' % num_faces)
        out.write('(\n')
        for i in xrange(num_faces):
          out.write('%i\n' % face_table[i][thing])
        out.write(')\n')

    # Write boundaries
    with open(os.path.join(name, 'boundary'), 'w') as out:
      out.write(header('polyBoundaryMesh', 'boundary'))
      out.write('%i\n' % len(topologysets))
      out.write('(\n')
      for ts_id, (start, end) in enumerate(zip(bnd_bnds[:-1], bnd_bnds[1:])):
        ts_name = topologysets[ts_id]
        out.write('%s\n' % ts_name)
        out.write('{\n')
        out.write('    type patch;\n')
        out.write('    nFaces %i;\n' % (end - start))
        out.write('    startFace %i;\n' % start)
        out.write('}\n')
      out.write(')\n')

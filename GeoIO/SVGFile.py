from GoTools import *
from GeoUtils.Factory import *
from math import *
import xml.etree.ElementTree as etree

def GetBezierRepresentation(curve):
    """ Compute a Bezier representation of a given spline curve. The input curve must be of order
        less than or equal to 4. The bezier representation is of order 4, and maximal knot multiplicity,
        i.e. consist of a C0-discretization with knot multiplicity 3.

        @param curve  : Spline curve
        @type  curve  : Curve
        @returns      : Bezier curve
        @rtype        : Curve
    """
    # error test input. Another way would be to approximate higher-order curves. Consider looking at Curve.Rebuild()
    if curve.GetOrder() > 4:
        raise RuntimeError('Bezier representation, only supported for curves of order 4 or less')

    bezier = curve.Clone()
    bezier.RaiseOrder(4-curve.GetOrder()) # make sure curve is cubic

    # make sure it is C0 everywhere (knot multiplicity of 3)
    knots        = bezier.GetKnots(True)
    thisKnot     = knots[0]
    multiplicity = 0
    newKnots     = []
    for i in range(len(knots)):
        if thisKnot == knots[i]:
            multiplicity += 1
        else:
            for j in range(3-multiplicity):
                newKnots.append(thisKnot)
            thisKnot     = knots[i]
            multiplicity = 1
    for k in newKnots:
        bezier.InsertKnot(k)

    return bezier

def WriteSVG_Curve(xmlRoot, curve, scale=1, center=[0,0], offset=[0,0], fill='none', stroke='#000000', width='2'):
    """ Writes a Curve to the xml tree. This will draw a single curve

        @param xmlRoot: XML Root node to append this surface to
        @type  xmlRoot: Etree.element
        @param curve  : The spline curve to write
        @type  curve  : Curve
        @param scale  : Scaling factor of image size divided by geometry size
        @type  scale  : Float
        @param center : Lower left corner of the bounding box, measured in geometry coordinates
        @type  center : Two floats
        @param offset : Lower left corner of drawing space, measured in image pixel coordinates
        @type  offset : Two floats
        @param fill   : Fill color in hex or none
        @type  fill   : String
        @param stroke : Line color written in hex, i.e. '#000000'
        @type  stroke : String
        @param width  : Line width, measured in pixels
        @type  width  : Int
        @returns: None
        @rtype  : NoneType
    """
    curveNode = etree.SubElement(xmlRoot, 'path')
    curveNode.attrib['style'] = 'fill:%s;stroke:%s;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' %(fill,stroke,width)
    bezier = GetBezierRepresentation(curve)
    bezier = Translate(bezier, [-x for x in center])
    bezier = Scale(    bezier, scale)
    bezier = Translate(bezier, offset)
    pathString = 'M %f,%f C ' % (bezier[0][0], bezier[0][1])

    for i in range(1,len(bezier)):
        pathString += '%f,%f ' % (bezier[i][0],bezier[i][1])

    curveNode.attrib['d'] = pathString

def WriteSVG_Surface(xmlRoot, surface, scale=1, center=[0,0], offset=[0,0], fill='#ffcc99'):
    """ Writes a Surface to the xml tree. This will draw the surface along with all knot lines

        @param xmlRoot: XML Root node to append this surface to
        @type  xmlRoot: Etree.element
        @param surface: The spline surface to write
        @type  surface: Surface
        @param scale  : Scaling factor of image size divided by geometry size
        @type  scale  : Float
        @param center : Lower left corner of the bounding box, measured in geometry coordinates
        @type  center : Two floats
        @param offset : Lower left corner of drawing space, measured in image pixel coordinates
        @type  offset : Two floats
        @param fill   : Surface color written in hex, i.e. '#ffcc99'
        @type  fill   : String
        @returns: None
        @rtype  : NoneType
    """
    # fetch boundary curves and create a connected, oriented bezier loop from it
    bndryCurves = surface.GetEdges()
    crv1        = CombineCurves(bndryCurves[0], bndryCurves[1])
    crv2        = CombineCurves(crv1,           bndryCurves[2].FlipParametrization())
    boundary    = CombineCurves(crv2,           bndryCurves[3].FlipParametrization())
    
    # fetch all meshlines (i.e. elements, also known as knot spans)
    knot = surface.GetKnots()
    knotlines = []
    for k in knot[0][1:-1]:
        knotlines.append(surface.GetConstParCurve(k, 0))
    for k in knot[1][1:-1]:
        knotlines.append(surface.GetConstParCurve(k, 1))

    # create a group node for all elements corresponding to this surface patch
    groupNode = etree.SubElement(xmlRoot, 'g')
    
    # fill interior with a peach color
    WriteSVG_Curve(groupNode, boundary, scale, center, offset, fill, width=2)

    # draw all meshlines
    for meshline in knotlines:
        WriteSVG_Curve(groupNode, meshline, scale, center, offset, width=1)


def WriteSVG(filename, entities, width=1000, height=1000, margin=0.05):
  """ Writes a list of planar curves and surfaces to vector graphics SVG file.
      The image will never be stretched, and the geometry will keep width/height ratio
      of original geometry, regardless of provided width/height ratio from arguments.

      @param filename: Filename to write results to
      @type  filename: String
      @param entities: Primitives to write
      @type  entities: List of Curves and/or Surfaces
      @param width   : Maximum image width in pixels
      @type  width   : Int
      @param height  : Maximum image height in pixels 
      @type  height  : Int
      @param margin  : White-space around all edges of image, given in percentage of total size (default 5%)
      @type  margin  : Float
      @returns: None
      @rtype  : NoneType
  """
    if GetDimension() != 2:
        raise RuntimeError('SVG files only applicable for 2D geometries')

    # in case of single-entry input, wrap it in a list so we can iterate on it
    if type(entities) is not list:
        entities = [entities]

    # compute the bounding box for all geometries 
    boundingbox = [1e10, 1e10, -1e10, -1e10]
    for entry in entities:
        for controlpoint in entry:
            boundingbox[0] = min(boundingbox[0], controlpoint[0])
            boundingbox[1] = min(boundingbox[1], controlpoint[1])
            boundingbox[2] = max(boundingbox[2], controlpoint[0])
            boundingbox[3] = max(boundingbox[3], controlpoint[1])

    # compute scaling factors by keeping aspect ratio, and never exceed width or height size (including margins)
    geometryRatio = (boundingbox[3]-boundingbox[1])/(boundingbox[2]-boundingbox[0])
    imageRatio    = 1.0*height / width
    if geometryRatio > imageRatio: # scale by y-coordinate
        marginPixels = height*margin
        scale        = height*(1-2*margin) / (boundingbox[3]-boundingbox[1])
        width        = height/geometryRatio
    else:                          # scale by x-coordinate
        marginPixels = width*margin
        scale        = width*(1-2*margin) / (boundingbox[2]-boundingbox[0])
        height       = width*geometryRatio
    center           = [boundingbox[0], boundingbox[1]]
    offset           = [marginPixels, marginPixels]


    # create xml root tag
    xmlRoot = etree.Element('svg',  {'xmlns':'http://www.w3.org/2000/svg', 'version':'1.1', 'width':str(width), 'height':str(height)})

    # populate tree with all curves and surfaces in entities
    for entry in entities:
        if type(entry) is Curve:
            WriteSVG_Curve(xmlRoot, entry, scale, center, offset)
        elif type(entry) is Surface:
            WriteSVG_Surface(xmlRoot, entry, scale, center, offset)

    # write to file
    xmlTree = etree.ElementTree(xmlRoot)
    xmlTree.write(filename, encoding='us-ascii', xml_declaration=True)

  

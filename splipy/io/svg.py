from __future__ import print_function

from splipy import Curve, Surface, SplineObject, BSplineBasis
import xml.etree.ElementTree as etree
from xml.dom import minidom
import numpy as np
import re
from .master import MasterIO

def read_number_and_unit(mystring):
    try:
        for i in range(1, len(mystring)):
            number=float(mystring[:i])
    except ValueError:
        unit=mystring[i-1:]
    return (number, unit)


def bezier_representation(curve):
    """ Compute a Bezier representation of a given spline curve. The input
        curve must be of order less than or equal to 4. The bezier
        representation is of order 4, and maximal knot multiplicity, i.e.
        consist of a C0-discretization with knot multiplicity 3.

        :param curve  : Spline curve
        :type  curve  : Curve
        :returns      : Bezier curve
        :rtype        : Curve
    """
    # error test input. Another way would be to approximate higher-order curves. Consider looking at Curve.rebuild()
    if curve.order(0) > 4 or curve.rational:
        raise RuntimeError('Bezier representation, only supported for non-rational curves of order 4 or less')

    bezier = curve.clone()
    bezier.raise_order(4-curve.order(0)) # make sure curve is cubic

    # make it non-periodic
    if bezier.periodic():
        bezier = bezier.split(bezier.start(0))

    # make sure it is C0 everywhere
    for k in bezier.knots(0):
        bezier.insert_knot( [k]*bezier.continuity(k) )

    return bezier


class SVG(MasterIO):

    namespace = '{http://www.w3.org/2000/svg}'

    def __init__(self, filename, width=1000, height=1000, margin=0.05):
        """  Constructor
          :param filename: Filename to write results to
          :type  filename: String
          :param width   : Maximum image width in pixels
          :type  width   : Int
          :param height  : Maximum image height in pixels
          :type  height  : Int
          :param margin  : White-space around all edges of image, given in percentage of total size (default 5%)
          :type  margin  : Float
        """
        if filename[-4:] != '.svg':
            filename += '.svg'
        self.filename = filename

        self.width  = width
        self.height = height
        self.margin = margin

        self.all_objects = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # in case something goes wrong, print error and abandon writing process
        if exc_type is not None:
            print(exc_type, exc_value, traceback)
            return False

        # compute the bounding box for all geometries
        boundingbox = [np.inf, np.inf, -np.inf, -np.inf]
        for entry in self.all_objects:
            bb = entry.bounding_box()
            boundingbox[0] = min(boundingbox[0], bb[0][0])
            boundingbox[1] = min(boundingbox[1], bb[1][0])
            boundingbox[2] = max(boundingbox[2], bb[0][1])
            boundingbox[3] = max(boundingbox[3], bb[1][1])

        # compute scaling factors by keeping aspect ratio, and never exceed width or height size (including margins)
        geometryRatio = float(boundingbox[3]-boundingbox[1])/(boundingbox[2]-boundingbox[0])
        imageRatio    = 1.0*self.height / self.width
        if geometryRatio > imageRatio: # scale by y-coordinate
            marginPixels = self.height*self.margin
            self.scale   = self.height*(1-2*self.margin) / (boundingbox[3]-boundingbox[1])
            self.width   = self.height/geometryRatio + 2*marginPixels
        else:                          # scale by x-coordinate
            marginPixels = self.width*self.margin
            self.scale   = self.width*(1-2*self.margin) / (boundingbox[2]-boundingbox[0])
            self.height  = self.width*geometryRatio + 2*marginPixels
        self.center      = [boundingbox[0], boundingbox[1]]
        self.offset      = [marginPixels, marginPixels]

        # create xml root tag
        self.xmlRoot = etree.Element('svg',  {'xmlns':'http://www.w3.org/2000/svg',
                                              'version':'1.1',
                                              'width':str(self.width),
                                              'height':str(self.height)})

        # populate tree with all curves and surfaces in entities
        for entry in self.all_objects:
            if isinstance(entry, Curve):
                self.write_curve(self.xmlRoot, entry)
            elif isinstance(entry, Surface):
                self.write_surface(entry)

        # if no objects are stored, then we've most likely only called read()
        if len(self.all_objects) > 0:
            rough_string = etree.tostring(self.xmlRoot) # entire xml-file on one line
            reparsed     = minidom.parseString(rough_string)
            result       = reparsed.toprettyxml(indent="  ") # adds newline and inline
            f = open(self.filename, 'w')
            f.write(result)

    def write_curve(self, xmlNode, curve, fill='none', stroke='#000000', width=2):
        """ Writes a Curve to the xml tree. This will draw a single curve

            :param xmlNode: Node in xml tree
            :type  xmlNode: Etree.element
            :param curve  : The spline curve to write
            :type  curve  : Curve
            :param fill   : Fill color in hex or none
            :type  fill   : String
            :param stroke : Line color written in hex, i.e. '#000000'
            :type  stroke : String
            :param width  : Line width, measured in pixels
            :type  width  : Int
            :returns: None
            :rtype  : NoneType
        """
        curveNode = etree.SubElement(xmlNode, 'path')
        curveNode.attrib['style'] = 'fill:%s;stroke:%s;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' %(fill,stroke,width)
        bezier  = bezier_representation(curve)
        bezier -= self.center
        bezier *= self.scale
        bezier += self.offset
        pathString = 'M %f,%f C ' % (bezier[0][0],  self.height + 2*self.margin - bezier[0][1])

        for i in range(1,len(bezier)):
            pathString += '%f,%f ' % (bezier[i][0], self.height + 2*self.margin - bezier[i][1])

        curveNode.attrib['d'] = pathString

    def write_surface(self, surface, fill='#ffcc99'):
        """ Writes a Surface to the xml tree. This will draw the surface along with all knot lines

            :param surface: The spline surface to write
            :type  surface: Surface
            :param fill   : Surface color written in hex, i.e. '#ffcc99'
            :type  fill   : String
            :returns: None
            :rtype  : NoneType
        """
        # fetch boundary curves and create a connected, oriented bezier loop from it
        bndry_curves = surface.edges()
        bndry_curves[0].reverse()
        bndry_curves[3].reverse()
        boundary      = bndry_curves[0]
        boundary.append(bndry_curves[2])
        boundary.append(bndry_curves[1])
        boundary.append(bndry_curves[3])

        # fetch all meshlines (i.e. elements, also known as knot spans)
        knot = surface.knots()
        knotlines = []
        for k in knot[0][1:-1]:
            knotlines.append(surface.const_par_curve(k, 0))
        for k in knot[1][1:-1]:
            knotlines.append(surface.const_par_curve(k, 1))

        # create a group node for all elements corresponding to this surface patch
        groupNode = etree.SubElement(self.xmlRoot, 'g')

        # fill interior with a peach color
        self.write_curve(groupNode, boundary, fill, width=2)

        # draw all meshlines
        for meshline in knotlines:
            self.write_curve(groupNode, meshline, width=1)


    def write(self, obj):
        """ Writes a list of planar curves and surfaces to vector graphics SVG file.
            The image will never be stretched, and the geometry will keep width/height ratio
            of original geometry, regardless of provided width/height ratio from arguments.

            :param obj: Primitives to write
            :type  obj: List of Curves and/or Surfaces
            :returns: None
            :rtype  : NoneType
        """
        # actually this is a dummy method. It will collect all geometries provided
        # and ton't actually write them to file until __exit__ is called

        if isinstance(obj[0], SplineObject): # input SplineModel or list
            for o in obj:
                self.write(o)
            return

        if obj.dimension != 2:
            raise RuntimeError('SVG files only applicable for 2D geometries')

        # have to clone stuff we put here, in case they change on the outside
        self.all_objects.append( obj.clone() )

    def read(self):
        tree = etree.parse(self.filename)
        root = tree.getroot()
        self.width,_  = read_number_and_unit(root.attrib['width'])
        self.height,_ = read_number_and_unit(root.attrib['height'])
        result = []
        for path in root.iter(SVG.namespace + 'path'):
            result.append(self.curve_from_path(path.attrib['d']))
        return result

    def curve_from_path(self, path):
        # see http://www.w3schools.com/svg/svg_path.asp for documentation
        # and also https://www.w3.org/TR/SVG/paths.html

        # figure out the largest polynomial order of this path
        if re.search('[cCsS]', path):
            order = 4
        elif re.search('[qQtTaA]', path):
            order = 3
        else:
            order = 2
        result = None

        # each 'piece' is an operator (M,C,Q,L etc) and accomponying list of argument points
        for piece in re.findall('[a-zA-Z][^a-zA-Z]*', path):

            # if not single-letter command (i.e. 'z')
            if len(piece)>1:
                # points is a (string-)list of (x,y)-coordinates for the given operator
                points = re.findall('-?\d+\.?\d*', piece[1:])

                # convert string-list to a list of numpy arrays (of size 2)
                np_pts = np.reshape(np.array(points).astype('float'), (int(len(points)/2),2))

            if piece[0] == 'm' or piece[0] == 'M':
                # I really hope it always start with a move command (think it does)
                startpoint = np_pts[0]
                if len(np_pts) > 1:
                    if piece[0] == 'M':
                        knot = [0] + list(range(len(np_pts))) + [len(np_pts)-1]
                        curve_piece = Curve(BSplineBasis(2, knot), np_pts)
                    elif piece[0] == 'm':
                        knot = [0] + list(range(len(np_pts))) + [len(np_pts)-1]
                        controlpoints = [startpoint]
                        for cp in np_pts[1:]:
                            controlpoints.append(cp + controlpoints[-1])
                        curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
                else:
                    continue
            elif piece[0] == 'c':
                # cubic spline, relatively positioned
                controlpoints = [startpoint]
                knot = list(range(int(len(np_pts)/3)+1)) * 3
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    startpoint = controlpoints[int((len(controlpoints)-1)/3)*3]
                    controlpoints.append(cp + startpoint)
                curve_piece = Curve(BSplineBasis(4, knot), controlpoints)
            elif piece[0] == 'C':
                # cubic spline, absolute position
                controlpoints = [startpoint]
                knot = list(range(int(len(np_pts)/3)+1)) * 3
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    controlpoints.append(cp)
                curve_piece = Curve(BSplineBasis(4, knot), controlpoints)
            elif piece[0] == 'q':
                # quadratic spline, relatively positioned
                controlpoints = [startpoint]
                knot = list(range(int(len(np_pts)/2)+1)) * 2
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    startpoint = controlpoints[int((len(controlpoints)-1)/2)*2]
                    controlpoints.append(cp + startpoint)
                curve_piece = Curve(BSplineBasis(3, knot), controlpoints)
            elif piece[0] == 'Q':
                # quadratic spline, absolute position
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)/2+1)) * 2
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    controlpoints.append(cp)
                curve_piece = Curve(BSplineBasis(3, knot), controlpoints)
            elif piece[0] == 'l':
                # linear spline, relatively positioned
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)+1))
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    startpoint = controlpoints[-1]
                    controlpoints.append(cp + startpoint)
                curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
            elif piece[0] == 'L':
                # linear spline, absolute position
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)+1))
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    controlpoints.append(cp)
                curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
            elif piece[0] == 'z' or piece[0] == 'Z':
                # periodic curve
                curve_piece = Curve(BSplineBasis(2), [startpoint, result[0]])
                result.append(curve_piece.make_periodic(0))
                continue
            else:
                raise RuntimeError('Unknown path parameter:' + piece)

            if result is None:
                result = curve_piece
            else:
                result.append(curve_piece)
            startpoint = controlpoints[-1]

        # invert y-values since these are image coordinates
        result[:,1] = self.height - result[:,1]
        return result

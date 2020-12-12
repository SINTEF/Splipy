import xml.etree.ElementTree as etree
from xml.dom import minidom
import re

import numpy as np

from ..curve import Curve
from ..surface import Surface
from ..splineobject import SplineObject
from ..basis import BSplineBasis
from .. import curve_factory, state

from .master import MasterIO


def read_number_and_unit(mystring):
    unit = ''
    try:
        for i in range(1, len(mystring)+1):
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
        parent_map = dict((c, p) for p in tree.iter() for c in p)
        if 'width' in root.attrib:
            self.width,_  = read_number_and_unit(root.attrib['width'])
            self.height,_ = read_number_and_unit(root.attrib['height'])
        result = []
        for path in root.iter(SVG.namespace + 'path'):
            crvs = self.curves_from_path(path.attrib['d'])
            parent = path
            while parent != root:
                if 'transform' in parent.attrib:
                    for crv in crvs:
                        self.transform(crv, parent.attrib['transform'])
                parent = parent_map[parent]

            # invert y-values since these are image coordinates
            for crv in crvs:
                crv *= [1,-1]
                crv += [0, self.height]
                result.append(crv)
        return result

    def transform(self, curve, operation):
        # intended input operation string: 'translate(-10,-20) scale(2) rotate(45) translate(5,10)'
        all_operations = re.findall(r'[^\)]*\)', operation.lower())
        all_operations.reverse()
        for one_operation in all_operations:
            parts = re.search(r'([a-z]*)\w*\((.*)\)', one_operation.strip())
            func = parts.group(1)
            args = [float(d) for d in parts.group(2).split(',')]
            if func == 'translate':
                if len(args)==1:
                    args.append(0)
                curve += args
            elif func == 'scale':
                if len(args)==1:
                    args.append(args[0])
                curve *= args
            elif func == 'rotate':
                curve.rotate(-args[0]/360*2*np.pi)
            elif func == 'matrix':
                M = np.array([[args[0], args[2], args[4]],
                               [args[1], args[3], args[5]],
                               [      0,       0,      1]])
                n = len(curve)
                if not curve.rational:
                    cp = np.ones((n, 3))  # pad with weights=1
                    cp[:, :-1] = np.reshape(curve.controlpoints, (n, 2))
                    cp = cp @ M.T
                    curve.controlpoints = np.reshape(np.array(cp[:,:-1]), curve.controlpoints.shape)
                else:
                    cp = np.reshape(curve.controlpoints, (n, 3))
                    cp = cp @ M.T
                    curve.controlpoints = np.reshape(np.array(cp), curve.controlpoints.shape)

    def curves_from_path(self, path):
        # see https://www.w3schools.com/graphics/svg_path.asp for documentation
        # and also https://www.w3.org/TR/SVG/paths.html

        # figure out the largest polynomial order of this path
        if re.search('[cCsS]', path):
            order = 4
        elif re.search('[qQtTaA]', path):
            order = 3
        else:
            order = 2
        last_curve = None
        result = []

        # each 'piece' is an operator (M,C,Q,L etc) and accomponying list of argument points
        for piece in re.findall('[a-zA-Z][^a-zA-Z]*', path):

            # if not single-letter command (i.e. 'z')
            if len(piece)>1:
                # points is a (string-)list of (x,y)-coordinates for the given operator
                points = re.findall(r'-?\d+\.?\d*', piece[1:])

                if piece[0].lower() != 'a' and piece[0].lower() != 'v' and piece[0].lower() != 'h':
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
            elif piece[0] == 's':
                # smooth cubic spline, relative position
                controlpoints = [startpoint]
                knot = list(range(int(len(np_pts)/2)+1)) * 3
                knot += [knot[0], knot[-1]]
                knot.sort()
                x0  = np.array(last_curve[-1])
                xn1 = np.array(last_curve[-2])
                controlpoints.append(2*x0 -xn1)
                startpoint = controlpoints[-1]
                for i, cp in enumerate(np_pts):
                    if i % 2 == 0 and i>0:
                        startpoint = controlpoints[-1]
                        controlpoints.append(2*controlpoints[-1] - controlpoints[-2])
                    controlpoints.append(cp + startpoint)
                curve_piece = Curve(BSplineBasis(4, knot), controlpoints)
            elif piece[0] == 'S':
                # smooth cubic spline, absolute position
                controlpoints = [startpoint]
                knot = list(range(int(len(np_pts)/2)+1)) * 3
                knot += [knot[0], knot[-1]]
                knot.sort()
                x0  = np.array(last_curve[-1])
                xn1 = np.array(last_curve[-2])
                controlpoints.append(2*x0 -xn1)
                for i,cp in enumerate(np_pts):
                    if i % 2 == 0 and i>0:
                        controlpoints.append(2*controlpoints[-1] - controlpoints[-2])
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
            elif piece[0] == 'h':
                # horizontal piece, relatively positioned
                np_pts = np.array(points).astype('float')
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)+1))
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    startpoint = controlpoints[-1]
                    controlpoints.append(np.array([cp, 0]) + startpoint)
                curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
            elif piece[0] == 'H':
                # horizontal piece, absolute position
                np_pts = np.array(points).astype('float')
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)+1))
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    controlpoints.append([cp, startpoint[1]])
                curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
            elif piece[0] == 'v':
                # vertical piece, relatively positioned
                np_pts = np.array(points).astype('float')
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)+1))
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    startpoint = controlpoints[-1]
                    controlpoints.append(np.array([0, cp]) + startpoint)
                curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
            elif piece[0] == 'V':
                # vertical piece, absolute position
                np_pts = np.array(points).astype('float')
                controlpoints = [startpoint]
                knot = list(range(len(np_pts)+1))
                knot += [knot[0], knot[-1]]
                knot.sort()
                for cp in np_pts:
                    controlpoints.append([startpoint[0], cp])
                curve_piece = Curve(BSplineBasis(2, knot), controlpoints)
            elif piece[0] == 'A' or piece[0] == 'a':
                np_pts = np.reshape(np.array(points).astype('float'), (int(len(points))))
                rx              = float(points[0])
                ry              = float(points[1])
                x_axis_rotation = float(points[2])
                large_arc_flag  = (points[3] != '0')
                sweep_flag      = (points[4] != '0')
                xend            = np.array([float(points[5]), float(points[6]) ])
                if piece[0] == 'a':
                    xend += startpoint

                R = np.array([[ np.cos(x_axis_rotation), np.sin(x_axis_rotation)],
                              [-np.sin(x_axis_rotation), np.cos(x_axis_rotation)]])
                xp = np.linalg.solve(R, (startpoint - xend)/2)
                if sweep_flag == large_arc_flag:
                    cprime = -(np.sqrt(abs(rx**2*ry**2 - rx**2*xp[1]**2 - ry**2*xp[0]**2) /
                                       (rx**2*xp[1]**2 + ry**2*xp[0]**2)) *
                                       np.array([rx*xp[1]/ry, -ry*xp[0]/rx]))
                else:
                    cprime = +(np.sqrt(abs(rx**2*ry**2 - rx**2*xp[1]**2 - ry**2*xp[0]**2) /
                                       (rx**2*xp[1]**2 + ry**2*xp[0]**2)) *
                                       np.array([rx*xp[1]/ry, -ry*xp[0]/rx]))
                center = np.linalg.solve(R.T, cprime) + (startpoint+xend)/2
                def arccos(vec1, vec2):
                    return (np.sign(vec1[0]*vec2[1] - vec1[1]*vec2[0]) *
                            np.arccos(vec1.dot(vec2)/np.linalg.norm(vec1)/np.linalg.norm(vec2)))
                tmp1 = np.divide( xp - cprime, [rx,ry])
                tmp2 = np.divide(-xp - cprime, [rx,ry])
                theta1 = arccos(np.array([1,0]), tmp1)
                delta_t= arccos(tmp1, tmp2) % (2*np.pi)
                if not sweep_flag and delta_t > 0:
                    delta_t -= 2*np.pi
                elif sweep_flag and delta_t < 0:
                    delta_t += 2*np.pi
                curve_piece = (curve_factory.circle_segment(delta_t)*[rx,ry]).rotate(theta1) + center
                # curve_piece = curve_factory.circle_segment(delta_t)

            elif piece[0] == 'z' or piece[0] == 'Z':
                # periodic curve
                # curve_piece = Curve(BSplineBasis(2), [startpoint, last_curve[0]])
                # curve_piece.reparam([0, curve_piece.length()])
                # last_curve.append(curve_piece).make_periodic(0)
                last_curve.make_periodic(0)
                result.append(last_curve)
                last_curve = None
                continue
            else:
                raise RuntimeError('Unknown path parameter:' + piece)

            if(curve_piece.length()>state.controlpoint_absolute_tolerance):
                curve_piece.reparam([0, curve_piece.length()])
                if last_curve is None:
                    last_curve = curve_piece
                else:
                    last_curve.append(curve_piece)
            startpoint = last_curve[-1,:2] # disregard rational weight (if any)

        if last_curve is not None:
            result.append(last_curve)
        return result

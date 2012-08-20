#include "curvefactory.h"

#include "curve.h"
#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "sislP.h"
// #include "GoTools/geometry/SISL_code.h"

extern "C" {
PyObject* CurveFactory_module;

PyDoc_STRVAR(generate_approximate_curve__doc__, "Construct a spline curve that approximatively interpolates a point cloud\n"
                                                "@param points: The point cloud\n"
                                                "@type points: List of (Point, list of floats or tuple of floats)\n"
                                                "@param parvals: The parameter values for the points\n"
                                                "@type parvals: List of floats\n"
                                                "@param order: (optional) The order of the spline curve (p+1)\n"
                                                "@type order: int >= 1\n"
                                                "@param maxiter: (optional) The maxium number of iterations in interpolation scheme\n"
                                                "@type maxiter: int >= 1\n"
                                                "@return: Spline interpolating the point cloud\n"
                                                "@rtype: Curve");
PyObject* Generate_ApproximateCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"points", "parvals", "order", "maxiter", NULL };
  PyObject* pointso;
  PyObject* parvals;
  int order=3;
  int max_iter=5;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO|ii",
                                   (char**)keyWords,&pointso,&parvals,
                                                    &order,&max_iter))
    return NULL;

  if (!PyObject_TypeCheck(pointso,&PyList_Type) || 
      !PyObject_TypeCheck(parvals,&PyList_Type))
    return NULL;

  // get points
  std::vector<double> points;
  for (int i=0; i < PyList_Size(pointso); ++i) {
    PyObject* entryo = PyList_GetItem(pointso,i);
    shared_ptr<Go::Point> entry = PyObject_AsGoPoint(entryo);
    if (entry) {
      points.push_back((*entry)[0]);
      points.push_back((*entry)[1]);
      if (modState.dim == 3)
        points.push_back((*entry)[2]);
    }
  }

  // get parameters
  std::vector<double> params;
  for (int i=0; i < PyList_Size(parvals); ++i)
    params.push_back(PyFloat_AsDouble(PyList_GetItem(parvals,i)));

  Go::ApproxCurve approx(points,params,modState.dim,modState.approxTolerance,
                         params.size(),order);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  double maxdist;
  double avdist;
  result->data = approx.getApproxCurve(maxdist,avdist,max_iter);

  return (PyObject*)result;
}


PyDoc_STRVAR(generate_circle__doc__, "Generate a circle\n"
                                     "@param center: The center of the circle\n"
                                     "@type center: Point, list of floats or tuple of floats\n"
                                     "@param radius: The radius of the circle\n"
                                     "@type radius: float\n"
                                     "@param normal: The normal of the circle\n"
                                     "@type normal: Point, list of floats or tuple of floats\n"
                                     "@return: The circle\n"
                                     "@rtype: Curve");
PyObject* Generate_Circle(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "radius", "normal", NULL };
  PyObject* centero;
  PyObject* normalo = 0;
  double radius;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od|O",
                                   (char**)keyWords,&centero,&radius,&normalo))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  if (!center)
    return NULL;

  Go::Point normal(0.0,0.0);
  if (modState.dim == 3) {
    shared_ptr<Go::Point> norm = PyObject_AsGoPoint(normalo);
    if (!norm)
      return NULL;
    normal = *norm;
  }

  Go::Point x_axis;
  if (modState.dim == 3) {
    Go::Point unit;
    if (fabs(normal[0]) < fabs(normal[1]) && fabs(normal[0]) < fabs(normal[2]))
      unit = Go::Point(1.0, 0.0, 0.0);
    else if (fabs(normal[1]) < fabs(normal[2]))
      unit = Go::Point(0.0, 1.0, 0.0);
    else
      unit = Go::Point(0.0, 0.0, 1.0);
    x_axis = normal % unit;
  } else
    x_axis = Go::Point(1.0,0.0);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::Circle(radius,*center,normal,x_axis));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_circle_segment__doc__, "Generate a circle segment\n"
                                             "@param center: The center of the circle\n"
                                             "@type center: Point, list of floats or tuple of floats\n"
                                             "@param start: The starting point\n"
                                             "@type start: Point, list of floats or tuple of floats\n"
                                             "@param angle: The angle of the circle segment in radians\n"
                                             "@type angle: float\n"
                                             "@param normal: The normal of the circle\n"
                                             "@type normal: Point, list of floats or tuple of floats\n"
                                             "@return: The circle segment\n"
                                             "@rtype: Curve");
PyObject* Generate_CircleSegment(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "start", "angle", "normal", NULL };
  PyObject* centero;
  PyObject* starto;
  PyObject* normalo = 0;
  double angle;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOd|O",
                            (char**)keyWords,&centero,&starto,&angle,&normalo))
      return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> start = PyObject_AsGoPoint(starto);
  if (!center || !start)
    return NULL;
  Go::Point normal(0.0,0.0);
  if (modState.dim == 3) {
    shared_ptr<Go::Point> norm = PyObject_AsGoPoint(normalo);
    if (!norm)
      return NULL;
    normal = *norm;
  }

  Go::Point x_axis = *start-*center;
  Go::Point normal2 = normal-((x_axis*normal)/x_axis.length2())*x_axis;
  double radius = x_axis.length();

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::Circle(radius,*center,normal2,x_axis));
  static_pointer_cast<Go::Circle>(result->data)->setParamBounds(0.0,angle);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_ellipse__doc__, "Generate an ellipse\n"
                                      "@param center: Center of ellipse\n"
                                      "@type center: Point, list of floats or tuple of floats\n"
                                      "@param semi_axis: The semi-axis of the ellipse\n"
                                      "@type semi_axis: Point, list of floats or tuple of floats\n"
                                      "@param radius1: The first radius of the ellipse\n"
                                      "@type radius1: float\n"
                                      "@param radius2: The second radius of the ellipse\n"
                                      "@type radius2: float\n"
                                      "@param normal: Normal of ellipse\n"
                                      "@type normal: Point, list of floats or tuple of floats\n"
                                      "@return: The ellipse\n"
                                      "@rtype: Curve");
PyObject* Generate_Ellipse(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "semi_axis", "radius1", "radius2", "normal", NULL };
  PyObject* centero;
  PyObject* axiso;
  PyObject* normalo = 0;
  double radius1, radius2;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOdd|O",
                                   (char**)keyWords,&centero,&axiso,
                                   &radius1,&radius2,&normalo))
      return NULL;

  if (centero->ob_type != &Point_Type || axiso->ob_type != &Point_Type)
    return NULL;
  if (modState.dim == 3 && normalo->ob_type != &Point_Type)
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;
  Go::Point normal(0.0,0.0);
  if (modState.dim == 3) {
    shared_ptr<Go::Point> norm = PyObject_AsGoPoint(normalo);
    if (!norm)
      return NULL;
    normal = *norm;
  }

  Go::Point normal2 = normal-((*axis*normal)/axis->length2())*(*axis);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::Ellipse(*center,*axis,normal2,radius1,radius2));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_elliptic_segment__doc__, "Generate an elliptic segment\n"
                                               "@param center: Center of ellipse\n"
                                               "@type center: Point, list of floats or tuple of floats\n"
                                               "@param semi_axis: The semi-axis of the ellipse\n"
                                               "@type semi_axis: Point, list of floats or tuple of floats\n"
                                               "@param radius1: The first radius of the ellipse\n"
                                               "@type radius1: float\n"
                                               "@param radius2: The second radius of the ellipse\n"
                                               "@type radius2: float\n"
                                               "@param startangle: The starting angle of the elliptic segment\n"
                                               "@type startangle: float\n"
                                               "@param endangle: The ending angle of the elliptic segment\n"
                                               "@type endangle: float\n"
                                               "@param normal: Normal of ellipse\n"
                                               "@type normal: Point, list of floats or tuple of floats\n"
                                               "@return: The ellipse\n"
                                               "@rtype: Curve");
PyObject* Generate_EllipticSegment(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "semi_axis", "radius1", "radius2",
                                   "startangle", "endangle", "normal", NULL };
  PyObject* centero;
  PyObject* axiso;
  PyObject* normalo = 0;
  double radius1, radius2, startangle, endangle;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOdddd|O",
                                   (char**)keyWords,&centero,&axiso,&radius1,
                                   &radius2,&startangle,&endangle,&normalo))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;
  Go::Point normal(0.0,0.0);
  if (modState.dim == 3) {
    shared_ptr<Go::Point> norm = PyObject_AsGoPoint(normalo);
    if (!norm)
      return NULL;
    normal = *norm;
  }

  Go::Point normal2 = normal-((*axis*normal)/axis->length2())*(*axis);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::Ellipse(*center,*axis,normal2,radius1,radius2));
  static_pointer_cast<Go::Ellipse>(result->data)->setParamBounds(startangle,endangle);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_helix__doc__, "Generate a helix\n"
                                    "@param center: Center of helix\n"
                                    "@type center: Point, list of floats or tuple of floats\n"
                                    "@param start: Starting point of helix\n"
                                    "@type start: Point, list of floats or tuple of floats\n"
                                    "@param axis: Center axis of helix\n"
                                    "@type axis: Point, list of floats or tuple of floats\n"
                                    "@param frequency: Frequency of helix\n"
                                    "@type frequency: float\n"
                                    "@param quadrants: Number of quadrants in helix\n"
                                    "@type quadrants: 1 <= int <= 4\n"
                                    "@return: The helix\n"
                                    "@rtype: Curve");
PyObject* Generate_Helix(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "start", "axis", "frequency", "quadrants", NULL };
  PyObject* centero;
  PyObject* starto;
  PyObject* axiso;
  double frequency;
  int quadrants;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOdi",
                                   (char**)keyWords,&centero,&starto,
                                   &axiso,&frequency,&quadrants))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  shared_ptr<Go::Point> start = PyObject_AsGoPoint(starto);

  if (!center || !axis || !start)
    return NULL;

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  SISLCurve* sisl_helix;
  int stat;
  s1012(start->begin(), center->begin(), axis->begin(),
        frequency, quadrants, 0, &sisl_helix, &stat);
  result->data.reset(Go::SISLCurve2Go(sisl_helix));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_interpolate_curve__doc__, "Construct a cubic spline that interpolates a point cloud\n"
                                                "@param points: The point cloud\n"
                                                "@type points: List of (Point, list of floats or tuple of floats)\n"
                                                "@param parvals: The parameter values for the points\n"
                                                "@type parvals: List of floats\n"
                                                "@return: Cubic spline interpolating the point cloud\n"
                                                "@rtype: Curve");
PyObject* Generate_InterpolateCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"points", "parvals", NULL };
  PyObject* pointso;
  PyObject* parvals;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&pointso,&parvals))
    return NULL;

  if (!PyObject_TypeCheck(pointso,&PyList_Type) || 
      !PyObject_TypeCheck(parvals,&PyList_Type))
    return NULL;

  // get points
  std::vector<double> points;
  for (int i=0; i < PyList_Size(pointso); ++i) {
    PyObject* entryo = PyList_GetItem(pointso,i);
    shared_ptr<Go::Point> entry = PyObject_AsGoPoint(entryo);
    if (entry) {
      points.push_back((*entry)[0]);
      points.push_back((*entry)[1]);
      if (modState.dim == 3)
        points.push_back((*entry)[2]);
    }
  }

  // get parameters
  std::vector<double> params;
  for (int i=0; i < PyList_Size(parvals); ++i)
    params.push_back(PyFloat_AsDouble(PyList_GetItem(parvals,i)));

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Go::SplineInterpolator interp;
  interp.setNaturalConditions();
  std::vector<double> coefs;
  interp.interpolate(params.size(),modState.dim,&params[0],&points[0],coefs);

  result->data.reset(new Go::SplineCurve(interp.basis(),coefs.begin(),
                                         modState.dim,false));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_intersect_cylinder__doc__,"Generate the intersection curve between a surface and the infinite extension of a cylinder\n"
                                                "Set the gap tolerance for more accurate results\n"
                                                "@param surface: a surface to intersect with\n"
                                                "@type surface: Surface\n"
                                                "@param cylinder: the cylinder to intersect with \n"
                                                "@type cylinder: Cylinder\n"
                                                "@return: The intersection curves and/or points\n"
                                                "@rtype: A tuple with a List of Curves and a List of Points");
PyObject* Generate_IntersectCylinder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "cylinder", NULL };
  PyObject* surfo;
  PyObject* cylindero;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,
                                   &surfo,&cylindero))
    return NULL;

  // get ParamSurfaces from Python objects
  shared_ptr<Go::ParamSurface> cyl_ps    = PyObject_AsGoSurface(cylindero);
  shared_ptr<Go::ParamSurface> param_srf = PyObject_AsGoSurface(surfo);
  if (!cyl_ps || !param_srf)
    return NULL;
 
  // check for Cylinder and SplineSurface from ParamSurfaces
  if (cyl_ps->instanceType() != Go::Class_Cylinder) {
    std::cerr << "argument not a cylinder\n";
    return NULL;
  }
  shared_ptr<Go::SplineSurface> srf = convertSplineSurface(param_srf);
  shared_ptr<Go::Cylinder>      cyl = static_pointer_cast<Go::Cylinder>(cyl_ps);
  if (!srf || !cyl)
    return NULL;

  // get SISL surface from SplineSurface
  SISLSurf* sisl_srf = GoSurf2SISL(*srf, true);
  if (!sisl_srf)
    return NULL;

  // SISL input arguments
  Go::Point x_ax, y_ax, z_ax, pt;
  pt   = cyl->getLocation();
  cyl->getCoordinateAxes(x_ax, y_ax, z_ax);
  double radius = cyl->getRadius();
  double center[] = {pt[0],   pt[1],   pt[2]  };
  double axis[]   = {z_ax[0], z_ax[1], z_ax[2]};
  int dim = cyl->dimension();

  // SISL output arguments
  int num_of_curves, num_of_pts; // note that intersections may be at single-points
  double      *pts;
  SISLIntcurve **curves;
  int stat;
  
  // SISL surface intersect cylinder (generating points and curves)
  s1853(sisl_srf, center, axis, radius, dim, 0, modState.gapTolerance, // input arguments
        &num_of_pts, &pts, &num_of_curves, &curves, &stat);            // output arguments
  if (stat < 0) {
    std::cerr << __FUNCTION__ << " ERROR: " << stat << std::endl;
    return NULL;
  } else if (stat > 0)
    std::cerr << __FUNCTION__ << " WARNING: " << stat << std::endl;

  PyObject* result    = PyTuple_New(2);
  PyObject* curveList = PyList_New(0);
  PyObject* pointList = PyList_New(0);

  // store all resulting points
  for (int i=0; i<num_of_pts; i++) {
    Go::Point *intersectPt = new Go::Point();;
    srf->point(*intersectPt, pts[2*i], pts[2*i+1]);
    Point* pyPoint = (Point*) Point_Type.tp_alloc(&Point_Type, 0);
    pyPoint->data = shared_ptr<Go::Point>(intersectPt);
    PyList_Append(pointList, (PyObject*) pyPoint);
  }

  // generate the cylinder intersection curves by SISL s1316
  for (int i=0; i<num_of_curves; i++) {
    double maxStep = 0; // "If maxstep <= gapTolerance, maxstep is neglected. mastep=0 is recommended"
    int makeCrv = 1;    // "1: Make only a geometric curve"
    int graphic = 0;    // "0: Don't draw the curve"
    s1316(sisl_srf, center, axis, radius, dim, 0, modState.gapTolerance, maxStep, // input arguments
          curves[i], makeCrv, graphic, &stat);
  
    if (stat == 3) {
      // ERROR
      std::cerr << "Iteration stopped due to singular point or\n"
                << "degenerate surface. A part of an intersec-\n"
                << "tion curve may have been traced out. If no\n"
                << "curve is traced out, the curve pointers in\n"
                << "the SISLIntcurve object point to NULL.\n";
      return NULL;
    } else if (stat < 0) {
      // ERROR
      std::cerr << __FUNCTION__ << " ERROR " << stat << std::endl;
      return NULL;
    }

    Go::SplineCurve *crv = Go::SISLCurve2Go(curves[i]->pgeom);

    Curve* pyCurve = (Curve*) Curve_Type.tp_alloc(&Curve_Type, 0);
    pyCurve->data = shared_ptr<Go::SplineCurve>(crv);
    PyList_Append(curveList, (PyObject*) pyCurve);
  }

  PyTuple_SetItem(result, 0, (PyObject*) curveList);
  PyTuple_SetItem(result, 1, (PyObject*) pointList);

  return result;
}

PyDoc_STRVAR(generate_line__doc__, "Generate an infinite line\n"
                                   "@param p0: The starting point of the line\n"
                                   "@type p0: Point, list of floats or tuple of floats\n"
                                   "@param direction: The direction of the line\n"
                                   "@type direction: Point, list of floats or tuple of floats\n"
                                   "@return: The line\n"
                                   "@rtype: Curve");
PyObject* Generate_Line(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"p0", "dir", NULL };
  PyObject* p0o;
  PyObject* diro;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&p0o,&diro))
    return NULL;

  shared_ptr<Go::Point> p0 = PyObject_AsGoPoint(p0o);
  shared_ptr<Go::Point> dir = PyObject_AsGoPoint(diro);
  if (!p0 || !dir)
    return NULL;

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::Line(*p0,*dir));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_line_segment__doc__, "Generate a line segment\n"
                                           "@param p0: The starting point of the line\n"
                                           "@type p0: Point, list of floats or tuple of floats\n"
                                           "@param p1: The ending point of the line\n"
                                           "@type p1: Point, list of floats or tuple of floats\n"
                                           "@param relative: (optional) p1 is relative to p0\n"
                                           "@type relative: bool\n"
                                           "@return: The line segment\n"
                                           "@rtype: Curve");
PyObject* Generate_LineSegment(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"p0", "p1", "relative", NULL };
  PyObject* p0o;
  PyObject* p1o;
  bool vector=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO|b",
                                   (char**)keyWords,&p0o,&p1o,&vector))
    return NULL;

  shared_ptr<Go::Point> p0 = PyObject_AsGoPoint(p0o);
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(p1o);
  if (!p0 || !p1)
    return NULL;

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  if (vector)
    result->data.reset(new Go::Line(*p0,*p1));
  else
    result->data.reset(new Go::Line(*p0,*p1-*p0));
  static_pointer_cast<Go::Line>(result->data)->setParamBounds(0.0,1.0);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_crvnonrational__doc__,"Generate a non-rational representation (approximation) of a rational spline curve\n"
                                            "@param original: The initial (rational) curve\n"
                                            "@type original: Curve\n"
                                            "@return: Non-rational B-spline representation\n"
                                            "@rtype: Curve");
PyObject* Generate_CrvNonRational(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"original", NULL };
  PyObject* originalo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&originalo))
    return NULL;

  shared_ptr<Go::ParamCurve> crv = PyObject_AsGoCurve(originalo);
  if (!crv)
    return NULL;

  shared_ptr<Go::SplineCurve> crv_base = convertSplineCurve(crv);

  if (!crv_base)
    return NULL;

  // if it's already B-spline, just return itself
  if(!crv_base->rational()) {
    Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
    result->data = crv_base;

    return (PyObject*)result;
  }

  // extract basis functions and evaluation points
  int dimension = crv_base->dimension();
  Go::BsplineBasis basis = crv_base->basis();
  std::vector<double> greville(basis.numCoefs());
  for(int i=0; i<basis.numCoefs(); i++)
    greville[i] = basis.grevilleParameter(i);

  // evaluate original spline at interpolation points
  std::vector<double> interpolationPoints;
  std::vector<double> weights(0);
  crv_base->gridEvaluator(interpolationPoints, greville);

  Go::SplineCurve *nonrational_crv = 
        Go::CurveInterpolator::regularInterpolation(basis,
                                                    greville,
                                                    interpolationPoints,
                                                    dimension,
                                                    false,
                                                    weights);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data = shared_ptr<Go::SplineCurve>(nonrational_crv);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_resample_curve__doc__,"Generate a surface by sampling another\n"
                                            "@param surface: The surface to sample\n"
                                            "@type surface: Surface\n"
                                            "@param n: The number of points to sample in the first direction\n"
                                            "@type n: integer\n"
                                            "@param p: The order of the curve\n"
                                            "@type p: integer\n"
                                            "@return: The resampled curve\n"
                                            "@rtype: Curve");
PyObject* Generate_ResampleCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curve", "n", "p", NULL };
  PyObject* curveo;
  int n, p;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Oii",
                                   (char**)keyWords,&curveo,&n,&p))
    return NULL;

  shared_ptr<Go::SplineCurve> curve = convertSplineCurve(PyObject_AsGoCurve(curveo));

  if (!curve || n < 1 || p < 1)
    return NULL;

  // create basis functions and evaluation points
  int dimension = curve->dimension();
  // repeated knot at start

  int k=1;
  std::vector<double> knot;
  for (int i=0;i<p;++i)
    knot.push_back(0);
  for (int i=0;i<n-p;++i)
    knot.push_back(k++);
  for (int i=0;i<p;++i)
    knot.push_back(k);

  Go::BsplineBasis basis(n,p,knot.begin());

  basis.rescale(curve->basis().startparam(), curve->basis().endparam());

  std::vector<double> greville(basis.numCoefs());
  for(int i=0; i<basis.numCoefs(); i++)
    greville[i] = basis.grevilleParameter(i);

  // evaluate original spline at interpolation points
  std::vector<double> interpolationPoints;
  std::vector<double> weights(0);
  curve->gridEvaluator(interpolationPoints, greville);

  Go::SplineCurve *res = 
        Go::CurveInterpolator::regularInterpolation(basis,
                                                    greville,
                                                    interpolationPoints,
                                                    dimension,
                                                    false,
                                                    weights);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(res);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_spline_curve__doc__, "Generate a spline curve\n"
                                           "@param params: The parameter values\n"
                                           "@type params: List of float\n"
                                           "@param knots: The knot values\n"
                                           "@type knots: List of float\n"
                                           "@param coeffs: The coefficients\n"
                                           "@type coeffs: List of float\n"
                                           "@param order: The order of the curve (1+p)\n"
                                           "@type order: integer\n"
                                           "@param rational: Is the curve rational?\n"
                                           "@type rational: Boolean\n"
                                           "@return: The spline curve\n"
                                           "@rtype: Curve");
PyObject* Generate_SplineCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"params", "knots", "coeffs", "order", "rational", NULL };
  PyObject* paramso = 0;
  PyObject* knotso  = 0;
  PyObject* coeffso = 0;
  bool rational=false;
  int order=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOi|b",
                                   (char**)keyWords,&paramso, &knotso,
                                                    &coeffso,&order,&rational))
    return NULL;

  if (!PyObject_TypeCheck(paramso,&PyList_Type) || 
      !PyObject_TypeCheck(knotso,&PyList_Type) ||
      !PyObject_TypeCheck(coeffso,&PyList_Type))
    return NULL;

  // get parameters
  std::vector<double> params;
  for (int i=0; i < PyList_Size(paramso); ++i)
    params.push_back(PyFloat_AsDouble(PyList_GetItem(paramso,i)));

  // get knots
  std::vector<double> knots;
  for (int i=0; i < PyList_Size(knotso); ++i)
    knots.push_back(PyFloat_AsDouble(PyList_GetItem(knotso,i)));

  // get coeffs
  std::vector<double> coeffs;
  for (int i=0; i < PyList_Size(knotso); ++i)
    coeffs.push_back(PyFloat_AsDouble(PyList_GetItem(coeffso,i)));

  if (params.size() != knots.size() ||
       knots.size() != coeffs.size()*(rational?1+modState.dim:
                                                 modState.dim)) {
    std::cerr << "Error constructing spline curve" << std::endl;
    return NULL;
  }

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::SplineCurve(params.size(), order,
                                         knots.begin(), coeffs.begin(),
                                         modState.dim, rational));

  return (PyObject*)result;
}


  PyMethodDef CurveFactory_methods[] = {
     {(char*)"ApproximateCurve",      (PyCFunction)Generate_ApproximateCurve,  METH_VARARGS|METH_KEYWORDS, generate_approximate_curve__doc__},
     {(char*)"Circle",                (PyCFunction)Generate_Circle,            METH_VARARGS|METH_KEYWORDS, generate_circle__doc__},
     {(char*)"CircleSegment",         (PyCFunction)Generate_CircleSegment,     METH_VARARGS|METH_KEYWORDS, generate_circle_segment__doc__},
     {(char*)"Ellipse",               (PyCFunction)Generate_Ellipse,           METH_VARARGS|METH_KEYWORDS, generate_ellipse__doc__},
     {(char*)"EllipticSegment",       (PyCFunction)Generate_EllipticSegment,   METH_VARARGS|METH_KEYWORDS, generate_elliptic_segment__doc__},
     {(char*)"Helix",                 (PyCFunction)Generate_Helix,             METH_VARARGS|METH_KEYWORDS, generate_helix__doc__},
     {(char*)"InterpolateCurve",      (PyCFunction)Generate_InterpolateCurve,  METH_VARARGS|METH_KEYWORDS, generate_interpolate_curve__doc__},
     {(char*)"IntersectCylinder",     (PyCFunction)Generate_IntersectCylinder, METH_VARARGS|METH_KEYWORDS, generate_intersect_cylinder__doc__},
     {(char*)"Line",                  (PyCFunction)Generate_Line,              METH_VARARGS|METH_KEYWORDS, generate_line__doc__},
     {(char*)"LineSegment",           (PyCFunction)Generate_LineSegment,       METH_VARARGS|METH_KEYWORDS, generate_line_segment__doc__},
     {(char*)"NonRationalCurve",      (PyCFunction)Generate_CrvNonRational,    METH_VARARGS|METH_KEYWORDS, generate_crvnonrational__doc__},
     {(char*)"ResampleCurve",         (PyCFunction)Generate_ResampleCurve,     METH_VARARGS|METH_KEYWORDS, generate_resample_curve__doc__},
     {(char*)"SplineCurve",           (PyCFunction)Generate_SplineCurve,       METH_VARARGS|METH_KEYWORDS, generate_spline_curve__doc__},
     {NULL,                           NULL,                                    0,                          NULL}
  };

PyDoc_STRVAR(curve_factory__doc__,"A module with methods for generating curves");

PyMODINIT_FUNC
init_CurveFactory_Module()
{
  CurveFactory_module = Py_InitModule3((char*)"GoTools.CurveFactory",CurveFactory_methods,curve_factory__doc__);
  modState.addInfo(CurveFactory_module);
}
}

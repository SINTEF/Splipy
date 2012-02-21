#include "curvefactory.h"

#include "curve.h"
#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

extern "C" {
PyObject* CurveFactory_module;

PyDoc_STRVAR(generate_circle__doc__, "Generate a circle\n"
                                     "@param center: The center of the circle\n"
                                     "@type center: Point, list of floats or tuple of floats\n"
                                     "@param radius: The radius of the circle\n"
                                     "@type radius: float\n"
                                     "@param normal: The normal of the circle\n"
                                     "@type normal: Point, list of floats or tuple of floats\n"
                                     "@return: Curve describing the circle");
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
                                             "@param radius: The radius of the circle\n"
                                             "@type radius: float\n"
                                             "@param angle: The angle of the circle segment in radians\n"
                                             "@type angle: float\n"
                                             "@param normal: The normal of the circle\n"
                                             "@type normal: Point, list of floats or tuple of floats\n"
                                             "@return: Curve describing the circle segment");
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
                                      "@return: Curve describing the ellipse");
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
                                               "@return: Curve describing the ellipse");
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
                                    "@return: Curve describing the helix");
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

PyDoc_STRVAR(generate_interpolate_curve__doc__, "Construct a spline curve from a point cloud\n"
                                                "@param points: The point cloud\n"
                                                "@type points: List of (Point, list of floats or tuple of floats)\n"
                                                "@param parvals: The parameter values for the points\n"
                                                "@type parvals: List of floats\n"
                                                "@param order: (optional) The order of the spline curve\n"
                                                "@type order: int >= 1\n"
                                                "@param maxiter: (optional) The maxium number of iterations in interpolation scheme\n"
                                                "@type maxiter: int >= 1\n"
                                                "@return: Spline curve interpolating the point cloud\n");
PyObject* Generate_InterpolateCurve(PyObject* self, PyObject* args, PyObject* kwds)
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

PyDoc_STRVAR(generate_line__doc__, "Generate an infinite line\n"
                                   "@param p0: The starting point of the line\n"
                                   "@type p0: Point, list of floats or tuple of floats\n"
                                   "@param direction: The direction of the line\n"
                                   "@type direction: Point, list of floats or tuple of floats\n"
                                   "@return: Curve describing the line");
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
                                           "@return: Curve describing the line segment");
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

  PyMethodDef CurveFactory_methods[] = {
     {(char*)"Circle",                (PyCFunction)Generate_Circle,           METH_VARARGS|METH_KEYWORDS, generate_circle__doc__},
     {(char*)"CircleSegment",         (PyCFunction)Generate_CircleSegment,    METH_VARARGS|METH_KEYWORDS, generate_circle_segment__doc__},
     {(char*)"Ellipse",               (PyCFunction)Generate_Ellipse,          METH_VARARGS|METH_KEYWORDS, generate_ellipse__doc__},
     {(char*)"EllipticSegment",       (PyCFunction)Generate_EllipticSegment,  METH_VARARGS|METH_KEYWORDS, generate_elliptic_segment__doc__},
     {(char*)"Helix",                 (PyCFunction)Generate_Helix,            METH_VARARGS|METH_KEYWORDS, generate_helix__doc__},
     {(char*)"InterpolateCurve",      (PyCFunction)Generate_InterpolateCurve, METH_VARARGS|METH_KEYWORDS, generate_interpolate_curve__doc__},
     {(char*)"Line",                  (PyCFunction)Generate_Line,             METH_VARARGS|METH_KEYWORDS, generate_line__doc__},
     {(char*)"LineSegment",           (PyCFunction)Generate_LineSegment,      METH_VARARGS|METH_KEYWORDS, generate_line_segment__doc__},
     {NULL,                           NULL,                                   0,                          NULL}
  };

PyDoc_STRVAR(curve_factory__doc__,"A module with methods for generating curves");

PyMODINIT_FUNC
init_CurveFactory_Module()
{
  CurveFactory_module = Py_InitModule3((char*)"GoTools.CurveFactory",CurveFactory_methods,curve_factory__doc__);
  modState.addInfo(CurveFactory_module);
}
}

#include "curvefactory.h"

#include "curve.h"
#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

#ifdef HAS_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL GEOMOD_ARRAY_API
#define NO_IMPORT_ARRAY
#include <arrayobject.h>
#endif

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

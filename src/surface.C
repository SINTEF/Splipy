//===========================================================================
//                                                                           
// File: surface.C
//                                                                           
// Created: Mon Feb 20 13:22:00 2012                                         
//                                                                           
// Author: Arne Morten Kvarving <arne.morten.kvarving@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description: GoTools.Surface python class implementation
//                                                                           
//===========================================================================

#include "surface.h"
#include "curve.h"
#include "pyutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineSurface.h"

#include <fstream>
#include <sstream>

extern "C"
{
PyTypeObject Surface_Type;

PyObject* Surface_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Surface* self;
  self = (Surface*)type->tp_alloc(type,0);

  // only generators allocate the paramsurface!
  return (PyObject*)self;
}

void Surface_Dealloc(Surface* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyDoc_STRVAR(surface_clone__doc__,"Clone a surface\n"
                                  "@return: New copy of surface\n");
PyObject* Surface_Clone(PyObject* self, PyObject* args, PyObject* kwds)
{
  Surface* res = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  shared_ptr<Go::ParamSurface> srf = PyObject_AsGoSurface(self);
  res->data.reset(srf->clone());
 
  return (PyObject*)res;
}

PyObject* Surface_Str(Surface* self)
{
  std::stringstream str;
  printSurfaceToStream(str,self->data);
  return PyString_FromString(str.str().c_str());
}

PyDoc_STRVAR(surface_flip_parametrization__doc__,"Flip surface parametrization\n"
                                                 "@param direction: The parametric direction to flip (0=u, 1=v)\n"
                                                 "@type direction: int\n"
                                                 "@return: None");
PyObject* Surface_FlipParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", NULL };
  PyObject* axiso;
  int direction = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&direction))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;

   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }

   Go::SplineSurface *ssurf = surface->asSplineSurface();
   ssurf->reverseParameterDirection(direction == 0);

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(surface_raise_order__doc__,"Raise order of a spline surface\n"
                                        "@param raise_u: Raise of order in u\n"
                                        "@type raise_u: int (>= 0)\n"
                                        "@param raise_v: Raise of order in v\n"
                                        "@type raise_v: int (>= 0)\n"
                                        "@returns: None");
PyObject* Surface_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"raise_u", "raise_v", NULL };
  int raise_u=0, raise_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords,&raise_u,&raise_v))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }
   static_pointer_cast<Go::SplineSurface>(surface)->raiseOrder(raise_u,raise_v);

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(surface_get_edges__doc__,"Return the four edge curves in (parametric) order: bottom, right, top, left\n"
                                      "@return: A list of the four edge curves");
PyObject* Surface_GetEdges(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
  if (!surface->isSpline()) {
    Surface* surf = (Surface*)self;
    surf->data = convertSplineSurface(surface);
    surface = surf->data;
  }

  Go::SplineSurface* ss = surface->asSplineSurface();
  Curve* curve0 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Curve* curve1 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Curve* curve2 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Curve* curve3 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  curve0->data = shared_ptr<Go::SplineCurve>(ss->edgeCurve(0));
  curve1->data = shared_ptr<Go::SplineCurve>(ss->edgeCurve(1));
  curve2->data = shared_ptr<Go::SplineCurve>(ss->edgeCurve(2));
  curve3->data = shared_ptr<Go::SplineCurve>(ss->edgeCurve(3));

  PyObject* result = PyList_New(0);
  PyList_Append(result, (PyObject*) curve0);
  PyList_Append(result, (PyObject*) curve1);
  PyList_Append(result, (PyObject*) curve2);
  PyList_Append(result, (PyObject*) curve3);

  return result;
}


PyDoc_STRVAR(surface_get_knots__doc__,"Return unique knots for a spline surface\n"
                                      "@return: Tuple with List of float");
PyObject* Surface_GetKnots(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
  if (!surface->isSpline()) {
    Surface* surf = (Surface*)self;
    surf->data = convertSplineSurface(surface);
    surface = surf->data;
  }

  PyObject* result = PyTuple_New(2);
  for (int i=0;i<2;++i) {
    PyObject* list = PyList_New(0);
    std::vector<double> knots;
    static_pointer_cast<Go::SplineSurface>(surface)->basis(i).knotsSimple(knots);
    for (std::vector<double>::iterator it  = knots.begin();
                                       it != knots.end();++it) {
      PyList_Append(list,Py_BuildValue((char*)"d",*it));
    }
    PyTuple_SetItem(result,i,list);
  }

  return result;
}

PyDoc_STRVAR(surface_insert_knot__doc__,"Insert a knot in a spline surface\n"
                                        "@param direction: Direction to insert knot in\n"
                                        "@type direction: int (0 or 1)\n"
                                        "@param knot: The knot to insert\n"
                                        "@type knot: float\n"
                                        "@return: None");
PyObject* Surface_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", "knot", NULL };
  int direction=0;
  double knot;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"id",
                                   (char**)keyWords,&direction,&knot))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }
   if (direction == 0)
     static_pointer_cast<Go::SplineSurface>(surface)->insertKnot_u(knot);
   else
     static_pointer_cast<Go::SplineSurface>(surface)->insertKnot_v(knot);

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(surface_project__doc__,"Project the surface onto an axis or plane parallel to the cartesian coordinate system\n"
                                    "@param axis: The axis or plane to project onto (\"X\",\"Y\",\"Z\" or a combination of these)\n"
                                    "@type axis: string\n"
                                    "@return: None");
PyObject* Surface_Project(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", NULL };
  shared_ptr<Go::ParamSurface> surf = PyObject_AsGoSurface(self);
  char *sAxis;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&sAxis) || !surf)
    return NULL;

  bool bAxis[3];
  bAxis[0] = true;
  bAxis[1] = true;
  bAxis[2] = true;
  while(*sAxis != 0) {
    if(*sAxis == 'x' || *sAxis == 'X')
      bAxis[0] = false;
    else if(*sAxis == 'y' || *sAxis == 'Y')
      bAxis[1] = false;
    else if(*sAxis == 'z' || *sAxis == 'Z')
      bAxis[2] = false;
    sAxis++;
  }

  Surface* srf = (Surface*)self;
  srf->data = convertSplineSurface(surf);

  shared_ptr<Go::SplineSurface> ssurf    = static_pointer_cast<Go::SplineSurface>(srf->data);
  bool                          rational = ssurf->rational();
  int                           dim      = ssurf->dimension();
  std::vector<double>::iterator coefs    = (rational) ? ssurf->rcoefs_begin() : ssurf->coefs_begin();
  std::vector<double>::iterator coefsEnd = (rational) ? ssurf->rcoefs_end()   : ssurf->coefs_end();
  while(coefs != coefsEnd) {
    if(bAxis[0] && dim>0)
      coefs[0] = 0.0;
    if(bAxis[1] && dim>1)
      coefs[1] = 0.0;
    if(bAxis[2] && dim>2)
      coefs[2] = 0.0;
    coefs += (dim+rational);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(surface_rotate__doc__,"Rotate a surface around an axis\n"
                                   "@param axis: The axis to rotate around\n"
                                   "@type axis: Point, list of floats or tuple of floats\n"
                                   "@param angle: Angle to rotate surface with in radians\n"
                                   "@type angle: float\n"
                                   "@return: None");
PyObject* Surface_Rotate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", "angle", NULL };
  PyObject* axiso;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  if (!surface || !axis)
    return NULL;

   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }

   Go::GeometryTools::rotateSplineSurf(*axis, angle,
                        *static_pointer_cast<Go::SplineSurface>(surface));

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(surface_swap_parametrization__doc__,"Swaps the two surface parameter directions\n"
                                                 "@return: None");
PyObject* Surface_SwapParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;

   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }

   Go::SplineSurface *ssurf = surface->asSplineSurface();
   ssurf->swapParameterDirection();

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(surface_translate__doc__,"Translate a surface along a given vector\n"
                                      "@param vector: The vector to translate along\n"
                                      "@type vector: Point, list of floats or tuple of floats\n"
                                      "@return: None");
PyObject* Surface_Translate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"vector", NULL };
  PyObject* veco;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&veco))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  shared_ptr<Go::Point>        vec     = PyObject_AsGoPoint(veco);
  if (!surface || !vec)
    return NULL;

   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }

   Go::GeometryTools::translateSplineSurf(*vec, 
                              *static_pointer_cast<Go::SplineSurface>(surface));

   Py_INCREF(Py_None);
   return Py_None;
}

PyObject* Surface_Add(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::ParamSurface> surf = PyObject_AsGoSurface(o1);
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(o2);
  Surface* result = NULL;
  if (surf && point) {
    result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    if (surf->isSpline())
      result->data = shared_ptr<Go::ParamSurface>(
          new Go::SplineSurface(*static_pointer_cast<Go::SplineSurface>(surf)));
    else
      result->data = convertSplineSurface(surf);

    Go::GeometryTools::translateSplineSurf(*point, 
                         *static_pointer_cast<Go::SplineSurface>(result->data));
  }

  return (PyObject*) result;
}

PyMethodDef Surface_methods[] = {
     {(char*)"Clone",               (PyCFunction)Surface_Clone,                 METH_VARARGS|METH_KEYWORDS, surface_clone__doc__},
     {(char*)"FlipParametrization", (PyCFunction)Surface_FlipParametrization,   METH_VARARGS|METH_KEYWORDS, surface_flip_parametrization__doc__},
     {(char*)"GetEdges",            (PyCFunction)Surface_GetEdges,              METH_VARARGS|METH_KEYWORDS, surface_get_edges__doc__},
     {(char*)"GetKnots",            (PyCFunction)Surface_GetKnots,              METH_VARARGS|METH_KEYWORDS, surface_get_knots__doc__},
     {(char*)"InsertKnot",          (PyCFunction)Surface_InsertKnot,            METH_VARARGS|METH_KEYWORDS, surface_insert_knot__doc__},
     {(char*)"Project",             (PyCFunction)Surface_Project,               METH_VARARGS|METH_KEYWORDS, surface_project__doc__},
     {(char*)"RaiseOrder",          (PyCFunction)Surface_RaiseOrder,            METH_VARARGS|METH_KEYWORDS, surface_raise_order__doc__},
     {(char*)"Rotate",              (PyCFunction)Surface_Rotate,                METH_VARARGS|METH_KEYWORDS, surface_rotate__doc__},
     {(char*)"SwapParametrization", (PyCFunction)Surface_SwapParametrization,   METH_VARARGS|METH_KEYWORDS, surface_swap_parametrization__doc__},
     {(char*)"Translate",           (PyCFunction)Surface_Translate,             METH_VARARGS|METH_KEYWORDS, surface_translate__doc__},
     {NULL,           NULL,                     0,            NULL}
   };

PyNumberMethods Surface_operators = {0};

PyDoc_STRVAR(surface__doc__, "A parametric description of a surface");
void init_Surface_Type()
{
  InitializeTypeObject(&Surface_Type);
  Surface_operators.nb_add = Surface_Add;
  Surface_Type.tp_name = "GoTools.Surface";
  Surface_Type.tp_basicsize = sizeof(Surface);
  Surface_Type.tp_dealloc = (destructor)Surface_Dealloc;
  Surface_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  Surface_Type.tp_doc = surface__doc__;
  Surface_Type.tp_methods = Surface_methods;
  Surface_Type.tp_base = 0;
  Surface_Type.tp_new = Surface_New;
  Surface_Type.tp_str = (reprfunc)Surface_Str;
  Surface_Type.tp_as_number = &Surface_operators;
  PyType_Ready(&Surface_Type);
}

}

//! \brief Check if a surface is a spline surface
//! \return The surface itself if it s a spline surface,
//          the surface converted to a spline surface otherwise.
shared_ptr<Go::SplineSurface> convertSplineSurface(shared_ptr<Go::ParamSurface> surface)
{
  if (surface->instanceType() == Go::Class_SplineSurface)
    return dynamic_pointer_cast<Go::SplineSurface, Go::ParamSurface>(surface);
  shared_ptr<Go::ElementarySurface> e_surface = dynamic_pointer_cast<Go::ElementarySurface, Go::ParamSurface>(surface);
  if (e_surface)
    return shared_ptr<Go::SplineSurface>(e_surface->geometrySurface());
  return shared_ptr<Go::SplineSurface>();
}

void printSurfaceToStream(std::ostream& str, shared_ptr<Go::ParamSurface> surf)
{
  if (surf) {
    if (surf->instanceType() == Go::Class_Plane)
      str << "Plane:" << std::endl;
    if (surf->instanceType() == Go::Class_Sphere)
      str << "Sphere surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Cylinder)
      str << "Cylinder surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Cone)
      str << "Cone surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Torus)
      str << "Torus surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Disc)
      str << "Circular disc:" << std::endl;
    str << *surf;
  } else
    str << "(empty)";
}

void WriteSurfaceG2(std::ofstream& g2_file, Surface* surface, bool convert)
{
  if (!surface->data)
    return;
  if (convert) {
    shared_ptr<Go::SplineSurface> srf = convertSplineSurface(surface->data);
    if (srf) {
      srf->writeStandardHeader(g2_file);
      srf->write(g2_file);
    }
  } else {
    surface->data->writeStandardHeader(g2_file);
    surface->data->write(g2_file);
  }
}

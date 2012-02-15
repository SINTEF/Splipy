#include "surfacefactory.h"

#include "curve.h"
#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"
#include "surface.h"

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Disc.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/Torus.h"

PyObject* Generate_Plane(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"p0", "normal", NULL };
  PyObject* p0o;
  PyObject* normalo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&p0o,&normalo))
    return NULL;

  shared_ptr<Go::Point> p0     = PyObject_AsGoPoint(p0o);
  shared_ptr<Go::Point> normal = PyObject_AsGoPoint(normalo);
  if (!p0o || !normal)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(new Go::Plane(*p0,*normal));

  return (PyObject*)result;
}

PyObject* Generate_SphereSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "radius", NULL };
  PyObject* centero;
  double radius;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&centero,&radius))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  if (!center)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Go::Point x_axis(1.0, 0.0, 0.0);
  Go::Point z_axis(0.0, 0.0, 1.0);

  result->data.reset(new Go::Sphere(radius,*center,x_axis,z_axis));

  return (PyObject*)result;
}

PyObject* Generate_CylinderSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "axis", "radius", "height", NULL };
  PyObject* centero;
  PyObject* axiso;
  double radius;
  double height = std::numeric_limits<double>::infinity();
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOd|d",
                                   (char**)keyWords,&centero,
                                   &axiso,&radius,&height))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis   = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(new Go::Cylinder(radius,*center,*axis,someNormal(*axis)));
  if (height != std::numeric_limits<double>::infinity())
    static_pointer_cast<Go::Cylinder>(result->data)->setParameterBounds(0.0,0.0,2.0*M_PI,height);

  return (PyObject*)result;
}

PyObject* Generate_ConeSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"apex", "axis", "angle", "radius", "height", NULL };
  PyObject* centero;
  PyObject* axiso;
  double apex;
  double angle;
  double height = std::numeric_limits<double>::infinity();
  double radius = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOd|dd",
                                   (char**)keyWords,&centero,
                                   &axiso,&apex,&angle,&height,&radius))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis   = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);

  result->data.reset(new Go::Cone(radius,*center,*axis,someNormal(*axis), angle));
  static_pointer_cast<Go::Cylinder>(result->data)->setParameterBounds(0.0,0.0,2.0*M_PI,height);

  return (PyObject*)result;
}

PyObject* Generate_TorusSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "axis", "major_radius", "minor_radius", NULL };
  PyObject* centero;
  PyObject* axiso;
  double major_radius;
  double minor_radius;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOdd",
                                   (char**)keyWords,&centero,&axiso,
                                   &major_radius,&minor_radius))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis   = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);

  result->data.reset(new Go::Torus(major_radius,minor_radius,
                                      *center,*axis,someNormal(*axis)));

  return (PyObject*)result;
}

PyObject* Generate_CircularDisc(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "boundarypoint", "normal", NULL };
  PyObject* centero;
  PyObject* bpointo;
  PyObject* normalo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO|O",
                                   (char**)keyWords,&centero,&bpointo,&normalo))
    return NULL;

  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> bpoint = PyObject_AsGoPoint(bpointo);
  if (!center || !bpoint)
    return NULL;

  Go::Point normal(0.0,0.0);
  if (modState.dim == 3) {
    shared_ptr<Go::Point> norm = PyObject_AsGoPoint(normalo);
    if (!norm)
      return NULL;
    normal = *norm;
  }

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);

  Go::Point x_axis = *bpoint-*center;
  Go::Point z_axis = normal-((normal*x_axis)/x_axis.length2())*x_axis;
  result->data.reset(new Go::Disc(*center,x_axis.length(),x_axis,z_axis));

  return (PyObject*)result;
}

PyObject* Generate_Rectangle(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"corner", "axis_x", "axis_y", "length_x", "length_y", NULL };
  PyObject* cornero;
  PyObject* axisxo;
  PyObject* axisyo;
  double length_x, length_y;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOdd",
                                   (char**)keyWords,&cornero,&axisxo,
                                   &axisyo,&length_x,&length_y))
    return NULL;

  shared_ptr<Go::Point> corner = PyObject_AsGoPoint(cornero);
  shared_ptr<Go::Point> axisx  = PyObject_AsGoPoint(axisxo);
  shared_ptr<Go::Point> axisy  = PyObject_AsGoPoint(axisyo);
  if (!corner || !axisx || !axisy)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);

  Go::Point dir_v_ortog = *axisy-(*axisx*(*axisy))/axisx->length2()*(*axisx);
  if (modState.dim == 3)
    *axisy = (*axisx) % dir_v_ortog;
  result->data.reset(new Go::Plane(*corner,*axisy, *axisx));
  static_pointer_cast<Go::Plane>(result->data)->setParameterBounds(0.0,0.0,length_x,length_y);

  return (PyObject*)result;
}

PyObject* Generate_SweepCurveRotational(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curve", "pos", "axis", "angle", NULL };
  PyObject* curveo;
  PyObject* poso;
  PyObject* axiso;
  double angle = 2*M_PI;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO|d",
                                   (char**)keyWords,&curveo,&poso,&axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(curveo);
  shared_ptr<Go::Point> pos   = PyObject_AsGoPoint(poso);
  shared_ptr<Go::Point> axis  = PyObject_AsGoPoint(axiso);

  if (!curve || !pos || !axis)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  shared_ptr<Go::SplineCurve> c = convertSplineCurve(curve);
  result->data.reset(Go::SweepSurfaceCreator::rotationalSweptSurface(*c, angle, *pos, *axis));

  return (PyObject*)result;
}

PyObject* Generate_TrimSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"original", "trim", NULL };
  PyObject* originalo;
  PyObject* trimo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&originalo,&trimo))
    return NULL;

  shared_ptr<Go::ParamSurface> surf1 = PyObject_AsGoSurface(originalo);
  shared_ptr<Go::ParamSurface> surf2 = PyObject_AsGoSurface(trimo);
  if (!surf1 || !surf2)
    return NULL;

  shared_ptr<Go::SplineSurface> surf_base = convertSplineSurface(surf1);
  shared_ptr<Go::SplineSurface> surf_remove = convertSplineSurface(surf2);

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  std::vector<std::shared_ptr<Go::BoundedSurface> > trimmed_sfs =
      Go::BoundedUtils::trimSurfWithSurf(surf_base, surf_remove, fabs(modState.gapTolerance));
  result->data = trimmed_sfs[0];

  return (PyObject*)result;
}

PyObject* Generate_AddLoop(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"original", "loop", NULL };
  PyObject* originalo;
  PyObject* loopo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&originalo,&loopo))
    return NULL;

  shared_ptr<Go::ParamSurface> par_surf = PyObject_AsGoSurface(originalo);
  if (!par_surf)
    return NULL;

  std::vector<shared_ptr<Go::ParamCurve> > loops;
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(loopo);
  if (curve)
    loops.push_back(curve);
  else if (PyObject_TypeCheck(loopo,&PyList_Type)) {
    for (int i=0;i<PyList_Size(loopo); ++i) {
      PyObject* obj = PyList_GetItem(loopo,i);
      shared_ptr<Go::ParamCurve> curv = PyObject_AsGoCurve(obj);
      if (curv)
        loops.push_back(curv);
    }
  }
  if (loops.size() < 1)
    return NULL;

  shared_ptr<Go::ParamSurface> underlying;   // The underlying surface for all curves
  std::vector<Go::CurveLoop> former_loops;        // The already existing loops (at least the outer boundary)
  std::vector<std::vector<shared_ptr<Go::CurveOnSurface> > > curves;   // The CurveOnSurface-instances used when making the new surfaces
  std::vector<double> tolerances;     // The tolerance used for each loop

  // Get underlying surface and already existing loops
  Go::ClassType surf_type = par_surf->instanceType();
  if (surf_type == Go::Class_BoundedSurface) {
    shared_ptr<Go::BoundedSurface> bs = dynamic_pointer_cast<Go::BoundedSurface, Go::ParamSurface>(par_surf);
    underlying = bs->underlyingSurface();
    former_loops = bs->allBoundaryLoops();
    for (int i = 0; i < former_loops.size(); ++i)
      tolerances.push_back(former_loops[i].getSpaceEpsilon());
  } else {
    underlying = par_surf;
    shared_ptr<Go::ParamSurface> surf_as_spline;
    if (surf_type == Go::Class_SplineSurface)
      surf_as_spline = par_surf;
    else {
      // If the original surface is an ElementarySurface, we must convert to SplineSurface
      // to get boundary curves, as outerBoundaryLoop() often returns nothing otherwise
      shared_ptr<Go::ElementarySurface> elementary = dynamic_pointer_cast<Go::ElementarySurface, Go::ParamSurface>(par_surf);
      surf_as_spline = shared_ptr<Go::ParamSurface>(elementary->geometrySurface());
    }
    former_loops.push_back(surf_as_spline->outerBoundaryLoop());
    tolerances.push_back(fabs(modState.gapTolerance));
  }

  // Store former loops in 'curves'
  for (int i = 0; i < former_loops.size(); ++i) {
    Go::CurveLoop loop = former_loops[i];
    std::vector<shared_ptr<Go::CurveOnSurface> > cos_loop;
    for (std::vector<shared_ptr<Go::ParamCurve> >::const_iterator it  = loop.begin();
                                                                  it != loop.end(); ++it) {
      shared_ptr<Go::CurveOnSurface> cos;
      if ((*it)->instanceType() == Go::Class_CurveOnSurface)
        cos = dynamic_pointer_cast<Go::CurveOnSurface, Go::ParamCurve>(*it);
      else
        cos = shared_ptr<Go::CurveOnSurface>(new Go::CurveOnSurface(underlying, *it, false));
      cos_loop.push_back(cos);
    }
    curves.push_back(cos_loop);
  }

  // Insert new loops
  for (size_t i = 0; i < loops.size(); ++i)
  {
    shared_ptr<Go::ParamCurve> par_curv = loops[i];

    // Currently assumes we have one periodic curve, split it into two curves
    double start_par = par_curv->startparam();
    double end_par = par_curv->endparam();
    double middle = 0.5 * (start_par + end_par);

    // Get the two curves
    shared_ptr<Go::ParamCurve> pcrv1(par_curv->subCurve(start_par, middle));
    shared_ptr<Go::ParamCurve> pcrv2(par_curv->subCurve(middle, end_par));

    // Get corresponing CurveOnSurface instances
    shared_ptr<Go::CurveOnSurface> cos1(new Go::CurveOnSurface(underlying, pcrv1, false));
    shared_ptr<Go::CurveOnSurface> cos2(new Go::CurveOnSurface(underlying, pcrv2, false));

    // Get parameter curves, needed to determine loop orientation
//    if (!(cos1->ensureParCrvExistence(gapTol)) || !(cos2->ensureParCrvExistence(gapTol)))
//      raiseError(cmd, "Could not create parameter curves for curve on surface, needed to determine orientation");

    std::vector<shared_ptr<Go::CurveOnSurface> > new_loop;
    new_loop.push_back(cos1);
    new_loop.push_back(cos2);

    // Reverse if curve is counter-clockwise in parameter space
    if (Go::LoopUtils::paramIsCCW(new_loop, fabs(modState.gapTolerance)))
    {
      cos1->reverseParameterDirection();
      cos2->reverseParameterDirection();
    }

    // Store new loop
    curves.push_back(new_loop);
    tolerances.push_back(fabs(modState.gapTolerance));
  }

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(new Go::BoundedSurface(underlying, curves, tolerances));
}

PyObject* Generate_LoftCurves(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curves", NULL };
  PyObject* curveso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&curveso))
    return NULL;

  if (!PyObject_TypeCheck(curveso,&PyList_Type))
    return NULL;

  std::vector<shared_ptr<Go::SplineCurve> > curves;
  for (int i=0; i < PyList_Size(curveso); ++i) {
    PyObject* curvo = PyList_GetItem(curveso,i);
    shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(curvo);
    if (!curve)
      continue;
    curves.push_back(convertSplineCurve(curve));
  }
  if (curves.size() < 2)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(Go::CoonsPatchGen::loftSurface(curves.begin(),curves.size()));

  return (PyObject*)result;
}

PyObject* Generate_SweepCurveLinear(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curve1", "curve2", "point", NULL };
  PyObject* curve1o;
  PyObject* curve2o;
  PyObject* pointo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO",
                                   (char**)keyWords,&curve1o,&curve2o,&pointo))
    return NULL;

  shared_ptr<Go::ParamCurve> curve1 = PyObject_AsGoCurve(curve1o);
  shared_ptr<Go::ParamCurve> curve2 = PyObject_AsGoCurve(curve2o);
  shared_ptr<Go::Point> point       = PyObject_AsGoPoint(pointo);
  if (!curve1 || !curve2 || !point)
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(Go::SweepSurfaceCreator::linearSweptSurface(*convertSplineCurve(curve1),
                                                                    *convertSplineCurve(curve2),*point));

  return (PyObject*)result;
}

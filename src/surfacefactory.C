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
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Disc.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/SurfaceInterpolator.h"


extern "C" {
PyObject* SurfaceFactory_module;

PyDoc_STRVAR(generate_addloop__doc__,"Add a loop to a surface\n"
                                     "@param original: The original surface\n"
                                     "@type original: Surface\n"
                                     "@param loop: The curve describing the loop\n"
                                     "@type loop: Curve\n"
                                     "@return: Surface with loop added");
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

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_circular_disc__doc__,"Generate a circular disc\n"
                                           "@param center: The center of the disc\n"
                                           "@type center: Point, list of floats or tuple of floats\n"
                                           "@param boundarypoint: A point on the boundary of the disc\n"
                                           "@type boundarypoint: Point, list of floats or tuple of floats\n"
                                           "@param normal: The normal of the disc\n"
                                           "@type normal: Point, list of floats or tuple of floats\n"
                                           "@return: Surface describing the disc");
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

PyDoc_STRVAR(generate_cone_surface__doc__,"Generate a cone surface\n"
                                          "@param apex: The apex of the cone\n"
                                          "@type apex: Point, list of floats or tuple of floats\n"
                                          "@param axis: The axis of the cone\n"
                                          "@type axis: Point, list of floats or tuple of floats\n"
                                          "@param angle: The angle of the cone in radians\n"
                                          "@type angle: float\n"
                                          "@param radius: (optional) The radius of the cone\n"
                                          "@type radius: float\n"
                                          "@param height: (optional) The height of the cone\n"
                                          "@type height: float\n"
                                          "@return: Surface describing the cone surface");
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

PyDoc_STRVAR(generate_coons_surface_patch__doc__,"Create a new SplineSurface representing the coons patch defined by a loop of four boundary curves\n"
                                                 "@param curves: There must be exactly four curves, they must all be of type 'SplineCurve', and they must be nonrational\n"
                                                 "@type curves: List of Curve\n"
                                                 "@return: Surface describing the interior of the four curves");
PyObject* Generate_CoonsSurfacePatch(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curves", NULL };
  PyObject* curveso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&curveso))
    return NULL;

  if (!PyObject_TypeCheck(curveso,&PyList_Type))
    return NULL;

  std::vector<shared_ptr<Go::ParamCurve> > curves;
  for (int i=0; i < PyList_Size(curveso); ++i) {
    PyObject* curvo = PyList_GetItem(curveso,i);
    shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(curvo);
    if (!curve)
      continue;
    curves.push_back(convertSplineCurve(curve));
  }
  if (curves.size() != 4)
    return NULL;

  Go::CurveLoop loop(curves, modState.gapTolerance);
  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(Go::CoonsPatchGen::createCoonsPatch(loop));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_cylinder_surface__doc__,"Generate a cylinder surface\n"
                                              "@param center: The center of the cylinder\n"
                                              "@type center: Point, list of floats or tuple of floats\n"
                                              "@param axis: The axis of the cylinder\n"
                                              "@type axis: Point, list of floats or tuple of floats\n"
                                              "@param radius: The radius of the cylinder\n"
                                              "@type radius: float\n"
                                              "@param height: (optional) The height of the cylinder\n"
                                              "@type height: float\n"
                                              "@return: Surface describing the cylinder surface");
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

PyDoc_STRVAR(generate_linear_curve_sweep__doc__,"Generate a surface by linearly sweeping a curve along a curve\n"
                                                "@param curve1: The curve to sweep\n"
                                                "@type curve1: Curve\n"
                                                "@param curve2: The curve to sweep along\n"
                                                "@type curve2: Curve\n"
                                                "@param point: The start point of the new surface\n"
                                                "@type point: Point\n"
                                                "@return: Surface describing the swept area");
PyObject* Generate_LinearCurveSweep(PyObject* self, PyObject* args, PyObject* kwds)
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

  Go::Point pt(*point);
  pt *= -1;
  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(Go::SweepSurfaceCreator::linearSweptSurface(*convertSplineCurve(curve1),
                                                                    *convertSplineCurve(curve2),pt));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_loft_curves__doc__,"Generate a surface by lofting curves\n"
                                         "@param curves: The curves to loft\n"
                                         "@type curves: List of Curve\n"
                                         "@return: Surface describing the lofted surface");
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

PyDoc_STRVAR(generate_nonrational__doc__,"Generate a non-rational representation (approximation) of a rational spline surface\n"
                                         "@param original: The surface to trim\n"
                                         "@type original: Surface\n"
                                         "@return: Non-rational B-spline representation of the surface");
PyObject* Generate_NonRational(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"original", NULL };
  PyObject* originalo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&originalo))
    return NULL;

  shared_ptr<Go::ParamSurface> surf = PyObject_AsGoSurface(originalo);
  if (!surf)
    return NULL;

  shared_ptr<Go::SplineSurface> surf_base = convertSplineSurface(surf);

  if (!surf_base)
    return NULL;

  // if it's already B-spline, just return itself
  if(!surf_base->rational()) {
    Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    result->data = surf_base;

    return (PyObject*)result;
  }

  // extract basis functions and evaluation points
  int dimension = surf_base->dimension();
  Go::BsplineBasis basis_u = surf_base->basis_u();
  Go::BsplineBasis basis_v = surf_base->basis_v();
  std::vector<double> greville_u(basis_u.numCoefs());
  std::vector<double> greville_v(basis_v.numCoefs());
  for(int i=0; i<basis_u.numCoefs(); i++)
    greville_u[i] = basis_u.grevilleParameter(i);
  for(int i=0; i<basis_v.numCoefs(); i++)
    greville_v[i] = basis_v.grevilleParameter(i);

  // evaluate original spline at interpolation points
  std::vector<double> interpolationPoints;
  std::vector<double> weights(0);
  surf_base->gridEvaluator(interpolationPoints, greville_u, greville_v);

  Go::SplineSurface *nonrational_surf = Go::SurfaceInterpolator::regularInterpolation(basis_u,
                                                                                      basis_v,
                                                                                      greville_u,
                                                                                      greville_v,
                                                                                      interpolationPoints,
                                                                                      dimension,
                                                                                      false,
                                                                                      weights);

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data = shared_ptr<Go::SplineSurface>(nonrational_surf);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_plane__doc__,"Generate an infinite plane\n"
                                   "@param p0: A point on the plane\n"
                                   "@type center: Point, list of floats or tuple of floats\n"
                                   "@param normal: The normal of the plane\n"
                                   "@type center: Point, list of floats or tuple of floats\n"
                                   "@return: Surface describing the plane");
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

PyDoc_STRVAR(generate_rectangle__doc__,"Generate a rectangle\n"
                                       "@param corner: A corner of the rectangle\n"
                                       "@type corner: Point, list of floats or tuple of floats\n"
                                       "@param u_axis: The first axis of the rectangle\n"
                                       "@type u_axis: Point, list of floats or tuple of floats\n"
                                       "@param v_axis: The second axis of the rectangle\n"
                                       "@type v_axis: Point, list of floats or tuple of floats\n"
                                       "@param u_length: Length of sides along u-axis\n"
                                       "@type u_length: float\n"
                                       "@param v_length: Length of sides along v-axis\n"
                                       "@type v_length: float\n"
                                       "@return: Surface describing the rectangle");
PyObject* Generate_Rectangle(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"corner", "u_axis", "v_axis", "u_length", "v_length", NULL };
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
  result->data.reset(new Go::Plane(*corner,*axisx % dir_v_ortog, *axisx));
  static_pointer_cast<Go::Plane>(result->data)->setParameterBounds(0.0,0.0,length_x,length_y);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_rotational_curve_sweep__doc__,"Generate a surface by rotating a curve\n"
                                                    "@param curve: The curve to rotate\n"
                                                    "@type curve: Curve\n"
                                                    "@param pos: The point to rotate the curve about\n"
                                                    "@type pos: Point, list of floats or tuple of floats\n"
                                                    "@param axis: The axis to rotate the curve about\n"
                                                    "@type axis: Point, list of floats or tuple of floats\n"
                                                    "@param angle: (optional )The angle to rotate the curve in radians (defaults to a full rotation)\n"
                                                    "@type angle: float\n"
                                                    "@return: Surface describing the swept surface");
PyObject* Generate_RotationalCurveSweep(PyObject* self, PyObject* args, PyObject* kwds)
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

PyDoc_STRVAR(generate_sphere_surface__doc__,"Generate a sphere surface\n"
                                            "@param center: The center of the sphere\n"
                                            "@type center: Point, list of floats or tuple of floats\n"
                                            "@param radius: The radius of the sphere\n"
                                            "@type radius: float\n"
                                            "@return: Surface describing the sphere surface");
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

PyDoc_STRVAR(generate_torus_surface__doc__,"Generate a torus surface\n"
                                           "@param center: The center of the torus\n"
                                           "@type center: Point, list of floats or tuple of floats\n"
                                           "@param axis: The center axis of the torus\n"
                                           "@type center: Point, list of floats or tuple of floats\n"
                                           "@param major_radius: The major radius of the torus\n"
                                           "@type major_radius: float\n"
                                           "@param minor_radius: The minor radius of the torus\n"
                                           "@type minor_radius: float\n"
                                           "@return: Surface describing the torus surface");
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

PyDoc_STRVAR(generate_trim_surface__doc__,"Generate a surface by trimming one surface with another\n"
                                          "@param original: The surface to trim\n"
                                          "@type original: Surface\n"
                                          "@param trim: The surface to trim with\n"
                                          "@type trim: Surface\n"
                                          "@return: Surface describing the trimmed surface");
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

PyMethodDef SurfaceFactory_methods[] = {
     {(char*)"AddLoop",               (PyCFunction)Generate_AddLoop,              METH_VARARGS|METH_KEYWORDS, generate_addloop__doc__},
     {(char*)"CircularDisc",          (PyCFunction)Generate_CircularDisc,         METH_VARARGS|METH_KEYWORDS, generate_circular_disc__doc__},
     {(char*)"ConeSurface",           (PyCFunction)Generate_ConeSurface,          METH_VARARGS|METH_KEYWORDS, generate_cone_surface__doc__},
     {(char*)"ConvertNonRational",    (PyCFunction)Generate_NonRational,          METH_VARARGS|METH_KEYWORDS, generate_nonrational__doc__},
     {(char*)"CoonsSurfacePatch",     (PyCFunction)Generate_CoonsSurfacePatch,    METH_VARARGS|METH_KEYWORDS, generate_coons_surface_patch__doc__},
     {(char*)"CylinderSurface",       (PyCFunction)Generate_CylinderSurface,      METH_VARARGS|METH_KEYWORDS, generate_cylinder_surface__doc__},
     {(char*)"LinearCurveSweep",      (PyCFunction)Generate_LinearCurveSweep,     METH_VARARGS|METH_KEYWORDS, generate_linear_curve_sweep__doc__},
     {(char*)"LoftCurves",            (PyCFunction)Generate_LoftCurves,           METH_VARARGS|METH_KEYWORDS, generate_loft_curves__doc__},
     {(char*)"Plane",                 (PyCFunction)Generate_Plane,                METH_VARARGS|METH_KEYWORDS, generate_plane__doc__},
     {(char*)"Rectangle",             (PyCFunction)Generate_Rectangle,            METH_VARARGS|METH_KEYWORDS, generate_rectangle__doc__},
     {(char*)"RotationalCurveSweep",  (PyCFunction)Generate_RotationalCurveSweep, METH_VARARGS|METH_KEYWORDS, generate_rotational_curve_sweep__doc__},
     {(char*)"SphereSurface",         (PyCFunction)Generate_SphereSurface,        METH_VARARGS|METH_KEYWORDS, generate_sphere_surface__doc__},
     {(char*)"TorusSurface",          (PyCFunction)Generate_TorusSurface,         METH_VARARGS|METH_KEYWORDS, generate_torus_surface__doc__},
     {(char*)"TrimSurface",           (PyCFunction)Generate_TrimSurface,          METH_VARARGS|METH_KEYWORDS, generate_trim_surface__doc__},
     {NULL,                           NULL,                                       0,                          NULL}
  };

PyDoc_STRVAR(surface_factory__doc__,"A module with methods for generating surfaces");

PyMODINIT_FUNC
init_SurfaceFactory_Module()
{
  SurfaceFactory_module = Py_InitModule3((char*)"GoTools.SurfaceFactory",SurfaceFactory_methods,surface_factory__doc__);
  modState.addInfo(SurfaceFactory_module);
}

}

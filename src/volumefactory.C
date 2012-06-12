#include "volumefactory.h"

#include "curve.h"
#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"
#include "surface.h"
#include "volume.h"

#include "GoTools/trivariate/ConeVolume.h"
#include "GoTools/trivariate/CylinderVolume.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/trivariate/Parallelepiped.h"
#include "GoTools/trivariate/SphereVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/TorusVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
#include "GoTools/geometry/GeometryTools.h"

extern "C" {
PyObject* VolumeFactory_module;

PyDoc_STRVAR(generate_box__doc__, "Generate a box\n"
                                  "@param corner: A corner of the box\n"
                                  "@type corner: Point, list of floats or tuple of floats\n"
                                  "@param u_axis: First axis of the box\n"
                                  "@type u_axis: Point, list of floats or tuple of floats\n"
                                  "@param v_axis: Second axis of the box\n"
                                  "@type v_axis: Point, list of floats or tuple of floats\n"
                                  "@param u_length: Length of the box along the first axis\n"
                                  "@type u_length: float\n"
                                  "@param v_length: Length of the box along the second axis\n"
                                  "@type v_length: float\n"
                                  "@param w_length: Length of the box along the third axis\n"
                                  "@type w_length: float\n"
                                  "@return: Volume describing the box");
PyObject* Generate_Box(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"corner", "u_axis", "v_axis", "u_length", "v_length", "w_length", NULL };
  PyObject* cornero;
  PyObject* axisxo;
  PyObject* axisyo;
  double length_x, length_y, length_z;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOddd",
                                   (char**)keyWords,&cornero,&axisxo,
                                   &axisyo,&length_x,&length_y,&length_z))
    return NULL;
       
  shared_ptr<Go::Point> corner = PyObject_AsGoPoint(cornero);
  shared_ptr<Go::Point> axisx  = PyObject_AsGoPoint(axisxo);
  shared_ptr<Go::Point> axisy  = PyObject_AsGoPoint(axisyo);
  if (!corner || !axisx || !axisy)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  Go::Point dir_v_ortog = *axisy-(*axisx*(*axisy))/axisx->length2()*(*axisx);
  result->data.reset(new Go::Parallelepiped(*corner,*axisx, dir_v_ortog,
                                            *axisx % dir_v_ortog, 
                                            length_x, length_y, length_z));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_cone__doc__, "Generate a cone\n"
                                   "@param apex: The apex of the cone\n"
                                   "@type apex: Point, list of floats or tuple of floats\n"
                                   "@param axis: The axis of the cone\n"
                                   "@type axis: Point, list of floats or tuple of floats\n"
                                   "@param angle: The angle of the cone in radians\n"
                                   "@type angle: float\n"
                                   "@param height: (optional) The height of the cone\n"
                                   "@type height: float\n"
                                   "@param radius: (optional) The radius of the cone\n"
                                   "@type radius: float\n"
                                   "@return: Volume describing the cone");
PyObject* Generate_Cone(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"apex", "axis", "angle", "height", "radius", NULL };
  PyObject* apexo;
  PyObject* axiso;
  double angle, height = std::numeric_limits<double>::infinity(), radius=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOd|dd",
                                   (char**)keyWords,&apexo,&axiso,
                                                    &angle,&height,&radius))
    return NULL;
       
  shared_ptr<Go::Point> apex = PyObject_AsGoPoint(apexo);
  shared_ptr<Go::Point> axis  = PyObject_AsGoPoint(axiso);
  if (!apex || !axis)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  result->data.reset(new Go::ConeVolume(radius, *apex, *axis,
                                        someNormal(*axis), angle));
  static_pointer_cast<Go::ConeVolume>(result->data)->setParameters(0.0,height,2);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_contract_surf_to__doc__,"Generate a volume by contracting a surface to a point\n"
                                              "@param surface: The surface\n"
                                              "@type surface: Surface\n"
                                              "@param point: The point to contract to\n"
                                              "@type point: Point, list of floats or tuple of floats\n"
                                              "@return: Volume with resulting geometry");
PyObject* Generate_ContractSurfaceTo(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "point", NULL };
  PyObject* surfaceo;
  PyObject* pointo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&surfaceo,&pointo))
    return NULL;

  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(PyObject_AsGoSurface(surfaceo));
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(pointo);
  if (!surf || !point)
    return NULL;

  // create "basis" for the point
  double knots[] = {0,0,1,1};

  int dim = surf->dimension()+surf->rational();
  int n1  = surf->numCoefs_u();
  int n2  = surf->numCoefs_v();
  int p1  = surf->order_u();
  int p2  = surf->order_v();

  std::vector<double> coefs(n1*n2*dim*2);
  std::copy(surf->ctrl_begin(),surf->ctrl_end(),coefs.begin());

  int j=n1*n2*dim;
  for (size_t i=0; i<n1*n2; i++) {
    coefs[j++] = (*point)[0];
    coefs[j++] = (*point)[1];
    if (surf->dimension() > 2)
      coefs[j++] = (*point)[2];
    if (surf->rational())
      coefs[j++] = 1; 
  }

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data.reset(new Go::SplineVolume(n1, n2, 2, p1, p2, 2,
                                          surf->basis_u().begin(), surf->basis_v().begin(), knots,
                                          coefs.begin(), surf->dimension(),
                                          surf->rational()));
  return (PyObject*)result;
}

PyDoc_STRVAR(generate_cylinder__doc__, "Generate a cylinder\n"
                                       "@param center: The center of the cylinder\n"
                                       "@type center: Point, list of floats or tuple of floats\n"
                                       "@param boundarypoint: A point on the boundary of the cylinder\n"
                                       "@type boundarypoint: Point, list of floats or tuple of floats\n"
                                       "@param normal: The normal of the cylinder\n"
                                       "@type center: Point, list of floats or tuple of floats\n"
                                       "@param height: (optional) The height of the cylinder\n"
                                       "@type height: float\n"
                                       "@return: Volume describing the cylinder");
PyObject* Generate_Cylinder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "boundarypoint", "normal", "height", NULL };
  PyObject* centero;
  PyObject* bpointo;
  PyObject* normalo;
  double height = -1.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO|d",
                                   (char**)keyWords,&centero,&bpointo,
                                                    &normalo,&height))
    return NULL;
       
  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> bpoint = PyObject_AsGoPoint(bpointo);
  shared_ptr<Go::Point> normal = PyObject_AsGoPoint(normalo);
  if (!center || !bpoint || !normal)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  Go::Point v = *bpoint-*center;
  Go::Point x_axis = v- ((v*(*normal))/normal->length2())*(*normal);
  double rad = x_axis.length();
  result->data.reset(new Go::CylinderVolume(*center,rad,*normal,x_axis));
  if (height > 0)
    static_pointer_cast<Go::CylinderVolume>(result->data)->setParameters(0,height,2);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_extrude_surface__doc__, "Generate a volume by extruding along a given direction\n"
                                              "@param surface: The surface to extrude\n"
                                              "@type surface: Surface\n"
                                              "@param direction: The direction to extrude\n"
                                              "@type direction: Point, list of floats or tuple of floats\n"
                                              "@param amount: (optional) Scaling factor of direction\n"
                                              "@type amount: float\n"
                                              "@return: Volume describing the extruded volume");
PyObject* Generate_ExtrudeSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "direction", "amount", NULL };
  PyObject* surfo;
  PyObject* pointo;
  double amount = 1.0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO|d",
                                   (char**)keyWords,&surfo,&pointo,&amount))
    return NULL;

  shared_ptr<Go::Point> point = PyObject_AsGoPoint(pointo);
  if (!point)
    return NULL;

  Volume* result = NULL;
  shared_ptr<Go::ParamSurface> s1 = PyObject_AsGoSurface(surfo);
  if (!s1)
    s1 = PyObject_AsGoSurface(surfo);

  
  Go::Point extrude = *point;
  extrude *= amount;
  shared_ptr<Go::ParamSurface> s2(s1->clone());
  Go::GeometryTools::translateSplineSurf(extrude,
                              *static_pointer_cast<Go::SplineSurface>(s2));
  

  if (s1 && s2) {
    std::vector<shared_ptr<Go::SplineSurface> > surfaces;
    surfaces.push_back(static_pointer_cast<Go::SplineSurface>(s1));
    surfaces.push_back(static_pointer_cast<Go::SplineSurface>(s2));
    result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
    result->data.reset(Go::LoftVolumeCreator::loftVolume(surfaces.begin(),surfaces.size()));
  }

  return (PyObject*)result;
}


PyDoc_STRVAR(generate_loft_surfaces__doc__, "Generate a volume by lofting surfaces\n"
                                            "@param surfaces: The surfaces to loft\n"
                                            "@type center: List of Surface\n"
                                            "@return: Volume describing the lofted surfaces");
PyObject* Generate_LoftSurfaces(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surfaces", NULL };
  PyObject* surfaceso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&surfaceso))
    return NULL;

  if (!PyObject_TypeCheck(surfaceso,&PyList_Type))
    return NULL;

  std::vector<shared_ptr<Go::SplineSurface> > surfaces;
  for (int i=0; i < PyList_Size(surfaceso); ++i) {
    PyObject* surfo = PyList_GetItem(surfaceso,i);
    shared_ptr<Go::ParamSurface> surf = PyObject_AsGoSurface(surfo);
    if (!surf)
      continue;
    surfaces.push_back(convertSplineSurface(surf));
  }
  if (surfaces.size() < 2)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data.reset(Go::LoftVolumeCreator::loftVolume(surfaces.begin(),surfaces.size()));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_linear_surface_sweep__doc__, "Generate a volume by sweeping a surface along a curve\n"
                                                   "@param surface: The surface to sweep along\n"
                                                   "@type surface: Surface\n"
                                                   "@param curve: The curve to sweep\n"
                                                   "@type curve: Curve\n"
                                                   "@param point: The starting point of the sweep\n"
                                                   "@type point: Point, list of floats or tuple of floats\n"
                                                   "@return: Volume describing the swept volume");
PyObject* Generate_LinearSurfaceSweep(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "curve", "point", NULL };
  PyObject* o1;
  PyObject* o2;
  PyObject* pointo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO",
                                   (char**)keyWords,&o1,&o2,&pointo))
    return NULL;

  shared_ptr<Go::Point> point = PyObject_AsGoPoint(pointo);
  if (!point)
    return NULL;

  Volume* result = NULL;
  shared_ptr<Go::ParamSurface> s = PyObject_AsGoSurface(o1);
  if (!s)
    s = PyObject_AsGoSurface(o2);
  shared_ptr<Go::ParamCurve> c = PyObject_AsGoCurve(o2);
  if (!c)
    c = PyObject_AsGoCurve(o1);

  if (s && c) {
    result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
    result->data.reset(Go::SweepVolumeCreator::linearSweptVolume(*convertSplineSurface(s),
                                                                   *convertSplineCurve(c),*point));
  }

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_volnonrational__doc__,"Generate a non-rational representation (approximation) of a rational spline volume\n"
                                            "@param original: The initial (rational) volume\n"
                                            "@type original: Volume\n"
                                            "@return: Non-rational B-spline representation of the volume");
PyObject* Generate_VolNonRational(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"original", NULL };
  PyObject* originalo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&originalo))
    return NULL;

  shared_ptr<Go::ParamVolume> vol = PyObject_AsGoVolume(originalo);
  if (!vol)
    return NULL;

  shared_ptr<Go::SplineVolume> vol_base = convertSplineVolume(vol);

  if (!vol_base)
    return NULL;

  // if it's already B-spline, just return itself
  if(!vol_base->rational()) {
    Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
    result->data = vol_base;

    return (PyObject*)result;
  }

  // extract basis functions and evaluation points
  int dimension = vol_base->dimension();
  Go::BsplineBasis basis_u = vol_base->basis(0);
  Go::BsplineBasis basis_v = vol_base->basis(1);
  Go::BsplineBasis basis_w = vol_base->basis(2);
  std::vector<double> greville_u(basis_u.numCoefs());
  std::vector<double> greville_v(basis_v.numCoefs());
  std::vector<double> greville_w(basis_w.numCoefs());
  for(int i=0; i<basis_u.numCoefs(); i++)
    greville_u[i] = basis_u.grevilleParameter(i);
  for(int i=0; i<basis_v.numCoefs(); i++)
    greville_v[i] = basis_v.grevilleParameter(i);
  for(int i=0; i<basis_w.numCoefs(); i++)
    greville_w[i] = basis_w.grevilleParameter(i);

  // evaluate original spline at interpolation points
  std::vector<double> interpolationPoints;
  std::vector<double> weights(0);
  vol_base->gridEvaluator(greville_u, greville_v, greville_w,
                          interpolationPoints);

  Go::SplineVolume *nonrational_vol = 
    Go::VolumeInterpolator::regularInterpolation(basis_u,
                                                 basis_v,
                                                 basis_w,
                                                 greville_u,
                                                 greville_v,
                                                 greville_w,
                                                 interpolationPoints,
                                                 dimension,
                                                 false,
                                                 weights);

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data = shared_ptr<Go::SplineVolume>(nonrational_vol);

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_rotational_surface_sweep__doc__, "Generate a volume by rotating a surface\n"
                                                       "@param surface: The surface to rotate\n"
                                                       "@type surface: Surface\n"
                                                       "@param position: The starting position of the surface\n"
                                                       "@type position: Point, list of floats or tuple of floats\n"
                                                       "@param axis: The axis to rotate the around\n"
                                                       "@type axis: Point, list of floats or tuple of floats\n"
                                                       "@param angle: (optional) The angle to rotate with in radians. Defaults to a full rotation\n"
                                                       "@type angle: float\n"
                                                       "@return: Volume describing the rotated surface");
PyObject* Generate_RotationalSurfaceSweep(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "position", "axis", "angle", NULL };
  PyObject* surfaceo;
  PyObject* poso;
  PyObject* axiso;
  double angle=2*M_PI;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO|d",
                                   (char**)keyWords,&surfaceo,&poso,
                                                    &axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(surfaceo);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  shared_ptr<Go::Point> pos = PyObject_AsGoPoint(poso);
  if (!surface || !axis || !pos)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data.reset(Go::SweepVolumeCreator::rotationalSweptVolume(*convertSplineSurface(surface),
                                                                   angle,*pos,*axis));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_parallelepiped__doc__, "Generate a parallelepiped\n"
                                             "@param corner: A corner of the parallelepiped\n"
                                             "@type corner: Point, list of floats or tuple of floats\n"
                                             "@param u_axis: First axis of the parallelepiped\n"
                                             "@type u_axis: Point, list of floats or tuple of floats\n"
                                             "@param v_axis: Second axis of the parallelepiped\n"
                                             "@type v_axis: Point, list of floats or tuple of floats\n"
                                             "@param v_axis: Third axis of the parallelepiped\n"
                                             "@type v_axis: Point, list of floats or tuple of floats\n"
                                             "@param u_length: Length of the box along the first axis\n"
                                             "@type u_length: float\n"
                                             "@param v_length: Length of the box along the second axis\n"
                                             "@type v_length: float\n"
                                             "@param w_length: Length of the box along the third axis\n"
                                             "@type w_length: float\n"
                                             "@return: Volume describing the parallellepiped");
PyObject* Generate_Parallelepiped(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"corner", "u_axis", "v_axis", "w_axis", "u_length", "v_length", "w_length", NULL };
  PyObject* cornero;
  PyObject* axisxo;
  PyObject* axisyo;
  PyObject* axiszo;
  double length_x, length_y, length_z;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOOddd",
                                   (char**)keyWords,&cornero,&axisxo,
                                                    &axisyo,&axiszo,&length_x,
                                                    &length_y,&length_z))
    return NULL;
       
  shared_ptr<Go::Point> corner = PyObject_AsGoPoint(cornero);
  shared_ptr<Go::Point> axisx  = PyObject_AsGoPoint(axisxo);
  shared_ptr<Go::Point> axisy  = PyObject_AsGoPoint(axisyo);
  shared_ptr<Go::Point> axisz  = PyObject_AsGoPoint(axiszo);
  if (!corner || !axisx || !axisy || !axisz)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  result->data.reset(new Go::Parallelepiped(*corner,*axisx, *axisy, *axisz, 
                                            length_x, length_y, length_z));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_sphere__doc__, "Generate a sphere\n"
                                     "@param center: Center of the sphere\n"
                                     "@type center: Point, list of floats or tuple of floats\n"
                                     "@param radius: Radius of the sphere\n"
                                     "@type radius: float\n"
                                     "@return: Volume describing the sphere");
PyObject* Generate_Sphere(PyObject* self, PyObject* args, PyObject* kwds)
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

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);

  Go::Point x_axis(1.f,0.f,0.f);
  Go::Point z_axis(0.f,0.f,1.f);
  
  result->data.reset(new Go::SphereVolume(radius,*center,x_axis,z_axis));

  return (PyObject*)result;
}

PyDoc_STRVAR(generate_torus__doc__, "Generate a torus\n"
                                    "@param center: Center of the torus\n"
                                    "@type center: Point, list of floats or tuple of floats\n"
                                    "@param axis: Center axis of the torus\n"
                                    "@type axis: Point, list of floats or tuple of floats\n"
                                    "@param major_radius: Major radius of the torus\n"
                                    "@type major_radius: float\n"
                                    "@param minor_radius: Minor radius of the torus\n"
                                    "@type minor_radius: float\n"
                                    "@return: Volume describing the torus");
PyObject* Generate_Torus(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "axis", "major_radius", "minor_radius", NULL };
  PyObject* centero;
  PyObject* axiso;
  double major_radius, minor_radius;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOdd",
                                   (char**)keyWords,&centero,&axiso,
                                                    &major_radius,&minor_radius))
    return NULL;
       
  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis   = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);

  result->data.reset(new Go::TorusVolume(major_radius,minor_radius, *center,
                                         *axis, someNormal(*axis)));

  return (PyObject*)result;
}

PyMethodDef VolumeFactory_methods[] = {
     {(char*)"Box",                   (PyCFunction)Generate_Box,                    METH_VARARGS|METH_KEYWORDS, generate_box__doc__},
     {(char*)"Cone",                  (PyCFunction)Generate_Cone,                   METH_VARARGS|METH_KEYWORDS, generate_cone__doc__},
     {(char*)"ContractSurfaceTo",     (PyCFunction)Generate_ContractSurfaceTo,      METH_VARARGS|METH_KEYWORDS, generate_contract_surf_to__doc__},
     {(char*)"Cylinder",              (PyCFunction)Generate_Cylinder,               METH_VARARGS|METH_KEYWORDS, generate_cylinder__doc__},
     {(char*)"ExtrudeSurface",        (PyCFunction)Generate_ExtrudeSurface,         METH_VARARGS|METH_KEYWORDS, generate_extrude_surface__doc__},
     {(char*)"LoftSurfaces",          (PyCFunction)Generate_LoftSurfaces,           METH_VARARGS|METH_KEYWORDS, generate_loft_surfaces__doc__},
     {(char*)"LinearSurfaceSweep",    (PyCFunction)Generate_LinearSurfaceSweep,     METH_VARARGS|METH_KEYWORDS, generate_linear_surface_sweep__doc__},
     {(char*)"NonRationalVolume",     (PyCFunction)Generate_VolNonRational,         METH_VARARGS|METH_KEYWORDS, generate_volnonrational__doc__},
     {(char*)"RotationalSurfaceSweep",(PyCFunction)Generate_RotationalSurfaceSweep, METH_VARARGS|METH_KEYWORDS, generate_rotational_surface_sweep__doc__},
     {(char*)"Parallelepiped",        (PyCFunction)Generate_Parallelepiped,         METH_VARARGS|METH_KEYWORDS, generate_parallelepiped__doc__},
     {(char*)"Sphere",                (PyCFunction)Generate_Sphere,                 METH_VARARGS|METH_KEYWORDS, generate_sphere__doc__},
     {(char*)"Torus",                 (PyCFunction)Generate_Torus,                  METH_VARARGS|METH_KEYWORDS, generate_torus__doc__},
     {NULL,                           NULL,                                         0,                          NULL}
   };

PyDoc_STRVAR(volume_factory__doc__,"A module with methods for generating volumes");

PyMODINIT_FUNC
init_VolumeFactory_Module()
{
  VolumeFactory_module = Py_InitModule3((char*)"GoTools.VolumeFactory",VolumeFactory_methods,volume_factory__doc__);
  modState.addInfo(VolumeFactory_module);
}

}

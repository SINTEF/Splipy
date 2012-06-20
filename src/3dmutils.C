#ifdef ENABLE_OPENNURBS

#include "curve.h"
#include "surface.h"
#include "pyutils.h"

#include <GoTools/geometry/SplineSurface.h>
#include <opennurbs.h>

Curve* ONCurveToGoCurve(const ON_Curve* curve)
{
  ON_NurbsCurve* n_curve = curve->NurbsCurve();
  if (!n_curve)
    return NULL;

  ON_Interval domain = n_curve->Domain();
  if (n_curve->IsPeriodic()) {
    if (!n_curve->InsertKnot(domain.Min(), n_curve->Order()-1) ||
        !n_curve->InsertKnot(domain.Max(), n_curve->Order()-1))
      return NULL;
  }
  double knot_vector[n_curve->KnotCount()+2]; // g2 is operating with knot multiplicity Order() at the start and end knot, openNURBS is not
  double *coefs = NULL; // g2 is storing controlpoints as P1_x*w,P1_y*w,P1_z*w,w P2_x*w,... that is each controlpoint is premultiplied by its weight. openNURBS is not
  if (n_curve->IsRational()) {
    coefs = new double[n_curve->CVCount() * n_curve->CVSize()];
    int k=0;
    for (int j=0; j<n_curve->CVCount(); j++) {
      for (int d=0; d<n_curve->Dimension(); d++)
        coefs[k++] = n_curve->CV(j)[d] * n_curve->Weight(j);
      coefs[k++] = n_curve->Weight(j);
    }
  } else
    coefs = n_curve->m_cv;

  knot_vector[0] = n_curve->SuperfluousKnot(0);
  std::copy(n_curve->Knot(), n_curve->Knot() + n_curve->KnotCount(), knot_vector+1);
  knot_vector[n_curve->KnotCount()+1] = n_curve->SuperfluousKnot(1);

  Curve* result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  result->data.reset(new Go::SplineCurve(n_curve->CVCount(), n_curve->Order(),
                                         knot_vector, coefs, 
                                         n_curve->Dimension(), 
                                         n_curve->IsRational()));

  if (n_curve->IsRational())
    delete[] coefs;

  return result;
}

ON_NurbsCurve* GoCurveToONCurve(PyObject* curve)
{
  shared_ptr<Go::SplineCurve> crv = 
        convertSplineCurve(PyObject_AsGoCurve(curve));
  if (!crv)
    return NULL;

  ON_NurbsCurve* result = new ON_NurbsCurve(crv->dimension(), crv->rational(),
                                            crv->order(), crv->numCoefs());

  for (int j=1;j<crv->numCoefs()+crv->order()-1;++j)
    result->SetKnot(j-1,(*(crv->knotsBegin()+j)));

  // WARNING: assuming dim=3 atm
  if (crv->rational()) {
    std::vector<double>::const_iterator cit = crv->rcoefs_begin();
    for (int i=0;i<crv->numCoefs();++i) {
      double x = *cit++;
      double y = *cit++;
      double z = *cit++;
      double w = *cit++;
      ON_4dPoint p(x/w, y/w, z/w, w);
      result->SetCV(i,p);
      cit++;
    }
  } else {
    std::vector<double>::const_iterator cit = crv->coefs_begin();
    for (int i=0;i<crv->numCoefs();++i) {
      double x=*cit++;
      double y=*cit++;
      double z=*cit++;
      ON_3dPoint p(x,y,z);
      result->SetCV(i,p);
    }
  }

  return result;
}

Surface* ONSurfaceToGoSurface(const ON_Surface* surf)
{
  ON_NurbsSurface* n_surf= NULL;
  if (surf->HasNurbForm())
    n_surf = surf->NurbsSurface();
  else {
    n_surf = ON_NurbsSurface::New();
    surf->GetNurbForm(*n_surf);
  }
  if (!n_surf)
    return NULL;

  // openNURBS is storing major index first, goTools minor index
  double x_comp[n_surf->m_cv_count[0]][n_surf->m_cv_count[1]];
  double y_comp[n_surf->m_cv_count[0]][n_surf->m_cv_count[1]];
  double z_comp[n_surf->m_cv_count[0]][n_surf->m_cv_count[1]];
  double w_comp[n_surf->m_cv_count[0]][n_surf->m_cv_count[1]];
              
  int k=0;
  for(int i=0; i<n_surf->m_cv_count[0]; i++) {
    for(int j=0; j<n_surf->m_cv_count[1]; j++) {
      x_comp[i][j] = n_surf->m_cv[k++];
      y_comp[i][j] = n_surf->m_cv[k++];
      z_comp[i][j] = n_surf->m_cv[k++];
      if (n_surf->m_is_rat)
        w_comp[i][j] = n_surf->m_cv[k++];
    }
  }

  // then interleave
  double comp[n_surf->m_cv_count[0]*n_surf->m_cv_count[1]*(n_surf->m_is_rat?4:3)];
  k = 0;
  for(int j=0; j<n_surf->m_cv_count[1]; j++) {
    for(int i=0; i<n_surf->m_cv_count[0]; i++) {
      comp[k++] = x_comp[i][j];
      comp[k++] = y_comp[i][j];
      comp[k++] = z_comp[i][j];
      if (n_surf->m_is_rat)
        comp[k++] = w_comp[i][j];
    }
  }

  // fix up knot vectors
  double u[n_surf->m_knot_capacity[0]+2];
  u[0] = n_surf->m_knot[0][0];
  std::copy(n_surf->m_knot[0],n_surf->m_knot[0]+n_surf->m_knot_capacity[0], u+1);
  u[n_surf->m_knot_capacity[0]+1] = n_surf->m_knot[0][n_surf->m_knot_capacity[0]-1];

  double v[n_surf->m_knot_capacity[1]+2];
  v[0] = n_surf->m_knot[1][0];
  std::copy(n_surf->m_knot[1],n_surf->m_knot[1]+n_surf->m_knot_capacity[1], v+1);
  v[n_surf->m_knot_capacity[1]+1] = n_surf->m_knot[1][n_surf->m_knot_capacity[1]-1];

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data.reset(new Go::SplineSurface(n_surf->m_cv_count[0],
                                           n_surf->m_cv_count[1],
                                           n_surf->m_order[0],
                                           n_surf->m_order[1],
                                           u,v,comp,3,n_surf->m_is_rat));

  if (!surf->HasNurbForm())
    delete n_surf;

  return result;
}

ON_NurbsSurface* GoSurfaceToONSurface(PyObject* surf)
{
  shared_ptr<Go::SplineSurface> srf = 
        convertSplineSurface(PyObject_AsGoSurface(surf));
  if (!srf)
    return NULL;

  ON_NurbsSurface* result = new ON_NurbsSurface(srf->dimension(),
                                                srf->rational(),
                                                srf->order_u(),
                                                srf->order_v(),
                                                srf->numCoefs_u(),
                                                srf->numCoefs_v());

  for (int j=1;j<srf->numCoefs_u()+srf->order_u()-1;++j)
    result->SetKnot(0,j-1,(*srf->basis_u().begin()+j));
  for (int j=1;j<srf->numCoefs_v()+srf->order_v()-1;++j)
    result->SetKnot(1,j-1,(*srf->basis_v().begin()+j));
  
  // WARNING: Assuming dim = 3
  std::vector<double>::iterator cit = srf->ctrl_begin();
  for (int i=0;i<srf->numCoefs_u();++i) {
    for (int j=0;j<srf->numCoefs_v();++j) {
      if (srf->rational()) {
        double x = *cit++;
        double y = *cit++;
        double z = *cit++;
        double w = *cit++;
        ON_4dPoint p(x/w, y/w, z/w, w);
        result->SetCV(i,j,p);
      } else {
        double x = *cit++;
        double y = *cit++;
        double z = *cit++;
        ON_3dPoint p(x,y,z);
        result->SetCV(i,j,p);
      }
    }
  }

  return result;
}

static void WriteEntity3DM(ONX_Model& model, PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&Curve_Type)) {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = GoCurveToONCurve(obj);
    mo.m_bDeleteObject = true;
    mo.m_attributes.m_layer_index = 0;
  }
  if (PyObject_TypeCheck(obj,&Surface_Type)) {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = GoSurfaceToONSurface(obj);
    mo.m_bDeleteObject = true;
    mo.m_attributes.m_layer_index = 0;
  }
}

void DoWrite3DM(const std::string& fname, PyObject* objectso)
{
  ONX_Model model; 
  model.m_properties.m_RevisionHistory.NewRevision();
  model.m_properties.m_Application.m_application_name = "GeoModeller";
  model.m_properties.m_Application.m_application_URL = "http://sintef.no";
  model.m_properties.m_Application.m_application_details = "The SINTEF ICT GoTools/OpenNURBS based python bindings";
  
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0;
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01;

  if (PyObject_TypeCheck(objectso,&PyList_Type)) {
    for (int i=0; i < PyList_Size(objectso); ++i) {
      PyObject* obj = PyList_GetItem(objectso,i);
      WriteEntity3DM(model,obj);
    }
  } else if (PyObject_TypeCheck(objectso,&PyTuple_Type)) {
    for (int i=0;i < PyTuple_Size(objectso); ++i) {
      PyObject* obj = PyTuple_GetItem(objectso,i);
      WriteEntity3DM(model,obj);
    }
    PyTuple_GetItem(objectso,0);
  } else
    WriteEntity3DM(model,objectso);

  FILE* fp = ON::OpenFile(fname.c_str(), "wb");
  ON_BinaryFile archive(ON::write3dm, fp);
  ON_TextLog error_log;
  model.Polish();
  model.Write(archive,4,"",&error_log);
  ON::CloseFile(fp);
}

PyObject* DoRead3DM(const std::string& fname, const std::string& type)
{
  FILE* archive_fp = ON::OpenFile(fname.c_str(),"rb");
  if (!archive_fp)
    return NULL;

  ON_BinaryFile archive(ON::read3dm, archive_fp);
  ONX_Model model;

  ON_TextLog dump_to_stdout;
  if (!model.Read(archive,&dump_to_stdout)) {
    std::cerr << "Failed to read 3DM file " << fname << std::endl;
    return NULL;
  }
  ON::CloseFile(archive_fp);

  bool surfaces(true);
  bool curves(true);
  if (!type.empty() && type == "curves")
    surfaces = false;
  if (type == "surfaces")
    curves = false;

  PyObject* result=NULL;
  PyObject* curr=NULL;
  for (int i=0;i<model.m_object_table.Count();++i) {
    if (curr) {
      if (!result)
        result = PyList_New(0);
      PyList_Append(result,curr);
      curr = NULL;
    }
    if (curves && model.m_object_table[i].m_object->ObjectType() == ON::curve_object)
      curr = (PyObject*)ONCurveToGoCurve((const ON_Curve*)model.m_object_table[i].m_object);
    if (surfaces && model.m_object_table[i].m_object->ObjectType() == ON::brep_object) {
      ON_Brep *brep_obj = (ON_Brep*) model.m_object_table[i].m_object;
      int sc=0;
      while (brep_obj->FaceIsSurface(sc)) {
        if (curr) {
          if (!result)
            result = PyList_New(0);
          PyList_Append(result,curr);
        }
        ON_Surface *surf = brep_obj->m_S[sc++];
        curr = (PyObject*)ONSurfaceToGoSurface(surf);
      }
    }
  }
  model.Destroy();
  if (!result)
    result = curr;
  if (PyObject_TypeCheck(result,&PyList_Type))
    PyList_Append(result,curr);

  return result;
}


#endif

#ifdef ENABLE_OPENNURBS

#include "curve.h"
#include "surface.h"
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
#endif

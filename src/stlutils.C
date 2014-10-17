#include "stlutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

namespace GeoModeller {

int WriteSurfaceSTL(std::ofstream& stl_file, Surface* surface,
                    bool ascii, int res[2])
{
  int nTriangles = 0;
  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(surface->data);

  std::vector<double> pts;
  std::vector<double> normals;
  std::vector<double> param_u;
  std::vector<double> param_v;
  surf->gridEvaluator(res[0], res[1], pts, normals, param_u, param_v, false);
  int ip = 0;
  // loop over all squares
  for(int j=0; j<res[1]-1; j++)
    for(int i=0; i<res[0]-1; i++, ip++)
    {
      int i1 =   j  *res[0] + ( i );
      int i2 =   j  *res[0] + (i+1);
      int i3 = (j+1)*res[0] + ( i );
      int i4 = (j+1)*res[0] + (i+1);
      Go::Point p1(pts.begin()+3*i1, pts.begin()+3*i1+3);
      Go::Point p2(pts.begin()+3*i2, pts.begin()+3*i2+3);
      Go::Point p3(pts.begin()+3*i3, pts.begin()+3*i3+3);
      Go::Point p4(pts.begin()+3*i4, pts.begin()+3*i4+3);

      // evaluate the normal of the first triangle as the average of the three corners
      double n1[] = { (normals[i1*3  ] + normals[i2*3  ] + normals[i3*3  ])/3.0,
                      (normals[i1*3+1] + normals[i2*3+1] + normals[i3*3+1])/3.0,
                      (normals[i1*3+2] + normals[i2*3+2] + normals[i3*3+2])/3.0};
      // and the second triangle
      double n2[] = { (normals[i2*3  ] + normals[i3*3  ] + normals[i4*3  ])/3.0,
                      (normals[i2*3+1] + normals[i3*3+1] + normals[i4*3+1])/3.0,
                      (normals[i2*3+2] + normals[i3*3+2] + normals[i4*3+2])/3.0};
      int t1[] = {i1,i3,i2}; // first triangle coordinates
      int t2[] = {i2,i3,i4}; // second triangle coordinates

      if(ascii) {
        if( ((p3-p1)%(p2-p1)).length()/2.0 > modState.gapTolerance ) {
          stl_file << "facet normal ";
          for(int k=0; k<3; k++)
            stl_file << n1[k] << " " ;
          stl_file << std::endl;

          stl_file << "outer loop" << std::endl;
          for(int t=0; t<3; t++) // for all triangle corners
          {
            stl_file << "vertex ";
            for(int k=0; k<3; k++) // for (x,y,z) coordinate
              stl_file << pts[ t1[t]*3 + k ] << " " ;
            stl_file << std::endl;
          }
          stl_file << "endloop" << std::endl << "endfacet" << std::endl;
          nTriangles++;
        }
        if( ((p4-p2)%(p3-p2)).length()/2.0 > modState.gapTolerance ) {
          stl_file << "facet normal ";
          for(int k=0; k<3; k++)
            stl_file << n2[k] << " " ;
          stl_file << std::endl;

          stl_file << "outer loop" << std::endl;
          for(int t=0; t<3; t++) // for all triangle corners
          {
            stl_file << "vertex ";
            for(int k=0; k<3; k++) // for (x,y,z) coordinate
              stl_file << pts[ t2[t]*3 + k ] << " " ;
            stl_file << std::endl;
          }
          stl_file << "endloop" << std::endl << "endfacet" << std::endl;
          nTriangles++;
        }

      }
      else  // Writing binary file output
      { 
        // storing stuff as float representation REAL32 in binary
        std::vector<float> fpts;
        float fn1[3];
        float fn2[3];
        fpts.resize(    pts.size()    );
        std::copy(pts.begin(), pts.end(), fpts.begin());
        std::copy(n1,          n1+3,      fn1 );
        std::copy(n2,          n2+3,      fn2 );
        unsigned short noAttribute = 0;

        if( ((p3-p1)%(p2-p1)).length()/2.0 > modState.gapTolerance ) {
          stl_file.write((char*) fn1, 3*sizeof(float)) ;
  
          for(int t=0; t<3; t++) // for all triangle corners
            stl_file.write((char*) &(fpts[ t1[t]*3 ]), 3*sizeof(float)) ;
          stl_file.write((char*) &noAttribute, sizeof(unsigned short) );
          nTriangles++;
        }
  
        if( ((p4-p2)%(p3-p2)).length()/2.0 > modState.gapTolerance ) {
          stl_file.write((char*) fn2, 3*sizeof(float)) ;
  
          for(int t=0; t<3; t++) // for all triangle corners
            stl_file.write((char*) &(fpts[ t2[t]*3 ]), 3*sizeof(float)) ;
          stl_file.write((char*) &noAttribute, sizeof(unsigned short) );
          nTriangles++;
        }
      }
  }
  return nTriangles;
}

int WriteVolumeSTL(std::ofstream& stl_file, Volume* volume,
                   bool ascii, int res[3])
{
  int nTriangles = 0;
  shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume->data);
  if (vol->isLeftHanded()) {
    vol = shared_ptr<Go::SplineVolume>(vol->clone());
    vol->reverseParameterDirection(2);
  }

  std::vector<shared_ptr<Go::ParamSurface> > edges = vol->getAllBoundarySurfaces();

  edges[1]->reverseParameterDirection(true);
  edges[2]->reverseParameterDirection(true);
  edges[5]->reverseParameterDirection(true);
  
  for(int i=0; i<6; i++) {
    Surface* surface = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    surface->data = shared_ptr<Go::ParamSurface>(edges.at(i));

    int surfRes[] = { res[i/2==0], res[2-(i/2==2)] };
    nTriangles += WriteSurfaceSTL(stl_file, surface, ascii, surfRes);
  }
  return nTriangles;
}

}

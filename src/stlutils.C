#include "stlutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"

#include <array>

namespace GeoModeller {

struct QuadTesselated {
  std::array<Go::Point, 4> point;
  std::array<std::array<int, 3>, 2> indices;
  std::array<std::array<double, 3>, 2> normal;
};

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
  // loop over all squares
  std::vector<QuadTesselated> tess((res[1]-1)*(res[0]-1));
#pragma omp parallel for schedule(static)
  for(int j=0; j<res[1]-1; j++) {
    int ip = j*(res[0]-1);
    for(int i=0; i<res[0]-1; i++, ip++) {
      int i1 =   j  *res[0] + ( i );
      int i2 =   j  *res[0] + (i+1);
      int i3 = (j+1)*res[0] + ( i );
      int i4 = (j+1)*res[0] + (i+1);

      tess[ip].point[0] = Go::Point(pts.begin()+3*i1, pts.begin()+3*i1+3);
      tess[ip].point[1] = Go::Point(pts.begin()+3*i2, pts.begin()+3*i2+3);
      tess[ip].point[2] = Go::Point(pts.begin()+3*i3, pts.begin()+3*i3+3);
      tess[ip].point[3] = Go::Point(pts.begin()+3*i4, pts.begin()+3*i4+3);

      tess[ip].normal[0] = { (normals[i1*3  ] + normals[i2*3  ] + normals[i3*3  ])/3.0,
                             (normals[i1*3+1] + normals[i2*3+1] + normals[i3*3+1])/3.0,
                             (normals[i1*3+2] + normals[i2*3+2] + normals[i3*3+2])/3.0 };
      // and the second triangle
      tess[ip].normal[1] = { (normals[i2*3  ] + normals[i3*3  ] + normals[i4*3  ])/3.0,
                             (normals[i2*3+1] + normals[i3*3+1] + normals[i4*3+1])/3.0,
                             (normals[i2*3+2] + normals[i3*3+2] + normals[i4*3+2])/3.0};
      tess[ip].indices[0] = {i1,i3,i2}; // first triangle coordinates
      tess[ip].indices[1] = {i2,i3,i4}; // second triangle coordinates
    }
  }

  std::vector<float> fpts;
  if (!ascii) {
    fpts.resize(    pts.size()    );
    std::copy(pts.begin(), pts.end(), fpts.begin());
  }

  int ip = 0;
  for(int j=0; j<res[1]-1; j++) {
    for(int i=0; i<res[0]-1; i++, ip++) {
      if(ascii) {
        if (((tess[ip].point[2]-tess[ip].point[0])%(tess[ip].point[1]-tess[ip].point[0])).length()/2.0 > modState.gapTolerance) {
          stl_file << "facet normal ";
          for(int k=0; k<3; k++)
            stl_file << tess[ip].normal[0][k] << " " ;
          stl_file << std::endl;

          stl_file << "outer loop" << std::endl;
          for(int t=0; t<3; t++) // for all triangle corners
          {
            stl_file << "vertex ";
            for(int k=0; k<3; k++) // for (x,y,z) coordinate
              stl_file << pts[ tess[ip].indices[0][t]*3 + k ] << " " ;
            stl_file << std::endl;
          }
          stl_file << "endloop" << std::endl << "endfacet" << std::endl;
          nTriangles++;
        }
        if (((tess[ip].point[3]-tess[ip].point[1])%(tess[ip].point[2]-tess[ip].point[1])).length()/2.0 > modState.gapTolerance) {
          stl_file << "facet normal ";
          for(int k=0; k<3; k++)
            stl_file << tess[ip].normal[1][k] << " " ;
          stl_file << std::endl;

          stl_file << "outer loop" << std::endl;
          for(int t=0; t<3; t++) { // for all triangle corners
            stl_file << "vertex ";
            for(int k=0; k<3; k++) // for (x,y,z) coordinate
              stl_file << pts[ tess[ip].indices[1][t]*3 + k ] << " " ;
            stl_file << std::endl;
          }
          stl_file << "endloop" << std::endl << "endfacet" << std::endl;
          nTriangles++;
        }
      } else {// Writing binary file output
        // storing stuff as float representation REAL32 in binary
        float fn1[3];
        float fn2[3];
        unsigned short noAttribute = 0;

        if (((tess[ip].point[2]-tess[ip].point[0])%(tess[ip].point[1]-tess[ip].point[0])).length()/2.0 > modState.gapTolerance) {
          std::copy(tess[ip].normal[0].begin(), tess[ip].normal[0].end(), fn1);
          stl_file.write((char*) fn1, 3*sizeof(float)) ;
  
          for(int t=0; t<3; t++) // for all triangle corners
            stl_file.write((char*) &(fpts[ tess[ip].indices[0][t]*3 ]), 3*sizeof(float)) ;
          stl_file.write((char*) &noAttribute, sizeof(unsigned short) );
          nTriangles++;
        }
  
        if (((tess[ip].point[3]-tess[ip].point[1])%(tess[ip].point[2]-tess[ip].point[1])).length()/2.0 > modState.gapTolerance) {
          std::copy(tess[ip].normal[1].begin(), tess[ip].normal[1].end(), fn2);
          stl_file.write((char*) fn2, 3*sizeof(float)) ;
  
          for(int t=0; t<3; t++) // for all triangle corners
            stl_file.write((char*) &(fpts[ tess[ip].indices[1][t]*3 ]), 3*sizeof(float)) ;
          stl_file.write((char*) &noAttribute, sizeof(unsigned short) );
          nTriangles++;
        }
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

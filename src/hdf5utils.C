#ifdef ENABLE_HDF5

#include "curve.h"
#include "geomodeller.h"
#include "surface.h"
#include "pyutils.h"

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <hdf5.h>

#include <sstream>
#include <sys/stat.h>

static bool checkGroupExistence(int parent, const char* path)
{
  bool result = false;
  // turn off errors to avoid cout spew
  H5E_BEGIN_TRY {
    result = H5Gget_objinfo((hid_t)parent,path,0,NULL) == 0;
  } H5E_END_TRY;

  return result;
}

PyObject* DoReadHDF5Field(const std::string& fname, const std::string& fldname,
                          int patch, int level)
{
  hid_t file = H5Fopen(fname.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  std::stringstream str;
  str << level << "/" << patch;
  hid_t group = H5Gopen2(file,str.str().c_str(),H5P_DEFAULT);

  hid_t set = H5Dopen2(group,fldname.c_str(),H5P_DEFAULT);
  hsize_t siz = H5Dget_storage_size(set) / 8;
  double* data = new double[siz];
  H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(set);
  H5Gclose(group);
  H5Fclose(file);

  PyObject* result = PyList_New(0); 
  for (int i=0;i<siz;++i)
    PyList_Append(result, PyFloat_FromDouble(data[i]));

  delete[] data;

  return result;
}

PyObject* DoReadHDF5Geometry(const std::string& fname, const std::string& geoname,
                             int patch, int level)
{
  hid_t file = H5Fopen(fname.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  std::stringstream str;
  str << level << "/basis/" << geoname;
  hid_t group = H5Gopen2(file,str.str().c_str(),H5P_DEFAULT);

  std::stringstream str2;
  str2 << patch;
  hid_t set = H5Dopen2(group,str2.str().c_str(),H5P_DEFAULT);
  hsize_t siz = H5Dget_storage_size(set);
  char* temp = new char[siz];
  H5Dread(set,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,temp);
  H5Dclose(set);
  H5Gclose(group);
  H5Fclose(file);
  std::stringstream f;
  f << temp;
  delete[] temp;

  Go::ObjectHeader header;
  f >> header;
  if (header.classType() == Go::SplineCurve::classType())
    return ReadG2<Curve, Go::SplineCurve>(f,Curve_Type);
  if (header.classType() == Go::SplineSurface::classType())
    return ReadG2<Surface, Go::SplineSurface>(f,Surface_Type);
  if (header.classType() == Go::SplineVolume::classType())
    return ReadG2<Volume, Go::SplineVolume>(f,Volume_Type);

  return NULL;
}

void DoWriteHDF5Field(const std::string& fname, const std::string& fldname,
                      int patch, int level, PyObject* data, bool trunc)
{
  struct stat fb;
  hid_t file;
  if (trunc || stat(fname.c_str(),&fb) == -1) {
    file = H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  } else
    file = H5Fopen(fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  std::stringstream str;
  str << level; 
  if (!checkGroupExistence(file, str.str().c_str())) {
    hid_t grp = H5Gcreate2(file,str.str().c_str(),0, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(grp);
  }
  str << "/" <<patch;
  hid_t group;
  if (checkGroupExistence(file, str.str().c_str()))
    group = H5Gopen2(file,str.str().c_str(),H5P_DEFAULT);
  else
    group = H5Gcreate2(file,str.str().c_str(),0, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t siz = PyList_Size(data);
  double* wrData = new double[siz];
  for (int i=0;i<siz;++i)
    wrData[i] = PyFloat_AsDouble(PyList_GetItem(data,i));

  hid_t space = H5Screate_simple(1,&siz,NULL);
  hid_t set = H5Dcreate2(group,fldname.c_str(),
                         H5T_NATIVE_DOUBLE,space,
                         H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t file_space = H5Dget_space(set);
  hsize_t stride = 1;
  hsize_t start = 0;
  H5Sselect_hyperslab(file_space,H5S_SELECT_SET,&start,&stride,&siz,NULL);
  hid_t mem_space = H5Screate_simple(1,&siz,NULL);
  H5Dwrite(set,H5T_NATIVE_DOUBLE,mem_space,file_space,H5P_DEFAULT,wrData);
  H5Sclose(mem_space);
  H5Sclose(file_space);
  H5Dclose(set);
  H5Sclose(space);
  H5Gclose(group);
  H5Fclose(file);
  delete[] wrData;
}

void DoWriteHDF5Geometry(const std::string& fname, const std::string& geoname,
                         int patch, int level, PyObject* data, bool trunc)
{
  struct stat fb;
  hid_t file;
  if (trunc || stat(fname.c_str(),&fb) == -1) {
    file = H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  } else
    file = H5Fopen(fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  std::stringstream str;
  str << level; 
  if (!checkGroupExistence(file, str.str().c_str())) {
    hid_t grp = H5Gcreate2(file,str.str().c_str(),0, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(grp);
  }
  str << "/basis";
  if (!checkGroupExistence(file, str.str().c_str())) {
    hid_t grp = H5Gcreate2(file,str.str().c_str(),0, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(grp);
  }
  str << "/" + geoname;
  hid_t group;
  if (checkGroupExistence(file, str.str().c_str()))
    group = H5Gopen2(file,str.str().c_str(),H5P_DEFAULT);
  else
    group = H5Gcreate2(file,str.str().c_str(),0, H5P_DEFAULT, H5P_DEFAULT);

  std::stringstream strd;
  if (PyObject_TypeCheck(data, &Curve_Type)) {
    ((Curve*)data)->data->writeStandardHeader(strd);
    ((Curve*)data)->data->write(strd);
  } else if (PyObject_TypeCheck(data, &Surface_Type)) {
    ((Surface*)data)->data->writeStandardHeader(strd);
    ((Surface*)data)->data->write(strd);
  } else if (PyObject_TypeCheck(data, &Volume_Type)) {
    ((Volume*)data)->data->writeStandardHeader(strd);
    ((Volume*)data)->data->write(strd);
  }

  hsize_t siz = strd.str().size();

  std::stringstream strp;
  strp << patch;

  hid_t space = H5Screate_simple(1,&siz,NULL);
  hid_t set = H5Dcreate2(group,strp.str().c_str(),
                         H5T_NATIVE_CHAR,space,
                         H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t file_space = H5Dget_space(set);
  hsize_t stride = 1;
  hsize_t start = 0;
  H5Sselect_hyperslab(file_space,H5S_SELECT_SET,&start,&stride,&siz,NULL);
  hid_t mem_space = H5Screate_simple(1,&siz,NULL);
  H5Dwrite(set,H5T_NATIVE_CHAR,mem_space,file_space,H5P_DEFAULT,strd.str().c_str());
  H5Sclose(mem_space);
  H5Sclose(file_space);
  H5Dclose(set);
  H5Sclose(space);
  H5Gclose(group);
  H5Fclose(file);
}

#endif

#define PY_SSIZE_T_CLEAN
#include "Python.h"

#include "scatmech.h"
#include "brdf.h"
#include "local.h"
#include "crossrcw.h"
#include "rcw.h"
//#include "getmodel.h"
#include "sphrscat.h"
#include <sstream>
#include <iostream>
using namespace std;

using namespace SCATMECH;

#define HERE() cerr<<"I am here: " << __LINE__ << endl;

std::map<int,Model*> mapModel;
std::map<int,BRDF_Model_Ptr> mapBRDF_Model; 
std::map<int,Local_BRDF_Model_Ptr> mapLocal_BRDF_Model;
std::map<int,Model_Ptr<RCW_Model> > mapRCW_Model;
std::map<int,Model_Ptr<CrossRCW_Model> > mapCrossRCW_Model;
std::map<int,Free_Space_Scatterer_Ptr> mapFSS_Model;
std::map<int,StackModel_Ptr> mapStackModel; 
std::map<int,Model_Ptr<Model> > mapLone_Model;
	       
int modelct = 0;

class SCATPY_exception: public std::runtime_error
{
public:
  SCATPY_exception(const std::string& _msg) : std::runtime_error(_msg) {}
};


// Convert std::complex to Py_complex...
static PyObject * PYCOMPLEX(const COMPLEX& x)
{
  Py_complex result;
  result.real = real(x);
  result.imag = imag(x);
  return PyComplex_FromCComplex(result);
}

// Convert a JonesMatrix to a PyObject...
static PyObject * PyJones(const JonesMatrix& jones)
{
  // This had to be done manually (not using Py_Builder), because someone reported
  // a segmentation fault using it with complex numbers.
    COMPLEX ss = jones.SS(), sp = jones.SP(), ps = jones.PS(), pp = jones.PP();

    PyObject *_ss = PYCOMPLEX(ss);
    PyObject *_sp = PYCOMPLEX(sp);
    PyObject *_ps = PYCOMPLEX(ps);
    PyObject *_pp = PYCOMPLEX(pp);
    
    PyObject *inner1 = PyTuple_New(2);
    PyTuple_SetItem(inner1, 0, _ss);
    PyTuple_SetItem(inner1, 1, _ps);
    PyObject *inner2 = PyTuple_New(2);
    PyTuple_SetItem(inner2, 0, _sp);
    PyTuple_SetItem(inner2, 1, _pp);

    PyObject *result = PyTuple_New(2);
    PyTuple_SetItem(result, 0, inner1);
    PyTuple_SetItem(result, 1, inner2);

    return result;
}  

// Convert a MuellerMatrix to a PyObject...
static PyObject * PyMueller(const MuellerMatrix& m)
{
    PyObject *result = PyTuple_New(4);
    for (int i=0;i<4;++i) {
      PyObject *row = PyTuple_New(4);
      for (int j=0;j<4;++j) {
	PyTuple_SetItem(row, j, PyFloat_FromDouble(m[i][j]));
      }
      PyTuple_SetItem(result, i, row);
    }
    return result;
}

// Convert a Vector to a PyObject...
static PyObject * PyVector(const Vector& v)
{
  PyObject *result = PyTuple_New(3);
  PyTuple_SetItem(result, 0, PyFloat_FromDouble(v.x));
  PyTuple_SetItem(result, 1, PyFloat_FromDouble(v.y));
  PyTuple_SetItem(result, 2, PyFloat_FromDouble(v.z));
  return result;
}

static PyObject * Get_Model(PyObject *self, PyObject *args)
{
  try { 
    const char *modelname;

    if (!PyArg_ParseTuple(args, "s", &modelname))
        return NULL;

    Model_Ptr<Model> newModel = std::string(modelname);
    mapLone_Model[++modelct] = newModel;
    mapModel[modelct] = mapLone_Model[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Free_Model(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapLone_Model.erase(handle);
    
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Get_BRDF_Model(PyObject *self, PyObject *args)
{
  try {
    const char *modelname;

    if (!PyArg_ParseTuple(args, "s", &modelname))
        return NULL;

    BRDF_Model_Ptr newModel = std::string(modelname);
    mapBRDF_Model[++modelct] = newModel;
    mapModel[modelct] = mapBRDF_Model[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Free_BRDF_Model(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapBRDF_Model.erase(handle);
    
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Get_Local_BRDF_Model(PyObject *self, PyObject *args)
{
  try {
    const char *modelname;

    if (!PyArg_ParseTuple(args, "s", &modelname))
        return NULL;

    Local_BRDF_Model_Ptr newModel = std::string(modelname);
    mapLocal_BRDF_Model[++modelct] = newModel;
    mapModel[modelct] = mapLocal_BRDF_Model[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Free_Local_BRDF_Model(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapLocal_BRDF_Model.erase(handle);
    
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Get_StackModel(PyObject *self, PyObject *args)
{
  try {
    const char *modelname;

    if (!PyArg_ParseTuple(args, "s", &modelname))
        return NULL;

    StackModel_Ptr newModel = std::string(modelname);
    mapStackModel[++modelct] = newModel;
    mapModel[modelct] = mapStackModel[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Free_StackModel(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapStackModel.erase(handle);
    
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Get_Free_Space_Scatterer(PyObject *self, PyObject *args)
{
  try {
    const char *modelname;

    if (!PyArg_ParseTuple(args, "s", &modelname))
        return NULL;

    Free_Space_Scatterer_Ptr newModel = std::string(modelname);
    mapFSS_Model[++modelct] = newModel;
    mapModel[modelct] = mapFSS_Model[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Free_Free_Space_Scatterer(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapFSS_Model.erase(handle);
    
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Get_RCW_Model(PyObject *self, PyObject *args)
{
  try {
    Model_Ptr<RCW_Model> newModel = std::string("RCW_Model");
    mapRCW_Model[++modelct] = newModel;
    mapModel[modelct] = mapRCW_Model[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
    return NULL;
  }    
}

static PyObject * Free_RCW_Model(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapRCW_Model.erase(handle);
      
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Get_CrossRCW_Model(PyObject *self, PyObject *args)
{
  try {
    Model_Ptr<CrossRCW_Model> newModel = std::string("CrossRCW_Model");
    mapCrossRCW_Model[++modelct] = newModel;
    mapModel[modelct] = mapCrossRCW_Model[modelct].get();
      
    return PyLong_FromLong(modelct);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Free_CrossRCW_Model(PyObject *self, PyObject *args)
{
  try {
    int handle;

    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    mapModel.erase(handle);
    mapCrossRCW_Model.erase(handle);
      
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * BRDF(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double thetai,thetas,phis;
    double rotation=0;
    const char *coords = "psps";
    
    if (!PyArg_ParseTuple(args, "iddd|ds", &handle, &thetai, &thetas, &phis, &rotation, &coords))
        return NULL;

    BRDF_Model::Coordinate_System _coords;
    if (std::string(coords)=="psps") _coords = BRDF_Model::psps;
    else if (std::string(coords)=="xyxy") _coords = BRDF_Model::xyxy;
    else if (std::string(coords)=="plane") _coords = BRDF_Model::plane;
    else return NULL;
    
    MuellerMatrix mueller = mapBRDF_Model[handle]->Mueller(thetai,thetas,phis,rotation,_coords);

    return PyMueller(mueller);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * BRDFJones(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double thetai,thetas,phis;
    double rotation=0;
    const char *coords = "psps";
    
    if (!PyArg_ParseTuple(args, "iddd|ds", &handle, &thetai, &thetas, &phis, &rotation, &coords))
        return NULL;

    BRDF_Model::Coordinate_System _coords;
    if (std::string(coords)=="psps") _coords = BRDF_Model::psps;
    else if (std::string(coords)=="xyxy") _coords = BRDF_Model::xyxy;
    else if (std::string(coords)=="plane") _coords = BRDF_Model::plane;
    else return NULL;
    
    JonesMatrix jones = mapBRDF_Model[handle]->Jones(thetai,thetas,phis,rotation,_coords);

    return PyJones(jones);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * VectoredBRDF(PyObject *self, PyObject *args)
{
  try {
    int handle;
    Vector source, viewer, normal, xaxis;
    const char *coords = "plane";
    
    if (!PyArg_ParseTuple(args, "idddddddddddd|s", &handle,
			  &source.x, &source.y, &source.z,
			  &viewer.x, &viewer.y, &viewer.z,
			  &normal.x, &normal.y, &normal.z,
			  &xaxis.x, &xaxis.y, &xaxis.z, &coords))
        return NULL;

    BRDF_Model::Coordinate_System _coords;
    if (std::string(coords)=="psps") _coords = BRDF_Model::psps;
    else if (std::string(coords)=="xyxy") _coords = BRDF_Model::xyxy;
    else if (std::string(coords)=="plane") _coords = BRDF_Model::plane;
    else return NULL;
    
    MuellerMatrix mueller = mapBRDF_Model[handle]->Mueller(source,viewer,normal,xaxis,_coords);

    return PyMueller(mueller);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * VectoredBRDFJones(PyObject *self, PyObject *args)
{
  try {
    int handle;
    Vector source, viewer, normal, xaxis;
    const char *coords = "plane";
    
    if (!PyArg_ParseTuple(args, "idddddddddddd|s", &handle,
			  &source.x, &source.y, &source.z,
			  &viewer.x, &viewer.y, &viewer.z,
			  &normal.x, &normal.y, &normal.z,
			  &xaxis.x, &xaxis.y, &xaxis.z, &coords))
        return NULL;

    BRDF_Model::Coordinate_System _coords;
    if (std::string(coords)=="psps") _coords = BRDF_Model::psps;
    else if (std::string(coords)=="xyxy") _coords = BRDF_Model::xyxy;
    else if (std::string(coords)=="plane") _coords = BRDF_Model::plane;
    else return NULL;
    
    JonesMatrix jones = mapBRDF_Model[handle]->Jones(source,viewer,normal,xaxis,_coords);

    return PyJones(jones);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * LocalDSC(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double thetai,thetas,phis;
    double rotation=0;
    const char *coords = "psps";
    
    if (!PyArg_ParseTuple(args, "iddd|ds", &handle, &thetai, &thetas, &phis, &rotation, &coords))
        return NULL;

    BRDF_Model::Coordinate_System _coords;
    if (std::string(coords)=="psps") _coords = BRDF_Model::psps;
    else if (std::string(coords)=="xyxy") _coords = BRDF_Model::xyxy;
    else if (std::string(coords)=="plane") _coords = BRDF_Model::plane;
    else return NULL;
    
    MuellerMatrix mueller = mapLocal_BRDF_Model[handle]->MuellerDSC(thetai,thetas,phis,rotation,_coords);

    return PyMueller(mueller);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * LocalDSCJones(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double thetai,thetas,phis;
    double rotation=0;
    const char *coords = "psps";
    
    if (!PyArg_ParseTuple(args, "iddd|ds", &handle, &thetai, &thetas, &phis, &rotation, &coords))
        return NULL;

    BRDF_Model::Coordinate_System _coords;
    if (std::string(coords)=="psps") _coords = BRDF_Model::psps;
    else if (std::string(coords)=="xyxy") _coords = BRDF_Model::xyxy;
    else if (std::string(coords)=="plane") _coords = BRDF_Model::plane;
    else return NULL;
    
    JonesMatrix jones = mapLocal_BRDF_Model[handle]->JonesDSC(thetai,thetas,phis,rotation,_coords);

    return PyJones(jones);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}


static PyObject * FSSjones(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double kix,kiy,kiz,ksx,ksy,ksz;

    if (!PyArg_ParseTuple(args, "idddddd", &handle, &kix, &kiy, &kiz, &ksx, &ksy, &ksz))
        return NULL;
    
    Free_Space_Scatterer_Ptr &model = mapFSS_Model[handle];

    double lambda = model->get_lambda();
    double medium = model->get_medium().n(lambda);
    double k = 2.*pi/lambda*medium;  
    JonesMatrix jones = model->jones(Vector(kix,kiy,kiz),Vector(ksx,ksy,ksz))/k;

    return PyJones(jones);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
} 

static PyObject * FSSext(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double kx,ky,kz;

    if (!PyArg_ParseTuple(args, "iddd", &handle, &kx, &ky, &kz))
        return NULL;

    Vector k(kx,ky,kz);
    Free_Space_Scatterer_Ptr &model = mapFSS_Model[handle];

    MuellerMatrix mueller = model->extinction(k);

    return PyMueller(mueller);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * RCWDiffractionEfficiency(PyObject *self, PyObject *args)
{
  try {
    int handle,i;
    
    if (!PyArg_ParseTuple(args, "ii", &handle, &i))
        return NULL;

    MuellerMatrix mueller = mapRCW_Model[handle]->GetIntensity(i);

    return PyMueller(mueller);
		  
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * RCWDiffractionAmplitude(PyObject *self, PyObject *args)
{
  try {
    int handle,i;
    
    if (!PyArg_ParseTuple(args, "ii", &handle, &i))
        return NULL;

    JonesMatrix jones = mapRCW_Model[handle]->GetAmplitude(i);

    return PyJones(jones);
		  
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * CrossRCWDiffractionEfficiency(PyObject *self, PyObject *args)
{
  try {
    int handle,i,j;
    
    if (!PyArg_ParseTuple(args, "iii", &handle, &i,&j))
        return NULL;

    MuellerMatrix mueller = mapCrossRCW_Model[handle]->GetIntensity(i,j);

    return PyMueller(mueller);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * CrossRCWDiffractionAmplitude(PyObject *self, PyObject *args)
{
  try {
    int handle,i,j;
    
    if (!PyArg_ParseTuple(args, "iii", &handle, &i,&j))
        return NULL;

    JonesMatrix jones = mapCrossRCW_Model[handle]->GetAmplitude(i,j);

    return PyJones(jones);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * RCWDirection(PyObject *self, PyObject *args)
{
  try {
    int handle,i;
    
    if (!PyArg_ParseTuple(args, "ii", &handle, &i))
        return NULL;

    Vector direction = mapRCW_Model[handle]->GetDirection(i);

    return PyVector(direction);
		  
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * GetGratingEpsilon(PyObject *self, PyObject *args)
{
  try {
    int handle;
    double x,z;
    int direction = 0;
    
    if (!PyArg_ParseTuple(args, "idd|i", &handle, &x, &z, &direction))
        return NULL;

    COMPLEX eps;
    
    const Grating_Ptr &grating = mapRCW_Model[handle]->get_grating(); 

    int level = grating->get_level(z);

    double lambda = mapRCW_Model[handle]->get_lambda();

    if (grating->get_lambda()!=lambda) {
      grating->set_lambda(lambda);
    }
    if (level==-1) {
      eps = grating->get_medium_i().epsilon(lambda);
    } else if (level==-2) eps = grating->get_medium_t().epsilon(lambda);
    else eps = grating->eps(x,level,direction);

    return PyComplex_FromDoubles(real(eps),imag(eps));
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * GetGratingDefinition(PyObject *self, PyObject *args)
{
  try {

    int handle;
    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    Model_Ptr<RCW_Model> &model = mapRCW_Model[handle];
    const Grating_Ptr &grating = model->get_grating();

    typedef std::vector<std::vector<double> > vvdouble;
    typedef std::vector<std::vector<COMPLEX> > vvcomplex;
    typedef std::vector<double> vdouble;
    typedef std::vector<COMPLEX> vcomplex;

    grating->set_lambda(model->get_lambda());
    
    const vvdouble & cposition = grating->get_position();
    const vvcomplex & cmaterialx = grating->get_materialx();
    const vvcomplex & cmaterialy = grating->get_materialy();
    const vvcomplex & cmaterialz = grating->get_materialz();
    const vvcomplex & cmaterialmux = grating->get_materialmux();
    const vvcomplex & cmaterialmuy = grating->get_materialmuy();
    const vvcomplex & cmaterialmuz = grating->get_materialmuz();
    const vdouble & cthickness = grating->get_thickness();

    PyObject *position = PyList_New(0);
    for (vvdouble::const_iterator p=cposition.begin(); p!=cposition.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vdouble::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PyFloat_FromDouble(*q));
      }
      PyList_Append(position, A);
    }

    PyObject *materialx = PyList_New(0);
    for (vvcomplex::const_iterator p=cmaterialx.begin(); p!=cmaterialx.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vcomplex::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PYCOMPLEX(*q));
      }
      PyList_Append(materialx, A);
    }
    
    PyObject *materialy = PyList_New(0);
    for (vvcomplex::const_iterator p=cmaterialy.begin(); p!=cmaterialy.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vcomplex::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PYCOMPLEX(*q));
      }
      PyList_Append(materialy, A);
    }
    
    PyObject *materialz = PyList_New(0);
    for (vvcomplex::const_iterator p=cmaterialz.begin(); p!=cmaterialz.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vcomplex::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PYCOMPLEX(*q));
      }
      PyList_Append(materialz, A);
    }

    PyObject *materialmux = PyList_New(0);
    for (vvcomplex::const_iterator p=cmaterialmux.begin(); p!=cmaterialmux.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vcomplex::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PYCOMPLEX(*q));
      }
      PyList_Append(materialmux, A);
    }
    
    PyObject *materialmuy = PyList_New(0);
    for (vvcomplex::const_iterator p=cmaterialmuy.begin(); p!=cmaterialmuy.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vcomplex::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PYCOMPLEX(*q));
      }
      PyList_Append(materialmuy, A);
    }
    
    PyObject *materialmuz = PyList_New(0);
    for (vvcomplex::const_iterator p=cmaterialmuz.begin(); p!=cmaterialmuz.end(); ++p) {
      PyObject *A = PyList_New(0);
      for (vcomplex::const_iterator q=p->begin(); q!=p->end(); ++q) {
	PyList_Append(A, PYCOMPLEX(*q));
      }
      PyList_Append(materialmuz, A);
    }
     
    PyObject *thickness = PyList_New(0);
    for (vdouble::const_iterator p=cthickness.begin(); p!=cthickness.end(); ++p) {
      PyList_Append(thickness, PyFloat_FromDouble(*p));
    }
    
    PyObject *dict = PyDict_New();

    PyDict_SetItemString(dict, "position", position);
    PyDict_SetItemString(dict, "materialx", materialx);
    PyDict_SetItemString(dict, "materialy", materialy);
    PyDict_SetItemString(dict, "materialz", materialz);
    PyDict_SetItemString(dict, "materialmux", materialx);
    PyDict_SetItemString(dict, "materialmuy", materialy);
    PyDict_SetItemString(dict, "materialmuz", materialz);
    PyDict_SetItemString(dict, "thickness", thickness);
    PyDict_SetItemString(dict, "medium_i", PYCOMPLEX(grating->get_medium_i().epsilon(grating->get_lambda())));
    PyDict_SetItemString(dict, "medium_t", PYCOMPLEX(grating->get_medium_t().epsilon(grating->get_lambda())));
    PyDict_SetItemString(dict, "period", PyFloat_FromDouble(grating->get_period()));

    cerr << ">>>>>>>" << grating->get_inheritance().get_name() << endl;
    if (grating->get_inheritance().get_name() == "Generic_Grating") {
      Generic_Grating* _grating = dynamic_cast<Generic_Grating*>(grating.get());

      Generic_Grating::segmentvector segs = _grating->Get_Boundaries();

      PyObject *seglist = PyList_New(0);
      
      for (Generic_Grating::segmentvector::iterator v = segs.begin(); v != segs.end(); ++v) {	
	PyObject *segment = PyDict_New();
	PyDict_SetItemString(segment, "x1", PyFloat_FromDouble(real(v->vertex1)));
	PyDict_SetItemString(segment, "y1", PyFloat_FromDouble(imag(v->vertex1)));
	PyDict_SetItemString(segment, "x2", PyFloat_FromDouble(real(v->vertex2)));
	PyDict_SetItemString(segment, "y2", PyFloat_FromDouble(imag(v->vertex2)));
	PyDict_SetItemString(segment, "mat1", PyUnicode_FromString(v->mat1.c_str()));
	PyDict_SetItemString(segment, "mat2", PyUnicode_FromString(v->mat2.c_str()));
	PyList_Append(seglist, segment);
      }
      PyDict_SetItemString(dict, "segments", seglist);      
    }
    
    return dict;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * CrossRCWDirection(PyObject *self, PyObject *args)
{
  try {
    int handle,i,j;
    
    if (!PyArg_ParseTuple(args, "iii", &handle, &i,&j))
        return NULL;

    Vector direction = mapCrossRCW_Model[handle]->GetDirection(i,j);

    return PyVector(direction);
		  
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * SetParameter(PyObject *self, PyObject *args)
{
  try {
    const char *parameter;
    const char *value=0;
    int handle;
    
    if (!PyArg_ParseTuple(args, "iss", &handle, &parameter, &value))
        return NULL;
 
    mapModel[handle]->set_parameter(parameter,value);
    
    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;  
  }     
}

static PyObject * GetParameter(PyObject *self, PyObject *args)
{
  try {
    const char *parameter;
    int handle;
    
    if (!PyArg_ParseTuple(args, "is", &handle, &parameter))
        return NULL;

    std::string value = mapModel[handle]->get_parameter(parameter);

    return PyUnicode_FromString(value.c_str());
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * PrintParameters(PyObject *self, PyObject *args)
{
  try {
    int handle;
    
    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    std::ostringstream os;
    mapModel[handle]->print_parameters(os);

    return PyUnicode_FromString(os.str().c_str());
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * AskParameters(PyObject *self, PyObject *args)
{
  try {
    int handle;
    
    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    Model* model = mapModel[handle];
    
    model->AskUser();

    Py_INCREF(Py_None);
    return Py_None;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * GetParameterDictionary(PyObject *self, PyObject *args)
{
  try {
    int handle;
    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    ModelParameterList mpl;

    typedef std::deque<std::string> StringList;
    StringList names;

    Model *model = mapModel[handle];
    model->get_parameter_names(names);

    PyObject *dict = PyDict_New();
    PyObject *dictsub;

    for (StringList::iterator q=names.begin();q!=names.end();++q) {

       ParameterInfo info = model->get_parameter_info((*q));
       std::string value = model->get_parameter((*q));

       dictsub = Py_BuildValue("{s:s,s:s,s:s,s:s}",
			       "value",value.c_str(),
			       "description",(info.description).c_str(),
			       "type", (info.type).c_str(),
			       "default", (info.defaultvalue).c_str());
       PyDict_SetItemString(dict,q->c_str(),dictsub);
       Py_DECREF(&dictsub);
    }
    return dict;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * GetModelName(PyObject *self, PyObject *args)
{
  try {
    int handle;
    
    if (!PyArg_ParseTuple(args, "i", &handle))
        return NULL;

    Model* model = mapModel[handle];
    
    
    return PyUnicode_FromString(model->get_inheritance().get_name().c_str());
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * ReflectionCoefficient(PyObject *self, PyObject *args)
{
  try {
    int handle;
    COMPLEX theta;
    double lambda;
    const char* n0;
    const char* nt;
    const char* type;
    typedef std::string string;
    
    if (!PyArg_ParseTuple(args, "iDdsss", &handle,&theta,&lambda,&n0,&nt,&type))
        return NULL;

    string stype(type);
    JonesMatrix result;
    if (stype==string("12")) result = mapStackModel[handle]->r12(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else if (stype==string("21")) result = mapStackModel[handle]->r21(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else if (stype==string("12i")) result = mapStackModel[handle]->r12i(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else if (stype==string("21i")) result = mapStackModel[handle]->r21i(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else throw SCATPY_exception("invalid type string: " + string(type));

    return PyJones(result);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * TransmissionCoefficient(PyObject *self, PyObject *args)
{
  try {
    int handle;
    COMPLEX theta;
    double lambda;
    const char* n0;
    const char* nt;
    const char* type;
    typedef std::string string;
    
    if (!PyArg_ParseTuple(args, "iDdsss", &handle,&theta,&lambda,&n0,&nt,&type))
        return NULL;

    string stype(type);
    JonesMatrix result;
    if (stype==string("12")) result = mapStackModel[handle]->t12(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else if (stype==string("21")) result = mapStackModel[handle]->t21(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else if (stype==string("12i")) result = mapStackModel[handle]->t12i(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else if (stype==string("21i")) result = mapStackModel[handle]->t21i(theta,lambda,dielectric_function(n0),dielectric_function(nt));
    else throw SCATPY_exception("invalid type string: " + string(type));
    
    return PyJones(result);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * GetOpticalConstant(PyObject *self, PyObject *args)
{
  try {
    double lambda;
    const char* material;
    
    if (!PyArg_ParseTuple(args, "ds", &lambda,&material)) return NULL;

    COMPLEX nk = dielectric_function(material).index(lambda);
    return  PYCOMPLEX(nk);
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * GetModelDictionary(PyObject *self, PyObject *args)
{
  try {
    const char* model;
    if (!PyArg_ParseTuple(args, "s", &model))
        return NULL;

    const Inheritance* modelinh = Model::inheritance.get_named_inheritance(model);
    if (!modelinh) SCATMECH_exception(std::string("There are no models under")+std::string(model));

    InheritanceList il = modelinh->get_progeny();


    PyObject *dict = PyDict_New();
    PyObject *dictsub;

    for (InheritanceList::iterator q=il.begin();q!=il.end();++q) {

       if (*q) {
           std::string name = (*q)->get_name();
           std::string desc = (*q)->get_description();
           const Inheritance* _parent = (*q)->get_parent();
           std::string parent = (_parent) ? _parent->get_name() : "Model";

           dictsub = Py_BuildValue("{s:s,s:s}", "description",desc.c_str(), "parent", parent.c_str());
           PyDict_SetItemString(dict,name.c_str(),dictsub);
           Py_DECREF(&dictsub);
       }
    }
    return dict;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * Lu_Chipman_Decomposition(PyObject *self, PyObject *args)
{
  try {

    MuellerMatrix m;
    if (!PyArg_ParseTuple(args, "dddddddddddddddd",
			  &m[0][0],&m[0][1],&m[0][2],&m[0][3],
			  &m[1][0],&m[1][1],&m[1][2],&m[1][3],
			  &m[2][0],&m[2][1],&m[2][2],&m[2][3],
			  &m[3][0],&m[3][1],&m[3][2],&m[3][3]
			  ))
        return NULL;

    MuellerMatrix depolarizer,retarder,diattenuator;
    m.Lu_Chipman_Decomposition(depolarizer, retarder, diattenuator);

    PyObject *result = PyTuple_New(3);
    PyTuple_SetItem(result, 0, PyMueller(depolarizer));
    PyTuple_SetItem(result, 1, PyMueller(retarder));
    PyTuple_SetItem(result, 2, PyMueller(diattenuator));
      
    return result;
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

extern "C" {

static PyMethodDef SCATPYMethods[] = {
     {"Get_Model",  Get_Model, METH_VARARGS, "Gets an instance of a Model."},
     {"Free_Model",  Free_Model, METH_VARARGS, "Frees an instance of a Model."},
     {"Get_BRDF_Model",  Get_BRDF_Model, METH_VARARGS, "Gets an instance of a BRDF_Model."},
     {"Free_BRDF_Model",  Free_BRDF_Model, METH_VARARGS, "Frees an instance of a BRDF_Model."},
     {"Get_Local_BRDF_Model",  Get_Local_BRDF_Model, METH_VARARGS, "Gets an instance of a Local_BRDF_Model."},
     {"Free_Local_BRDF_Model",  Free_Local_BRDF_Model, METH_VARARGS, "Frees an instance of a Local_BRDF_Model."},
     {"Get_StackModel",  Get_StackModel, METH_VARARGS, "Gets an instance of a StackModel."},
     {"Free_StackModel",  Free_StackModel, METH_VARARGS, "Frees an instance of a StackModel."},
     {"Get_Free_Space_Scatterer",  Get_Free_Space_Scatterer, METH_VARARGS, "Gets an instance of a Free_Space_Scatterer."},
     {"Free_Free_Space_Scatterer",  Free_Free_Space_Scatterer, METH_VARARGS, "Frees an instance of a Free_Space_Scatterer."},
     {"Get_RCW_Model",  Get_RCW_Model, METH_VARARGS, "Gets an instance of a RCW_Model."},
     {"Free_RCW_Model",  Free_RCW_Model, METH_VARARGS, "Frees an instance of a RCW_Model."},
     {"GetGratingEpsilon", GetGratingEpsilon, METH_VARARGS, "Gets epsilon at x,z"},
     {"GetGratingDefinition", GetGratingDefinition, METH_VARARGS, "Gets the grating definition in terms of transitions, dielectric constants, and thicknesses"},
     {"Get_CrossRCW_Model",  Get_CrossRCW_Model, METH_VARARGS, "Gets an instance of a CrossRCW_Model."},
     {"Free_CrossRCW_Model",  Free_CrossRCW_Model, METH_VARARGS, "Frees an instance of a CrossRCW_Model."},
     {"BRDF",  BRDF, METH_VARARGS, "Returns the Mueller matrix BRDF for a specified geometry"},
     {"VectoredBRDF", VectoredBRDF,  METH_VARARGS, "Returns the Mueller matrix BRDF for a specified geometry (expressed as vectors"},
     {"BRDFJones",  BRDFJones, METH_VARARGS, "Returns the Jones matrix BRDF for a specified geometry"},
     {"VectoredBRDFJones", VectoredBRDFJones,  METH_VARARGS, "Returns the Jones matrix BRDF for a specified geometry (expressed as vectors"},
     {"LocalDSC",  LocalDSC, METH_VARARGS, "Returns the Mueller matrix differential scattering cross section for a specified geometry"},
     {"LocalDSCJones",  LocalDSCJones, METH_VARARGS, "Returns the Jones matrix differential scattering cross section for a specified geometry"},
     {"FSSjones", FSSjones, METH_VARARGS, "Returns the jones scattering matrix for a specified geometry"},
     {"FSSext", FSSext, METH_VARARGS, "Returns the Mueller matrix extinction cross section"},
     {"RCWDiffractionEfficiency", RCWDiffractionEfficiency, METH_VARARGS, "Returns the Mueller matrix diffraction efficiency for a specific diffraction order (RCW_Model)"},
     {"RCWDiffractionAmplitude", RCWDiffractionAmplitude, METH_VARARGS, "Returns the Jones matrix diffraction amplitude for a specific diffraction order (RCW_Model)"},
     {"CrossRCWDiffractionEfficiency", CrossRCWDiffractionEfficiency, METH_VARARGS, "Returns the Mueller matrix diffraction efficiency for a specific diffraction order (CrossRCW_Model)"},
     {"RCWDirection", RCWDirection, METH_VARARGS, "Returns the direction [x,y,z] with unit length for a specific diffraction order (RCW_Model)"},
     {"CrossRCWDirection", CrossRCWDirection, METH_VARARGS, "Returns the direction [x,y,z] with unit length for a specific diffraction order (CrossRCW_Model)"},
     {"GetParameterDictionary", GetParameterDictionary, METH_VARARGS, "Gets a list of all the parameters for the Model"},
     {"SetParameter",  SetParameter, METH_VARARGS, "Sets a parameter for the Model"},
     {"GetParameter",  GetParameter, METH_VARARGS, "Gets a parameter for the Model"},
     {"PrintParameters",  PrintParameters, METH_VARARGS, "Gets a list of all the parameters for the Model"},
     {"AskParameters",  AskParameters, METH_VARARGS, "Queries the user for model parameters"}, 
     {"GetModelDictionary", GetModelDictionary, METH_VARARGS, "Gets a dictionary of all models"},
     {"GetModelName", GetModelName, METH_VARARGS, "Gets the name of a model"},
     {"ReflectionCoefficient", ReflectionCoefficient, METH_VARARGS, "Gets the reflection coeffient"},
     {"TransmissionCoefficient", TransmissionCoefficient, METH_VARARGS, "Gets the transmission coefficient"},     
     {"GetOpticalConstant", GetOpticalConstant, METH_VARARGS, "Gets the optical constant from a file"},
     {"Lu_Chipman_Decomposition",Lu_Chipman_Decomposition, METH_VARARGS, "Performs a Lu-Chipman decomposition of a Mueller matrix"},
     
    {NULL, NULL, 0, NULL}        /* Sentinel */  
};

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
} CustomObject;

static PyTypeObject SCATPYType = {
    PyVarObject_HEAD_INIT(NULL, 0)
};

static PyModuleDef SCATPYmodule;
  

PyMODINIT_FUNC
PyInit_SCATPY(void)
{
    PyObject *m;

    Register((Model*)0);
    //Register_Special_Models();

    ///    if (PyType_Ready(&SCATPYType) < 0)
    ///    return NULL;

    SCATPYmodule.m_base.m_init = NULL;
    SCATPYmodule.m_base.m_index = 0;
    SCATPYmodule.m_base.m_copy = NULL;
    
    SCATPYmodule.m_name = "SCATPY";
    SCATPYmodule.m_doc = "Interface to SCATMECH C++ codes.";
    SCATPYmodule.m_size = -1;
    SCATPYmodule.m_methods = SCATPYMethods;

    
    m = PyModule_Create(&SCATPYmodule);
    if (m == NULL)
        return NULL;

    ///Py_INCREF(&SCATPYType);
    ///PyModule_AddObject(m, "SCATPY", (PyObject *) &SCATPYType);

    return m;
}
}

int main(int argv,char** argc)
{
  Register((Model*)0);
}

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

using namespace SCATMECH;

std::map<int,Model*> mapModel;
std::map<int,BRDF_Model_Ptr> mapBRDF_Model; 
std::map<int,Local_BRDF_Model_Ptr> mapLocal_BRDF_Model;
std::map<int,Model_Ptr<RCW_Model> > mapRCW_Model;
std::map<int,Model_Ptr<CrossRCW_Model> > mapCrossRCW_Model;
std::map<int,Free_Space_Scatterer_Ptr> mapFSS_Model;
std::map<int,StackModel_Ptr> mapStackModel;
std::map<int,Model_Ptr<Model> > mapLone_Model;
	       
int modelct = 0;

class SCATPY_exception: public std::exception
{
public:
  SCATPY_exception(const std::string& _msg) : msg(_msg) {}
#ifndef _MSC_VER
  virtual ~SCATPY_exception() _GLIBCXX_USE_NOEXCEPT {}
#endif
  virtual const char *what() const throw() {
    return msg.c_str();
  }
private:
  std::string msg;
};

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
    const char *modelname;

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

    return Py_BuildValue("[[d,d,d,d],[d,d,d,d],[d,d,d,d],[d,d,d,d]]",
		  mueller[0][0],mueller[0][1],mueller[0][2],mueller[0][3],
		  mueller[1][0],mueller[1][1],mueller[1][2],mueller[1][3],
		  mueller[2][0],mueller[2][1],mueller[2][2],mueller[2][3],
		  mueller[3][0],mueller[3][1],mueller[3][2],mueller[3][3]);
		  
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

    return Py_BuildValue("[[d,d,d,d],[d,d,d,d],[d,d,d,d],[d,d,d,d]]",
		  mueller[0][0],mueller[0][1],mueller[0][2],mueller[0][3],
		  mueller[1][0],mueller[1][1],mueller[1][2],mueller[1][3],
		  mueller[2][0],mueller[2][1],mueller[2][2],mueller[2][3],
		  mueller[3][0],mueller[3][1],mueller[3][2],mueller[3][3]);
		  
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

    return Py_BuildValue("[[D,D],[D,D]]",
			 jones.SS(),jones.SP(),jones.PS(),jones.PP());
		  
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

    return Py_BuildValue("[[d,d,d,d],[d,d,d,d],[d,d,d,d],[d,d,d,d]]",
		  mueller[0][0],mueller[0][1],mueller[0][2],mueller[0][3],
		  mueller[1][0],mueller[1][1],mueller[1][2],mueller[1][3],
		  mueller[2][0],mueller[2][1],mueller[2][2],mueller[2][3],
		  mueller[3][0],mueller[3][1],mueller[3][2],mueller[3][3]);
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

    return Py_BuildValue("[[d,d,d,d],[d,d,d,d],[d,d,d,d],[d,d,d,d]]",
		  mueller[0][0],mueller[0][1],mueller[0][2],mueller[0][3],
		  mueller[1][0],mueller[1][1],mueller[1][2],mueller[1][3],
		  mueller[2][0],mueller[2][1],mueller[2][2],mueller[2][3],
		  mueller[3][0],mueller[3][1],mueller[3][2],mueller[3][3]);
		  
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

    return Py_BuildValue("[[d,d,d,d],[d,d,d,d],[d,d,d,d],[d,d,d,d]]",
		  mueller[0][0],mueller[0][1],mueller[0][2],mueller[0][3],
		  mueller[1][0],mueller[1][1],mueller[1][2],mueller[1][3],
		  mueller[2][0],mueller[2][1],mueller[2][2],mueller[2][3],
		  mueller[3][0],mueller[3][1],mueller[3][2],mueller[3][3]);
		  
  }
  catch (std::exception& e) {
    PyErr_SetString(PyExc_Exception,(format("At %s, ",__func__) + e.what()).c_str());
    return NULL;
  }    
}

static PyObject * RCWDirection(PyObject *self, PyObject *args)
{
  try {
    int handle,i,j;
    
    if (!PyArg_ParseTuple(args, "ii", &handle, &i))
        return NULL;

    Vector direction = mapRCW_Model[handle]->GetDirection(i);

    return Py_BuildValue("[d,d,d]", direction.x,direction.y,direction.z);
		  
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
    if (grating->get_lambda()!=lambda) grating->set_lambda(lambda);
    if (level==-1) eps = grating->get_medium_i().epsilon(lambda);
    else if (level==-2) eps = grating->get_medium_t().epsilon(lambda);
    else eps = grating->eps(x,level,direction);

    return PyComplex_FromDoubles(real(eps),imag(eps));
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

    return Py_BuildValue("[d,d,d]", direction.x,direction.y,direction.z);
		  
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


    return Py_BuildValue("[[D,D],[D,D]]",
			 result[1],result[2],result[3],result[0]);
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
    
    return Py_BuildValue("[[D,D],[D,D]]",
			 result[1],result[2],result[3],result[0]);
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
    return Py_BuildValue("D",nk);
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
     {"Get_CrossRCW_Model",  Get_CrossRCW_Model, METH_VARARGS, "Gets an instance of a CrossRCW_Model."},
     {"Free_CrossRCW_Model",  Free_CrossRCW_Model, METH_VARARGS, "Frees an instance of a CrossRCW_Model."},
     {"BRDF",  BRDF, METH_VARARGS, "Returns the Mueller matrix BRDF for a specified geometry"},
     {"LocalDSC",  LocalDSC, METH_VARARGS, "Returns the Mueller matrix differential scattering cross section for a specified geometry"},
     {"FSSjones", FSSjones, METH_VARARGS, "Returns the jones scattering matrix for a specified geometry"},
     {"FSSext", FSSext, METH_VARARGS, "Returns the Mueller matrix extinction cross section"},
     {"RCWDiffractionEfficiency", RCWDiffractionEfficiency, METH_VARARGS, "Returns the Mueller matrix diffraction efficiency for a specific diffraction order (RCW_Model)"},
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

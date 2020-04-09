//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: local.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include <complex>
#include "scatmech.h"
#include "brdf.h"
#include "local.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {


    //
    // The following function returns the jones matrix for scattering using a
    // specific coordinate system, presumably, but not necessarily different,
    // than the {s,p} system.
    //
    // This function was added in Version 3 (TAG: 18 JAN 2001)
    //
    JonesMatrix
    Local_BRDF_Model::
    JonesDSC(double thetai,double thetas,double phis,double rotation, Coordinate_System cs)
    {
        set_geometry(thetai,thetas,phis,rotation);
        JonesMatrix j = jonesDSC();
        convert(j,cs);
        return j;
    }

    //
    // The following function returns the mueller matrix for scattering using a
    // specific coordinate system, presumably, but not necessarily different,
    // than the {s,p} system.
    //
    // This function was added in Version 3 (TAG: 18 JAN 2001)
    //
    MuellerMatrix
    Local_BRDF_Model::
    MuellerDSC(double thetai,double thetas,double phis,double rotation,Coordinate_System cs)
    {
        set_geometry(thetai,thetas,phis,rotation);
        MuellerMatrix m = muellerDSC();
        convert(m,cs);
        return m;
    }

    void
    Local_BRDF_Model::
    setup()
    {
        BRDF_Model::setup();
    }

    JonesMatrix
    Local_BRDF_Model::
    JonesDSC(const Vector& source, const Vector& viewer, const Vector& normal, const Vector& xaxis,
             Coordinate_System cs)
    {
        if (set_geometry(source,viewer,normal,xaxis)) {
            return JonesDSC(thetai,thetas,phis,rotation,cs);
        } else {
            return JonesZero();
        }
    }

    MuellerMatrix
    Local_BRDF_Model::
    MuellerDSC(const Vector& source, const Vector& viewer, const Vector& normal, const Vector& xaxis,
               Coordinate_System cs)
    {
        if (set_geometry(source,viewer,normal,xaxis)) {
            return MuellerDSC(thetai,thetas,phis,rotation,cs);
        } else {
            return MuellerZero();
        }
    }

    DEFINE_VIRTUAL_MODEL(Local_BRDF_Model,BRDF_Model,
                         "All models for discrete defects.");

    DEFINE_PARAMETER(Local_BRDF_Model,double,density,"Density [um^-2]","1",0xFF);


} // namespace SCATMECH


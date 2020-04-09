//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: lambert.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "lambert.h"

using namespace std;


namespace SCATMECH {


    MuellerMatrix
    Lambertian_BRDF_Model::
    mueller()
    {
        SETUP();

        return MuellerDepolarizer(reflectance->Get_Reflectance(lambda)/pi);
    }

    DEFINE_MODEL(Lambertian_BRDF_Model,BRDF_Model,
                 "A totally diffuse and depolarizing BRDF.");

    DEFINE_PTRPARAMETER(Lambertian_BRDF_Model,Reflectance_Ptr,reflectance,"Reflectance","Table_Reflectance",0xFF);


} // namespace SCATMECH



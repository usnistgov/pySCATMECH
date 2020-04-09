//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reg_lambert.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "lambert.h"
#include "diffuse.h"


namespace SCATMECH {

    void Register(const Lambertian_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Lambertian_BRDF_Model);
            Register_Model(Diffuse_Subsurface_BRDF_Model);
        }
    }


} // namespace SCATMECH {


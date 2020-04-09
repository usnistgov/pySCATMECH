//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reg_local.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "sphprt.h"
#include "sphdfct.h"
#include "onelayer.h"
#include "bobvlieg.h"
#include "subsphere.h"
#include "axipart.h"
#include "raystack.h"
#include "subbobvlieg.h"

namespace SCATMECH {


    void Register(const Local_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Local_BRDF_Model);
            Register_Model(Rayleigh_Defect_BRDF_Model);
            Register_Model(Rayleigh_Stack_BRDF_Model);
            Register_Model(OneLayer_BRDF_Model);
            Register_Model(Double_Interaction_BRDF_Model);
            Register_Model(Subsurface_Particle_BRDF_Model);
            Register_Model(Bobbert_Vlieger_BRDF_Model);
            Register_Model(Axisymmetric_Particle_BRDF_Model);
            Register_Model(Subsurface_Bobbert_Vlieger_BRDF_Model);
            Register_Model(Subsurface_Axisymmetric_Particle_BRDF_Model);
            Register((Axisymmetric_Shape*)0);
        }
    }


} // namespace SCATMECH {


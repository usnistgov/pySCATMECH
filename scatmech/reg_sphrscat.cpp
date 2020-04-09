//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reg_brdf.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "rayscat.h"
#include "raygscat.h"
#include "miescat.h"
#include "coatedmie.h"
#include "axifree.h"
#include "nsphere.h"
#include "oasphere.h"

namespace SCATMECH {


    void Register(const Free_Space_Scatterer* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;
            Register_Model(Free_Space_Scatterer);
            Register_Model(TMatrix_Axisymmetric_Scatterer);
            Register_Model(Optically_Active_Sphere_Scatterer);
            Register((SphericalScatterer*)NULL);
        }
    }

    void Register(const SphericalScatterer* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(SphericalScatterer);
            Register_Model(RayleighScatterer);
            Register_Model(RayleighGansSphereScatterer);
            Register_Model(MieScatterer);
            Register_Model(CoatedMieScatterer);
            Register_Model(MultilayerCoatedMieScatterer);
        }
    }



} // namespace SCATMECH {


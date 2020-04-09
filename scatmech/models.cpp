//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: models.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "brdf.h"
#include "sphrscat.h"
#include "facet.h"
#include "rough.h"
#include "torrspar.h"
#include "rcw.h"
#include "crossrcw.h"
#include "reflectance.h"
#include "filmtran.h"

using namespace std;

namespace SCATMECH {

    void Register(const Model* x)
    {
        static int Models_Registered = 0;
        if (!Models_Registered) {
            Models_Registered=1;

            Register((BRDF_Model*)NULL);
            Register((Free_Space_Scatterer*)NULL);
            Register((Slope_Distribution_Function*)NULL);
            Register((PSD_Function*)NULL);
            Register((Shadow_Function*)NULL);
            Register((Reflectance*)NULL);
            Register((Grating*)NULL);
            Register((CrossGrating*)NULL);
            Register_Model(RCW_Model);
            Register_Model(CrossRCW_Model);
			Register((StackModel*)NULL);
        }
    }

} // namespace SCATMECH


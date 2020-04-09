//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: roughnes.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_ROUGHNESS_H
#define SCATMECH_ROUGHNESS_H

#include "scatmech.h"
#include "urough.h"
#include "filmtran.h"


namespace SCATMECH {



    //
    // Roughness_Stack_BRDF_Model is a Roughness_BRDF_Model with a dielectric_stack which
    // calculates the scatter from a single rough interface in a dielectric stack.
    //
    class Roughness_Stack_BRDF_Model :
        public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();

            // The rough layer index...(0 = most buried)...
            DECLARE_PARAMETER(int,this_layer);

            // The films
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:

            // The Jones matrix for scattering...
            virtual JonesMatrix jones();

    };

    void Register(const Roughness_Stack_BRDF_Model* x);


} // namespace SCATMECH



#endif


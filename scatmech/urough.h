//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: urough.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_UROUGH_H
#define SCATMECH_UROUGH_H

#include "scatmech.h"
#include "rough.h"


namespace SCATMECH {


    //
    // Microroughness_BRDF_Model is a Roughness_BRDF_Model that corresponds to
    // Rayleigh-Rice Theory (First-order vector perturbation theory.
    //
    class Microroughness_BRDF_Model : public Roughness_BRDF_Model
    {
        public:
            DECLARE_MODEL();
        protected:
            // Routines that return the Mueller and Jones Matrices...
            JonesMatrix jones();
    };


} // namespace SCATMECH


#endif


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crough.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_CROUGH_H
#define SCATMECH_CROUGH_H

#include "rough.h"
#include "allrough.h"


namespace SCATMECH {

    //
    // Correlated_Roughness_BRDF_Model is a Roughness_BRDF_Model that calculates
    // the scattering from a single thin film with correlated roughness.
    //
    class Correlated_Roughness_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            // The thickness of the film (in units of length)
            DECLARE_PARAMETER(double,thickness);
            // Dielectric function of the film
            DECLARE_PARAMETER(dielectric_function,film);

        protected:

            void setup();

            // Routines that return the Mueller and Jones Matrices...
            JonesMatrix jones();

        private:
            Correlated_Roughness_Stack_BRDF_Model model;

    };

} // namespace SCATMECH


#endif


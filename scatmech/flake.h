//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: flake.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_FLAKE_H
#define SCATMECH_FLAKE_H

#include "facet.h"


namespace SCATMECH {

    class Subsurface_Facet_BRDF_Model:
        public Facet_BRDF_Model
    {
        public:
            // Functions which overload the Facet_BRDF_Model definitions for
            // local_angle and local_slope...
            virtual double local_angle(double thetai,double thetas,double phis);
            virtual double local_slope(double thetai,double thetas,double phis);

            DECLARE_MODEL();

            // The optical properties of the overcoat...
            DECLARE_PARAMETER(dielectric_function,overcoat);

            // Any films which exist above the overcoat...
            DECLARE_PARAMETER(StackModel_Ptr,overcoat_films);

        protected:

            // Routine to perform one time calculations...
            void setup();

            // Jones matrix for scattering...
            JonesMatrix jones();
    };


} // namespace SCATMECH



#endif

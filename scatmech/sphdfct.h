//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sphdfct.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_SPHDFCT_H
#define SCATMECH_SPHDFCT_H

#include "scatmech.h"
#include "filmtran.h"
#include "local.h"


namespace SCATMECH {



    //
    // Spherical_Defect_BRDF_Model is a BRDF_Model which handles a
    // spherically-symetric scatterer beneath a surface.
    //
    class Rayleigh_Defect_BRDF_Model :
        public Local_BRDF_Model
    {
        public:

            DECLARE_MODEL();

            // The radius of the sphere in units of length...
            DECLARE_PARAMETER(double,radius);

            // The distance beneath the surface in units of length...
            DECLARE_PARAMETER(double,distance);

            // The dielectric function of the particle
            DECLARE_PARAMETER(dielectric_function,defect);

            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:


            // Routine to carry out one-time calculations...
            void setup();

            // The spherical scattering model (Mie, Rayleigh-Gans, or Rayleigh)...
            // SphericalScatterer *scatterer;

            // The distance beneath the surface of the center of the sphere
            // (*2*pi/lambda)...
            double kd;

            // The radius of the sphere (*2*pi/lambda)...
            double kr;

        protected:
            // The Jones Matrix for scattering...
            virtual JonesMatrix jonesDSC();

    };


} // namespace SCATMECH



#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: onelayer.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_ONELAYER_H
#define SCATMECH_ONELAYER_H

#include "scatmech.h"
#include "local.h"


namespace SCATMECH {



    //
    // OneLayer_BRDF_Model is a BRDF_Model that calculates the scatter
    // from a very small spherical defect (in the Rayleigh approximation)
    // in a single dielectric film.
    //
    class OneLayer_BRDF_Model : public Local_BRDF_Model
    {
        public:

            DECLARE_MODEL();

            // The radius of the defect...
            DECLARE_PARAMETER(double,radius);

            // The optical properties of the defect...
            DECLARE_PARAMETER(dielectric_function,defect);

            // The optical properties of the film...
            DECLARE_PARAMETER(dielectric_function,film);

            // The thickness of the film...
            DECLARE_PARAMETER(double,tau);

            // The depth of the defect in the film
            // (depth is positive and measured from exposed interface)...
            DECLARE_PARAMETER(double,depth);

        protected:
            // Routine which returns the Jones matrix for scattering...
            JonesMatrix jonesDSC();

    };


} // namespace SCATMECH


#endif


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rayinst.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_INSTSIG_H
#define SCATMECH_INSTSIG_H

#include "scatmech.h"
#include "instrument.h"
#include "filmtran.h"


namespace SCATMECH {


    //
    // Rayleigh_Instrument_BRDF_Model is an Instrument_BRDF_Model for the
    // Rayleigh scatter in the air surrounding a sample.
    //
    class Rayleigh_Instrument_BRDF_Model : public Instrument_BRDF_Model
    {
        public:

            DECLARE_MODEL();

            // The field of view of the detector [um]...
            DECLARE_PARAMETER(double,field_of_view);

            // The index of refraction of the air minus 1 ...
            DECLARE_PARAMETER(dielectric_function,air);

            // The number density of the air [um^-3]
            DECLARE_PARAMETER(double,number_density);

            // Any dielectric films on the substrate...
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:
            // The Mueller matrix BRDF for scattering...
            MuellerMatrix mueller();

        private:
    };


} // namespace SCATMECH



#endif

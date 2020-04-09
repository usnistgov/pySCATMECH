//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: twoface.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_TWOFACE_H
#define SCATMECH_TWOFACE_H

#include "scatmech.h"
#include "rough.h"
#include "roughnes.h"

namespace SCATMECH {



    //
    // Two_Face_BRDF_Model is a Roughness_BRDF_Model that calculates the scatter from
    // one of the two interfaces of a dielectric film.
    //
    class Two_Face_BRDF_Model: public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();

            // The interface number... (1 = buried, 2 = exposed)
            DECLARE_PARAMETER(int,face);

            // The film thickness...
            DECLARE_PARAMETER(double,thickness);

            // The dielectric constant of the film...
            DECLARE_PARAMETER(dielectric_function,film);

        protected:

            void setup();

            // Routines which returns the scattering Jones matrix...
            virtual JonesMatrix jones();

        private:
            Roughness_Stack_BRDF_Model model;

    };


} // namespace SCATMECH


#endif

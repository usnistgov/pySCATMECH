//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rayscat.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_RAYSCAT_H
#define SCATMECH_RAYSCAT_H

#include "scatmech.h"
#include "sphrscat.h"


namespace SCATMECH {



    //
    // RayleighScatterer is a SphericalScatterer in the Rayleigh approximation.
    //
    class RayleighScatterer: public SphericalScatterer {
        public:

            // The scattering matrix elements are as follows...
            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            // The scattering, extinction, and backscattering cross sections...
            double Csca();
            double Cext();
            double Cback();

        public:
            DECLARE_MODEL();

        private:
            // The polarizability is its only parameters...
            COMPLEX polarizability;

            virtual void setup();
    };


} // namespace SCATMECH


#endif

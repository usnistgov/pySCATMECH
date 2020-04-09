//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: raygscat.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_RAYGSCAT_H
#define SCATMECH_RAYGSCAT_H

#include "scatmech.h"
#include "sphrscat.h"


namespace SCATMECH {



    class RayleighGansSphereScatterer: public SphericalScatterer {
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

            virtual void setup();
    };


} // namespace SCATMECH


#endif

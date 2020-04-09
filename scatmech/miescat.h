//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: miescat.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_MIESCAT_H
#define SCATMECH_MIESCAT_H

#include "sphrscat.h"
#include "dielfunc.h"
#include <vector>


namespace SCATMECH {



    //
    // MieScatterer is a SphericalScatterer which returns the Mie
    // solution for the free space scattering of a homogeneous sphere.
    //
    class MieScatterer: public SphericalScatterer {
        public:

            MieScatterer();

            // The scattering matrix elements...
            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();

        private:

            std::vector<COMPLEX> A,B;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;

            virtual void setup();
    };


} // namespace SCATMECH


#endif

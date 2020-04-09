//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: coated.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_COATED_H
#define SCATMECH_COATED_H

#include "sphrscat.h"
#include "dielfunc.h"
#include <vector>


namespace SCATMECH {


    //
    // CoatedMieScatterer is a SphericalScatterer which returns the Mie
    // solution for the free space scattering of a homogeneous sphere with
    // a single coating.
    //
    class CoatedMieScatterer: public SphericalScatterer {
        public:
            CoatedMieScatterer();

            // The scattering matrix elements...
            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,thickness);
            DECLARE_PARAMETER(dielectric_function,coating);
            DECLARE_PARAMETER(int,nmax);

        private:

            std::vector<COMPLEX> A,B;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;

            void setup();
    };


} // namespace SCATMECH


#endif

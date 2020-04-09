//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: nsphere.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_NSPHERE_H
#define SCATMECH_NSPHERE_H

#include "sphrscat.h"
#include "dielfunc.h"
#include "filmtran.h"
#include <vector>

namespace SCATMECH {
    //
    // MultilayerCoatedMieScatterer is a SphericalScatterer that returns the
    // solution for the free space scattering of a sphere with
    // any number of coatings.
    //
    class MultilayerCoatedMieScatterer: public SphericalScatterer {
        public:
            MultilayerCoatedMieScatterer();

            // The scattering matrix elements...
            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        private:

            std::vector<COMPLEX> A,B;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;

            void setup();
    };
} // namespace SCATMECH


#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: finiteaperture.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_FINITEAPERTURE_H
#define SCATMECH_FINITEAPERTURE_H

#include "scatmech.h"
#include "instrument.h"
#include "vector3d.h"


namespace SCATMECH {


    class Finite_Aperture_Instrument_BRDF_Model : public Instrument_BRDF_Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(int,integralmode)
            DECLARE_PARAMETER(BRDF_Model_Ptr,model);

        public:
            Finite_Aperture_Instrument_BRDF_Model();
        protected:

            virtual MuellerMatrix mueller();
            static Vector four_angles(double theta,double phi,double alpha,double beta);
    };


}


#endif



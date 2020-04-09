//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: two_source.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef TWO_SOURCE_H
#define TWO_SOURCE_H

#include "brdf.h"

namespace SCATMECH {

    class Two_Source_BRDF_Model : public BRDF_Model
    {
        public:
            MuellerMatrix mueller();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,factor1);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source1);
            DECLARE_PARAMETER(double,factor2);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source2);
            DECLARE_PARAMETER(double,correlation);
    };

    class Three_Source_BRDF_Model : public BRDF_Model
    {
        public:
            MuellerMatrix mueller();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,factor1);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source1);
            DECLARE_PARAMETER(double,factor2);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source2);
            DECLARE_PARAMETER(double,factor3);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source3);
    };

    class Four_Source_BRDF_Model : public BRDF_Model
    {
        public:
            MuellerMatrix mueller();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,factor1);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source1);
            DECLARE_PARAMETER(double,factor2);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source2);
            DECLARE_PARAMETER(double,factor3);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source3);
            DECLARE_PARAMETER(double,factor4);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source4);
    };

}

#endif


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: transmit.h
//**
//** Thomas A. Germer
//** Sensor Science Division,
//** National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//**
//** Version 7.00 (January 2015)
//**
//******************************************************************************

#ifndef SCATMECH_TRANSMIT_H

#include "brdf.h"
#include "filmtran.h"

namespace SCATMECH {

    class Transmit_BRDF_Model : public BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(BRDF_Model_Ptr,model);
            DECLARE_PARAMETER(StackModel_Ptr,films);

        protected:
            virtual MuellerMatrix mueller();
    };

}

#endif


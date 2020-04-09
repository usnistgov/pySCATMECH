//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: instrument.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_INSTRUMENT_H
#define SCATMECH_INSTRUMENT_H

#include "brdf.h"


namespace SCATMECH {



    class Instrument_BRDF_Model : public BRDF_Model
    {
        public:
            DECLARE_MODEL();
    };

    void Register(const Instrument_BRDF_Model* x);


}


#endif

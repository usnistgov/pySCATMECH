//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reg_instrument.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "instrument.h"
#include "rayinst.h"
#include "finiteaperture.h"
#include "focussedbeam.h"


namespace SCATMECH {


    void Register(const Instrument_BRDF_Model* x)
    {
        static bool regd=false;
        if (!regd) {
            regd=true;

            Register_Model(Instrument_BRDF_Model);
            Register_Model(Rayleigh_Instrument_BRDF_Model);
            Register_Model(Finite_Aperture_Instrument_BRDF_Model);
            Register_Model(Focussed_Beam_Instrument_BRDF_Model);
        }
    }


}


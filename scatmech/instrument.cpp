//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: instrument.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "instrument.h"


namespace SCATMECH {


    DEFINE_VIRTUAL_MODEL(Instrument_BRDF_Model,BRDF_Model,
                         "A virtual class for effective BRDF models for instrument functions");


}


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: gaussianbeam.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** This software was developed at the National Institute of Standards and Technology 
//** by employees of the Federal Government in the course of their official duties. 
//** Pursuant to Title 17 Section 105 of the United States Code this software is not 
//** subject to copyright protection and is in the public domain. SCATMECH is an experimental 
//** system. NIST assumes no responsibility whatsoever for its use by other parties, and 
//** makes no guarantees, expressed or implied, about its quality, reliability, or any other 
//** characteristic. We would appreciate acknowledgment if the software is used. This software 
//** can be redistributed and/or modified freely provided that any derivative works bear some 
//** notice that they are derived from it, and any modified versions bear some notice 
//** that they have been modified.
//******************************************************************************

#ifndef SCATMECH_GAUSSIANBEAM_H
#define SCATMECH_GAUSSIANBEAM_H

#include "scatmech.h"
#include "instrument.h"
#include "inherit.h"
#include "vector3d.h"
#include "focussedbeam.h"
#include "local.h"

namespace SCATMECH {

    class Gaussian_Beam_Instrument_BRDF_Model : public Instrument_BRDF_Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(int,integralmode)
			DECLARE_PARAMETER(double, NA);
			DECLARE_PARAMETER(double, xNA);
			DECLARE_PARAMETER(double, focus_x);
			DECLARE_PARAMETER(double, focus_y);
			DECLARE_PARAMETER(double, focus_z);
			DECLARE_PARAMETER(double, offset_theta);
			DECLARE_PARAMETER(double, offset_phi);
			DECLARE_PARAMETER(Local_BRDF_Model_Ptr,model);

        public:
			Gaussian_Beam_Instrument_BRDF_Model();
        protected:
            virtual JonesMatrix jones();
    };

} // namespace SCATMECH {


#endif



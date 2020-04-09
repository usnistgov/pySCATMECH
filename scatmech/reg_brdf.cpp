//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reg_brdf.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "rough.h"
#include "facet.h"
#include "lambert.h"
#include "local.h"
#include "instrument.h"
#include "firstdiffuse.h"
#include "two_source.h"
#include "transmit.h"
#include "rcw.h"
#include "crossrcw.h"
#include "zernikeexpansion.h"
#include "polydisperse.h"

namespace SCATMECH {


	void Register(const BRDF_Model* x)
	{
		static bool Models_Registered = false;
		if (!Models_Registered) {
			Models_Registered = true;

			Register_Model(BRDF_Model);
			Register((Roughness_BRDF_Model*)0);
			Register((Facet_BRDF_Model*)0);
			Register((Lambertian_BRDF_Model*)0);
			Register((Local_BRDF_Model*)0);
			Register((Instrument_BRDF_Model*)0);
			Register((First_Diffuse_BRDF_Model*)0);
			Register_Model(Two_Source_BRDF_Model);
			Register_Model(Three_Source_BRDF_Model);
			Register_Model(Four_Source_BRDF_Model);
			Register_Model(Transmit_BRDF_Model);
			Register_Model(RCW_BRDF_Model);
			Register_Model(CrossRCW_BRDF_Model);
			Register_Model(ZernikeExpansion_BRDF_Model);
			Register_Model(Polydisperse_Sphere_BRDF_Model);

		}
	}



} // namespace SCATMECH {


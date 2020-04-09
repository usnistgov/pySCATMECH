//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: polydisperse.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_POLYDISPERSE_H
#define SCATMECH_POLYDISPERSE_H

#include "bobvlieg.h"
#include "firstdiffuse.h"
#include <vector>
#include "sphrscat.h"
#include "oasphere.h"
#include "sphprt.h"

namespace SCATMECH {
	
	///
	/// @brief Model for polydisperse spheres on a surface.
	///
	/// Polydisperse_Sphere_BRDF_Model uses the Double_Interaction_BRDF_Model with MieScatterer and
	/// SurfaceParticleSizeDistribution to estimate the scatter from a surface contaminated with particles.
	///
	class Polydisperse_Sphere_BRDF_Model : public BRDF_Model {
	public:             		
		DECLARE_MODEL();

		/// The optical constants of the particle
		DECLARE_PARAMETER(dielectric_function, particle);

		/// Any dielectric layers on the surface
		DECLARE_PARAMETER(StackModel_Ptr, stack);

		/// The distribution of particles diameters
		DECLARE_PARAMETER(SurfaceParticleSizeDistribution_Ptr,distribution);

		/// The starting diameter for integration 
		DECLARE_PARAMETER(double, Dstart);

		/// The ending diameter for integration
		DECLARE_PARAMETER(double, Dend);

		/// The fractional step for integration
		DECLARE_PARAMETER(double, Dstep);

		/// If non-zero, the fractional area coverage of spheres.
		/// If zero, use distribution to determine coverage.
		DECLARE_PARAMETER(double, fractional_coverage);

		/// A parameter that requires minimum absorption for the particle
		/// to eliminate rainbow effects. The additional absorption is diameter dependent, 
		/// so that there is 1/e absorption over the diameter of the particle when
		/// antirainbow = 1. Recommended value = 10 to completely eliminate rainbows.
		DECLARE_PARAMETER(double, antirainbow);
		
	protected:
		MuellerMatrix mueller();
		
	private:
		void setup();

		typedef Model_Ptr<Double_Interaction_BRDF_Model> Double_Interaction_BRDF_Model_Ptr;
		typedef std::list<Double_Interaction_BRDF_Model_Ptr> ModelList;
		typedef ModelList::iterator ModelListIterator;

		ModelList submodels;

		double distributioncoverage;
	};


} // namespace SCATMECH;

#endif

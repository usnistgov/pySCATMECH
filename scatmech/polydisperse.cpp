//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: polydisperse.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "polydisperse.h"
#include "miescat.h"

using namespace std;

namespace SCATMECH {

void
Polydisperse_Sphere_BRDF_Model::
setup()
{
    BRDF_Model::setup();

	// Clear the list of sub-models...
	submodels.clear();

	// Get the fractional area coverage...
	double sumdistribution = 0;
	double sumdistributionA = 0;

	// Then start creating a list of models, each with a specific diameter...
	for (double diameter = Dstart; diameter <= Dend; diameter *= 1. + Dstep) {

		// Push one back...
		submodels.push_back(new Double_Interaction_BRDF_Model);

		// 
		submodels.back()->set_lambda(lambda);
		submodels.back()->set_substrate(substrate);
		submodels.back()->set_type(type);

		// The diameter steps are logarithmic, so the density of spheres at this 
		// diameter is given by...
		double density = distribution->surfacedensity(diameter) * diameter * Dstep;
		submodels.back()->set_density(density);

		sumdistribution += density;
		sumdistributionA += density*pi*sqr(diameter / 2);

		double radius = diameter / 2;
		submodels.back()->set_stack(stack);
		submodels.back()->set_distance(radius);

		MieScatterer *mie = new MieScatterer;
		mie->set_medium(optical_constant(1, 0));
		mie->set_lambda(lambda);
		mie->set_radius(radius);

		// Antirainbow calculation...
		// First, get the optical constants of the particle...
		optical_constant nsphere = particle.index(lambda);

		// Then, calculate the minimum absorption...
		double absmin = antirainbow / diameter;

		// and what it means in terms of coefficient (k)...
		double kmin = absmin / 4. / pi*lambda;

		// If the absorption coefficient is less than this value, set it...
		if (nsphere.k < kmin) nsphere.k = kmin;
		mie->set_sphere(nsphere);

		submodels.back()->set_scatterer(mie);

		submodels.back()->set_alpha(0.);
		submodels.back()->set_beta(0.);
	}

	distributioncoverage = sumdistributionA;
}

MuellerMatrix 
Polydisperse_Sphere_BRDF_Model::
mueller()
{
	SETUP();

	// Zero the result...
	MuellerMatrix result = MuellerZero();

	// Sum up all the models...
	for (ModelListIterator m = submodels.begin(); m != submodels.end(); ++m) {
		result += (*m)->Mueller(thetai, thetas, phis, rotation, model_cs);
	}

	if (fractional_coverage != 0) result *= fractional_coverage / distributioncoverage;

	return result;
}

DEFINE_MODEL(Polydisperse_Sphere_BRDF_Model,BRDF_Model,"Double interaction theory with a Mie scatterer and polydispersity.");
DEFINE_PTRPARAMETER(Polydisperse_Sphere_BRDF_Model,SurfaceParticleSizeDistribution_Ptr,distribution,"Size distribution","Regular_SurfaceParticleSizeDistribution",0xFF);
DEFINE_PTRPARAMETER(Polydisperse_Sphere_BRDF_Model,StackModel_Ptr,stack,"Films on substrate","No_StackModel",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,dielectric_function,particle,"Optical properties of the particle","(1.5,0.0)",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,double,Dstart,"Starting diameter [um]","1",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,double,Dend,"Ending diameter [um]","1000",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,double,Dstep,"Fractional diameter step","0.01",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model, double, fractional_coverage, "Fractional area coverage (0 = use distribution's coverage)", "0", 0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model, double, antirainbow, "Anti-rainbow factor", "0", 0xFF);

}
//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: firstdiffuse.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "firstdiffuse.h"

#include "askuser.h"
#include "matrix3d.h"

using namespace std;


namespace SCATMECH {


	///
	/// First_Diffuse_BRDF_Model calculates the net singly scattered light
	/// in a multiply scattering medium with a given polarized phase function.  
	/// This model should generally be combined with a Lambertian or subsurface 
	/// diffuse scattering model to account for the multiple scattering.
	/// This class was revised in 2018 to account for an arbitrary polarized 
	/// phase function, rather than a Rayleigh-based phased function.
	///
	/// Note that anisotropic absorption or scattering is not currently implemented.
	///
	MuellerMatrix
		First_Diffuse_BRDF_Model::
		mueller()
	{
		SETUP();

		throw_backward();
		throw_transmission();

		double n = substrate.n(lambda);

		// Angles of refraction internal to the material... 
		double sin_thetai_internal = sin(thetai) / n;
		double cos_thetai_internal = sqrt(1. - sqr(sin_thetai_internal));
		double sin_thetas_internal = sin(thetas) / n;
		double cos_thetas_internal = sqrt(1. - sqr(sin_thetas_internal));

		// The incident and scattering khat vectors outside the material...
		Vector inKo(sin(thetai), 0., -cos(thetai));
		Vector outKo(sin(thetas)*cos(phis), sin(thetas)*sin(phis), cos(thetas));

		// The incident and scattering khat vectors inside the material...
		// (one of my compilers not have a COMPLEX sin(COMPLEX)!)
		Vector inKi(sin_thetai_internal, 0., -cos_thetai_internal);
		Vector outKi(sin_thetas_internal*cos(phis),
			sin_thetas_internal*sin(phis),
			cos_thetas_internal);

		// Use the absorption coefficient (k) to determine the absorption rate of the host material and 
		// add the absorption coeeficient of the scatterer...
		double k_absorb_host = 4 * pi * substrate.k(lambda) / lambda;
		double k_absorb_i = k_absorb_host + phase_function->absorption_coefficient(inKi)[0][0];
		double k_absorb_o = k_absorb_i;

		double k_scat_i = phase_function->scattering_coefficient(inKi);
		double k_scat_o = k_scat_i;

		// The attenuation coefficient is...
		double k_atten_i = k_absorb_i + k_scat_i;
		double k_atten_o = k_atten_i;

		// The surface normal...
		Vector normal(0., 0., 1.);

		// The following definitions make {s,p,k} right handed...
		Vector cnormal = normal;
		Vector inSi = perpto(inKi, cnormal);
		Vector inPi = perpto(inKi, inSi);
		Vector outSi = perpto(outKi, cnormal);
		Vector outPi = perpto(outKi, outSi);
		Vector inSo = perpto(inKo, cnormal);
		Vector inPo = perpto(inKo, inSo);
		Vector outSo = perpto(outKo, cnormal);
		Vector outPo = perpto(outKo, outSo);
		Vector Perp, inPar, outPar;
		GetBasisVectorsParPerp(inKi, outKi, Perp, inPar, outPar);
		MuellerMatrix inrot = GetJonesRotator(inPar, Perp, inSi, inPi);
		MuellerMatrix outrot = GetJonesRotator(outSi, outPi, outPar, Perp);

		// The transmittance for incident light...
		double field_power_i = n*cos_thetai_internal / cos(thetai);
		MuellerMatrix intrans = MuellerMatrix(stack->t12(thetai, lambda, vacuum, substrate))*field_power_i;

		// The transmittance for scattered light...
		double field_power_s = n*cos_thetas_internal / cos(thetas);
		MuellerMatrix outtrans = MuellerMatrix(stack->t12(thetas, lambda, vacuum, substrate))*field_power_s;

		// The Jacobian converting internal solid angle to external solid angle...
		double jacobian = cos(thetas) / cos_thetas_internal / sqr(n);

		// The cosine of the internal, local scattering angle...
		double cos_scatterangle = inKi*outKi;

		// The probability of scattering that angle per solid angle, if the scattering were scalar...
		MuellerMatrix scatter = outtrans*outrot*phase_function->f(inKi, outKi)*inrot*intrans;

		// The probability that a ray is scattered, and that that scattered ray reaches 
		// the surface, is given by...
		double first = k_scat_i*cos_thetas_internal / (k_atten_o*cos_thetai_internal + k_atten_i*cos_thetas_internal);

		// Return the whole value with factors common to all elements...
		MuellerMatrix result = scatter*first*jacobian / cos(thetas);

		MuellerMatrix indepol = MuellerZero();
		MuellerMatrix outdepol = MuellerZero();
		indepol[0][0] = 1;
		indepol[1][1] = indepol[2][2] = 1. - depoll*cos_thetas_internal / (cos_thetai_internal + cos_thetas_internal);
		indepol[3][3] = 1. - depolc*cos_thetas_internal / (cos_thetai_internal + cos_thetas_internal);
		outdepol[0][0] = 1;
		outdepol[1][1] = outdepol[2][2] = 1. - depoll*cos_thetai_internal / (cos_thetai_internal + cos_thetas_internal);
		outdepol[3][3] = 1. - depolc*cos_thetai_internal / (cos_thetai_internal + cos_thetas_internal);

		return outdepol*result*indepol;
	}

	void Register(const First_Diffuse_BRDF_Model* x)
	{
		static bool Models_Registered = false;
		if (!Models_Registered) {
			Models_Registered = true;

			Register_Model(First_Diffuse_BRDF_Model);
			Register((Polarized_Phase_Function*)0);
		}
	}
	DEFINE_MODEL(First_Diffuse_BRDF_Model, BRDF_Model,"Model for single scattering in a diffuse material under a smooth interface.");
	DEFINE_PARAMETER(First_Diffuse_BRDF_Model, double, depoll, "Linear depolarization coefficient", "0", 0xFF);
	DEFINE_PARAMETER(First_Diffuse_BRDF_Model, double, depolc, "Circular depolarization coefficient", "0", 0xFF);
	DEFINE_PTRPARAMETER(First_Diffuse_BRDF_Model, Polarized_Phase_Function_Ptr, phase_function, "Phase Function", "Unpolarized_Phase_Function", 0xFF);
	DEFINE_PTRPARAMETER(First_Diffuse_BRDF_Model, StackModel_Ptr, stack, "Film stack on substrate", "No_StackModel", 0xFF);
	
} // namespace SCATMECH

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sizedistribution.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_SIZEDISTRIBUTION_H
#define SCATMECH_SIZEDISTRIBUTION_H

#include "inherit.h"

namespace SCATMECH {


	///
	/// @brief Abstract class for a distribution.
	///
	class Distribution : public Model {
	public:

		DECLARE_MODEL();

		virtual double pdf(double d) = 0; ///< Returns the normalized probability density function
	};

	typedef Model_Ptr<Distribution> Distribution_Ptr;
	void Register(const Distribution* x);

	///
	/// @brief Class for a volumetric particle diameter distribution
	///
	class VolumeParticleSizeDistribution : public Model {
	public:
		/// Returns the number of particles having diameter between D and D+dD, divided by dD.
		virtual double volumedensity(double D) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<VolumeParticleSizeDistribution> VolumeParticleSizeDistribution_Ptr;
	void Register(const VolumeParticleSizeDistribution* x);

	///
	/// @brief Volumetric particle diameter distribution that takes ad Distribution_Ptr
	/// 
	class Regular_VolumeParticleSizeDistribution : public Model {
	public:
		/// Returns the number of particles having diameter between D and D+dD, divided by dD.
		virtual double volumedensity(double D) {
			SETUP();
			return distribution->pdf(D)*numberdensity;
		};
		DECLARE_MODEL();
		DECLARE_PARAMETER(Distribution_Ptr, distribution); ///< The normalized size diameter distribution
		DECLARE_PARAMETER(double, numberdensity);          ///< The total number of particles divided by volume
	};

	///
	/// @brief Class for a particle diameter distribution on a surface
	///
	class SurfaceParticleSizeDistribution : public Model {
	public:
		/// Returns the number of particles having diameter between D and D+dD, divided by dD.
		virtual double surfacedensity(double d) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<SurfaceParticleSizeDistribution> SurfaceParticleSizeDistribution_Ptr;
	void Register(const SurfaceParticleSizeDistribution* x);

	///
	/// @brief SurfaceParticleSizeDistribution that takes a Distribution_Ptr
	///
	class Regular_SurfaceParticleSizeDistribution : public SurfaceParticleSizeDistribution {
	public:
		/// Returns the number of particles having diameter between D and D+dD, divided by dD.
		virtual double surfacedensity(double d) {
			SETUP();
			return distribution->pdf(d)*numberdensity;
		}

		DECLARE_MODEL();
		DECLARE_PARAMETER(Distribution_Ptr, distribution); ///< The normalized size diameter distribution
		DECLARE_PARAMETER(double, numberdensity);          ///< The total number of particles divided by area
	};

	///
	/// @brief A log-normal distribution
	///
	class Log_Normal_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, sigma);  ///< The width parameter
		DECLARE_PARAMETER(double, median); ///< The median
	};

	///
	/// @brief A normal distribution
	///
	class Normal_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, mean);  ///< The mean
		DECLARE_PARAMETER(double, sigma); ///< The standard deviation
	};

	///
	/// @brief Distribution defined by W. Weibull, J. Appl. Mech.-Trans. ASME 18, 293–297 (1951).
	///
	class Weibull_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, exponent); 
		DECLARE_PARAMETER(double, scale);
	};

	///
	/// @brief Particle size distribution defined by D. Deirmendjian, Electromagnetic Scattering on Spherical Polydispersions (Elsevier, 1969)
	///
	class Modified_Gamma_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, Lambda);
		DECLARE_PARAMETER(double, mu); 
		DECLARE_PARAMETER(double, gamma); 
	};

	///
	/// @brief Adds two distributions
	///
	class Bimodal_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Distribution_Ptr, distributionA); ///< The first distribution
		DECLARE_PARAMETER(Distribution_Ptr, distributionB); ///< The second distribution
		DECLARE_PARAMETER(double, fractionA);		        ///< The fraction of A
	};

	///
	/// @brief The surface distribution defined by IEST-STD-CC1246E
	///
	class CC1246E_SurfaceParticleSizeDistribution : public SurfaceParticleSizeDistribution {
	public:
		virtual double surfacedensity(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, cleanliness);  ///< The cleanliness parameter
		DECLARE_PARAMETER(double, slope);        ///< Defined as 0.926 by the standard
	protected:
		void setup();
	private:
		double coverage(double slope, double level);
		double cov926;
		double cov;
	};

} // namespace SCATMECH

#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: phasefunction.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_PHASEFUNCTION_H
#define SCATMECH_PHASEFUNCTION_H

#include "inherit.h"
#include "sphrscat.h"
#include "sizedistribution.h"
#include "miescat.h"

namespace SCATMECH {

	///
	/// Phase_Fucntion is an abstract class for handling a scalar phase function.
	///
	class Phase_Function : public Model
	{
	public:
		/// f(double theta) returns the phase function as a function of angle in radians.
		/// The result is normalized to 1 when integrated over solid angle.
		/// That is, the integral of 2*pi*sin(theta)*f(theta) from 0 to pi = 1.
		virtual double f(double theta) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<Phase_Function> Phase_Function_Ptr;
	void Register(const Phase_Function* x);

	///
	/// Polarized_Phase_Function is an abstract class for handling a polarized phase function.
	/// In particular it provides the Mueller matrix phase function, the absorption 
	/// coefficient, and the scattering coefficient as a function of direction.
	///
	class Polarized_Phase_Function : public Model {
	public:
		/// The virtual function f(const Vector& in,const Vector& out) provides
		/// the Mueller matrix phase function (fraction scattered per solid angle) 
		/// from direction "in" to direction "out". The returned value is assumed 
		/// to be normalized so that the integral of the [0][0] element over
		/// solid angle is unity. 
		virtual MuellerMatrix f(const Vector& in, const Vector& out) = 0;

		/// The virtual function absorption_coefficient(const Vector& in) returns
		/// the Naperian (1/e) absorption coeffient, with units of inverse length. 
		virtual MuellerMatrix absorption_coefficient(const Vector& in) = 0;

		/// The virtual function scattering_coefficient(const Vector& in) returns
		/// the Naperian (1/e) scattering coeffient, with units of inverse length.
		virtual double scattering_coefficient(const Vector& in) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<Polarized_Phase_Function> Polarized_Phase_Function_Ptr;
	void Register(const Polarized_Phase_Function* x);

	//
	// Definitions of classes inheriting Phase_Function...
	//

	///
	/// Henyey_Greenstein_Phase_Function handles the phase function defined
	/// by L.G. Henyey and J.L. Greenstein, J. L., Astrophysical Journal, 93, 70 - 83 (1941).
	///
	class Henyey_Greenstein_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g); 
	};

	///
	/// Double_Henyey_Greenstein_Phase_Function handles the double sided 
	/// Henyey-Greenstein phase function.
	///
	class Double_Henyey_Greenstein_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g1);
		DECLARE_PARAMETER(double, g2);
		DECLARE_PARAMETER(double, c);
	};

	///
	/// Isotropic_Phase_Function provides an isotropic phase function.
	///
	class Isotropic_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta) {
			return 0.25 / pi;
		}
		DECLARE_MODEL();
	};

	///
	/// Rayleigh_Phase_Function provides the phase function appropriate
	/// for Rayleigh scattering, but does not handle polarization.
	///
	class Rayleigh_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta) {
			return 0.75*(1. + sqr(cos(theta)));
		}
		DECLARE_MODEL();
	};

	///
	/// Table_Phase_Function uses a tabulated phase function. The
	/// table should be normalized and be a function of angle in degrees.
	///
	class Table_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(Table, table);
	};

	///
	/// Reynolds_McCormick_Phase_Function is the two-parameter phase function
	/// defined by L.0. Reynolds and N.J. McCormick, J. Opt. Soc. Am. 70, 1206 (1980).
	///
	class Reynolds_McCormick_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g);
		DECLARE_PARAMETER(double, alpha);
	protected:
		void setup();
	private:
		double K;
	};

	///
	/// Double_Reynolds_McCormick_Phase_Function is a three parameter Reynolds-McCormick 
	/// phase function with both (1+c)*RM(g,alpha)/2 and (1-c)*RM(-g,alpha)/2, 
	/// where RM(g,alpha) is the regular Reynolds-McCormick phase function.
	///
	class Double_Reynolds_McCormick_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g);
		DECLARE_PARAMETER(double, alpha);
		DECLARE_PARAMETER(double, c);
	protected:
		void setup();
	private:
		double K;
	};

	///
	/// Legendre_Phase_Function is the sum of Legendre polynomials.
	///
	class Legendre_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, c0);
		DECLARE_PARAMETER(double, c1);
		DECLARE_PARAMETER(double, c2);
		DECLARE_PARAMETER(double, c3);
		DECLARE_PARAMETER(double, c4);
		DECLARE_PARAMETER(double, c5);
	};

	///
	/// Gaussian_Phase_Function is a Gaussian in angle.
	///
	class Gaussian_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, width);
		DECLARE_PARAMETER(double, background);
	};

	// ******************************************************************
	// * Definitions of classes inheriting Unpolarized_Phase_Function   *
	// ******************************************************************

	///
	/// Unpolarized_Phase_Function forwards a scalar phase function to 
	/// the Polarized_Phase_Function class. It is a totally depolarizing
	/// phase function.
	///
	class Unpolarized_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(double, log_absorption);
		DECLARE_PARAMETER(double, log_scattering);
	};

	///
	/// CodeWhitney_Phase_Function provides the phase function defined in
	/// A.D. Code and B.A. Whitney, "Polarization from scattering in blobs,"
	/// The Astrophysical Journal 441, 400-407 (1995).
	/// It is adapted to allow for any scalar phase function, not just 
	/// Henyey-Greenstein.
	///
	class CodeWhitney_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(double, log_absorption);
		DECLARE_PARAMETER(double, log_scattering);
		DECLARE_PARAMETER(double, plmax);
		DECLARE_PARAMETER(double, pcmax);
		DECLARE_PARAMETER(double, s);
	};

	///
	/// Constrained_Phase_Function provides the phase function described 
	/// in T.A. Germer, “Polarized single-scattering phase function determined 
	/// for a common reflectance standard from bidirectional reflectance 
	/// measurements,” Proc. SPIE 10655, 1065504 (2018). It uses a cubic spline
	/// interpolation, defining the elements at 0 degrees, 60 degrees, 120 degrees,
	/// and 180 degrees, and accounting for some required symmetries.
	///
	class Constrained_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(double, log_absorption);
		DECLARE_PARAMETER(double, log_scattering);
		DECLARE_PARAMETER(double, m11_0);
		DECLARE_PARAMETER(double, m11_60);
		DECLARE_PARAMETER(double, m11_120);
		DECLARE_PARAMETER(double, m11_180);
		DECLARE_PARAMETER(double, m22_60);
		DECLARE_PARAMETER(double, m22_120);
		DECLARE_PARAMETER(double, m33_0);
		DECLARE_PARAMETER(double, m33_60);
		DECLARE_PARAMETER(double, m33_120);
		DECLARE_PARAMETER(double, m33_180);
		DECLARE_PARAMETER(double, m10_60);
		DECLARE_PARAMETER(double, m10_120);
		DECLARE_PARAMETER(double, m23_60);
		DECLARE_PARAMETER(double, m23_120);
	};

	///
	/// SphericalScatterer_Phase_Function uses a SphericalScatterer to define
	/// a phase function.
	///
	class SphericalScatterer_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);
	protected:
		void setup();

	private:
		double csection;
		double absorption;
		double scattering;
		double k;

		DECLARE_MODEL();
		DECLARE_PARAMETER(SphericalScatterer_Ptr, scatterer);
		DECLARE_PARAMETER(double, volume_density);
	};

	///
	/// Distributed_Mie_Phase_Function averages the Mie scattering phase
	/// function over a distribution.
	///
	class Distributed_Mie_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);
	protected:
		void setup();

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, lambda);
		DECLARE_PARAMETER(dielectric_function, medium);
		DECLARE_PARAMETER(dielectric_function, sphere);
		DECLARE_PARAMETER(VolumeParticleSizeDistribution_Ptr, distribution);
		DECLARE_PARAMETER(double, integralStart);
		DECLARE_PARAMETER(double, integralEnd);
		DECLARE_PARAMETER(double, integralStep);

	private:
		typedef Model_Ptr<MieScatterer> MieScatterer_Ptr;
		typedef std::list<MieScatterer_Ptr> MieList;
		typedef MieList::iterator MieListIterator;
		typedef std::list<double> DensityList;
		typedef DensityList::iterator DensityListIterator;

		MieList mielist;
		DensityList densitylist;
		double norm;
		double scattering;
		double absorption;
	};

}
#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: phasefunction.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "phasefunction.h"
#include <math.h>

namespace SCATMECH {

double
Henyey_Greenstein_Phase_Function::
f(double theta)
{
	SETUP();
	return (1. - sqr(g)) / pow(1. + sqr(g) - 2.*g*cos(theta), 1.5) / 4. / pi;
}

double
Double_Henyey_Greenstein_Phase_Function::
f(double theta)
{
	SETUP();
	double q1 = (1. - sqr(g1)) / pow(1. + sqr(g1) - 2.*g1*cos(theta), 1.5) / 4. / pi;
	double q2 = (1. - sqr(g2)) / pow(1. + sqr(g2) - 2.*g2*cos(theta), 1.5) / 4. / pi;
	return fabs((1 + c) / 2.*q1 + (1 - c) / 2.*q2);
}

double
Table_Phase_Function::
f(double theta)
{
	SETUP();
	return table.value(theta / deg);
}

void
Reynolds_McCormick_Phase_Function::
setup()
{
	Phase_Function::setup();
	double temp1 = pow(1 - sqr(g), 2 * alpha);
	double temp2 = pow(1 + sqr(g), 2 * alpha);
	K = alpha*g / pi*temp1 / (temp2 - temp1);
}

double
Reynolds_McCormick_Phase_Function::
f(double theta)
{
	SETUP();
	return K*pow(1 + sqr(g) - 2 * g*cos(theta), -(alpha + 1));
}

void
Double_Reynolds_McCormick_Phase_Function::
setup()
{
	Phase_Function::setup();
	double temp1 = pow(1 - sqr(g), 2 * alpha);
	double temp2 = pow(1 + sqr(g), 2 * alpha);
	K = alpha*g / pi*temp1 / (temp2 - temp1);
}

double
Double_Reynolds_McCormick_Phase_Function::
f(double theta)
{
	SETUP();
	double q1 = K*pow(1 + sqr(g) - 2 * g*cos(theta), -(alpha + 1));
	double q2 = K*pow(1 + sqr(g) + 2 * g*cos(theta), -(alpha + 1));
	return fabs((1 + c) / 2.*q1 + (1 - c) / 2.*q2);
}

namespace {
	double LegendreP(int i, double x) {
		if (i < 0) return 0;
		if (i == 0) return 1;
		if (i == 1) return x;
		double P1 = x;
		double P2 = 1;
		for (int j = 2;j <= i;++j) {
			double P = ((2.*j - 1.)*x*P1 - (j - 1.)*P2) / j;
			P2 = P1;
			P1 = P;
		}
		return P1;
	}
}

double
Legendre_Phase_Function::
f(double theta)
{
	SETUP();
	double result = 0;
	double costheta = cos(theta);
	if (c0 != 0) result += c0*LegendreP(0, costheta);
	if (c1 != 0) result += 3 * c1*LegendreP(1, costheta);
	if (c2 != 0) result += 5 * c2*LegendreP(2, costheta);
	if (c3 != 0) result += 7 * c3*LegendreP(3, costheta);
	if (c4 != 0) result += 9 * c4*LegendreP(4, costheta);
	if (c5 != 0) result += 11 * c5*LegendreP(5, costheta);
	return result;
}

double Gaussian_Phase_Function::f(double theta)
{
	SETUP();
	double wrad = width*deg;
	double result = (1 - background)*exp(-sqr(theta / wrad) / 2) + background;
	// The following is an approximation good for small widths ...
	double wrad2 = sqr(wrad);
	double wrad4 = sqr(wrad2);
	double wrad6 = wrad4*wrad2;
	double norm = 2 * background + (1 - background) * wrad2 + (background - 1) / 3 * wrad4 + (1 - background) / 15 * wrad6;
	return result / norm / 2 / pi;
}


void Register(const Phase_Function* x)
{
	static bool Models_Registered = false;
	if (!Models_Registered) {
		Models_Registered = true;

		Register_Model(Phase_Function);
		Register_Model(Henyey_Greenstein_Phase_Function);
		Register_Model(Double_Henyey_Greenstein_Phase_Function);
		Register_Model(Reynolds_McCormick_Phase_Function);
		Register_Model(Double_Reynolds_McCormick_Phase_Function);
		Register_Model(Isotropic_Phase_Function);
		Register_Model(Legendre_Phase_Function);
		Register_Model(Rayleigh_Phase_Function);
		Register_Model(Gaussian_Phase_Function);
		Register_Model(Table_Phase_Function);
	}
}


DEFINE_VIRTUAL_MODEL(Phase_Function, Model,
	"General scalar phase function for volume scattering");

DEFINE_MODEL(Henyey_Greenstein_Phase_Function, Phase_Function,
	"Henyey-Greenstein phase function");

DEFINE_MODEL(Double_Henyey_Greenstein_Phase_Function, Phase_Function,
	"Double Henyey-Greenstein phase function with retroreflection");

DEFINE_MODEL(Isotropic_Phase_Function, Phase_Function,
	"Scattering which is isotropic.");

DEFINE_MODEL(Rayleigh_Phase_Function, Phase_Function,
	"Phase function appropriate for Rayleigh scattering");

DEFINE_MODEL(Table_Phase_Function, Phase_Function,
	"Phase function provided by a table of values versus angle");

DEFINE_MODEL(Reynolds_McCormick_Phase_Function, Phase_Function,
	"Phase function described by Reynolds and McCormick, J. Opt. Soc. Am. vol. 70, pp 1206-1212 (1980)");

DEFINE_MODEL(Double_Reynolds_McCormick_Phase_Function, Phase_Function,
	"A Reynolds-McCormick phase function with a opposition effect.");

DEFINE_MODEL(Legendre_Phase_Function, Phase_Function,
	"Phase function defined by a series of up to six Legendre polynomials (l=0 to l=5)");

DEFINE_MODEL(Gaussian_Phase_Function, Phase_Function,
	"Phase function defined by a Gaussian");

DEFINE_PARAMETER(Henyey_Greenstein_Phase_Function, double, g, "Asymmetry parameter (g)", "0.01", 0xFF);

DEFINE_PARAMETER(Double_Henyey_Greenstein_Phase_Function, double, g1, "Asymmetry parameter (g1)", "0.01", 0xFF);
DEFINE_PARAMETER(Double_Henyey_Greenstein_Phase_Function, double, g2, "Asymmetry parameter (g2)", "0.01", 0xFF);
DEFINE_PARAMETER(Double_Henyey_Greenstein_Phase_Function, double, c, "Parameter c (-1<c<1)", "0.1", 0xFF);

DEFINE_PARAMETER(Table_Phase_Function, Table, table, "Tabulated phase function", "0.07957747", 0xFF);

DEFINE_PARAMETER(Reynolds_McCormick_Phase_Function, double, g, "Asymmetry parameter (g)", "0.01", 0xFF);
DEFINE_PARAMETER(Reynolds_McCormick_Phase_Function, double, alpha, "Peaking factor (alpha)", "0.5", 0xFF);

DEFINE_PARAMETER(Double_Reynolds_McCormick_Phase_Function, double, g, "Asymmetry parameter (g)", "0.01", 0xFF);
DEFINE_PARAMETER(Double_Reynolds_McCormick_Phase_Function, double, alpha, "Peaking factor (alpha)", "0.5", 0xFF);
DEFINE_PARAMETER(Double_Reynolds_McCormick_Phase_Function, double, c, "Mixing Factor (c)", "0.1", 0xFF);

DEFINE_PARAMETER(Legendre_Phase_Function, double, c0, "Legendre coefficient (l=0)", "1", 0xFF);
DEFINE_PARAMETER(Legendre_Phase_Function, double, c1, "Legendre coefficient (l=1)", "0.65", 0xFF);
DEFINE_PARAMETER(Legendre_Phase_Function, double, c2, "Legendre coefficient (l=2)", "0.42", 0xFF);
DEFINE_PARAMETER(Legendre_Phase_Function, double, c3, "Legendre coefficient (l=3)", "0", 0xFF);
DEFINE_PARAMETER(Legendre_Phase_Function, double, c4, "Legendre coefficient (l=4)", "0", 0xFF);
DEFINE_PARAMETER(Legendre_Phase_Function, double, c5, "Legendre coefficient (l=5)", "0", 0xFF);

DEFINE_PARAMETER(Gaussian_Phase_Function, double, width, "Width [deg] parameter", "20", 0xFF);
DEFINE_PARAMETER(Gaussian_Phase_Function, double, background, "Background parameter", "0.1", 0xFF);


void Register(const Polarized_Phase_Function* x)
{
	static bool regd = false;
	if (!regd) {
		regd = true;

		Register((Phase_Function*)0);
		Register((Distribution*)0);
		Register_Model(Polarized_Phase_Function);
		Register_Model(Unpolarized_Phase_Function);
		Register_Model(CodeWhitney_Phase_Function);
		Register_Model(Constrained_Phase_Function);
		Register_Model(SphericalScatterer_Phase_Function);
		Register_Model(Distributed_Mie_Phase_Function);
	}
}

MuellerMatrix Unpolarized_Phase_Function::f(const Vector& in, const Vector& out)
{
	SETUP();
	double costheta = unit(in)*unit(out);
	if (costheta > 1.) costheta = 1.;
	if (costheta < -1.) costheta = -1.;
	double theta = acos(costheta);

	MuellerMatrix result = MuellerZero();
	result[0][0] = phase_function->f(theta);
	return result;
}
MuellerMatrix Unpolarized_Phase_Function::absorption_coefficient(const Vector& in)
{
	MuellerMatrix result = MuellerZero();
	result[0][0] = result[1][1] = result[2][2] = result[3][3] = exp(log_absorption);
	return result;
}
double Unpolarized_Phase_Function::scattering_coefficient(const Vector& in)
{
	return exp(log_scattering);
}

MuellerMatrix CodeWhitney_Phase_Function::f(const Vector& in, const Vector& out)
{
	SETUP();
	//
	// Phase function taken from A.D. Code and B.A. Whitney, "Polarization from scattering in blobs,"
	// The Astrophysical Journal 441, 400-407 (1995).
	// Adapted to allow for any phase function, not just Henyey-Greenstein.
	//
	double costheta = unit(in)*unit(out);
	if (costheta > 1.) costheta = 1.;
	if (costheta < -1.) costheta = -1.;
	double theta = acos(costheta);

	MuellerMatrix result = MuellerZero();
	double pf = phase_function->f(theta);

	double thetaf = theta*(1 + 3.13*s*exp(-7.*theta / pi));
	double costhetaf = cos(thetaf);
	double a = (1 + costheta*costheta);
	double af = (1 + costhetaf*costhetaf);

	result[0][0] = result[1][1] = 1.;
	result[0][1] = result[1][0] = plmax*(costheta*costheta - 1) / a;
	result[2][2] = result[3][3] = 2 * costheta / a;
	result[2][3] = -(result[3][2] = pcmax*(1 - costhetaf*costhetaf) / af);
	return result*pf;
}

MuellerMatrix CodeWhitney_Phase_Function::absorption_coefficient(const Vector& in)
{
	MuellerMatrix result = MuellerZero();
	result[0][0] = result[1][1] = result[2][2] = result[3][3] = exp(log_absorption);
	return result;
}

double CodeWhitney_Phase_Function::scattering_coefficient(const Vector& in)
{
	return exp(log_scattering);
}

namespace {
	double spline(double y0, double y1, double y2, double y3, double x)
	{
		static double pi2 = pi*pi;
		static double pi3 = pi*pi2;
		static double pi4 = pi*pi3;
		static double pi5 = pi*pi4;
		double x2 = x*x;
		double x3 = x2*x;

		if (x <= pi / 3) {
			return (5 * pi3*y0 + 27 * x3*(6 * y0 - 9 * y1 + 4 * y2 - y3) +
				9 * pi*x2*(-11 * y0 + 14 * y1 - 4 * y2 + y3)) / (5.*pi3);
		}
		if (x <= 2 * pi / 3) {
			return (9 * pi*x2*(16 * y0 - 34 * y1 + 29 * y2 - 11 * y3) +
				pi3*(14 * y0 - 16 * y1 + 11 * y2 - 4 * y3)
				+ 27 * x3*(-3 * y0 + 7 * y1 - 7 * y2 + 3 * y3) +
				9 * pi2*x*(-9 * y0 + 16 * y1 - 11 * y2 + 4 * y3)) / (5.*pi3);
		}
		return (-9 * sqr(pi - x)*(-3 * x*(y0 - 4 * y1 + 9 * y2) +
			pi*(2 * y0 - 8 * y1 + 13 * y2)) +
			(2 * pi - 3 * x)*(34 * pi2 - 93 * pi*x + 54 * x2)*y3) / (5.*pi3);
	}
}

MuellerMatrix Constrained_Phase_Function::f(const Vector& in, const Vector& out)
{
	SETUP();

	double costheta = unit(in)*unit(out);
	if (costheta>1.) costheta = 1.;
	if (costheta<-1.) costheta = -1.;
	double theta = acos(costheta);

	double m11a = m11_0;
	double m22a = m11_0;
	double m33a = m33_0;

	double m22d = -m11_180;
	double m33d = m33_180;

	MuellerMatrix result = MuellerZero();
	result[0][0] = 1.;
	result[1][0] = result[0][1] = spline(0., m10_60, m10_120, 0, theta);
	result[1][1] = spline(m11a, m11_60, m11_120, m11_180, theta);

	result[2][3] = spline(0., m23_60, m23_120, 0, theta);
	result[3][2] = -result[2][3];
	result[2][2] = spline(m22a, m22_60, m22_120, m22d, theta);
	result[3][3] = spline(m33a, m33_60, m33_120, m33d, theta);

	double pf = phase_function->f(theta);

	return result*pf;
}

MuellerMatrix Constrained_Phase_Function::absorption_coefficient(const Vector& in)
{
	MuellerMatrix result = MuellerZero();
	result[0][0] = result[1][1] = result[2][2] = result[3][3] = exp(log_absorption);
	return result;
}

double Constrained_Phase_Function::scattering_coefficient(const Vector& in)
{
	return exp(log_scattering);
}


MuellerMatrix SphericalScatterer_Phase_Function::f(const Vector& in, const Vector& out)
{
	SETUP();

	return MuellerMatrix(scatterer->jones(in, out) / k) / csection;
}

MuellerMatrix SphericalScatterer_Phase_Function::absorption_coefficient(const Vector& in)
{
	SETUP();
	MuellerMatrix result = MuellerZero();
	result[0][0] = result[1][1] = result[2][2] = result[3][3] = absorption;
	return result;
}
double SphericalScatterer_Phase_Function::scattering_coefficient(const Vector& in)
{
	SETUP();
	return scattering;
}

void SphericalScatterer_Phase_Function::setup()
{
	Polarized_Phase_Function::setup();

	csection = 0.;
	double dt = pi / 1000.;
	k = 2 * pi / scatterer->get_lambda();
	for (double t = 0;t <= pi;t += dt) {
		Vector kin(0, 0, 1);
		Vector kout(0, sin(t), cos(t));
		csection += MuellerMatrix(scatterer->jones(kin, kout) / k)[0][0] * sin(t);
	}
	csection *= 2 * pi*dt;
	scattering = csection*volume_density;
	double extinction = 4 * pi / sqr(k)*volume_density*real(scatterer->jones(Vector(1, 0, 0), Vector(1, 0, 0))[0]);
	absorption = extinction - scattering;
}

void Distributed_Mie_Phase_Function::setup()
{
	Polarized_Phase_Function::setup();

	// Clear the list of sub-models...
	mielist.clear();
	densitylist.clear();

	for (double diameter = integralStart; diameter <= integralEnd; diameter *= 1. + integralStep) {
		mielist.push_back(new MieScatterer);
		mielist.back()->set_lambda(lambda);
		mielist.back()->set_medium(medium);
		mielist.back()->set_sphere(sphere);
		mielist.back()->set_radius(diameter / 2);

		double density = distribution->volumedensity(diameter) * diameter * integralStep;
		densitylist.push_back(density);
	}

	DensityListIterator density = densitylist.begin();
	MieListIterator mie = mielist.begin();
	norm = 0;
	absorption = 0;
	scattering = 0;
	for (; mie != mielist.end();++mie, ++density) {
		norm += mie->get()->Csca()*(*density);
		absorption += mie->get()->Cabs()*(*density);
		scattering += mie->get()->Csca()*(*density);
	}
}

MuellerMatrix
Distributed_Mie_Phase_Function::
f(const Vector& in, const Vector& out)
{
	SETUP();

	double k = 2 * pi / lambda;
	DensityListIterator density = densitylist.begin();
	MieListIterator mie = mielist.begin();
	MuellerMatrix result = MuellerZero();
	for (; mie != mielist.end();++mie, ++density) {
		result += MuellerMatrix(mie->get()->jones(in, out) / k)*(*density);
	}
	return result / norm / sqr(medium.n(lambda));
}

MuellerMatrix
Distributed_Mie_Phase_Function::
absorption_coefficient(const Vector& in)
{
	SETUP();
	MuellerMatrix result = MuellerZero();
	result[0][0] = result[1][1] = result[2][2] = result[3][3] = absorption;
	return result / sqr(medium.n(lambda));
}

double
Distributed_Mie_Phase_Function::
scattering_coefficient(const Vector& in)
{
	SETUP();
	return scattering / sqr(medium.n(lambda));
}


DEFINE_VIRTUAL_MODEL(Polarized_Phase_Function, Model, "Base class for phase functions with polarization");

DEFINE_MODEL(Unpolarized_Phase_Function, Polarized_Phase_Function, "A plug-in of a Phase_Function into the Polarized_Phase_Function class");
DEFINE_PTRPARAMETER(Unpolarized_Phase_Function, Phase_Function_Ptr, phase_function, "Phase function", "Henyey_Greenstein_Phase_Function", 0xFF);
DEFINE_PARAMETER(Unpolarized_Phase_Function, double, log_absorption, "Natural logarithm of absorption coefficient [1/um]", "-20", 0xFF);
DEFINE_PARAMETER(Unpolarized_Phase_Function, double, log_scattering, "Natural logarithm of scattering coefficient [1/um]", "0", 0xFF);

DEFINE_MODEL(CodeWhitney_Phase_Function, Polarized_Phase_Function, "A phase function described by Code and Whitney");
DEFINE_PTRPARAMETER(CodeWhitney_Phase_Function, Phase_Function_Ptr, phase_function, "Phase function", "Henyey_Greenstein_Phase_Function", 0xFF);
DEFINE_PARAMETER(CodeWhitney_Phase_Function, double, log_absorption, "Natural logarithm of absorption coefficient [1/um]", "-20", 0xFF);
DEFINE_PARAMETER(CodeWhitney_Phase_Function, double, log_scattering, "Natural logarithm of scattering coefficient [1/um]", "0", 0xFF);
DEFINE_PARAMETER(CodeWhitney_Phase_Function, double, plmax, "Maximum degree of linear polarization", "1", 0xFF);
DEFINE_PARAMETER(CodeWhitney_Phase_Function, double, pcmax, "Maximum degree of circular polarization", "0", 0xFF);
DEFINE_PARAMETER(CodeWhitney_Phase_Function, double, s, "Circular polarization asymmetry parameter", "0", 0xFF);

DEFINE_MODEL(Constrained_Phase_Function, Polarized_Phase_Function, "A general phase function for an isotropic medium that uses cubic splines and values set at 0 degrees, 60 degrees, 120 degrees, and 180 degrees.");
DEFINE_PTRPARAMETER(Constrained_Phase_Function, Phase_Function_Ptr, phase_function, "Phase function", "Henyey_Greenstein_Phase_Function", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, log_absorption, "Natural logarithm of absorption coefficient [1/um]", "-20", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, log_scattering, "Natural logarithm of scattering coefficient [1/um]", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m11_0, "m11 at 0 degrees", "1", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m11_60, "m11 at 60 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m11_120, "m11 at 120 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m11_180, "m11 (and m22) at 180 degrees", "0", 0xFF);

//DEFINE_PARAMETER(Constrained_Phase_Function, double, m22_0, "m22 at 0 degrees", "1", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m22_60, "m22 at 60 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m22_120, "m22 at 120 degrees", "0", 0xFF);
//DEFINE_PARAMETER(Constrained_Phase_Function, double, m22_180, "m22 at 180 degrees", "0", 0xFF);

DEFINE_PARAMETER(Constrained_Phase_Function, double, m33_0, "m33 at 0 degrees", "1", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m33_60, "m33 at 60 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m33_120, "m33 at 120 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m33_180, "m33 at 180 degrees", "0", 0xFF);

DEFINE_PARAMETER(Constrained_Phase_Function, double, m10_60, "m10 at 60 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m10_120, "m10 at 120 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m23_60, "m23 at 60 degrees", "0", 0xFF);
DEFINE_PARAMETER(Constrained_Phase_Function, double, m23_120, "m23 at 120 degrees", "0", 0xFF);


DEFINE_MODEL(SphericalScatterer_Phase_Function, Polarized_Phase_Function, "A plug-in of a SphericalScatterer into the Polarized_Phase_Function class which gives a Rayleigh-like Mueller matrix");
DEFINE_PTRPARAMETER(SphericalScatterer_Phase_Function, SphericalScatterer_Ptr, scatterer, "Scatterer", "MieScatterer", 0xFF);
DEFINE_PARAMETER(SphericalScatterer_Phase_Function, double, volume_density, "Volume number density [um^-3]", "1E-4", 0xFF);

DEFINE_MODEL(Distributed_Mie_Phase_Function, Polarized_Phase_Function, "A phase function corresponding to size distribution of spherical scatterers");
DEFINE_PARAMETER(Distributed_Mie_Phase_Function, double, lambda, "Wavelength [um]", "0.532", 0xFF);
DEFINE_PARAMETER(Distributed_Mie_Phase_Function, dielectric_function, medium, "Medium surrounding spheres", "(1.33,0)", 0xFF);
DEFINE_PARAMETER(Distributed_Mie_Phase_Function, dielectric_function, sphere, "Sphere medium", "(1,0)", 0xFF);
DEFINE_PTRPARAMETER(Distributed_Mie_Phase_Function, VolumeParticleSizeDistribution_Ptr, distribution, "Size distribution", "Regular_VolumeParticleSizeDistribution", 0xFF);
DEFINE_PARAMETER(Distributed_Mie_Phase_Function, double, integralStart, "Start diameter for integration [um]", "0.1", 0xFF);
DEFINE_PARAMETER(Distributed_Mie_Phase_Function, double, integralEnd, "End diameter for integration [um]", "100", 0xFF);
DEFINE_PARAMETER(Distributed_Mie_Phase_Function, double, integralStep, "Fractional step diameter for integration", "0.01", 0xFF);

} // namespace SCATMECH

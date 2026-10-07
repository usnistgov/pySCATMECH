//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: gaussianbeam.cpp
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

#include "scatmech.h"
#include "gaussianbeam.h"
#include "askuser.h"
#include "matrix3d.h"
#include <cstdlib>

using namespace std;

namespace SCATMECH {

    //
    // Constructor...
    //
    Gaussian_Beam_Instrument_BRDF_Model::
    Gaussian_Beam_Instrument_BRDF_Model()
    {
        model_cs = xyxy;
    }

    JonesMatrix
    Gaussian_Beam_Instrument_BRDF_Model::
    jones()
    {
        SETUP();

        if (lambda!=model->get_lambda()) error("lambda!=model.lambda");
        if (type!=model->get_type()) error("type!=model.type");
        if (substrate.index(lambda)!=model->get_substrate().index(lambda)) error("substrate!=model.substrate");

       JonesMatrix j=JonesZero();

        Vector x(1,0,0);
        Vector z(0,0,1);

		// Rotation matrix for offset polar angle (rotates about y)...
		Matrix Moffset_theta = Matrix(cos(offset_theta*deg), 0, sin(offset_theta*deg), 0, 1, 0, -sin(offset_theta*deg), 0, cos(offset_theta*deg));
		// Rotation matrix for offset azimuth angle (rotates about z)...
		Matrix Moffset_phi = Matrix(cos(offset_phi*deg), sin(offset_phi*deg), 0, -sin(offset_phi*deg), cos(offset_phi*deg), 0,	0, 0, 1);
		// Net rotation matrix for offset...
		Matrix Moffset = Moffset_phi*Moffset_theta;
		
		// Scattering direction...
		Vector vs(sin(thetas)*cos(phis-rotation), sin(thetas)*sin(phis-rotation), cos(thetas));
		if (type == 1 || type == 2) vs.z = -vs.z;

		// Rotation matrix about y for incident angle thetai...
		Matrix M_thetai = Matrix(cos(thetai), 0, -sin(thetai), 0, 1, 0, sin(thetai), 0, cos(thetai));
		// Rotation matrix about z for incident azimuth angle...
		Matrix M_phii = Matrix(cos(-rotation), sin(-rotation), 0,-sin(-rotation), cos(-rotation), 0, 0, 0, 1);
		// Net rotation matrix for incident direction...
		Matrix Minc = M_phii*M_thetai;

		// The focal point...
		Vector focalpoint(focus_x, focus_y, focus_z);
		double k = 2 * pi / lambda;

		if (integralmode > 1) {

			// Uses a circle integral for integration...
			Circle_Integral ci(integralmode);

			// Perform integration...
			for (int i = 0; i < ci.n(); ++i) {
				double alpha = ci.theta(i);
				double r = ci.r(i);
				double w = ci.w(i);

				// Rotation matrix about z for the integration variable azimuth...
				Matrix M_phi = Matrix(cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha), 0, 0, 0, 1);
				// Rotation matrix about y for integration variable polar...
				Matrix M_theta = Matrix(cos(r*xNA*NA), 0, sin(r*xNA*NA), 0, 1, 0, -sin(r*xNA*NA), 0, cos(r*xNA*NA));
				// Net rotation matrix for integration...
				Matrix Mint = M_phi*M_theta;

				// The incident direction for this ray...
				Vector vi = Moffset*Minc*Mint*Vector(0, 0, 1);
				if (type == 2 || type == 3) vi.z = -vi.z;

				// The Gaussian ...
				double gauss = sqrt(2*pi)*sqr(xNA)*exp(-sqr(r*xNA))*NA/lambda;

				// The propagation phase...
				COMPLEX phase = exp(COMPLEX(0, focalpoint*vi*k));

				// The incident and scattering angles for this ray...
				double thetai_ = acos(vi.z);
				double thetas_ = acos(vs.z);

				// A couple weighting factors, including the integration weight...
				double weight = w;

				j += model->JonesDSC(vi, vs, z, x, BRDF_Model::xyxy)*weight*gauss*phase;
			}
			return (j / sqrt(cos(thetai)*cos(thetas)));
		}
		else { // integralmode==1

			// Incident direction of propagation...
			Vector kinhat(sin(thetai)*cos(rotation), sin(thetai)*sin(rotation), -cos(thetai));
			Vector sinhat(sin(rotation), cos(rotation), 0.);
			Vector pinhat(cos(thetai)*cos(rotation), -cos(thetai)*sin(rotation), sin(thetai));

			kinhat = Moffset*kinhat;
			sinhat = Moffset*sinhat;
			pinhat = Moffset*pinhat;

			// Focus position...
			Vector focus(focus_x, focus_y, focus_z);

			// Location of scatterer along beam horizontally...
			double zz = focus*kinhat;
			// ... and transversely...
			double r = sqrt(sqr(focus*sinhat) + sqr(focus*pinhat));

			double w0 = lambda / pi / NA;  // Beam waist at focus
			double zR = pi*sqr(w0)/lambda; // Rayleigh range

			double w = w0*sqrt(1. + sqr(zz / zR)); // Beam waist at zz
			double R = zz + sqr(zR) / zz;          // Radius of wavefront at zz
			double psi = atan(zz / zR);            // 

			COMPLEX j(0, 1);
			// Equation 17.1 from A.E. Siegman, _Lasers_ ...
			COMPLEX u = (zz != 0) ? sqrt(2 / pi)*exp(-j*(k*zz + psi)) / w*exp(-sqr(r / w) - j*k*sqr(r) / (2 * R)) :
								   sqrt(2 / pi) / w*exp(-sqr(r / w));

			// The incident direction for this ray...
			Vector vi = Moffset*Minc*Vector(0, 0, 1);

			return model->JonesDSC(vi, vs, z, x, BRDF_Model::xyxy)*u/sqrt(cos(thetai)*cos(thetas));
		}
    }

    DEFINE_MODEL(Gaussian_Beam_Instrument_BRDF_Model,Instrument_BRDF_Model,"A BRDF model measured with a Gaussian focussed beam.");
    DEFINE_PTRPARAMETER(Gaussian_Beam_Instrument_BRDF_Model,Local_BRDF_Model_Ptr,model,"Model to be integrated","Rayleigh_Defect_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model,int,integralmode,"Order of Gauss-Zernike integral (1-25)","10",0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, xNA, "Extent of integration", "3", 0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, NA, "Gaussian numerical aperture", "0.05", 0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, focus_x, "Position of sample in x [um]", "0", 0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, focus_y, "Position of sample in y [um]", "0", 0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, focus_z, "Position of sample in z [um]", "0", 0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, offset_theta, "Offset polar angle [deg]", "0", 0xFF);
	DEFINE_PARAMETER(Gaussian_Beam_Instrument_BRDF_Model, double, offset_phi, "Offset azimuthal angle [deg]", "0", 0xFF);
}


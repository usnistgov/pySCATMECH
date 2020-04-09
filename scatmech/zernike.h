//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: zernike.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_ZERNIKE_H
#define SCATMECH_ZERNIKE_H

namespace SCATMECH {

	double Radial_Zernike_Polynomial(int n,int m,double rho);

	double Zernike_Polynomial(int n,int m,double x,double y);
	double Sin_Zernike_Polynomial(int n,int m,double x,double y);
	double Cos_Zernike_Polynomial(int n,int m,double x,double y);

};

	
#endif

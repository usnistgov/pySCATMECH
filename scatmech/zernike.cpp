//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: zernike.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "zernike.h"
#include "bobvlieg.h"

using namespace std;

namespace SCATMECH { 

	using namespace BobVlieg_Supp;

	double Radial_Zernike_Polynomial(int n,int l,double rho)
	{
		int m = abs(l);
		if ((n-m)%2==1) return 0.;

		double sum = 0;
		for (int s=0;s<=(n-m)/2;++s) {
			sum += mpow(s)*Fact(n-s)/Fact(s)/Fact((n+m)/2-s)/Fact((n-m)/2-s)*pow(rho,(double)(n-2*s));		
		}
		return sum;
	};

	double Zernike_Polynomial(int n,int l,double x,double y)
	{
		double r2 = sqr(x)+sqr(y);
		if (r2>1) return 0.;
		
		double r = sqrt(r2);
		double theta = atan2(y,x);

		if (l>=0) return Radial_Zernike_Polynomial(n,l,r)*cos(l*theta);
		else return Radial_Zernike_Polynomial(n,l,r)*sin(l*theta);
	}

	double Sin_Zernike_Polynomial(int n,int l,double x,double y)
	{
		double r2 = sqr(x)+sqr(y);
		if (r2>1) return 0.;
		
		double r = sqrt(r2);
		double theta = atan2(y,x);

		return Radial_Zernike_Polynomial(n,l,r)*sin(l*theta);
	}

	double Cos_Zernike_Polynomial(int n,int l,double x,double y)
	{
		double r2 = sqr(x)+sqr(y);
		if (r2>1) return 0.;
		
		double r = sqrt(r2);
		double theta = atan2(y,x);

		return Radial_Zernike_Polynomial(n,l,r)*cos(l*theta);
	}

};

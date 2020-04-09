//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: zernikeexpansion.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef ZERNIKEEXPANSION_H
#define ZERNIKEEXPANSION_H

#include "brdf.h"
#include <list>

namespace SCATMECH {

	///
	/// The class ZernikeExpansion_BRDF_Model implements the theory of 
	/// J.J. Koenderink and A.J. van Doorn, "Phenomenological description of bidirectional surface reflection," J. Opt. Soc. Am. A vol. 15, pp. 2903-2912 (1998)
	/// as extended by
	/// T.A. Germer, "Full four-dimensional and reciprocal Mueller matrix bidirectional reflectance distribution function of sintered polytetrafluoroethylene,"
	/// Appl. Opt., to be submitted (2017)
 	///
	class ZernikeExpansion_BRDF_Model : public BRDF_Model {
	public:

		// List of the parameters ...
		DECLARE_MODEL();
		DECLARE_PARAMETER(std::string,coefficientfile); /// File containing indices and coefficients
		DECLARE_PARAMETER(double,scale);                /// An overall scale factor (if the coefficients correspond to unit reflectance)

		// Function evaluating Mueller matrix
		MuellerMatrix mueller();

	protected:

		// Function performing one-time code, namely reading the coefficient file
		void setup();

		// Internal structure holding indices and coefficient 
		struct Coefficient {
			Coefficient(short _i,short _j,short _n,short _m,short _k,short _l,short _p, double _coeff) : i(_i),j(_j),n(_n),m(_m),k(_k),l(_l),p(_p),coeff(_coeff) {}
			short i,j,n,m,k,l,p;
			double coeff;
		};

		// The list of coefficients
		std::list<Coefficient> coeff;
	};

}

#endif

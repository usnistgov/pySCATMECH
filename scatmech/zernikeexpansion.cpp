//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: zernikeexpansion.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "zernike.h"
#include "zernikeexpansion.h"
#include "bobvlieg.h"

using namespace std;
using namespace SCATMECH;

namespace SCATMECH {

// The azimuthal function 
static double az(int k,double phi)
{
	if (k>0) return cos(k*phi);
	if (k<0) return -sin(k*phi);
	else  return 1./sqrt(2.0);
}

void ZernikeExpansion_BRDF_Model::setup()
{
	// Standard procedure is for setup() to call the parent class setup()
	BRDF_Model::setup();

	// The model uses the xyxy basis, rather than the default psps
	model_cs = xyxy; 

	// Open the coefficient file
	ifstream_with_comments cfile(coefficientfile.c_str());

	// Throw exception if the coefficient file doesn't open
	if (!cfile) error("Cannot open coefficient file.");

	// Clear the list of coefficients
	coeff.clear();

	// Read the coefficients until the end of the file
	while (!cfile.eof()) {
		short i,j,n,m,k,l,p;
		double c;

		// Read a line
		string line = cfile.getstrline();

		// If the line isn't blank...
		if (line.size()!=0) {

			// Create a stream from the line
			istringstream sline(line.c_str());

			// The following allows for comma or whitespace delimited data...
			sline >> i;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> j;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> n;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> m;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> k;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> l;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> p;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> c;

			// If there was a failure, then throw an exception
			if (sline.fail()) error("Error reading a line in coefficient file");

			// Push the coefficient onto the list of coefficients
			coeff.push_back(Coefficient(i,j,n,m,k,l,p,c));
		}
	}
	return;
}

MuellerMatrix ZernikeExpansion_BRDF_Model::mueller()
{
	// It is standard procedure for all Models to call SETUP() at beginning of 
	// routines that use 
	SETUP();

	// This code does not accept anything for transmission mode or for
	// radiation upwelling in the material
	throw_transmission();
	throw_backward();

	// Start the sum with zero
	MuellerMatrix result = MuellerZero();
	double brdf = 0;

	// Incident and scattering rho's and phi's
	double rhoi = sqrt(2.)*sin(thetai/2.);
	double rhor = sqrt(2.)*sin(thetas/2.);
	double phii = pi-rotation;
	double phir = phis-rotation;

	// Iterate through all the coefficients...
	for (list<Coefficient>::iterator q=coeff.begin();q!=coeff.end();++q) {
		short i = q->i;
		short j = q->j;
		short n = q->n;
		short m = q->m;
		short k = q->k;
		short l = q->l;
		short p = q->p;
		double coeff = q->coeff;

		// The following just evaluates lambda^p
		double power = p==0?1.:p==1?lambda:p==2?sqr(lambda):p==3?sqr(lambda):pow(lambda,p);

		if (i==1&&j==1) {
			// The basis set for the unpolarized BRDF...
			double prefact = 0.5/pi*sqrt((n+1.)*(m+1.)/((n==0||(n==m&&l==0))?4.:((n==m||l==0)?2.:1.)))*power;
			brdf += coeff * prefact * (Radial_Zernike_Polynomial(n,l,rhoi)*Radial_Zernike_Polynomial(m,l,rhor) +
						 Radial_Zernike_Polynomial(m,l,rhoi)*Radial_Zernike_Polynomial(n,l,rhor)) * cos(l*(phir-phii));
		} else if (i==j) {
			// The basis set for the normalized diagonal elements...
			double prefact = 1./sqrt(2.+2.*(n==m?1:0)*(k==l?1:0))*sqrt((n+1.)*(m+1.))/pi*power;
			result[i-1][j-1] += coeff * prefact * (Radial_Zernike_Polynomial(n,k,rhoi)*Radial_Zernike_Polynomial(m,l,rhor)*az(k,phii)*az(l,phir) +
												 Radial_Zernike_Polynomial(m,l,rhoi)*Radial_Zernike_Polynomial(n,k,rhor)*az(k,phir)*az(l,phii));
		} else if (i<j) {		
			// The baiss set for the off-diagonal elements...
			double prefact = sqrt((n+1.)*(m+1.))/pi*power;
			result[i-1][j-1] += coeff * prefact * Radial_Zernike_Polynomial(n,k,rhoi)*Radial_Zernike_Polynomial(m,l,rhor)*az(k,phii)*az(l,phir);
			result[j-1][i-1] += coeff * prefact * Radial_Zernike_Polynomial(n,k,rhor)*Radial_Zernike_Polynomial(m,l,rhoi)*az(k,phir)*az(l,phii);
		} else throw("i>j in the coefficient file");
	}

	// The normalized Mueller matrix has 1 here...
	result[0][0] = 1.; 
	// The following three elements have sign changes upon transpose...
	result[2][0] = -result[2][0];
	result[2][1] = -result[2][1];
	result[3][2] = -result[3][2];

	// Change to un-normalized Mueller matrix...
	result *= brdf*scale;

	// And return!
	return result;
}

//
// These are required for Models in SCATMECH...
//
DEFINE_MODEL(ZernikeExpansion_BRDF_Model,BRDF_Model,"Model using the model of Koenderink and van Doorn, extended by Germer to include the Mueller matrix");
DEFINE_PARAMETER(ZernikeExpansion_BRDF_Model,string,coefficientfile,"File containing coefficients","",0xFF);
DEFINE_PARAMETER(ZernikeExpansion_BRDF_Model,double,scale,"Scale factor","1",0xFF);


} // namespace SCATMECH;

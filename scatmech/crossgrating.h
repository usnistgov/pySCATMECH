//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crossgrating.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef CROSSGRATING_H
#define CROSSGRATING_H

#include <vector>
#include <algorithm>
#include "matrixmath.h"
#include "inherit.h"
#include "dielfunc.h"
#include "grating.h"
#include "scattabl.h"

namespace SCATMECH {

    class CrossGrating : public Model {
        public:

            CrossGrating() : order1(1),order2(1),lambda(0.5) {}

            // Return the fourier decomposition Toeplitz matrices for specific order of a specific level.
            // Returns CFARRAY(2*order1+1,2*order2+1,2*order1+1,2*order2+1,levels);
            // level = 0 is the closest level to the incident direction.
            
			// type = 0 returns the inverse of the Fourier matrix of eps_z ([[eps_z]]^-1) ...
			CFARRAY get_EPS0() {
                SETUP();
                return EPS0;
            }
            // type = 11 returns the inverse of the Fourier matrix of 1/eps_1 ([[1/eps_1]]^-1) ...
            CFARRAY get_EPS11() {
                SETUP();
                return EPS11;
            }
			// type = 11 returns the inverse of the Fourier matrix of 1/eps_2 ([[1/eps_2]]^-1) ...
			CFARRAY get_EPS12() {
				SETUP();
				return EPS12;
			}

            // type = 2 returns Eq. (8) of LF Li, with eps -> eps_1
            CFARRAY get_EPS2() {
                SETUP();
                return EPS2;
            }
            // type = 3 returns Eq. (9) of LF Li, with eps -> eps_2
            CFARRAY get_EPS3() {
                SETUP();
                return EPS3;
            }

			// type = 0 returns the inverse of the Fourier matrix of eps_z ([[eps_z]]^-1) ...
			CFARRAY get_MU0() {
				SETUP();
				return MU0;
			}
			// type = 11 returns the inverse of the Fourier matrix of 1/eps_1 ([[1/eps_1]]^-1) ...
			CFARRAY get_MU11() {
				SETUP();
				return MU11;
			}
			// type = 11 returns the inverse of the Fourier matrix of 1/eps_2 ([[1/eps_2]]^-1) ...
			CFARRAY get_MU12() {
				SETUP();
				return MU12;
			}

			// type = 2 returns Eq. (8) of LF Li, with eps -> eps_1
			CFARRAY get_MU2() {
				SETUP();
				return MU2;
			}
			// type = 3 returns Eq. (9) of LF Li, with eps -> eps_2
			CFARRAY get_MU3() {
				SETUP();
				return MU3;
			}


            // Return the thickness of a specific level (counting from TOP)...
            double get_thick(int level) {
                SETUP();
                return thick(level);
            }

            // Get the number of levels in the structure (used by CrossRCW_Model)
            int get_levels() {
                SETUP();
                return levels;
            }

            // Set and get the wavelength (used by CrossRCW_Model)
            void set_lambda(double _lambda) {
                lambda = _lambda;
                set_recalc(0xFF);
            }
            double get_lambda() const {
                return lambda;
            }

            // Set and get the order in the 1st coordinate (used by CrossRCW_Model)
            void set_order1(int _order1) {
                order1 = _order1;
                set_recalc(0xFF);
            }
            int get_order1() {
                return order1;
            }

            // Get the order in the 2nd coordinate (used by CrossRCW_Model)
            void set_order2(int _order2) {
                order2 = _order2;
                set_recalc(0xFF);
            }
            int get_order2() {
                return order2;
            }

            // Get the parallelogram angle (in degrees) (used by CrossRCW_Model)
            double get_zeta() {
                SETUP();
                return zeta;
            }

            // Get the lattice constant in the 1st coordinate (used by CrossRCW_Model)
            double get_d1() {
                SETUP();
                return d1;
            }

            // Get the lattice constant in the 2nd coordinate (used by CrossRCW_Model)
            double get_d2() {
                SETUP();
                return d2;
            }

        protected:

            void setup();

            // Overloads the Model function to add further functionality...
            virtual void set_parameter_base(
                const STRING& parameter, ///< The parameter name
                const STRING& value      ///< String represention of a value
            );

            // Child class must allocate and fill the Fourier factorization of the grating
            // Each of these is dimensioned (2*order1+1,2*order2+1,2*order1+1,2*order2+1,levels)
            CFARRAY EPS0;  // EPS0(i,j,k,l,level) is the inverse of the Fourier matrix of eps ([[eps_3]]^-1) ...
            CFARRAY EPS11;  // EPS11(i,j,k,l,level) is the inverse of the Fourier matrix of 1/eps_1 ([[1/eps_1]]^-1)
			CFARRAY EPS12;  // EPS12(i,j,k,l,level) is the inverse of the Fourier matrix of 1/eps_2 ([[1/eps_2]]^-1)
			CFARRAY EPS2;  // EPS2(i,j,k,l,level) is Eq. (8) of LF Li with eps --> eps_1
            CFARRAY EPS3;  // EPS3(i,j,k,l,level) is Eq. (9) of LF Li with eps --> eps_2
			CFARRAY MU0;  // MU0(i,j,k,l,level) is the inverse of the Fourier matrix of mu ([[mu_3]]^-1) ...
			CFARRAY MU11;  // MU11(i,j,k,l,level) is the inverse of the Fourier matrix of 1/mu_1 ([[1/mu_1]]^-1)
			CFARRAY MU12;  // MU12(i,j,k,l,level) is the inverse of the Fourier matrix of 1/mu_2 ([[1/mu_2]]^-1)
			CFARRAY MU2;  // MU2(i,j,k,l,level) is Eq. (8) of LF Li with eps --> mu_1
			CFARRAY MU3;  // MU3(i,j,k,l,level) is Eq. (9) of LF Li with eps --> mu_2

            // It is the responsibility of any child class to set these values...
            int levels;	   // <-- The number of levels in the structure
            double zeta;   // <-- The angle that the second coordinate of the unit cell makes with the y axis
            double d1;     // <-- The lattice constant along the first coordinate, which is along the x axis.
            double d2;     // <-- The lattice constant along the second coordinate.
            DFARRAY thick;  // Child class must allocate and set this as thick(nlevels). thick(1) is nearest the substrate.

            // It is the responsibility of CrossRCW_Model to set these values before setup()...
            double lambda;
            int order1;
            int order2;

        protected:

            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,medium_i);
            DECLARE_PARAMETER(dielectric_function,medium_t);
    };

    //
    // A Gridded_CrossGrating is a CrossGrating for which child classes define the real space
    // dielectric constant on a grid and then perform a Fourier factorization.
    //
    class Gridded_CrossGrating : public CrossGrating {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,zeta);
            DECLARE_PARAMETER(double,d1);
            DECLARE_PARAMETER(double,d2);
            DECLARE_PARAMETER(int,grid1);
            DECLARE_PARAMETER(int,grid2);

		public: 
			Gridded_CrossGrating() : nonmagnetic(true), isotropic(true) {}
        protected:
            void setup() {
                CrossGrating::setup();

                CrossGrating::zeta = zeta;
                CrossGrating::d1 = d1;
                CrossGrating::d2 = d2;

                coszeta=cos(zeta*deg);
                sinzeta=sin(zeta*deg);
            }

            // Overloads the Model function to add further functionality...
            virtual void set_parameter_base(
                const STRING& parameter, ///< The parameter name
                const STRING& value      ///< String represention of a value
            );

            // It is up to the child class to allocate and set this [with eps.allocate(grid1,grid2,levels)]
			// if isotropic is set:
            CFARRAY eps;

			// It is up to the child class to allocate and set these [with, e.g., eps1.allocate(grid1,grid2,levels)] 
			// if the isotropic is unset or nonmagnetic is unset;
			CFARRAY eps1, eps2, eps3;
			CFARRAY mu;
			CFARRAY mu1, mu2, mu3;

            // Child class must call FourierFactorize at end of setup(), after setting eps
            void FourierFactorize();

			// The following is called by FourierFactorize...
			void FourierFactorize(CFARRAY& eps1, CFARRAY& eps2, CFARRAY& eps3, CFARRAY& _EPS0, CFARRAY& _EPS11, CFARRAY& _EPS12, CFARRAY& _EPS2, CFARRAY& _EPS3);

            // epsilon returns the dielectric function at lambda in form N+iK
            COMPLEX epsilon(const dielectric_function& e) {
                return (COMPLEX)(e.epsilon(lambda));
            }

            double coszeta,sinzeta;
            // The following gives the x and y coordinates of an index (i,j) in eps.
            void getxy(int i,int j,double &x, double &y);
			
			bool isotropic;
			bool nonmagnetic;
    };

    void Register(const CrossGrating* x);
    typedef Model_Ptr<CrossGrating> CrossGrating_Ptr;

}

#endif

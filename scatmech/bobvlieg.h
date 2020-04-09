//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: bobvlieg.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_BOBVLIEG_H
#define SCATMECH_BOBVLIEG_H

#include <stdio.h>
#include <vector>
#include "local.h"
#include "scatmatrix.h"
#include "filmtran.h"

namespace SCATMECH {

    class Bobbert_Vlieger_BRDF_Model: public Local_BRDF_Model
    {
        public:

            COMPLEX Epp(double thetai,double thetas,double phis);
            COMPLEX Eps(double thetai,double thetas,double phis);
            COMPLEX Esp(double thetai,double thetas,double phis);
            COMPLEX Ess(double thetai,double thetas,double phis);

            // The electric field in the region D (see B&V, Fig. 1)
            CVector EField(double thetai,const JonesVector& inpol,const Vector r);

            DECLARE_MODEL();

            // Index of particle...
            DECLARE_PARAMETER(dielectric_function,sphere);

            // Radius of particle...
            DECLARE_PARAMETER(double,radius);
            // Coatings on sphere...
            DECLARE_PARAMETER(StackModel_Ptr,spherecoat);

            // Substrate coatings...
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            // Distance between particle and substrate...
            DECLARE_PARAMETER(double,delta);
            // Scattering order (-1 means exact)...
            DECLARE_PARAMETER(int,order);
            // Maximum chosen l value (negative value signifies add to nominal)...
            DECLARE_PARAMETER(int,lmax);
            // Norm_Inc_Approx =0 for full solution, =1 to assume that r(theta)=r(0)
            DECLARE_PARAMETER(int,Norm_Inc_Approx);
            // Number of iterative improvement steps
            DECLARE_PARAMETER(int,improve);

        public:
            Bobbert_Vlieger_BRDF_Model();

            // The real part of the following are the partial extinction cross
            // sections for S or P polarization, respectively, where the "partial"
            // means that when type=0, the downward reflectance is reduced by the
            // result, when type=1, the downward transmittance is reduced by the result,
            // when type=2, the upward reflectance is reduced by the result, and
            // when type=3, the upward transmittance is reduced by the result.
            // The total extinction cross section can be obtained by summing these
            // values for type=0 and type=1 or for type=2 and type=3.
            COMPLEX PartialExtinctionS(double theta);
            COMPLEX PartialExtinctionP(double theta);

        protected:
            void setup();

            JonesMatrix jonesDSC();

        private:

            COMPLEX N0;          // Index of substrate
            COMPLEX N2;          // Index of particle
            COMPLEX N3;          // Index of particle film

            //Private bookkeeping...
            double  old_thetai;  // Last value of incident angle
            double  old_thetas;  // Last value of scattering angle
            double  old_phis;    // Last value of scattering azimuth angle
            int     LMAX;        // Maximum order of spherical harmonic
            int     BH_LMAX;     // Bohren and Huffman recommended LMAX (used for output)
            double  q;           // 2*PI/lambda*particle radius
            double  qq;          // 2*PI/lambda*particle-surface separation
            double  k;           // 2*PI/lambda
            double  d;           // Sum of thickness of all coatings on sphere
            int     sqrsize;     // Dimension of matrices and vectors
            double	r0;			 // Radius of inner sphere...

            //Data arrays...
            //The following are kept after setup()...

            ScatterTMatrix ScatMatrix;          // The scattering matrix
            ScatterTMatrix ScatMatrixInverse;   // The inverse of the scattering matrix
            ScatterTMatrix A;				// The reflection matrix
            std::vector<COMPLEX> Wp,Ws;     // The vectors W (depends upon thetai)
            std::vector<COMPLEX> Zp,Zs;     // The vectors Z (depends upon thetas)
            std::vector<COMPLEX> Vp,Vs;     // The vectors V (scattering vector)
            std::vector<COMPLEX> VSRp,VSRs; // The vectors V^SR (scattered then reflected)
            std::vector<COMPLEX> eIP;       // The vector Phi (depends upon phis)

            //std::vector<double> lvector;    // sqrt((2.*ll+1.)/(ll*(ll+1.)))

            // Function which returns a unique number for the set (l,m,f)
            int index(int l,int m, int f) const {
                return 2*(l*l+l-1+m)+f;
            }
            int index(int l,int m) const {
                return l*l+l-1+m;
            }

            //Private Functions...
            COMPLEX VIp(int l,int m,int f,double thetai) const;
            COMPLEX VIRp(int l,int m,int f,double thetai) const;
            COMPLEX VITp(int l,int m,int f,double thetai) const;
            COMPLEX VIs(int l,int m,int f,double thetai) const;
            COMPLEX VIRs(int l,int m,int f,double thetai) const;
            COMPLEX VITs(int l,int m,int f,double thetai) const;

            COMPLEX Bmatrix(int l_,int m_,int f_,int l,int m,int f) const;
            void    Amatrix(ScatterTMatrix& A);
            void    invert_block_diagonal(ScatterTMatrix& AA);
            void    set_geometry(double thetai,double thetas,double phis);

            COMPLEX E(std::vector<COMPLEX>& W,std::vector<COMPLEX>& Z);
            void    calculate_W(double thetai);
            void    calculate_Z(double thetas);
            void    calculate_eIP(double phis);

            COMPLEX	PartialExtinction(double theta,int pol);

            void    iterative_improvement(ScatterTMatrix& Ainv, ScatterTMatrix& A,
                                          std::vector<COMPLEX>& b,std::vector<COMPLEX>& x);
    };

    class Gauss_Laguerre_Integration { 
        public:
            static double* weights[];
            static double* zeros[];
    };


    namespace BobVlieg_Supp {

        static const int efield=0;
        static const int hfield=1;

        double  Fact(int j);
        double  sqrtFact(int j);
        double  mpow(int m);  // (-1)**m
        COMPLEX ipow(int m);  // (i)**m
        COMPLEX LegendreP(int l,int m,COMPLEX x);
        COMPLEX Legendre(int l,int m,COMPLEX x);
        COMPLEX Ptilde(int l,int m,COMPLEX x);
        COMPLEX Jhalf(int l,COMPLEX rho);
        COMPLEX Yhalf(int l,COMPLEX rho);
        COMPLEX j(int l,COMPLEX rho);
        COMPLEX y(int l,COMPLEX rho);
        COMPLEX h(int l,COMPLEX rho);
        COMPLEX psi(int l,COMPLEX rho);
        COMPLEX zeta(int l,COMPLEX rho);
        COMPLEX psi_(int l,COMPLEX rho);
        COMPLEX zeta_(int l,COMPLEX rho);
        COMPLEX chi(int l,COMPLEX rho);
        COMPLEX chi_(int l,COMPLEX rho);
        COMPLEX dLegendreCosx_dx(int l,int m,COMPLEX cosx);
        COMPLEX P_tilde(int l,int m,COMPLEX x);
        COMPLEX d(int l,int m,int m_,COMPLEX cosalpha2,COMPLEX sinalpha2);
        COMPLEX dplus(int l,int m,COMPLEX cosalpha2,COMPLEX sinalpha2);
        COMPLEX dminus(int l,int m,COMPLEX cosalpha2,COMPLEX sinalpha2);
        COMPLEX arccosine(const COMPLEX& a);
        COMPLEX arcsine(const COMPLEX& a);
        void cvert(COMPLEX** V, int N);
        COMPLEX ContinuedFraction_Jratio(double nu, const COMPLEX& x);

        extern double lvector[];
        static const COMPLEX cI(0.,1.);
    }

} // namespace SCATMECH



#endif

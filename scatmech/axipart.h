//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: axipart.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_AXIPART_H
#define SCATMECH_AXIPART_H

#include <stdio.h>
#include <vector>
#include "local.h"
#include "bobvlieg.h"
#include "tmatrix.h"


namespace SCATMECH {


    class Axisymmetric_Particle_BRDF_Model: public Local_BRDF_Model
    {
        public:

            COMPLEX Epp(double thetai,double thetas,double phis);
            COMPLEX Eps(double thetai,double thetas,double phis);
            COMPLEX Esp(double thetai,double thetas,double phis);
            COMPLEX Ess(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(Axisymmetric_Shape_Ptr,Shape);
            DECLARE_PARAMETER(dielectric_function,particle);
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,delta);
            DECLARE_PARAMETER(int,lmax);
            DECLARE_PARAMETER(int,mmax);
            DECLARE_PARAMETER(int,order);
            DECLARE_PARAMETER(int,Norm_Inc_Approx);
            DECLARE_PARAMETER(int,improve);

        public:
            Axisymmetric_Particle_BRDF_Model();

            MuellerMatrix Specular(double theta);

            //COMPLEX PartialExtinctionS(double theta);
            //COMPLEX PartialExtinctionP(double theta);

        protected:
            void setup();

            JonesMatrix jonesDSC();

        private:
            //Primary parameters (accessible through set_ and get_ functions)...
            dielectric_function n3;

            COMPLEX N0;          // Index of substrate
            COMPLEX N1;          // Index of substrate film
            COMPLEX N2;          // Index of particle
            COMPLEX N3;          // Index of particle film

            //Private bookkeeping...
            double  old_thetai;  // Last value of incident angle
            double  old_thetas;  // Last value of scattering angle
            double  old_phis;    // Last value of scattering azimuth angle
            int     LMAX;         // Maximum order of spherical harmonic
            int     MMAX;        // Maximum degree of spherical harmonic
            int     BH_LMAX;     // Bohren and Huffman recommended lmax (used for output)
            double  qq;          // 2*PI/lambda*particle-surface separation
            double  k;           // 2*PI/lambda
            int     sqrsize;     // Dimension of matrices and vectors
            int     old_LMAX;

            // Function which returns a unique number for the set (l,m,f)
            int index(int l,int m, int f) const {
                return 2*(l*l+l-1+m)+f;
            }
            int index(int l,int m) const {
                return l*l+l-1+m;
            }

            //Data arrays...
            //The following are kept after setup()...

            ScatterTMatrix ScatMatrix;          // The scattering matrix
            ScatterTMatrix ScatMatrixInverse;   // The inverse of the scattering matrix
            std::vector<COMPLEX> Wp,Ws;     // The vectors W (depends upon thetai)
            std::vector<COMPLEX> Zp,Zs;     // The vectors Z (depends upon thetas)
            std::vector<COMPLEX> Vp,Vs;     // The vectors V (scattering vector)
            std::vector<COMPLEX> eIP;       // The vector Phi (depends upon phis)

            //std::vector<double> lvector;    // sqrt((2.*ll+1.)/(ll*(ll+1.)))

            //Private Functions...
            COMPLEX VIp(int l,int m,int f,double thetai) const;
            COMPLEX VIRp(int l,int m,int f,double thetai) const;
            COMPLEX VITp(int l,int m,int f,double thetai) const;
            COMPLEX VIs(int l,int m,int f,double thetai) const;
            COMPLEX VIRs(int l,int m,int f,double thetai) const;
            COMPLEX VITs(int l,int m,int f,double thetai) const;

            void    Bmatrix(ScatterTMatrix& T);
            void    Amatrix(ScatterTMatrix& A);
            void    invert_block_diagonal(ScatterTMatrix& AA);
            void    set_geometry(double thetai,double thetas,double phis);

            COMPLEX E(std::vector<COMPLEX>& W,std::vector<COMPLEX>& Z);
            void    calculate_W(double thetai);
            void    calculate_Z(double thetas);
            void    calculate_eIP(double phis);


            void    iterative_improvement(ScatterTMatrix& Ainv, ScatterTMatrix& A,
                                          std::vector<COMPLEX>& b,std::vector<COMPLEX>& x);

    };


} // namespace SCATMECH



#endif

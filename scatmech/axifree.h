//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: axifree.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_AXIFREE_H
#define SCATMECH_AXIFREE_H

#include <stdio.h>
#include <vector>
#include "tmatrix.h"
#include "scatmatrix.h"


namespace SCATMECH {

    class TMatrix_Axisymmetric_Scatterer : public Free_Space_Scatterer
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(Axisymmetric_Shape_Ptr,Shape);
            DECLARE_PARAMETER(dielectric_function,particle);
            DECLARE_PARAMETER(int,lmax);
            DECLARE_PARAMETER(int,mmax);

        public:
            TMatrix_Axisymmetric_Scatterer();

            virtual JonesMatrix jones(const Vector& kin,const Vector& kout);

        protected:
            void setup();

        private:

            COMPLEX N2;          // Index of particle
            double N0;           // Index of surrounding medium

            double  a;           // Vertical radius of particle

            //Private bookkeeping...
            double  old_thetai;  // Last value of incident angle
            double  old_thetas;  // Last value of scattering angle
            double  old_phis;    // Last value of scattering azimuth angle

            int     LMAX;          // Maximum order of spherical harmonic
            int     MMAX;        // Maximum degree of spherical harmonic
            int     BH_LMAX;     // Bohren and Huffman recommended lmax (used for output)
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

            ScatterTMatrix ScatMatrix;      // The scattering matrix
            std::vector<COMPLEX> Wp,Ws;     // The vectors W (depends upon thetai)
            std::vector<COMPLEX> Zp,Zs;     // The vectors Z (depends upon thetas)
            std::vector<COMPLEX> Vp,Vs;     // The vectors V (scattering vector)
            std::vector<COMPLEX> eIP;       // The vector Phi (depends upon phis)

            std::vector<double> lvector;    // sqrt((2.*ll+1.)/(ll*(ll+1.)))

            //Private Functions...
            COMPLEX VIp(int l,int m,int f,double thetai) const;
            COMPLEX VIs(int l,int m,int f,double thetai) const;

            void    Bmatrix(ScatterTMatrix& T);

            void    set_geometry(double thetai,double thetas,double phis);

            COMPLEX E(std::vector<COMPLEX>& W,std::vector<COMPLEX>& Z);
            void    calculate_W(double thetai);
            void    calculate_Z(double thetas);
            void    calculate_eIP(double phis);
    };


} // namespace SCATMECH



#endif

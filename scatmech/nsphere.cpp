//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: nsphere.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "nsphere.h"
#include "bobvlieg.h"
#include "matrixmath.h"

using namespace SCATMECH;
using namespace std;
using namespace SCATMECH::BobVlieg_Supp;

//
// This solution to the multilayer scattering problem is based upon
// J. Sinzig and M. Quinten, "Scattering and absorption by spherical multilayer particles,"
// Applied Physics A 58 (2), 157-162 (1994).
// There are some errors in the manuscript, namely in Eqs. (9a) and (9b). The superscript on
// S and T should be r-1, not r. This is pretty obvious by looking at Eqs. (4.56), (4.57), and (8.2)
// of C.F. Bohren and D.R.Huffman, "Absorption and Scattering of Light by Small Particles," (Wiley, New York, 1983).
//
namespace SCATMECH {

    MultilayerCoatedMieScatterer::
    MultilayerCoatedMieScatterer()
    {
    }

    void
    MultilayerCoatedMieScatterer::
    setup()
    {
        SphericalScatterer::setup();

        int nlayers = stack->get_n();
        int r = nlayers+1;
        CFARRAY m(r,1);
        CFARRAY n(r+1,1);
        CFARRAY k(r+1,1);
        DFARRAY R(r,1);
        CFARRAY y(r,1);
        DFARRAY x(r,1);
        for (int s=1; s<=r+1; ++s) {
            if (s==1) {
                n(s) = sphere.index(lambda);
            } else if (s==r+1) {
                n(s) = medium.index(lambda);
            } else {
                n(s) = stack->get_e()[s-2].index(lambda);
            }
        }
        for (int s=1; s<=r; ++s) {
            m(s) = n(s+1)/n(s);
        }
        for (int s=1; s<=r+1; ++s) {
            k(s) = 2*pi*n(s)/lambda;
        }

        for (int s=1; s<=r; ++s) {
            if (s==1) {
                R(s)=SphericalScatterer::radius;
                for (int t=0; t<nlayers; ++t) {
                    R(s)-=stack->get_t()[t];
                }
                if (R(s)<0.) throw SCATMECH_exception("Thickness of stack exceeds particle radius.");
            } else {
                R(s)=R(s-1)+stack->get_t()[s-2];
            }
            y(s) = k(s)*R(s);
            x(s) = 2*pi/lambda*R(s);
        }

        //C
        //C*** Series expansion terminated after NSTOP terms
        //C    Logarithmic derivatives calculated from NMX on down
        double XSTOP=x(r)+4.*pow(x(r),0.3333)+2.;

        NMX=(int)(XSTOP+15);
        NSTOP=(int)(XSTOP);

        A.resize(NMX);
        B.resize(NMX);
        E.resize(NMX);
        F.resize(NMX);

        // This loop is done in reverse, since the spherical Bessel functions use a
        // lookup table and it is more efficient if the highest order is calculated,
        // since all orders less than or equal to that order are calculated.
        //for (int n=1; n<=NMX; ++n) {
        for (int n=NMX; n>=1; --n) {

            COMPLEX T = 0.;
            COMPLEX S = 0.;

            for (int s=1; s<=r-1; ++s) {
                // Eq. (8a)...
                T = -(m(s)*psi(n,m(s)*y(s))*(psi_(n,y(s))+T*chi_(n,y(s)))-psi_(n,m(s)*y(s))*(psi(n,y(s))+T*chi(n,y(s))))/
                    (m(s)*chi(n,m(s)*y(s))*(psi_(n,y(s))+T*chi_(n,y(s)))-chi_(n,m(s)*y(s))*(psi(n,y(s))+T*chi(n,y(s))));

                // Eq. (8b)...
                S = -(psi(n,m(s)*y(s))*(psi_(n,y(s))+S*chi_(n,y(s)))-m(s)*psi_(n,m(s)*y(s))*(psi(n,y(s))+S*chi(n,y(s))))/
                    (chi(n,m(s)*y(s))*(psi_(n,y(s))+S*chi_(n,y(s)))-m(s)*chi_(n,m(s)*y(s))*(psi(n,y(s))+S*chi(n,y(s))));
            }

            // Eq. (9a)...
            A[n-1] = -(m(r)*psi(n,m(r)*y(r))*(psi_(n,y(r))+T*chi_(n,y(r)))-psi_(n,m(r)*y(r))*(psi(n,y(r))+T*chi(n,y(r))))/
                     (m(r)*zeta(n,m(r)*y(r))*(psi_(n,y(r))+T*chi_(n,y(r)))-zeta_(n,m(r)*y(r))*(psi(n,y(r))+T*chi(n,y(r))));

            // Eq. (9b)...
            B[n-1] = -(psi(n,m(r)*y(r))*(psi_(n,y(r))+S*chi_(n,y(r)))-m(r)*psi_(n,m(r)*y(r))*(psi(n,y(r))+S*chi(n,y(r))))/
                     (zeta(n,m(r)*y(r))*(psi_(n,y(r))+S*chi_(n,y(r)))-m(r)*zeta_(n,m(r)*y(r))*(psi(n,y(r))+S*chi(n,y(r))));

            // The following two lines ensure the A and B coefficients
            // match the conventions of MieScatterer and CoatedMieScatterer...
            A[n-1] = -A[n-1];
            B[n-1] = -B[n-1];

            E[n-1]=n;
            F[n-1]=(2.*E[n-1]+1.)/(E[n-1]*(E[n-1]+1.));
        }
    }

    //
    // Routine to calculate s1(angle)
    //
    COMPLEX
    MultilayerCoatedMieScatterer::
    s1(double theta)
    {
        SETUP();

        COMPLEX S1=0.;
        double PI0=0.;
        double PI1=1.;
        double AMU=cos(theta);
        int N;
        for (N=0; N<NSTOP; ++N) {
            double PPI=PI1;
            double TAU=E[N]*AMU*PPI-(E[N]+1.)*PI0;
            S1+=F[N]*(A[N]*PPI+B[N]*TAU);
            PI1=((2.*E[N]+1.)*AMU*PPI-(E[N]+1.)*PI0)/E[N];
            PI0=PPI;
        }
        return S1;
    }

    //
    // Routine to calculate s2(angle)
    //
    COMPLEX
    MultilayerCoatedMieScatterer::
    s2(double theta)
    {
        SETUP();

        COMPLEX S2=0.;
        double PI0=0.;
        double PI1=1.;
        double AMU=cos(theta);
        int N;
        for (N=0; N<NSTOP; ++N) {
            double PPI=PI1;
            double TAU=E[N]*AMU*PPI-(E[N]+1.)*PI0;
            S2+=F[N]*(A[N]*TAU+B[N]*PPI);
            PI1=((2.*E[N]+1.)*AMU*PPI-(E[N]+1.)*PI0)/E[N];
            PI0=PPI;
        }
        return S2;
    }

    //
    // Routine to calculate scattering cross section
    //
    double
    MultilayerCoatedMieScatterer::
    Csca()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Eq. (4.61) of Bohren and Huffman
            result += (2.*n+1.)*(norm(A[i])+norm(B[i]));
        }
        result *= 2.*pi/sqr(k);
        return result;
    }

    //
    // Routine to calculate extinction cross section
    //
    double
    MultilayerCoatedMieScatterer::
    Cext()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Eq. (4.62) of Bohren and Huffman
            result += (2.*n+1.)*real(A[i]+B[i]);
        }
        result *= 2.*pi/sqr(k);
        return result;
    }

    //
    // Routine to calculate backscatter cross section
    //
    double
    MultilayerCoatedMieScatterer::
    Cback()
    {
        SETUP();

        COMPLEX result1=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Section (4.6) of Bohren and Huffman
            result1 += ((2.*n+1.)*pow(-1.,(double)n))*(A[i]-B[i]);
        }
        double result = pi/sqr(k)*norm(result1);
        return result;
    }

    DEFINE_MODEL(MultilayerCoatedMieScatterer,SphericalScatterer,
                 "Scattering from a sphere with any number of coatings.");

    DEFINE_PTRPARAMETER(MultilayerCoatedMieScatterer,StackModel_Ptr,stack,"Coating stack on core sphere","No_StackModel",0xFF);

} // namespace SCATMECH;



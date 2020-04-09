//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: miescat.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "miescat.h"

//
// These routines are loosely based upon the program BHMIE found
// in C. F. Bohren and D. R. Huffman, "Absorption and Scattering of Light by
// Small Particles," (Wiley, New York, 1983), as modified by
// B.T.Draine.
//
// This version was further modified to place it in the context of a
// SphericalScatterer C++ object class.  This change necessitated the
// following modifications:
// * Calculate the coeffs. A[i], B[i], E[i], and F[i] at class construction.
// * Place calculation of so s1(angle) and s2(angle) in their own routines.
// * Arrays changed to scalars, and vice versa, to accomodate these changes.
//
// T.A.Germer (15 NOV 99)
//

using namespace std;


namespace SCATMECH {


    MieScatterer::
    MieScatterer()
    {
        NSTOP=0;
        NMX=0;
    }

    void
    MieScatterer::
    setup()
    {
        SphericalScatterer::setup();

        //COMPLEX REFREL=refrel;
        COMPLEX REFREL = m;

        COMPLEX Y=x*REFREL;
        double ymod=abs(Y);
        //C
        //C*** Series expansion terminated after NSTOP terms
        //C    Logarithmic derivatives calculated from NMX on down
        double XSTOP=x+4.*pow(x,0.3333)+2.;
        // Typecaste to int added to prevent warning (TAG 22 JAN 2001)
        NMX=(int)(((XSTOP>ymod)? XSTOP : ymod )+15);

        // Typecaste to int added to prevent warning (TAG 22 JAN 2001)
        NSTOP=(int)(XSTOP);
        // The following were changed to use the SCATMEM functions...
        //A = new COMPLEX[NMX];
        //B = new COMPLEX[NMX];
        //E = new double[NMX];
        //F = new double[NMX];
        A.resize(NMX);
        B.resize(NMX);
        E.resize(NMX);
        F.resize(NMX);

        //Logarithmic derivative D(J) calculated by downward recurrence
        //beginning with initial value (0.,0.) at J=NMX
        vector<COMPLEX> D;
        D.resize(NMX);
        D[NMX-1]=COMPLEX(0.,0.);
        int NN=NMX-1;
        int N;

        for (N=0; N<NN; ++N) {
            E[N]=NMX-N;
            D[NMX-N-2]=(E[N]/Y)-(1./(D[NMX-N-1]+E[N]/Y));
        }

        //Riccati-Bessel functions with real argument X
        //calculated by upward recurrence
        double  PSI0=cos(x);
        double  PSI1=sin(x);
        double  CHI0=-sin(x);
        double  CHI1=cos(x);
        COMPLEX XI1=COMPLEX(PSI1,-CHI1);

        for (N=0; N<NSTOP; ++N) {
            E[N]=N+1;
            F[N]=(2.*E[N]+1.)/(E[N]*(E[N]+1.));
            //C for given N, PSI  = psi_n        CHI  = chi_n
            //C              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
            //C              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
            //C Calculate psi_n and chi_n
            double PSI=(2.*E[N]-1.)*PSI1/x-PSI0;
            double CHI=(2.*E[N]-1.)*CHI1/x-CHI0;
            COMPLEX XI=COMPLEX(PSI,-CHI);

            A[N]=(D[N]/REFREL+E[N]/x)*PSI-PSI1;
            A[N]=A[N]/((D[N]/REFREL+E[N]/x)*XI-XI1);
            B[N]=(REFREL*D[N]+E[N]/x)*PSI-PSI1;
            B[N]=B[N]/((REFREL*D[N]+E[N]/x)*XI-XI1);

            PSI0=PSI1;
            PSI1=PSI;
            CHI0=CHI1;
            CHI1=CHI;
            XI1=COMPLEX(PSI1,-CHI1);
        }
    }

    //
    // Routine to calculate s1(angle)
    //
    COMPLEX
    MieScatterer::
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
    MieScatterer::
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
    MieScatterer::
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
    MieScatterer::
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
    MieScatterer::
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

    DEFINE_MODEL(MieScatterer,SphericalScatterer,
                 "Scattering by a homogeneous sphere.");



} // namespace SCATMECH




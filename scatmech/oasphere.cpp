//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: oasphere.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "scatmech.h"
#include "oasphere.h"
#include "bobvlieg.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {

    Optically_Active_Sphere_Scatterer::
    Optically_Active_Sphere_Scatterer()
    {
        NSTOP=0;
        NMX=0;
    }

    void
    Optically_Active_Sphere_Scatterer::
    setup()
    {
        Free_Space_Scatterer::setup();

        int N;
        COMPLEX cI(0.,1.);

        if (imag(COMPLEX(medium.index(lambda)))!=0) error("medium must be non-absorbing");

        double nmed = real(COMPLEX(medium.index(lambda)));
        COMPLEX nL = COMPLEX(materialleft.index(lambda));
        COMPLEX nR = COMPLEX(materialright.index(lambda));
        COMPLEX mL = nL/nmed;
        COMPLEX mR = nR/nmed;
        COMPLEX m = 2./(1./mR+1./mL);

        k = 2.*pi*nmed/lambda;
        x = k*radius;

        double absmL = abs(mL);
        double absmR = abs(mR);

        double ymod=x*(absmL > absmR ? absmL : absmR);

        //C
        //C*** Series expansion terminated after NSTOP terms
        //C
        double XSTOP=x+4.*pow(x,0.3333)+2.;

        NMX=(int)(((XSTOP>ymod)? XSTOP : ymod )+50);

        NSTOP=(int)(XSTOP);

        COMPLEX mLx = mL*x;
        COMPLEX mRx = mR*x;

        // The logarithmic derivatives...
        vector<COMPLEX> DmLx(NMX);
        vector<COMPLEX> DmRx(NMX);
        for (N=NMX-1; N>=0; --N) {
            DmLx[N] = BobVlieg_Supp::psi_(N,mLx)/BobVlieg_Supp::psi(N,mLx);
            DmRx[N] = BobVlieg_Supp::psi_(N,mRx)/BobVlieg_Supp::psi(N,mRx);
        }

        // The Ricatti-Bessel functions....
        vector<COMPLEX> psix(NMX);
        vector<COMPLEX> xix(NMX);
        for (N=NMX-1; N>=0; --N) {
            psix[N] = BobVlieg_Supp::psi(N,x);
            xix[N] = BobVlieg_Supp::zeta(N,x);
        }

        a.resize(NSTOP);
        b.resize(NSTOP);
        c.resize(NSTOP);
        E.resize(NSTOP);
        F.resize(NSTOP);

        for (N=0; N<NSTOP; ++N) {
            E[N]=N+1;
            F[N]=(2.*E[N]+1.)/(E[N]*(E[N]+1.));

            int n=N+1;
            double nn=n;

            // The following three expressions were obtained by taking the expressions
            // for an, bn, and cn from p. 188 of Bohren and Huffman, making the substitutions
            // given on p. 127.
            COMPLEX denom =
                (-2.*m*sqr(x)*sqr(xix[n-1]) + x*(4.*m*nn + (1. + sqr(m))*x*(DmLx[n] + DmRx[n]))*xix[n-1]*
                 xix[n] - (x*DmLx[n]*(nn + sqr(m)*nn + 2.*m*x*DmRx[n]) + nn*(2.*m*nn + (1. + sqr(m))*x*DmRx[n]))*sqr(xix[n]));

            COMPLEX an = (psix[n]*(x*(2.*m*nn + x*(DmLx[n] + DmRx[n]))*xix[n-1] - (nn*(2.*m*nn + (1. + sqr(m))*x*DmLx[n]) +
                                   x*(nn + sqr(m)*nn + 2.*m*x*DmLx[n])*DmRx[n])*xix[n]) + m*x*psix[n-1]*(-2.*x*xix[n-1] + (2.*nn + m*x*(DmLx[n] + DmRx[n]))*
                                           xix[n]))/denom;

            COMPLEX bn = (psix[n]*(m*x*(2.*nn + m*x*(DmLx[n] + DmRx[n]))*xix[n-1] - (nn*(2.*m*nn + (1. + sqr(m))*x*DmLx[n]) +
                                   x*(nn + sqr(m)*nn + 2.*m*x*DmLx[n])*DmRx[n])*xix[n]) + x*psix[n-1]*(-2.*m*x*xix[n-1] + (2.*m*nn + x*(DmLx[n] + DmRx[n]))*
                                           xix[n]))/denom;

            COMPLEX cn =(cI*m*sqr(x)*(DmLx[n] - DmRx[n])*(psix[n]*xix[n-1] - psix[n-1]*xix[n]))/denom;

            a[N]=an;
            b[N]=bn;
            c[N]=cn;
        }
    }

    JonesMatrix Optically_Active_Sphere_Scatterer::jones(const Vector& kin,const Vector& kout)
    {
        SETUP();

        double costheta=kin*kout/Norm(kin)/Norm(kout);
        if (costheta>1.) costheta=1;
        if (costheta<-1.) costheta=-1;
        double theta = acos(costheta);
        COMPLEX S1=0.;
        COMPLEX S2=0.;
        COMPLEX S3=0.;
        double PI0=0.;
        double PI1=1.;
        double AMU=cos(theta);
        int N;
        for (N=0; N<NSTOP; ++N) {
            double PPI=PI1;
            double TAU=E[N]*AMU*PPI-(E[N]+1.)*PI0;
            S1+=F[N]*(a[N]*PPI+b[N]*TAU);
            S2+=F[N]*(a[N]*TAU+b[N]*PPI);
            S3+=F[N]*c[N]*(PPI+TAU);
            PI1=((2.*E[N]+1.)*AMU*PPI-(E[N]+1.)*PI0)/E[N];
            PI0=PPI;
        }
        return JonesMatrix(S1,S2,S3,-S3);
    }

    //
    // Routine to calculate scattering cross section
    //
    double
    Optically_Active_Sphere_Scatterer::
    CscaL()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Eq. (8.17) of Bohren and Huffman
            result += (2.*n+1.)*(norm(a[i])+norm(b[i])+2*norm(c[i])-2*imag((a[i]+b[i])*conj(c[i])));
        }
        result *= 2.*pi/sqr(k);
        return result;
    }

    double
    Optically_Active_Sphere_Scatterer::
    CscaR()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Eq. (8.17) of Bohren and Huffman
            result += (2.*n+1.)*(norm(a[i])+norm(b[i])+2*norm(c[i])+2*imag((a[i]+b[i])*conj(c[i])));
        }
        result *= 2.*pi/sqr(k);
        return result;
    }

    //
    // Routine to calculate extinction cross section
    //
    double
    Optically_Active_Sphere_Scatterer::
    CextL()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Eq. (8.17) of Bohren and Huffman
            result += (2.*n+1.)*real(a[i]+b[i]-2.*COMPLEX(0,1)*c[i]);
        }
        result *= 2.*pi/sqr(k);
        return result;
    }

    double
    Optically_Active_Sphere_Scatterer::
    CextR()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;
            // Eq. (8.17) of Bohren and Huffman
            result += (2.*n+1.)*real(a[i]+b[i]+2.*COMPLEX(0,1)*c[i]);
        }
        result *= 2.*pi/sqr(k);
        return result;
    }

    //
    // Routine to calculate backscatter cross section
    //
    double
    Optically_Active_Sphere_Scatterer::
    CbackL()
    {
        MuellerMatrix m = jones(Vector(0,0,1),Vector(0,0,-1));
        StokesVector s = m*StokesVector(1,0,0,1);
        return s[0];
    }

    double
    Optically_Active_Sphere_Scatterer::
    CbackR()
    {
        MuellerMatrix m = jones(Vector(0,0,1),Vector(0,0,-1));
        StokesVector s = m*StokesVector(1,0,0,1);
        return s[0];
    }

    DEFINE_MODEL(Optically_Active_Sphere_Scatterer,Free_Space_Scatterer,"Scattering from a homogeneous optically active or circularly dichroic sphere.");
    DEFINE_PARAMETER(Optically_Active_Sphere_Scatterer,dielectric_function,materialleft,"Optical properties of sphere for left circular polarization","(1.59,0)",0xFF);
    DEFINE_PARAMETER(Optically_Active_Sphere_Scatterer,dielectric_function,materialright,"Optical properties of sphere for right circular polarization","(1.59,0)",0xFF);
    DEFINE_PARAMETER(Optically_Active_Sphere_Scatterer,double,radius,"Radius of sphere [um]","1",0xFF);


} // namespace SCATMECH




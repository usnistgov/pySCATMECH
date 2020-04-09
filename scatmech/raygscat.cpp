//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: raygscat.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "scatmech.h"
#include "raygscat.h"

using namespace std;



namespace SCATMECH {


    const double EulerGamma = 0.5772156649015329;

    // Cosine integral...
    double
    Ci(double x)
    {
        double x2 = sqr(x);
        double x4 = sqr(x2);
        double x6 = sqr(x4);
        double x8 = sqr(x6);

        if (x<1.0) {
            // From Abramowitz and Stegun (5.2.16)
            return EulerGamma - x2/4. + x4/96. - x6/4320. + x8/322560. + log(x);
        } else {
            double g;
            {
                // From Abramowitz and Stegun (5.2.39)
                const double a1 = 42.242855;
                const double a2 = 302.757865;
                const double a3 = 352.018498;
                const double a4 = 21.821899;
                const double b1 = 48.196927;
                const double b2 = 482.485984;
                const double b3 = 1114.978885;
                const double b4 = 449.690326;
                g = (x8+a1*x6+a2*x4+a3*x2+a4)/(x8+b1*x6+b2*x4+b3*x2+b4)/x2;
            }
            double f;
            {
                // From Abramowitz and Stegun (5.2.38)
                const double a1 = 38.027264;
                const double a2 = 265.187033;
                const double a3 = 335.677320;
                const double a4 = 38.102495;
                const double b1 = 40.021433;
                const double b2 = 322.624911;
                const double b3 = 570.236280;
                const double b4 = 157.105423;
                f = (x8+a1*x6+a2*x4+a3*x2+a4)/(x8+b1*x6+b2*x4+b3*x2+b4)/x;
            }
            // From Abramowitz and Stegun (5.2.9)
            return f*sin(x)-g*cos(x);
        }
    }

    void
    RayleighGansSphereScatterer::
    setup()
    {
        SphericalScatterer::setup();
    }

    COMPLEX
    RayleighGansSphereScatterer::
    s1(double angle)
    {
        SETUP();

        double c = cos(angle);
        if (c==1.) { // Check for singularity...
            return COMPLEX(0.,2./3.)*(m-1.)*cube(x);
        }
        if (c==-1.) { // Check for other singularity...
            return COMPLEX(0.,0.25)*(m-1.)*(sin(2.*x)-2.*x*cos(2.*x));
        }
        COMPLEX u = 2.*x*sin(angle/2.);
        COMPLEX factor = COMPLEX(0.,2.)*cube(x/u)*
                         (m-1.)*(sin(u)-u*cos(u));

        return factor;
    }

    COMPLEX
    RayleighGansSphereScatterer::
    s2(double angle)
    {
        return s1(angle)*cos(angle);
    }

    //
    // Routine to calculate scattering cross section
    //
    double
    RayleighGansSphereScatterer::
    Csca()
    {
        SETUP();

        // From van de Hulst, Sec. 7.22:
        double phi = 2.5+2.*sqr(x)-sin(4.*x)/(4.*x)-7./(16*sqr(x))*(1-cos(4*x))
                     + (1./(2.*sqr(x))-2.)*(EulerGamma+log(4*x)-Ci(4*x));
        return  pi*sqr(radius)*norm(m-1.) *phi;
    }

    //
    // Routine to calculate extinction cross section
    //
    double
    RayleighGansSphereScatterer::
    Cext()
    {
        SETUP();

        // From van de Hulst 7.12 ...
        return pi*sqr(radius)*8.*x/3.*imag(m-1.)+Csca();
    }

    //
    // Routine to calculate backscatter cross section
    //
    double
    RayleighGansSphereScatterer::
    Cback()
    {
        SETUP();

        // Backscattering cross section is defined as follows...
        return 4.*pi*MuellerMatrix(s(pi))[0][0]/sqr(k);
    }

    DEFINE_MODEL(RayleighGansSphereScatterer,SphericalScatterer,
                 "Scattering by a sphere in the limit of small index contrast.");


} // namespace SCATMECH



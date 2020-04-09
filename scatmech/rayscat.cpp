//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rayscat.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "scatmech.h"
#include "rayscat.h"

using namespace std;


namespace SCATMECH {


    void
    RayleighScatterer::
    setup()
    {
        SphericalScatterer::setup();
        COMPLEX n2=sqr((COMPLEX)m);
        polarizability = (n2-1.)/(n2+2.)*cube(x);
    }


    COMPLEX
    RayleighScatterer::
    s1(double angle)
    {
        SETUP();
        return polarizability;
    }

    COMPLEX
    RayleighScatterer::
    s2(double angle)
    {
        SETUP();
        return polarizability*cos(angle);
    }

    double
    RayleighScatterer::
    Csca()
    {
        SETUP();
        COMPLEX temp1 = (sqr(m)-1.)/(sqr(m)+2.);

        double Q=  8./3.*sqr(sqr(x))*norm(temp1);
        return pi*sqr(radius)*Q;
    }

    double
    RayleighScatterer::
    Cext()
    {
        SETUP();

        COMPLEX temp1 = (sqr(m)-1.)/(sqr(m)+2.);
        COMPLEX temp2 = (sqr(sqr(m))+27.*sqr(m)+38.)/(2.*sqr(m)+3.);

        double Q = 4.*x*imag(temp1*(1.+sqr(x)/15.*temp1*temp2))+
                   8./3.*sqr(sqr(x))*real(sqr(temp1));
        return pi*sqr(radius)*Q;
    }

    double
    RayleighScatterer::
    Cback()
    {
        SETUP();

        COMPLEX temp1 = (sqr(m)-1.)/(sqr(m)+2.);

        double Q=4.*sqr(sqr(x))*norm(temp1);
        return pi*sqr(radius)*Q;
    }

    DEFINE_MODEL(RayleighScatterer,SphericalScatterer,
                 "Scattering by a sphere in the limit of small size.");



} // namespace SCATMECH


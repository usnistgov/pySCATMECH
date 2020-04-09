//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: fresnel.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "fresnel.h"

using namespace std;


namespace SCATMECH {

    //
    // Fresnel reflection coefficient for s-polarized light
    //
    COMPLEX
    rs(dielectric_constant epsilon, double theta)
    {
        COMPLEX e=epsilon;
        return (cos(theta)-sqrt(e-sqr(sin(theta))))/
               (cos(theta)+sqrt(e-sqr(sin(theta))));
    }

    //
    // Fresnel reflection coefficient for p-polarized light
    //
    COMPLEX
    rp(dielectric_constant epsilon, double theta)
    {
        COMPLEX e=epsilon;
        return (e*cos(theta)-sqrt(e-sqr(sin(theta))))/
               (e*cos(theta)+sqrt(e-sqr(sin(theta))));
    }

    //
    // Fresnel transmission coefficient for s-polarized light
    //
    COMPLEX
    ts(dielectric_constant epsilon, double theta)
    {
        return 2.*cos(theta)/
               (cos(theta)+sqrt((COMPLEX)epsilon-sqr(sin(theta))));
    }

    //
    // Fresnel transmission coefficient for p-polarized light
    //
    COMPLEX
    tp(dielectric_constant epsilon, double theta)
    {
        return 2.*sqrt((COMPLEX)epsilon)*cos(theta)/
               ((COMPLEX)epsilon*cos(theta)+sqrt((COMPLEX)epsilon-sqr(sin(theta))));
    }

    //
    // Fresnel tranmittance (intensity) for s-polarized light
    //
    double
    Ts(dielectric_constant epsilon, double theta)
    {
        return abs(sqrt((COMPLEX)epsilon-sqr(sin(theta)))/
                   cos(theta)*sqr(abs(ts(epsilon,theta))));
    }

    //
    // Fresnel tranmittance (intensity) for p-polarized light
    //
    double
    Tp(dielectric_constant epsilon, double theta)
    {
        return abs(sqrt((COMPLEX)epsilon-sqr(sin(theta)))/
                   cos(theta)*sqr(abs(tp(epsilon,theta))));
    }

    //
    // Snell's law...
    //

    //
    // Changed the following two functions so that the incident angle
    // would be complex.  (Version 3.02, TAG 7 AUG 2002)
    //
    COMPLEX
    sin_internal_angle(optical_constant n,COMPLEX thetai)
    {
        return sin(thetai)/(COMPLEX)n;
    }

    COMPLEX
    cos_internal_angle(optical_constant n,COMPLEX thetai)
    {
        return sqrt(COMPLEX(1,0)-sqr(sin(thetai)/(COMPLEX)(n)));
    }


} // namespace SCATMECH



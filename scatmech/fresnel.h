//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: fresnel.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef FRESNEL_H
#define FRESNEL_H

#include "optconst.h"


namespace SCATMECH {

    //
    // Fresnel reflection coefficients for s- and p-polarized light
    //
    COMPLEX rs( dielectric_constant epsilon, double theta );
    COMPLEX rp( dielectric_constant epsilon, double theta );

    //
    // Fresnel reflectances (intensity) for s- and p-polarized light
    //
    inline double Rs( dielectric_constant epsilon, double theta )
    {
        using namespace std;
        COMPLEX _rs=rs(epsilon,theta);
        return abs(_rs*conj(_rs));
    }

    inline double Rp( dielectric_constant epsilon, double theta )
    {
        using namespace std;
        COMPLEX _rp=rp(epsilon,theta);
        return abs(_rp*conj(_rp));
    }

    //
    // Fresnel transmission coefficients for s- and p-polarized light
    //
    COMPLEX ts( dielectric_constant epsilon, double theta );
    COMPLEX tp( dielectric_constant epsilon, double theta );

    //
    // Fresnel transmission (intensity) for s- and p-polarized light
    //
    double Ts( dielectric_constant epsilon, double thetas );
    double Tp( dielectric_constant epsilon, double thetas );

    //
    // Snell's law of refraction...
    //
    COMPLEX sin_internal_angle( optical_constant n, COMPLEX thetai );
    COMPLEX cos_internal_angle( optical_constant n, COMPLEX thetai );

} // namespace SCATMECH


#endif


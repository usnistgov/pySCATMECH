//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rayinst.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "rayinst.h"
#include "filmtran.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {


    static double subtend(double theta1,double phi1,double theta2, double phi2)
    {
        // This function returns the angle subtended by two vectors defined by
        // (theta1,phi1) and (theta2,phi2)...
        double cosreturnvalue = cos(theta1)*cos(theta2) +
                                cos(phi1)*cos(phi2)*sin(theta1)*sin(theta2) +
                                sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2);
        if (cosreturnvalue>1) return 0.;
        if (cosreturnvalue<-1) return pi;
        else return acos(cosreturnvalue);
    }

    //
    // Mueller Matrix for scattering...
    //
    MuellerMatrix
    Rayleigh_Instrument_BRDF_Model::
    mueller()
    {
        // All jones and mueller routines should call setup() ...
        SETUP();

        throw_backward();
        throw_transmission();

        if (is_reflection()) {
            // This routine uses the derivation found in Thomas A. Germer and Clara C. Asmail,
            // "A goniometric optical scatter instrument for bidirectional reflectance distribution
            // function measurements with out-of-plane and polarimetry capabilities," in
            // Scattering and Surface Roughness, Z.-H. Gu and A. A. Maradudin, Eds. Proc. SPIE
            // 3141, pp. 220-231 (1997).


            double n_minus_1 = air.n(lambda)-1.;

            // Reflection coefficients ...
            complex<double> rss = stack->rs12(thetas,lambda,vacuum,substrate);
            complex<double> rsi = stack->rs12(thetai,lambda,vacuum,substrate);
            complex<double> rps = stack->rp12(thetas,lambda,vacuum,substrate);
            complex<double> rpi = stack->rp12(thetai,lambda,vacuum,substrate);

            // Angle subtended by viewing and scattering directions...
            double subtendi = subtend(-thetai,0.,thetas,phis);
            double sini = fabs(sin(subtendi));
            double subtendr = subtend(thetai,0.,thetas,phis);
            double sinr = fabs(sin(subtendr));

            // Equation 42 of the paper...
            JonesMatrix term1(
                (complex<double>)(sin(thetai)*sin(thetas)-cos(thetas)*cos(thetai)*cos(phis)),
                (complex<double>)cos(phis),
                (complex<double>)(-cos(thetai)*sin(phis)),
                (complex<double>)(-cos(thetas)*sin(phis)));
            // Equation 43 of the paper...
            JonesMatrix term2(
                (-rpi*sin(thetai)*sin(thetas)-rpi*cos(thetas)*cos(thetai)*cos(phis)),
                rsi*cos(phis),
                -rpi*cos(thetai)*sin(phis),
                -rsi*cos(thetas)*sin(phis));
            // Equation 44 of the paper...
            JonesMatrix term3(
                (-rps*sin(thetai)*sin(thetas)-rps*cos(thetas)*cos(thetai)*cos(phis)),
                -rss*cos(phis),
                rss*cos(thetai)*sin(phis),
                -rps*cos(thetas)*sin(phis));
            // Equation 45 of the paper...
            JonesMatrix term4(
                (rps*rpi*sin(thetai)*sin(thetas)-rps*rpi*cos(thetas)*cos(thetai)*cos(phis)),
                -rss*rsi*cos(phis),
                rss*rpi*cos(thetai)*sin(phis),
                -rps*rsi*cos(thetas)*sin(phis));

            // Convert the terms to Mueller matrices for incoherent addition...
            MuellerMatrix sterm1 = term1;
            MuellerMatrix sterm2 = term2;
            MuellerMatrix sterm3 = term3;
            MuellerMatrix sterm4 = term4;

            // Rayleigh scattering terms (Eq. 36)...
            double i0=4.*sqr(pi)/pow(lambda,4.)*sqr(n_minus_1)/number_density;

            // Other geometric terms (Eq. 46)...
            double factor = i0/2./cos(thetas)*field_of_view;
            MuellerMatrix result=((sterm1/sini)+(sterm2/sinr)+(sterm3/sinr)+(sterm4/sini))*factor;
            return result;
        }
        return MuellerZero();
    }

    DEFINE_MODEL(Rayleigh_Instrument_BRDF_Model,Instrument_BRDF_Model,
                 "Effective BRDF due to Rayleigh scatter by the air surrounding a smooth sample.");

    DEFINE_PTRPARAMETER(Rayleigh_Instrument_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_PARAMETER(Rayleigh_Instrument_BRDF_Model,double,field_of_view,"Field of view of detector [um]","1000",0xFF);

    DEFINE_PARAMETER(Rayleigh_Instrument_BRDF_Model,dielectric_function,air,"Index of refraction of air","1.0002748",0xFF);

    DEFINE_PARAMETER(Rayleigh_Instrument_BRDF_Model,double,number_density,"Number density of air [um^-3]","2.51E7",0xFF);



} // namespace SCATMECH


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: onelayer.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "onelayer.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {


    //
    // Routine which returns the scattering matrix for
    // a small defect (Rayleigh approximation) in a
    // dielectric film.
    //
    // The scattering matrix is described by Thomas A. Germer in "Polarized light
    // scattering by microroughness and small defects in dielectric layers,"
    // J. Opt. Soc. Am. A 18(6), 1279-1288 (2001).
    //
    JonesMatrix
    OneLayer_BRDF_Model::
    jonesDSC()
    {
        // All BRDF_Model::jones(,,) should call setup() if recalc is set...
        SETUP();

        throw_backward();
        throw_transmission();

        COMPLEX i(0,1);

        COMPLEX eps1 = substrate.epsilon(lambda);
        COMPLEX eps2 = film.epsilon(lambda);
        COMPLEX epss = defect.epsilon(lambda);

        double k = 2*pi/lambda;
        double ki = k*sin(thetai);
        double kr = k*sin(thetas);

        COMPLEX qi3 = k*cos(thetai);
        COMPLEX qr3 = k*cos(thetas);
        COMPLEX qi2 = sqrt(sqr(k)*eps2-sqr(ki));
        COMPLEX qr2 = sqrt(sqr(k)*eps2-sqr(kr));
        COMPLEX qi1 = sqrt(sqr(k)*eps1-sqr(ki));
        COMPLEX qr1 = sqrt(sqr(k)*eps1-sqr(kr));



        // Changed power of eps2 to 1 - TAG 31 DEC 2002...
        COMPLEX S = 4.*(epss-eps2)/(epss+2.*eps2)*cube(radius)*
                    exp(i*(qi2+qr2)*depth-i*(qi3+qr3)*tau)*eps2*qi3*qr3;

        COMPLEX Ki_plus = exp(2.*i*qi2*tau)+1.;
        COMPLEX Kr_plus = exp(2.*i*qr2*tau)+1.;
        COMPLEX Ki_minus = Ki_plus-2.;
        COMPLEX Kr_minus = Kr_plus-2.;

        COMPLEX Li_plus = exp(2.*i*qi2*(tau-depth))+1.;
        COMPLEX Lr_plus = exp(2.*i*qr2*(tau-depth))+1.;
        COMPLEX Li_minus = Li_plus-2.;
        COMPLEX Lr_minus = Lr_plus-2.;

        COMPLEX F_pi_plus  = eps2*Ki_minus*qi1- eps1*Ki_plus*qi2;
        COMPLEX F_pr_plus  = eps2*Kr_minus*qr1- eps1*Kr_plus*qr2;
        COMPLEX F_pi_minus = eps2*Ki_plus*qi1 - eps1*Ki_minus*qi2;
        COMPLEX F_pr_minus = eps2*Kr_plus*qr1 - eps1*Kr_minus*qr2;

        COMPLEX F_si_plus  = Ki_minus*qi1 - Ki_plus*qi2;
        COMPLEX F_sr_plus  = Kr_minus*qr1 - Kr_plus*qr2;
        COMPLEX F_si_minus = Ki_plus*qi1  - Ki_minus*qi2;
        COMPLEX F_sr_minus = Kr_plus*qr1  - Kr_minus*qr2;

        COMPLEX G_pi_plus  = eps2*Li_minus*qi1- eps1*Li_plus*qi2;
        COMPLEX G_pr_plus  = eps2*Lr_minus*qr1- eps1*Lr_plus*qr2;
        COMPLEX G_pi_minus = eps2*Li_plus*qi1 - eps1*Li_minus*qi2;
        COMPLEX G_pr_minus = eps2*Lr_plus*qr1 - eps1*Lr_minus*qr2;

        COMPLEX G_si_plus  = Li_minus*qi1 - Li_plus*qi2;
        COMPLEX G_sr_plus  = Lr_minus*qr1 - Lr_plus*qr2;
        COMPLEX G_si_minus = Li_plus*qi1  - Li_minus*qi2;
        COMPLEX G_sr_minus = Lr_plus*qr1  - Lr_minus*qr2;

        COMPLEX Gammapi = eps2*F_pi_plus*qi3 - F_pi_minus*qi2;
        COMPLEX Gammapr = eps2*F_pr_plus*qr3 - F_pr_minus*qr2;
        COMPLEX Gammasi =      F_si_plus*qi3 - F_si_minus*qi2;
        COMPLEX Gammasr =      F_sr_plus*qr3 - F_sr_minus*qr2;

        JonesMatrix s;

        s.PP() = (G_pi_plus*G_pr_plus*ki*kr-
                  G_pi_minus*G_pr_minus*qi2*qr2*cos(phis))/Gammapi/Gammapr;
        s.PS() = G_pi_minus*G_sr_plus*k*qi2*sin(phis)/Gammapi/Gammasr;
        s.SP() = G_si_plus*G_pr_minus*k*qr2*sin(phis)/Gammasi/Gammapr;
        s.SS() =  G_si_plus*G_sr_plus*k*k*cos(phis)/Gammasi/Gammasr;

        s = s*S;

        return s;
    }

    DEFINE_MODEL(OneLayer_BRDF_Model,Local_BRDF_Model,
                 "Scattering by a Rayleigh defect in a single dielectric layer.");

    DEFINE_PARAMETER(OneLayer_BRDF_Model,double,radius,"Radius [um]","0.01",0xFF);

    DEFINE_PARAMETER(OneLayer_BRDF_Model,dielectric_function,defect,"Particle","(1,0)",0xFF);

    DEFINE_PARAMETER(OneLayer_BRDF_Model,dielectric_function,film,"Film","(1.59,0)",0xFF);

    DEFINE_PARAMETER(OneLayer_BRDF_Model,double,tau,"Thickness [um]","0.05",0xFF);

    DEFINE_PARAMETER(OneLayer_BRDF_Model,double,depth,"Depth [um]","0",0xFF);


} // namespace SCATMECH



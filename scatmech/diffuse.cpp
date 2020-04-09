//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: diffuse.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "diffuse.h"
#include "fresnel.h"

using namespace std;


namespace SCATMECH {


    static
    COMPLEX
    ArcSin(const COMPLEX& a)
    {
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }

    Diffuse_Subsurface_BRDF_Model::
    Diffuse_Subsurface_BRDF_Model()
    {
    }

    MuellerMatrix Diffuse_Subsurface_BRDF_Model::
    mueller()
    {
        SETUP();

        throw_transmission();
        throw_backward();

        optical_constant vacuum(1.,0.);

        double s = abs(sqrt((COMPLEX)substrate.epsilon(lambda)-sqr(sin(thetai)))/cos(thetai)
                       * sqrt((COMPLEX)substrate.epsilon(lambda)-sqr(sin(thetas)))/cos(thetas));

        // The following line was added for Version 4.00...
        s *= feedback/sqr(substrate.n(lambda));

        MuellerMatrix m = Lambertian_BRDF_Model::mueller();

        return s*(((MuellerMatrix)(stack->t12(thetas,lambda,vacuum,substrate))*m)*
                  (MuellerMatrix)(stack->t12(thetai,lambda,vacuum,substrate)));
    }

    void
    Diffuse_Subsurface_BRDF_Model::
    setup()
    {
        Lambertian_BRDF_Model::setup();

        // The following code was added for Version 4.00...

        substrate.force_nonabsorbing();

        double rho=0;
        double normalization=0;
        for (double theta=0; theta<=pi/2; theta+=deg/5.) {
            // Need angle outside of material...
            COMPLEX ti = ArcSin((COMPLEX)(sin(theta)*substrate.n(lambda)));

            // Calculate reflectance...
            JonesMatrix r = stack->r21(ti,lambda,substrate,vacuum);
            double rr=0.5*(norm(r[0])+norm(r[1]));

            // Integrate the reflectance...
            double dOmega = sin(2.*theta);
            rho += rr*dOmega;
            normalization += dOmega;
        }
        rho = rho/normalization;

        // The amount of light which escapes the material is
        // increased by this factor...
        double R = reflectance->Get_Reflectance(lambda);
        feedback = 1./(1.-rho*R);
    }

    DEFINE_MODEL(Diffuse_Subsurface_BRDF_Model,Lambertian_BRDF_Model,
                 "Totally diffuse scattering under a smooth surface.");

    DEFINE_PTRPARAMETER(Diffuse_Subsurface_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);


} // namespace SCATMECH





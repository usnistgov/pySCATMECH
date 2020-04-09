//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: transmit.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "transmit.h"
#include "filmtran.h"

using namespace std;

namespace SCATMECH {

    MuellerMatrix
    Transmit_BRDF_Model::
    mueller()
    {
        if (lambda!=model->get_lambda()) error("lambda!=model.lambda");
        if ((type)!=(model->get_type())) error("type!=model.type");
        if (substrate.index(lambda)!=model->get_substrate().index(lambda)) error("substrate!=model.substrate");

        SETUP();

        if (type!=0)
            if (substrate.k(lambda)!=0) error("Substrate must be non-absorbing for type " + to_string(type));

        switch (type) {
            case 0:
                return model->Mueller(thetai,thetas,phis,rotation);
                break;
            case 1:
            {
                double thetas_internal = asin(sin(thetas)/substrate.n(lambda));

                MuellerMatrix T = films->t12(thetas,lambda,vacuum,substrate);
                double transmittancefactorout = substrate.n(lambda)*cos(thetas_internal)/cos(thetas);
                double jacobian = cos(thetas)/cos(thetas_internal)/sqr(substrate.n(lambda));
                double ratiocosines = cos(thetas_internal)/cos(thetas);

                T = T*(jacobian*transmittancefactorout*ratiocosines);

                return T * model->Mueller(thetai,thetas_internal,phis,rotation);
            }
            break;
            case 2:
            {
                double thetas_internal = asin(sin(thetas)/substrate.n(lambda));
                double thetai_internal = asin(sin(thetai)/substrate.n(lambda));
                MuellerMatrix Ti = films->t12(thetai,lambda,vacuum,substrate);
                MuellerMatrix Ts = films->t12(thetas,lambda,vacuum,substrate);
                double transmittancefactorin = cos(thetai_internal)/cos(thetai)*substrate.n(lambda);
                double transmittancefactorout = cos(thetas_internal)/cos(thetas)*substrate.n(lambda);
                double jacobian = cos(thetas)/cos(thetas_internal)/sqr(substrate.n(lambda));
                double ratiocosines = cos(thetas_internal)/cos(thetas);
                double factor = transmittancefactorin*transmittancefactorout*jacobian*ratiocosines;

                MuellerMatrix internalBRDF = model->Mueller(thetai_internal,thetas_internal,phis,rotation);
                return Ts*internalBRDF*Ti*factor;
            }
            break;
            case 3:
            {
                double thetai_internal = asin(sin(thetai)/substrate.n(lambda));
                MuellerMatrix Ti = films->t12(thetai,lambda,vacuum,substrate);
                double transmittancefactorin = cos(thetai_internal)/cos(thetai)*substrate.n(lambda);
                double factor = transmittancefactorin;

                MuellerMatrix internalBRDF = model->Mueller(thetai_internal,thetas,phis,rotation);
                return internalBRDF*Ti*factor;
            }
            break;
            default:
                error("Invalid type = " + to_string(type));
        }
        return MuellerZero();
    }

    DEFINE_MODEL(Transmit_BRDF_Model,
                 BRDF_Model,
                 "Model which evaluates another BRDF_Model in transmission, outside of a material.");

    DEFINE_PTRPARAMETER(Transmit_BRDF_Model,BRDF_Model_Ptr,model,"The model","Microroughness_BRDF_Model",0xFF);
    DEFINE_PTRPARAMETER(Transmit_BRDF_Model,StackModel_Ptr,films,"Films on the bottom side","No_StackModel",0xFF);

} // namespace SCATMECH

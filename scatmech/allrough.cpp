//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: allrough.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "allrough.h"

using namespace std;

namespace SCATMECH {

    void
    Correlated_Roughness_Stack_BRDF_Model::
    setup()
    {
        Roughness_BRDF_Model::setup();

        model.set_lambda(lambda);
        model.set_type(type);
        model.set_substrate(substrate);
        model.set_psd(psd);
        model.set_stack(stack);
    }

    JonesMatrix
    Correlated_Roughness_Stack_BRDF_Model::
    jones()
    {
        SETUP();

        JonesMatrix J=JonesZero();

        for (int layer=0; layer<=stack->get_n(); ++layer) {
            model.set_this_layer(layer);
            J = J + model.Jones(thetai,thetas,phis,rotation);
        }
        return J;
    }

    void
    Uncorrelated_Roughness_Stack_BRDF_Model::
    setup()
    {
        Roughness_BRDF_Model::setup();

        model.set_lambda(lambda);
        model.set_type(type);
        model.set_substrate(substrate);
        model.set_psd(psd);
        model.set_stack(stack);
    }

    MuellerMatrix
    Uncorrelated_Roughness_Stack_BRDF_Model::
    mueller()
    {
        SETUP();

        MuellerMatrix M=MuellerZero();
        for (int layer=0; layer<=stack->get_n(); ++layer) {
            model.set_this_layer(layer);
            M = M + model.Mueller(thetai,thetas,phis,rotation);
        }
        return M;
    }

    void
    Growth_Roughness_Stack_BRDF_Model::
    setup()
    {
        Roughness_BRDF_Model::setup();

        model.set_lambda(lambda);
        model.set_type(type);
        model.set_substrate(substrate);
        model.set_psd(psd);
        model.set_stack(stack);
    }
    //
    // Growth_Roughness_Stack_BRDF_Model uses the film growth model described in
    // E. Spiller, D. Stearns, M. Krumrey, J. Appl. Phys. 74(1), 107-118 (1993),
    // with some modifications:
    // (1) the intrinsic infinite-thickness roughness of all
    // coating interfaces are the same. This is treated as such because the
    // software would be too complicated, as there is no high index, low index assumptions
    // in dielectric_stack.
    // (2) the replication factor has an extra parameter, the spatial frequency exponent, which
    // described in the paper is 2.
    //
    // The substrate PSD function is given by the Roughness_Stack_BRDF_Model::psd, while
    // the intrinsic infinite-thickness roughness is given by intrinsic.
    //
    MuellerMatrix
    Growth_Roughness_Stack_BRDF_Model::
    mueller()
    {
        SETUP();

        int N = stack->get_n();
        double fx,fy;
        Bragg_Frequency(fx,fy);
        double qpow = pow(4*sqr(pi)*(sqr(fx)+sqr(fy)),exponent/2.);
        double nu = pow(relaxation,exponent-1.);
        double fourpisqr = 4.*sqr(pi);
        double PSDsub = psd->psd(fx,fy);
        double FTsub = sqrt(PSDsub);
        double PSDint0 = intrinsic->psd(fx,fy);
        double FTint0 = sqrt(PSDint0);

        vector<double> a(N);
        vector<double> PSDint(N);
        for (int i=0; i<N; ++i) {
			// The following was changed 2/27/2019...
			// double argument = fourpisqr*nu*qpow*stack->get_t()[i];
			double argument = nu*qpow*stack->get_t()[i];

            // The replication factor between the layers...
            a[i] = exp(-argument);
            // The intrinsic roughness of the layer...
            PSDint[i] = PSDint0*(1.-sqr(a[i]));
        }

        vector<JonesMatrix> J(N+1);
        for (int i=0; i<=N; ++i) {
            model.set_this_layer(i);
            // The Jones matrix for scattering by the layer, per unit surface roughness...
            J[i] = model.Jones(thetai,thetas,phis,rotation)/FTsub;
        }

        // Cross power spectral density matrix...
        vector<vector<double> > covar(N+1,vector<double>(N+1));
        // First the diagonal elements (the power spectra)...
        for (int i=0; i<=N; ++i) {
            covar[i][i] = 0.;
            for (int j=0; j<=i; ++j) {
                double aprod=1;
                for (int k=j; k<i; ++k) {
                    aprod *= sqr(a[k]);
                }
                covar[i][i] += aprod*(j==0 ? PSDsub : PSDint[j-1]);
            }
        }
        // Then the off diagonal elements...
        for (int i=0; i<=N; ++i) {
            for (int j=i+1; j<=N; ++j) {
                covar[j][i] = covar[i][j] = a[j-1]*covar[i][j-1];
            }
        }

        // Calculate the Mueller matrix from the covariance matrix...
        MuellerMatrix result = MuellerZero();
        for (int i=0; i<=N; ++i) {
            for (int j=0; j<=N; ++j) {
                result += covar[i][j]*ReCrossMueller(J[i],J[j]);
            }
        }
        return result;
    }

    DEFINE_MODEL(Correlated_Roughness_Stack_BRDF_Model,Roughness_BRDF_Model,
                 "All films in a stack equally rough and correlated");
	DEFINE_PTRPARAMETER(Correlated_Roughness_Stack_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);


    DEFINE_MODEL(Uncorrelated_Roughness_Stack_BRDF_Model,Roughness_BRDF_Model,
                 "All films in a stack equally rough but uncorrelated");
    DEFINE_PTRPARAMETER(Uncorrelated_Roughness_Stack_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_MODEL(Growth_Roughness_Stack_BRDF_Model,Roughness_BRDF_Model,
                 "All films in a stack having intrinsic roughness, but also correlation with roughness from the stack below.");
    DEFINE_PTRPARAMETER(Growth_Roughness_Stack_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);
    DEFINE_PARAMETER(Growth_Roughness_Stack_BRDF_Model,double,relaxation,"Correlation relaxation parameter [um]","0",0xFF);
    DEFINE_PARAMETER(Growth_Roughness_Stack_BRDF_Model,double,exponent,"Correlation spatial frequency exponent parameter","2",0xFF);
    DEFINE_PTRPARAMETER(Growth_Roughness_Stack_BRDF_Model,PSD_Function_Ptr,intrinsic,"Intrinsic PSD of the film interfaces","ABC_PSD_Function",0xFF);


} // namespace SCATMECH;


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: roughnes.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "roughnes.h"
#include "allrough.h"
#include <vector>

using namespace std;


namespace SCATMECH {


    //
    // Jones scattering matrix for scattering from a single rough interface...
    //
    JonesMatrix
    Roughness_Stack_BRDF_Model::
    jones()
    {
        double theta0=thetai;
        double theta=thetas;
        double phi=phis;
        // This program follows the derivation outlined in:
        //      J. M. Elson, "Multilayer-coated optics: guided-wave coupling
        //                    and scattering by means of interface random
        //                    roughness," J. Opt. Soc. Am. A 12(4) 729 (1995).
        //		J. M. Elson, "Theory and Software for Light Scattering From
        //					  Multilayer Optical Components with Interfacial
        //				      Roughness," Naval Air Warfare Center Weapons
        //                    Division (NAWCWPNS), Tech. Pub. 8084 (1992).
        //
        // Define variables with names in line with those used by Elson.
        // Note that I use C/C++ conventions and start counting at zero,
        // instead of 1.
        //
        // Note use of class JonesMatrix for doing 2x2 matrix calculations.
        // These matrices have nothing to do with the Jones formalism, but
        // rather with transfer matrices, and are used only for convenience.
        //
        // Note that I assume omega_0/c = 1 and I(ntensity) = 1.
        // (Use I to be sqrt(-1).)
        //
        SETUP();

        const int M11=1,M12=2,M21=3,M22=0;

        int L=stack->get_n();
        int m = is_forward() ? this_layer : L-this_layer;
        int n;

        if (this_layer<0) error("this_layer<0");
        if (this_layer>L) error("this_layer>number of layers");

        COMPLEX I(0,1);
        vector<double> d(L+1);
        vector<COMPLEX> eps(L+2);
        vector<COMPLEX> q0(L+2);
        vector<COMPLEX> q1(L+2);

        //
        // The parameter theta is defined differently for
        // reflection than for transmission...
        //
        double sintheta = sin(theta);
        COMPLEX scatter_medium = is_reflection() ? COMPLEX(1.,0.) : (COMPLEX)(substrate.index(lambda));
        COMPLEX k = sin(theta)*scatter_medium;

        double k0 = sin(theta0);

        double _lambda = lambda;
        //
        // Define variables as per Figure 1 of paper.
        //
        if (is_forward()) {
            eps[0]=(COMPLEX)(substrate.epsilon(lambda));
            for (n=1; n<L+1; ++n) {
                eps[n]=(COMPLEX)(stack->get_e()[n-1].epsilon(lambda));
            }
            eps[L+1]=1;

            d[0]=0;
            for (n=1; n<L+1; ++n) d[n]=d[n-1]+2*pi*stack->get_t()[n-1]/lambda;
        } else { // is_backward()
            _lambda = lambda/substrate.n(lambda);
            scatter_medium = is_reflection() ? 1. : 1/substrate.n(lambda);
            k = sin(theta)*scatter_medium;

            eps[0]=1./(COMPLEX)(substrate.epsilon(lambda));
            for (n=1; n<L+1; ++n) {
                eps[n]=(COMPLEX)(stack->get_e()[L-n].epsilon(lambda))/(COMPLEX)(substrate.epsilon(lambda));
            }
            eps[L+1]=1.;

            d[0]=0;
            for (n=1; n<L+1; ++n) d[n]=d[n-1]+2*pi*stack->get_t()[L-n]/_lambda;
        }

        //
        // Symbols defined in the paper:
        // alpha_i, q^(0)_i, q^(1)_i, beta_i, mu_i, p_i
        //
        for (n=0; n<L+2; ++n) {
            q0[n]=sqrt(eps[n]-k0*k0);
            q1[n]=sqrt(eps[n]-k*k);
        }

        //
        // We will be needing a bunch of different matrix products:
        //
        JonesMatrix unitm(1.,1.,0.,0.);
        JonesMatrix P01L1 = unitm, S01L1 = unitm, P01m1 = unitm, S01m1 = unitm,
                    P11L1 = unitm, S11L1 = unitm, P11m1 = unitm, S11m1 = unitm;

        int j;
        for (j=0; j<L+1; ++j) {
            COMPLEX denom0 = 2.*q0[j]*sqrt(eps[j]*eps[j+1]);
            JonesMatrix p0(
                (eps[j]*q0[j+1]+eps[j+1]*q0[j])/denom0*exp(-I*(q0[j+1]-q0[j])*d[j]),
                (eps[j]*q0[j+1]+eps[j+1]*q0[j])/denom0*exp(+I*(q0[j+1]-q0[j])*d[j]),
                (eps[j]*q0[j+1]-eps[j+1]*q0[j])/denom0*exp(-I*(q0[j+1]+q0[j])*d[j]),
                (eps[j]*q0[j+1]-eps[j+1]*q0[j])/denom0*exp(+I*(q0[j+1]+q0[j])*d[j]));
            JonesMatrix s0(
                (q0[j]+q0[j+1])/(2.*q0[j])*exp(-I*(q0[j+1]-q0[j])*d[j]),
                (q0[j]+q0[j+1])/(2.*q0[j])*exp(+I*(q0[j+1]-q0[j])*d[j]),
                (q0[j]-q0[j+1])/(2.*q0[j])*exp(-I*(q0[j+1]+q0[j])*d[j]),
                (q0[j]-q0[j+1])/(2.*q0[j])*exp(+I*(q0[j+1]+q0[j])*d[j]));
            COMPLEX denom1 = 2.*q1[j]*sqrt(eps[j]*eps[j+1]);
            JonesMatrix p1(
                (eps[j]*q1[j+1]+eps[j+1]*q1[j])/denom1*exp(-I*(q1[j+1]-q1[j])*d[j]),
                (eps[j]*q1[j+1]+eps[j+1]*q1[j])/denom1*exp(+I*(q1[j+1]-q1[j])*d[j]),
                (eps[j]*q1[j+1]-eps[j+1]*q1[j])/denom1*exp(-I*(q1[j+1]+q1[j])*d[j]),
                (eps[j]*q1[j+1]-eps[j+1]*q1[j])/denom1*exp(+I*(q1[j+1]+q1[j])*d[j]));
            JonesMatrix s1(
                (q1[j]+q1[j+1])/(2.*q1[j])*exp(-I*(q1[j+1]-q1[j])*d[j]),
                (q1[j]+q1[j+1])/(2.*q1[j])*exp(+I*(q1[j+1]-q1[j])*d[j]),
                (q1[j]-q1[j+1])/(2.*q1[j])*exp(-I*(q1[j+1]+q1[j])*d[j]),
                (q1[j]-q1[j+1])/(2.*q1[j])*exp(+I*(q1[j+1]+q1[j])*d[j]));

            P01L1 *= p0;
            S01L1 *= s0;
            P11L1 *= p1;
            S11L1 *= s1;
            if (j<m) {
                P01m1 *= p0;
                S01m1 *= s0;
                P11m1 *= p1;
                S11m1 *= s1;
            }
        }

        COMPLEX am0 = -q0[L+1]*P01m1[M12]/(q0[m]*P01L1[M11]);
        COMPLEX bm0 =  q0[L+1]*P01m1[M11]/(q0[m]*P01L1[M11]);
        COMPLEX gm0 = -q0[L+1]*S01m1[M12]/(q0[m]*S01L1[M11]);
        COMPLEX fm0 =  q0[L+1]*S01m1[M11]/(q0[m]*S01L1[M11]);

        COMPLEX e0mx = q0[m]/sqrt(eps[m])*(am0*exp(I*q0[m]*d[m])
                                           +bm0*exp(-I*q0[m]*d[m]));
        COMPLEX e0mz = -k0/sqrt(eps[m])*  (am0*exp(I*q0[m]*d[m])
                                           -bm0*exp(-I*q0[m]*d[m]));
        COMPLEX e0my = (gm0*exp(I*q0[m]*d[m])+fm0*exp(-I*q0[m]*d[m]));

        COMPLEX eta1 = P11m1[M11]*exp(-I*q1[m]*d[m])+P11m1[M12]*exp(I*q1[m]*d[m]);
        COMPLEX eta2 = P11m1[M11]*exp(-I*q1[m]*d[m])-P11m1[M12]*exp(I*q1[m]*d[m]);

        COMPLEX xi1  = S11m1[M11]*exp(-I*q1[m]*d[m])-S11m1[M12]*exp(I*q1[m]*d[m]);
        COMPLEX xi2  = S11m1[M21]*exp(-I*q1[m]*d[m])-S11m1[M22]*exp(I*q1[m]*d[m]);

        COMPLEX Xs  = (eps[m+1]-eps[m])*(-e0my*cos(phi));
        COMPLEX Xp  = (eps[m+1]-eps[m])*(e0mx*sin(phi));
        COMPLEX Y1p = (eps[m+1]-eps[m])/eps[m+1]*e0mz;
        COMPLEX Y1s = 0;
        COMPLEX Y2s = (eps[m+1]-eps[m])*(e0my*sin(phi));
        COMPLEX Y2p = (eps[m+1]-eps[m])*(e0mx*cos(phi));

        JonesMatrix J;

        if (is_reflection()) {

            COMPLEX mu_pp = (k*eps[m]*Y1p*eta1-q1[m]*Y2p*eta2)/
                            (q1[m]*sqrt(eps[m]));
            COMPLEX mu_sp = Xp*xi1/q1[m];
            COMPLEX mu_ps = (k*eps[m]*Y1s*eta1-q1[m]*Y2s*eta2)/
                            (q1[m]*sqrt(eps[m]));
            COMPLEX mu_ss = Xs*xi1/q1[m];

            double fx,fy;
            Bragg_Frequency(fx,fy);

            COMPLEX prefactor = sqrt(sqr(pi)/pow(_lambda,4.)*
                                     q1[L+1]*eps[L+1]/q0[L+1] * psd->psd(fx,fy));

            J.PP() =  mu_pp/P11L1[M11];
            J.SS() =  -mu_ss/S11L1[M11];
            J.SP() =  (mu_ps)/P11L1[M11];
            J.PS() =  (-mu_sp)/S11L1[M11];

            J = J*prefactor;

        } else { // is_transmission()
            COMPLEX eta3 = P11m1[M21]*exp(-I*q1[m]*d[m])+P11m1[M22]*exp(I*q1[m]*d[m]);
            COMPLEX eta4 = P11m1[M21]*exp(-I*q1[m]*d[m])-P11m1[M22]*exp(I*q1[m]*d[m]);

            COMPLEX alpha_pp = -k*sqrt(eps[m])*Y1p/q1[m]*(eta1*P11L1[M21]-eta3*P11L1[M11])
                               + Y2p/sqrt(eps[m])*(eta2*P11L1[M21]-eta4*P11L1[M11]);
            COMPLEX alpha_sp =  Xp/q1[m]*(xi2*S11L1[M11]-xi1*S11L1[M21]);
            COMPLEX alpha_ps = k*sqrt(eps[m])*Y1s/q1[m]*(eta1*P11L1[M21]-eta3*P11L1[M11])
                               + Y2s/sqrt(eps[m])*(eta2*P11L1[M21]-eta4*P11L1[M11]);
            COMPLEX alpha_ss = Xs/q1[m]*(xi2*S11L1[M11]-xi1*S11L1[M21]);

            double fx,fy;
            Bragg_Frequency(fx,fy);

            COMPLEX prefactor = sqrt(
                                    pow(2.*pi/_lambda,4.)*q1[0]*eps[0]/(16.*sqr(pi)*q0[L+1])
                                    * psd->psd(fx,fy)
                                );

            J.PP() =  alpha_pp/P11L1[M11];
            J.SS() =  alpha_ss/S11L1[M11];
            J.SP() =  alpha_ps/P11L1[M11];
            J.PS() =  alpha_sp/S11L1[M11];

            J = J*prefactor;
        }

        if (is_forward()) {
            return J;
        } else {
            return JonesMatrix(J[0],J[1],-J[2],-J[3]);
        }
    }

    void Register(const Roughness_Stack_BRDF_Model* x)
    {
        static bool Models_Registered = false;

        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Roughness_Stack_BRDF_Model);
            Register_Model(Correlated_Roughness_Stack_BRDF_Model);
            Register_Model(Uncorrelated_Roughness_Stack_BRDF_Model);
            Register_Model(Growth_Roughness_Stack_BRDF_Model);
        }
    }

    DEFINE_MODEL(Roughness_Stack_BRDF_Model,Roughness_BRDF_Model,
                 "Scattering by a single rough interface in a stack of films.");

    DEFINE_PTRPARAMETER(Roughness_Stack_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_PARAMETER(Roughness_Stack_BRDF_Model,int,this_layer,"Rough interface","0",0xFF);


} // namespace SCATMECH


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: urough.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include <string>
#include "scatmech.h"
#include "urough.h"

using namespace std;


namespace SCATMECH {


    //
    // Jones Matrix for scattering...
    //
    JonesMatrix
    Microroughness_BRDF_Model::
    jones()
    {
        // All routines should call setup()...
        SETUP();

        if (is_reflection()) {
            // This routine follows the Rayleigh-Rice first-order vector
            // perturbation theory. For details, see Thomas A. Germer,
            // "Angular dependence and polarization of out-of-plane optical
            // scattering from particulate contamination, subsurface defects,
            // and surface microroughness,"
            // Applied Optics vol. 36 no. 33, pp 8798-8805 (1997).

            // An abreviation...
            COMPLEX e = (COMPLEX)(substrate.epsilon(lambda));
            double _lambda = lambda;

            if (is_backward()) {
                e = 1./e;
                _lambda = _lambda/substrate.n(lambda);
            }

            // p-->p scattering
            COMPLEX pp = (e-1.)*
                         (
                             e*sin(thetai)*sin(thetas)-
                             sqrt(e-sqr(sin(thetas)))*sqrt(e-sqr(sin(thetai)))*cos(phis)
                         )/
                         (e*cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (e*cos(thetas)+sqrt(e-sqr(sin(thetas))));

            // s-->s scattering
            COMPLEX ss = (e-1.)*cos(phis)/
                         (cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (cos(thetas)+sqrt(e-sqr(sin(thetas))));

            // p-->s scattering
            COMPLEX ps = -(e-1.)*sqrt(e-sqr(sin(thetai)))*sin(phis)/
                         (e*cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (cos(thetas)+sqrt(e-sqr(sin(thetas))));
            // s-->p scattering
            COMPLEX sp = -(e-1.)*sqrt(e-sqr(sin(thetas)))*sin(phis)/
                         (cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (e*cos(thetas)+sqrt(e-sqr(sin(thetas))));

            // Calculate spatial frequency...
            double fx,fy;
            Bragg_Frequency(fx,fy);

            // Various factors common to all the scattering elements...
            COMPLEX factor=sqrt(
                               16.*sqr(pi)/sqr(sqr(_lambda))*
                               cos(thetai)*cos(thetas)*
                               psd->psd(fx,fy)
                           );

            // Return the Jones matrix...
            JonesMatrix result = JonesMatrix(pp,ss,ps,sp)*factor;
            if (is_forward()) {
                return result;
            } else {
                return JonesMatrix(result[0],result[1],-result[2],-result[3]);
            }

        } else { /* is_transmission() */

            // The solution for transmissive scatter was obtained from
            //      J. M. Elson, "Multilayer-coated optics: guided-wave coupling
            //                    and scattering by means of interface random
            //                    roughness," J. Opt. Soc. Am. A 12(4) 729 (1995).
            //
            // The case for L=0 simplifies considerably to the following...
            //
            COMPLEX q0[2];
            COMPLEX q1[2];

            COMPLEX eps=(COMPLEX)(substrate.epsilon(lambda));
            double _lambda = lambda;

            if (is_backward()) {
                eps = 1./eps;
                _lambda = lambda/substrate.n(lambda);
            }

            COMPLEX index = sqrt(eps);
            COMPLEX k = sin(thetas)*index;
            double k0 = sin(thetai);

            q0[0]=sqrt(eps-k0*k0);
            q1[0]=sqrt(eps-k*k);
            q0[1]=sqrt(1.-k0*k0);
            q1[1]=sqrt(1.-k*k);

            double fx,fy;
            Bragg_Frequency(fx,fy);

            double PSD = psd->psd(fx,fy);

            JonesMatrix J;
            J.PP() =  index*(k*k0+cos(phis)*q0[0]*q1[1])/(q0[0]+eps*q0[1])/(q1[0]+eps*q1[1]);
            J.SS() =  cos(phis)/(q0[0]+q0[1])/(q1[0]+q1[1]);
            J.SP() =  index*q1[1]*sin(phis)/(q0[0]+q0[1])/(q1[0]+eps*q1[1]);
            J.PS() =  -q0[0]*sin(phis)/(q0[0]+eps*q0[1])/(q1[0]+q1[1]);

            COMPLEX common = 4.*pi*(1.-eps)*index*sqrt(q1[0]*q0[1]*PSD)/sqr(_lambda);

            if (is_forward()) {
                return J*common;
            } else {
                return JonesMatrix(J[0],J[1],-J[2],-J[3])*common;
            }
        }
    }

    DEFINE_MODEL(Microroughness_BRDF_Model,
                 Roughness_BRDF_Model,
                 "Scattering by a single rough surface in the smooth-surface limit (no films).");


} // namespace SCATMECH


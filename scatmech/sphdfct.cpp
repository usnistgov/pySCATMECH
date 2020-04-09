//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sphdfct.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "sphdfct.h"
#include "vector3d.h"
#include "matrix3d.h"
#include "fresnel.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {


    static
    COMPLEX
    ArcSin(const COMPLEX& a)
    {
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }

    //
    // Routine for calculating jones matrix from arbitrary spherical defect...
    //
    JonesMatrix
    Rayleigh_Defect_BRDF_Model::
    jonesDSC()
    {
        SETUP();

        if (is_forward()) {
            if (is_reflection()) {
                // This routine uses the methodology outlined in Thomas A. Germer,
                // "Angular dependence and polarization of out-of-plane optical
                // scattering from particulate contamination, subsurface defects,
                // and surface microroughness," Appl. Opt. 36, 8798--8805 (1997).
                //
                // The code follows the derivation, except that the transmission
                // coefficients are chosen to be appropriate for a film.
                //
                // Also, there was an error in the paper regarding the Eq. (8)
                // and Appendix A.
                //

                COMPLEX n = substrate.index(lambda);
                COMPLEX e = substrate.epsilon(lambda);
                COMPLEX e_part=defect.epsilon(lambda);
                dielectric_constant epsilon=substrate.epsilon(lambda);

                COMPLEX x = kr*n;
                double k=2.*pi/lambda;

                // Angles of refraction internal to the material...
                COMPLEX sin_thetai_internal = sin(thetai)/n; //sin_internal_angle(epsilon,thetai);
                COMPLEX cos_thetai_internal = sqrt(1.0-sqr(sin_thetai_internal)); //cos_internal_angle(epsilon,thetai);
                COMPLEX sin_thetas_internal = sin(thetas)/n; //sin_internal_angle(epsilon,thetas);
                COMPLEX cos_thetas_internal = sqrt(1.0-sqr(sin_thetas_internal)); //cos_internal_angle(epsilon,thetas);

                // The incident and scattering khat vectors outside the material...
                CVector inKo((COMPLEX)(sin(thetai)),(COMPLEX)(0.),(COMPLEX)(-cos(thetai)));
                CVector outKo((COMPLEX)(sin(thetas)*cos(phis)),
                              (COMPLEX)(sin(thetas)*sin(phis)),(COMPLEX)(cos(thetas)));

                // The incident and scattering khat vectors inside the material...
                // (one of my compilers not have a COMPLEX sin(COMPLEX)!)
                CVector inKi(sin_thetai_internal,(COMPLEX)(0.),-cos_thetai_internal);
                CVector outKi(sin_thetas_internal*cos(phis),
                              sin_thetas_internal*sin(phis),
                              cos_thetas_internal);

                // The surface normal...
                Vector normal(0.,0.,-1.);

                // The following definitions make {s,p,k} right handed...
                CVector cnormal = normal;
                CVector inSi =  perpto(inKi,cnormal);
                CVector inPi =  perpto(inKi,inSi);
                CVector outSi = perpto(outKi,cnormal);
                CVector outPi = perpto(outKi,outSi);
                CVector inSo =  perpto(inKo,cnormal);
                CVector inPo =  perpto(inKo,inSo);
                CVector outSo = perpto(outKo,cnormal);
                CVector outPo = perpto(outKo,outSo);

                // Transmission coefficients...
                const dielectric_constant vacuum(1.);

                COMPLEX tsi = stack->ts12(thetai,lambda,vacuum,substrate);
                COMPLEX tpi = stack->tp12(thetai,lambda,vacuum,substrate);
                COMPLEX tss = stack->ts21(thetas,lambda,substrate,vacuum);
                COMPLEX tps = stack->tp21(thetas,lambda,substrate,vacuum);

                // Transmission transfer matrices...
                CMatrix ti = (tsi*outer(inSi,inSo))+(tpi*outer(inPi,inPo));
                // Eq. 8 of the paper is corrected here...
                // Factor of n corrected...TAG 31 DEC 2002...
                CMatrix ts = cos(thetas)/cos_thetas_internal/n*
                             ((tss*outer(outSo,outSi))+(tps*outer(outPo,outPi)));

                COMPLEX rayleigh = (e_part-e)/(e_part+2.*e)*cube(x);
                CMatrix scatterlocal = (outer(outPi,outPi)+outer(outSi,outSi))*rayleigh;

                // Loss and phase factors...
                COMPLEX I(0,1);
                // Error found in following two lines in SCATMECH Version 1.0
                // COMPLEX alpha = exp(2.*I*kd*cos(thetai));
                // COMPLEX beta  = exp(2.*I*kd*cos(thetas));
                // Replaced with following two lines (TAG 11 April 2000):
                COMPLEX alpha = exp(I*n*kd*cos_thetai_internal);
                COMPLEX beta  = exp(I*n*kd*cos_thetas_internal);

                // The scatter matrix globally includes transmission to/from defect...
                CMatrix scatterglobal = ts*(scatterlocal*ti)*alpha*beta;

                COMPLEX pp = outPo*(scatterglobal*inPo);
                COMPLEX sp = outPo*(scatterglobal*inSo);
                COMPLEX ps = outSo*(scatterglobal*inPo);
                COMPLEX ss = outSo*(scatterglobal*inSo);

                COMPLEX common=1./(k*n);
                // Return the whole value with factors common to all elements...
                return JonesMatrix(pp,ss,ps,sp)*common;

            } else { // is_transmission()

                COMPLEX n = substrate.index(lambda);
                COMPLEX e = substrate.epsilon(lambda);
                COMPLEX e_part=defect.epsilon(lambda);
                dielectric_constant epsilon=substrate.epsilon(lambda);

                COMPLEX x = kr*n;
                double k=2.*pi/lambda;

                COMPLEX sin_thetai_internal = sin(thetai)/n;
                COMPLEX cos_thetai_internal = sqrt(1.-sqr(sin_thetai_internal));

                // The incident and scattering khat vectors at the defect...
                CVector inK(sin_thetai_internal,0.,-cos_thetai_internal);
                CVector outK(sin(thetas)*cos(phis),sin(thetas)*sin(phis),-cos(thetas));
                CVector outKr(outK.x,outK.y,-outK.z);

                // The surface normal...
                CVector normal(0.,0.,-1.);

                // The following definitions make {s,p,k} right handed...
                CVector inS =  perpto(inK,normal);
                CVector inP =  perpto(inK,inS);
                CVector outS =  perpto(outK,normal);
                CVector outP =  perpto(outK,outS);
                CVector outSr = outS;
                CVector outPr = perpto(outKr,outSr);

                // Reflection and transmission coefficients...

                COMPLEX ts = stack->ts12(thetai,lambda,vacuum,substrate)*sqrt(n);
                COMPLEX tp = stack->tp12(thetai,lambda,vacuum,substrate)*sqrt(n);
                COMPLEX rs = stack->rs21i(fabs(thetas),lambda,substrate,vacuum);
                COMPLEX rp = stack->rp21i(fabs(thetas),lambda,substrate,vacuum);

                COMPLEX rayleigh = (e_part-e)/(e_part+2.*e)*cube(x);

                // Loss and phase factors...
                COMPLEX I(0,1);
                COMPLEX alpha = exp(2.*I*n*kd*cos(thetas));

                COMPLEX pp = tp*(outP*inP + outPr*inP*(alpha*rp));
                COMPLEX sp = ts*(outP*inS + outPr*inS*(alpha*rp));
                COMPLEX ps = tp*(outS*inP + outSr*inP*(alpha*rs));
                COMPLEX ss = ts*(outS*inS + outSr*inS*(alpha*rs));

                COMPLEX common=rayleigh/(k*n);
                // Return the whole value with factors common to all elements...
                return JonesMatrix(pp,ss,ps,sp)*common;
            }
        } else { // is_backward()
            if (is_reflection()) {

                COMPLEX n = substrate.index(lambda);
                COMPLEX e = substrate.epsilon(lambda);
                COMPLEX e_part=defect.epsilon(lambda);
                dielectric_constant epsilon=substrate.epsilon(lambda);

                COMPLEX x = kr*n;
                double k=2.*pi/lambda;

                // The incident and scattering khat vectors...
                Vector inK(sin(thetai),0.,cos(thetai));
                Vector outK(sin(thetas)*cos(phis),sin(thetas)*sin(phis),-cos(thetas));
                Vector inKr(inK.x,inK.y,-inK.z);
                Vector outKr(outK.x,outK.y,-outK.z);

                // The surface normal...
                Vector normal(0.,0.,-1.);

                // The following definitions make {s,p,k} right handed...
                Vector inS =  perpto(inK,normal);
                Vector inP =  perpto(inK,inS);
                Vector inSr =  inS;
                Vector inPr =  perpto(inKr,inSr);
                Vector outS =  perpto(outK,normal);
                Vector outP =  perpto(outK,outS);
                Vector outSr = outS;
                Vector outPr =  perpto(outKr,outSr);

                // Reflection coefficients...

                COMPLEX rsi = stack->rs21i(fabs(thetai),lambda,substrate,vacuum);
                COMPLEX rpi = stack->rp21i(fabs(thetai),lambda,substrate,vacuum);
                COMPLEX rss = stack->rs21i(fabs(thetas),lambda,substrate,vacuum);
                COMPLEX rps = stack->rp21i(fabs(thetas),lambda,substrate,vacuum);

                COMPLEX rayleigh = (e_part-e)/(e_part+2.*e)*cube(x);

                // Loss and phase factors...
                COMPLEX I(0,1);
                COMPLEX alpha = exp(2.*I*n*kd*cos(thetai));
                COMPLEX beta  = exp(2.*I*n*kd*cos(thetas));

                COMPLEX pp = outP*inP + outP*inPr*(alpha*rpi) + outPr*inP*(beta*rps) + outPr*inPr*(alpha*beta*rps*rpi);
                COMPLEX sp = outP*inS + outP*inSr*(alpha*rsi) + outPr*inS*(beta*rps) + outPr*inSr*(alpha*beta*rps*rsi);
                COMPLEX ps = outS*inP + outS*inPr*(alpha*rpi) + outSr*inP*(beta*rss) + outSr*inPr*(alpha*beta*rss*rpi);
                COMPLEX ss = outS*inS + outS*inSr*(alpha*rsi) + outSr*inS*(beta*rss) + outSr*inSr*(alpha*beta*rss*rsi);

                COMPLEX common=rayleigh/(k*n);
                // Return the whole value with factors common to all elements...
                return JonesMatrix(pp,ss,ps,sp)*common;

            } else { // is_transmission()

                COMPLEX n = substrate.index(lambda);
                COMPLEX e = substrate.epsilon(lambda);
                COMPLEX e_part=defect.epsilon(lambda);
                dielectric_constant epsilon=substrate.epsilon(lambda);

                COMPLEX x = kr*n;
                double k=2.*pi/lambda;

                COMPLEX sin_thetas_internal = sin(thetas)/n;
                COMPLEX cos_thetas_internal = sqrt(1.-sqr(sin_thetas_internal));

                // The incident and scattering khat vectors at the defect...
                CVector inK(sin(thetai),0.,cos(thetai));
                CVector outK(sin_thetas_internal*cos(phis),sin_thetas_internal*sin(phis),cos_thetas_internal);
                CVector inKr(inK.x,inK.y,-inK.z);

                // The surface normal...
                CVector normal(0.,0.,-1.);

                // The following definitions make {s,p,k} right handed...
                CVector inS =  perpto(inK,normal);
                CVector inP =  perpto(inK,inS);
                CVector inSr =  inS;
                CVector inPr =  perpto(inKr,inSr);
                CVector outS =  perpto(outK,normal);
                CVector outP =  perpto(outK,outS);

                // Reflection and transmission coefficients...

                COMPLEX rsi = stack->rs21i(fabs(thetai),lambda,substrate,vacuum);
                COMPLEX rpi = stack->rp21i(fabs(thetai),lambda,substrate,vacuum);
                COMPLEX tss = stack->ts21(thetas,lambda,substrate,vacuum)*sqrt(cos(thetas)/cos_thetas_internal/n);
                COMPLEX tps = stack->tp21(thetas,lambda,substrate,vacuum)*sqrt(cos(thetas)/cos_thetas_internal/n);

                COMPLEX rayleigh = (e_part-e)/(e_part+2.*e)*cube(x);

                // Loss and phase factors...
                COMPLEX I(0,1);
                COMPLEX alpha = exp(2.*I*n*kd*cos(thetai));

                COMPLEX pp = tps*(outP*inP + outP*inPr*(alpha*rpi));
                COMPLEX sp = tps*(outP*inS + outP*inSr*(alpha*rsi));
                COMPLEX ps = tss*(outS*inP + outS*inPr*(alpha*rpi));
                COMPLEX ss = tss*(outS*inS + outS*inSr*(alpha*rsi));

                COMPLEX jacobian = cos(thetas)/sqr(n)/cos_thetas_internal;
                COMPLEX common=rayleigh/(k*n)*sqrt(jacobian);
                // Return the whole value with factors common to all elements...
                return JonesMatrix(pp,ss,ps,sp)*common;
            }
        }
    }

    //
    // Routines to carry out one-time calculations...
    //
    void
    Rayleigh_Defect_BRDF_Model::
    setup()
    {
        // Call parent's setup()...
        Local_BRDF_Model::setup();

        double k = 2.*pi/lambda;
        kd = k*distance;
        kr = k*radius;
    }

    DEFINE_MODEL(Rayleigh_Defect_BRDF_Model,Local_BRDF_Model,
                 "Scattering by a subsurface defect in the Rayleigh limit.");

    DEFINE_PTRPARAMETER(Rayleigh_Defect_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_PARAMETER(Rayleigh_Defect_BRDF_Model,double,radius,"Defect radius [um]","0.001",0xFF);

    DEFINE_PARAMETER(Rayleigh_Defect_BRDF_Model,double,distance,"Distance from center to surface [um]","0",0xFF);

    DEFINE_PARAMETER(Rayleigh_Defect_BRDF_Model,dielectric_function,defect,"Defect","(1,0)",0xFF);



} // namespace SCATMECH


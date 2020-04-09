//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: raystack.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "raystack.h"

using namespace std;
using namespace SCATMECH;

static
complex<double>
arcsin(const complex<double>& x)
{
    complex<double> i(0.,1.);
    return -i*log(i*x+sqrt(1.-x*x));
}

//
// Jones matrix for scattering...
//
JonesMatrix Rayleigh_Stack_BRDF_Model::jonesDSC()
{
    SETUP();

    if (is_reflection()) {

        // This routine uses the results of Thomas A. Germer, "Application of
        // bidirectional ellipsometry to the characterization of roughness
        // and defects in dielectric layers," in Flatness, Roughness, and
        // Discrete Defect Characterization for Computer Disks, Wafers, and
        // Flat Panel Displays II," J.C. Stover, Editor,
        // Proc. SPIE 3275,121--131 (1998).
        //
        // In particular, it uses the result in Section 2.B. Rayleigh defect
        // in a multilayer system.
        //

        // Some abbreviations...
        complex<double> I(0,1);
        double k = 2*pi/lambda_eff;

        COMPLEX n3 = material3.index(lambda_eff);
        COMPLEX n2 = material2.index(lambda_eff);
        COMPLEX n1 = material1.index(lambda_eff);
        COMPLEX e3 = material3.epsilon(lambda_eff);
        COMPLEX e2 = material2.epsilon(lambda_eff);
        COMPLEX e1 = material1.epsilon(lambda_eff);

        COMPLEX q3i = k*sqrt(e3-sqr(sin(thetai)));
        COMPLEX q2i = k*sqrt(e2-sqr(sin(thetai)));
        COMPLEX q1i = k*sqrt(e1-sqr(sin(thetai)));
        COMPLEX q3s = k*sqrt(e3-sqr(sin(thetas)));
        COMPLEX q2s = k*sqrt(e2-sqr(sin(thetas)));
        COMPLEX q1s = k*sqrt(e1-sqr(sin(thetas)));

        // The total phase is not mentioned in the text, but must be kept
        // in mind if one is interested in coherent scattering between two
        // defects...
        complex<double> total_phase(0,0);
        for (int i=0; i<stack->get_n(); ++i) {
            complex<double> ee=stack->get_e()[i].epsilon(lambda_eff);
            double tt = k*stack->get_t()[i];
            ee = 1.;
            total_phase += sqrt(ee-sqr(sin(thetai)))*tt;
            total_phase += sqrt(ee-sqr(sin(thetas)))*tt;
        }

        // The internal angle of propagation...
        complex<double> sinthetaii = sin(thetai)/(complex<double>)(n2);
        complex<double> sinthetasi = sin(thetas)/(complex<double>)(n2);

        complex<double> costhetaii = sqrt(1.-sqr(sinthetaii));
        complex<double> costhetasi = sqrt(1.-sqr(sinthetasi));

        //complex<double> thetaii = arcsin(sinthetaii);
        //complex<double> thetasi = arcsin(sinthetasi);

        // Calculate the different transmission coefficients...
        complex<double> ts32i = above.ts12(thetai,lambda_eff,material3,material2);
        complex<double> tp32i = above.tp12(thetai,lambda_eff,material3,material2);
        complex<double> ts23s = above.ts21(thetas,lambda_eff,material2,material3);
        complex<double> tp23s = above.tp21(thetas,lambda_eff,material2,material3);

        // Calculate the different reflection coeff. for interface above
        // the defect...
        complex<double> rs23i = above.rs21(thetai,lambda_eff,material2,material3);
        complex<double> rp23i = above.rp21(thetai,lambda_eff,material2,material3);
        complex<double> rs23s = above.rs21(thetas,lambda_eff,material2,material3);
        complex<double> rp23s = above.rp21(thetas,lambda_eff,material2,material3);

        // Calculate the different reflection coeff. for interface below
        // the defect...
        complex<double> rs21i = below.rs12(thetai,lambda_eff,material2,material1);
        complex<double> rp21i = below.rp12(thetai,lambda_eff,material2,material1);
        complex<double> rs21s = below.rs12(thetas,lambda_eff,material2,material1);
        complex<double> rp21s = below.rp12(thetas,lambda_eff,material2,material1);


        // The phases between the defect and the layers...
        complex<double> p1 = q2i*(tau-d);
        complex<double> p2 = q2i*d;
        complex<double> p12 = q2i*tau;
        complex<double> p3 = q2s*(tau-d);
        complex<double> p4 = q2s*d;
        complex<double> p34 = q2s*tau;

        // The phase factors...
        complex<double> exp2ip12 = exp(2.*I*p12);
        complex<double> exp2ip34 = exp(2.*I*p34);
        complex<double> exp2ip1 = exp(2.*I*p1);
        complex<double> exp2ip2 = exp(2.*I*p2);
        complex<double> exp2ip3 = exp(2.*I*p3);
        complex<double> exp2ip4 = exp(2.*I*p4);

        // The various C factors...
        complex<double> Csi = 1. - rs21i * rs23i * exp2ip12;
        complex<double> Cpi = 1. - rp21i * rp23i * exp2ip12;
        complex<double> Css = 1. - rs21s * rs23s * exp2ip34;
        complex<double> Cps = 1. - rp21s * rp23s * exp2ip34;

        // The various B parameters...
        complex<double> Bs_plus_i = 1. + rs21i*exp2ip2;
        complex<double> Bp_plus_i = 1. + rp21i*exp2ip2;
        complex<double> Bp_minus_i = 1. - rp21i*exp2ip2;
        complex<double> Bs_plus_s = 1. + rs21s*exp2ip4;
        complex<double> Bp_plus_s = 1. + rp21s*exp2ip4;
        complex<double> Bp_minus_s = 1. - rp21s*exp2ip4;

        JonesMatrix j;
        j[1] = ts32i*ts23s*Bs_plus_s*Bs_plus_i*cos(phis)/Css/Csi;
        j[2] = -tp32i*ts23s*Bs_plus_s*Bp_minus_i*costhetaii*sin(phis)/Css/Cpi;
        j[3] = -ts32i*tp23s*Bp_minus_s*Bs_plus_i*costhetasi*sin(phis)/Cps/Csi;
        j[0] = tp32i*tp23s*(
                   Bp_plus_s*Bp_plus_i*sinthetaii*sinthetasi
                   -Bp_minus_s*Bp_minus_i*costhetaii*costhetasi*cos(phis))/Cps/Cpi;

        // Return the matrix with the overall phase factor...
        //j = j*exp(I*(p1+p3))*exp(-I*total_phase);
        j = j*exp(I*(p1+p3));

        COMPLEX ns = materialS.index(lambda_eff);

        COMPLEX polarizability = (sqr(ns)-sqr(n2))/(sqr(ns)+2.*sqr(n2))*cube(radius);

        COMPLEX csect2=polarizability*sqr(k*n2);

        COMPLEX jacobian = 1./sqr((COMPLEX)n2);

        j = j*cos(thetas)/costhetasi*csect2*sqrt(jacobian);

        if (is_backward()) {
            return JonesMatrix(j[0],j[1],-j[2],-j[3]);
        } else {
            return j;
        }
    } else {

        complex<double> I(0,1);
        double k = 2*pi/lambda_eff;
        COMPLEX n3 = material3.index(lambda_eff);
        COMPLEX n2 = material2.index(lambda_eff);
        COMPLEX n1 = material1.index(lambda_eff);
        COMPLEX e3 = material3.epsilon(lambda_eff);
        COMPLEX e2 = material2.epsilon(lambda_eff);
        COMPLEX e1 = material1.epsilon(lambda_eff);

        COMPLEX q3i = k*sqrt(e3-sqr(sin(thetai)));
        COMPLEX q2i = k*sqrt(e2-sqr(sin(thetai)));
        COMPLEX q1i = k*sqrt(e1-sqr(sin(thetai)));
        COMPLEX q3s = k*sqrt(e3-sqr(sin(thetas)*n1));
        COMPLEX q2s = k*sqrt(e2-sqr(sin(thetas)*n1));
        COMPLEX q1s = k*sqrt(e1-sqr(sin(thetas)*n1));

        // External angle on scattering...
        COMPLEX sinthetasext = sin(thetas)*n1;
        COMPLEX thetasext = arcsin(sinthetasext);

        // The total phase is not mentioned in the text, but must be kept in mind if one
        // is interested in coherent scattering between two defects...
        complex<double> total_phase(0,0);
        for (int i=0; i<stack->get_n(); ++i) {
            complex<double> ee=stack->get_e()[i].epsilon(lambda_eff);
            double tt = k*stack->get_t()[i];
            total_phase += sqrt(ee-sqr(sin(thetai)))*tt;
            total_phase -= sqrt(ee-sqr(sinthetasext))*tt;
        }

        // Calculate the different transmission coefficients...
        complex<double> ts32i = above.ts12(thetai,lambda_eff,material3,material2);
        complex<double> tp32i = above.tp12(thetai,lambda_eff,material3,material2);
        complex<double> ts21s = below.ts12(thetasext,lambda_eff,material2,material1);
        complex<double> tp21s = below.tp12(thetasext,lambda_eff,material2,material1);

        // Calculate the different reflection coeff. for interface above
        // the defect...
        complex<double> rs23i = above.rs21(thetai,lambda_eff,material2,material3);
        complex<double> rp23i = above.rp21(thetai,lambda_eff,material2,material3);
        complex<double> rs23s = above.rs21(thetasext,lambda_eff,material2,material3);
        complex<double> rp23s = above.rp21(thetasext,lambda_eff,material2,material3);

        // Calculate the different reflection coeff. for interface below
        // the defect...
        complex<double> rs21i = below.rs12(thetai,lambda_eff,material2,material1);
        complex<double> rp21i = below.rp12(thetai,lambda_eff,material2,material1);
        complex<double> rs21s = below.rs12(thetasext,lambda_eff,material2,material1);
        complex<double> rp21s = below.rp12(thetasext,lambda_eff,material2,material1);

        // The internal angle of propagation...
        complex<double> thetaii = arcsin(sin(thetai)/n2);
        complex<double> thetasi = arcsin(sinthetasext/n2);

        // The phases between the defect and the layers...
        complex<double> p1 = d*q2i;
        complex<double> p2 = (tau-d)*q2i;
        complex<double> p3 = d*q2s;
        complex<double> p4 = (tau-d)*q2s;

        // The phase factors...
        complex<double> exp2ip12 = exp(2.*I*(p1+p2));
        complex<double> exp2ip34 = exp(2.*I*(p3+p4));
        complex<double> exp2ip1 = exp(2.*I*p1);
        complex<double> exp2ip2 = exp(2.*I*p2);
        complex<double> exp2ip3 = exp(2.*I*p3);
        complex<double> exp2ip4 = exp(2.*I*p4);

        // The various C factors...
        complex<double> Csi = 1. - rs21i * rs23i * exp2ip12;
        complex<double> Cpi = 1. - rp21i * rp23i * exp2ip12;
        complex<double> Css = 1. - rs21s * rs23s * exp2ip34;
        complex<double> Cps = 1. - rp21s * rp23s * exp2ip34;

        // The various D factors...
        complex<double> Ds_plus  = 1. + rs23s * exp2ip4;
        complex<double> Dp_plus  = 1. + rp23s * exp2ip4;
        complex<double> Dp_minus = 1. - rp23s * exp2ip4;

        // The various B factors...
        complex<double> Bs_plus  = 1. + rs21i * exp2ip1;
        complex<double> Bp_plus  = 1. + rp21i * exp2ip1;
        complex<double> Bp_minus = 1. - rp21i * exp2ip1;

        // The matrix elements...
        JonesMatrix j;
        j.SS() = ts32i*ts21s*Ds_plus*Bs_plus*cos(phis)/Css/Csi;
        j.SP() = ts32i*tp21s*Dp_minus*Bs_plus*cos(thetasi)*sin(phis)/Cps/Csi;
        j.PS() = -tp32i*ts21s*Ds_plus*Bp_minus*cos(thetaii)*sin(phis)/Css/Cpi;
        j.PP() = tp32i*tp21s*(Dp_plus*Bp_plus*sin(thetaii)*sin(thetasi)+Dp_minus*Bp_minus*cos(thetaii)*cos(thetasi)*cos(phis))/Cps/Cpi;
        // Return the matrix with the overall phase factor...

        COMPLEX ns = materialS.index(lambda_eff);
        COMPLEX polarizability = (sqr(ns)-sqr(n2))/(sqr(ns)+2.*sqr(n2))*cube(radius);
        //COMPLEX polarizability = (sqr(ns)-sqrt(n2))/(4.*pi);
        COMPLEX csect2=polarizability*sqr(k*n2);
        COMPLEX jacobian = sqr(n1/n2);

        j = j*(cos(thetas)/cos(thetasi)*csect2*sqrt(jacobian*n1));
        j = j*exp(I*(p2+p3));

        if (is_backward()) {
            return JonesMatrix(j[0],j[1],-j[2],-j[3]);
        } else {
            return j;
        }

    }
}

//
// Routine to perform housekeeping which only is neccesary when
// the defect location changes...
//
void Rayleigh_Stack_BRDF_Model::setup()

{
    // Call parent's setup()...
    Local_BRDF_Model::setup();

    double bscale = substrate.n(lambda);

    // Remove all films...
    above.wash();
    below.wash();

    int layer_number;

    double stackthickness = stack->get_total_thickness();

    if (depth>=0 && depth<stackthickness) {
        layer_number=0;

        double _depth = depth;
        for (int ii=stack->get_n()-1; ii>=0; --ii) {
            double t = stack->get_t()[ii];
            if (_depth<=t&&_depth>=0) {
                layer_number = ii;
                if (is_forward()) {
                    d = t-_depth;
                } else {
                    d = _depth;
                }
            }
            _depth -= t;
        }

        // Grow films below the defect...
        if (is_forward()) {
            for (int ii=0; ii<layer_number; ++ii) {
                below.grow(stack->get_e()[ii].index(lambda),stack->get_t()[ii]);
            }
        } else {
            for (int ii=layer_number-1; ii>=0; --ii) {
                above.grow(optical_constant(COMPLEX(stack->get_e()[ii].index(lambda))/bscale),stack->get_t()[ii]);
            }
        }

        // Get defect's layer information...
        int i=layer_number;
        if (is_forward()) {
            material2 = stack->get_e()[i].index(lambda);
        } else {
            material2 = optical_constant(COMPLEX(stack->get_e()[i].index(lambda))/bscale);
        }
        tau = stack->get_t()[i];

        // Grow films above the defect...
        if (is_forward()) {
            for (int i=layer_number+1; i<stack->get_n(); ++i) {
                above.grow(stack->get_e()[i].index(lambda),stack->get_t()[i]);
            }
        } else {
            for (int i=stack->get_n()-1; i>=layer_number+1; --i) {
                below.grow(optical_constant(COMPLEX(stack->get_e()[i].index(lambda))/bscale),stack->get_t()[i]);
            }
        }
    } else if (depth<0 && is_forward()) {
        // Grow films below the defect...
        for (int ii=0; ii<stack->get_n(); ++ii) {
            below.grow(stack->get_e()[ii].index(lambda),stack->get_t()[ii]);
        }

        // Get defect's layer information...
        material2 = optical_constant(1.0);
        tau = -depth;
        d = tau;
    } else if (depth<0 && is_backward()) {

        // Get defect's layer information...
        material2 = optical_constant(1./bscale);
        tau = -depth;
        d = 0;
        for (int ii=stack->get_n()-1; ii>=0; --ii) {
            above.grow(optical_constant(COMPLEX(stack->get_e()[ii].index(lambda))/bscale),stack->get_t()[ii]);
        }

    } else if (depth>stackthickness && is_forward()) {
        // Get defect's layer information...
        material2 = substrate.index(lambda);
        tau = depth-stackthickness;
        d = 0;

        // Grow films above the defect...
        for (int i=0; i<stack->get_n(); ++i) {
            above.grow(stack->get_e()[i].index(lambda),stack->get_t()[i]);
        }

    } else {
        // Get defect's layer information...
        material2 = optical_constant(1.);
        tau = depth-stackthickness;
        d = tau;

        // Grow films above the defect...
        for (int i=stack->get_n()-1; i>=0; --i) {
            below.grow(optical_constant(COMPLEX(stack->get_e()[i].index(lambda))/bscale),stack->get_t()[i]);
        }
    }

    if (is_forward()) {
        material3 = optical_constant(1.);
        material1 = substrate.index(lambda);
        materialS = sphere.index(lambda);
        lambda_eff = lambda;
    } else {
        material3 = optical_constant(1.);
        material1 = optical_constant(1./bscale);
        materialS = optical_constant(COMPLEX(sphere.index(lambda))/bscale);
        lambda_eff = lambda/bscale;
    }
}

DEFINE_MODEL(Rayleigh_Stack_BRDF_Model,Local_BRDF_Model,"A Rayleigh defect in a stack of layers");
DEFINE_PTRPARAMETER(Rayleigh_Stack_BRDF_Model,StackModel_Ptr,stack,"Dielectric stack","No_StackModel",0xFF);
DEFINE_PARAMETER(Rayleigh_Stack_BRDF_Model,double,depth,"Depth into stack [um]","0",0xFF);
DEFINE_PARAMETER(Rayleigh_Stack_BRDF_Model,double,radius,"Sphere radius [um]","0.01",0xFF);
DEFINE_PARAMETER(Rayleigh_Stack_BRDF_Model,dielectric_function,sphere,"Sphere optical properties","(1,0)",0xFF);

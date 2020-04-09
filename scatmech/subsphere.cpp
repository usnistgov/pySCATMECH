//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: subsphere.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "subsphere.h"
#include "askuser.h"
#include "rayscat.h"

using namespace std;


namespace SCATMECH {


    static
    COMPLEX
    ArcSin(const COMPLEX& a)
    {
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }

    void
    Subsurface_Particle_BRDF_Model::
    setup()
    {
        Local_BRDF_Model::setup();

        if (substrate.k(lambda)!=0)
            error("Substrate cannot be absorbing");
        if (scatterer->get_lambda()!=lambda) error("scatterer.lambda!=lambda");
        if (scatterer->get_medium().n(lambda)!=substrate.n(lambda)) error("scatterer.medium!=substrate");

        Euler = Matrix(cos(alpha)*cos(beta),  sin(alpha), -cos(alpha)*sin(beta),
                       -cos(beta)*sin(alpha), cos(alpha), sin(alpha)*sin(beta),
                       sin(beta),                 0,            cos(beta));
    }

    JonesMatrix
    Subsurface_Particle_BRDF_Model::
    jonesDSC()
    {
        SETUP();

        double n = substrate.n(lambda);
        double k=2.*pi/lambda;

        Matrix _Euler = Euler*Matrix(cos(rotation),sin(rotation),0,
                                     -sin(rotation),cos(rotation),0,
                                     0,0,1);

        if (fabs(sin(thetai)-sin(thetas)*cos(phis))<1E-5) thetai+=2.5E-5;
        if (fabs(sin(thetai)+sin(thetas)*cos(phis))<1E-5) thetai+=2.5E-5;
        if (fabs(thetai)<1E-5) thetai=1E-5;
        if (fabs(thetas)<1E-5) thetas=1E-5;

        if (is_forward()) {

            // Angle of incidence in the vicinity of the particle...
            double sin_thetai_internal = sin(thetai)/n;
            double cos_thetai_internal = sqrt(1.-sqr(sin_thetai_internal));

            // Scatterer does not handle evanescent waves...
            if (fabs(sin_thetai_internal)>1) return JonesZero();

            // The incident direction local to the particle...
            Vector inKi(sin_thetai_internal,0.,-cos_thetai_internal);

            // The local polarization basis... {s,p,k} must be right handed...
            Vector normal(0.,0.,1.);
            Vector inSi =  perpto(inKi,normal);
            Vector inPi =  cross(inKi,inSi);

            if (is_reflection()) {

                // Scattering direction local to the particle...
                double sin_thetas_internal = sin(thetas)/n;

                // Scatterer doesn't handle evanescent waves...
                if (fabs(sin_thetas_internal)>1) return JonesZero();

                double cos_thetas_internal = sqrt(1.-sqr(sin_thetas_internal));

                // Scatterering direction, local to particle...
                Vector outKi(sin_thetas_internal*cos(phis),sin_thetas_internal*sin(phis),cos_thetas_internal);

                // Polarization basis set for scattering direction make {s,p,k} right handed...
                Vector outSi = perpto(outKi,normal);
                Vector outPi = cross(outKi,outSi);

                // Transmission coefficients into and out of the material...
                JonesMatrix ti = stack->t12(thetai,lambda,vacuum,substrate)*(COMPLEX)sqrt(n);
                JonesMatrix ts = stack->t12(thetas,lambda,vacuum,substrate)*(COMPLEX)sqrt(n*cos_thetas_internal/cos(thetas));

                // Polarization basis local to the scattering plane... {par*,perp,k} are right handed...
                Vector perp = perpto(inKi,outKi);
                Vector pari = cross(perp,inKi);
                Vector pars = cross(perp,outKi);

                // Conversion between {par*,perp,k} and {s,p,k} ...
                JonesMatrix matrixin = GetJonesRotator(pari,perp,inSi,inPi);
                JonesMatrix matrixout = GetJonesRotator(outSi,outPi,pars,perp);

                // Scattering matrix in {par*,perp,k} basis...
                JonesMatrix scatter = scatterer->jones(_Euler*inKi,_Euler*outKi);

                // Convert to {s,p,k} basis...
                scatter=matrixout*scatter*matrixin;

                // The Jacobian handles the change in the solid angle on refraction...
                double jacobian = cos(thetas)/cos_thetas_internal/sqr(n);

                // The scatter matrix globally includes transmission to/from defect...
                scatter = ts*scatter*ti;

                COMPLEX common=jacobian/sqr(k*n);

                // Return the whole value with factors common to all elements...
                return scatter*sqrt(common);

            } else { // is_transmission()

                if (abs(thetas)<1E-8) thetas=1E-8;

                // Sines and cosines of incident angles at the particle...
                double sin_thetai = sin(thetai)/n;
                double cos_thetai = sqrt(1.-sqr(sin_thetai));

                // Scatterer does not handle evanescent fields...
                if (fabs(sin_thetai)>1) return JonesZero();

                // Incident direction local to particle...
                Vector inKi(sin_thetai,0,-cos_thetai);

                // Direction of scattering and polarization basis vectors...
                Vector outKi(sin(thetas)*cos(phis),sin(thetas)*sin(phis),-cos(thetas));
                Vector outSi = perpto(outKi,normal);
                Vector outPi = cross(outKi,outSi);

                // Direction of scattering which is reflected into scattering direction...
                Vector outKr(sin(thetas)*cos(phis),sin(thetas)*sin(phis),cos(thetas));
                Vector outSr = outSi;
                Vector outPr = cross(outKr,outSr);

                // Calculate the direct scattered waves...
                JonesMatrix scatter_direct;
                {
                    // Local polarization basis vectors....
                    Vector perp = perpto(inKi,outKi);
                    Vector pari = cross(perp,inKi);
                    Vector pars = cross(perp,outKi);

                    // Rotation matrices...
                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inSi,inPi);
                    JonesMatrix matrixout = GetJonesRotator(outSi,outPi,pars,perp);

                    // Scattering matrix...
                    JonesMatrix scatter = scatterer->jones(_Euler*inKi,_Euler*outKi);

                    // Scattering matrix in global coordinates...
                    scatter_direct=matrixout*scatter*matrixin;
                }

                // Calculate the scattered wave reflected in interface...
                JonesMatrix scatter_indirect;
                {
                    // Local polarization basis vectors...
                    Vector perp = perpto(inKi,outKr);
                    Vector pari = cross(perp,inKi);
                    Vector pars = cross(perp,outKr);

                    // The reflection coefficients...
                    JonesMatrix r = stack->r21i(fabs(thetas),lambda,substrate,vacuum);

                    // Rotation matrices...
                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inSi,inPi);
                    JonesMatrix matrixout = GetJonesRotator(outSr,outPr,pars,perp);;

                    // Scattering matrix...
                    JonesMatrix scatter = scatterer->jones(_Euler*inKi,_Euler*outKr);

                    // Phase difference due to extra path length...
                    COMPLEX phase = exp(COMPLEX(0,2)*cos(thetas)*k*depth*n);

                    // Scattering matrix in global coordinates, accounting for phase and reflection
                    scatter_indirect=r*matrixout*scatter*matrixin*phase;
                }

                // Transmission into substrate...
                JonesMatrix ti = stack->t12(thetai,lambda,vacuum,substrate)*(COMPLEX)sqrt(n);

                JonesMatrix scatter = (scatter_direct+scatter_indirect)*ti/k/n;

                return scatter;
            }
        } else { // is_backward()
            if (is_reflection()) {

                double n = substrate.n(lambda);

                double k = 2.*pi/lambda;
                double kd = k*depth*n;

                // The incident and scattering khat vectors...
                Vector inK(sin(thetai),0.,cos(thetai));
                Vector outK(sin(thetas)*cos(phis),sin(thetas)*sin(phis),-cos(thetas));

                // The surface normal ...
                Vector normal(0.,0.,-1.);

                // The reflected khats...
                Vector inKr(inK.x,inK.y,-inK.z);
                Vector outKr(outK.x,outK.y,-outK.z);

                // The polarization vectors...
                // The following definitions make {s,p,k} right handed...
                Vector inS =   perpto(inK,normal);
                Vector inP =   perpto(inK,inS);
                Vector outS =  perpto(outK,normal);
                Vector outP =  perpto(outK,outS);

                // The reflected light polarization vectors...
                Vector inSr =   inS;
                Vector inPr =   perpto(inKr,inSr);
                Vector outSr =  outS;
                Vector outPr =  perpto(outKr,outSr);

                // Basis vectors appropriate for the scattering plane frame of reference
                // "perp" vectors...
                Vector b1in = perpto(inK,outK),   b1out=b1in;
                Vector b2in = perpto(inKr,outK),  b2out=b2in;
                Vector b3in = perpto(inK,outKr),  b3out=b3in;
                Vector b4in = perpto(inKr,outKr), b4out=b4in;
                // "par" vectors...
                Vector a1in = cross(b1in,inK),   a1out = cross(b1out,outK);
                Vector a2in = cross(b2in,inKr),  a2out = cross(b2out,outK);
                Vector a3in = cross(b3in,inK),   a3out = cross(b3out,outKr);
                Vector a4in = cross(b4in,inKr),  a4out = cross(b4out,outKr);

                // The reflection coefficients...
                JonesMatrix ri = stack->r21i(fabs(thetai),lambda,substrate,vacuum);
                JonesMatrix rs = stack->r21i(fabs(thetas),lambda,substrate,vacuum);

                // The scattering matrices...
                JonesMatrix S1 = scatterer->jones(_Euler*inK,_Euler*outK);
                JonesMatrix S2 = scatterer->jones(_Euler*inKr,_Euler*outK);
                JonesMatrix S3 = scatterer->jones(_Euler*inK,_Euler*outKr);
                JonesMatrix S4 = scatterer->jones(_Euler*inKr,_Euler*outKr);

                // Phase factors between the different waves...
                COMPLEX I(0,1);
                COMPLEX alpha = exp(2.*I*kd*cos(thetai));
                COMPLEX beta  = exp(2.*I*kd*cos(thetas));

                // Basis rotation matrices...
                JonesMatrix RotIn1  = GetJonesRotator(a1in,b1in,inS,inP);
                JonesMatrix RotIn2  = GetJonesRotator(a2in,b2in,inSr,inPr);
                JonesMatrix RotIn3  = GetJonesRotator(a3in,b3in,inS,inP);
                JonesMatrix RotIn4  = GetJonesRotator(a4in,b4in,inSr,inPr);

                JonesMatrix RotOut1 = GetJonesRotator(outS,outP,a1out,b1out);
                JonesMatrix RotOut2 = GetJonesRotator(outS,outP,a2out,b2out);
                JonesMatrix RotOut3 = GetJonesRotator(outSr,outPr,a3out,b3out);
                JonesMatrix RotOut4 = GetJonesRotator(outSr,outPr,a4out,b4out);

                S1 =    RotOut1*S1*RotIn1;
                S2 =    RotOut2*S2*RotIn2*alpha*ri;
                S3 = rs*RotOut3*S3*RotIn3*beta;
                S4 = rs*RotOut4*S4*RotIn4*alpha*beta*ri;

                COMPLEX common=1./k/n;

                return (S1+S2+S3+S4)*common;

            } else { // is_transmission()

                double n = substrate.n(lambda);
                double k=2.*pi/lambda;

                Vector normal(0,0,1);

                // Incident directions on particle
                Vector inK(sin(thetai),0.,cos(thetai));
                Vector inKr(inK.x,inK.y,-inK.z);

                // Polarization bases for incident waves...
                Vector inS = perpto(inK,normal);
                Vector inP = cross(inK,inS);
                Vector inSr = inS;
                Vector inPr = cross(inKr,inSr);

                // Sine of scattering angle local to particle
                double sin_thetas_inside = sin(thetas)/n;

                // This model does not handle evanescent waves at the particle...
                if (fabs(sin_thetas_inside)>1.) return JonesZero();

                double thetas_inside = asin(sin_thetas_inside);
                double cos_thetas_inside = sqrt(1.-sqr(sin_thetas_inside));

                // Direction of waves scattered from particle and polarization basis...
                Vector outK(sin_thetas_inside*cos(phis),sin_thetas_inside*sin(phis),cos_thetas_inside);
                Vector outS = perpto(outK,normal);
                Vector outP = cross(outK,outS);

                // Illumination directly from source
                JonesMatrix scatter_direct;
                {
                    // Local polarization basis...
                    Vector perp = perpto(inK,outK);
                    Vector pari = cross(perp,inK);
                    Vector pars = cross(perp,outK);

                    // Basis rotation matrices...
                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inS,inP);
                    JonesMatrix matrixout = GetJonesRotator(outS,outP,pars,perp);

                    // Scattering matrix in local basis...
                    JonesMatrix scatter = scatterer->jones(_Euler*inK,_Euler*outK);

                    // Scattering matrix in global basis...
                    scatter_direct = matrixout*scatter*matrixin;
                }

                // Illumination through reflection...
                JonesMatrix scatter_indirect;
                {
                    // Local polarization basis...
                    Vector perp = perpto(inKr,outK);
                    Vector pari = cross(perp,inKr);
                    Vector pars = cross(perp,outK);

                    // Basis rotation matrices...
                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inSr,inPr);
                    JonesMatrix matrixout = GetJonesRotator(outS,outP,pars,perp);

                    // Scattering matrix in local basis...
                    JonesMatrix scatter = scatterer->jones(_Euler*inKr,_Euler*outK);

                    // Phase and reflection coefficients...
                    COMPLEX phase = exp(COMPLEX(0,2)*cos(thetai)*k*n*depth);
                    JonesMatrix r = stack->r21i(fabs(thetai),lambda,substrate,vacuum);

                    // Scattering matrix in global basis...
                    scatter_indirect = phase*matrixout*scatter*matrixin*r;
                }

                // Transmission out of material...
                JonesMatrix t = stack->t21(thetas,lambda,substrate,vacuum)*
                                sqrt(cos(thetas)/cos(thetas_inside)/n);

                JonesMatrix scatter = t*(scatter_direct+scatter_indirect)/k/n;

                // Jacobian accounts for refraction of solid angle through interface...
                double jacobian = cos(thetas)/cos(thetas_inside)/sqr(n);

                return scatter*sqrt(jacobian);

            }
        }
    }

    DEFINE_MODEL(Subsurface_Particle_BRDF_Model,Local_BRDF_Model,
                 "Scattering by a particle beneath an interface of a nonabsorbing media");

    DEFINE_PTRPARAMETER(Subsurface_Particle_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);
    DEFINE_PARAMETER(Subsurface_Particle_BRDF_Model,double,depth,"Depth of center of particle [um]","0",0xFF);
    DEFINE_PTRPARAMETER(Subsurface_Particle_BRDF_Model,Free_Space_Scatterer_Ptr,scatterer,"The scattering function","MieScatterer",0xFF);
    DEFINE_PARAMETER(Subsurface_Particle_BRDF_Model,double,alpha,"Rotation of particle about z-axis [rad]","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Particle_BRDF_Model,double,beta,"Rotation of particle about y-axis [rad]","0",0xFF);



} //namespace SCATMECH



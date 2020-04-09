//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: facet.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "facet.h"

using namespace std;


namespace SCATMECH {

    double
    Facet_BRDF_Model::
    jacobian(double ti,double tt,double pt)
    {
        if (is_transmission()) {
            double n = is_forward() ? substrate.n(lambda) : 1./substrate.n(lambda);
            double j = (sqr(n)*(-(n*cos(ti)*(3 + cos(2*tt))) +
                                cos(pt)*sin(2*ti)*sin(tt) +
                                cos(tt)*sqr(cos(ti)) +
                                cos(tt)*(1 - 2*n*cos(pt)*sin(ti)*sin(tt) +
                                         2*sqr(n) - sqr(sin(ti)))))/
                       (2.*(-cos(ti) + n*cos(tt))*
                        pow(1 - 2*n*cos(ti)*cos(tt) -
                            2*n*cos(pt)*sin(ti)*sin(tt) + sqr(n),1.5));
            return j;
        } else { // is_reflection()
            return 1./sqrt(8.*(1. + cos(ti)*cos(tt) - cos(pt)*sin(ti)*sin(tt)));
        }
    }

    //
    // The following returns the local incident angle for a specific in/out
    // direction... (From Eq. 2 of Barrick)
    //
    double
    Facet_BRDF_Model::
    local_angle(double thetai,double thetas,double phis)
    {
        Vector ki(sin(thetai),0.,-cos(thetai));

        double cos_angle = -nhat(thetai,thetas,phis)*ki;
        if (cos_angle > 1.) cos_angle = 1.;

        return acos(fabs(cos_angle));
    }

    //
    // The following return the surface slope for a specific in/out direction...
    // (from Eq. 9 of Barrick)
    //
    double
    Facet_BRDF_Model::
    local_slope(double thetai,double thetas,double phis)
    {
        double nz = nhat(thetai,thetas,phis).z;
        return sqrt(1.-sqr(nz))/nz; // = tan(acos(nz));
    }

    double
    Facet_BRDF_Model::
    evaluate_sdf(const Vector& _nhat)
    {
        double cosa = cos(rotation);
        double sina = sin(rotation);
        double sx = _nhat.x/_nhat.z;
        double sy = _nhat.y/_nhat.z;
        return sdf->f(cosa*sx+sina*sy,-sina*sx+cosa*sy);
    }

    Vector
    Facet_BRDF_Model::
    nhat(double ti,double tt,double pt)
    {
        if (is_transmission()) {
            double n = is_forward() ? substrate.n(lambda) : 1./substrate.n(lambda);

            double denom = sqrt(1 - 2*n*cos(ti)*cos(tt) -
                                2*n*cos(pt)*sin(ti)*sin(tt) + sqr(n));

            Vector result((sin(ti) - n*cos(pt)*sin(tt))/denom,
                          -(n*sin(pt)*sin(tt))/denom,
                          (-cos(ti) + n*cos(tt))/denom);

            if (result.z>=0.)   return result;
            else                return -result;
        } else { // is_reflection()
            Vector ki=polar(1.,pi-ti,0.);
            Vector kt=polar(1.,tt,pt);
            return unit(kt-ki);
        }

    }

    //
    // Calculate Jones matrix for scattering...
    //
    JonesMatrix
    Facet_BRDF_Model::
    jones()
    {
        SETUP();

        if (is_transmission()) {

            Vector _nhat = nhat(thetai,thetas,phis);
            Vector ki=polar(1.,pi-thetai,0.);
            Vector kt=polar(1.,pi-thetas,phis);

            double cos_iota_i = -ki*_nhat;
            if (cos_iota_i > 1.) cos_iota_i = 1.;
            double cos_iota_t = -kt*_nhat;
            if (cos_iota_i < 0.) return JonesZero();
            if (cos_iota_t < 0.) return JonesZero();

            double tn = acos(_nhat.z);
            double pn = atan2(_nhat.y,_nhat.x);

            double ti = thetai;
            double tt = thetas;
            double pt = phis;

            double iota = acos(cos_iota_i);
            COMPLEX ts,tp;
            if (is_forward()) {
                ts = stack->ts12(iota,lambda,vacuum,substrate);
                tp = stack->tp12(iota,lambda,vacuum,substrate);
            } else {
                ts = stack->ts21i(iota,lambda,substrate,vacuum);
                tp = stack->tp21i(iota,lambda,substrate,vacuum);
            }

            JonesMatrix j;
            if (fabs(iota) > 1E-7) {

                double a1 = cos(pn - pt)*cos(tt)*sin(tn)*(cos(tn)*sin(ti) + cos(pn)*cos(ti)*sin(tn))+
                            cos(tn)*(cos(tn)*sin(ti) + cos(pn)*cos(ti)*sin(tn))*sin(tt);
                double a2 = sin(pn)*sin(pn - pt)*sqr(sin(tn));
                double a3 = sin(tn)*(cos(tn)*sin(pn)*sin(tt)+cos(pn - pt)*cos(tt)*sin(pn)*sin(tn));
                double a4 = sin(tn)*(sin(pn - pt)*(cos(tn)*sin(ti) + cos(pn)*cos(ti)*sin(tn)));
                double factor1 = sqrt(2*cos(pn)*cos(ti)*cos(tn)*sin(ti)*sin(tn) +
                                      sqr(cos(tn)*sin(ti)) + (sqr(cos(pn)*cos(ti)) + sqr(sin(pn)))*
                                      sqr(sin(tn)))*sqrt(sqr(sin(pn - pt)*sin(tn)* sin(tt)) +
                                              sqr(cos(pn)*cos(tt)*sin(tn) + cos(pt)*cos(tn)*sin(tt)) +
                                              sqr(cos(tt)*sin(pn)*sin(tn) + cos(tn)*sin(pt)*sin(tt)));

                j.SS() = (tp*a2 + ts*a1)/factor1;
                j.PS() = (tp*a4 - ts*a3)/factor1;
                j.SP() = (tp*a3 - ts*a4)/factor1;
                j.PP() = (tp*a1 + ts*a2)/factor1;
            } else {
                j.SS()=ts;
                j.PP()=tp;
                j.SP()=0.;
                j.PS()=0.;
            }

            // Angle distribution function...
            // double s=sqr(1.+sqr(slope))* sdf->f(slope); // = sec(acos(nhat.z))^4 *sdf->(slope)
            double s = 1./pow(_nhat.z,4.) * evaluate_sdf(_nhat);

            // Transmission Coefficient to Transmittance ratio...
            if (is_forward()) {
                s *= substrate.n(lambda)*cos_iota_t/cos_iota_i;
            } else {
                s *= cos_iota_t/cos_iota_i/substrate.n(lambda);
                j.SP()=-j.SP();
                j.PS()=-j.PS();
            }

            // BTDF requires this factor...
            s *= 1./cos(thetas);

            // Probability of striking this facet is proportional to
            // the exposed area of the facet [cos(iota)] normalized by
            // the average of that exposed area [approx. cos(thetai)]...
            s *= cos_iota_i/cos(thetai);

            // Conversion from dOmega_n to dOmega_t...
            s *= fabs(jacobian(thetai,thetas,phis));

            return j*sqrt(s);

        } else { // REFLECTION...
            Vector _nhat = nhat(thetai,thetas,phis);
            Vector ki=polar(1.,pi-thetai,0.);

            double cos_iota = -ki*_nhat;
            if (cos_iota > 1.) cos_iota = 1.;
            double iota = acos(cos_iota);

            COMPLEX rsiota,rpiota;
            if (is_forward()) {
                rsiota = stack->rs12(iota,lambda,vacuum,substrate);
                rpiota = stack->rp12(iota,lambda,vacuum,substrate);
            } else {
                rsiota = stack->rs21i(iota,lambda,substrate,vacuum);
                rpiota = stack->rp21i(iota,lambda,substrate,vacuum);
            }

            double a1 = sqr(sin(2*iota));
            double a2 = cos(thetai)*sin(thetas)+sin(thetai)*cos(thetas)*cos(phis);
            double a3 = sin(thetai)*cos(thetas)+cos(thetai)*sin(thetas)*cos(phis);

            JonesMatrix j;
            if (fabs(iota) > 1E-10) {
                j[0] = (sin(thetai)*sin(thetas)*sqr(sin(phis))*rsiota+a2*a3*rpiota)/a1;
                j[1] = (a2*a3*rsiota+sin(thetai)*sin(thetas)*sqr(sin(phis))*rpiota)/a1;
                j[2] = -sin(phis)*(sin(thetas)*a2*rsiota-sin(thetai)*a3*rpiota)/a1;
                j[3] = -sin(phis)*(sin(thetai)*a3*rsiota-sin(thetas)*a2*rpiota)/a1;
            } else {
                j[0] = rpiota;
                j[1] = rsiota;
                j[2] = 0.;
                j[3] = 0.;
            }

            if (is_backward()) {
                j.SP()=-j.SP();
                j.PS()=-j.PS();
            }

            // Angle distribution function...
            // double s=sqr(1.+sqr(slope))* sdf->f(slope); // = sec(acos(nhat.z))^4 *sdf->(slope)
            double s = 1./pow(_nhat.z,4.) * evaluate_sdf(_nhat);

            // BTDF requires this factor...
            s *= 1./cos(thetas);

            // Probability of striking this facet is proportional to
            // the exposed area of the facet [cos(iota)] normalized by
            // the average of that exposed area [approx. cos(thetai)]...
            s *= cos_iota/cos(thetai);

            // Conversion from dOmega_n to dOmega_t...
            s *= fabs(jacobian(thetai,thetas,phis));

            return j*sqrt(s);

        }
    }

    double
    Gaussian_Slope_Distribution_Function::
    f(double slope)
    {
        SETUP();
        return  1./sqr(s)/pi*exp(-sqr(slope/s));
    }

    double
    Exponential_Slope_Distribution_Function::
    f(double slope)
    {
        SETUP();
        const double sqrt6 = 2.44948974278;
        const double threeoverpi = 0.954929658551;
        return threeoverpi/sqr(s)*exp(-sqrt6*slope/s);
    }

    double
    Table_Slope_Distribution_Function::
    f(double slope)
    {
        SETUP();
        return T.value(slope);
    }

    double
    Anisotropic_Exponential_Slope_Distribution_Function::
    f(double slopex,double slopey)
    {
        SETUP();
        const double sqrt6 = 2.44948974278;
        const double threeoverpi = 0.954929658551;
        return threeoverpi/fabs(sx*sy)*exp(-sqrt6*sqrt(sqr(slopex/sx)+sqr(slopey/sy)));
    };

    //
    // Bumpy_Slope_Distribution_Function definitions and code...
    //

    double
    Ellipsoid_Slope_Distribution_Function::
    f(double slopex,double slopey)
    {
        SETUP();
        return density*sqr(rx*ry*rz/(sqr(rz)+sqr(rx*slopex)+sqr(ry*slopey)));
    }

    //
    // Two_Slope_Distribution_Function definitions and code...
    //

    double
    Two_Slope_Distribution_Function::
    f(double slopex,double slopey)
    {
        SETUP();
        double sdf1 = s1->f(slopex,slopey);
        double sdf2 = s2->f(slopex,slopey);
        return sdf1*(1.-fract)+sdf2*fract;
    }

    //
    // Angle_Distribution_Function definitions and code...
    //
    double
    Angle_Distribution_Function::
    f(double slope)
    {
        SETUP();
        // Conversion from angle distribution to slope distribution...
        return P(atan(slope))/sqr(1.+sqr(slope));
    }


    //
    // Gaussian_Angle_Distribution_Function definitions and code...
    //

    void
    Gaussian_Angle_Distribution_Function::
    setup()
    {
        Angle_Distribution_Function::setup();
        sigmarad = sigma*deg;
        norm = 1.-exp(-sqr(pi/sigmarad/2.));
    }

    double
    Gaussian_Angle_Distribution_Function::
    P(double thetan)
    {
        SETUP();
        return  1./sqr(sigmarad)/pi*exp(-sqr(thetan/sigmarad))/norm;
    }



    void
    Exponential_Angle_Distribution_Function::
    setup()
    {
        Angle_Distribution_Function::setup();
        const double a = 3.84764949048;
        sigmarad = sigma*deg;
        norm = (sigmarad-exp(-a/sigmarad)*(a+sigmarad))/sigmarad;
    }

    double
    Exponential_Angle_Distribution_Function::
    P(double thetan)
    {
        SETUP();
        const double sqrt6 = 2.44948974278;
        const double threeoverpi = 0.954929658551;
        return threeoverpi/sqr(sigmarad)*exp(-sqrt6*thetan/sigmarad)/norm;
    }


    double
    Table_Angle_Distribution_Function::
    P(double thetan)
    {
        SETUP();
        return T.value(thetan*deg);
    }


    void Register(const Slope_Distribution_Function* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Slope_Distribution_Function);
            Register_Model(Unit_Slope_Distribution_Function);
            Register_Model(Exponential_Slope_Distribution_Function);
            Register_Model(Gaussian_Slope_Distribution_Function);
            Register_Model(Table_Slope_Distribution_Function);
            Register_Model(Anisotropic_Exponential_Slope_Distribution_Function);
            Register_Model(Ellipsoid_Slope_Distribution_Function);
            Register_Model(Angle_Distribution_Function);
            Register_Model(Gaussian_Angle_Distribution_Function);
            Register_Model(Exponential_Angle_Distribution_Function);
            Register_Model(Table_Angle_Distribution_Function);
            Register_Model(Two_Slope_Distribution_Function);
        }
    }

    DEFINE_MODEL(Facet_BRDF_Model,BRDF_Model,
                 "Scattering by a rough surface in the facet approximation.");

    DEFINE_VIRTUAL_MODEL(Slope_Distribution_Function,Model,
                         "All slope distribution functions");

    DEFINE_MODEL(Unit_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "Unit slope distribution function");

    DEFINE_MODEL(Exponential_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "Slope distribution proportional to exp(-slope/s)");

    DEFINE_MODEL(Gaussian_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "Slope distribution proportional to exp(-sqr(slope/s))");

    DEFINE_MODEL(Table_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "Slope distribution given by a table of values.");

    DEFINE_MODEL(Anisotropic_Exponential_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "An anisotropic SDF, with exponential behavior");



    DEFINE_PTRPARAMETER(Facet_BRDF_Model,Slope_Distribution_Function_Ptr,sdf,"Slope distribution function","Exponential_Slope_Distribution_Function",0xFF);

    DEFINE_PTRPARAMETER(Facet_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_PARAMETER(Exponential_Slope_Distribution_Function,double,s,"RMS slope","0.1",0xFF);

    DEFINE_PARAMETER(Gaussian_Slope_Distribution_Function,double,s,"RMS slope","0.1",0xFF);

    DEFINE_PARAMETER(Table_Slope_Distribution_Function,Table,T,"Slope distribution function","1",0xFF);

    DEFINE_PARAMETER(Anisotropic_Exponential_Slope_Distribution_Function,double,sx,"RMS slope in x direction","0.1",0xFF);
    DEFINE_PARAMETER(Anisotropic_Exponential_Slope_Distribution_Function,double,sy,"RMS slope in y direction","0.01",0xFF);

    DEFINE_MODEL(Ellipsoid_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "Ellipsoidal bumps or dips on surface");
    DEFINE_PARAMETER(Ellipsoid_Slope_Distribution_Function,double,rx,"Semiaxis of ellipsoid in x-direction [um]","10",0xFF);
    DEFINE_PARAMETER(Ellipsoid_Slope_Distribution_Function,double,ry,"Semiaxis of ellipsoid in y-direction [um]","10",0xFF);
    DEFINE_PARAMETER(Ellipsoid_Slope_Distribution_Function,double,rz,"Height or depth of bumps [um]","1",0xFF);
    DEFINE_PARAMETER(Ellipsoid_Slope_Distribution_Function,double,density,"Number density of ellipsoids [um^-2]","0.0001",0xFF);

    DEFINE_MODEL(Two_Slope_Distribution_Function,
                 Slope_Distribution_Function,
                 "Two-component slope distribution function");

    DEFINE_PTRPARAMETER(Two_Slope_Distribution_Function,Slope_Distribution_Function_Ptr,s1,"First slope distribution function","Exponential_Slope_Distribution_Function",0xFF);
    DEFINE_PTRPARAMETER(Two_Slope_Distribution_Function,Slope_Distribution_Function_Ptr,s2,"Second slope distribution function","Exponential_Slope_Distribution_Function",0xFF);
    DEFINE_PARAMETER(Two_Slope_Distribution_Function,double,fract,"Relative fraction of second to first (0<=f<=1)","0.5",0xFF);

    DEFINE_VIRTUAL_MODEL(Angle_Distribution_Function,
                         Slope_Distribution_Function,
                         "Slope distribution function which defines the distribution of surface normal angles");

    DEFINE_MODEL(Gaussian_Angle_Distribution_Function,
                 Angle_Distribution_Function,
                 "An orientation distribution which is normal in angle");
    DEFINE_PARAMETER(Gaussian_Angle_Distribution_Function,double,sigma,"Angle standard deviation [deg]","5",0xFF);

    DEFINE_MODEL(Exponential_Angle_Distribution_Function,
                 Angle_Distribution_Function,
                 "An orientation distribution which is exponential in angle");
    DEFINE_PARAMETER(Exponential_Angle_Distribution_Function,double,sigma,"Angle standard deviation [deg]","5",0xFF);

    DEFINE_MODEL(Table_Angle_Distribution_Function,
                 Angle_Distribution_Function,
                 "An orientation distribution defined by a table");
    DEFINE_PARAMETER(Table_Angle_Distribution_Function,Table,T,"Angle distribution function","1",0xFF);


} // namespace SCATMECH



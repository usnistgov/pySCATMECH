//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: torrspar.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "torrspar.h"
#include "vector3d.h"

using namespace std;


namespace SCATMECH {


    //
    // Scattering matrix for Shadowed_Facet_BRDF_Model...
    //
    JonesMatrix
    Shadowed_Facet_BRDF_Model::
    jones()
    {
        // This program implements the model outlined in
        // K.E.Torrance and E.M.Sparrow, "Theory of off-specular reflection
        // from roughened surfaces,"
        // J. Opt. Soc. Am. vol. 57, no. 9, pp. 1105-1114 (1967).

        // All BRDF_Model::jones(,,) should call setup if recalc is set...
        SETUP();

        double s = shadow->f(thetai,thetas,phis);

        // The scattering matrix is the Facet_BRDF_Model scattering
        // matrix multiplied by the T&S shadow function...
        return Facet_BRDF_Model::jones()*
               (COMPLEX)(sqrt(s));
    }

    static
    double min3(double x,double y, double z)
    {
        if (x<y) {
            if (x<z) return x;
            else return z;
        } else {
            if (y<z) return y;
            else return z;
        }
    }

    double
    Torrance_Sparrow_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        // The article by Torrance and Sparrow gives a more complicated
        // version of the shadowing function.  It is significantly
        // simplified in R.L. Cook and K.E. Torrance, "A Reflectance
        // Model for Computer Graphics,"
        // Computer Graphics vol. 15, no. 3, pp. 307-316 (1981)
        //
        // In SCATMECH Version 3.00, we use the simplified expressions...
        //
        double cos_thetas = cos(thetas);
        double cos_thetai = cos(thetai);
        double sin_thetas = sin(thetas);
        double sin_thetai = sin(thetai);
        double cos_phis = cos(phis);
        double sin_phis = sin(phis);

        Vector in(-sin_thetai,0.,cos_thetai);
        Vector out(sin_thetas*cos_phis,sin_thetas*sin_phis,cos_thetas);
        Vector normal = unit(in+out);
        double cos_delta = normal.z;
        double cos_beta = normal*in;
        double temp=2.*cos_delta/cos_beta;
        return min3(1.,temp*cos_thetai,temp*cos_thetas);
    }

    double
    Table_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        return T.value(thetai)*T.value(thetas);
    }

    double
    Smith_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        SETUP();

        return S(thetai)*S(thetas);
    }

    double
    Smith_Shadow_Function::
    S(double theta)
    {
        if (theta!=0.) {
            const double sqrt2 = 1.4142135623730951;
            const double sqrttowoverpi = 0.7978845608028654;
            double tant = tan(fabs(theta));
            double mu = 1./tant;
            double Lambda = 0.5*(sqrttowoverpi*w/mu*exp(-sqr(mu/w)/2)-erfc(mu/w/sqrt2));
            return (1.-0.5*erfc(mu/w/sqrt2))/(Lambda+1.);
        } else {
            return 1.;
        }
    }


    double
    Maxwell_Beard_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        SETUP();

        double temp = sqr(sin(thetai))-2.*sin(thetai)*sin(thetas)*cos(phis)+
                      sqr(sin(thetas));

        double localslope = (temp>0.) ? sqrt(temp)/(cos(thetai)+cos(thetas)) : 0.;

        double thetan = atan(localslope);

        double cos_angle = sqrt((1.-sin(thetai)*sin(thetas)*cos(phis)+
                                 cos(thetai)*cos(thetas))/2.);
        double beta = (cos_angle<1.) ? acos(cos_angle) : 0.;

        return (1.+(thetan/Omega)*exp(-2*beta/tau))/(1.+(thetan/Omega));
    }

    void Register(const Shadow_Function* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Shadow_Function);
            Register_Model(Unit_Shadow_Function);
            Register_Model(Torrance_Sparrow_Shadow_Function);
            Register_Model(Smith_Shadow_Function);
            Register_Model(Maxwell_Beard_Shadow_Function);
            Register_Model(Table_Shadow_Function);
        }
    }

    DEFINE_VIRTUAL_MODEL(Shadow_Function,Model,
                         "Shadow Function");

    DEFINE_MODEL(Unit_Shadow_Function,Shadow_Function,
                 "Unit Shadow Function");

    DEFINE_MODEL(Torrance_Sparrow_Shadow_Function,Shadow_Function,
                 "Shadow function appropriate for V-grooves");

    DEFINE_MODEL(Shadowed_Facet_BRDF_Model,Facet_BRDF_Model,
                 "The facet scattering model with a shadow function.");

    DEFINE_MODEL(Smith_Shadow_Function,Shadow_Function,
                 "Shadow function developed by Smith appropriate for gaussian rough surfaces");
    DEFINE_PARAMETER(Smith_Shadow_Function,double,w,"RMS surface slope","0.35",0xFF);

    DEFINE_MODEL(Maxwell_Beard_Shadow_Function,Shadow_Function,
                 "Empirical shadow function developed by Maxwell-Beard");
    DEFINE_PARAMETER(Maxwell_Beard_Shadow_Function,double,Omega,"Omega [rad]","0.7",0xFF);
    DEFINE_PARAMETER(Maxwell_Beard_Shadow_Function,double,tau,"tau [rad]","0.25",0xFF);

    DEFINE_MODEL(Table_Shadow_Function,Shadow_Function,
                 "Shadow function from a table.");

    DEFINE_PTRPARAMETER(Shadowed_Facet_BRDF_Model,Shadow_Function_Ptr,shadow,"Shadow function","Torrance_Sparrow_Shadow_Function",0xFF);

    DEFINE_PARAMETER(Table_Shadow_Function,Table,T,"Tabulated shadow function","1",0xFF);


} // namespace SCATMECH


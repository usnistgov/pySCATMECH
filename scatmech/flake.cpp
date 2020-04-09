//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: flake.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "flake.h"
#include "vector3d.h"
#include "matrix3d.h"
#include "fresnel.h"

using namespace std;


namespace SCATMECH {


    double
    Subsurface_Facet_BRDF_Model::
    local_angle(double thetai,double thetas,double phis)
    {
        double _thetai = asin(sin(thetai)/overcoat.n(lambda));
        double _thetas = asin(sin(thetas)/overcoat.n(lambda));

        return Facet_BRDF_Model::local_angle(_thetai,_thetas,phis);
    }

    double
    Subsurface_Facet_BRDF_Model::
    local_slope(double thetai,double thetas,double phis)
    {
        double _thetai = asin(sin(thetai)/overcoat.n(lambda));
        double _thetas = asin(sin(thetas)/overcoat.n(lambda));

        return Facet_BRDF_Model::local_slope(_thetai,_thetas,phis);
    }


    JonesMatrix
    Subsurface_Facet_BRDF_Model::
    jones()
    {
        SETUP();

        throw_transmission();
        throw_backward();

        if (thetai<0) {
            thetai=-thetai;
            phis = pi+phis;
        }

        if (thetas<0) {
            thetas=-thetas;
            phis=pi+phis;
        }

        dielectric_constant epsilon = overcoat.epsilon(lambda);
        double e = overcoat.e1(lambda);
        double n = overcoat.n(lambda);

        // Other useful terms...
        double cos_phis = cos(phis);
        double sin_phis = sin(phis);
        double cos_thetai = cos(thetai);
        double sin_thetai = sin(thetai);
        double cos_thetas = cos(thetas);
        double sin_thetas = sin(thetas);

        // Angles of refraction internal to the material...
        double sin_thetai_internal = sin_thetai/n;
        double cos_thetai_internal = sqrt(1.0-sqr(sin_thetai_internal));
        double sin_thetas_internal = sin_thetas/n;
        double cos_thetas_internal = sqrt(1.0-sqr(sin_thetas_internal));
        double thetai_internal = acos(cos_thetai_internal);
        double thetas_internal = acos(cos_thetas_internal);

        // Internal angle of incidence on flake...
        double iota = local_angle(thetai,thetas,phis);

        // Facet scattering coefficients...
        double a1 = sqr(sin(2*iota));
        double a2 = cos_thetai_internal*sin_thetas_internal+
                    sin_thetai_internal*cos_thetas_internal*cos_phis;
        double a3 = sin_thetai_internal*cos_thetas_internal+
                    cos_thetai_internal*sin_thetas_internal*cos_phis;

        // We need the vacuum "external" angle for the flake reflection
        // coeficients...
        double iota_external = asin(sin(iota)*n);

        // The reflection coefficients are obtained from the
        // dielectric stack...
        COMPLEX rsiota = stack->rs12(iota_external,lambda,overcoat,substrate);
        COMPLEX rpiota = stack->rp12(iota_external,lambda,overcoat,substrate);

        // The local surface normal...
        Vector _nhat = nhat(thetai_internal,thetas_internal,phis);

        // The slope of the flake accessed by the geometry...
        double slope = sqrt(sqr(_nhat.x)+sqr(_nhat.y))/_nhat.z;

        // Get the slope distribution function and prefactor...
        COMPLEX s=pi*sqr(1.+sqr(slope))* evaluate_sdf(_nhat);

        // Radar cross section to BRDF conversion factor...
        s *= 1./(4.*pi*cos_thetai*cos_thetas);

        // The polarimetric part of the calculation...
        JonesMatrix j;
        if (fabs(iota) > 1E-5) {
            j[0] = (sin_thetai_internal*sin_thetas_internal*sqr(sin_phis)*
                    rsiota+ a2*a3*rpiota)/a1;
            j[1] = (a2*a3*rsiota+ sin_thetai_internal*sin_thetas_internal*
                    sqr(sin_phis)*rpiota)/a1;
            j[2] = -sin_phis*(sin_thetas_internal*a2*rsiota
                              - sin_thetai_internal*a3*rpiota)/a1;
            j[3] = -sin_phis*(sin_thetai_internal*a3*rsiota
                              -sin_thetas_internal*a2*rpiota)/a1;
        } else {
            j[0] = rpiota;
            j[1] = rsiota;
            j[2] = 0.;
            j[3] = 0.;
        }

        // The transmission matrices into and out of the overcoat...
        JonesMatrix jti(overcoat_films->tp12(thetai,lambda,vacuum,epsilon),
                        overcoat_films->ts12(thetai,lambda,vacuum,epsilon),
                        (COMPLEX)(0.),(COMPLEX)(0.));

        JonesMatrix jts(overcoat_films->tp12(thetas,lambda,vacuum,epsilon),
                        overcoat_films->ts12(thetas,lambda,vacuum,epsilon),
                        (COMPLEX)(0.),(COMPLEX)(0.));

        // Total Jones matrix
        j = jts*(j*jti);

        j = j*sqrt(s);

        return j;
    }

    void
    Subsurface_Facet_BRDF_Model::
    setup()
    {
        Facet_BRDF_Model::setup();
        overcoat.force_nonabsorbing();
    }


    DEFINE_MODEL(Subsurface_Facet_BRDF_Model,Facet_BRDF_Model,
                 "Facet scattering model applied to an interface under a smooth overlayer.");

    DEFINE_PARAMETER(Subsurface_Facet_BRDF_Model,dielectric_function,overcoat,"Overcoat","(1.59,0)",0xFF);

	DEFINE_PTRPARAMETER(Subsurface_Facet_BRDF_Model,StackModel_Ptr,overcoat_films,"Overcoat films","No_StackModel",0xFF);


} // namespace SCATMECH



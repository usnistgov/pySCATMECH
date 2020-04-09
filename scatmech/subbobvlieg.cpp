//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: subbobvlieg.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "subbobvlieg.h"

using namespace std;

namespace SCATMECH {

    void
    Subsurface_Bobbert_Vlieger_BRDF_Model::
    setup()
    {
        Local_BRDF_Model::setup();

        substrate.force_nonabsorbing();

        double n = substrate.n(lambda);

        model.set_lambda(lambda/n);
        int backtype;
        switch (type) {
            case 0:
                backtype = 2;
                break;
            case 1:
                backtype = 3;
                break;
            case 2:
                backtype = 0;
                break;
            case 3:
                backtype = 1;
                break;
            default:
                error("Invalid type");
        }
        model.set_type(backtype);
        model.set_substrate((optical_constant)(1./n));
        model.set_density(density);
        model.set_sphere((optical_constant)((COMPLEX)(sphere.index(lambda))/n));
        model.set_radius(radius);
        dielectric_stack temp;
        for (int i=0; i<spherecoat->get_n(); ++i) {
            int j = i;
            COMPLEX index = (COMPLEX)(spherecoat->get_e()[j].index(lambda))/n;
            double thickness = spherecoat->get_t()[j];
            temp.grow(optical_constant(index),thickness);
        }
		Stack_StackModel temp0;
		temp0.set_stack(temp);
        model.set_spherecoat(temp0.clone());
        temp.wash();
        for (int i=0; i<stack->get_n(); ++i) {
            // Need to reverse order of films...
            int j = stack->get_n()-i-1;
            COMPLEX index = (COMPLEX)(stack->get_e()[j].index(lambda))/n;
            double thickness = stack->get_t()[j];
            temp.grow(optical_constant(index),thickness);
        }
		temp0.set_stack(temp);
        model.set_stack(temp0.clone());
        model.set_delta(delta);
        model.set_lmax(lmax);
        model.set_order(order);
        model.set_Norm_Inc_Approx(Norm_Inc_Approx);
        model.set_improve(improve);
    }

    JonesMatrix
    Subsurface_Bobbert_Vlieger_BRDF_Model::
    jonesDSC()
    {
        SETUP();

        JonesMatrix j = model.JonesDSC(thetai,thetas,phis,rotation);
        return JonesMatrix(j[0],j[1],-j[2],-j[3]);
    }

    void
    Subsurface_Axisymmetric_Particle_BRDF_Model::
    setup()
    {
        Local_BRDF_Model::setup();

        substrate.force_nonabsorbing();

        double n = substrate.n(lambda);

        model.set_lambda(lambda/n);
        int backtype;
        switch (type) {
            case 0:
                backtype = 2;
                break;
            case 1:
                backtype = 3;
                break;
            case 2:
                backtype = 0;
                break;
            case 3:
                backtype = 1;
                break;
            default:
                error("Invalid type");
        }
        model.set_type(backtype);
        model.set_substrate((optical_constant)(1./n));
        model.set_density(density);
        model.set_particle((optical_constant)((COMPLEX)(particle.index(lambda))/n));
        model.set_Shape(Shape);
        dielectric_stack temp;
        for (int i=0; i<stack->get_n(); ++i) {
            // Need to reverse order of films...
            int j = stack->get_n()-i-1;
            COMPLEX index = (COMPLEX)(stack->get_e()[j].index(lambda))/n;
            double thickness = stack->get_t()[j];
            temp.grow(optical_constant(index),thickness);
        }
		Stack_StackModel temp0;
		temp0.set_stack(temp);
        model.set_stack(temp0.clone());
        model.set_delta(delta);
        model.set_lmax(lmax);
        model.set_mmax(mmax);
        model.set_order(order);
        model.set_Norm_Inc_Approx(Norm_Inc_Approx);
        model.set_improve(improve);
    }

    JonesMatrix
    Subsurface_Axisymmetric_Particle_BRDF_Model::
    jonesDSC()
    {
        SETUP();

        JonesMatrix j = model.JonesDSC(thetai,thetas,phis,rotation);
        return JonesMatrix(j[0],j[1],-j[2],-j[3]);
    }


    DEFINE_MODEL(Subsurface_Bobbert_Vlieger_BRDF_Model,Local_BRDF_Model,
                 "Theory for scattering by a sphere in a substrate.");

    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,dielectric_function,sphere,"Sphere optical properties","(1.59,0)",0xFF);
    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,double,radius,"Particle radius [um]","0.05",0xFF);
    DEFINE_PTRPARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,StackModel_Ptr,spherecoat,"Coatings on the sphere","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,StackModel_Ptr,stack,"Substrate films","No_StackModel",0xFF);
    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,double,delta,"Separation of particle from substrate [um] (in contact: 0)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,int,lmax,"Maximum spherical harmonic order (use Bohren & Huffman estimate: 0)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,int,order,"Perturbation order (exact: -1)","-1",0xFF);
    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,int,Norm_Inc_Approx,"Normal Incidence Approximation (exact: 0)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Bobbert_Vlieger_BRDF_Model,int,improve,"Iterative improvement steps (recommend: 3)","3",0xFF);

    DEFINE_MODEL(Subsurface_Axisymmetric_Particle_BRDF_Model,Local_BRDF_Model,
                 "Axisymmetric particle on a substrate.");

    DEFINE_PTRPARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,Axisymmetric_Shape_Ptr,Shape,"Particle Shape","Ellipsoid_Axisymmetric_Shape",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,dielectric_function,particle,"Particle optical properties","(1.59,0)",0xFF);
    DEFINE_PTRPARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,StackModel_Ptr,stack,"Substrate films","No_StackModel",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,double,delta,"Separation of particle from substrate [um] (in contact: 0)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,int,lmax,"Maximum polar order (lmax)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,int,mmax,"Maximum azimuthal order (mmax)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,int,order,"Perturbation order (exact: -1)","-1",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,int,Norm_Inc_Approx,"Normal incidence approximation (exact: 0)","0",0xFF);
    DEFINE_PARAMETER(Subsurface_Axisymmetric_Particle_BRDF_Model,int,improve,"Iterative improvement steps (recommend: 3)","3",0xFF);


}

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crough.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "crough.h"

using namespace std;


namespace SCATMECH {

    void
    Correlated_Roughness_BRDF_Model::
    setup()
    {
        Roughness_BRDF_Model::setup();

        model.set_lambda(lambda);
        model.set_type(type);
        model.set_substrate(substrate);
        model.set_psd(psd);
        SingleFilm_StackModel stack;
		stack.set_material(film);
		stack.set_thickness(thickness);
        model.set_stack(stack.clone());
    }
    //
    // Jones matrix for scattering...
    //
    JonesMatrix
    Correlated_Roughness_BRDF_Model::
    jones()
    {
        SETUP();

        return model.Jones(thetai,thetas,phis,rotation);
    }

    DEFINE_MODEL(Correlated_Roughness_BRDF_Model,Roughness_BRDF_Model,
                 "Scattering from a single dielectric layer on a substrate with equal and correlated roughness.");

    DEFINE_PARAMETER(Correlated_Roughness_BRDF_Model,double,thickness,"Film thickness [um]","0.05",0xFF);

    DEFINE_PARAMETER(Correlated_Roughness_BRDF_Model,dielectric_function,film,"Overlayer","(1.46,0)",0xFF);


} // namespace SCATMECH




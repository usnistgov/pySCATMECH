//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: two_source.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "two_source.h"

using namespace std;

namespace SCATMECH {

    MuellerMatrix Two_Source_BRDF_Model::mueller()
    {
        SETUP();

        if (lambda!=source1->get_lambda()) error("lambda!=source1.lambda");
        if (type!=source1->get_type()) error("type!=source1.type");
        if (substrate.index(lambda)!=source1->get_substrate().index(lambda)) error("substrate!=source1.substrate");

        if (lambda!=source2->get_lambda()) error("lambda!=source2.lambda");
        if (type!=source2->get_type()) error("type!=source2.type");
        if (substrate.index(lambda)!=source2->get_substrate().index(lambda)) error("substrate!=source2.substrate");

        if (factor1<0.) error("factor1<0");
        if (factor2<0.) error("factor2<0");

        if (correlation==0) {
            return source1->Mueller(thetai,thetas,phis,rotation)*factor1+
                   source2->Mueller(thetai,thetas,phis,rotation)*factor2;
        }

        JonesMatrix j1 = source1->Jones(thetai,thetas,phis,rotation)*sqrt(factor1);
        JonesMatrix j2 = source2->Jones(thetai,thetas,phis,rotation)*sqrt(factor2);

        return (1-correlation)*(MuellerMatrix(j1)+MuellerMatrix(j2)) +correlation*(MuellerMatrix(j1+j2));
    }

    DEFINE_MODEL(Two_Source_BRDF_Model,BRDF_Model,"Sum of two BRDF models");
    DEFINE_PARAMETER(Two_Source_BRDF_Model,double,factor1,"Intensity scale factor for first scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Two_Source_BRDF_Model,BRDF_Model_Ptr,source1,"First scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Two_Source_BRDF_Model,double,factor2,"Intensity scale factor for second scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Two_Source_BRDF_Model,BRDF_Model_Ptr,source2,"Second scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Two_Source_BRDF_Model,double,correlation,"Correlation coefficient","0",0xFF);

    MuellerMatrix Three_Source_BRDF_Model::mueller()
    {
        SETUP();

        if (lambda!=source1->get_lambda()) error("lambda!=source1.lambda");
        if (type!=source1->get_type()) error("type!=source1.type");
        if (substrate.index(lambda)!=source1->get_substrate().index(lambda)) error("substrate!=source1.substrate");

        if (lambda!=source2->get_lambda()) error("lambda!=source2.lambda");
        if (type!=source2->get_type()) error("type!=source2.type");
        if (substrate.index(lambda)!=source2->get_substrate().index(lambda)) error("substrate!=source2.substrate");

        if (lambda!=source3->get_lambda()) error("lambda!=source3.lambda");
        if (type!=source3->get_type()) error("type!=source3.type");
        if (substrate.index(lambda)!=source3->get_substrate().index(lambda)) error("substrate!=source3.substrate");

        if (factor1<0.) error("factor1<0");
        if (factor2<0.) error("factor2<0");
        if (factor3<0.) error("factor3<0");

        return source1->Mueller(thetai,thetas,phis,rotation)*factor1+
               source2->Mueller(thetai,thetas,phis,rotation)*factor2+
               source3->Mueller(thetai,thetas,phis,rotation)*factor3;
    }

    DEFINE_MODEL(Three_Source_BRDF_Model,BRDF_Model,"Sum of three BRDF models");
    DEFINE_PARAMETER(Three_Source_BRDF_Model,double,factor1,"Intensity scale factor for first scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Three_Source_BRDF_Model,BRDF_Model_Ptr,source1,"First scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Three_Source_BRDF_Model,double,factor2,"Intensity scale factor for second scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Three_Source_BRDF_Model,BRDF_Model_Ptr,source2,"Second scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Three_Source_BRDF_Model,double,factor3,"Intensity scale factor for third scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Three_Source_BRDF_Model,BRDF_Model_Ptr,source3,"Third scattering source","Microroughness_BRDF_Model",0xFF);

    MuellerMatrix Four_Source_BRDF_Model::mueller()
    {
        SETUP();

        if (lambda!=source1->get_lambda()) error("lambda!=source1.lambda");
        if (type!=source1->get_type()) error("type!=source1.type");
        if (substrate.index(lambda)!=source1->get_substrate().index(lambda)) error("substrate!=source1.substrate");

        if (lambda!=source2->get_lambda()) error("lambda!=source2.lambda");
        if (type!=source2->get_type()) error("type!=source2.type");
        if (substrate.index(lambda)!=source2->get_substrate().index(lambda)) error("substrate!=source2.substrate");

        if (lambda!=source3->get_lambda()) error("lambda!=source3.lambda");
        if (type!=source3->get_type()) error("type!=source3.type");
        if (substrate.index(lambda)!=source3->get_substrate().index(lambda)) error("substrate!=source3.substrate");

        if (lambda!=source4->get_lambda()) error("lambda!=source4.lambda");
        if (type!=source4->get_type()) error("type!=source4.type");
        if (substrate.index(lambda)!=source4->get_substrate().index(lambda)) error("substrate!=source4.substrate");

        if (factor1<0.) error("factor1<0");
        if (factor2<0.) error("factor2<0");
        if (factor3<0.) error("factor3<0");
        if (factor4<0.) error("factor4<0");

        return source1->Mueller(thetai,thetas,phis,rotation)*factor1+
               source2->Mueller(thetai,thetas,phis,rotation)*factor2+
               source3->Mueller(thetai,thetas,phis,rotation)*factor3+
               source4->Mueller(thetai,thetas,phis,rotation)*factor4;
    }

    DEFINE_MODEL(Four_Source_BRDF_Model,BRDF_Model,"Sum of four BRDF models");
    DEFINE_PARAMETER(Four_Source_BRDF_Model,double,factor1,"Intensity scale factor for first scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Four_Source_BRDF_Model,BRDF_Model_Ptr,source1,"First scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Four_Source_BRDF_Model,double,factor2,"Intensity scale factor for second scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Four_Source_BRDF_Model,BRDF_Model_Ptr,source2,"Second scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Four_Source_BRDF_Model,double,factor3,"Intensity scale factor for third scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Four_Source_BRDF_Model,BRDF_Model_Ptr,source3,"Third scattering source","Microroughness_BRDF_Model",0xFF);
    DEFINE_PARAMETER(Four_Source_BRDF_Model,double,factor4,"Intensity scale factor for fourth scattering source","1",0xFF);
    DEFINE_PTRPARAMETER(Four_Source_BRDF_Model,BRDF_Model_Ptr,source4,"Fourth scattering source","Microroughness_BRDF_Model",0xFF);

}

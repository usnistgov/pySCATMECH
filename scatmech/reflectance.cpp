//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reflectance.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "reflectance.h"
#include "scateval.h"

using namespace std;

namespace SCATMECH {

    void
    Kubelka_Munk_Reflectance::
    setup()
    {
        //Evaluator::VMAP vars;
        //vars["x"]=lambda;
        double s = absorption.value(lambda);
        double r = scattering.value(lambda);

        //double s = Evaluator(absorbance,vars);
        //double r = Evaluator(scattering,vars);
        double Hprime = substrate->Get_Reflectance(lambda);
        if (s==0) {
            double temp1 = (1-Hprime)*r*thickness;
            reflectance = (temp1+Hprime)/(temp1+1);
        } else if (r == 0) {
            reflectance = Hprime*exp(-2*s*thickness);
        } else {
            double a = 1+s/r;
            double Hinf = a - sqrt(a*a-1.);
            double temp1 = 1/Hinf-Hinf;
            double temp2 = Hprime-Hinf;
            double temp3 = Hprime-1/Hinf;
            double temp4 = exp(r*thickness*temp1);
            double H = (temp2/Hinf-Hinf*temp3*temp4)/
                       (temp2-temp3*temp4);
            reflectance = H;
        }
    };

    void
    Table_Reflectance::
    setup()
    {
        reflectance = table.value(lambda);
    }

    void
    Equation_Reflectance::
    setup()
    {
        Evaluator::VMAP vars;
        vars["x"] = lambda;
        reflectance = Evaluator(expression,vars);
    }

    void
    Register(const Reflectance* x)
    {
        static bool regd=false;
        if (!regd) {
            regd=true;

            Register_Model(Reflectance);
            Register_Model(Table_Reflectance);
            Register_Model(Equation_Reflectance);
            Register_Model(Kubelka_Munk_Reflectance);
        }
    }

    DEFINE_VIRTUAL_MODEL(Reflectance,Model,"Generalized model for reflectance");

    DEFINE_MODEL(Kubelka_Munk_Reflectance,Reflectance,"Kubelka-Munk reflectance for a paint with pigments");
    DEFINE_PARAMETER(Kubelka_Munk_Reflectance,double,thickness,"Thickness of the coating [um]","1",0xFF);
    DEFINE_PARAMETER(Kubelka_Munk_Reflectance,Table,absorption,"Absorption coefficient function [1/um]","1",0xFF);
    DEFINE_PARAMETER(Kubelka_Munk_Reflectance,Table,scattering,"Scattering coefficient function [1/um]","1",0xFF);
    DEFINE_PTRPARAMETER(Kubelka_Munk_Reflectance,Reflectance_Ptr,substrate,"Reflectance of the substrate","Table_Reflectance",0xFF);

    DEFINE_MODEL(Table_Reflectance,Reflectance,"Reflectance as a function of wavelength taken from a table");
    DEFINE_PARAMETER(Table_Reflectance,Table,table,"Table of reflectance","1",0xFF);

    DEFINE_MODEL(Equation_Reflectance,Reflectance,"Reflectance as a function of wavelength from an expression");
    DEFINE_PARAMETER(Equation_Reflectance,string,expression,"Expression for reflectance (function of x)","1",0xFF);

}

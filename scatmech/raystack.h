//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: raystack.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_STACK_BR_H
#define SCATMECH_STACK_BR_H

#include "scatmech.h"
#include "filmtran.h"
#include "local.h"

namespace SCATMECH {

    //
    // Stack_BRDF_Model is a BRDF_Model for a defect in a dielectric stack.
    //
    class Rayleigh_Stack_BRDF_Model : public Local_BRDF_Model
    {
        public:
            // Routine to perform housekeeping...
            virtual void setup();

            // Jones matrix for scattering...
            SCATMECH::JonesMatrix jonesDSC();

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,depth);
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(dielectric_function,sphere);

        protected:

            // Normalized distance of defect from upper surface...
            //double kd;
            double d;
            double lambda_eff;

            // The optical constant of the layer containg the defect...
            dielectric_function material1;
            dielectric_function material2;
            dielectric_function material3;
            dielectric_function materialS;

            // Normalized distance of defect from defect to upper layer
            //double ktau;
            double tau;

            // The dielectric stack below the defect...
            dielectric_stack below;

            // The dielectric stack above the defect...
            dielectric_stack above;


    };

}

#endif

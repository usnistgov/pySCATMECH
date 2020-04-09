//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: subsphere.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_SUBSURFACE_PARTICLE_H
#define SCATMECH_SUBSURFACE_PARTICLE_H

#include "local.h"
#include "sphrscat.h"
#include "filmtran.h"


namespace SCATMECH {


    class Subsurface_Particle_BRDF_Model : public Local_BRDF_Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,depth);
            DECLARE_PARAMETER(Free_Space_Scatterer_Ptr,scatterer);
            DECLARE_PARAMETER(StackModel_Ptr,stack);

            // Particle can be rotated from its natural orientation...
            // Rotation of particle about z-axis...
            DECLARE_PARAMETER(double,alpha);
            // Rotation of particle about y-axis...
            DECLARE_PARAMETER(double,beta);
            // Rotation of particle about z-axis again is performed by brdf variable "rotation"


        protected:
            virtual void setup();
            virtual JonesMatrix jonesDSC();
            Matrix Euler;
    };


} // namespace SCATMECH


#endif


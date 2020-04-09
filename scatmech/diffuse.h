//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: diffuse.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_DIFFUSE_H
#define SCATMECH_DIFFUSE_H

#include "filmtran.h"
#include "lambert.h"


namespace SCATMECH {

    class Diffuse_Subsurface_BRDF_Model: public Lambertian_BRDF_Model
    {
        public:

            Diffuse_Subsurface_BRDF_Model();

        public:

            DECLARE_MODEL();

            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:
            MuellerMatrix mueller();
            void setup();

        private:
            double feedback;
    };


} // namespace SCATMECH




#endif

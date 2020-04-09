//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: lambert.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_LAMBERT_H
#define SCATMECH_LAMBERT_H

#include "brdf.h"
#include "scattabl.h"
#include "reflectance.h"

namespace SCATMECH {



    class Lambertian_BRDF_Model : public BRDF_Model {
        public:
            DECLARE_MODEL();

            // The total reflectance of the sample
            //DECLARE_PARAMETER(Table,reflectance);
            DECLARE_PARAMETER(Reflectance_Ptr,reflectance);
        protected:

            MuellerMatrix mueller();
    };

    typedef Model_Ptr<Lambertian_BRDF_Model> Lambertian_BRDF_Model_Ptr;

    void Register(const Lambertian_BRDF_Model* x);


} // namespace SCATMECH


#endif

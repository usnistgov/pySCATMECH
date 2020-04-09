//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: local.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_LOCAL_H
#define SCATMECH_LOCAL_H

#include <complex>
#include "brdf.h"


namespace SCATMECH {


    class Local_BRDF_Model : public BRDF_Model {
        public:
            MuellerMatrix MuellerDSC(double thetai,
                                     double thetas,
                                     double phis,
                                     double rotation,
                                     Coordinate_System cs = psps);

            JonesMatrix JonesDSC(double thetai,
                                 double thetas,
                                 double phis,
                                 double rotation,
                                 Coordinate_System cs = psps);

            // Jones matrix for scattering where the directions to the source, the viewer,
            // and the surface normal are specified by Vectors...
            JonesMatrix JonesDSC(const Vector& source,
                                 const Vector& viewer,
                                 const Vector& normal,
                                 const Vector& xaxis,
                                 Coordinate_System cs=plane);

            // Mueller matrix for scattering where the directions to the source, the viewer,
            // and the surface normal are specified by Vectors...
            MuellerMatrix MuellerDSC(const Vector& source,
                                     const Vector& viewer,
                                     const Vector& normal,
                                     const Vector& xaxis,
                                     Coordinate_System cs=plane);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,density);

        protected:

            virtual void setup();

            virtual JonesMatrix jonesDSC() {
                error("Attempt to call jonesDSC for a depolarizing model.");
                return JonesMatrix();
            }

            virtual MuellerMatrix muellerDSC() {
                return this->jonesDSC();
            }

        private:

            JonesMatrix jones() {
                return this->jonesDSC()*sqrt(density/(fabs(cos(thetai)*cos(thetas))));
            }
            MuellerMatrix mueller() {
                return this->muellerDSC()*(density/(fabs(cos(thetai)*cos(thetas))));
            }
    };

    typedef Model_Ptr<Local_BRDF_Model> Local_BRDF_Model_Ptr;

    void Register(const Local_BRDF_Model* x);


} // namespace SCATMECH


#endif

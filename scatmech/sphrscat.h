//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sphrscat.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_SPHRSCAT_H
#define SCATMECH_SPHRSCAT_H

#include "scatmech.h"
#include "optconst.h"
#include "mueller.h"
#include "dielfunc.h"
#include "inherit.h"
#include "vector3d.h"
#include "matrix3d.h"


namespace SCATMECH {


    class Free_Space_Scatterer : public Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,lambda);
            DECLARE_PARAMETER(dielectric_function,medium);

        public:
            //
            // Generalized scattering matrix...
            //
            // The coordinate system for the polarization is parallel and perpendicular to the plane of scattering.
            //
            // The Vectors kin and kout are unit vectors in the incident and scattering directions, respectively.
            //
            virtual JonesMatrix jones(const Vector& kin,const Vector& kout) = 0;

			// 
			// Exctinction cross section...
			//
			MuellerMatrix extinction(const Vector& k);
    };

    typedef Model_Ptr<Free_Space_Scatterer> Free_Space_Scatterer_Ptr;

    //
    // SphericalScatterer is a virtual base class for the scattering
    // from any spherically symmetric scatterer.
    //
    class SphericalScatterer : public Free_Space_Scatterer
    {
        public:
            // Two COMPLEX functions characterize a spherically-symmetric
            // scatterer...
            virtual COMPLEX s1(double angle)=0;
            virtual COMPLEX s2(double angle)=0;

            // Result can also be returned as a Jones scattering matrix...
            JonesMatrix s(double angle) {
                return JonesMatrix(s1(angle),s2(angle),(COMPLEX)0.,(COMPLEX)0.);
            }

            JonesMatrix jones(const Vector& kin,const Vector& kout) {
                double dprod = kin*kout;
                if (dprod>1.) dprod = 1.;
                if (dprod<-1.) dprod = -1.;
                double theta = acos(dprod);
                return s(theta);
            }

            virtual double Csca()=0;
            virtual double Cext()=0;
            virtual double Cback()=0;
            double Cabs() {
                return Cext()-Csca();
            }

            double Qsca() {
                return Csca()/(pi*sqr(radius));
            }
            double Qext() {
                return Cext()/(pi*sqr(radius));
            }
            double Qback() {
                return Cback()/(pi*sqr(radius));
            }
            double Qabs() {
                return Cabs()/(pi*sqr(radius));
            }

        public:
            DECLARE_MODEL();

            DECLARE_PARAMETER(dielectric_function,sphere);
            DECLARE_PARAMETER(double,radius);

        protected:

            void setup() {
                Free_Space_Scatterer::setup();
                nmed = medium.n(lambda);
                nsphere = sphere.index(lambda);
                k = 2.*pi*nmed/lambda;
                x = k*radius;
                m = nsphere/nmed;
            }

            double nmed;
            COMPLEX nsphere;
            COMPLEX m;
            double x;
            double k;
    };

    typedef Model_Ptr<SphericalScatterer> SphericalScatterer_Ptr;

    void Register(const Free_Space_Scatterer* x);
    void Register(const SphericalScatterer* x);


} // namespace SCATMECH


#endif

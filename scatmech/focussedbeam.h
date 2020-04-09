//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: focussedbeam.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_FOCUSSEDBEAM_H
#define SCATMECH_FOCUSSEDBEAM_H

#include "scatmech.h"
#include "instrument.h"
#include "inherit.h"
#include "vector3d.h"


namespace SCATMECH {


    class Focussed_Beam_Instrument_BRDF_Model : public Instrument_BRDF_Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(int,integralmode)
            DECLARE_PARAMETER(double,focal_point);
            DECLARE_PARAMETER(BRDF_Model_Ptr,model);

        public:
            Focussed_Beam_Instrument_BRDF_Model();
        protected:
            virtual MuellerMatrix mueller();

            static Vector four_angles(double theta,double phi,double alpha,double beta);
    };

    class Circle_Integral {
        public:
            Circle_Integral(int _order=1) {
                set_order(_order);
            }
            void set_order(int _order);
            double x(int i) const {
                return cos(theta(i))*r(i);
            }
            double y(int i) const {
                return sin(theta(i))*r(i);
            }
            double w(int i) const;
            int n() const;
            double r(int i) const;
            double theta(int i) const;
            int get_maxorder() const {
                return maxorder;
            }

        private:

            void check(int i) const;
            int order;
            int offset;
            int nn;
            static double table[];
            static int maxorder;
    };

} // namespace SCATMECH {


#endif



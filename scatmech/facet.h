//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: facet.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_FACET_H
#define SCATMECH_FACET_H

#include "brdf.h"
#include "filmtran.h"
#include "askuser.h"


namespace SCATMECH {


    class Slope_Distribution_Function :
        public Model
    {
        public:
            virtual double f(double slopex,double slopey) {
                return f(sqrt(sqr(slopex)+sqr(slopey)));
            }

        protected:

            virtual double f(double slope) {
                error("Single argument slope distribution not defined");
                return 0;
            }
        public:

            DECLARE_MODEL();
    };

    typedef Model_Ptr<Slope_Distribution_Function> Slope_Distribution_Function_Ptr;

    class Facet_BRDF_Model : public BRDF_Model
    {
        public:
            virtual double local_angle(double thetai,double thetas,double phis);
            virtual double local_slope(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(Slope_Distribution_Function_Ptr,sdf);
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:
            JonesMatrix jones();

            double evaluate_sdf(const Vector& _nhat);
            Vector nhat(double thetai,double thetas,double phis);

        private:
            double jacobian(double thetai,double thetas,double phis);
    };

    typedef Model_Ptr<Facet_BRDF_Model> Facet_BRDF_Model_Ptr;

    class Unit_Slope_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            double f(double slope) {
                return 1.;
            }
            DECLARE_MODEL();
    };

    class Exponential_Slope_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            double f(double slope);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,s);
    };

    class Gaussian_Slope_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            double f(double slope);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,s);
    };

    class Table_Slope_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            double f(double slope);
            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,T);
    };

    class Anisotropic_Exponential_Slope_Distribution_Function :
        public Slope_Distribution_Function
    {
        public:
            double f(double slopex,double slopey);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,sx);
            DECLARE_PARAMETER(double,sy);
    };

    class Ellipsoid_Slope_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            double f(double slopex,double slopey);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,rx);
            DECLARE_PARAMETER(double,ry);
            DECLARE_PARAMETER(double,rz);
            DECLARE_PARAMETER(double,density);
    };

    class Two_Slope_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            double f(double slopex,double slopey);

            DECLARE_MODEL();
            DECLARE_PARAMETER(Slope_Distribution_Function_Ptr,s1);
            DECLARE_PARAMETER(Slope_Distribution_Function_Ptr,s2);
            DECLARE_PARAMETER(double,fract);
    };

    class Angle_Distribution_Function:
        public Slope_Distribution_Function
    {
        public:
            // P() is the distribution of angles, which is assumed to be normalized so that
            // the integral of 2 pi theta P(theta) dtheta from 0 to pi/2 is 1...
            virtual double P(double thetax)=0;

            DECLARE_MODEL();

        private:
            // The following function is now made private...
            double f(double slope);
    };

    class Gaussian_Angle_Distribution_Function:
        public Angle_Distribution_Function
    {
        public:
            double P(double thetan);
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,sigma);

        private:
            double sigmarad;
            double norm;
    };

    class Exponential_Angle_Distribution_Function:
        public Angle_Distribution_Function
    {
        public:
            double P(double thetan);
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,sigma);

        private:
            double sigmarad;
            double norm;
    };

    class Table_Angle_Distribution_Function:
        public Angle_Distribution_Function
    {
        public:
            double P(double thetan);

            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,T);
    };


    void Register(const Facet_BRDF_Model* x);
    void Register(const Slope_Distribution_Function* x);


} // namespace SCATMECH


#endif

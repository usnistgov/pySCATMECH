//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rough.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_ROUGH_H
#define SCATMECH_ROUGH_H

#include "scatmech.h"
#include "brdf.h"
#include "askuser.h"


namespace SCATMECH {


    class PSD_Function : public Model {
        public:
            virtual double psd(double fx,double fy) {
                return psd(sqrt(fx*fx+fy*fy));
            }
        protected:
            virtual double psd(double /* f */) {
                error("Attempt to evaluate f(double)");
                return 0;
            }
        public:
            DECLARE_MODEL();
    };

    typedef Model_Ptr<PSD_Function> PSD_Function_Ptr;

    class Unit_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double /* f */) {
                return 1.;
            }
            DECLARE_MODEL();
    };

    class ABC_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double f);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,A);
            DECLARE_PARAMETER(double,B);
            DECLARE_PARAMETER(double,C);
    };

    class Table_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double f);

            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,T);
    };

    class Fractal_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double f);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,A);
            DECLARE_PARAMETER(double,exponent);
    };

    class Gaussian_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double f);
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,length);
            DECLARE_PARAMETER(double,sigma);
        private:
            double A,B;
    };

    class Elliptical_Mesa_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,axisx);
            DECLARE_PARAMETER(double,axisy);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,density);
    };

    class Rectangular_Mesa_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,lengthx);
            DECLARE_PARAMETER(double,lengthy);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,density);
    };

    class Triangular_Mesa_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,side);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,density);
    };

    class Rectangular_Pyramid_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,lengthx);
            DECLARE_PARAMETER(double,lengthy);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,density);
    };

    class Triangular_Pyramid_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,length);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,density);
    };

    class Parabolic_Dimple_PSD_Function:
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,axisx);
            DECLARE_PARAMETER(double,axisy);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,density);
    };

    class Double_PSD_Function :
        public PSD_Function
    {
        public:
            double psd(double fx,double fy);

            DECLARE_MODEL()
            DECLARE_PARAMETER(PSD_Function_Ptr,psd1);
            DECLARE_PARAMETER(PSD_Function_Ptr,psd2);
    };

    //
    // Roughness_BRDF_Model is a virtual base class for all BRDF_Model's
    // that require a power spectrum.
    //
    class Roughness_BRDF_Model : public BRDF_Model
    {
        public:
            // Function to return the spatial frequency for a given scattering
            // geometry...
            void Bragg_Frequency(double& fx,double& fy);

            DECLARE_MODEL();
            DECLARE_PARAMETER(PSD_Function_Ptr,psd);
    };

    typedef Model_Ptr<Roughness_BRDF_Model> Roughness_BRDF_Model_Ptr;

    void Register(const Roughness_BRDF_Model* x);
    void Register(const PSD_Function* x);


} // namespace SCATMECH



#endif

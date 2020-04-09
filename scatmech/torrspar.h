//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: torrspar.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_TORRSPAR_H
#define SCATMECH_TORRSPAR_H

#include "scatmech.h"
#include "facet.h"


namespace SCATMECH {


    class Shadow_Function : public Model {
        public:
            virtual double f(double thetai,double thetas,double phis) = 0;
            DECLARE_MODEL();
    };

    class Unit_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis) {
                return 1.;
            }
            DECLARE_MODEL();
    };

    class Torrance_Sparrow_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);
            DECLARE_MODEL();
    };

    //
    // A shadow function developed for Gaussian rough surfaces by
    // Bruce G. Smith, "Geometrical Shadowing of a Random Rough Surfaces,"
    // IEEE Trans. Antennas and Propagation, AP-15, 668 (1967)
    //
    class Smith_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,w);

        private:
            double S(double theta);
    };

    //
    // An empirical shadow function developed by
    // J.R. Maxwell, J. Beard, S. Weiner, D. Ladd, and S. Ladd,
    // "Bidirectional Reflectance Model Validation and Utilization"
    // Environmental Research Institute of Michigan, Ann Arbor, Oct. 1973
    // Document available from National Technical Information Service
    // Report AD913816
    //
    class Maxwell_Beard_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,tau);
            DECLARE_PARAMETER(double,Omega);
    };

    class Table_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);
            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,T);
    };

    typedef Model_Ptr<Shadow_Function> Shadow_Function_Ptr;

    //
    // Shadowed_Facet_BRDF_Model is a Facet_BRDF_Model which
    // implements the shadowing functions described by K.E.Torrance and
    // E.M.Sparrow in "Theory of off-specular reflection from roughened surfaces,"
    // J. Opt. Soc. Am. vol. 57, no. 9, pp. 1105-1114 (1967).
    //
    class Shadowed_Facet_BRDF_Model: public Facet_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(Shadow_Function_Ptr,shadow);

        protected:
            // Routine which returns the Jones matrix for scattering...
            JonesMatrix jones();
    };

    void Register(const Shadow_Function* x);


} // namespace SCATMECH


#endif

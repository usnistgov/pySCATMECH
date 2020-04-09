//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: reflectance.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_REFLECTANCE_H
#define SCATMECH_REFLECTANCE_H

#include "inherit.h"
#include "scattabl.h"
#include <string>

namespace SCATMECH {

    ///
    /// Reflectance is an abstract class for models describing the
    /// diffuse reflectance of a coating or material.
    ///
    class Reflectance : public Model {
        public:
            ///
            /// Get_Reflectance returns the reflectance as a
            /// function of wavelength.
            ///
            virtual double Get_Reflectance(double _lambda) {
                lambda = _lambda;
                if (lambda!=old_lambda) {
                    set_recalc(1);
                    old_lambda=lambda;
                }
                SETUP();
                return reflectance;
            }
        protected:

            /// The wavelength at which to evaluate the reflectance
            double lambda;

            /// The previous wavelength used.
            double old_lambda;

            /// The reflectance of the material.  This parameter must
            /// be set by an child class' setup().
            double reflectance;

            DECLARE_MODEL();
    };

    typedef Model_Ptr<Reflectance> Reflectance_Ptr;
    void Register(const Reflectance* x);

    ///
    /// Kubelka-Munk model for reflectance from a coating, described
    /// in Z. Tech. Phys. 12, 593 (1931).  It is a 1-D solution to the
    /// radiative transfer problem of a film of given thickness, absorption
    /// coefficient, scattering coefficient, and underlying substrate
    /// reflectance.  It is commonly used to describe paint coatings.
    ///
    class Kubelka_Munk_Reflectance : public Reflectance {
        protected:
            /// setup() calculates the reflectance.
            void setup();

            DECLARE_MODEL();

            /// Thickness of the coating [um]
            DECLARE_PARAMETER(double,thickness);

            /// Absorption coefficient of the coating [1/um]
            DECLARE_PARAMETER(Table,absorption);

            /// Scattering coefficient of the coating [1/um]
            DECLARE_PARAMETER(Table,scattering);

            /// Reflectance of the underlying substrate
            DECLARE_PARAMETER(Reflectance_Ptr,substrate);
    };

    ///
    /// A reflectance whose wavelength dependence is given by a
    /// Table.
    ///
    class Table_Reflectance : public Reflectance {
        protected:
            /// setup() calculates the reflectance.
            void setup();

            DECLARE_MODEL();

            /// The table containing the reflectance as a function of wavelength
            DECLARE_PARAMETER(Table,table);
    };

    ///
    /// A reflectance given by an expression.
    ///
    class Equation_Reflectance : public Reflectance {
        protected:
            /// setup() calculates the reflectance.
            void setup();

            DECLARE_MODEL();

            /// The string containing the expression for the reflectance, where
            /// the variable x is the wavelength.
            DECLARE_PARAMETER(std::string,expression);
    };
}

#endif

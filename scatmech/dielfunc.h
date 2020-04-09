//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: dielfunc.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_DIELFUNC_H
#define SCATMECH_DIELFUNC_H

#include "scatmech.h"
#include "optconst.h"
#include "scattabl.h"
#include "inherit.h"

namespace SCATMECH {

    ///
    /// The class dielectric_function stores the optical properties (n & k) as a function of wavelength
    ///
    ///
    class dielectric_function {
        public:
            /// @brief Constructor from a string
            ///
            /// Constructor from a string. If the string can be read as a complex number [e.g., "(2.3,0.2)"], then it will
            /// be interpreted as constant optical constants. If not, it attempts to read a file containing at least three
            /// columns (wavelength, n, and k).
            dielectric_function(const std::string& filename);
            dielectric_function(const char* filename);

            /// Constructor with an optical constant
            dielectric_function(const optical_constant& nn=optical_constant(COMPLEX(1,0)));
            /// Constructor with a dielectric constant
            dielectric_function(const dielectric_constant& ee);

            /// Set the dielectric function to a value, function, or file
            void set(const std::string& value);

            dielectric_constant epsilon(double lambda) const {
                update(lambda);
                return (dielectric_constant)(last.y);
            }

            optical_constant index(double lambda) const {
                update(lambda);
                return last.y;
            }

            dielectric_constant epsilon() const {
                return (dielectric_constant)(last.y);
            }
            optical_constant index() const {
                return last.y;
            }

            double e1(double lambda) const {
                update(lambda);
                return ((dielectric_constant)(last.y)).e1;
            }

            double e2(double lambda) const {
                update(lambda);
                return ((dielectric_constant)(last.y)).e2;
            }

            double n(double lambda) const {
                update(lambda);
                return last.y.n;
            }

            double k(double lambda) const {
                update(lambda);
                return last.y.k;
            }

            static dielectric_function AskUser(const std::string& query, const std::string& deflt);
            static dielectric_function AskUser(const std::string& query,const optical_constant& deflt);

            void force_nonabsorbing();

            const std::string& get_name() const {
                return name;
            }
            friend std::ostream& operator<<(std::ostream& os,const dielectric_function& df) {
                return os << df.name;
            }

            void set_parameter(const std::string& param, const std::string& value);
            double get_parameter(const std::string& param);

            void set_interpolationx(Table::InterpolationMode m) {
                n1.set_interpolationx(m);
                n2.set_interpolationx(m);
            }
            Table::InterpolationMode get_interpolationx() {
                return n1.get_interpolationx();
            }

            void set_interpolationy(Table::InterpolationMode m) {
                n1.set_interpolationy(m);
                n2.set_interpolationy(m);
            }
            Table::InterpolationMode get_interpolationy() {
                return n1.get_interpolationy();
            }

            void set_n(const Table& n);
            void set_k(const Table& k);

			const Table& get_n() {return n1;}
			const Table& get_k() {return n2;}

        private:
            Table n1;
            Table n2;

            typedef last_value<double,optical_constant> LAST;
            mutable LAST last;

            int _AskUser(const std::string& query,const std::string& deflt);

            void update(double lambda) const {
                if (last.x!=lambda) {
                    last.x=lambda;
                    last.y=optical_constant(n1.value(lambda),n2.value(lambda));
                }
            }

            std::string name;
    };

    template <>
    void ModelParameterSet(dielectric_function& variable,const std::string& subparameter,const std::string& value);
    template <>
    std::string ModelParameterGet(dielectric_function& variable,const std::string& subparameter);
    template <>
    void ModelParameterAskUser(dielectric_function& variable,const std::string& prompt);

    extern dielectric_function vacuum;

} // namespace SCATMECH


#endif

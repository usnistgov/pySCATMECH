//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: scateval.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_SCATEVAL_H
#define SCATMECH_SCATEVAL_H

#include "scatmech.h"
#include <string>
#include <stack>
#include <vector>
#include <sstream>
#include <map>

namespace SCATMECH {

    /// @brief A class for evaluating expressions
    ///
    /// A class that takes a string expression (e.g., "3+4*x") and evaluates
    /// it using an optional list of variables.  Results can be multiple
    /// values separated by commas and surrounded by parentheses.
    class Evaluator {
        public:

            /// A vector of doubles
            typedef std::vector<double> VDOUBLE;

            /// A string
            typedef std::string STRING;

            /// A map from a string to a double
            typedef std::map<STRING,double> VMAP;

            /// Constructor
            Evaluator(const STRING& str,                 ///< Expression to evaluate
                      const VMAP& _variables = VMAP(),   ///< Map of variables
                      bool _top=true                     ///< Normally true, but false for recursive calls
                     ) : input(str), variables(_variables), top(_top) {
				if (top) {
					// If top level, provide a meaningful error message...
					try {
						evaluate();
					}
					catch (std::exception& e) {
						throw SCATMECH_exception("Error during evaluation of expression \"" + str + "\": " + e.what());
					}
				} else {
					evaluate();
				}
            }

            /// Get the result, forcing it to have only one value.
            operator double() const {
                if (result.size()!=1) error("Call to double() with more than one value");
                return result[0];
            }

            /// Get the result as a character string
            STRING ResultString() const;

            /// Get the i-th result
            double Result(int i) {
                return result[i];
            }

            /// Get the number of results
            int NResult() {
                return result.size();
            }

        private:

            const VMAP& variables;

            std::istringstream input;
            std::stack<double> val_stack;
            std::stack<int> op_stack;
            std::stack<int> prec_stack;
            VDOUBLE result;

            void evaluate();
            void operate();
            double call_function(const STRING& s,const VDOUBLE& args);
            void error(const std::string& message) const;
            bool lower_prec(int prec);
            void get_value();
            void get_operator();
            VDOUBLE get_paren(bool allow_empty=false);

            bool top;
    };

    void Register_Evaluator_Function(double (*f)(const std::vector<double>& args),const std::string& name,int nargs);
    bool Evaluate_Function(double *value,const std::string& name, const std::vector<double>& args);
    std::vector<double> OneLineFunction(const std::string& expression,double x,Evaluator::VMAP variables=Evaluator::VMAP());
}

#endif

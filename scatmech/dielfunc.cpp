//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: dielfunc.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include <ctype.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>
#include "scatmech.h"
#include "dielfunc.h"
#include "askuser.h"

using namespace std;

namespace SCATMECH {

    dielectric_function::
    dielectric_function(const optical_constant& nn)
    {
        last.x = -1;
        last.y = nn;
        n1=nn.n;
        n2=nn.k;
        ostringstream temp;
        temp<<(COMPLEX)nn;
        name=temp.str();
    }

    dielectric_function::
    dielectric_function(const dielectric_constant& ee)
    {
        optical_constant nn = ee;
        last.x = -1;
        last.y = nn;
        n1=nn.n;
        n2=nn.k;
        ostringstream temp;
        temp<<(COMPLEX)nn;
        name=temp.str();
    }

    dielectric_function::
    dielectric_function(const std::string& filename)
    {
        set(filename);
    }

    dielectric_function::
    dielectric_function(const char* filename)
    {
        set(filename);
    }

    int
    dielectric_function::
    _AskUser(const string& query,const string& deflt)
    {
        last.x=-1;

        string response = SCATMECH::AskUser(query + " [n, (n,k), or filename]",deflt);

        istringstream rr(response);

        // First, try interpreting the response as a double...
        double result1;
        rr >> result1;

        if (rr.fail()) {
            rr.clear();
            COMPLEX result3;
            rr >> result3;
            if (rr.fail()) {
                // Finally, try to interpret the response as a filename...
                string filename;
                rr.clear();
                rr >> filename;
                set(filename);
            } else {
                *this = dielectric_function(optical_constant(result3));
            }
        } else {
            double result2 = SCATMECH::AskUser(query + " (k)",0.0);

            *this = dielectric_function(optical_constant(result1,result2));
        }
        return 1;
    }

    dielectric_function
    dielectric_function::
    AskUser(const std::string& query, const std::string& deflt)
    {
        dielectric_function result;
        bool success;
        do {
            try {
                result._AskUser(query,deflt);
                success = true;
            } catch (SCATMECH_exception& e) {
                SCATMECH_output << e.what() << endl;
                success = false;
            }
        } while (!success);
        return result;
    }

    dielectric_function
    dielectric_function::
    AskUser(const std::string& query,const optical_constant& deflt)
    {
        ostringstream gg;
        gg << (COMPLEX)deflt;
        return AskUser(query,gg.str());
    }

    void
    dielectric_function::
    force_nonabsorbing()
    {
        n2=0.;
        last.x=-1;
    }

    void
    dielectric_function::
    set(const std::string& value)
    {
        last.x=-1;

        istringstream rr(value);

        // First try interpreting it as a double...
        double result1;
        rr >> result1;

        if (rr.fail()) {
            // Next, try interpreting it as a complex...
            rr.clear();
            COMPLEX result2;
            rr >> result2;
            if (rr.fail()) {
                // Finally, try it as a filename...
                string filename;
                rr.clear();
                rr >> filename;
                n1.set(filename,2);
                n2.set(filename,3);
                name=filename;
                //read(filename);
            } else {
                *this = dielectric_function(optical_constant(result2));
            }
        } else {
            *this = dielectric_function(optical_constant(result1));
        }
    }

    void
    dielectric_function::
    set_n(const Table& n)
    {
        n1 = n;
        last.x=-1.;
		index(1.);
    }

    void
    dielectric_function::
    set_k(const Table& k)
    {
        n2 = k;
        last.x=-1.;
		index(1.);
    }

    void
    dielectric_function::
    set_parameter(const string& param,const string& value)
    {
        n1.set_parameter(param,value);
        n2.set_parameter(param,value);
        update(last.x);
    }

    double
    dielectric_function::
    get_parameter(const string& param)
    {
        return from_string<double>(n1.get_parameter(param));
    }

    template <>
    void
    ModelParameterSet(dielectric_function& variable,const string& parameter,const string& value)
    {
        if (parameter.empty()) {
            variable = value;
        } else if (parameter == "n") {
            Table n;
            n.set(value);
            variable.set_n(n);
        } else if (parameter == "k") {
            Table k;
            k.set(value);
            variable.set_k(k);
        } else {
            variable.set_parameter(parameter,value);
        }
    }

    template <>
    string
    ModelParameterGet(dielectric_function& variable,const string& parameter)
    {
        if (parameter.empty()) {
            ostringstream oss;
            oss << variable;
            return oss.str();
        } else if (parameter == "n") {
            ostringstream oss;
            oss << variable.index().n;
            return oss.str();
        } else if (parameter == "k") {
            ostringstream oss;
            oss << variable.index().k;
            return oss.str();
        } else {
            return to_string(variable.get_parameter(parameter));
        }
    }


    template <>
    void
    ModelParameterAskUser(dielectric_function& value,const string& prompt)
    {
        value = dielectric_function::AskUser(prompt,value.index());
    }

} // namespace SCATMECH


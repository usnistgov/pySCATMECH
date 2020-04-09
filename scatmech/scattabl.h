//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: scattabl.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_SCATTABL_H
#define SCATMECH_SCATTABL_H

#include <vector>
#include <utility>
#include <string>
#include "inherit.h"
#include "scateval.h"

namespace SCATMECH {

    template <class X,class Y>
    struct last_value {
        X x;
        Y y;
    };

    class FormulaFile;
    class TableFile;

    class Table
    {
        public:
            typedef last_value<double,double> LAST;
            typedef std::pair<double,double> PAIR;
            typedef std::vector<PAIR> VECTOR;

            enum InterpolationMode {
                LIN, LOG, SPLINE, LOGSPLINE
            };

            Table(const std::string& filename,int col=2);
            Table(const char* filename,int col=2);
            Table(double singlevalue=0.);
            Table(double *l,double *v,int nn);
            Table(const std::vector<double>& x,const std::vector<double>& y);
            Table(const VECTOR& v);

            double value(double l) const ;
            double value() const {
                return last.y;
            };

            int set(const std::string& filename, int ycol=2);

            static Table AskUser(const std::string& prompt,const std::string& deflt);
            static Table AskUser(const std::string& prompt,double deflt);

            void set_parameter(const std::string& param,const std::string& value);
            std::string get_parameter(const std::string& param) const;

            const std::string& get_name() const {
                return name;
            }
            friend std::ostream& operator<<(std::ostream& os,const Table& t) {
                return os << t.name;
            }

            void set_interpolationx(InterpolationMode m) {
                interpolationx=m;
            }
            InterpolationMode get_interpolationx() const {
                return interpolationx;
            }

            void set_interpolationy(InterpolationMode m) {
                interpolationy=m;
            }
            InterpolationMode get_interpolationy() const {
                return interpolationy;
            }

            const VECTOR& get_table() const {
                return values;
            }

        private:

            int _AskUser(const std::string& prompt,const std::string& deflt);

            bool ReadFunctionFormat(const std::string& fname);

            mutable LAST last;
            VECTOR values;
            std::string name;
            InterpolationMode interpolationx,interpolationy;

            int icol;
            Evaluator::VMAP params;
            FormulaFile *formulafile;
            TableFile *tablefile;
            std::string formulastring;
    };

    void ClearTableCache();

    template <>
    void ModelParameterSet(Table& variable,const std::string& subparameter,const std::string& value);

    template <>
    std::string ModelParameterGet(Table& variable,const std::string& subparameter);

    template <>
    void ModelParameterAskUser(Table& variable,const std::string& prompt);

} // namespace SCATMECH

#endif

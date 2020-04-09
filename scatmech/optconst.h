//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: optconst.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_OPTCONST_H
#define SCATMECH_OPTCONST_H

#include <iostream>
#include <complex>
#include "scatmech.h"

namespace SCATMECH {

    class dielectric_constant;

    class optical_constant {
        public:
            double n;
            double k;
            optical_constant(): n(1.), k(0.) {}
            optical_constant(const COMPLEX& x): n(std::real(x)),k(std::imag(x)) {}
            optical_constant(double x, double y=0.): n(x),k(y) {}
            optical_constant(const dielectric_constant& e);
            optical_constant(const optical_constant& oc): n(oc.n), k(oc.k) {}
            bool operator==(const optical_constant& a) const {
                return (n==a.n&&k==a.k);
            }
            bool operator!=(const optical_constant& a) const {
                return (n!=a.n||k!=a.k);
            }

            operator COMPLEX() const {
                return COMPLEX(n,k);
            }
            friend std::istream& operator>>(std::istream& is,optical_constant& nn)
            {
                COMPLEX temp;
                is >> temp;
                nn=temp;
                return is;
            }
            friend std::ostream& operator<<(std::ostream& os,const optical_constant& nn)
            {
                return os << '(' << nn.n << ',' << nn.k << ')';
            }
    };

    class dielectric_constant {
        public:
            double e1;
            double e2;
            dielectric_constant(): e1(1.),e2(0.) {}
            dielectric_constant( const COMPLEX& x): e1(std::real(x)),e2(std::imag(x)) {}
            dielectric_constant( double x, double y=0.): e1(x),e2(y) {}
            dielectric_constant( const optical_constant& n);
            dielectric_constant( const dielectric_constant& dc) : e1(dc.e1),e2(dc.e2) {}
            bool operator==(const dielectric_constant& a) const {
                return (e1==a.e1&&e2==a.e2);
            }
            bool operator!=(const dielectric_constant& a) const {
                return (e1!=a.e1||e2!=a.e2);
            }
            operator COMPLEX() const {
                return COMPLEX(e1,e2);
            }
            friend std::istream& operator>>(std::istream& is,dielectric_constant& e)
            {
                COMPLEX temp;
                is >> temp;
                e=temp;
                return is;
            }
            friend std::ostream& operator<<(std::ostream& os,const dielectric_constant& ee)
            {
                return os << '(' << ee.e1 << ',' << ee.e2 << ')';
            }

    };

    inline
    dielectric_constant::
    dielectric_constant(const optical_constant& n)
    {
        *this = dielectric_constant((COMPLEX)n*(COMPLEX)n);
    }

    inline
    optical_constant::
    optical_constant(const dielectric_constant& e)
    {
        *this =  optical_constant(std::sqrt((COMPLEX)e));
    }

    optical_constant AskUser(const std::string& query,const optical_constant& deflt);

} // namespace SCATMECH


#endif


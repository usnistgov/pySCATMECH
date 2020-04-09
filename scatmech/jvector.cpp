//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: jvector.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "mueller.h"

using namespace std;


namespace SCATMECH {


    JonesVector::
    JonesVector(const StokesVector& x)
    {
        StokesVector xx = x.pol_part();
        // Next line corrected for Version 4 (TAG 5 MAR 2003)...
        double arg = (xx[3]==0.&&xx[2]==0.) ? 0. : atan2(xx[3],xx[2]);
        double mag = sqrt((xx[0]-xx[1])/2.);
        (*this)[0] = sqrt((xx[0]+xx[1])/2.);
        (*this)[1] = COMPLEX(mag*cos(arg),mag*sin(arg));
    }

    double
    JonesVector::
    DOCP() const
    {
        return ((StokesVector)(*this)).DOCP();
    }

    double
    JonesVector::
    e() const
    {
        return ((StokesVector)(*this)).e();
    }

    double
    JonesVector::
    DOLP() const
    {
        return ((StokesVector)(*this)).DOLP();
    }

    double
    JonesVector::
    eta() const
    {
        return ((StokesVector)(*this)).eta();
    }


    ostream& operator<<(ostream& os,const JonesVector& j)
    {
        return os << '(' << j[0] << ',' << j[1] << ')';
    }

    istream& operator>>(istream& is,JonesVector& j)
    {
        COMPLEX js,jp;
        char c;
        is >> c;
        if (c=='(') {
            is >> js;
            is >> c;
            if (c==',') {
                is >> jp;
                is >> c;
                if (c==')') {
                    j = JonesVector(js,jp);
                    return is;
                }
            }
        }
        is.setstate(ios::failbit);
        return is;
    }

    JonesVector
    JonesVector::
    rotate(const double angle) const
    {
        return JonesRotator(angle)*(*this);
    }


} // namespace SCATMECH


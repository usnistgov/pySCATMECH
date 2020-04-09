//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: stokes.cpp
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


    double
    StokesVector::
    psi() const
    {
        return ((JonesVector)((*this).pol_part())).psi();
    }

    double
    StokesVector::
    delta() const
    {
        return ((JonesVector)((*this).pol_part())).delta();
    }

    StokesVector
    StokesVector::
    operator*(const MuellerMatrix& matrix) const
    {
        StokesVector out;
        int i,j;
        for (i=0; i<4; ++i) {
            out.s[i]=0;
            for (j=0; j<4; ++j)
                out.s[i]+=s[j]*matrix.m[j][i];
        }
        return out;
    }

    StokesVector
    StokesVector::
    rotate(double angle) const
    {
        StokesVector temp(*this);
        return MuellerMatrix(JonesRotator(angle))*temp;
    }

    StokesVector
    StokesVector::
    pol_part() const
    {
        StokesVector ss;
        double dop = DOP();
        ss.s[0] = s[0]*dop;
        ss.s[1] = s[1];
        ss.s[2] = s[2];
        ss.s[3] = s[3];
        return ss;
    }

    StokesVector
    StokesVector::
    unpol_part() const
    {
        StokesVector ss;
        double dop = DOP();
        ss.s[0] = s[0]*(1.-dop);
        ss.s[1] = 0;
        ss.s[2] = 0;
        ss.s[3] = 0;
        return ss;
    }

    ostream& operator<<(ostream& os,const StokesVector& s)
    {
        return os << '(' << s[0] << ',' << s[1] << ',' << s[2] << ',' << s[3] << ')';
    }

    istream& operator>>(istream& is, StokesVector& s)
    {
        double s1,s2,s3,s4;
        char c;
        is >> c;
        if (c=='(') {
            is >> s1;
            is >> c;
            if (c==',') {
                is >> s2;
                is >> c;
                if (c==',') {
                    is >> s3;
                    is >> c;
                    if (c==',') {
                        if (c==',') {
                            is >> s4;
                            is >> c;
                            if (c==')') {
                                s = StokesVector(s1,s2,s3,s4);
                                return is;
                            }
                        }
                    }
                }
            }
        }
        is.setstate(ios::failbit);
        return is;
    }

} // namespace SCATMECH


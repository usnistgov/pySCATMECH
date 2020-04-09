//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: jmatrix.cpp
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
#include "vector3d.h"

using namespace std;


namespace SCATMECH {



    //******************************************************************************
    //**
    //** This file contains the executables for the class JonesMatrix.
    //**
    //******************************************************************************


    JonesMatrix::
    JonesMatrix(const MuellerMatrix& mm)
    {
        //
        // Transformation from MuellerMatrix to JonesMatrix from
        // R.A. Chipman, "Polarimetry" in _Handbook_of_Optics_Vol._II_
        // (McGraw-Hill, New York, 1995).
        //
        // Modified by TAG to use atan2(,)
        //
        j[0] = sqrt(COMPLEX((mm.m[0][0]-mm.m[0][1]-mm.m[1][0]+mm.m[1][1])/2.))*
               exp(COMPLEX(0,1)*atan2((mm.m[3][2]-mm.m[2][3]),(mm.m[2][2]+mm.m[3][3])));
        j[1] = sqrt(COMPLEX((mm.m[0][0]+mm.m[0][1]+mm.m[1][0]+mm.m[1][1])/2.));
        j[2] = sqrt(COMPLEX((mm.m[0][0]-mm.m[0][1]+mm.m[1][0]-mm.m[1][1])/2.))*
               exp(COMPLEX(0,1)*atan2((-mm.m[0][3]-mm.m[1][3]),(mm.m[0][2]+mm.m[1][2])));
        j[3] = sqrt(COMPLEX((mm.m[0][0]+mm.m[0][1]-mm.m[1][0]-mm.m[1][1])/2.))*
               exp(COMPLEX(0,1)*atan2((mm.m[3][0]+mm.m[3][1]),(mm.m[2][0]+mm.m[2][1])));
    }

    //******************************************************************************
    //**
    //** The following are the stream i/o functions ...
    //**
    //******************************************************************************

    ostream& operator<<(ostream& os, const JonesMatrix& j)
    {
        return os << '(' << j[0] << ',' << j[1] << ',' << j[2] << ',' << j[3] << ')';
    }

    istream& operator>>(istream& is, JonesMatrix& j)
    {
        COMPLEX j1,j2,j3,j4;
        char c;
        is >> c;
        if (c=='(') {
            is >> j1;
            is >> c;
            if (c==',') {
                is >> j2;
                is >> c;
                if (c==',') {
                    is >> j3;
                    is >> c;
                    if (c==',') {
                        is >> j4;
                        is >> c;
                        if (c==')') {
                            j = JonesMatrix(j1,j2,j3,j4);
                            return is;
                        }
                    }
                }
            }
        }
        is.setstate(ios::failbit);
        return is;
    }

    //******************************************************************************
    //**
    //** The following are various optimum polarization modifying elements...
    //**
    //******************************************************************************

    JonesMatrix
    JonesRotator(double angle)
    {
        double c = cos(angle);
        double s = sin(angle);
        return JonesMatrix(c,c,s,-s);
    }

    JonesMatrix
    JonesLinearRetarder(double phase,double angle)
    {
        JonesMatrix out;
        COMPLEX phaser=exp(-COMPLEX(0,1)*phase/2.);
        out[0]=phaser;       // advanced by phase/2
        out[1]=1./phaser;     // delayed by phase/2
        out[2]=out[3]=0;
        return out.rotate(angle);
    }

    JonesMatrix
    JonesCircularRetarder(double phase)
    {
        JonesMatrix out;
        out[0]=cos(phase/2);
        out[1]=out[0];
        out[2]=sin(phase/2);
        out[3]=-out[2];
        return out;
    }

    JonesMatrix
    JonesLinearPolarizer(double angle,double diattenuation)
    {
        JonesMatrix out(1,sqrt(-(diattenuation-1.)/(diattenuation+1.)),0,0);
        return out.rotate(angle);
    }

    JonesMatrix
    JonesCircularPolarizer(double diattenuation)
    {
        double e = sqrt((diattenuation-1.)/(diattenuation+1.));
        double x = (1.+e)/2.;
        COMPLEX y = COMPLEX(0.,(1.-e)/2.);
        JonesMatrix out(x,x,y,-y);
        return out;
    }

    //******************************************************************************
    //**
    //** The following function returns a JonesMatrix rotated...
    //**
    //******************************************************************************

    JonesMatrix
    JonesMatrix::
    rotate(const double angle) const
    {
        JonesMatrix rot = JonesRotator(angle);
        return rot*((*this)*rot.transpose());
    }

    //******************************************************************************
    //**
    //** The following returns the JonesMatrix given two Jones eigenvectors and
    //** their respective eigenvalues (transmittances).
    //**
    //******************************************************************************

    JonesMatrix
    JonesGeneralized(const JonesVector& a, const JonesVector& b,
                     const COMPLEX& ma, const COMPLEX& mb)
    {
        COMPLEX det = a[0]*b[1]-a[1]*b[0];

        return JonesMatrix(a[0]*b[1]*mb-a[1]*b[0]*ma,
                           a[0]*b[1]*ma-a[1]*b[0]*mb,
                           a[0]*b[0]*(mb-ma),
                           a[1]*b[1]*(ma-mb))/det;
    }

    void GetBasisVectorsSP(const Vector& k,Vector& s,Vector& p)
    {
        static Vector zhat(0,0,1);
        if (fabs(k.x)>1E-10||fabs(k.y)>1E-10) {
            s = perpto(zhat,k);
            p = perpto(k,s);
        } else {
            s = Vector(0,1,0);
            p = perpto(k,s);
        }
    }

    void GetBasisVectorsXY(const Vector& k,Vector& x,Vector& y)
    {
        Vector s,p;
        GetBasisVectorsSP(k,s,p);
        if (fabs(k.x)>1E-10||fabs(k.y)>1E-10) {
            double norm = sqrt(sqr(k.x)+sqr(k.y));
            double kx = k.x/norm;
            double ky = k.y/norm;
			if (k.z>0) 
                x = -p*kx-s*ky;
            else
                x = p*kx-s*ky;
        } else {
            x = Vector(1,0,0);
        }
        y = perpto(x,k); // {y,x,k} is right-handed.
    }

    void GetBasisVectorsParPerp(const Vector& kin,const Vector& kout, Vector& perp, Vector& parin, Vector& parout)
    {
        perp = perpto(kin,kout);
        parin = perpto(kin,perp);
        parout = perpto(kout,perp);
    }

    JonesMatrix GetJonesRotator(const Vector& xo, const Vector& yo, const Vector& xi, const Vector& yi)
    {
        JonesMatrix result;
        result.PP() = yo*yi;
        result.SS() = xo*xi;
        result.SP() = yo*xi;
        result.PS() = xo*yi;
        return result;
    }



} // namespace SCATMECH


//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: vector.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "vector3d.h"

using namespace std;


namespace SCATMECH {



    Vector3D<double>
    perpto(const Vector3D<double>& a,const Vector3D<double>& b)
    {
        // Modified by TAG 14 MAR 2005 by adding these two lines, and comparing
        // _norm to a small number, rather than zero...
        Vector3D<double> _a = unit(a);
        Vector3D<double> _b = unit(b);
		
        Vector3D<double> temp=cross(_a,_b);
        double _norm=Norm(temp);
        if (_norm>1E-10) return temp/_norm;
        else {
            static Vector3D<double> x(1.,0.,0.);
            Vector3D<double> temp=cross(_a,x);
            double _norm = Norm(temp);
            if (_norm>1E-10) return temp/_norm;
            else {
                static Vector3D<double> y(0.,1.,0.);
                Vector3D<double> temp=cross(_a,y);
                double _norm = Norm(temp);

                return temp/_norm;
            }
        }
    }

    Vector3D<COMPLEX >
    perpto(const Vector3D<COMPLEX >& a,const Vector3D<COMPLEX >& b)
    {
        Vector3D<COMPLEX > temp=cross(a,b);
        //double _norm=Norm(temp);
        COMPLEX _norm = sqrt(sqr(temp.x)+sqr(temp.y)+sqr(temp.z));
        if (abs(_norm)>1E-10) return temp/_norm;
        else {
            Vector3D<COMPLEX > x((COMPLEX)(1.),(COMPLEX)(0.),(COMPLEX)(0.));
            Vector3D<COMPLEX > temp=cross(a,x);
            //double _norm = Norm(temp);
            COMPLEX _norm = sqrt(sqr(temp.x)+sqr(temp.y)+sqr(temp.z));
            if (abs(_norm)>1E-10) return temp/_norm;
            else {
                Vector3D<COMPLEX > y((COMPLEX)(0.),(COMPLEX)(1.),(COMPLEX)(0.));
                Vector3D<COMPLEX > temp=cross(a,y);
                //double _norm = Norm(temp);
                COMPLEX _norm = sqrt(sqr(temp.x)+sqr(temp.y)+sqr(temp.z));
                return temp/_norm;
            }
        }
    }


} // namespace SCATMECH




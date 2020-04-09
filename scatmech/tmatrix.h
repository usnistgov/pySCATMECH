//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: tmatrix.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef TMATRIX_H
#define TMATRIX_H

#include "inherit.h"
#include "sphrscat.h"
#include "scatmatrix.h"

namespace SCATMECH {

    class Axisymmetric_Shape;
    typedef Model_Ptr<Axisymmetric_Shape> Axisymmetric_Shape_Ptr;
    void Register(const Axisymmetric_Shape* x);

    class shape_rec {
        public:
            double shape; // The radial coordinate
            double dshape; // The derivative of the radial coordinate
            double abs;      // The angle at which the radial coordinate is evaluated
            double wt;      // The weighting factor for integration
    };

    typedef std::vector<shape_rec> shape_vec;

    class Axisymmetric_Shape : public Model {
        public:
            // Get_Volume() returns the volume of the particle...
            double Get_Volume();

            // Get_MaxRadius() returns the largest distance of the
            // particle surface from the origin of the particle...
            double Get_MaxRadius();

            // Get_Base_Length() returns the distance from the origin of
            // the particle to the lowest vertical point of the particle...
            double Get_Base_Length();

            // Get_Breast() returns the largest horizontal distance of the
            // particle surface from the axis...
            double Get_Breast();

            // Get_TMatrix() finds the T-Matrix ...
            void Get_TMatrix(ScatterTMatrix& T, int lmax, int mmax, double k, COMPLEX& index);

            // Write() creates a file containing the particle shape...
            void Write(const std::string& filename);

        protected:
            void Insert_Points(double start,double end);
            void Flip();
            void Renormalize(double newV);

            // The following overrides the same named function of Model...
            virtual void set_parameter_base(const STRING& parameter, const STRING& value);

            DECLARE_MODEL();
            DECLARE_PARAMETER(int,npoints);

        private:
            shape_vec shape;

        protected:
            virtual void setup();
            virtual double f(double theta) const = 0;
            virtual double df(double theta) const;
    };

    class Ellipsoid_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;
            double df(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,vertical);
            DECLARE_PARAMETER(double,horizontal);
            DECLARE_PARAMETER(double,offset);
    };

    class Cylinder_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;
            double df(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(double,length);
            DECLARE_PARAMETER(double,corner);
            DECLARE_PARAMETER(int,renorm);
        private:
            double theta0,theta1,theta2,c;
    };

    class Table_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,table);
            DECLARE_PARAMETER(double,scale);
    };

    class Conical_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;
            double df(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,offset);
            DECLARE_PARAMETER(double,corner_base);
            DECLARE_PARAMETER(double,corner_apex);
            DECLARE_PARAMETER(int,updown);
            DECLARE_PARAMETER(int,renorm);
        private:
            double z0,beta,a0,theta0,theta1,theta2,theta3,c0,g;
    };

    class Indented_Ellipsoid_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,vertical);  // vertical minor axis
            DECLARE_PARAMETER(double,horizontal);  // horizontal minor axis
            DECLARE_PARAMETER(double,indent);  // indentation length
            DECLARE_PARAMETER(int,updown);
            DECLARE_PARAMETER(int,renorm);
        private:
            double theta0,B;
    };

    class Double_Conical_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(double,corner_apex);
            DECLARE_PARAMETER(double,corner_waist);
            DECLARE_PARAMETER(int,renorm);
        private:
            double theta0,theta1,alpha,beta,c1,c2;
    };

    class Chebyshev_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(int,n);        // Number of half-oscillations
            DECLARE_PARAMETER(double,start);// Start angle of oscillations
            DECLARE_PARAMETER(double,end);  // End angle of oscillations
            DECLARE_PARAMETER(double,amplitude);   // Relative amplitude of oscilations
            DECLARE_PARAMETER(double,vertical);   // Vertical mean radius
            DECLARE_PARAMETER(double,horizontal);   // Horizontal mean radius
            DECLARE_PARAMETER(int,renorm);
        private:
            double a;
    };

    class Gauss_Legendre_Integration
    {
        public:
            Gauss_Legendre_Integration(int order,double x1,double x2);
            double weight(int i) const {
                return w[i];
            }
            double abscissa(int i) const {
                return x[i];
            }
            int n() const {
                return x.size();
            }
        private:
            std::vector<double> w;
            std::vector<double> x;
    };


}
#endif

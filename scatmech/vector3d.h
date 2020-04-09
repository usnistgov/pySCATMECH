//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: vector3d.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_VECTOR3D_H
#define SCATMECH_VECTOR3D_H

#include "scatmech.h"
#include <ios>

namespace SCATMECH {

    ///
    /// template class Vector3D is used to handle 3-dimensional vectors
    ///
    template <class TYPE>
    class Vector3D {
        public:
            ///
            /// The components of the vector
            ///
            TYPE x,y,z;

            /// Default class constructor
            Vector3D<TYPE>() {}

            /// Class constructor from three values
            Vector3D<TYPE>(const TYPE& _x,const TYPE& _y,const TYPE& _z)
            {
                x=_x;
                y=_y;
                z=_z;
            }

            /// Type-converting copy constructor
            template <class TYPE2>
            Vector3D<TYPE>(const Vector3D<TYPE2>& a)
            {
                x=a.x;
                y=a.y;
                z=a.z;
            }

            /// Assignment operator
            Vector3D<TYPE>& operator=(const Vector3D<TYPE>& a)
            {
                x=a.x;
                y=a.y;
                z=a.z;
                return *this;
            }

            /// Addition of vectors
            Vector3D<TYPE> operator+(const Vector3D<TYPE>& a) const
            {
                return Vector3D<TYPE>(x+a.x,y+a.y,z+a.z);
            }

            /// Addition of vectors
            const Vector3D<TYPE>& operator+=(const Vector3D<TYPE>& a)
            {
                return *this = *this+a;
            }

            /// Subtraction of vectors
            Vector3D<TYPE> operator-(const Vector3D<TYPE>& a) const
            {
                return Vector3D<TYPE>(x-a.x,y-a.y,z-a.z);
            }

            /// Subtraction of vectors
            const Vector3D<TYPE>& operator-=(const Vector3D<TYPE>& a)
            {
                return *this = *this-a;
            }

            /// Unary minus sign
            Vector3D<TYPE> operator-() const
            {
                return Vector3D<TYPE>(-x,-y,-z);
            }

            /// Scalar product of two vectors
            TYPE operator*(const Vector3D<TYPE>& a) const
            {
                return x*a.x+y*a.y+z*a.z;
            }

            /// Product of vector with scalar
            Vector3D<TYPE> operator*(const TYPE& b) const
            {
                return Vector3D<TYPE>(x*b,y*b,z*b);
            }

            /// Product of a scalar with a vector
            friend Vector3D<TYPE> operator*(const TYPE& b,
                                            const Vector3D<TYPE>& a)
            {
                return Vector3D<TYPE>(a.x*b,a.y*b,a.z*b);
            }

            /// Division of a vector by a scalar
            Vector3D<TYPE> operator/(const TYPE& a) const
            {
                return Vector3D<TYPE>(x/a,y/a,z/a);
            }

            /// Division of vector by a scalar
            const Vector3D<TYPE>& operator/=(const TYPE& a)
            {
                return *this = *this/a;
            }

            /// The cross product of two vectors
            friend Vector3D<TYPE> cross(const Vector3D<TYPE>& a,
                                        const Vector3D<TYPE>& b)
            {   return Vector3D<TYPE>(a.y*b.z-a.z*b.y,
                                      a.z*b.x-a.x*b.z,
                                      a.x*b.y-a.y*b.x);
            }

            /// Multiplication of vector by scalar
            const Vector3D<TYPE>& operator*=(const TYPE& a)
            {
                return *this = *this*a;
            }

            /// Comparison
            bool operator==(const Vector3D<TYPE>& b) const
            {
                return (x==b.x && y==b.y && z==b.z);
            }

            /// Comparison
            bool operator!=(const Vector3D<TYPE>& b) const
            {
                return (x!=b.x || y!=b.y || z!=b.z);
            }

            /// Send to output stream
            friend std::ostream& operator<<(std::ostream& os,const Vector3D<TYPE>& v)
            {
                return os << '(' << v.x << ',' << v.y << ',' << v.z << ')';
            }

            /// Extract from an input stream
            friend std::istream& operator>>(std::istream& is,Vector3D<TYPE>& v)
            {   TYPE vx,vy,vz;
                char c;
                is >> c;
                if (c=='(') {
                    is >> vx;
                    is >> c;
                    if (c==',') {
                        is >> vy;
                        is >> c;
                        if (c==',') {
                            is >> vz;
                            is >> c;
                            if (c==')') {
                                v.x=vx, v.y=vy, v.z=vz;
                                return is;
                            }
                        }
                    }
                }
                // If it gets here, there was a problem
                is.setstate(std::ios::failbit);
                return is;
            }
            const TYPE& operator[](int i) const {
                return *(&x+i);
            }
            TYPE& operator[](int i) {
                return *(&x+i);
            }

    };

    /// The norm of a real vector:
    inline
    double
    Norm(const Vector3D<double>& a)
    {
        return sqrt(sqr(a.x)+sqr(a.y)+sqr(a.z));
    }

    /// A vector of unit length:
    inline
    Vector3D<double>
    unit(const Vector3D<double>& a)
    {
        return a/Norm(a);

    }

    /// A function returning a vector having specific polar coordinates:
    inline
    Vector3D<double>
    polar(double r,double theta,double phi)
    {
        return Vector3D<double>(r*sin(theta)*cos(phi),
                                r*sin(theta)*sin(phi),
                                r*cos(theta));
    }

    /// Complex conjugate of a complex vector
    inline
    Vector3D<COMPLEX>
    Conj(const Vector3D<COMPLEX>& a)
    {
        using std::conj;
        return Vector3D<COMPLEX>(conj(a.x),conj(a.y),conj(a.z));
    }

    /// The norm of a complex vector:
    inline
    double
    Norm(const Vector3D<COMPLEX >& a)
    {
        using std::norm;
        return sqrt(norm(a.x)+norm(a.y)+norm(a.z));
    }

    /// A vector of unit length:
    inline
    Vector3D<COMPLEX >
    unit(const Vector3D<COMPLEX >& a)
    {
        return a/Norm(a);
    }

    /// A vector perpendicular to two vectors:
    Vector3D<double> perpto(const Vector3D<double>& a,
                            const Vector3D<double>& b);
    Vector3D<COMPLEX > perpto(const Vector3D<COMPLEX >& a,
                              const Vector3D<COMPLEX >& b);

    /// A real vector
    typedef Vector3D<double> Vector;
    /// A complex vector
    typedef Vector3D<COMPLEX > CVector;


} // namespace SCATMECH

#endif // SCATMECH_VECTOR3D_H

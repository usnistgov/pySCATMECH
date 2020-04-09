//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: matrix3d.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_MATRIX3D_H
#define SCATMECH_MATRIX3D_H

#include "vector3d.h"


namespace SCATMECH {


    ///
    /// Template class Matrix3D handles 3D matrices
    ///
    template <class TYPE>
    class Matrix3D {
        private:
            ///The elements of the matrix
            TYPE x[3][3];

        public:
            /// Access to the elements
            const TYPE* operator[](int i) const {
                return x[i];
            }
            TYPE* operator[](int i) {
                return x[i];
            }

            /// Default constructor
            Matrix3D() {};

            /// Copy constructor
            template <class TYPE2>
            Matrix3D(const Matrix3D<TYPE2>& a)
            {
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        x[i][j]=a[i][j];
            }

            /// Element by element constructor
            Matrix3D(const TYPE& a,const TYPE& b,const TYPE& c,
                     const TYPE& d,const TYPE& e,const TYPE& f,
                     const TYPE& g,const TYPE& h,const TYPE& i)
            {
                x[0][0]=a;
                x[0][1]=b;
                x[0][2]=c;
                x[1][0]=d;
                x[1][1]=e;
                x[1][2]=f;
                x[2][0]=g;
                x[2][1]=h;
                x[2][2]=i;
            }

            /// Constructor that specifies three vectors
            Matrix3D(const Vector3D<TYPE>& a,const Vector3D<TYPE>& b,const Vector3D<TYPE>& c)
            {
                x[0][0]=a.x;
                x[0][1]=a.y;
                x[0][2]=a.z;
                x[1][0]=b.x;
                x[1][1]=b.y;
                x[1][2]=b.z;
                x[2][0]=c.x;
                x[2][1]=c.y;
                x[2][2]=c.z;
            };

            /// Assignment operator
            Matrix3D<TYPE>& operator=(const Matrix3D<TYPE>& a)
            {
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        x[i][j]=a.x[i][j];
                return *this;
            };

            /// Addition of two matrices
            Matrix3D<TYPE> operator+(const Matrix3D<TYPE>& a) const
            {
                Matrix3D<TYPE> result;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result.x[i][j]=x[i][j]+a.x[i][j];
                return result;
            };
            Matrix3D<TYPE> operator+=(const Matrix3D<TYPE>& a) {
                return (*this = *this + a);
            }

            /// Subtraction between two matrices
            Matrix3D<TYPE> operator-(const Matrix3D<TYPE>& a) const
            {
                Matrix3D<TYPE> result;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result.x[i][j]=x[i][j]-a.x[i][j];
                return result;
            };
            Matrix3D<TYPE> operator-=(const Matrix3D<TYPE>& a) {
                return (*this = *this - a);
            }

            /// Unary minus sign
            Matrix3D<TYPE> operator-() const
            {
                Matrix3D<TYPE> result;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result.x[i][j]=-x[i][j];
                return result;
            };

            /// Multiplication of two matrices
            Matrix3D<TYPE> operator*(const Matrix3D<TYPE>& a) const
            {
                Matrix3D<TYPE> result;
                int i,j;
                for (i=0; i<3; ++i)
                    for (j=0; j<3; ++j)
                        result.x[i][j] = 0.;
                for (i=0; i<3; ++i)
                    for (j=0; j<3; ++j)
                        for (int k=0; k<3; ++k)
                            result.x[i][k]+=x[i][j]*a.x[j][k];
                return result;
            }
            Matrix3D<TYPE> operator*=(const Matrix3D<TYPE>& a) {
                return (*this = *this * a);
            }


            /// Left multiplication of a vector by a matrix
            Vector3D<TYPE> operator*(const Vector3D<TYPE>& a) const
            {
                return Vector3D<TYPE>(x[0][0]*a.x+x[0][1]*a.y+x[0][2]*a.z,
                                      x[1][0]*a.x+x[1][1]*a.y+x[1][2]*a.z,
                                      x[2][0]*a.x+x[2][1]*a.y+x[2][2]*a.z);
            };

            /// Right multiplication of a vector by a matrix
            friend Vector3D<TYPE> operator*(const Vector3D<TYPE>& b,const Matrix3D<TYPE>& a)
            {
                return Vector3D<TYPE>(a.x[0][0]*b.x + a.x[1][0]*b.y + a.x[2][0]*b.z,
                                      a.x[0][1]*b.x + a.x[1][1]*b.y + a.x[2][1]*b.z,
                                      a.x[0][2]*b.x + a.x[1][2]*b.y + a.x[2][2]*b.z);
            };

            /// Multiplication of a matrix by a scalar:
            Matrix3D<TYPE> operator*(TYPE b) const
            {
                Matrix3D<TYPE> result;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result.x[i][j]=x[i][j]*b;
                return result;
            };
            Matrix3D<TYPE> operator*=(const TYPE& a) {
                return (*this = *this * a);
            }

            /// Multiplication of a matrix by a scalar:
            friend Matrix3D<TYPE> operator*(TYPE b,const Matrix3D<TYPE>& a)
            {
                return a*b;
            };

            /// Division of matrix by a vector:
            Matrix3D<TYPE> operator/(TYPE b) const {
                return (*this)*(1./b);
            };
            Matrix3D<TYPE> operator/=(const TYPE& b) {
                return (*this = *this / b);
            }

            Matrix3D<TYPE> transpose() const {
                Matrix3D<TYPE> a;
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        a[j][i]=x[i][j];
                    }
                }
                return a;
            }

            TYPE determinant() const {
                return -x[0][2]*x[1][1]*x[2][0] + x[0][1]*x[1][2]*x[2][0] +
                       x[0][2]*x[1][0]*x[2][1] - x[0][0]*x[1][2]*x[2][1] -
                       x[0][1]*x[1][0]*x[2][2] + x[0][0]*x[1][1]*x[2][2];
            }

            Matrix3D<TYPE> inverse() const {
                Matrix3D<TYPE> a;
                a[0][0] = x[1][1]*x[2][2] - x[1][2]*x[2][1];
                a[0][1] = x[0][2]*x[2][1] - x[0][1]*x[2][2];
                a[0][2] = x[0][1]*x[1][2] - x[0][2]*x[1][1];

                a[1][0] = x[1][2]*x[2][0] - x[1][0]*x[2][2];
                a[1][1] = x[0][0]*x[2][2] - x[0][2]*x[2][0];
                a[1][2] = x[0][2]*x[1][0] - x[0][0]*x[1][2];

                a[2][0] = x[1][0]*x[2][1] - x[1][1]*x[2][0];
                a[2][1] = x[0][1]*x[2][0] - x[0][0]*x[2][1];
                a[2][2] = x[0][0]*x[1][1] - x[0][1]*x[1][0];
                return a/determinant();
            }
    };

    /// The outer product of two vectors is a matrix:
    template <class TYPE>
    Matrix3D<TYPE> outer(const Vector3D<TYPE>& a,const Vector3D<TYPE>& b)
    {
        Matrix3D<TYPE> result;
        result[0][0]= a.x*b.x;
        result[0][1]= a.x*b.y;
        result[0][2]= a.x*b.z;
        result[1][0]= a.y*b.x;
        result[1][1]= a.y*b.y;
        result[1][2]= a.y*b.z;
        result[2][0]= a.z*b.x;
        result[2][1]= a.z*b.y;
        result[2][2]= a.z*b.z;
        return result;
    };

    /// A real matrix
    typedef Matrix3D<double> Matrix;
    /// A complex matrix
    typedef Matrix3D< COMPLEX > CMatrix;

} // namespace SCATMECH


#endif // SCATMECH_MATRIX3D_H

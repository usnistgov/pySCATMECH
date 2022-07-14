//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: mueller.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_MUELLER_HPP
#define SCATMECH_MUELLER_HPP

#include "scatmech.h"
#include "optconst.h"
#include "vector3d.h"

//*****************************************************************************
//**
//** Some attempt was made to include most of the terminology found in
//**    R. A. Chipman, "Polarimetry", Chapter 22 of _Handbook_of_Optics_,
//**    Volume II, (McGraw-Hill, New York, 1995)
//**
//*****************************************************************************
//
// Class declarations...
//

namespace SCATMECH {



    class JonesMatrix;
    class JonesVector;
    class MuellerMatrix;
    class StokesVector;

    /**
     * @brief Class to handle Jones vectors
     *
     * Class to handle data storage and the various operations
     *           associated with Jones vectors
     */
    class JonesVector {
        public:

            ///
            /// Copy constructor
            ///
            JonesVector() {}
            JonesVector(const JonesVector& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
            }
            ///
            /// Constructor taking two COMPLEX values
            ///
            JonesVector(const COMPLEX& s,const COMPLEX& p)
            {
                j[0]=s;
                j[1]=p;
            }
            ///
            /// Conversion constructor from StokesVector. The conversion from
            /// StokesVector to JonesVector uses only the polarized part of the
            /// JonesVector and returns a real s-component.
            ///
            explicit JonesVector(const StokesVector& x);

            ///
            /// Assignment operator
            ///
            JonesVector& operator=(const JonesVector& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
                return *this;
            }

            ///
            /// Multiplication by a constant
            ///
            JonesVector operator*(const COMPLEX& x) const
            {
                return JonesVector(j[0]*x,j[1]*x);
            }
            friend JonesVector operator*(const COMPLEX& x,const JonesVector& y)
            {
                return y*x;
            }
            JonesVector& operator*=(const COMPLEX& x)
            {
                return *this = x*(*this);
            }
            ///
            /// Division by a constant
            ///
            JonesVector operator/(const COMPLEX& x) const
            {
                return (1./x)*(*this);
            }
            JonesVector& operator/=(const COMPLEX& x)
            {
                return *this = (1./x)*(*this);
            }

            ///
            /// Product of two JonesVectors. This would normally
            /// only make physical sense for conj(J1)*J2.
            ///
            COMPLEX operator*(const JonesVector& v) const
            {
                return j[0]*v[0]+j[1]*v[1];
            }
            ///
            /// Complex conjugate of a JonesVector
            ///
            friend JonesVector conj(const JonesVector& v)
            {
                return JonesVector(conj(v[0]),conj(v[1]));
            }

            ///
            /// Addition and subtraction operators...
            ///
            JonesVector operator+(const JonesVector& a) const
            {
                return JonesVector(j[0]+a.j[0],j[1]+a.j[1]);
            }
            JonesVector operator-(const JonesVector& a) const
            {
                return JonesVector(j[0]-a.j[0],j[1]-a.j[1]);
            }
            JonesVector operator+=(const JonesVector& a)
            {
                return *this = *this + a;
            }
            JonesVector operator-=(const JonesVector& a)
            {
                return *this = *this - a;
            }
            JonesVector operator-() const {
                return JonesVector(-j[0],-j[1]);
            }
            const JonesVector& operator+() const {
                return *this;
            }

            ///
            /// Element indexing
            ///
            COMPLEX& operator[](int i) {
                return j[i];
            }
            const COMPLEX& operator[](int i) const {
                return j[i];
            }

            ///
            /// Return intensity
            ///
            double intensity() const {
                return std::norm(j[0])+std::norm(j[1]);
            }

            ///
            /// arctan of the component ratio
            ///
            double psi() const
            {
                return atan(std::abs(j[1])/std::abs(j[0]));
            }

            ///
            /// phase between two components
            ///
            double delta() const
            {
                return std::arg(j[1])-std::arg(j[0]);
            }

            ///
            /// principal angle of polarization
            ///
            double eta() const;

            ///
            /// degree of linear polarization
            ///
            double DOLP() const;

            ///
            /// degree of polarization (always 1 for Jones!)
            ///
            double DOP() const
            {
                return 1.;
            }

            ///
            /// degree of circular polarization
            ///
            double DOCP() const;

            /// ellipticity (ratio of minor to major axes)
            double e() const;

            ///
            /// eccentricity [sqrt(1-e^2)]
            ///
            double epsilon() const
            {
                double t = e();
                return sqrt(1.-t*t);
            }

            ///
            /// Vector rotated by angle...
            ///
            JonesVector rotate(const double angle) const;

            ///
            /// S component of wave
            ///
            COMPLEX& S() {
                return j[0];
            }
            const COMPLEX& S() const {
                return j[0];
            }

            ///
            /// P component of wave
            ///
            COMPLEX& P() {
                return j[1];
            }
            const COMPLEX& P() const {
                return j[1];
            }

            friend class JonesMatrix;
            friend class MuellerMatrix;
            friend class StokesVector;

        private:
            ///
            /// The elements of the vector
            ///
            COMPLEX j[2];
    };

    /**
     * @brief Class to handle Jones Matrices
     *
     * Class to store the elements and handle operations of a Jones matrix.  The elements are stored
     * in the order {pp,ss,ps,sp}, where ps and sp are, counter to normal matrix notation, p->s and s->p, respectively.
     */
    class JonesMatrix {
        public:
            ///
            /// Default constructor leaves elements unassigned
            ///
            JonesMatrix() {}

            ///
            /// Copy constructor
            ///
            JonesMatrix(const JonesMatrix& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
                j[2]=x.j[2];
                j[3]=x.j[3];
            }

            ///
            /// Constructor from four values
            ///
            JonesMatrix(const COMPLEX& pp,const COMPLEX& ss,
                        const COMPLEX& ps,const COMPLEX& sp)
            {
                j[0]=pp;
                j[1]=ss;
                j[2]=ps,
                j[3]=sp;
            }

            ///
            /// @brief Conversion from a Jones matrix to a Mueller matrix.
            ///
            /// Conversion from a Jones matrix to a Mueller matrix using the
            /// expression given in R.A. Chipman, "Polarimetry" in _Handbook_of_Optics_Vol._II_
            /// (McGraw-Hill, New York, 1995)
            ///
            explicit JonesMatrix(const MuellerMatrix& x);

            ///
            /// Assignment operator
            ///
            JonesMatrix& operator=(const JonesMatrix& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
                j[2]=x.j[2];
                j[3]=x.j[3];
                return *this;
            }

            ///
            /// Multiplication of two Jones matrices
            ///
            JonesMatrix operator*(const JonesMatrix& matrix) const
            {
                return JonesMatrix(j[3]*matrix.j[2]+j[0]*matrix.j[0],
                                   j[1]*matrix.j[1]+j[2]*matrix.j[3],
                                   j[1]*matrix.j[2]+j[2]*matrix.j[0],
                                   j[3]*matrix.j[1]+j[0]*matrix.j[3]);
            }
            JonesMatrix& operator*=(const JonesMatrix& a)
            {
                return *this = *this * a;
            }

            ///
            /// Multiplication of a Jones matrix by a constant
            ///
            JonesMatrix operator*(const COMPLEX& x) const
            {
                return JonesMatrix(j[0]*x,j[1]*x,j[2]*x,j[3]*x);
            }
            friend JonesMatrix operator*(const COMPLEX& x,const JonesMatrix& y)
            {
                return y*x;
            }
            JonesMatrix& operator*=(const COMPLEX& a) {
                return *this = *this * a;
            }

            ///
            /// Division of a Jones matrix by a constant
            ///
            JonesMatrix operator/(const COMPLEX& x) const
            {
                return JonesMatrix(j[0]/x,j[1]/x,j[2]/x,j[3]/x);
            }
            JonesMatrix& operator/=(const COMPLEX& a) {
                return *this = *this / a;
            }

            ///
            /// Multiplication of a Jones matrix and a Jones vector
            ///
            JonesVector operator*(const JonesVector& a) const
            {
                return JonesVector(j[1]*a.j[0]+j[2]*a.j[1],j[3]*a.j[0]+j[0]*a.j[1]);
            }

            ///
            /// Addition of two Jones matrices
            ///
            JonesMatrix operator+(const JonesMatrix& a) const
            {
                return JonesMatrix(j[0]+a.j[0],j[1]+a.j[1],j[2]+a.j[2],j[3]+a.j[3]);
            }
            JonesMatrix operator-(const JonesMatrix& a) const
            {
                return JonesMatrix(j[0]-a.j[0],j[1]-a.j[1],j[2]-a.j[2],j[3]-a.j[3]);
            }
            JonesMatrix& operator+=(const JonesMatrix& a)
            {
                return *this = *this + a;
            }
            JonesMatrix& operator-=(const JonesMatrix& a)
            {
                return *this = *this - a;
            }
            JonesMatrix operator-() const {
                return JonesMatrix(-j[0],-j[1],-j[2],-j[3]);
            }
            const JonesMatrix& operator+() const {
                return *this;
            }

            ///
            /// Element indexing
            ///
            COMPLEX& operator[](int i) {
                return j[i];
            };
            const COMPLEX& operator[](int i) const {
                return j[i];
            };

            ///
            /// Return matrix rotated by angle
            ///
            JonesMatrix rotate(double angle) const;

            ///
            /// Return matrix transpose
            ///
            JonesMatrix transpose() const {
                return JonesMatrix(j[0],j[1],j[3],j[2]);
            }
            friend JonesMatrix transpose(const JonesMatrix& j) {
                return JonesMatrix(j[0],j[1],j[3],j[2]);
            }

            ///
            /// Return the Hermitian transpose
            ///
            JonesMatrix hermitian() const {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[3]), std::conj(j[2]));
            }
            friend JonesMatrix hermitian(const JonesMatrix& j) {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[3]), std::conj(j[2]));
            }

            ///
            /// Return the complex conjugate
            ///
            JonesMatrix conj() const {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[2]), std::conj(j[3]));
            }
            friend JonesMatrix conj(const JonesMatrix& j) {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[2]), std::conj(j[3]));
            }

            ///
            /// Public access to p->p matrix element
            ///
            COMPLEX& PP() {
                return j[0];
            }
            const COMPLEX& PP() const {
                return j[0];
            }

            ///
            /// Public access to s->s matrix element
            ///
            COMPLEX& SS() {
                return j[1];
            }
            const COMPLEX& SS() const {
                return j[1];
            }

            ///
            /// Public access to p->s matrix element
            ///
            COMPLEX& PS() {
                return j[2];
            }
            const COMPLEX& PS() const {
                return j[2];
            }

            ///
            /// Public access to s->p matrix element
            ///
            COMPLEX& SP() {
                return j[3];
            }
            const COMPLEX& SP() const {
                return j[3];
            }


            friend class JonesVector;
            friend class MuellerMatrix;
            friend class StokesVector;
        private:

            /// @brief Array of matrix elements
            ///
            /// Array of matrix elements. Note that the elements are p->p, s->s, p->s, and s->p, in order.
            COMPLEX j[4];
    };

    ///
    /// @brief Class to handle Mueller matrices
    ///
    /// Class to handle the storage of the elements and operations of a Mueller matrix
    ///
    class MuellerMatrix {

        public:

            ///
            /// Default constructor. Values are left unassigned.
            ///
            MuellerMatrix() {}

            ///
            /// Constructor from a Jones matrix
            ///
            MuellerMatrix(const JonesMatrix& j);

            ///
            /// Copy constructor
            ///
            MuellerMatrix(const MuellerMatrix& x) {
                *this=x;
            }

            ///
            /// @brief Constructor from four Stokes vectors.
            ///
            /// Constructor from four Stokes vectors. The first row of the matrix is s1, the second row is s2, etc.
            ///
            MuellerMatrix(const StokesVector& s1,const StokesVector& s2,const StokesVector& s3, const StokesVector& s4);

            ///
            /// Assignment operator
            ///
            MuellerMatrix& operator=(const MuellerMatrix& x);

            //
            // Public access to the elements of the Mueller matrix
            //
            double* operator[](int i) {
                return m[i];
            }
            const double* operator[](int i) const {
                return m[i];
            }

            ///
            /// Maximum "transmission"
            ///
            double Tmax() const
            {
                return m[0][0]+sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[0][3]*m[0][3]);
            }
            ///
            /// Minimum "transmission"
            ///
            double Tmin() const
            {
                return m[0][0]-sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[0][3]*m[0][3]);
            }
            ///
            /// Diattenuation
            ///
            double diattenuation() const
            {
                return sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[0][3]*m[0][3])/m[0][0];
            }
            ///
            /// Linear diattenuation
            ///
            double linear_diattenuation() const
            {
                return sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2])/m[0][0];
            }
            ///
            /// Polarization dependent loss
            ///
            double polarization_dependent_loss() const
            {
                return 10.*std::log(Tmax()/Tmin())/std::log(10.);
            }
            ///
            /// Polarizance
            ///
            double polarizance() const
            {
				return sqrt(m[1][0] * m[1][0] + m[2][0] * m[2][0] + m[3][0] * m[3][0]) / m[0][0];
            }
            ///
            /// Depolarization index
            ///
            double depolarization_index() const
            {
                MuellerMatrix mmt = (*this)*(this->transpose());
                return sqrt((mmt[0][0]+mmt[1][1]+mmt[2][2]+mmt[3][3]-sqr(m[0][0]))/(3.*sqr(m[0][0])));
            }
            ///
            /// Extinction ratio
            ///
            double extinction_ratio() const {
                return Tmax()/Tmin();
            }

            ///
            /// Return matrix rotated by angle
            ///
            MuellerMatrix rotate(const double angle) const;

            ///
            /// Return scattering matrix for a parity conversion
            ///
            MuellerMatrix parity() const;

            ///
            /// Return transpose of matrix
            ///
            MuellerMatrix transpose() const;

            ///
            /// Return the inverse of the matrix
            ///
            MuellerMatrix inverse() const;

            ///
            /// Returns true for a Mueller matrix that is a mapping of the space of Stokes vectors into itself.
            ///
            bool valid() const;

            /// @brief Return true if matrix is physically valid, otherwise false.
            ///
            /// Return true if matrix is physically valid, otherwise false.
            /// physically_valid() uses the criterion that the Mueller matrix be a
            /// convex sum of Jones-Mueller matrices. This is a more restrictive
            /// condition than valid()...
            ///
            bool physically_valid() const;

            ///
            /// The matrix logarithm of a Mueller matrix
            ///
            MuellerMatrix log() const;

            ///
            /// The matrix exponential of a Mueller matrix
            ///
            MuellerMatrix exp() const;

            ///
            /// The Cloude decomposition, which expresses the Mueller matrix as the sum of four non-depolarizing Mueller matrices
            ///
            void Cloude_Decomposition(MuellerMatrix& M1,MuellerMatrix& M2,MuellerMatrix& M3,MuellerMatrix& M4) const;

            ///
            /// Lu-Chipman decomposition, which expresses the Mueller matrix as a product of a depolarizer, a retarder, and a diattenuator
            ///
            void Lu_Chipman_Decomposition(MuellerMatrix& depolarizer,MuellerMatrix& retarder, MuellerMatrix& diattenuator) const;

            ///
            /// Returns the closest Mueller matrix that is non-depolarizing
            ///
            MuellerMatrix Closest_NonDepolarizing(int rank=1) const;

            ///
            /// The Mueller matrix normalized, M/M[0][0]
            ///
            MuellerMatrix normalized() const {
                return *this/m[0][0];
            }

			///
			/// The polarimetric entropy, defined by Cloude and Pottier
			///
			double entropy() const;

            ///
            /// Multiplication of two Mueller matrices
            ///
            MuellerMatrix operator*(const MuellerMatrix& matrix) const;
            MuellerMatrix& operator*=(const MuellerMatrix& matrix)
            {
                return *this = *this * matrix;
            }

            ///
            /// Multiplication by a StokesVector
            ///
            StokesVector operator*(const StokesVector &s) const;

            ///
            /// Multiplcation by a constant
            ///
            MuellerMatrix operator*(double d) const;
            friend MuellerMatrix operator*(double d,const MuellerMatrix &v)
            {
                return v*d;
            }
            MuellerMatrix& operator*=(const double d) {
                return *this = *this * d;
            }

            ///
            /// Division by a constant
            ///
            MuellerMatrix operator/(double d) const;
            MuellerMatrix& operator/=(double d) {
                return *this = *this * (1./d);
            }

            ///
            /// Addition of two Mueller matrices
            ///
            MuellerMatrix operator+(const MuellerMatrix& a) const;
            MuellerMatrix& operator+=(const MuellerMatrix& matrix)
            {
                return *this = *this + matrix;
            }

            ///
            /// Subtraction of two Mueller matrices
            ///
            MuellerMatrix operator-(const MuellerMatrix& a) const;
            MuellerMatrix& operator-=(const MuellerMatrix& matrix)
            {
                return *this = *this - matrix;
            }

            MuellerMatrix operator-() const;
            const MuellerMatrix& operator+() const {
                return *this;
            }

            friend class JonesVector;
            friend class JonesMatrix;
            friend class StokesVector;

        private:

            ///
            /// The elements of the Mueller matrix...
            ///
            double m[4][4];
    };

    ///
    /// @brief Class to handle Stokes vectors
    ///
    /// Class to store the elements and handle the operations of a Stokes vector
    ///
    class StokesVector {

        public:

            ///
            /// Default constructor leaves elements unassigned.
            ///
            StokesVector() {}

            ///
            /// Constructor from four values
            ///
            StokesVector(double I,double Q,double U,double V)
            {
                s[0]=I;
                s[1]=Q;
                s[2]=U;
                s[3]=V;
            }

            ///
            /// Copy constructor
            ///
            StokesVector(const StokesVector& x)
            {
                s[0]=x.s[0];
                s[1]=x.s[1];
                s[2]=x.s[2];
                s[3]=x.s[3];
            }

            ///
            /// Conversion from a Jones vector to a Stokes vector
            ///
            StokesVector(const JonesVector& j)
            {
                s[0] = std::norm(j.j[0])+std::norm(j.j[1]);
                s[1] = std::norm(j.j[0])-std::norm(j.j[1]);
                s[2] = 2*std::real(std::conj(j.j[0])*j.j[1]);  // (s+p) minus (s-p)
                s[3] = 2*std::imag(std::conj(j.j[0])*j.j[1]);  // Left minus Right
            }

            ///
            /// Assignment operator
            ///
            StokesVector& operator=(const StokesVector& x)
            {
                s[0]=x.s[0];
                s[1]=x.s[1];
                s[2]=x.s[2];
                s[3]=x.s[3];
                return *this;
            }

            ///
            /// The first element (intensity) of the Stokes vector
            ///
            double& I() {
                return s[0];
            }
            double I() const {
                return s[0];
            }

            ///
            /// The second element (Is - Ip) of the Stokes vector
            ///
            double& Q() {
                return s[1];
            }
            double Q() const {
                return s[1];
            }

            ///
            /// The third element [I(s+p) - I(s-p)] of the Stokes vector
            ///
            double& U() {
                return s[2];
            }
            double U() const {
                return s[2];
            }

            ///
            /// The fourth element [I(lcp) - I(rcp)] of the Stokes vector
            ///
            double& V() {
                return s[3];
            }
            double V() const {
                return s[3];
            }

            ///
            /// The normalized Q, Q/I
            ///
            double q() const {
                return Q()/I();
            }

            ///
            /// The normalized U, U/I
            ///
            double u() const {
                return U()/I();
            }

            ///
            /// The normalized V, V/I
            ///
            double v() const {
                return V()/I();
            }

            ///
            /// Public access to the Stokes vector elements
            ///
            double& operator[](int i) {
                return s[i];
            }
            double operator[](int i) const {
                return s[i];
            }

            ///
            /// Right multiplication by a Mueller matrix
            ///
            StokesVector operator*(const MuellerMatrix& matrix) const;

            ///
            /// Inner product with another Stokes vector
            ///
            double operator*(const StokesVector& a) const
            {
                return a.s[0]*s[0] + a.s[1]*s[1] + a.s[2]*s[2] + a.s[3]*s[3];
            }

            ///
            /// Addition of Stokes vectors
            ///
            StokesVector operator+(const StokesVector& a) const
            {
                return StokesVector(s[0]+a.s[0],s[1]+a.s[1],s[2]+a.s[2],s[3]+a.s[3]);
            }
            StokesVector operator+=(const StokesVector& a)
            {
                return *this = *this+a;
            }

            ///
            /// Subtraction of Stokes vectors
            ///
            StokesVector operator-(const StokesVector& a) const
            {
                return StokesVector(s[0]-a.s[0],s[1]-a.s[1],s[2]-a.s[2],s[3]-a.s[3]);
            }
            StokesVector operator-=(const StokesVector& a)
            {
                return *this = *this-a;
            }

            ///
            /// Negation of a Stokes vector
            ///
            StokesVector operator-() const
            {
                return StokesVector(-s[0],-s[1],-s[2],-s[3]);
            }
            const StokesVector& operator+() const {
                return *this;
            }

            ///
            /// Multiplication by a constant
            ///
            StokesVector operator*(double d) const
            {
                return StokesVector(s[0]*d,s[1]*d,s[2]*d,s[3]*d);
            }
            friend StokesVector operator*(double d,const StokesVector& s)
            {
                return StokesVector(s[0]*d,s[1]*d,s[2]*d,s[3]*d);
            }
            StokesVector& operator*=(double d) {
                return *this = *this *  d;
            }

            ///
            /// Division by a constant
            ///
            StokesVector operator/(double d) const
            {
                return StokesVector(s[0]/d,s[1]/d,s[2]/d,s[3]/d);
            }
            StokesVector& operator/=(double d) {
                return *this = *this /  d;
            }

            ///
            /// Return vector rotated by angle...
            ///
            StokesVector rotate(double angle) const;

            ///
            /// Principle angle of polarization
            ///
            double eta() const
            {
                return atan2(s[2],s[1])/2.;
            }

            ///
            /// Intensity
            ///
            double intensity() const
            {
                return s[0];
            }

            ///
            /// Degree of linear polarization
            ///
            double DOLP() const
            {
                return sqrt(sqr(s[2])+sqr(s[1]))/s[0];
            }

            ///
            /// Degree of polarization
            ///
            double DOP() const
            {
                return sqrt(sqr(s[2])+sqr(s[1])+sqr(s[3]))/s[0];
            }

            ///
            /// Degree of circular polarization
            ///
            double DOCP() const
            {
                return s[3]/s[0];
            }

            ///
            /// Ellipticity (ratio of minor to major axes)
            ///
            double e() const
            {
                return s[3]/(s[0]+sqrt(sqr(s[1])+sqr(s[2])));
            }

            ///
            /// Phase between two components
            ///
            double delta() const;

            ///
            /// Arctangent of the component ratio
            ///
            double psi() const;

            ///
            /// eccentricity [sqrt(1-e^2)]
            ///
            double eccentricity() const
            {
                return sqrt(1.-sqr(e()));
            }

            ///
            /// Returns true if Stokes vector is valid, false otherwise...
            ///
            bool valid() const
            {
                return s[0]>=sqrt(sqr(s[1])+sqr(s[2])+sqr(s[3]));
            }

            ///
            /// Polarized part of Stokes Vector
            ///
            StokesVector pol_part() const;

            ///
            /// Unpolarized part of Stokes Vector
            ///
            StokesVector unpol_part() const;

            friend class JonesVector;
            friend class MuellerMatrix;
            friend class JonesMatrix;

        private:

            ///
            /// The elements of the Stokes vector
            ///
            double s[4];
    };

    ///
    /// Stream insertion and extraction operators
    ///
    std::ostream& operator<<(std::ostream& os,const JonesVector& j);
    std::istream& operator>>(std::istream& is,JonesVector& j);
    std::ostream& operator<<(std::ostream& os,const StokesVector& j);
    std::istream& operator>>(std::istream& is, StokesVector& j);
    std::ostream& operator<<(std::ostream& os,const MuellerMatrix& mm);
    std::istream& operator>>(std::istream& is, MuellerMatrix& mm);
    std::ostream& operator<<(std::ostream& os,const JonesMatrix& j);
    std::istream& operator>>(std::istream& is, JonesMatrix& j);

    ///
    /// Returns a zero Jones matrix
    ///
    inline JonesMatrix JonesZero() {
        return JonesMatrix(0.,0.,0.,0.);
    }
    ///
    /// Returns a unit Jones matrix
    ///
    inline JonesMatrix JonesUnit() {
        return JonesMatrix(1.,1.,0.,0.);
    }

    ///
    /// Returns a Jones rotation matrix
    ///
    JonesMatrix JonesRotator(double angle);
    ///
    /// Returns a Jones matrix for a linear retarder
    ///
    JonesMatrix JonesLinearRetarder(double phase, double angle=0);
    ///
    /// Returns a Jones matrix for a circular retarder
    ///
    JonesMatrix JonesCircularRetarder(double phase);
    ///
    /// Returns a Jones matrix for a linear polarizer
    ///
    JonesMatrix JonesLinearPolarizer(double angle=0,double diattenuation=1);
    ///
    /// Returns a Jones matrix for a circular polarizer
    ///
    JonesMatrix JonesCircularPolarizer(double diattenuation=1);

    ///
    /// Returns Jones matrix given eigenvectors and eigenvalues
    ///
    JonesMatrix JonesGeneralized(
        const JonesVector& a, ///< The first eigenvector
        const JonesVector& b, ///< The second eigenvector
        const COMPLEX& ma,    ///< The eigenvalue associated with eigenvector a
        const COMPLEX& mb     ///< The eigenvalue associated with eigenvector b
    );

    ///
    /// Returns a zero Mueller matrix
    ///
    MuellerMatrix MuellerZero();
    ///
    /// Returns a unit Mueller matrix
    ///
    MuellerMatrix MuellerUnit(
        double attenuation=1.  ///< An overall scaling factor
    );
    ///
    /// Returns a simple Mueller matrix depolarizer
    ///
    MuellerMatrix MuellerDepolarizer(
        double attenuation=1.,  ///< An overal scaling factor
        double depolarization=1. ///< The amount of depolarization, 0<=depolarization<=1
    );
	// 
	// Returns a diagonal Mueller matrix
	//
	MuellerMatrix MuellerDiagonal(double m00,double m11, double m22, double m33);
    ///
    /// Returns a partial linear polarizer
    ///
    MuellerMatrix MuellerPartialLinearPolarizer(
        double tmax,   ///< The maximum transmittance
        double tmin,   ///< The minimum transmittance
        double angle   ///< The angle of the maximum transmittance
    );

    ///
    /// Returns a zero Stokes vector
    ///
    inline StokesVector StokesZero() {
        return StokesVector(0,0,0,0);
    }

    ///
    /// Returns a unit, unpolarized Stokes vector
    ///
    inline StokesVector StokesVectorUnitUnpol() {
        return StokesVector(1,0,0,0);
    }
    ///
    /// Returns a linearly s-polarized Stokes vector
    ///
    inline StokesVector StokesVectorUnitS() {
        return StokesVector(1,1,0,0);
    }
    ///
    /// Returns a linearly p-polarized Stokes vector
    ///
    inline StokesVector StokesVectorUnitP() {
        return StokesVector(1,-1,0,0);
    }
    ///
    /// Returns a linear polarized Stokes vector along angle eta
    ///
    inline StokesVector StokesVectorUnitLinear(double eta) {
        return StokesVector(1,cos(2*eta),sin(2*eta),0);
    }
    ///
    /// Returns a right-circularly-polarized Stokes vector
    ///
    inline StokesVector StokesVectorUnitRCP() {
        return StokesVector(1,0,0,-1);
    }
    ///
    /// Returns a left circularly polarized Stokes vector
    ///
    inline StokesVector StokesVectorUnitLCP() {
        return StokesVector(1,0,0,1);
    }
    ///
    /// Returns a generalized unit Stokes vector
    ///
    inline StokesVector StokesVectorUnitGeneral(
        double eta,     ///< The principal angle of the polarization
        double DOCP=0., ///< The degree of circular polarization, -1<=DOCP<=1
        double DOP=1.   ///< The degree of polarization, 0<=DOP<=1
    )
    {
        double temp=sqrt(1.-sqr(DOCP));
        return StokesVector(1,
                            DOP*temp*cos(2*eta),
                            DOP*temp*sin(2*eta),
                            DOP*DOCP);
    }

    ///
    /// Returns a zero Jones vector
    ///
    inline JonesVector JonesVectorZero() {
        return JonesVector(0.,0.);
    }
    ///
    /// Returns an s-polarized Jones vector
    ///
    inline JonesVector JonesVectorUnitS() {
        return JonesVector(1.,0.);
    }
    ///
    /// Returns a p-polarized Jones vector
    ///
    inline JonesVector JonesVectorUnitP() {
        return JonesVector(0.,1.);
    }
    ///
    /// Returns a linearly polarized Jones vector aligned along an angle
    ///
    inline JonesVector JonesVectorUnitLinear(double eta) {
        return JonesVector(cos(eta),sin(eta));
    }
    ///
    /// Returns a right-circularly-polarized Jones vector
    ///
    inline JonesVector JonesVectorUnitRCP() {
        return JonesVector(1.,COMPLEX(0,-1));
    }
    ///
    /// Returns a left-circularly-polarized Jones vector
    ///
    inline JonesVector JonesVectorUnitLCP() {
        return JonesVector(1.,COMPLEX(0,1));
    }
    ///
    /// Returns a generalized Jones vector
    ///
    inline JonesVector JonesVectorUnitGeneral(
        double eta, ///< Principal angle of the polarization
        double DOCP ///< Degree of circular polarization, -1<=DOCP<=1
    )
    {
        return JonesVector(StokesVectorUnitGeneral(eta,DOCP));
    }
    ///
    /// @brief Given a direction of propagation k, retrieves s and p directions
    ///
    /// Given a direction of propagation k, retrieves s and p directions. The s direction is
    /// perpendicular to the z axis and k.
    void GetBasisVectorsSP(const Vector& k,Vector& s,Vector& p);
    ///
    /// @brief Given a direction of propagation k, retrieves x and y directions
    ///
    /// Given a direction of propagation k, retrieves x and y directions. The x and y vectors are perpendicular
    /// to the direction of propagation. If k has no y component and a positive x component, x=p and y=s. The x-y basis
    /// set has no discontinuity when k is propagating along the z axis.
    void GetBasisVectorsXY(const Vector& k,Vector& x,Vector& y);
    ///
    /// @brief Retrieves basis vector appropriate for the scattering plane
    ///
    /// Given an incident direction kin and scattering direction kout, retrieves the vector perp perpendicular to the
    /// plane containing kin and kout, and vectors parin and parout perpendicular to perp and kin and kout, respectively.
    ///
    void GetBasisVectorsParPerp(const Vector& kin,const Vector& kout, Vector& perp, Vector& parin, Vector& parout);
    ///
    /// Returns a Jones matrix rotator that rotates from basis xi-yi to basis xo-yo
    ///
    JonesMatrix GetJonesRotator(const Vector& xo, const Vector& yo, const Vector& xi, const Vector& yi);
    ///
    /// @brief Returns the real part of the cross Mueller matrix
    ///
    /// Returns the real part of the cross Mueller matrix defined in Eq. (43)
    /// of T.A.Germer, "Measuring Interfacial Roughness by Polarized
    /// Optical Scattering," in "Light Scattering and Nanoscale Surface Roughness," Ed. by
    /// A.A. Maradudin, (Springer,New York, 2007).
    ///
    /// It is the polarimetric equivalent of conj(j1)*j2, which one gets from the product
    /// of (j1 + j2)*conj(j1 + j2) = j1*conj(j1) + j1*conj(j2) + j2*conj(j1) + j2*conj(j2),
    /// and use for similar expressions, when there is partial coherence between j1 and j2.
    ///
    /// ReCrossMueller(j1,j1) is the same as MuellerMatrix(j1).
    ///
    MuellerMatrix ReCrossMueller(const JonesMatrix& j1,const JonesMatrix& j2);
    ///
    /// @brief Returns the imaginary part of the cross Mueller matrix
    ///
    /// ImCrossMueller(j1,j1) is the zero matrix.
    ///
    MuellerMatrix ImCrossMueller(const JonesMatrix& j1,const JonesMatrix& j2);


} // namespace SCATMECH

#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: mueller.cpp
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
#include "matrixmath.h"
#include "vector3d.h"
#include "matrix3d.h"
#include <cmath>
#include <limits>

using namespace std;


namespace SCATMECH {


    //
    // Conversion from JonesMatrix to MuellerMatrix
    //
    MuellerMatrix::
    MuellerMatrix(const JonesMatrix& jj)
    {
        //
        // from
        // Bohren and Huffman
        // _Absorption_and_Scattering_of_Light_by_Small_Particles_,
        // (Wiley, New York, 1983).
        //
        //init();
        double norm0 = norm(jj.j[0]);
        double norm1 = norm(jj.j[1]);
        double norm2 = norm(jj.j[2]);
        double norm3 = norm(jj.j[3]);

        m[0][0] = 0.5*(norm0+norm1+norm2+norm3);
        m[0][1] = 0.5*(norm1-norm0+norm3-norm2);
        m[0][2] = real(jj.j[1]*conj(jj.j[2])+jj.j[0]*conj(jj.j[3]));
        m[0][3] = imag(jj.j[1]*conj(jj.j[2])-jj.j[0]*conj(jj.j[3]));
        m[1][0] = 0.5*(norm1-norm0-norm3+norm2);
        m[1][1] = 0.5*(norm1+norm0-norm3-norm2);
        m[1][2] = real(jj.j[1]*conj(jj.j[2])-jj.j[0]*conj(jj.j[3]));
        m[1][3] = imag(jj.j[1]*conj(jj.j[2])+jj.j[0]*conj(jj.j[3]));
        m[2][0] = real(jj.j[1]*conj(jj.j[3])+jj.j[0]*conj(jj.j[2]));
        m[2][1] = real(jj.j[1]*conj(jj.j[3])-jj.j[0]*conj(jj.j[2]));
        m[2][2] = real(jj.j[0]*conj(jj.j[1])+jj.j[2]*conj(jj.j[3]));
        m[2][3] = imag(jj.j[1]*conj(jj.j[0])+jj.j[3]*conj(jj.j[2]));
        m[3][0] = imag(jj.j[3]*conj(jj.j[1])+jj.j[0]*conj(jj.j[2]));
        m[3][1] = imag(jj.j[3]*conj(jj.j[1])-jj.j[0]*conj(jj.j[2]));
        m[3][2] = imag(jj.j[0]*conj(jj.j[1])-jj.j[2]*conj(jj.j[3]));
        m[3][3] = real(jj.j[0]*conj(jj.j[1])-jj.j[2]*conj(jj.j[3]));

    }

    //
    // The following two functions return the real and imaginary parts of Mueller product
    // defined in Eq. (43) of T.A.Germer, "Measuring Interfacial Roughness by Polarized
    // Optical Scattering," in "Light Scattering and Nanoscale Surface Roughness," Ed. by
    // A.A. Maradudin, (Springer,New York, 2007).
    //
    // It is the polarimetric equivalent of conj(j1)*j2, which one gets from the product
    // of (j1 + j2)*conj(j1 + j2) = j1*conj(j1) + j1*conj(j2) + j2*conj(j1) + j2*conj(j2),
    // and use for similar expressions, when there is partial coherence between j1 and j2.
    //
    // ReCrossMueller(j1,j1) is the same as MuellerMatrix(j1).
    //
    MuellerMatrix ReCrossMueller(const JonesMatrix& j1,const JonesMatrix& j2)
    {
        COMPLEX j2PPj1PP = j2.PP()*conj(j1.PP());
        COMPLEX j2SSj1SS = j2.SS()*conj(j1.SS());
        COMPLEX j2PSj1PS = j2.PS()*conj(j1.PS());
        COMPLEX j2SPj1SP = j2.SP()*conj(j1.SP());
        COMPLEX j2PSj1PP = j2.PS()*conj(j1.PP());
        COMPLEX j2PPj1PS = j2.PP()*conj(j1.PS());
        COMPLEX j2SSj1SP = j2.SS()*conj(j1.SP());
        COMPLEX j2SPj1SS = j2.SP()*conj(j1.SS());
        COMPLEX j2SPj1PP = j2.SP()*conj(j1.PP());
        COMPLEX j2SSj1PS = j2.SS()*conj(j1.PS());
        COMPLEX j2PPj1SP = j2.PP()*conj(j1.SP());
        COMPLEX j2PSj1SS = j2.PS()*conj(j1.SS());
        COMPLEX j2SSj1PP = j2.SS()*conj(j1.PP());
        COMPLEX j2SPj1PS = j2.SP()*conj(j1.PS());
        COMPLEX j2PSj1SP = j2.PS()*conj(j1.SP());
        COMPLEX j2PPj1SS = j2.PP()*conj(j1.SS());

        MuellerMatrix result;

        result[0][0] = ( real(j2PPj1PP) + real(j2PSj1PS) + real(j2SPj1SP) + real(j2SSj1SS))/2.;
        result[1][0] = (-real(j2PPj1PP) + real(j2PSj1PS) - real(j2SPj1SP) + real(j2SSj1SS))/2.;
        result[2][0] = ( real(j2PSj1PP) + real(j2PPj1PS) + real(j2SSj1SP) + real(j2SPj1SS))/2.;
        result[3][0] = (-imag(j2PSj1PP) + imag(j2PPj1PS) - imag(j2SSj1SP) + imag(j2SPj1SS))/2.;
        result[0][1] = (-real(j2PPj1PP) - real(j2PSj1PS) + real(j2SPj1SP) + real(j2SSj1SS))/2.;
        result[1][1] = ( real(j2PPj1PP) - real(j2PSj1PS) - real(j2SPj1SP) + real(j2SSj1SS))/2.;
        result[2][1] = (-real(j2PSj1PP) - real(j2PPj1PS) + real(j2SSj1SP) + real(j2SPj1SS))/2.;
        result[3][1] = ( imag(j2PSj1PP) - imag(j2PPj1PS) - imag(j2SSj1SP) + imag(j2SPj1SS))/2.;
        result[0][2] = ( real(j2SPj1PP) + real(j2SSj1PS) + real(j2PPj1SP) + real(j2PSj1SS))/2.;
        result[1][2] = (-real(j2SPj1PP) + real(j2SSj1PS) - real(j2PPj1SP) + real(j2PSj1SS))/2.;
        result[2][2] = ( real(j2SSj1PP) + real(j2SPj1PS) + real(j2PSj1SP) + real(j2PPj1SS))/2.;
        result[3][2] = (-imag(j2SSj1PP) + imag(j2SPj1PS) - imag(j2PSj1SP) + imag(j2PPj1SS))/2.;
        result[0][3] = ( imag(j2SPj1PP) + imag(j2SSj1PS) - imag(j2PPj1SP) - imag(j2PSj1SS))/2.;
        result[1][3] = (-imag(j2SPj1PP) + imag(j2SSj1PS) + imag(j2PPj1SP) - imag(j2PSj1SS))/2.;
        result[2][3] = ( imag(j2SSj1PP) + imag(j2SPj1PS) - imag(j2PSj1SP) - imag(j2PPj1SS))/2.;
        result[3][3] = ( real(j2SSj1PP) - real(j2SPj1PS) - real(j2PSj1SP) + real(j2PPj1SS))/2.;
        
        return result;
    }

    MuellerMatrix ImCrossMueller(const JonesMatrix& j1,const JonesMatrix& j2)
    {
        COMPLEX j2PPj1PP = j2.PP()*conj(j1.PP());
        COMPLEX j2SSj1SS = j2.SS()*conj(j1.SS());
        COMPLEX j2PSj1PS = j2.PS()*conj(j1.PS());
        COMPLEX j2SPj1SP = j2.SP()*conj(j1.SP());
        COMPLEX j2PSj1PP = j2.PS()*conj(j1.PP());
        COMPLEX j2PPj1PS = j2.PP()*conj(j1.PS());
        COMPLEX j2SSj1SP = j2.SS()*conj(j1.SP());
        COMPLEX j2SPj1SS = j2.SP()*conj(j1.SS());
        COMPLEX j2SPj1PP = j2.SP()*conj(j1.PP());
        COMPLEX j2SSj1PS = j2.SS()*conj(j1.PS());
        COMPLEX j2PPj1SP = j2.PP()*conj(j1.SP());
        COMPLEX j2PSj1SS = j2.PS()*conj(j1.SS());
        COMPLEX j2SSj1PP = j2.SS()*conj(j1.PP());
        COMPLEX j2SPj1PS = j2.SP()*conj(j1.PS());
        COMPLEX j2PSj1SP = j2.PS()*conj(j1.SP());
        COMPLEX j2PPj1SS = j2.PP()*conj(j1.SS());

        MuellerMatrix result;

        result[0][0] = ( imag(j2PPj1PP) + imag(j2PSj1PS) + imag(j2SPj1SP) + imag(j2SSj1SS))/2.;
        result[1][0] = (-imag(j2PPj1PP) + imag(j2PSj1PS) - imag(j2SPj1SP) + imag(j2SSj1SS))/2.;
        result[2][0] = ( imag(j2PSj1PP) + imag(j2PPj1PS) + imag(j2SSj1SP) + imag(j2SPj1SS))/2.;
        result[3][0] = ( real(j2PSj1PP) - real(j2PPj1PS) + real(j2SSj1SP) - real(j2SPj1SS))/2.;
        result[0][1] = (-imag(j2PPj1PP) - imag(j2PSj1PS) + imag(j2SPj1SP) + imag(j2SSj1SS))/2.;
        result[1][1] = ( imag(j2PPj1PP) - imag(j2PSj1PS) - imag(j2SPj1SP) + imag(j2SSj1SS))/2.;
        result[2][1] = (-imag(j2PSj1PP) - imag(j2PPj1PS) + imag(j2SSj1SP) + imag(j2SPj1SS))/2.;
        result[3][1] = (-real(j2PSj1PP) + real(j2PPj1PS) + real(j2SSj1SP) - real(j2SPj1SS))/2.;
        result[0][2] = ( imag(j2SPj1PP) + imag(j2SSj1PS) + imag(j2PPj1SP) + imag(j2PSj1SS))/2.;
        result[1][2] = (-imag(j2SPj1PP) + imag(j2SSj1PS) - imag(j2PPj1SP) + imag(j2PSj1SS))/2.;
        result[2][2] = ( imag(j2SSj1PP) + imag(j2SPj1PS) + imag(j2PSj1SP) + imag(j2PPj1SS))/2.;
        result[3][2] = ( real(j2SSj1PP) - real(j2SPj1PS) + real(j2PSj1SP) - real(j2PPj1SS))/2.;
        result[0][3] = (-real(j2SPj1PP) - real(j2SSj1PS) + real(j2PPj1SP) + real(j2PSj1SS))/2.;
        result[1][3] = ( real(j2SPj1PP) - real(j2SSj1PS) - real(j2PPj1SP) + real(j2PSj1SS))/2.;
        result[2][3] = (-real(j2SSj1PP) - real(j2SPj1PS) + real(j2PSj1SP) + real(j2PPj1SS))/2.;
        result[3][3] = ( imag(j2SSj1PP) - imag(j2SPj1PS) - imag(j2PSj1SP) + imag(j2PPj1SS))/2.;

        return result;
    }

    MuellerMatrix::
    MuellerMatrix(const StokesVector& s1,const StokesVector& s2,const StokesVector& s3, const StokesVector& s4)
    {
        for (int i=0; i<4; ++i) {
            m[0][i]=s1[i];
            m[1][i]=s2[i];
            m[2][i]=s3[i];
            m[3][i]=s4[i];
        }
    }

    //
    // Routine to change Mueller matrix to one for a parity conversion...
    //
    MuellerMatrix
    MuellerMatrix::
    parity() const
    {
        MuellerMatrix mm = *this;
        mm[0][2] = -mm[0][2];
        mm[0][3] = -mm[0][3];
        mm[1][2] = -mm[1][2];
        mm[1][3] = -mm[1][3];
        mm[2][0] = -mm[2][0];
        mm[2][1] = -mm[2][1];
        mm[3][0] = -mm[3][0];
        mm[3][1] = -mm[3][1];
        return mm;
    }

    MuellerMatrix&
    MuellerMatrix::
    operator=(const MuellerMatrix& x)
    {
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                m[i][j]=x.m[i][j];
            }
        }
        return *this;
    }

    //
    // Define multiplication between two Mueller matrices
    //
    MuellerMatrix
    MuellerMatrix::
    operator*(const MuellerMatrix& matrix) const
    {
        MuellerMatrix out;
        int i,j,k;
        for (i=0; i<4; ++i) {
            for (j=0; j<4; ++j) {
                out.m[i][j]=0;
                for (k=0; k<4; ++k) {
                    out.m[i][j]+=m[i][k]*matrix.m[k][j];
                }
            }
        }
        return out;
    }

    //
    // Define multiplication between a Mueller matrix and a constant...
    //
    MuellerMatrix
    MuellerMatrix::
    operator*(double d) const
    {
        MuellerMatrix out;

        int i,j;
        for (i=0; i<4; ++i)
            for (j=0; j<4; ++j) {
                out.m[i][j]=m[i][j]*d;
            }
        return out;
    }

    MuellerMatrix
    MuellerMatrix::
    operator/(double d) const
    {
        MuellerMatrix out;

        int i,j;
        for (i=0; i<4; ++i)
            for (j=0; j<4; ++j) {
                out.m[i][j]=m[i][j]/d;
            }
        return out;
    }

    //
    // Multiplication with a StokesVector...
    //
    StokesVector
    MuellerMatrix::
    operator*(const StokesVector& in) const
    {
        StokesVector out;
        int i,j;
        for (i=0; i<4; ++i) {
            out.s[i]=0;
            for (j=0; j<4; ++j)
                out.s[i]+=m[i][j]*in.s[j];
        }
        return out;
    }

    //
    // Define addition of two Mueller matrices...
    //
    MuellerMatrix
    MuellerMatrix::
    operator+(const MuellerMatrix& a) const
    {
        MuellerMatrix c;

        int i,j;
        for (i=0; i<4; ++i) {
            for (j=0; j<4; ++j) {
                c.m[i][j] = m[i][j] + a.m[i][j];
            }
        }
        return c;
    }

    //
    // Define subtraction of two Mueller matrices...
    //
    MuellerMatrix
    MuellerMatrix::
    operator-(const MuellerMatrix& a) const
    {
        MuellerMatrix c;

        int i,j;
        for (i=0; i<4; ++i) {
            for (j=0; j<4; ++j) {
                c.m[i][j] = m[i][j] - a.m[i][j];
            }
        }
        return c;
    }

    MuellerMatrix
    MuellerMatrix::
    operator-() const
    {
        MuellerMatrix c;

        int i,j;
        for (i=0; i<4; ++i) {
            for (j=0; j<4; ++j) {
                c.m[i][j] = -m[i][j];
            }
        }
        return c;
    }

    MuellerMatrix
    MuellerMatrix::
    rotate(const double angle) const
    {
        MuellerMatrix temp(*this);
        return (MuellerMatrix(JonesRotator(angle))*temp)*
               MuellerMatrix(JonesRotator(-angle));
    }


    MuellerMatrix
    MuellerMatrix::
    transpose() const
    {
        MuellerMatrix mm;
        int i,j;
        for (i=0; i<4; ++i)
            for (j=0; j<4; ++j)
                mm[i][j]=m[j][i];
        return mm;
    }


    MuellerMatrix
    MuellerZero()
    {
        MuellerMatrix m;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                m[i][j]=0.;
        return m;
    }

    MuellerMatrix
    MuellerUnit(double transmittance)
    {
        MuellerMatrix m=MuellerZero();
        m[0][0]=m[1][1]=m[2][2]=m[3][3]=transmittance;
        return m;
    }

    MuellerMatrix
    MuellerDepolarizer(double transmittance,double depolarization)
    {
        MuellerMatrix m=MuellerZero();
        m[0][0]=transmittance;
        m[1][1]=m[2][2]=m[3][3]=transmittance*(1-depolarization);
        return m;
    }

    MuellerMatrix
    MuellerDiagonal(double m00,double m11, double m22, double m33)
    {
        MuellerMatrix m=MuellerZero();
		m[0][0]=m00;
		m[1][1]=m11;
		m[2][2]=m22;
		m[3][3]=m33;
        return m;
    }

    MuellerMatrix
    MuellerPartialLinearPolarizer(double tmax, double tmin, double angle)
    {
        MuellerMatrix m;

        double tpt=0.5*(tmax+tmin);
        double tmt=0.5*(tmax-tmin);
        double tt= sqrt(tmax*tmin);
        double cos2a = cos(2*angle);
        double sin2a = sin(2*angle);

        m[0][0]=tpt;
        m[0][1]=tmt*cos2a;
        m[0][2]=tmt*sin2a;
        m[0][3]=0;
        m[1][0]=m[0][1];
        m[1][1]=tpt*sqr(cos2a)+tt*sqr(sin2a); // Correction: TAG 23 JUN 2004
        m[1][2]=(tpt-tt)*cos2a*sin2a;
        m[1][3]=0;
        m[2][0]=m[0][2];
        m[2][1]=m[1][2]; // Correction: TAG 5 DEC 2003
        m[2][2]=tt*sqr(cos2a)+tpt*sqr(sin2a);
        m[2][3]=0;
        m[3][0]=m[3][1]=m[3][2]=0;
        m[3][3]=tt;

        return m;
    }

    ostream& operator<<(ostream& os,const MuellerMatrix& m)
    {
        os << '(';
        for (int i=0; i<4; ++i) {
            os << '(';
            for (int j=0; j<4; ++j) {
                os << m[i][j];
                if (j!=3) os << ',';
            }
            os << ')';
            if (i!=3) os << ',';
        }
        os << ')';
        return os;
    }

    istream& operator>>(istream& is,MuellerMatrix& m)
    {
        StokesVector s1,s2,s3,s4;
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
                                m = MuellerMatrix(s1,s2,s3,s4);
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

    ///
    /// The function valid() uses the condition described in
    /// Clark R. Givens and Alexander B. Kostinski, "A simple necessary and
    /// sufficient condition on physically realizable Mueller matrices,"
    /// Journal of Modern Optics 40, 471-481 (1993).
    ///
    /// Returns true if the matrix is valid and false if it is not.
    ///
    bool
    MuellerMatrix::valid() const
    {
		// Expanding the Mueller matrix slightly helps to ensure that numerical errors don't cause 
		// legitimate Mueller matrices to appear invalid.
		MuellerMatrix expanded = *this + sqrt(numeric_limits<double>::epsilon()) * MuellerDepolarizer(m[0][0], 1);
        CFARRAY Q(4,1);
        CFARRAY W(4,4);
        CFARRAY M(4,4);
        for (int i=1; i<=4; ++i) {
            double Gi = i==1 ? 1 : -1;
            for (int j=1; j<=4; ++j) {
                M(i,j) = 0;
                for (int k=1; k<=4; ++k) {
                    double Gk = k==1 ? 1 : -1;
                    M(i,j) += Gi*expanded.m[k-1][i-1]*Gk*expanded.m[k-1][j-1];
                }
            }
        }
        eigen(M,Q,W,4);
        int imax=0;
        double max = -numeric_limits<double>::max();
		double small = (norm(Q(1)) + norm(Q(2)) + norm(Q(3)) + norm(Q(4)))*numeric_limits<double>::epsilon()*10.;
        for (int i=1; i<=4; ++i) {
            if (fabs(imag(Q(i)))>small) return false;
            if (abs(Q(i))>max) {
                max = abs(Q(i));
                imax = i;
            }
        }
		// Check that the vector with the largest eigenvalue is a valid Stokes vector.  
		// After normalizing to the 1st element, the others must be real and the 
		// sum-of-squares must be less than 1.
		COMPLEX W1 = W(1, imax);
		COMPLEX W2 = W(2, imax) / W1;
		COMPLEX W3 = W(3, imax) / W1;
		COMPLEX W4 = W(4, imax) / W1;
		if (fabs(imag(W2)) > small) return false;
		if (fabs(imag(W3)) > small) return false;
		if (fabs(imag(W4)) > small) return false;
		if (1. < sqr(real(W2)) + sqr(real(W3)) + sqr(real(W4))) return false;
        return true;
    }

    void MuellerToHermitian(const MuellerMatrix& m, CFARRAY& M)
    {
        M.allocate(4,4);

        M(1,1)=(m[0][0]+m[1][1]+m[0][1]+m[1][0])/2.;
        M(1,2)=COMPLEX(m[0][2]+m[1][2],m[0][3]+m[1][3])/2.;
        M(1,3)=COMPLEX(m[2][0]+m[2][1],-m[3][0]-m[3][1])/2.;
        M(1,4)=COMPLEX(m[2][2]+m[3][3],m[2][3]-m[3][2])/2.;
        M(2,2)=(m[0][0]-m[1][1]-m[0][1]+m[1][0])/2.;
        M(2,3)=COMPLEX(m[2][2]-m[3][3],-m[2][3]-m[3][2])/2.;
        M(2,4)=COMPLEX(m[2][0]-m[2][1],-m[3][0]+m[3][1])/2.;
        M(3,3)=(m[0][0]-m[1][1]+m[0][1]-m[1][0])/2.;
        M(3,4)=COMPLEX(m[0][2]-m[1][2],m[0][3]-m[1][3])/2.;
        M(4,4)=(m[0][0]+m[1][1]-m[0][1]-m[1][0])/2.;
        M(2,1)=conj(M(1,2));
        M(3,1)=conj(M(1,3));
        M(4,1)=conj(M(1,4));
        M(3,2)=conj(M(2,3));
        M(4,2)=conj(M(2,4));
        M(4,3)=conj(M(3,4));
    }

    void HermitianToMueller(CFARRAY& M, MuellerMatrix& m)
    {
        m[0][0] = real(M(1,1)+M(2,2)+M(3,3)+M(4,4))/2.;
        m[0][1] = real(M(1,1)-M(2,2)+M(3,3)-M(4,4))/2.;
        m[0][2] = real(M(1,2)+M(2,1)+M(3,4)+M(4,3))/2.;
        m[0][3] = imag(M(1,2)-M(2,1)+M(3,4)-M(4,3))/2.;
        m[1][0] = real(M(1,1)+M(2,2)-M(3,3)-M(4,4))/2.;
        m[1][1] = real(M(1,1)-M(2,2)-M(3,3)+M(4,4))/2.;
        m[1][2] = real(M(1,2)+M(2,1)-M(3,4)-M(4,3))/2.;
        m[1][3] = imag(M(1,2)-M(2,1)-M(3,4)+M(4,3))/2.;
        m[2][0] = real(M(1,3)+M(3,1)+M(2,4)+M(4,2))/2.;
        m[2][1] = real(M(1,3)+M(3,1)-M(2,4)-M(4,2))/2.;
        m[2][2] = real(M(1,4)+M(4,1)+M(2,3)+M(3,2))/2.;
        m[2][3] = imag(M(1,4)-M(4,1)-M(2,3)+M(3,2))/2.;
        m[3][0] = imag(M(3,1)-M(1,3)-M(2,4)+M(4,2))/2.;
        m[3][1] = imag(M(3,1)-M(1,3)+M(2,4)-M(4,2))/2.;
        m[3][2] = imag(M(4,1)-M(1,4)+M(3,2)-M(2,3))/2.;
        m[3][3] = real(M(1,4)+M(4,1)-M(2,3)-M(3,2))/2.;
    }

    ///
    /// The function physically_valid() calculates the Mueller matrix coherency matrix, which is
    /// the coherency matrix for the four Jones-Mueller matrices generated by the Pauli matrices.
    /// This coherency matrix must be a valid coherency matrix (positive semi-definite) if the matrix
    /// is the convex sum of Jones-Mueller matrices.  This test is a more stringent test than valid()
    /// for the validity of the Mueller matrix.
    /// See B.N. Simon, et al, "A complete characterization of pre-Mueller and Mueller matrices
    /// in polarization optics," J. Opt. Soc. Am. A 27(2), 188-199 (2010).
    ///
    bool
    MuellerMatrix::physically_valid() const
    {
        CFARRAY Q(4,1);
        CFARRAY W(4,4);
        CFARRAY M(4,4);

        MuellerToHermitian(*this,M);
        eigen(M,Q,W,4);
		double small = (norm(Q(1)) + norm(Q(2)) + norm(Q(3)) + norm(Q(4))) * numeric_limits<double>::epsilon() * 10.;
        for (int i=1; i<=4; ++i) {
			// Eigenvalues must be real and positive...
            if (real(Q(i)) < -small) return false;
			if (fabs(imag(Q(i))) > small) return false;
        }
        return true;
    }

	///
	/// The function entropy() calculates the polarization entropy defined by Cloude and Pottier in 
	/// S.R. Cloude and E. Pottier, "Concept of polarization entropy in optical scattering," Opt. Eng. 34(6) 1599-1610 (1995)
	///
	double
	MuellerMatrix::entropy() const
	{
		CFARRAY Q(4, 1);
		CFARRAY W(4, 4);
		CFARRAY M(4, 4);

		MuellerToHermitian(*this, M);
		eigen(M, Q, W, 4);

		double sumQ = real(Q(1)) + real(Q(2)) + real(Q(3)) + real(Q(4));
		double result = 0.;
		for (int i = 1; i <= 4; ++i) {
			if (real(Q(i))<0.) return std::numeric_limits<double>::quiet_NaN(); // Invalid Mueller matrix
			double Pi = real(Q(i)) / sumQ;
			if (Pi > 0) {
				result += Pi*std::log(Pi);
			}
		}
		return -result/std::log(4.);
	}

    ///
    /// The Cloude decomposition expresses the Mueller matrix as the sum of four non-depolarizing (Jones-Mueller)
    /// matrices.
    ///
    void MuellerMatrix::Cloude_Decomposition(MuellerMatrix& M1,
            MuellerMatrix& M2,
            MuellerMatrix& M3,
            MuellerMatrix& M4) const
    {
        CFARRAY Q(4,1);
        CFARRAY W(4,4);
        CFARRAY M(4,4);
        MuellerToHermitian(*this,M);
        eigen(M,Q,W,4);
        CFARRAY Wn(4,4);

        for (int i=1; i<=4; ++i) {
            for (int j=1; j<=4; ++j) {
                for (int k=1; k<=4; ++k) {
                    Wn(j,k) = Q(i)*W(j,i)*conj(W(k,i));
                }
            }
            switch (i) {
                case 1:
                    HermitianToMueller(Wn,M1);
                    break;
                case 2:
                    HermitianToMueller(Wn,M2);
                    break;
                case 3:
                    HermitianToMueller(Wn,M3);
                    break;
                case 4:
                    HermitianToMueller(Wn,M4);
                    break;
                default:
                    break;
            }
        }
    }

    MuellerMatrix MuellerMatrix::Closest_NonDepolarizing(int rank) const
    {
        MuellerMatrix M1,M2,M3,M4;
        Cloude_Decomposition(M1,M2,M3,M4);
        MuellerMatrix &m1=M1,&m2=M2,&m3=M3,&m4=M4;
        if (m4[0][0]>m3[0][0]) swap(m3,m4);
        if (m3[0][0]>m2[0][0]) swap(m2,m3);
        if (m2[0][0]>m1[0][0]) swap(m1,m2);
		if (rank==1) return m1;
        if (m4[0][0]>m3[0][0]) swap(m3,m4);
        if (m3[0][0]>m2[0][0]) swap(m2,m3);
		if (rank==2) return m2;
        if (m4[0][0]>m3[0][0]) swap(m3,m4);
		if (rank==3) return m3;
		return m4;
    }

    MuellerMatrix MuellerMatrix::log() const
    {
        CFARRAY Q(4,1);
        CFARRAY W(4,4);
        CFARRAY Wi(4,4);
        CFARRAY M(4,4);
        for (int i=1; i<=4; ++i) {
            for (int j=1; j<=4; ++j) {
                M(i,j)=m[i-1][j-1];
            }
        }

        eigen(M,Q,W,4);
        for (int i=1; i<=16; ++i) Wi(i)=W(i);
        Inverse(Wi,4);

        for (int i=1; i<=4; ++i) Q(i)=std::log(Q(i));
        MuellerMatrix result;
        for (int i=1; i<=4; ++i) {
            for (int j=1; j<=4; ++j) {
                result[i-1][j-1]=0;
                for (int k=1; k<=4; ++k) {
                    result[i-1][j-1] += real(W(i,k)*Q(k)*Wi(k,j));
                }
            }
        }
        return result;
    }

    MuellerMatrix MuellerMatrix::exp() const
    {
        CFARRAY Q(4,1);
        CFARRAY W(4,4);
        CFARRAY Wi(4,4);
        CFARRAY M(4,4);
        for (int i=1; i<=4; ++i) {
            for (int j=1; j<=4; ++j) {
                M(i,j)=m[i-1][j-1];
            }
        }

        eigen(M,Q,W,4);
        for (int i=1; i<=16; ++i) Wi(i)=W(i);
        Inverse(Wi,4);

        for (int i=1; i<=4; ++i) Q(i)=std::exp(Q(i));
        MuellerMatrix result;
        for (int i=1; i<=4; ++i) {
            for (int j=1; j<=4; ++j) {
                result[i-1][j-1]=0;
                for (int k=1; k<=4; ++k) {
                    result[i-1][j-1] += real(W(i,k)*Q(k)*Wi(k,j));
                }
            }
        }
        return result;
    }

    inline int LeviCivita(int i,int j,int k)
    {
        int ijk = i*100+j*10+k;
        switch (ijk) {
            case 123:
            case 231:
            case 312:
            case  12:
            case 120:
            case 201:
                return 1;
            case 132:
            case 213:
            case 321:
            case 102:
            case 210:
            case  21:
                return -1;
            default:
                return 0;
        }
    }

    ///
    /// The Lu-Chipman decomposition decomposes a Mueller matrix into the product of a depolarizer, a retarder,
    /// and a diattenuator, as described in S.-Y. Lu and R.A. Chipman, "Interpretation of Mueller matrices based
    /// upon a polar decomposition," J. Opt. Soc. Am. A, 13(5), 1106-1113 (1996). (Referred to as L&C.)
    ///
    void
    MuellerMatrix::
    Lu_Chipman_Decomposition(MuellerMatrix& depolarizer,MuellerMatrix& retarder, MuellerMatrix& diattenuator) const
    {
        // Use M as a pseudonym for *this...
        const MuellerMatrix &M = *this;

        // The diattenuator takes the net transmittance of the matrix. See L&C, Eq. (18)
        double Tu = M[0][0];

        // If the transmittance is zero, then the decomposition is trivial.
        if (Tu==0.) {
            retarder = MuellerUnit();
            depolarizer = MuellerUnit();
            diattenuator=MuellerZero();
            return;
        }

        // From L&C, Eq. (36), Darrow and Parrow are the diattenuation and polarizance, respectively.
        Vector Darrow = Vector(M[0][1]/Tu,M[0][2]/Tu,M[0][3]/Tu);
        Vector Parrow = Vector(M[1][0]/Tu,M[2][0]/Tu,M[3][0]/Tu);
        double D = sqrt(sqr(Darrow.x)+sqr(Darrow.y)+sqr(Darrow.z));

        // Thus, the diattenuator is given by L&C, Eq. (18),
        if (D>10.*numeric_limits<double>::epsilon()) {
            diattenuator[0][0] = 1.;
            diattenuator[1][0] = diattenuator[0][1] = Darrow.x;
            diattenuator[2][0] = diattenuator[0][2] = Darrow.y;
            diattenuator[3][0] = diattenuator[0][3] = Darrow.z;
            for (int i=1; i<4; ++i) {
                for (int j=1; j<4; ++j) {
                    diattenuator[i][j] = (i==j)? sqrt(1.-sqr(D)) : 0;
                    diattenuator[i][j] += (1.-sqrt(1.-sqr(D)))*Darrow[i-1]*Darrow[j-1]/sqr(D);
                }
            }
        } else {
            diattenuator = MuellerUnit();
        }
        diattenuator = diattenuator*Tu;

        // If the matrix is non-depolarizing, then use
        if (1-M.depolarization_index() < sqrt(numeric_limits<double>::epsilon())) {
            depolarizer=MuellerUnit();
            // If the diattenuation is unity (or very close!)...
            if (1-D < sqrt(numeric_limits<double>::epsilon())) {
                // Get the retardance vector direction...
                Vector Rarrow = cross(Parrow,Darrow);
                double normRarrow = Norm(Rarrow);
                Vector Rhat = (normRarrow!=0) ? unit(Rarrow) : perpto(Parrow,Vector(0,0,1));
                double cosR = Parrow*Darrow/Norm(Parrow)/Norm(Darrow);
                double sinR = (fabs(cosR)<1) ? sqrt(1-sqr(cosR)) : 0.;

                // The retarder is determined from L&C, Eqs. (14) and (15)...
                retarder[0][0]=1.;
                retarder[0][1]=retarder[0][2]=retarder[0][3]=0;
                retarder[1][0]=retarder[2][0]=retarder[3][0]=0;
                for (int i=1; i<4; ++i) {
                    for (int j=1; j<4; ++j) {
                        retarder[i][j] = ((i==j) ? cosR : 0) + Rhat[i-1]*Rhat[j-1]*(1-cosR);
                        for (int k=1; k<=3; ++k) {
                            retarder[i][j] += LeviCivita(i-1,j-1,k-1)*Rhat[k-1]*sinR;
                        }
                    }
                }
            } else {
                // If the matrix does not have unit diattenutation, then
                // just use L&C, Eq. (35)...
                retarder = M*diattenuator.inverse();
            }
            return;
        }

        // L&C, Eq. (47)...
        MuellerMatrix Mprime = M*diattenuator.inverse();

        // L&C, Eq. (48)...
        Matrix mprime;
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                mprime[i][j] = Mprime[i+1][j+1];
            }
        }
        // Also from L&C, Eq. (48)...
        Vector PDelta(Mprime[1][0],Mprime[2][0],Mprime[3][0]);

        // Get the determinant of mprime...
        double mprimedeterminant = mprime.determinant();

        // If mprime is not singular...
        if (mprimedeterminant!=0) {
            // To evaluate L&C, Eq. (52), we need the eigenvalues of
            // m'(m')^T ...
            Matrix mmT = mprime*mprime.transpose();

            // Set up matrices and calculate eigenvalues...
            CFARRAY Q(3,1);
            CFARRAY W(3,3);
            CFARRAY MMT(3,3);
            for (int i=1; i<=3; ++i) {
                for (int j=1; j<=3; ++j) {
                    MMT(i,j) = mmT[i-1][j-1];
                }
            }
            eigen(MMT,Q,W,3);

            // These are the square roots of the eigenvalues
            double lambda1 = sqrt(real(Q(1)));
            double lambda2 = sqrt(real(Q(2)));
            double lambda3 = sqrt(real(Q(3)));

            Matrix ident(1,0,0,0,1,0,0,0,1);
            // Finally, this is L&C, Eq. (52)...
            Matrix mDelta = (mmT+(lambda1*lambda2+lambda2*lambda3+lambda3*lambda1)*ident).inverse()*((lambda1+lambda2+lambda3)*mmT+lambda1*lambda2*lambda3*ident);

            // Use the right sign, using the determinant...
            if (mprimedeterminant<0) {
                mDelta = -mDelta;
            }

            // Thus, the depolarizer is given by L&C, Eq. (48)...
            depolarizer[0][0] = 1;
            depolarizer[0][1] = depolarizer[0][2] = depolarizer[0][3] = 0;
            depolarizer[1][0] = PDelta.x;
            depolarizer[2][0] = PDelta.y;
            depolarizer[3][0] = PDelta.z;
            for (int i=1; i<4; ++i) {
                for (int j=1; j<4; ++j) {
                    depolarizer[i][j] = mDelta[i-1][j-1];
                }
            }

            // And the retarder is determined by L&C, Eq. (53)...
            retarder = depolarizer.inverse()*Mprime;

        } else {
            // If the matrix is singular, then we have to use the
            // formalism described in L&C, Appendix B...

            // First, perform the singular value decomposition
            // of mprime...
            CFARRAY X(3,3);
            CFARRAY S(3,1);
            CFARRAY E(3,1);
            CFARRAY U(3,3);
            CFARRAY V(3,3);
            CFARRAY WORK(3,1);
            int info;
            for (int i=1; i<=3; ++i) {
                for (int j=1; j<=3; ++j) {
                    X(i,j) = mprime[i-1][j-1];
                }
            }
            CMLIB::CSVDC(X,3,3,3,S,E,U,3,V,3,WORK,11,info);

            // The singular values are...
            double lambda1 = real(S(1));
            double lambda2 = real(S(2));
            double lambda3 = real(S(3));

            Matrix mR;
            // How one determines the retardance depends upon how many of the singular
            // values are zero...
            if (fabs(lambda2/lambda1)>10.*numeric_limits<double>::epsilon()) {
                // One of them is zero...
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        // L&C, Eq. (B3)...
                        mR[i][j] = real(V(i+1,1)*U(j+1,1)+V(i+1,2)*U(j+1,2)+V(i+1,3)*U(j+1,3));
                    }
                }
            } else if (lambda1!=0) {
                // Two of them are zero...
                Vector v(real(V(1,1)),real(V(2,1)),real(V(3,1)));
                Vector u(real(U(1,1)),real(U(2,1)),real(U(3,1)));
                Vector pv = perpto(v,Vector(1,0,0));
                Vector ppv = perpto(pv,v);
                Vector pu = perpto(u,Vector(1,0,0));
                Vector ppu = perpto(pu,u);
                // Essentially, L&C, Eq. (B9)...
                mR = outer(v,u)+outer(pv,pu)+outer(ppv,ppu);
            } else {
                // The case of all of them zero is trivial...
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        mR[i][j] = i==j ? 1.:0.;
                    }
                }
            }

            // The full retarder Mueller matrix is...
            for (int i=0; i<4; ++i) {
                for (int j=0; j<4; ++j) {
                    if (i==0&&j==0) retarder[i][j]=1;
                    else if (i==0||j==0) retarder[i][j]=0;
                    else retarder[i][j] = mR[i-1][j-1];
                }
            }

            // The depolarizer is determined from...
            // (The matrix mR is always invertible.)
            Matrix mDelta = mprime*mR.inverse();

            depolarizer[0][0] = 1;
            depolarizer[0][1] = depolarizer[0][2] = depolarizer[0][3] = 0;
            depolarizer[1][0] = PDelta.x;
            depolarizer[2][0] = PDelta.y;
            depolarizer[3][0] = PDelta.z;
            for (int i=1; i<4; ++i) {
                for (int j=1; j<4; ++j) {
                    depolarizer[i][j] = mDelta[i-1][j-1];
                }
            }
        }
		//for (int j=0;j<4;++j) {
		//	if (depolarizer[j][j]<0) {
		//		for (int i=0;i<4;++i) depolarizer[i][j] = -depolarizer[i][j];
		//		for (int i=0;i<4;++i) retarder[j][i] = -retarder[j][i];
		//	}
		//}
    }

    MuellerMatrix
    MuellerMatrix::
    inverse() const
    {
        MuellerMatrix n;
        n[0][0]= -m[1][3]*m[2][2]*m[3][1] + m[1][2]*m[2][3]*m[3][1] + m[1][3]*m[2][1]*m[3][2] - m[1][1]*m[2][3]*m[3][2] -
                 m[1][2]*m[2][1]*m[3][3] + m[1][1]*m[2][2]*m[3][3];
        n[0][1]=  m[0][3]*m[2][2]*m[3][1] - m[0][2]*m[2][3]*m[3][1] - m[0][3]*m[2][1]*m[3][2] + m[0][1]*m[2][3]*m[3][2] +
                  m[0][2]*m[2][1]*m[3][3] - m[0][1]*m[2][2]*m[3][3];
        n[0][2]= -m[0][3]*m[1][2]*m[3][1] + m[0][2]*m[1][3]*m[3][1] + m[0][3]*m[1][1]*m[3][2] - m[0][1]*m[1][3]*m[3][2] -
                 m[0][2]*m[1][1]*m[3][3] + m[0][1]*m[1][2]*m[3][3];
        n[0][3]=  m[0][3]*m[1][2]*m[2][1] - m[0][2]*m[1][3]*m[2][1] - m[0][3]*m[1][1]*m[2][2] + m[0][1]*m[1][3]*m[2][2] +
                  m[0][2]*m[1][1]*m[2][3] - m[0][1]*m[1][2]*m[2][3];
        n[1][0]=  m[1][3]*m[2][2]*m[3][0] - m[1][2]*m[2][3]*m[3][0] - m[1][3]*m[2][0]*m[3][2] + m[1][0]*m[2][3]*m[3][2] +
                  m[1][2]*m[2][0]*m[3][3] - m[1][0]*m[2][2]*m[3][3];
        n[1][1]= -m[0][3]*m[2][2]*m[3][0] + m[0][2]*m[2][3]*m[3][0] + m[0][3]*m[2][0]*m[3][2] - m[0][0]*m[2][3]*m[3][2] -
                 m[0][2]*m[2][0]*m[3][3] + m[0][0]*m[2][2]*m[3][3];
        n[1][2]=  m[0][3]*m[1][2]*m[3][0] - m[0][2]*m[1][3]*m[3][0] - m[0][3]*m[1][0]*m[3][2] + m[0][0]*m[1][3]*m[3][2] +
                  m[0][2]*m[1][0]*m[3][3] - m[0][0]*m[1][2]*m[3][3];
        n[1][3]= -m[0][3]*m[1][2]*m[2][0] + m[0][2]*m[1][3]*m[2][0] + m[0][3]*m[1][0]*m[2][2] - m[0][0]*m[1][3]*m[2][2] -
                 m[0][2]*m[1][0]*m[2][3] + m[0][0]*m[1][2]*m[2][3];
        n[2][0]= -m[1][3]*m[2][1]*m[3][0] + m[1][1]*m[2][3]*m[3][0] + m[1][3]*m[2][0]*m[3][1] - m[1][0]*m[2][3]*m[3][1] -
                 m[1][1]*m[2][0]*m[3][3] + m[1][0]*m[2][1]*m[3][3];
        n[2][1]=  m[0][3]*m[2][1]*m[3][0] - m[0][1]*m[2][3]*m[3][0] - m[0][3]*m[2][0]*m[3][1] + m[0][0]*m[2][3]*m[3][1] +
                  m[0][1]*m[2][0]*m[3][3] - m[0][0]*m[2][1]*m[3][3];
        n[2][2]= -m[0][3]*m[1][1]*m[3][0] + m[0][1]*m[1][3]*m[3][0] + m[0][3]*m[1][0]*m[3][1] - m[0][0]*m[1][3]*m[3][1] -
                 m[0][1]*m[1][0]*m[3][3] + m[0][0]*m[1][1]*m[3][3];
        n[2][3]=  m[0][3]*m[1][1]*m[2][0] - m[0][1]*m[1][3]*m[2][0] - m[0][3]*m[1][0]*m[2][1] + m[0][0]*m[1][3]*m[2][1] +
                  m[0][1]*m[1][0]*m[2][3] - m[0][0]*m[1][1]*m[2][3];
        n[3][0]=  m[1][2]*m[2][1]*m[3][0] - m[1][1]*m[2][2]*m[3][0] - m[1][2]*m[2][0]*m[3][1] + m[1][0]*m[2][2]*m[3][1] +
                  m[1][1]*m[2][0]*m[3][2] - m[1][0]*m[2][1]*m[3][2];
        n[3][1]= -m[0][2]*m[2][1]*m[3][0] + m[0][1]*m[2][2]*m[3][0] + m[0][2]*m[2][0]*m[3][1] - m[0][0]*m[2][2]*m[3][1] -
                 m[0][1]*m[2][0]*m[3][2] + m[0][0]*m[2][1]*m[3][2];
        n[3][2]=  m[0][2]*m[1][1]*m[3][0] - m[0][1]*m[1][2]*m[3][0] - m[0][2]*m[1][0]*m[3][1] + m[0][0]*m[1][2]*m[3][1] +
                  m[0][1]*m[1][0]*m[3][2] - m[0][0]*m[1][1]*m[3][2];
        n[3][3]= -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] -
                 m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2];

        double denom = m[0][3]*m[1][2]*m[2][1]*m[3][0] - m[0][2]*m[1][3]*m[2][1]*m[3][0] - m[0][3]*m[1][1]*m[2][2]*m[3][0] +
                       m[0][1]*m[1][3]*m[2][2]*m[3][0] + m[0][2]*m[1][1]*m[2][3]*m[3][0] - m[0][1]*m[1][2]*m[2][3]*m[3][0] -
                       m[0][3]*m[1][2]*m[2][0]*m[3][1] + m[0][2]*m[1][3]*m[2][0]*m[3][1] + m[0][3]*m[1][0]*m[2][2]*m[3][1] -
                       m[0][0]*m[1][3]*m[2][2]*m[3][1] - m[0][2]*m[1][0]*m[2][3]*m[3][1] + m[0][0]*m[1][2]*m[2][3]*m[3][1] +
                       m[0][3]*m[1][1]*m[2][0]*m[3][2] - m[0][1]*m[1][3]*m[2][0]*m[3][2] - m[0][3]*m[1][0]*m[2][1]*m[3][2] +
                       m[0][0]*m[1][3]*m[2][1]*m[3][2] + m[0][1]*m[1][0]*m[2][3]*m[3][2] - m[0][0]*m[1][1]*m[2][3]*m[3][2] -
                       m[0][2]*m[1][1]*m[2][0]*m[3][3] + m[0][1]*m[1][2]*m[2][0]*m[3][3] + m[0][2]*m[1][0]*m[2][1]*m[3][3] -
                       m[0][0]*m[1][2]*m[2][1]*m[3][3] - m[0][1]*m[1][0]*m[2][2]*m[3][3] + m[0][0]*m[1][1]*m[2][2]*m[3][3];

        return n/denom;
    }


} // namespace SCATMECH

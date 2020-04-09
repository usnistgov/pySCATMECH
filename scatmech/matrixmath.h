//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: matrixmath.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_MATRIXMATH_H
#define SCATMECH_MATRIXMATH_H

#include <vector>
#include <valarray>

#include "scatmech.h"

namespace SCATMECH {

    int eigen(std::vector<COMPLEX>& A,
              std::vector<COMPLEX>& Q,
              std::vector<COMPLEX>& W,int nmat);

    void LUdecompose(std::vector<COMPLEX>& matrix,int nmat,std::vector<int>& pivot);
    void LUbacksubstitute(std::vector<COMPLEX>& matrix,int nmat,std::vector<int>& pivot,std::vector<COMPLEX>& b);
    void LUImprove(std::vector<COMPLEX>& a, std::vector<COMPLEX>& alud, int n, std::vector<int>& indx, std::vector<COMPLEX>& b, std::vector<COMPLEX>& x);

    void Inverse(std::vector<COMPLEX>& matrix, int nmat);

    //
    // The template class FARRAY is a Fortran-like array.
    //
    template <class T>
    class FARRAY
    {
        public:
            FARRAY() {
                p = NULL;
                array(0);
                owner=false;
            }
            FARRAY(T* t) {
                p = t;
                owner=false;
            }
            FARRAY(T& t) {
                p = &t;
                array(0);
                owner=false;
            }
            FARRAY(std::vector<T>& t) {
                p = &t[0];
                array(0);
                owner=false;
            }
            FARRAY(std::valarray<T>& t) {
                p = &t[0];
                array(0);
                owner=false;
            }

            FARRAY(const FARRAY& a) {
                p = (T*)a.p;
                int i;
                for (i=0; i<8; ++i) step[i] = a.step[i];
                #ifdef _DEBUG
                for (i=0; i<8; ++i) dims[i] = a.dims[i];
                #endif
                owner=false;
            }

            const FARRAY& operator=(const FARRAY& a) {
                p = (T*)a.p;
                int i;
                for (i=0; i<8; ++i) step[i] = a.step[i];
                #ifdef _DEBUG
                for (i=0; i<8; ++i) dims[i] = a.dims[i];
                #endif
                owner=false;
                return *this;
            }

            //FARRAY(const FARRAY& a) {p = (T*)(a.p); step1 = a.step1; owner=false;}
            FARRAY(int i,int j) {
                p = NULL;
                allocate(i,j);
            }
            FARRAY(int i,int j,int k) {
                p = NULL;
                allocate(i,j,k);
            }
            FARRAY(int i,int j,int k,int l) {
                p = NULL;
                allocate(i,j,k,l);
            }
            FARRAY(int i,int j,int k,int l,int I) {
                p = NULL;
                allocate(i,j,k,l,I);
            }
            FARRAY(int i,int j,int k,int l,int I,int J) {
                p = NULL;
                allocate(i,j,k,l,I,J);
            }
            FARRAY(int i,int j,int k,int l,int I,int J,int K) {
                p = NULL;
                allocate(i,j,k,l,I,J,K);
            }
            FARRAY(int i,int j,int k,int l,int I,int J,int K,int L) {
                p = NULL;
                allocate(i,j,k,l,I,J,K,L);
            }

            ~FARRAY() {
                deallocate();
            }

            T& operator()(int i)
            {
                #ifdef _DEBUG
                //checkbounds(i);
                #endif
                return p[i-1];
            }
            const T& operator()(int i) const
            {
                #ifdef _DEBUG
                //checkbounds(i);
                #endif
                return p[i-1];
            }
            T& operator()(int i,int j)
            {
                #ifdef _DEBUG
                checkbounds(i,j);
                #endif
                return p[(i-1)+step[0]*(j-1)];
            }
            const T& operator()(int i,int j) const
            {
                #ifdef _DEBUG
                checkbounds(i,j);
                #endif
                return p[(i-1)+step[0]*(j-1)];
            }
            T& operator()(int i,int j,int k)
            {
                #ifdef _DEBUG
                checkbounds(i,j,k);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)];
            }
            const T& operator()(int i,int j,int k) const
            {
                #ifdef _DEBUG
                checkbounds(i,j,k);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)];
            }
            T& operator()(int i,int j,int k,int l)
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)];
            }
            const T& operator()(int i,int j,int k,int l) const
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)];
            }
            T& operator()(int i,int j,int k,int l,int I)
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l,I);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)];
            }
            const T& operator()(int i,int j,int k,int l,int I) const
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l,I);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)];
            }
            T& operator()(int i,int j,int k,int l,int I,int J)
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,k,I,J);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)+step[4]*(J-1)];
            }
            const T& operator()(int i,int j,int k,int l,int I,int J) const
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,k,I,J);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)+step[4]*(J-1)];
            }
            T& operator()(int i,int j,int k,int l,int I,int J,int K)
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l,I,J,K);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)+step[4]*(J-1)+step[5]*(K-1)];
            }
            const T& operator()(int i,int j,int k,int l,int I,int J,int K) const
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l,I,J,K);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)+step[4]*(J-1)+step[5]*(K-1)];
            }

            T& operator()(int i,int j,int k,int l,int I,int J,int K,int L)
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l,I,J,K,L);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)+step[4]*(J-1)+step[5]*(K-1)+step[6]*(L-1)];
            }
            const T& operator()(int i,int j,int k,int l,int I,int J,int K,int L) const
            {
                #ifdef _DEBUG
                checkbounds(i,j,k,l,I,J,K,L);
                #endif
                return p[(i-1)+step[0]*(j-1)+step[1]*(k-1)+step[2]*(l-1)+step[3]*(I-1)+step[4]*(J-1)+step[5]*(K-1)+step[6]*(L-1)];
            }

            T& operator[](int i) {
                return p[i];
            }
            const T& operator[](int i) const {
                return p[i];
            }

            #ifdef _DEBUG
            void checkbounds(int i,int j=1,int k=1,int l=1,int I=1,int J=1,int K=1,int L=1) const {
                if (i>dims[0] || j>dims[1] || k>dims[2] || l>dims[3] || I>dims[4] || J>dims[5] || K>dims[6] || L>dims[7])
                    throw SCATMECH_exception("Array out of bounds in FARRAY");
            }
            #endif

            #ifdef _DEBUG
            void set_dims(int i,int j=1,int k=1,int l=1,int I=1,int J=1,int K=1,int L=1) {
                dims[0]=i;
                dims[1]=j;
                dims[2]=k;
                dims[3]=l;
                dims[4]=I;
                dims[5]=J;
                dims[6]=K;
                dims[7]=L;
            }
            #endif
            void array(int i) {
                step[0] = i;
                #ifdef _DEBUG
                set_dims(i);
                #endif
            }
            void array(int i,int j) {
                step[1] = j*(step[0] = i);
                #ifdef _DEBUG
                set_dims(i,j);
                #endif
            }
            void array(int i,int j,int k) {
                step[2] = k*(step[1] = j*(step[0] = i) );
                #ifdef _DEBUG
                set_dims(i,j,k);
                #endif
            }
            void array(int i,int j,int k,int l) {
                step[3] = l*(step[2] = k*(step[1] = j*(step[0] = i) ) );
                #ifdef _DEBUG
                set_dims(i,j,k,l);
                #endif
            }
            void array(int i,int j,int k,int l,int I) {
                step[4] = I*(step[3] = l*(step[2] = k*(step[1] = j*(step[0] = i) ) ) );
                #ifdef _DEBUG
                set_dims(i,j,k,l,I);
                #endif
            }
            void array(int i,int j,int k,int l,int I,int J) {
                step[5] = J*(step[4] = I*(step[3] = l*(step[2] = k*(step[1] = j*(step[0] = i) ) ) ) );
                #ifdef _DEBUG
                set_dims(i,j,k,l,I,J);
                #endif
            }
            void array(int i,int j,int k,int l,int I,int J,int K) {
                step[6] = K*(step[5] = J*(step[4] = I*(step[3] = l*(step[2] = k*(step[1] = j*(step[0] = i) ) ) ) ) );
                #ifdef _DEBUG
                set_dims(i,j,k,l,I,J,K);
                #endif
            }
            void array(int i,int j,int k,int l,int I,int J,int K,int L) {
                step[7] = L*(step[6] = K*(step[5] = J*(step[4] = I*(step[3] = l*(step[2] = k*(step[1] = j*(step[0] = i) ) ) ) ) ) );
                #ifdef _DEBUG
                set_dims(i,j,k,l,I,J,K,L);
                #endif
            }

            void allocate(int i) {
                deallocate();
                array(i);
                p = new T[i];
                owner=true;
            }
            void allocate(int i,int j) {
                deallocate();
                unsigned size = i*j;
                array(i,j);
                p = new T[size];
                owner=true;
            }
            void allocate(int i,int j,int k) {
                deallocate();
                unsigned size = i*j*k;
                array(i,j,k);
                p = new T[size];
                owner=true;
            }
            void allocate(int i,int j,int k,int l) {
                deallocate();
                unsigned size = i*j*k*l;
                array(i,j,k,l);
                p = new T[size];
                owner=true;
            }
            void allocate(int i,int j,int k,int l,int I) {
                deallocate();
                unsigned size = i*j*k*l*I;
                array(i,j,k,l,I);
                p = new T[size];
                owner=true;
            }
            void allocate(int i,int j,int k,int l,int I,int J) {
                deallocate();
                unsigned size = i*j*k*l*I*J;
                array(i,j,k,l,I,J);
                p = new T[size];
                owner=true;
            }
            void allocate(int i,int j,int k,int l,int I,int J,int K) {
                deallocate();
                unsigned size = i*j*k*l*I*J*K;
                array(i,j,k,l,I,J,K);
                p = new T[size];
                owner=true;
            }
            void allocate(int i,int j,int k,int l,int I,int J,int K,int L) {
                deallocate();
                unsigned size = i*j*k*l*I*J*K*L;
                p = new T[size];
                owner=true;
            }

            void deallocate() {
                if (owner && p!=NULL) delete[] p;
                p=0;
                #ifdef _DEBUG
                for (int i=0; i<8; ++i) dims[i] = 0;
                #endif
            }

            std::string show() {
                std::ostringstream o;
                o << p << ' ';
                for (int i=0; i<8; ++i) o << step[i] << ' ';
                return o.str();
            }

            T& operator=(const T& v) {
                *p = v;
                return *p;
            }

            operator T() const {
                return *p;
            }

        protected:
            bool owner;
            int step[8];
            #ifdef _DEBUG
            int dims[8];
            #endif
            T* p;
    };

    typedef FARRAY<double> DFARRAY;
    typedef FARRAY<COMPLEX> CFARRAY;
    typedef FARRAY<int> IFARRAY;

    int eigen(CFARRAY A,CFARRAY Q, CFARRAY W, int nmat);
    void LUdecompose(CFARRAY matrix,int nmat,IFARRAY pivot);
    void LUbacksubstitute(CFARRAY matrix,int nmat,IFARRAY pivot,CFARRAY b);
    void LUdecompose(DFARRAY matrix,int nmat,IFARRAY pivot);
    void LUbacksubstitute(DFARRAY matrix,int nmat,IFARRAY pivot,DFARRAY b);
    void LUImprove(CFARRAY a, CFARRAY alud, int n, IFARRAY indx, CFARRAY b, CFARRAY x);
    void Inverse(CFARRAY matrix, int nmat);
    void Inverse(DFARRAY matrix, int nmat);

    namespace CMLIB {

        void CGEEV(DFARRAY A,int LDA, int N, DFARRAY E0,DFARRAY V, int LDV, DFARRAY WORK, int JOB, int& INFO);
        void CGEFA(CFARRAY A, int LDA, int N, IFARRAY IPVT, int& INFO);
        void CGESL(CFARRAY A, int LDA, int N, IFARRAY IPVT, CFARRAY B, int JOB);
        void CGEDI(CFARRAY A, int LDA, int N, IFARRAY IPVT, CFARRAY DET, CFARRAY WORK, int JOB);

        void SCOPY(int N,DFARRAY SX,int INCX,DFARRAY SY,int INCY);
        void ZCOPY(int N,CFARRAY ZX,int INCX,CFARRAY ZY,int INCY);
        void CBABK2(int& NM,int& N,int& LOW,int& IGH,DFARRAY SCALE,int& M,DFARRAY ZR,DFARRAY ZI);
        void CBAL(int& NM,int& N,DFARRAY AR,DFARRAY AI,int& LOW,int& IGH,DFARRAY SCALE);
        void COMQR(int& NM,int& N,int& LOW,int& IGH,DFARRAY HR,DFARRAY HI,DFARRAY WR,DFARRAY WI,int& IERR);
        void COMQR2(int& NM,int& N,int& LOW,int& IGH,DFARRAY ORTR,DFARRAY ORTI,
                    DFARRAY HR,DFARRAY HI,DFARRAY WR,DFARRAY WI,DFARRAY ZR,DFARRAY ZI,int& IERR);
        void CORTH(int& NM,int& N,int& LOW,int& IGH,DFARRAY AR,DFARRAY AI,DFARRAY ORTR,DFARRAY ORTI);
        void CSROOT(const double& XR, const double& XI, DFARRAY YR,DFARRAY YI);
        double PYTHAG(const double& A,const double& B);
        void CDIV(const double& AR,const double& AI,const double& BR,const double& BI,DFARRAY CR,DFARRAY CI);
        void CAXPY(int N, COMPLEX CA, FARRAY<COMPLEX> CX, int INCX, FARRAY<COMPLEX> CY, int INCY);
        void CSCAL(int N, COMPLEX CA, FARRAY<COMPLEX> CX, int INCX);
        int ICAMAX(int N, FARRAY<COMPLEX> CX, int INCX);
        COMPLEX CDOTC(int N, FARRAY<COMPLEX> CX, int INCX, FARRAY<COMPLEX> CY, int INCY);
        void CSWAP(int N,CFARRAY CX, int INCX, CFARRAY CY, int INCY);
        void ZGEMV(const char *TRANS,int M,int N,COMPLEX ALPHA,CFARRAY A,int LDA,CFARRAY X,int INCX,COMPLEX BETA,CFARRAY Y,int INCY);
        void ZGEMM(const char* TRANSA,const char* TRANSB,int M,int N,int K,COMPLEX ALPHA,CFARRAY A,int LDA,CFARRAY B,int LDB,COMPLEX BETA,CFARRAY C,int LDC);

        void DGEFA(FARRAY<double> A,int LDA,int N,FARRAY<int> IPVT,int& INFO);
        void DAXPY(int N,double DA,FARRAY<double> DX,int INCX,FARRAY<double> DY,int INCY);
        void DSCAL(int N,double DA,FARRAY<double> DX,int INCX);
        int IDAMAX(int N,FARRAY<double> DX,int INCX);

        void DGESL(DFARRAY A,int LDA,int N,IFARRAY IPVT,DFARRAY B,int JOB);
        double DDOT(int N,DFARRAY DX,int INCX,DFARRAY DY,int INCY);

        void DSWAP(int N,FARRAY<double> DX,int INCX,FARRAY<double> DY,int INCY);
        void DGEDI(FARRAY<double> A,int LDA,int N,FARRAY<int> IPVT,FARRAY<double> DET,FARRAY<double> WORK,int JOB);

        /// *** Added for SVD...

        void SROTG(double& SA,double& SB,double& SC,double& SS);
        double SCNRM2(const int& N,CFARRAY CX,const int& INCX);
        void CSROT(const int& N,CFARRAY CX,const int& INCX,CFARRAY CY,const int& INCY,const double& C,const double& S);
        void CSVDC(CFARRAY X,const int& LDX,const int& N,const int& P,CFARRAY S,CFARRAY E,CFARRAY U,const int& LDU,CFARRAY V,const int& LDV,CFARRAY WORK,const int& JOB,int& INFO);

    }

}

//
// Use #define USING_BLAS to use optimized BLAS routines
//
#ifdef USING_BLAS
extern "C" {
    struct cmplx {
        double re;
        double im;
    };
    extern void dcopy_(int *N,double *SX,int *INCX,double *SY,int *INCY);
    extern void zaxpy_(int *N, cmplx *CA, cmplx *CX, int *INCX, cmplx *CY, int *INCY);
    extern void zscal_(int *N, cmplx *CA, cmplx *CX, int *INCX);
    extern int izamax_(int *N, cmplx *CX, int *INCX);
    extern cmplx zdotc_(int *N, cmplx *CX, int *INCX, cmplx *CY, int *INCY);
    extern void zswap_(int *N, cmplx *CX, int *INCX, cmplx *CY, int *INCY);
    extern void zgemv_(const char *TRANS,int* M, int* N, cmplx *ALPHA, cmplx *A, int* LDA, cmplx *X, int* INCX, cmplx* BETA, cmplx* Y, int* INCY);
    extern void zgemm_(const char *TRANSA,const char*  TRANSB, int* M, int* N, int* K, cmplx* ALPHA, cmplx* A, int* LDA, cmplx* B, int* LDB, cmplx* BEAT, cmplx* C, int *LDC);

    extern void dgefa_(double* A,int* LDA,int* N,int* IPVT,int* INFO);
    extern void daxpy_(int* N,double* DA,double* DX,int* INCX,double* DY,int* INCY);
    extern void dscal_(int* N,double* DA,double* DX,int* INCX);
    extern int idamax_(int* N,double* DX,int* INCX);
    extern void dswap_(int* N,double* DX,int* INCX,double* DY,int* INCY);
    extern void dgedi_(double* A,int* LDA,int* N,int* IPVT,double* DET,double* WORK,int* JOB);
    extern void dgesl_(double& A,int* LDA,int* N,int* IPVT,double* B,int* JOB);
    extern double ddot_(int* N,double* DX,int* INCX,double* DY,int* INCY);
}
#endif


#endif

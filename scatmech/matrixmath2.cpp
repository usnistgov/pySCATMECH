//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: matrixmath2.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include <vector>
#include <cmath>
#include <limits>
#include "scatmech.h"

#include "matrixmath.h"

using namespace std;

#define THROW(x) throw SCATMECH_exception(x);

namespace SCATMECH {

    namespace CMLIB {

        //template <class T> inline T mIn(const T a,const T b) {
        //    return (a>b)? b : a;
        //}
        //template <class T> inline T mAx(const T a,const T b) {
        //    return (a<b)? b : a;
        //}


        inline double CABS1(COMPLEX ZDUM) {
            return abs(real(ZDUM)) + abs(imag(ZDUM));
        }
        inline COMPLEX CSIGN(const COMPLEX& ZDUM1,const COMPLEX& ZDUM2) {
            return abs(ZDUM1)*(ZDUM2/abs(ZDUM2));
        }
        inline int MIN0(int a,int b) {
            return a>b? b : a;
        }
        inline int MAX0(int a,int b) {
            return a>b? a : b;
        }
        template<class T> T max(const T& a,const T& b) {
            return a>b ? a : b;
        }
        inline double AMAX1(double a,double b,double c,double d,double e) {
            return max(a,max(b,max(c,max(d,e))));
        }


        void CSVDC(CFARRAY X,const int& LDX,const int& N,const int& P,CFARRAY S,CFARRAY E,CFARRAY U,const int& LDU,CFARRAY V,const int& LDV,CFARRAY WORK,const int& JOB,int& INFO)
        {

            X.array(LDX,P);
            if (JOB/10==1) U.array(LDU,N);
            else if (JOB/10==2) U.array(LDU,N<P?N:P);
            else U.array(LDU);
            V.array(LDV,P);

            //    INTEGER LDX,N,P,LDU,LDV,JOB,INFO
            //    COMPLEX X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
            //C***BEGIN PROLOGUE  CSVDC
            //C***DATE WRITTEN   790319   (YYMMDD)
            //C***REVISION DATE  820801   (YYMMDD)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C***CATEGORY NO.  D6
            //C***KEYWORDS  COMPLEX,LINEAR ALGEBRA,LINPACK,MATRIX,
            //C             SINGULAR VALUE DECOMPOSITION
            //C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
            //C***PURPOSE  Perform the singular value decomposition of a COMPLEX NXP
            //C            matrix.
            //C***DESCRIPTION
            //C
            //C     CSVDC is a subroutine to reduce a complex NxP matrix X by
            //C     unitary transformations U and V to diagonal form.  The
            //C     diagonal elements S(I) are the singular values of X.  The
            //C     columns of U are the corresponding left singular vectors,
            //C     and the columns of V the right singular vectors.
            //C
            //C     On Entry
            //C
            //C         X         COMPLEX(LDX,P), where LDX .GE. N.
            //C                   X contains the matrix whose singular value
            //C                   decomposition is to be computed.  X is
            //C                   destroyed by CSVDC.
            //C
            //C         LDX       INTEGER.
            //C                   LDX is the leading dimension of the array X.
            //C
            //C         N         INTEGER.
            //C                   N is the number of rows of the matrix X.
            //C
            //C         P         INTEGER.
            //C                   P is the number of columns of the matrix X.
            //C
            //C         LDU       INTEGER.
            //C                   LDU is the leading dimension of the array U
            //C                   (see below).
            //C
            //C         LDV       INTEGER.
            //C                   LDV is the leading dimension of the array V
            //C                   (see below).
            //C
            //C         WORK      COMPLEX(N).
            //C                   WORK is a scratch array.
            //C
            //C         JOB       INTEGER.
            //C                   JOB controls the computation of the singular
            //C                   vectors.  It has the decimal expansion AB
            //C                   with the following meaning
            //C
            //C                        A .EQ. 0    Do not compute the left singular
            //C                                    vectors.
            //C                        A .EQ. 1    Return the N left singular vectors
            //C                                    in U.
            //C                        A .GE. 2    Return the first MIN(N,P)
            //C                                    left singular vectors in U.
            //C                        B .EQ. 0    Do not compute the right singular
            //C                                    vectors.
            //C                        B .EQ. 1    Return the right singular vectors
            //C                                    in V.
            //C
            //C     On Return
            //C
            //C         S         COMPLEX(MM), where MM = MIN(N+1,P).
            //C                   The first MIN(N,P) entries of S contain the
            //C                   singular values of X arranged in descending
            //C                   order of magnitude.
            //C
            //C         E         COMPLEX(P).
            //C                   E ordinarily contains zeros.  However see the
            //C                   discussion of INFO for exceptions.
            //C
            //C         U         COMPLEX(LDU,K), where LDU .GE. N.  If JOBA .EQ. 1
            //C                                   then K .EQ. N.  If JOBA .GE. 2 then
            //C                                   K .EQ. MIN(N,P).
            //C                   U contains the matrix of right singular vectors.
            //C                   U is not referenced if JOBA .EQ. 0.  If N .LE. P
            //C                   or if JOBA .GT. 2, then U may be identified with X
            //C                   in the subroutine call.
            //C
            //C         V         COMPLEX(LDV,P), where LDV .GE. P.
            //C                   V contains the matrix of right singular vectors.
            //C                   V is not referenced if JOB .EQ. 0.  If P .LE. N,
            //C                   then V may be identified with X in the
            //C                   subroutine call.
            //C
            //C         INFO      INTEGER.
            //C                   The singular values (and their corresponding
            //C                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
            //C                   are correct (here M=MIN(N,P)).  Thus if
            //C                   INFO.EQ. 0, all the singular values and their
            //C                   vectors are correct.  In any event, the matrix
            //C                   B = CTRANS(U)*X*V is the bidiagonal matrix
            //C                   with the elements of S on its diagonal and the
            //C                   elements of E on its super-diagonal (CTRANS(U)
            //C                   is the conjugate-transpose of U).  Thus the
            //C                   singular values of X and B are the same.
            //C
            //C     LINPACK.  This version dated 03/19/79 .
            //C     Stewart, G. W., University of Maryland, Argonne National Lab.
            //C
            //C     CSVDC uses the following functions and subprograms.
            //C
            //C     External CSROT
            //C     BLAS CAXPY,CDOTC,CSCAL,CSWAP,SCNRM2,SROTG
            //C     Fortran ABS,AIMAG,AMAX1,CABS,CMPLX
            //C     Fortran CONJG,MAX0,MIN0,MOD,REAL,SQRT
            //C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
            //C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
            //C***ROUTINES CALLED  CAXPY,CDOTC,CSCAL,CSROT,CSWAP,SCNRM2,SROTG
            //C***END PROLOGUE  CSVDC
            //C
            //C
            using namespace SCATMECH::CMLIB;
            int I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,
                MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1;
            COMPLEX T,R;
            double B,C,CS,EL,EMM1,F,G,SCALE,SHIFT,SL,SM,SN,SMM1,T1,TEST,ZTEST;
            bool WANTU,WANTV;
            COMPLEX ZDUM,ZDUM1,ZDUM2;
            //C
            //C     SET THE MAXIMUM NUMBER OF ITERATIONS.
            //C
	    int MAXIT = 30;
	    
            //C***FIRST EXECUTABLE STATEMENT  CSVDC
            //C
            //C     DETERMINE WHAT IS TO BE COMPUTED.
            //C
            WANTU = false;
            WANTV = false;
            JOBU = (JOB%100)/10;
            NCU = N;
            if (JOBU > 1) NCU = MIN0(N,P);
            if (JOBU != 0) WANTU = true;
            if ((JOB%10) != 0) WANTV = true;
            //C
            //C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
            //C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
            //C
            INFO = 0;
            NCT = MIN0(N-1,P);
            NRT = MAX0(0,MIN0(P-2,N));
            LU = MAX0(NCT,NRT);
            if (LU >= 1) {
                for (L=1; L<=LU; ++L) {
                    LP1 = L + 1;
                    if (L <= NCT) {
                        //C
                        //C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
                        //C           PLACE THE L-TH DIAGONAL IN S(L).
                        //C
                        S(L) = COMPLEX(SCNRM2(N-L+1,X(L,L),1),0.0E0);
                        if (CABS1(S(L)) != 0.0E0) {
                            if (CABS1(X(L,L)) != 0.0E0) S(L) = CSIGN(S(L),X(L,L));
                            CSCAL(N-L+1,1.0E0/S(L),X(L,L),1);
                            X(L,L) = COMPLEX(1.0E0,0.0E0) + X(L,L);
                        }
                        S(L) = -S(L);
                    }
                    if (P >= LP1) {
                        for (J = LP1; J<=P; ++J) {
                            if (L <= NCT) {
                                if (CABS1(S(L)) != 0.0E0) {
                                    //C
                                    //C              APPLY THE TRANSFORMATION.
                                    //C
                                    T = -CDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L);
                                    CAXPY(N-L+1,T,X(L,L),1,X(L,J),1);
                                }
                            }
                            //C           PLACE THE L-TH ROW OF X INTO  E FOR THE
                            //C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
                            //C
                            E(J) = conj(X(L,J));
                        }
                    }
                    if (!(!WANTU || L > NCT)) {
                        //C
                        //C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
                        //C           MULTIPLICATION.
                        //C
                        for (I=L; I<=N; ++I) {
                            U(I,L) = X(I,L);
                        }
                    }
                    if (L <= NRT) {
                        //C
                        //C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
                        //C           L-TH SUPER-DIAGONAL IN E(L).
                        //C
                        E(L) = COMPLEX(SCNRM2(P-L,E(LP1),1),0.0E0);
                        if (CABS1(E(L)) != 0.0E0) {
                            if (CABS1(E(LP1)) != 0.0E0) E(L) = CSIGN(E(L),E(LP1));
                            CSCAL(P-L,1.0E0/E(L),E(LP1),1);
                            E(LP1) = COMPLEX(1.0E0,0.0E0) + E(LP1);
                        }
                        E(L) = -conj(E(L));
                        if (!(LP1 > N || CABS1(E(L)) == 0.0E0)) {
                            //C
                            //C              APPLY THE TRANSFORMATION.
                            //C
                            for (I=LP1; I<=N; ++I) {
                                WORK(I) = COMPLEX(0.0E0,0.0E0);
                            }
                            for (J=LP1; J<=P; ++J) {
                                CAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1);
                            }
                            for (J=LP1; J<=P; ++J) {
                                COMPLEX temp = conj(-E(J)/E(LP1));
                                CAXPY(N-L,temp,WORK(LP1),1,X(LP1,J),1);
                            }
                        }
                        if (WANTV) {
                            //C
                            //C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
                            //C              BACK MULTIPLICATION.
                            //C
                            for (I=LP1; I<=P; ++I) {
                                V(I,L) = E(I);
                            }
                        }
                    }
                }
            }
            //C
            //C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
            //C
            M = MIN0(P,N+1);
            NCTP1 = NCT + 1;
            NRTP1 = NRT + 1;
            if (NCT < P) S(NCTP1) = X(NCTP1,NCTP1);
            if (N < M) S(M) = COMPLEX(0.0E0,0.0E0);
            if (NRTP1 < M) E(NRTP1) = X(NRTP1,M);
            E(M) = COMPLEX(0.0E0,0.0E0);
            //C
            //C     IF REQUIRED, GENERATE U.
            //C
            if (WANTU) {
                if (NCU >= NCTP1) {
                    for (J=NCTP1; J<=NCU; ++J) {
                        for (I=1; I<=N; ++I) {
                            U(I,J) = COMPLEX(0.0E0,0.0E0);
                        }
                        U(J,J) = COMPLEX(1.0E0,0.0E0);
                    }
                }
                if (NCT >= 1) {
                    for (LL = 1; LL<=NCT; ++LL) {
                        L = NCT - LL + 1;
                        if (CABS1(S(L)) != 0.0E0) {
                            LP1 = L + 1;
                            if (NCU >= LP1) {
                                for (J=LP1; J<=NCU; ++J) {
                                    T = -CDOTC(N-L+1,U(L,L),1,U(L,J),1)/U(L,L);
                                    CAXPY(N-L+1,T,U(L,L),1,U(L,J),1);
                                }
                            }
                            CSCAL(N-L+1,COMPLEX(-1.0E0,0.0E0),U(L,L),1);
                            U(L,L) = COMPLEX(1.0E0,0.0E0) + U(L,L);
                            LM1 = L - 1;
                            if (LM1 >= 1) {
                                for (I = 1; I<=LM1; ++I) {
                                    U(I,L) = COMPLEX(0.0E0,0.0E0);
                                }
                            }
                        } else {
                            for (I=1; I<=N; ++I) {
                                U(I,L) = COMPLEX(0.0E0,0.0E0);
                            }
                            U(L,L) = COMPLEX(1.0E0,0.0E0);
                        }
                    }
                }
            }
            //C
            //C     IF IT IS REQUIRED, GENERATE V.
            //C
            if (WANTV) {
                for (LL=1; LL<=P; ++LL) {
                    L = P - LL + 1;
                    LP1 = L + 1;
                    if (L <= NRT) {
                        if (CABS1(E(L)) != 0.0E0) {
                            for (J = LP1; J<=P; ++J) {
                                T = -CDOTC(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L);
                                CAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1);
                            }
                        }
                    }
                    for (I=1; I<=P; ++I) {
                        V(I,L) = COMPLEX(0.0E0,0.0E0);
                    }
                    V(L,L) = COMPLEX(1.0E0,0.0E0);
                }
            }
            //C
            //C     TRANSFORM S AND E SO THAT THEY ARE REAL.
            //C
            for (I=1; I<=M; ++I) {
                if (CABS1(S(I)) != 0.0E0) {
                    T = COMPLEX(abs(S(I)),0.0E0);
                    R = S(I)/T;
                    S(I) = T;
                    if (I < M) E(I) = E(I)/R;
                    if (WANTU) CSCAL(N,R,U(1,I),1);
                }
                //C     ...EXIT
                if (I == M) goto Line390;
                if (CABS1(E(I)) != 0.0E0) {
                    T = COMPLEX(abs(E(I)),0.0E0);
                    R = T/E(I);
                    E(I) = T;
                    S(I+1) = S(I+1)*R;
                    if (WANTV) CSCAL(P,R,V(1,I+1),1);
                }
            }

Line390:
            ;
            //C
            //C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
            //C
            MM = M;
            ITER = 0;
            //C
            //C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
            //C
            while (M != 0) {
                //C
                //C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
                //C        FLAG AND RETURN.
                //C
                if (ITER >= MAXIT) {
                    INFO = M;
                    //C     ......EXIT
                    return;
                }
                //C
                //C        THIS SECTION OF THE PROGRAM INSPECTS FOR
                //C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
                //C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
                //C
                //C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
                //C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
                //C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
                //C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
                //C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
                //C
                for (LL=1; LL<=M; ++LL) {
                    L = M - LL;
                    //C        ...EXIT
                    if (L == 0) goto Line440;
                    TEST = abs(S(L)) + abs(S(L+1));
                    ZTEST = TEST + abs(E(L));
                    if (ZTEST == TEST) {
                        E(L) = COMPLEX(0.0E0,0.0E0);
                        //C        ......EXIT
                        goto Line440;
                    }
                }
Line440:
                if (L == M - 1) {
                    KASE = 4;
                } else {
                    LP1 = L + 1;
                    MP1 = M + 1;
                    for (LLS = LP1; LLS<=MP1; ++LLS) {
                        LS = M - LLS + LP1;
                        //C           ...EXIT
                        if (LS == L) goto Line480;
                        TEST = 0.0E0;
                        if (LS != M) TEST = TEST + abs(E(LS));
                        if (LS != L + 1) TEST = TEST + abs(E(LS-1));
                        ZTEST = TEST + abs(S(LS));
                        if (ZTEST == TEST) {
                            S(LS) = COMPLEX(0.0E0,0.0E0);
                            //C           ......EXIT
                            goto Line480;
                        }
                    }
Line480:
                    if (LS == L) {
                        KASE = 3;
                    } else if (LS == M) {
                        KASE = 1;
                    } else {
                        KASE = 2;
                        L = LS;
                    }
                }
                L = L + 1;
                //C
                //C        PERFORM THE TASK INDICATED BY KASE.
                //C
                switch (KASE) {
                    //C
                    //C        DEFLATE NEGLIGIBLE S(M).
                    //C
                    case 1:
                        MM1 = M - 1;
                        F = real(E(M-1));
                        E(M-1) = COMPLEX(0.0E0,0.0E0);
                        for (KK=L; KK<=MM1; ++KK) {
                            K = MM1 - KK + L;
                            T1 = real(S(K));
                            SROTG(T1,F,CS,SN);
                            S(K) = COMPLEX(T1,0.0E0);
                            if (K != L) {
                                F = -SN*real(E(K-1));
                                E(K-1) = CS*E(K-1);
                            }
                            if (WANTV) CSROT(P,V(1,K),1,V(1,M),1,CS,SN);
                        }
                        break;
                    //C
                    //C        SPLIT AT NEGLIGIBLE S(L).
                    //C
                    case 2:
                        F = real(E(L-1));
                        E(L-1) = COMPLEX(0.0E0,0.0E0);
                        for (K=L; K<=M; ++K) {
                            T1 = real(S(K));
                            SROTG(T1,F,CS,SN);
                            S(K) = COMPLEX(T1,0.0E0);
                            F = -SN*real(E(K));
                            E(K) = CS*E(K);
                            if (WANTU) CSROT(N,U(1,K),1,U(1,L-1),1,CS,SN);
                        }
                        break;
                    //C
                    //C        PERFORM ONE QR STEP.
                    //C
                    case 3:
                        //C
                        //C           CALCULATE THE SHIFT.
                        //C
                        SCALE = AMAX1(abs(S(M)),abs(S(M-1)),abs(E(M-1)),abs(S(L)),abs(E(L)));
                        SM = real(S(M))/SCALE;
                        SMM1 = real(S(M-1))/SCALE;
                        EMM1 = real(E(M-1))/SCALE;
                        SL = real(S(L))/SCALE;
                        EL = real(E(L))/SCALE;
                        B = ((SMM1 + SM)*(SMM1 - SM) + sqr(EMM1))/2.0E0;
                        C = sqr(SM*EMM1);
                        SHIFT = 0.0E0;
                        if (!(B == 0.0E0 && C == 0.0E0)) {
                            SHIFT = sqrt(sqr(B)+C);
                            if (B < 0.0E0) SHIFT = -SHIFT;
                            SHIFT = C/(B + SHIFT);
                        }
                        F = (SL + SM)*(SL - SM) - SHIFT;
                        G = SL*EL;
                        //C
                        //C           CHASE ZEROS.
                        //C
                        MM1 = M - 1;
                        for (K=L; K<=MM1; ++K) {
                            SROTG(F,G,CS,SN);
                            if (K != L) E(K-1) = COMPLEX(F,0.0E0);
                            F = CS*real(S(K)) + SN*real(E(K));
                            E(K) = CS*E(K) - SN*S(K);
                            G = SN*real(S(K+1));
                            S(K+1) = CS*S(K+1);
                            if (WANTV) CSROT(P,V(1,K),1,V(1,K+1),1,CS,SN);
                            SROTG(F,G,CS,SN);
                            S(K) = COMPLEX(F,0.0E0);
                            F = CS*real(E(K)) + SN*real(S(K+1));
                            S(K+1) = -SN*E(K) + CS*S(K+1);
                            G = SN*real(E(K+1));
                            E(K+1) = CS*E(K+1);
                            if (WANTU && K < N) CSROT(N,U(1,K),1,U(1,K+1),1,CS,SN);
                        }
                        E(M-1) = COMPLEX(F,0.0E0);
                        ITER = ITER + 1;
                        break;
                    //C
                    //C        CONVERGENCE.
                    //C
                    case 4:
                        //C
                        //C           MAKE THE SINGULAR VALUE  POSITIVE
                        //C
                        if (real(S(L)) < 0.0E0) {
                            S(L) = -S(L);
                            if (WANTV) CSCAL(P,COMPLEX(-1.0E0,0.0E0),V(1,L),1);
                        }
                        //C
                        //C           ORDER THE SINGULAR VALUE.
                        //C

                        while (L != MM) {
                            //C           ...EXIT
                            if (real(S(L)) >= real(S(L+1))) goto Line640;
                            T = S(L);
                            S(L) = S(L+1);
                            S(L+1) = T;
                            if (WANTV && L < P) CSWAP(P,V(1,L),1,V(1,L+1),1);
                            if (WANTU && L < N) CSWAP(N,U(1,L),1,U(1,L+1),1);
                            L = L + 1;
                        }
Line640:
                        ITER = 0;
                        M = M - 1;
                        break;
                }
            }
            return;
        }

        void CSROT(const int& N,CFARRAY CX,const int& INCX,CFARRAY CY,const int& INCY,const double& C,const double& S)
        {
            COMPLEX CTEMP;
            int I,IX,IY;
            //C***BEGIN PROLOGUE  CSROT
            //C***DATE WRITTEN   810223   (YYMMDD)
            //C***REVISION DATE  820801   (YYMMDD)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C***CATEGORY NO.  D1B10
            //C***KEYWORDS  BLAS,COMPLEX,LINEAR ALGEBRA,PLANE ROTATION,VECTOR
            //C***AUTHOR  DONGARRA, J., (ANL)
            //C***PURPOSE  Applies a plane rotation to complex vectors.
            //C***DESCRIPTION
            //C
            //C     CSROT applies the complex Givens rotation
            //C
            //C          (X)   ( C S)(X)
            //C          (Y) = (-S C)(Y)
            //C
            //C     N times where for I = 0,...,N-1
            //C
            //C          X = CX(1+I*INCX)
            //C          Y = CY(1+I*INCY)
            //C
            //C     Argument Description
            //C
            //C        N      (integer)  number of elements in each vector
            //C
            //C        CX     (complex array)  beginning of one vector
            //C
            //C        INCX   (integer)  memory spacing of successive elements
            //C               of vector CX
            //C
            //C        CY     (complex array)  beginning of the other vector
            //C
            //C        INCY   (integer)  memory spacing of successive elements
            //C               of vector CY
            //C
            //C        C      (real)  cosine term of the rotation
            //C
            //C        S      (real)  sine term of the rotation.
            //C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
            //C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
            //C***ROUTINES CALLED  (NONE)
            //C***END PROLOGUE  CSROT
            //C
            //C***FIRST EXECUTABLE STATEMENT  CSROT
            if (N<=0) return;
            if (INCX==1&&INCY==1) goto Line20;
            //C
            //C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
            //C         TO 1
            //C
            IX = 1;
            IY = 1;
            if (INCX<0) IX = (-N+1)*INCX + 1;
            if (INCY<0) IY = (-N+1)*INCY + 1;

            for (I=1; I<=N; ++I) {
                CTEMP = C*CX(IX) + S*CY(IY);
                CY(IY) = C*CY(IY) - S*CX(IX);
                CX(IX) = CTEMP;
                IX = IX + INCX;
                IY = IY + INCY;
            }
            return;
            //C
            //C       CODE FOR BOTH INCREMENTS EQUAL TO 1
            //C
Line20:
            for (I=1; I<=N; ++I) {
                CTEMP = C*CX(I) + S*CY(I);
                CY(I) = C*CY(I) - S*CX(I);
                CX(I) = CTEMP;
            }
            return;
        }

        double SCNRM2(const int& N,CFARRAY CX,const int& INCX)
        {
            //C***BEGIN PROLOGUE  SCNRM2
            //C***DATE WRITTEN   791001   (YYMMDD)
            //C***REVISION DATE  820801   (YYMMDD)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***CATEGORY NO.  D1A3B
            //C***KEYWORDS  BLAS,COMPLEX,LINEAR ALGEBRA,NORM,UNITARY,VECTOR
            //C***AUTHOR  LAWSON, C. L., (JPL)
            //C           HANSON, R. J., (SNLA)
            //C           KINCAID, D. R., (U. OF TEXAS)
            //C           KROGH, F. T., (JPL)
            //C***PURPOSE  Unitary norm of complex vector
            //C***DESCRIPTION
            //C
            //C                B L A S  Subprogram
            //C    Description of Parameters
            //C
            //C     --Input--
            //C        N  number of elements in input vector(s)
            //C       CX  complex vector with N elements
            //C     INCX  storage spacing between elements of CX
            //C
            //C     --Output--
            //C   SCNRM2  single precision result (zero if N .LE. 0)
            //C
            //C     unitary norm of the complex N-vector stored in CX() with storage
            //C     increment INCX .
            //C     If N .LE. 0, return with result = 0.
            //C     If N .GE. 1, then INCX must be .GE. 1
            //C
            //C           C. L. Lawson, 1978 Jan 08
            //C
            //C     Four phase method     using two built-in constants that are
            //C     hopefully applicable to all machines.
            //C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
            //C         CUTHI = minimum of  SQRT(V)      over all known machines.
            //C     where
            //C         EPS = smallest no. such that EPS + 1. .GT. 1.
            //C         U   = smallest positive no.   (underflow limit)
            //C         V   = largest  no.            (overflow  limit)
            //C
            //C     Brief outline of algorithm..
            //C
            //C     Phase 1    scans zero components.
            //C     Move to phase 2 when a component is nonzero and .LE. CUTLO
            //C     Move to phase 3 when a component is .GT. CUTLO
            //C     Move to phase 4 when a component is .GE. CUTHI/M
            //C     where M = N for X() real and M = 2*N for complex.
            //C
            //C     Values for CUTLO and CUTHI..
            //C     From the environmental parameters listed in the IMSL converter
            //C     document the limiting values are as follows..
            //C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
            //C                   Univac and DEC at 2**(-103)
            //C                   Thus CUTLO = 2**(-51) = 4.44089E-16
            //C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
            //C                   Thus CUTHI = 2**(63.5) = 1.30438E19
            //C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
            //C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
            //C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
            //C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
            //C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
            //C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
            //C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
            //C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
            //C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
            //C***ROUTINES CALLED  (NONE)
            //C***END PROLOGUE  SCNRM2
            bool IMAG, SCALE;
            int          NEXT,NN,I;
            double       HITEST, SUM, XMAX, ABSX, _SCNRM2;
            const double ONE = 1.;
            const double ZERO = 0.;
            const double CUTLO = sqrt(std::numeric_limits<double>::min()/std::numeric_limits<double>::epsilon());
            const double CUTHI = sqrt(std::numeric_limits<double>::max());
            //C***FIRST EXECUTABLE STATEMENT  SCNRM2
            if (N > 0) goto Line10;
            _SCNRM2 = ZERO;
            goto Line300;
            //C
Line10:
            NEXT = 30;
            SUM = ZERO;
            NN = N * INCX;
            //C                                                 BEGIN MAIN LOOP
            for (I=1; I<=NN; I+=INCX) { // DO 210 I=1,NN,INCX
                ABSX = fabs(real(CX(I)));
                IMAG = false;
                if (NEXT==30) goto Line30;
                if (NEXT==50) goto Line50;
                if (NEXT==70) goto Line70;
                if (NEXT==90) goto Line90;
                if (NEXT==110) goto Line110;
                throw;
Line30:
                if ( ABSX > CUTLO) goto Line85;
                NEXT = 50;
                SCALE = false;
                //C
                //C                        PHASE 1.  SUM IS ZERO
                //C
Line50:
                if ( ABSX == ZERO) goto Line200;
                if ( ABSX > CUTLO) goto Line85;
                //C
                //C                                PREPARE FOR PHASE 2.
                NEXT = 70;
                goto Line105;
                //C
                //C                                PREPARE FOR PHASE 4.
                //C
Line100:
                NEXT = 110;
                SUM = (SUM / ABSX) / ABSX;
Line105:
                SCALE = true;
                XMAX = ABSX;
                goto Line115;
                //C
                //C                   PHASE 2.  SUM IS SMALL.
                //C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
                //C
Line70:
                if ( ABSX > CUTLO ) goto Line75;
                //C
                //C                     COMMON CODE FOR PHASES 2 AND 4.
                //C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
                //C
Line110:
                if ( ABSX <= XMAX ) goto Line115;
                SUM = ONE + SUM * sqr(XMAX / ABSX);
                XMAX = ABSX;
                goto Line200;
                //C
Line115:
                SUM = SUM + sqr(ABSX/XMAX);
                goto Line200;
                //C
                //C
                //C                  PREPARE FOR PHASE 3.
                //C
Line75:
                SUM = (SUM * XMAX) * XMAX;
                //C
Line85:
                NEXT = 90;
                SCALE = false;
                //C
                //C     FOR REAL OR D.P. SET HITEST = CUTHI/N
                //C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
                //C
                HITEST = CUTHI/double( N );
                //C
                //C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
                //C
Line90:
                if (ABSX >= HITEST) goto Line100;
                SUM = SUM + sqr(ABSX);
Line200:
                ;
                //C                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS.
                //C
                if (IMAG) goto Line210;
                ABSX = fabs(imag(CX(I)));
                IMAG = true;
                if (NEXT==50) goto Line50;
                if (NEXT==70) goto Line70;
                if (NEXT==90) goto Line90;
                if (NEXT==110) goto Line110;
                throw;
Line210:
                ;
            }
            //C
            //C              END OF MAIN LOOP.
            //C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
            //C
            _SCNRM2 = sqrt(SUM);
            if (SCALE) _SCNRM2 = _SCNRM2 * XMAX;
Line300:
            return (_SCNRM2);
        }

        void SROTG(double& SA,double& SB,double& SC,double& SS)
        {
            //C***BEGIN PROLOGUE  SROTG
            //C***DATE WRITTEN   791001   (YYMMDD)
            //C***REVISION DATE  820801   (YYMMDD)
            //C***CATEGORY NO.  D1B10
            //C***KEYWORDS  BLAS,GIVENS ROTATION,LINEAR ALGEBRA,VECTOR
            //C***AUTHOR  LAWSON, C. L., (JPL)
            //C           HANSON, R. J., (SNLA)
            //C           KINCAID, D. R., (U. OF TEXAS)
            //C           KROGH, F. T., (JPL)
            //C***PURPOSE  Construct s.p. plane Givens rotation
            //C***DESCRIPTION
            //C
            //C                B L A S  Subprogram
            //C    Description of Parameters
            //C
            //C     --Input--
            //C       SA  single precision scalar
            //C       SB  single precision scalar
            //C
            //C     --Output--
            //C       SA  single precision result R
            //C       SB  single precision result Z
            //C       SC  single precision result
            //C       SS  single precision result
            //C
            //C     Designed by C. L. Lawson, JPL, 1977 Sept 08
            //C
            //C
            //C     Construct the Givens transformation
            //C
            //C         ( SC  SS )
            //C     G = (        ) ,    SC**2 + SS**2 = 1 ,
            //C         (-SS  SC )
            //C
            //C     which zeros the second entry of the 2-vector  (SA,SB)**T.
            //C
            //C     The quantity R = (+/-)SQRT(SA**2 + SB**2) overwrites SA in
            //C     storage.  The value of SB is overwritten by a value Z which
            //C     allows SC and SS to be recovered by the following algorithm:
            //C
            //C           If Z=1  set  SC=0.  and  SS=1.
            //C           If ABS(Z) .LT. 1  set  SC=SQRT(1-Z**2)  and  SS=Z
            //C           If ABS(Z) .GT. 1  set  SC=1/Z  and  SS=SQRT(1-SC**2)
            //C
            //C     Normally, the subprogram SROT(N,SX,INCX,SY,INCY,SC,SS) will
            //C     next be called to apply the transformation to a 2 by N matrix.
            //C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
            //C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
            //C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
            //C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
            //C***ROUTINES CALLED  (NONE)
            //C***END PROLOGUE  SROTG
            //C
            //C***FIRST EXECUTABLE STATEMENT  SROTG
            double U,V,R;
            if (fabs(SA) <= fabs(SB)) goto Line10;
            //C
            //C *** HERE ABS(SA) .GT. ABS(SB) ***
            //C
            U = SA + SA;
            V = SB / U;
            //C
            //C     NOTE THAT U AND R HAVE THE SIGN OF SA
            //C
            R = sqrt(.25 + sqr(V)) * U;
            //C
            //C     NOTE THAT SC IS POSITIVE
            //C
            SC = SA / R;
            SS = V * (SC + SC);
            SB = SS;
            SA = R;
            return;
            //C
            //C *** HERE ABS(SA) .LE. ABS(SB) ***
            //C
Line10:
            if (SB == 0.) goto Line20;
            U = SB + SB;
            V = SA / U;
            //C
            //C     NOTE THAT U AND R HAVE THE SIGN OF SB
            //C     (R IS IMMEDIATELY STORED IN SA)
            //C
            SA = sqrt(.25 + sqr(V)) * U;
            //C
            //C     NOTE THAT SS IS POSITIVE
            //C
            SS = SB / SA;
            SC = V * (SS + SS);
            if (SC == 0.) goto Line15;
            SB = 1. / SC;
            return;
Line15:
            SB = 1.;
            return;
            //C
            //C *** HERE SA = SB = 0. ***
            //C
Line20:
            SC = 1.;
            SS = 0.;
            return;
            //C
        }
    } //namespace CMLIB
} // namespace SCATMECH


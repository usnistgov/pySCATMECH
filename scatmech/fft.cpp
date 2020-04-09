//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: fft.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "matrixmath.h"
#include "fft.h"
#include <iostream>

using namespace std;
using namespace SCATMECH;


/// All of these routines have been converted to C++ for SCATMECH by T.A. Germer 3/23/2010

namespace SCATMECH {

    namespace {
        // The following map stores work arrays for CFFTB and CFFTF...
        map<int,vector<double> > fftwork;
    }

    void fft1d(std::valarray<COMPLEX>& data,int isign)
    {
        int n=data.size();

        // If this size has not been done before, we need to create a work array.
        // This only needs to be done once for each size...
        if (fftwork.find(n)==fftwork.end()) {
            fftwork[n] = vector<double>(4*n+15);
            vector<double> &work = fftwork[n];
            CMLIB::CFFTI(n,(double*)(&(work[0])));
        }

        vector<double> &work = fftwork[n];

        if (isign==1) {
            CMLIB::CFFTB(n,(double*)(&data[0]),(double*)(&(work[0])));
        } else {
            CMLIB::CFFTF(n,(double*)(&data[0]),(double*)(&(work[0])));
        }
    }

    void fft1d(CFARRAY data,int N,int isign)
    {
        int n=N;

        // If this size has not been done before, we need to create a work array.
        // This only needs to be done once for each size...
        if (fftwork.find(n)==fftwork.end()) {
            fftwork[n] = vector<double>(4*n+15);
            vector<double> &work = fftwork[n];
            CMLIB::CFFTI(n,(double*)(&(work[0])));
        }

        vector<double> &work = fftwork[n];

        if (isign==1) {
            CMLIB::CFFTB(n,(double*)(&data(1)),(double*)(&(work[0])));
        } else {
            CMLIB::CFFTF(n,(double*)(&data(1)),(double*)(&(work[0])));
        }
    }

    namespace CMLIB {

        void CFFTF(int& N,DFARRAY C,DFARRAY WSAVE)
        {
            //C***BEGIN PROLOGUE  CFFTF
            //C***DATE WRITTEN   790601   (YYMMDD)
            //C***REVISION DATE  800626   (YYMMDD)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***CATEGORY NO.  J1A2
            //C***KEYWORDS  FOURIER TRANSFORM
            //C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
            //C***PURPOSE  Forward transform of a complex, periodic sequence.
            //C***DESCRIPTION
            //C
            //C  Subroutine CFFTF computes the forward complex discrete Fourier
            //C  transform (the Fourier analysis).  Equivalently, CFFTF computes
            //C  the Fourier coefficients of a complex periodic sequence.
            //C  The transform is defined below at output parameter C.
            //C
            //C  The transform is not normalized.  To obtain a normalized transform
            //C  the output must be divided by N.  Otherwise a call of CFFTF
            //C  followed by a call of CFFTB will multiply the sequence by N.
            //C
            //C  The array WSAVE which is used by subroutine CFFTF must be
            //C  initialized by calling subroutine CFFTI(N,WSAVE).
            //C
            //C  Input Parameters
            //C
            //C
            //C  N      the length of the complex sequence C.  The method is
            //C         more efficient when N is the product of small primes.
            //C
            //C  C      a complex array of length N which contains the sequence
            //C
            //C  WSAVE   a real work array which must be dimensioned at least 4*N+15
            //C          in the program that calls CFFTF.  The WSAVE array must be
            //C          initialized by calling subroutine CFFTI(N,WSAVE), and a
            //C          different WSAVE array must be used for each different
            //C          value of N.  This initialization does not have to be
            //C          repeated so long as N remains unchanged.  Thus subsequent
            //C          transforms can be obtained faster than the first.
            //C          The same WSAVE array can be used by CFFTF and CFFTB.
            //C
            //C  Output Parameters
            //C
            //C  C      for J=1,...,N
            //C
            //C             C(J)=the sum from K=1,...,N of
            //C
            //C                   C(K)*EXP(-I*J*K*2*PI/N)
            //C
            //C                         where I=SQRT(-1)
            //C
            //C  WSAVE   contains initialization calculations which must not be
            //C          destroyed between calls of subroutine CFFTF or CFFTB
            //C***REFERENCES  (NONE)
            //C***ROUTINES CALLED  CFFTF1
            //C***END PROLOGUE  CFFTF
            //C***FIRST EXECUTABLE STATEMENT  CFFTF
            if (N==1) return;
            int IW1 = N+N+1;
            int IW2 = IW1+N+N;
            CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2));
            return;
        }

        void CFFTB(int& N,DFARRAY C,DFARRAY WSAVE)
        {
            //C***BEGIN PROLOGUE  CFFTB
            //C***DATE WRITTEN   790601   (YYMMDD)
            //C***REVISION DATE  830401   (YYMMDD)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***CATEGORY NO.  J1A2
            //C***KEYWORDS  FOURIER TRANSFORM
            //C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
            //C***PURPOSE  Unnormalized inverse of CFFTF.
            //C***DESCRIPTION
            //C
            //C  Subroutine CFFTB computes the backward complex discrete Fourier
            //C  transform (the Fourier synthesis).  Equivalently, CFFTB computes
            //C  a complex periodic sequence from its Fourier coefficients.
            //C  The transform is defined below at output parameter C.
            //C
            //C  A call of CFFTF followed by a call of CFFTB will multiply the
            //C  sequence by N.
            //C
            //C  The array WSAVE which is used by subroutine CFFTB must be
            //C  initialized by calling subroutine CFFTI(N,WSAVE).
            //C
            //C  Input Parameters
            //C
            //C
            //C  N      the length of the complex sequence C.  The method is
            //C         more efficient when N is the product of small primes.
            //C
            //C  C      a complex array of length N which contains the sequence
            //C
            //C  WSAVE   a real work array which must be dimensioned at least 4*N+15
            //C          in the program that calls CFFTB.  The WSAVE array must be
            //C          initialized by calling subroutine CFFTI(N,WSAVE), and a
            //C          different WSAVE array must be used for each different
            //C          value of N.  This initialization does not have to be
            //C          repeated so long as N remains unchanged.  Thus subsequent
            //C          transforms can be obtained faster than the first.
            //C          The same WSAVE array can be used by CFFTF and CFFTB.
            //C
            //C  Output Parameters
            //C
            //C  C      For J=1,...,N
            //C
            //C             C(J)=the sum from K=1,...,N of
            //C
            //C                   C(K)*EXP(I*J*K*2*PI/N)
            //C
            //C                         where I=SQRT(-1)
            //C
            //C  WSAVE   contains initialization calculations which must not be
            //C          destroyed between calls of subroutine CFFTF or CFFTB
            //C***REFERENCES  (NONE)
            //C***ROUTINES CALLED  CFFTB1
            //C***END PROLOGUE  CFFTB
            //C***FIRST EXECUTABLE STATEMENT  CFFTB
            if (N == 1) return;
            int IW1 = N+N+1;
            int IW2 = IW1+N+N;
            CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2));
            return;
        }

        void CFFTF1(int& N,DFARRAY C,DFARRAY CH,DFARRAY WA,DFARRAY IFAC)
        {
            //C***BEGIN PROLOGUE  CFFTF1
            //C***REFER TO  CFFTF
            //C***ROUTINES CALLED  PASSF,PASSF2,PASSF3,PASSF4,PASSF5
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  CFFTF1
            //C***FIRST EXECUTABLE STATEMENT  CFFTF1
            int NF,NA,L1,IW,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,N2,I,K1,NAC;
            NF = (int)IFAC(2);
            NA = 0;
            L1 = 1;
            IW = 1;
            for (K1=1; K1<=NF; ++K1) {
                IP = (int)IFAC(K1+2);
                L2 = IP*L1;
                IDO = N/L2;
                IDOT = IDO+IDO;
                IDL1 = IDOT*L1;
                if (IP != 4) goto L103;
                IX2 = IW+IDOT;
                IX3 = IX2+IDOT;
                if (NA != 0) goto L101;
                PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3));
                goto L102;
L101:
                PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3));
L102:
                NA = 1-NA;
                goto L115;
L103:
                if (IP != 2) goto L106;
                if (NA != 0) goto L104;
                PASSF2 (IDOT,L1,C,CH,WA(IW));
                goto L105;
L104:
                PASSF2 (IDOT,L1,CH,C,WA(IW));
L105:
                NA = 1-NA;
                goto L115;
L106:
                if (IP != 3) goto L109;
                IX2 = IW+IDOT;
                if (NA != 0) goto L107;
                PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2));
                goto L108;
L107:
                PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2));
L108:
                NA = 1-NA;
                goto L115;
L109:
                if (IP != 5) goto L112;
                IX2 = IW+IDOT;
                IX3 = IX2+IDOT;
                IX4 = IX3+IDOT;
                if (NA != 0) goto L110;
                PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4));
                goto L111;
L110:
                PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4));
L111:
                NA = 1-NA;
                goto L115;
L112:
                if (NA != 0) goto L113;
                PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW));
                goto L114;
L113:
                PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW));
L114:
                if (NAC != 0) NA = 1-NA;
L115:
                L1 = L2;
                IW = IW+(IP-1)*IDOT;
            }
            if (NA == 0) return;
            N2 = N+N;
            for (I=1; I<=N2; ++I) {
                C(I) = CH(I);
            }
            return;
        }

        void CFFTB1(int& N,DFARRAY C,DFARRAY CH,DFARRAY WA,DFARRAY IFAC)
        {
            //C***BEGIN PROLOGUE  CFFTB1
            //C***REFER TO  CFFTB
            //C***ROUTINES CALLED  PASSB,PASSB2,PASSB3,PASSB4,PASSB5
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  CFFTB1
            //C***FIRST EXECUTABLE STATEMENT  CFFTB1
            int NF,NA,L1,IW,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,N2,I,K1,NAC;
            NF = (int)IFAC(2);
            NA = 0;
            L1 = 1;
            IW = 1;
            for (K1=1; K1<=NF; ++K1) {
                IP = (int)IFAC(K1+2);
                L2 = IP*L1;
                IDO = N/L2;
                IDOT = IDO+IDO;
                IDL1 = IDOT*L1;
                if (IP != 4) goto L103;
                IX2 = IW+IDOT;
                IX3 = IX2+IDOT;
                if (NA != 0) goto L101;
                PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3));
                goto L102;
L101:
                PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3));
L102:
                NA = 1-NA;
                goto L115;
L103:
                if (IP != 2) goto L106;
                if (NA != 0) goto L104;
                PASSB2 (IDOT,L1,C,CH,WA(IW));
                goto L105;
L104:
                PASSB2 (IDOT,L1,CH,C,WA(IW));
L105:
                NA = 1-NA;
                goto L115;
L106:
                if (IP != 3) goto L109;
                IX2 = IW+IDOT;
                if (NA != 0) goto L107;
                PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2));
                goto L108;
L107:
                PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2));
L108:
                NA = 1-NA;
                goto L115;
L109:
                if (IP != 5) goto L112;
                IX2 = IW+IDOT;
                IX3 = IX2+IDOT;
                IX4 = IX3+IDOT;
                if (NA != 0) goto L110;
                PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4));
                goto L111;
L110:
                PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4));
L111:
                NA = 1-NA;
                goto L115;
L112:
                if (NA != 0) goto L113;
                PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW));
                goto L114;
L113:
                PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW));
L114:
                if (NAC != 0) NA = 1-NA;
L115:
                L1 = L2;
                IW = IW+(IP-1)*IDOT;
            }
            if (NA == 0) return;
            N2 = N+N;
            for (I=1; I<=N2; ++I) {
                C(I) = CH(I);
            }
            return;
        }

        void CFFTI(int& N,DFARRAY WSAVE)
        {
            //C***BEGIN PROLOGUE  CFFTI
            //C***DATE WRITTEN   790601   (YYMMDD)
            //C***REVISION DATE  830401   (YYMMDD)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***CATEGORY NO.  J1A2
            //C***KEYWORDS  FOURIER TRANSFORM
            //C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
            //C***PURPOSE  Initialize for CFFTF and CFFTB.
            //C***DESCRIPTION
            //C
            //C  Subroutine CFFTI initializes the array WSAVE which is used in
            //C  both CFFTF and CFFTB.  The prime factorization of N together with
            //C  a tabulation of the trigonometric functions are computed and
            //C  stored in WSAVE.
            //C
            //C  Input Parameter
            //C
            //C  N       the length of the sequence to be transformed
            //C
            //C  Output Parameter
            //C
            //C  WSAVE   a work array which must be dimensioned at least 4*N+15.
            //C          The same work array can be used for both CFFTF and CFFTB
            //C          as long as N remains unchanged.  Different WSAVE arrays
            //C          are required for different values of N.  The contents of
            //C          WSAVE must not be changed between calls of CFFTF or CFFTB.
            //C***REFERENCES  (NONE)
            //C***ROUTINES CALLED  CFFTI1
            //C***END PROLOGUE  CFFTI
            //C***FIRST EXECUTABLE STATEMENT  CFFTI
            if (N == 1) return;
            int IW1 = N+N+1;
            int IW2 = IW1+N+N;
            CFFTI1 (N,WSAVE(IW1),WSAVE(IW2));
            return;
        }

        void CFFTI1(int& N,DFARRAY WA,DFARRAY IFAC)
        {
            //C***BEGIN PROLOGUE  CFFTI1
            //C***REFER TO  CFFTI
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  CFFTI1
            const int NTRYH[5] = {0,3,4,2,5};
            int NL,NF,J,NQ,NR,NTRY,IB,I,IP,LD,L2,IDO,IDOT,IPM,I1,II,L1,K1;
            double TPI,ARGH,FI,ARGLD,ARG;
            //C***FIRST EXECUTABLE STATEMENT  CFFTI1
            NL = N;
            NF = 0;
            J = 0;
L101:
            J = J+1;
            if (J-4>0) goto L103;
            NTRY = NTRYH[J];
            goto L104;
L103:
            NTRY = NTRY+2;
L104:
            NQ = NL/NTRY;
            NR = NL-NTRY*NQ;
            if (NR != 0) goto L101;
            NF = NF+1;
            IFAC(NF+2) = NTRY;
            NL = NQ;
            if (NTRY != 2) goto L107;
            if (NF == 1) goto L107;
            for (I=1; I<=NF; ++I) {
                IB = NF-I+2;
                IFAC(IB+2) = IFAC(IB+1);
            }
            IFAC(3) = 2;
L107:
            if (NL != 1) goto L104;
            IFAC(1) = N;
            IFAC(2) = NF;
            TPI = 6.28318530717959;
            ARGH = TPI/(double)(N);
            I = 2;
            L1 = 1;
            for (K1=1; K1<=NF; ++K1) {
                IP = (int)IFAC(K1+2);
                LD = 0;
                L2 = L1*IP;
                IDO = N/L2;
                IDOT = IDO+IDO+2;
                IPM = IP-1;
                for (J=1; J<=IPM; ++J) {
                    I1 = I;
                    WA(I-1) = 1.;
                    WA(I) = 0.;
                    LD = LD+L1;
                    FI = 0.;
                    ARGLD = (double)(LD)*ARGH;
                    for (II=4; II<=IDOT; II+=2) {
                        I = I+2;
                        FI = FI+1.;
                        ARG = FI*ARGLD;
                        WA(I-1) = cos(ARG);
                        WA(I) = sin(ARG);
                    }
                    if (IP > 5) {
                        WA(I1-1) = WA(I-1);
                        WA(I1) = WA(I);
                    }
                }
                L1 = L2;
            }
            return;
        }


        void PASSF(int& NAC,int& IDO,int& IP,int& L1,int& IDL1,DFARRAY CC,DFARRAY C1,DFARRAY C2,DFARRAY CH,DFARRAY CH2,DFARRAY WA)
        {
            //C***BEGIN PROLOGUE  PASSF
            //C***REFER TO  CFFTF
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSF
            CH.array(IDO,L1,IP);
            CC.array(IDO,IP,L1);
            C1.array(IDO,L1,IP);
            C2.array(IDL1,IP);
            CH2.array(IDL1,IP);
            int IDOT,NT,IPP2,IPPH,IDP,J,K,JC,I,INC,IDL,IDLJ,L,LC,IK,IDIJ,IDJ;
            double WAR,WAI;
            //C***FIRST EXECUTABLE STATEMENT  PASSF
            IDOT = IDO/2;
            NT = IP*IDL1;
            IPP2 = IP+2;
            IPPH = (IP+1)/2;
            IDP = IP*IDO;
            //C
            if (IDO < L1) goto L106;
            for (J=2; J<=IPPH; ++J) {
                JC = IPP2-J;
                for (K=1; K<=L1; ++K) {
                    #pragma omp parallel for
                    for (I=1; I<=IDO; ++I) {
                        CH(I,K,J) = CC(I,J,K)+CC(I,JC,K);
                        CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K);
                    }
                }
            }
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=1; I<=IDO; ++I) {
                    CH(I,K,1) = CC(I,1,K);
                }
            }
            goto L112;
L106:
            for (J=2; J<=IPPH; ++J) {
                JC = IPP2-J;
                for (I=1; I<=IDO; ++I) {
                    #pragma omp parallel for
                    for (K=1; K<=L1; ++K) {
                        CH(I,K,J) = CC(I,J,K)+CC(I,JC,K);
                        CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K);
                    }
                }
            }
            for (I=1; I<=IDO; ++I) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    CH(I,K,1) = CC(I,1,K);
                }
            }
L112:
            IDL = 2-IDO;
            INC = 0;
            for (L=2; L<=IPPH; ++L) {
                LC = IPP2-L;
                IDL = IDL+IDO;
                #pragma omp parallel for
                for (IK=1; IK<=IDL1; ++IK) {
                    C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2);
                    C2(IK,LC) = -WA(IDL)*CH2(IK,IP);
                }
                IDLJ = IDL;
                INC = INC+IDO;
                for (J=3; J<=IPPH; ++J) {
                    JC = IPP2-J;
                    IDLJ = IDLJ+INC;
                    if (IDLJ > IDP) IDLJ = IDLJ-IDP;
                    WAR = WA(IDLJ-1);
                    WAI = WA(IDLJ);
                    #pragma omp parallel for
                    for (IK=1; IK<=IDL1; ++IK) {
                        C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J);
                        C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC);
                    }
                }
            }
            for (J=2; J<=IPPH; ++J) {
                #pragma omp parallel for
                for (IK=1; IK<=IDL1; ++IK) {
                    CH2(IK,1) = CH2(IK,1)+CH2(IK,J);
                }
            }
            for (J=2; J<=IPPH; ++J) {
                JC = IPP2-J;
                #pragma omp parallel for
                for (IK=2; IK<=IDL1; IK+=2) {
                    CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC);
                    CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC);
                    CH2(IK,J) = C2(IK,J)+C2(IK-1,JC);
                    CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC);
                }
            }
            NAC = 1;
            if (IDO == 2) return;
            NAC = 0;
            #pragma omp parallel for
            for (IK=1; IK<=IDL1; ++IK) {
                C2(IK,1) = CH2(IK,1);
            }
            for (J=2; J<=IP; ++J) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    C1(1,K,J) = CH(1,K,J);
                    C1(2,K,J) = CH(2,K,J);
                }
            }
            if (IDOT > L1) goto L127;
            IDIJ = 0;
            for (J=2; J<=IP; ++J) {
                IDIJ = IDIJ+2;
                for (I=4; I<=IDO; I+=2) {
                    IDIJ = IDIJ+2;
                    #pragma omp parallel for
                    for (K=1; K<=L1; ++K) {
                        C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J);
                        C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J);
                    }
                }
            }
            return;
L127:
            IDJ = 2-IDO;
            for (J=2; J<=IP; ++J) {
                IDJ = IDJ+IDO;
                for (K=1; K<=L1; ++K) {
                    IDIJ = IDJ;
                    #pragma omp parallel for
                    for (I=4; I<=IDO; I+=2) {
                        IDIJ = IDIJ+2;
                        C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J);
                        C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J);
                    }
                }
            }
            return;
        }



        void PASSF2(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1)
        {
            //C***BEGIN PROLOGUE  PASSF2
            //C***REFER TO  CFFTF
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSF2
            CC.array(IDO,2,L1);
            CH.array(IDO,L1,2);
            int K,I;
            double TR2,TI2;
            //C***FIRST EXECUTABLE STATEMENT  PASSF2
            if (IDO > 2) goto L102;
            for (K=1; K<=L1; ++K) {
                CH(1,K,1) = CC(1,1,K)+CC(1,2,K);
                CH(1,K,2) = CC(1,1,K)-CC(1,2,K);
                CH(2,K,1) = CC(2,1,K)+CC(2,2,K);
                CH(2,K,2) = CC(2,1,K)-CC(2,2,K);
            }
            return;
L102:
            if (IDO/2<L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K);
                    TR2 = CC(I-1,1,K)-CC(I-1,2,K);
                    CH(I,K,1) = CC(I,1,K)+CC(I,2,K);
                    TI2 = CC(I,1,K)-CC(I,2,K);
                    CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2;
                    CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K);
                    TR2 = CC(I-1,1,K)-CC(I-1,2,K);
                    CH(I,K,1) = CC(I,1,K)+CC(I,2,K);
                    TI2 = CC(I,1,K)-CC(I,2,K);
                    CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2;
                    CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2;
                }
            }
            return;
        }


        void PASSF3(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2)
        {
            //C***BEGIN PROLOGUE  PASSF3
            //C***REFER TO  CFFTF
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSF3
            CC.array(IDO,3,L1);
            CH.array(IDO,L1,3);
            const double TAUR = -0.5;
            const double TAUI = -0.866025403784439;
            int K,I;
            double TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3;
            //C***FIRST EXECUTABLE STATEMENT  PASSF3
            if (IDO != 2) goto L102;
            for (K=1; K<=L1; ++K) {
                TR2 = CC(1,2,K)+CC(1,3,K);
                CR2 = CC(1,1,K)+TAUR*TR2;
                CH(1,K,1) = CC(1,1,K)+TR2;
                TI2 = CC(2,2,K)+CC(2,3,K);
                CI2 = CC(2,1,K)+TAUR*TI2;
                CH(2,K,1) = CC(2,1,K)+TI2;
                CR3 = TAUI*(CC(1,2,K)-CC(1,3,K));
                CI3 = TAUI*(CC(2,2,K)-CC(2,3,K));
                CH(1,K,2) = CR2-CI3;
                CH(1,K,3) = CR2+CI3;
                CH(2,K,2) = CI2+CR3;
                CH(2,K,3) = CI2-CR3;
            }
            return;
L102:
            if (IDO/2<L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    TR2 = CC(I-1,2,K)+CC(I-1,3,K);
                    CR2 = CC(I-1,1,K)+TAUR*TR2;
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2;
                    TI2 = CC(I,2,K)+CC(I,3,K);
                    CI2 = CC(I,1,K)+TAUR*TI2;
                    CH(I,K,1) = CC(I,1,K)+TI2;
                    CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K));
                    CI3 = TAUI*(CC(I,2,K)-CC(I,3,K));
                    DR2 = CR2-CI3;
                    DR3 = CR2+CI3;
                    DI2 = CI2+CR3;
                    DI3 = CI2-CR3;
                    CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2;
                    CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2;
                    CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3;
                    CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    TR2 = CC(I-1,2,K)+CC(I-1,3,K);
                    CR2 = CC(I-1,1,K)+TAUR*TR2;
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2;
                    TI2 = CC(I,2,K)+CC(I,3,K);
                    CI2 = CC(I,1,K)+TAUR*TI2;
                    CH(I,K,1) = CC(I,1,K)+TI2;
                    CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K));
                    CI3 = TAUI*(CC(I,2,K)-CC(I,3,K));
                    DR2 = CR2-CI3;
                    DR3 = CR2+CI3;
                    DI2 = CI2+CR3;
                    DI3 = CI2-CR3;
                    CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2;
                    CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2;
                    CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3;
                    CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3;
                }
            }
            return;
        }


        void PASSF4(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3)
        {
            //C***BEGIN PROLOGUE  PASSF4
            //C***REFER TO  CFFTF
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSF4
            CC.array(IDO,4,L1);
            CH.array(IDO,L1,4);
            int K,I;
            double TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4;
            //C***FIRST EXECUTABLE STATEMENT  PASSF4
            if (IDO != 2) goto L102;
            for (K=1; K<=L1; ++K) {
                TI1 = CC(2,1,K)-CC(2,3,K);
                TI2 = CC(2,1,K)+CC(2,3,K);
                TR4 = CC(2,2,K)-CC(2,4,K);
                TI3 = CC(2,2,K)+CC(2,4,K);
                TR1 = CC(1,1,K)-CC(1,3,K);
                TR2 = CC(1,1,K)+CC(1,3,K);
                TI4 = CC(1,4,K)-CC(1,2,K);
                TR3 = CC(1,2,K)+CC(1,4,K);
                CH(1,K,1) = TR2+TR3;
                CH(1,K,3) = TR2-TR3;
                CH(2,K,1) = TI2+TI3;
                CH(2,K,3) = TI2-TI3;
                CH(1,K,2) = TR1+TR4;
                CH(1,K,4) = TR1-TR4;
                CH(2,K,2) = TI1+TI4;
                CH(2,K,4) = TI1-TI4;
            }
            return;
L102:
            if (IDO/2<L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    TI1 = CC(I,1,K)-CC(I,3,K);
                    TI2 = CC(I,1,K)+CC(I,3,K);
                    TI3 = CC(I,2,K)+CC(I,4,K);
                    TR4 = CC(I,2,K)-CC(I,4,K);
                    TR1 = CC(I-1,1,K)-CC(I-1,3,K);
                    TR2 = CC(I-1,1,K)+CC(I-1,3,K);
                    TI4 = CC(I-1,4,K)-CC(I-1,2,K);
                    TR3 = CC(I-1,2,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = TR2+TR3;
                    CR3 = TR2-TR3;
                    CH(I,K,1) = TI2+TI3;
                    CI3 = TI2-TI3;
                    CR2 = TR1+TR4;
                    CR4 = TR1-TR4;
                    CI2 = TI1+TI4;
                    CI4 = TI1-TI4;
                    CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2;
                    CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2;
                    CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3;
                    CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3;
                    CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4;
                    CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    TI1 = CC(I,1,K)-CC(I,3,K);
                    TI2 = CC(I,1,K)+CC(I,3,K);
                    TI3 = CC(I,2,K)+CC(I,4,K);
                    TR4 = CC(I,2,K)-CC(I,4,K);
                    TR1 = CC(I-1,1,K)-CC(I-1,3,K);
                    TR2 = CC(I-1,1,K)+CC(I-1,3,K);
                    TI4 = CC(I-1,4,K)-CC(I-1,2,K);
                    TR3 = CC(I-1,2,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = TR2+TR3;
                    CR3 = TR2-TR3;
                    CH(I,K,1) = TI2+TI3;
                    CI3 = TI2-TI3;
                    CR2 = TR1+TR4;
                    CR4 = TR1-TR4;
                    CI2 = TI1+TI4;
                    CI4 = TI1-TI4;
                    CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2;
                    CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2;
                    CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3;
                    CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3;
                    CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4;
                    CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4;
                }
            }
            return;
        }

        void PASSF5(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3,DFARRAY WA4)
        {
            //C***BEGIN PROLOGUE  PASSF5
            //C***REFER TO  CFFTF
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSF5
            CC.array(IDO,5,L1);
            CH.array(IDO,L1,5);
            const double TR11 = .309016994374947;
            const double TI11 = -.951056516295154;
            const double TR12 = -.809016994374947;
            const double TI12 = -.587785252292473;
            int K,I;
            double TI5,TI2,TI4,TI3,TR5,TR2,TR4,TR3,CR2,CI2,CR3,CI3,CR5,CI5,CR4,CI4;
            double DR3,DR4,DI3,DI4,DR5,DR2,DI5,DI2;
            //C***FIRST EXECUTABLE STATEMENT  PASSF5
            if (IDO != 2) goto L102;
            for (K=1; K<=L1; ++K) {
                TI5 = CC(2,2,K)-CC(2,5,K);
                TI2 = CC(2,2,K)+CC(2,5,K);
                TI4 = CC(2,3,K)-CC(2,4,K);
                TI3 = CC(2,3,K)+CC(2,4,K);
                TR5 = CC(1,2,K)-CC(1,5,K);
                TR2 = CC(1,2,K)+CC(1,5,K);
                TR4 = CC(1,3,K)-CC(1,4,K);
                TR3 = CC(1,3,K)+CC(1,4,K);
                CH(1,K,1) = CC(1,1,K)+TR2+TR3;
                CH(2,K,1) = CC(2,1,K)+TI2+TI3;
                CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3;
                CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3;
                CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3;
                CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3;
                CR5 = TI11*TR5+TI12*TR4;
                CI5 = TI11*TI5+TI12*TI4;
                CR4 = TI12*TR5-TI11*TR4;
                CI4 = TI12*TI5-TI11*TI4;
                CH(1,K,2) = CR2-CI5;
                CH(1,K,5) = CR2+CI5;
                CH(2,K,2) = CI2+CR5;
                CH(2,K,3) = CI3+CR4;
                CH(1,K,3) = CR3-CI4;
                CH(1,K,4) = CR3+CI4;
                CH(2,K,4) = CI3-CR4;
                CH(2,K,5) = CI2-CR5;
            }
            return;
L102:
            if (IDO/2 < L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    TI5 = CC(I,2,K)-CC(I,5,K);
                    TI2 = CC(I,2,K)+CC(I,5,K);
                    TI4 = CC(I,3,K)-CC(I,4,K);
                    TI3 = CC(I,3,K)+CC(I,4,K);
                    TR5 = CC(I-1,2,K)-CC(I-1,5,K);
                    TR2 = CC(I-1,2,K)+CC(I-1,5,K);
                    TR4 = CC(I-1,3,K)-CC(I-1,4,K);
                    TR3 = CC(I-1,3,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3;
                    CH(I,K,1) = CC(I,1,K)+TI2+TI3;
                    CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3;
                    CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3;
                    CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3;
                    CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3;
                    CR5 = TI11*TR5+TI12*TR4;
                    CI5 = TI11*TI5+TI12*TI4;
                    CR4 = TI12*TR5-TI11*TR4;
                    CI4 = TI12*TI5-TI11*TI4;
                    DR3 = CR3-CI4;
                    DR4 = CR3+CI4;
                    DI3 = CI3+CR4;
                    DI4 = CI3-CR4;
                    DR5 = CR2+CI5;
                    DR2 = CR2-CI5;
                    DI5 = CI2-CR5;
                    DI2 = CI2+CR5;
                    CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2;
                    CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2;
                    CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3;
                    CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3;
                    CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4;
                    CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4;
                    CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5;
                    CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    TI5 = CC(I,2,K)-CC(I,5,K);
                    TI2 = CC(I,2,K)+CC(I,5,K);
                    TI4 = CC(I,3,K)-CC(I,4,K);
                    TI3 = CC(I,3,K)+CC(I,4,K);
                    TR5 = CC(I-1,2,K)-CC(I-1,5,K);
                    TR2 = CC(I-1,2,K)+CC(I-1,5,K);
                    TR4 = CC(I-1,3,K)-CC(I-1,4,K);
                    TR3 = CC(I-1,3,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3;
                    CH(I,K,1) = CC(I,1,K)+TI2+TI3;
                    CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3;
                    CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3;
                    CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3;
                    CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3;
                    CR5 = TI11*TR5+TI12*TR4;
                    CI5 = TI11*TI5+TI12*TI4;
                    CR4 = TI12*TR5-TI11*TR4;
                    CI4 = TI12*TI5-TI11*TI4;
                    DR3 = CR3-CI4;
                    DR4 = CR3+CI4;
                    DI3 = CI3+CR4;
                    DI4 = CI3-CR4;
                    DR5 = CR2+CI5;
                    DR2 = CR2-CI5;
                    DI5 = CI2-CR5;
                    DI2 = CI2+CR5;
                    CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2;
                    CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2;
                    CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3;
                    CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3;
                    CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4;
                    CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4;
                    CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5;
                    CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5;
                }
            }
            return;
        }

        void PASSB(int& NAC,int& IDO,int& IP,int& L1,int& IDL1,DFARRAY CC,DFARRAY C1,DFARRAY C2,DFARRAY CH,DFARRAY CH2,DFARRAY WA)
        {
            //C***BEGIN PROLOGUE  PASSB
            //C***REFER TO  CFFTB
            //C***ROUTINES CALLED  (NONE)
            //C***END PROLOGUE  PASSB
            CH.array(IDO,L1,IP);
            CC.array(IDO,IP,L1);
            C1.array(IDO,L1,IP);
            C2.array(IDL1,IP);
            CH2.array(IDL1,IP);
            int IDOT,NT,IPP2,IPPH,IDP,J,K,JC,I,INC,IDL,IDLJ,L,LC,IK,IDIJ,IDJ;
            double WAR,WAI;
            //C***FIRST EXECUTABLE STATEMENT  PASSB
            IDOT = IDO/2;
            NT = IP*IDL1;
            IPP2 = IP+2;
            IPPH = (IP+1)/2;
            IDP = IP*IDO;
            //C
            if (IDO < L1) goto L106;
            for (J=2; J<=IPPH; ++J) {
                JC = IPP2-J;
                for (K=1; K<=L1; ++K) {
                    #pragma omp parallel for
                    for (I=1; I<=IDO; ++I) {
                        CH(I,K,J) = CC(I,J,K)+CC(I,JC,K);
                        CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K);
                    }
                }
            }
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=1; I<=IDO; ++I) {
                    CH(I,K,1) = CC(I,1,K);
                }
            }
            goto L112;
L106:
            for (J=2; J<=IPPH; ++J) {
                JC = IPP2-J;
                for (I=1; I<=IDO; ++I) {
                    #pragma omp parallel for
                    for (K=1; K<=L1; ++K) {
                        CH(I,K,J) = CC(I,J,K)+CC(I,JC,K);
                        CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K);
                    }
                }
            }
            for (I=1; I<=IDO; ++I) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    CH(I,K,1) = CC(I,1,K);
                }
            }
L112:
            IDL = 2-IDO;
            INC = 0;
            for (L=2; L<=IPPH; ++L) {
                LC = IPP2-L;
                IDL = IDL+IDO;
                #pragma omp parallel for
                for (IK=1; IK<=IDL1; ++IK) {
                    C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2);
                    C2(IK,LC) = WA(IDL)*CH2(IK,IP);
                }
                IDLJ = IDL;
                INC = INC+IDO;
                for (J=3; J<=IPPH; ++J) {
                    JC = IPP2-J;
                    IDLJ = IDLJ+INC;
                    if (IDLJ > IDP) IDLJ = IDLJ-IDP;
                    WAR = WA(IDLJ-1);
                    WAI = WA(IDLJ);
                    #pragma omp parallel for
                    for (IK=1; IK<=IDL1; ++IK) {
                        C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J);
                        C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC);
                    }
                }
            }
            for (J=2; J<=IPPH; ++J) {
                #pragma omp parallel for
                for (IK=1; IK<=IDL1; ++IK) {
                    CH2(IK,1) = CH2(IK,1)+CH2(IK,J);
                }
            }
            for (J=2; J<=IPPH; ++J) {
                JC = IPP2-J;
                #pragma omp parallel for
                for (IK=2; IK<=IDL1; IK+=2) {
                    CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC);
                    CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC);
                    CH2(IK,J) = C2(IK,J)+C2(IK-1,JC);
                    CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC);
                }
            }
            NAC = 1;
            if (IDO == 2) return;
            NAC = 0;
            for (IK=1; IK<=IDL1; ++IK) {
                C2(IK,1) = CH2(IK,1);
            }
            for (J=2; J<=IP; ++J) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    C1(1,K,J) = CH(1,K,J);
                    C1(2,K,J) = CH(2,K,J);
                }
            }
            if (IDOT > L1) goto L127;
            IDIJ = 0;
            for (J=2; J<=IP; ++J) {
                IDIJ = IDIJ+2;
                for (I=4; I<=IDO; I+=2) {
                    IDIJ = IDIJ+2;
                    #pragma omp parallel for
                    for (K=1; K<=L1; ++K) {
                        C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J);
                        C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J);
                    }
                }
            }
            return;
L127:
            IDJ = 2-IDO;
            for (J=2; J<=IP; ++J) {
                IDJ = IDJ+IDO;
                for (K=1; K<=L1; ++K) {
                    IDIJ = IDJ;
                    #pragma omp parallel for
                    for (I=4; I<=IDO; I+=2) {
                        IDIJ = IDIJ+2;
                        C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J);
                        C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J);
                    }
                }
            }
            return;
        }


        void PASSB2(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1)
        {
            //C***BEGIN PROLOGUE  PASSB2
            //C***REFER TO  CFFTB
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSB2
            CC.array(IDO,2,L1);
            CH.array(IDO,L1,2);
            int K,I;
            double TR2,TI2;
            //C***FIRST EXECUTABLE STATEMENT  PASSB2
            if (IDO > 2) goto L102;
            for (K=1; K<=L1; ++K) {
                CH(1,K,1) = CC(1,1,K)+CC(1,2,K);
                CH(1,K,2) = CC(1,1,K)-CC(1,2,K);
                CH(2,K,1) = CC(2,1,K)+CC(2,2,K);
                CH(2,K,2) = CC(2,1,K)-CC(2,2,K);
            }
            return;
L102:
            if (IDO/2 < L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K);
                    TR2 = CC(I-1,1,K)-CC(I-1,2,K);
                    CH(I,K,1) = CC(I,1,K)+CC(I,2,K);
                    TI2 = CC(I,1,K)-CC(I,2,K);
                    CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2;
                    CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K);
                    TR2 = CC(I-1,1,K)-CC(I-1,2,K);
                    CH(I,K,1) = CC(I,1,K)+CC(I,2,K);
                    TI2 = CC(I,1,K)-CC(I,2,K);
                    CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2;
                    CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2;
                }
            }
            return;
        }

        void PASSB3(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2)
        {
            //C***BEGIN PROLOGUE  PASSB3
            //C***REFER TO  CFFTB
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSB3
            CC.array(IDO,3,L1);
            CH.array(IDO,L1,3);
            const double TAUR = -0.5;
            const double TAUI = 0.866025403784439;
            int K,I;
            double TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3;
            //C***FIRST EXECUTABLE STATEMENT  PASSB3
            if (IDO != 2) goto L102;
            for (K=1; K<=L1; ++K) {
                TR2 = CC(1,2,K)+CC(1,3,K);
                CR2 = CC(1,1,K)+TAUR*TR2;
                CH(1,K,1) = CC(1,1,K)+TR2;
                TI2 = CC(2,2,K)+CC(2,3,K);
                CI2 = CC(2,1,K)+TAUR*TI2;
                CH(2,K,1) = CC(2,1,K)+TI2;
                CR3 = TAUI*(CC(1,2,K)-CC(1,3,K));
                CI3 = TAUI*(CC(2,2,K)-CC(2,3,K));
                CH(1,K,2) = CR2-CI3;
                CH(1,K,3) = CR2+CI3;
                CH(2,K,2) = CI2+CR3;
                CH(2,K,3) = CI2-CR3;
            }
            return;
L102:
            if (IDO/2 < L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    TR2 = CC(I-1,2,K)+CC(I-1,3,K);
                    CR2 = CC(I-1,1,K)+TAUR*TR2;
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2;
                    TI2 = CC(I,2,K)+CC(I,3,K);
                    CI2 = CC(I,1,K)+TAUR*TI2;
                    CH(I,K,1) = CC(I,1,K)+TI2;
                    CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K));
                    CI3 = TAUI*(CC(I,2,K)-CC(I,3,K));
                    DR2 = CR2-CI3;
                    DR3 = CR2+CI3;
                    DI2 = CI2+CR3;
                    DI3 = CI2-CR3;
                    CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2;
                    CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2;
                    CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3;
                    CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    TR2 = CC(I-1,2,K)+CC(I-1,3,K);
                    CR2 = CC(I-1,1,K)+TAUR*TR2;
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2;
                    TI2 = CC(I,2,K)+CC(I,3,K);
                    CI2 = CC(I,1,K)+TAUR*TI2;
                    CH(I,K,1) = CC(I,1,K)+TI2;
                    CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K));
                    CI3 = TAUI*(CC(I,2,K)-CC(I,3,K));
                    DR2 = CR2-CI3;
                    DR3 = CR2+CI3;
                    DI2 = CI2+CR3;
                    DI3 = CI2-CR3;
                    CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2;
                    CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2;
                    CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3;
                    CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3;
                }
            }
            return;
        }

        void PASSB4(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3)
        {
            //C***BEGIN PROLOGUE  PASSB4
            //C***REFER TO  CFFTB
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSB4
            CC.array(IDO,4,L1);
            CH.array(IDO,L1,4);
            int K,I;
            double TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4;
            //C***FIRST EXECUTABLE STATEMENT  PASSB4
            if (IDO != 2) goto L102;
            for (K=1; K<=L1; ++K) {
                TI1 = CC(2,1,K)-CC(2,3,K);
                TI2 = CC(2,1,K)+CC(2,3,K);
                TR4 = CC(2,4,K)-CC(2,2,K);
                TI3 = CC(2,2,K)+CC(2,4,K);
                TR1 = CC(1,1,K)-CC(1,3,K);
                TR2 = CC(1,1,K)+CC(1,3,K);
                TI4 = CC(1,2,K)-CC(1,4,K);
                TR3 = CC(1,2,K)+CC(1,4,K);
                CH(1,K,1) = TR2+TR3;
                CH(1,K,3) = TR2-TR3;
                CH(2,K,1) = TI2+TI3;
                CH(2,K,3) = TI2-TI3;
                CH(1,K,2) = TR1+TR4;
                CH(1,K,4) = TR1-TR4;
                CH(2,K,2) = TI1+TI4;
                CH(2,K,4) = TI1-TI4;
            }
            return;
L102:
            if (IDO/2 < L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    TI1 = CC(I,1,K)-CC(I,3,K);
                    TI2 = CC(I,1,K)+CC(I,3,K);
                    TI3 = CC(I,2,K)+CC(I,4,K);
                    TR4 = CC(I,4,K)-CC(I,2,K);
                    TR1 = CC(I-1,1,K)-CC(I-1,3,K);
                    TR2 = CC(I-1,1,K)+CC(I-1,3,K);
                    TI4 = CC(I-1,2,K)-CC(I-1,4,K);
                    TR3 = CC(I-1,2,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = TR2+TR3;
                    CR3 = TR2-TR3;
                    CH(I,K,1) = TI2+TI3;
                    CI3 = TI2-TI3;
                    CR2 = TR1+TR4;
                    CR4 = TR1-TR4;
                    CI2 = TI1+TI4;
                    CI4 = TI1-TI4;
                    CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2;
                    CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2;
                    CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3;
                    CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3;
                    CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4;
                    CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    TI1 = CC(I,1,K)-CC(I,3,K);
                    TI2 = CC(I,1,K)+CC(I,3,K);
                    TI3 = CC(I,2,K)+CC(I,4,K);
                    TR4 = CC(I,4,K)-CC(I,2,K);
                    TR1 = CC(I-1,1,K)-CC(I-1,3,K);
                    TR2 = CC(I-1,1,K)+CC(I-1,3,K);
                    TI4 = CC(I-1,2,K)-CC(I-1,4,K);
                    TR3 = CC(I-1,2,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = TR2+TR3;
                    CR3 = TR2-TR3;
                    CH(I,K,1) = TI2+TI3;
                    CI3 = TI2-TI3;
                    CR2 = TR1+TR4;
                    CR4 = TR1-TR4;
                    CI2 = TI1+TI4;
                    CI4 = TI1-TI4;
                    CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2;
                    CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2;
                    CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3;
                    CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3;
                    CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4;
                    CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4;
                }
            }
            return;
        }

        void PASSB5(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3,DFARRAY WA4)
        {
            //C***BEGIN PROLOGUE  PASSB5
            //C***REFER TO  CFFTB
            //C***ROUTINES CALLED  (NONE)
            //C***REVISION HISTORY  (YYMMDD)
            //C   000330  Modified array declarations.  (JEC)
            //C
            //C***END PROLOGUE  PASSB5
            CC.array(IDO,5,L1);
            CH.array(IDO,L1,5);
            const double TR11 = 0.309016994374947;
            const double TI11 = 0.951056516295154;
            const double TR12 = -0.809016994374947;
            const double TI12 = 0.587785252292473;
            int K,I;
            double TI5,TI2,TI4,TI3,TR5,TR2,TR4,TR3,CR2,CI2,CR3,CI3,CR5,CI5,CR4,CI4;
            double DR3,DR4,DI3,DI4,DR5,DR2,DI5,DI2;
            //C***FIRST EXECUTABLE STATEMENT  PASSB5
            if (IDO != 2) goto L102;
            for (K=1; K<=L1; ++K) {
                TI5 = CC(2,2,K)-CC(2,5,K);
                TI2 = CC(2,2,K)+CC(2,5,K);
                TI4 = CC(2,3,K)-CC(2,4,K);
                TI3 = CC(2,3,K)+CC(2,4,K);
                TR5 = CC(1,2,K)-CC(1,5,K);
                TR2 = CC(1,2,K)+CC(1,5,K);
                TR4 = CC(1,3,K)-CC(1,4,K);
                TR3 = CC(1,3,K)+CC(1,4,K);
                CH(1,K,1) = CC(1,1,K)+TR2+TR3;
                CH(2,K,1) = CC(2,1,K)+TI2+TI3;
                CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3;
                CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3;
                CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3;
                CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3;
                CR5 = TI11*TR5+TI12*TR4;
                CI5 = TI11*TI5+TI12*TI4;
                CR4 = TI12*TR5-TI11*TR4;
                CI4 = TI12*TI5-TI11*TI4;
                CH(1,K,2) = CR2-CI5;
                CH(1,K,5) = CR2+CI5;
                CH(2,K,2) = CI2+CR5;
                CH(2,K,3) = CI3+CR4;
                CH(1,K,3) = CR3-CI4;
                CH(1,K,4) = CR3+CI4;
                CH(2,K,4) = CI3-CR4;
                CH(2,K,5) = CI2-CR5;
            }
            return;
L102:
            if (IDO/2 < L1) goto L105;
            for (K=1; K<=L1; ++K) {
                #pragma omp parallel for
                for (I=2; I<=IDO; I+=2) {
                    TI5 = CC(I,2,K)-CC(I,5,K);
                    TI2 = CC(I,2,K)+CC(I,5,K);
                    TI4 = CC(I,3,K)-CC(I,4,K);
                    TI3 = CC(I,3,K)+CC(I,4,K);
                    TR5 = CC(I-1,2,K)-CC(I-1,5,K);
                    TR2 = CC(I-1,2,K)+CC(I-1,5,K);
                    TR4 = CC(I-1,3,K)-CC(I-1,4,K);
                    TR3 = CC(I-1,3,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3;
                    CH(I,K,1) = CC(I,1,K)+TI2+TI3;
                    CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3;
                    CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3;
                    CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3;
                    CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3;
                    CR5 = TI11*TR5+TI12*TR4;
                    CI5 = TI11*TI5+TI12*TI4;
                    CR4 = TI12*TR5-TI11*TR4;
                    CI4 = TI12*TI5-TI11*TI4;
                    DR3 = CR3-CI4;
                    DR4 = CR3+CI4;
                    DI3 = CI3+CR4;
                    DI4 = CI3-CR4;
                    DR5 = CR2+CI5;
                    DR2 = CR2-CI5;
                    DI5 = CI2-CR5;
                    DI2 = CI2+CR5;
                    CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2;
                    CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2;
                    CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3;
                    CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3;
                    CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4;
                    CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4;
                    CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5;
                    CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5;
                }
            }
            return;
L105:
            for (I=2; I<=IDO; I+=2) {
                #pragma omp parallel for
                for (K=1; K<=L1; ++K) {
                    TI5 = CC(I,2,K)-CC(I,5,K);
                    TI2 = CC(I,2,K)+CC(I,5,K);
                    TI4 = CC(I,3,K)-CC(I,4,K);
                    TI3 = CC(I,3,K)+CC(I,4,K);
                    TR5 = CC(I-1,2,K)-CC(I-1,5,K);
                    TR2 = CC(I-1,2,K)+CC(I-1,5,K);
                    TR4 = CC(I-1,3,K)-CC(I-1,4,K);
                    TR3 = CC(I-1,3,K)+CC(I-1,4,K);
                    CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3;
                    CH(I,K,1) = CC(I,1,K)+TI2+TI3;
                    CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3;
                    CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3;
                    CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3;
                    CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3;
                    CR5 = TI11*TR5+TI12*TR4;
                    CI5 = TI11*TI5+TI12*TI4;
                    CR4 = TI12*TR5-TI11*TR4;
                    CI4 = TI12*TI5-TI11*TI4;
                    DR3 = CR3-CI4;
                    DR4 = CR3+CI4;
                    DI3 = CI3+CR4;
                    DI4 = CI3-CR4;
                    DR5 = CR2+CI5;
                    DR2 = CR2-CI5;
                    DI5 = CI2-CR5;
                    DI2 = CI2+CR5;
                    CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2;
                    CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2;
                    CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3;
                    CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3;
                    CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4;
                    CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4;
                    CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5;
                    CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5;
                }
            }
            return;
        }

    } // namespace CMLIB
} // namespace SCATMECH




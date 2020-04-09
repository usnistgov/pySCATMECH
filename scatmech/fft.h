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
#ifndef SCATMECH_FFT_H
#define SCATMECH_FFT_H

#include "scatmech.h"
#include "matrixmath.h"
#include <map>
#include <valarray>

#include <iostream>

namespace SCATMECH {

    namespace CMLIB {

        void CFFTI(int& N,DFARRAY WSAVE);
        void CFFTI1(int& N,DFARRAY WA,DFARRAY IFAC);

        void CFFTF(int& N,DFARRAY C,DFARRAY WSAVE);
        void CFFTF1(int& N,DFARRAY C,DFARRAY CH,DFARRAY WA,DFARRAY IFAC);
        void PASSF(int& NAC,int& IDO,int& IP,int& L1,int& IDL1,DFARRAY CC,DFARRAY C1,DFARRAY C2,DFARRAY CH,DFARRAY CH2,DFARRAY WA);
        void PASSF2(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1);
        void PASSF3(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2);
        void PASSF4(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3);
        void PASSF5(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3,DFARRAY WA4);

        void CFFTB(int& N,DFARRAY C,DFARRAY WSAVE);
        void CFFTB1(int& N,DFARRAY C,DFARRAY CH,DFARRAY WA,DFARRAY IFAC);
        void PASSB(int& NAC,int& IDO,int& IP,int& L1,int& IDL1,DFARRAY CC,DFARRAY C1,DFARRAY C2,DFARRAY CH,DFARRAY CH2,DFARRAY WA);
        void PASSB2(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1);
        void PASSB3(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2);
        void PASSB4(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3);
        void PASSB5(int& IDO,int& L1,DFARRAY CC,DFARRAY CH,DFARRAY WA1,DFARRAY WA2,DFARRAY WA3,DFARRAY WA4);

    } // CMLIB

    void fft1d(std::valarray<COMPLEX>& data,int isign);
    void fft1d(CFARRAY data,int N,int isign);

} // SCATMECH

#endif

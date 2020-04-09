//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: scatmatrix.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_SCATMATRIX_H
#define SCATMECH_SCATMATRIX_H

#include "scatmech.h"

namespace SCATMECH {

    class ScatterTMatrix {
        public:
            ScatterTMatrix(int lmax=0);
            ScatterTMatrix(const ScatterTMatrix& matrix2);
            ScatterTMatrix& operator=(const ScatterTMatrix& matrix2);
            ~ScatterTMatrix();

            void resize(int lmax);
            void free();

            COMPLEX** operator[](int m) {
                return matrix[m];
            }
        private:
            COMPLEX ***matrix;
            int LMAX;
    };

}

#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: scatmatrix.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include <cstdlib>
#include "scatmatrix.h"

using namespace std;

namespace SCATMECH {

    static
    int
    index(int l,int m, int f) {
        return 2*(l*l+l-1+m)+f;
    }

    ScatterTMatrix::
    ScatterTMatrix(int lmax)
    {
        matrix=0;
        LMAX=0;
        resize(lmax);
    }

    ScatterTMatrix::
    ScatterTMatrix(const ScatterTMatrix& matrix2) {
        matrix=0;
        LMAX=0;
        *this = matrix2;
    }

    ScatterTMatrix&
    ScatterTMatrix::
    operator=(const ScatterTMatrix& matrix2) {
        if (matrix2.LMAX!=LMAX) {
            resize(matrix2.LMAX);
        }
        for (int m=-LMAX; m<=LMAX; ++m) {
            int beginl = (m==0) ? 1 : abs(m);
            for (int l=beginl; l<=LMAX; ++l) {
                for (int f=0; f<=1; ++f) {
                    int i = index(l,m,f);
                    for (int l_=beginl; l_<=LMAX; ++l_) {
                        for (int f_=0; f_<=1; ++f_) {
                            int i_ = index(l_,m,f_);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            matrix[mm][ll_][ll]=matrix2.matrix[mm][ll_][ll];
                        }
                    }
                }
            }
        }
        return *this;
    }

    void
    ScatterTMatrix::
    free()
    {
        if (matrix) {
            for (int m=-LMAX; m<=LMAX; ++m) {
                int n= 2*((m==0) ? LMAX : LMAX-abs(m)+1);
                for (int l=0; l<n; ++l) {
                    if (matrix[m+LMAX][l]) delete[] matrix[m+LMAX][l];
                }
                if (matrix[m+LMAX]) delete[] matrix[m+LMAX];
            }
            delete[] matrix;
        }
        LMAX=0;
    }

    ScatterTMatrix::
    ~ScatterTMatrix()
    {
        free();
    }

    void
    ScatterTMatrix::
    resize(int lmax)
    {
        try {
            if (lmax!=LMAX)    free();
            if (lmax!=0) {
                matrix = new COMPLEX**[2*lmax+1];
                for (int m=-lmax; m<=lmax; ++m) {
                    int n= 2*((m==0) ? lmax : lmax-abs(m)+1);
                    matrix[lmax+m] = new COMPLEX*[n];
                    for (int j=0; j<n; ++j) matrix[lmax+m][j] = new COMPLEX[n];
                }
            } else {
                matrix=0;
            }

            LMAX = lmax;

        } catch (bad_alloc) {
            matrix=0;
            throw SCATMECH_exception("Error allocating memory in ScatterTMatrix::resize()");
        }
    }
}

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: tmatrix.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "axipart.h"
#include <fstream>
#include <algorithm>

using namespace std;

namespace SCATMECH {

    using namespace BobVlieg_Supp;

    static void
    solve(vector<vector<COMPLEX> >& a,vector<vector<COMPLEX> >& b,int n)
    {
        //c     ...............................................................
        //c     .                      -1                                     .
        //c     .  calculate [T]' = [A]' * [B]' using Gauss-Jordan reduction  .
        //c     ...............................................................
        vector<int> ls(n);
        //c     ............................
        //c     .  start reduction of [a]  .
        //c     ............................
        for (int i=0; i<n; ++i) {
            //c     ..........................................................
            //c     .  search for the maximum element in the ith row of [a]  .
            //c     ..........................................................
            COMPLEX aijmax = a[i][0]; // a(i,1)
            int jmax = 0;
            for (int j=1; j<n; ++j) {
                if(abs(a[i][j])>abs(aijmax)) {
                    aijmax = a[i][j];
                    jmax = j;
                }
            }
            //c     ....................................................
            //c     .  normalize the ith row of [a] and [b] by aijmax  .
            //c     ....................................................
            for (int jj=0; jj<n; ++jj) {
                a[i][jj] /= aijmax;
                b[i][jj] /= aijmax;
            }
            //c     .......................................................
            //c     .  use row transformations to obtain zeros above and  .
            //c     .    below the jmax element of the ith row of [a] -   .
            //c     .    apply the same row transformations to [b]        .
            //c     .......................................................
            for (int k=0; k<n; ++k) {
                if(k!=i) {
                    COMPLEX arat = -a[k][jmax];
                    for (int j=0; j<n; ++j) {
                        if(abs(a[i][j])>0.0) {
                            a[k][j] += arat*a[i][j];
                        }
                    }
                    a[k][jmax] = 0.0;
                    for (int jj=0; jj<n; ++jj) {
                        if(abs(b[i][jj])>0.0) {
                            b[k][jj] += arat*b[i][jj];
                        }
                    }
                }
                //c     ........................................................
                //c     .  store row counter (i) in array ls(*), such that     .
                //c     .    ls(*) contains the location of the pivot (unity)  .
                //c     .    element of each column (after reduction)          .
                //c     ........................................................
                ls[jmax] = i;
            }
        }
        //c     ....................................................
        //c     .  the reduction of [a] is complete - perform      .
        //c     .    row interchanges as indicated in array ls(*)  .
        //c     ....................................................
        for (int ii=0; ii<n; ++ii) {
            int k = ii;
            //c     .........................................................
            //c     .  put the integer value in ls(k) into k                .
            //c     .                                                       .
            //c     .  if k is less than i, then that row has already been  .
            //c     .    involved in an interchange so iterate k = ls(k)    .
            //c     .    until a value of k greater than i (corresponding   .
            //c     .    to a row stored above the ith row) is obtained     .
            //c     .........................................................
            do {
                k = ls[k];
            } while (k<ii);
            if (k>0) {
                for (int j=0; j<n; ++j) {
                    swap(b[ii][j],b[k][j]);
                }
            }
        }
        //c     ........................................
        //c     .  [T]' is stored in the b(*,*) array  .
        //c     ........................................
        return;
    }

    void
    Axisymmetric_Shape::
    Get_TMatrix(ScatterTMatrix& T, int nrank, int mmax, double k, COMPLEX& index)
    {
        //         This routine was adapted from "Light Scattering by Particles:
        //         Computational Methods," by P.W. Barber and S.C. Hill, pp. 164-183.
        //         Copyright (c) 1990 by World Scientific Publishing Co Pte Ltd.
        //         See: http://www.worldscibooks.com/engineering/0784.html
        //         Permission to incorporate this work into SCATMECH was granted.

        //           The code was modified from FORTRAN to C++, and several changes
        //         were made:
        //              (1) Dynamic allocation of memory as needed
        //              (2) Arbitrary axisymmetric shapes can be handled through
        //                  the elements of a shape_vec.  The Gauss-Legendre integration
        //                  is effectively handled by the Axisymmetric_Shape class, rather
        //                  than by this routine.
        //              (3) No convergence testing is performed.
        //              (4) The simplification that Barber and Hill incorporated for
        //                  mirror symmetric particles was removed.
        //
        //c     ..................................................................
        //c     .  Light Scattering by Particles: Computational Methods          .
        //c     .  by P.W. Barber and S.C. Hill                                  .
        //c     .  copyright (c) 1990 by World Scientific Publishing Co Pte Ltd  .
        //c     .                                                                .
        //c     .  equation numbers in columns 73-80 are references to the text  .
        //c     ..................................................................
        //c     .................................................................
        //c     .  calculate the scattering by axisymmetric dielectric objects  .
        //c     .  convergence test over three parameters:                      .
        //c     .    nrank : number of terms (matrix order)                     .
        //c     .    ntheta: number of integration steps                        .
        //c     .    nm    : number of azimuthal modes m                        .
        //c     .  inputs:  nrank  = matrix order                               .
        //c     .           ntheta = integration steps                          .
        //c     .           ic = convergence case                               .
        //c     .                ntheta ( = 0 ), nrank ( = 1 ), or nm ( = 2 )   .
        //c     .           x = size parameter (ka)                             .
        //c     .           aovrb = a/b ratio                                   .
        //c     .           cm = complex index of refraction, (real,imag)       .
        //c     .                (imag is positive for absorption)              .
        //c     .                                                               .
        //c     .  calculate [T] = [B] * [A]   =  { [A]' * [B]'}'               .
        //c     .                                                               .
        //c     .    use transposed matrices to permit the efficient            .
        //c     .       -1                                      -1              .
        //c     .     [ ]  * [ ]  operation rather than [ ] * [ ]               .
        //c     .                                                               .
        //c     .    (1)  obtain the transposed matrices [A]' and [B]'          .
        //c     .                             -1                                .
        //c     .    (2)  calculate [T]' = [A]' * [B]'                          .
        //c     .                                                               .
        //c     .    (3)  transpose [T]' to obtain [T]                          .
        //c     .                                                               .
        //c     .  two simplifications are incorporated in the program to       .
        //c     .    speed up the computation:                                  .
        //c     .                                                               .
        //c     .  (1)  [Simplification has been removed for the SCATMECH  .
        //c     .        version.]                                              .
        //c     .                                                               .
        //c     .  (2)  all matrix elements and incident and scattered          .
        //c     .       field coefficients are zero when n < m by virtue        .
        //c     .       of the behavior of the associated Legendre              .
        //c     .       functions.  all matrices and coefficient arrays         .
        //c     .       are compressed to take advantage of this behavior       .
        //c     .................................................................

        vector<double> pnmllg;
        vector<COMPLEX> hankel;
        vector<COMPLEX> bslcmp;

        int isize=2*(nrank+1);

        vector<vector<COMPLEX> > a(isize,vector<COMPLEX>(isize));
        vector<vector<COMPLEX> > b(isize,vector<COMPLEX>(isize));

        //c     ...........................
        //c     .  set program constants  .
        //c     ...........................
        double pi2 = pi/2.0;

        int nranki = nrank+1;
        int nr2 = 2*nrank;

        //c     .........................................
        //c     .  set the complex index of refraction  .
        //c     .    for an exp(-iwt) time variation    .
        //c     .........................................
        COMPLEX cm = index;
        COMPLEX dcn = sqr(cm);
        //c     .................................................................
        //c     .  convergence over m (ic = 2) requires nrank azimuthal modes,  .
        //c     .    a nonsymmetric orientation, and a file for the matrices    .
        //c     .................................................................

        int nm = (mmax<nrank) ? mmax : nrank;

        //c     .....................................................
        //c     .  set the integration points and weighting values  .
        //c     .....................................................

        shape_vec shape2 = shape;

        int ntheta = shape2.size();

        for (int i=0; i<ntheta; ++i) {
            shape2[i].shape *= k;
            shape2[i].dshape *= k;
        }

        //c     ............................................
        //c     .  enter a loop for each azimuthal mode m  .
        //c     ............................................
        for (int im = 1; im<=nm+1; ++im) {
            ostringstream msg;
            msg << "Calculating T matrix: m = " << im-1 << ret;
            message(msg.str());

            //c     ...............................
            //c     .  set m-dependent variables  .
            //c     ...............................
            int kmv = im-1;
            if (nm==1) kmv = 1;

            double cmv = kmv;
            double cm2 = sqr(cmv);
            double prodm = 1.0;

            double em;
            if(kmv==0) {
                em = 1.0;
            } else {
                em = 2.0;
                double quanm = cmv;
                for (int i=1; i<=kmv; ++i) {
                    quanm = quanm+1.0;
                    prodm = quanm*prodm/2.0;
                }
            }
            double qem = 2.0/em;
            double twm = 2.0*cmv;
            //c     ...................................................
            //c     .  set indices for matrix compression when n < m  .
            //c     .    note: useful only when m > 1                 .
            //c     ...................................................
            int ij = kmv-1;
            if (ij<0) ij = 0;
            int ijt = 2*ij;
            int ns = nrank-ij;
            int ns2 = 2*ns;
            //c     ............................................
            //c     .  initialize all matrix elements to zero  .
            //c     ............................................
            for (int i=0; i<isize; ++i) {
                for (int j=0; j<isize; ++j) {
                    a[i][j] = b[i][j] = 0.;
                }
            }
            //c     ................................................
            //c     .  enter a loop to integrate over the surface  .
            //c     .    (theta is the integration variable)       .
            //c     ................................................
            for (int ithta=0; ithta<ntheta; ++ithta) {

                double theta = shape2[ithta].abs;
                double sinth = sin(theta);
                double costh = cos(theta);

                //c     .................................................
                //c     .  calculate the associated Legendre functions  .
                //c     .    at each integration point                  .
                //c     .................................................
                pnmllg.resize(nranki+1);

                for (int n=0; n<=nranki; ++n) {
                    if (n>=kmv) {
                        pnmllg[n]=real(Legendre(n,kmv,costh))/sinth;
                    } else {
                        pnmllg[n]=0.;
                    }
                }
                pnmllg[nranki]=0.;

                //c     ......................................
                //c     .  calculate kr and it's derivative  .
                //c     .    at each integration point       .
                //c     ......................................
                double ckr = shape2[ithta].shape;
                double dckr = shape2[ithta].dshape;
                //c     ........................................................
                //c     .  calculate the Hankel functions (real argument) and  .
                //c     .    Bessel functions (complex argument) at each       .
                //c     .    integration point                                 .
                //c     ........................................................
                COMPLEX ckpr = cm*ckr;
                hankel.resize(nranki+1);
                bslcmp.resize(nranki+1);
                for (int ii=0; ii<nranki; ++ii) {
                    hankel[ii] = h(ii,ckr);
                    bslcmp[ii] = j(ii,ckpr);
                }
                hankel[nranki]=bslcmp[nranki]=0.;

                double d = dckr*sinth;
                double wtsin = shape2[ithta].wt*sinth;
                //c     .............................................
                //c     .  enter a loop for each row of the matrix  .
                //c     .............................................
                for (int irow=1; irow<=nrank; ++irow) {

                    if(irow>ij) {
                        int irow1 = irow+nrank;
                        double crow = irow;
                        double crowm = cmv+crow;
                        double crow1 = crow+1.0;
                        double br = real(hankel[irow-1])/real(hankel[irow+1-1]);
                        COMPLEX h = hankel[irow-1]/hankel[irow+1-1];
                        //c     ................................................
                        //c     .  enter a loop for each column of the matrix  .
                        //c     ................................................
                        for (int icol=1; icol<=nrank; ++icol) {

                            if (icol>ij) {

                                int icol1 = icol+nrank;

                                double ccol = icol;
                                double ccolm = cmv+ccol;
                                double ccol1 = ccol+1.0;
                                //c     .....................................
                                //c     .  calculate variable combinations  .
                                //c     .....................................
                                double crij = crow+ccol;
                                double crssij = crow*ccol;
                                double cmcrco = cm2+qem*crssij*sqr(costh);
                                double pnr0c0 = pnmllg[irow-1]*pnmllg[icol-1];
                                double pnr0c1 = pnmllg[irow-1]*pnmllg[icol+1-1];
                                double pnr1c0 = pnmllg[irow+1-1]*pnmllg[icol-1];
                                double pnr1c1 = pnmllg[irow+1-1]*pnmllg[icol+1-1];
                                double b1a = crow*costh*pnr1c1-crowm*pnr0c1;
                                double b1b = ccol*costh*pnr1c1-ccolm*pnr1c0;
                                COMPLEX bk = cm*bslcmp[icol-1]/bslcmp[icol+1-1];
                                COMPLEX bbk = real(hankel[irow+1-1])*bslcmp[icol+1-1]*wtsin;
                                COMPLEX hbk = hankel[irow+1-1]*bslcmp[icol+1-1]*wtsin;
                                //c     .....................................................
                                //c     .  calculate the [A]' and  [B]' submatrices, e.g.,  .
                                //c     .            _                             _        .
                                //c     .           |               |               |       .
                                //c     .           |             ,               , |       .
                                //c     .           | (K + cm * J)  |-(I + cm * L)  |       .
                                //c     .     ,     |                               |       .
                                //c     .  [A]  =   |- - - - - - - -|- - - - - - - -|       .
                                //c     .           |             ,               , |       .
                                //c     .           | (L + cm * I)  | (J + cm * K)  |       .
                                //c     .           |                               |       .
                                //c     .           |_              |              _|       .
                                //c     .                                                   .
                                //c     .....................................................

                                //c     ....................................................
                                //c     .  increment the -(I + cm * L)' submatrix element  .
                                //c     ....................................................
                                if(kmv != 0) {
                                    double b1 = b1a+b1b;
                                    COMPLEX sa = pnr1c1*
                                                 (crow*crow1*bk+ccol*ccol1*h-crssij*(crij+2.0)/ckr)*d;
                                    sa = (ckr*(1.0+h*bk)-ccol*h-crow*bk+crssij/ckr)*b1*ckr+sa;
                                    // eq 3.23 ...
                                    a[icol-ij-1][irow1-ijt-1] += cmv*sa*hbk;
                                    COMPLEX sb = sa+(br-h)*(ccol*ccol1*pnr1c1*d+(bk*ckr-ccol)*b1*ckr);
                                    b[icol-ij-1][irow1-ijt-1] += +cmv*sb*bbk;
                                    //c     ...................................................
                                    //c     .  increment the (L + cm * I)' submatrix element  .
                                    //c     ...................................................
                                    COMPLEX c = (dcn-1.0)*b1*sqr(ckr);
                                    // eq 3.24 ...
                                    a[icol1-ijt-1][irow-ij-1] -= cmv*(sa+c)*hbk/cm;
                                    b[icol1-ijt-1][irow-ij-1] -= cmv*(sb+c)*bbk/cm;
                                }
                                //c     ...................................................
                                //c     .  increment the (J + cm * K)' submatrix element  .
                                //c     ...................................................
                                double a12 = cmcrco*pnr1c1-qem*(crow*ccolm*costh*pnr1c0+
                                                                ccol*crowm*costh*pnr0c1-crowm*ccolm*pnr0c0);
                                double a22 = a12*sqr(ckr);
                                b1a = ccol*ccol1*b1a;
                                b1b = crow*crow1*b1b;
                                COMPLEX sa = (ckr*(bk-dcn*h)+dcn*crow-ccol)*a12*ckr+(b1a-dcn*b1b)*qem*d;
                                //  eq 3.25 ...
                                a[icol1-ijt-1][irow1-ijt-1] += sa*hbk/cm;
                                b[icol1-ijt-1][irow1-ijt-1] += (sa-(br-h)*dcn*a22)*bbk/cm;
                                //c     ...................................................
                                //c     .  increment the (K + cm * J)' submatrix element  .
                                //c     ...................................................
                                sa = (ckr*(bk-h)+crow-ccol)*a12*ckr+(b1a-b1b)*qem*d;
                                //  eq 3.26 ...
                                a[icol-ij-1][irow-ij-1] += sa*hbk;
                                b[icol-ij-1][irow-ij-1] += (sa-(br-h)*a22)*bbk;

                            } // matches if (icol>ij)
                        } // matches for (int icol=1;icol<=nrank;++icol)
                    } // matches if(irow>ij)
                } // matches for (int irow=1;irow<=nrank;++irow)
            } // matches for (int ithta=0;ithta<ntheta;++ithta)
            //c     ..........................................
            //c     .                      -1                .
            //c     .  calculate [T]' = [A]' * [B]'          .
            //c     ..........................................
            solve(a,b,ns2);
            //c     ........................................
            //c     .  [T]' is stored in the b(*,*) array  .
            //c     .    transpose [T]' to obtain [T]      .
            //c     ........................................
            for (int ir=0; ir<ns2; ++ir) {
                for (int jr=0; jr<ns2; ++jr) {
                    a[ir][jr]=b[jr][ir];
                }
            }

            //c     ..................................................
            //c     .  for convergence over nm (ic = 2) save the     .
            //c     .    [T] matrix and normalization constants      .
            //c     .    for later use in scattering programs        .
            //c     .  the [T] matrix is stored in the a(*,*) array  .
            //c     ..................................................

            int lmax=nrank;
            int m = kmv;

            //
            // The following copies the matrix a (for the mth azimuthal order)
            // into the the returned T-matrix.
            //
            int mm = lmax+m;
            int mm2 = lmax-m;
            int beginl = (m==0) ? 1 : abs(m);
            for (int l=beginl; l<=lmax; ++l) {
                for (int f=0; f<=1; ++f) {
                    int ll = 2*(lmax-l)+f;
                    for (int l_=beginl; l_<=lmax; ++l_) {
                        for (int f_=0; f_<=1; ++f_) {
                            int ll_ = 2*(lmax-l_)+f_;

                            int nn = ((f==0) ? ns : 0 ) + (l-beginl);
                            int nn_ = ((f_==0) ? ns : 0 ) + (l_-beginl);

                            T[mm][ll_][ll]=-a[nn_][nn];
                            T[mm2][ll_][ll]= (f==f_) ? -a[nn_][nn] : a[nn_][nn];
                        }
                    }
                }
            }
        } // matches for (int im = 1; im<=nm+1; ++im)
    }

} // namespace SCATMECH

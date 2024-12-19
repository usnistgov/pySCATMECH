//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: axipart1.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
//#include <float.h>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "scatmech.h"
#include "bobvlieg.h"
#include "axipart.h"
#include "dielfunc.h"
#include "askuser.h"

using namespace std;

namespace SCATMECH {

    using namespace BobVlieg_Supp;

    // The following code implements the theory described in
    //   P.A. Bobbert and J. Vlieger, Physica 137A (1986) 209-242.
    // Some details of the theory can be found also in
    //   P.A. Bobbert, J. Vlieger, and  R. Grief, Physica 137A (1986) 243-257
    // These articles will be referred to as B&V and BV&G, respectively.
    //

    static complex<float> cfI(0.,1.);

    static
    COMPLEX
    D(int n,int m,int f)
    {
        int mm = (m==0 ? 1 : 2);

        double d = mm*(2*n+1)*Fact(n-abs(m))/(4*n*(n+1)*Fact(n+abs(m)));

        double e = sqrt((2*n+1)*Fact(n-abs(m))/Fact(n+abs(m)));

        COMPLEX ff = (f==0) ? 1. : -1.;

        return 1./(ff*(double)(e/d));
    }

    //
    // Bmatrix() calculates the particle's scattering matrix found
    // in Eq. (4.2) of B&V.
    //
    void
    Axisymmetric_Particle_BRDF_Model::
    Bmatrix(ScatterTMatrix& T)
    {
        //float xdummy;
        //long ldummy;
        complex<float> dummy;

        int m;
        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (int l=beginl; l<=LMAX; ++l) {
                for (int f=0; f<=1; ++f) {
                    //int i = index(l,m,f);
                    int ll = 2*(LMAX-l)+f;
                    for (int l_=beginl; l_<=LMAX; ++l_) {
                        for (int f_=0; f_<=1; ++f_) {
                            //int i_ = index(l_,m,f_);
                            int ll_ = 2*(LMAX-l_)+f_;
                            T[mm][ll_][ll]=0;
                        }
                    }
                }
            }
        }

        Shape->Get_TMatrix(T, LMAX, MMAX, k, N2);

        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (int l=beginl; l<=LMAX; ++l) {
                for (int f=0; f<=1; ++f) {
                    int ll = 2*(LMAX-l)+f;
                    for (int l_=beginl; l_<=LMAX; ++l_) {
                        for (int f_=0; f_<=1; ++f_) {
                            int ll_ = 2*(LMAX-l_)+f_;
                            T[mm][ll_][ll]=T[mm][ll_][ll]*D(l_,m,f_)/D(l,m,f);
                        }
                    }
                }
            }
        }


    }


    //
    // Amatrix() calculates the matrix A given by Eq. (8.17) of B&V
    //
    void
    Axisymmetric_Particle_BRDF_Model::
    Amatrix(ScatterTMatrix& A)
    {
        // Choose number of points in Gaussian integration, depending upon
        // LMAX, and round up to nearest 10.
        int npts=((LMAX/10)+1)*10;
        if (npts>100) npts=100;
        int ipt = npts/10-1;

        // Calculate the coefficients in front of the
        // Legendre polynomials found in Eqs. (8.5a) and (8.5b) of B&V.
        vector<COMPLEX> U1(index(LMAX,LMAX)+1);
        vector<COMPLEX> U2(index(LMAX,LMAX)+1);
        vector<COMPLEX> U3(index(LMAX,LMAX)+1);
        vector<COMPLEX> U4(index(LMAX,LMAX)+1);
        vector<COMPLEX> V1(index(LMAX,LMAX)+1);
        vector<COMPLEX> V2(index(LMAX,LMAX)+1);
        for (int l=1; l<=LMAX; ++l) {
            for (int m=0; m<=l; ++m) { // We don't need negative m values
                int i=index(l,m);

                double temp1=(2.*l-1.)*(2.*l+1.);
                double temp2=(2.*l+3.)*(2.*l+1.);

                U1[i] = (l-1.)/2.*ipow(-(l-1))*sqrt((l+m-1.)*(l+m)/temp1);
                U2[i] = (l-1.)/2.*ipow(-(l-1))*sqrt((l-m-1.)*(l-m)/temp1);
                U3[i] =-(l+2.)/2.*ipow(-(l+1))*sqrt((l-m+1.)*(l-m+2.)/temp2);
                U4[i] =-(l+2.)/2.*ipow(-(l+1))*sqrt((l+m+1.)*(l+m+2.)/temp2);

                V1[i] = 0.5*ipow(-(l-1))*sqrt((l-m+1.)*(l+m));
                V2[i] =-0.5*ipow(-(l-1))*sqrt((l+m+1.)*(l-m));
            }
        }

        // Normal incidence reflection coefficients...
        COMPLEX rsnormal = stack->rs12(0.,lambda,vacuum,substrate);
        COMPLEX rpnormal = stack->rp12(0.,lambda,vacuum,substrate);

        vector<COMPLEX> reflect_s(npts);
        vector<COMPLEX> reflect_p(npts);
        vector<vector<COMPLEX> > Umatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));
        vector<vector<COMPLEX> > Vmatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));
        vector<vector<COMPLEX> > dminusmatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));
        vector<vector<COMPLEX> > dplusmatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));

        for (int j=0; j<npts; ++j) {
            ostringstream msg;
            msg << "Creating lookup tables: j = " << j << '/' << npts
                << ret;
            message(msg.str());

            COMPLEX cosa = 1.0 - Gauss_Laguerre_Integration::zeros[ipt][j]/(2.*cI*qq);
            COMPLEX alpha = arccosine(cosa);
            COMPLEX cosa2 = cos((pi-alpha)/2.);
            COMPLEX sina2 = sin((pi-alpha)/2.);

            if (Norm_Inc_Approx) { // Normal incidence approximation
                reflect_s[j] = rsnormal;
                reflect_p[j] = rpnormal;
            } else {               // Exact solution
                reflect_s[j] = stack->rs12(alpha,lambda,vacuum,substrate);
                reflect_p[j] = stack->rp12(alpha,lambda,vacuum,substrate);
            }

            // Calculate U, V, d+, and d- for each angle and (l,m) combination...
            for (int l=1; l<=LMAX; ++l) {
                int _mmax = (MMAX>l) ? l : MMAX;
                for (int m=0; m<=_mmax; ++m) { // We don't need negative m values
                    int i = index(l,m);
                    Umatrix[j][i] = U1[i]*Ptilde(l-1,m-1,cosa)+
                                    U2[i]*Ptilde(l-1,m+1,cosa)+
                                    U3[i]*Ptilde(l+1,m-1,cosa)+
                                    U4[i]*Ptilde(l+1,m+1,cosa);
                    Vmatrix[j][i] = V1[i]*Ptilde(l,m-1,cosa)+
                                    V2[i]*Ptilde(l,m+1,cosa);

                    dminusmatrix[j][i] = dminus(l,m,cosa2,sina2);
                    dplusmatrix[j][i] = dplus(l,m,cosa2,sina2);
                }
            }
        }

        // Constant from change of integration variable (see Sec. 4 of BV&G)...
        COMPLEX expx = exp(2.*cI*qq)/(2.*cI*qq);

        // Construct A matrix...
        for (int m=0; m<=MMAX; ++m) {
            ostringstream msg;
            msg << "Constructing A matrix: m = " << m << '/' << MMAX << ret;
            message(msg.str());

            int beginl = (m==0) ? 1 : m;
            for (int l=beginl; l<=LMAX; ++l) {
                for (int l_=beginl; l_<=LMAX; ++l_) {
                    if (l>old_LMAX || l_>old_LMAX) {
                        COMPLEX Aee=0,Ahe=0,Aeh=0,Ahh=0;

                        int ii = index(l,m);
                        int ii_ = index(l_,m);

                        // Integration by Gaussian method [see Eq. (4.2) of BV&G]
                        for (int kk=0; kk<npts; ++kk) {
                            double weight = Gauss_Laguerre_Integration::weights[ipt][kk];
                            Aee += weight*(reflect_p[kk]*Vmatrix[kk][ii]*dminusmatrix[kk][ii_]
                                           +reflect_s[kk]*Umatrix[kk][ii]*dplusmatrix[kk][ii_]);
                            Ahe += weight*(reflect_p[kk]*Vmatrix[kk][ii]*dplusmatrix[kk][ii_]
                                           +reflect_s[kk]*Umatrix[kk][ii]*dminusmatrix[kk][ii_]);
                            Aeh += weight*(reflect_p[kk]*Umatrix[kk][ii]*dminusmatrix[kk][ii_]
                                           +reflect_s[kk]*Vmatrix[kk][ii]*dplusmatrix[kk][ii_]);
                            Ahh += weight*(reflect_p[kk]*Umatrix[kk][ii]*dplusmatrix[kk][ii_]
                                           +reflect_s[kk]*Vmatrix[kk][ii]*dminusmatrix[kk][ii_]);
                        }

                        // Constant value in front of Eq. (8.17) of BV...
                        COMPLEX aa = expx*mpow(m-1)*ipow(l_-1)*lvector[l_];

                        Aee *= aa;
                        Ahe *= -cI*aa;
                        Aeh *= cI*aa;
                        Ahh *= aa;

                        // Put the results into the A matrix...

                        A[LMAX+m][2*(LMAX-l_)][2*(LMAX-l)] = Aee;
                        A[LMAX+m][2*(LMAX-l_)+1][2*(LMAX-l)] = Ahe;
                        A[LMAX+m][2*(LMAX-l_)+1][2*(LMAX-l)+1] = Ahh;
                        A[LMAX+m][2*(LMAX-l_)][2*(LMAX-l)+1] = Aeh;

                        // Use the symmetry relations described in
                        // Eqs. (4.5) of BV&G...

                        A[LMAX-m][2*(LMAX-l_)][2*(LMAX-l)] = Aee;
                        A[LMAX-m][2*(LMAX-l_)+1][2*(LMAX-l)] = -Ahe;
                        A[LMAX-m][2*(LMAX-l_)+1][2*(LMAX-l)+1] = Ahh;
                        A[LMAX-m][2*(LMAX-l_)][2*(LMAX-l)+1] = -Aeh;
                    }
                }
            }
        }
    }

    //
    // Matrix 1-BA in Eq. 5.3 of B&V is block diagonal.
    // The following routine inverts that matrix.
    //
    void
    Axisymmetric_Particle_BRDF_Model::
    invert_block_diagonal(ScatterTMatrix& AA)
    {
        int m;
        for (m=-MMAX; m<=MMAX; ++m) {
            ostringstream msg;
            msg << "Inverting matrix (m = " << m << ") " << ret;
            message(msg.str());

            int n= 2*((m==0) ? LMAX : LMAX-abs(m)+1);

            cvert(AA[LMAX+m],n);
        }
        message("Done Inverting Matrices");
    }

    //
    // setup() is called whenever the model needs to be "recalculated"
    //
    void
    Axisymmetric_Particle_BRDF_Model::
    setup()
    {
        Local_BRDF_Model::setup();

        Shape->Write("shape.dat");
        double length = Shape->Get_Base_Length();

        ostringstream msg;
        msg << "Distance of origin from surface: " << length << endl;
        message(msg.str());

        k   = 2.0*pi/lambda;
        double q = k * length;
        qq  = q+k*delta;
        double qqq = Shape->Get_MaxRadius();
        // Bohren and Huffman recommended lmax...
        BH_LMAX = (int)(qqq*k+4.05*pow(qqq*k,0.3333)+2.);

        if (lmax==0) LMAX = BH_LMAX;
        if (lmax>0) LMAX = lmax;
        if (lmax<0) LMAX = BH_LMAX+abs(lmax);

        if (mmax==0) MMAX = LMAX;
        if (mmax>0) MMAX = mmax;
        if (mmax<0) MMAX = LMAX+abs(mmax);

        if (MMAX>LMAX) MMAX = LMAX;

        N0  = COMPLEX(substrate.index(lambda));
        N2 = COMPLEX(particle.index(lambda));

        ScatterTMatrix BMatrix(LMAX);

        Bmatrix(BMatrix);

        N2 = particle.index(lambda);
        N3 = n3.index(lambda);

        int m,l,l_,f,f_,l__,f__;

        if (delta<0.)
            error("delta cannot be less than zero.");


        // The Bohren and Huffman recommended LMAX...
        // BH_LMAX = (int)(q+4.*pow(q,0.3333)+2.);
        BH_LMAX = LMAX;

        sqrsize = index(LMAX,LMAX,1)+1;

        message("Allocating Memory");

        //Non-temporary-arrays...
        ScatMatrix.resize(LMAX);
        ScatMatrixInverse.resize(LMAX);

        Zp.resize(sqrsize);
        Zs.resize(sqrsize);
        eIP.resize(sqrsize);
        Wp.resize(sqrsize);
        Ws.resize(sqrsize);
        Vp.resize(sqrsize);
        Vs.resize(sqrsize);

        vector<COMPLEX> B(sqrsize);

        //
        // Calculate A matrix or read it from file and
        // write it to a file...
        //
        if (order!=0) {
            //int retvalue=1;
            old_LMAX=0;

            if (old_LMAX<LMAX) {
                message("Calculating A Matrix");
                Amatrix(ScatMatrix);
            }
        }

        //
        // Matrix A becomes B^(-1)*(1-B*A)...
        //

        message("Calculating matrix BA");
        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (l=beginl; l<=LMAX; ++l) {
                for (f=0; f<=1; ++f) {
                    //int i = index(l,m,f);
                    int ll = 2*(LMAX-l)+f;
                    for (l_=beginl; l_<=LMAX; ++l_) {
                        for (f_=0; f_<=1; ++f_) {
                            //int i_ = index(l_,m,f_);
                            int ll_ = 2*(LMAX-l_)+f_;
                            ScatMatrixInverse[mm][ll_][ll]=0.;
                            if (order!=0) {
                                for (l__=beginl; l__<=LMAX; ++l__) {
                                    for (f__=0; f__<=1; ++f__) {
                                        //int i__ = index(l__,m,f__);
                                        int ll__ = 2*(LMAX-l__)+f__;

                                        ScatMatrixInverse[mm][ll_][ll] += BMatrix[mm][ll_][ll__]*ScatMatrix[mm][ll__][ll];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        message("Calculating matrix 1-BA");

        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (l=beginl; l<=LMAX; ++l) {
                for (f=0; f<=1; ++f) {
                    //int i = index(l,m,f);
                    int ll = 2*(LMAX-l)+f;
                    for (l_=beginl; l_<=LMAX; ++l_) {
                        for (f_=0; f_<=1; ++f_) {
                            //int i_ = index(l_,m,f_);
                            int ll_ = 2*(LMAX-l_)+f_;
                            ScatMatrixInverse[mm][ll_][ll]=(double)(ll_==ll)-ScatMatrixInverse[mm][ll_][ll];
                        }
                    }
                }
            }
        }

        message("Inverting B");
        invert_block_diagonal(BMatrix);

        message("Creating B^(-1)(1-BA)");
        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (l=beginl; l<=LMAX; ++l) {
                for (f=0; f<=1; ++f) {
                    //int i = index(l,m,f);
                    int ll = 2*(LMAX-l)+f;
                    for (l_=beginl; l_<=LMAX; ++l_) {
                        for (f_=0; f_<=1; ++f_) {
                            //int i_ = index(l_,m,f_);
                            int ll_ = 2*(LMAX-l_)+f_;
                            ScatMatrix[mm][ll_][ll] =0.;
                            for (l__=beginl; l__<=LMAX; ++l__) {
                                for (f__=0; f__<=1; ++f__) {
                                    //int i__ = index(l__,m,f__);
                                    int ll__ = 2*(LMAX-l__)+f__;

                                    ScatMatrix[mm][ll_][ll] += BMatrix[mm][ll_][ll__]*ScatMatrixInverse[mm][ll__][ll];
                                }
                            }
                        }
                    }
                }
            }
        }

        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (l=beginl; l<=LMAX; ++l) {
                for (f=0; f<=1; ++f) {
                    //int i = index(l,m,f);
                    int ll = 2*(LMAX-l)+f;
                    for (l_=beginl; l_<=LMAX; ++l_) {
                        for (f_=0; f_<=1; ++f_) {
                            //int i_ = index(l_,m,f_);
                            int ll_ = 2*(LMAX-l_)+f_;
                            ScatMatrixInverse[mm][ll_][ll]=ScatMatrix[mm][ll_][ll];
                        }
                    }
                }
            }
        }

        //
        // Invert the matrix B^(-1)*(1-B*A)
        //
        invert_block_diagonal(ScatMatrix);

        // Reset previous geometry...
        old_thetai=-1000;
        old_thetas=-1000;
        old_phis=-1000;
    }

    //
    // The following function must be called before any attempt to evaluate
    // the scattering field.  It sets the geometry for the measurement...
    //
    void
    Axisymmetric_Particle_BRDF_Model::
    set_geometry(double _thetai,double _thetas,double _phis)
    {
        SETUP();

        if (abs(_thetai)<1E-6) _thetai = 1E-6;
        if (abs(_thetas)<1E-6) _thetas = 1E-6;

        if (thetai < 0) {
            _thetai = -_thetai;
            _phis = pi + _phis;
        }

        if (_thetas < 0) {
            _thetas = -_thetas;
            _phis = pi + _phis;
        }

        if (_thetai != old_thetai) {
            calculate_W(_thetai);
            old_thetai = _thetai;
        }

        if (_thetas != old_thetas) {
            calculate_Z(pi - _thetas);
            old_thetas = _thetas;
        }

        if (_phis != old_phis) {
            calculate_eIP(_phis);
            old_phis = _phis;
        }
    }

    // The electric field version of the p-->p
    // differential scattering cross-section...
    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    Epp(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai, _thetas, _phis);
        return E(Wp,Zp);
    }

    // The electric field version of the p-->s
    // differential scattering cross-section...
    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    Eps(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai, _thetas, _phis);
        return E(Wp,Zs);
    }

    // The electric field version of the s-->p
    // differential scattering cross-section...
    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    Esp(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai, _thetas, _phis);
        return E(Ws,Zp);
    }

    // The electric field version of the s-->s
    // differential scattering cross-section...
    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    Ess(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai, _thetas, _phis);
        return E(Ws,Zs);
    }

    //
    // Jones matrix differential scattering cross section...
    //
    JonesMatrix
    Axisymmetric_Particle_BRDF_Model::
    jonesDSC()
    {
        set_geometry(thetai,thetas,phis);

        COMPLEX pp=E(Wp,Zp);
        COMPLEX ps=E(Wp,Zs);
        COMPLEX sp=E(Ws,Zp);
        COMPLEX ss=E(Ws,Zs);

        return JonesMatrix(pp,ss,ps,sp);
    }

    //
    // Constructor for class Axisymmetric_Particle_BRDF_Model...
    //
    Axisymmetric_Particle_BRDF_Model::
    Axisymmetric_Particle_BRDF_Model()
    {
        n3=vacuum;
        // Initialize array addresses to zero ...
        ScatMatrix=0;
        ScatMatrixInverse=0;

        // Initialize old angles ...
        old_thetai=-1000.;
        old_thetas=-1000.;
        old_phis=-1000;
        LMAX = 0;
        MMAX = 0;
    }

    MuellerMatrix
    Axisymmetric_Particle_BRDF_Model::
    Specular(double _theta)
    {
        SETUP();

        // This function returns the Mueller matrix specular reflectance (for type==0 or 
        // type==2) or the regular transmittance (for type==1 or type==3) for an incident 
        // angle of theta (in radians). The result is derived from the optical theorem 
        // and is only valid when density is suffiently low that multiple scattering 
        // between spheres can be neglected.

        switch (type) {
        case 0:
        {
            JonesMatrix X = JonesDSC(_theta, _theta, 0., 0.);
            JonesMatrix r = stack->r12(_theta, lambda, vacuum, substrate);
            r *= exp(2. * cI * qq * cos(_theta));
            MuellerMatrix sigma = (4. * pi / k) * ReCrossMueller(X, r);

            return MuellerMatrix(r) - sigma * (density / cos(_theta));
        }
        break;
        case 1:
        {
            double index = substrate.n(lambda);
            double sint = sin(_theta) / index;
            if (sint >= 1. || substrate.k(lambda) != 0) return MuellerZero();
            double thetat = asin(sint);
            JonesMatrix X = JonesDSC(_theta, thetat, 0., 0.);
            COMPLEX phase = exp(cI * qq * (cos(_theta) - index * cos(thetat)));
            //double factor = cos(thetat)/cos(_theta);
            JonesMatrix t = stack->t12(_theta, lambda, vacuum, substrate);
            t *= phase;
            MuellerMatrix sigma = (4. * pi / k / sqrt(index)) * ReCrossMueller(X, t);

            return MuellerMatrix(t) * (cos(thetat) / cos(_theta) * index) - sigma * (density / cos(_theta));
        }
        break;
        case 2:
        {
            double index = substrate.n(lambda);
            JonesMatrix X = JonesDSC(_theta, _theta, 0., 0.);
            JonesMatrix r = stack->r21i(_theta, lambda, substrate, vacuum);
            r *= exp(-2. * cI * index * qq * cos(_theta));
            MuellerMatrix sigma = (4. * pi / k / index) * ReCrossMueller(X, r);

            return MuellerMatrix(r) - sigma * (density / cos(_theta));
        }
        break;
        case 3:
        {
            double index = substrate.n(lambda);
            double sint = sin(_theta) * index;
            COMPLEX cost = sqrt(1. - sqr(sint));
            if (imag(cost) < 0) cost = -cost;
            if (sint >= 1.) return MuellerZero();
            double thetat = asin(sint);
            JonesMatrix X = JonesDSC(_theta, thetat, 0., 0.);
            COMPLEX phase = exp(cI * qq * (cost - index * cos(_theta)));
            //double factor = cos(thetat)/cos(_theta);
            JonesMatrix t = stack->t21i(_theta, lambda, substrate, vacuum);
            t *= phase;
            MuellerMatrix sigma = (4. * pi / k / sqrt(index)) * ReCrossMueller(X, t);

            return MuellerMatrix(t) * (cos(_theta) / index / cos(thetat)) - sigma * (density / cos(_theta));
        }
        break;
        default:
            error("Invalid type = " + to_string(type));
        }
        return MuellerZero();
    }

    DEFINE_MODEL(Axisymmetric_Particle_BRDF_Model,Local_BRDF_Model,
                 "Axisymmetric particle on a substrate.");

    DEFINE_PTRPARAMETER(Axisymmetric_Particle_BRDF_Model,
                        Axisymmetric_Shape_Ptr,
                        Shape,
                        "Particle Shape",
                        "Ellipsoid_Axisymmetric_Shape",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,dielectric_function,particle,"Particle optical properties","(1.59,0)",0xFF);
    DEFINE_PTRPARAMETER(Axisymmetric_Particle_BRDF_Model,StackModel_Ptr,stack,"Substrate films","No_StackModel",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,double,delta,"Separation of particle from substrate [um] (in contact: 0)","0",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,int,lmax,"Maximum polar order (lmax)","0",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,int,mmax,"Maximum azimuthal order (mmax)","0",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,int,order,"Perturbation order (exact: -1)","-1",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,int,Norm_Inc_Approx,"Normal incidence approximation (exact: 0)","0",0xFF);
    DEFINE_PARAMETER(Axisymmetric_Particle_BRDF_Model,int,improve,"Iterative improvement steps (recommend: 3)","3",0xFF);


} // namespace SCATMECH




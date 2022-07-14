//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: axifree.cpp
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
#include "axifree.h"
#include "dielfunc.h"
#include "askuser.h"
#include "matrix3d.h"
#include "float.h"

using namespace std;

namespace SCATMECH {

    using namespace BobVlieg_Supp;

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
    // Bmatrix() calculates the Mie scattering matrix found
    // in Eq. (4.2) of B&V.
    //
    void
    TMatrix_Axisymmetric_Scatterer::
    Bmatrix(ScatterTMatrix& T)
    {
        int m;
        for (m=-MMAX; m<=MMAX; ++m) {
            int mm = LMAX+m;
            int beginl = (m==0) ? 1 : abs(m);
            for (int l=beginl; l<=LMAX; ++l) {
                for (int f=0; f<=1; ++f) {
                    int i = index(l,m,f);
                    int ll = 2*(LMAX-l)+f;
                    for (int l_=beginl; l_<=LMAX; ++l_) {
                        for (int f_=0; f_<=1; ++f_) {
                            int i_ = index(l_,m,f_);
                            int ll_ = 2*(LMAX-l_)+f_;
                            T[mm][ll_][ll]=0;
                        }
                    }
                }
            }
        }

        COMPLEX relN2 = N2/N0;
        Shape->Get_TMatrix(T, LMAX, MMAX, k*N0, relN2);

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
    // The complex arccosine is not part of the standard C++ library (!)...
    //
    static COMPLEX arccosine(const COMPLEX& a)
    {
        double x=real(a);
        double y=imag(a);

        return pi/2. - arg(sqrt(1. - sqr(x + cI*y)) +
                           cI*(x + cI*y)) +
               cI*log(abs(sqrt(1. - sqr(x + cI*y)) +
                          cI*(x + cI*y)));
    }

    //
    // setup() is called whenever the model needs to be "recalculated"
    //
    void
    TMatrix_Axisymmetric_Scatterer::
    setup()
    {
        Free_Space_Scatterer::setup();

        Shape->Write("shape.dat");
        double length = Shape->Get_Base_Length();

        k = 2.0*pi/lambda;
        double q = k * Shape->Get_MaxRadius();

        // Bohren and Huffman suggested lmax...
        BH_LMAX = (int)(q+4.05*pow(q,0.3333)+2.);

        if (lmax==0) LMAX = BH_LMAX;
        if (lmax>0) LMAX = lmax;
        if (lmax<0) LMAX = BH_LMAX+abs(lmax);

        if (mmax==0) MMAX = LMAX;
        if (mmax>0) MMAX = mmax;
        if (mmax<0) MMAX = LMAX-abs(mmax);

        if (MMAX>LMAX) MMAX = LMAX;

        if (medium.k(lambda)!=0) error("Surrounding medium must be non-absorbing.");

        N2 = COMPLEX(particle.index(lambda));
        N0 = medium.n(lambda);

        ScatMatrix.resize(LMAX);
        Bmatrix(ScatMatrix);

        // The Bohren and Huffman recommended LMAX...
        // BH_LMAX = (int)(q+4.*pow(q,0.3333)+2.);
        BH_LMAX = LMAX;

        sqrsize = index(LMAX,LMAX,1)+1;

        lvector.resize(LMAX+1);
        for (int ll=1; ll<=LMAX; ++ll) {
            lvector[ll]=sqrt((2.*ll+1.)/(ll*(ll+1.)));
        }

        Zp.resize(sqrsize);
        Zs.resize(sqrsize);
        eIP.resize(sqrsize);
        Wp.resize(sqrsize);
        Ws.resize(sqrsize);
        Vp.resize(sqrsize);
        Vs.resize(sqrsize);

        vector<COMPLEX> B(sqrsize);

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
    TMatrix_Axisymmetric_Scatterer::
    set_geometry(double thetai,double thetas,double phis)
    {
        SETUP();

        /*
        if (thetai<0) {
            thetai=-thetai;
            phis = pi+phis;
        }

        if (thetas<0) {
            thetas=-thetas;
            phis=pi+phis;
        }
        */

        if (thetai>=pi) {
            thetai-=2*pi;
        }
        if (thetai<-pi) {
            thetai+=2*pi;
        }

        if (thetas>=pi) {
            thetas-=2*pi;
        }

        if (thetas<-pi) {
            thetas+=2*pi;
        }

        if (thetai<0) {
            thetai=-thetai;
            phis = pi+phis;
        }

        if (thetas<0) {
            thetas=-thetas;
            phis=pi+phis;
        }

        if (thetai!=old_thetai) {
            calculate_W(thetai);
            old_thetai=thetai;
        }

        if (thetas!=old_thetas) {
            calculate_Z(pi-thetas);
            old_thetas=thetas;
        }

        if (phis!=old_phis) {
            calculate_eIP(phis);
            old_phis=phis;
        }
    }

    JonesMatrix
    TMatrix_Axisymmetric_Scatterer::
        jones(const Vector& kin, const Vector& kout)
    {
        SETUP();

        double inphi = 0;  atan2(kin.y, kin.x);
        double cosinphi = cos(inphi);
        double sininphi = sin(inphi);

        Vector _kin(kin.x*cosinphi+kin.y*sininphi,-kin.x*sininphi+kin.y*cosinphi,kin.z);
        Vector _kout(kout.x*cosinphi+kout.y*sininphi,-kout.x*sininphi+kout.y*cosinphi,kout.z);

        _kin = -_kin;

        if (fabs(_kin.y)<1E-5 && fabs(_kin.x)<1E-5) {
            if (_kin.z>0) _kin=polar(sqrt(1-sqr(1E-5)),1E-5,0);
            else _kin=polar(-sqrt(1-sqr(1E-5)),1E-5,0);
        }
        if (fabs(_kout.y)<1E-5 && fabs(_kout.x)<1E-5) {
            if (_kout.z>0) _kout=polar(sqrt(1-sqr(1E-5)),1E-5,0);
            else _kout=polar(-sqrt(1-sqr(1E-5)),1E-5,0);
        }

        double _thetai = acos(_kin.z);
        double _thetas = acos(_kout.z);

        double _phis =  atan2(_kout.y,_kout.x);

        set_geometry(_thetai,_thetas,_phis);

        COMPLEX pp=E(Wp,Zp);
        COMPLEX ps=E(Wp,Zs);
        COMPLEX sp=E(Ws,Zp);
        COMPLEX ss=E(Ws,Zs);

        JonesMatrix scatter = JonesMatrix(pp,ss,ps,sp)*k*COMPLEX(0,1);

        Vector sr = perpto(_kout,Vector(0,0,-1));
        Vector si = perpto(-_kin,Vector(0,0,1));
        Vector pr = perpto(_kout,sr);
        Vector pi = perpto(-_kin,si);

        Vector ai = perpto(-_kin,_kout);
        Vector ar = ai;
        Vector bi = -perpto(-_kin,ai);
        Vector br = -perpto(_kout,ar);

        JonesMatrix rotatein = GetJonesRotator(bi,ai,si,pi).transpose();
        JonesMatrix rotateout = GetJonesRotator(br,ar,sr,pr);

		return rotateout*scatter*rotatein;
    }

    //
    // Constructor for class TMatrix_Axisymmetric_Scatterer...
    //
    TMatrix_Axisymmetric_Scatterer::
    TMatrix_Axisymmetric_Scatterer()
    {
    }

    //
    // The element of the incident plane wave vector for p-polarized light
    // from Eq. (2.14) of BV&G...
    //
    COMPLEX
    TMatrix_Axisymmetric_Scatterer::
    VIp(int l,int m,int f,double thetai) const
    {
        COMPLEX result;
        COMPLEX costheta2=cos(thetai/2.);
        COMPLEX sintheta2=sin(thetai/2.);
        if (f==efield) {
            // Eq. (2.14a) of BV&G
            result =  1./k*ipow(l-1)*lvector[l]*mpow(m-1)*
                      dminus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.14b) of BV&G
            result =  -1./k*ipow(l)*lvector[l]*mpow(m-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result;
    }


    //
    // The element of the incident plane wave vector for s-polarized light
    // from Eq. (2.16) of BV&G...
    //
    COMPLEX
    TMatrix_Axisymmetric_Scatterer::
    VIs(int l,int m,int f,double thetai) const
    {
        COMPLEX result;
        COMPLEX costheta2=cos(thetai/2.);
        COMPLEX sintheta2=sin(thetai/2.);

        if (f==efield) {
            // Eq. (2.16a) of BV&G
            result =  -cI/k*
                      ipow(l-1)*lvector[l]*mpow(m-1)*
                      dplus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.16b) of BV&G
            result = cI/k*
                     ipow(l)*lvector[l]*mpow(m-1)*
                     dminus(l,m,costheta2,sintheta2);
        }
        return result;
    }

    //
    // Routine to calculate the scattered wave vector for a specific incident
    // angle.  This is evaluation of Eq. (5.3) of B&V...
    //
    void
    TMatrix_Axisymmetric_Scatterer::
    calculate_W(double thetai)
    {
        int i,m,l,f,l_,f_;

        // First calculate the V vector of the incident wave...
        for (l=1; l<=LMAX; ++l) {
            for (m=-l; m<=l; ++m) {
                for (f=0; f<=1; ++f) {
                    i=index(l,m,f);
                    Vp[i]=VIp(l,m,f,thetai);
                    Vs[i]=VIs(l,m,f,thetai);
                }
            }
        }

        // Then multiply the incident vector V by the scattering matrix
        for (m=-MMAX; m<=MMAX; ++m) {
            int beginl = (m==0) ? 1 : abs(m);
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    Wp[i_]=0;
                    Ws[i_]=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            Wp[i_]+=ScatMatrix[mm][ll_][ll]*Vp[i];
                            Ws[i_]+=ScatMatrix[mm][ll_][ll]*Vs[i];
                        }
                    }
                }
            }
        }
    }

    //
    // The vector Z is that part of the scattering function which depends
    // upon the scattering angle thetas.  It must be evaluated anytime
    // that thetas changes...
    //
    void
    TMatrix_Axisymmetric_Scatterer::
    calculate_Z(double thetas)
    {
        double d_cost = cos(thetas);
        if (d_cost > 0.999999) d_cost = 0.999999;
        if (d_cost < -0.999999) d_cost = -0.999999;

        COMPLEX cost= d_cost;
        double sint = real(sqrt(1.-sqr(cost)));

        for (int i=0,l=1; l<=LMAX; ++l) {
            if (l<=BH_LMAX) {
                int _mmax = (MMAX>l) ? l : MMAX;
                for (int m=-_mmax; m<=_mmax; ++m) {
                    COMPLEX temp1 = Ptilde(l,m,cost)*(double)m/sint;
                    COMPLEX temp2 = P_tilde(l,m,cost);
                    for (int f=0; f<=1; ++f,++i) {
                        if (f==efield) {
                            Zp[i] = ipow(l)*mpow(l)*temp2;
                            Zs[i] = -ipow(l+1)*mpow(l+1)*temp1;
                        } else {
                            Zp[i] = -ipow(l+1)*mpow(l+1)*temp1;
                            Zs[i] = -ipow(l)*mpow(l)*temp2;
                        }
                    }
                }
            } else {
                int _mmax = (MMAX>l) ? l : MMAX;
                for (int m=-_mmax; m<=_mmax; ++m) {
                    for (int f=0; f<=1; ++f,++i) {
                        Zp[i] = 0.;
                        Zs[i] = 0.;
                        Zp[i] = 0.;
                        Zs[i] = 0.;
                    }
                }
            }
        }
    }

    //
    // The vector eIP is that part of the scattering function which depends
    // upon the azimuthal scattering angle phis.  It must be evaluated anytime
    // that phis changes...
    //
    // This function basically calculates exp(i m phis)...
    //
    void
    TMatrix_Axisymmetric_Scatterer::
    calculate_eIP(double phis)
    {
        vector<COMPLEX> expmphi(2*LMAX+2);

        COMPLEX expphi=exp(cI*phis);
        expmphi[LMAX]=1.;

        for (int m=1; m<=LMAX; ++m) {
            expmphi[LMAX+m]=expmphi[LMAX+m-1]*expphi;
            expmphi[LMAX-m]=expmphi[LMAX-m+1]/expphi;
        }

        for (int i=0,l=1; l<=BH_LMAX; ++l) {
            int _mmax = (MMAX>l) ? l : MMAX;
            for (int m=-_mmax; m<=_mmax; ++m) {
                for (int f=0; f<=1; ++f,++i) {
                    eIP[i] = expmphi[m+LMAX];
                }
            }
        }
    }

    //
    // The total scattered electric field depends upon the incident vector W and
    // the scattered vector Z.
    //
    COMPLEX
    TMatrix_Axisymmetric_Scatterer::
    E(vector<COMPLEX>& W,vector<COMPLEX>& Z)
    {
        COMPLEX E(0.);

        for (int i=0,l=1; l<=BH_LMAX; ++l) {
            int _mmax = (MMAX>l) ? l : MMAX;
            for (int m=-_mmax; m<=_mmax; ++m) {
                for (int f=0; f<=1; ++f,++i) {
                    E += W[i]*Z[i]*eIP[i];
                }
            }
        }

        return E;
    }


    DEFINE_MODEL(TMatrix_Axisymmetric_Scatterer,Free_Space_Scatterer,
                 "Axisymmetric particle solved using T-Matrix method");

    DEFINE_PTRPARAMETER(TMatrix_Axisymmetric_Scatterer,
                        Axisymmetric_Shape_Ptr,Shape,
                        "Particle Shape",
                        "Ellipsoid_Axisymmetric_Shape",0xFF);
    DEFINE_PARAMETER(TMatrix_Axisymmetric_Scatterer,dielectric_function,particle,"Particle","(1.59,0)",0xFF);
    DEFINE_PARAMETER(TMatrix_Axisymmetric_Scatterer,int,lmax,"Maximum polar order","0",0xFF);
    DEFINE_PARAMETER(TMatrix_Axisymmetric_Scatterer,int,mmax,"Maximum azimuthal order","0",0xFF);


} // namespace SCATMECH




//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rcw.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "rcw.h"
#include "matrixmath.h"

using namespace std;

//
// Equation numbers refer to Moharam et al. J. Opt. Soc. Am. A 12, 1068 (1995) or
// Moharam et al. J. Opt. Soc. Am. A 12, 1077 (1995)

namespace SCATMECH {

    using namespace CMLIB;

    namespace {
        template <class T> const T& min(const T& a,const T& b) {
            return (a<b) ? a : b;
        }
        template <class T> const T& max(const T& a,const T& b) {
            return (a>b) ? a : b;
        }
        static const COMPLEX cI(0,1);
    }

    CVector RCW_Model::GetEField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        SETUP();

        CVector result(0.,0.,0.);

        if (type==0 && pos.z<0) return result;
        if (type==1 && pos.z>-totalthickness) return result;
        if (type==2 && pos.z>-totalthickness) return result;
        if (type==3 && pos.z<0) return result;

        Vector R = pos;
        if (type==1) R.z += totalthickness;
        if (type==2) R.z += totalthickness;

        if (type==0 && incident) {
            Vector kin(nI*k0*sin(thetai*deg),0.,-nI*k0*cos(thetai*deg));
            Vector s(0.,1.,0.);
            Vector p = cross(kin,s)/(k0*nI);
            COMPLEX phase = exp(cI*COMPLEX(kin*R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        if (type==2 && incident) {
            CVector kin(conj(nII)*k0*sin(thetai*deg),0.,conj(nII)*k0*cos(thetai*deg));
            CVector s(0.,1.,0.);
            CVector p = cross(kin,s)/(k0*conj(nII));
            COMPLEX phase = exp(cI*COMPLEX(kin*(CVector)R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        for (int i=0; i<nmat; ++i) {
            JonesVector out = r[i]*inpol;
            double x = real(kvector[i].x);
            double y = real(kvector[i].y);
            double rho = sqrt(x*x+y*y);
            COMPLEX k = k0 * (type==0 || type==3 ? nI : conj(nII));
            CVector s = rotation!=0 ? CVector(-y/rho,x/rho,0) : CVector(0,1,0);
            CVector p = cross(kvector[i],s)/k;
            COMPLEX phase = exp(cI*(kvector[i]*(CVector)R));
            result += (out.S()*s+out.P()*p)*phase;
        }
        return result;
    }

    CVector RCW_Model::GetHField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        SETUP();

        CVector result(0.,0.,0.);

        if (type==0 && pos.z<0) return result;
        if (type==1 && pos.z>-totalthickness) return result;
        if (type==2 && pos.z>-totalthickness) return result;
        if (type==3 && pos.z<0) return result;

        Vector R = pos;
        if (type==1) R.z += totalthickness;
        if (type==2) R.z += totalthickness;

        if (type==0 && incident) {
            Vector kin(nI*k0*sin(thetai*deg),0.,-nI*k0*cos(thetai*deg));
            Vector p(0.,-1,0.);
            Vector s = cross(p,kin)/(k0*nI);
            COMPLEX phase = nI*exp(cI*COMPLEX(kin*R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        if (type==2 && incident) {
            CVector kin(conj(nII)*k0*sin(thetai*deg),0.,conj(nII)*k0*cos(thetai*deg));
            CVector p(0.,-1.,0.);
            CVector s = cross(p,kin)/(k0*conj(nII));
            COMPLEX phase = conj(nII)*exp(cI*COMPLEX(kin*(CVector)R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        for (int i=0; i<nmat; ++i) {
            JonesVector out = r[i]*inpol;
            double x = real(kvector[i].x);
            double y = real(kvector[i].y);
            COMPLEX N = (type==0 || type==3 ? nI : conj(nII));
            double rho = sqrt(x*x+y*y);
            COMPLEX k = k0 * N;
            CVector p = rotation!=0 ? CVector(y/rho,-x/rho,0) : CVector(0,-1,0);
            CVector s = cross(p,kvector[i])/k;
            COMPLEX phase = N*exp(cI*(kvector[i]*(CVector)R));
            result += (out.S()*s+out.P()*p)*phase;
        }
        return result;
    }

    CVector RCW_Model::GetBField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        return GetHField(inpol,pos,incident);
    }

    CVector RCW_Model::GetDField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        return GetEField(inpol,pos,incident) * (type==0 || type==3 ? sqr(nI) : sqr(conj(nII)));
    }


    void
    RCW_Model::setup()
    {
        int i;

        Model::setup();

        if (order<0) error("order < 0");

        if (grating->get_lambda()!=lambda) grating->set_lambda(lambda);

        // Exchange some information with grating...
        period = grating->get_period();

        totalthickness=0;
        for (i=0; i<grating->get_levels(); ++i) {
            totalthickness+=grating->get_thickness(i);
        }

        dielectric_function medium_i = grating->get_medium_i();
        dielectric_function medium_t = grating->get_medium_t();

        k0 = 2*pi/lambda;
        n = order;
        nmat = 2*n+1;

        if (medium_i.k(lambda)!=0.)
            error("medium_i must be non-absorbing");

        nI = medium_i.n(lambda);
        eI = sqr(nI);

        // Note: Moharam, et al., use n-ik convention...
        nII = conj((COMPLEX)medium_t.index(lambda));
        eII = sqr(nII);

        kxi.resize(nmat);
        Kx.resize(nmat);
        kIzi.resize(nmat);
        kIIzi.resize(nmat);
        YI.resize(nmat);
        YII.resize(nmat);
        ZI.resize(nmat);
        ZII.resize(nmat);
        V.resize(nmat);
        kvector.resize(nmat);
        r.resize(nmat);
        R.resize(nmat);

        if (type==1 || type==2 || type==3) {
            if (imag(nII)!=0) {
                error("medium_t must be non-absorbing for type == 1, 2, or 3");
            }
        }

        // Incident wavevector, with coordinates aligned along grating.
        // Note that SCATMECH generally assumes the z axis points out
        // of the grating.  rotation represents a counterclockwise rotation
        // of the grating, looking down at the surface.
        if (type==0 || type==1) {
            inckx = k0*nI*sin(thetai*deg)*cos(rotation*deg);
            incky = -k0*nI*sin(thetai*deg)*sin(rotation*deg);
            inckz = -k0*nI*cos(thetai*deg);
            // This program has a problem if kx of incident light is zero...
            if (fabs(inckx/k0/nI)<1E-7) {
                inckx = k0*nI*1E-7;
                inckz = -k0*sqrt(eI - sqr(inckx/k0) - sqr(incky/k0));
            }
        } else {
            inckx = k0*real(nII)*sin(thetai*deg)*cos(rotation*deg);
            incky = -k0*real(nII)*sin(thetai*deg)*sin(rotation*deg);
            inckz = -k0*real(nII)*cos(thetai*deg);
            // This program has a problem if kx of incident light is zero...
            if (fabs(inckx/k0/real(nII))<1E-7) {
                inckx = k0*real(nII)*1E-7;
                inckz = -k0*sqrt(real(eII) - sqr(inckx/k0) - sqr(incky/k0));
            }
        }

        // Eq. (51) ...
        // Moharam et al. assume z axis points into grating, so y(Moharam) = -y(SCATMECH)
        ky = -incky;

        // Kvector (magnitude) associated with grating...
        double Kgrating = 2.*pi/period;

        for (i=0; i<nmat; ++i) {
            // Eqs. (6), (50) ...
            kxi[i] = inckx -(i-n)*Kgrating;

            Kx[i] = kxi[i]/k0;
            COMPLEX kk = sqr(Kx[i]);
            double mu = 1.;

            // Eq. (7), (52) ... (written differently, because nII is not necc. real)...
            COMPLEX aa = k0*sqrt(eI*mu - kk - sqr(ky/k0));
            if (imag(aa)>0) aa = conj(aa);
            kIzi[i]  = aa;
            COMPLEX bb = k0*sqrt(eII*mu - kk - sqr(ky/k0));
            if (imag(bb)>0) bb = conj(bb);
            kIIzi[i] = bb;

            // Defined after Eq. (24)...
            YI[i] = kIzi[i]/k0/mu;
            YII[i] = kIIzi[i]/k0/mu;

            // Defined after Eq. (44)...
            ZI[i] = kIzi[i]/k0/sqr(nI);
            ZII[i] = kIIzi[i]/k0/sqr(nII);
        }

        // Min and max orders that propagate...
        min_order=order+1;
        max_order=-order-1;

        // Rotation matrix elements to rotate from sample coordinates to
        // incident wave coordinates...
        double cost = cos(rotation*deg);
        double sint = sin(rotation*deg);

        for (i=-order; i<=order; ++i) {
            int j=i+n;
            R[j] = JonesZero();
            V[j] = Vector(0,0,0);
            double kx = inckx-i*Kgrating;
            double ky = incky;

            if (type==0 || type==3) {
                COMPLEX eps = medium_i.epsilon(lambda);
                COMPLEX kz = sqrt(sqr(k0)*eps-sqr(kx)-sqr(ky));

                // Diffracted wavevector in incident wave coordinates...
                kvector[j] = CVector(kx*cost-ky*sint,kx*sint+ky*cost,kz);

                // If it is propagating...
                if (imag(kz)==0.) {
                    V[j] = Vector(real(kvector[j].x),real(kvector[j].y),real(kvector[j].z))/(k0*nI);
                    if (i<min_order) min_order=i;
                    if (i>max_order) max_order=i;
                } else {
                    V[j] = Vector(0.,0.,0.);
                }
            } else { // type==1 || type==2
                COMPLEX eps = medium_t.epsilon(lambda);
                COMPLEX kz = -sqrt(sqr(k0)*eps-sqr(kx)-sqr(ky));

                // Diffracted wavevector in incident wave coordinates...
                kvector[j] = CVector(kx*cost-ky*sint,kx*sint+ky*cost,kz);

                // If it is propagating...
                if (imag(kz)==0) {
                    V[j] = Vector(real(kvector[j].x),real(kvector[j].y),-real(kvector[j].z))/(k0*real(nII));
                    if (i<min_order) min_order=i;
                    if (i>max_order) max_order=i;
                } else {
                    V[j] = Vector(0.,0.,0.);
                }
            }
        }

        if (rotation==0) {
            switch (type) {
                case 0:
                    InPlaneReflection(false);
                    break;
                case 1:
                    InPlaneTransmission(false);
                    break;
                case 2:
                    InPlaneReflection(true);
                    break;
                case 3:
                    InPlaneTransmission(true);
                    break;
                default:
                    error("Invalid type (must be 0, 1, 2, or 3)");
            }
        } else {
            switch (type) {
                case 0:
                    if (grating->is_anisotropic()) {
                        ConicalAnisoReflection(false);
                    } else {
                        ConicalReflection(false);
                    }
                    break;
                case 1:
                    if (grating->is_anisotropic()) {
                        ConicalAnisoTransmission(false);
                    } else {
                        ConicalTransmission(false);
                    }
                    break;
                case 2:
                    if (grating->is_anisotropic()) {
                        ConicalAnisoReflection(true);
                    } else {
                        ConicalReflection(true);
                    }
                    break;
                case 3:
                    if (grating->is_anisotropic()) {
                        ConicalAnisoTransmission(true);
                    } else {
                        ConicalTransmission(true);
                    }
                    break;
                default:
                    error("Invalid type (must be 0, 1, 2, or 3)");
            }
        }

        kxi.resize(0);
        Kx.resize(0);
        kIzi.resize(0);
        YI.resize(0);
        YII.resize(0);
        ZI.resize(0);
        ZII.resize(0);
    }

    //***********************************************************************
    //**                                                                   **
    //** In-Plane Reflection                                               **
    //**                                                                   **
    //***********************************************************************
    void RCW_Model::InPlaneReflection(bool backward)
    {
        int i,j,k;

        int nnmat = 2*nmat;

        vector<COMPLEX> epsinvx(2*nmat+1),epsy(2*nmat+1),epsz(2*nmat+1);
        vector<COMPLEX> muinvx(2*nmat+1),muy(2*nmat+1),muz(2*nmat+1);
        vector<COMPLEX> Ex(sqr(nmat)),Ey(sqr(nmat)),Ez(sqr(nmat));
        vector<COMPLEX> Exinv(sqr(nmat)),Eyinv(sqr(nmat)),Ezinv(sqr(nmat));
        vector<COMPLEX> Mx(sqr(nmat)),My(sqr(nmat)),Mz(sqr(nmat));
        vector<COMPLEX> Mxinv(sqr(nmat)),Myinv(sqr(nmat)),Mzinv(sqr(nmat));

        vector<COMPLEX> A(sqr(nmat));
        vector<COMPLEX> B(sqr(nmat));
        vector<COMPLEX> EB(sqr(nmat));
        vector<COMPLEX> MA(sqr(nmat));
        vector<COMPLEX> Ws(sqr(nmat)),Qs(nmat),Xs(nmat),Vs(sqr(nmat));
        vector<COMPLEX> Wp(sqr(nmat)),Qp(nmat),Xp(nmat),Vp(sqr(nmat));
        vector<COMPLEX> Fs(sqr(nnmat)),Fp(sqr(nnmat));
        vector<int> indexs(nnmat),indexp(nnmat);
        vector<COMPLEX> WXVXs(nnmat);
        vector<COMPLEX> WXVXp(nnmat);
        vector<COMPLEX> as(2*nmat*nmat),ap(2*nmat*nmat);
        vector<COMPLEX> fs(nmat*nmat),gs(nmat*nmat);
        vector<COMPLEX> fp(nmat*nmat),gp(nmat*nmat);

        for (i=0; i<nmat; ++i) {
            for (j=0; j<nmat; ++j) {
                int ij = i+nmat*j;
                fs[ij] = (i==j) ? 1. : 0.;
                fp[ij] = (i==j) ? 1. : 0.;
                if (backward) {
                    gs[ij] = (i==j) ? cI*YI[i] : 0.;
                    gp[ij] = (i==j) ? cI*ZI[i] : 0.;
                } else {
                    gs[ij] = (i==j) ? cI*YII[i] : 0.;
                    gp[ij] = (i==j) ? cI*ZII[i] : 0.;
                }
            }
        }

        forloopint layerloop(0,grating->get_levels()-1,backward ? +1 : -1 );
        for (forloopint::iterator layer = layerloop.begin(); layer!=layerloop.end(); ++layer) {

            string layerstr = format("(layer = %d)",(int)layer);

            double d = grating->get_thickness(layer);
            if (d<0.) error("Thickness of layer " + to_string((int)layer) + " = " + to_string(d) + " less than zero");

            //
            // Get the Fourier components for this layer...
            //
            message(layerstr + "Creating E and M matrices");
            for (i=-nmat; i<=nmat; ++i) {
                if (grating->is_anisotropic()) {
                    epsinvx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    epsy[i+nmat]=conj(grating->fouriery(i,layer,0));
                    epsz[i+nmat]=conj(grating->fourierz(i,layer,0));
                } else {
                    epsinvx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    epsy[i+nmat]=conj(grating->fourierx(i,layer,0));
                    epsz[i+nmat]=epsy[i+nmat];
                }
                if (grating->is_magnetic()) {
                    muinvx[i+nmat]=conj(grating->fouriermux(i,layer,1));
                    muy[i+nmat]=conj(grating->fouriermuy(i,layer,0));
                    muz[i+nmat]=conj(grating->fouriermuz(i,layer,0));
                } else {
                    muinvx[i+nmat]= (i==0 ? 1.: 0); 
                    muy[i+nmat]= (i==0 ? 1.: 0);
                    muz[i+nmat]= (i==0 ? 1.: 0);
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Ex[i+nmat*j] = Exinv[i+nmat*j] = epsinvx[i-j+nmat];
                    Ey[i+nmat*j] = Eyinv[i+nmat*j] = epsy[i-j+nmat];
                    Ez[i+nmat*j] = Ezinv[i+nmat*j] = epsz[i-j+nmat];
                    Mx[i+nmat*j] = Mxinv[i+nmat*j] = muinvx[i-j+nmat];
                    My[i+nmat*j] = Myinv[i+nmat*j] = muy[i-j+nmat];
                    Mz[i+nmat*j] = Mzinv[i+nmat*j] = muz[i-j+nmat];
                }
            }

            message(layerstr + "Inverting E and M matrices");
            Inverse(Ex,nmat);
            Inverse(Eyinv,nmat);
            Inverse(Ezinv,nmat);
            if (grating->is_magnetic()) {
                Inverse(Mx,nmat);
                Inverse(Myinv,nmat);
                Inverse(Mzinv,nmat);
            }

            //
            // Building matrix A...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    A[i+nmat*j] = Kx[i]*Mzinv[i+nmat*j]*Kx[j] - Ey[i+nmat*j];
                    //A[i+nmat*j] = (i==j ? sqr(Kx[i]) : 0) - Ey[i+nmat*j];
                }
            }

            //
            // Building matrix B...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    // Eq. (36) ...
                    //B[i+nmat*j] = Kx[i] * Ezinv[i+nmat*j] * Kx[j] - (i==j ? 1. : 0.);
                    B[i+nmat*j] = Kx[i] * Ezinv[i+nmat*j] * Kx[j] - My[i+nmat*j];
                }
            }

            //
            // Building matrix M*A...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    COMPLEX& _MA = MA[i+nmat*j];
                    _MA = 0;
                    for (int k=0; k<nmat; ++k) {
                        _MA += Mx[i+nmat*k]*A[k+nmat*j];
                    }
                }
            }

            //
            // Building matrix E*B...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    COMPLEX& _EB = EB[i+nmat*j];
                    _EB = 0;
                    for (int k=0; k<nmat; ++k) {
                        _EB += Ex[i+nmat*k]*B[k+nmat*j];
                    }
                }
            }

            //
            // Solving the (S) eigenproblem for this layer...
            //
            message(layerstr + "Computing s-pol eigenvalues and eigenvectors");
            eigen(MA,Qs,Ws,nmat);

            //
            // Creating the Q, X, and V matrices for S-pol...
            //
            for (i=0; i<nmat; ++i) {
                Qs[i] = sqrt(Qs[i]);
                Xs[i] = exp(-k0*Qs[i]*d);
            }
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    COMPLEX& _Vs = Vs[i+nmat*j];
                    _Vs = 0;
                    for (k=0; k<nmat; ++k) {
                        _Vs += Mxinv[i+nmat*k]*Ws[k+nmat*j];
                    }
                    _Vs *= Qs[j];
                }
            }

            //
            // Solving the (P) eigenproblem for this layer...
            //
            message(layerstr + "Computing p-pol eigenvalues and eigenvectors");
            eigen(EB,Qp,Wp,nmat);

            //
            // Creating the Q, X, and V matrices for P-pol...
            //
            for (i=0; i<nmat; ++i) {
                Qp[i] = sqrt(Qp[i]);
                Xp[i] = exp(-k0*Qp[i]*d);
            }
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    COMPLEX& _Vp = Vp[i+nmat*j];
                    _Vp = 0;
                    for (k=0; k<nmat; ++k) {
                        _Vp += Exinv[i+nmat*k]*Wp[k+nmat*j];
                    }
                    _Vp *= Qp[j];
                }
            }

            //
            // Building large matrices...
            //
            message(layerstr + "Creating large matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    Fs[  i   +nnmat*   j    ] = -Ws[ij];
                    Fs[  i   +nnmat*(j+nmat)] =  fs[ij];
                    Fs[i+nmat+nnmat*   j    ] =  Vs[ij];
                    Fs[i+nmat+nnmat*(j+nmat)] =  gs[ij];
                    Fp[  i   +nnmat*   j    ] = -Wp[ij];
                    Fp[  i   +nnmat*(j+nmat)] =  fp[ij];
                    Fp[i+nmat+nnmat*   j    ] =  Vp[ij];
                    Fp[i+nmat+nnmat*(j+nmat)] =  gp[ij];
                }
            }

            //
            // LU decomposing large matrices...
            //
            message(layerstr + "LU decomposing large matrices");
            LUdecompose(Fs,nnmat,indexs);
            LUdecompose(Fp,nnmat,indexp);

            message(layerstr + "LU backsubstitution");
            for (i=0; i<nmat; ++i) {
                // Do one column at a time (the i-th column)...
                for (j=0; j<nmat; ++j) {
                    WXVXs[  j   ] = Ws[j+nmat*i]*Xs[i];
                    WXVXs[j+nmat] = Vs[j+nmat*i]*Xs[i];
                    WXVXp[  j   ] = Wp[j+nmat*i]*Xp[i];
                    WXVXp[j+nmat] = Vp[j+nmat*i]*Xp[i];
                }

                LUbacksubstitute(Fs,nnmat,indexs,WXVXs);
                LUbacksubstitute(Fp,nnmat,indexp,WXVXp);
                for (j=0; j<nnmat; ++j) {
                    as[j+nnmat*i] = WXVXs[j];
                    ap[j+nnmat*i] = WXVXp[j];
                }
            }

            message(layerstr + "Calculating f and g matrices");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    fs[ij] = 0;
                    gs[ij] = 0;
                    fp[ij] = 0;
                    gp[ij] = 0;
                    for (k=0; k<nmat; ++k) {
                        COMPLEX I = (j==k) ? 1. : 0.;
                        COMPLEX Xsa = Xs[k]*as[k+nnmat*j];
                        COMPLEX Xpa = Xp[k]*ap[k+nnmat*j];
                        int ik = i+nmat*k;
                        fs[ij] += Ws[ik]*(I+Xsa);
                        gs[ij] += Vs[ik]*(I-Xsa);
                        fp[ij] += Wp[ik]*(I+Xpa);
                        gp[ij] += Vp[ik]*(I-Xpa);
                    }
                }
            }
        } // ...End of loop on layer

        message("Creating final matrices");
        vector<COMPLEX> Gs(nnmat*nnmat),Gp(nnmat*nnmat);
        for (i=0; i<nmat; ++i) {
            for (j=0; j<nmat; ++j) {
                COMPLEX I = (i==j) ? 1. : 0.;
                COMPLEX iZ,iY;
                if (backward) {
                    iZ = (i==j) ? cI*ZII[i] : 0.;
                    iY = (i==j) ? cI*YII[i] : 0.;
                } else {
                    iZ = (i==j) ? cI*ZI[i] : 0.;
                    iY = (i==j) ? cI*YI[i] : 0.;
                }
                int ij=i+nmat*j;
                Gs[  i   +nnmat*   j    ] = I;
                Gs[  i   +nnmat*(j+nmat)] = fs[ij];
                Gs[i+nmat+nnmat*   j    ] = -iY;
                Gs[i+nmat+nnmat*(j+nmat)] = gs[ij];
                Gp[  i   +nnmat*   j    ] = I;
                Gp[  i   +nnmat*(j+nmat)] = fp[ij];
                Gp[i+nmat+nnmat*   j    ] = -iZ;
                Gp[i+nmat+nnmat*(j+nmat)] = gp[ij];
            }
        }

        vector<COMPLEX> RRs(nnmat,0.);
        vector<COMPLEX> RRp(nnmat,0.);
        RRs[n] = -1.;
        //RRs[n+nmat] = COMPLEX(0.,-cos(thetai*deg)*nI);
        RRs[n+nmat] = COMPLEX(0.,inckz/k0);
        RRp[n] = -1.;
        if (backward) {
            RRp[n+nmat] = COMPLEX(0.,inckz/k0/sqr(real(nII)));
        } else {
            RRp[n+nmat] = COMPLEX(0.,inckz/k0/sqr(nI));
        }

        message("Solving for reflected fields");
        {
            vector<int> indexs(nnmat),indexp(nnmat);
            message("LU decomposition of Gs");
            LUdecompose(Gs,nnmat,indexs);
            message("LU decomposition of Gs");
            LUbacksubstitute(Gs,nnmat,indexs,RRs);
            message("LU decomposition of Gp");
            LUdecompose(Gp,nnmat,indexp);
            message("LU decomposition of Gp");
            LUbacksubstitute(Gp,nnmat,indexp,RRp);
        }

        message("Finishing");
        for (i=-order; i<=order; ++i) {
            j= n+i;
            r[j].PS() = 0.;
            r[j].PP() = conj(RRp[j]);
            r[j].SP() = 0.;
            r[j].SS() = conj(RRs[j]);
            if (backward) {
                R[j]=MuellerMatrix(r[j])*fabs(real(kIIzi[j])/inckz);
            } else {
                R[j]=MuellerMatrix(r[j])*fabs(real(kIzi[j])/inckz);
            }
        }
    }

    //***********************************************************************
    //**                                                                   **
    //** In-Plane Transmission                                             **
    //**                                                                   **
    //***********************************************************************
    void RCW_Model::InPlaneTransmission(bool backward)
    {
        int i,j,k;

        int nnmat = 2*nmat;

        vector<COMPLEX> epsinvx(2*nmat+1),epsy(2*nmat+1),epsz(2*nmat+1);
        vector<COMPLEX> muinvx(2*nmat+1),muy(2*nmat+1),muz(2*nmat+1);
        vector<COMPLEX> Ex(sqr(nmat)),Ey(sqr(nmat)),Ez(sqr(nmat));
        vector<COMPLEX> Exinv(sqr(nmat)),Eyinv(sqr(nmat)),Ezinv(sqr(nmat));
        vector<COMPLEX> Mx(sqr(nmat)),My(sqr(nmat)),Mz(sqr(nmat));
        vector<COMPLEX> Mxinv(sqr(nmat)),Myinv(sqr(nmat)),Mzinv(sqr(nmat));
        vector<COMPLEX> A(sqr(nmat));
        vector<COMPLEX> B(sqr(nmat));
        vector<COMPLEX> EB(sqr(nmat));
        vector<COMPLEX> MA(sqr(nmat));
        vector<COMPLEX> Ws(sqr(nmat)),Qs(nmat),Xs(nmat),Vs(sqr(nmat));
        vector<COMPLEX> Wp(sqr(nmat)),Qp(nmat),Xp(nmat),Vp(sqr(nmat));
        vector<COMPLEX> Fs(sqr(nnmat)),Fp(sqr(nnmat));
        vector<int> indexs(nnmat),indexp(nnmat);
        vector<COMPLEX> WXVXs(nnmat);
        vector<COMPLEX> WXVXp(nnmat);
        vector<COMPLEX> as(2*nmat*nmat),ap(2*nmat*nmat),as_(2*nmat*nmat),ap_(2*nmat*nmat);

        vector<COMPLEX> fs(nmat*nmat),gs(nmat*nmat),ss(nmat*nmat),ts(nmat*nmat);
        vector<COMPLEX> fp(nmat*nmat),gp(nmat*nmat),sp(nmat*nmat),tp(nmat*nmat);
        double costhetai = backward ? -inckz/k0/real(nII) : -inckz/k0/nI;

        for (i=0; i<nmat; ++i) {
            for (j=0; j<nmat; ++j) {
                int ij = i+nmat*j;
                if (backward) {
                    fs[ij] = (i==j) ? 1. : 0.;
                    gs[ij] = (i==j) ? cI*nII*costhetai : 0.;
                    ss[ij] = (i==j) ? 1. : 0.;
                    ts[ij] = (i==j) ? -cI*YII[i] : 0.;
                    fp[ij] = (i==j) ? 1. : 0.;
                    gp[ij] = (i==j) ? cI*costhetai/nII : 0.;
                    sp[ij] = (i==j) ? 1. : 0.;
                    tp[ij] = (i==j) ? -cI*ZII[i] : 0.;
                } else {
                    fs[ij] = (i==j) ? 1. : 0.;
                    gs[ij] = (i==j) ? cI*nI*costhetai : 0.;
                    ss[ij] = (i==j) ? 1. : 0.;
                    ts[ij] = (i==j) ? -cI*YI[i] : 0.;
                    fp[ij] = (i==j) ? 1. : 0.;
                    gp[ij] = (i==j) ? cI*costhetai/nI : 0.;
                    sp[ij] = (i==j) ? 1. : 0.;
                    tp[ij] = (i==j) ? -cI*ZI[i] : 0.;
                }
            }
        }

        forloopint layerloop(0,grating->get_levels()-1,backward ? -1 : +1 );
        for (forloopint::iterator layer = layerloop.begin(); layer!=layerloop.end(); ++layer) {

            string layerstr = format("(layer = %d)",(int)layer);

            double d = grating->get_thickness(layer);
            if (d<0.) error("Thickness of layer " + to_string((int)layer) + " = " + to_string(d) + " less than zero");

            //
            // Get the Fourier components for this layer...
            //
            message(layerstr + "Creating E and M matrices");
            for (i=-nmat; i<=nmat; ++i) {
                if (grating->is_anisotropic()) {
                    epsinvx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    epsy[i+nmat]=conj(grating->fouriery(i,layer,0));
                    epsz[i+nmat]=conj(grating->fourierz(i,layer,0));
                } else {
                    epsinvx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    epsy[i+nmat]=conj(grating->fourierx(i,layer,0));
                    epsz[i+nmat]=epsy[i+nmat];
                }
                if (grating->is_magnetic()) {
                    muinvx[i+nmat]=conj(grating->fouriermux(i,layer,1));
                    muy[i+nmat]=conj(grating->fouriermuy(i,layer,0));
                    muz[i+nmat]=conj(grating->fouriermuz(i,layer,0));
                } else {
                    muinvx[i+nmat]=(i==0 ? 1.:0.);
                    muy[i+nmat]=(i==0 ? 1.:0.);
                    muz[i+nmat]=(i==0 ? 1.:0.);
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Ex[i+nmat*j] = Exinv[i+nmat*j] = epsinvx[i-j+nmat];
                    Ey[i+nmat*j] = Eyinv[i+nmat*j] = epsy[i-j+nmat];
                    Ez[i+nmat*j] = Ezinv[i+nmat*j] = epsz[i-j+nmat];
                    Mx[i+nmat*j] = Mxinv[i+nmat*j] = muinvx[i-j+nmat];
                    My[i+nmat*j] = Myinv[i+nmat*j] = muy[i-j+nmat];
                    Mz[i+nmat*j] = Mzinv[i+nmat*j] = muz[i-j+nmat];
                }
            }

            message(layerstr + "Inverting E and M matrices");
            Inverse(Ex,nmat);
            Inverse(Eyinv,nmat);
            Inverse(Ezinv,nmat);
            if (grating->is_magnetic()) {
                Inverse(Mx,nmat);
                Inverse(Myinv,nmat);
                Inverse(Mzinv,nmat);
            }

            //
            // Building matrix A...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    A[i+nmat*j] = Kx[i] * Mzinv[i+nmat*j] * Kx[j] - Ey[i+nmat*j];
                }
            }

            //
            // Building matrix B...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    // Eq. (36) ...
                    B[i+nmat*j] = Kx[i] * Ezinv[i+nmat*j] * Kx[j] - My[i+nmat*j];
                }
            }

            //
            // Building matrix M*A...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    COMPLEX& _MA = MA[i+nmat*j];
                    _MA = 0;
                    for (int k=0; k<nmat; ++k) {
                        _MA += Mx[i+nmat*k]*A[k+nmat*j];
                    }
                }
            }

            //
            // Building matrix E*B...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    COMPLEX& _EB = EB[i+nmat*j];
                    _EB = 0;
                    for (int k=0; k<nmat; ++k) {
                        _EB += Ex[i+nmat*k]*B[k+nmat*j];
                    }
                }
            }

            //
            // Solving the (S) eigenproblem for this layer...
            //
            message(layerstr + "Computing s-pol eigenvalues and eigenvectors");
            eigen(MA,Qs,Ws,nmat);

            //
            // Creating the Q, X, and V matrices for S-pol...
            //
            for (i=0; i<nmat; ++i) {
                Qs[i] = sqrt(Qs[i]);
                Xs[i] = exp(-k0*Qs[i]*d);
            }
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Vs[i+nmat*j] = 0;
                    for (k=0; k<nmat; ++k) {
                        Vs[i+nmat*j] += Mxinv[i+nmat*k]*Ws[k+nmat*j];
                    }
                    Vs[i+nmat*j] *= Qs[j];
                }
            }

            //
            // Solving the (P) eigenproblem for this layer...
            //
            message(layerstr + "Computing p-pol eigenvalues and eigenvectors");
            eigen(EB,Qp,Wp,nmat);

            //
            // Creating the Q, X, and V matrices for P-pol...
            //
            for (i=0; i<nmat; ++i) {
                Qp[i] = sqrt(Qp[i]);
                Xp[i] = exp(-k0*Qp[i]*d);
            }
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Vp[i+nmat*j] = 0;
                    for (k=0; k<nmat; ++k) {
                        //if (oldmethod) Vp[i+nmat*j] += Einv[i+nmat*k]*Wp[k+nmat*j];
                        Vp[i+nmat*j] += Exinv[i+nmat*k]*Wp[k+nmat*j];
                    }
                    Vp[i+nmat*j] *= Qp[j];
                }
            }

            //
            // Building large matrices...
            //
            message(layerstr + "Creating large matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    Fs[  i   +nnmat*   j    ] = -ss[ij];
                    Fs[  i   +nnmat*(j+nmat)] =  Ws[ij];
                    Fs[i+nmat+nnmat*   j    ] = -ts[ij];
                    Fs[i+nmat+nnmat*(j+nmat)] =  Vs[ij];
                    Fp[  i   +nnmat*   j    ] = -sp[ij];
                    Fp[  i   +nnmat*(j+nmat)] =  Wp[ij];
                    Fp[i+nmat+nnmat*   j    ] = -tp[ij];
                    Fp[i+nmat+nnmat*(j+nmat)] =  Vp[ij];
                }
            }

            //
            // LU decomposing large matrices...
            //
            message(layerstr + "LU decomposing large matrices");

            LUdecompose(Fs,nnmat,indexs);
            LUdecompose(Fp,nnmat,indexp);

            message(layerstr + "LU backsubstitution");
            for (i=0; i<nmat; ++i) {
                // Do one column at a time (the i-th column)...
                for (j=0; j<nmat; ++j) {
                    int ji=j+nmat*i;
                    WXVXs[  j   ] = fs[ji];
                    WXVXs[j+nmat] = gs[ji];
                    WXVXp[  j   ] = fp[ji];
                    WXVXp[j+nmat] = gp[ji];
                }

                LUbacksubstitute(Fs,nnmat,indexs,WXVXs);
                LUbacksubstitute(Fp,nnmat,indexp,WXVXp);

                for (j=0; j<nnmat; ++j) {
                    as[j+nnmat*i] = WXVXs[j];
                    ap[j+nnmat*i] = WXVXp[j];
                }

                for (j=0; j<nmat; ++j) {
                    int ji=j+nmat*i;
                    WXVXs[  j   ] = -Ws[ji]*Xs[i];
                    WXVXs[j+nmat] =  Vs[ji]*Xs[i];
                    WXVXp[  j   ] = -Wp[ji]*Xp[i];
                    WXVXp[j+nmat] =  Vp[ji]*Xp[i];
                }

                LUbacksubstitute(Fs,nnmat,indexs,WXVXs);
                LUbacksubstitute(Fp,nnmat,indexp,WXVXp);

                for (j=0; j<nnmat; ++j) {
                    as_[j+nnmat*i] = WXVXs[j];
                    ap_[j+nnmat*i] = WXVXp[j];
                }
            }

            message(layerstr + "Calculating f, g, s, and t matrices");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    fs[ij] =  0;
                    gs[ij] =  0;
                    ss[ij] =  Ws[ij];
                    ts[ij] = -Vs[ij];
                    fp[ij] =  0;
                    gp[ij] =  0;
                    sp[ij] =  Wp[ij];
                    tp[ij] = -Vp[ij];
                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        fs[ij] += Ws[ik]*Xs[k]*as[k+nmat+nnmat*j];
                        gs[ij] += Vs[ik]*Xs[k]*as[k+nmat+nnmat*j];
                        ss[ij] += Ws[ik]*Xs[k]*as_[k+nmat+nnmat*j];
                        ts[ij] += Vs[ik]*Xs[k]*as_[k+nmat+nnmat*j];
                        fp[ij] += Wp[ik]*Xp[k]*ap[k+nmat+nnmat*j];
                        gp[ij] += Vp[ik]*Xp[k]*ap[k+nmat+nnmat*j];
                        sp[ij] += Wp[ik]*Xp[k]*ap_[k+nmat+nnmat*j];
                        tp[ij] += Vp[ik]*Xp[k]*ap_[k+nmat+nnmat*j];
                    }
                }
            }
        } // ...End of loop on layer

        message("Creating final matrices");
        vector<COMPLEX> Gs(nnmat*nnmat),Gp(nnmat*nnmat);
        for (i=0; i<nmat; ++i) {
            for (j=0; j<nmat; ++j) {
                int ij=i+nmat*j;
                COMPLEX I = (i==j) ? 1. : 0.;
                COMPLEX iZ,iY;
                if (backward) {
                    iZ = (i==j) ? cI*ZI[i] : 0.;
                    iY = (i==j) ? cI*YI[i] : 0.;
                } else {
                    iZ = (i==j) ? cI*ZII[i] : 0.;
                    iY = (i==j) ? cI*YII[i] : 0.;
                }

                Gs[  i   +nnmat*   j    ] = -ss[ij];
                Gs[  i   +nnmat*(j+nmat)] = I;
                Gs[i+nmat+nnmat*   j    ] = -ts[ij];
                Gs[i+nmat+nnmat*(j+nmat)] = iY;
                Gp[  i   +nnmat*   j    ] = -sp[ij];
                Gp[  i   +nnmat*(j+nmat)] = I;
                Gp[i+nmat+nnmat*   j    ] = -tp[ij];
                Gp[i+nmat+nnmat*(j+nmat)] = iZ;
            }
        }

        vector<COMPLEX> Ts(nnmat,0.);
        vector<COMPLEX> Tp(nnmat,0.);

        for (i=0; i<nmat; ++i) {
            int in = i+nmat*n;
            Ts[i]=fs[in];
            Ts[i+nmat]=gs[in];
            Tp[i]=fp[in];
            Tp[i+nmat]=gp[in];
        }

        message("Solving for transmitted fields");
        {
            vector<int> indexs(nnmat),indexp(nnmat);

            LUdecompose(Gs,nnmat,indexs);
            LUbacksubstitute(Gs,nnmat,indexs,Ts);

            LUdecompose(Gp,nnmat,indexp);
            LUbacksubstitute(Gp,nnmat,indexp,Tp);
        }

        message("Finishing");
        for (i=-order; i<=order; ++i) {
            j= n+i;
            if (backward) {
                r[j].PS() = 0;
                r[j].PP() = conj(Tp[j+nmat]/nI*nII);
                r[j].SP() = 0.;
                r[j].SS() = conj(Ts[j+nmat]);
                R[j]=MuellerMatrix(r[j])*fabs(real(kIzi[j])/inckz);
            } else {
                r[j].PP() = conj(Tp[j+nmat]/nII*nI);
                r[j].SP() = 0.;
                r[j].SS() = conj(Ts[j+nmat]);
                R[j]=MuellerMatrix(r[j])*fabs(real(kIIzi[j])/inckz);
            }
        }
    }

    //***********************************************************************
    //**                                                                   **
    //** Conical Reflection for Isotropic Media                            **
    //**                                                                   **
    //***********************************************************************
    void RCW_Model::ConicalReflection(bool backward)
    {
        int i,j,k;

        int nnmat = 4*nmat;
        int nmat2 = 2*nmat;
        int nmat4 = 4*nmat;

        vector<COMPLEX> EEE(2*nmat+1),EEEE(2*nmat+1);
        vector<COMPLEX> E(sqr(nmat)), Einv(sqr(nmat));
        vector<COMPLEX> EE(sqr(nmat)), EEinv(sqr(nmat));
        vector<COMPLEX> MMM(2*nmat+1),MMMM(2*nmat+1);
        vector<COMPLEX> M(sqr(nmat)), Minv(sqr(nmat));
        vector<COMPLEX> MM(sqr(nmat)), MMinv(sqr(nmat));
        vector<COMPLEX> A(sqr(nmat)), Ainv(sqr(nmat));
        vector<COMPLEX> B(sqr(nmat)), Binv(sqr(nmat));
        vector<COMPLEX> Matrix1(sqr(nmat)),Matrix2(sqr(nmat));
        vector<COMPLEX> W1(sqr(nmat)),Q1(nmat),X1(nmat);
        vector<COMPLEX> W2(sqr(nmat)),Q2(nmat),X2(nmat);
        vector<COMPLEX> V11(sqr(nmat)),V12(sqr(nmat)),V21(sqr(nmat)),V22(sqr(nmat)),v21(sqr(nmat)),v12(sqr(nmat));
        vector<COMPLEX> Vsp(sqr(nmat)),Wsp(sqr(nmat)),Wpp(sqr(nmat)),Vpp(sqr(nmat)),
               Vss(sqr(nmat)),Wss(sqr(nmat)),Wps(sqr(nmat)),Vps(sqr(nmat));
        vector<COMPLEX> BigMat(sqr(nnmat)),BigMatLU(sqr(nnmat));
        vector<int>     index(nnmat);
        vector<COMPLEX> WXVXl(nnmat),WXVXr(nnmat),WXVXlx(nnmat),WXVXrx(nnmat);
        vector<COMPLEX> ab(nmat2*nmat4), fg(nmat2*nmat4);
        vector<double>  Fc(nmat),Fs(nmat);

        //
        // Set the initial values for the fg matrix...
        //
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                fg[i+nmat4*j]=0;
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                fg[    i     +nnmat*   i    ]=1;
                fg[  i+nmat  +nnmat*   i    ]=cI*YI[i];
                fg[(i+2*nmat)+nnmat*(i+nmat)]=1;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=cI*ZI[i];
            } else {
                fg[    i     +nnmat*   i    ]=1;
                fg[  i+nmat  +nnmat*   i    ]=cI*YII[i];
                fg[(i+2*nmat)+nnmat*(i+nmat)]=1;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=cI*ZII[i];
            }
        }

        //
        // Iterate through each layer...
        //
        forloopint layerloop(0,grating->get_levels()-1,backward ? +1 : -1 );
        for (forloopint::iterator layer = layerloop.begin(); layer!=layerloop.end(); ++layer) {


            string layerstr = format("(layer = %d)",(int)layer);

            double d = grating->get_thickness(layer);
            if (d<0.) error("Thickness of layer " + to_string((int)layer) + " = " + to_string(d) + " less than zero");

            //
            // Get the Fourier components for this layer...
            //
            for (i=-nmat; i<=nmat; ++i) {
                EEE[i+nmat]=conj(grating->fourierx(i,layer,0));
                EEEE[i+nmat]=conj(grating->fourierx(i,layer,1));
                MMM[i+nmat]=conj(grating->fouriermux(i,layer,0));
                MMMM[i+nmat]=conj(grating->fouriermux(i,layer,1));
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Einv[i+nmat*j] = E[i+nmat*j] = EEE[i-j+nmat];
                    EEinv[i+nmat*j] = EE[i+nmat*j] = EEEE[i-j+nmat];
                    Minv[i+nmat*j] = M[i+nmat*j] = MMM[i-j+nmat];
                    MMinv[i+nmat*j] = MM[i+nmat*j] = MMMM[i-j+nmat];
                }
            }

            message(layerstr + "Inverting E matrix");
            Inverse(Einv,nmat);
            Inverse(EEinv,nmat);
            if (grating->is_magnetic()) {
                Inverse(Minv,nmat);
                Inverse(MMinv,nmat);
            }

            //
            // Building matrix A...
            //
            message(layerstr + "Building and inverting A matrix");
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    Ainv[ij] = A[ij] = Kx[i]*Minv[ij]*Kx[j] - E[ij];
                }
            }
            Inverse(Ainv,nmat);

            //
            // Building matrix B...
            //
            message(layerstr + "Building and inverting B matrix");
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    // Eq. (36) ...
                    int ij = i+nmat*j;
                    Binv[ij] = B[ij] = Kx[i]*Einv[ij]*Kx[j] - M[ij];
                }
            }
            Inverse(Binv,nmat);

            //
            // Building matrices ky^2 I + A  and ky^2 I + BE ...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    Matrix1[ij] = 0;
                    Matrix2[ij] = 0;
                    for (int k=0; k<nmat; ++k) {
                        int ik=i+nmat*k;
                        int kj=k+nmat*j;
                        Matrix1[ij] += A[ik]*MMinv[kj];
                        Matrix2[ij] += B[ik]*EEinv[kj];
                    }
                }
                int ii=i+nmat*i;
                Matrix1[ii] += sqr(ky/k0);
                Matrix2[ii] += sqr(ky/k0);
            }

            //
            // Calculate eigenvalues and eigenvectors...
            //
            message(layerstr + "Computing eigenvalues and eigenvectors");
            eigen(Matrix1,Q1,W1,nmat);
            eigen(Matrix2,Q2,W2,nmat);

            for (i=0; i<nmat; ++i) {
                Q1[i] = sqrt(Q1[i]);
                X1[i] = exp(-k0*Q1[i]*d);
                Q2[i] = sqrt(Q2[i]);
                X2[i] = exp(-k0*Q2[i]*d);
            }


            message(layerstr + "Making V matrices");
            // Eqs. (65)...
            for (j=0; j<nmat; ++j) {
                for (i=0; i<nmat; ++i) {
                    int ij = i+nmat*j;
                    V11[ij] = 0;
                    v12[ij] = 0;
                    v21[ij] = 0;
                    V22[ij] = 0;
                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        int kj = k+nmat*j;
                        V11[ij] += Ainv[ik]*W1[kj];
                        v12[ij] += Ainv[ik]*Kx[k]*Minv[kj];
                        v21[ij] += Binv[ik]*Kx[k]*Einv[kj];
                        V22[ij] += Binv[ik]*W2[kj];
                    }
                    V11[ij] *= Q1[j];
                    v12[ij] *= ky/k0;
                    v21[ij] *= ky/k0;
                    V22[ij] *= Q2[j];
                }
            }
            for (j=0; j<nmat; ++j) {
                for (i=0; i<nmat; ++i) {
                    int ij = i+nmat*j;
                    V12[ij] = 0;
                    V21[ij] = 0;
                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        int kj = k+nmat*j;
                        V12[ij] += v12[ik]*W2[kj];
                        V21[ij] += v21[ik]*W1[kj];
                    }
                }
            }

            for (i=0; i<nmat; ++i) {
                double phii = atan2(ky,kxi[i]); //atan(ky/kxi[i]);
                Fc[i] = cos(phii);
                Fs[i] = sin(phii);
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    Vsp[ij] = Fc[i]*V12[ij]-Fs[i]*W2[ij];
                    Wsp[ij] = Fs[i]*V22[ij];
                    Wpp[ij] = Fc[i]*V22[ij];
                    Vpp[ij] = Fs[i]*V12[ij]+Fc[i]*W2[ij];
                    Vss[ij] = Fc[i]*V11[ij];
                    Wss[ij] = Fs[i]*V21[ij]+Fc[i]*W1[ij];
                    Wps[ij] = Fc[i]*V21[ij]-Fs[i]*W1[ij];
                    Vps[ij] = Fs[i]*V11[ij];
                }
            }

            //
            // Building large matrices...
            //
            message(layerstr + "Building large matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    BigMat[    i     +nnmat*     j    ] = -Vss[ij];
                    BigMat[    i     +nnmat* (j+nmat) ] = -Vsp[ij];
                    BigMat[ (i+nmat) +nnmat*     j    ] =  Wss[ij];
                    BigMat[ (i+nmat) +nnmat* (j+nmat) ] =  Wsp[ij];
                    BigMat[(i+2*nmat)+nnmat*     j    ] =  Wps[ij];
                    BigMat[(i+2*nmat)+nnmat* (j+nmat) ] =  Wpp[ij];
                    BigMat[(i+3*nmat)+nnmat*     j    ] = -Vps[ij];
                    BigMat[(i+3*nmat)+nnmat* (j+nmat) ] = -Vpp[ij];
                }
            }
            for (i=0; i<nmat4; ++i) {
                for (j=0; j<nmat2; ++j) {
                    BigMat[i+nnmat*(j+2*nmat)] = fg[i+nnmat*j];
                }
            }

            //
            // LU decomposing large matrices...
            //
            message(layerstr + "LU decomposing large matrix");
            LUdecompose(BigMat,nnmat,index);

            message(layerstr + "Creating ab matrix by LU backsubstitution");
            for (j=0; j<nmat; ++j) {
                // Do one column at a time (the j-th column)...
                for (i=0; i<nmat; ++i) {
                    int ij = i+nmat*j;
                    WXVXl[   i    ] = Vss[ij]*X1[j];
                    WXVXl[ i+nmat ] = Wss[ij]*X1[j];
                    WXVXl[i+nmat*2] = Wps[ij]*X1[j];
                    WXVXl[i+nmat*3] = Vps[ij]*X1[j];

                    WXVXr[   i    ] = Vsp[ij]*X2[j];
                    WXVXr[ i+nmat ] = Wsp[ij]*X2[j];
                    WXVXr[i+nmat*2] = Wpp[ij]*X2[j];
                    WXVXr[i+nmat*3] = Vpp[ij]*X2[j];
                }
                LUbacksubstitute(BigMat,nnmat,index,WXVXl);
                LUbacksubstitute(BigMat,nnmat,index,WXVXr);
                for (i=0; i<nmat4; ++i) {
                    ab[i+nmat4*   j    ] = WXVXl[i];
                    ab[i+nmat4*(j+nmat)] = WXVXr[i];
                }
            }

            //
            // Calculating new fg matrix...
            //
            message(layerstr + "Calculating new fg matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    fg[     i    +nnmat*   j    ] = Vss[ij];
                    fg[  i+nmat  +nnmat*   j    ] = Wss[ij];
                    fg[(i+nmat*2)+nnmat*   j    ] = Wps[ij];
                    fg[(i+nmat*3)+nnmat*   j    ] = Vps[ij];
                    fg[     i    +nnmat*(nmat+j)] = Vsp[ij];
                    fg[  i+nmat  +nnmat*(nmat+j)] = Wsp[ij];
                    fg[(i+nmat*2)+nnmat*(nmat+j)] = Wpp[ij];
                    fg[(i+nmat*3)+nnmat*(nmat+j)] = Vpp[ij];

                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        COMPLEX a11 = ab[  k    +nnmat*   j    ];
                        COMPLEX a12 = ab[  k    +nnmat*(j+nmat)];
                        COMPLEX a21 = ab[k+nmat +nnmat*   j    ];
                        COMPLEX a22 = ab[k+nmat +nnmat*(j+nmat)];

                        fg[     i    +nnmat*   j    ] +=  Vss[ik]*X1[k]*a11+Vsp[ik]*X2[k]*a21;
                        fg[  i+nmat  +nnmat*   j    ] += -Wss[ik]*X1[k]*a11-Wsp[ik]*X2[k]*a21;
                        fg[(i+nmat*2)+nnmat*   j    ] += -Wps[ik]*X1[k]*a11-Wpp[ik]*X2[k]*a21;
                        fg[(i+nmat*3)+nnmat*   j    ] +=  Vps[ik]*X1[k]*a11+Vpp[ik]*X2[k]*a21;
                        fg[     i    +nnmat*(nmat+j)] +=  Vss[ik]*X1[k]*a12+Vsp[ik]*X2[k]*a22;
                        fg[  i+nmat  +nnmat*(nmat+j)] += -Wss[ik]*X1[k]*a12-Wsp[ik]*X2[k]*a22;
                        fg[(i+nmat*2)+nnmat*(nmat+j)] += -Wps[ik]*X1[k]*a12-Wpp[ik]*X2[k]*a22;
                        fg[(i+nmat*3)+nnmat*(nmat+j)] +=  Vps[ik]*X1[k]*a12+Vpp[ik]*X2[k]*a22;
                    }
                }
            }
        } // ...End of loop on layer

        message("Creating final large matrix");
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                BigMat[i+nnmat*    j     ] = 0;
                BigMat[i+nnmat*(j+2*nmat)] = fg[i+nnmat*j];
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                BigMat[    i      +nnmat*   i    ]=-1;
                BigMat[ (i+nmat)  +nnmat*   i    ]=cI*YII[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat)]=-1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat)]=cI*ZII[i];
            } else {
                BigMat[    i      +nnmat*   i    ]=-1;
                BigMat[ (i+nmat)  +nnmat*   i    ]=cI*YI[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat)]=-1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat)]=cI*ZI[i];
            }
        }

        message("Performing final LU decomposition");
        LUdecompose(BigMat,nnmat,index);

        vector<COMPLEX> RRs(nnmat,0.);
        vector<COMPLEX> RRp(nnmat,0.);

        if (backward) {
            double costhetai = -inckz/k0/real(nII);
            RRs[   n    ] = -1.;
            RRs[ n+nmat ] = -cI*nII*costhetai;
            RRp[n+2*nmat] = -cI*nII;
            RRp[n+3*nmat] = costhetai;
        } else {
            double costhetai = -inckz/k0/nI;
            RRs[   n    ] = -1.;
            RRs[ n+nmat ] = -cI*nI*costhetai;
            RRp[n+2*nmat] = -cI*nI;
            RRp[n+3*nmat] = costhetai;
        }

        message("Solving for fields");
        LUbacksubstitute(BigMat,nnmat,index,RRs);
        LUbacksubstitute(BigMat,nnmat,index,RRp);

        message("Finishing");
        for (i=-order; i<=order; ++i) {
            j= n+i;
            // Conventions for incident polarizations are handled above.
            // Conventions for reflected polarizations need factor of -1 for S
            // and sqrt(-1) for P.  Conjugates are needed for conversion to -iwt from iwt.
            if (backward) {
                r[j].PS() = conj(   RRp[j] );
                r[j].PP() = conj(  cI*RRp[j+nmat]/nII );
                r[j].SP() = conj( -cI*RRs[j+nmat]/nII );
                r[j].SS() = conj(   -RRs[j] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIIzi[j])/inckz);
            } else {
                r[j].PS() = conj(   -RRp[j] );
                r[j].PP() = conj( cI*RRp[j+nmat]/nI );
                r[j].SP() = conj( cI*RRs[j+nmat]/nI );
                r[j].SS() = conj(   -RRs[j] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIzi[j])/inckz);
            }
        }
    }

    //***********************************************************************
    //**                                                                   **
    //** Conical Transmission for Isotropic Media                          **
    //**                                                                   **
    //***********************************************************************
    void RCW_Model::ConicalTransmission(bool backward)
    {
        int i,j,k;

        int nnmat = 4*nmat;
        int nmat2 = 2*nmat;
        int nmat4 = 4*nmat;

        vector<COMPLEX> EEE(2*nmat+1),EEEE(2*nmat+1);
        vector<COMPLEX> E(sqr(nmat)), Einv(sqr(nmat));
        vector<COMPLEX> EE(sqr(nmat)), EEinv(sqr(nmat));
        vector<COMPLEX> MMM(2*nmat+1),MMMM(2*nmat+1);
        vector<COMPLEX> M(sqr(nmat)), Minv(sqr(nmat));
        vector<COMPLEX> MM(sqr(nmat)), MMinv(sqr(nmat));
        vector<COMPLEX> A(sqr(nmat)), Ainv(sqr(nmat));
        vector<COMPLEX> B(sqr(nmat)), Binv(sqr(nmat));
        vector<COMPLEX> Matrix1(sqr(nmat)),Matrix2(sqr(nmat));
        vector<COMPLEX> W1(sqr(nmat)),Q1(nmat),X1(nmat);
        vector<COMPLEX> W2(sqr(nmat)),Q2(nmat),X2(nmat);
        vector<COMPLEX> V11(sqr(nmat)),V12(sqr(nmat)),V21(sqr(nmat)),V22(sqr(nmat)),v21(sqr(nmat)),v12(sqr(nmat));
        vector<COMPLEX> Vsp(sqr(nmat)),Wsp(sqr(nmat)),Wpp(sqr(nmat)),Vpp(sqr(nmat)),
               Vss(sqr(nmat)),Wss(sqr(nmat)),Wps(sqr(nmat)),Vps(sqr(nmat));
        vector<COMPLEX> BigMat(sqr(nnmat)),BigMat2(sqr(nnmat));
        vector<int>     index(nnmat);
        vector<COMPLEX> WXVX(nnmat);
        vector<COMPLEX> ab(nmat4*nmat4), fg(nmat2*nmat4),st(nmat2*nmat4);
        vector<double>  Fc(nmat),Fs(nmat);

        //
        // Set the initial values for the fg matrix...
        //
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                fg[i+nmat4*j]=0;
                st[i+nmat4*j]=0;
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                double costhetai = -inckz/k0/real(nII);
                fg[    i     +nnmat*   i    ]=-1;
                fg[  i+nmat  +nnmat*   i    ]=-cI*nII*costhetai;
                fg[(i+2*nmat)+nnmat*(i+nmat)]=-cI*nII;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=costhetai;
                st[    i     +nnmat*   i    ]=1;
                st[  i+nmat  +nnmat*   i    ]=-cI*YII[i];
                st[(i+2*nmat)+nnmat*(i+nmat)]=1;
                st[(i+3*nmat)+nnmat*(i+nmat)]=-cI*ZII[i];
            } else {
                double costhetai = -inckz/k0/nI;
                fg[    i     +nnmat*   i    ]=-1;
                fg[  i+nmat  +nnmat*   i    ]=-cI*nI*costhetai;
                fg[(i+2*nmat)+nnmat*(i+nmat)]=-cI*nI;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=costhetai;
                st[    i     +nnmat*   i    ]=1;
                st[  i+nmat  +nnmat*   i    ]=-cI*YI[i];
                st[(i+2*nmat)+nnmat*(i+nmat)]=1;
                st[(i+3*nmat)+nnmat*(i+nmat)]=-cI*ZI[i];
            }
        }

        //
        // Iterate through each layer...
        //
        forloopint layerloop(0,grating->get_levels()-1,backward ? -1 : +1 );
        for (forloopint::iterator layer = layerloop.begin(); layer!=layerloop.end(); ++layer) {

            string layerstr = format("(layer = %d)",(int)layer);

            double d = grating->get_thickness(layer);
            if (d<0.) error("Thickness of layer " + to_string((int)layer) + " = " + to_string(d) + " less than zero");

            //
            // Get the Fourier components for this layer...
            //
            for (i=-nmat; i<=nmat; ++i) {
                EEE[i+nmat]=conj(grating->fourierx(i,layer,0));
                EEEE[i+nmat]=conj(grating->fourierx(i,layer,1));
            }
            if (grating->is_magnetic()) {
                for (i=-nmat; i<=nmat; ++i) {
                    MMM[i+nmat]=conj(grating->fouriermux(i,layer,0));
                    MMMM[i+nmat]=conj(grating->fouriermux(i,layer,1));
                }
            } else {
                for (i=-nmat; i<=nmat; ++i) {
                    MMM[i+nmat]= (i==0 ? 1. : 0.);
                    MMMM[i+nmat]= (i==0 ? 1. : 0.);
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Einv[i+nmat*j] = E[i+nmat*j] = EEE[i-j+nmat];
                    EEinv[i+nmat*j] = EE[i+nmat*j] = EEEE[i-j+nmat];
                    Minv[i+nmat*j] = M[i+nmat*j] = MMM[i-j+nmat];
                    MMinv[i+nmat*j] = MM[i+nmat*j] = MMMM[i-j+nmat];
                }
            }

            message(layerstr + "Inverting E matrix");
            Inverse(Einv,nmat);
            Inverse(EEinv,nmat);

            if (grating->is_magnetic()) {
                Inverse(Minv,nmat);
                Inverse(MMinv,nmat);
            }

            //
            // Building matrix A...
            //
            message(layerstr + "Building and inverting A matrix");
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    Ainv[ij] = A[ij] = Kx[i]*Minv[ij]*Kx[j] - E[ij];
                }
            }
            Inverse(Ainv,nmat);

            //
            // Building matrix B...
            //
            message(layerstr + "Building and inverting B matrix");
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    // Eq. (36) ...
                    int ij = i+nmat*j;
                    Binv[ij] = B[ij] = Kx[i]*Einv[ij]*Kx[j] - M[ij];
                }
            }
            Inverse(Binv,nmat);

            //
            // Building matrices ky^2 I + A  and ky^2 I + BE ...
            //
            for (i=0; i<nmat; ++i) {
                for (int j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    Matrix1[ij] = 0;
                    Matrix2[ij] = 0;
                    for (int k=0; k<nmat; ++k) {
                        int ik=i+nmat*k;
                        int kj=k+nmat*j;
                        Matrix1[ij] += A[ik]*MMinv[kj];
                        Matrix2[ij] += B[ik]*EEinv[kj];
                    }
                }
                int ii=i+nmat*i;
                Matrix1[ii] += sqr(ky/k0);
                Matrix2[ii] += sqr(ky/k0);
            }

            //
            // Calculate eigenvalues and eigenvectors...
            //
            message(layerstr + "Computing eigenvalues and eigenvectors");
            eigen(Matrix1,Q1,W1,nmat);
            eigen(Matrix2,Q2,W2,nmat);

            for (i=0; i<nmat; ++i) {
                Q1[i] = sqrt(Q1[i]);
                X1[i] = exp(-k0*Q1[i]*d);
                Q2[i] = sqrt(Q2[i]);
                X2[i] = exp(-k0*Q2[i]*d);
            }

            message(layerstr + "Making V matrices");
            // Eqs. (65)...
            for (j=0; j<nmat; ++j) {
                for (i=0; i<nmat; ++i) {
                    int ij = i+nmat*j;
                    V11[ij] = 0;
                    v12[ij] = 0;
                    v21[ij] = 0;
                    V22[ij] = 0;
                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        int kj = k+nmat*j;
                        V11[ij] += Ainv[ik]*W1[kj];
                        v12[ij] += Ainv[ik]*Kx[k]*Minv[kj];
                        v21[ij] += Binv[ik]*Kx[k]*Einv[kj];
                        V22[ij] += Binv[ik]*W2[kj];
                    }
                    V11[ij] *= Q1[j];
                    v12[ij] *= ky/k0;
                    v21[ij] *= ky/k0;
                    V22[ij] *= Q2[j];
                }
            }
            for (j=0; j<nmat; ++j) {
                for (i=0; i<nmat; ++i) {
                    int ij = i+nmat*j;
                    V12[ij] = 0;
                    V21[ij] = 0;
                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        int kj = k+nmat*j;
                        V12[ij] += v12[ik]*W2[kj];
                        V21[ij] += v21[ik]*W1[kj];
                    }
                }
            }

            for (i=0; i<nmat; ++i) {
                double phii = atan2(ky,kxi[i]); //atan(ky/kxi[i]);
                Fc[i] = cos(phii);
                Fs[i] = sin(phii);
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    Vsp[ij] = Fc[i]*V12[ij]-Fs[i]*W2[ij];
                    Wsp[ij] = Fs[i]*V22[ij];
                    Wpp[ij] = Fc[i]*V22[ij];
                    Vpp[ij] = Fs[i]*V12[ij]+Fc[i]*W2[ij];
                    Vss[ij] = Fc[i]*V11[ij];
                    Wss[ij] = Fs[i]*V21[ij]+Fc[i]*W1[ij];
                    Wps[ij] = Fc[i]*V21[ij]-Fs[i]*W1[ij];
                    Vps[ij] = Fs[i]*V11[ij];
                }
            }

            //
            // Building large matrices...
            //
            message(layerstr + "Building large matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    BigMat[    i     +nnmat* (j+2*nmat)] = Vss[ij];
                    BigMat[    i     +nnmat* (j+3*nmat)] = Vsp[ij];
                    BigMat[ (i+nmat) +nnmat* (j+2*nmat)] = Wss[ij];
                    BigMat[ (i+nmat) +nnmat* (j+3*nmat)] = Wsp[ij];
                    BigMat[(i+2*nmat)+nnmat* (j+2*nmat)] = Wps[ij];
                    BigMat[(i+2*nmat)+nnmat* (j+3*nmat)] = Wpp[ij];
                    BigMat[(i+3*nmat)+nnmat* (j+2*nmat)] = Vps[ij];
                    BigMat[(i+3*nmat)+nnmat* (j+3*nmat)] = Vpp[ij];
                    BigMat2[    i     +nnmat* (j+2*nmat)] = -Vss[ij]*X1[j];
                    BigMat2[    i     +nnmat* (j+3*nmat)] = -Vsp[ij]*X2[j];
                    BigMat2[ (i+nmat) +nnmat* (j+2*nmat)] = Wss[ij]*X1[j];
                    BigMat2[ (i+nmat) +nnmat* (j+3*nmat)] = Wsp[ij]*X2[j];
                    BigMat2[(i+2*nmat)+nnmat* (j+2*nmat)] = Wps[ij]*X1[j];
                    BigMat2[(i+2*nmat)+nnmat* (j+3*nmat)] = Wpp[ij]*X2[j];
                    BigMat2[(i+3*nmat)+nnmat* (j+2*nmat)] = -Vps[ij]*X1[j];
                    BigMat2[(i+3*nmat)+nnmat* (j+3*nmat)] = -Vpp[ij]*X2[j];
                }
            }
            for (i=0; i<nmat4; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij = i+nnmat*j;
                    BigMat[ij] = -st[ij];
                    BigMat2[ij] = fg[ij];
                }
            }

            //
            // LU decomposing large matrix...
            //
            message(layerstr + "LU decomposing large matrix");
            LUdecompose(BigMat,nnmat,index);

            message(layerstr + "Creating ab matrix by LU backsubstitution");
            for (j=0; j<nnmat; ++j) {
                // Do one column at a time (the j-th column)...
                for (i=0; i<nnmat; ++i) {
                    int ij = i+nnmat*j;
                    WXVX[i] = BigMat2[ij];
                }
                LUbacksubstitute(BigMat,nnmat,index,WXVX);
                for (i=0; i<nnmat; ++i) {
                    int ij = i+nnmat*j;
                    ab[ij] = WXVX[i];
                }
            }

            //
            // Calculating new fg matrix...
            //
            message(layerstr + "Calculating new fg matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    fg[     i    +nnmat*   j    ] = 0;
                    fg[  i+nmat  +nnmat*   j    ] = 0;
                    fg[(i+nmat*2)+nnmat*   j    ] = 0;
                    fg[(i+nmat*3)+nnmat*   j    ] = 0;
                    fg[     i    +nnmat*(nmat+j)] = 0;
                    fg[  i+nmat  +nnmat*(nmat+j)] = 0;
                    fg[(i+nmat*2)+nnmat*(nmat+j)] = 0;
                    fg[(i+nmat*3)+nnmat*(nmat+j)] = 0;
                    st[     i    +nnmat*   j    ] =  Vss[ij];
                    st[  i+nmat  +nnmat*   j    ] = -Wss[ij];
                    st[(i+nmat*2)+nnmat*   j    ] = -Wps[ij];
                    st[(i+nmat*3)+nnmat*   j    ] =  Vps[ij];
                    st[     i    +nnmat*(nmat+j)] =  Vsp[ij];
                    st[  i+nmat  +nnmat*(nmat+j)] = -Wsp[ij];
                    st[(i+nmat*2)+nnmat*(nmat+j)] = -Wpp[ij];
                    st[(i+nmat*3)+nnmat*(nmat+j)] =  Vpp[ij];

                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        COMPLEX a31 = ab[(k+nmat*2)+nnmat*    j     ];
                        COMPLEX a32 = ab[(k+nmat*2)+nnmat* (j+nmat) ];
                        COMPLEX a33 = ab[(k+nmat*2)+nnmat*(j+nmat*2)];
                        COMPLEX a34 = ab[(k+nmat*2)+nnmat*(j+nmat*3)];
                        COMPLEX a41 = ab[(k+nmat*3)+nnmat*    j     ];
                        COMPLEX a42 = ab[(k+nmat*3)+nnmat* (j+nmat) ];
                        COMPLEX a43 = ab[(k+nmat*3)+nnmat*(j+nmat*2)];
                        COMPLEX a44 = ab[(k+nmat*3)+nnmat*(j+nmat*3)];

                        fg[     i    +nnmat*   j    ] += Vss[ik]*X1[k]*a31 + Vsp[ik]*X2[k]*a41;
                        fg[  i+nmat  +nnmat*   j    ] += Wss[ik]*X1[k]*a31 + Wsp[ik]*X2[k]*a41;
                        fg[(i+nmat*2)+nnmat*   j    ] += Wps[ik]*X1[k]*a31 + Wpp[ik]*X2[k]*a41;
                        fg[(i+nmat*3)+nnmat*   j    ] += Vps[ik]*X1[k]*a31 + Vpp[ik]*X2[k]*a41;
                        fg[     i    +nnmat*(nmat+j)] += Vss[ik]*X1[k]*a32 + Vsp[ik]*X2[k]*a42;
                        fg[  i+nmat  +nnmat*(nmat+j)] += Wss[ik]*X1[k]*a32 + Wsp[ik]*X2[k]*a42;
                        fg[(i+nmat*2)+nnmat*(nmat+j)] += Wps[ik]*X1[k]*a32 + Wpp[ik]*X2[k]*a42;
                        fg[(i+nmat*3)+nnmat*(nmat+j)] += Vps[ik]*X1[k]*a32 + Vpp[ik]*X2[k]*a42;
                        st[     i    +nnmat*   j    ] += Vss[ik]*X1[k]*a33 + Vsp[ik]*X2[k]*a43;
                        st[  i+nmat  +nnmat*   j    ] += Wss[ik]*X1[k]*a33 + Wsp[ik]*X2[k]*a43;
                        st[(i+nmat*2)+nnmat*   j    ] += Wps[ik]*X1[k]*a33 + Wpp[ik]*X2[k]*a43;
                        st[(i+nmat*3)+nnmat*   j    ] += Vps[ik]*X1[k]*a33 + Vpp[ik]*X2[k]*a43;
                        st[     i    +nnmat*(nmat+j)] += Vss[ik]*X1[k]*a34 + Vsp[ik]*X2[k]*a44;
                        st[  i+nmat  +nnmat*(nmat+j)] += Wss[ik]*X1[k]*a34 + Wsp[ik]*X2[k]*a44;
                        st[(i+nmat*2)+nnmat*(nmat+j)] += Wps[ik]*X1[k]*a34 + Wpp[ik]*X2[k]*a44;
                        st[(i+nmat*3)+nnmat*(nmat+j)] += Vps[ik]*X1[k]*a34 + Vpp[ik]*X2[k]*a44;
                    }
                }
            }
        } // ...End of loop on layer

        message("Creating final large matrix");
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                BigMat[i+nnmat*j] = -st[i+nnmat*j];
            }
            for (j=nmat2; j<nmat4; ++j) {
                BigMat[i+nnmat*j] = 0;
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                BigMat[    i      +nnmat*(i+nmat*2)] = 1;
                BigMat[ (i+nmat)  +nnmat*(i+nmat*2)] = cI*YI[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat*3)] = 1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat*3)] = cI*ZI[i];
            } else {
                BigMat[    i      +nnmat*(i+nmat*2)] = 1;
                BigMat[ (i+nmat)  +nnmat*(i+nmat*2)] = cI*YII[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat*3)] = 1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat*3)] = cI*ZII[i];
            }
        }

        message("Performing final LU decomposition");
        LUdecompose(BigMat,nnmat,index);

        vector<COMPLEX> Ts(nnmat);
        vector<COMPLEX> Tp(nnmat);

        for (i=0; i<nnmat; ++i) {
            int ins = i+nnmat*n;
            Ts[i]=fg[ins];
            int inp = i+nnmat*(n+nmat);
            Tp[i]=fg[inp];
        }

        message("Solving for fields");
        LUbacksubstitute(BigMat,nnmat,index,Ts);
        LUbacksubstitute(BigMat,nnmat,index,Tp);

        message("Finishing");
        for (i=-order; i<=order; ++i) {
            j= n+i;
            // Conventions for incident polarizations are handled above.
            // Conventions for reflected polarizations need factor of -1 for S
            // and sqrt(-1) for P.  Conjugates are needed for conversion to -iwt from iwt.
            if (backward) {
                r[j].PS() = conj(    Tp[j+2*nmat] );
                r[j].PP() = conj( cI*Tp[j+3*nmat]/nI );
                r[j].SP() = conj(-cI*Ts[j+3*nmat]/nI );
                r[j].SS() = conj(   -Ts[j+2*nmat] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIzi[j])/inckz);
            } else {
                r[j].PS() = conj(   -Tp[j+2*nmat] );
                r[j].PP() = conj( cI*Tp[j+3*nmat]/nII );
                r[j].SP() = conj( cI*Ts[j+3*nmat]/nII );
                r[j].SS() = conj(   -Ts[j+2*nmat] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIIzi[j])/inckz);
            }
        }
    }

    //***********************************************************************
    //**                                                                   **
    //** Conical Reflection for Anisotropic Media                          **
    //**                                                                   **
    //***********************************************************************
    void RCW_Model::ConicalAnisoReflection(bool backward)
    {
        int i,j,k;

        int nnmat = 4*nmat;
        int nmat2 = 2*nmat;
        int nmat4 = 4*nmat;

        vector<COMPLEX> EEEx(2*nmat+1),EEEEx(2*nmat+1);
        vector<COMPLEX> Ex(sqr(nmat)), Exinv(sqr(nmat));
        vector<COMPLEX> EEx(sqr(nmat)), EExinv(sqr(nmat));
        vector<COMPLEX> EEEy(2*nmat+1),EEEEy(2*nmat+1);
        vector<COMPLEX> Ey(sqr(nmat)), Eyinv(sqr(nmat));
        vector<COMPLEX> EEy(sqr(nmat)), EEyinv(sqr(nmat));
        vector<COMPLEX> EEEz(2*nmat+1),EEEEz(2*nmat+1);
        vector<COMPLEX> Ez(sqr(nmat)), Ezinv(sqr(nmat));
        vector<COMPLEX> EEz(sqr(nmat)), EEzinv(sqr(nmat));

        vector<COMPLEX> MMMx(2*nmat+1),MMMMx(2*nmat+1);
        vector<COMPLEX> Mx(sqr(nmat)), Mxinv(sqr(nmat));
        vector<COMPLEX> MMx(sqr(nmat)), MMxinv(sqr(nmat));
        vector<COMPLEX> MMMy(2*nmat+1),MMMMy(2*nmat+1);
        vector<COMPLEX> My(sqr(nmat)), Myinv(sqr(nmat));
        vector<COMPLEX> MMy(sqr(nmat)), MMyinv(sqr(nmat));
        vector<COMPLEX> MMMz(2*nmat+1),MMMMz(2*nmat+1);
        vector<COMPLEX> Mz(sqr(nmat)), Mzinv(sqr(nmat));
        vector<COMPLEX> MMz(sqr(nmat)), MMzinv(sqr(nmat));


        vector<COMPLEX> F(sqr(nmat2)), G(sqr(nmat2));
        vector<COMPLEX> FG(sqr(nmat2));
        vector<COMPLEX> W(sqr(nmat2)),Q(nmat2),v(sqr(nmat2));
        vector<COMPLEX> X(nmat2);

        vector<COMPLEX> Ws(nmat*nmat2);
        vector<COMPLEX> Wp(nmat*nmat2);
        vector<COMPLEX> Vs(nmat*nmat2);
        vector<COMPLEX> Vp(nmat*nmat2);

        vector<COMPLEX> X1(nmat);
        vector<COMPLEX> X2(nmat);
        vector<COMPLEX> Vs1(sqr(nmat)),Vs2(sqr(nmat)),Vp1(sqr(nmat)),Vp2(sqr(nmat)),
               Ws1(sqr(nmat)),Ws2(sqr(nmat)),Wp1(sqr(nmat)),Wp2(sqr(nmat));
        vector<COMPLEX> BigMat(sqr(nnmat)),BigMatLU(sqr(nnmat));
        vector<int>     index(nnmat);
        vector<COMPLEX> WXVXl(nnmat),WXVXr(nnmat),WXVXlx(nnmat),WXVXrx(nnmat);
        vector<COMPLEX> ab(nmat2*nmat4), fg(nmat2*nmat4);
        vector<double>  Fc(nmat),Fs(nmat);

        //
        // Set the initial values for the fg matrix...
        //
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                fg[i+nmat4*j]=0;
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                fg[    i     +nnmat*   i    ]=1;
                fg[  i+nmat  +nnmat*   i    ]=cI*YI[i];
                fg[(i+2*nmat)+nnmat*(i+nmat)]=1;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=cI*ZI[i];
            } else {
                fg[    i     +nnmat*   i    ]=1;
                fg[  i+nmat  +nnmat*   i    ]=cI*YII[i];
                fg[(i+2*nmat)+nnmat*(i+nmat)]=1;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=cI*ZII[i];
            }
        }

        //
        // Iterate through each layer...
        //
        forloopint layerloop(0,grating->get_levels()-1,backward ? +1 : -1 );
        for (forloopint::iterator layer = layerloop.begin(); layer!=layerloop.end(); ++layer) {

            string layerstr = format("(layer = %d)",(int)layer);

            double d = grating->get_thickness(layer);
            if (d<0.) error("Thickness of layer " + to_string((int)layer) + " = " + to_string(d) + " less than zero");

            //
            // Get the Fourier components for this layer...
            //
            for (i=-nmat; i<=nmat; ++i) {
                if (grating->is_anisotropic()) {
                    EEEx[i+nmat]=conj(grating->fourierx(i,layer,0));
                    EEEEx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    EEEy[i+nmat]=conj(grating->fouriery(i,layer,0));
                    EEEEy[i+nmat]=conj(grating->fouriery(i,layer,1));
                    EEEz[i+nmat]=conj(grating->fourierz(i,layer,0));
                    EEEEz[i+nmat]=conj(grating->fourierz(i,layer,1));
                } else {
                    EEEx[i+nmat]=conj(grating->fourierx(i,layer,0));
                    EEEEx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    EEEy[i+nmat]=EEEx[i+nmat];
                    EEEEy[i+nmat]=EEEEx[i+nmat];
                    EEEz[i+nmat]=EEEx[i+nmat];
                    EEEEz[i+nmat]=EEEEx[i+nmat];
                }
                if (grating->is_magnetic()) {
                    MMMx[i+nmat]=conj(grating->fouriermux(i,layer,0));
                    MMMMx[i+nmat]=conj(grating->fouriermux(i,layer,1));
                    MMMy[i+nmat]=conj(grating->fouriermuy(i,layer,0));
                    MMMMy[i+nmat]=conj(grating->fouriermuy(i,layer,1));
                    MMMz[i+nmat]=conj(grating->fouriermuz(i,layer,0));
                    MMMMz[i+nmat]=conj(grating->fouriermuz(i,layer,1));
                } else {
                    MMMx[i+nmat]=(i==0 ? 1.: 0.);
                    MMMMx[i+nmat]=(i==0 ? 1.: 0.);
                    MMMy[i+nmat]=(i==0 ? 1.: 0.);
                    MMMMy[i+nmat]=(i==0 ? 1.: 0.);
                    MMMz[i+nmat]=(i==0 ? 1.: 0.);
                    MMMMz[i+nmat]=(i==0 ? 1.: 0.);
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Exinv[i+nmat*j] = Ex[i+nmat*j] = EEEx[i-j+nmat];
                    EExinv[i+nmat*j] = EEx[i+nmat*j] = EEEEx[i-j+nmat];
                    Eyinv[i+nmat*j] = Ey[i+nmat*j] = EEEy[i-j+nmat];
                    EEyinv[i+nmat*j] = EEy[i+nmat*j] = EEEEy[i-j+nmat];
                    Ezinv[i+nmat*j] = Ez[i+nmat*j] = EEEz[i-j+nmat];
                    EEzinv[i+nmat*j] = EEz[i+nmat*j] = EEEEz[i-j+nmat];
                    Mxinv[i+nmat*j] = Mx[i+nmat*j] = MMMx[i-j+nmat];
                    MMxinv[i+nmat*j] = MMx[i+nmat*j] = MMMMx[i-j+nmat];
                    Myinv[i+nmat*j] = My[i+nmat*j] = MMMy[i-j+nmat];
                    MMyinv[i+nmat*j] = MMy[i+nmat*j] = MMMMy[i-j+nmat];
                    Mzinv[i+nmat*j] = Mz[i+nmat*j] = MMMz[i-j+nmat];
                    MMzinv[i+nmat*j] = MMz[i+nmat*j] = MMMMz[i-j+nmat];
                }
            }

            message(layerstr + "Inverting E and M matrices");
            Inverse(EExinv,nmat);
            Inverse(Eyinv,nmat);
            Inverse(Ezinv,nmat);

            if (grating->is_magnetic()) {
                Inverse(MMxinv,nmat);
                Inverse(Myinv,nmat);
                Inverse(Mzinv,nmat);
            }

            //
            // Building matrix R and S (the upper right and lower left off-diagonals of Eq. (57) of M&G)...
            //
            message(layerstr + "Building matrices F and G ");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    int ija = i+nmat2*j;
                    int ijb = i+nmat2*(j+nmat);
                    int ijc = (i+nmat)+nmat2*j;
                    int ijd = (i+nmat)+nmat2*(j+nmat);
                    F[ija] = (ky/k0)*Ezinv[ij]*Kx[j];
                    F[ijb] = MMxinv[ij] - sqr(ky/k0)*Ezinv[ij];
                    F[ijc] = Kx[i]*Ezinv[ij]*Kx[j] - My[ij];
                    F[ijd] = -Kx[i]*Ezinv[ij]*(ky/k0);
                    G[ija] = (ky/k0)*Mzinv[ij]*Kx[j];
                    G[ijb] = EExinv[ij] - sqr(ky/k0)*Mzinv[ij];
                    G[ijc] = Kx[i]*Mzinv[ij]*Kx[j] - Ey[ij];
                    G[ijd] = -Kx[i]*Mzinv[ij]*(ky/k0);
                }
            }

            message(layerstr + "Multiply F and G");
            for (i=0; i<nmat2; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij=i+nmat2*j;
                    FG[ij]=0;
                    for (k=0; k<nmat2; ++k) {
                        int ik=i+nmat2*k;
                        int kj=k+nmat2*j;
                        FG[ij] += F[ik]*G[kj];
                    }
                }
            }

            message(layerstr + "Invert F");
            Inverse(F,nmat2);

            message(layerstr + "Finding eigenvalues of FG");
            eigen(FG,Q,W,nmat2);

            for (i=0; i<nmat2; ++i) {
                Q[i] = sqrt(Q[i]);
                X[i] = exp(-k0*Q[i]*d);
            }
            for (i=0; i<nmat; ++i) {
                X1[i] = X[i];
                X2[i] = X[i+nmat];
            }

            message(layerstr + "Multiply inverse(F)*W*Q");
            for (i=0; i<nmat2; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij=i+nmat2*j;
                    v[ij]= 0;
                    for (k=0; k<nmat2; ++k) {
                        int ik=i+nmat2*k;
                        int kj=k+nmat2*j;
                        v[ij] += F[ik]*W[kj]*Q[j];
                        //v[ij] += G[ik]*W[kj]/Q[j];
                    }
                }
            }

            for (i=0; i<nmat; ++i) {
                double phii = atan2(ky,kxi[i]); //atan(ky/kxi[i]);
                Fc[i] = cos(phii);
                Fs[i] = sin(phii);
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij=i+nmat*j;
                    int yij = i+nmat2*j;
                    int xij = (i+nmat)+nmat2*j;
                    Ws[ij] = Fc[i]*W[yij]-Fs[i]*W[xij];
                    Wp[ij] = Fs[i]*W[yij]+Fc[i]*W[xij];
                    Vp[ij] = Fc[i]*v[yij]-Fs[i]*v[xij];
                    Vs[ij] = Fs[i]*v[yij]+Fc[i]*v[xij];
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    int ij1=i+nmat*j;
                    int ij2=i+nmat*(j+nmat);
                    Vs1[ij] = Vs[ij1];
                    Vs2[ij] = Vs[ij2];
                    Vp1[ij] = Vp[ij1];
                    Vp2[ij] = Vp[ij2];
                    Wp1[ij] = Wp[ij1];
                    Wp2[ij] = Wp[ij2];
                    Ws1[ij] = Ws[ij1];
                    Ws2[ij] = Ws[ij2];
                }
            }
            //
            // Building large matrices...
            //
            message(layerstr + "Building large matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    BigMat[    i     +nnmat*     j    ] = -Ws1[ij];
                    BigMat[    i     +nnmat* (j+nmat) ] = -Ws2[ij];
                    BigMat[ (i+nmat) +nnmat*     j    ] =  Vs1[ij];
                    BigMat[ (i+nmat) +nnmat* (j+nmat) ] =  Vs2[ij];
                    BigMat[(i+2*nmat)+nnmat*     j    ] =  Vp1[ij];
                    BigMat[(i+2*nmat)+nnmat* (j+nmat) ] =  Vp2[ij];
                    BigMat[(i+3*nmat)+nnmat*     j    ] =  -Wp1[ij];
                    BigMat[(i+3*nmat)+nnmat* (j+nmat) ] =  -Wp2[ij];
                }
            }
            for (i=0; i<nmat4; ++i) {
                for (j=0; j<nmat2; ++j) {
                    BigMat[i+nnmat*(j+2*nmat)] = fg[i+nnmat*j];
                }
            }


            //
            // LU decomposing large matrices...
            //
            message(layerstr + "LU decomposing large matrix");
            LUdecompose(BigMat,nnmat,index);

            message(layerstr + "Creating ab matrix by LU backsubstitution");
            for (j=0; j<nmat; ++j) {
                // Do one column at a time (the j-th column)...
                for (i=0; i<nmat; ++i) {
                    int ij = i+nmat*j;
                    WXVXl[   i    ] = Ws1[ij]*X1[j];
                    WXVXl[ i+nmat ] = Vs1[ij]*X1[j];
                    WXVXl[i+nmat*2] = Vp1[ij]*X1[j];
                    WXVXl[i+nmat*3] = Wp1[ij]*X1[j];

                    WXVXr[   i    ] = Ws2[ij]*X2[j];
                    WXVXr[ i+nmat ] = Vs2[ij]*X2[j];
                    WXVXr[i+nmat*2] = Vp2[ij]*X2[j];
                    WXVXr[i+nmat*3] = Wp2[ij]*X2[j];
                }
                LUbacksubstitute(BigMat,nnmat,index,WXVXl);
                LUbacksubstitute(BigMat,nnmat,index,WXVXr);
                for (i=0; i<nmat4; ++i) {
                    ab[i+nmat4*   j    ] = WXVXl[i];
                    ab[i+nmat4*(j+nmat)] = WXVXr[i];
                }
            }

            //
            // Calculating new fg matrix...
            //
            message(layerstr + "Calculating new fg matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    fg[     i    +nnmat*   j    ] = Ws1[ij];
                    fg[  i+nmat  +nnmat*   j    ] = Vs1[ij];
                    fg[(i+nmat*2)+nnmat*   j    ] = Vp1[ij];
                    fg[(i+nmat*3)+nnmat*   j    ] = Wp1[ij];
                    fg[     i    +nnmat*(nmat+j)] = Ws2[ij];
                    fg[  i+nmat  +nnmat*(nmat+j)] = Vs2[ij];
                    fg[(i+nmat*2)+nnmat*(nmat+j)] = Vp2[ij];
                    fg[(i+nmat*3)+nnmat*(nmat+j)] = Wp2[ij];

                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        COMPLEX a11 = ab[  k    +nnmat*   j    ];
                        COMPLEX a12 = ab[  k    +nnmat*(j+nmat)];
                        COMPLEX a21 = ab[k+nmat +nnmat*   j    ];
                        COMPLEX a22 = ab[k+nmat +nnmat*(j+nmat)];

                        fg[     i    +nnmat*   j    ] +=  Ws1[ik]*X1[k]*a11+Ws2[ik]*X2[k]*a21;
                        fg[  i+nmat  +nnmat*   j    ] += -Vs1[ik]*X1[k]*a11-Vs2[ik]*X2[k]*a21;
                        fg[(i+nmat*2)+nnmat*   j    ] += -Vp1[ik]*X1[k]*a11-Vp2[ik]*X2[k]*a21;
                        fg[(i+nmat*3)+nnmat*   j    ] +=  Wp1[ik]*X1[k]*a11+Wp2[ik]*X2[k]*a21;
                        fg[     i    +nnmat*(nmat+j)] +=  Ws1[ik]*X1[k]*a12+Ws2[ik]*X2[k]*a22;
                        fg[  i+nmat  +nnmat*(nmat+j)] += -Vs1[ik]*X1[k]*a12-Vs2[ik]*X2[k]*a22;
                        fg[(i+nmat*2)+nnmat*(nmat+j)] += -Vp1[ik]*X1[k]*a12-Vp2[ik]*X2[k]*a22;
                        fg[(i+nmat*3)+nnmat*(nmat+j)] +=  Wp1[ik]*X1[k]*a12+Wp2[ik]*X2[k]*a22;
                    }
                }
            }
        } // ...End of loop on layer

        message("Creating final large matrix");
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                BigMat[i+nnmat*    j     ] = 0;
                BigMat[i+nnmat*(j+2*nmat)] = fg[i+nnmat*j];
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                BigMat[    i      +nnmat*   i    ]=-1;
                BigMat[ (i+nmat)  +nnmat*   i    ]=cI*YII[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat)]=-1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat)]=cI*ZII[i];
            } else {
                BigMat[    i      +nnmat*   i    ]=-1;
                BigMat[ (i+nmat)  +nnmat*   i    ]=cI*YI[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat)]=-1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat)]=cI*ZI[i];
            }
        }

        message("Performing final LU decomposition");
        LUdecompose(BigMat,nnmat,index);

        vector<COMPLEX> RRs(nnmat,0.);
        vector<COMPLEX> RRp(nnmat,0.);

        if (backward) {
            double costhetai = -inckz/k0/real(nII);
            RRs[   n    ] = -1.;
            RRs[ n+nmat ] = -cI*nII*costhetai;
            RRp[n+2*nmat] = -cI*nII;
            RRp[n+3*nmat] = costhetai;
        } else {
            double costhetai = -inckz/k0/nI;
            RRs[   n    ] = -1.;
            RRs[ n+nmat ] = -cI*nI*costhetai;
            RRp[n+2*nmat] = -cI*nI;
            RRp[n+3*nmat] = costhetai;
        }

        message("Solving for fields");
        LUbacksubstitute(BigMat,nnmat,index,RRs);
        LUbacksubstitute(BigMat,nnmat,index,RRp);

        message("Finishing");
        for (i=-order; i<=order; ++i) {
            j= n+i;
            // Conventions for incident polarizations are handled above.
            // Conventions for reflected polarizations need factor of -1 for S
            // and sqrt(-1) for P.  Conjugates are needed for conversion to -iwt from iwt.
            if (backward) {
                r[j].PS() = conj(    RRp[j] );
                r[j].PP() = conj( cI*RRp[j+nmat]/nII );
                r[j].SP() = conj(-cI*RRs[j+nmat]/nII );
                r[j].SS() = conj(   -RRs[j] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIIzi[j])/inckz);
            } else {
                r[j].PS() = conj(   -RRp[j] );
                r[j].PP() = conj( cI*RRp[j+nmat]/nI );
                r[j].SP() = conj( cI*RRs[j+nmat]/nI );
                r[j].SS() = conj(   -RRs[j] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIzi[j])/inckz);
            }
        }
    }


    //***********************************************************************
    //**                                                                   **
    //** Conical Transmission for Anisotropic Media                        **
    //**                                                                   **
    //***********************************************************************
    void RCW_Model::ConicalAnisoTransmission(bool backward)
    {
        int i,j,k;

        int nnmat = 4*nmat;
        int nmat2 = 2*nmat;
        int nmat4 = 4*nmat;

        vector<COMPLEX> EEEx(2*nmat+1),EEEEx(2*nmat+1);
        vector<COMPLEX> Ex(sqr(nmat)), Exinv(sqr(nmat));
        vector<COMPLEX> EEx(sqr(nmat)), EExinv(sqr(nmat));
        vector<COMPLEX> EEEy(2*nmat+1),EEEEy(2*nmat+1);
        vector<COMPLEX> Ey(sqr(nmat)), Eyinv(sqr(nmat));
        vector<COMPLEX> EEy(sqr(nmat)), EEyinv(sqr(nmat));
        vector<COMPLEX> EEEz(2*nmat+1),EEEEz(2*nmat+1);
        vector<COMPLEX> Ez(sqr(nmat)), Ezinv(sqr(nmat));
        vector<COMPLEX> EEz(sqr(nmat)), EEzinv(sqr(nmat));

        vector<COMPLEX> MMMx(2*nmat+1),MMMMx(2*nmat+1);
        vector<COMPLEX> Mx(sqr(nmat)), Mxinv(sqr(nmat));
        vector<COMPLEX> MMx(sqr(nmat)), MMxinv(sqr(nmat));
        vector<COMPLEX> MMMy(2*nmat+1),MMMMy(2*nmat+1);
        vector<COMPLEX> My(sqr(nmat)), Myinv(sqr(nmat));
        vector<COMPLEX> MMy(sqr(nmat)), MMyinv(sqr(nmat));
        vector<COMPLEX> MMMz(2*nmat+1),MMMMz(2*nmat+1);
        vector<COMPLEX> Mz(sqr(nmat)), Mzinv(sqr(nmat));
        vector<COMPLEX> MMz(sqr(nmat)), MMzinv(sqr(nmat));


        vector<COMPLEX> F(sqr(nmat2)), G(sqr(nmat2));
        vector<COMPLEX> FG(sqr(nmat2));
        vector<COMPLEX> W(sqr(nmat2)),Q(nmat2),v(sqr(nmat2));
        vector<COMPLEX> X(nmat2);

        vector<COMPLEX> Ws(nmat*nmat2);
        vector<COMPLEX> Wp(nmat*nmat2);
        vector<COMPLEX> Vs(nmat*nmat2);
        vector<COMPLEX> Vp(nmat*nmat2);

        vector<COMPLEX> X1(nmat);
        vector<COMPLEX> X2(nmat);
        vector<COMPLEX> Vs1(sqr(nmat)),Vs2(sqr(nmat)),Vp1(sqr(nmat)),Vp2(sqr(nmat)),
               Ws1(sqr(nmat)),Ws2(sqr(nmat)),Wp1(sqr(nmat)),Wp2(sqr(nmat));

        //vector<COMPLEX> Matrix1(sqr(nmat)),Matrix2(sqr(nmat));
        vector<COMPLEX> BigMat(sqr(nnmat)),BigMat2(sqr(nnmat));
        vector<int>     index(nnmat);
        vector<COMPLEX> WXVX(nnmat);
        vector<COMPLEX> ab(nmat4*nmat4), fg(nmat2*nmat4),st(nmat2*nmat4);
        vector<double>  Fc(nmat),Fs(nmat);

        //
        // Set the initial values for the fg matrix...
        //
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                fg[i+nmat4*j]=0;
                st[i+nmat4*j]=0;
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                double costhetai = -inckz/k0/real(nII);
                fg[    i     +nnmat*   i    ]=-1;
                fg[  i+nmat  +nnmat*   i    ]=-cI*nII*costhetai;
                fg[(i+2*nmat)+nnmat*(i+nmat)]=-cI*nII;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=costhetai;
                st[    i     +nnmat*   i    ]=1;
                st[  i+nmat  +nnmat*   i    ]=-cI*YII[i];
                st[(i+2*nmat)+nnmat*(i+nmat)]=1;
                st[(i+3*nmat)+nnmat*(i+nmat)]=-cI*ZII[i];
            } else {
                double costhetai = -inckz/k0/nI;
                fg[    i     +nnmat*   i    ]=-1;
                fg[  i+nmat  +nnmat*   i    ]=-cI*nI*costhetai;
                fg[(i+2*nmat)+nnmat*(i+nmat)]=-cI*nI;
                fg[(i+3*nmat)+nnmat*(i+nmat)]=costhetai;
                st[    i     +nnmat*   i    ]=1;
                st[  i+nmat  +nnmat*   i    ]=-cI*YI[i];
                st[(i+2*nmat)+nnmat*(i+nmat)]=1;
                st[(i+3*nmat)+nnmat*(i+nmat)]=-cI*ZI[i];
            }
        }

        //
        // Iterate through each layer...
        //
        forloopint layerloop(0,grating->get_levels()-1,backward ? -1 : +1 );
        for (forloopint::iterator layer = layerloop.begin(); layer!=layerloop.end(); ++layer) {

            string layerstr = format("(layer = %d)",(int)layer);

            double d = grating->get_thickness(layer);
            if (d<0.) error("Thickness of layer " + to_string((int)layer) + " = " + to_string(d) + " less than zero");

            //
            // Get the Fourier components for this layer...
            //
            for (i=-nmat; i<=nmat; ++i) {
                if (grating->is_anisotropic()) {
                    EEEx[i+nmat]=conj(grating->fourierx(i,layer,0));
                    EEEEx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    EEEy[i+nmat]=conj(grating->fouriery(i,layer,0));
                    EEEEy[i+nmat]=conj(grating->fouriery(i,layer,1));
                    EEEz[i+nmat]=conj(grating->fourierz(i,layer,0));
                    EEEEz[i+nmat]=conj(grating->fourierz(i,layer,1));
                } else {
                    EEEx[i+nmat]=conj(grating->fourierx(i,layer,0));
                    EEEEx[i+nmat]=conj(grating->fourierx(i,layer,1));
                    EEEy[i+nmat]=EEEx[i+nmat];
                    EEEEy[i+nmat]=EEEEx[i+nmat];
                    EEEz[i+nmat]=EEEx[i+nmat];
                    EEEEz[i+nmat]=EEEEx[i+nmat];
                }
                if (grating->is_magnetic()) {
                    MMMx[i+nmat]=conj(grating->fouriermux(i,layer,0));
                    MMMMx[i+nmat]=conj(grating->fouriermux(i,layer,1));
                    MMMy[i+nmat]=conj(grating->fouriermuy(i,layer,0));
                    MMMMy[i+nmat]=conj(grating->fouriermuy(i,layer,1));
                    MMMz[i+nmat]=conj(grating->fouriermuz(i,layer,0));
                    MMMMz[i+nmat]=conj(grating->fouriermuz(i,layer,1));
                } else {
                    MMMx[i+nmat]=(i==0 ? 1.: 0.);
                    MMMMx[i+nmat]=(i==0 ? 1.: 0.);
                    MMMy[i+nmat]=(i==0 ? 1.: 0.);
                    MMMMy[i+nmat]=(i==0 ? 1.: 0.);
                    MMMz[i+nmat]=(i==0 ? 1.: 0.);
                    MMMMz[i+nmat]=(i==0 ? 1.: 0.);
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    Exinv[i+nmat*j] = Ex[i+nmat*j] = EEEx[i-j+nmat];
                    EExinv[i+nmat*j] = EEx[i+nmat*j] = EEEEx[i-j+nmat];
                    Eyinv[i+nmat*j] = Ey[i+nmat*j] = EEEy[i-j+nmat];
                    EEyinv[i+nmat*j] = EEy[i+nmat*j] = EEEEy[i-j+nmat];
                    Ezinv[i+nmat*j] = Ez[i+nmat*j] = EEEz[i-j+nmat];
                    EEzinv[i+nmat*j] = EEz[i+nmat*j] = EEEEz[i-j+nmat];
                    Mxinv[i+nmat*j] = Mx[i+nmat*j] = MMMx[i-j+nmat];
                    MMxinv[i+nmat*j] = MMx[i+nmat*j] = MMMMx[i-j+nmat];
                    Myinv[i+nmat*j] = My[i+nmat*j] = MMMy[i-j+nmat];
                    MMyinv[i+nmat*j] = MMy[i+nmat*j] = MMMMy[i-j+nmat];
                    Mzinv[i+nmat*j] = Mz[i+nmat*j] = MMMz[i-j+nmat];
                    MMzinv[i+nmat*j] = MMz[i+nmat*j] = MMMMz[i-j+nmat];
                }
            }

            message(layerstr + "Inverting E and M matrices");
            Inverse(EExinv,nmat);
            Inverse(Eyinv,nmat);
            Inverse(Ezinv,nmat);

            if (grating->is_magnetic()) {
                Inverse(MMxinv,nmat);
                Inverse(Myinv,nmat);
                Inverse(Mzinv,nmat);
            }

            //
            // Building matrix R and S (the upper right and lower left off-diagonals of Eq. (57) of M&G)...
            //
            message(layerstr + "Building matrices F and G ");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    int ija = i+nmat2*j;
                    int ijb = i+nmat2*(j+nmat);
                    int ijc = (i+nmat)+nmat2*j;
                    int ijd = (i+nmat)+nmat2*(j+nmat);
                    F[ija] = (ky/k0)*Ezinv[ij]*Kx[j];
                    F[ijb] = MMxinv[ij] - sqr(ky/k0)*Ezinv[ij];
                    F[ijc] = Kx[i]*Ezinv[ij]*Kx[j] - My[ij];
                    F[ijd] = -Kx[i]*Ezinv[ij]*(ky/k0);
                    G[ija] = (ky/k0)*Mzinv[ij]*Kx[j];
                    G[ijb] = EExinv[ij] - sqr(ky/k0)*Mzinv[ij];
                    G[ijc] = Kx[i]*Mzinv[ij]*Kx[j] - Ey[ij];
                    G[ijd] = -Kx[i]*Mzinv[ij]*(ky/k0);
                }
            }

            message(layerstr + "Multiply F and G");
            for (i=0; i<nmat2; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij=i+nmat2*j;
                    FG[ij]=0;
                    for (k=0; k<nmat2; ++k) {
                        int ik=i+nmat2*k;
                        int kj=k+nmat2*j;
                        FG[ij] += F[ik]*G[kj];
                    }
                }
            }

            message(layerstr + "Invert F");
            Inverse(F,nmat2);

            message(layerstr + "Finding eigenvalues of FG");
            eigen(FG,Q,W,nmat2);

            for (i=0; i<nmat2; ++i) {
                Q[i] = sqrt(Q[i]);
                X[i] = exp(-k0*Q[i]*d);
            }
            for (i=0; i<nmat; ++i) {
                X1[i] = X[i];
                X2[i] = X[i+nmat];
            }

            message(layerstr + "Multiply inverse(F)*W*Q");
            for (i=0; i<nmat2; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij=i+nmat2*j;
                    v[ij]= 0;
                    for (k=0; k<nmat2; ++k) {
                        int ik=i+nmat2*k;
                        int kj=k+nmat2*j;
                        v[ij] += F[ik]*W[kj]*Q[j];
                        //v[ij] += G[ik]*W[kj]/Q[j];
                    }
                }
            }

            for (i=0; i<nmat; ++i) {
                double phii = atan2(ky,kxi[i]); //atan(ky/kxi[i]);
                Fc[i] = cos(phii);
                Fs[i] = sin(phii);
            }


            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij=i+nmat*j;
                    int yij = i+nmat2*j;
                    int xij = (i+nmat)+nmat2*j;
                    Ws[ij] = Fc[i]*W[yij]-Fs[i]*W[xij];
                    Wp[ij] = Fs[i]*W[yij]+Fc[i]*W[xij];
                    Vp[ij] = Fc[i]*v[yij]-Fs[i]*v[xij];
                    Vs[ij] = Fs[i]*v[yij]+Fc[i]*v[xij];
                }
            }

            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij=i+nmat*j;
                    int ij1=i+nmat*j;
                    int ij2=i+nmat*(j+nmat);
                    Vs1[ij] = Vs[ij1];
                    Vs2[ij] = Vs[ij2];
                    Vp1[ij] = Vp[ij1];
                    Vp2[ij] = Vp[ij2];
                    Wp1[ij] = Wp[ij1];
                    Wp2[ij] = Wp[ij2];
                    Ws1[ij] = Ws[ij1];
                    Ws2[ij] = Ws[ij2];
                }
            }

            //
            // Building large matrices...
            //
            message(layerstr + "Building large matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    BigMat[    i     +nnmat* (j+2*nmat)] = Ws1[ij];
                    BigMat[    i     +nnmat* (j+3*nmat)] = Ws2[ij];
                    BigMat[ (i+nmat) +nnmat* (j+2*nmat)] = Vs1[ij];
                    BigMat[ (i+nmat) +nnmat* (j+3*nmat)] = Vs2[ij];
                    BigMat[(i+2*nmat)+nnmat* (j+2*nmat)] = Vp1[ij];
                    BigMat[(i+2*nmat)+nnmat* (j+3*nmat)] = Vp2[ij];
                    BigMat[(i+3*nmat)+nnmat* (j+2*nmat)] = Wp1[ij];
                    BigMat[(i+3*nmat)+nnmat* (j+3*nmat)] = Wp2[ij];
                    BigMat2[    i     +nnmat* (j+2*nmat)] = -Ws1[ij]*X1[j];
                    BigMat2[    i     +nnmat* (j+3*nmat)] = -Ws2[ij]*X2[j];
                    BigMat2[ (i+nmat) +nnmat* (j+2*nmat)] = Vs1[ij]*X1[j];
                    BigMat2[ (i+nmat) +nnmat* (j+3*nmat)] = Vs2[ij]*X2[j];
                    BigMat2[(i+2*nmat)+nnmat* (j+2*nmat)] = Vp1[ij]*X1[j];
                    BigMat2[(i+2*nmat)+nnmat* (j+3*nmat)] = Vp2[ij]*X2[j];
                    BigMat2[(i+3*nmat)+nnmat* (j+2*nmat)] = -Wp1[ij]*X1[j];
                    BigMat2[(i+3*nmat)+nnmat* (j+3*nmat)] = -Wp2[ij]*X2[j];
                }
            }
            for (i=0; i<nmat4; ++i) {
                for (j=0; j<nmat2; ++j) {
                    int ij = i+nnmat*j;
                    BigMat[ij] = -st[ij];
                    BigMat2[ij] = fg[ij];
                }
            }

            //
            // LU decomposing large matrix...
            //
            message(layerstr + "LU decomposing large matrix");
            LUdecompose(BigMat,nnmat,index);

            message(layerstr + "Creating ab matrix by LU backsubstitution");
            for (j=0; j<nnmat; ++j) {
                // Do one column at a time (the j-th column)...
                for (i=0; i<nnmat; ++i) {
                    int ij = i+nnmat*j;
                    WXVX[i] = BigMat2[ij];
                }
                LUbacksubstitute(BigMat,nnmat,index,WXVX);
                for (i=0; i<nnmat; ++i) {
                    int ij = i+nnmat*j;
                    ab[ij] = WXVX[i];
                }
            }

            //
            // Calculating new fg matrix...
            //
            message(layerstr + "Calculating new fg matrix");
            for (i=0; i<nmat; ++i) {
                for (j=0; j<nmat; ++j) {
                    int ij = i+nmat*j;
                    fg[     i    +nnmat*   j    ] = 0;
                    fg[  i+nmat  +nnmat*   j    ] = 0;
                    fg[(i+nmat*2)+nnmat*   j    ] = 0;
                    fg[(i+nmat*3)+nnmat*   j    ] = 0;
                    fg[     i    +nnmat*(nmat+j)] = 0;
                    fg[  i+nmat  +nnmat*(nmat+j)] = 0;
                    fg[(i+nmat*2)+nnmat*(nmat+j)] = 0;
                    fg[(i+nmat*3)+nnmat*(nmat+j)] = 0;
                    st[     i    +nnmat*   j    ] =  Ws1[ij];
                    st[  i+nmat  +nnmat*   j    ] = -Vs1[ij];
                    st[(i+nmat*2)+nnmat*   j    ] = -Vp1[ij];
                    st[(i+nmat*3)+nnmat*   j    ] =  Wp1[ij];
                    st[     i    +nnmat*(nmat+j)] =  Ws2[ij];
                    st[  i+nmat  +nnmat*(nmat+j)] = -Vs2[ij];
                    st[(i+nmat*2)+nnmat*(nmat+j)] = -Vp2[ij];
                    st[(i+nmat*3)+nnmat*(nmat+j)] =  Wp2[ij];

                    for (k=0; k<nmat; ++k) {
                        int ik = i+nmat*k;
                        COMPLEX a31 = ab[(k+nmat*2)+nnmat*    j     ];
                        COMPLEX a32 = ab[(k+nmat*2)+nnmat* (j+nmat) ];
                        COMPLEX a33 = ab[(k+nmat*2)+nnmat*(j+nmat*2)];
                        COMPLEX a34 = ab[(k+nmat*2)+nnmat*(j+nmat*3)];
                        COMPLEX a41 = ab[(k+nmat*3)+nnmat*    j     ];
                        COMPLEX a42 = ab[(k+nmat*3)+nnmat* (j+nmat) ];
                        COMPLEX a43 = ab[(k+nmat*3)+nnmat*(j+nmat*2)];
                        COMPLEX a44 = ab[(k+nmat*3)+nnmat*(j+nmat*3)];

                        fg[     i    +nnmat*   j    ] += Ws1[ik]*X1[k]*a31 + Ws2[ik]*X2[k]*a41;
                        fg[  i+nmat  +nnmat*   j    ] += Vs1[ik]*X1[k]*a31 + Vs2[ik]*X2[k]*a41;
                        fg[(i+nmat*2)+nnmat*   j    ] += Vp1[ik]*X1[k]*a31 + Vp2[ik]*X2[k]*a41;
                        fg[(i+nmat*3)+nnmat*   j    ] += Wp1[ik]*X1[k]*a31 + Wp2[ik]*X2[k]*a41;
                        fg[     i    +nnmat*(nmat+j)] += Ws1[ik]*X1[k]*a32 + Ws2[ik]*X2[k]*a42;
                        fg[  i+nmat  +nnmat*(nmat+j)] += Vs1[ik]*X1[k]*a32 + Vs2[ik]*X2[k]*a42;
                        fg[(i+nmat*2)+nnmat*(nmat+j)] += Vp1[ik]*X1[k]*a32 + Vp2[ik]*X2[k]*a42;
                        fg[(i+nmat*3)+nnmat*(nmat+j)] += Wp1[ik]*X1[k]*a32 + Wp2[ik]*X2[k]*a42;
                        st[     i    +nnmat*   j    ] += Ws1[ik]*X1[k]*a33 + Ws2[ik]*X2[k]*a43;
                        st[  i+nmat  +nnmat*   j    ] += Vs1[ik]*X1[k]*a33 + Vs2[ik]*X2[k]*a43;
                        st[(i+nmat*2)+nnmat*   j    ] += Vp1[ik]*X1[k]*a33 + Vp2[ik]*X2[k]*a43;
                        st[(i+nmat*3)+nnmat*   j    ] += Wp1[ik]*X1[k]*a33 + Wp2[ik]*X2[k]*a43;
                        st[     i    +nnmat*(nmat+j)] += Ws1[ik]*X1[k]*a34 + Ws2[ik]*X2[k]*a44;
                        st[  i+nmat  +nnmat*(nmat+j)] += Vs1[ik]*X1[k]*a34 + Vs2[ik]*X2[k]*a44;
                        st[(i+nmat*2)+nnmat*(nmat+j)] += Vp1[ik]*X1[k]*a34 + Vp2[ik]*X2[k]*a44;
                        st[(i+nmat*3)+nnmat*(nmat+j)] += Wp1[ik]*X1[k]*a34 + Wp2[ik]*X2[k]*a44;
                    }
                }
            }
        } // ...End of loop on layer

        message("Creating final large matrix");
        for (i=0; i<nmat4; ++i) {
            for (j=0; j<nmat2; ++j) {
                BigMat[i+nnmat*j] = -st[i+nnmat*j];
            }
            for (j=nmat2; j<nmat4; ++j) {
                BigMat[i+nnmat*j] = 0;
            }
        }
        for (i=0; i<nmat; ++i) {
            if (backward) {
                BigMat[    i      +nnmat*(i+nmat*2)] = 1;
                BigMat[ (i+nmat)  +nnmat*(i+nmat*2)] = cI*YI[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat*3)] = 1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat*3)] = cI*ZI[i];
            } else {
                BigMat[    i      +nnmat*(i+nmat*2)] = 1;
                BigMat[ (i+nmat)  +nnmat*(i+nmat*2)] = cI*YII[i];
                BigMat[(i+2*nmat) +nnmat*(i+nmat*3)] = 1;
                BigMat[(i+3*nmat) +nnmat*(i+nmat*3)] = cI*ZII[i];
            }
        }

        message("Performing final LU decomposition");
        LUdecompose(BigMat,nnmat,index);

        vector<COMPLEX> Ts(nnmat);
        vector<COMPLEX> Tp(nnmat);

        for (i=0; i<nnmat; ++i) {
            int ins = i+nnmat*n;
            Ts[i]=fg[ins];
            int inp = i+nnmat*(n+nmat);
            Tp[i]=fg[inp];
        }

        message("Solving for fields");
        LUbacksubstitute(BigMat,nnmat,index,Ts);
        LUbacksubstitute(BigMat,nnmat,index,Tp);

        message("Finishing");
        for (i=-order; i<=order; ++i) {
            j= n+i;
            // Conventions for incident polarizations are handled above.
            // Conventions for reflected polarizations need factor of -1 for S
            // and sqrt(-1) for P.  Conjugates are needed for conversion to -iwt from iwt.
            if (backward) {
                r[j].PS() = conj(    Tp[j+2*nmat] );
                r[j].PP() = conj( cI*Tp[j+3*nmat]/nI );
                r[j].SP() = conj(-cI*Ts[j+3*nmat]/nI );
                r[j].SS() = conj(   -Ts[j+2*nmat] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIzi[j])/inckz);
            } else {
                r[j].PS() = conj(   -Tp[j+2*nmat] );
                r[j].PP() = conj( cI*Tp[j+3*nmat]/nII );
                r[j].SP() = conj( cI*Ts[j+3*nmat]/nII );
                r[j].SS() = conj(   -Ts[j+2*nmat] );
                R[j]=MuellerMatrix(r[j])*fabs(real(kIIzi[j])/inckz);
            }
        }
    }

    void RCW_BRDF_Model::setup()
    {
        BRDF_Model::setup();

        //
        // The following only need to be done if the RCW-related parameters change, or
        // if the last incident direction changed...
        //
        if (grating->get_recalc()) {
            grating->SETUP();
            RCW.set_grating(grating);
        }

        if (order<0) error("order < 0");
        if (grating->get_medium_t().index(lambda)!=substrate.index(lambda)) error("grating->medium_t != substrate");

        if (order!=RCW.get_order()) RCW.set_order(order);
        if (type!=RCW.get_type()) RCW.set_type(type);
        if (lambda!=RCW.get_lambda()) RCW.set_lambda(lambda);
        if (thetai/deg!=RCW.get_thetai()) RCW.set_thetai(thetai/deg);
        if (rotation/deg!=RCW.get_rotation()) RCW.set_rotation(rotation/deg);

        V.resize(2*order+1);
        R.resize(2*order+1);

        for (int i=-order; i<=order; ++i) {
            int j=i+order;

            V[j]=RCW.GetDirection(i);
            R[j]=RCW.GetIntensity(i);
        }

        cosalpha = cos(alpha);
        Omega = pi*sqr(alpha);

        low = RCW.GetMinimumPropagatingOrder();
        high = RCW.GetMaximumPropagatingOrder();
    }

    MuellerMatrix RCW_BRDF_Model::mueller()
    {
        SETUP();

        Vector v = polar(1.,thetas,phis);
        MuellerMatrix result = MuellerZero();

        for (int i=low; i<=high; ++i) {
            if (V[i+order]*v>=cosalpha) {
                result += R[i+order]/(Omega*v.z);
            }
        }
        return result;
    }

    DEFINE_MODEL(RCW_Model,Model,"Rigorous coupled wave theory for a grating");
    DEFINE_PARAMETER(RCW_Model,int,order,"Order of calculation","25",0xFF);
    DEFINE_PARAMETER(RCW_Model,int,type,"(0) for Forward/Reflection, (1) for Forward/Transmission, (2) for Backward/Reflection, or (3) for Backward/Transmission","0",0xFF);
    DEFINE_PARAMETER(RCW_Model,double,lambda,"Wavelength (vacuum)","0.532",0xFF);
    DEFINE_PARAMETER(RCW_Model,double,thetai,"Incident angle [deg]","0",0xFF);
    DEFINE_PARAMETER(RCW_Model,double,rotation,"Azimuthal rotation of sample [deg]","0",0xFF);
    DEFINE_PTRPARAMETER(RCW_Model,Grating_Ptr,grating,"Grating","Single_Line_Grating",0xFF);

    DEFINE_MODEL(RCW_BRDF_Model,BRDF_Model,"Rigorous coupled wave theory for a grating, form fitted into a BRDF_Model");
    DEFINE_PARAMETER(RCW_BRDF_Model,double,alpha,"Half angle of diffraction cone [rad]","0.0175",0x01);
    DEFINE_PARAMETER(RCW_BRDF_Model,int,order,"Maximum order","25",0x02);
    DEFINE_PTRPARAMETER(RCW_BRDF_Model,Grating_Ptr,grating,"Grating","Single_Line_Grating",0x02);
}

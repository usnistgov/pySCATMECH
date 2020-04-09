//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: rough.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#include "scatmech.h"
#include "rough.h"

using namespace std;


namespace SCATMECH {

    //
    // Routine to get the isotropic spatial frequency...
    //
    void
    Roughness_BRDF_Model::
    Bragg_Frequency(double& fx,double& fy)
    {
        double _fx,_fy;
        switch (type) {
            case Type_DOWNUP:
            {
                _fx = sin(thetas)*cos(phis)-sin(thetai);
                _fy = sin(thetas)*sin(phis);
            }
            break;
            case Type_DOWNDOWN:
            {
                double m = substrate.n(lambda);
                _fx = sin(thetas)*cos(phis)*m-sin(thetai);
                _fy = sin(thetas)*sin(phis)*m;
            }
            break;
            case Type_UPDOWN:
            {
                double m = substrate.n(lambda);
                _fx = sin(thetas)*cos(phis)*m-sin(thetai)*m;
                _fy = sin(thetas)*sin(phis)*m;
            }
            break;
            case Type_UPUP:
            {
                double m = substrate.n(lambda);
                _fx = sin(thetas)*cos(phis)-sin(thetai)*m;
                _fy = sin(thetas)*sin(phis);
            }
            break;
            default:
                error("Invalid value for BRDF_Model::type = " + to_string(type));
        }

        _fx = _fx/lambda;
        _fy = _fy/lambda;

        double cosrot = cos(rotation);
        double sinrot = sin(rotation);

        fx =  cosrot*_fx + sinrot*_fy;
        fy = -sinrot*_fx + cosrot*_fy;
    }

    double
    ABC_PSD_Function::
    psd(double f)
    {
        SETUP();
        return A/pow(1+sqr(B*f),C/2);
    }

    double
    Table_PSD_Function::
    psd(double f)
    {
        SETUP();
        return T.value(f);
    }

    double
    Fractal_PSD_Function::
    psd(double f)
    {
        SETUP();
        return A/pow(f,exponent);
    }

    void
    Gaussian_PSD_Function::
    setup()
    {
        PSD_Function::setup();
        A = pi*sqr(sigma*length);
        B = pi*length;
    }

    double
    Gaussian_PSD_Function::
    psd(double f)
    {
        SETUP();
        return A*exp(-sqr(B*f));
    }

    /*  THIS SECTION REMOVED SINCE COMPILERS NOW INCLUDE BESSEL FUNCTIONS j0, j1, and jn ...
    static double poly(double x,double* coeff,int n) {
        double accum=0.;
        for (int i=0; i<n; ++i) accum = coeff[i]+x*accum;
        return accum;
    }

    double j1pos(double x)
    {
        static double coeff1[] = {1.,376.9991397,99447.43394,18583304.74,2300535178.0,144725228443.0};
        static double coeff2[] = {-30.16036606,15704.48260,-2972611.439,242396853.1,-7895059235.0,72362614232.0};
        static double coeff3[] = {-0.240337019e-6,0.2457520174e-5,-0.3516396496e-4,0.183105e-2,1.0};
        static double coeff4[] = {0.105787412e-6,-0.88228987e-6,0.8449199096e-5,-0.2002690873e-3,0.04687499995};

        if (fabs(x) < 8.) {
            return x*poly(x*x,coeff2,6)/poly(x*x,coeff1,6);
        } else {
            return sqrt(0.636619772/x)*(cos(x-2.356194491)*poly(sqr(8./x),coeff3,5)-
                                        8./x*sin(x-2.356194491)*poly(sqr(8./x),coeff4,5));
        }
    }

    double j1(double x) {
    	return jn(1,x);
        if (x>0) return j1pos(x);
        else return -j1pos(-x);
    }
    */
    #ifdef _MSC_VER
#define J1(x) _j1(x)
#define Jn(n,x) _jn(n,x)
    #else
#define J1(x) j1(x)
#define Jn(n,x) jn(n,x)
    #endif

    double
    Elliptical_Mesa_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();
        double f = sqrt(sqr(fx*axisx/2.)+sqr(fy*axisy/2.));
        double rx_ry = axisx*axisy/4.;
        //double temp =  (f!=0) ? _j1(2*pi*f)/f : pi;
        double temp = (f!=0) ? J1(2*pi*f)/f : pi;
        return sqr(rx_ry*temp*height)*density;
    }

    static double sinc(double x)
    {
        if (x==0.) return 1;
        else return sin(x)/x;
    }

    double
    Rectangular_Mesa_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();
        return sqr(lengthx*lengthy*sinc(fx*lengthx*pi)*sinc(fy*lengthy*pi)*height)*density;
    }


    double
    Triangular_Mesa_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();

        static const double sqrt3=sqrt(3.);
        static const COMPLEX i(0,1);
        COMPLEX s;
        double eps = 1./side*0.000001;
        if (fabs(fy)<eps) {
            s = (sqrt3 - 3.*i*fx*pi*side - sqrt3*cos(sqrt3*fx*pi*side) +
                 i*sqrt3*sin(sqrt3*fx*pi*side))/(6.*sqr(fx*pi));
            if (fabs(fx) < eps) {
                s = sqrt3*sqr(side)/4.;
            }
        } else if (fabs(sqrt3*fx - fy)<eps) {
            (-i/4.*sqrt3*(exp(i*fy*pi*side)*fy*pi*side - sin(fy*pi*side)))/sqr(fy*pi);
        } else {
            s = (sqrt3*fy*(cos(sqrt3*fx*side*pi)-cos(fy*side*pi)-i*sin(sqrt3*fx*side*pi)) +
                 3.*i*fx*sin(fy*side*pi))/(2.*(-3.*sqr(fx)*fy+cube(fy))*sqr(pi));
        }
        return density * sqr(height) * norm(s);
    }

    double
    Rectangular_Pyramid_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();
        double denom = (cube(fx*pi)*fy*sqr(lengthx) - fx*cube(fy*pi)*sqr(lengthy));
        double s;

        if (fabs(denom*lengthx*lengthy) < 1E-6) {
            if (fabs(fx*lengthx) < 1E-6) {
                if (fabs(fy*lengthy) < 1E-6) {
                    s = (height*lengthx*lengthy)/3.;
                } else {
                    s = (height*lengthx*(-(fy*lengthy*pi*cos(fy*lengthy*pi)) + sin(fy*lengthy*pi)))/(cube(fy*pi)*sqr(lengthy));
                }
            } else if (fabs(fy*lengthy) < 1E-6) {
                s = (height*lengthy*(-(fx*lengthx*pi*cos(fx*lengthx*pi)) + sin(fx*lengthx*pi)))/(cube(fx*pi)*sqr(lengthx));
            } else {
                s = -height*lengthx*(-2.*fy*lengthy*pi+sin(2*fy*lengthy*pi))/(4*cube(fy*pi)*sqr(lengthy));
            }
        } else {
            s = (fy*height*lengthy*cos(fy*lengthy*pi)*sin(fx*lengthx*pi) -
                 fx*height*lengthx*cos(fx*lengthx*pi)*sin(fy*lengthy*pi))/
                denom;
        }
        return density * sqr(s);
    }

    double Triangular_Pyramid_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();
        COMPLEX I(0,1);
        double h = height;
        double L = length;
        double denom = (4.*fx*fy*(3.*pow(fx,4) - 10.*pow(fx,2)*pow(fy,2) + 3.*pow(fy,4))*L*pow(pi,3));

        COMPLEX s = 0;
        if (fabs(denom*sqr(L)*cube(L))<1E-6) {
            if (fabs(fx*L)<1E-6) {
                if (fabs(fy*L)<1E-6) {
                    s = (h*sqr(L))/(4.*sqrt(3.));
                } else {
                    s=(h*(COMPLEX(0,-9) + (COMPLEX(0,1)*(8. + exp(I*sqrt(3.)*fy*L*pi) +
                                                         COMPLEX(0,2)*sqrt(3.)*fy*L*pi))/exp((I*fy*L*pi)/sqrt(3.))))/
                      (4.*cube(fy)*L*cube(pi));
                }
            } else if (fabs(fy*L)<1E-6) {
                s = (sqrt(3.)*h*(fx*L*pi - sin(fx*L*pi)))/(2.*cube(fx)*L*cube(pi));
            } else if (fabs((3*sqr(fx) - sqr(fy))*sqr(L))<1E-6) {
                s = (3.*h*(2.*sqrt(3.)*fy*L*pi - 3.*sin((2.*fy*L*pi)/sqrt(3.))))/(16.*cube(fy)*L*cube(pi));
            } else if (fabs((3*sqr(fx) - sqr(fy))*sqr(L))<1E-6) {
                s = (3.*sqrt(3.)*h*(COMPLEX(0,9) + (COMPLEX(0,-1) - 4.*exp(2.*I*fx*L*pi)*(COMPLEX(0,2) + fx*L*pi))/
                                    exp((4.*I*fx*L*pi)/3.)))/(32.*cube(fx)*L*cube(pi));
            } else {
                s=0.;
            }
        } else {
            s = (3.*h*(COMPLEX(0,1)*fx*(3.*pow(fx,2) - 9.*pow(fy,2) +
                                        exp((2.*I*fy*L*pi)/sqrt(3.))*(-3.*pow(fx,2) + pow(fy,2))) +
                       (COMPLEX(0,8)*fx*pow(fy,2)*cos(fx*L*pi) -
                        2.*sqrt(3.)*fy*(pow(fx,2) + pow(fy,2))*sin(fx*L*pi))/exp((I*fy*L*pi)/sqrt(3.)))) /
                denom;
        }
        return density*norm(s);
    }

    double Parabolic_Dimple_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();
        double fn = sqrt(sqr(fx*axisx/2.)+sqr(fy*axisy/2.));
        double rx_ry = axisx*axisy/4.;

        double temp = (fn!=0) ? Jn(2,2*pi*fn)/sqr(fn*pi) : 0.5;
        return sqr(rx_ry*pi*temp*height)*density;
    }

    double
    Double_PSD_Function::
    psd(double fx,double fy)
    {
        SETUP();
        return psd1->psd(fx,fy) + psd2->psd(fx,fy);
    }

    void Register(const PSD_Function* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(PSD_Function);
            Register_Model(Unit_PSD_Function);
            Register_Model(ABC_PSD_Function);
            Register_Model(Table_PSD_Function);
            Register_Model(Fractal_PSD_Function);
            Register_Model(Gaussian_PSD_Function);
            Register_Model(Elliptical_Mesa_PSD_Function);
            Register_Model(Rectangular_Mesa_PSD_Function);
            Register_Model(Triangular_Mesa_PSD_Function);
            Register_Model(Rectangular_Pyramid_PSD_Function);
            Register_Model(Triangular_Pyramid_PSD_Function);
            Register_Model(Parabolic_Dimple_PSD_Function);
            Register_Model(Double_PSD_Function);
        }
    }

    DEFINE_VIRTUAL_MODEL(Roughness_BRDF_Model,BRDF_Model,
                         "Scattering by a rough surface in the smooth surface limit.");
    DEFINE_PTRPARAMETER(Roughness_BRDF_Model,PSD_Function_Ptr,psd,"Power spectral density function","ABC_PSD_Function",0xFF);


    DEFINE_VIRTUAL_MODEL(PSD_Function,Model,
                         "Virtual class for Power Spectral Density functions");

    DEFINE_MODEL(Unit_PSD_Function,PSD_Function,
                 "Unit power spectrum");

    DEFINE_MODEL(ABC_PSD_Function,PSD_Function,
                 "PSD = A/pow(1+sqr(B*f),C/2)");
    DEFINE_PARAMETER(ABC_PSD_Function,double,A,"Power Spectrum A Parameter [um^4]","0.01",0xFF);
    DEFINE_PARAMETER(ABC_PSD_Function,double,B,"Power Spectrum B Parameter [um]","362",0xFF);
    DEFINE_PARAMETER(ABC_PSD_Function,double,C,"Power Spectrum C Parameter","2.5",0xFF);

    DEFINE_MODEL(Table_PSD_Function,PSD_Function,
                 "PSD given by a table of values.");
    DEFINE_PARAMETER(Table_PSD_Function,Table,T,"Table of PSD versus spatial frequency","1",0xFF);

    DEFINE_MODEL(Fractal_PSD_Function,
                 PSD_Function,
                 "A PSD with a fractal behavior");
    DEFINE_PARAMETER(Fractal_PSD_Function,double,A,"Power spectrum amplitude parameter","0.01",0xFF);
    DEFINE_PARAMETER(Fractal_PSD_Function,double,exponent,"Power spectrum exponent","2.5",0xFF);

    DEFINE_MODEL(Gaussian_PSD_Function,
                 PSD_Function,
                 "A PSD with gaussian statistics");
    DEFINE_PARAMETER(Gaussian_PSD_Function,double,sigma,"Standard deviation of height [um]","0.05",0xFF);
    DEFINE_PARAMETER(Gaussian_PSD_Function,double,length,"Correlation length [um]","5",0xFF);

    DEFINE_MODEL(Elliptical_Mesa_PSD_Function,
                 PSD_Function,
                 "PSD for elliptical mesas or depressions");
    DEFINE_PARAMETER(Elliptical_Mesa_PSD_Function,double,axisx,"Axis of mesas along x-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Elliptical_Mesa_PSD_Function,double,axisy,"Axis of mesas along y-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Elliptical_Mesa_PSD_Function,double,height,"Height or depth  of mesas [um]","0.01",0xFF);
    DEFINE_PARAMETER(Elliptical_Mesa_PSD_Function,double,density,"Surface number density of mesas [um^-2]","0.01",0xFF);

    DEFINE_MODEL(Rectangular_Mesa_PSD_Function,
                 PSD_Function,
                 "PSD for rectangular mesas or depressions");
    DEFINE_PARAMETER(Rectangular_Mesa_PSD_Function,double,lengthx,"Length of mesas along x-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Rectangular_Mesa_PSD_Function,double,lengthy,"Length of mesas along y-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Rectangular_Mesa_PSD_Function,double,height,"Height or depth  of mesas [um]","0.01",0xFF);
    DEFINE_PARAMETER(Rectangular_Mesa_PSD_Function,double,density,"Surface number density of mesas [um^-2]","0.01",0xFF);

    DEFINE_MODEL(Triangular_Mesa_PSD_Function,
                 PSD_Function,
                 "PSD for isosceles mesas or depressions");
    DEFINE_PARAMETER(Triangular_Mesa_PSD_Function,double,side,"Length of side of triangle [um]","1.5",0xFF);
    DEFINE_PARAMETER(Triangular_Mesa_PSD_Function,double,height,"Height or depth of mesas [um]","0.01",0xFF);
    DEFINE_PARAMETER(Triangular_Mesa_PSD_Function,double,density,"Surface number density of mesas [um^-2]","0.01",0xFF);

    DEFINE_MODEL(Rectangular_Pyramid_PSD_Function,
                 PSD_Function,
                 "PSD for rectangular pyramids.");
    DEFINE_PARAMETER(Rectangular_Pyramid_PSD_Function,double,lengthx,"Length of side of pyramid in x-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Rectangular_Pyramid_PSD_Function,double,lengthy,"Length of side of pyramid in y-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Rectangular_Pyramid_PSD_Function,double,height,"Height or depth of pyramid [um]","0.01",0xFF);
    DEFINE_PARAMETER(Rectangular_Pyramid_PSD_Function,double,density,"Surface number density of pyramids [um^-2]","0.01",0xFF);

    DEFINE_MODEL(Triangular_Pyramid_PSD_Function,
                 PSD_Function,
                 "PSD for triangular pyramids.");
    DEFINE_PARAMETER(Triangular_Pyramid_PSD_Function,double,length,"Length of side of base of pyramid [um]","1.5",0xFF);
    DEFINE_PARAMETER(Triangular_Pyramid_PSD_Function,double,height,"Height or depth of pyramid [um]","0.01",0xFF);
    DEFINE_PARAMETER(Triangular_Pyramid_PSD_Function,double,density,"Surface number density of pyramids [um^-2]","0.01",0xFF);

    DEFINE_MODEL(Parabolic_Dimple_PSD_Function,
                 PSD_Function,
                 "PSD for parabolic bumps or dimples");
    DEFINE_PARAMETER(Parabolic_Dimple_PSD_Function,double,axisx,"Axis of feature along x-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Parabolic_Dimple_PSD_Function,double,axisy,"Axis of feature along y-direction [um]","1.5",0xFF);
    DEFINE_PARAMETER(Parabolic_Dimple_PSD_Function,double,height,"Height or depth of features [um]","0.01",0xFF);
    DEFINE_PARAMETER(Parabolic_Dimple_PSD_Function,double,density,"Surface number density of features [um^-2]","0.01",0xFF);

    DEFINE_MODEL(Double_PSD_Function,
                 PSD_Function,
                 "Sum of two PSD functions");
    DEFINE_PTRPARAMETER(Double_PSD_Function,PSD_Function_Ptr,psd1,"First PSD Function","ABC_PSD_Function",0xFF);
    DEFINE_PTRPARAMETER(Double_PSD_Function,PSD_Function_Ptr,psd2,"Second PSD Function","ABC_PSD_Function",0xFF);


} // namespace SCATMECH




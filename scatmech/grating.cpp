//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: grating.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "grating.h"
#include <fstream>
#include <algorithm>
#include <valarray>
#include <list>
#include <limits>
#include "scateval.h"

using namespace std;

//#define HERE    cerr << __FILE__ << " " << __LINE__ << endl;

namespace SCATMECH {

    static inline double mymax(double x,double y) {
        return (x>y)?x:y;
    }

    void Register(const Grating* x)
    {
        static bool regd=false;
        if (!regd) {
            regd=true;

            Register_Model(Grating);

            // Future releases...
            Register_Model(Generic_Grating);
            Register_Model(Single_Line_Grating);
            Register_Model(Sinusoidal_Relief_Grating);
            Register_Model(Sinusoidal_Volume_Grating);
            Register_Model(Triangular_Grating);
            Register_Model(Corner_Rounded_Grating);
            Register_Model(Dielectric_Stack_Grating);
            Register_Model(Overlaid_Grating);
        }
    }

    Fourier_Component_Calculator::
    Fourier_Component_Calculator(double _period,int _order,int _recip)
        : period(_period), order(_order), sum(0.), begin(true), recip(_recip) {}

    void Fourier_Component_Calculator::enter(double x,COMPLEX y)
    {
        if (begin) {
            begin = false;
        } else {
            COMPLEX yy = recip ? 1./y : y;
            double x1 = last_x;
            double x2 = x;

            if (order!=0) {
                sum += COMPLEX(0.,0.5)*
                       (exp(COMPLEX(0.,2*order*pi*x1/period)) -
                        exp(COMPLEX(0.,2*order*pi*x2/period)))*yy/(order*pi);
            } else {
                sum += yy*(x2 - x1)/period;
            }
        }
        last_y = y;
        last_x = x;
    }

    COMPLEX Grating::epsilon(const dielectric_function& e)
    {
        return (COMPLEX)e.epsilon(lambda);
    }

    void Single_Line_Grating::setup()
    {
        Grating::setup();

        COMPLEX eline = epsilon(material);
        COMPLEX espace = epsilon(space);

        position.clear();
        materialx.clear();
        thickness.clear();

        for (int level=0; level<nlevels; ++level) {
            position.push_back(vector<double>());
            materialx.push_back(vector<COMPLEX>());
            materialmux.push_back(vector<COMPLEX>());
            thickness.push_back(height/nlevels);

            double h = (0.5+level)/(double)nlevels;
            double x1 = -topwidth/2.-(bottomwidth-topwidth)*h/2.-offset*h;
            double x2 = topwidth/2.+(bottomwidth-topwidth)*h/2.-offset*h;

            position[level].push_back(-period/2);
            position[level].push_back(x1);
            position[level].push_back(x2);
            position[level].push_back(period/2);
            materialx[level].push_back(espace);
            materialx[level].push_back(espace);
            materialx[level].push_back(eline);
            materialx[level].push_back(espace);
            materialmux[level].push_back(1.);
            materialmux[level].push_back(1.);
            materialmux[level].push_back(1.);
            materialmux[level].push_back(1.);
        }
        materialz = materialy = materialx;
        materialmuz = materialmuy = materialmux;
    }

    void Corner_Rounded_Grating::setup()
    {
        Grating::setup();

        const char no_straight_section[]=
            "No straight section of sidewall";
        const char no_flat_section_top[]=
            "No flat section on top of line";
        const char no_flat_section_bottom[]=
            "No flat section on bottom of trench";

        COMPLEX eps_line=epsilon(material);
        COMPLEX vacuum = epsilon(medium_i);

        vector<double> y(nlevels);
        vector<double> x(nlevels);

        position.resize(nlevels);
        materialx.resize(nlevels);
        materialmux.resize(nlevels);
        thickness.resize(nlevels);

        double theta=sidewall*deg;
        double cottheta = sidewall!=90. ? 1./tan(theta) : 0.;
        double sintheta = sin(theta);
        double costheta = cos(theta);
        double r1=radiust;
        double r2=radiusb;

        if (sidewall<=90) {
            // See TAG Notes 3 April 2006...
            double a2=r1*sintheta;
            double a4=r2*sintheta;
            double a5=r2*(1-costheta)*cottheta;
            double a3=(height-(r1+r2)*(1.-costheta))*cottheta;
            double a1=width/2-a3-a2-a5;
            double a6=a2+a3+a4;

            if (a3<0) error(no_straight_section);
            if (a1<0) error(no_flat_section_top);
            if (a6+a1>period/2.) error(no_flat_section_bottom);

            for (int i=0; i<nlevels; ++i) {
                double xx = (i+0.5)/(double)(nlevels)*a6;
                double yy;
                if (xx<a2) {
                    yy = height-r1+sqrt(r1*r1-xx*xx);
                } else if (xx<a2+a3) {
                    double xxx = a6-a4+a5-xx;
                    yy = xxx*tan(theta);
                } else {
                    double xxx = a6-xx;
                    yy = r2-sqrt(r2*r2-xxx*xxx);
                }
                x[i]=a1+xx;
                y[i]=yy;
            }
        } else {
            // See TAG Notes 25 January 2007...
            double b1 = width/2+r2*((costheta-1)*cottheta-1+sintheta);
            double b2 = width/2+r2*(costheta-1)*cottheta;
            double b3 = width/2+r2*((costheta-1)*cottheta+sintheta);
            double w2 = width-2*height*cottheta;
            double c3 = w2/2-r1*((costheta-1)*cottheta-1+sintheta);
            double c2 = w2/2-r1*(costheta-1)*cottheta;
            double c1 = w2/2-r1*((costheta-1)*cottheta+sintheta);

            if (c2<b2) error(no_straight_section);
            if (c1<0) error(no_flat_section_top);
            if (b3>period/2.) error(no_flat_section_bottom);

            double tspan=c3-c1;
            double bspan=b3-b1;
            double sspan=c3-b1;
            double span=tspan+bspan+sspan;

            for (int i=0; i<nlevels; ++i) {
                double xx = (i+0.5)/(double)(nlevels)*span;
                double yy;
                if (xx<tspan) {
                    double x=xx;
                    yy = height-r1+sqrt(r1*r1-x*x);
                    xx = c1+x;
                } else if (xx<tspan+c3-c2) {
                    double x = xx-tspan;
                    yy = height-r1-sqrt(sqr(r1)-sqr(r1-x));
                    xx = c3-x;
                } else if (xx<tspan+c3-b2) {
                    double x=xx-tspan-c3+c2;
                    double y1 = r2*(1-costheta);
                    double y2 = r1*(1-costheta);
                    yy = height-y2-x/(c2-b2)*(height-y1-y2);
                    xx = c2-x;
                } else if (xx<tspan+c3-b1) {
                    double x=xx-tspan-c3+b2;
                    x = b2-b1-x;
                    yy = r2+sqrt(sqr(r2)-sqr(r2-x));
                    xx = b1+x;
                } else if (xx<tspan+c3-b1+b3-b1) {
                    double x=xx-tspan-c3+b1;
                    yy = r2-sqrt(sqr(r2)-sqr(r2-x));
                    xx = b1+x;
                }
                x[i]=xx;
                y[i]=yy;
            }
        }
        for (int i=0; i<nlevels; ++i) {
            position[i].resize(4);
            position[i][0]=-period/2;
            position[i][1]=-x[i];
            position[i][2]=x[i];
            position[i][3]=period/2;
            materialx[i].resize(4);
            materialx[i][0]=vacuum;
            materialx[i][1]=vacuum;
            materialx[i][2]=eps_line;
            materialx[i][3]=vacuum;
            materialmux[i].resize(4);
            materialmux[i][0]=1.;
            materialmux[i][1]=1.;
            materialmux[i][2]=1.;
            materialmux[i][3]=1.;

            if (i==0) thickness[i] = height-(y[0]+y[1])/2.;
            else if (i==nlevels-1) thickness[i] = (y[nlevels-2]+y[nlevels-1])/2.;
            else thickness[i] = (y[i-1]+y[i])/2-(y[i]+y[i+1])/2.;
        }
        materialz = materialy = materialx;
        materialmuz = materialmuy = materialmux;
    }

    namespace {
        bool tdiff(const dielectric_stack& a, const dielectric_stack& b) {
            if (a.get_n()!=b.get_n()) return true;
            for (int i=0; i<a.get_n(); ++i) {
                if (a.get_t()[i]!=b.get_t()[i]) return true;
            }
            return false;
        }
        bool tdiff(StackModel& a, StackModel& b) {
			return tdiff(a.get_stack(),b.get_stack());
		}
	}

    void Dielectric_Stack_Grating::setup()
    {
        Grating::setup();

        position.clear();
        materialx.clear();
        materialy.clear();
        materialz.clear();
        materialmux.clear();
        materialmuy.clear();
        materialmuz.clear();
        thickness.clear();

        int n = stackepsx->get_n();

        if (stackepsx->get_n()==0 && stackepsy->get_n()==0 && stackepsz->get_n()==0 &&
                stackmux->get_n()==0 && stackmuy->get_n()==0 && stackmuz->get_n()==0) { // This is an empty stack
            return;
        }

        if (stackepsy->get_n()==0 && stackepsz->get_n()==0 &&
                stackmux->get_n()==0 && stackmuy->get_n()==0 && stackmuz->get_n()==0) { // This is an isotropic stack
            materialx.resize(stackepsx->get_n());
            thickness.resize(stackepsx->get_n());
            position.resize(stackepsx->get_n());
            for (int i=0; i<stackepsx->get_n(); ++i) {
                int ii = n-i-1;
                thickness[i]=stackepsx->get_t()[ii];
                position[i].resize(2);
                materialx[i].resize(2);
                position[i][0] = 0;
                position[i][1] = period;
                materialx[i][0] = epsilon(stackepsx->get_e()[ii]);
                materialx[i][1] = materialx[i][0];
            }
            materialz=materialy=materialx;
            return;
        }
        if (stackepsy->get_n()==0 && stackepsz->get_n()==0 &&
                stackmuy->get_n()==0 && stackmuz->get_n()==0) { // This is a isotropic magnetic material
            if (tdiff(*stackepsx,*stackmux)) error("Layers inconsistent");
            magnetic = true;
            // TODO: code for isotropic magnetic material
            materialx.resize(stackepsx->get_n());
            materialmux.resize(stackepsx->get_n());
            thickness.resize(stackepsx->get_n());
            position.resize(stackepsx->get_n());
            for (int i=0; i<stackepsx->get_n(); ++i) {
                int ii = n-i-1;
                thickness[i]=stackepsx->get_t()[ii];
                position[i].resize(2);
                materialx[i].resize(2);
                materialmux[i].resize(2);
                position[i][0] = 0;
                position[i][1] = period;
                materialx[i][0] = epsilon(stackepsx->get_e()[ii]);
                materialx[i][1] = materialx[i][0];
                materialmux[i][0] = epsilon(stackmux->get_e()[ii]);
                materialmux[i][1] = materialmux[i][0];
            }
            materialz=materialy=materialx;
            materialmuz=materialmuy=materialmux;
            return;
        }

        if (stackmux->get_n()==0 && stackmuy->get_n()==0 && stackmuz->get_n()==0) { // This is an anisotropic, but nonmagnetic, stack
            if (tdiff(*stackepsx,*stackepsy) || tdiff(*stackepsx,*stackepsz)) error("Layers inconsistent");
            anisotropic = true;
            materialx.resize(stackepsx->get_n());
            materialy.resize(stackepsx->get_n());
            materialz.resize(stackepsx->get_n());
            thickness.resize(stackepsx->get_n());
            position.resize(stackepsx->get_n());
            for (int i=0; i<stackepsx->get_n(); ++i) {
                int ii = n-i-1;
                thickness[i]=stackepsx->get_t()[ii];
                position[i].resize(2);
                materialx[i].resize(2);
                materialy[i].resize(2);
                materialz[i].resize(2);
                position[i][0] = 0;
                position[i][1] = period;
                materialx[i][0] = epsilon(stackepsx->get_e()[ii]);
                materialx[i][1] = materialx[i][0];
                materialy[i][0] = epsilon(stackepsy->get_e()[ii]);
                materialy[i][1] = materialy[i][0];
                materialz[i][0] = epsilon(stackepsz->get_e()[ii]);
                materialz[i][1] = materialz[i][0];
            }
            return;
        }

        if (tdiff(*stackepsx,*stackepsy) || tdiff(*stackepsx,*stackepsz) || tdiff(*stackepsx,*stackmux) || tdiff(*stackepsx,*stackmuy) || tdiff(*stackepsx,*stackmuz)) error("Layers inconsistent");

        magnetic = true;
        anisotropic = true;
        materialx.resize(stackepsx->get_n());
        materialy.resize(stackepsx->get_n());
        materialz.resize(stackepsx->get_n());
        materialmux.resize(stackepsx->get_n());
        materialmuy.resize(stackepsx->get_n());
        materialmuz.resize(stackepsx->get_n());
        thickness.resize(stackepsx->get_n());
        position.resize(stackepsx->get_n());
        for (int i=0; i<stackepsx->get_n(); ++i) {
            int ii = n-i-1;
            thickness[i]=stackepsx->get_t()[ii];
            position[i].resize(2);
            materialx[i].resize(2);
            materialy[i].resize(2);
            materialz[i].resize(2);
            materialmux[i].resize(2);
            materialmuy[i].resize(2);
            materialmuz[i].resize(2);
            position[i][0] = 0;
            position[i][1] = period;
            materialx[i][0] = epsilon(stackepsx->get_e()[ii]);
            materialx[i][1] = materialx[i][0];
            materialy[i][0] = epsilon(stackepsy->get_e()[ii]);
            materialy[i][1] = materialy[i][0];
            materialz[i][0] = epsilon(stackepsz->get_e()[ii]);
            materialz[i][1] = materialz[i][0];
            materialmux[i][0] = epsilon(stackmux->get_e()[ii]);
            materialmux[i][1] = materialmux[i][0];
            materialmuy[i][0] = epsilon(stackmuy->get_e()[ii]);
            materialmuy[i][1] = materialmuy[i][0];
            materialmuz[i][0] = epsilon(stackmuz->get_e()[ii]);
            materialmuz[i][1] = materialmuz[i][0];
        }
        return;

    }

    void Overlaid_Grating::setup()
    {
        Grating::setup();

        bottom->set_lambda(lambda);
        top->set_lambda(lambda);

        if (period!=top->get_period()) error("period != top.period");
        if (period!=bottom->get_period()) error("period != bottom.period");
        if (medium_i.index(lambda)!=top->get_medium_i().index(lambda)) error("medium_i != top.medium_i");
        if (medium_t.index(lambda)!=bottom->get_medium_t().index(lambda)) error ("medium_t != bottom.medium_t");
        if (top->get_medium_t().index(lambda)!=bottom->get_medium_i().index(lambda)) error ("top.medium_t != bottom.medium_i");
        if (separation<0) error("separation<0");

        thickness = top->get_thickness();
        thickness.push_back(separation);
        thickness.insert(thickness.end(),bottom->get_thickness().begin(),bottom->get_thickness().end());

        if (bottom->is_anisotropic()||top->is_anisotropic()) anisotropic = true;
        if (bottom->is_magnetic() || top->is_magnetic()) magnetic = true;
    }
    COMPLEX Overlaid_Grating::fourierx(int order,int level,int recip) {
        SETUP();
        if (level<top->get_levels()) return top->fourierx(order,level,recip)*exp(COMPLEX(0,2*pi/period*order*overlay));
        if (level==top->get_levels()) {
            if (order==0) {
                if (recip==0) return epsilon(top->get_medium_t());
                else return 1./epsilon(top->get_medium_t());
            } else return 0.;
        }
        return bottom->fourierx(order,level-top->get_levels()-1,recip);
    }

    COMPLEX Overlaid_Grating::fouriery(int order,int level,int recip) {
        SETUP();
        if (level<top->get_levels()) return top->fouriery(order,level,recip)*exp(COMPLEX(0,2*pi/period*order*overlay));
        if (level==top->get_levels()) {
            if (order==0) {
                if (recip==0) return epsilon(top->get_medium_t());
                else return 1./epsilon(top->get_medium_t());
            } else return 0.;
        }
        return bottom->fouriery(order,level-top->get_levels()-1,recip);
    }

    COMPLEX Overlaid_Grating::fourierz(int order,int level,int recip) {
        SETUP();
        if (level<top->get_levels()) return top->fourierz(order,level,recip)*exp(COMPLEX(0,2*pi/period*order*overlay));
        if (level==top->get_levels()) {
            if (order==0) {
                if (recip==0) return epsilon(top->get_medium_t());
                else return 1./epsilon(top->get_medium_t());
            } else return 0.;
        }
        return bottom->fourierz(order,level-top->get_levels()-1,recip);
    }

    COMPLEX Overlaid_Grating::fouriermux(int order,int level,int recip) {
        SETUP();
        if (level<top->get_levels()) return top->fouriermux(order,level,recip)*exp(COMPLEX(0,2*pi/period*order*overlay));
        if (level==top->get_levels()) return (order==0) ? 1. : 0.;
        return bottom->fouriermux(order,level-top->get_levels()-1,recip);
    }

    COMPLEX Overlaid_Grating::fouriermuy(int order,int level,int recip) {
        SETUP();
        if (level<top->get_levels()) return top->fouriermuy(order,level,recip)*exp(COMPLEX(0,2*pi/period*order*overlay));
        if (level==top->get_levels()) return (order==0) ? 1. : 0.;
        return bottom->fouriermuy(order,level-top->get_levels()-1,recip);
    }

    COMPLEX Overlaid_Grating::fouriermuz(int order,int level,int recip) {
        SETUP();
        if (level<top->get_levels()) return top->fouriermuz(order,level,recip)*exp(COMPLEX(0,2*pi/period*order*overlay));
        if (level==top->get_levels()) return (order==0) ? 1. : 0.;
        return bottom->fouriermuz(order,level-top->get_levels()-1,recip);
    }

    void Triangular_Grating::setup()
    {
        Grating::setup();

        COMPLEX vacuum = epsilon(medium_i);
        COMPLEX m=epsilon(material);

        position.resize(nlevels);
        materialx.resize(nlevels);
        materialmux.resize(nlevels);
        thickness.resize(nlevels);

        double apex1=period*aspect;
        double apex2=period*(1.-aspect);

        for (int level=0; level<nlevels; ++level) {
            double x1 = -period/2.+apex1*(nlevels-level-0.5)/nlevels;
            double x2 = period/2.-apex2*(nlevels-level-0.5)/nlevels;

            position[level].resize(4);
            position[level][0]=-period/2;
            position[level][1]=x1;
            position[level][2]=x2;
            position[level][3]=period/2;
            materialx[level].resize(4);
            materialx[level][0]=vacuum;
            materialx[level][1]=vacuum;
            materialx[level][2]=m;
            materialx[level][3]=vacuum;
            materialmux[level].resize(4);
            materialmux[level][0]=1.;
            materialmux[level][1]=1.;
            materialmux[level][2]=1.;
            materialmux[level][3]=1.;
            thickness[level] = amplitude/nlevels;
        }
        materialz = materialy = materialx;
        materialmuz = materialmuy = materialmux;
    }

    void Sinusoidal_Relief_Grating::setup()
    {
        Grating::setup();

        int _nlevels = base!=0 ? nlevels-1 : nlevels;

        COMPLEX vacuum=epsilon(medium_i);
        COMPLEX e=epsilon(material);

        position.resize(nlevels);
        materialx.resize(nlevels);
        materialmux.resize(nlevels);
        thickness.resize(nlevels);

        for (int level=0; level<nlevels; ++level) {
            position[level].resize(0);
            materialx[level].resize(0);
            materialmux[level].resize(0);
            if (level==_nlevels) {
                position[level].push_back(-period/2);
                position[level].push_back(period/2);
                materialx[level].push_back(e);
                materialx[level].push_back(e);
                materialmux[level].push_back(1.);
                materialmux[level].push_back(1.);
                thickness[level] = base;
            } else {
                double xh,tk;
                switch (option) {
                    case 0:
                    {
                        xh = period/2.*(level+0.5)/_nlevels;
                        double x1 = period/2.*(double)level/(double)_nlevels;
                        double x0 = period/2.*(double)(level+1)/(double)_nlevels;
                        tk = amplitude/2.*(cos(2*pi/period*x1)-cos(2*pi/period*x0));
                    }
                    break;
                    case 1:
                    {
                        double z = (level+0.5)/_nlevels;
                        xh = acos(1.-z*2.)*period/pi/2.;
                        tk = amplitude/_nlevels;
                    }
                    break;
                    default:
                        error("Invalid option");
                }
                position[level].push_back(-period/2);
                position[level].push_back(-xh);
                position[level].push_back(xh);
                position[level].push_back(period/2);
                materialx[level].push_back(vacuum);
                materialx[level].push_back(vacuum);
                materialx[level].push_back(e);
                materialx[level].push_back(vacuum);
                materialmux[level].push_back(1.);
                materialmux[level].push_back(1.);
                materialmux[level].push_back(1.);
                materialmux[level].push_back(1.);

                thickness[level] = tk;
            }
        }

        materialz = materialy = materialx;
        materialmuz = materialmuy = materialmux;
    }

    void Sinusoidal_Volume_Grating::setup()
    {
        Grating::setup();

        position.clear();
        materialx.clear();
        thickness.clear();
        for (int i=0; i<nlevels; ++i) {
            thickness.push_back(thick/nlevels);
        }
    }

    COMPLEX Sinusoidal_Volume_Grating::fourierx(int order,int level,int recip)
    {
        SETUP();
        COMPLEX a =((COMPLEX)maximum.epsilon(lambda)+(COMPLEX)minimum.epsilon(lambda))/2.;
        COMPLEX b =((COMPLEX)maximum.epsilon(lambda)-(COMPLEX)minimum.epsilon(lambda))/2.;
        double h = (level+0.5)*thick/nlevels;
        COMPLEX phase = exp(COMPLEX(0,-2*pi*order*h*tan(tilt*deg)/period+pi*order/2));

        if (recip==0) {
            // The following is trivially the Fourier expansion coefficients for a sinusoidal grating...
            if (order==0) return a;
            if (order==1) return phase*COMPLEX(0,0.5)*b;
            if (order==-1) return phase*COMPLEX(0,-0.5)*b;
            else return 0.;
        } else {
            // The solution to the inverse (1/epsilon) is a little trickier.  Assume that
            // epsilon = a + b*sin(2*pi*x/period). We expand 1/epsilon in powers of b, which
            // gives terms proportional to sin^j, each of which can be "easily" Fourier
            // expanded. The following code gives the resulting expansion...
            COMPLEX result = 0.;
            COMPLEX factor = COMPLEX(0,-1)/a;
            for (int j=0; j<20; ++j) {
                bool iszero = true;
                COMPLEX result0 = factor;
                for (int i=-j; i<=j; i+=2) {
                    if (i!=order) {
                        result0 /= double(order-i);
                    } else {
                        result0 *= COMPLEX(0,1);
                        iszero = false;
                    }
                }
                factor *= COMPLEX(0,-1)*b*(double)(j+1)/a;
                if (!iszero) {
                    result += result0;
                }
                if (abs(result0/result) < numeric_limits<double>::epsilon() && j>5) break;
                if (j==20) error("1/epsilon series did not converge");
            }
            return phase*result;
        }
    }

    COMPLEX Sinusoidal_Volume_Grating::eps(double x, int level, int direction)
    {
        COMPLEX a = ((COMPLEX)maximum.epsilon(lambda) + (COMPLEX)minimum.epsilon(lambda)) / 2.;
        COMPLEX b = ((COMPLEX)maximum.epsilon(lambda) - (COMPLEX)minimum.epsilon(lambda)) / 2.;
        double h = (level + 0.5) * thick / nlevels;
        COMPLEX result = a + b * cos(2 * pi * h * tan(tilt * deg) / period + 2 * pi * x / period);
        return result;
    }

    namespace {

        string read_pair(istream& is)
        {
            is >> ws;

            if (is.peek()!='(') {
                string result;
                is >> result;
                return result;
            }


            // Get the contents of the parentheses and place them in the string contents...
            string result;
            result += is.get();
            int level=1;
            while (level>0) {
                if (is.eof()) throw SCATMECH_exception("Mismatched parentheses");
                char next = is.get();
                if (next==')') --level;
                if (next=='(') ++level;

                if (level>=0) result += next;
            }
            return result;
        }

    }

    void Generic_Grating::setup()
    {
        Grating::setup();


        if (filename.size()==0) error("Empty grating description file name");
        string fname = find_file(filename);
        ifstream_with_comments file(fname.c_str());
        if (!file) error("Cannot open grating description file (" + fname + ")");

        int i,j;

        epsx.clear();
        epsy.clear();
        epsz.clear();
        mux.clear();
        muy.clear();
        muz.clear();
        position.clear();
        thickness.clear();
        materialx.clear();
        materialy.clear();
        materialz.clear();
        materialmux.clear();
        materialmuy.clear();
        materialmuz.clear();
        vars.clear();
        boundaries.clear();

        string stemp;

        //*****************************************************************
        //*
        //* First block: read the parameters ...
        //*
        //*****************************************************************

        // Read and check the section heading...
        file >> stemp;
        if (stemp!="PARAMETERS") error("Expected PARAMETERS label");

        // Evaluate the values passed through the parameter pstring...
        Evaluator pe(pstring);

        vars["period"] = period;
        vars["lambda"] = lambda;
        vars["mediumin"] = medium_i.n(lambda);
        vars["mediumik"] = medium_i.k(lambda);
        vars["mediumtn"] = medium_t.n(lambda);
        vars["mediumtk"] = medium_t.k(lambda);

        i=0;
        while (1) {
            file >> stemp;
            if (stemp=="END") break;
            if (file.fail()) error("Error reading file in PARAMETERS section");

            varsmap::iterator a = vars.find(stemp);
            if (a!=vars.end()) error("Duplicate parameter label: " + stemp);

            if (i>pe.NResult())
                error("Number of values in pstring (" + to_string(pe.NResult()) +
                      ") less than number of parameters (>" + to_string(i) + ")");
            double v = pe.Result(i);

            varsmap::iterator p = override.find(stemp);
            if (p!=override.end()) v = p->second;

            vars.insert(a,varsmap::value_type(stemp,v));
            parameters[stemp]=v;

            i++;
        }
        if (i!=pe.NResult())
            error("Number of values in pstring (" + to_string(pe.NResult()) +
                  ") greater than number of parameters (" + to_string(i) + ")");

        //*****************************************************************
        //*
        //* Read in working variables...
        //*
        //*****************************************************************
        file >> stemp;
        if (stemp=="WORKING") {
            while (1) {
                file >> stemp;
                if (stemp=="END") break;

                string expression = file.getquoted();

                if (file.fail()) gg_error("Error reading file in WORKING section:");

                double value = (double)(Evaluator(expression,vars));

                vars[stemp]=value;
            }
            file >> stemp;
        }

        //*****************************************************************
        //*
        //* Second block: read in materials...
        //*
        //*****************************************************************

        if (stemp!="MATERIALS") gg_error("Expected MATERIALS label");

        epsx["medium_i"] = epsilon(medium_i);
        epsy["medium_i"] = epsilon(medium_i);
        epsz["medium_i"] = epsilon(medium_i);
        epsx["medium_t"] = epsilon(medium_t);
        epsy["medium_t"] = epsilon(medium_t);
        epsz["medium_t"] = epsilon(medium_t);

        mux["medium_i"] = 1.;
        muy["medium_i"] = 1.;
        muz["medium_i"] = 1.;
        mux["medium_t"] = 1.;
        muy["medium_t"] = 1.;
        muz["medium_t"] = 1.;

        while (1) {
            // Read material index...
            file >> stemp;
            if (stemp=="END") break;

            if (stemp=="ANISO") {
                anisotropic=true;
                file >> stemp;

                // Read material dielectric functions...
                string smaterialx = read_pair(file);
                string smaterialy = read_pair(file);
                string smaterialz = read_pair(file);

                if (file.fail()) gg_error("Error reading file in MATERIALS section for ANISO material");

                string smaterial2x;
                try {
                    smaterial2x = Evaluator(smaterialx,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2x = smaterialx;
                }
                dielectric_function dfx(smaterial2x);

                string smaterial2y;
                try {
                    smaterial2y = Evaluator(smaterialy,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2y = smaterialy;
                }
                dielectric_function dfy(smaterial2y);

                string smaterial2z;
                try {
                    smaterial2z = Evaluator(smaterialz,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2z = smaterialz;
                }
                dielectric_function dfz(smaterial2z);

                epsmap_t::iterator a = epsx.find(stemp);
                if (a!=epsx.end()) gg_error("Duplicate material index: " + stemp);
                epsx.insert(a,epsmap_t::value_type(stemp,epsilon(dfx)));
                epsy[stemp] = epsilon(dfy);
                epsz[stemp] = epsilon(dfz);
                mux[stemp] = 1.;
                muy[stemp] = 1.;
                muz[stemp] = 1.;
            } else if (stemp=="ANISOMAGNETIC") {
                anisotropic=true;
                magnetic=true;
                file >> stemp;

                // Read material dielectric functions...
                string smaterialx = read_pair(file);
                string smaterialy = read_pair(file);
                string smaterialz = read_pair(file);
                string smaterialmux = read_pair(file);
                string smaterialmuy = read_pair(file);
                string smaterialmuz = read_pair(file);

                if (file.fail()) gg_error("Error reading file in MATERIALS section for ANISOMAGNETIC material");

                string smaterial2x;
                try {
                    smaterial2x = Evaluator(smaterialx,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2x = smaterialx;
                }
                dielectric_function dfx(smaterial2x);

                string smaterial2y;
                try {
                    smaterial2y = Evaluator(smaterialy,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2y = smaterialy;
                }
                dielectric_function dfy(smaterial2y);

                string smaterial2z;
                try {
                    smaterial2z = Evaluator(smaterialz,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2z = smaterialz;
                }
                dielectric_function dfz(smaterial2z);

                string smaterial2mux;
                try {
                    smaterial2mux = Evaluator(smaterialmux,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2mux = smaterialmux;
                }
                dielectric_function dfmux(smaterial2mux);

                string smaterial2muy;
                try {
                    smaterial2muy = Evaluator(smaterialmuy,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2muy = smaterialmuy;
                }
                dielectric_function dfmuy(smaterial2muy);

                string smaterial2muz;
                try {
                    smaterial2muz = Evaluator(smaterialmuz,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2muz = smaterialmuz;
                }
                dielectric_function dfmuz(smaterial2muz);

                epsmap_t::iterator a = epsx.find(stemp);
                if (a!=epsx.end()) gg_error("Duplicate material index: " + stemp);
                epsx.insert(a,epsmap_t::value_type(stemp,epsilon(dfx)));
                epsy[stemp] = epsilon(dfy);
                epsz[stemp] = epsilon(dfz);
                mux[stemp] = epsilon(dfmux);
                muy[stemp] = epsilon(dfmuy);
                muz[stemp] = epsilon(dfmuz);

            } else if (stemp=="MAGNETIC") {
                magnetic=true;

                file >> stemp;

                // Read material dielectric function...
                string smaterial = read_pair(file);
                string smaterialmu = read_pair(file);

                if (file.fail()) gg_error("Error reading file in MATERIALS section for MAGNETIC material");

                string smaterial2;
                try {
                    smaterial2 = Evaluator(smaterial,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2 = smaterial;
                }
                dielectric_function df(smaterial2);

                string smaterial2mu;
                try {
                    smaterial2mu = Evaluator(smaterialmu,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2mu = smaterialmu;
                }
                dielectric_function dfmu(smaterial2mu);

                epsmap_t::iterator a = epsx.find(stemp);
                if (a!=epsx.end()) gg_error("Duplicate material index: " + stemp);
                COMPLEX e = epsilon(df);
                COMPLEX m = epsilon(dfmu);
                epsx.insert(a,epsmap_t::value_type(stemp,e));
                epsy[stemp] = e;
                epsz[stemp] = e;
                mux[stemp] = m;
                muy[stemp] = m;
                muz[stemp] = m;
            } else {
                // Read material dielectric function...
                string smaterial = read_pair(file);

                if (file.fail()) gg_error("Error reading file in MATERIALS section");

                string smaterial2;
                try {
                    smaterial2 = Evaluator(smaterial,vars).ResultString();
                } catch (SCATMECH_exception&) {
                    smaterial2 = smaterial;
                }
                dielectric_function df(smaterial2);

                epsmap_t::iterator a = epsx.find(stemp);
                if (a!=epsx.end()) gg_error("Duplicate material index: " + stemp);
                COMPLEX e = epsilon(df);
                epsx.insert(a,epsmap_t::value_type(stemp,e));
                epsy[stemp] = e;
                epsz[stemp] = e;
                mux[stemp] = 1.;
                muy[stemp] = 1.;
                muz[stemp] = 1.;
            }
        }

        //*****************************************************************
        //*
        //* Third block: read vertices...
        //*
        //*****************************************************************

        file >> stemp;
        if (stemp!="VERTICES") gg_error("Expected MATERIALS label");

        typedef map<string,COMPLEX> verticesmap;
        verticesmap vertices;

        while (1) {
            file >> stemp;
            if (stemp=="END") break;
            COMPLEX vertex = get_complex_value(file);
            if (file.fail()) gg_error("Error reading file in VERTICES section");

            verticesmap::iterator a = vertices.find(stemp);

            if (a!=vertices.end()) gg_error("Duplicate vertex label: " + stemp);

            vertices.insert(a,verticesmap::value_type(stemp,vertex));
        }

        //*****************************************************************
        //*
        //* Fourth block: read boundaries...
        //*
        //*****************************************************************
        file >> stemp;
        if (stemp!="BOUNDARIES") gg_error("Expected BOUNDARIES label");

        message("Boundaries:\n");

        while (!file.eof()) {
            string vertex1,vertex2;
            string imat1,imat2;

            // Read first vertex...
            file >> vertex1;
            if (vertex1=="END") break;
            // If end of file, then end of block...
            if (!file.eof()) {
                if (file.fail()) gg_error("Error reading file at vertex #1");
                verticesmap::iterator v1 = vertices.find(vertex1);
                if (v1==vertices.end()) gg_error("Unknown vertex label: " + vertex1);

                // Read second vertex...
                file >> vertex2;
                if (file.fail()) gg_error("Error reading file at vertex #2");
                verticesmap::iterator v2 = vertices.find(vertex2);
                if (v2==vertices.end()) gg_error("Unknown vertex label: " + vertex2);

                // Read left material index...
                file >> imat1;
                if (file.fail()) gg_error("Error reading file at material #1");
                epsmap_t::iterator m1 = epsx.find(imat1);
                if (m1==epsx.end()) gg_error ("Unknown material label: " + imat1);

                // Read right material index...
                file >> imat2;
                if (file.fail()) gg_error("Error reading file at material #2");
                epsmap_t::iterator m2 = epsx.find(imat2);
                if (m2==epsx.end()) gg_error ("Unknown material: " + imat2);

                // Load into a segment...
                segment seg;
                seg.vertex1=v1->second;
                seg.vertex2=v2->second;
                seg.mat1=imat1;
                seg.mat2=imat2;

                message(imat1 + tab + imat2 + tab + format("%g",real(seg.vertex1)) + tab + format("%g",imag(seg.vertex1)) + '\n' +
                        imat1 + tab + imat2 + tab + format("%g",real(seg.vertex2)) + tab + format("%g",imag(seg.vertex2)) + "\n\n");

                // Store in list of boundaries...
                boundaries.push_back(seg);
            }
        }

        //*****************************************************************
        //*
        //* Get a list of all vertical and horizontal positions...
        //*
        //*****************************************************************
        list<double> X,Y;
        for (i=0; i<(int)boundaries.size(); ++i) {
            X.push_back(real(boundaries[i].vertex1));
            X.push_back(real(boundaries[i].vertex2));
            Y.push_back(imag(boundaries[i].vertex1));
            Y.push_back(imag(boundaries[i].vertex2));
        }

        // Sort them and throw away duplicates
        X.sort();
        X.unique();
        Y.sort();
        Y.unique();

        // Total period of the grating...
        double _period = X.back()-X.front();

        // Total height of grating...
        double height = Y.back()-Y.front();

        double firstx = X.front();
        double prevx = firstx;
        double firsty = Y.front();
        double prevy = firsty;

        list<double>::iterator p;

        vector<double> td;
        double total_td=0;

        //*****************************************************************
        //*
        //* Check for consistency along slices in the vertical direction
        //*
        //*****************************************************************
        list<bounds> Yslice;
        for (p=++(X.begin()); p!=X.end(); ++p) {
            // Get the horizontal position between the two x-positions
            double width = (*p-prevx);
            double x0 = prevx+width/2.;
            prevx = *p;
            //cerr << "Checking at x = " << x0 << endl;

            Yslice.clear();
            for (j=0; j<(int)boundaries.size(); ++j) {
                double x1=real(boundaries[j].vertex1);
                double y1=imag(boundaries[j].vertex1);
                double x2=real(boundaries[j].vertex2);
                double y2=imag(boundaries[j].vertex2);

                // If segment crosses horizontal position...
                if ((x1>=x0 && x2<x0) || (x1<=x0 && x2>x0)) {
                    bounds b;
                    // Get vertical position of crossing...
                    b.x = (x0*y1 - x0*y2 + y2*x1 - y1*x2)/(x1 - x2);

                    // Store positions and material indices...
                    if (x2>x1) {
                        b.mat1 = boundaries[j].mat2;
                        b.mat2 = boundaries[j].mat1;
                    } else {
                        b.mat1 = boundaries[j].mat1;
                        b.mat2 = boundaries[j].mat2;
                    }
                    Yslice.push_back(b);
                }
            }
            Yslice.sort();
            // Store positions and materials in microlayer
            string imat="medium_t";
            for (list<bounds>::iterator q=Yslice.begin(); q!=Yslice.end(); ++q) {
                if (q->mat1!=imat)
                    gg_error("Inconsistent material index discovered along vertical line at x = " + to_string(x0) + ": " + q->mat1 + "!=" + imat);
                imat = q->mat2;
            }
            if (imat!="medium_i")
                gg_error("Inconsistent material index discovered along vertical line at x = " + to_string(x0) + ": " + string("medium_i") + "!=" + imat);
        }

        //*****************************************************************
        //*
        //* Find maximum traversed distance for each sublayer...
        //*
        //* Note: Sublayers refer to regions between different given
        //*       y-values, which are further divided into microlayers.
        //*
        //*****************************************************************
        for (p=++(Y.begin()); p!=Y.end(); ++p) {
            // Get a vertical position in the middle of the sublayer...
            double subheight=(*p-prevy);
            double h=prevy+subheight/2.;

            double max=0;
            for (j=0; j<(int)boundaries.size(); ++j) {
                double x1=real(boundaries[j].vertex1);
                double y1=imag(boundaries[j].vertex1);
                double x2=real(boundaries[j].vertex2);
                double y2=imag(boundaries[j].vertex2);

                // If segment crosses vertical position...
                if ((y1>=h && y2<h) || (y1<=h && y2>h)) {
                    // Calculate the traversed distance...
                    double tdist = fabs((x2-x1)/(y2-y1)*subheight);
                    // Check to see if the materials on both sides are the same...
                    if (boundaries[j].mat1==boundaries[j].mat2) tdist=0;
                    // If this is a bigger distance, save it...
                    if (max<tdist) max=tdist;
                }
            }
            // Store maximum traversed distance...
            td.push_back(max);
            total_td += max;
            prevy=*p;
        }

        //*****************************************************************
        //*
        //* Calculate and store boundaries for each microlayer...
        //*
        //*****************************************************************
        int k=0;
        int kk=0;
        prevy=firsty;
        for (p=++(Y.begin()); p!=Y.end(); ++p) {
            double subheight=(*p-prevy);
            // Scale number of microlayers in sublayer by maximum traversed distance....
            int mlayers = (total_td==0) ? 1 : (int)(td[kk]/total_td*nlayers+0.5);
            // Make sure there is at least one microlayer...
            if (mlayers<1) mlayers = 1;
            for (i=0; i<mlayers; ++i) {
                // Calculate vertical position in middle of microlayer...
                double h=prevy+subheight*(i+0.5)/mlayers;
                // Start new microlayer...
                position.push_back(vector<double>());
                materialx.push_back(vector<COMPLEX>());
                materialy.push_back(vector<COMPLEX>());
                materialz.push_back(vector<COMPLEX>());
                materialmux.push_back(vector<COMPLEX>());
                materialmuy.push_back(vector<COMPLEX>());
                materialmuz.push_back(vector<COMPLEX>());
                thickness.push_back(subheight/mlayers);

                // Get locations of boundaries...
                list<bounds> X;
                for (j=0; j<(int)boundaries.size(); ++j) {
                    double x1=real(boundaries[j].vertex1);
                    double y1=imag(boundaries[j].vertex1);
                    double x2=real(boundaries[j].vertex2);
                    double y2=imag(boundaries[j].vertex2);

                    // If segment crosses vertical position...
                    if ((y1>=h && y2<h) || (y1<=h && y2>h)) {
                        bounds b;
                        // Get horizontal position of crossing...
                        b.x = (h*x1 - h*x2 + x2*y1 - x1*y2)/(y1 - y2);
                        // Store positions and material indices...
                        if (y2>y1) {
                            b.mat1 = boundaries[j].mat1;
                            b.mat2 = boundaries[j].mat2;
                        } else {
                            b.mat1 = boundaries[j].mat2;
                            b.mat2 = boundaries[j].mat1;
                        }
                        X.push_back(b);
                    }
                }
                // If there are no boundaries in the microlayer...
                if (X.size()==0) {
                    // Then it is a homogeneous layer, and we need to make one...
                    string mat0="medium_t";
                    for (list<bounds>::iterator p=Yslice.begin(); p!=Yslice.end(); ++p) {
                        if (p->x<h) {
                            mat0=p->mat2;
                        }
                    }
                    bounds b;
                    b.x=0;
                    b.mat1=mat0;
                    b.mat2=mat0;
                    X.push_back(b);
                }

                X.sort();
                // Store positions and materials in microlayer
                string imat=X.begin()->mat1;
                for (list<bounds>::iterator q=X.begin(); q!=X.end(); ++q) {
                    if (q->mat1!=imat)
                        gg_error("Inconsistent material index discovered along horizontal line at x = " +
                                 to_string(q->x) + ", y = " +
                                 to_string(h) + ": " + q->mat1 + "!=" + imat);
                    imat = q->mat2;
                    position[k].push_back(q->x);
                    materialx[k].push_back(epsx[q->mat1]);
                    materialy[k].push_back(epsy[q->mat1]);
                    materialz[k].push_back(epsz[q->mat1]);
                    materialmux[k].push_back(mux[q->mat1]);
                    materialmuy[k].push_back(muy[q->mat1]);
                    materialmuz[k].push_back(muz[q->mat1]);
                }
                if (X.begin()->mat1!=imat)
                    gg_error("Inconsistent material index discovered along horizontal line at x = " +
                             to_string(X.begin()->x) + ", y = " +
                             to_string(h) + ": " + X.begin()->mat1 + "!=" + imat);
                position[k].push_back(X.begin()->x+_period);
                materialx[k].push_back(epsx[X.begin()->mat1]);
                materialy[k].push_back(epsy[X.begin()->mat1]);
                materialz[k].push_back(epsz[X.begin()->mat1]);
                materialmux[k].push_back(mux[X.begin()->mat1]);
                materialmuy[k].push_back(muy[X.begin()->mat1]);
                materialmuz[k].push_back(muz[X.begin()->mat1]);
                // Increment microlayer number...
                k++;
            }
            prevy=*p;
            // Increment sublayer number...
            ++kk;
        }

        period = _period;

        // Grating wants information stored from the top microlayer...
        reverse(position.begin(),position.end());
        reverse(materialx.begin(),materialx.end());
        reverse(materialy.begin(),materialy.end());
        reverse(materialz.begin(),materialz.end());
        reverse(materialmux.begin(),materialmux.end());
        reverse(materialmuy.begin(),materialmuy.end());
        reverse(materialmuz.begin(),materialmuz.end());
        reverse(thickness.begin(),thickness.end());
    }

    void Generic_Grating::gg_error(const string& message) const
    {
        ofstream file("Generic_Grating_Error.txt");

        print_ggparameters(file);
        file << endl;
        print_variables(file);
        file << endl;
        print_materials(file);
        file << endl;
        print_boundaries(file);

        Grating::error(message + "\n - written error file to Generic_Grating_Error.txt");
    }

    void Generic_Grating::print_ggparameters(ostream& os) const
    {
        os << "Parameters:" << endl;
        for (varsmap::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            os << p->first << tab << p->second << endl;
        }
    }

    void Generic_Grating::print_variables(ostream& os) const
    {
        os << "Variables: " << endl;
        for (varsmap::const_iterator p=vars.begin(); p!=vars.end(); ++p) {
            os << p->first << tab << p->second << endl;
        }
    }
    void Generic_Grating::print_materials(ostream& os) const
    {
        os << "Materials:" << endl;
        os << "Material" << tab << "epsx" << tab << "epsy" << tab << "epsz" << tab << "mux" << tab << "muy" << tab << "muz" << endl;
        for (epsmap_t::const_iterator px=epsx.begin(),py=epsy.begin(),pz=epsz.begin(),qx=mux.begin(),qy=muy.begin(),qz=muz.begin();
                px!=epsx.end(); ++px,++py,++pz,++qx,++qy,++qz) {
            os << px->first << tab << px->second << tab << py->second << tab << pz->second << tab
               << qx->second << tab << qy->second << tab << qz->second << endl;
        }
    }

    void Generic_Grating::print_boundaries(ostream& os) const
    {
        os << "Boundaries:" << endl
           << "mat1" << tab << "mat2" << tab << "x1" << tab << "z1" << endl
           << "mat1" << tab << "mat2" << tab << "x2" << tab << "z2" << endl << endl;

        for (segmentvector::const_iterator p=boundaries.begin(); p!=boundaries.end(); ++p) {
            os << p->mat1 << tab << p->mat2 << tab << real(p->vertex1) << tab << imag(p->vertex1) << endl
               << p->mat1 << tab << p->mat2 << tab << real(p->vertex2) << tab << imag(p->vertex2) << endl << endl;
        }
    }

    COMPLEX
    Generic_Grating::
    get_complex_value(istream& is) const
    {
        string a = read_pair(is);

        if (is.eof()) error("Premature end of file while reading vertex");

        string b = Evaluator(a,vars).ResultString();

        istringstream c(b);
        COMPLEX d;
        c >> d;
        return d;
    }

    void Generic_Grating::error(const string& message) const
    {
        Grating::error("(input file = " + filename + "): " + message);
    }

    COMPLEX Grating::fourierx(int i,int level,int recip)
    {
        SETUP();

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        Fourier_Component_Calculator z(period,i,recip);
        for (int j=0; j<(int)position[level].size(); ++j) {
            double x=position[level][j];
            z.enter(x,materialx[level][j]);
        }
        COMPLEX result = z.result();
        return result;
    }

    COMPLEX Grating::fouriery(int i,int level,int recip)
    {
        SETUP();

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        Fourier_Component_Calculator z(period,i,recip);
        for (int j=0; j<(int)position[level].size(); ++j) {
            double x=position[level][j];
            z.enter(x,materialy[level][j]);
        }
        COMPLEX result = z.result();
        return result;
    }

    COMPLEX Grating::fourierz(int i,int level,int recip)
    {
        SETUP();

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        Fourier_Component_Calculator z(period,i,recip);
        for (int j=0; j<(int)position[level].size(); ++j) {
            double x=position[level][j];
            z.enter(x,materialz[level][j]);
        }
        COMPLEX result = z.result();
        return result;
    }

    COMPLEX Grating::fouriermux(int i,int level,int recip)
    {
        SETUP();

        if (!magnetic) return i==0 ? 1 : 0;

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        Fourier_Component_Calculator z(period,i,recip);
        for (int j=0; j<(int)position[level].size(); ++j) {
            double x=position[level][j];
            z.enter(x,materialmux[level][j]);
        }
        COMPLEX result = z.result();
        return result;
    }

    COMPLEX Grating::fouriermuy(int i,int level,int recip)
    {
        SETUP();

        if (!magnetic) return i==0 ? 1 : 0;

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        Fourier_Component_Calculator z(period,i,recip);
        for (int j=0; j<(int)position[level].size(); ++j) {
            double x=position[level][j];
            z.enter(x,materialmuy[level][j]);
        }
        COMPLEX result = z.result();
        return result;
    }

    COMPLEX Grating::fouriermuz(int i,int level,int recip)
    {
        SETUP();

        if (!magnetic) return i==0 ? 1 : 0;

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        Fourier_Component_Calculator z(period,i,recip);
        for (int j=0; j<(int)position[level].size(); ++j) {
            double x=position[level][j];
            z.enter(x,materialmuz[level][j]);
        }
        COMPLEX result = z.result();
        return result;
    }

    COMPLEX Grating::eps(double x,int level,int direction)
    {
        SETUP();

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));

        epsvector_t& m = direction==2 ? materialz : (direction==1 ? materialy : materialx);

        double x0=position[level][0];
        x -= period*(int)(x/period);
        while (x<x0) x += period;
        while (x>x0+period) x -= period;

        int j;
        for (j=0; j<(int)position[level].size(); ++j) {
            if (position[level][j]>x) {
                return m[level][j];
            }
        }
        error("Invalid position in eps()");
        return 0.;
    }

    COMPLEX Grating::mu(double x,int level,int direction)
    {
        SETUP();

        if (level<0 || level>=(int)thickness.size()) error("Invalid level: " + to_string(level));
        if (!magnetic) return COMPLEX(1,0);

        epsvector_t& m = direction==2 ? materialmuz : (direction==1 ? materialmuy : materialmux);

        double x0=position[level][0];
        x -= period*(int)(x/period);
        while (x<x0) x += period;
        while (x>x0+period) x -= period;

        int j;
        for (j=0; j<(int)position[level].size(); ++j) {
            if (position[level][j]>x) {
                return m[level][j];
            }
        }
        error("Invalid position in mu()");
        return 0.;
    }

    STRING Generic_Grating::get_parameter_base(const STRING& parameter) const
    {
        if (parameter.substr(0,6)=="param.") {
            const_cast<Generic_Grating*>(this)->SETUP();
            varsmap::const_iterator p = parameters.find(parameter.substr(6));
            if (p!=parameters.end()) return format("%16g",p->second);
        }
        return Grating::get_parameter_base(parameter);
    }

    void Generic_Grating::set_parameter_base(const STRING& parameter, const STRING& value)
    {
        if (parameter.substr(0,6)=="param.") {
            override[parameter.substr(6)] = from_string<double>(value);
            set_recalc(1);
        } else {
            Grating::set_parameter_base(parameter,value);
        }
    }

    void  Generic_Grating::print_parameters(std::ostream& os,const std::string& prefix) const
    {
        Grating::print_parameters(os);
        for (varsmap::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            if (prefix.size()!=0) os << prefix << '.';
            os << "param." << p->first << " = " << p->second << " (double) " << endl;
        }
    }


    DEFINE_VIRTUAL_MODEL(Grating,Model,"Generalized description of a grating");
    DEFINE_PARAMETER(Grating,double,period,"Period of the grating [um]","1",0xFF);
    DEFINE_PARAMETER(Grating,dielectric_function,medium_i,"Incident medium","(1,0)",0xFF);
    DEFINE_PARAMETER(Grating,dielectric_function,medium_t,"Transmitted medium","(4.05,0.05)",0xFF);

    DEFINE_MODEL(Single_Line_Grating,Grating,"Model for a single line with adjustable sidewall angles");
    DEFINE_PARAMETER(Single_Line_Grating,dielectric_function,material,"Line material","(4.05,0.05)",0xFF);
    DEFINE_PARAMETER(Single_Line_Grating,dielectric_function,space,"Space between lines","(1,0)",0xFF);
    DEFINE_PARAMETER(Single_Line_Grating,double,height,"Height of grating","0.2",0xFF);
    DEFINE_PARAMETER(Single_Line_Grating,double,topwidth,"Width of top of line [um]","0.2",0xFF);
    DEFINE_PARAMETER(Single_Line_Grating,double,bottomwidth,"Width of bottom of line [um]","0.2",0xFF);
    DEFINE_PARAMETER(Single_Line_Grating,double,offset,"Shift of bottom relative to top [um]","0",0xFF);
    DEFINE_PARAMETER(Single_Line_Grating,int,nlevels,"Number of levels","10",0xFF);

    DEFINE_MODEL(Corner_Rounded_Grating,Grating,"A single line grating with corner rounding and a sidewall angle.");
    DEFINE_PARAMETER(Corner_Rounded_Grating,dielectric_function,material,"Line optical properties","(4.05,0.05)",0xFF);
    DEFINE_PARAMETER(Corner_Rounded_Grating,double,height,"Height of line [um]","0.1",0xFF);
    DEFINE_PARAMETER(Corner_Rounded_Grating,double,width,"Bottom width of line [um]","0.1",0xFF);
    DEFINE_PARAMETER(Corner_Rounded_Grating,double,sidewall,"Sidewall angle [deg]","88",0xFF);
    DEFINE_PARAMETER(Corner_Rounded_Grating,double,radiusb,"Bottom radius [um]","0.010",0xFF);
    DEFINE_PARAMETER(Corner_Rounded_Grating,double,radiust,"Top radius [um]","0.001",0xFF);
    DEFINE_PARAMETER(Corner_Rounded_Grating,int,nlevels,"Number of levels","10",0xFF);

    DEFINE_MODEL(Sinusoidal_Relief_Grating,Grating,"A grating with a sinusoidal relief");
    DEFINE_PARAMETER(Sinusoidal_Relief_Grating,dielectric_function,material,"Line material","(4.05,0.05)",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Relief_Grating,double,amplitude,"Amplitude of sinusoid [um]","0.4",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Relief_Grating,double,base,"Base of sinusoid [um]","0.",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Relief_Grating,int,option,"(0) Divide grating evenly horizontally, or (1) divide grating evenly vertically","0",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Relief_Grating,int,nlevels,"Number of levels","20",0xFF);

    DEFINE_MODEL(Sinusoidal_Volume_Grating,Grating,"A volume grating with a sinusoidal variation in dielectric function");
    DEFINE_PARAMETER(Sinusoidal_Volume_Grating,dielectric_function,minimum,"Maximum index of refraction","(1.52,0)",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Volume_Grating,dielectric_function,maximum,"Minimum index of refraction","(1.50,0)",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Volume_Grating,double,thick,"Thickness of layer [um]","0.1",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Volume_Grating,double,tilt,"Tilt angle of modulation [deg]","0",0xFF);
    DEFINE_PARAMETER(Sinusoidal_Volume_Grating,int,nlevels,"Number of levels","1",0xFF);

    DEFINE_MODEL(Triangular_Grating,Grating,"A grating having a triangular profile");
    DEFINE_PARAMETER(Triangular_Grating,dielectric_function,material,"Triangle medium","(4.05,0.05)",0xFF);
    DEFINE_PARAMETER(Triangular_Grating,double,amplitude,"Amplitude of grating [um]","0.4",0xFF);
    DEFINE_PARAMETER(Triangular_Grating,double,aspect,"Aspect ratio of triangle (location of apex relative to period)","0.5",0xFF);
    DEFINE_PARAMETER(Triangular_Grating,int,nlevels,"Number of levels","20",0xFF);

    DEFINE_MODEL(Generic_Grating,Grating,"A grating described by a file specification");
    DEFINE_PARAMETER(Generic_Grating,string,filename,"Filename","",0xFF);
    DEFINE_PARAMETER(Generic_Grating,string,pstring,"Parameter string","",0xFF);
    DEFINE_PARAMETER(Generic_Grating,int,nlayers,"Approximate number of levels","20",0xFF);

    DEFINE_MODEL(Dielectric_Stack_Grating,Grating,"A grating with zero layers.");
    DEFINE_PTRPARAMETER(Dielectric_Stack_Grating,StackModel_Ptr,stackepsx,"Stack of films","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Dielectric_Stack_Grating,StackModel_Ptr,stackepsy,"Stack of films","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Dielectric_Stack_Grating,StackModel_Ptr,stackepsz,"Stack of films","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Dielectric_Stack_Grating,StackModel_Ptr,stackmux,"Stack of films","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Dielectric_Stack_Grating,StackModel_Ptr,stackmuy,"Stack of films","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Dielectric_Stack_Grating,StackModel_Ptr,stackmuz,"Stack of films","No_StackModel",0xFF);

    DEFINE_MODEL(Overlaid_Grating,Grating,"One grating on top of another, with a specified offset");
    DEFINE_PTRPARAMETER(Overlaid_Grating,Grating_Ptr,top,"Top grating","Single_Line_Grating",0xFF);
    DEFINE_PTRPARAMETER(Overlaid_Grating,Grating_Ptr,bottom,"Bottom grating","Single_Line_Grating",0xFF);
    DEFINE_PARAMETER(Overlaid_Grating,double,overlay,"Overlay distance [um]","0",0xFF);
    DEFINE_PARAMETER(Overlaid_Grating,double,separation,"Thickness of medium between gratings [um]","0",0xFF);

}


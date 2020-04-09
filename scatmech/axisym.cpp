//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: axisym.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "scatmech.h"
#include "axipart.h"
#include "bobvlieg.h"
#include "askuser.h"
#include <vector>
#include <fstream>

using namespace std;


namespace SCATMECH {

    using namespace BobVlieg_Supp;

    Gauss_Legendre_Integration::
    Gauss_Legendre_Integration(int order,double x1,double x2)
    {
        x.resize(order);
        w.resize(order);
        const double epsilon = 6.0e-16;
        int m=(order+1)/2;
        for (int i=0; i<m; ++i) {
            double z=cos(pi*(i+0.75)/(order+0.5));
            double change;
            double deriv;
            do {
                double p = real(Legendre(order,0,z));
                double p1 = real(Legendre(order-1,0,z));
                deriv=order*(z*p-p1)/(z*z-1.0);
                change = p/deriv;
                z -= change;
            } while (fabs(change) > epsilon);
            x[i] = - (x[order-1-i] = z);
            w[i]= w[order-1-i]= 2.0/((1.0-z*z)*sqr(deriv));
        }

        for (int j=0; j<order; ++j) {
            x[j] = (x1 + x2 + x[j]*(x2-x1))/2.;
            w[j] = w[j]*(x2-x1)/2.;
        }
    }

    void
    Axisymmetric_Shape::
    Insert_Points(double angle1,double angle2)
    {
        int npts = (int)(npoints*fabs(angle2-angle1)/pi);

        shape_vec shape1(npts);

        Gauss_Legendre_Integration gli(npts,angle1,angle2);

        for (int i=0; i<gli.n(); ++i) {
            double theta = gli.abscissa(i);
            shape1[i].wt = gli.weight(i);
            shape1[i].abs =    theta;
            shape1[i].shape = f(theta);
            shape1[i].dshape = df(theta);
        }

        shape.insert(shape.end(),shape1.begin(),shape1.end());
    }

    void
    Axisymmetric_Shape::
    setup()
    {
        Model::setup();
        shape.clear();
    }

    double
    Axisymmetric_Shape::
    Get_Volume()
    {
        SETUP();

        double Volume=0;
        for (unsigned int i=0; i<shape.size(); ++i) {
            double theta = shape[i].abs;
            double r = shape[i].shape;
            double dr = shape[i].dshape;
            double dtheta = shape[i].wt;
            Volume+=pi*sqr(r*sin(theta))*(r*sin(theta)-cos(theta)*dr)*dtheta;
        }
        return Volume;
    }

    void
    Axisymmetric_Shape::
    Renormalize(double newV)
    {
        double V = Get_Volume();
        double factor = pow(newV/V,1./3.);
        for (unsigned int i=0; i<shape.size(); ++i) {
            shape[i].shape=shape[i].shape*factor;
            shape[i].dshape=shape[i].dshape*factor;
        }
    }

    double
    Axisymmetric_Shape::
    Get_Base_Length()
    {
        SETUP();

        double zmax=0;
        for (shape_vec::const_iterator p=shape.begin(); p!=shape.end(); ++p) {
            double z = p->shape*cos(p->abs);
            if (z>zmax) {
                zmax = z;
            }
        }
        return zmax;
    }

    double
    Axisymmetric_Shape::
    Get_Breast()
    {
        SETUP();

        double xmax=0;
        for (shape_vec::const_iterator p=shape.begin(); p!=shape.end(); ++p) {
            double x = p->shape*sin(p->abs);
            if (x>xmax) {
                xmax = x;
            }
        }
        return xmax;
    }

    double
    Axisymmetric_Shape::
    Get_MaxRadius()
    {
        SETUP();

        if (shape.size()==0) return 0;
        double maxr = shape[0].shape;
        for (shape_vec::const_iterator p=shape.begin(); p!=shape.end(); ++p) {
            if (p->shape>maxr) maxr = p->shape;
        }
        return maxr;
    }

    void
    Axisymmetric_Shape::
    Flip()
    {
        for (unsigned int i=0; i<shape.size(); ++i) {
            shape[i].abs=pi-shape[i].abs;
            shape[i].dshape=-shape[i].dshape;
        }
    }

    double
    Axisymmetric_Shape::
    df(double theta) const
    {
        const double dtheta = 0.000001;
        return (f(theta+dtheta)-f(theta-dtheta))/2./dtheta;
    }

    void
    Axisymmetric_Shape::
    Write(const string& filename)
    {
        SETUP();

        ofstream shapefile(filename.c_str());
        shapefile << "theta\tr\tDr\tx\tz\tweight\n";
        int i;
        double _length = 0;//Get_Base_Length();

        for (i=0; i<(int)shape.size(); ++i) {
            double theta = shape[i].abs;
            shapefile << theta/deg << tab
                      << shape[i].shape << tab
                      << shape[i].dshape << tab
                      << shape[i].shape*sin(theta) << tab
                      << -shape[i].shape*cos(theta)+_length << tab
                      << shape[i].wt << endl;
        }

        for (i=(int)shape.size()-1; i>=0; --i) {
            double theta = 2.*pi-shape[i].abs;
            shapefile << theta/deg << tab
                      << shape[i].shape << tab
                      << -shape[i].dshape << tab
                      << shape[i].shape*sin(theta) << tab
                      << -shape[i].shape*cos(theta)+_length << tab
                      << shape[i].wt << endl;
        }
    }

    void
    Ellipsoid_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();

        if (vertical<=0 || horizontal<=0 || fabs(offset)>vertical)
            error("Parameters in out of range");
        Insert_Points(0.,pi);
    }

    double
    Ellipsoid_Axisymmetric_Shape::
    f(double theta) const
    {
        double costh=cos(theta);
        double sinth=sin(theta);
        return (horizontal*(-(horizontal*offset*costh) + vertical*sqrt(sqr(horizontal)*sqr(costh) +
                            (sqr(vertical) - sqr(offset))*sqr(sinth))))/
               (sqr(horizontal)*sqr(costh) + sqr(vertical)*sqr(sinth));
    }

    double
    Ellipsoid_Axisymmetric_Shape::
    df(double theta) const
    {
        double a = vertical;
        double b = horizontal;
        double d = offset;
        return
            (b*sin(theta)*(-(a*(pow(a,4.) - 3*pow(b,4.) + 5*sqr(b)*sqr(d) +
                                sqr(a)*(2*sqr(b) - sqr(d)))*cos(theta)) +
                           a*(sqr(a) - sqr(b))*(sqr(a) - sqr(b) - sqr(d))*cos(3*theta) -
                           2*b*d*(-3*sqr(a) + sqr(b) + (-sqr(a) + sqr(b))*cos(2*theta))*
                           sqrt(sqr(b)*sqr(cos(theta)) +
                                (sqr(a) - sqr(d))*sqr(sin(theta)))))/
            (4.*sqr(sqr(b)*sqr(cos(theta)) + sqr(a)*sqr(sin(theta)))*
             sqrt(sqr(b)*sqr(cos(theta)) + (sqr(a) - sqr(d))*sqr(sin(theta)))
            );
    }

    double
    Table_Axisymmetric_Shape::
    f(double theta) const
    {
        double r = const_cast<Table&>(table).value(theta/deg);
        if (scale>0) r*=scale;
        if (r<=0) error("radius <= 0");
        return r;
    }

    void
    Table_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();

        Insert_Points(0.,pi);
        if (scale<0) Renormalize(4.*pi*cube(scale)/3.);
    }

    void
    Conical_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();

        z0 = height/3.-offset;
        beta = atan(radius/height);
        a0 = radius-corner_base*tan(beta)-corner_base/cos(beta);

        theta0 = atan2(a0,z0);
        theta1 = atan2(a0,(z0-corner_base));
        theta2 = atan2((a0+corner_base*cos(beta)),(z0-corner_base-corner_base*sin(beta)));
        c0 = sqrt(sqr(a0)+sqr(z0-corner_base));

        g = height-corner_apex/sin(beta)-z0;
        theta3 = pi - atan(corner_apex*cos(beta)/(g+corner_apex*sin(beta)));

        if (a0<0||g<0||theta0>theta1||theta1>theta2||theta2>theta3||offset>height/3.||offset<-2.*height/3.)
            error("Invalid parameters");

        Insert_Points(0.,theta0);
        Insert_Points(theta0,theta2);
        Insert_Points(theta2,theta3);
        Insert_Points(theta3,pi);
        if (updown) Flip();
        if (renorm) Renormalize(pi*radius*radius*height/3.);
    }

    double
    Conical_Axisymmetric_Shape::
    f(double theta) const
    {
        if (theta<theta0) {
            return z0/cos(theta);
        } else if (theta<theta2) {
            double alpha = theta - theta1;
            return (2*c0*cos(alpha) + 2*sqrt(-sqr(c0) + sqr(corner_base) + sqr(c0)*sqr(cos(alpha))))/2.;
        } else if (theta<theta3) {
            return sin(beta)*(height-z0)/sin(theta-beta);
        } else {
            double alpha = pi-theta;
            return (2*g*cos(alpha) + 2*sqrt(-sqr(g) + sqr(corner_apex) + sqr(g)*sqr(cos(alpha))))/2.;
        }
    }

    double
    Conical_Axisymmetric_Shape::
    df(double theta) const
    {
        if (theta<theta0) {
            return z0*tan(theta)/cos(theta);
        } else if (theta<theta2) {
            double alpha = theta - theta1;
            return c0*(-1. - (c0*cos(alpha))/sqrt(-sqr(c0) + sqr(corner_base) + sqr(c0)*sqr(cos(alpha))))*sin(alpha);
        } else if (theta<theta3) {
            return -((height - z0)/tan(beta - theta)/sin(beta - theta)*sin(beta));
        } else {
            return g*(1. - (1.*g*cos(theta))/sqrt(-sqr(g) + sqr(corner_apex) + sqr(g)*sqr(cos(theta))))*sin(theta);
        }
    }

    double
    Double_Conical_Axisymmetric_Shape::
    f(double theta) const
    {
        if (theta>pi/2.) theta = pi-theta;

        if (theta<theta0) {
            return (2*c1*cos(theta) + 2*sqrt(-sqr(c1) + sqr(corner_apex) + sqr(c1)*sqr(cos(theta))))/2.;
        } else if (theta<theta1) {
            return (height/2.)*sin(alpha)/sin(pi-alpha-theta);
        } else {
            return (2*c2*cos(pi/2-theta) + 2*sqrt(-sqr(c2) + sqr(corner_waist) + sqr(c2)*sqr(cos(pi/2-theta))))/2.;
        }
    }

    void
    Double_Conical_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();

        alpha = atan(radius/(height/2.));
        beta = pi/2 - alpha;
        c1 = height/2. - corner_apex/sin(alpha);
        theta0 = atan((corner_apex*cos(alpha))/(c1+corner_apex*sin(alpha)));

        c2 = radius - corner_waist/sin(beta);
        theta1 = pi/2 - atan((corner_waist*cos(beta))/(c2+corner_waist*sin(beta)));

        if (c1<0 || c2<0 || theta1<theta0)
            error("Invalid parameters");

        Insert_Points(0.,theta0);
        Insert_Points(theta0,theta1);
        Insert_Points(theta1,pi-theta1);
        Insert_Points(pi-theta1,pi-theta0);
        Insert_Points(pi-theta0,pi);

        if (renorm) Renormalize(pi*radius*radius*height/3.);
    }

    double
    Indented_Ellipsoid_Axisymmetric_Shape::
    f(double theta) const
    {
        if (theta>theta0) {
            double theta1=atan(B/vertical);
            double sqrtantheta = sqr(tan(theta));
            double a2 = sqr(vertical);
            double b2 = sqr(B);
            double c2 = sqr(indent);
            double sint = sin(theta);
            double cost = cos(theta);

            double d= (B*(B*indent*cost + vertical*sqrt(4*b2*sqr(cost) +
                          (4*a2 - c2)*sqr(sint))))/
                      (2.*(b2*sqr(cost) + a2*sqr(sint)));
            return d;
        } else {
            return (vertical-indent/2)/cos(theta);
        }
    }

    void
    Indented_Ellipsoid_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();
        B = (horizontal==0) ? vertical : horizontal;
        theta0 = atan(fabs((2.*B*sqrt(indent)*sqrt(2. - indent/vertical))/(sqrt(vertical)*(-2.*vertical + indent))));

        Insert_Points(0.,theta0);
        Insert_Points(theta0,pi);
        if (updown) Flip();
        if (renorm) Renormalize(4.*pi/3.*vertical*B*B);
    }

    void
    Cylinder_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();

        theta0 = atan2(2.*(radius-corner),length);
        theta1 = atan2((radius-corner),(length/2.-corner));
        theta2 = atan2(radius,(length/2.-corner));
        c = sqrt(sqr(length/2.-corner)+sqr(radius-corner));

        if (corner>radius||corner>(length/2.)||corner<0||radius<0||length<0||theta0>theta1 || theta1>theta2)
            error("Invalid parameters");

        Insert_Points(0.,theta0);
        Insert_Points(theta0,theta2);
        Insert_Points(theta2,pi-theta2);
        Insert_Points(pi-theta2,pi-theta0);
        Insert_Points(pi-theta0,pi);
        if (renorm) Renormalize(pi*radius*radius*length);
    }

    double
    Cylinder_Axisymmetric_Shape::
    f(double theta) const
    {
        if (theta>pi/2.) theta = pi-theta;
        if (theta<theta0) {
            return length/2./cos(theta);
        } else if (theta<theta2) {
            double alpha = theta-theta1;
            return (2*c*cos(alpha) + 2*sqrt(sqr(corner) - sqr(c) + sqr(c)*sqr(cos(alpha))))/2.;
        } else {
            return radius/sin(theta);
        }
    }


    double
    Cylinder_Axisymmetric_Shape::
    df(double theta) const
    {
        bool otherside = false;
        double result;
        if (theta>pi/2.) {
            otherside=true;
            theta = pi-theta;
        }
        if (theta<theta0) {
            result = length/2.*tan(theta)/cos(theta);
        } else if (theta<theta2) {
            double alpha = theta-theta1;
            result = c*(-1. - (c*cos(alpha))/sqrt(sqr(corner) - sqr(c) + sqr(c)*sqr(cos(alpha))))*sin(alpha);
        } else {
            result = -radius/tan(theta)/sin(theta);
        }
        if (otherside) return -result;
        else           return result;
    }

    void
    Chebyshev_Axisymmetric_Shape::
    setup()
    {
        Axisymmetric_Shape::setup();

        a = n*pi/(start-end);

        Insert_Points(0.,start);
        Insert_Points(start,end);
        Insert_Points(end,pi);

        if (renorm) Renormalize(4.*pi/3.*vertical*horizontal*horizontal);
    }

    double
    Chebyshev_Axisymmetric_Shape::
    f(double theta) const
    {
        double costh=cos(theta);
        double sinth=sin(theta);

        double r = vertical*horizontal/sqrt(sqr(horizontal*costh)+sqr(vertical*sinth));
        double c;

        if (theta<=start) {
            c = 1.+amplitude;
        } else if (theta<=end) {
            c = 1.+amplitude*(cos(a*(theta-start)));
        } else {
            c = (n%2==0) ? 1.+amplitude : 1.-amplitude;
        }

        return r*c;
    }

    void Axisymmetric_Shape::set_parameter_base(const STRING& parameter, const STRING& value)
    {
        if (parameter=="WriteShape") {
            Write(value);
        } else {
            Model::set_parameter_base(parameter,value);
        }
    }


    void Register(const Axisymmetric_Shape* x)
    {
        static bool regd=false;
        if (!regd) {
            regd=true;

            Register_Model(Axisymmetric_Shape);
            Register_Model(Ellipsoid_Axisymmetric_Shape);
            Register_Model(Conical_Axisymmetric_Shape);
            Register_Model(Indented_Ellipsoid_Axisymmetric_Shape);
            Register_Model(Double_Conical_Axisymmetric_Shape);
            Register_Model(Table_Axisymmetric_Shape);
            Register_Model(Cylinder_Axisymmetric_Shape);
            Register_Model(Chebyshev_Axisymmetric_Shape);
            //Register_Model(Meniscus_Axisymmetric_Shape);
        }
    }

    DEFINE_VIRTUAL_MODEL(Axisymmetric_Shape,Model,"Axisymmetric shape.");
    DEFINE_PARAMETER(Axisymmetric_Shape,int,npoints,"Number of integration points for T matrix","100",0xFF);

    DEFINE_MODEL(Table_Axisymmetric_Shape,Axisymmetric_Shape,"Shape taken from a table");
    DEFINE_PARAMETER(Table_Axisymmetric_Shape,Table,table,"Table of values","0.05",0xFF);
    DEFINE_PARAMETER(Table_Axisymmetric_Shape,double,scale,"Volume-equivalent sphere radius [um] if negative, scaling factor if positive","1",0xFF);

    DEFINE_MODEL(Ellipsoid_Axisymmetric_Shape,Axisymmetric_Shape,"An ellipsoid with a vertical axis of rotation");
    DEFINE_PARAMETER(Ellipsoid_Axisymmetric_Shape,double,vertical,"Vertical radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(Ellipsoid_Axisymmetric_Shape,double,horizontal,"Horizontal radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(Ellipsoid_Axisymmetric_Shape,double,offset,"Origin offset [um]","0",0xFF);

    DEFINE_MODEL(Conical_Axisymmetric_Shape,Axisymmetric_Shape,"Cone");
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,double,radius,"Base radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,double,height,"Vertical height [um]","0.1",0xFF);
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,double,offset,"Offset of origin from default position [um]","0",0xFF);
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,double,corner_base,"Radius of curvature at base corner [um]","0",0xFF);
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,double,corner_apex,"Radius of curvature at apex [um]","0",0xFF);
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,int,updown,"Direction (0) apex up, (1) apex down","0",0xFF);
    DEFINE_PARAMETER(Conical_Axisymmetric_Shape,int,renorm,"Renormalize volume to that of cone (0=no,1=yes)","0",0xFF);

    DEFINE_MODEL(Indented_Ellipsoid_Axisymmetric_Shape,Axisymmetric_Shape,"Elliptical Bump");
    DEFINE_PARAMETER(Indented_Ellipsoid_Axisymmetric_Shape,double,vertical,"Vertical minor axis [um]","0.05",0xFF);
    DEFINE_PARAMETER(Indented_Ellipsoid_Axisymmetric_Shape,double,horizontal,"Horizontal minor axis (0 indicates sphere) [um]","0",0xFF);
    DEFINE_PARAMETER(Indented_Ellipsoid_Axisymmetric_Shape,double,indent,"Indentation [um]","0.01",0xFF);
    DEFINE_PARAMETER(Indented_Ellipsoid_Axisymmetric_Shape,int,updown,"Direction (0: flat towards surface, 1: flat away from surface)","0",0xFF);
    DEFINE_PARAMETER(Indented_Ellipsoid_Axisymmetric_Shape,int,renorm,"Renormalize volume (0=no,1=yes)","0",0xFF);

    DEFINE_MODEL(Double_Conical_Axisymmetric_Shape,Axisymmetric_Shape,"Saucer-shaped particle");
    DEFINE_PARAMETER(Double_Conical_Axisymmetric_Shape,double,height,"Vertical height [um]","0.1",0xFF);
    DEFINE_PARAMETER(Double_Conical_Axisymmetric_Shape,double,radius,"Horizontal radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(Double_Conical_Axisymmetric_Shape,double,corner_apex,"Radius of curvature at apex [um]","0",0xFF);
    DEFINE_PARAMETER(Double_Conical_Axisymmetric_Shape,double,corner_waist,"Radius of curvature at waist [um]","0",0xFF);
    DEFINE_PARAMETER(Double_Conical_Axisymmetric_Shape,int,renorm,"Renormalize volume (0=no,1=yes)","0",0xFF);

    DEFINE_MODEL(Cylinder_Axisymmetric_Shape,Axisymmetric_Shape,"Cylinder-shaped particle");
    DEFINE_PARAMETER(Cylinder_Axisymmetric_Shape,double,radius,"Radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(Cylinder_Axisymmetric_Shape,double,length,"Length [um]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_Axisymmetric_Shape,double,corner,"Corner radius [um]","0",0xFF);
    DEFINE_PARAMETER(Cylinder_Axisymmetric_Shape,int,renorm,"Renormalize volume to that of cylinder (0=no,1=yes)","0",0xFF);

    DEFINE_MODEL(Chebyshev_Axisymmetric_Shape,Axisymmetric_Shape,"A particle based loosely on Chebyshev polynomial");
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,double,vertical,"Mean vertical axis [um]","0.05",0xFF);
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,double,horizontal,"Mean horizontal axis [um]","0.05",0xFF);
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,double,amplitude,"Fractional amplitude of oscillations","0.01",0xFF);
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,int,n,"Number of half oscillations","8",0xFF);
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,double,start,"Start angle [rad]","0",0xFF);
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,double,end,"End angle [rad]","3.14159265",0xFF);
    DEFINE_PARAMETER(Chebyshev_Axisymmetric_Shape,int,renorm,"Renormalize volume to that of ellipsoid (0=no,1=yes)","0",0xFF);



} // namespace SCATMECH





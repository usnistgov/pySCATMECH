//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: grating.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_GRATING_H
#define SCATMECH_GRATING_H

#include "inherit.h"
#include "dielfunc.h"
#include "filmtran.h"
#include <map>

namespace SCATMECH {

    class RCW_Model;

    class Grating : public Model {
        public:
            Grating() {
                anisotropic = false;
                magnetic = false;
                lambda=1.;
            }

            virtual COMPLEX fourier(int order,int level,int recip=0)
            {
                SETUP();
                if (anisotropic) error("Calling fourier with an anisotropic grating");
                return fourierx(order,level,recip);
            }

            // Return the fourier component for epsilon for specific order of a specific level.
            // level = 0 is the closest level to the incident direction.
            // recip == 0 returns the fourier component for epsilon.
            // recip == 1 returns the fourier component for 1/epsilon.
            virtual COMPLEX fourierx(int order,int level,int recip=0);
            virtual COMPLEX fouriery(int order,int level,int recip=0);
            virtual COMPLEX fourierz(int order,int level,int recip=0);

            // Return the fourier component for mu for specific order of a specific level.
            // level = 0 is the closest level to the incident direction.
            // recip == 0 returns the fourier component for mu.
            // recip == 1 returns the fourier component for 1/mu.
            virtual COMPLEX fouriermux(int order,int level,int recip=0);
            virtual COMPLEX fouriermuy(int order,int level,int recip=0);
            virtual COMPLEX fouriermuz(int order,int level,int recip=0);

            // Return the number of levels...
            int get_levels() {
                SETUP();
                return thickness.size();
            }

            // Return the thickness of a specific level (counting from TOP)...
            double get_thickness(int level) {
                SETUP();
                return thickness[level];
            }

            // Return the dielectric constant at position x at level level.  The
            // electric field is direction = 0 (x), 1 (y), and 2 (z)...
            virtual COMPLEX eps(double x,int level,int direction);

            // Return the magnetic constant at position x at level level.  The
            // magnetic field is direction = 0 (x), 1 (y), and 2 (z)...
            COMPLEX mu(double x,int level,int direction);

            // Return the total thickness of the grating...
            double get_total_thickness() {
                SETUP();
                double result=0;
                for (int i=0; i<(int)thickness.size(); ++i) result += thickness[i];
                return result;
            }

            // Return the level number at position z. The top of the grating is 0.
            // If z is inside the range [-get_total_thickness(),0] the function
            // returns the level number.
            // If the position is above the grating, the function returns -1.
            // If the function is below the grating, the function returns -2.
            int get_level(double z) {
                SETUP();
                if (z>0) return -1;
                for (int i=0; i<(int)thickness.size(); ++i) {
                    z += thickness[i];
                    if (z>0) return i;
                }
                return -2;
            }

            void set_lambda(double _lambda) {
                lambda = _lambda;
                set_recalc(0x01);
            }
            double get_lambda() const {
                return lambda;
            }

            bool is_anisotropic() {
                SETUP();
                return anisotropic;
            }
            bool is_magnetic() {
                SETUP();
                return magnetic;
            }

            friend class RCW_Model;

            const std::vector<std::vector<double> > & get_position() {
                SETUP();
                return position;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialx() {
                SETUP();
                return materialx;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialy() {
                SETUP();
                return materialy;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialz() {
                SETUP();
                return materialz;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialmux() {
                SETUP();
                return materialmux;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialmuy() {
                SETUP();
                return materialmuy;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialmuz() {
                SETUP();
                return materialmuz;
            }

            const std::vector<double> & get_thickness() {
                SETUP();
                return thickness;
            }

        protected:

            bool anisotropic;
            bool magnetic;

            typedef std::vector<std::vector<COMPLEX> > epsvector_t;

            // It is the responsibility of every child class to fill the following...
            std::vector<std::vector<double> > position;
            epsvector_t materialx,materialy,materialz;
            epsvector_t materialmux,materialmuy,materialmuz;
            std::vector<double> thickness;

            // epsilon returns the dielectric function at lambda in form N+iK
            COMPLEX epsilon(const dielectric_function& e);

            double lambda;

        protected:

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,period);
            DECLARE_PARAMETER(dielectric_function,medium_i);
            DECLARE_PARAMETER(dielectric_function,medium_t);
    };

    void Register(const Grating* x);
    typedef Model_Ptr<Grating> Grating_Ptr;

    class Fourier_Component_Calculator {
        public:
            Fourier_Component_Calculator(double _period,int _order,int _recip);
            void enter(double x,COMPLEX y);
            COMPLEX result() {
                return sum;
            }

        private:
            double period,last_x;
            int order,recip;
            COMPLEX sum,last_y;
            bool begin;
    };

    class Generic_Grating : public Grating {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(std::string,filename);
            DECLARE_PARAMETER(std::string,pstring);
            DECLARE_PARAMETER(int,nlayers);

        public:

            struct segment {
                COMPLEX vertex1,vertex2;
                std::string mat1,mat2;
            };
            typedef std::vector<segment> segmentvector;

            struct bounds {
                double x;
                std::string mat1,mat2;

                bool operator>(const bounds& b) const {
                    return x>b.x;
                }
                bool operator<(const bounds& b) const {
                    return x<b.x;
                }
            };

            const segmentvector& Get_Boundaries() {
                SETUP();
                return boundaries;
            }

            void print_parameters(std::ostream& os,const std::string& prefix="") const;
            void print_ggparameters(std::ostream& os) const;
            void print_variables(std::ostream& os) const;
            void print_materials(std::ostream& os) const;
            void print_boundaries(std::ostream& os) const;

        protected:
            void setup();
            void error(const std::string& message) const;

            // The following three override the same named function of Model...
            virtual STRING get_parameter_base(const STRING& parameter) const;
            virtual void set_parameter_base(const STRING& parameter, const STRING& value);

        private:

            typedef std::map<std::string,COMPLEX> epsmap_t;
            epsmap_t epsx,epsy,epsz;
            epsmap_t mux,muy,muz;

            std::string filecontents;
            std::string old_filename;

            typedef std::map<std::string,double> varsmap;
            varsmap vars;

            varsmap parameters;
            varsmap override;

            segmentvector boundaries;

            void gg_error(const std::string& message) const;

            COMPLEX get_complex_value(std::istream& is) const;
    };

    class Dielectric_Stack_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stackepsx);
            DECLARE_PARAMETER(StackModel_Ptr,stackepsy);
            DECLARE_PARAMETER(StackModel_Ptr,stackepsz);
            DECLARE_PARAMETER(StackModel_Ptr,stackmux);
            DECLARE_PARAMETER(StackModel_Ptr,stackmuy);
            DECLARE_PARAMETER(StackModel_Ptr,stackmuz);
    };

    class Single_Line_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(dielectric_function,space);
            DECLARE_PARAMETER(double,topwidth);
            DECLARE_PARAMETER(double,bottomwidth);
            DECLARE_PARAMETER(double,offset);
            DECLARE_PARAMETER(int,nlevels);
            DECLARE_PARAMETER(double,height);
    };

    class Corner_Rounded_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,width);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(double,sidewall);
            DECLARE_PARAMETER(double,radiusb);
            DECLARE_PARAMETER(double,radiust);
            DECLARE_PARAMETER(int,nlevels);

    };

    class Triangular_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(double,amplitude);
            DECLARE_PARAMETER(double,aspect);
            DECLARE_PARAMETER(int,nlevels);
    };

    class Sinusoidal_Relief_Grating : public Grating {
        protected:
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,amplitude);
            DECLARE_PARAMETER(double,base);
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(int,nlevels);
            DECLARE_PARAMETER(int,option);

    };

    class Sinusoidal_Volume_Grating : public Grating {
        protected:
            void setup();
        public:
            virtual COMPLEX fourierx(int order,int level,int recip=0);
            virtual COMPLEX eps(double x, int level, int direction);

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,thick);
            DECLARE_PARAMETER(double,tilt);
            DECLARE_PARAMETER(dielectric_function,minimum);
            DECLARE_PARAMETER(dielectric_function,maximum);
            DECLARE_PARAMETER(int,nlevels);
    };


    class Overlaid_Grating : public Grating {
        public:
            virtual COMPLEX fourierx(int order,int level,int recip=0);
            virtual COMPLEX fouriery(int order,int level,int recip=0);
            virtual COMPLEX fourierz(int order,int level,int recip=0);
            virtual COMPLEX fouriermux(int order,int level,int recip=0);
            virtual COMPLEX fouriermuy(int order,int level,int recip=0);
            virtual COMPLEX fouriermuz(int order,int level,int recip=0);

        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(Grating_Ptr,bottom);
            DECLARE_PARAMETER(Grating_Ptr,top);
            DECLARE_PARAMETER(double,overlay);
            DECLARE_PARAMETER(double,separation);
    };


}
#endif

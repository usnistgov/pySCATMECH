//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: gcross.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_GCROSS_H
#define SCATMECH_GCROSS_H

#include "vector3d.h"
#include "matrix3d.h"
#include "crossgrating.h"
#include "crossgrating2.h"

namespace SCATMECH {

    class Generic_CrossGrating : public Gridded_CrossGrating {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(std::string,filename);
            DECLARE_PARAMETER(std::string,pstring);
            DECLARE_PARAMETER(int,nlayers);

        public:

			struct material {
				material() {}
				material(const COMPLEX& _eps1, const COMPLEX& _eps2, const COMPLEX& _eps3, const COMPLEX& _mu1, const COMPLEX& _mu2, const COMPLEX& _mu3) :
					eps1(_eps1), eps2(_eps2), eps3(_eps3), mu1(_mu1), mu2(_mu2), mu3(_mu3) {}
				bool operator==(const material& m) const {
					if (eps1 != m.eps1) return false;
					if (eps2 != m.eps2) return false;
					if (eps3 != m.eps3) return false;
					if (mu1 != m.mu1) return false;
					if (mu2 != m.mu2) return false;
					if (mu3 != m.mu3) return false;
					return true;
				}
				bool operator!=(const material& m) const {
					if (eps1 != m.eps1) return true;
					if (eps2 != m.eps2) return true;
					if (eps3 != m.eps3) return true;
					if (mu1 != m.mu1) return true;
					if (mu2 != m.mu2) return true;
					if (mu3 != m.mu3) return true;
					return false;
				}

				COMPLEX eps1, eps2, eps3, mu1, mu2, mu3;
			};

            struct boundary {
                int n;
                std::string *text;
                Vector vertex1,vertex2,vertex3;
                material mat1,mat2;
            };
            typedef std::vector<boundary> boundary_vector;

            struct bounds {
                double x;
                material mat1,mat2;
                int n; // The line number of this boundary
                std::string *text;

                bool operator<(const bounds& b) {
                    if (x==b.x) return mat2==b.mat1;
                    else return x<b.x;
                }
                bool operator==(const bounds& b) {
                    if (x!=b.x) return false;
					if (mat1 != b.mat1) return false;
					if (mat2 != b.mat2) return false;
                    return true;
                }

            };

            const boundary_vector& Get_Boundaries() {
                SETUP();
                return boundaries;
            }

            void print_parameters(std::ostream& os,const std::string& prefix) const;

        protected:
            void setup();

            // The following three override the same named function of Model...
            virtual STRING get_parameter_base(const STRING& parameter) const;
            virtual void set_parameter_base(const STRING& parameter, const STRING& value);

        private:
            //std::vector<std::vector<std::string> > material;
            //std::vector<std::vector<double> > position;
            std::vector<double> thickness;

            std::string filecontents;
            std::string old_filename;

            typedef std::map<std::string,material> Materialmap;
            Materialmap matmap;

            typedef std::map<std::string,double> varsmap;
            varsmap vars;

            typedef std::map<std::string,Vector> verticesmap;
            verticesmap vertices;


            material mat_t;
            material mat_i;

            varsmap parameters;
            varsmap override;

            boundary_vector boundaries;

            material find_epsilon(const Vector& s,bool throwerror);
            Vector get_triplet_value(std::istream& is) const;
			COMPLEX stringToEpsilon(const std::string& expression, varsmap& vars);
            void error(const std::string& message) const;
    };

}

#endif

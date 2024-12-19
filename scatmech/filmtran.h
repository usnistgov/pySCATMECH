//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: filmtran.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_FILMTRAN_H
#define SCATMECH_FILMTRAN_H

//
// In all of functions in this file, the angles were changed to be complex.
// (Version 3.02, TAG 7 AUG 2002)
//

#include "mueller.h"
#include "dielfunc.h"

namespace SCATMECH {


    class Thin_Film_Transfer_Matrix;
	class StackModel;

    //
    // dielectric_stack is a class which holds all the relavant information about
    // a dielctric stack.
    //
    class dielectric_stack
    {
        public:
            // Routines to read a file containing the film stack information...
            // The following member function was removed 10 OCT 2002
            // void read_stack_file();
            void read_stack_file(const std::string& filename);

            // Routine to ask user for a stack file...
            // Added 10 OCT 2002
            static dielectric_stack AskUser(const std::string& query, const std::string& deflt="");

            // Routine to print state of film...
            void print(std::ostream& os) const;

            // Constructors...
            dielectric_stack() {
                n=0;
                e.resize(0);
                t.resize(0);
            }
            //dielectric_stack(const dielectric_stack& ds);

            // Routine to remove all the films...
            void wash();

            // Routine to add films (one at a time)...
            void grow(const dielectric_function& epsilon, double thickness);

            const std::string& get_name() const {
                return name;
            }
            friend std::ostream& operator<<(std::ostream& os,const dielectric_stack& ds);

            int get_n() const {
                return n;
            }
            const std::vector<dielectric_function>& get_e() const {
                return e;
            }
            const std::vector<double>& get_t() const {
                return t;
            }

            double get_total_thickness() const {
                double result=0;
                for (int i=0; i<n; ++i) result +=t[i];
                return result;
            }

			// The following reverses the order of the stack...
			void reverse();

            // The reflection and transmission coefficients ...
            // for light propagating from the top layer to the bottom layer...
            // NOTE: theta is the "vacuum" angle: theta=asin(sin(internalangle)*n0);
            COMPLEX rp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX rs12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX tp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX ts12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            // for light propagating from the bottom layer to the top layer
            // NOTE: theta is the "vacuum" angle: theta=asin(sin(internalangle)*n0);
            COMPLEX rp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX rs21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX tp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX ts21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;

            // Jones Matrix versions of reflection coefficients...
            // NOTE: theta is the "vacuum" angle: theta=asin(sin(internalangle)*n0);
            JonesMatrix r12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            JonesMatrix t12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            JonesMatrix r21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            JonesMatrix t21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;

            // The reflection and transmission coefficients ...
            // for light propagating from the top layer to the bottom layer...
            // NOTE: angle is the internal angle of incidence
            COMPLEX rp12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX rs12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX tp12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX ts12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            // for light propagating from the bottom layer to the top layer
            // NOTE: angle is the internal angle of incidence
            COMPLEX rp21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX rs21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX tp21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            COMPLEX ts21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;

            // Jones Matrix versions of reflection coefficients...
            // NOTE: angle is the internal angle of incidence
            JonesMatrix r12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            JonesMatrix t12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            JonesMatrix r21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;
            JonesMatrix t21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const;

        protected:
            // The number of layers...
            int n;

            // An array containing the list of dielectric constants...
            std::vector<dielectric_function> e;

            // An array containing the list of thicknesses ...
            std::vector<double> t;

            // A routine used by the constructors and copiers...
            void init();

            std::string name;

        private:
            // The thin film transfer matrices for the stack...
            Thin_Film_Transfer_Matrix matrix_forward(COMPLEX theta,double lambda) const;
            Thin_Film_Transfer_Matrix matrix_backward(COMPLEX theta,double lambda) const;
    };

    //
    // Thin_Film_Transfer_Matrix is a private subclass which handles the
    // reflection and transmission coefficients for a
    // dielectric stack of films.
    //
    class Thin_Film_Transfer_Matrix {
        public:
            // The constuctor requires the incident angle...
            Thin_Film_Transfer_Matrix();

            // The matrix multiplication operator...
            Thin_Film_Transfer_Matrix operator*(Thin_Film_Transfer_Matrix a);

            // The transfer matrix for a single layer...

            // The reflection and transmission coefficients...
            // n0 and nt are the indices on the incident and tranmitting
            // sides of the film.
            COMPLEX rs(COMPLEX theta,const optical_constant& n0,const optical_constant& nt);
            COMPLEX ts(COMPLEX theta,const optical_constant& n0,const optical_constant& nt);
            COMPLEX rp(COMPLEX theta,const optical_constant& n0,const optical_constant& nt);
            COMPLEX tp(COMPLEX theta,const optical_constant& n0,const optical_constant& nt);
            friend Thin_Film_Transfer_Matrix
            Transfer_Matrix(const optical_constant& n,double thickness,
                            COMPLEX theta, double lambda);
        private:
            // There are two matrices {A,B,C,D} for s and p polarized light.
            COMPLEX As;
            COMPLEX Ap;
            COMPLEX Bs;
            COMPLEX Bp;
            COMPLEX Cs;
            COMPLEX Cp;
            COMPLEX Ds;
            COMPLEX Dp;

    };

    inline
    Thin_Film_Transfer_Matrix::
    Thin_Film_Transfer_Matrix()
    {
        As = 1;
        Bs = 0;
        Cs = 0;
        Ds = 1;
        Ap = 1;
        Bp = 0;
        Cp = 0;
        Dp = 1;
    }

    //
    // grow() adds a new layer above the current layers...
    //
    inline
    void
    dielectric_stack::
    grow(const dielectric_function& epsilon,double thickness)
    {
        e.push_back(epsilon);
        t.push_back(thickness);
        ++n;
    }

    //
    // wash() removes all the films...
    //
    inline
    void
    dielectric_stack::
    wash()
    {
        n=0;
        e.resize(0);
        t.resize(0);
    }

    template <>
    void ModelParameterSet(dielectric_stack& variable,const std::string& subparameter,const std::string& value);
    template <>
    std::string ModelParameterGet(dielectric_stack& variable,const std::string& subparameter);
    template <>
    void ModelParameterAskUser(dielectric_stack& variable,const std::string& prompt);

	class StackModel : public Model, protected dielectric_stack
	{
		DECLARE_MODEL();
	public:

		const dielectric_stack& get_stack() {SETUP(); return *this;}

        // Routine to remove all the films...
		void wash() {SETUP(); dielectric_stack::wash();}

        // Routine to add films (one at a time)...
		void grow(const dielectric_function& epsilon, double thickness) {
			SETUP(); dielectric_stack::grow(epsilon,thickness);
		}

		int get_n() {
			SETUP();
			return dielectric_stack::get_n();
		}
        
		const std::vector<dielectric_function>& get_e() {
			SETUP(); return dielectric_stack::get_e();
        }

		const std::vector<double>& get_t() {
			SETUP();
			return dielectric_stack::get_t();
		}
       
		double get_total_thickness() {
			SETUP();
			return dielectric_stack::get_total_thickness();
		}

		// The following reverses the order of the stack...
		void reverse() {SETUP(); dielectric_stack::reverse();}

        // The reflection and transmission coefficients ...
        // for light propagating from the top layer to the bottom layer...
        // NOTE: theta is the "vacuum" angle: theta=asin(sin(internalangle)*n0);
        COMPLEX rp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rp12(theta,lambda,n0,nt);}
        COMPLEX rs12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rs12(theta,lambda,n0,nt);}
		COMPLEX tp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::tp12(theta,lambda,n0,nt);}
        COMPLEX ts12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::ts12(theta,lambda,n0,nt);}

        // for light propagating from the bottom layer to the top layer
        // NOTE: theta is the "vacuum" angle: theta=asin(sin(internalangle)*n0);
        COMPLEX rp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rp21(theta,lambda,n0,nt);}
		COMPLEX rs21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rs21(theta,lambda,n0,nt);}
        COMPLEX tp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::tp21(theta,lambda,n0,nt);}
        COMPLEX ts21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::ts21(theta,lambda,n0,nt);}

        // Jones Matrix versions of reflection coefficients...
        // NOTE: theta is the "vacuum" angle: theta=asin(sin(internalangle)*n0);
        JonesMatrix r12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::r12(theta,lambda,n0,nt);}
        JonesMatrix t12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::t12(theta,lambda,n0,nt);}
        JonesMatrix r21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::r21(theta,lambda,n0,nt);}
        JonesMatrix t21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::t21(theta,lambda,n0,nt);}

        // The reflection and transmission coefficients ...
        // for light propagating from the top layer to the bottom layer...
        // NOTE: angle is the internal angle of incidence
        COMPLEX rp12i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rp12i(theta,lambda,n0,nt);}
        COMPLEX rs12i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rs12i(theta,lambda,n0,nt);}
        COMPLEX tp12i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::tp12i(theta,lambda,n0,nt);}
        COMPLEX ts12i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::ts12i(theta,lambda,n0,nt);}

		// for light propagating from the bottom layer to the top layer
        // NOTE: angle is the internal angle of incidence
        COMPLEX rp21i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rp21i(theta,lambda,n0,nt);}
        COMPLEX rs21i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::rs21i(theta,lambda,n0,nt);}
        COMPLEX tp21i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::tp21i(theta,lambda,n0,nt);}
        COMPLEX ts21i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::ts21i(theta,lambda,n0,nt);}

        // Jones Matrix versions of reflection coefficients...
        // NOTE: angle is the internal angle of incidence  
		JonesMatrix r12i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::r12i(theta,lambda,n0,nt);}
        JonesMatrix t12i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::t12i(theta,lambda,n0,nt);}
        JonesMatrix r21i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::r21i(theta,lambda,n0,nt);}
		JonesMatrix t21i(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) 
			{SETUP(); return dielectric_stack::t21i(theta,lambda,n0,nt);}
	};

    typedef Model_Ptr<StackModel> StackModel_Ptr;

    void Register(const StackModel* x);

	class No_StackModel : public StackModel
	{
		DECLARE_MODEL();
	protected:
		void setup() {StackModel::setup(); wash();}
	};

	class Stack_StackModel : public StackModel
	{
		DECLARE_MODEL();
		DECLARE_PARAMETER(dielectric_stack,stack);
	protected:
		void setup() {
			StackModel::setup(); 
			dielectric_stack::operator=(stack);
		}
	};

	class SingleFilm_StackModel : public StackModel
	{
		DECLARE_MODEL();
		DECLARE_PARAMETER(dielectric_function,material);
		DECLARE_PARAMETER(double,thickness);
	protected:
		void setup() {
			StackModel::setup();
			wash();
			grow(material,thickness);
		}
	};

	class DoubleFilm_StackModel : public StackModel
	{
		DECLARE_MODEL();
		DECLARE_PARAMETER(dielectric_function,material1);
		DECLARE_PARAMETER(double,thickness1);
		DECLARE_PARAMETER(dielectric_function,material2);
		DECLARE_PARAMETER(double,thickness2);
	protected:
		void setup() {
			StackModel::setup();
			wash();
			grow(material1,thickness1);
			grow(material2,thickness2);
		}
	};

	class GradedFilm_StackModel : public StackModel
	{
		DECLARE_MODEL();
		DECLARE_PARAMETER(dielectric_function,start);
		DECLARE_PARAMETER(dielectric_function,end);
		DECLARE_PARAMETER(double,thickness);
		DECLARE_PARAMETER(int,steps);
	protected:
		void setup() {
			StackModel::setup();
			wash();
			if (steps<2) error("steps<2");
			const Table::VECTOR& startT = start.get_n().get_table();
			const Table::VECTOR& endT = end.get_n().get_table();

			int startn = startT.size();
			int endn = endT.size();
			double start1 = startT[0].first;
			double start2 = startT[startn-1].first;
			double end1 = endT[0].first;
			double end2 = endT[endn-1].first;
			int nlambda = startn>endn ? startn : endn;
			double lambdamin = start1 < end1 ? start1 : end1;
			double lambdamax = start2 > end2 ? start2 : end2;
			for (int i=0;i<steps;++i) {
				double f = (i+0.5)/steps;
				Table::VECTOR nn,kk;
				for (int j=0;j<nlambda;++j) {
					double x = double(j)/(nlambda-1.);
					double f1 = (1-f);
					double f2 = f;
					double lambda = lambdamin*(1-x)+lambdamax*x;
					COMPLEX e1 = start.epsilon(lambda);
					COMPLEX e2 = end.epsilon(lambda);
					COMPLEX eps = e1*f1+e2*f2;
					COMPLEX nk = sqrt(eps);
					nn.push_back(Table::PAIR(lambda,real(nk)));
					kk.push_back(Table::PAIR(lambda,imag(nk)));
				}
				dielectric_function eps;
				eps.set_n(nn);
				eps.set_k(kk);
				grow(eps,thickness/steps);
			}
		}
	};

} // namespace SCATMECH



#endif

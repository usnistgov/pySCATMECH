//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: filmtran.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

//*****************************************************************************
//*
//* Operations needed to calculate reflectance from and transmittance through
//* a film stack.
//*
//* See G. R. Fowles, _Introduction_to_Modern_Optics_, Second Edition,
//*     (Holt, Rinehart and Winston, New York, 1975)
//*     Section 4.4, Theory of Multilayer Films
//*     Note that some of the work is left to the reader in Problem 4.10.
//*
//*****************************************************************************

//
// In all of functions in this file, the angles were changed to be complex.
// (Version 3.02, TAG 7 AUG 2002)
//

//*****************************************************************************
// Major changes were made in SCATMECH version 3 to the classes
// dielectric_stack and Thin_Film_Transfer_Matrix so that the wavelength
// could more easily be changed.   (TAG 19 JAN 2001)
//*****************************************************************************
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "scatmech.h"
#include "filmtran.h"
#include "fresnel.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {



    // Largest number of layers envisioned...
    // static const int MAXSTACK = 64;

    // The square root of -1...
    static const COMPLEX I = COMPLEX(0,1);


    //
    // The following gives the transfer matrices for a
    // layer of given thickness (normalized by 2*pi/lambda),
    // given index, and given external angle of incidence.
    //
    Thin_Film_Transfer_Matrix
    Transfer_Matrix(const optical_constant& n,double thickness,COMPLEX theta,double lambda)
    {
        Thin_Film_Transfer_Matrix tftm;
        COMPLEX cos_thetai=cos_internal_angle(n,theta);
        COMPLEX klcostheta = (COMPLEX)n*2.*pi*thickness/lambda*cos_thetai;
        COMPLEX p = cos_thetai*(COMPLEX)n;
        COMPLEX q = cos_thetai/(COMPLEX)n;
        COMPLEX sinklcostheta = sin(klcostheta);
        tftm.As = tftm.Ap = cos(klcostheta);
        tftm.Bs = -I/p*sinklcostheta;
        tftm.Cs = -I*p*sinklcostheta;
        tftm.Bp = -I/q*sinklcostheta;
        tftm.Cp = -I*q*sinklcostheta;
        tftm.Ds = tftm.Dp = tftm.As;
        return tftm;
    }

    //
    // The following gives the reflectance and transmittance coefficients for
    // s- or p-polarized light with indices n0 and nt surrounding the layers.
    // The light is assumed to be incident on the n0 side.  That is nt is the
    // index of the substrate.
    //
    COMPLEX
    Thin_Film_Transfer_Matrix::
    rs(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX p0 = (COMPLEX)n0*cos_internal_angle(n0,theta);
        COMPLEX pt = (COMPLEX)nt*cos_internal_angle(nt,theta);
        return ((As+Bs*pt)*p0-(Cs+Ds*pt))/((As+Bs*pt)*p0+(Cs+Ds*pt));
    }

    COMPLEX
    Thin_Film_Transfer_Matrix::
    rp(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX q0 = cos_internal_angle(n0,theta)/(COMPLEX)n0;
        COMPLEX qt = cos_internal_angle(nt,theta)/(COMPLEX)nt;
        return ((Ap+Bp*qt)*q0-(Cp+Dp*qt))/((Ap+Bp*qt)*q0+(Cp+Dp*qt));
    }

    COMPLEX
    Thin_Film_Transfer_Matrix::
    ts(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX p0 = (COMPLEX)n0*cos_internal_angle(n0,theta);
        COMPLEX pt = (COMPLEX)nt*cos_internal_angle(nt,theta);
        return (2.*p0)/((As+Bs*pt)*p0+(Cs+Ds*pt));
    }

    COMPLEX
    Thin_Film_Transfer_Matrix::
    tp(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX p0 = cos_internal_angle(n0,theta)/(COMPLEX)n0;
        COMPLEX pt = cos_internal_angle(nt,theta)/(COMPLEX)nt;
        return (2.*p0)/((Ap+Bp*pt)*p0+(Cp+Dp*pt))*(COMPLEX)n0/(COMPLEX)nt;
    }

    //
    // Multiple layers can be accounted for my multiplication of the transfer
    // matrices for each layer...
    //
    Thin_Film_Transfer_Matrix
    Thin_Film_Transfer_Matrix::
    operator*(Thin_Film_Transfer_Matrix a)
    {
        Thin_Film_Transfer_Matrix t;
        t.As = As*a.As+Bs*a.Cs;
        t.Ap = Ap*a.Ap+Bp*a.Cp;
        t.Bs = As*a.Bs+Bs*a.Ds;
        t.Bp = Ap*a.Bp+Bp*a.Dp;
        t.Cs = Cs*a.As+Ds*a.Cs;
        t.Cp = Cp*a.Ap+Dp*a.Cp;
        t.Ds = Cs*a.Bs+Ds*a.Ds;
        t.Dp = Cp*a.Bp+Dp*a.Dp;
        return t;
    }

    //
    // The dielectric_stack class holds information about a stack of films...
    //

    //
    // Copy constructor...
    //
    /*dielectric_stack::
    dielectric_stack(const dielectric_stack& ds)
    {
        n=ds.n;
        e.resize(0);
        t.resize(0);
        int i;
        for (i=0;i<n;++i) {
            e.push_back(ds.e[i]);
            t.push_back(ds.t[i]);
        }
    }*/


    //
    // The net transfer matrix for light from the "air" side of the layer.
    //Thin_Film_Transfer_Matrix
    Thin_Film_Transfer_Matrix
    dielectric_stack::
    matrix_forward(COMPLEX theta,double lambda) const
    {
        Thin_Film_Transfer_Matrix M;
        int i;
        for (i=n-1; i>=0; --i)
            M = M*Transfer_Matrix(e[i].index(lambda),t[i],theta,lambda);
        return M;
    }

    //
    // The net transfer matrix for light from the "substrate" side of the layer.
    //
    Thin_Film_Transfer_Matrix
    dielectric_stack::
    matrix_backward(COMPLEX theta,double lambda) const
    {
        Thin_Film_Transfer_Matrix M;
        int i;
        for (i=0; i<n; ++i)
            M = M*Transfer_Matrix(e[i].index(lambda),t[i],theta,lambda);
        return M;
    }

	void dielectric_stack::reverse()
	{
		for (int i=0;i<n/2;++i) {
			swap(e[i],e[n-i-1]);
			swap(t[i],t[n-i-1]);
		}
	}

	//
    // The following give the reflectances of the dielectric stack...
    //   (r,t) reflectance, transmittance
    //   (s,p) polarization
    //   (12,21) forward vs. backward reflection or transmittance
    //   The first argument is always the index of where the light is from
    //   and the second is always the index of where the light is going.
    //
    COMPLEX
    dielectric_stack::
    rs12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rs21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    JonesMatrix
    dielectric_stack::
    r12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(rp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(tp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    r21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(rp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(tp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    static
    COMPLEX
    ArcSin(COMPLEX a)
    {
        if (imag(a)==0.) a += COMPLEX(0,-1E-10); // Needed to ensure the right branch cut...
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }

    COMPLEX
    dielectric_stack::
    rs12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rs21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    JonesMatrix
    dielectric_stack::
    r12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(rp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(tp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    r21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(rp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(tp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    //
    // The following two routines read a file containing data for a thin film
    // File consists of COMPLEX indices followed by thicknesses in order of growth.
    // For example:
    //      (1.8,0.) 0.5
    //      (1.3,0.) 0.6
    //      silicon.dat 0.3
    //      (1.4,0.) 1.0
    //      ...
    // A file can be blank, indicating that no films have been grown.
    //

    dielectric_stack
    dielectric_stack::
    AskUser(const string& prompt,const string& deflt)
    {
        // Define and initialize the dielectric stack...
        dielectric_stack stack;

        while (1) {
            stack.wash();

            try {
                // Create and display the full prompt...
                string prompt2 = prompt
                                 + "\n(':filename' or index and thickness) <"
                                 + deflt + "> :";
                SCATMECH_output << prompt2;

                // Get response from user...
                string response = SCATMECH_input.getstrline();


                if (response.size()==0) {
                    // If response is blank, return the default...
                    stack.read_stack_file(deflt);
                } else if (response[0]==':') {
                    // If response begins with a colon, use the next word as a filename...
                    istringstream responsestr(response);
                    responsestr.get();
                    string filename;
                    responsestr >> filename;
                    stack.read_stack_file(filename);
                } else {
                    int i=1;
                    bool finished=false;

                    string dielfunc;
                    double thickness=-1;

                    while (!finished) {

                        // Extract the dielectric function and thickness of the layer...
                        istringstream sresponse(response);
                        dielfunc.clear();
                        thickness=-1;
                        sresponse >> dielfunc >> thickness;

                        // If no dielectric function, assume finished...
                        if (dielfunc.size()==0) finished=true;
                        // Otherwise...
                        else {
                            // If thickness is less than zero, or no thickness given...
                            if (thickness<0.) throw SCATMECH_exception("Thickness must be greater than zero.");

                            // If both the dielectric function and the thickness is good, then grow film...
                            stack.grow(dielectric_function(dielfunc),thickness);

                            // Prompt user for next layer...
                            ostringstream message;
                            message << "Layer #" << ++i << ":";
                            SCATMECH_output << message.str() << endl;

                            // Get next response...
                            response = getstrline(sresponse);
                            string next;
                            istringstream temp(response);
                            temp >> next;
                            if (next.size()==0) {
                                response = SCATMECH_input.getstrline();
                            }
                        }
                    }
                }
                ostringstream message;
                stack.print(message);
                SCATMECH_output << message.str();
                return stack;
            }
            catch(const exception& e) {
                SCATMECH_output << e.what() << endl;
                //instack_unwind();
            }
        }
    }

    void
    dielectric_stack::
    read_stack_file(const string& film_file_name)
    {
        // If string pointer is NULL.. then assume no film...
        if (film_file_name.size()==0) {
            wash();
            return;
        }

        // Open file stream...
        string fname = find_file(film_file_name);
        ifstream_with_comments film_file(fname.c_str());
        // and complain if it doesn't exist...

        if (!film_file) {
            throw SCATMECH_exception("Cannot open file: " + film_file_name);
        }

        do {
            // Read values...
            string material;
            film_file >> material;

            // A blank line exits the routine...
            if (material.size()==0) return;

            // Get the film thickness
            double thickness=-1.;
            film_file >> thickness;
            if (thickness<0.) throw SCATMECH_exception("Thickness must be greater than zero.");

            // Deposit that film...
            grow(dielectric_function(material),thickness);

        } while (!film_file.eof());
    }

    //
    // The following prints out the status of the film...
    //
    void
    dielectric_stack::
    print(ostream& os) const
    {
        os << endl << "#  index                 thickness " << endl
           << "-----------------------------------" << endl;
        if (n==0) {
            os << "No layers" << endl;
        } else {
            for (int i=0; i<n; ++i) {
                string newname = e[i].get_name().substr(0,20);
                while (newname.size()<20) newname += ' ';
                os << i << ": "
                   << newname
                   << "  " << t[i] << endl;
            }
        }
    }

    ostream&
    operator<<(std::ostream& os,const dielectric_stack& ds)
    {
        if (ds.n==0) {
            os << "(empty)";
        } else {
            for (int i=0; i<ds.n; ++i) {
                string name = ds.e[i].get_name();
                os << name << ' ' << ds.t[i];
                if (i!=ds.n-1) os << ' ';
            }
        }
        return os;
    }

    template <>
    void
    ModelParameterSet(dielectric_stack& variable,const string& parameter,const string& value)
    {
        if (parameter.empty()) {
            variable.wash();
			istringstream iss(value);
            string material;
            iss >> material;
            if (iss.fail()) return;
            double thickness;
            iss >> thickness;
            if (iss.fail()) {
                variable.read_stack_file(material);
                return;
            }
            variable.grow(material,thickness);
            while (1) {
                iss >> material;
                if (iss.fail()) return;
                iss >> thickness;
                if (iss.fail()) throw SCATMECH_exception("Invalid thickness parameter for dielectric_stack");
                variable.grow(material,thickness);
            }
        } else if (parameter == "grow") {
            istringstream iss(value);
            string material;
            double thickness;
            iss >> material >> thickness;
            if (iss.fail()) throw SCATMECH_exception("Invalid parameter for dielectric_stack");
            variable.grow(material,thickness);
        } else if (parameter == "wash") {
            variable.wash();
        }
    }

    template <>
    string
    ModelParameterGet(dielectric_stack& variable,const string& parameter)
    {
        ostringstream oss;
        oss << variable;
        return oss.str();
    }

    template <>
    void
    ModelParameterAskUser(dielectric_stack& value,const string& prompt)
    {
        value = dielectric_stack::AskUser(prompt);
    }

    void Register(const StackModel* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(StackModel);
			Register_Model(No_StackModel);
            Register_Model(Stack_StackModel);
            Register_Model(SingleFilm_StackModel);
			Register_Model(DoubleFilm_StackModel);
            Register_Model(GradedFilm_StackModel);

        }
    }

	DEFINE_VIRTUAL_MODEL(StackModel,Model,"Abstract model for film stacks, allowing parametric variation of dielectric stacks");
	
	DEFINE_MODEL(No_StackModel,StackModel,"An empty film stack");

	DEFINE_MODEL(Stack_StackModel,StackModel,"StackModel that takes a dielectric stack");
	DEFINE_PARAMETER(Stack_StackModel,dielectric_stack,stack,"Film stack","",0xFF);

	DEFINE_MODEL(SingleFilm_StackModel,StackModel,"A single film");
	DEFINE_PARAMETER(SingleFilm_StackModel,dielectric_function,material,"Material","(1.5,0)",0xFF);
	DEFINE_PARAMETER(SingleFilm_StackModel,double,thickness,"Thickness [um]","0.1",0xFF);

	DEFINE_MODEL(DoubleFilm_StackModel,StackModel,"A double film");
	DEFINE_PARAMETER(DoubleFilm_StackModel,dielectric_function,material1,"Material deposited first","(1.5,0)",0xFF);
	DEFINE_PARAMETER(DoubleFilm_StackModel,double,thickness1,"Thickness of first film [um]","0.1",0xFF);
	DEFINE_PARAMETER(DoubleFilm_StackModel,dielectric_function,material2,"Material deposited second","(1.5,0)",0xFF);
	DEFINE_PARAMETER(DoubleFilm_StackModel,double,thickness2,"Thickness of second film[um]","0.1",0xFF);

	DEFINE_MODEL(GradedFilm_StackModel,StackModel,"A graded index film");
	DEFINE_PARAMETER(GradedFilm_StackModel,dielectric_function,start,"Material closest to substrate","(1.5,0)",0xFF);
	DEFINE_PARAMETER(GradedFilm_StackModel,dielectric_function,end,"Material farthest from substrate","(1,0)",0xFF);
	DEFINE_PARAMETER(GradedFilm_StackModel,double,thickness,"Thickness [um]","0.1",0xFF);
	DEFINE_PARAMETER(GradedFilm_StackModel,int,steps,"Number of steps","5",0xFF);

} // namespace SCATMECH




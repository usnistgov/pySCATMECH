//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: gcross.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "gcross.h"
#include "random.h"
#include <list>
#include <algorithm>

using namespace std;

namespace SCATMECH {


    ///
    /// The following are only available to this file...
    ///
    namespace {

        ///
        /// read_paren reads a token in the stream.  A token is a character string, possibly with
        /// parentheses, that has either no whitespace, or has matched parentheses.
        ///
        string read_paren(istream& is)
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

        ///
        /// The equal struct defines a comparison between two numbers that yields true
        /// if the numbers differ by a tolerance and false otherwise...
        ///
        struct equal  {
            public:
                equal(double tolerance) {
                    tol=tolerance;
                }
                bool operator()(const double& a,const double& b) const {
                    if (fabs(a-b)>tol) return false;
                    return true;
                }
            private:
                double tol;
        };

    }

    ///
    /// The following function returns true if the line segment defined by endpoints
    /// s1 and s2 intersects the triangle defined by corners t1, t2,
    /// and t3.  If the intersection exists, then intersection will be the point of intersection.
    ///
    bool TriangleLineIntersection(const Vector& t1, const Vector& t2, const Vector& t3,
                                  const Vector& s1,const Vector& s2,
                                  Vector& intersection)
    {
        double denom = s1.y*t1.z*t2.x - s2.y*t1.z*t2.x - s1.x*t1.z*t2.y + s2.x*t1.z*t2.y - s1.y*t1.x*t2.z +
                       s2.y*t1.x*t2.z + s1.x*t1.y*t2.z - s2.x*t1.y*t2.z - s1.y*t1.z*t3.x + s2.y*t1.z*t3.x +
                       s1.y*t2.z*t3.x - s2.y*t2.z*t3.x + s1.x*t1.z*t3.y - s2.x*t1.z*t3.y - s1.x*t2.z*t3.y +
                       s2.x*t2.z*t3.y +
                       s2.z*(t1.y*t2.x - t1.x*t2.y - t1.y*t3.x + t2.y*t3.x + t1.x*t3.y - t2.x*t3.y) +
                       s1.z*(t1.x*t2.y - t2.y*t3.x + t1.y*(-t2.x + t3.x) - t1.x*t3.y +
                             t2.x*t3.y) + ((s1.y - s2.y)*(t1.x - t2.x) - (s1.x - s2.x)*(t1.y - t2.y))*t3.z;

        if (denom==0.) return false;

        double qa = s1.x*s2.z*t1.y - s1.x*s2.y*t1.z - s2.z*t1.y*t3.x + s2.y*t1.z*t3.x -
                    s1.x*s2.z*t3.y + s2.z*t1.x*t3.y + s1.x*t1.z*t3.y - s2.x*t1.z*t3.y +
                    s1.z*(-(s2.x*t1.y) + s2.y*(t1.x - t3.x) + t1.y*t3.x + s2.x*t3.y - t1.x*t3.y) +
                    s1.x*s2.y*t3.z - s2.y*t1.x*t3.z - s1.x*t1.y*t3.z + s2.x*t1.y*t3.z +
                    s1.y*(-(s2.z*t1.x) + s2.x*t1.z + s2.z*t3.x - t1.z*t3.x - s2.x*t3.z +
                          t1.x*t3.z);
        qa /= denom;

        double qb = -(s1.x*s2.z*t1.y) + s1.x*s2.y*t1.z + s2.z*t1.y*t2.x -
                    s2.y*t1.z*t2.x + s1.x*s2.z*t2.y - s2.z*t1.x*t2.y - s1.x*t1.z*t2.y +
                    s2.x*t1.z*t2.y +
                    s1.z*(s2.x*t1.y - t1.y*t2.x + s2.y*(-t1.x + t2.x) - s2.x*t2.y +
                          t1.x*t2.y) + (-(s1.x*s2.y) + s2.y*t1.x + s1.x*t1.y - s2.x*t1.y)*t2.z +
                    s1.y*(s2.z*t1.x - s2.x*t1.z - s2.z*t2.x + t1.z*t2.x + s2.x*t2.z -
                          t1.x*t2.z);

        qb /= denom;

        double ra =  -(s1.x*t1.z*t2.y) + s1.x*t1.y*t2.z + t1.z*t2.y*t3.x -
                     t1.y*t2.z*t3.x + s1.x*t1.z*t3.y - t1.z*t2.x*t3.y - s1.x*t2.z*t3.y +
                     t1.x*t2.z*t3.y +
                     s1.z*(t1.x*t2.y - t2.y*t3.x + t1.y*(-t2.x + t3.x) - t1.x*t3.y +
                           t2.x*t3.y) + (-(s1.x*t1.y) + t1.y*t2.x + s1.x*t2.y - t1.x*t2.y)*t3.z +
                     s1.y*(t1.z*t2.x - t1.x*t2.z - t1.z*t3.x + t2.z*t3.x + t1.x*t3.z - t2.x*t3.z);
        ra /= denom;

        intersection = s1 + ra*(s2-s1);

        const double small = 1E-10;
        if (qa>= -small && qb>= -small && qa<1.+small && qb<1.+small && qa+qb<1.+small) return true;
        else return false;
    }

    ///
    /// The following function returns true if the line segment passing through s
    /// intersects the triangle defined by corners t1, t2,
    /// and t3.  If the intersection exists, then intersection will be the point of intersection.
    /// The line is in the (direction==1) x, (direction==2) y, or (direction==3) z direction.
    ///
    bool TriangleLineIntersection(const Vector& t1, const Vector& t2, const Vector& t3,
                                  const Vector& s,int direction,
                                  Vector& intersection)
    {
        double ra,rb;
        switch (direction) {
            case 1:
                ra = (s.z*(t1.y - t3.y) + t1.z*t3.y - t1.y*t3.z +
                      s.y*(-t1.z + t3.z))/(-(t2.z*t3.y) + t1.z*(-t2.y + t3.y) +
                                           t1.y*(t2.z - t3.z) + t2.y*t3.z);
                rb = (s.z*t1.y - s.y*t1.z - s.z*t2.y +
                      t1.z*t2.y + s.y*t2.z - t1.y*t2.z)/(t1.z*t2.y - t1.y*t2.z - t1.z*t3.y +
                              t2.z*t3.y + t1.y*t3.z - t2.y*t3.z);
                break;
            case 2:
                ra = (s.z*(t1.x - t3.x) + t1.z*t3.x - t1.x*t3.z +
                      s.x*(-t1.z + t3.z))/(-(t2.z*t3.x) + t1.z*(-t2.x + t3.x) +
                                           t1.x*(t2.z - t3.z) + t2.x*t3.z);
                rb = (s.z*t1.x - s.x*t1.z - s.z*t2.x +
                      t1.z*t2.x + s.x*t2.z - t1.x*t2.z)/(t1.z*t2.x - t1.x*t2.z - t1.z*t3.x +
                              t2.z*t3.x + t1.x*t3.z - t2.x*t3.z);
                break;
            case 3:
                ra = (s.y*(t1.x - t3.x) + t1.y*t3.x - t1.x*t3.y +
                      s.x*(-t1.y + t3.y))/(-(t2.y*t3.x) + t1.y*(-t2.x + t3.x) +
                                           t1.x*(t2.y - t3.y) + t2.x*t3.y);
                rb = (s.y*t1.x - s.x*t1.y - s.y*t2.x +
                      t1.y*t2.x + s.x*t2.y - t1.x*t2.y)/(t1.y*t2.x - t1.x*t2.y - t1.y*t3.x +
                              t2.y*t3.x + t1.x*t3.y - t2.x*t3.y);
                break;
            default:
                throw SCATMECH_exception("Invalid direction in TriangleLineIntersection");
        }
        intersection = t1+ra*(t2-t1)+rb*(t3-t1);

        const double small = 1E-10;
        if (ra>= -small && rb>= -small && ra<1.+small && rb<1.+small && ra+rb<1.+small) return true;
        else return false;
    }

    ///
    /// find_epsilon is a function which returns the dielectric constant at location s.
    /// If throwerror is true and the function finds an inconsistency, then an exception is thrown.
    ///
	Generic_CrossGrating::material 
		Generic_CrossGrating::find_epsilon(const Vector& s,bool throwerror)
    {
        //
        // Make lists of boundaries passed in the x, y, and z directions, respectively...
        //
        list<bounds> boundslistx;
        list<bounds> boundslisty;
        list<bounds> boundslistz;

        for (int j=0; j<(int)boundaries.size(); ++j) {
            Vector &v1=boundaries[j].vertex1;
            Vector &v2=boundaries[j].vertex2;
            Vector &v3=boundaries[j].vertex3;
            Vector intersection;
            if ( TriangleLineIntersection(v1,v2,v3,s,1,intersection) ) {
                bounds b;
                b.x = intersection.x;
                b.n = boundaries[j].n;
                b.text = boundaries[j].text;

                if (cross(v2-v1,v3-v2).x>=0) {
                    b.mat1 = boundaries[j].mat1;
                    b.mat2 = boundaries[j].mat2;
                } else {
                    b.mat1 = boundaries[j].mat2;
                    b.mat2 = boundaries[j].mat1;
                }
                boundslistx.push_back(b);
            }
            if ( TriangleLineIntersection(v1,v2,v3,s,2,intersection) ) {
                bounds b;
                b.x = intersection.y;
                b.n = boundaries[j].n;
                b.text = boundaries[j].text;

                if (cross(v2-v1,v3-v2).y>=0) {
                    b.mat1 = boundaries[j].mat1;
                    b.mat2 = boundaries[j].mat2;
                } else {
                    b.mat1 = boundaries[j].mat2;
                    b.mat2 = boundaries[j].mat1;
                }
                boundslisty.push_back(b);
            }
            if ( TriangleLineIntersection(v1,v2,v3,s,3,intersection) ) {
                bounds b;
                b.x = intersection.z;
                b.n = boundaries[j].n;
                b.text = boundaries[j].text;

                if (cross(v2-v1,v3-v2).z>=0) {
                    b.mat1 = boundaries[j].mat1;
                    b.mat2 = boundaries[j].mat2;
                } else {
                    b.mat1 = boundaries[j].mat2;
                    b.mat2 = boundaries[j].mat1;
                }
                boundslistz.push_back(b);
            }
        }

        // Sort uniquely all the boundary points...
        boundslistx.sort();
        boundslisty.sort();
        boundslistz.sort();
        boundslistx.unique();
        boundslisty.unique();
        boundslistz.unique();

        // matxa and matxb are the materials at the two ends of the line segment going through s in the x-direction...
        material matxa(0.,0.,0.,0.,0.,0.),matxb(0.,0.,0.,0.,0.,0.);
        if (!boundslistx.empty()) {
            matxa=boundslistx.back().mat2;
            matxb=boundslistx.front().mat1;
        }
        // Find the material at s going in the increasing x direction...
        for (list<bounds>::iterator p=boundslistx.begin(); p!=boundslistx.end(); ++p) {
            if (p->x <= s.x) matxa = p->mat2;
        }
        // Find the material at s going in the decreasing x direction...
        for (list<bounds>::reverse_iterator _p=boundslistx.rbegin(); _p!=boundslistx.rend(); ++_p) {
            if (_p->x >= s.x) matxb = _p->mat1;
        }


        // Repeat in the y direction...
        material matya(0.,0.,0.,0.,0.,0.),matyb(0., 0., 0., 0., 0., 0.);
        if (!boundslisty.empty()) {
            matya=boundslisty.back().mat2;
            matyb=boundslisty.front().mat1;
        }
        for (list<bounds>::iterator q=boundslisty.begin(); q!=boundslisty.end(); ++q) {
            if (q->x <= s.y) matya = q->mat2;
        }
        for (list<bounds>::reverse_iterator _q=boundslisty.rbegin(); _q!=boundslisty.rend(); ++_q) {
            if (_q->x >= s.y) matyb = _q->mat1;
        }

        // repeat in the z direction...
        material matza=mat_t,matzb=mat_i;
        if (boundslistz.empty()) error("Cannot find boundary from medium_i to medium_t at " + to_string(s));

        for (list<bounds>::iterator r=boundslistz.begin(); r!=boundslistz.end(); ++r) {
            if (r->x <= s.z) matza = r->mat2;
        }
        for (list<bounds>::reverse_iterator _r=boundslistz.rbegin(); _r!=boundslistz.rend(); ++_r) {
            if (_r->x >= s.z) matzb = _r->mat1;
        }

        if (throwerror) {
            if (matxa!=matxb) {
                ostringstream err;
                err << "Inconsistent optical properties found at (x," << s.y << ',' << s.z << ")." << endl
                    << "The boundaries found were:" << endl;
                for (list<bounds>::iterator p=boundslistx.begin(); p!=boundslistx.end(); ++p) {
                    err << "x=" << p->x << " boundary #" << p->n
                        << " (defined as " << *(p->text) << ") from eps1 " << p->mat1.eps1 << " to eps1 " << p->mat2.eps2 << endl;
                }
                err << "matxa = " << matxa.eps1 << ", matxb = " << matxb.eps1 << endl;
                err << "matya = " << matxa.eps1 << ", matyb = " << matxb.eps1 << endl;
                err << "matza = " << matxa.eps1 << ", matzb = " << matxb.eps1 << endl;
                err << endl;
                error(err.str());
            }

            if (matya!=matyb) {
                ostringstream err;
                err << "Inconsistent optical properties found at (" << s.x << ",y," << s.z << ")." << endl
                    << "The boundaries found were:" << endl;
                for (list<bounds>::iterator p=boundslisty.begin(); p!=boundslisty.end(); ++p) {
                    err << "y=" << p->x << " boundary #" << p->n
                        << " (defined as " << *(p->text) << ") from eps1 " << p->mat1.eps1 << " to eps1 " << p->mat2.eps1 << endl;
                }
                err << "matxa = " << matxa.eps1 << ", matxb = " << matxb.eps1 << endl;
                err << "matya = " << matxa.eps1 << ", matyb = " << matxb.eps1 << endl;
                err << "matza = " << matxa.eps1 << ", matzb = " << matxb.eps1 << endl;
                err << endl;
                error(err.str());
            }

            if (boundslistz.size()>1 && matza!=matzb) {
                ostringstream err;
                err << "Inconsistent optical properties found at (" << s.x << ',' << s.y << ",z)." << endl
                    << "The boundaries found were:" << endl;
                for (list<bounds>::iterator p=boundslistz.begin(); p!=boundslistz.end(); ++p) {
                    err << "z=" << p->x << " boundary #" << p->n
                        << " (defined as " << *(p->text) << ") from eps1 " << p->mat1.eps1 << " to eps1 " << p->mat2.eps1 << endl;
                }
                err << "matxa = " << matxa.eps1 << ", matxb = " << matxb.eps1 << endl;
                err << "matya = " << matxa.eps1 << ", matyb = " << matxb.eps1 << endl;
                err << "matza = " << matxa.eps1 << ", matzb = " << matxb.eps1 << endl;

                err << endl;
                error(err.str());
            }
            if (!boundslistx.empty() && matxa!=matza) {
                ostringstream err;
                err << "Inconsistent optical properties found at " << s << " in x and z scans. " << endl
                    << "x scan said the material should be " << matxa.eps1 << endl
                    << "z scan said the material should be " << matza.eps1 << endl;

                err << "matxa = " << matxa.eps1 << ", matxb = " << matxb.eps1 << endl;
                err << "matya = " << matxa.eps1 << ", matyb = " << matxb.eps1 << endl;
                err << "matza = " << matxa.eps1 << ", matzb = " << matxb.eps1 << endl;

                err << endl;
                error(err.str());
            }
            if (!boundslisty.empty() && matya!=matza) {
                ostringstream err;
                err << "Inconsistent optical properties found at " << s << " in y and z scans. " << endl
                    << "y scan said the material should be " << matya.eps1 << endl
                    << "z scan said the material should be " << matza.eps1 << endl;
                err << "matxa = " << matxa.eps1 << ", matxb = " << matxb.eps1 << endl;
                err << "matya = " << matxa.eps1 << ", matyb = " << matxb.eps1 << endl;
                err << "matza = " << matxa.eps1 << ", matzb = " << matxb.eps1 << endl;

                err << endl;
                error(err.str());
            }
        }
        return matza;
    }
	
	COMPLEX Generic_CrossGrating::stringToEpsilon(const string& expression, Generic_CrossGrating::varsmap& vars)
	{
		string expressiona;
		try {
			expressiona = Evaluator(expression, vars).ResultString();
		}
		catch (SCATMECH_exception&) {
			expressiona = expression;
		}
		dielectric_function df(expressiona);
		return epsilon(df);
	}


    void Generic_CrossGrating::setup()
    {
        CrossGrating::setup();
		
		isotropic = true;
		nonmagnetic = true;

        string fname = find_file(filename);
        ifstream_with_comments file(fname.c_str());
        if (!file) error("Cannot open file");

        int i,j;
		matmap.clear();

		//position.clear();
        thickness.clear();
        //material.clear();
        vars.clear();
        boundaries.clear();

        string stemp;

        PointerCollector<string> StringTrash;

        //*****************************************************************
        //*
        //* First block: read the parameters ...
        //*
        //*****************************************************************

        file >> stemp;
        if (stemp!="CROSSGRATING") error("Expected CROSSGRATING label");

        // Read and check the section heading...
        file >> stemp;
        if (stemp!="PARAMETERS") error("Expected PARAMETERS label");

        // Evaluate the values passed through the parameter pstring...
        Evaluator pe(pstring);

        vars["d1"] = d1;
        vars["d2"] = d2;
        vars["grid1"] = grid1;
        vars["grid2"] = grid2;
        vars["lambda"] = lambda;
        vars["zeta"] = zeta;

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

                if (file.fail()) error("Error reading file in WORKING section");

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

        if (stemp!="MATERIALS") error("Expected MATERIALS label");

		matmap["medium_i"] = mat_i = material(epsilon(medium_i), epsilon(medium_i), epsilon(medium_i), 1., 1., 1.);
		matmap["medium_t"] = mat_t = material(epsilon(medium_t), epsilon(medium_t), epsilon(medium_t), 1., 1., 1.);

        while (1) {
            // Read material index...
            file >> stemp;
            if (stemp=="END") break;

			if (stemp == "ANISO") {
				isotropic = false;
				file >> stemp;
				if (file.fail()) error("Error reading file in MATERIALS section for ANISO material");
				string smaterial1 = read_paren(file);
				string smaterial2 = read_paren(file);
				string smaterial3 = read_paren(file);
				if (file.fail()) error("Error reading file in MATERIALS section for ANISO material");

				COMPLEX e1 = stringToEpsilon(smaterial1, vars);
				COMPLEX e2 = stringToEpsilon(smaterial2, vars);
				COMPLEX e3 = stringToEpsilon(smaterial3, vars);

				matmap[stemp] = material(e1, e2, e3, 1., 1., 1.);

			} else if (stemp == "MAGNETIC") {
				nonmagnetic = false;

				file >> stemp;
				if (file.fail()) error("Error reading file in MATERIALS section for MAGNETIC material");
				string smaterial1 = read_paren(file);
				string smaterial2 = read_paren(file);
				if (file.fail()) error("Error reading file in MATERIALS section for MAGNETIC material");

				COMPLEX e = stringToEpsilon(smaterial1, vars);
				COMPLEX m = stringToEpsilon(smaterial2, vars);

				matmap[stemp] = material(e,e,e,m,m,m);

			} else if (stemp == "ANISOMAGNETIC") {
				isotropic = false;
				nonmagnetic = false;

				file >> stemp;
				if (file.fail()) error("Error reading file in MATERIALS section for ANISOMAGNETIC material");
				string smaterial1 = read_paren(file);
				string smaterial2 = read_paren(file);
				string smaterial3 = read_paren(file);
				string smaterial4 = read_paren(file);
				string smaterial5 = read_paren(file);
				string smaterial6 = read_paren(file);
				if (file.fail()) error("Error reading file in MATERIALS section for ANISOMAGNETIC material");

				COMPLEX e1 = stringToEpsilon(smaterial1, vars);
				COMPLEX e2 = stringToEpsilon(smaterial2, vars);
				COMPLEX e3 = stringToEpsilon(smaterial3, vars);
				COMPLEX m1 = stringToEpsilon(smaterial4, vars);
				COMPLEX m2 = stringToEpsilon(smaterial5, vars);
				COMPLEX m3 = stringToEpsilon(smaterial6, vars);
	
				matmap[stemp] = material(e1, e2, e3, m1, m2, m3);

			} else { // Isotropic material...
				// Read material dielectric function...
				string smaterial = read_paren(file);
	
				if (file.fail()) error("Error reading file in MATERIALS section");
	
				COMPLEX e = stringToEpsilon(smaterial, vars);
	
				matmap[stemp] = material(e, e, e, 1., 1., 1.);
			}
        }


        //*****************************************************************
        //*
        //* Third block: read vertices...
        //*
        //*****************************************************************

        file >> stemp;
        if (stemp!="VERTICES") error("Expected MATERIALS label");


        vertices.clear();

        zeta = vars["zeta"];
        coszeta = cos(zeta*deg);
        sinzeta = sin(zeta*deg);

        Vector x1(1,0,0);
        Vector x2(sin(zeta*deg),cos(zeta*deg),0);
        Vector x3(0,0,1);
        Vector x_1(1,-tan(zeta*deg),0);
        Vector x_2(0,1./cos(zeta*deg),0);
        Vector x_3(0,0,1);

        Matrix convert = outer(Vector(1,0,0),x_1)+outer(Vector(0,1,0),x_2)+outer(Vector(0,0,1),x_3);
        Matrix convert_back = outer(x1,Vector(0,0,1)) + outer(x2,Vector(0,1,0)) + outer(x3,Vector(0,0,1));

        while (1) {
            file >> stemp;
            if (stemp=="END") break;
            Vector vertex = get_triplet_value(file);
            if (file.fail()) error("Error reading file in VERTICES section");

            verticesmap::iterator a = vertices.find(stemp);

            if (a!=vertices.end()) error("Duplicate vertex label: " + stemp);

            vertex = convert*vertex;

            vertices.insert(a,verticesmap::value_type(stemp,vertex));
        }

        //*****************************************************************
        //*
        //* Fourth block: read boundaries...
        //*
        //*****************************************************************
        file >> stemp;
        if (stemp!="BOUNDARIES") error("Expected BOUNDARIES label");

        int ibound=1;
        while (!file.eof()) {
            string vertex1,vertex2,vertex3,vertex4;
            string imat1,imat2;

            // Read first vertex...
            file >> vertex1;

            if (file.fail()) error("Error reading file in BOUNDARIES section");

            if (vertex1=="END") break;

            if (vertex1=="QUAD") {
                file >> vertex1;
                if (file.fail()) error("Error reading file at vertex #1");
                verticesmap::iterator v1 = vertices.find(vertex1);
                if (v1==vertices.end()) error("Unknown vertex label: " + vertex1);

                // Read second vertex...
                file >> vertex2;
                if (file.fail()) error("Error reading file at vertex #2");
                verticesmap::iterator v2 = vertices.find(vertex2);
                if (v2==vertices.end()) error("Unknown vertex label: " + vertex2);

                // Read third vertex...
                file >> vertex3;
                if (file.fail()) error("Error reading file at vertex #3");
                verticesmap::iterator v3 = vertices.find(vertex3);
                if (v3==vertices.end()) error("Unknown vertex label: " + vertex3);

                // Read fouth vertex...
                file >> vertex4;
                if (file.fail()) error("Error reading file at vertex #4");
                verticesmap::iterator v4 = vertices.find(vertex4);
                if (v4==vertices.end()) error("Unknown vertex label: " + vertex4);

                // Read left material index...
                file >> imat1;
                if (file.fail()) error("Error reading file at material #1");
                Materialmap::iterator m1 = matmap.find(imat1);
                if (m1==matmap.end()) error ("Unknown material label: " + imat1);

                // Read right material index...
                file >> imat2;
                if (file.fail()) error("Error reading file at material #2");
                Materialmap::iterator m2 = matmap.find(imat2);
                if (m2==matmap.end()) error ("Unknown material: " + imat2);

                string *text = StringTrash.New();
                *text = "QUAD " + vertex1 + " " + vertex2 + " " + vertex3 + " " + vertex4 + " " + imat1 + " " + imat2;

                // Load into a triangle boundary...
                boundary tb1;
                tb1.n=ibound;
                tb1.text=text;
                tb1.vertex1=v1->second;
                tb1.vertex2=v2->second;
                tb1.vertex3=v3->second;
				tb1.mat1 = matmap[imat1];
				tb1.mat2 = matmap[imat2];

				// Load into a triangle boundary...
                boundary tb2;
                tb2.n=ibound;
                tb2.text=text;
                tb2.vertex1=v1->second;
                tb2.vertex2=v3->second;
                tb2.vertex3=v4->second;
				tb2.mat1 = matmap[imat1];
				tb2.mat2 = matmap[imat2];

				// Store in list of boundaries...
                boundaries.push_back(tb1);
                boundaries.push_back(tb2);
                ++ibound;

            } else {

                if (file.fail()) error("Error reading file at vertex #1");
                verticesmap::iterator v1 = vertices.find(vertex1);
                if (v1==vertices.end()) error("Unknown vertex label: " + vertex1);

                // Read second vertex...
                file >> vertex2;
                if (file.fail()) error("Error reading file at vertex #2");
                verticesmap::iterator v2 = vertices.find(vertex2);
                if (v2==vertices.end()) error("Unknown vertex label: " + vertex2);

                // Read third vertex...
                file >> vertex3;
                if (file.fail()) error("Error reading file at vertex #3");
                verticesmap::iterator v3 = vertices.find(vertex3);
                if (v3==vertices.end()) error("Unknown vertex label: " + vertex3);

                // Read left material index...
                file >> imat1;
                if (file.fail()) error("Error reading file at material #1");
                Materialmap::iterator m1 = matmap.find(imat1);
                if (m1==matmap.end()) error ("Unknown material label: " + imat1);

                // Read right material index...
                file >> imat2;
                if (file.fail()) error("Error reading file at material #2");
                Materialmap::iterator m2 = matmap.find(imat2);
                if (m2==matmap.end()) error ("Unknown material: " + imat2);

                string *text = StringTrash.New();
                *text = vertex1 + " " + vertex2 + " " + vertex3 + " " + imat1 + " " + imat2;

                // Load into a triangle boundary...
                boundary tb;
                tb.n=ibound;
                tb.text=text;
                tb.vertex1=v1->second;
                tb.vertex2=v2->second;
                tb.vertex3=v3->second;
				tb.mat1 = matmap[imat1];
				tb.mat2 = matmap[imat2];

                // Store in list of boundaries...
                boundaries.push_back(tb);
                ++ibound;
            }
        }

        //*****************************************************************
        //*
        //* Get a list of all vertical and horizontal positions...
        //*
        //*****************************************************************
        list<double> X,Y,Z;
        for (i=0; i<(int)boundaries.size(); ++i) {
            X.push_back(boundaries[i].vertex1.x);
            X.push_back(boundaries[i].vertex2.x);
            X.push_back(boundaries[i].vertex3.x);
            Y.push_back(boundaries[i].vertex1.y);
            Y.push_back(boundaries[i].vertex2.y);
            Y.push_back(boundaries[i].vertex3.y);
            Z.push_back(boundaries[i].vertex1.z);
            Z.push_back(boundaries[i].vertex2.z);
            Z.push_back(boundaries[i].vertex3.z);
        }

        // Sort them ...

        X.sort();
        Y.sort();
        Z.sort();

        // Total period of the grating...
        d1 = X.back()-X.front();
        d2 = Y.back()-Y.front();

        grid1 = (int)(vars["grid1"]);
        grid2 = (int)(vars["grid2"]);

        // Total height of grating...
        double height = Z.back()-Z.front();

        double firstx = X.front();
        double prevx = firstx;
        double firsty = Y.front();
        double prevy = firsty;
        double firstz = Z.front();
        double prevz = firstz;

        // Throw away duplicates...
        list<double>::iterator it;
        it = unique(X.begin(),X.end(),equal(d1*1E-10));
        X.erase(it,X.end());

        it = unique(Y.begin(),Y.end(),equal(d2*1E-10));
        Y.erase(it,Y.end());

        it = unique(Z.begin(),Z.end(),equal(height*1E-10));
        Z.erase(it,Z.end());

        list<double>::iterator p;

        vector<double> td;
        double total_td=0;

        //*****************************************************************
        //*
        //* Find maximum traversed distance for each sublayer...
        //*
        //* Note: Sublayers refer to regions between different given
        //*       y-values, which are further divided into microlayers.
        //*
        //*****************************************************************
        for (p=++(Z.begin()); p!=Z.end(); ++p) {
            // Get a vertical position in the middle of the sublayer...
            double subheight=(*p-prevz);
            double h=prevz+subheight/2.;

            double max=0;
            for (j=0; j<(int)boundaries.size(); ++j) {
                Vector v1=boundaries[j].vertex1;
                Vector v2=boundaries[j].vertex2;
                Vector v3=boundaries[j].vertex3;

                // If segment crosses vertical position...
                if ((v1.z>=h && v2.z<h) || (v1.z<=h && v2.z>h) ||
                        (v2.z>=h && v3.z<h) || (v2.z<=h && v3.z>h) ||
                        (v3.z>=h && v1.z<h) || (v3.z<=h && v1.z>h)) {
                    Vector v2v1 = v2-v1;
                    Vector v3v1 = v3-v1;
                    Vector v = perpto(v2v1,v3v1);

                    // Calculate the traversed distance...
                    double tdist = fabs(tan(asin(v.z)));

                    tdist = subheight*tdist;

                    // Check to see if the materials on both sides are the same...
                    if (boundaries[j].mat1==boundaries[j].mat2) tdist=0;
                    // If this is a bigger distance, save it...
                    if (max<tdist) max=tdist;
                }
            }
            // Store maximum traversed distance...
            td.push_back(max);
            total_td += max;
            prevz=*p;
        }

        //*****************************************************************
        //*
        //* Check for consistency in grating...
        //*
        //*****************************************************************

        {
            int itry =0;
            bool good = false;
            do {
                try {
                    Random_Number random = Random_Number::uniform(0.1,0.9);
                    random();
                    double x=X.front();
                    for (list<double>::iterator p=++(X.begin()); p!=X.end(); ++p) {
                        double rr = random();
                        double xx1=(*p)*rr+x*(1-rr);
                        double xx2=(*p)*(1-rr)+x*rr;

                        double y=Y.front();
                        for (list<double>::iterator q=++(Y.begin()); q!=Y.end(); ++q) {
                            double rr = random();
                            double yy1=(*q)*rr+y*(1-rr);
                            double yy2=(*q)*(1-rr)+y*rr;

                            double z=Z.front();
                            for (list<double>::iterator r=++(Z.begin()); r!=Z.end(); ++r) {
                                double rr = random();
                                double zz1 = (*r)*rr+z*(1-rr);
                                double zz2 = (*r)*(1-rr)+z*rr;

                                bool throwerror = true;
                                find_epsilon(Vector(xx1,yy1,zz1),throwerror);
                                find_epsilon(Vector(xx1,yy1,zz2),throwerror);
                                find_epsilon(Vector(xx1,yy2,zz1),throwerror);
                                find_epsilon(Vector(xx1,yy2,zz2),throwerror);
                                find_epsilon(Vector(xx2,yy1,zz1),throwerror);
                                find_epsilon(Vector(xx2,yy1,zz2),throwerror);
                                find_epsilon(Vector(xx2,yy2,zz1),throwerror);
                                find_epsilon(Vector(xx2,yy2,zz2),throwerror);
                                z=*r;
                            }
                            y=*q;
                        }
                        x=*p;
                    }
                    good=true;
                }
                catch(exception&) {
                    ++itry;
                    if (itry>3) throw;
                }
            } while (good==false);
        }

        //*****************************************************************
        //*
        //* Count number of microlayers, so that we can allocate memory...
        //*
        //*****************************************************************
        levels = 0;
        int kk=0;
        int nvert = 0;
        for (p=++(Z.begin()); p!=Z.end(); ++p) {
            if (td[kk]<1E-10) ++nvert;
            int mlayers = (total_td==0) ? 1 : (int)(td[kk]/total_td*nlayers+0.5);
            // Make sure there is at least one microlayer...
            if (mlayers<1) mlayers = 1;
            levels += mlayers;
            // Increment sublayer number...
            ++kk;
        }

		if (isotropic) {
			eps.allocate(grid1, grid2, levels);
		} else {
			eps1.allocate(grid1, grid2, levels);
			eps2.allocate(grid1, grid2, levels);
			eps3.allocate(grid1, grid2, levels);
		}

		if (!nonmagnetic) {
			if (isotropic) {
				mu.allocate(grid1, grid2, levels);
			}
			else {
				mu1.allocate(grid1, grid2, levels);
				mu2.allocate(grid1, grid2, levels);
				mu3.allocate(grid1, grid2, levels);
			}
		}
        message("Number of levels: " + to_string(levels) + "\n");

        //*****************************************************************
        //*
        //* Calculate and store epsilon for grid ...
        //*
        //*****************************************************************
        int k=0;
        kk=0;
        prevz=firstz;
        for (p=++(Z.begin()); p!=Z.end(); ++p) {
            double subheight=(*p-prevz);
            // Scale number of microlayers in sublayer by maximum traversed distance....
            int mlayers = (total_td==0) ? 1 : (int)(td[kk]/total_td*(nlayers)+0.5);
            // Make sure there is at least one microlayer...
            if (mlayers<1) mlayers = 1;
            for (i=0; i<mlayers; ++i) {
                // Calculate vertical position in middle of microlayer...
                double h=prevz+subheight*(i+0.5)/mlayers;

                for (int jj=0; jj<grid1; ++jj) {

                    list<bounds> boundslist;

                    double x = firstx+(jj+0.5)*d1/grid1;
                    Vector point(x,0,h);

                    for (j=0; j<(int)boundaries.size(); ++j) {
                        Vector v1=boundaries[j].vertex1;
                        Vector v2=boundaries[j].vertex2;
                        Vector v3=boundaries[j].vertex3;
                        Vector intersection;
                        bool ts = TriangleLineIntersection(v1,v2,v3,point,2,intersection);
                        if (ts) {
                            bounds b;
                            b.x = intersection.y;
                            b.n = boundaries[j].n;
                            b.text = boundaries[j].text;

                            if (cross(v2-v1,v3-v1).y>=0) {
                                b.mat1 = boundaries[j].mat1;
                                b.mat2 = boundaries[j].mat2;
                            } else {
                                b.mat1 = boundaries[j].mat2;
                                b.mat2 = boundaries[j].mat1;
                            }
                            boundslist.push_back(b);
                        }
                    }
                    boundslist.sort();
                    boundslist.unique();
                    if (boundslist.empty()) {
                        material e = find_epsilon(Vector(x,Y.front()*0.1111111+Y.back()*0.8888889,h),false);
                        for (int m=0; m<grid2; ++m) {
							if (isotropic) {
								eps(jj + 1, m + 1, (k + i + 1)) = e.eps1;
							} else {
								eps1(jj + 1, m + 1, (k + i + 1)) = e.eps1;
								eps2(jj + 1, m + 1, (k + i + 1)) = e.eps2;
								eps3(jj + 1, m + 1, (k + i + 1)) = e.eps3;
							}
							if (!nonmagnetic) {
								if (isotropic) {
									mu(jj + 1, m + 1, (k + i + 1)) = e.mu1;
								}
								else {
									mu1(jj + 1, m + 1, (k + i + 1)) = e.mu1;
									mu2(jj + 1, m + 1, (k + i + 1)) = e.mu2;
									mu3(jj + 1, m + 1, (k + i + 1)) = e.mu3;
								}
							}
						}
                    } else {
                        int j=0;

                        material mat = boundslist.back().mat2;

                        for (list<bounds>::iterator q=boundslist.begin(); q!=boundslist.end(); ++q) {

                            int knext = int((q->x - firsty)/d2*grid2);

                            for (int m=j; m<knext; ++m) {
								if (isotropic) {
									eps(jj + 1, m + 1, (k + i + 1)) = mat.eps1;
								}
								else {
									eps1(jj + 1, m + 1, (k + i + 1)) = mat.eps1;
									eps2(jj + 1, m + 1, (k + i + 1)) = mat.eps2;
									eps3(jj + 1, m + 1, (k + i + 1)) = mat.eps3;
								}
								if (!nonmagnetic) {
									if (isotropic) {
										mu(jj + 1, m + 1, (k + i + 1)) = mat.mu1;
									}
									else {
										mu1(jj + 1, m + 1, (k + i + 1)) = mat.mu1;
										mu2(jj + 1, m + 1, (k + i + 1)) = mat.mu2;
										mu3(jj + 1, m + 1, (k + i + 1)) = mat.mu3;
									}
								}
                            }

                            mat = q->mat2;
                            j=knext;
                        }
                        for (int m=j; m<grid2; ++m) {
							if (isotropic) {
								eps(jj + 1, m + 1, (k + i + 1)) = mat.eps1;
							}
							else {
								eps1(jj + 1, m + 1, (k + i + 1)) = mat.eps1;
								eps2(jj + 1, m + 1, (k + i + 1)) = mat.eps2;
								eps3(jj + 1, m + 1, (k + i + 1)) = mat.eps3;
							}
							if (!nonmagnetic) {
								if (isotropic) {
									mu(jj + 1, m + 1, (k + i + 1)) = mat.mu1;
								}
								else {
									mu1(jj + 1, m + 1, (k + i + 1)) = mat.mu1;
									mu2(jj + 1, m + 1, (k + i + 1)) = mat.mu2;
									mu3(jj + 1, m + 1, (k + i + 1)) = mat.mu3;
								}
							}
						}
                    }
                }
                thickness.push_back(subheight/mlayers);
            }

            prevz=*p;
            // Increment sublayer number...
            ++kk;
            k+=mlayers;
        }

        // Grating wants information stored from the top microlayer...
        //reverse(thickness.begin(),thickness.end());
        thick.allocate(levels);
        for (i=1; i<=levels; ++i) thick(i) = thickness[i-1];

        CrossGrating::zeta = zeta;
        CrossGrating::d1 = d1;
        CrossGrating::d2 = d2;

        FourierFactorize();
    }

    Vector
    Generic_CrossGrating::
    get_triplet_value(istream& is) const
    {
        string a;
        is >> a;

        if (is.eof()) error("Premature end of file while reading vertex");

        Evaluator eval(a,vars);

        if (eval.NResult()!=3) error("Vertices require three values");

        Vector result(eval.Result(0),eval.Result(1),eval.Result(2));
        return result;
    }

    void Generic_CrossGrating::error(const string& message) const
    {
        ofstream file("bounds.dat");
        for (boundary_vector::const_iterator q=boundaries.begin(); q!=boundaries.end(); ++q) {
            file << q->mat1.eps1 << tab << q->mat2.eps1 << tab << q->vertex1.x << tab << q->vertex1.y << tab << q->vertex1.z << endl
                 << q->mat1.eps1 << tab << q->mat2.eps1 << tab << q->vertex2.x << tab << q->vertex2.y << tab << q->vertex2.z << endl
                 << q->mat1.eps1 << tab << q->mat2.eps1 << tab << q->vertex3.x << tab << q->vertex3.y << tab << q->vertex3.z << endl
                 << q->mat1.eps1 << tab << q->mat2.eps1 << tab << q->vertex1.x << tab << q->vertex1.y << tab << q->vertex1.z << endl << endl;
        }
        ofstream varfile("vardump.dat");
        varfile << "Variables: " << endl << endl;
        for (varsmap::const_iterator p=vars.begin(); p!=vars.end(); ++p) {
            varfile << p->first << " = " << p->second << endl;
        }
        varfile << endl << "Vertices: " << endl << endl;
        for (verticesmap::const_iterator r=vertices.begin(); r!=vertices.end(); ++r) {
            varfile << r->first << " = " << r->second << endl;
        }

        varfile << endl << "Boundaries: " << endl << endl;
        for (boundary_vector::const_iterator s=boundaries.begin(); s!=boundaries.end(); ++s) {
            varfile << s->n << ": " << *(s->text) << " : " << s->mat1.eps1 << tab << s->mat2.eps1 << tab << s->vertex1 << tab << s->vertex2 << tab << s->vertex3  << endl;
        }

        Model::error("(input file = " + filename + "): " + message + " (addtional information sent to file vardump.dat)");
    }

    STRING Generic_CrossGrating::get_parameter_base(const STRING& parameter) const
    {
        if (parameter.substr(0,6)=="param.") {
            const_cast<Generic_CrossGrating*>(this)->SETUP();
            varsmap::const_iterator p = parameters.find(parameter.substr(6));
            if (p!=parameters.end()) return format("%16g",p->second);
        }
        return CrossGrating::get_parameter_base(parameter);
    }

    void Generic_CrossGrating::set_parameter_base(const STRING& parameter, const STRING& value)
    {
        if (parameter.substr(0,6)=="param.") {
            override[parameter.substr(6)] = from_string<double>(value);
            set_recalc(1);
        } else if (parameter=="SaveBoundaries") {
            SETUP();

            ofstream file(value.c_str());
            if (!file) error("Cannot open file: " + value);
            for (boundary_vector::iterator p=boundaries.begin(); p!=boundaries.end(); ++p) {
                file << p->mat1.eps1 << tab << p->mat2.eps1 << tab << p->vertex1.x << tab << p->vertex1.y << tab << p->vertex1.z << endl
                     << p->mat1.eps1 << tab << p->mat2.eps1 << tab << p->vertex2.x << tab << p->vertex2.y << tab << p->vertex2.z << endl
                     << p->mat1.eps1 << tab << p->mat2.eps1 << tab << p->vertex3.x << tab << p->vertex3.y << tab << p->vertex3.z << endl
                     << p->mat1.eps1 << tab << p->mat2.eps1 << tab << p->vertex1.x << tab << p->vertex1.y << tab << p->vertex1.z << endl << endl;
            }
            set_recalc(0);
        } else {
            Gridded_CrossGrating::set_parameter_base(parameter,value);
        }
    }

    void  Generic_CrossGrating::print_parameters(std::ostream& os, const std::string& prefix) const
    {
        Gridded_CrossGrating::print_parameters(os);
        for (varsmap::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            if (prefix.size()!=0) os << prefix << '.';
            os << "param." << p->first << " = " << p->second << " (double) " << endl;
        }
    }

    DEFINE_MODEL(Generic_CrossGrating,Gridded_CrossGrating,"Generic file-based specification of a crossed grating");
    DEFINE_PARAMETER(Generic_CrossGrating,string,filename,"File containing crossed grating description","",0xFF);
    DEFINE_PARAMETER(Generic_CrossGrating,string,pstring,"String containing parameters","",0xFF);
    DEFINE_PARAMETER(Generic_CrossGrating,int,nlayers,"Number of layers","10",0xFF);

}

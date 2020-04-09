//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: scattabl.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include <float.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <list>
#include <map>
#include <limits>
#include "askuser.h"
#include "scatmech.h"
#include "scattabl.h"
#include "scateval.h"
#include <iostream>

using namespace std;

namespace SCATMECH {

    ///
    /// Class TableFile stores the data taken from a file.  The purpose of the class
    /// is to prevent rereading of a file everytime a new Table is created from the same
    /// file.
    ///
    class TableFile {
        public:
            typedef std::vector<double> vectord;
            typedef std::vector<vectord> vvectord;

            TableFile(const std::string& filename);     ///< The constructor
            TableFile() {};
            const std::string& get_name() {
                return name;
            }
            const vectord& operator[](unsigned int i) {
                if (i>=columns.size()) throw SCATMECH_exception("Out of bounds error in TableFile::operator[] (i = " + to_string(i) + ",size = " + to_string(columns.size())+")");
                return columns[i];
            }
            int ncolumns() {
                return columns.size();
            }
        private:
            vvectord columns; ///< The data in the columns
            std::string name; ///< The name of the file
    };

    class FormulaFile {
        public:
            typedef std::vector<double> vectord;
            typedef std::vector<vectord> vvectord;

            FormulaFile(const std::string& filename);  ///< The constructor
            FormulaFile() {};

            const std::string& get_name() {
                return name;   ///< The name of the file
            }
            std::vector<double> get_column(int i, const Evaluator::VMAP& params);

            Evaluator::VMAP get_parameter_map() const; ///< Returns the parameters and their default values

        private:
            std::string name;
            std::string filecontents;
            vvectord lastcolumns;
            Evaluator::VMAP lastparams;

            void error(const std::string& msg) const;
    };

    //
    // Case-independent string compare...
    //
    bool icompare(const std::string& a,const std::string& b)
    {
        if (a.size()!=b.size()) return false;
        for (string::const_iterator p=a.begin(),q=b.begin(); p!=a.end(); ++p,++q) {
            if (tolower(*p)!=tolower(*q)) return false;
        }
        return true;
    }

    TableFile::
    TableFile(const string& filename)
    {
        typedef list<double> listd;
        typedef vector<listd> vlistd;

        vlistd values;

        bool done;
        name = filename;

        string fname = find_file(filename);

        //ifstream file(fname.c_str(),ios::binary);
        ifstream_with_comments file(fname.c_str());

        if (!file) throw SCATMECH_exception("Problem opening file " + fname);

        ifstream::pos_type top=file.tellg();

        // Skip header.  The end of the header is the beginning of any line
        // containing only numbers.  The first row containing only numbers
        // defines the number of columns in the data.
        int ncol;
        for (done=false; !done;) {
            top = file.tellg();
            string line = getstrline(file);

            if (file.fail()) throw SCATMECH_exception("No data in file " + fname);

            istringstream linestream(line);

            bool isheader=false;
            bool isblankline=true;
            ncol=0;
            while (!linestream.fail()) {
                string s;
                linestream >> s;
                if (!linestream.fail()) {
                    istringstream charstream(s);
                    double x;
                    charstream >> x;
                    if (charstream.fail()) isheader=true;
					if (charstream.peek() == ',') charstream.get();
					isblankline=false;
                    ncol++;
                }
            }
            isheader = isheader | isblankline;
            if (!isheader) done=true;
        }
        // Go to the top of the data...
        file.seekg(top,ios::beg);

        // Size arrays...
        columns.resize(ncol);
        values.resize(ncol);
        for (int i=0; i<ncol; ++i) {
            values[i].resize(0);
        }

        // Temporary array for reading columns...
        vector<double> tempvals(ncol);

        // Read the data...
        for (done=false; !done;) {
            string line = getstrline(file);
            if (!file.fail()) {
                istringstream linestream(line);
                int i;
                // Try reading all of the columns...
                for (i=0; i<ncol; ++i)  {
                    double x;
                    linestream >> x;
                    if (!linestream.fail()) {
                        tempvals[i] = x;
                    } else i=ncol;
					if (linestream.peek() == ',') linestream.get();
                }
                // If all of the columns are present, store the row...
                if (i==ncol) {
                    for (i=0; i<ncol; ++i) {
                        values[i].push_back(tempvals[i]);
                    }
                }
            } else done = true;
        }

        // Copy data into arrays...
        for (int i=0; i<ncol; ++i) {
            columns[i].resize(values[i].size());
            listd::iterator q = values[i].begin();
            vectord::iterator p = columns[i].begin();
            for(; p!=columns[i].end(); ++p,++q) {
                *p = *q;
            }
        }
    }

    typedef map<string,TableFile> TableFileMap;
    TableFileMap tablefiles;
    mutex_t maps_mutex;

    typedef map<string,FormulaFile> FormulaFileMap;
    FormulaFileMap formulafiles;

    FormulaFile::FormulaFile(const string& filename)
    {
        string fname = find_file(filename);

        ifstream_with_comments file(fname.c_str());
        if (!file) error("Cannot open file " + filename);

        ostringstream stream;
        stream << file.rdbuf();

        filecontents = stream.str();
        name = filename;
    }

    Evaluator::VMAP FormulaFile::get_parameter_map() const
    {
        Evaluator::VMAP result;
        istringstream file(filecontents);

        string word;
        file >> word;
        if (file.fail()) error("Failure to read file");
        if (word!="PARAMETERS") error("Expecting PARAMETERS");

        while (1) {
            file >> word;
            if (word=="END") break;
            double value;
            file >> value;
            if (file.fail()) error("Expecting value in PARAMETERS section");

            result[word]=value;
        }
        return result;
    }

    vector<double> FormulaFile::get_column(int i, const Evaluator::VMAP &params)
    {
        if (params==lastparams && lastcolumns.size()!=0) {
            if (i<0 || i>(int)lastcolumns.size()) error("Column number " + to_string(i) + " out of range");
            return lastcolumns[i-1];
        }

        lastparams = params;
        lastcolumns.clear();

        istringstream file(filecontents);
        Evaluator::VMAP vars;

        string word;
        file >> word;
        if (file.fail()) error("Failure to read file");
        if (word!="PARAMETERS") error("Cannot find PARAMETERS section");

        while (1) {
            file >> word;
            if (file.fail()) error("Failure to read file");
            if (word=="END") break;
            double value;
            file >> value;
            if (file.fail()) error("Expecting value in PARAMETERS section");

            Evaluator::VMAP::const_iterator it =  params.find(word);
            if (it==params.end()) {
                vars[word]=value;
            } else {
                vars[word]=it->second;
            }
        }

        string variable;
        double begin,end,step;
        file >> word;
        if (word=="FOR") {
            file >> variable;
            if (file.fail()) error("Cannot read FOR parameter");
            file >> word;
            if (file.fail()||word!="FROM") error("FROM expected");
            file >> begin;
            if (file.fail()) error("Cannot read FROM value");
            file >> word;
            if (file.fail()||word!="TO") error("TO expected");
            file >> end;
            if (file.fail()) error("Cannot read TO value");
            file >> word;
            if (file.fail()||word!="BY") error("BY expected");
            file >> step;
            if (file.fail()) error("Cannot read BY value");
        } else error("FOR or AT expected");

        ifstream::pos_type backto=file.tellg();

        vector<double> values;

        for (double x=begin; x<=end; x+=step) {

            vars[variable]=x;
            file.seekg(backto,ios::beg);

            file >> word;

            if (file.fail()||word!="WORKING") error("WORKING expected");

            while (1) {
                file >> word;
                if (file.fail()) error("Failure to read file");
                if (word=="END") break;
                string expression;
                file >> expression;
                if (file.fail()) error("Failure to read file");
                vars[word]=Evaluator(expression,vars);
            }

            file >> word;
            if (word!="VALUES") error("VALUES expected");

            int icol=0;
            while (true) {
                file >> word;
                if (file.fail()) error("Cannot read expression");
                if (word=="END") break;
                double value = Evaluator(word,vars);
                if ((int)lastcolumns.size()<icol+1) lastcolumns.push_back(vectord());
                lastcolumns[icol].push_back(value);
                ++icol;
            }
            if (icol<i) error("Not enough columns in file");
        }
        if (i<0 || i>(int)lastcolumns.size()) error("Column number " + to_string(i) + " out of range");
        return lastcolumns[i-1];
    }

    void FormulaFile::error(const string& msg) const
    {
        throw SCATMECH_exception("Error in FormulaFile: " + msg);
    }

    void ClearTableCache()
    {
        maps_mutex.lock();
        tablefiles.clear();
        formulafiles.clear();
        maps_mutex.unlock();
    }


    Table::
    Table(double singlevalue)
    {
        values.resize(1);
        values[0]=PAIR(0.,singlevalue);
        last.x = 0;
        last.y= singlevalue;
        formulafile=NULL;
        tablefile=NULL;
        formulastring.clear();
        interpolationx=LIN;
        interpolationy=LIN;
    }

    namespace {

        /// The following returns the slope at x2 of a quadratic that goes through
        /// the points (x1,y1), (x2,y2), and (x3,y3)...
        double dquad(double x1,double y1,double x2,double y2,double x3,double y3) {
            return (sqr(x3)*(y1 - y2) - 2*x2*(x3*(y1 - y2) + x1*(y2 - y3)) +
                    sqr(x2)*(y1 - y3) + sqr(x1)*(y2 - y3))/((x1 - x2)*(x1 - x3)*(x2 - x3));
        }

        /// The following returns a cubic, evaluated at x, that goes through
        /// (x1,y1) and (x2,y2), has a slope dy1 at x1, and a slope of dy2 at x2...
        double cubic(double x, double x1, double y1, double dy1, double x2, double y2, double dy2) {
            return (dy2*sqr(x - x1)*(x - x2)*(x1 - x2) + dy1*(x - x1)*sqr(x - x2)*(x1 - x2) -
                    sqr(x - x2)*(2*x - 3*x1 + x2)*y1 + sqr(x - x1)*(2*x + x1 - 3*x2)*y2)/cube(x1 - x2);
        }

        /// The following returns a spline, valid between x2 and x3, ...
        double spline(double x,double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
            return cubic(x,x2,y2,dquad(x1,y1,x2,y2,x3,y3),x3,y3,dquad(x2,y2,x3,y3,x4,y4));
        }
        double leftspline(double x,double x2, double y2, double x3, double y3, double x4, double y4) {
            return cubic(x,x2,y2,0.,x3,y3,dquad(x2,y2,x3,y3,x4,y4));
        }
        double rightspline(double x,double x1, double y1, double x2, double y2, double x3, double y3) {
            return cubic(x,x2,y2,dquad(x1,y1,x2,y2,x3,y3),x3,y3,0.);
        }
    }

    double Table::value(double lambda) const
    {
        if (last.x==lambda) return last.y;

        // Note that this function is declared const, which it is by contract. However, internally, the
        // class keeps track of the previously returned value and this routine can update its data and
        // that previous value. The following two statements serve to break that contract...
        VECTOR &V = const_cast<VECTOR&>(values);

        // If we are using a formula table file, and a parameter has changed, we need to recalculate
        // the table...
        if (formulafile!=NULL && last.x==numeric_limits<double>::max()) {
            // Get the two columns...
            vector<double> x = formulafile->get_column(1,params);
            vector<double> y = formulafile->get_column(icol,params);

            // Store them in values...
            V.resize(x.size());
            for (int i=0; i<(int)x.size(); ++i) {
                V[i]=PAIR(x[i],y[i]);
            }
            sort(V.begin(),V.end());

            // Reset the previous used values...
            last.x = x[0];
            last.y = y[0];
        }

        // if we are using a formula string, then calculate the value
        if (formulastring.size()!=0) {
            if (icol==1) {
                last.x = lambda;
                last.y = lambda;
                return last.y;
            }
            vector<double> fresult = OneLineFunction(formulastring,lambda,params);
            if ((int)fresult.size()<icol-1) throw SCATMECH_exception("No column " + to_string(icol) + " in formula");
            last.x = lambda;
            last.y = fresult[icol-2];
            return last.y;
        }

        if (V.size()==1) {
            last.x = lambda;
            last.y = V[0].second;
            return last.y;
        }

        last.x = lambda;

        int bound1=0;
        int bound2=V.size()-1;

        if (lambda < V[bound1].first) {
            last.y = V[bound1].second;
            return last.y;
        }
        if (lambda > V[bound2].first) {
            last.y = V[bound2].second;
            return last.y;
        }

        while (bound2-bound1>1) {
            int split = (bound1+bound2)/2;
            if (lambda< V[split].first) bound2=split;
            else bound1=split;
        }

        if (interpolationy==SPLINE && V.size()>2) {
            if (interpolationx==LIN) {
                if (bound1==0) {
                    return last.y = leftspline(lambda,V[0].first,V[0].second,V[1].first,V[1].second,V[2].first,V[2].second);
                }
                if (bound2==V.size()-1) {
                    int N = V.size()-1;
                    return last.y = rightspline(lambda,V[N-2].first,V[N-2].second,V[N-1].first,V[N-1].second,V[N].first,V[N].second);
                }
                return last.y = spline(lambda,V[bound1-1].first,V[bound1-1].second,V[bound1].first,V[bound1].second,
                                       V[bound2].first,V[bound2].second,V[bound2+1].first,V[bound2+1].second);
            } else {
                if (bound1==0) {
                    return last.y = leftspline(log(lambda),log(V[0].first),V[0].second,log(V[1].first),V[1].second,log(V[2].first),V[2].second);
                }
                if (bound2==V.size()-1) {
                    int N = V.size()-1;
                    return last.y = rightspline(log(lambda),log(V[N-2].first),V[N-2].second,log(V[N-1].first),V[N-1].second,log(V[N].first),V[N].second);
                }
                return last.y = spline(log(lambda),log(V[bound1-1].first),V[bound1-1].second,log(V[bound1].first),V[bound1].second,
                                       log(V[bound2].first),V[bound2].second,log(V[bound2+1].first),V[bound2+1].second);
            }
        }

        if (interpolationy==LOGSPLINE && V.size()>2) {
            if (interpolationx==LIN) {
                if (bound1==0) {
                    return last.y = exp(leftspline(lambda,V[0].first,log(V[0].second),V[1].first,log(V[1].second),V[2].first,log(V[2].second)));
                }
                if (bound2==V.size()-1) {
                    int N = V.size()-1;
                    return last.y = exp(rightspline(lambda,V[N-2].first,log(V[N-2].second),V[N-1].first,log(V[N-1].second),V[N].first,log(V[N].second)));
                }
                return last.y = exp(spline(lambda,V[bound1-1].first,log(V[bound1-1].second),V[bound1].first,log(V[bound1].second),
                                           V[bound2].first,log(V[bound2].second),V[bound2+1].first,log(V[bound2+1].second)));
            } else {
                if (bound1==0) {
                    return last.y = exp(leftspline(log(lambda),log(V[0].first),log(V[0].second),log(V[1].first),log(V[1].second),log(V[2].first),log(V[2].second)));
                }
                if (bound2==V.size()-1) {
                    int N = V.size()-1;
                    return last.y = exp(rightspline(log(lambda),log(V[N-2].first),log(V[N-2].second),log(V[N-1].first),log(V[N-1].second),log(V[N].first),log(V[N].second)));
                }
                return last.y = exp(spline(log(lambda),log(V[bound1-1].first),log(V[bound1-1].second),log(V[bound1].first),log(V[bound1].second),
                                           log(V[bound2].first),log(V[bound2].second),log(V[bound2+1].first),log(V[bound2+1].second)));
            }
        }

        double x = lambda;
        double x1 = V[bound1].first;
        double x2 = V[bound2].first;
        double y1 = V[bound1].second;
        double y2 = V[bound2].second;

        if (interpolationy==LOG) {
            y1 = log(y1);
            y2 = log(y2);
        }
        if (interpolationx==LOG) {
            x1 = log(x1);
            x2 = log(x2);
            x = log(x);
        }

        last.y = (y2-y1)*(x-x1)/(x2-x1)+y1;

        if (interpolationy==LOG) {
            last.y = exp(last.y);
        }

        return last.y;
    }

    Table::
    Table(double *l,double *v,int nn)
    {
        values.resize(nn);
        for (int i=0; i<nn; ++i) {
            values[i]=PAIR(l[i],v[i]);
        }

        std::sort(values.begin(),values.end());

        last.x=l[0];
        last.y=v[0];
        name = "(internal array)";
        formulafile=NULL;
        tablefile=NULL;
        formulastring.clear();
        interpolationx=LIN;
        interpolationy=LIN;
    }

    Table::
    Table(const std::vector<double>& x,const std::vector<double>& y)
    {
        if (x.size()!=y.size()) throw SCATMECH_exception("Size mismatch in Table::Table()");

        values.resize(x.size());
        for (int i=0; i<(int)x.size(); ++i) {
            values[i]=PAIR(x[i],y[i]);
        }

        std::sort(values.begin(),values.end());

        last.x=x[0];
        last.y=y[0];
        name = "(internal array)";
        formulafile=NULL;
        tablefile=NULL;
        formulastring.clear();
        interpolationx=LIN;
        interpolationy=LIN;
    }

    Table::
    Table(const Table::VECTOR& v)
    {
        values = v;

        std::sort(values.begin(),values.end());

        last.x=v[0].first;
        last.y=v[0].second;
        name = "(internal array)";
        formulafile=NULL;
        tablefile=NULL;
        formulastring.clear();
        interpolationx=LIN;
        interpolationy=LIN;
    }

    Table::
    Table(const string& filename, int col)
    {
        formulafile=NULL;
        tablefile=NULL;
        formulastring.clear();
        set(filename,col);
        interpolationx=LIN;
        interpolationy=LIN;
    }

    Table::
    Table(const char* filename, int col)
    {
        formulafile=NULL;
        tablefile=NULL;
        formulastring.clear();
        set(filename,col);
        interpolationx=LIN;
        interpolationy=LIN;
    }

    bool
    Table::
    ReadFunctionFormat(const string& fname)
    {
        ifstream_with_comments file(fname.c_str());
        if (!file) throw SCATMECH_exception("Cannot open file " + fname);
        string word;
        file >> word;
        if (word.size()==0) throw SCATMECH_exception("File " + fname + " is empty");
        if (word=="PARAMETERS") {
            while (1) {
                file >> word;
                if (file.fail()) throw SCATMECH_exception("Failure to read file");
                if (word=="END") break;
                double value;
                file >> value;
                if (file.fail()) throw SCATMECH_exception("Expecting value in PARAMETERS section");

                Evaluator::VMAP::const_iterator it =  params.find(word);
                if (it==params.end()) {
                    params[word]=value;
                }
            }
            file >> word;
            if (file.fail()) throw SCATMECH_exception("Missing formula or FOR statement in " + fname);
            if (word=="FOR") return false;
            if (word[0]!='@') throw SCATMECH_exception("Missing formula or FOR statement in " + fname);
            formulastring=word;
            while (!file.eof()) {
                file >> word;
                if (!file.eof()) formulastring+=word;
            }
            return true;
        }
        return false;
    }

    int
    Table::
    set(const string& filename, int col)
    {
        formulastring.clear();
        last.x=numeric_limits<double>::max();

        // First, try to read the value as a single value
        try {
            Evaluator eval(filename,params);
            if (eval.NResult()>=col-1) {
                double result = eval.Result(col-2);
                last.y=result;
                *this=result;
                name = to_string(result);
                return 0;
            }
        } catch (exception&) {}

        vector<double> x,y;

        // Then see if it is a simple functional form: "@x:..."
        if (filename.size()>0) {
            if (filename[0]=='@') {
                formulastring = filename;
                icol = col;
                last.x = numeric_limits<double>::max();
                last.y = numeric_limits<double>::max();
                ostringstream str;
                str << filename << " (col=" << col << ")";
                name = str.str();
                return 0;
            }
        }

        string filename2 = filename;
        size_t found = filename.find('!');
        if (found!=string::npos) {
            if (filename.size()==found+1) throw SCATMECH_exception("No characters beyond ! in filename");
            filename2 = filename.substr(0,found);
            string postfix = filename.substr(found+1);
            istringstream postfixstrm(postfix);
            int newcol;
            postfixstrm >> newcol;
            if (!postfixstrm.fail()) {
                col = newcol+col-2;
            } else throw SCATMECH_exception("No column number following ! in filename");
        }

        maps_mutex.lock();
        TableFileMap::iterator t =  tablefiles.find(filename2);
        FormulaFileMap::iterator f =  formulafiles.find(filename2);
		TableFileMap::iterator tend = tablefiles.end();
		FormulaFileMap::iterator fend = formulafiles.end();
		maps_mutex.unlock();

        if (t!=tend) {
            x = (t->second)[0];
            y = (t->second)[col-1];
            icol=col;
            params.clear();
            tablefile=&(t->second);
        } else if (f!=fend) {
            x = f->second.get_column(1,params);
            y = f->second.get_column(col,params);
            params = f->second.get_parameter_map();
            formulafile = &(f->second);
            icol=col;
        } else {
            string fname = find_file(filename2);
            if (ReadFunctionFormat(fname)) {
                last.x = numeric_limits<double>::max();
                last.y = numeric_limits<double>::max();
                icol=col;
            } else {
                string word;
                {
                    ifstream_with_comments file(fname.c_str());
                    if (!file) throw SCATMECH_exception("Cannot open file " + filename);
                    file >> word;
                    if (word.size()==0) throw SCATMECH_exception("File " + filename + " is empty");
                }
                if (word=="PARAMETERS") {
			        maps_mutex.lock();
                    formulafiles[filename2]=FormulaFile(filename2);
                    f = formulafiles.find(filename2);
					maps_mutex.unlock();
                    x = f->second.get_column(1,params);
                    y = f->second.get_column(col,params);
                    params = f->second.get_parameter_map();
                    formulafile = &(f->second);
                    icol=col;
                } else {
					maps_mutex.lock();
					tablefiles[filename2]=TableFile(filename2);
                    t=tablefiles.find(filename2);
					maps_mutex.unlock();
                    x = (t->second)[0];
                    y = (t->second)[col-1];
                    params.clear();
                    icol=col;
                    tablefile = &(t->second);
                }
            }
        }

		if (x.size()>0) {
            values.resize(x.size());

            for (int i=0; i<(int)x.size(); ++i) {
                values[i]=PAIR(x[i],y[i]);
            }

            sort(values.begin(),values.end());

            last.x = x[0];
            last.y = y[0];
        }

        ostringstream str;
        str << filename2 << " (col=" << col << ")";
        name = str.str();

        return 0;
    }

    int
    Table::
    _AskUser(const string& query, const string& deflt)
    {
        string response = SCATMECH::AskUser(query + " (value or filename)",deflt);

        set(response);
        return 1;
    }

    Table
    Table::
    AskUser(const std::string& prompt,const std::string& deflt)
    {
        Table result;
        bool success;
        do {
            try {
                result._AskUser(prompt,deflt);
                success = true;
            } catch (SCATMECH_exception& e) {
                SCATMECH_output << e.what() << endl;
                success = false;
            }
        } while (!success);

        return result;
    }

    Table
    Table::
    AskUser(const std::string& prompt,double deflt)
    {
        ostringstream ww3;
        ww3 << deflt << ends;
        return AskUser(prompt,ww3.str());
    }

    void Table::set_parameter(const std::string& param,const std::string& value)
    {
        if (param=="") {
            set(value);
            return;
        }
        if (icompare(param,"interpolationx")) {
            if (value=="lin") interpolationx = LIN;
            else if (value=="log") interpolationx = LOG;
            else throw SCATMECH_exception("Invalid Table interpolation mode - must be lin or log");
            last.x = last.y = numeric_limits<double>::max();

            return;
        }
        if (icompare(param,"interpolationy")) {
            if (value=="lin") interpolationy = LIN;
            else if (value=="log") interpolationy = LOG;
            else if (value=="spline") interpolationy = SPLINE;
            else if (value=="logspline") interpolationy = LOGSPLINE;
            else throw SCATMECH_exception("Invalid Table interpolation mode - must be lin, log, spline, or logspline");
            last.x = last.y = numeric_limits<double>::max();

            return;
        }
        params[param]=from_string<double>(value);
		last.x=numeric_limits<double>::max();
    }

    string Table::get_parameter(const std::string& param) const
    {
        if (param=="") {
            ostringstream oss;
            oss << *this;
            return oss.str();
        } else if (icompare(param,"interpolationx")) {
            if (get_interpolationx() == Table::LIN) return string("lin");
            if (get_interpolationx() == Table::LOG) return string("log");
            return string("invalid");
        } else if (icompare(param,"interpolationy")) {
            if (get_interpolationy() == Table::LIN) return string("lin");
            if (get_interpolationy() == Table::LOG) return string("log");
            if (get_interpolationy() == Table::SPLINE) return string("spline");
            if (get_interpolationy() == Table::LOGSPLINE) return string("logspline");
            return string("invalid");
        } else {
            //if (formulafile!=NULL) {
            Evaluator::VMAP::const_iterator it = params.find(param);
            if (it!=params.end()) return to_string(it->second);
            //}
            throw SCATMECH_exception("No parameter " + param + " in Table");
        }
    }

    template <>
    void
    ModelParameterSet(Table& variable,const string& parameter,const string& value)
    {
        variable.set_parameter(parameter,value);
    }

    template <>
    string
    ModelParameterGet(Table& variable,const string& parameter)
    {
        return variable.get_parameter(parameter);
    }

    template <>
    void
    ModelParameterAskUser(Table& value,const string& prompt)
    {
        value = Table::AskUser(prompt,value.value());
    }

} // namespace SCATMECH


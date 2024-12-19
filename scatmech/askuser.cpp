//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: askuser.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <stack>
#include <fstream>
#include "scatmech.h"
#include "optconst.h"
#include "askuser.h"
#include "dielfunc.h"

using namespace std;


namespace SCATMECH {

    string
    Get_SCATMECH_Version()
    {
        return std::string("SCATMECH 7.23 (Build: ") + std::string(__DATE__) + std::string(")");
    }

    optical_constant default_substrate(4.15,0.05);
    optical_constant default_particle(1.59,0);
    double default_lambda=0.532;
    dielectric_function vacuum(optical_constant(1,0));

    streambuf_with_echo SCATMECH_echo(cin.rdbuf());
    istream_with_comments SCATMECH_input(&SCATMECH_echo);
    std::ostream SCATMECH_output(cerr.rdbuf());

    class SCATMECH_output_starter {
        public:
            SCATMECH_output_starter() {
                SCATMECH_output.setf(ios::unitbuf);
            }
    } _SCATMECH_output_starter;

    //char SCATMECH_comment=';';

    stack<ifstream*> instack;

    // Environment variables can be set...
    string SCATMECH_DIRECTORY; // setenv SCATMECH_DIRECTORY

    void
    check_SCATMECH_DIRECTORY()
    {
        static bool firsttime=true;
        if (firsttime) {

            #if _MSC_VER > 1000
            char *buffer=0;
            size_t size=0;
            SCATMECH_DIRECTORY=std::string();
            if (_dupenv_s(&buffer,&size,"SCATMECH")==0) {
                if (buffer!=NULL) {
                    SCATMECH_DIRECTORY=std::string(buffer);
                    free(buffer);
                }
            }
            #else
            const char* scatmech_directory = getenv("SCATMECH");
            if (scatmech_directory == NULL) {
                SCATMECH_DIRECTORY=std::string();
            } else {
                SCATMECH_DIRECTORY=std::string(scatmech_directory);
            }
            #endif
            firsttime=false;
        }
    }

    string
    find_file(const string& name)
    {
        check_SCATMECH_DIRECTORY();
        string pname = SCATMECH_DIRECTORY;

        ifstream file(name.c_str());

        #ifdef _MSC_VER
        const char delim = '\\';
        const char separator = ';';
        #else
        const char delim = '/';
        const char separator = ':';
        #endif

        if (!file) {
            unsigned int begin=0;
            while (begin<pname.size()) {
                pname = pname.substr(begin);
                int end = pname.find_first_of(separator);
                if (end==string::npos) end = pname.size();
                string directory = pname.substr(0,end);
                if (directory[directory.size()-1] != delim) directory += delim;
				string filename=directory+name;
				ifstream _file(filename.c_str());
                if (_file) return filename;
                begin = end+1;
            }
            throw SCATMECH_exception("Cannot open file \"" + name + "\"");
        } else {
            return name;
        }
    }

    void
    set_istream(istream &is)
    {
        SCATMECH_input = is;
    }

    void
    set_ostream(ostream &os)
    {
        SCATMECH_output.rdbuf(os.rdbuf());
    }

    void
    streambuf_with_comments::
    clear_comment()
    {
        int c;
        do {
            c = istr->sbumpc();
        } while (c!='\n');
    }

    int
    streambuf_with_comments::
    underflow()
    {
        int c = istr->sbumpc();
        if (c==comment_delimiter) {
            clear_comment();
            c='\n';
        }
        istr->sputbackc((char)c);
        return c;
    }

    int
    streambuf_with_comments::
    uflow()
    {
        int c = istr->sbumpc();
        if (c==comment_delimiter) {
            clear_comment();
            return '\n';
        }
        return c;
    }

    int
    streambuf_with_comments::
    pbackfail(char c)
    {
        return (c!=EOF) ? istr->sputbackc(c) : EOF;
    }

    streambuf_with_comments::
    streambuf_with_comments(std::streambuf *s)
    {
        istr = s;
        setp(0, 0);
        setg(0,0,0);
        comment_delimiter=';';
    }

    streambuf_with_comments&
    streambuf_with_comments::
    operator=(std::streambuf *s)
    {
        istr = s;
        setp(0, 0);
        setg(0,0,0);
        comment_delimiter=';';
        return *this;
    }

    streambuf_with_comments::
    ~streambuf_with_comments()
    {
        //delete istr;
    }

    streambuf_with_echo::
    streambuf_with_echo(std::streambuf *s,const std::string& name)
    {
        strbuf = s;
        echo = false;
        set_echo(name);
    }

    streambuf_with_echo&
    streambuf_with_echo::
    operator=(std::streambuf *s)
    {
        strbuf = s;
        if (echo) echostream.close();
        echo = false;
        return *this;
    }

    streambuf_with_echo::~streambuf_with_echo()
    {}

    void
    streambuf_with_echo::
    set_echo(const std::string& name)
    {
        if (echo) echostream.close();
        if (name.size()!=0) {
            echostream.open(name.c_str());
            if (!echostream) throw SCATMECH_exception("Cannot open file for streambuf_with_echo: " + name);
            echo = true;
        } else {
            echo = false;
        }
    }

    int
    streambuf_with_echo::
    underflow()
    {
        int c = strbuf->sbumpc();
        if (echo) {
            echostream.put((char)c);
            echostream.flush();
        }
        strbuf->sputbackc((char)c);
        return c;
    }

    int
    streambuf_with_echo::
    uflow()
    {
        int c = strbuf->sbumpc();
        if (echo) {
            echostream.put((char)c);
            echostream.flush();
        }
        return c;
    }

    int
    streambuf_with_echo::
    pbackfail(char c)
    {
        return (c!=EOF) ? strbuf->sputbackc(c) : EOF;
    }

    void
    Set_SCATMECH_Input_Echo(const std::string& name)
    {
        SCATMECH_echo.set_echo(name);
    }

    void
    SCATMECH_Flush_Echo()
    {
        SCATMECH_echo.flush_echo();
    }

    string
    getstrline(istream& str)
    {
        int c;
        string result;
        while ((c=str.get())!='\n') {
            if (str.eof()) {
                //if (result.size()>0) str.clear();
                return result;
            }
            result += (char)c;
        }
        //if (result.size()>0) str.clear();
        return result;
    }

    string
    istream_with_comments::
    getstrline()
    {
        return SCATMECH::getstrline(*this);
    }

    string
    istream_with_comments::
    getquoted()
    {
        const char quote='\"';
        string result;
        ws(*this);
        if (peek()==quote) {
            get();
            int c;
            while ((c=get())!=quote&&!fail()) result += (char)c;
            if (fail()) throw SCATMECH_exception("Failed to find closing quotation mark");
        } else *this >> result;
        return result;
    }


} // namespace SCATMECH




//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: askuser.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_ASKUSER_H
#define SCATMECH_ASKUSER_H

#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <stack>
#include <list>
#include <cctype>     // std::tolower
#include <cstddef>    // std::size_t
#include <algorithm>

namespace SCATMECH {

    extern std::string SCATMECH_DIRECTORY;
    void check_SCATMECH_DIRECTORY();
    std::string find_file(const std::string& name);

    /// @brief An input streambuf with some added features
    ///
    /// The streambuf_with_comments class is a stream buffer that removes comments
    /// (remainder of any line starting with a semicolon) and allows
    /// for redirection of the input (any <filename).
    ///
    class streambuf_with_comments : public std::streambuf {
        private:
            std::streambuf* istr;
            char comment_delimiter;

        protected:
            void clear_comment();
            int underflow();
            int uflow();
            int pbackfail(char c);

            std::streampos seekpos(std::streampos sp, std::ios_base::openmode which) {
                return istr->pubseekpos(sp,which);
            }
            std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,std::ios_base::openmode which) {
                return istr->pubseekoff(off,way,which);
            }
            int sync() {
                return istr->pubsync();
            }

        public:
            streambuf_with_comments(std::streambuf *s);
            streambuf_with_comments& operator=(std::streambuf *s);
            virtual ~streambuf_with_comments();
            void set_comment_delimiter(char c) {
                comment_delimiter=c;
            }
    };

    /// @brief An input streambuf with some added features
    ///
    /// The streambuf_with_echo class is a stream buffer that, for input,
    /// echos the stream to a file.
    ///
    class streambuf_with_echo : public std::streambuf {
        private:
            std::streambuf* strbuf;
            std::ofstream echostream;
            bool echo;

        protected:
            int underflow();
            int uflow();
            int pbackfail(char c);

        public:
            streambuf_with_echo(std::streambuf *s,const std::string& name="");
            streambuf_with_echo& operator=(std::streambuf *s);
            virtual ~streambuf_with_echo();
            void set_echo(const std::string& name);
            void flush_echo() {
                echostream.flush();
            }
    };

    /// @brief An input stream with some added features
    ///
    /// The istream_with_comments uses streambuf_with_comments to ignore
    /// comments and allow for redirection.
    ///
    class istream_with_comments : public std::istream {
        public:
            typedef std::istream istream;
            istream_with_comments(std::streambuf* _buf) : buf(_buf), istream(&buf) {}

            istream_with_comments& operator=(std::istream& _str) {
                buf = _str.rdbuf();
                return *this;
            }

            /// Get an entire line and return as a string
            std::string getstrline();

            /// Get the next token, which may be surrounded by double quotes
            std::string getquoted();

        protected:
            streambuf_with_comments buf;
    };

    /// @brief An input file stream with some added features
    ///
    /// The ifstream_with_comments uses streambuf_with_comments to ignore
    /// comments and allow for redirection.
    ///
    class ifstream_with_comments : public istream_with_comments {
        public:
            ifstream_with_comments(const char* name) : fbuf(), istream_with_comments(&fbuf) {
                using namespace std;
                std::string newname = find_file(name);
                open(newname.c_str());
            }
            ifstream_with_comments() : istream_with_comments(&fbuf) {}
            void open(const char* name) {
                using namespace std;
                fbuf.open(name,ios::in);
                if (!fbuf.is_open()) setstate(failbit);
            }
            void close() {
                fbuf.close();
            }
        protected:
            std::filebuf fbuf;
    };

    /// The input stream used by the SCATMECH library
    extern istream_with_comments SCATMECH_input;

    /// The output stream used by the SCATMECH library
    extern std::ostream SCATMECH_output;

    /// @brief Template function to query user for a value
    ///
    /// Queries the user for a value.  It sends the prompt to the SCATMECH output
    /// stream, followed by the default value in angle brackets <> and waits
    /// for user input.  If something throws an exception, it asks again.
    template <class TYPE>
    TYPE AskUser(
        const std::string& prompt,  ///< The prompt
        const TYPE& deflt
    )
    {
        while (1) {
            try {
                SCATMECH_output << prompt << " <" << deflt << "> :";

                std::istringstream line(SCATMECH_input.getstrline());

                TYPE result;
                line >> result;

                if (line.fail()) result = deflt;

                return result;
            }
            catch (std::exception& e) {
                SCATMECH_output << e.what() << std::endl;
                if (SCATMECH_input.eof()) throw e;
            }
        }
    }

    /// Function to query the user for a string
    inline std::string AskUser(const std::string& query, const char* deflt)
    {
        return AskUser(query,std::string(deflt));
    }

    /// Set the standard SCATMECH input stream
    void set_istream(std::istream &is);

    /// Set the standard SCATMECH output stream
    void set_ostream(std::ostream &os);

    /// @brief Template function for C style formatting
    ///
    /// This function is useful for output streams when one wants
    /// to use tradtional C-style printf formats.  For example:
    ///      cout << format("%16.12E",x) << endl;
    template <class TYPE>
    std::string format(const std::string& fmt, const TYPE& value) {
        char buffer[256];
        #ifdef _MSC_VER
        sprintf_s(buffer,256,fmt.c_str(),value);
        #else
        snprintf(buffer,256,fmt.c_str(),value);
        #endif
        return std::string(buffer);
    }

    /// A string version of getline()...
    std::string getstrline(std::istream& str);

    /// Get a token from an input stream
    //std::string get_token(std::istream &is,const char* delim);

	class CommandLineParser {
	private:
		typedef std::list<std::string> TokenList;
	public:
		CommandLineParser(const int &argc, const char **argv) {
			for (int i = 1; i < argc; ++i) 
				tokens.push_back(argv[i]);
		}
		std::string getOption(const std::string &option) const {
			TokenList::const_iterator itr;
			itr = std::find(tokens.begin(), tokens.end(), option);
			if (itr != tokens.end() && ++itr != tokens.end()) {
				return *itr;
			}
			else {
				return std::string("");
			}
		}
		bool optionExists(const std::string &option) const {
			return (std::find(tokens.begin(), tokens.end(), option)!= tokens.end());
		}
	private:
		TokenList tokens;
	};



} // namespace SCATMECH


#endif

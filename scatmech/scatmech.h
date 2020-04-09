//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: scatmech.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_SCATMECH_H
#define SCATMECH_SCATMECH_H

#include <exception>
#include <complex>
#include <string>
#include <list>
#include <iomanip>
#include "askuser.h"

#ifdef USE_PTHREADS
#include <pthread.h>
#endif


namespace SCATMECH {

    typedef std::string STRING;
    typedef std::complex<double> COMPLEX;

    //
    // Some useful values...
    //
    const double pi = 4.*atan(1.);
    const double deg = pi/180.;

    const char tab   = '\t';
    const char space = ' ';
    const char ret   = '\r';

    //
    // Some useful inline functions...
    //
    template <class T>
    inline T sqr(const T& x) {
        return x*x;
    }
    template <class T>
    inline T cube(const T& x) {
        return x*x*x;
    }

    //
    // The following are useful for error handling...
    //

    class SCATMECH_exception : public std::exception
    {
        public:

            SCATMECH_exception(const std::string m)
            {
                message = "SCATMECH: ";
                message += std::string(m);
            }

            ~SCATMECH_exception() throw() {}

            virtual const char *what() const throw()
            {
                return message.c_str();
            }

        private:
            std::string message;
    };

    template <class T>
    std::string to_string(const T& t,int precision=-1)
    {
        std::ostringstream oss;
        if (precision<0) oss << std::setprecision(precision);
        oss << t;
        return oss.str();
    }

    template <class T>
    T from_string(const std::string& s)
    {
        T result;
        std::istringstream iss(s);
        iss >> result;
        //if (iss.fail()) throw SCATMECH_exception("Unrecognized value in from_string<T>()");
        return result;
    }

    template <>
    inline std::string from_string<std::string>(const std::string& s)
    {
        return s;
    }

    std::string Get_SCATMECH_Version();

    void Set_SCATMECH_Input_Echo(const std::string& name);
    void SCATMECH_Flush_Echo();

    #ifdef USE_PTHREAD
    class mutex_t {
        public:
            mutex_t() {
                pthread_mutex_init(&it,NULL);
            }
            void lock() {
                pthread_mutex_lock(&it);
            }
            void unlock() {
                pthread_mutex_unlock(&it);
            }
        private:
            pthread_mutex_t it;
    };

    #else
    class mutex_t {
        public:
            void lock() {}
            void unlock() {}
    };
    #endif

    ///
    /// PointerCollector<T> stores pointers to T and deletes them
    /// when it has been decided they are no longer needed, or
    /// at the end of the program.
    ///
    template <class T>
    class PointerCollector {
        private:
            typedef std::list<T*> PointerList;
            typedef typename PointerList::iterator PointerListIter;
            PointerList trash;

        public:
            ~PointerCollector() {
                empty();
            }

            T* New() {
                T* t = new T;
                trash.push_back(t);
                return t;
            }

            T* New(T* t) {
                trash.push_back(t);
                return t;
            }

            T* New(const T& tt) {
                T* t = new T(tt);
                trash.push_back(t);
                return t;
            }

            void empty() {
                for (PointerListIter q = trash.begin(); q!=trash.end(); ++q) {
                    delete (*q);
                }
                trash.clear();
            }
    };

    //
    // This is a simple class for a loop int that iterates from one value to another
    // in a specified direction.
    //
    class forloopint {
        public:
            /// iterator class for looper, which looks like an integer...
            class iterator {
                public:
                    iterator(int i,bool _downward) : here(i), downward(_downward) {}
                    int operator ++() {
                        return (downward ? --here : ++here);
                    }
                    bool operator!=(iterator& i) {
                        return (i.here!=here || i.downward!=downward);
                    }
                    operator const int&() {
                        return here;
                    }
                private:
                    int here; /// It's current value
                    bool downward; /// Whether it iterates downward (true) or upward (false)
            };

            forloopint(int _bottom,int _top,int _direction = +1) {
                downward = _direction<0 ? true : false;
                if (downward) {
                    first = _top;
                    last = _bottom;
                } else {
                    first = _bottom;
                    last = _top;
                }
            }
            iterator begin() {
                if (downward && first<last) return end();
                if (!downward && first>last) return end();
                return iterator(first,downward);
            }
            iterator end() {
                return (downward ? iterator(last-1,downward) : iterator(last+1,downward));
            }

        private:
            int first;
            int last;
            bool downward;
    };


} // namespace SCATMECH


#endif

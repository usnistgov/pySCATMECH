//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: random.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef RANDOM_H
#define RANDOM_H

#include <ctime>
#include <vector>
#include "scattabl.h"
#include <cmath>

namespace SCATMECH {

    namespace CMLIB {
        ///
        /// class UNIrand gives a uniform random number on the interval [0,1)
        /// adapted from the CMLIB library...
        ///
        class UNIrand {
            public:
                UNIrand();
                void seed(int JD);
                double operator()();

            private:
                static int seed_add;
                static int seed_adder;
                int M[17],I,J,M1,M2;
                bool seeded;
        };
    }
    ///
    /// class Random_Number gives a random number from an arbitrary distribution.
    ///
    class Random_Number {
        public:
            ///
            /// The default distribution is a uniform distribution on [0,1);
            ///
            Random_Number();

            ///
            /// Other distributions can be set with the following constructor...
            ///
            Random_Number(const SCATMECH::Table& distrib);

            ///
            /// The following are some 'canned' distributions...
            ///

            /// Uniform distribution on [xmin,xmax)...
            static Table uniform(double xmin=0.,double xmax=1.);
            /// Normal distribution with standard deviation stdev.  The number of
            /// points in the curve is npoints, and the curve extends out nstdev*stdev.
            static Table normal(double stdev = 1., int npoints=400, double nstdev=5.);
            /// Sin-squared distribution on [0,2*pi).  The number of points in the
            /// distribution is set by npoints.
            static Table sinesqr(int npoints=400);
            /// Exponential distribution with 1/e constant decay.  The number of points
            /// in the distribution is npoints, and it extends ndecay*decay.
            static Table exponential(double decay, int npoints=400,double ndecay=10.);
            /// Log-Normal distribution with mode mode, and standard deviation stdev.  The
            /// number of points in the distribution is npoints, and it extends nstdev*stdev/2, logarithmically
            /// in each direction from the mode.
            static Table log_normal(double mode, double stdev = 1., int npoints=400, double nstdev=5.);

            ///
            /// The following returns a random number with the given distribution...
            ///
            double operator()();

            ///
            /// The following returns a random number with the given distribution, with
            /// a uniform background of offset.
            ///
            double operator()(double offset);

            ///
            /// The following sets the seed for the random number generator...
            ///
            void set_seed(long seed);

        private:

            void set_distribution(const SCATMECH::Table& pdf);
            void initialize();
            SCATMECH::Table lookup;
            CMLIB::UNIrand rand;
            static int seed_add;
            static int seed_adder;
    };

    ///
    /// Multivariate_Random_Number returns a set of N multivariate normally distributed random numbers
    /// having a specific covariance matrix.
    ///
    class Multivariate_Random_Number {
        public:
            ///
            /// The default constructor initializes the class with N=0.
            ///
            Multivariate_Random_Number();
            ///
            /// The following constructor initializes the class for n random numbers,
            /// with nxn covariance matrix covmatrix;
            ///
            Multivariate_Random_Number(int n, const std::vector<double>& covmatrix); // elements are covmatrix[i+n*j]

            ///
            /// The following returns (in the vector values) the set of N random numbers.
            ///
            void get(std::vector<double>& values);

            ///
            /// The following sets the seed for the random number generator...
            ///
            void set_seed(long seed) {
                random.set_seed(seed);
            }

        private:
            Random_Number random;
            std::vector<COMPLEX> eigenvalues,eigenvectors;
            std::vector<double> temp;
            int N;
    };

    ///
    /// The following returns a SCATMECH::Table of values for a self-affine random curve in the
    /// interval [0,period).  The correlation function is exp(-pow(sqr(z/corrlength),exponent))
    ///
    Table Random_Curve(double rms,        //< The root-mean-square amplitude of the fluctuations
                       double corrlength, //< The correlation length
                       double exponent,   //< The Hurst exponent
                       double period,     //< The length over which the function is defined (and periodic)
                       int n,             //< The number of points in the curve
                       int seed=0         //< The random number seed.
                      );
}

#endif

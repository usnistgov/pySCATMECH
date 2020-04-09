//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: oasphere.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_OASPHERE_H
#define SCATMECH_OASPHERE_H

#include "sphrscat.h"
#include "dielfunc.h"
#include <vector>

namespace SCATMECH {

    ///
    /// Optically_Active_Sphere_Scatterer is a Free_Space_Scatterer which
    /// returns the solution for the free space scattering of a homogeneous
    /// sphere consisting of an optically active or circularly dichroic medium.
    /// The calculations follow the description given in Section 8.3 of
    /// Bohren and Huffman, Absorption and Scattering of Light by Small
    /// Particles (Wiley, New York, 1983).
    ///
    class Optically_Active_Sphere_Scatterer: public Free_Space_Scatterer {
        public:
            Optically_Active_Sphere_Scatterer();

            virtual JonesMatrix jones(const Vector& kin,const Vector& kout);

            double CscaL(); /// Scattering cross section for left-circular polarization (LCP)
            double CextL(); /// Extinction cross section for LCP
            double CbackL();/// Backscattering cross section for LCP
            double CabsL()  /// Absorption cross section for LCP
            {
                return CextL()-CscaL();
            }
            double CscaR(); /// Scattering cross section for right-circular polarization (RCP)
            double CextR(); /// Extinction cross section for RCP
            double CbackR();/// Backscattering cross section for RCP
            double CabsR()  /// Absorption cross section for RCP
            {
                return CextR()-CscaR();
            }

            /// Scattering efficiency for LCP
            double QscaL() {
                return CscaL()/(pi*sqr(radius));
            }
            /// Extinction efficiency for LCP
            double QextL() {
                return CextL()/(pi*sqr(radius));
            }
            /// Backscattering efficiency for LCP
            double QbackL() {
                return CbackL()/(pi*sqr(radius));
            }
            /// Absorption efficiency for LCP
            double QabsL() {
                return CabsL()/(pi*sqr(radius));
            }

            /// Scattering efficiency for RCP
            double QscaR() {
                return CscaR()/(pi*sqr(radius));
            }
            /// Extinction efficiency for RCP
            double QextR() {
                return CextL()/(pi*sqr(radius));
            }
            /// Backscattering efficiency for RCP
            double QbackR() {
                return CbackR()/(pi*sqr(radius));
            }
            /// Absorption efficiency for LCP
            double QabsR() {
                return CabsR()/(pi*sqr(radius));
            }

        public:

            DECLARE_MODEL();
            /// Radius of the sphere
            DECLARE_PARAMETER(double,radius);
            /// Optical properties for left-circularly polarized light
            DECLARE_PARAMETER(dielectric_function,materialleft);
            /// Optical properties for right-circularly polarized light
            DECLARE_PARAMETER(dielectric_function,materialright);

        protected:

            std::vector<COMPLEX> a,b,c;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;
            double k;
            double x;
            void setup();
    };

} // namespace SCATMECH


#endif

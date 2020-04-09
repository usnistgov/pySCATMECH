//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sphprt.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_DIPARTICLE_H
#define SCATMECH_DIPARTICLE_H

#include "scatmech.h"
#include "filmtran.h"
#include "sphrscat.h"
#include "local.h"
#include "matrix3d.h"


namespace SCATMECH {



    //
    // Double_Interaction_BRDF_Model is a base class inheriting BRDF_Model and
    // dielectric_stack which solves for the scattering from a defect above a
    // dielectric stack (or none at all).
    //
    // The dielectric_stack properties of the class was changed from member
    // status to inherited status by TAG 31.11.2000
    // The class Spherical_Particle_BRDF_Model was changed to
    // Double_Interaction_BRDF_Model by TAG 22.12.2002

    class Double_Interaction_BRDF_Model :
        public Local_BRDF_Model
    {
        public:

            // The following line was removed by TAG 31.11.2000:
            // void set_stack_file() {stack.read_stack_file(lambda); recalc=1;}
            // The following line was added by TAG 31.11.2000:
            // void set_stack_file() {read_stack_file(); recalc=1;}

        public:
            DECLARE_MODEL();

            // The distance of the center of the sphere to the surface in
            // units of length...
            DECLARE_PARAMETER(double,distance);

            // The particular scatterer's class...
            DECLARE_PARAMETER(Free_Space_Scatterer_Ptr,scatterer);

            // Any dielectric layers on the surface...
            DECLARE_PARAMETER(StackModel_Ptr,stack);

            // The parameters for this virtual model are the Euler angles for the rotation of the
            // scatterer...
            //    First: Rotate alpha about the z axis
            //    Second: Rotate beta about the x axis
            //    Third: Rotate gamma about the z axis.
            //
            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(double,beta);

        protected:
            // The dielectric stack below the particle...
            // The following line was removed by TAG 31.11.2000:
            // dielectric_stack stack;


            // Routine to carry out once-only routines...
            void setup();

            // Jones matrix for scattering...
            virtual JonesMatrix jonesDSC();

            // Euler rotation matrix for scatterer
            Matrix Euler;
    };


} // namespace SCATMECH


#endif

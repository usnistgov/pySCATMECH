//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: brdf.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

#ifndef SCATMECH_BRDF_H
#define SCATMECH_BRDF_H

#include <limits>
#include "scatmech.h"
#include "mueller.h"
#include "dielfunc.h"
#include "inherit.h"
#include "vector3d.h"

namespace SCATMECH {

    ///
    /// BRDF_Model is the virtual base class for all surface scattering models...
    ///
    class BRDF_Model : public Model {
        public:

            /// Coordinate_System is used for defining how the basis set for
            /// polarization is defined...
            enum Coordinate_System {
                psps,     ///< p (TM) and s (TE)
                xyxy,     ///< x and y (rotated from p and s by phis, to remove singularity at surface normal)
                plane,    ///< parallel and perpendicular to the scattering plane
                undefined ///< Don't use this
            };

            /// Type is used for defining whether the calculation
            /// is done in reflection or transmission and forward versus backward
            enum Type {
                //REFLECTION = 0, // removed
                //TRANSMISSION = 1, // removed
                Type_DOWNUP = 0,    ///< Incident down, scattering up (reflection)
                Type_DOWNDOWN = 1,  ///< Incident down, scattering down (transmission)
                Type_UPDOWN = 2,    ///< Incident up, scattering down (reflection from substrate)
                Type_UPUP = 3       ///< Incident up, scattering up (transmission from substrate)
            };


            /// The constructor for BRDF_Model...
            BRDF_Model() {
                model_cs = psps; /// The default polarization coordinate system is psps.
                last_thetai = last_thetas = last_phis = last_rotation = std::numeric_limits<double>::infinity();
                recalc_on_incident_change = recalc_on_scatter_change = 0x00;
            }

            ///
            /// Public interface for the Jones matrix for scattering
            ///
            /// When typecast to MuellerMueller matrix, Jones returns the Mueller matrix BRDF
            ///
            JonesMatrix Jones(double thetai,				///< Incident polar angle
                              double thetas,				///< Scattering polar angle
                              double phis,					///< Scattering azimuthal, out-of-plane, angle
                              double rotation,				///< Rotation of the sample
                              Coordinate_System cs = psps); ///< Basis set for polarization
            ///
            /// Public interface for the Mueller matrix BRDF
            ///
            MuellerMatrix Mueller(double thetai,			    ///< Incident polar angle
                                  double thetas,			    ///< Scattering polar angle
                                  double phis,				    ///< Scattering azimuthal, out-of-plane, angle
                                  double rotation,			    ///< Rotation of the sample
                                  Coordinate_System cs = psps); ///< Basis set for polarization


            /// Jones matrix for scattering where the directions to the source, the viewer,
            /// the surface normal, and the x-axis on the sample are specified by Vectors...
            JonesMatrix Jones(const Vector& source,		            ///< Vector pointing from the sample to the source
                              const Vector& viewer,		            ///< Vector pointing from the sample to the viewer
                              const Vector& normal,		            ///< Vector pointing away from the surface
                              const Vector& xaxis=Vector(0,0,0),	///< Vector pointing along the sample's fiducial x-axis
                              Coordinate_System cs=plane);			///< Basis set for polarization

            /// Mueller matrix for scattering where the directions to the source, the viewer,
            /// the surface normal, and the x-axis on the sample are specified by Vectors...
            MuellerMatrix Mueller(const Vector& source,		            ///< Vector pointing from the sample to the source
                                  const Vector& viewer,		            ///< Vector pointing from the sample to the viewer
                                  const Vector& normal,		            ///< Vector pointing away from the surface
                                  const Vector& xaxis=Vector(0,0,0),	///< Vector pointing along the sample's fiducial x-axis
                                  Coordinate_System cs=plane);			///< Basis set for polarization


            /// The following returns the coordinate system used by the model...
            Coordinate_System get_model_cs() const {
                return model_cs;
            }


            DECLARE_MODEL();

            ///
            /// The wavelength of the light in vacuum [&mu;m]
            ///
            DECLARE_PARAMETER(double,lambda);

            ///
            /// The optical constants of the substrate, expressed as a complex number (n,k) or, optionally, as a function of wavelength.
            ///
            DECLARE_PARAMETER(dielectric_function,substrate);

            ///  	Indicates whether the light is incident from above the substrate or from within the substrate
            ///	    and whether the scattering is evaluated in reflection or transmission. The choices are: <br>
            /// 0 : Light is incident from the above the substrate, and scattering is evaluated in reflection.<br>
            /// 1 : Light is incident from the above the substrate, and scattering is evaluated in transmission.<br>
            /// 2 : Light is incident from the within the substrate, and scattering is evaluated in reflection.<br>
            /// 3 : Light is incident from the within the substrate, and scattering is evaluated in transmission.<br>
            ///     For 1, 2, and 3, the substrate must be non-absorbing.<br>
            ///     Not all modes are supported for all models inheriting BRDF_Model.
            DECLARE_PARAMETER(int,type);

        protected:

            /// Routine which does housekeeping if (recalc!=0) ...
            virtual void setup();

            /// Coordinate system for which the model is defined
            /// (default: BRDF_Model::psps)
            /// Note: there is no set_model_cs() function since it is up to the
            /// model to define this, if it is not psps.
            Coordinate_System model_cs;

            /// Set the scattering geometry by angles
            /// The returned value = 1 if the geometry is valid, 0 if not.
            int set_geometry(double _thetai,double _thetas,double _phis, double _rotation);
            /// Set the scattering geometry by vectors
            /// The returned value = 1 if the geometry is valid, 0 if not.
            int set_geometry(const Vector& source,const Vector& viewer,const Vector& normal,const Vector& xaxis);

            /// Convert the Jones matrix from that for the model to another basis.
            void convert(JonesMatrix& J,Coordinate_System to) const;
            /// Convert the Mueller matrix from that for the model to another basis.
            void convert(MuellerMatrix& M,Coordinate_System to) const;

            /// The incident angle (in radians). Normally set by set_geometry()
            double thetai;
            /// The scattering polar angle (in radians). Normally set by set_geometry()
            double thetas;
            /// The scattering azimuthal angle (in radians). Normally set by set_geometry()
            double phis;
            /// The rotation of the sample (in radians). Normally set by set_geometry()
            double rotation;

            ///
            /// All BRDF_Model must define ONE of the virtual functions jones() or mueller().
            ///

            /// The Jones matrix for scattering. This function, when typecast to MuellerMatrix must return the
            /// Mueller matrix BRDF.
            virtual JonesMatrix jones() {
                //  This default function should never be called.
                error("Attempt to call jones() for depolarizing model.");
                return JonesMatrix();
            }

            ///
            /// The Matrix matrix BRDF for scattering
            ///
            virtual MuellerMatrix mueller()
            {
                // The default is to typecast from jones()
                return this->jones();
            }

            ///
            /// This routine should be called by the class constructor if setup() should be called when
            /// either the incident direction changes or the scattering direction changes.
            /// The
            ///
            void set_recalc_on_direction_change(int incident_change,int scatter_change) {
                recalc_on_incident_change = incident_change;
                recalc_on_scatter_change = scatter_change;
            };

        public:
            /// Test if type is reflection
            bool is_reflection() const {
                return type==Type_DOWNUP||type==Type_UPDOWN;
            }
            /// Test if type is transmission
            bool is_transmission() const  {
                return type==Type_DOWNDOWN||type==Type_UPUP;
            }
            /// Test if incident light is from vacuum
            bool is_forward() const {
                return type==Type_DOWNDOWN||type==Type_DOWNUP;
            }
            /// Test if incident light is from substrate
            bool is_backward() const  {
                return type==Type_UPDOWN||type==Type_UPUP;
            }

            /// Test if normal reflection
            bool is_down_to_up() const  {
                return type==Type_DOWNUP;
            }
            /// Test if normal transmission
            bool is_down_to_down() const  {
                return type==Type_DOWNDOWN;
            }
            /// Test if reverse reflection
            bool is_up_to_down() const  {
                return type==Type_UPDOWN;
            }
            /// Test if reverse transmission
            bool is_up_to_up() const  {
                return type==Type_UPUP;
            }

            /// Test if the incident light is upward propagating
            bool is_inup() const  {
                return is_backward();
            }
            /// Test if the incident light is downward propagating
            bool is_indown() const  {
                return is_forward();
            }
            /// Test if the scattered light is upward propagating
            bool is_outup() const  {
                return type==Type_DOWNUP||type==Type_UPUP;
            }
            /// Test if the scattered light is downward propagating
            bool is_outdown() {
                return type==Type_DOWNDOWN||type==Type_UPDOWN;
            }

        protected:
            /// Throw an exception if the incident light is downward propagating
            void throw_forward() {
                if (is_forward()) error("No forward (type==0 or 1) code for model");
            }
            /// Throw an exception if the incident light is upward propagating
            void throw_backward() {
                if (is_backward()) error("No backward (type==2 or 3) code for model");
            }
            /// Throw an exception if the scattering is in transmission
            void throw_transmission() {
                if (is_transmission()) error("No transmission (type==1 or 3) code for model");
            }
            /// Throw an exception if the scattering is in reflection
            void throw_reflection() {
                if (is_reflection()) error("No reflection (type==0 or 2) code for model");
            }

        private:
            double last_thetai;
            double last_thetas;
            double last_phis;
            double last_rotation;
            int recalc_on_incident_change;
            int recalc_on_scatter_change;
    };

    ///
    /// A smart pointer to a BRDF_Model
    ///
    typedef Model_Ptr<BRDF_Model> BRDF_Model_Ptr;

    void Register(const BRDF_Model* x);


} // namespace SCATMECH


#endif


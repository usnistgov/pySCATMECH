//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crossgrating2.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef CROSSGRATING2_H
#define CROSSGRATING2_H

#include "crossgrating.h"

namespace SCATMECH {

    class Cylinder_CrossGrating : public Gridded_CrossGrating {
        protected:

            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,rtop);
            DECLARE_PARAMETER(Table,rbottom);
            DECLARE_PARAMETER(double,thickness);
            DECLARE_PARAMETER(int,nlevels);
            DECLARE_PARAMETER(dielectric_function,inside);
            DECLARE_PARAMETER(dielectric_function,outside);
    };

	class Rectangle_CrossGrating : public Gridded_CrossGrating {
	protected:

		void setup();

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, length1);
		DECLARE_PARAMETER(double, length2);
		DECLARE_PARAMETER(double, zetaa);
		DECLARE_PARAMETER(double, thickness);
		DECLARE_PARAMETER(dielectric_function, inside);
		DECLARE_PARAMETER(dielectric_function, outside);
	};

	class OneD_CrossGrating : public CrossGrating {
        protected:

            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(Grating_Ptr,grating);
            DECLARE_PARAMETER(double,d2);
            DECLARE_PARAMETER(double,zeta);
    };

    class Overlaid_CrossGrating : public CrossGrating {
        protected:
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(CrossGrating_Ptr,top);
            DECLARE_PARAMETER(CrossGrating_Ptr,bottom);
            DECLARE_PARAMETER(double,overlay1);
            DECLARE_PARAMETER(double,overlay2);
            DECLARE_PARAMETER(double,separation);
    };

    class Overlaid_1D_CrossGrating : public CrossGrating {
        protected:
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(Grating_Ptr,top);
            DECLARE_PARAMETER(Grating_Ptr,bottom);
            DECLARE_PARAMETER(double,angle);
            DECLARE_PARAMETER(double,separation);
    };

    class Null_CrossGrating : public CrossGrating {
        protected:
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,d1);
            DECLARE_PARAMETER(double,d2);
            DECLARE_PARAMETER(double,zeta);
    };

    class Sphere_CrossGrating : public Gridded_CrossGrating {
        protected:

            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,diameter);
            DECLARE_PARAMETER(double,above);
            DECLARE_PARAMETER(double,below);
            DECLARE_PARAMETER(int,nlevels);
            DECLARE_PARAMETER(dielectric_function,sphere);
            DECLARE_PARAMETER(dielectric_function,surrounding);
    };

    class Pyramidal_Pit_CrossGrating : public Gridded_CrossGrating {
        protected:
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,side);
            DECLARE_PARAMETER(double,depth);
            DECLARE_PARAMETER(int,nlevels);
    };

}
#endif

//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crossgrating2.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "crossgrating2.h"

using namespace std;


namespace SCATMECH {

    void Cylinder_CrossGrating::setup()
    {
        Gridded_CrossGrating::setup();

        levels = nlevels;
        eps.allocate(grid1,grid2,levels);
        thick.allocate(levels);

        COMPLEX einside = epsilon(inside);
        COMPLEX eoutside = epsilon(outside);

        for (int n=1; n<=levels; ++n) {
            thick(n) = thickness/levels;
        }
        for (int i=1; i<=grid1; ++i) {
            for (int j=1; j<=grid2; ++j) {
                double x,y;
                getxy(i,j,x,y);

                double theta=atan2(y,x);
                double r = sqrt(x*x+y*y);

                double _rtop = rtop.value(theta/deg);
                double _rbottom = rbottom.value(theta/deg);

                for (int n=1; n<=levels; ++n) {
                    double radius = (n-0.5)/levels*(_rtop-_rbottom)+_rbottom;

                    if (r<=radius) {
                        eps(i,j,n) = einside;
                    } else {
                        eps(i,j,n) = eoutside;
                    }
                }
            }
        }

        FourierFactorize();
    }

	void Rectangle_CrossGrating::setup()
	{
		Gridded_CrossGrating::setup();

		levels = 1;
		eps.allocate(grid1, grid2, levels);
		thick.allocate(levels);

		COMPLEX einside = epsilon(inside);
		COMPLEX eoutside = epsilon(outside);

		thick(1) = thickness;

		double r1 = length1 / 2;
		double r2 = length2 / 2;

		double tanzetaa = tan(zetaa*deg);
		double seczetaa = 1. / cos(zetaa*deg);

		for (int i = 1; i <= grid1; ++i) {
			for (int j = 1; j <= grid2; ++j) {
				double x, y;
				getxy(i, j, x, y);

				double x1a = x - y*tanzetaa;
				double x2a = y*seczetaa;

				if (x1a<r1 && x1a>-r1 && x2a<r2 && x2a>-r2) {
					eps(i, j, 1) = einside;
				}
				else
				{
					eps(i, j, 1) = eoutside;
				}			
			}
		}

		FourierFactorize();
	}

	void Sphere_CrossGrating::setup()
    {
        Gridded_CrossGrating::setup();

        if (above<0) error("abovelevels<0");
        if (below<0) error("belowlevels<0");

        int abovelevels = (above!=0. ? 1 : 0);
        int belowlevels = (below!=0. ? 1 : 0);
        levels =  nlevels + abovelevels + belowlevels ;

        eps.allocate(grid1,grid2,levels);
        thick.allocate(levels);

        COMPLEX esurround = (COMPLEX)surrounding.epsilon(lambda);
        COMPLEX esphere = (COMPLEX)sphere.epsilon(lambda);
        for (int i=1; i<=grid1*grid2*levels; ++i) eps(i)=esurround;

        if (below!=0.) thick(1) = below;
        if (above!=0.) thick(levels) = above;


        double radius = diameter/2.;
        double t = diameter/nlevels;

        for (int n=1; n<=nlevels; ++n) {

            double r = (n<=(nlevels+1)/2) ? (2*n-1)*diameter/nlevels/2 : (2*(nlevels-n+1)-1)*diameter/nlevels/2;

            if (n==1 || n==nlevels) {
                thick(n+belowlevels)=radius-sqrt(sqr(radius)-sqr(r+radius/nlevels));
            } else if ((double)n==(nlevels+1.)/2.) {
                thick(n+belowlevels)=2*sqrt(sqr(radius)-sqr(r-radius/nlevels));
            } else {
                thick(n+belowlevels)=sqrt(sqr(radius)-sqr(r-radius/nlevels))-sqrt(sqr(radius)-sqr(r+radius/nlevels));
            }

            for (int i=1; i<=grid1; ++i) {
                for (int j=1; j<=grid2; ++j) {
                    double x,y;
                    getxy(i,j,x,y);

                    if (sqr(x)+sqr(y)<sqr(r)) {
                        eps(i,j,n+belowlevels) = esphere;
                    }
                }
            }
        }

        FourierFactorize();
    }

    void Pyramidal_Pit_CrossGrating::setup()
    {
        Gridded_CrossGrating::setup();

        levels = nlevels;

        eps.allocate(grid1,grid2,nlevels);
        thick.allocate(levels);

        COMPLEX medium_t_eps = (COMPLEX)medium_t.epsilon(lambda);
        COMPLEX medium_i_eps = (COMPLEX)medium_i.epsilon(lambda);
        for (int i=1; i<=grid1*grid2*levels; ++i) eps(i)=medium_t_eps;

        for (int n=1; n<=levels; ++n) {
            thick(n) = depth/levels;
            double hside = side*(n-0.5)/levels/2.;
            for (int i=1; i<=grid1; ++i) {
                for (int j=1; j<=grid2; ++j) {
                    double x,y;
                    getxy(i,j,x,y);

                    if (fabs(x)<hside && fabs(y)<hside) eps(i,j,n) = medium_i_eps;
                }
            }
        }
        FourierFactorize();
    }

    //
    // The routine FillE0E1E2E3Matrices fills the E-matrices using 1D Grating data.
    //
    // "level2d" is the level number of the 2D grating.
    // "level1d" is the level number of the 1D grating.
    // "flip" specifies
    // (if false) the grating direction is aligned along the "1" axis, or
    // (if true) the grating direction is aligned along the "2" axis.
    //
    //
	void FillE0E1E2E3Matrices(CFARRAY& EPS0, CFARRAY& EPS11, CFARRAY& EPS12, CFARRAY& EPS2, CFARRAY& EPS3, 
					   		CFARRAY& MU0, CFARRAY& MU11, CFARRAY& MU12, CFARRAY& MU2, CFARRAY& MU3, 
							int M1, int M2, int level2d, int level1d, Grating_Ptr& grating, bool flip)
	{
        int M = flip ? M2 : M1;

        CFARRAY ematrix(M,M);
		CFARRAY mmatrix(M, M);
		for (int i=1; i<=M; ++i) {
            for (int k=1; k<=M; ++k) {
				ematrix(i, k) = (grating->fourierz(i - k, level1d, 0));
				mmatrix(i, k) = (grating->fouriermuz(i - k, level1d, 0));
			}
        }
        Inverse(ematrix,M);
		Inverse(mmatrix, M);
		for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                for (int k=1; k<=M1; ++k) {
                    for (int l=1; l<=M2; ++l) {
                        if (flip) {
							EPS0(i, j, k, l, level2d) = (i == k) ? ematrix(j, l) : 0.;
							MU0(i, j, k, l, level2d) = (i == k) ? mmatrix(j, l) : 0.;
						} else {
							EPS0(i, j, k, l, level2d) = (j == l) ? ematrix(i, k) : 0.;
							MU0(i, j, k, l, level2d) = (j == l) ? mmatrix(i, k) : 0.;
						}
                    }
                }
            }
        }

        for (int i=1; i<=M; ++i) {
            for (int j=1; j<=M; ++j) {
				ematrix(i, j) = (grating->fourierx(i - j, level1d, 1));
				mmatrix(i, j) = (grating->fouriermux(i - j, level1d, 1));
			}
        }
        Inverse(ematrix,M);
		Inverse(mmatrix, M);
		for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                for (int k=1; k<=M1; ++k) {
                    for (int l=1; l<=M2; ++l) {
                        if (flip) {
							EPS12(i, j, k, l, level2d) = (i == k) ? ematrix(j, l) : 0.;
							EPS3(i, j, k, l, level2d) = (i == k) ? ematrix(j, l) : 0.;
							MU12(i, j, k, l, level2d) = (i == k) ? mmatrix(j, l) : 0.;
							MU3(i, j, k, l, level2d) = (i == k) ? mmatrix(j, l) : 0.;
						} else {
							EPS11(i, j, k, l, level2d) = (j == l) ? ematrix(i, k) : 0.;
							EPS2(i, j, k, l, level2d) = (j == l) ? ematrix(i, k) : 0.;
							MU11(i, j, k, l, level2d) = (j == l) ? mmatrix(i, k) : 0.;
							MU2(i, j, k, l, level2d) = (j == l) ? mmatrix(i, k) : 0.;
						}
                    }
                }
            }
        }

        for (int i=1; i<=M; ++i) {
            for (int k=1; k<=M; ++k) {
				ematrix(i, k) = (grating->fouriery(i - k, level1d, 0));
				mmatrix(i, k) = (grating->fouriermuy(i - k, level1d, 0));
			}
        }
        for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                for (int k=1; k<=M1; ++k) {
                    for (int l=1; l<=M2; ++l) {
                        if (flip) {
							EPS11(i, j, k, l, level2d) = (i == k) ? ematrix(j, l) : 0.;
							EPS2(i, j, k, l, level2d) = (i == k) ? ematrix(j, l) : 0.;
							MU11(i, j, k, l, level2d) = (i == k) ? mmatrix(j, l) : 0.;
							MU2(i, j, k, l, level2d) = (i == k) ? mmatrix(j, l) : 0.;
						} else {
							EPS12(i, j, k, l, level2d) = (j == l) ? ematrix(i, k) : 0.;
							EPS3(i, j, k, l, level2d) = (j == l) ? ematrix(i, k) : 0.;
							MU12(i, j, k, l, level2d) = (j == l) ? mmatrix(i, k) : 0.;
							MU3(i, j, k, l, level2d) = (j == l) ? mmatrix(i, k) : 0.;
						}
                    }
                }
            }
        }
    }


    void OneD_CrossGrating::setup()
    {
        CrossGrating::setup();

        if (grating->get_lambda()!=lambda) grating->set_lambda(lambda);
        if ((COMPLEX)grating->get_medium_i().index(lambda)!=(COMPLEX)medium_i.index(lambda)) error("grating.medium_i!=medium_i");
        if ((COMPLEX)grating->get_medium_t().index(lambda)!=(COMPLEX)medium_t.index(lambda)) error("grating.medium_t!=medium_t");

        levels = grating->get_levels();

        double period = grating->get_period();
        d1 = period/cos(zeta*deg);
        CrossGrating::d2 = d2;
        CrossGrating::zeta = zeta;

        thick.allocate(levels);

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        EPS0.allocate(M1,M2,M1,M2,levels);
        EPS11.allocate(M1,M2,M1,M2,levels);
		EPS12.allocate(M1, M2, M1, M2, levels);
		EPS2.allocate(M1,M2,M1,M2,levels);
        EPS3.allocate(M1,M2,M1,M2,levels);
		MU0.allocate(M1, M2, M1, M2, levels);
		MU11.allocate(M1, M2, M1, M2, levels);
		MU12.allocate(M1, M2, M1, M2, levels);
		MU2.allocate(M1, M2, M1, M2, levels);
		MU3.allocate(M1, M2, M1, M2, levels);

        for (int level=1; level<=levels; ++level) {
            thick(level) = grating->get_thickness(levels-level);
            FillE0E1E2E3Matrices(EPS0,EPS11,EPS12,EPS2,EPS3,
				MU0, MU11, MU12, MU2, MU3,
				M1,M2,level,levels-level,grating,false);
        }
    }

    void Overlaid_CrossGrating::setup()
    {
        CrossGrating::setup();

        top->set_lambda(lambda);
        top->set_order1(order1);
        top->set_order2(order2);

        bottom->set_lambda(lambda);
        bottom->set_order1(order1);
        bottom->set_order2(order2);

        if (top->get_zeta()!=bottom->get_zeta()) error("top.zeta != bottom.zeta");
        if (top->get_d1()!=bottom->get_d1()) error("top.d1 != bottom.d1");
        if (top->get_d2()!=bottom->get_d2()) error("top.d2 != bottom.d2");
        if (top->get_medium_i().index(lambda)!=medium_i.index(lambda)) error("top.medium_i!=medium_i");
        if (bottom->get_medium_t().index(lambda)!=medium_t.index(lambda)) error("bottom.medium_t!=medium_t");
        if (top->get_medium_t().index(lambda)!=bottom->get_medium_i().index(lambda)) error("top.medium_t!=bottom.medium_i");

        int tlevels = top->get_levels();
        int blevels = bottom->get_levels();

        levels = tlevels + blevels;
        if (separation>0.) ++levels;

        d1 = top->get_d1();
        d2 = top->get_d2();
        zeta = top->get_zeta();

        CFARRAY tE0 = top->get_EPS0();
        CFARRAY tE11 = top->get_EPS11();
		CFARRAY tE12 = top->get_EPS12();
		CFARRAY tE2 = top->get_EPS2();
        CFARRAY tE3 = top->get_EPS3();
		CFARRAY tMU0 = top->get_MU0();
		CFARRAY tMU11 = top->get_MU11();
		CFARRAY tMU12 = top->get_MU12();
		CFARRAY tMU2 = top->get_MU2();
		CFARRAY tMU3 = top->get_MU3();

        CFARRAY bE0 = bottom->get_EPS0();
        CFARRAY bE11 = bottom->get_EPS11();
		CFARRAY bE12 = bottom->get_EPS12();
		CFARRAY bE2 = bottom->get_EPS2();
        CFARRAY bE3 = bottom->get_EPS3();
		CFARRAY bMU0 = bottom->get_MU0();
		CFARRAY bMU11 = bottom->get_MU11();
		CFARRAY bMU12 = bottom->get_MU12();
		CFARRAY bMU2 = bottom->get_MU2();
		CFARRAY bMU3 = bottom->get_MU3();

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        EPS0.allocate(M1,M2,M1,M2,levels);
        EPS11.allocate(M1,M2,M1,M2,levels);
		EPS12.allocate(M1,M2,M1,M2,levels);
		EPS2.allocate(M1,M2,M1,M2,levels);
        EPS3.allocate(M1,M2,M1,M2,levels);
		MU0.allocate(M1, M2, M1, M2, levels);
		MU11.allocate(M1, M2, M1, M2, levels);
		MU12.allocate(M1, M2, M1, M2, levels);
		MU2.allocate(M1, M2, M1, M2, levels);
		MU3.allocate(M1, M2, M1, M2, levels);
		thick.allocate(levels);

        double k1 = 2*pi/d1;
        double k2 = 2*pi/d2;
        COMPLEX cI(0,1);

        for (int level=1; level<=levels; ++level) {
            if (level<=blevels) {
                int blevel = level;
                thick(level) = bottom->get_thick(blevel);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                EPS0(i,j,k,l,level) = bE0(i,j,k,l,blevel);
                                EPS11(i,j,k,l,level) = bE11(i,j,k,l,blevel);
								EPS12(i, j, k, l, level) = bE12(i,j,k,l,blevel);
								EPS2(i,j,k,l,level) = bE2(i,j,k,l,blevel);
                                EPS3(i,j,k,l,level) = bE3(i,j,k,l,blevel);
								MU0(i, j, k, l, level) = bMU0(i, j, k, l, blevel);
								MU11(i, j, k, l, level) = bMU11(i, j, k, l, blevel);
								MU12(i, j, k, l, level) = bMU12(i, j, k, l, blevel);
								MU2(i, j, k, l, level) = bMU2(i, j, k, l, blevel);
								MU3(i, j, k, l, level) = bMU3(i, j, k, l, blevel);
							}
                        }
                    }
                }
            }
            if (level==blevels+1 && separation>0.) {
                thick(level) = separation;
                COMPLEX smedium = top->get_medium_t().epsilon(lambda);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                EPS0(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                EPS11(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
								EPS12(i,j,k,l,level) = (i == k && j == l) ? smedium : 0.;
								EPS2(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                EPS3(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
								MU0(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU11(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU12(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU2(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU3(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
							}
                        }
                    }
                }
            }
            int iseparation = separation>0.? 1 : 0;
            if (level>blevels + iseparation) {
                int tlevel = level - blevels - iseparation;
                thick(level) = top->get_thick(tlevel);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                COMPLEX phase = exp(-cI*(k1*overlay1*(i-k)+k2*overlay2*(j-l)));
                                EPS0(i,j,k,l,level) = tE0(i,j,k,l,tlevel)*phase;
                                EPS11(i,j,k,l,level) = tE11(i,j,k,l,tlevel)*phase;
								EPS12(i,j,k,l,level) = tE12(i,j,k,l,tlevel)*phase;
								EPS2(i,j,k,l,level) = tE2(i,j,k,l,tlevel)*phase;
                                EPS3(i,j,k,l,level) = tE3(i,j,k,l,tlevel)*phase;
								MU0(i, j, k, l, level) = tMU0(i, j, k, l, tlevel)*phase;
								MU11(i, j, k, l, level) = tMU11(i, j, k, l, tlevel)*phase;
								MU12(i, j, k, l, level) = tMU12(i, j, k, l, tlevel)*phase;
								MU2(i, j, k, l, level) = tMU2(i, j, k, l, tlevel)*phase;
								MU3(i, j, k, l, level) = tMU3(i, j, k, l, tlevel)*phase;
							}
                        }
                    }
                }
            }
        }
    }


    void Overlaid_1D_CrossGrating::setup()
    {
        CrossGrating::setup();

        top->set_lambda(lambda);

        bottom->set_lambda(lambda);

        if (top->get_medium_t().index(lambda)!=bottom->get_medium_i().index(lambda)) error("top.medium_t!=bottom.medium_i");
        if (top->get_medium_i().index(lambda)!=medium_i.index(lambda)) error("top.medium_i!=medium_i");
        if (bottom->get_medium_t().index(lambda)!=medium_t.index(lambda)) error("bottom.medium_t!=medium_t");
        if (angle == 0.) error("angle=0 not supported");

        zeta = 90.-angle;
        d1 = top->get_period()/sin(angle*deg);
        d2 = bottom->get_period()/sin(angle*deg);

        int tlevels = top->get_levels();
        int blevels = bottom->get_levels();

        levels = tlevels + blevels;
        if (separation>0.) ++levels;

        thick.allocate(levels);

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        EPS0.allocate(M1,M2,M1,M2,levels);
        EPS11.allocate(M1,M2,M1,M2,levels);
		EPS12.allocate(M1,M2,M1,M2,levels);
		EPS2.allocate(M1,M2,M1,M2,levels);
        EPS3.allocate(M1,M2,M1,M2,levels);
		MU0.allocate(M1, M2, M1, M2, levels);
		MU11.allocate(M1, M2, M1, M2, levels);
		MU12.allocate(M1, M2, M1, M2, levels);
		MU2.allocate(M1, M2, M1, M2, levels);
		MU3.allocate(M1, M2, M1, M2, levels);

        int iseparation = separation>0.? 1 : 0;

        for (int level=1; level<=levels; ++level) {
            if (level<=blevels) {
                int blevel = blevels-level;
                thick(level) = bottom->get_thickness(blevel);
                FillE0E1E2E3Matrices(EPS0,EPS11,EPS12,EPS2,EPS3,
					MU0, MU11, MU12, MU2, MU3,
					M1,M2,level,blevel,bottom,true);
			} else if (level==blevels+1 && separation>0.) {
                thick(level) = separation;
                COMPLEX smedium = top->get_medium_t().epsilon(lambda);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                EPS0(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                EPS11(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
								EPS12(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
								EPS2(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                EPS3(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
								MU0(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU11(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU12(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU2(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU3(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
							}
                        }
                    }
                }
            } else if (level>blevels + iseparation) {
                int tlevel = blevels-level+blevels+iseparation;
                thick(level) = top->get_thickness(tlevel);
                FillE0E1E2E3Matrices(EPS0,EPS11,EPS12,EPS2,EPS3,
					MU0, MU11, MU12, MU2, MU3,
					M1,M2,level,tlevel,top,false);
				for (int i = 1; i <= M1; ++i) {
					for (int j = 1; j <= M2; ++j) {
						for (int k = 1; k <= M1; ++k) {
							for (int l = 1; l <= M2; ++l) {
								MU0(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU11(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU12(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU2(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
								MU3(i, j, k, l, level) = (i == k && j == l) ? 1. : 0.;
							}
						}
					}
				}

            }
        }
    }


    void Null_CrossGrating::setup()
    {
        CrossGrating::setup();

        levels = 0;

        CrossGrating::d1 = d1;
        CrossGrating::d2 = d2;
        CrossGrating::zeta = zeta;

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        EPS0.allocate(M1,M2,M1,M2,0);
        EPS11.allocate(M1,M2,M1,M2,0);
		EPS12.allocate(M1,M2,M1,M2,0);
		EPS2.allocate(M1,M2,M1,M2,0);
        EPS3.allocate(M1,M2,M1,M2,0);
		MU0.allocate(M1, M2, M1, M2, 0);
		MU11.allocate(M1, M2, M1, M2, 0);
		MU12.allocate(M1, M2, M1, M2, 0);
		MU2.allocate(M1, M2, M1, M2, 0);
		MU3.allocate(M1, M2, M1, M2, 0);
		thick.allocate(levels);
    }



    DEFINE_MODEL(Cylinder_CrossGrating,Gridded_CrossGrating,"Contact holes in a 2-d array");
    DEFINE_PARAMETER(Cylinder_CrossGrating,Table,rtop,"Radius of top of holes [um] as function of angle [deg]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,Table,rbottom,"Radius of bottom of holes [um] as function of angle [deg]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,double,thickness,"Thickness of grating [um]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,int,nlevels,"Number of levels in grating","1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,dielectric_function,inside,"Medium inside holes","(1,0)",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,dielectric_function,outside,"Medium outside holes","(1.5,0)",0xFF);

	DEFINE_MODEL(Rectangle_CrossGrating, Gridded_CrossGrating, "Rectangles in a 2-d array");
	DEFINE_PARAMETER(Rectangle_CrossGrating, double, length1, "Length of rectangle in first direction [um]", "0.1", 0xFF);
	DEFINE_PARAMETER(Rectangle_CrossGrating, double, length2, "Length of rectangle in second direction [um]", "0.1", 0xFF);
	DEFINE_PARAMETER(Rectangle_CrossGrating, double, zetaa, "Skew angle [deg]", "0", 0xFF);
	DEFINE_PARAMETER(Rectangle_CrossGrating, double, thickness, "Thickness of grating [um]", "0.1", 0xFF);
	DEFINE_PARAMETER(Rectangle_CrossGrating, dielectric_function, inside, "Medium inside rectangles", "(1,0)", 0xFF);
	DEFINE_PARAMETER(Rectangle_CrossGrating, dielectric_function, outside, "Medium outside rectangles", "(1.5,0)", 0xFF);

	DEFINE_MODEL(OneD_CrossGrating,CrossGrating,"One dimensional grating");
    DEFINE_PARAMETER(OneD_CrossGrating,double,d2,"Lattice constant #2 [um]","0.5",0xFF);
    DEFINE_PARAMETER(OneD_CrossGrating,double,zeta,"Angle of lattice vectors from perpendicular [deg]","0",0xFF);
    DEFINE_PTRPARAMETER(OneD_CrossGrating,Grating_Ptr,grating,"Grating","Single_Line_Grating",0xFF);

    DEFINE_MODEL(Overlaid_CrossGrating,CrossGrating,"A crossed grating on top of another");
    DEFINE_PTRPARAMETER(Overlaid_CrossGrating,CrossGrating_Ptr,top,"Top grating","OneD_CrossGrating",0xFF);
    DEFINE_PTRPARAMETER(Overlaid_CrossGrating,CrossGrating_Ptr,bottom,"Bottom grating","OneD_CrossGrating",0xFF);
    DEFINE_PARAMETER(Overlaid_CrossGrating,double,overlay1,"Overlay along 1st coordinate [um]","0",0xFF);
    DEFINE_PARAMETER(Overlaid_CrossGrating,double,overlay2,"Overlay along 2nd coordinate [um]","0",0xFF);
    DEFINE_PARAMETER(Overlaid_CrossGrating,double,separation,"Vertical separation between gratings [um]","0",0xFF);

    DEFINE_MODEL(Overlaid_1D_CrossGrating,CrossGrating,"One 1D grating on top of another at some angle");
    DEFINE_PTRPARAMETER(Overlaid_1D_CrossGrating,Grating_Ptr,top,"Top grating","Single_Line_Grating",0xFF);
    DEFINE_PTRPARAMETER(Overlaid_1D_CrossGrating,Grating_Ptr,bottom,"Bottom grating","Single_Line_Grating",0xFF);
    DEFINE_PARAMETER(Overlaid_1D_CrossGrating,double,angle,"Angle between gratings [deg]","90",0xFF);
    DEFINE_PARAMETER(Overlaid_1D_CrossGrating,double,separation,"Vertical separation between gratings [um]","0",0xFF);

    DEFINE_MODEL(Null_CrossGrating,CrossGrating,"A trivial crossed grating with no layers");
    DEFINE_PARAMETER(Null_CrossGrating,double,d1,"Lattice constant #1 [um]","0.5",0xFF);
    DEFINE_PARAMETER(Null_CrossGrating,double,d2,"Lattice constant #2 [um]","0.5",0xFF);
    DEFINE_PARAMETER(Null_CrossGrating,double,zeta,"Angle of lattice vectors from perpendicular [deg]","0",0xFF);

    DEFINE_MODEL(Sphere_CrossGrating,Gridded_CrossGrating,"Grating consisting of spheres embedded in a medium");
    DEFINE_PARAMETER(Sphere_CrossGrating,double,diameter,"Diameter of sphere [um]","0.05",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,double,above,"Distance of sphere from top of layer [um]","0.05",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,double,below,"Distance of sphere from bottom of layer [um]","0.05",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,int,nlevels,"Number of levels in structure","7",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,dielectric_function,sphere,"Sphere","(1,0)",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,dielectric_function,surrounding,"Medium around sphere","(1.5,0)",0xFF);

    DEFINE_MODEL(Pyramidal_Pit_CrossGrating,Gridded_CrossGrating,"Grating consisting of pyramidal pits");
    DEFINE_PARAMETER(Pyramidal_Pit_CrossGrating,double,side,"Length of side of base of pyramid [um]","0.05",0xFF);
    DEFINE_PARAMETER(Pyramidal_Pit_CrossGrating,double,depth,"Depth of pit [um]","0.05",0xFF);
    DEFINE_PARAMETER(Pyramidal_Pit_CrossGrating,int,nlevels,"Number of levels","10",0xFF);


}

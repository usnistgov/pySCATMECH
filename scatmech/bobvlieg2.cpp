//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: bobvlieg2.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include <float.h>
#include <assert.h>
#include <vector>
#include <cstdlib>
#include "scatmech.h"
#include "bobvlieg.h"

using namespace std;

namespace SCATMECH {

    using namespace BobVlieg_Supp;

    //
    // factorial_list is a lookup table for n!
    //
    static double factorial_list[171]=
    {   1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3.6288e6,3.99168e7,4.790016e8,6.2270208e9,
        8.71782912e10,1.307674368e12,2.0922789888e13,3.55687428096e14,6.402373705728e15,
        1.21645100408832e17,
        2.43290200817664e18,5.109094217170944e19,1.1240007277776077e21,2.585201673888498e22,
        6.204484017332394e23,1.5511210043330986e25,4.0329146112660565e26,1.0888869450418352e28,
        3.0488834461171387e29,8.841761993739702e30,2.6525285981219107e32,8.222838654177922e33,
        2.631308369336935e35,8.683317618811886e36,2.9523279903960416e38,1.0333147966386145e40,
        3.7199332678990125e41,1.3763753091226346e43,5.230226174666011e44,2.0397882081197444e46,
        8.159152832478977e47,3.345252661316381e49,1.40500611775288e51,6.041526306337383e52,
        2.658271574788449e54,1.1962222086548019e56,5.502622159812089e57,2.5862324151116818e59,
        1.2413915592536073e61,6.082818640342675e62,3.0414093201713376e64,1.5511187532873822e66,
        8.065817517094388e67,4.2748832840600255e69,2.308436973392414e71,1.2696403353658276e73,
        7.109985878048635e74,4.0526919504877214e76,2.3505613312828785e78,1.3868311854568984e80,
        8.32098711274139e81,5.075802138772248e83,3.146997326038794e85,1.98260831540444e87,
        1.2688693218588417e89,8.247650592082472e90,5.443449390774431e92,3.647111091818868e94,
        2.4800355424368305e96,1.711224524281413e98,1.1978571669969892e100,8.504785885678623e101,
        6.1234458376886085e103,4.4701154615126844e105,3.307885441519386e107,2.48091408113954e109,
        1.8854947016660504e111,1.4518309202828587e113,1.1324281178206297e115,8.946182130782976e116,
        7.156945704626381e118,5.797126020747368e120,4.753643337012842e122,3.945523969720659e124,
        3.314240134565353e126,2.81710411438055e128,2.4227095383672734e130,2.107757298379528e132,
        1.8548264225739844e134,1.650795516090846e136,1.4857159644817615e138,1.352001527678403e140,
        1.2438414054641308e142,1.1567725070816416e144,1.087366156656743e146,1.032997848823906e148,
        9.916779348709496e149,9.619275968248212e151,9.426890448883248e153,9.332621544394415e155,
        9.332621544394415e157,9.42594775983836e159,9.614466715035127e161,9.90290071648618e163,
        1.0299016745145628e166,1.081396758240291e168,1.1462805637347084e170,1.226520203196138e172,
        1.324641819451829e174,1.4438595832024937e176,1.588245541522743e178,1.7629525510902446e180,
        1.974506857221074e182,2.2311927486598138e184,2.5435597334721877e186,2.925093693493016e188,
        3.393108684451898e190,3.969937160808721e192,4.684525849754291e194,5.574585761207606e196,
        6.689502913449127e198,8.094298525273444e200,9.875044200833601e202,1.214630436702533e205,
        1.506141741511141e207,1.882677176888926e209,2.372173242880047e211,3.0126600184576594e213,
        3.856204823625804e215,4.974504222477287e217,6.466855489220474e219,8.47158069087882e221,
        1.1182486511960043e224,1.4872707060906857e226,1.9929427461615188e228,2.6904727073180504e230,
        3.659042881952549e232,5.012888748274992e234,6.917786472619489e236,9.615723196941089e238,
        1.3462012475717526e241,1.898143759076171e243,2.695364137888163e245,3.854370717180073e247,
        5.5502938327393044e249,8.047926057471992e251,1.1749972043909107e254,1.727245890454639e256,
        2.5563239178728654e258,3.80892263763057e260,5.713383956445855e262,8.62720977423324e264,
        1.3113358856834524e267,2.0063439050956823e269,3.0897696138473508e271,4.789142901463394e273,
        7.471062926282894e275,1.1729568794264145e278,1.853271869493735e280,2.9467022724950384e282,
        4.7147236359920616e284,7.590705053947219e286,1.2296942187394494e289,2.0044015765453026e291,
        3.287218585534296e293,5.423910666131589e295,9.003691705778438e297,1.503616514864999e300,
        2.5260757449731984e302,4.269068009004705e304,7.257415615307999e306
    };

    //
    // sqrtfactorial_list is a lookup table for sqrt(n!)
    //
    static double sqrtfactorial_list[171] =
    {   1.,1.,1.4142135623730951,2.449489742783178,4.898979485566356,10.954451150103322,26.832815729997478,
        70.9929573971954,200.79840636817812,602.3952191045344,1904.9409439665053,6317.974358922328,
        21886.105181141756,78911.47445080469,295259.70128007646,1.143535905863913e6,4.574143623455652e6,
        1.8859677306253146e7,8.001483428544985e7,3.487765766344294e8,1.5597762686284978e9,7.147792818185865e9,
        3.352612008237171e10,1.6078562354540588e11,7.876854713229384e11,3.9384273566146914e12,
        2.008211794424596e13,1.0434974580907397e14,5.521669535672285e14,2.973510046012911e15,
        1.6286585271694956e16,9.067986906793549e16,5.129628026803635e17,2.9467469553410734e18,
        1.7182339742875652e19,1.016520927791757e20,6.099125566750543e20,3.7099532465014094e21,
        2.28696877430935e22,1.4282115417961528e23,9.032802905233225e23,5.783815921445271e24,
        3.748341123420973e25,2.4579516484946126e26,1.6304206741784308e27,1.0937194378152021e28,
        7.417966136220959e28,5.085501366740237e29,3.523338699662023e30,2.466337089763416e31,
        1.743963680863606e32,1.2454391808865588e33,8.980989654316716e33,6.538259159791715e34,
        4.804619624270389e35,3.56320127885842e36,2.666455677120592e37,2.013129889124823e38,
        1.5331540468207618e39,1.1776379687564844e40,9.121944481710788e40,7.124466393192017e41,
        5.609810447812647e42,4.452649004137245e43,3.562119203309796e44,2.8718723147247465e45,
        2.3331200978034607e46,1.909741105966688e47,1.574812859496909e48,1.3081378078327272e49,
        1.094466613011557e50,9.222139602976429e50,7.825244940376378e51,6.685892207860283e52,
        5.751421947239992e53,4.980877514193197e54,4.3422283469044446e55,3.810289910601106e56,
        3.3651569321810684e57,2.991016905800263e58,2.6752468492881884e59,2.4077221643593693e60,
        2.1802851503903892e61,1.986334304622628e62,1.820505461284133e63,1.678423103505356e64,
        1.5565055535934566e65,1.4518117296604018e66,1.361920123419132e67,1.2848328747704294e68,
        1.2188994890809338e69,1.1627560052213888e70,1.1152763807523813e71,1.0755335917960172e72,
        1.0427685057848376e73,1.0163650175128551e74,9.958302741285531e74,9.807790764615756e75,
        9.709217501366034e76,9.66054943799493e77,9.66054943799493e78,9.708732028353835e79,
        9.805338706559364e80,9.951331929187258e81,1.0148407138632952e83,1.039902283024848e84,
        1.0706449288791818e85,1.1074837259283487e86,1.1509308491181515e87,1.2016070835354183e88,
        1.2602561412358772e89,1.327762234396748e90,1.4051714689748985e91,1.4937177607097713e92,
        1.5948541417547208e93,1.7102905289724946e94,1.8420392733196266e95,1.992470115411702e96,
        2.164376549899368e97,2.3610560690520687e98,2.586407337108586e99,2.845048070819445e100,
        3.142458305345292e101,3.485154855530143e102,3.880904200712948e103,4.33898280347932e104,
        4.870496117317051e105,5.488770370909735e106,6.209834799433722e107,7.053016533709026e108,
        8.041676124552936e109,9.20411901861271e110,1.0574727661722567e112,1.2195370868041225e113,
        1.4117162413748447e114,1.6402660477245912e115,1.9128624838060233e116,2.238948134342328e117,
        2.6301685255168514e118,3.1009229588851587e119,3.669061525201986e120,4.356769168863748e121,
        5.1916896458553485e122,6.208357848239801e123,7.45002941788776e124,8.971023385028037e125,
        1.083972879914858e127,1.3142472714275037e128,1.598850811637179e129,1.9516461353510197e130,
        2.390268595042376e131,2.937211223973046e132,3.621237199747419e133,4.479223040992358e134,
        5.558569612631789e135,6.920363358569689e136,8.643531064491464e137,1.083031338155279e139,
        1.3613492826948325e140,1.7165961296982576e141,2.1713414369905209e142,2.7551234190045313e143,
        3.506699614651146e144,4.477054362575132e145,5.733427060261861e146,7.364720406187589e147,
        9.488778480804807e148,1.2262204185483942e150,1.5893633143410597e151,2.0661723086433777e152,
        2.6939590968142033e153
    };

    namespace BobVlieg_Supp {

        //
        // The following table containes lvector[L]=sqrt((2.*L+1.)/(L*(L+1.)))
        //
        double lvector[] = {0,1.224744871391589,0.9128709291752768,0.7637626158259733,
                            0.6708203932499369,0.6055300708194983,0.5563486402641868,0.5175491695067657,
                            0.485912657903775,0.45946829173634074,0.43693144875265144,0.4174235549683609,
                            0.40032038451271784,0.38516444325982163,0.37161167647860327,0.3593976442141304,
                            0.34831527300961795,0.3381997707972616,0.32891812735531006,0.3203616377585937,
                            0.3124404705204619,0.30507965037608303,0.29821603968282906,0.2917960375608824,
                            0.2857738033247041,0.28010986855435577,0.27477004112270953,0.26972453123756235,
                            0.26494724821174376,0.26041522988109395,0.2561081760691415,0.25200806438709267,
                            0.24809883172443659,0.2443661085521326,0.2407969959889677,0.2373798777259908,
                            0.234104260543897,0.230960638422895,0.22794037622744648,0.22503560971771921,
                            0.22223915924615767,0.21954445497885577,0.2169454718656561,0.21443667289146487,
                            0.2120129593904461,0.20966962740703488,0.2074023292527436,0.2052070395430288,
                            0.20308002510990453,0.20101781827814696,0.19901719306948054,0.1970751439629669,
                            0.19518886689325063,0.19335574221320762,0.19157331938539016,0.1898393031986818,
                            0.18815154133374953,0.18650801312401472,0.18490681937862027,0.18334617315079021,
                            0.18182439134950373,0.1803398871049201,0.17889116280878825,0.1774768037604235,
                            0.17609547235694445,0.17474590277351437,0.17342689608547893,0.17213731578966118,
                            0.17087608368677457,0.16964217609103785,0.16843462033669862,0.16725249155436342,
                            0.16609490969284757,0.1649610367647458,0.16385007429612855,0.16276126096272245,
                            0.16169387039666994,0.16064720914950814,0.15962061479838402,0.15861345418375222,
                            0.1576251217679012,0.15665503810463674,0.15570264841133516,0.1547674212353693,
                            0.15384884720762404,0.15294643787645767,0.1520597246160449,0.15118825760355725,
                            0.15033160486010963,0.1494893513508266,0.14866109813976924,0.14784646159581327,
                            0.14704507264588745,0.1462565760722691,0.1454806298508984,0.14471690452791044,
                            0.1439650826318039,0.14322485811886332,0.14249593584963455,0.14177803109441922,
                            0.14107086906590566,0.1403741844771942,0.1396877211236017,0.13901123148674815,
                            0.13834447635953573,0.13768722449072923,0.1370392522479392,0.13640034329789139,
                            0.13577028830294485,0.13514888463289143,0.1345359360911355
                           };
        //
        // The complex arccosine is not part of the standard C++ library (!)...
        //
        COMPLEX arcsine(const COMPLEX& a)
        {
            // Old version
            //COMPLEX q = sqrt(1. - sqr(a)) + COMPLEX(0.,1.)*a;
            //return COMPLEX(arg(q),-log(abs(q)));
            COMPLEX i(0.,1.);
            return -i*log(i*a+sqrt(1.-a*a));
        }
        //
        // The complex arccosine is not part of the standard C++ library (!)...
        //
        COMPLEX
        arccosine(const COMPLEX& a)
        {
            double x=real(a);
            double y=imag(a);

            return pi/2. - arg(sqrt(1. - sqr(x + cI*y)) +
                               cI*(x + cI*y)) +
                   cI*log(abs(sqrt(1. - sqr(x + cI*y)) +
                              cI*(x + cI*y)));
        }

        //
        // Routine to return the factorial of j
        //
        double
        Fact(int j)
        {
            if (j <= 170 && j>=0) return factorial_list[j];

            throw SCATMECH_exception("Fact("+to_string(j)+") out of range.");
            return 0.;
        }

        //
        // Routine to return the square root of the factorial of j
        //
        double
        sqrtFact(int j)
        {
            if (j <= 170 && j>=0) return sqrtfactorial_list[j];
            throw SCATMECH_exception("sqrtFact("+to_string(j)+") out of range.");
            return 0.;
        }

        //
        // mpow returns (-1)**m
        //
        double
        mpow(int m)
        {
            int result;
            if (m<0) m=-m;
            if (m & 0x1) result = -1;
            else result = 1;
            return result;
        }

        //
        // ipow returns pow(sqrt(-1.),m)...
        //
        COMPLEX
        ipow(int m)
        {
            COMPLEX result;
            if (m>=0) {
                switch (m%4) {
                    case (0):
                        result= COMPLEX(1.,0.);
                        break;
                    case (1):
                        result= COMPLEX(0.,1.);
                        break;
                    case (2):
                        result= COMPLEX(-1.,0.);
                        break;
                    case (3):
                        result= COMPLEX(0.,-1.);
                        break;
                    default:
                        result= COMPLEX(1.,0.);
                        break;
                }
            } else {
                m=-m;
                switch (m%4) {
                    case (0):
                        result= COMPLEX(1.,0.);
                        break;
                    case (1):
                        result= COMPLEX(0.,-1.);
                        break;
                    case (2):
                        result= COMPLEX(-1.,0.);
                        break;
                    case (3):
                        result= COMPLEX(0.,1.);
                        break;
                    default:
                        result= COMPLEX(1.,0.);
                        break;
                }
            }
            return result;
        }

        //
        // Routine to calculate the associated Legendre polynomials
        //
        COMPLEX
        LegendreP(int l,int m,COMPLEX x)
        {
            COMPLEX temp1,temp2,temp3,temp4,result;
            double temp5;
            int i,ll;

            if (m>l) return 0.;

            if (m<0) {
                result = mpow(-m)*Fact(l+m)/Fact(l-m)*LegendreP(l,-m,x);
            } else {
                temp3=1.0;

                if (m>0) {
                    temp1=sqrt(1.0-sqr(x));
                    temp5 = 1.0;
                    for (i=1; i<=m; ++i) {
                        temp3 *= -temp5*temp1;
                        temp5 += 2.0;
                    }
                }
                if (l==m) {
                    result = temp3;
                } else {
                    temp4=x*(2.*m+1.)*temp3;
                    if (l==(m+1)) {
                        result = temp4;
                    } else {
                        for (ll=(m+2); ll<=l; ++ll) {
                            temp2 = (x*(2.*ll-1.)*temp4-(ll+m-1.)*temp3)/(double)(ll-m);
                            temp3=temp4;
                            temp4=temp2;
                        }
                        result = temp2;
                    }
                }
            }
            return result;
        }

        //
        // Routine to calculate the associated Legendre polynomials with the
        // convention (-1)^m used by B&V...
        //
        COMPLEX
        Legendre(int l,int m,COMPLEX x)
        {
            COMPLEX result;
            if (m>l||m<-l) result = 0;
            else if (m>=0) result = LegendreP(l,m,x);
            else result = mpow(-m)*Fact(l+m)/Fact(l-m)*LegendreP(l,-m,x);
            result *=mpow(m);
            return result;
        }

        //
        // Routine to calculate the normalized associated Legendre polynomials...
        //
        COMPLEX
        Ptilde(int l,int m,COMPLEX x)
        {
            COMPLEX result;
            if (m>l || m<-l) {
                result = 0.;
            } else {
                double y=(double)(2.*l+1.)*Fact(l-m)/Fact(l+m);
                result = sqrt(y) * Legendre(l,m,x);
            }
            return result;
        }

        //
        // Routine to calculate the derivative of the associated Legendre polynomials...
        //
        COMPLEX
        dLegendreCosx_dx(int l,int m,COMPLEX cosx)
        {
            COMPLEX sinx = sqrt(1.-sqr(cosx));
            COMPLEX cscx = 1./sinx;

            COMPLEX result = -(cscx*((double)(l + m)*Legendre(l-1,m,cosx) -
                                     (double)(l)*cosx*Legendre(l,m,cosx)));
            return result;
        }

        //
        // Routine to calculate the derivative of the normalized
        // associated Legendre polynomials...
        //
        COMPLEX
        P_tilde(int l,int m,COMPLEX x)
        {
            COMPLEX result;
            if (m>l || m<-l) result = 0.;
            else result = sqrt((2.*l+1.)*Fact(l-m)/Fact(l+m))*dLegendreCosx_dx(l,m,x);
            return result;
        }

        //
        // The following is the ratio J_{\nu-1}(z)/J_\nu(z) using the continued fraction
        // method described in W.J.Lentz, "Generating Bessel functions in Mie Scattering
        // calculations using continued fractions," Appl. Opt. 15(3), 668-671 (1976).
        // This routine implements Eq. (5) of Lentz.
        //
        COMPLEX
        ContinuedFraction_Jratio(double nu, const COMPLEX& x) {
            COMPLEX f0 = 2.*nu/x;
            COMPLEX C0 = f0;
            COMPLEX D0 = 0.;
            int j=0;
            COMPLEX xx = x*x;
            COMPLEX fj;
            for (int j=1; true; ++j) {
                double mpow = j%2==0 ? 1. : -1.;
                COMPLEX bj = mpow*2.*(nu+j)/x;
                COMPLEX aj = 1.;

                COMPLEX Dj = bj + aj*D0;
                if (Dj == 0.) Dj = sqr(std::numeric_limits<double>::epsilon());
                COMPLEX Cj = bj + aj/C0;
                if (Cj == 0.) Cj = sqr(std::numeric_limits<double>::epsilon());
                Dj = 1./Dj;
                COMPLEX Deltaj = Cj*Dj;
                fj = f0*Deltaj;
                f0 = fj;
                C0 = Cj;
                D0 = Dj;
                if (abs(Deltaj-1.)<=std::numeric_limits<double>::epsilon()) break;
            }

            return f0;
        }

        //
        // This routine implements Eq. (3) of Lentz.
        //
        COMPLEX
        ContinuedFraction_A(double n, const COMPLEX& x) {
            COMPLEX f0 = ContinuedFraction_Jratio(n+0.5,x);
            return -n/x+f0;
        }

        //
        // The following structure is used for a map between (l,rho) and
        // values.
        //
        struct lrhopair {
            lrhopair(int _l,COMPLEX& _rho) : l(_l), rho(_rho) {}
            bool operator<(const lrhopair& p) const {
                if (l<p.l) return true;
                if (l>p.l) return false;
                if (rho.real()<p.rho.real()) return true;
                if (rho.real()>p.rho.real()) return false;
                if (rho.imag()<p.rho.imag()) return true;
                return false;
            }
            int l;
            COMPLEX rho;
        };

        // This is the map datatype...
        typedef map<lrhopair,COMPLEX> lrhopairmap;

        //
        // The spherical Bessel function j(l,rho)...
        //
        COMPLEX
        j(int l,COMPLEX rho)
        {
            COMPLEX result;

            if (rho==0.) throw SCATMECH_exception("Encountered j(l,0.)");

            static lrhopairmap previous;
            static mutex_t mutex;
            mutex.lock();
            // First look in the tables to see if it has already been determined...
            lrhopairmap::iterator p = previous.find(lrhopair(l,rho));
			if (p!=previous.end()) {
				mutex.unlock();
				return p->second;
			}
            COMPLEX zm = sin(rho)/rho;
            previous[lrhopair(0,rho)]=zm;
			if (l==0) {
				mutex.unlock();
				return zm;
			}

            COMPLEX zmm = zm;

            zm = (zm-cos(rho))/rho;
            previous[lrhopair(1,rho)]=zm;
			if (l==1){
				mutex.unlock();
				return zm;
			}

            for (int i=2; i<=l; ++i) {
                double nu = i+0.5;
                COMPLEX cf = ContinuedFraction_Jratio(nu,rho);
                COMPLEX zmmm = zm;
                if (norm(cf)<1E-8) {
                    zm = (double)(2*i-1)*zm/rho - zmm;
                } else {
                    zm /= cf;
                }
                previous[lrhopair(i,rho)]=zm;
                zmm = zmmm;
            }
            result = zm;
            mutex.unlock();
            return result;
        }

        //
        // The spherical Bessel function y(l,rho)...
        //
        COMPLEX
        y(int l,COMPLEX rho)
        {
            COMPLEX result;

            static lrhopairmap previous;
            static mutex_t mutex;
            mutex.lock();
            lrhopairmap::iterator p = previous.find(lrhopair(l,rho));
			if (p!=previous.end()) {
				mutex.unlock();
				return p->second;
			}

            if (rho==0.)  throw SCATMECH_exception("Encountered y(l,0.)");

            COMPLEX z0 = -cos(rho)/rho;
            if (l==0) {
                previous[lrhopair(0,rho)]=z0;
                result = z0;
            } else {
                COMPLEX z1 = (z0-sin(rho))/rho;
                previous[lrhopair(1,rho)]=z1;
                if (l==1) {
                    result = z1;
                } else {
                    COMPLEX z2;
                    int n;
                    for (n=1; n<l; ++n) {
                        z2 = (double)(2*n+1)/rho*z1-z0;
                        previous[lrhopair(n+1,rho)]=z2;
                        z0=z1;
                        z1=z2;
                    }
                    result = z2;
                }
            }

            mutex.unlock();
            return result;
        }

        //
        // The spherical Bessel function of the third kind, h(l,rho)...
        //
        COMPLEX
        h(int l,COMPLEX rho)
        {
            COMPLEX result = j(l,rho)+cI*y(l,rho);
            return result;
        }

        //
        // The Ricatti-Bessel function psi(l,rho)
        //
        COMPLEX
        psi(int l,COMPLEX rho)
        {
            COMPLEX result = rho*j(l,rho);
            return result;
        }

        //
        // The Ricatti-Bessel function zeta(l,rho)
        //
        COMPLEX
        zeta(int l,COMPLEX rho)
        {
            COMPLEX result = rho*h(l,rho);
            return result;
        }

        //
        // The derivative of the Ricatti-Bessel function psi
        //
        COMPLEX
        psi_(int l,COMPLEX rho)
        {
            static lrhopairmap previous;

            mutex_t mutex;
            mutex.lock();
            lrhopairmap::iterator p = previous.find(lrhopair(l,rho));
			if (p!=previous.end()) {
				mutex.unlock();
				return p->second;
			}
            COMPLEX result = psi(l,rho)*ContinuedFraction_A(l,rho);

            previous[lrhopair(l,rho)]=result;
            mutex.unlock();
            return result;
        }

        //
        // The derivative of the Ricatti-Bessel function zeta
        //
        COMPLEX
        zeta_(int l,COMPLEX rho)
        {
            COMPLEX result = h(l,rho)+rho/(2.*l+1.)*
                             ((double)(l)*h(l-1,rho)-(l+1.)*h(l+1,rho));

            return result;
        }

        //
        // The Ricatti-Bessel function chi(l,rho)
        //
        COMPLEX
        chi(int l,COMPLEX rho)
        {
            COMPLEX result = -rho*y(l,rho);
            return result;
        }

        //
        // The derivative of the Ricatti-Bessel function chi
        //
        COMPLEX
        chi_(int l,COMPLEX rho)
        {
            COMPLEX result = -(y(l,rho)+rho/(2.*l+1.)*
                               ((double)(l)*y(l-1,rho)-(l+1.)*y(l+1,rho)));
            return result;
        }

        //
        // Routine which returns the rotation matrix element d^l_mm'(alpha)
        // Eq. (3.5) of B&V...
        //
        COMPLEX
        d(int l,int m,int m_,COMPLEX cosa2,COMPLEX sina2)
        // cosa2 and sina2 are cosine and sine of alpha/2
        {
            int kmin = (m_-m>0) ? m_-m : 0;
            int kmax = (l+m_<l-m) ? l+m_ : l-m;

            COMPLEX sum = 0;

            int k;
            for (k=kmin; k<=kmax; ++k) {
                sum+=mpow(k)*pow(cosa2,2*l+m_-m-2*k)*pow(-sina2,m-m_+2*k)/
                     (Fact(l-m-k)*Fact(l+m_-k)*Fact(k+m-m_)*Fact(k));
            }
            COMPLEX result = sqrtFact(l+m_)*sqrtFact(l-m_)*
                             sqrtFact(l+m)*sqrtFact(l-m)*sum;

            return result;

        }

        //
        // Eq. (8.16) of B&V for d^l'_m,+
        //
        COMPLEX
        dplus(int l_,int m,COMPLEX cosalpha2,COMPLEX sinalpha2)
        {
            COMPLEX result = 0.5*(d(l_,m,1,cosalpha2,sinalpha2)+
                                  d(l_,m,-1,cosalpha2,sinalpha2));
            return result;
        }

        //
        // Eq. (8.16) of B&V for d^l'_m,-
        //
        COMPLEX
        dminus(int l_,int m,COMPLEX cosalpha2,COMPLEX sinalpha2)
        {
            COMPLEX result = 0.5*(d(l_,m,1,cosalpha2,sinalpha2)-
                                  d(l_,m,-1,cosalpha2,sinalpha2));
            return result;
        }

    } // namespace BobVlieg_Supp;

    //
    // The element of the incident plane wave vector for p-polarized light
    // from Eq. (2.14) of BV&G...
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    VIp(int l,int m,int f,double thetai) const
    {
        COMPLEX result;
        COMPLEX costheta2=cos(thetai/2.);
        COMPLEX sintheta2=sin(thetai/2.);
        if (f==efield) {
            // Eq. (2.14a) of BV&G
            result =  1./k*ipow(l-1)*lvector[l]*mpow(m-1)*
                      dminus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.14b) of BV&G
            result =  -1./k*ipow(l)*lvector[l]*mpow(m-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result;
    }

    //
    // The element of the incident plane wave vector for p-polarized light
    // reflected from the surface from Eq. (2.15) of BV&G...
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    VIRp(int l,int m,int f,double thetai) const
    {
        COMPLEX result;
        COMPLEX cost = cos(thetai);
        COMPLEX phase = exp(2.*cI*qq*cost);
        COMPLEX _rp = stack->rp12(thetai,lambda,vacuum,substrate);
        COMPLEX temp = phase*_rp;
        COMPLEX costheta2=cos(thetai/2.);
        COMPLEX sintheta2=sin(thetai/2.);

        if (f==efield) {
            // Eq. (2.15a) of BV&G
            result = 1./k*phase*_rp*
                     ipow(l-1)*lvector[l]*mpow(l)*
                     dminus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.15b) of BV&G
            result =  -1./k*phase*_rp*
                      ipow(l)*lvector[l]*mpow(l-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result;
    }

    //
    // The element of the incident plane wave vector for p-polarized light
    // incident from the material...
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    VITp(int l,int m,int f,double thetai) const
    {
        COMPLEX sint = sin(thetai)*substrate.n(lambda);
        COMPLEX cost = sqrt(1.-sqr(sint));
        if (imag(cost)<0) cost = -cost;
        COMPLEX _thetai = arcsine(sint);
        if (imag(_thetai)>0) _thetai = conj(_thetai);
        COMPLEX sintheta2=sqrt(1.-cost)/sqrt(2.);
        if (real(sintheta2)<0) sintheta2 = -sintheta2;
        COMPLEX costheta2=sint/2./sintheta2;
        COMPLEX _tp = stack->tp12(_thetai,lambda,vacuum,substrate);
        COMPLEX phase = exp(cI*qq*(cost-substrate.n(lambda)*cos(thetai)));

        COMPLEX result;
        if (f==efield) {
            // Eq. (2.14a) of BV&G
            result =  _tp/k*
                      ipow(l-1)*lvector[l]*mpow(l)*
                      dminus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.14b) of BV&G
            result =  -_tp/k*
                      ipow(l)*lvector[l]*mpow(l-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result*phase;
    }

    //
    // The element of the incident plane wave vector for s-polarized light
    // from Eq. (2.16) of BV&G...
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    VIs(int l,int m,int f,double thetai) const
    {
        COMPLEX result;
        COMPLEX costheta2=cos(thetai/2.);
        COMPLEX sintheta2=sin(thetai/2.);

        if (f==efield) {
            // Eq. (2.16a) of BV&G
            result =  -cI/k*
                      ipow(l-1)*lvector[l]*mpow(m-1)*
                      dplus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.16b) of BV&G
            result = cI/k*
                     ipow(l)*lvector[l]*mpow(m-1)*
                     dminus(l,m,costheta2,sintheta2);
        }
        return result;
    }

    //
    // The element of the incident plane wave vector for s-polarized light
    // reflected from the surface from Eq. (2.17) of BV&G...
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    VIRs(int l,int m,int f,double thetai) const
    {
        COMPLEX result;
        COMPLEX cost = cos(thetai);
        COMPLEX phase = exp(2.*cI*qq*cost);
        COMPLEX _rs = stack->rs12(thetai,lambda,vacuum,substrate);
        COMPLEX costheta2=cos(thetai/2.);
        COMPLEX sintheta2=sin(thetai/2.);

        if (f==efield) {
            // Eq. (2.17a) of BV&G
            result =  -cI/k*phase*_rs*
                      ipow(l-1)*lvector[l]*
                      mpow(l-1)*dplus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.17b) of BV&G
            result =  cI/k*phase*_rs*
                      ipow(l)*lvector[l]*mpow(l)*
                      dminus(l,m,costheta2,sintheta2);
        }
        return result;
    }

    //
    // The element of the incident plane wave vector for s-polarized light
    // incident from the material...
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    VITs(int l,int m,int f,double thetai) const
    {
        COMPLEX sint = sin(thetai)*substrate.n(lambda);
        COMPLEX cost = sqrt(1.-sqr(sint));
        if (imag(cost)<0) cost = -cost;
        COMPLEX _thetai = arcsine(sint);
        if (imag(_thetai)>0) _thetai = conj(_thetai);
        COMPLEX sintheta2=sqrt(1.-cost)/sqrt(2.);
        if (real(sintheta2)<0) sintheta2 = -sintheta2;
        COMPLEX costheta2=sint/2./sintheta2;
        COMPLEX _ts = stack->ts12(_thetai,lambda,vacuum,substrate);
        COMPLEX phase = exp(cI*qq*(cost-substrate.n(lambda)*cos(thetai)));

        COMPLEX result;
        if (f==efield) {
            // Eq. (2.16a) of BV&G
            result =  -cI/k*_ts*
                      ipow(l-1)*lvector[l]*
                      mpow(l-1)*dplus(l,m,costheta2,sintheta2);
        } else {
            // Eq. (2.16b) of BV&G
            result = cI/k*_ts*
                     ipow(l)*lvector[l]*mpow(l)*
                     dminus(l,m,costheta2,sintheta2);
        }
        return result*phase;
    }

    //
    // The following should iteratively improve the solution to A.x=b
    //
    void
    Bobbert_Vlieger_BRDF_Model::
    iterative_improvement(ScatterTMatrix& Ainv, ScatterTMatrix& A, vector<COMPLEX>& b, vector<COMPLEX>& x)
    {
        vector<COMPLEX> temp1(sqrsize);
        int m,l,f,l_,f_;

        for (m=-LMAX; m<=LMAX; ++m) {
            COMPLEX temp2;
            int beginl = (m==0) ? 1 : abs(m);
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    temp1[i_]=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            temp1[i_]+=A[mm][ll_][ll]*x[i];
                        }
                    }
                    temp1[i_]-=b[i_];
                }
            }
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    temp2=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            temp2+=Ainv[mm][ll_][ll]*temp1[i];
                        }
                    }
                    x[i_]-=temp2;
                }
            }
        }
    }

    //
    // Routine to calculate the scattered wave vector for a specific incident
    // angle.  This is evaluation of Eq. (5.3) of B&V...
    //
    void
    Bobbert_Vlieger_BRDF_Model::
    calculate_W(double thetai)
    {
        int i,m,l,f,l_,f_;

        // First calculate the V vector of the incident wave...
        if (is_down_to_up()||is_down_to_down()) {
            for (l=1; l<=LMAX; ++l) {
                for (m=-l; m<=l; ++m) {
                    for (f=0; f<=1; ++f) {
                        i=index(l,m,f);
                        Vp[i]=VIp(l,m,f,thetai)+VIRp(l,m,f,thetai);
                        Vs[i]=VIs(l,m,f,thetai)+VIRs(l,m,f,thetai);
                    }
                }
            }
        } else { // is_backward()
            COMPLEX sint = sin(thetai)*substrate.n(lambda);
            COMPLEX cost = sqrt(1.-sqr(sint));
            // The following factor accounts for the backwards transmission coefficient, compared
            // to the forward transmission coefficient, and the Poynting vector across the
            // interface...
            double factor = abs(COMPLEX(cos(thetai))/cost*COMPLEX(sqrt(substrate.n(lambda))));

            for (l=1; l<=LMAX; ++l) {
                for (m=-l; m<=l; ++m) {
                    for (f=0; f<=1; ++f) {
                        i=index(l,m,f);
                        Vp[i]=VITp(l,m,f,thetai)*factor;
                        Vs[i]=VITs(l,m,f,thetai)*factor;
                    }
                }
            }
        }


        // Then multiply the incident vector V by the scattering matrix
        for (m=-LMAX; m<=LMAX; ++m) {
            int beginl = (m==0) ? 1 : abs(m);
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    Wp[i_]=0;
                    Ws[i_]=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            Wp[i_]+=ScatMatrix[mm][ll_][ll]*Vp[i];
                            Ws[i_]+=ScatMatrix[mm][ll_][ll]*Vs[i];
                        }
                    }
                }
            }
        }

        // Iteratively improve the solution, since the matrix
        // inversion isn't perfect...
        if (order<0) {
            for (i=0; i<improve; ++i) {
                iterative_improvement(ScatMatrix,ScatMatrixInverse,Vp,Wp);
                iterative_improvement(ScatMatrix,ScatMatrixInverse,Vs,Ws);
            }
        }

        // Then multiply the reflection matrix A and the vector W to get VSR...
        for (m=-LMAX; m<=LMAX; ++m) {
            int beginl = (m==0) ? 1 : abs(m);
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    VSRp[i_]=0;
                    VSRs[i_]=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            VSRp[i_]+=A[mm][ll_][ll]*Wp[i];
                            VSRs[i_]+=A[mm][ll_][ll]*Ws[i];
                        }
                    }
                }
            }
        }
    }

    //
    // The vector Z is that part of the scattering function which depends
    // upon the scattering angle thetas.  It must be evaluated anytime
    // that thetas changes...
    //
    void
    Bobbert_Vlieger_BRDF_Model::
    calculate_Z(double thetas)
    {
        if (is_down_to_up()||is_up_to_up()) {
            double d_cost = cos(thetas);
            COMPLEX cost= d_cost;
            double sint = real(sqrt(1.-sqr(cost)));

            COMPLEX _cost=cos(pi-thetas);
            COMPLEX rp = stack->rp12(pi-thetas,lambda,vacuum,substrate);
            COMPLEX rs = stack->rs12(pi-thetas,lambda,vacuum,substrate);
			COMPLEX phase = exp(2.*cI*qq*_cost);
            COMPLEX rpphase = rp*phase;
            COMPLEX rsphase = rs*phase;

            for (int i=0,l=1; l<=LMAX; ++l) {
                if (l<=BH_LMAX) {
                    for (int m=-l; m<=l; ++m) {
                        COMPLEX temp1 = Ptilde(l,m,cost)*(double)m/sint;
                        COMPLEX temp2 = P_tilde(l,m,cost);
                        for (int f=0; f<=1; ++f,++i) {
                            if (f==efield) {
                                Zp[i] = ipow(l)*mpow(l)*temp2;
                                Zs[i] = -ipow(l+1)*mpow(l+1)*temp1;
                                Zp[i] = Zp[i]*(1.+mpow(l-m+1)*rpphase);
                                Zs[i] = Zs[i]*(1.+mpow(l-m)*rsphase);
                            } else {
                                Zp[i] = -ipow(l+1)*mpow(l+1)*temp1;
                                Zs[i] = -ipow(l)*mpow(l)*temp2;
                                Zp[i] = Zp[i]*(1.+mpow(l-m)*rpphase);
                                Zs[i] = Zs[i]*(1.+mpow(l-m+1)*rsphase);
                            }
                        }
                    }
                } else {
                    for (int m=-l; m<=l; ++m) {
                        for (int f=0; f<=1; ++f,++i) {
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                        }
                    }
                }
            }
        } else { // is_down_to_down()||is_up_to_down()
            double index = substrate.n(lambda);
            COMPLEX sint = (COMPLEX)(sin(pi-thetas)*index);
            COMPLEX _thetas = arcsine(sint);
            if (imag(_thetas)>0) _thetas = conj(_thetas);
            COMPLEX cost= sqrt(1.-sqr(sint));
            if (imag(cost)<0) cost = -cost;
            COMPLEX tp = stack->tp12(_thetas,lambda,vacuum,substrate);
            COMPLEX ts = stack->ts12(_thetas,lambda,vacuum,substrate);
			COMPLEX phase = exp(cI*qq*(cost-substrate.n(lambda)*cos(pi-thetas)));
            // The following factor accounts for the transmittance across the interface and
            // the Jacobian as the solid angle across the interface changes...
            double factor = abs(COMPLEX(cos(pi-thetas))/cost*COMPLEX(sqrt(cube(index))));
            if (imag(cost)!=0.) phase = -phase;
            COMPLEX tpphase = tp*phase*factor;
            COMPLEX tsphase = ts*phase*factor;

            for (int i=0,l=1; l<=LMAX; ++l) {
                if (l<=BH_LMAX) {
                    for (int m=-l; m<=l; ++m) {
                        COMPLEX temp1 = Ptilde(l,m,cost)*(double)m/sint;
                        COMPLEX temp2 = P_tilde(l,m,cost);
                        for (int f=0; f<=1; ++f,++i) {
                            if (f==efield) {
                                Zp[i] = ipow(l)*mpow(l)*temp2*tpphase;
                                Zs[i] = -ipow(l+1)*mpow(l+1)*temp1*tsphase;
                            } else {
                                Zp[i] = -ipow(l+1)*mpow(l+1)*temp1*tpphase;
                                Zs[i] = -ipow(l)*mpow(l)*temp2*tsphase;
                            }
                        }
                    }
                } else {
                    for (int m=-l; m<=l; ++m) {
                        for (int f=0; f<=1; ++f,++i) {
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                        }
                    }
                }
            }
        }
    }

    //
    // The vector eIP is that part of the scattering function which depends
    // upon the azimuthal scattering angle phis.  It must be evaluated anytime
    // that phis changes...
    //
    // This function basically calculates exp(i m phis)...
    //
    void
    Bobbert_Vlieger_BRDF_Model::
    calculate_eIP(double phis)
    {
        vector<COMPLEX> expmphi(2*LMAX+2);

        COMPLEX expphi=exp(cI*phis);
        expmphi[LMAX]=1.;

        for (int m=1; m<=LMAX; ++m) {
            expmphi[LMAX+m]=expmphi[LMAX+m-1]*expphi;
            expmphi[LMAX-m]=expmphi[LMAX-m+1]/expphi;
        }

        for (int i=0,l=1; l<=BH_LMAX; ++l) {
            for (int m=-l; m<=l; ++m) {
                for (int f=0; f<=1; ++f,++i) {
                    eIP[i] = expmphi[m+LMAX];
                }
            }
        }
    }

    //
    // The total scattered electric field depends upon the incident vector W and
    // the scattered vector Z.
    //
    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    E(vector<COMPLEX>& W,vector<COMPLEX>& Z)
    {
        COMPLEX E(0.);

        for (int i=0,l=1; l<=BH_LMAX; ++l) {
            for (int m=-l; m<=l; ++m) {
                for (int f=0; f<=1; ++f,++i) {
                    E += W[i]*Z[i]*eIP[i];
                }
            }
        }

        return E;
    }

} // namespace SCATMECH





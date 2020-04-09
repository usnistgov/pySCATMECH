//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: bobvlieg3.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************

//
// THIS FILE CONTAINS ROUTINES USED TO CALCULATE THE INTEGRALS OF THE
// BOBBERT & VLIEGER MODEL
//

#include "scatmech.h"
#include "bobvlieg.h"

using namespace std;

namespace SCATMECH {

    static double zeros10[]=
    {   0.13779347054049243083,0.72945454950317049816,1.8083429017403160482,
        3.4014336978548995145,5.5524961400638036324,8.3301527467644967002,
        11.843785837900065565,16.2792578313781021,21.996585811980761951,
        29.92069701227389156
    };

    static double weights10[]=
    {   0.30844111576502014155,0.40111992915527355152,0.21806828761180942159,
        0.062087456098677747393,0.0095015169751811005538,
        0.00075300838858753877546,0.000028259233495995655674,
        4.2493139849626863726e-7,1.8395648239796307809e-9,
        9.9118272196090085584e-13
    };

    static double weights20[]=
    {   0.16874680185111386215,0.29125436200606828172,0.26668610286700128855,
        0.16600245326950684003,0.07482606466879237054,0.024964417309283221073,
        0.0062025508445722368474,0.001144962386476908242,
        0.00015574177302781197478,0.000015401440865224915689,
        1.0864863665179823515e-6,5.3301209095567147509e-8,
        1.7579811790505820036e-9,3.7255024025123208726e-11,
        4.7675292515781905245e-13,3.3728443433624384124e-15,
        1.155014339500398831e-17,1.5395221405823435535e-20,
        5.2864427255691578288e-24,1.6564566124990232959e-28
    };

    static double zeros20[]=
    {   0.070539889691988753367,0.37212681800161144379,
        0.91658210248327356467,1.7073065310283438807,2.7491992553094321296,
        4.0489253138508869224,5.6151749708616165141,7.4590174536710633098,
        9.5943928695810967725,12.03880254696431631,14.814293442630739979,
        17.948895520519376017,21.478788240285010976,25.451702793186905504,
        29.932554631700612007,35.013434240479000006,40.833057056728571062,
        47.61999404734650214,55.810795750063898891,66.524416525615753819
    };

    static double weights30[]=
    {   0.11604408602039325562,0.22085112475069602168,0.24139982758787346417,
        0.19463676844641672701,0.12372841596688099223,0.06367878036898826934,
        0.026860475273380519411,0.0093380708816042351496,
        0.0026806968913369005385,0.0006351291219408776464,
        0.00012390745990688661704,0.000019828788438952961056,
        2.5893509291314845837e-6,2.7409428405360851638e-7,
        2.3328311650257961682e-8,1.5807455747783781037e-9,
        8.4274791230570478545e-11,3.485161234907977146e-12,
        1.0990180597534727279e-13,2.5883126649592354134e-15,
        4.4378380598403008662e-17,5.3659183082123539536e-19,
        4.3939468922917157783e-21,2.3114097943886493589e-23,
        7.2745884982925408323e-26,1.2391497014482743994e-28,
        9.8323750831056357108e-32,2.8423235534027969144e-35,
        1.8786080317495715678e-39,8.7459804404651875591e-45
    };

    static double zeros30[]=
    {   0.047407180540804851462,0.24992391675316022399,
        0.61483345439276828461,1.1431958256661007983,1.8364545546225722915,
        2.6965218745572151958,3.7258145077795089493,4.9272937658498824097,
        6.3045155909650745228,7.8616932933702604688,9.6037759854792620798,
        11.536546597956139701,13.666744693064236295,16.002221188981066255,
        18.552134840143150124,21.327204321783128928,24.340035764532693401,
        27.605554796780961028,31.141586701111235818,34.969652008249069544,
        39.116084949067889122,43.613652908484827807,48.503986163804200427,
        53.841385406507505617,59.699121859235495477,66.180617794438489652,
        73.44123859555988224,81.736810506727685722,91.556466522536838256,
        104.15752443105889451
    };

    static double weights40[]=
    {   0.08841210619034244094,0.1768147390957222956,0.21136311701596243103,
        0.19408119531860179966,0.14643428242412511441,0.093326798435770880507,
        0.050932204361044237026,0.02397619301568484184,
        0.0097746252467144596189,0.0034579399930184868613,
        0.001062246893896871935,0.00028327168532432471583,
        0.000065509405003246292798,0.000013116069073267784125,
        2.2684528787793650545e-6,3.3796264822006792108e-7,
        4.3228213222820885689e-8,4.7284937709907793279e-9,
        4.4031741042328488129e-10,3.4724414848038224856e-11,
        2.3053815449168221616e-12,1.2797725976766356072e-13,
        5.8941771723511529447e-15,2.2322175799045774184e-16,
        6.8803364842843023409e-18,1.7056037368180867485e-19,
        3.3537119406661829355e-21,5.1461995601366791408e-23,
        6.044762511587663289e-25,5.3105847773213307528e-27,
        3.3925280532805218961e-29,1.5217354931814569975e-31,
        4.5852916145026869176e-34,8.762158657486248561e-37,
        9.8274157251479333061e-40,5.8011520191697791085e-43,
        1.5309086846066868536e-46,1.3819863056493280997e-50,
        2.5666336050123721838e-55,2.7003609402170336406e-61
    };

    static double zeros40[]=
    {   0.035700394308888385122,0.188162283158698516,0.46269428131457645356,
        0.85977296397293492226,1.3800108205273371865,2.0242091359228267334,
        2.7933693535068164577,3.6887026779082702096,4.7116411465549726936,
        5.8638508783437181143,7.1472479081022882507,8.5640170175861637627,
        10.116634048451939407,11.807892294004584843,13.640933712537087228,
        15.619285893339073837,17.746905950095663043,20.02823283457489053,
        22.468249983498418351,25.072560772426203794,27.847480009168862721,
        30.800145739445462701,33.938657084913719609,37.272245880476004328,
        40.811492823886920466,44.568603175334462707,48.557763533059992281,
        52.795611187216932969,57.301863323393627495,62.100179072775111612,
        67.219370927126998799,72.695158847612462118,78.572802911571309281,
        84.911231135704984543,91.789874671236376992,99.32080871744680825,
        107.67244063938827252,117.12230951269068881,128.20184198825565119,
        142.28004446915999789
    };

    static double weights50[] =
    {   0.071404726135189883536,0.14714860696458836553,
        0.18567162757483134625,0.18438538252735394717,0.15420116860635562195,
        0.11168536990226878894,0.0710528854901958562,0.04002027691150833111,
        0.020050623080071713675,0.0089608512036462806293,
        0.0035781124153156601908,0.0012776171567890498762,
        0.00040803024498371894306,0.00011652883223097239333,
        0.000029741704936941654528,6.7778425265420277603e-6,
        1.3774795031713603562e-6,2.4928861817200917845e-7,
        4.0103543504278268139e-8,5.7233317481414250756e-9,
        7.2294342491826646715e-10,8.0617101422817794832e-11,
        7.9133930999437228417e-12,6.8157366176767798825e-13,
        5.1324267165894902132e-14,3.3656247624378144081e-15,
        1.9134763269650353765e-16,9.385589781827253322e-18,
        3.9500699645034111309e-19,1.4177495178275120254e-20,
        4.3099702762921753006e-22,1.1012575198455480135e-23,
        2.3446177556089867146e-25,4.118544154638230108e-27,
        5.9022467635964484321e-29,6.8120089165530655136e-31,
        6.2374494988121017565e-33,4.4524405796833773756e-35,
        2.4268623522504871339e-37,9.8529714810496861343e-40,
        2.8910788723184281378e-42,5.9061627081123604884e-45,
        8.0128745975039700879e-48,6.78957542439641676e-51,
        3.3081730108492522889e-54,8.250964876440456021e-58,
        8.8487281282980176955e-62,3.0648948898444166018e-66,
        1.9887082293307516128e-71,6.0495671522387830948e-78
    };

    static double zeros50[] =
    {   0.028630518339379081948,0.15088293567693373238,
        0.37094878153489643847,0.68909069988104788096,1.1056250235399134233,
        1.6209617511025014149,2.2356103759151801252,2.9501833666418351204,
        3.7653997744057823296,4.6820893875592845835,5.7011975747848903101,
        6.8237909097945513773,8.051063669390792161,9.3843453082584069218,
        10.825109031549148043,12.374981608757462086,14.035754599829908268,
        15.809397197844667382,17.698070933350249298,19.704146535461564638,
        21.830223306578255095,24.079151444411499486,26.454057841252983897,
        28.958376011937380403,31.595880956622864625,34.370729963090445532,
        37.287510610550486553,40.351297573586065598,43.567720269995021221,
        46.943043991603038084,50.484267963129922656,54.199244880168622273,
        58.096828017248526501,62.187054175688905655,66.481373878444817064,
        70.992944826619490451,75.737011547727312131,80.73140480247768871,
        85.997211136463229075,91.559690412533876382,97.449565614850555658,
        103.70489123669226289,110.37385880764028209,117.51919820311117549,
        125.22547013347343967,133.61202792272869499,142.85832548925411939,
        153.2603719726035866,165.38564331668254038,180.69834370921451684
    };

    static double weights60[] =
    {   0.059883611523733380661,0.12591096707540107193,
        0.16473078908210712807,0.17239118732674750403,0.15442926800152203448,
        0.12180351302605122744,0.085807976879846703278,0.054435367226483746819,
        0.031252789752223371974,0.016290259047635155781,
        0.0077247456605088567732,0.0033366894943076820164,
        0.0013138658634496515358,0.00047178872610582552603,
        0.00015450073634605151289,0.000046133630291513734031,
        0.000012555541586433208077,3.1126536116671386255e-6,
        7.023892316922956353e-7,1.441394233872203205e-7,
        2.6871155382478738257e-8,4.5453240465711692205e-9,
        6.9667460870539940245e-10,9.6611190516491905786e-11,
        1.2101372902012551474e-11,1.3666508971830879161e-12,
        1.3887583568740821736e-13,1.2670440927349749065e-14,
        1.0354183850146325428e-15,7.5590720583370492079e-17,
        4.9160556832367860941e-18,2.8393315957761979668e-19,
        1.4514379644950183886e-20,6.5427350200929260331e-22,
        2.5902519096833059549e-23,8.9664278435414917268e-25,
        2.7006780092750217553e-26,7.0398894915621409495e-28,
        1.5787547853764471914e-29,3.0258776894584834948e-31,
        4.9201673552566388837e-33,6.7316641160505118187e-35,
        7.6780965388276276462e-37,7.2246694240103543426e-39,
        5.5415003983611333256e-41,3.4176601279081430909e-43,
        1.6681495220378452232e-45,6.3255732716007591197e-48,
        1.8231396385814367186e-50,3.8906596692228002314e-53,
        5.9551615457676981661e-56,6.2854492261473107477e-59,
        4.3523952400430152023e-62,1.853356484986910249e-65,
        4.4482734830374030633e-69,5.3225663149557690379e-73,
        2.641206780522461118e-77,4.0165058425505468822e-82,
        1.0516941039201472127e-87,1.0909419486248200726e-94
    };

    static double zeros60[] =
    {   0.023897977262724994782,0.12593471888169076076,
        0.30957893432678988075,0.57499554209280526609,0.92236948211666379157,
        1.3519383600081679397,1.8639963442992054802,2.458895843822428585,
        3.137049009785895931,3.8989293872049917811,4.7450738001258887621,
        5.676084508246916826,6.6926316627865745904,7.7954560890310120109,
        8.9853724256576560911,10.263272655037909547,11.630130063841871782,
        13.0870036793502451,14.635043234018346591,16.275494719209406879,
        18.009706598857114263,19.839136765434035716,21.765360334373534885,
        23.790078389494180642,25.915127811604900409,28.142492346079812161,
        30.474315093739509195,32.912912644080369003,35.460791112322411216,
        38.120664393927127573,40.895475014812933693,43.788418035940638635,
        46.802968571856479942,49.942913610317752243,53.212388982588307206,
        56.615922542696984717,60.158484884500430204,63.845549279532241991,
        67.683162987059545208,71.678032714447406955,75.837627854657063291,
        80.170306292607892439,84.685469194509277238,89.39375349025279188,
        94.307274066118867443,99.439932542889869648,104.80781680747746147,
        110.42972668651628922,116.32787889753132982,122.52887338413981661,
        129.06505218529826502,135.97646860411319811,143.31384526024606358,
        151.14321669561511669,159.55362523885103398,168.67080654892222049,
        178.6839250131463791,189.90524696213377268,202.93398795040067732,
        219.31811577379970553
    };


    static double weights70[] =
    {   0.051563243012164270229,0.10998950162574020169,
        0.14768212415343261318,0.1604596898024627034,0.15098695443069266183,
        0.1265690432180260599,0.095895700825128814593,0.066215373316473835375,
        0.041883348541682807904,0.024350102511921056885,
        0.013040962538543840992,0.0064435427862992088916,
        0.0029402504437194517251,0.0012398485702619616405,
        0.0004833241913894170636,0.00017420528057267625014,
        0.000058053494376544513075,0.000017884107821123416248,
        5.0915545259831160314e-6,1.339062360432574834e-6,
        3.2515951374757588397e-7,7.2857300502117805672e-8,
        1.5053095939896640226e-8,2.8655660408396162473e-9,
        5.0216147142098852899e-10,8.0928576281928769829e-11,
        1.1981828053238699934e-11,1.6278069543179481076e-12,
        2.0267240842502463373e-13,2.3094413655179291787e-14,
        2.4049252084517906854e-15,2.2850231074703926892e-16,
        1.97757612701751013e-17,1.5560806329624218181e-18,
        1.1110387717731329791e-19,7.1829197288354060023e-21,
        4.1952030335548165159e-22,2.2080731117912345824e-23,
        1.0445433763390937251e-24,4.4283690388074087953e-26,
        1.6773185796577407372e-27,5.6569146501250765273e-29,
        1.6925807723996771783e-30,4.4750675084987417923e-32,
        1.0409975602882768174e-33,2.1205320042553443222e-35,
        3.7629887947219422552e-37,5.7841309737942407093e-39,
        7.6529496704905794344e-41,8.6552529401353799781e-43,
        8.3028189806637512776e-45,6.697380650711144377e-47,
        4.4987478877910064173e-49,2.4889052050787134716e-51,
        1.1200025797450440813e-53,4.0410218782554424425e-56,
        1.1497889564466691143e-58,2.5303745238710871401e-61,
        4.2097318455251804551e-64,5.1515627440700920985e-67,
        4.4853453555808463381e-70,2.6665794510821039574e-73,
        1.0275572305094832968e-76,2.3985560329365738937e-80,
        3.0959719429444488907e-84,1.9433197501433467897e-88,
        4.891505804009283045e-93,3.5953048497449364636e-98,
        4.2133030593353622796e-104,1.6795619696557684288e-111
    };

    static double zeros70[] =
    {   0.020508076855477421325,0.10806702025631410371,
        0.26563794660339549732,0.49333422708016144832,0.79127365982113265697,
        1.159606217703598357,1.5985170394817861328,2.1082274245705460469,
        2.6889955806317441561,3.3411173958832110116,4.064927300484995867,
        4.860799239590724741,5.7291477704892946164,6.6704292938635368482,
        7.6851434291307363944,8.7738345446358861885,9.9370934546970471203,
        11.175559297014786833,12.489921605721043349,13.880922597365373641,
        15.349359689447448566,16.896088273746363141,18.522024769724042349,
        20.22814998675812811,22.015512827969596426,23.885234373049034234,
        25.838512382870400714,27.876626274954639259,30.000942626181155446,
        32.212921267755075169,34.514122047582961748,36.906212347210986335,
        39.390975454736450165,41.970319912113620428,44.646289975661884062,
        47.421077353134833121,50.297034410420869247,53.276689077093000734,
        56.362761724244586683,59.558184342457255409,62.866122415109115968,
        66.289999966175402748,69.833528367009889907,73.500739619754960224,
        77.296025004691687492,81.224180196844650247,85.290458239813578278,
        89.500632134864984794,93.861069292904766126,98.378820752042406786,
        103.0617289508173172,107.9185590653046012,112.95916061514784235,
        118.19466844535009009,123.63775565084377488,129.3029560944188626,
        135.20708180192250383,141.36977226484451322,147.81423126299372723,
        154.56823716079665802,161.66556400425667395,169.14804149570311696,
        177.06865045639524215,185.49638276051802978,194.52430038011315135,
        204.28387229494048068,214.97299389279797307,226.91855381008155889,
        240.74830983405439524,258.08550153049492097
    };

    static double weights80[] =
    {   0.045272641464027452755,0.097622691129370068013,
        0.13365852798017004242,0.14937645827101128701,0.14584706978660375008,
        0.12798046766667039213,0.10240364873547614395,0.075344057336652620107,
        0.051240835821065727379,0.032323138479013890286,0.01895666011098740817,
        0.010353138932155252397,0.0052716047497087076864,
        0.0025045007398191357512,0.0011108140479285224387,
        0.00046009904063076981317,0.00017800427319655150323,
        0.000064328305181779756454,0.000021714152386661386056,
        6.8452422068302032728e-6,2.0148387768652888554e-6,
        5.5356096108309897758e-7,1.4190662350227729006e-7,
        3.3928319886067029298e-8,7.5618521942413492988e-9,
        1.5702103797210433426e-9,3.0358719151260974984e-10,
        5.4614960265651831348e-11,9.1353117403076832507e-12,
        1.4196250412935834672e-12,2.0478152050918212355e-13,
        2.7395267234969947724e-14,3.3954734164377129687e-15,
        3.8950030226439436306e-16,4.130565233562671889e-17,
        4.0446805024741155105e-18,3.6523826707208946996e-19,
        3.037339021429840488e-20,2.32275748066594179e-21,
        1.6309293123767971508e-22,1.049710736450737773e-23,
        6.182186966674321993e-25,3.3253390304932169409e-26,
        1.6303421327356487754e-27,7.2700624371570635928e-29,
        2.941821059548549614e-30,1.0775663295033320128e-31,
        3.5634814454894239535e-33,1.0609024804179843587e-34,
        2.8347893331015768975e-36,6.7761153137539824132e-38,
        1.4438133638965327344e-39,2.7317419619723957224e-41,
        4.5703808506247258942e-43,6.7309359032141144705e-45,
        8.6827148453882270607e-47,9.7574467664989971107e-49,
        9.4957661365042295333e-51,7.9503372413637032776e-53,
        5.6852252212220386817e-55,3.4443689191523066463e-57,
        1.7520838300614889115e-59,7.4077279619841976241e-62,
        2.5735398256614593841e-64,7.2517056812244499503e-67,
        1.6328052016335641917e-69,2.8875122894699307129e-72,
        3.9306497186688057597e-75,4.0218963695440024828e-78,
        3.0065656681535619953e-81,1.5862629302776148181e-84,
        5.6593984430409922697e-88,1.2934555058728238527e-91,
        1.7649977712530878877e-95,1.3078621805548315574e-99,
        4.6038641952371485236e-104,6.2980356469014355246e-109,
        2.405692653032077939e-114,1.3648444753359407889e-120,
        2.2905062537183813302e-128
    };

    static double zeros80[] =
    {   0.017960423300698365554,0.094639912994353988811,
        0.23262286812586756921,0.43199254780238748026,0.69282886135202183991,
        1.0152325561894714374,1.3993276878428727741,1.8452623038358451381,
        2.3532088716092615245,2.9233646865554263248,3.5559523140461340594,
        4.2512200823098780832,5.0094426336201647724,5.8309215386087190198,
        6.7159859778513171116,7.6649934948917730607,8.6783308251677010954,
        9.7564148057429307132,10.899693371287855377,12.108646642365699901,
        13.38378811277864737,14.725665943508585539,16.134864371662466579,
        17.61200524381443786,19.157749684241247922,20.772799909792096092,
        22.457901204540458311,24.213844068958647377,26.041466560165586693,
        27.941656841859465556,29.915355964900985501,31.963560902208920711,
        34.087327864726189875,36.287775928781454459,38.566091009292210458,
        40.9235302180312672,43.361426651731230296,45.881194661278886346,
        48.484335660833189136,51.172444544607010596,53.947216789554447121,
        56.810456334636223134,59.764084342109954943,62.810148963926477204,
        65.950836257456057343,69.188482420236277374,72.525587544263345359,
        75.964831127864174827,79.509089629088836962,83.161456401053689663,
        86.925264419615623448,90.804112300940755952,94.801894215947433207,
        98.922834446940579165,103.17152750803913023,107.55298497753990633,
        112.07269048412833362,116.73666467350366632,121.55154249095262557,
        126.52466579651554034,131.66419525212031087,136.97924668693697395,
        142.48005891216160193,148.17820245500444182,154.08684228179869786,
        160.22107287009571594,166.59835193405391874,173.23907133424950383,
        180.16732304903231798,187.41194967696377239,195.00802244153299145,
        202.99898419507493782,211.43987049483646669,220.40236815173573965,
        229.98320607568000435,240.31908705584154042,251.61587933049961117,
        264.2138238831991021,278.76673304600456365,296.96651199565134576
    };

    static double weights90[] =
    {   0.040349880733194453759,0.08774516590825644438,
        0.12197433007105024567,0.13934609737241291507,0.14002407789787965059,
        0.12732257003199616194,0.10629530202232139509,0.082164435226579932903,
        0.059116519378218446025,0.039729310560359240659,
        0.025000007759737647874,0.014755008821014049611,
        0.0081779335686123848974,0.0042602882246365180656,
        0.0020873891781625752096,0.00096234803399016944956,
        0.00041759572853942397377,0.00017059029683733717916,
        0.000065609060328746861807,0.000023756674643305893932,
        8.0982321811968663321e-6,2.5984893542075113491e-6,
        7.8468797635866041887e-7,2.2295410886798467612e-7,
        5.9587194750250550937e-8,1.4974990265515540693e-8,
        3.5374876257632836445e-9,7.8516100786321949079e-10,
        1.6366693239823142499e-10,3.2024852008286152608e-11,
        5.8790260158163194098e-12,1.0119666541227608527e-12,
        1.6323085751560917997e-13,2.4656306846273059641e-14,
        3.4852928272245770461e-15,4.6069395026860506574e-16,
        5.6898954252073074319e-17,6.5607032051167745985e-18,
        7.0561038488652820431e-19,7.0719545376382126601e-20,
        6.5984438666303071938e-21,5.725504161134784268e-22,
        4.6150050612209146816e-23,3.4514766652443002241e-24,
        2.3920505589556434351e-25,1.534246492043559792e-26,
        9.0943633731774701234e-28,4.9746077325609883178e-29,
        2.507108806381208327e-30,1.1622307736554326377e-31,
        4.9470799412827950519e-33,1.9298670030029676829e-34,
        6.8858947597425135659e-36,2.2424583737455556441e-37,
        6.6502020909542701992e-39,1.7915986619289898062e-40,
        4.3733929167498123721e-42,9.646398032256621866e-44,
        1.9168609218601826179e-45,3.4206494390742623471e-47,
        5.4629440099275337349e-49,7.7792042047128157282e-51,
        9.8376857185584135865e-53,1.1000512960261402582e-54,
        1.0825401124380840695e-56,9.327176410932607406e-59,
        6.9965900818267877287e-61,4.5411901751331572635e-63,
        2.5330558620365779701e-65,1.2051559929332908266e-67,
        4.8499140408727960365e-70,1.635546328896279327e-72,
        4.5738095372088182927e-75,1.0481860876578609484e-77,
        1.9422117925967768062e-80,2.8651770738024518238e-83,
        3.3056452572549452321e-86,2.9212348809053524052e-89,
        1.9293103837423878478e-92,9.2466473586023640843e-96,
        3.1030190283934115114e-99,6.9755934331239771853e-103,
        9.9321849493059713348e-107,8.3275601451342939583e-111,
        3.726831966377303718e-115,7.7504028259845558834e-120,
        6.0819838356019673209e-125,1.2781040029433570607e-130,
        3.7348584764408337619e-137,2.8365669716227786138e-145
    };

    static double zeros90[] =
    {   0.015975805570894463949,0.084180854862784609057,
        0.20690855381710546619,0.384223173487157175,0.6161816061342283208,
        0.90285482691855056406,1.2443300203706891562,1.640711100508238261,
        2.0921189850921238893,2.5986918332633257869,3.1605852935990622287,
        3.777972776540531172,4.4510457566246297243,5.180014107417716422,
        5.9651064712612580077,6.8065706657395314121,7.7046741288063573142,
        8.6597044046355533803,9.6719696724416088845,10.741799320735060605,
        11.869544569724794885,13.055579144855660791,14.300300004775465633,
        15.604128127363143753,16.967509357823105897,18.390915323263748485,
        19.874844418635857179,21.419822869415019443,23.026405876977915634,
        24.695178853253273935,26.426758751933337913,28.221795504321295022,
        30.080973568776281711,32.005013603714283144,33.994674275246735839,
        36.050754211807881456,38.174094119559089823,40.365579073989486211,
        42.626141004987886394,44.956761394777370139,47.358474210523538424,
        49.83236909620132654,52.379594851493598977,55.001363229169831744,
        57.698953086641581708,60.473714932317374961,63.327075913109049981,
        66.26054529612723316,69.275720505431591365,72.374293783899453749,
        75.558059561124090547,78.828922621096812465,82.188907178693730494,
        85.640166992213089865,89.184996661062968203,92.825844384027322991,
        96.565325685413664176,100.40624045517775367,104.35159009659028073,
        108.40459863343635551,112.56873610109279215,116.84774543601675053,
        121.24567339137983417,125.76690624971792181,130.41621128588601157,
        135.19878516805646921,140.12031078855197512,145.18702441447198724,
        150.40579557497050415,155.78422280722739213,161.33074933860639362,
        167.05480409443988088,172.96697524883302836,179.07922612392762777,
        185.40516697258039774,191.96040166011915458,198.76297649278026387,
        205.83397111132544527,213.19829141957099734,220.88575726600725318,
        228.93263306367919776,237.38384754538360584,246.29633092916886249,
        255.74425697622637107,265.82774031615169948,276.68831753604922606,
        288.53922079113514634,301.73302167704964291,316.94758454351505233,
        335.93816518930042563
    };

    static double zeros100[] =
    {   0.01438614699541966946443603,0.07580361202335712464299317,
        0.18631410205718717371146,0.3459691809914290908053783,
        0.5548109375809155095983409,0.8128912841156688450380813,
        1.120273835007540148567503,1.477034329923827069718562,
        1.883260826342394705801825,2.339053849646034171844423,
        2.844526542755359066517119,3.399804827445711944287159,
        4.0050275817586520174649,4.660346835568908459504028,
        5.365927985585117014883531,6.121950030804019789061374,
        6.928605829376173055020579,7.786102377862517434336654,
        8.694661113922167989204456,9.654518243555080727805574,
        10.6659250941216755505345,11.72914849447222507161972,
        12.84447118364103057156564,14.01219224969427464896413,
        15.23262760046669784302042,16.50611046808198959418502,
        17.83299194932638742984255,19.2136415841360666849015,
        20.64844797466834951965445,22.13781944765670440655272,
        23.68218476300237610017189,25.28199387183404153877519,
        26.93771872757426494355442,28.64985415389129217109679,
        30.41891877379094291959571,32.24545600452066584227049,
        34.13003512342151647128728,36.07325241037997322620358,
        38.07573237310709319104534,40.13812906211554576323368,
        42.26112748298479417841149,44.44544511431180607747884,
        46.69183354065153941911601,49.0010802107724369560282,
        51.37401033270399452135985,53.81148891835566043829179,
        56.31442299196170823768638,58.88376397828209029370608,
        61.52051028839614264679926,64.2257101231015601669677,
        67.00046451641931159676425,69.84593064455837428644683,
        72.76332542897458625633852,75.75392946593993619270919,
        78.8190913194114721995345,81.96023221906012400440603,
        85.17885121121890019965819,88.47653081739461122166173,
        91.85494326304930893124588,95.31585734883172060708188,
        98.86114604761352771688585,102.4927949239165609633515,
        106.2129114880468162055319,110.023735616030916960441,
        113.927651188971623711111,117.9271991325732268356693,
        122.0250920704416210770071,126.2242308447503875748128,
        130.5277232067994132655524,134.9389050402274075771759,
        139.4613645542401567648374,144.0989699772127241580661,
        148.8559013977582494611375,153.7366875479730306110451,
        158.7462485117131044339031,163.8899455825872328315617,
        169.1736398100030245805137,174.6037611823766267145259,
        180.1873909402456961946854,185.932360239666971431614,
        191.8473693722483291246716,197.9421331021432574060554,
        204.227559567030507720801,210.7159728615769433706683,
        217.4213932720014810961163,224.3598947888746081453609,
        231.550068025172484218833,239.0136297513149244679889,
        246.7762409672484904573225,254.8686292570474302786384,
        263.3281684691578931098124,272.2011700240925368258755,
        281.546328283897388799148,291.4401336163771072606972,
        301.9858552516391536657452,313.3295340040755243411836,
        325.6912634370265200020353,339.435101923449616535205,
        355.2613118885341324724827,374.984112834342678704884
    };

    static double weights100[] =
    {   0.036392605883401356537,0.07967674621295139855,
        0.11211510334248694468,0.13035661297514618374,0.1340433397284623804,
        0.12540709078066374997,0.10831411209726027355,0.087096638469959342035,
        0.065551009312310614325,0.046340133582644259874,
        0.030846308627681445601,0.01936782811397878911,0.011485442360179691617,
        0.0064389510016104295905,0.003414979989692662593,
        0.0017143197401822081619,0.00081487159158783785329,
        0.00036685483659948815078,0.00015645207417810679472,
        0.000063210870528885807631,0.000024195752265189292644,
        8.7743097637554872419e-6,3.0142674860009475187e-6,
        9.8083358993452597679e-7,3.022638743532254308e-7,
        8.8200583952959386115e-8,2.4364258562006733674e-8,
        6.3697113739017570248e-9,1.575600320459668008e-9,
        3.6863292013461326689e-10,8.1547989242246112119e-11,
        1.7050625568260654077e-11,3.3682141708667259378e-12,
        6.283524955366657366e-13,1.1064980159833061707e-13,
        1.8383501754549013832e-14,2.8801150572305857339e-15,
        4.2526128973372414067e-16,5.9144324672522439838e-17,
        7.7430891025001089332e-18,9.5362494490747949743e-19,
        1.1040945169869171454e-19,1.2008469704304970569e-20,
        1.2260070074904774033e-21,1.1740217146467280969e-22,
        1.0535940646087311407e-23,8.8532317684742975097e-25,
        6.9591987754834461444e-26,5.1123698154658325241e-27,
        3.5062721488171387421e-28,2.2426462726939223396e-29,
        1.3362094234457496265e-30,7.4074289575757510893e-32,
        3.8158601371363446901e-33,1.8242001976828543039e-34,
        8.0816323360641805778e-36,3.3130638018623846737e-37,
        1.254832997687677994e-38,4.3837909801617165881e-40,
        1.410139229217347485e-41,4.1688639996685860322e-43,
        1.1304819044771411726e-44,2.8060374967871239466e-46,
        6.3612705709821942695e-48,1.3139811959781289076e-49,
        2.4668106921703939992e-51,4.1977470228717237504e-53,
        6.4562691049120031533e-55,8.9473372486217303203e-57,
        1.1135669746997689184e-58,1.2402350422552966639e-60,
        1.2313765726083698151e-62,1.0853674489051309034e-64,
        8.4549564909552701838e-67,5.7926307566097949606e-69,
        3.4718547763506601612e-71,1.8098595335622393581e-73,
        8.1537526084108297035e-76,3.152472548793288043e-78,
        1.0379005899152634183e-80,2.8848755233073495343e-83,
        6.7048012259685912147e-86,1.2889637464360916588e-88,
        2.0248483451529085824e-91,2.5633922199660732721e-94,
        2.5739615075159617883e-97,2.0126547874800776647e-100,
        1.1994844193823276437e-103,5.3121327315398034324e-107,
        1.6959692594467450795e-110,3.7620972986344511541e-114,
        5.5396417544496093738e-118,5.1106404770917605505e-122,
        2.7399654694003410753e-126,7.7136114926382004229e-131,
        9.8824946009588260463e-136,4.6468630072942033152e-141,
        5.6260372950198530067e-147,8.9050314058891380744e-154,
        3.2465651634358090752e-162
    };

    double* Gauss_Laguerre_Integration::zeros[] =
    {   zeros10,zeros20,zeros30,zeros40,zeros50,
        zeros60,zeros70,zeros80,zeros90,zeros100
    };

    double* Gauss_Laguerre_Integration::weights[] =
    {   weights10,weights20,weights30,weights40,weights50,
        weights60,weights70,weights80,weights90,weights100
    };

} // namespace SCATMECH





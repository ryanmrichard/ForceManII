//This file is autogenerated from MakeTest.py

//If you value your sanity do not try to manually edit it!!!

#include "peptide.hpp"
#include<ForceManII/FManII.hpp>
const FManII::IVector peptide_FF_types={
230,
236,
177,
178,
233,
233,
233,
85,
80,
85,
85,
85,
180,
165,
177,
178,
183,
85,
85,
180,
166,
177,
178,
183,
85,
82,
80,
80,
85,
85,
85,
85,
85,
85,
85,
180,
166,
177,
178,
183,
85,
81,
82,
80,
80,
85,
85,
85,
85,
85,
85,
85,
85,
85,
180,
166,
177,
178,
183,
85,
82,
80,
81,
80,
85,
85,
85,
85,
85,
85,
85,
85,
85,
181,
188,
177,
178,
85,
81,
81,
187,
85,
85,
85,
85,
85,
85,
180,
166,
177,
178,
183,
85,
99,
96,
85,
85,
97,
180,
166,
177,
178,
183,
85,
100,
96,
80,
85,
97,
85,
85,
85,
180,
166,
177,
178,
183,
85,
148,
142,
85,
85,
146,
180,
166,
177,
178,
183,
85,
94,
90,
90,
90,
90,
90,
90,
85,
85,
91,
91,
91,
91,
91,
180,
166,
177,
178,
183,
85,
94,
90,
90,
90,
90,
90,
108,
109,
85,
85,
91,
91,
91,
91,
110,
180,
166,
177,
178,
183,
85,
446,
451,
453,
451,
450,
453,
85,
85,
454,
91,
91,
454,
180,
166,
177,
178,
183,
85,
81,
441,
455,
442,
444,
443,
90,
90,
90,
90,
85,
85,
91,
445,
91,
91,
91,
91,
180,
166,
177,
178,
183,
85,
216,
213,
214,
214,
85,
85,
180,
166,
177,
178,
183,
85,
81,
177,
178,
179,
85,
85,
182,
182,
180,
166,
177,
178,
183,
85,
81,
216,
213,
214,
214,
85,
85,
85,
85,
180,
166,
177,
178,
183,
85,
81,
81,
177,
178,
179,
85,
85,
85,
85,
182,
182,
180,
166,
177,
178,
183,
85,
81,
152,
144,
151,
85,
85,
85,
85,
85,
85,
85,
180,
166,
177,
178,
183,
85,
81,
81,
81,
235,
230,
85,
85,
85,
85,
85,
85,
85,
85,
233,
233,
233,
180,
225,
213,
214,
183,
85,
81,
251,
250,
246,
245,
243,
243,
85,
85,
85,
85,
85,
85,
247,
244,
244,
244,
244,
214,
};
const std::map<FManII::FFTerm_t,double> peptide_egys={
{{FManII::Model_t::ELECTROSTATICS,FManII::IntCoord_t::PAIR},-1.8260189166827088},
{{FManII::Model_t::LENNARD_JONES,FManII::IntCoord_t::PAIR14},0.2546045510698036},
{{FManII::Model_t::HARMONICOSCILLATOR,FManII::IntCoord_t::BOND},0.1148760114586295},
{{FManII::Model_t::FOURIERSERIES,FManII::IntCoord_t::IMPTORSION},0.0013598198338320242},
{{FManII::Model_t::FOURIERSERIES,FManII::IntCoord_t::TORSION},0.10448095136711852},
{{FManII::Model_t::HARMONICOSCILLATOR,FManII::IntCoord_t::ANGLE},0.05160271651620947},
{{FManII::Model_t::ELECTROSTATICS,FManII::IntCoord_t::PAIR14},0.7656098010293386},
{{FManII::Model_t::LENNARD_JONES,FManII::IntCoord_t::PAIR},145.1203020639047},
};
const FManII::ConnData peptide_conns={
{1,4,5,6,},
{0,2,7,8,},
{1,3,12,},
{2,},
{0,},
{0,},
{0,},
{1,},
{1,9,10,11,},
{8,},
{8,},
{8,},
{2,13,16,},
{12,14,17,18,},
{13,15,19,},
{14,},
{12,},
{13,},
{13,},
{14,20,23,},
{19,21,24,25,},
{20,22,35,},
{21,},
{19,},
{20,},
{20,26,27,28,},
{25,29,30,31,},
{25,32,33,34,},
{25,},
{26,},
{26,},
{26,},
{27,},
{27,},
{27,},
{21,36,39,},
{35,37,40,41,},
{36,38,54,},
{37,},
{35,},
{36,},
{36,42,45,46,},
{41,43,44,47,},
{42,48,49,50,},
{42,51,52,53,},
{41,},
{41,},
{42,},
{43,},
{43,},
{43,},
{44,},
{44,},
{44,},
{37,55,58,},
{54,56,59,60,},
{55,57,73,},
{56,},
{54,},
{55,},
{55,61,62,64,},
{60,63,65,66,},
{60,67,68,69,},
{61,70,71,72,},
{60,},
{61,},
{61,},
{62,},
{62,},
{62,},
{63,},
{63,},
{63,},
{56,74,80,},
{73,75,77,78,},
{74,76,87,},
{75,},
{74,},
{74,79,81,82,},
{78,80,83,84,},
{73,79,85,86,},
{78,},
{78,},
{79,},
{79,},
{80,},
{80,},
{75,88,91,},
{87,89,92,93,},
{88,90,98,},
{89,},
{87,},
{88,},
{88,94,95,96,},
{93,97,},
{93,},
{93,},
{94,},
{89,99,102,},
{98,100,103,104,},
{99,101,112,},
{100,},
{98,},
{99,},
{99,105,106,107,},
{104,108,},
{104,109,110,111,},
{104,},
{105,},
{106,},
{106,},
{106,},
{100,113,116,},
{112,114,117,118,},
{113,115,123,},
{114,},
{112,},
{113,},
{113,119,120,121,},
{118,122,},
{118,},
{118,},
{119,},
{114,124,127,},
{123,125,128,129,},
{124,126,143,},
{125,},
{123,},
{124,},
{124,130,136,137,},
{129,131,132,},
{130,133,138,},
{130,134,139,},
{131,135,140,},
{132,135,141,},
{133,134,142,},
{129,},
{129,},
{131,},
{132,},
{133,},
{134,},
{135,},
{125,144,147,},
{143,145,148,149,},
{144,146,164,},
{145,},
{143,},
{144,},
{144,150,157,158,},
{149,151,152,},
{150,153,159,},
{150,154,160,},
{151,155,161,},
{152,155,162,},
{153,154,156,},
{155,163,},
{149,},
{149,},
{151,},
{152,},
{153,},
{154,},
{156,},
{145,165,168,},
{164,166,169,170,},
{165,167,182,},
{166,},
{164,},
{165,},
{165,171,176,177,},
{170,172,173,},
{171,174,178,},
{171,175,179,},
{172,175,180,},
{173,174,181,},
{170,},
{170,},
{172,},
{173,},
{174,},
{175,},
{166,183,186,},
{182,184,187,188,},
{183,185,206,},
{184,},
{182,},
{183,},
{183,189,198,199,},
{188,190,191,},
{189,192,200,},
{189,193,194,},
{190,193,201,},
{191,192,195,},
{191,196,202,},
{193,197,203,},
{194,197,204,},
{195,196,205,},
{188,},
{188,},
{190,},
{192,},
{194,},
{195,},
{196,},
{197,},
{184,207,210,},
{206,208,211,212,},
{207,209,218,},
{208,},
{206,},
{207,},
{207,213,216,217,},
{212,214,215,},
{213,},
{213,},
{212,},
{212,},
{208,219,222,},
{218,220,223,224,},
{219,221,232,},
{220,},
{218,},
{219,},
{219,225,228,229,},
{224,226,227,},
{225,},
{225,230,231,},
{224,},
{224,},
{227,},
{227,},
{220,233,236,},
{232,234,237,238,},
{233,235,247,},
{234,},
{232,},
{233,},
{233,239,243,244,},
{238,240,245,246,},
{239,241,242,},
{240,},
{240,},
{238,},
{238,},
{239,},
{239,},
{234,248,251,},
{247,249,252,253,},
{248,250,264,},
{249,},
{247,},
{248,},
{248,254,258,259,},
{253,255,260,261,},
{254,256,257,},
{255,},
{255,262,263,},
{253,},
{253,},
{254,},
{254,},
{257,},
{257,},
{249,265,268,},
{264,266,269,270,},
{265,267,281,},
{266,},
{264,},
{265,},
{265,271,274,275,},
{270,272,276,277,},
{271,273,},
{272,278,279,280,},
{270,},
{270,},
{271,},
{271,},
{273,},
{273,},
{273,},
{266,282,285,},
{281,283,286,287,},
{282,284,303,},
{283,},
{281,},
{282,},
{282,288,292,293,},
{287,289,294,295,},
{288,290,296,297,},
{289,291,298,299,},
{290,300,301,302,},
{287,},
{287,},
{288,},
{288,},
{289,},
{289,},
{290,},
{290,},
{291,},
{291,},
{291,},
{283,304,307,},
{303,305,308,309,},
{304,306,327,},
{305,},
{303,},
{304,},
{304,310,316,317,},
{309,311,318,319,},
{310,312,320,321,},
{311,313,322,},
{312,314,315,},
{313,323,324,},
{313,325,326,},
{309,},
{309,},
{310,},
{310,},
{311,},
{311,},
{312,},
{314,},
{314,},
{315,},
{315,},
{305,},
};
const FManII::Vector peptide={
0.0,0.0,0.0,
0.0,0.0,2.779786929819,
2.681515509213033,0.0,3.8198788860826776,
3.222311623567086,1.287610156028886,5.675540674404963,
0.899573821447626,1.558107423820335,-0.637111668561405,
-1.799145753169263,-0.002216648585097,-0.637111668561405,
0.897653859842802,-1.559214803249889,-0.637111668561405,
-1.004179603568721,1.661819365548633,3.467361841190661,
-1.520957300602584,-2.223400576055697,3.824661782560836,
-1.510807582315665,-2.208564337316058,5.884383742079297,
-0.7134055759153021,-4.013728867760286,3.203883015722358,
-3.4822453230800248,-2.1196546192595966,3.203883015722358,
4.359558172377231,-1.449973523277777,2.617302620106813,
6.971613023412591,-1.6125542090154028,3.4226320270310313,
8.202006055946535,0.98500077450636,3.317041697669667,
9.662085833813483,1.662758559365166,4.991166958858668,
3.8121158917198654,-2.47064665527849,1.100275949561349,
7.053198163535688,-2.334168754626921,5.350167653618943,
7.992376862260775,-2.9075645320172128,2.188221437044473,
7.607205352278839,2.437206064177146,1.341756474791703,
8.652700035145068,4.948043411789688,1.024826749724535,
7.931518236785031,6.6083869013409124,3.259881265954395,
9.470165262370644,8.09397986978535,4.1651223756430324,
6.395324075533139,1.784914226746104,0.01936969138725,
10.704647661944783,4.826732451925833,0.892243574336295,
7.7173007883979805,6.085985379763785,-1.46087157305634,
8.746335957900017,8.74542510997332,-1.9267513903024769,
8.119784078165145,4.323904712142768,-3.7151257053344398,
5.6588902987158844,6.078229944304929,-1.536150697554144,
8.042426255079441,9.496275606360657,-3.710992874596497,
10.803884732531131,8.72388223369872,-2.020634867161986,
8.155445097303563,10.027821511998555,-0.427055396802132,
7.432770415412217,5.193995028710016,-5.4511376029372025,
7.096659971556699,2.556064159707279,-3.4491581107386238,
10.121773018993668,3.910281487670448,-3.9675930974648397,
5.5837642420231886,6.351554242175922,4.146744790400007,
4.668620968152192,7.823588426649285,6.266516572670922,
6.275558911528287,7.307542053573165,8.595411102062544,
6.859616522948517,9.039961582070838,10.027710018165203,
4.396033663416909,5.100406349104791,3.330182852197173,
4.756735111567284,9.828100151673123,5.800736910902202,
1.89116218075164,7.198314001682976,6.759810644839482,
0.7905139551704581,8.685626620777404,8.979057603171357,
-1.986942942504105,8.060352195811094,9.472353565065905,
2.415347603662383,8.499951593728209,11.361077999291737,
1.6400100379095839,5.173750394189859,7.044140706870412,
0.766193181692028,7.489125823856175,5.059030798575636,
1.041666098012514,10.71019211799651,8.694727541140429,
-2.707650419640903,9.152330358554746,11.063247734795391,
-2.214586894043001,6.063711179531442,9.924440272318332,
-3.138566526638562,8.51200615581204,7.825370438256912,
1.575916201540671,9.580073615616872,12.901038384439703,
4.300693207819881,9.254868780443873,11.017055272720276,
2.578079467479129,6.542293734875637,11.980649231223243,
6.988393790194911,4.916539789827069,8.96879072187312,
8.520880309056384,4.181974952012934,11.115814242631403,
11.03318408403639,5.581167758514303,11.061894690987268,
11.967646026692934,6.394830857050011,13.026228951758975,
6.454249511322138,3.566124590553768,7.730279426490432,
7.521143451287801,4.630879360699884,12.859885821577251,
8.92303856595344,1.32080508275166,11.092062276675662,
10.529972729877558,0.40222061785469404,13.310742317210838,
9.880819727314245,0.3505441709595,8.54454723734262,
10.932130986774615,-2.458947361680591,13.286990351255097,
7.116736420463835,0.332056981609113,11.04276499480062,
9.712403567722587,1.007707721990184,15.1017622486274,
12.336274875367165,1.39097060872323,13.360039599085882,
8.693406622674118,1.0407155658400469,7.009458565794294,
9.846838674580047,-1.70825560227633,8.490100452147551,
11.803763319748887,1.0004096002206662,8.194201487611966,
12.082868289420231,-3.045792318934596,14.891422517969778,
11.904279734829785,-3.0433016600810943,11.567619380711468,
9.12977836860202,-3.448566626504067,13.40932932205697,
12.123601333113125,5.8582563905553835,8.80387433481309,
14.516810796814374,7.154096303266365,8.501950923824571,
14.294138714352536,9.869008939882995,9.424981913973653,
16.025726763415047,10.848716591894165,10.62317023620306,
15.969764417976805,6.1762386929984245,9.586154753849476,
15.332168979562217,7.060912025022786,5.731555932170901,
13.441592613867165,5.654285036480691,4.059522597651723,
11.167061721727094,4.661816064865803,5.539423711417293,
17.21304160230375,6.240332529367338,5.553358550860179,
15.672000293890072,8.959404452885758,5.008354016728624,
14.369661613502922,4.1248031819417275,3.038664272504088,
12.828622194815235,6.843873215734158,2.4936616280985207,
11.14252551948592,2.602162135482945,5.5419767312284325,
9.43878447259731,5.60830989289431,4.939431930457815,
12.154288593448497,11.08127761045644,8.862945279503242,
11.713075369536776,13.67129447700808,9.634278515885349,
11.814085003100804,13.879708576882914,12.501102445305712,
12.837054152178153,15.703554938454397,13.511546490527968,
10.791541042370982,10.18443444306292,7.8722734373038925,
13.151689749894675,14.890201754980884,8.805163127937588,
9.163133922033857,14.57012563614203,8.61521920987924,
8.741980139495373,17.0904928719171,9.370211985634468,
7.63725199550394,13.377905068585887,9.317411151775818,
9.117933566102966,14.451526433072392,6.5593315043924605,
7.161711899002035,17.618091139963976,8.727215710343316,
10.748271986400848,12.003633968427451,13.80832228792245,
10.720069715741014,11.981665903805325,16.546302809686804,
13.410599217921575,12.000733239034334,17.562673035610565,
13.954881876725334,13.27919712127646,19.423631725880018,
9.932900575571084,10.54866023018278,12.880385569105924,
9.7124356930644,13.64015368243535,17.236863488121084,
9.269228808234285,9.666721331664501,17.486868567287818,
9.246712723075351,9.65323435728101,20.151253710410565,
10.384254069139779,7.183001551994109,16.51917390338872,
7.333510668114003,9.685144270331262,16.78299476926103,
8.344893907152793,8.214022823510586,20.703014014952796,
9.305432178731547,5.576607849702846,17.225166084249175,
12.331057341911535,6.962323130450667,17.15481536513068,
10.337780037892303,7.117044445800042,14.460954276031513,
15.092400546077895,10.577319814531952,16.333901387975178,
17.71302341474738,10.43495163797267,17.11485143991531,
18.9161911649578,13.045606202050182,17.01898564049334,
20.38420401252657,13.725015607423373,18.685485857316724,
14.541792974388953,9.56305096390796,14.813726433932043,
17.819025704374347,9.699133911827829,19.035818378965455,
19.133568235198464,8.572199597739656,15.423465421542783,
22.394692940853616,8.43813676690203,16.446330635690735,
19.054702410773537,9.170834555439033,13.45415173279006,
18.296903391924644,6.693936686588931,15.54554738961015,
23.143124587975034,6.7529981826491365,14.721762918677358,
18.289231104409307,14.507120281942782,15.060585119327092,
19.306510288533776,17.030839891622335,14.75433990500173,
18.588550472984984,18.66623143921481,17.00874521535895,
20.120145041179647,20.160228019132326,17.912115496318474,
17.072266464753195,13.85293682879675,13.743814601561946,
21.358316185884313,16.93139684090319,14.602708291644369,
18.337609979453696,18.178657567614948,12.286076194195369,
19.290067782703485,20.820891342694637,11.782275245527968,
21.539546577407382,21.182155929189737,10.437425289292271,
17.923661722111266,22.90937682149167,12.660026511254623,
22.42261742179307,23.63190977393384,9.97032659878323,
18.80673445622294,25.359130666235774,12.192925931019591,
21.05621325092684,25.72039525273087,10.848075974783894,
16.278400135678254,18.169173032876156,12.237682201343066,
18.80010907690551,16.96461444433783,10.68764768694774,
22.593630064893613,19.571037911925973,9.760304342091759,
16.18835091285042,22.630685701811917,13.697482299763644,
24.157930120779906,23.910599003887608,8.932870810274208,
17.75265096873671,26.970248683499538,12.870048767946093,
21.737440572701452,27.61020438967439,10.487743023201375,
16.251488547868902,18.378654827934774,17.914214981892254,
15.34045542802398,19.824688272525485,20.0535662151352,
16.973174903068003,19.306755952912344,22.364039695586047,
17.552466625543975,21.03374629864362,23.804806473845417,
15.069207939014897,17.121951140455984,17.09827854299178,
15.404218562344818,21.83357471321386,19.602913809908426,
12.573940043825727,19.16742078656541,20.56636792868821,
11.4351269121707,20.57793872091885,22.76998496336506,
10.246596979471073,22.908624710548047,22.376825581627624,
11.567806463584379,19.555014925265244,25.206234826191775,
9.190744708459137,24.216390683975604,24.419914172990907,
10.511954192572444,20.8627808986928,27.249323417555058,
9.323422370146826,23.193466888321996,26.85616403581762,
8.294719792418853,24.467603298049312,28.846716354960716,
12.348900134699672,17.14068399418494,20.85701912415234,
11.427490529449152,19.479482577484912,18.883795592328376,
10.144243750728867,23.697738378760658,20.497432614335473,
12.484673167557345,17.75705592001707,25.509530178248305,
8.273876114760183,26.014349689223778,24.116618820934374,
10.614305531588661,20.073667230480194,29.128716384847205,
8.490206276802935,23.609473057266445,30.400566226127832,
17.713359785973424,16.920258795924024,22.712172576183587,
19.27221743101139,16.184495871832866,24.83971491172127,
21.769713313141594,17.60948721753606,24.774721565781594,
22.713370662460612,18.417188230576475,26.737120747140565,
17.181851674855306,15.57423154039718,21.467765335715217,
18.283565705894283,16.609565166620555,26.596035694527814,
19.70295968321806,13.327829219603377,24.78976189492804,
21.306736115014562,12.363005049481568,26.928644476125715,
20.25376378640986,11.602089872476817,29.192014646574748,
23.85831595494994,12.529590064589886,27.33388865556281,
22.077133937032112,11.307921896947155,30.932626137234738,
24.249128297009044,11.862936309642444,29.82656384669507,
20.529478587656936,12.747418779341917,22.994552440993864,
17.906204433424882,12.320837703763022,24.769110969320252,
18.36333670906794,11.686144884467538,29.441142892882585,
25.305558824173612,13.53941427968775,26.30869151954443,
21.97571423292847,9.766306111832822,32.26620387794605,
25.645743517261415,12.692544916073333,30.82858916263633,
22.837147714740112,17.91529590515796,22.509438992631686,
25.214394663012275,19.23767701938445,22.196693120904165,
24.972598553541758,21.942880243869595,23.143118918797068,
26.704819660810585,22.930696710099568,24.333712432276684,
21.978509137666205,17.211928553696232,20.956745633769838,
26.686719865287944,18.26615377063164,23.26026303285518,
26.005966864736575,19.17446190560042,19.418565370503504,
28.43992260021264,20.474633070278188,18.812894962947084,
30.74021234730877,19.36299608807097,18.835059559072068,
28.853402205509795,23.105108970254324,18.0978547740513,
32.53500039234739,21.140829511536293,18.180354541553072,
31.454944510868348,23.46769636494172,17.71141525793074,
27.30754209582417,25.23029670720971,17.73285041982397,
32.48321056189286,25.799971614127664,17.002577149730854,
28.382085646593307,27.54164135134149,17.024095459567594,
30.901131963902063,27.890966649038084,16.65035679181911,
24.50916026227341,19.923131658196453,18.217795682573126,
26.070221327814554,17.23646097648543,18.72372257297815,
30.93153576533909,17.40288536488468,18.299658612416607,
34.19114680047706,20.586278081708308,17.410618123631668,
25.28948481491118,25.082021247208814,17.9988331322277,
34.49690446549725,25.986102065140205,16.727472729978224,
27.1396777362213,29.138939802447695,16.75886674783348,
31.56097136693317,29.501886245073003,15.585071008670086,
22.81572578092881,23.137736979180396,22.609592580322698,
22.355314830694837,25.717025871292485,23.405273046347084,
22.479724941180653,25.9040029195481,26.2726544449342,
23.493215224697174,27.73016797490809,27.288448304431316,
21.45332215998129,22.23489582149779,21.623907725008323,
23.774138110769936,26.95692233582511,22.57317945024066,
19.787479567802073,26.597887625997032,22.41588999727223,
19.24191567477778,29.305364190846948,23.218537883018058,
17.187656913093512,30.23811215203147,22.518687311621854,
20.889340997468196,30.39161837356992,24.516932705266168,
18.28387940040846,25.36370002367718,23.093201806797634,
19.71273145630718,26.46130390069008,20.361980834895935,
21.444546272488374,24.007033610832266,27.574326052047233,
21.440933116397407,23.963346925418563,30.312189410800265,
24.140068430731873,24.001755606144986,31.304911605697686,
24.6879812531605,25.271067545422383,33.17106326298491,
20.6356528423489,22.551169811646776,26.642124221673534,
20.422822452837774,25.606066830396372,31.024720043582665,
20.021918973997415,21.626506098243194,31.247154020569887,
19.981501514544682,21.520668214777274,34.121083319708895,
18.92137279361964,19.76934764693762,35.217842489204706,
21.020017439333536,23.471318969662633,35.33800260594126,
20.883019864034992,19.907577323580995,30.50784974997534,
18.090533975569908,21.612684642359646,30.53137305908641,
21.777935510919733,24.91418901294776,34.34476073639687,
21.230941095321757,23.403128207349567,37.233709577696445,
25.82512953621922,22.60514416534459,30.05012410037174,
28.453899013627286,22.48326439795805,30.80678172472327,
29.629796128090444,25.106663274105376,30.72098438537069,
31.105660787143513,25.787866029442128,32.379800976322805,
25.27126596665123,21.597273926919353,28.526880231322465,
28.584399710975646,21.733543957712133,32.72085305363956,
29.878034313457466,20.64832156537916,29.088236564436837,
32.675325775112576,20.420295889190484,29.775078262124765,
34.024671389476104,18.580001352610747,28.024394531943408,
36.33081361796016,18.235438164406432,28.40224335371797,
32.72475344808086,17.552230409795346,26.34103795929814,
29.65486334333453,21.180509977853315,27.110927445024654,
28.974291756478102,18.79737086394745,29.09149634176786,
32.89764069936248,19.801051579855073,31.72697434643688,
33.59954162345676,22.258874554572213,29.68433739958496,
28.970739071618784,26.577026940804462,28.779821944676105,
29.959770631659683,29.113260315983172,28.4846127301265,
29.24545609754367,30.723582758605566,30.75814017776641,
30.76989427341799,32.22592436849051,31.659771910816055,
27.748759099187865,25.92085927536998,27.4686960384401,
32.011041736555335,29.03583257303588,28.31411031304299,
28.95742784175227,30.27047559931908,26.034154017588495,
29.92435716662582,32.957246436507525,25.59269134984623,
28.886205958952846,34.041223499413746,23.139212855177806,
29.459529926755557,36.18939597357334,22.467987853062972,
27.24773226827232,32.60106332319685,21.872019298734074,
26.89857137673677,30.221983340715354,25.994910077974936,
29.40014717727722,29.06369658274368,24.42466305440126,
31.983268433695,32.99101583993095,25.54244542552471,
29.403299240226872,34.18560234442533,27.16190280357793,
26.761611595958005,30.888021044990378,22.559051858746894,
26.661059276083314,33.134532969891545,20.135912914831863,
26.919560563296592,30.405226290416703,31.681879815161366,
26.01307601390719,31.825147501291415,33.84056483499777,
27.671423952814035,31.305737415954876,36.13237727130725,
28.246082067164977,33.02718519536042,37.58160813227135,
25.74279039542651,29.14301783113196,30.866474389263807,
26.05257695625526,33.83804583025444,33.405289130143494,
23.257995361668378,31.135849159817784,34.37252459062726,
22.16754409752384,32.56571788829859,36.634154323440434,
18.906901271995878,31.8261811814074,37.35558178617904,
18.480374998470666,33.88849773470876,40.050778911552435,
23.05267285351155,29.102409509354338,34.62894339980066,
22.097193378405347,31.430134298358766,32.696602200022724,
22.370606493397823,34.59923879697957,36.37659411976968,
23.32695335273298,32.26940696349741,38.31068520581343,
16.55903011422272,33.694812159192196,40.76750662575639,
18.739989334291465,35.84647873670545,39.46623753027902,
19.772652777692382,33.33072432206953,41.55467009674835,
28.438690498867814,28.92422157667356,36.45506121201092,
30.02375297419524,28.18768008547493,38.56287858767542,
32.50608380526372,29.63840571969633,38.48714970839223,
33.45898002494296,30.440250471896835,40.44748530897122,
27.909720289752922,27.582654074480757,35.20476936481281,
29.046534091311585,28.588933163623256,40.331161335992334,
30.48279765253917,25.33600986876036,38.48651098100795,
32.147255075568445,24.41590795389423,40.661729213522115,
32.60629975391237,21.564237737179656,40.58536160685465,
34.27075717694164,20.644137712039512,42.76057983936881,
34.67163797907014,17.89665564800247,42.626947866056675,
31.280216666473432,24.775246909606533,36.67199986582213,
28.696124090897303,24.311574292315548,38.4531932220959,
31.34983417190819,24.976670913048054,42.47624032870794,
33.93392863721031,25.44034353033904,40.69504697243418,
33.40372065757262,21.003476667751816,38.77085049166883,
30.81962619227051,20.539802160734844,40.55204384794259,
33.42875952692687,21.128601094661477,44.576927768215945,
36.09184642829518,21.60644854659595,42.741283847295136,
35.77452364022229,17.342172248310085,44.082637810451175,
35.52512894255001,17.456812475432773,40.97742305779248,
32.9887406098402,17.001700316871954,42.72573329213166,
33.55018387035006,29.972754938930102,36.21503388606611,
35.91110736552925,31.321595996454565,35.89181704263354,
35.65049336466027,34.01671564306241,36.86173024484372,
37.48574044312339,34.910990672836874,38.049976718645,
32.68489345779488,29.27288924972598,34.66445135113397,
37.402552815361624,30.356766157154787,36.93463587349934,
36.67852508966215,31.28816863343514,33.1064025088193,
39.16235825289188,32.68100062225754,32.617358430399996,
39.92977597702478,32.647573259238115,29.831943896585756,
42.31076836220714,33.995327724318905,29.434670581370263,
43.40072073869059,34.244319910355536,27.16261900972776,
45.585083355265525,35.50532083473331,26.93759232868363,
42.3063086088731,33.232311172014384,25.115596008855377,
35.171797425756736,32.05056157616329,31.92686389237533,
36.776626434929106,29.355275633668416,32.40136086013131,
40.66908591679729,31.91860956925538,33.79689704684396,
39.06425690762492,34.613893622024264,33.322400079087984,
38.45916664594712,33.49013215787162,28.661394706671462,
40.11305294179394,30.712454162256346,29.150395321393013,
43.156116496400415,34.767503778492106,30.961750594177204,
46.40661527281942,35.692993302152885,25.225077284028092,
46.40997709535386,36.268103061015196,28.480498796648448,
43.12783863670102,33.419985529159945,23.40308096419984,
40.659884840956856,32.28185647831293,25.285206475272087,
33.6209049757624,35.132224673821085,36.396761275524284,
};

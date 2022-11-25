#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TXMLEngine.h"
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TH3.h>
#include <TObjArray.h>
#include <TString.h>
#include <TTree.h>
#include <stdio.h>
#include <stdlib.h>

#include <Riostream.h>

void write_tree(TString namefile, int HICnumber, bool verbose);

bool read_PosHIC(int halfcone, int halfdisk, int halflayer, int HICposition,
                 int HICnumber, bool verbose);
void read_allPosHICdiskonly(int halfcone, int disk, int halflayer);

void draw_MFTgeom_misaligned(bool verbose);
void draw_MFTgeom_ideal(bool verbose);
void draw_HICgeom_misaligned(int HICnumber, bool verbose, bool write);
void write_HICgeom_misaligned();
void write_HICgeom_ideal();

void xmlreadfile(
    const char *filename = "geofiles/Geometry.xml",
    const char *nameIdealPosChip = "data_txt_files/idealpositions.txt",
    const char *nameIdealPosDisk = "data_txt_files/idealpositionsdisk.txt");

void ReadMFTSurveyChip(
    const char *outputname = "data_txt_files/survey_true_padpositions.txt");

#endif

// some structures
typedef struct {
  Float_t x, y, z;
} POINT;

typedef struct {
  UInt_t nsensor, nladder;
} CHIPN;

typedef struct {
  UInt_t halfid, diskid, ladderid, chipid, HICid;
} GLOBALID;

POINT point;
CHIPN chipn;
GLOBALID sensor;

POINT pointsurvey;
GLOBALID sensorsurvey;

bool actif = false;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 4) {
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

const int NmaxSides = 2;
const int NmaxLadder00 = 12;
const int NmaxLadder01 = 12;
const int NmaxLadder02 = 13;
const int NmaxLadder03 = 16;
const int NmaxLadder04 = 17;

//______________________________________________________________________
//    BOTTOM MFT
//______________________________________________________________________
const Double_t BottomDISK00[NmaxSides][NmaxLadder00] = {
    //        {2026,3061,3067,3076,3112,3077,3114,3085,3090,3099,2032,2028},
    //        {2030,3101,3102,3105,3118,3127,3121,3106,3108,3111,2040,2036},
    //        {2048,3148,3167,3103,3164,3161,3165,3163,3157,3094,2046,2047},
    //        {2049,3168,3171,3176,3170,3205,3174,3178,3172,3169,2010,2044},
    {2062, 3255, 3248, 3254, 3242, 3243, 3229, 3244, 3250, 3252, 2063,
     2041}, // 00-01_07   Bottom Front
    {2066, 3256, 3257, 3258, 3259, 3253, 3261, 3246, 3251, 3260, 2065,
     2067} // Bottom Back

};

const Double_t BottomDISK01[NmaxSides][NmaxLadder01] = {
    {2071, 3289, 3294, 3279, 3281, 3292, 3283, 3296, 3277, 3291, 2076,
     2073}, // 00-01_04  Bottom Front
    {2074, 3295, 3297, 3284, 3290, 3287, 3293, 3299, 3300, 3285, 2075,
     2077} // Bottom Back
};

const Double_t BottomDISK02[NmaxSides][NmaxLadder02] = {
    {2037, 3116, 3117, 4013, 4012, 3132, 3131, 3130, 4014, 4016, 3150, 3119,
     2042}, // 00-02_04  Bottom Front
    {2043, 3123, 3153, 4020, 4029, 3136, 3129, 3147, 4030, 4031, 3149, 3124,
     2045}, // Bottom Back
    //{2068,3086,3262,4158,4153,3267,3266,3265,4159,4155,3263,3078,2069},
    //{2070,3272,3275,4156,4167,3270,3271,3278,4166,4161,3273,3274,2072},
};
const Double_t BottomDISK03[NmaxSides][NmaxLadder03] = {
    {3216, 3215, 4125, 4126, 4128, 4129, 4130, 4144, 4133, 4132, 4131, 4095,
     4097, 3224, 3226, 3227}, // 00-03_01  Bottom Front
    {3239, 3240, 4141, 4136, 4152, 4138, 4139, 4140, 4142, 4143, 4107, 4151,
     4127, 3245, 3247, 3249} // Bottom Back
    
};
const Double_t BottomDISK04[NmaxSides][NmaxLadder04] = {
    {3207, 3175, 4098, 4096, 5004, 5009, 4089, 4099, 4056, 4045, 4094, 5006,
     5017, 4032, 4010, 3203, 3204}, // 00-04_01  Bottom Front
    {3159, 3162, 4100, 4101, 5018, 5021, 4108, 4112, 4110, 4113, 4106, 5023,
     5022, 4102, 4104, 3209, 3210} // Bottom Back
};


//______________________________________________________________________
//    TOP MFT
//______________________________________________________________________
const Double_t TopDISK00[NmaxSides][NmaxLadder00] = {
    {2029, 3053, 3060, 3050, 3058, 3097, 3095, 3100, 3074, 3075, 2035,
     2038}, // 00-01_01  Top Front
    {2023, 3070, 3007, 3109, 3049, 3083, 3081, 3092, 3073, 3072, 2033,
     2022}, //  Top Back
    //        {2026,3061,3067,3076,3112,3077,3114,3085,3090,3099,2032,2028},
    //        {2030,3101,3102,3105,3118,3127,3121,3106,3108,3111,2040,2036},
    //        {2048,3148,3167,3103,3164,3161,3165,3163,3157,3094,2046,2047},
    //        {2049,3168,3171,3176,3170,3205,3174,3178,3172,3169,2010,2044},

};

const Double_t TopDISK01[NmaxSides][NmaxLadder01] = {
    {2060, 3158, 3196, 3182, 3166, 3191, 3195, 3201, 3194, 3197, 2058,
     2059}, // 00-01_06   Top Front
    {2056, 3200, 3199, 3192, 3186, 3183, 3202, 3190, 3193, 3198, 2054,
     2055}, //  Top Back
};

const Double_t TopDISK02[NmaxSides][NmaxLadder02] = {
    //{2068,3086,3262,4158,4153,3267,3266,3265,4159,4155,3263,3078,2069},
    //{2070,3272,3275,4156,4167,3270,3271,3278,4166,4161,3273,3274,2072},
    {2050, 3238, 3214, 4119, 4116, 3213, 3206, 3208, 4111, 4114, 3211, 3212,
     2061}, // 00-02_02   Top Front
    {2057, 3237, 3235, 4120, 4121, 3236, 3223, 3230, 4122, 4124, 3232, 3234,
     2053}, //  Top Back
};
const Double_t TopDISK03[NmaxSides][NmaxLadder03] = {
    {3139, 3141, 4075, 4023, 4092, 4036, 4067, 4084, 4079, 4080, 4086, 4076,
     4053, 3140, 3142, 3143}, // 00-03_03   Top Front
    {3144, 3145, 4054, 4066, 4059, 4060, 4061, 4062, 4063, 4069, 4083, 4085,
     4068, 3146, 3152, 3110}, //  Top Back
};
const Double_t TopDISK04[NmaxSides][NmaxLadder04] = {
    {3133, 3113, 4007, 4008, 5007, 5013, 4040, 4027, 4026, 4038, 4041, 5014,
     5010, 4019, 4047, 3156, 3134}, // 00-04_02   Top Front
    {3135, 3128, 4022, 4024, 5011, 5016, 4011, 4033, 4048, 4051, 4039, 5020,
     5012, 4034, 4006, 3151, 3137}, //  Top Back
};

double function_nchipPerHalfDisk(int half_MFT = 0, TString namedisk = "00", int half_disk = 0) {
  Int_t n = 0;
  if (half_MFT == 0) {
      if (namedisk == "00") {
        for (int i = 0; i < NmaxLadder00; i++) {
          n = n + Int_t(BottomDISK00[half_disk][i] / 1000.0);
        }
      }
      if (namedisk == "01") {
        for (int i = 0; i < NmaxLadder01; i++) {
          n = n + Int_t(BottomDISK01[half_disk][i] / 1000.0);
        }
      }
      if (namedisk == "02") {
        for (int i = 0; i < NmaxLadder02; i++) {
          n = n + Int_t(BottomDISK02[half_disk][i] / 1000.0);
        }
      }
      if (namedisk == "03") {
        for (int i = 0; i < NmaxLadder03; i++) {
          n = n + Int_t(BottomDISK03[half_disk][i] / 1000.0);
        }
      }
      if (namedisk == "04") {
        for (int i = 0; i < NmaxLadder04; i++) {
          n = n + Int_t(BottomDISK04[half_disk][i] / 1000.0);
        }
      }
  }
    if (half_MFT == 1) {
        if (namedisk == "00") {
          for (int i = 0; i < NmaxLadder00; i++) {
            n = n + Int_t(TopDISK00[half_disk][i] / 1000.0);
          }
        }
        if (namedisk == "01") {
          for (int i = 0; i < NmaxLadder01; i++) {
            n = n + Int_t(TopDISK01[half_disk][i] / 1000.0);
          }
        }
        if (namedisk == "02") {
          for (int i = 0; i < NmaxLadder02; i++) {
            n = n + Int_t(TopDISK02[half_disk][i] / 1000.0);
          }
        }
        if (namedisk == "03") {
          for (int i = 0; i < NmaxLadder03; i++) {
            n = n + Int_t(TopDISK03[half_disk][i] / 1000.0);
          }
        }
        if (namedisk == "04") {
          for (int i = 0; i < NmaxLadder04; i++) {
            n = n + Int_t(TopDISK04[half_disk][i] / 1000.0);
          }
        }
    }
  return n;
};

double function_nchip(int half_MFT = 0, TString namedisk = "00", int half_disk = 0,
                      int ladder = 0) {
  Int_t n = 0;
    if (half_MFT == 0) {
        if (namedisk == "00") {
          n = Int_t(BottomDISK00[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "01") {
          n = Int_t(BottomDISK01[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "02") {
          n = Int_t(BottomDISK02[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "03") {
          n = Int_t(BottomDISK03[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "04") {
          n = Int_t(BottomDISK04[half_disk][ladder] / 1000.0);
        }
    }
    if (half_MFT == 1) {
        if (namedisk == "00") {
          n = Int_t(TopDISK00[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "01") {
          n = Int_t(TopDISK01[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "02") {
          n = Int_t(TopDISK02[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "03") {
          n = Int_t(TopDISK03[half_disk][ladder] / 1000.0);
        }
        if (namedisk == "04") {
          n = Int_t(TopDISK04[half_disk][ladder] / 1000.0);
        }
    }

  return n;
};


const int NmaxHIC = 280;
//______________________________________________________________________
//    ALL HICs in MFT following the geometry ids
//______________________________________________________________________
const Double_t MFT_HICnumberGeo[NmaxHIC] = {
    // BOTTOM
    // layer 0 : DISK00-01/Disk_MD00-01_01/Front_MD00-01_01
    2029, 3053, 3060, 3050, 3058, 3097, 3095, 3100, 3074, 3075, 2035,
     2038,
    // layer 1 : DISK00-01/Disk_MD00-01_01/Back_MD00-01_01
    2023, 3070, 3007, 3109, 3049, 3083, 3081, 3092, 3073, 3072, 2033,
     2022,
    // layer 2 : DISK00-01/Disk_MD00-01_06/Front_MD00-01_06
    2060, 3158, 3196, 3182, 3166, 3191, 3195, 3201, 3194, 3197, 2058,
     2059,
    // layer 3 : DISK00-01/Disk_MD00-01_06/Back_MD00-01_06
    2056, 3200, 3199, 3192, 3186, 3183, 3202, 3190, 3193, 3198, 2054,
     2055,
    // layer 4 : DISK02/Disk_MD02_02/Front_MD02_02
    2037, 3116, 3117, 4013, 4012, 3132, 3131, 3130, 4014, 4016, 3150, 3119,
     2042,
    // layer 5 : DISK02/Disk_MD02_02/Back_MD02_02
    2043, 3123, 3153, 4020, 4029, 3136, 3129, 3147, 4030, 4031, 3149, 3124,
     2045,
    // layer 6 : DISK03/Disk_MD03_03/Front_MD03_03
    3139, 3141, 4075, 4023, 4092, 4036, 4067, 4084, 4079, 4080, 4086, 4076,
     4053, 3140, 3142, 3143,
    // layer 7 : DISK03/Disk_MD03_03/Back_MD03_03
    3110, 3152, 3146, 4068, 4085, 4083, 4069, 4063, 4062, 4061, 4060, 4059, 4066,
    4054, 3145, 3144,
    // layer 8 : DISK04/Disk_MD04_02/Front_MD04_02
    3133, 3113, 4007, 4008, 5007, 5013, 4040, 4027, 4026, 4038, 4041, 5014,
     5010, 4019, 4047, 3156, 3134,
    // layer 9 : DISK04/Disk_MD04_02/Back_MD04_02
    3137, 3151, 4006, 4034, 5012, 5020, 4039, 4051, 4048, 4033, 4011, 5016, 5011,
    4024, 4022, 3128, 3135,
    
    // TOP
    // layer 0 : DISK00-01/Disk_MD00-01_07/Front_MD00-01_07
    2062, 3255, 3248, 3254, 3242, 3243, 3229, 3244, 3250, 3252, 2063,
    2041,
   // layer 1 : DISK00-01/Disk_MD00-01_07/Back_MD00-01_07
   2066, 3256, 3257, 3258, 3259, 3253, 3261, 3246, 3251, 3260, 2065,
    2067,
   // layer 2 : DISK00-01/Disk_MD00-01_04/Front_MD00-01_04
   //2071, 3289, 3294, 3279, 3281, 3292, 3283, 3296, 3277, 3291, 2076, 2073,
   2048, 3148, 3167, 3103, 3164, 3161, 3165, 3163, 3157, 3094, 2046,
    2047,
   // layer 3 : DISK00-01/Disk_MD00-01_04/Back_MD00-01_04
   //2074, 3295, 3297, 3284, 3290, 3287, 3293, 3299, 3300, 3285, 2075, 2077,
   2049, 3168, 3171, 3176, 3170, 3205, 3174, 3178, 3172, 3169, 2010,
    2044,
   // layer 4 : DISK02/Disk_MD02_03/Front_MD02_03
    //2050, 3238, 3214, 4119, 4116, 3213, 3206, 3208, 4111, 4114, 3211, 3212, 2061,
    2068, 3086, 3262, 4158, 4153, 3267, 3266, 3265, 4159, 4155, 3263, 3078,
     2069,
   // layer 5 : DISK02/Disk_MD02_03/Back_MD02_03
    //2053, 3234, 3232, 4124, 4122, 3230, 3223, 3236, 4121, 4120, 3235, 3237, 2057,
    2070, 3272, 3275, 4156, 4167, 3270, 3271, 3278, 4166, 4161, 3273, 3274,
     2072,
   // layer 6 : DISK03/Disk_MD03_01/Front_MD03_01
   3216, 3215, 4125, 4126, 4128, 4129, 4130, 4144, 4133, 4132, 4131, 4095,
    4097, 3224, 3226, 3227,
   // layer 7 : DISK03/Disk_MD03_01/Back_MD03_01
   3239, 3240, 4141, 4136, 4152, 4138, 4139, 4140, 4142, 4143, 4107, 4151,
    4127, 3245, 3247, 3249,
   // layer 8 : DISK04/Disk_MD04_01/Front_MD04_01
   3207, 3175, 4098, 4096, 5004, 5009, 4089, 4099, 4056, 4045, 4094, 5006,
    5017, 4032, 4010, 3203, 3204,
   // layer 9 : DISK04/Disk_MD04_01/Back_MD04_01
   //3210, 3209, 4104, 4102, 5022, 5023, 4106, 4113, 4110, 4112, 4108, 5021, 5018, 4101, 4100, 3162, 3159
    3159, 3162, 4100, 4101, 5018, 5021,  4108,  4112, 4110, 4113, 4106, 5023, 5022, 4102, 4104, 3209, 3210
};

//______________________________________________________________________
// True DISK numbers installed
// OLD one
const Int_t DISK_POS_OLD[10] = {1, 6, 2, 3, 2,
                            7, 4, 4, 1, 1};
// New one GOOD one
const Int_t DISK_POS[10] = {1, 6, 2, 3, 2,  // BOTTOM h0
                            7, 4, 3, 1, 1};   // TOP h1
//______________________________________________________________________

const char *nameDISKintalled[20] = {
    // BOTTOM
    "DISK00-01/Disk_MD00-01_01/Front_MD00-01_01",
    "DISK00-01/Disk_MD00-01_01/Back_MD00-01_01",
    "DISK00-01/Disk_MD00-01_06/Front_MD00-01_06",
    "DISK00-01/Disk_MD00-01_06/Back_MD00-01_06",
    "DISK02/Disk_MD02_02/Front_MD02_02",
    "DISK02/Disk_MD02_02/Back_MD02_02",
    "DISK03/Disk_MD03_03/Front_MD03_03",
    "DISK03/Disk_MD03_03/Back_MD03_03",
    "DISK04/Disk_MD04_02/Front_MD04_02",
    "DISK04/Disk_MD04_02/Back_MD04_02",
    // TOP
    "DISK00-01/Disk_MD00-01_07/Front_MD00-01_07",
    "DISK00-01/Disk_MD00-01_07/Back_MD00-01_07",
    "DISK00-01/Disk_MD00-01_04/Front_MD00-01_04",
    "DISK00-01/Disk_MD00-01_04/Back_MD00-01_04",
    "DISK02/Disk_MD02_03/Front_MD02_03",
    "DISK02/Disk_MD02_03/Back_MD02_03",
    //"DISK02/Disk_MD02_04/Back_MD02_04",
    "DISK03/Disk_MD03_01/Front_MD03_01",
    "DISK03/Disk_MD03_01/Back_MD03_01",
    "DISK04/Disk_MD04_01/Front_MD04_01",
    "DISK04/Disk_MD04_01/Back_MD04_01"
};


Int_t Nbinx = 800;
Int_t Nbiny = 800;
Int_t Nbinz = 300;
Double_t Xmin = -18.0;
Double_t Xmax = 18.0;
Double_t Ymin = -18.0;
Double_t Ymax = 18.0;
Double_t Zmin = -100.0;
Double_t Zmax = -0.0;

Int_t NbinPosChip =
    6400; // it give resolution of 0.001 cm for positionning the chip
Double_t PosChipmin = -3.2;
Double_t PosChipmax = 3.2;
Double_t PosChippas = (PosChipmax - PosChipmin) / NbinPosChip;

Double_t xpas = (Xmax - Xmin) / Nbinx;
Double_t ypas = (Ymax - Ymin) / Nbiny;
Double_t zpas = (Zmax - Zmin) / Nbinz;

TH3F *hgeomHIC_misaligned =
    new TH3F("hgeomHIC_misaligned", "x:y:z distribution", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax, Nbinz, Zmin, Zmax);
TH3F *hgeomHIC_misaligned2 =
    new TH3F("hgeomHIC_misaligned2", "x:y:z distribution", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax, Nbinz, Zmin, Zmax);

TH2D *hprofileHIC_posx = new TH2D("hprofileHIC_posx", "", 5000, 1000, 6000,
                                  NbinPosChip, PosChipmin, PosChipmax);
TH2D *hprofileHIC_posy = new TH2D("hprofileHIC_posy", "", 5000, 1000, 6000,
                                  NbinPosChip, PosChipmin, PosChipmax);
TH2D *hprofileHIC_posz = new TH2D("hprofileHIC_posz", "", 5000, 1000, 6000,
                                  NbinPosChip, PosChipmin, PosChipmax);

TH2D *hprofileHIC_idealposxy_disk00f = new TH2D(
    "hprofileHIC_idealposxy_disk00f", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk00b = new TH2D(
    "hprofileHIC_idealposxy_disk00b", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk01f = new TH2D(
    "hprofileHIC_idealposxy_disk01f", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk01b = new TH2D(
    "hprofileHIC_idealposxy_disk01b", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk02f = new TH2D(
    "hprofileHIC_idealposxy_disk02f", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk02b = new TH2D(
    "hprofileHIC_idealposxy_disk02b", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk03f = new TH2D(
    "hprofileHIC_idealposxy_disk03f", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk03b = new TH2D(
    "hprofileHIC_idealposxy_disk03b", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk04f = new TH2D(
    "hprofileHIC_idealposxy_disk04f", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_idealposxy_disk04b = new TH2D(
    "hprofileHIC_idealposxy_disk04b", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);

TH2D *hprofileHIC_ideal2posxy_disk00f_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk00f_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk00b_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk00b_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk01f_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk01f_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk01b_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk01b_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk02f_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk02f_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk02b_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk02b_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk03f_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk03f_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk03b_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk03b_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk04f_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk04f_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk04b_bottom =
    new TH2D("hprofileHIC_ideal2posxy_disk04b_bottom", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);

TH2D *hprofileHIC_posxy_disk00f_bottom =
    new TH2D("hprofileHIC_posxy_disk00f_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk00b_bottom =
    new TH2D("hprofileHIC_posxy_disk00b_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk01f_bottom =
    new TH2D("hprofileHIC_posxy_disk01f_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk01b_bottom =
    new TH2D("hprofileHIC_posxy_disk01b_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk02f_bottom =
    new TH2D("hprofileHIC_posxy_disk02f_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk02b_bottom =
    new TH2D("hprofileHIC_posxy_disk02b_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk03f_bottom =
    new TH2D("hprofileHIC_posxy_disk03f_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk03b_bottom =
    new TH2D("hprofileHIC_posxy_disk03b_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk04f_bottom =
    new TH2D("hprofileHIC_posxy_disk04f_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk04b_bottom =
    new TH2D("hprofileHIC_posxy_disk04b_bottom", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax);

TH2D *hprofileHIC_ideal2posxy_disk00f_top =
    new TH2D("hprofileHIC_ideal2posxy_disk00f_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk00b_top =
    new TH2D("hprofileHIC_ideal2posxy_disk00b_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk01f_top =
    new TH2D("hprofileHIC_ideal2posxy_disk01f_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk01b_top =
    new TH2D("hprofileHIC_ideal2posxy_disk01b_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk02f_top =
    new TH2D("hprofileHIC_ideal2posxy_disk02f_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk02b_top =
    new TH2D("hprofileHIC_ideal2posxy_disk02b_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk03f_top =
    new TH2D("hprofileHIC_ideal2posxy_disk03f_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk03b_top =
    new TH2D("hprofileHIC_ideal2posxy_disk03b_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk04f_top =
    new TH2D("hprofileHIC_ideal2posxy_disk04f_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_ideal2posxy_disk04b_top =
    new TH2D("hprofileHIC_ideal2posxy_disk04b_top", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax);

TH2D *hprofileHIC_posxy_disk00f_top = new TH2D(
    "hprofileHIC_posxy_disk00f_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk00b_top = new TH2D(
    "hprofileHIC_posxy_disk00b_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk01f_top = new TH2D(
    "hprofileHIC_posxy_disk01f_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk01b_top = new TH2D(
    "hprofileHIC_posxy_disk01b_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk02f_top = new TH2D(
    "hprofileHIC_posxy_disk02f_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk02b_top = new TH2D(
    "hprofileHIC_posxy_disk02b_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk03f_top = new TH2D(
    "hprofileHIC_posxy_disk03f_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk03b_top = new TH2D(
    "hprofileHIC_posxy_disk03b_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk04f_top = new TH2D(
    "hprofileHIC_posxy_disk04f_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);
TH2D *hprofileHIC_posxy_disk04b_top = new TH2D(
    "hprofileHIC_posxy_disk04b_top", "", Nbinx, Xmin, Xmax, Nbiny, Ymin, Ymax);

TH3F *hprofileHIC_posxy_disk04b_ladder0 =
    new TH3F("hprofileHIC_posxy_disk04b_ladder0", "", Nbinx, Xmin, Xmax, Nbiny,
             Ymin, Ymax, Nbinz, Zmin, Zmax);
TH3F *hprofileHIC_ideal2posxy_disk04b_ladder0 =
    new TH3F("hprofileHIC_ideal2posxy_disk04b_ladder0", "", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax, Nbinz, Zmin, Zmax);

// TNtuple *ntupleHICpos = new TNtuple("ntupleHICpos","data from txt
// file","nhic:x:y:z");

TH3D *hgeom_ideal_ladder =
    new TH3D("hgeom_ideal_ladder", "x:y:z distribution", Nbinx, Xmin, Xmax,
             Nbiny, Ymin, Ymax, Nbinz, Zmin, Zmax);

Float_t x, y, z = 0;
Float_t xdisk, ydisk, zdisk = 0;
Float_t xladder, yladder, zladder = 0;
Float_t xchip, ychip, zchip = 0;
Float_t xmft, xhalf, ymft, yhalf, zmft, zhalf = 0;

TTree *treeChipPosition = new TTree("treeChipPos", "treeChipPos");
TTree *treeChipPositionIdeal =
    new TTree("treeChipPosIdeal", "treeChipPosIdeal");

ofstream myfile_idealposdisk;
ofstream myfile_idealpos;
ofstream myfile_surveypos;

Double_t dx_chip = 3.0;
Double_t dy_chip = 1.5;
Double_t dz_chip = 0.0;

Double_t dx_pad_sensor = 0.02;
Double_t dy_pad_sensor = 0.01;
Double_t dz_pad_sensor = 0.0;

Double_t dx_sensor_sensor = 0.01500;
Double_t dy_sensor_sensor = 0.00768 * 2;
Double_t dz_sensor_sensor = 0.0;
Double_t dy_width_sensitive = 0.017532;

Double_t factor_micrometer_centimeter = 0.0001;

Double_t dx_pad_AM_TL = 175.32 * factor_micrometer_centimeter;
Double_t dy_pad_AM_TL = 14852.8 * factor_micrometer_centimeter;

Double_t dx_pad_AM_TR = 29824.68 * factor_micrometer_centimeter;
Double_t dy_pad_AM_TR = 14852.8 * factor_micrometer_centimeter;

Double_t dx_pad_AM_BL = 175.32 * factor_micrometer_centimeter;
Double_t dy_pad_AM_BL = 76.8 * factor_micrometer_centimeter;

Double_t dx_pad_AM_BR = 29824.68 * factor_micrometer_centimeter;
Double_t dy_pad_AM_BR = 76.8 * factor_micrometer_centimeter;

// ideal pad positions in local frame
const Double_t D_AM[4][2] = {{dx_pad_AM_TR, dy_pad_AM_BL},
                             {dx_pad_AM_TR, dy_pad_AM_TR},
                             {dx_pad_AM_BL, dy_pad_AM_TL},
                             {dx_pad_AM_TL, dy_pad_AM_BR}};

// ideal pad positions in local frame
const Double_t Dglo_AM[4][2] = {{dy_pad_AM_BL, dx_pad_AM_TR},
                                {dy_pad_AM_TR, dx_pad_AM_TR},
                                {dy_pad_AM_TR, dx_pad_AM_TL},
                                {dy_pad_AM_BL, dx_pad_AM_TL}};


int nladder = 0;
int k_disk = 0;
int half_cone = 0;

//______________________________________________________________
// Run the workflow (currently really not optimized)
//
void ReadMFTSurvey() {
  xmlreadfile();
  ReadMFTSurveyChip();
}

void ReadMFTSurveyChip(const char *outputname) {
  myfile_surveypos.open(outputname);
  read_allPosHICdiskonly(0, 0, 0);
  read_allPosHICdiskonly(0, 0, 1);

  read_allPosHICdiskonly(0, 1, 2);
  read_allPosHICdiskonly(0, 1, 3);

  read_allPosHICdiskonly(0, 2, 4);
  read_allPosHICdiskonly(0, 2, 5);

  read_allPosHICdiskonly(0, 3, 6);
  read_allPosHICdiskonly(0, 3, 7);

  read_allPosHICdiskonly(0, 4, 8);
  read_allPosHICdiskonly(0, 4, 9);

    
  read_allPosHICdiskonly(1, 0, 10);
  read_allPosHICdiskonly(1, 0, 11);

  read_allPosHICdiskonly(1, 1, 12);
  read_allPosHICdiskonly(1, 1, 13);

  read_allPosHICdiskonly(1, 2, 14);
  read_allPosHICdiskonly(1, 2, 15);

  read_allPosHICdiskonly(1, 3, 16);
  read_allPosHICdiskonly(1, 3, 17);

  read_allPosHICdiskonly(1, 4, 18);
  read_allPosHICdiskonly(1, 4, 19);

  myfile_surveypos.close();
}

void read_allPosHICdiskonly(int halfcone = 0, int halfdisk = 0, int halflayer = 0) {
//  for (int side = 0; side < 2; side++) {
//    for (int pos = 1; pos < 35; pos++) {
//      for (int hic = 0; hic < NmaxHIC; hic++) {
//        if (read_PosHIC(halfcone, halfdisk, side, halflayer, pos, MFT_HICnumberGeo[hic],
//                        true) == false)
//          continue;
//      }
//    }
//  }
    for (int pos = 1; pos < 35; pos++) {
      for (int hic = 0; hic < NmaxHIC; hic++) {
        if (read_PosHIC(halfcone, halfdisk, halflayer, pos, MFT_HICnumberGeo[hic],
                        true) == false)
          continue;
      }
    }
}

bool read_PosHIC(int halfcone = 0, int halfdisk = 0, int halflayer = 0,
                 int HICposition = 0, int HICnumber = 0, bool verbose = true) {
  ifstream in;
  const char *nameSide[2] = {"Front", "Back"};
  const char *nameDisk[5] = {"00-01", "00-01", "02", "03", "04"};
  TString nameHICposition = to_string_with_precision(HICposition);
  TString nameHICnumber = to_string_with_precision(HICnumber);
    
  int halflayer2 = halflayer;
  TString namefile =
      Form("mftalign/positionsHIC/%s/Pos%02d_HIC_%4d.txt",
           nameDISKintalled[halflayer], HICposition, HICnumber);
  // cout<<namefile<<endl;

  in.open(namefile);
    
  if (!in.is_open()) {
    return false;
  }

  printf("\n----> Read the survey measurments for the file:\n----> "
         "%s\n",
         Form(" %s  -- Pos%02d_HIC_%4d.txt", nameDISKintalled[halflayer], HICposition, HICnumber));

  half_cone = halfcone;
  //
  //    if( (halflayer==7) && (halfdisk==0) && (side==0) ) half_cone=0; k_halfdisk= 0;  //
  //    half 0 halfdisk 0 front (bottom) if( (halflayer==7) && (halfdisk==0) && (side==1) )
  //    half_cone=0; k_halfdisk= 0;  // half 0 halfdisk 0 back (bottom) if( (halflayer==4)
  //    && (halfdisk==1) && (side==0) ) half_cone=0; k_halfdisk= 1;  // half 0 halfdisk 1
  //    front (bottom) if( (halflayer==4) && (halfdisk==1) && (side==1) ) half_cone=0;
  //    k_halfdisk= 1;  // half 0 halfdisk 1 back (bottom) if( (halflayer==4) && (halfdisk==2)
  //    && (side==0) ) half_cone=0; k_halfdisk= 2;  // half 0 halfdisk 2 front (bottom)
  //    if( (halflayer==4) && (halfdisk==2) && (side==1) ) half_cone=0; k_halfdisk= 2;  //
  //    half 0 halfdisk 2 back (bottom) if( (halflayer==1) && (halfdisk==3) && (side==0) )
  //    half_cone=0; k_halfdisk= 3;  // half 0 halfdisk 3 front (bottom) if( (halflayer==1)
  //    && (halfdisk==3) && (side==1) ) half_cone=0; k_halfdisk= 3;  // half 0 halfdisk 3
  //    back (bottom) if( (halflayer==1) && (halfdisk==4) && (side==0) ) half_cone=0;
  //    k_halfdisk= 4;  // half 0 halfdisk 4 front (bottom) if( (halflayer==1) && (halfdisk==4)
  //    && (side==1) ) half_cone=0; k_halfdisk= 4;  // half 0 halfdisk 4 back (bottom)
  //
  //    if( (halflayer==1) && (halfdisk==0) && (side==0) ) half_cone=1; k_halfdisk= 0;  //
  //    half 1 halfdisk 0 front (top) if( (halflayer==1) && (halfdisk==0) && (side==1) )
  //    half_cone=1; k_halfdisk= 0;  // half 1 halfdisk 0 back (top) if( (halflayer==6) &&
  //    (halfdisk==1) && (side==0) ) half_cone=1; k_halfdisk= 1;  // half 1 halfdisk 1 front
  //    (top) if( (halflayer==6) && (halfdisk==1) && (side==1) ) half_cone=1; k_halfdisk= 1;
  //    // half 1 halfdisk 1 back (top) if( (halflayer==2) && (halfdisk==2) && (side==0) )
  //    half_cone=1; k_halfdisk= 2;  // half 1 halfdisk 2 front (top) if( (halflayer==2) &&
  //    (halfdisk==2) && (side==1) ) half_cone=1; k_halfdisk= 2;  // half 1 halfdisk 2 back
  //    (top) if( (halflayer==3) && (halfdisk==3) && (side==0) ) half_cone=1; k_halfdisk= 3;
  //    // half 1 halfdisk 3 front (top) if( (halflayer==3) && (halfdisk==3) && (side==1) )
  //    half_cone=1; k_halfdisk= 3;  // half 1 halfdisk 3 back (top) if( (halflayer==2) &&
  //    (halfdisk==4) && (side==0) ) half_cone=1; k_halfdisk= 4;  // half 1 halfdisk 4 front
  //    (top) if( (halflayer==2) && (halfdisk==4) && (side==1) ) half_cone=1; k_halfdisk= 4;
  //    // half 1 halfdisk 4 back (top)

  //    Float_t nPos;
  //    Float_t x,y,z;
  //    Int_t nlines = 0;
  //    //TFile *f = new TFile("MFT_PosHICdata.root","RECREATE");
  //
  //    //TH3F *hgeomHIC_misaligned = new TH3F("hgeomHIC_misaligned","x:y:z
  //    distribution",-30,30,100,-30,30,100,0,-100,100);
  //    //TNtuple *ntuple = new TNtuple("ntuple","data fromascii file","x:y:z");
  //
  //    char line[127];
  //    while (1) {
  //        //in >> x ;
  //
  //        //if (!nPos==nIndice++) break;
  //        //cout<<"Number of line: "<<nlines<<endl;
  //
  //        if (nlines >19 ){
  //            in >> nPos >> x >> y >> z;
  //            //if (in.fail()) break;
  //            //sscanf(&line[0],"%f %f %f %f",&nPos,&x,&y,&z);
  //
  //            //printf("nPos=%d x=%f y=%f z=%f\n",nPos,x,y,z);
  //            //sscanf(&line[0],"%f %f %f %f",&nPos,&x,&y,&z);
  //
  //            printf("nPos=%f x=%f y=%f z=%f\n",nPos,x,y,z);
  //            if (!in.good()) break;
  //
  //            //hgeomHIC_misaligned->Fill(x,y,z);
  //            //ntuple->Fill(x,y,z);
  //        }
  //        nlines++;
  //        if (nlines >40 ) break;
  //    }
  //    printf(" found %d lines\n",nlines);
  //
  //    //ntuple->Write();
  //    in.close();
  //
  //    //f->Write();
  //    //f->Close();

  sensorsurvey.halfid = half_cone;
  sensorsurvey.diskid = halfdisk;
  sensorsurvey.ladderid = nladder;
  sensorsurvey.HICid = HICnumber;

  if ((halfdisk == 0) && (nladder < 24))
    nladder++;
  if ((halfdisk == 1) && (nladder < 24))
    nladder++;
  if ((halfdisk == 2) && (nladder < 26))
    nladder++;
  if ((halfdisk == 3) && (nladder < 32))
    nladder++;
  if ((halfdisk == 4) && (nladder < 34))
    nladder++;

  //    if( (k_halfdisk ==0) &&  (HICposition ==23) ) sensorsurvey.ladderid =
  //    HICposition; if( (k_halfdisk ==1) &&  (HICposition ==23) )
  //    sensorsurvey.ladderid = HICposition; if( (k_halfdisk ==2) &&  (HICposition
  //    ==25) ) sensorsurvey.ladderid = HICposition; if( (k_halfdisk ==3) &&
  //    (HICposition ==31) ) sensorsurvey.ladderid = HICposition; if( (k_halfdisk
  //    ==4) &&  (HICposition ==33) ) sensorsurvey.ladderid = HICposition;

  if ((halfdisk == 0) && (nladder == 24))
    nladder = 0;
  if ((halfdisk == 1) && (nladder == 24))
    nladder = 0;
  if ((halfdisk == 2) && (nladder == 26))
    nladder = 0;
  if ((halfdisk == 3) && (nladder == 32))
    nladder = 0;
  if ((halfdisk == 4) && (nladder == 34))
    nladder = 0;

  sensorsurvey.chipid = 0;

  // cout<<"  halfdisk "<<halfdisk<<"   "<<nladder<<"   pos "<<HICposition<<endl;

  write_tree(namefile, HICnumber, true);
  draw_HICgeom_misaligned(HICnumber, true, false);
  write_HICgeom_misaligned();

  return true;
}

void write_tree(TString namefile, int HICnumber = 0, bool verbose = true) {
  TTree *tree = new TTree("tree", "tree");
  tree->ReadFile(namefile, "bin:x:y:z");

  float bin, x, y, z;
  tree->SetBranchAddress("bin", &bin);
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  // tree->SetBranchAddress("nhic", &HICnumber);
  int nlines = 0;

  TFile *f = new TFile("temp_files/MFT_PosHICtree.root", "RECREATE");

  // Get the bin edges and values out of the tree
  std::vector<float> binEdges, valuesX, valuesY, valuesZ;
  for (int i = 20; i < tree->GetEntries(); i++) {

    tree->GetEntry(i);

    binEdges.push_back(bin);
    valuesX.push_back(x);
    valuesY.push_back(y);
    valuesZ.push_back(z);
    printf("i=%d x=%f y=%f z=%f\n", i, x, y, z);
  }
  tree->ResetBranchAddresses();
  tree->Write();

  f->Write();
  f->Close();
}

void draw_HICgeom_misaligned(int HICnumber = 0, bool verbose = true,
                             bool write = false) {

  TFile *f = new TFile("temp_files/MFT_PosHICtree.root");
  TTree *t1 = (TTree *)f->Get("tree");

  Float_t bin, x, y, z;
  t1->SetBranchAddress("x", &x);
  t1->SetBranchAddress("y", &y);
  t1->SetBranchAddress("z", &z);
  t1->SetBranchAddress("bin", &bin);

  Float_t conversion_cm = 0.1;
  Double_t z0_coneMFT = -46.0;

  Double_t x0_ref = 0.0;
  Double_t y0_ref = 0.0;
  Double_t z0_ref = 0.0;

  Double_t new_x = 0.0;
  Double_t new_y = 0.0;
  Double_t new_z = 0.0;

  Double_t par_limit = 0.5;
  Double_t par_limity = 0.3;

  // read all entries and fill the histograms
  Int_t nentries = (Int_t)t1->GetEntries();
  Int_t n0 = (nentries - 1);

  Float_t new_x_correct = 0.0;
  Float_t new_y_correct = 0.0;
  Float_t new_z_correct = 0.0;

  //    treeChipPosition->Branch("x",&new_x_correct);
  //    treeChipPosition->Branch("y",&new_y_correct);
  //    treeChipPosition->Branch("z",&new_z_correct);
  // treeChipPosition->SetBranchAddress("bin",&bin);
  int compt_pad = 0;

  t1->GetEntry(n0);

  // reference position
  x0_ref = x * conversion_cm;
  y0_ref = y * conversion_cm;
  z0_ref = z * conversion_cm;

  for (Int_t i = 0; i < nentries; i++) {

    t1->GetEntry(i);
    // hgeomHIC_misaligned->Fill(x*conversion_cm,y*conversion_cm,z0_coneMFT+z*conversion_cm);

    // cout<<" x: "<<x*conversion_cm<<"  y: "<<y*conversion_cm<<"  z: "<<
    // z*conversion_cm<<endl;

    for (int ichip = 0; ichip < 6; ichip++) {

      //            if( i==(n0-(ichip*4 -1))  ){
      //                new_x =  x*conversion_cm - dx_pad_sensor;
      //                new_y =  y*conversion_cm + dy_pad_sensor;
      //                new_z =  z*conversion_cm - dz_pad_sensor;
      //            }
      //            if( i==(n0-ichip*4)  ){
      //                new_x =  x*conversion_cm - dx_pad_sensor;
      //                new_y =  y*conversion_cm - dy_pad_sensor;
      //                new_z =  z*conversion_cm - dz_pad_sensor;
      //            }
      //            if( i==(n0-(ichip*4 +1))  ){
      //                new_x =  x*conversion_cm + dx_pad_sensor;
      //                new_y =  y*conversion_cm - dy_pad_sensor;
      //                new_z =  z*conversion_cm - dz_pad_sensor;
      //            }
      //            if( i==(n0-(ichip*4 +2)) ){
      //               new_x =  x*conversion_cm + dx_pad_sensor;
      //               new_y =  y*conversion_cm + dy_pad_sensor;
      //               new_z =  z*conversion_cm - dz_pad_sensor;
      //            }
      new_x = x * conversion_cm;
      new_y = y * conversion_cm;
      new_z = z * conversion_cm;

      new_x_correct = (x0_ref - new_x);
      new_y_correct = (y0_ref - new_y);
      new_z_correct = (z0_ref - new_z);

      // cout<<i<<" x "<<new_x_correct<<"   y "<<new_y_correct<<endl;

      if (i < n0) {

        if ((i == (n0 - 2)) || (i == (n0 - 3))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 1 * abs(D_AM[0][0] - D_AM[2][0]);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 1 * abs(D_AM[0][0] - D_AM[2][0]);
        }
        if ((i == (n0 - 6)) || (i == (n0 - 7))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 1 * abs(D_AM[0][0] - D_AM[2][0]) -
                            dx_chip - dx_sensor_sensor;
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 1 * abs(D_AM[0][0] - D_AM[2][0]) +
                            dx_chip + dx_sensor_sensor;
        }
        if ((i == (n0 - 10)) || (i == (n0 - 11))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 1 * abs(D_AM[0][0] - D_AM[2][0]) -
                            2 * (dx_chip + dx_sensor_sensor);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 1 * abs(D_AM[0][0] - D_AM[2][0]) +
                            2 * (dx_chip + dx_sensor_sensor);
        }
        if ((i == (n0 - 14)) || (i == (n0 - 15))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 1 * abs(D_AM[0][0] - D_AM[2][0]) -
                            3 * (dx_chip + dx_sensor_sensor);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 1 * abs(D_AM[0][0] - D_AM[2][0]) +
                            3 * (dx_chip + dx_sensor_sensor);
        }

        if ((i == (n0 - 4)) || (i == (n0 - 5))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 1 * (dx_chip + dx_sensor_sensor);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 1 * (dx_chip + dx_sensor_sensor);
        }
        if ((i == (n0 - 8)) || (i == (n0 - 9))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 2 * (dx_chip + dx_sensor_sensor);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 2 * (dx_chip + dx_sensor_sensor);
        }
        if ((i == (n0 - 12)) || (i == (n0 - 13))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 3 * (dx_chip + dx_sensor_sensor);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 3 * (dx_chip + dx_sensor_sensor);
        }
        if ((i == (n0 - 16)) || (i == (n0 - 17))) {
          if (new_x_correct > 0.0)
            new_x_correct = new_x_correct - 4 * (dx_chip + dx_sensor_sensor);
          if (new_x_correct < 0.0)
            new_x_correct = new_x_correct + 4 * (dx_chip + dx_sensor_sensor);
        }

        //                if( (i==(n0-18)) || (i==(n0-19))  ){
        //                    if(new_x_correct>0.0) new_x_correct =
        //                    new_x_correct - 5*(dx_chip + dx_sensor_sensor) ;
        //                    if(new_x_correct<0.0) new_x_correct =
        //                    new_x_correct + 5*(dx_chip + dx_sensor_sensor) ;
        //                }

        // if( abs(new_x_correct) < (ichip*(dx_chip + dx_sensor_sensor) +
        // par_limit) && abs(new_x_correct) > (ichip*dx_chip + dx_sensor_sensor
        // - par_limit) ){
        // }

        //                if( abs(new_y_correct) < ( (dy_chip -
        //                dy_width_sensitive) + par_limity) &&
        //                abs(new_y_correct) > ( (dy_chip - dy_width_sensitive)
        //                - par_limity) ){
        //                    if(new_y_correct>0.0) new_y_correct =
        //                    new_y_correct - (dy_chip - dy_width_sensitive);
        //                    if(new_y_correct<0.0) new_y_correct =
        //                    new_y_correct + (dy_chip - dy_width_sensitive);
        //                }
      }

      if (abs(new_y_correct) > (abs(D_AM[2][1] - D_AM[3][1]) - par_limity)) {
        if (new_y_correct > 0.0)
          new_y_correct = new_y_correct - abs(D_AM[2][1] - D_AM[3][1]);
        if (new_y_correct < 0.0)
          new_y_correct = new_y_correct + abs(D_AM[2][1] - D_AM[3][1]);
      }

      //            if( i == (n0+1) ){
      //                if( abs(new_y_correct) < ( (dy_chip -
      //                dy_width_sensitive) + par_limity) && abs(new_y_correct)
      //                > ( (dy_chip - dy_width_sensitive) - par_limity) ){
      //                    if(new_y_correct>0.0) new_y_correct = new_y_correct
      //                    - (dy_chip - dy_width_sensitive);
      //                    if(new_y_correct<0.0) new_y_correct = new_y_correct
      //                    + (dy_chip - dy_width_sensitive);
      //                }
      //            }
    }
    float aa = abs(D_AM[2][1] - D_AM[3][1]) - par_limity;

    // if(new_y_correct>1.0) cout<<i<<"  "<<n0<<"  "<<new_y_correct<<"  "
    // <<aa<<endl;

    if (new_y_correct > 1.0)
      new_y_correct = new_y_correct - abs(D_AM[2][1] - D_AM[3][1]);
    if (new_y_correct < -1.0)
      new_y_correct = new_y_correct + abs(D_AM[2][1] - D_AM[3][1]);

    if (new_y_correct > 1.0)
      cout << i << "  " << n0 << "  " << new_y_correct << "  " << aa << endl;

    //        new_x_correct=(x0_ref - new_x) ;
    //        new_y_correct=(y0_ref - new_y) ;
    //        new_z_correct=(z0_ref - new_z) ;

    //        new_x_correct = new_x ;
    //        new_y_correct = new_y ;
    //        new_z_correct = new_z ;

    hprofileHIC_posx->Fill(HICnumber, new_x_correct);
    hprofileHIC_posy->Fill(HICnumber, new_y_correct);
    hprofileHIC_posz->Fill(HICnumber, new_z_correct);

    // if( (x==0.000) || (y==0.000) ) continue;

    if (compt_pad == 4) {
      compt_pad = 0;
      sensorsurvey.chipid++;
    }

    pointsurvey.x = new_y_correct;
    pointsurvey.y = new_x_correct;
    pointsurvey.z = new_z_correct;

    myfile_surveypos << sensorsurvey.halfid << " " << sensorsurvey.diskid << " "
                     << sensorsurvey.ladderid << " " << sensorsurvey.chipid
                     << " ";
    if (compt_pad == 0)
      myfile_surveypos << compt_pad << " " << pointsurvey.x << " "
                       << pointsurvey.y << " " << pointsurvey.z
                       << " "; // write to file
    if (compt_pad == 1)
      myfile_surveypos << compt_pad << " " << pointsurvey.x << " "
                       << pointsurvey.y << " " << pointsurvey.z
                       << " "; // write to file
    if (compt_pad == 2)
      myfile_surveypos << compt_pad << " " << pointsurvey.x << " "
                       << pointsurvey.y << " " << pointsurvey.z
                       << " "; // write to file
    if (compt_pad == 3)
      myfile_surveypos << compt_pad << " " << pointsurvey.x << " "
                       << pointsurvey.y << " " << pointsurvey.z
                       << " "; // write to file
      
    //myfile_surveypos << HICnumber << " ";  // (!!)  CHECK the HIC number
    myfile_surveypos << "\n";

    compt_pad++;
    //        treeChipPosition->Fill();

    // cout<<" x0_ref - new_x "<<x0_ref - new_x<<" y0_ref - new_y "<<y0_ref -
    // new_y<<endl;
  }

  if (write == true) {
    TFile *fileGeom =
        new TFile("temp_files/MFT_PositionsChip.root", "recreate");
    // hgeomHIC_misaligned->Write();

    hprofileHIC_posx->Write();
    hprofileHIC_posy->Write();
    hprofileHIC_posz->Write();

    // fileGeom->Write();
    fileGeom->Close();
  }

  f->Close();

  // hgeomHIC_misaligned->Draw();
}
void write_HICgeom_misaligned() {
  TFile *fileGeom = new TFile("temp_files/MFT_PositionsChip.root", "recreate");
  // hgeomHIC_misaligned->Write();
  hprofileHIC_posx->Write();
  hprofileHIC_posy->Write();
  hprofileHIC_posz->Write();
  fileGeom->Close();
}

Double_t zPosDisk00f = -45.3;
Double_t zPosDisk00b = -46.7;

Double_t zPosDisk01f = -48.6;
Double_t zPosDisk01b = -50.0;

Double_t zPosDisk02f = -52.4;
Double_t zPosDisk02b = -53.8;

Double_t zPosDisk03f = -67.7;
Double_t zPosDisk03b = -69.1;
Double_t zPosDisk04f = -76.1;
Double_t zPosDisk04b = -77.5;
Int_t compt_chip_disk00 = 0;
Int_t compt_chip_disk01 = 0;
Int_t compt_chip_disk02 = 0;
Int_t compt_chip_disk03 = 0;
Int_t compt_chip_disk04 = 0;

void draw_MFTgeom_misaligned(bool verbose = true) {

  TFile *fideal = new TFile("temp_files/dataMFT_PosHIC_ideal.root");

  TH3D *hgeom_idealgeo = (TH3D *)fideal->Get("hgeom_ideal_ladder");

  TH2D *hgeom_ideal_disk00f =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk00f");
  TH2D *hgeom_ideal_disk00b =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk00b");
  TH2D *hgeom_ideal_disk01f =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk01f");
  TH2D *hgeom_ideal_disk01b =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk01b");
  TH2D *hgeom_ideal_disk02f =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk02f");
  TH2D *hgeom_ideal_disk02b =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk02b");
  TH2D *hgeom_ideal_disk03f =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk03f");
  TH2D *hgeom_ideal_disk03b =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk03b");
  TH2D *hgeom_ideal_disk04f =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk04f");
  TH2D *hgeom_ideal_disk04b =
      (TH2D *)fideal->Get("hprofileHIC_idealposxy_disk04b");

  Float_t conversion_cm = 0.01;

  Int_t nentriestot = (Int_t)hgeom_idealgeo->GetEntries();

  // read all entries and fill the histograms
  Int_t nentries0 = (Int_t)hgeom_ideal_disk00f->GetEntries();
  Int_t nentries1 = (Int_t)hgeom_ideal_disk00b->GetEntries();
  Int_t nentries2 = (Int_t)hgeom_ideal_disk01f->GetEntries();
  Int_t nentries3 = (Int_t)hgeom_ideal_disk01b->GetEntries();
  Int_t nentries4 = (Int_t)hgeom_ideal_disk02f->GetEntries();
  Int_t nentries5 = (Int_t)hgeom_ideal_disk02b->GetEntries();
  Int_t nentries6 = (Int_t)hgeom_ideal_disk03f->GetEntries();
  Int_t nentries7 = (Int_t)hgeom_ideal_disk03b->GetEntries();
  Int_t nentries8 = (Int_t)hgeom_ideal_disk04f->GetEntries();
  Int_t nentries9 = (Int_t)hgeom_ideal_disk04b->GetEntries();

  Int_t nentries = nentries0 + nentries1 + nentries2 + nentries3 + nentries4 +
                   nentries5 + nentries6 + nentries7 + nentries8 + nentries9;

  Float_t bin, x, y, z = 0.0;

  Int_t value = 0.0;
  Int_t value2 = 0.0;

  TFile *fileGeom = new TFile("temp_files/MFT_PositionsChip.root");
  TH3F *hMFTgeomHIC_misaligned = (TH3F *)fileGeom->Get("hgeomHIC_misaligned");
  TH2D *profileHIC_posx = (TH2D *)fileGeom->Get("hprofileHIC_posx");
  TH2D *profileHIC_posy = (TH2D *)fileGeom->Get("hprofileHIC_posy");
  TH2D *profileHIC_posz = (TH2D *)fileGeom->Get("hprofileHIC_posz");

  Float_t xposition = 0.0;
  Float_t yposition = 0.0;
  Float_t zposition = 0.0;
  Float_t z0_coneMFT = -46.0;

  Float_t xPosHICIdeal = 0.0;
  Float_t yPosHICIdeal = 0.0;
  Float_t zPosHICIdeal = 0.0;

  Double_t binProfile = 0;

  int nHIC_00fb, nHIC_01fb, nHIC_02fb, nHIC_03fb, nHIC_04fb = 0;
  int nHIC_00bb, nHIC_01bb, nHIC_02bb, nHIC_03bb, nHIC_04bb = 0;
  int nHIC_00ft, nHIC_01ft, nHIC_02ft, nHIC_03ft, nHIC_04ft = 0;
  int nHIC_00bt, nHIC_01bt, nHIC_02bt, nHIC_03bt, nHIC_04bt = 0;

  int comptnHIC_00fb, comptnHIC_01fb, comptnHIC_02fb, comptnHIC_03fb,
      comptnHIC_04fb = 0;
  int comptnHIC_00bb, comptnHIC_01bb, comptnHIC_02bb, comptnHIC_03bb,
      comptnHIC_04bb = 0;
  int comptnHIC_00ft, comptnHIC_01ft, comptnHIC_02ft, comptnHIC_03ft,
      comptnHIC_04ft = 0;
  int comptnHIC_00bt, comptnHIC_01bt, comptnHIC_02bt, comptnHIC_03bt,
      comptnHIC_04bt = 0;

  Float_t parlimit = 0.6;
  Int_t binLadder = 0;
  Int_t compt_total = 0;
  Int_t mdisk = 0;

  nHIC_03bb = 0;
  nHIC_03bt = 0;
  nHIC_03fb = 0;
  nHIC_03ft = 0;
  nHIC_02bb = 0;
  nHIC_02bt = 0;
  nHIC_02fb = 0;
  nHIC_02ft = 0;
  nHIC_01bb = 0;
  nHIC_01bt = 0;
  nHIC_01fb = 0;
  nHIC_01ft = 0;
  nHIC_00bb = 0;
  nHIC_00bt = 0;
  nHIC_00fb = 0;
  nHIC_00ft = 0;

  for (Int_t k = 0; k < Nbinz; k++) {
    for (Int_t j = 0; j < Nbiny; j++) {
      for (Int_t i = 0; i < Nbinx; i++) {

        value = hgeom_idealgeo->GetBinContent(i, j, k);

        if (value == 0)
          continue;

        xPosHICIdeal = Xmin + i * xpas;
        yPosHICIdeal = Ymin + j * ypas;
        zPosHICIdeal = Zmin + k * zpas;
        // if(  compt_total >900) continue;

        if (((zPosDisk00f + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk00f - parlimit))) { // front  disk00
          if (0.000 < yPosHICIdeal) {
            if (nHIC_00fb < NmaxLadder00) {
              // binProfile = profileHIC_posx->FindBin(DISK00[2][nHIC_00fb]); //
              // bottom MFT
              binLadder = BottomDISK00[0][nHIC_00fb];
              comptnHIC_00fb++;
              if (comptnHIC_00fb == function_nchip(0,"00", 0, nHIC_00fb))
                nHIC_00fb++;
              mdisk = 0;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_00ft < NmaxLadder00) {
              // binProfile = profileHIC_posx->FindBin(DISK00[0][nHIC_00ft]);
              // //top MFT
              binLadder = TopDISK00[0][nHIC_00ft];
              comptnHIC_00ft++;
              if (comptnHIC_00ft == function_nchip(1,"00", 0, nHIC_00ft))
                nHIC_00ft++;
              mdisk = 1;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk00++
                 << " on - DISK 00 FRONT -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }

        if (((zPosDisk00b + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk00b - parlimit))) { // back  disk00
          if (0.000 < yPosHICIdeal) {
            if (nHIC_00bb < NmaxLadder00) {
              // binProfile = profileHIC_posx->FindBin(DISK00[3][nHIC_00bb]);
              // //bottom MFT
              binLadder = BottomDISK00[1][nHIC_00bb];
              comptnHIC_00bb++;
              if (comptnHIC_00bb == function_nchip(0,"00", 1, nHIC_00bb))
                nHIC_00bb++;
              mdisk = 2;
              cout << "----> " << mdisk << endl;
            }

          } else {
            if (nHIC_00bt < NmaxLadder00) {
              // binProfile = profileHIC_posx->FindBin(DISK00[1][nHIC_00bt]);
              // //top MFT
              binLadder = TopDISK00[1][nHIC_00bt];
              comptnHIC_00bt++;
              if (comptnHIC_00bt == function_nchip(1,"00", 1, nHIC_00bt))
                nHIC_00bt++;
              mdisk = 3;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk00++
                 << " on - DISK 00 BACK -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }

        if (((zPosDisk01f + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk01f - parlimit))) { // front disk01
          if (0.000 < yPosHICIdeal) {
            if (nHIC_01fb < NmaxLadder01) {
              //                             //binProfile =
              //                             profileHIC_posx->FindBin(DISK01[0][nHIC_01fb]);
              //                             //bottom MFT
              binLadder = BottomDISK01[0][nHIC_01fb];
              comptnHIC_01fb++;
              if (comptnHIC_01fb == function_nchip(0,"01", 0, nHIC_01fb))
                nHIC_01fb++;
              mdisk = 4;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_01ft < NmaxLadder01) {
              // binProfile = profileHIC_posx->FindBin(DISK01[2][nHIC_01ft]);
              // //top MFT
              binLadder = TopDISK01[0][nHIC_01ft];
              comptnHIC_01ft++;
              if (comptnHIC_01ft == function_nchip(1,"01", 0, nHIC_01ft))
                nHIC_01ft++;
              mdisk = 5;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk01++
                 << " on -- DISK 01 FRONT -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }
        if (((zPosDisk01b + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk01b - parlimit))) { // back disk01
          if (0.000 < yPosHICIdeal) {
            if (nHIC_01bb < NmaxLadder01) {
              // binProfile = profileHIC_posx->FindBin(DISK01[1][nHIC_01bb]);
              // //bottom MFT
              binLadder = BottomDISK01[1][nHIC_01bb];
              comptnHIC_01bb++;
              if (comptnHIC_01bb == function_nchip(0,"01", 1, nHIC_01bb))
                nHIC_01bb++;
              mdisk = 6;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_01bt < NmaxLadder01) {
              //                            //binProfile =
              //                            profileHIC_posx->FindBin(DISK01[3][nHIC_01bt]);
              //                            //top MFT
              binLadder = TopDISK01[1][nHIC_01bt];
              comptnHIC_01bt++;
              if (comptnHIC_01bt == function_nchip(1,"02", 1, nHIC_01bt))
                nHIC_01bt++;
              mdisk = 7;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk01++
                 << " on -- DISK 01 BACK -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }

        if (((zPosDisk02f + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk02f - parlimit))) { // front disk02
          if (0.0000 < yPosHICIdeal) {
            if (nHIC_02fb < NmaxLadder02) {
              // binProfile = profileHIC_posx->FindBin(DISK02[0][nHIC_02fb]);
              // //bottom MFT
              binLadder = BottomDISK02[0][nHIC_02fb];
              comptnHIC_02fb++;
              if (comptnHIC_02fb == function_nchip(0,"02", 0, nHIC_02fb))
                nHIC_02fb++;
              mdisk = 8;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_02ft < NmaxLadder02) {
              // binProfile = profileHIC_posx->FindBin(DISK02[2][nHIC_02ft]);
              // //top MFT
              binLadder = TopDISK02[0][nHIC_02ft];
              comptnHIC_02ft++;
              if (comptnHIC_02ft == function_nchip(1,"02", 0, nHIC_02ft))
                nHIC_02ft++;
              mdisk = 9;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk02++
                 << " on --- DISK 02 FRONT -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }
        if (((zPosDisk02b + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk02b - parlimit))) { // back disk02
          if (0.0000 < yPosHICIdeal) {
            if (nHIC_02bb < NmaxLadder02) {
              // binProfile = profileHIC_posx->FindBin(DISK02[1][nHIC_02bb]);
              // //bottom MFT
              binLadder = BottomDISK02[1][nHIC_02bb];
              comptnHIC_02bb++;
              if (comptnHIC_02bb == function_nchip(0,"02", 1, nHIC_02bb))
                nHIC_02bb++;
              mdisk = 10;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_02bt < NmaxLadder02) {
              // binProfile = profileHIC_posx->FindBin(DISK02[3][nHIC_02bt]);
              // //top MFT
              binLadder = TopDISK02[1][nHIC_02bt];
              comptnHIC_02bt++;
              if (comptnHIC_02bt == function_nchip(1,"02", 1, nHIC_02bt))
                nHIC_02bt++;
              mdisk = 11;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk02++
                 << " on --- DISK 02 BACK -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }

        if (((zPosDisk03f + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk03f - parlimit))) { // front disk03
          if (0.000 < yPosHICIdeal) {
            if (nHIC_03fb < NmaxLadder03) {
              // binProfile = profileHIC_posx->FindBin(DISK03[0][nHIC_03fb]);
              // //bottom MFT
              binLadder = BottomDISK03[0][nHIC_03fb];
              comptnHIC_03fb++;
              if (comptnHIC_03fb == function_nchip(0,"03", 0, nHIC_03fb))
                nHIC_03fb++;
              mdisk = 12;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_03ft < NmaxLadder03) {
              // binProfile = profileHIC_posx->FindBin(DISK03[2][nHIC_03ft]);
              // //top MFT
              binLadder = TopDISK03[0][nHIC_03ft];
              comptnHIC_03ft++;
              if (comptnHIC_03ft == function_nchip(1,"03", 0, nHIC_03ft))
                nHIC_03ft++;
              mdisk = 13;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk03++
                 << " on ---- DISK 03 FRONT -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }
        if (((zPosDisk03b + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk03b - parlimit))) { // back disk03
          if (0.000 < yPosHICIdeal) {
            // cout<<"----> y "<<yPosHICIdeal<<"----> z "<<zPosHICIdeal<<"---->
            // nHIC_03bb "<<nHIC_03bb<<endl;

            if (nHIC_03bb < NmaxLadder03) {
              //                          //binProfile =
              //                          profileHIC_posx->FindBin(DISK03[1][nHIC_03bb]);
              //                          //bottom MFT
              binLadder = BottomDISK03[1][nHIC_03bb];
              comptnHIC_03bb++;
              if (comptnHIC_03bb == function_nchip(0,"03", 1, nHIC_03bb))
                nHIC_03bb++;
              mdisk = 14;
              cout << "----> " << mdisk << endl;
            }
          } else {

            if (nHIC_03bt < NmaxLadder03) {
              // binProfile = profileHIC_posx->FindBin(DISK03[3][nHIC_03bt]);
              // //top MFT
              binLadder = TopDISK03[1][nHIC_03bt];
              comptnHIC_03bt++;
              if (comptnHIC_03bt == function_nchip(1,"03", 1, nHIC_03bt))
                nHIC_03bt++;
              mdisk = 15;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk03++
                 << " on ---- DISK 03 BACK -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }
        if (((zPosDisk04f + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk04f - parlimit))) { // front disk04
          if (0.000 < yPosHICIdeal) {
            if (nHIC_04fb < NmaxLadder04) {
              // binProfile = profileHIC_posx->FindBin(DISK04[0][nHIC_04fb]);
              // //bottom MFT
              binLadder = BottomDISK04[0][nHIC_04fb];
              comptnHIC_04fb++;
              if (comptnHIC_04fb == function_nchip(0,"04", 0, nHIC_04fb))
                nHIC_04fb++;
              mdisk = 16;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_04ft < NmaxLadder04) {
              // binProfile = profileHIC_posx->FindBin(DISK04[2][nHIC_04ft]);
              // //top MFT
              binLadder = TopDISK04[0][nHIC_04ft];
              comptnHIC_04ft++;
              if (comptnHIC_04ft == function_nchip(1,"04", 0, nHIC_04ft))
                nHIC_04ft++;
              mdisk = 17;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk04++
                 << " on ----- DISK 04 FRONT -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }
        if (((zPosDisk04b + parlimit) > zPosHICIdeal) &&
            (zPosHICIdeal > (zPosDisk04b - parlimit))) { // back disk04
          if (0.000 < yPosHICIdeal) {
            if (nHIC_04bb < NmaxLadder04) {
              // binProfile = profileHIC_posx->FindBin(DISK04[1][nHIC_04bb]);
              // //bottom MFT
              binLadder = BottomDISK04[1][nHIC_04bb];
              comptnHIC_04bb++;
              if (comptnHIC_04bb == function_nchip(0,"04", 1, nHIC_04bb))
                nHIC_04bb++;
              mdisk = 18;
              cout << "----> " << mdisk << endl;
            }
          } else {
            if (nHIC_04bt < NmaxLadder04) {
              // binProfile = profileHIC_posx->FindBin(DISK04[3][nHIC_04bt]);
              // //top MFT
              binLadder = TopDISK04[1][nHIC_04bt];
              comptnHIC_04bt++;
              if (comptnHIC_04bt == function_nchip(1,"04", 1, nHIC_04bt))
                nHIC_04bt++;
              mdisk = 19;
              cout << "----> " << mdisk << endl;
            }
          }
          if (verbose)
            cout << compt_total << " chip: " << compt_chip_disk04++
                 << " on ----- DISK 04 BACK -----> x " << zPosHICIdeal << endl;
          compt_total++;
        }

        Int_t compt_chip = 0;
        xposition = 0;
        yposition = 0;
        zposition = 0;
        Float_t length_chip = 3.0;
        Float_t width_chip = 1.5;

        for (Int_t h = 0; h < NbinPosChip; h++) {
          if (compt_chip > 3)
            continue;

          if (profileHIC_posx->FindBin(binLadder, h) != 0.000) {
            // xposition = (PosChipmin + h*PosChippas) + xPosHICIdeal;
            if ((mdisk == 18) || (mdisk == 16) || (mdisk == 14) ||
                (mdisk == 12) || (mdisk == 10) || (mdisk == 8) ||
                (mdisk == 6) || (mdisk == 4) || (mdisk == 2) || (mdisk == 0)) {
              length_chip = -3.0;
              width_chip = 1.5;
            }
            if ((mdisk == 19) || (mdisk == 17) || (mdisk == 15) ||
                (mdisk == 13) || (mdisk == 11) || (mdisk == 9) ||
                (mdisk == 7) || (mdisk == 5) || (mdisk == 3) || (mdisk == 1)) {
              length_chip = 3.0;
              width_chip = 1.5;
            }
            if (compt_chip == 0)
              xposition = (0.0 + h * PosChippas) + xPosHICIdeal + width_chip;
            if (compt_chip == 1)
              xposition = (0.0 + h * PosChippas) + xPosHICIdeal;
            if (compt_chip == 2)
              xposition = (0.0 + h * PosChippas) + xPosHICIdeal;
            if (compt_chip == 3)
              xposition = (0.0 + h * PosChippas) + xPosHICIdeal + width_chip;
          }
          if (profileHIC_posy->FindBin(binLadder, h) != 0.000) {
            // yposition = (PosChipmin + h*PosChippas) + yPosHICIdeal;
            if ((mdisk == 18) || (mdisk == 16) || (mdisk == 14) ||
                (mdisk == 12) || (mdisk == 10) || (mdisk == 8) ||
                (mdisk == 6) || (mdisk == 4) || (mdisk == 2) || (mdisk == 0)) {
              length_chip = -3.0;
              width_chip = 1.5;
            }
            if ((mdisk == 19) || (mdisk == 17) || (mdisk == 15) ||
                (mdisk == 13) || (mdisk == 11) || (mdisk == 9) ||
                (mdisk == 7) || (mdisk == 5) || (mdisk == 3) || (mdisk == 1)) {
              length_chip = 3.0;
              width_chip = 1.5;
            }
            if (compt_chip == 0)
              yposition = (0.0 + h * PosChippas) + yPosHICIdeal + length_chip;
            if (compt_chip == 1)
              yposition = (0.0 + h * PosChippas) + yPosHICIdeal + length_chip;
            if (compt_chip == 2)
              yposition = (0.0 + h * PosChippas) + yPosHICIdeal;
            if (compt_chip == 3)
              yposition = (0.0 + h * PosChippas) + yPosHICIdeal;
          }
          if (profileHIC_posz->FindBin(binLadder, h) != 0.000) {
            // zposition = (PosChipmin + h*PosChippas) + zPosHICIdeal;

            if (compt_chip == 0)
              zposition = (0.0 + h * PosChippas) + zPosHICIdeal;
            if (compt_chip == 1)
              zposition = (0.0 + h * PosChippas) + zPosHICIdeal;
            if (compt_chip == 2)
              zposition = (0.0 + h * PosChippas) + zPosHICIdeal;
            if (compt_chip == 3)
              zposition = (0.0 + h * PosChippas) + zPosHICIdeal;
          }

          if (xposition * yposition * zposition == 0.00000)
            continue;

          // cout<<binLadder<<"   h "<<h<<"   dhx
          // "<<(PosChipmin+h*PosChippas)<<" --- xPosHICIdeal:
          // "<<xPosHICIdeal<<endl;
          if (verbose)
            cout << binLadder << "   x " << xposition << "   y " << yposition
                 << "  z: " << zposition << endl;

          // hgeomHIC_misaligned2->Fill(xposition,yposition,zposition);
          // cout<<"   -- profileHIC_posx->GetBinContent(h)
          // "<<profileHIC_posx->GetBinContent(h)<<endl;

          if (((zPosDisk00f + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk00f - parlimit))) { // front  disk00
            if (mdisk == 0)
              hprofileHIC_posxy_disk00f_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 1)
              hprofileHIC_posxy_disk00f_bottom->Fill(xposition - width_chip,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 0))
              hprofileHIC_ideal2posxy_disk00f_top->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 1))
              hprofileHIC_ideal2posxy_disk00f_bottom->Fill(xPosHICIdeal,
                                                           yPosHICIdeal);
          }
          if (((zPosDisk00b + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk00b - parlimit))) { // back  disk00
            if (mdisk == 2)
              hprofileHIC_posxy_disk00b_top->Fill(xposition - width_chip,
                                                  yposition - length_chip);
            if (mdisk == 3)
              hprofileHIC_posxy_disk00b_bottom->Fill(xposition,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 2))
              hprofileHIC_ideal2posxy_disk00b_top->Fill(xPosHICIdeal,
                                                        yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 3))
              hprofileHIC_ideal2posxy_disk00b_bottom->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
          }
          if (((zPosDisk01f + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk01f - parlimit))) { // front disk01
            if (mdisk == 4)
              hprofileHIC_posxy_disk01f_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 5)
              hprofileHIC_posxy_disk01f_bottom->Fill(xposition - width_chip,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 4))
              hprofileHIC_ideal2posxy_disk01f_top->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 5))
              hprofileHIC_ideal2posxy_disk01f_bottom->Fill(xPosHICIdeal,
                                                           yPosHICIdeal);
          }
          if (((zPosDisk01b + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk01b - parlimit))) { // back disk01
            if (mdisk == 6)
              hprofileHIC_posxy_disk01b_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 8)
              hprofileHIC_posxy_disk01b_bottom->Fill(xposition,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 6))
              hprofileHIC_ideal2posxy_disk01b_top->Fill(xPosHICIdeal,
                                                        yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 7))
              hprofileHIC_ideal2posxy_disk01b_bottom->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
          }
          if (((zPosDisk02f + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk02f - parlimit))) { // front disk02
            if (mdisk == 8)
              hprofileHIC_posxy_disk02f_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 9)
              hprofileHIC_posxy_disk02f_bottom->Fill(xposition - width_chip,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 8))
              hprofileHIC_ideal2posxy_disk02f_top->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 9))
              hprofileHIC_ideal2posxy_disk02f_bottom->Fill(xPosHICIdeal,
                                                           yPosHICIdeal);
          }
          if (((zPosDisk02b + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk02b - parlimit))) { // back disk02
            if (mdisk == 10)
              hprofileHIC_posxy_disk02b_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 11)
              hprofileHIC_posxy_disk02b_bottom->Fill(xposition,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 10))
              hprofileHIC_ideal2posxy_disk02b_top->Fill(xPosHICIdeal,
                                                        yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 11))
              hprofileHIC_ideal2posxy_disk02b_bottom->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
          }
          if (((zPosDisk03f + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk03f - parlimit))) { // front disk03
            if (mdisk == 12)
              hprofileHIC_posxy_disk03f_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 13)
              hprofileHIC_posxy_disk03f_bottom->Fill(xposition - width_chip,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 12))
              hprofileHIC_ideal2posxy_disk03f_top->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 13))
              hprofileHIC_ideal2posxy_disk03f_bottom->Fill(xPosHICIdeal,
                                                           yPosHICIdeal);
          }
          if (((zPosDisk03b + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk03b - parlimit))) { // back disk03
            if (mdisk == 14)
              hprofileHIC_posxy_disk03b_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 15)
              hprofileHIC_posxy_disk03b_bottom->Fill(xposition,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 14))
              hprofileHIC_ideal2posxy_disk03b_top->Fill(xPosHICIdeal,
                                                        yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 15))
              hprofileHIC_ideal2posxy_disk03b_bottom->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);
          }
          if (((zPosDisk04f + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk04f - parlimit))) { // front disk04
            if (mdisk == 16)
              hprofileHIC_posxy_disk04f_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 17)
              hprofileHIC_posxy_disk04f_bottom->Fill(xposition - width_chip,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 16))
              hprofileHIC_ideal2posxy_disk04f_top->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal + length_chip);
            if ((compt_chip == 0) && (mdisk == 17))
              hprofileHIC_ideal2posxy_disk04f_bottom->Fill(xPosHICIdeal,
                                                           yPosHICIdeal);

            // if( compt_chip == 1 )
            // hprofileHIC_ideal2posxy_disk04f->Fill(xPosHICIdeal,yPosHICIdeal+
            // length_chip); if( compt_chip == 2 )
            // hprofileHIC_ideal2posxy_disk04f->Fill(xPosHICIdeal,yPosHICIdeal);
            // if( compt_chip == 3 )
            // hprofileHIC_ideal2posxy_disk04f->Fill(xPosHICIdeal+
            // width_chip,yPosHICIdeal);
          }
          if (((zPosDisk04b + parlimit) > zPosHICIdeal) &&
              (zPosHICIdeal > (zPosDisk04b - parlimit))) { // back disk04
            if (mdisk == 18)
              hprofileHIC_posxy_disk04b_top->Fill(xposition,
                                                  yposition - length_chip);
            if (mdisk == 19)
              hprofileHIC_posxy_disk04b_bottom->Fill(xposition,
                                                     yposition - length_chip);
            if ((compt_chip == 0) && (mdisk == 18))
              hprofileHIC_ideal2posxy_disk04b_top->Fill(xPosHICIdeal,
                                                        yPosHICIdeal);
            if ((compt_chip == 0) && (mdisk == 19))
              hprofileHIC_ideal2posxy_disk04b_bottom->Fill(
                  xPosHICIdeal + width_chip, yPosHICIdeal);

            // if( (compt_chip == 0) && (mdisk == 18) )
            // hprofileHIC_ideal2posxy_disk04b->Fill(xPosHICIdeal+ width_chip,
            // yPosHICIdeal+ length_chip); if( compt_chip == 1 )
            // hprofileHIC_ideal2posxy_disk04b->Fill(xPosHICIdeal,yPosHICIdeal+
            // length_chip); if( compt_chip == 2 )
            // hprofileHIC_ideal2posxy_disk04b->Fill(xPosHICIdeal,yPosHICIdeal);
            // if( compt_chip == 3 )
            // hprofileHIC_ideal2posxy_disk04b->Fill(xPosHICIdeal+
            // width_chip,yPosHICIdeal);

            //                          if( binLadder == 3210 ){
            //                              hprofileHIC_posxy_disk04b_ladder0->Fill(xposition,yposition,zposition);
            //                              //if( (function_nchip("04",1,0)*4  )
            //                              ==
            //                              hprofileHIC_posxy_disk04b->GetEntries()
            //                              ) break;
            //                          }

            if ((binLadder == 5021) && (mdisk == 18))
              hprofileHIC_posxy_disk04b_ladder0->Fill(
                  xposition, yposition - length_chip, zposition);
            if ((compt_chip == 0) && (binLadder == 5021) && (mdisk == 18))
              hprofileHIC_ideal2posxy_disk04b_ladder0->Fill(
                  xPosHICIdeal, yPosHICIdeal, zPosHICIdeal);
          }

          compt_chip++;

        } // loop for

      } // loop z axis ideal
    }   // loop y axis ideal
  }     // loop x axis ideal

  // fideal->Close();

  // fileGeom->Write();
  fileGeom->Close();

  hgeomHIC_misaligned2->SetLineColor(2);

  TFile *fileGeomMisaligned =
      new TFile("temp_files/MFT_drawGeomMisaligned.root", "recreate");
  hprofileHIC_posxy_disk00f_bottom->Write();
  hprofileHIC_ideal2posxy_disk00f_bottom->Write();
  hprofileHIC_posxy_disk00f_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk00f_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk00b_bottom->Write();
  hprofileHIC_ideal2posxy_disk00b_bottom->Write();
  hprofileHIC_posxy_disk00b_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk00b_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk01f_bottom->Write();
  hprofileHIC_ideal2posxy_disk01f_bottom->Write();
  hprofileHIC_posxy_disk01f_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk01f_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk01b_bottom->Write();
  hprofileHIC_ideal2posxy_disk01b_bottom->Write();
  hprofileHIC_posxy_disk01b_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk01b_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk02f_bottom->Write();
  hprofileHIC_ideal2posxy_disk02f_bottom->Write();
  hprofileHIC_posxy_disk02f_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk02f_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk02b_bottom->Write();
  hprofileHIC_ideal2posxy_disk02b_bottom->Write();
  hprofileHIC_posxy_disk02b_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk02b_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk03f_bottom->Write();
  hprofileHIC_ideal2posxy_disk03f_bottom->Write();
  hprofileHIC_posxy_disk03f_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk03f_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk03b_bottom->Write();
  hprofileHIC_ideal2posxy_disk03b_bottom->Write();
  hprofileHIC_posxy_disk03b_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk03b_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk04f_bottom->Write();
  hprofileHIC_ideal2posxy_disk04f_bottom->Write();
  hprofileHIC_posxy_disk04f_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk04f_bottom->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk04b_bottom->Write();
  hprofileHIC_ideal2posxy_disk04b_bottom->Write();
  hprofileHIC_posxy_disk04b_bottom->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk04b_bottom->GetYaxis()->SetTitle("Y (cm)");

  hprofileHIC_posxy_disk00f_top->Write();
  hprofileHIC_ideal2posxy_disk00f_top->Write();
  hprofileHIC_posxy_disk00f_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk00f_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk00b_top->Write();
  hprofileHIC_ideal2posxy_disk00b_top->Write();
  hprofileHIC_posxy_disk00b_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk00b_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk01f_top->Write();
  hprofileHIC_ideal2posxy_disk01f_top->Write();
  hprofileHIC_posxy_disk01f_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk01f_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk01b_top->Write();
  hprofileHIC_ideal2posxy_disk01b_top->Write();
  hprofileHIC_posxy_disk01b_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk01b_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk02f_top->Write();
  hprofileHIC_ideal2posxy_disk02f_top->Write();
  hprofileHIC_posxy_disk02f_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk02f_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk02b_top->Write();
  hprofileHIC_ideal2posxy_disk02b_top->Write();
  hprofileHIC_posxy_disk02b_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk02b_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk03f_top->Write();
  hprofileHIC_ideal2posxy_disk03f_top->Write();
  hprofileHIC_posxy_disk03f_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk03f_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk03b_top->Write();
  hprofileHIC_ideal2posxy_disk03b_top->Write();
  hprofileHIC_posxy_disk03b_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk03b_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk04f_top->Write();
  hprofileHIC_ideal2posxy_disk04f_top->Write();
  hprofileHIC_posxy_disk04f_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk04f_top->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk04b_top->Write();
  hprofileHIC_ideal2posxy_disk04b_top->Write();
  hprofileHIC_posxy_disk04b_top->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk04b_top->GetYaxis()->SetTitle("Y (cm)");

  hgeom_idealgeo->Write();
  hgeomHIC_misaligned2->Write();
  hprofileHIC_posxy_disk04b_ladder0->Write();
  hprofileHIC_posxy_disk04b_ladder0->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_posxy_disk04b_ladder0->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_posxy_disk04b_ladder0->GetZaxis()->SetTitle("Z (cm)");
  hprofileHIC_ideal2posxy_disk04b_ladder0->Write();
  hprofileHIC_ideal2posxy_disk04b_ladder0->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_ideal2posxy_disk04b_ladder0->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_ideal2posxy_disk04b_ladder0->GetZaxis()->SetTitle("Z (cm)");

  hprofileHIC_posxy_disk00f_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk00b_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk01f_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk01b_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk02f_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk02b_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk03f_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk03b_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk04f_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk04b_top->SetMarkerStyle(2);
  hprofileHIC_posxy_disk00f_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk00b_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk01f_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk01b_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk02f_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk02b_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk03f_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk03b_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk04f_bottom->SetMarkerStyle(2);
  hprofileHIC_posxy_disk04b_bottom->SetMarkerStyle(2);

  //    fileGeomMisaligned->Write();
  fileGeomMisaligned->Close();
}

void draw_MFTgeom_ideal(bool verbose = true) {

  TFile *f = new TFile("temp_files/dataMFT_PosHIC_ideal.root");

  // TFile *fileGeom = new TFile("MFT_PosHIC.root","RECREATE");

  TH3D *hgeom_idealgeo = (TH3D *)f->Get("hgeom_ideal_ladder");

  Float_t conversion_cm = 0.01;

  // read all entries and fill the histograms
  Int_t nentries = (Int_t)hgeom_idealgeo->GetEntries();

  Int_t value = 0.0;

  for (Int_t i = 0; i < Nbinx; i++) {
    for (Int_t j = 0; j < Nbiny; j++) {
      for (Int_t k = 0; k < Nbinz; k++) {
        value = hgeom_idealgeo->GetBinContent(i, j, k);
        // if( value !=0) cout<<value<<endl;
        for (Int_t m = 0; m < value; m++) {
          hgeomHIC_misaligned->Fill(Xmin + i * xpas, Ymin + j * ypas,
                                    Zmin + k * zpas);
        }
      }
    }
  }

  f->Close();

  hgeomHIC_misaligned->Draw();
}

bool isCorrectValue(Float_t value) {
  if ((value < -100.0) && (value > 100.0)) {
    return false;
  }
  return true;
}

void DisplayNode(TXMLEngine *xml, XMLNodePointer_t node, Int_t level);

void xmlreadfile(const char *filename, const char *nameIdealPosChip,
                 const char *nameIdealPosDisk) {
  // First create engine
  TXMLEngine *xml = new TXMLEngine;
  // Now try to parse xml file
  // Only file with restricted xml syntax are supported
  XMLDocPointer_t xmldoc = xml->ParseFile(filename);
  if (xmldoc == 0) {
    delete xml;
    return;
  }

  myfile_idealpos.open(nameIdealPosChip);
  myfile_idealposdisk.open(nameIdealPosDisk);

  // take access to main node
  XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
  // display recursively all nodes and subnodes
  DisplayNode(xml, mainnode, 1);

  write_HICgeom_ideal();

  //    for(int i=0; i<valuesX.size(); i++){
  //        point.x = valuesX.at(i);
  //         point.y = valuesY.at(i);
  //         point.z = valuesZ.at(i);
  //        cout<<i<<endl;
  //         chipn.nsensor = i;
  //         chipn.nladder = binChip;
  //        tree.Fill();
  //    }
  //    tree.Write();
  //    f->Write();
  //    f->Close();

  myfile_idealpos.close();
  myfile_idealposdisk.close();

  // Release memory before exit
  xml->FreeDoc(xmldoc);

  delete xml;
}
int kpad = 0;

void DisplayNode(TXMLEngine *xml, XMLNodePointer_t node, Int_t level) {

  // this function display all accessible information about xml node and its
  // children
  // printf("%*c node: %s\n",level,' ', xml->GetNodeName(node));

  // TH3F *hgeom_ideal = new TH3F("hgeom_ideal_ladder","x:y:z
  // distribution",100,-30,30,100,-30,30,100,-100,0);

  // display namespace
  XMLNsPointer_t ns = xml->GetNS(node);
  if (ns != 0) {
    printf("%*c namespace: %s refer: %s\n", level + 2, ' ', xml->GetNSName(ns),
           xml->GetNSReference(ns));
  }

  Int_t binChip = 0;
  Int_t binladder = 0;

  const TString namemft = "MFT";
  const TString namehalf = "half";
  const TString namedisk = "disk";
  const TString nameladder = "ladder";
  const TString namechip = "chip";

  const TString nameposx = "xpos";
  const TString nameposy = "ypos";
  const TString nameposz = "zpos";

  Double_t xPosHICIdeal = 0.0;
  Double_t yPosHICIdeal = 0.0;
  Double_t zPosHICIdeal = 0.0;

  const TString nameihalf = "top";
  const TString nameidisk = "idisk";
  const TString nameiladder = "iladder";
  const TString nameichip = "ichip";

  int halfid = 0;
  int diskid = 0;
  int ladderid = 0;
  int chipid = 0;

  // display attributes
  XMLAttrPointer_t attr = xml->GetFirstAttr(node);
  while (attr != 0) {
    // printf("%*c attr: %s value: %s\n",level+2,' ', xml->GetAttrName(attr),
    // xml->GetAttrValue(attr));

    if (xml->GetNodeName(node) == namemft) {
      // printf("%*c attr: %s value: %s\n",level+2,' ', xml->GetAttrName(attr),
      // xml->GetAttrValue(attr));

      if (xml->GetAttrName(attr) == nameposx) {
        // xmft= atof( xml->GetAttrValue(attr) );
        xmft = 0.0;
      }

      if (xml->GetAttrName(attr) == nameposy) {
        // ymft= atof( xml->GetAttrValue(attr) );
        ymft = 0.0;
      }
      if (xml->GetAttrName(attr) == nameposz) {
        // zmft= atof( xml->GetAttrValue(attr) );
        zmft = 0.0;
      }
    }

    if (xml->GetNodeName(node) == namehalf) {
      if (xml->GetAttrName(attr) == nameihalf) {
        sensor.halfid = atof(xml->GetAttrValue(attr));
      }
      if (xml->GetAttrName(attr) == nameposx) {
        // xhalf= atof( xml->GetAttrValue(attr) );
        xhalf = 0.0;
      }

      if (xml->GetAttrName(attr) == nameposy) {
        // yhalf= atof( xml->GetAttrValue(attr) );
        yhalf = 0.0;
      }
      if (xml->GetAttrName(attr) == nameposz) {
        // zhalf= atof( xml->GetAttrValue(attr) );
        zhalf = -46.0;
      }
    }

    if (xml->GetNodeName(node) == namedisk) {
      if (xml->GetAttrName(attr) == nameidisk) {
        sensor.diskid = atof(xml->GetAttrValue(attr));
      }

      if (xml->GetAttrName(attr) == nameposx) {
        xdisk = atof(xml->GetAttrValue(attr));
      }

      if (xml->GetAttrName(attr) == nameposy) {
        ydisk = atof(xml->GetAttrValue(attr));
      }
      if (xml->GetAttrName(attr) == nameposz) {
        zdisk = atof(xml->GetAttrValue(attr));
      }

      Double_t xd = xmft + xhalf + xdisk;
      Double_t yd = ymft + yhalf + ydisk;
      Double_t zd = zmft + zhalf + zdisk;

      myfile_idealposdisk << sensor.halfid << " " << sensor.diskid << " " << xd
                          << " " << yd << " " << zd << " \n";
    }

    if (xml->GetNodeName(node) == nameladder) {
      if (xml->GetAttrName(attr) == nameiladder) {
        sensor.ladderid = atof(xml->GetAttrValue(attr));
      }

      if (xml->GetAttrName(attr) == nameposx) {
        xladder = atof(xml->GetAttrValue(attr));
      }

      if (xml->GetAttrName(attr) == nameposy) {
        yladder = atof(xml->GetAttrValue(attr));
      }
      if (xml->GetAttrName(attr) == nameposz) {
        zladder = atof(xml->GetAttrValue(attr));
        binladder++;
      }
    }

    if (xml->GetNodeName(node) == namechip) {

      if (xml->GetAttrName(attr) == nameichip) {
        sensor.chipid = atof(xml->GetAttrValue(attr));
      }

      if (xml->GetAttrName(attr) == nameposx) {
        xchip = atof(xml->GetAttrValue(attr));
      }

      if (xml->GetAttrName(attr) == nameposy) {
        ychip = atof(xml->GetAttrValue(attr));
      }
      if (xml->GetAttrName(attr) == nameposz) {
        zchip = atof(xml->GetAttrValue(attr));

        //               if(isCorrectValue(xmft)==false) xmft=0;
        //               if(isCorrectValue(xhalf)==false) xhalf=0;
        //               if(isCorrectValue(xdisk)==false) xdisk=0;
        //               if(isCorrectValue(xladder)==false) xladder=0;
        //               if(isCorrectValue(xchip)==false) xchip=0;

        if (actif == false) {
          point.x = xmft + xhalf + xdisk + xladder + ychip - dy_chip / 2.0;
          point.y = ymft + yhalf + ydisk + yladder - xchip - dx_chip;
          point.z = zmft + zhalf + zdisk + zladder + zchip;
        }
        if (actif == true) {
          point.x = xmft + xhalf + xdisk + xladder + ychip - dy_chip / 2.0;
          point.y = ymft + yhalf + ydisk - yladder - (-xchip);
          point.z = zmft + zhalf + zdisk - zladder - zchip;
        }

        //               if(actif==true){
        //                   point.x = xmft + xhalf + xdisk ;
        //                   point.y = ymft + yhalf + ydisk ;
        //                   point.z = zmft + zhalf + zdisk ;
        //
        //               }
        //               if(actif==false){
        //                   point.x = xmft + xhalf + xdisk ;
        //                   point.y = ymft + yhalf + ydisk ;
        //                   point.z = zmft + zhalf + zdisk ;
        //
        //               }

        chipn.nsensor = binChip++;
        chipn.nladder = binladder;

        //
        //                 binEdges.push_back(binChip++);
        //                 valuesX.push_back(x);
        //                 valuesY.push_back(y);
        //                 valuesZ.push_back(z);

        // cout<<" bin : "<<bin<<" binladder : "<<binladder<<endl;
        // cout<<"            point.x : "<<point.x<<" point.y : "<<point.y<<"
        // point.z "<<point.z<<" bin: "<<binChip<<"  binladder:
        // "<<binladder<<endl;

        // hgeom_ideal->Fill(point.x,point.y,point.z);
        // tree.Fill();
        for (int ipad = 0; ipad < 4; ipad++) {
          myfile_idealpos << sensor.halfid << " " << sensor.diskid << " "
                          << sensor.ladderid << " " << sensor.chipid << " "
                          << ipad << " ";
          //                   if(ipad==0) myfile_idealpos<< point.x +
          //                   D_AM[ipad][0] <<" "<< point.y + D_AM[ipad][1] <<"
          //                   "<< point.z << " " ; //write to file if(ipad==1)
          //                   myfile_idealpos<< point.x + dx_chip +
          //                   D_AM[ipad][0]<<" "<< point.y + D_AM[ipad][1]<<"
          //                   "<< point.z << " " ; //write to file if(ipad==2)
          //                   myfile_idealpos<< point.x + dx_chip +
          //                   D_AM[ipad][0]<<" "<< point.y + dy_chip +
          //                   D_AM[ipad][1] <<" "<< point.z << " " ; //write to
          //                   file if(ipad==3) myfile_idealpos<< point.x +
          //                   D_AM[ipad][0]<<" "<< point.y + dy_chip +
          //                   D_AM[ipad][1] <<" "<< point.z << " " ; //write to
          //                   file
          if (ipad == 0)
            myfile_idealpos << point.x + Dglo_AM[ipad][0] << " "
                            << point.y + Dglo_AM[ipad][1] << " " << point.z
                            << " "; // write to file
          if (ipad == 1)
            myfile_idealpos << point.x + Dglo_AM[ipad][0] << " "
                            << point.y + Dglo_AM[ipad][1] << " " << point.z
                            << " "; // write to file
          if (ipad == 2)
            myfile_idealpos << point.x + Dglo_AM[ipad][0] << " "
                            << point.y + Dglo_AM[ipad][1] << " " << point.z
                            << " "; // write to file
          if (ipad == 3)
            myfile_idealpos << point.x + Dglo_AM[ipad][0] << " "
                            << point.y + Dglo_AM[ipad][1] << " " << point.z
                            << " "; // write to file
          myfile_idealpos << "\n";
        }
      }
    }

    // cout<<" ---  "<<halfid<<" "<<diskid<<" "<<ladderid<<"  "<<chipid<<endl;
    attr = xml->GetNextAttr(attr);
  }

  // display content (if exists)
  const char *content = xml->GetNodeContent(node);
  if (content != 0)
    printf("%*c cont: %s\n", level + 2, ' ', content);

  // display all child nodes
  XMLNodePointer_t child = xml->GetChild(node);
  while (child != 0) {
    DisplayNode(xml, child, level + 2);
    child = xml->GetNext(child);
    // cout<<" -------level "<<level<<endl;
    // cout<<" -- level "<<level<<"         point.x : "<<point.x<<" point.y :
    // "<<point.y<<" point.z "<<point.z<<" chipn.nsensor: "<<chipn.nsensor<<"
    // chipn.nladder: "<<chipn.nladder<<endl;

    if (level == 1)
      actif = true;
    Double_t pointX = point.x;
    Double_t pointY = point.y;
    Double_t pointZ = point.z;

    if (level == 7)
      hgeom_ideal_ladder->Fill(pointX, pointY, pointZ);

    if (level == 7) {

      Float_t parlimit = 0.25;
      //           Double_t zPosDisk00f = -45.3;
      //           Double_t zPosDisk00b = -46.7;
      //
      //           Double_t zPosDisk01f = -48.6;
      //           Double_t zPosDisk01b = -50.0;
      //
      //           Double_t zPosDisk02f = -52.4;
      //           Double_t zPosDisk02b = -53.8;
      //
      //           Double_t zPosDisk03f = -67.7;
      //           Double_t zPosDisk03b = -69.1;
      //           Double_t zPosDisk04f = -76.1;
      //           Double_t zPosDisk04b = -77.5;

      if (((zPosDisk00f + parlimit) > point.z) &&
          (point.z > (zPosDisk00f - parlimit))) { // front  disk00
        hprofileHIC_idealposxy_disk00f->Fill(point.x, point.y);
      }

      if (((zPosDisk00b + parlimit) > point.z) &&
          (point.z > (zPosDisk00b - parlimit))) { // back  disk00
        hprofileHIC_idealposxy_disk00b->Fill(point.x, point.y);
      }

      if (((zPosDisk01f + parlimit) > point.z) &&
          (point.z > (zPosDisk01f - parlimit))) { // front disk01
        hprofileHIC_idealposxy_disk01f->Fill(point.x, point.y);
      }
      if (((zPosDisk01b + parlimit) > point.z) &&
          (point.z > (zPosDisk01b - parlimit))) { // back disk01
        hprofileHIC_idealposxy_disk01b->Fill(point.x, point.y);
      }

      if (((zPosDisk02f + parlimit) > point.z) &&
          (point.z > (zPosDisk02f - parlimit))) { // front disk02
        hprofileHIC_idealposxy_disk02f->Fill(point.x, point.y);
      }
      if (((zPosDisk02b + parlimit) > point.z) &&
          (point.z > (zPosDisk02b - parlimit))) { // back disk02
        hprofileHIC_idealposxy_disk02b->Fill(point.x, point.y);
      }

      if (((zPosDisk03f + parlimit) > point.z) &&
          (point.z > (zPosDisk03f - parlimit))) { // front disk03
        hprofileHIC_idealposxy_disk03f->Fill(point.x, point.y);
      }
      if (((zPosDisk03b + parlimit) > point.z) &&
          (point.z > (zPosDisk03b - parlimit))) { // back disk03
        hprofileHIC_idealposxy_disk03b->Fill(point.x, point.y);
      }
      if (((zPosDisk04f + parlimit) > point.z) &&
          (point.z > (zPosDisk04f - parlimit))) { // front disk04
        hprofileHIC_idealposxy_disk04f->Fill(point.x, point.y);
      }
      if (((zPosDisk04b + parlimit) > point.z) &&
          (zPosHICIdeal > (zPosDisk04b - parlimit))) { // back disk04
        hprofileHIC_idealposxy_disk04b->Fill(point.x, point.y);
      }
    }

    if (level == 3)
      cout << " --- level " << level << "         point.x : " << point.x
           << " point.y : " << point.y << " point.z " << point.z << endl;
    if (level == 1)
      cout << "  level " << level << "         point.x : " << point.x
           << " point.y : " << point.y << " point.z " << point.z << endl;
  }
}

void write_HICgeom_ideal() {

  TFile *f = new TFile("temp_files/dataMFT_PosHIC_ideal.root", "RECREATE");
  hgeom_ideal_ladder->Write();
  hgeom_ideal_ladder->GetXaxis()->SetTitle("X (cm)");
  hgeom_ideal_ladder->GetYaxis()->SetTitle("Y (cm)");
  hgeom_ideal_ladder->GetZaxis()->SetTitle("Z (cm)");

  hprofileHIC_idealposxy_disk00f->Write();
  hprofileHIC_idealposxy_disk00f->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk00f->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk00b->Write();
  hprofileHIC_idealposxy_disk00b->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk00b->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk01f->Write();
  hprofileHIC_idealposxy_disk01f->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk01f->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk01b->Write();
  hprofileHIC_idealposxy_disk01b->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk01b->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk02f->Write();
  hprofileHIC_idealposxy_disk02f->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk02f->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk02b->Write();
  hprofileHIC_idealposxy_disk02b->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk02b->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk03f->Write();
  hprofileHIC_idealposxy_disk03f->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk03f->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk03b->Write();
  hprofileHIC_idealposxy_disk03b->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk03b->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk04f->Write();
  hprofileHIC_idealposxy_disk04f->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk04f->GetYaxis()->SetTitle("Y (cm)");
  hprofileHIC_idealposxy_disk04b->Write();
  hprofileHIC_idealposxy_disk04b->GetXaxis()->SetTitle("X (cm)");
  hprofileHIC_idealposxy_disk04b->GetYaxis()->SetTitle("Y (cm)");

  f->Close();
}

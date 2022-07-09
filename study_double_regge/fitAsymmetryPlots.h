#ifndef FITASYMMETRYPLOTS
#define FITASYMMETRYPLOTS

#include <limits>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <chrono>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TTree.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TMath.h>
#include <TSystem.h>

using namespace std;


float fltmin=-1*std::numeric_limits<float>::max();
float fltmax=std::numeric_limits<float>::max();

struct criteria{
    // Can specify up to 2 regions that can be unattached.
    //   regions will be defined by [min,max] and should be merged with an or
    //   If you only care about region1 then we can choose a region2 which is always false, i.e. [fltmin,fltmin] .
    //   [fltmin,fltmin] should always be false since we do not use greater/less than or equal to sign. We use {<, >}  not {<=, >=}
    string variable;
    float minval1;
    float maxval1;
    float minval2;
    float maxval2;
};

struct defaultSelections{
    vector<criteria> selections={
        {"unusedEnergy",fltmin,0.01,fltmin,fltmin},
        {"chiSq",fltmin,16,fltmin,fltmin},
        {"photonTheta1",2.5,10.3,11.9,fltmax},
        {"photonTheta2",2.5,10.3,11.9,fltmax},
        {"photonTheta3",2.5,10.3,11.9,fltmax},
        {"photonTheta4",2.5,10.3,11.9,fltmax},
        {"photonE1",0.1,fltmax,fltmin,fltmin},
        {"photonE2",0.1,fltmax,fltmin,fltmin},
        {"photonE3",0.1,fltmax,fltmin,fltmin},
        {"photonE4",0.1,fltmax,fltmin,fltmin},
        {"proton_momentum",0.3,fltmax,fltmin,fltmin},
        {"proton_z",52,78,fltmin,fltmin},
        {"mmsq",-0.05,0.05,fltmin,fltmin},
    };
};

struct constructFitArgs{
    string weightVar;
    bool freePhase;
    vector<double> fluxRatios_90_0;
    vector<double> fluxRatios_45_135;
    defaultSelections defSels;
};

void constructAndFit(
    unordered_map<string, float*> &variable_map, 
    unordered_map<string, vector<float>> &array_variable_map,  
    unordered_map<string, vector<int>> &array_BeamAngles_map, 
    vector<criteria> selections,
    constructFitArgs args,
    string fitSaveLoc,
    ofstream* saveCsv
    );

map<string,float> extractAsymmetries(
        string sbTag, 
        vector<float> sbRegion, 
        constructFitArgs args,
        map<string,vector<pair<int,criteria>>> eventSelects
        );



#endif

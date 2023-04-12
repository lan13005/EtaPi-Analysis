#include "fitAsymmetryPlots.h"
#include <ctime>
#include <chrono>

// This one is for Asymmetry vs t_recoil which is equal to u3 in vincent/colin language.
//#include "/d/grid13/gluex/gluex_top/gluex_style.C"

double degToRad=TMath::Pi()/180;
// par[3] is used to shift phase by the para or perp orientation, either 0 for para or 90 for perp. 0/-45 is para and 45/90 is perp. 
//   THIS EQUATION IS ALSO USED FOR INSTRUMENTAL ASYMMETRY CALCULATION - MAYBE NEEDS A FLIP OF SIGN IN FRONT OF PAR[1]
int numDOFsig_sc = 3;
Double_t shiftedCos(Double_t *x, Double_t *par){
	return par[0]*(1.0 - par[1]*TMath::Cos(2*degToRad*(x[0]-par[2])));
}
int numDOFsig_flat = 1;
Double_t flat(Double_t *x, Double_t *par){
	return par[0];
}

// par0 and par1 = polarization fraction of perp and para components respectively 
int numDOFsig_asym=4;
Double_t asymmetry(Double_t *x, Double_t *par){
	return ((par[0]+par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3]))/(2+(par[0]-par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3]))));
}

// Number of t1 bins to consider, have to manually enter the names
int nt1bins=5; 
string mint1s[5] = {"01", "03", "05", "07", "09"};
//int nt1bins=1; 
//string mint1s[1] = {"05"};

int minEntriesInPhiHists=50; // set a minimum number entries in each phi histogram
int minEntriesForAsymPhiHists=30; // set a minimum number entries in each phi histogram
// WL for weighted log-likelihood but probably doenst work for asymmetry histogram since the entries are not counts
// M for improved MINUIT results
// E for error estimation - Minos
string fitOption = "E S Q"; 
bool verbose=true;
bool fitPhi=false; // (True) exit program if fitting phi (shifted cosine fits) do not converge after maxiter times 
int maxiter=200; // maximum reinitializations for asymmetry fits before we give up and exit, manually need to address this
string fileType="png"; // linux program 'montage' has trouble with pdfs...
string folder="/d/grid17/ln16/dselector_v3/study_double_regge/rootFiles/";
string outputFolder="fitAsymmetryPlots_results/";


static const int nDataSets = 3;
string dataSetTag[nDataSets] = { "2017_1" , "2018_1", "2018_8" };
string dataSetTag2[nDataSets+1] = { "2017_1" , "2018_1", "2018_8", "total" };
string polnames[7] = {"phi000","phi045","phi090","phi135","instAsym_000","instAsym_045","phiAMO"};
// We dont actually have a pol offset for the total dataset so we will just choose the largest dataset, 2018_1 
//    These numbers are from Alex's rho SDME analysis
float Phi0_0_90[nDataSets+1] = {3.1,4.2,3.1,4.2}; 
float Phi0_45_135[nDataSets+1] = {-45+3.2,-45+2.8,-45+3.1,-45+2.8};
// Last 2 elements is or {0,90} and {45/-45} summed which should be ~flat. Tegan/Will uses the asymmetry offsets. These -999 will be updated in the dataset loop
float poloffset[6] = {0, 45, 90, -45, -999, -999}; 

// From mike dugger's script
vector<vector<float>>  polfractions = { //2017, 2018_1, 2018_8, total
    {0.3537,0.3484,0.3472,0.3512}, // 0, 45, 90, -45 
    {0.3458,0.3416,0.3448,0.3597},
    {0.3563,0.3403,0.3430,0.3523},
    {0.3458,0.3416,0.3448,0.3597} // copy the 2018_1 fractions for the totals
};
vector<vector<float>>  polfractions_err = {
    {0.0100,0.0102,0.0098,0.0102},  
    {0.0055,0.0056,0.0055,0.0057},
    {0.0070,0.0075,0.0074,0.0077},
    {0.0055,0.0056,0.0055,0.0057} // copy the 2018_1 fractions for the totals
};

// ########################
// Current sideband parameters
// ########################
float pi0_peak=0.136;
float pi0_std=0.006425;
float eta_peak=0.548;
float eta_std=0.0123;
// ########################
// Sidebands that was used for older PWA
// ########################
//float pi0_peak=0.135881;
//float pi0_std=0.0076;
//float eta_peak=0.548625;
//float eta_std=0.0191;

//void fitAsymmetryPlots(int single){
int main(int argc, char* argv[]){
    if (argc!=2){
        cout << "ONLY TAKES IN 1 ARGUMENT!" << endl;
        return 1;
    }
    int single=stoi(argv[1]);

    std::chrono::steady_clock::time_point begin_time = std::chrono::steady_clock::now();
    cout << "*************** Running over systematic " << single << " ***************" << endl;

    map<string, float> percentages;
    vector< vector<float> > sbRegions = {
        //{3.0, 1.0, 2.0, 3.0, 1.0, 2.0}, // OLD PWA WIDTHS
        {3.0, 1.0, 3.0, 3.0, 2.0, 6.0}, // Nominal
    }; 
    bool freePhase=false;
    float fluxNormScaleFactor=1; // Scale flux normalization by +/- X% for systematics
    map<string,vector<pair<int,criteria>>> eventSelects={
        {"default", {make_pair(0,criteria{"default",0,0,0,0})} }, // Including a default where we do not make any selection modifications to the default
    };
    bool showSummary=false; // This is only used for event selection systematics section which is ran in a loop so we could do comparisons

    /////////////////////
    // SELECTION SYSTEMATICS + STANDARD RESULTS
    //   To check if  the selections is working in general we can look at outputFolder/diagnostics/ histograms
    //   where the "defaultSelections" are applied to the event selection variables. This is available for all systematic studies
    //   but is particually useful for this one where we are varying the selections. It does appear to work
    /////////////////////
    vector<string> sbTags={"_ASBS"}; 
    //Will overwrite the branch defined by weightVar with AccWeight*weightBS if weightVar!='AccWeight'
    //   and where weightBS is manullay defined around line 1030. Currently used to overwrite weightASBS branch to w.e. we want. i.e. if we do not want 
    //   to use pi0 sidebands
    string weightVar="weightASBS"; 
    showSummary=true;
    // L=Looser, T=tighter bounds
    eventSelects={
//        {"default", {make_pair(0,criteria{"default",0,0,0,0})} }, // Including a default where we do not make any selection modifications to the default
//        {"_evtSel_ueL1", {make_pair(0,criteria{"unusedEnergy",fltmin,0.13,fltmin,fltmin})} },
//        {"_evtSel_ueL2", {make_pair(0,criteria{"unusedEnergy",fltmin,0.17,fltmin,fltmin})} },
//        {"_evtSel_chiT", {make_pair(1,criteria{"chiSq",fltmin,11.5,fltmin,fltmin})} },
        {"_evtSel_chiL", {make_pair(1,criteria{"chiSq",fltmin,16,fltmin,fltmin})} },
//        {"_evtSel_gThetaBeamL", {   
//                                make_pair(2,criteria{"photonTheta1",2.1,10.3,11.9,fltmax}),
//                                make_pair(3,criteria{"photonTheta2",2.1,10.3,11.9,fltmax}),
//                                make_pair(4,criteria{"photonTheta3",2.1,10.3,11.9,fltmax}),
//                                make_pair(5,criteria{"photonTheta4",2.1,10.3,11.9,fltmax}),
//                            } },
//        {"_evtSel_gThetaBeamT", {
//                                make_pair(2,criteria{"photonTheta1",2.8,10.3,11.9,fltmax}),
//                                make_pair(3,criteria{"photonTheta2",2.8,10.3,11.9,fltmax}),
//                                make_pair(4,criteria{"photonTheta3",2.8,10.3,11.9,fltmax}),
//                                make_pair(5,criteria{"photonTheta4",2.8,10.3,11.9,fltmax}),
//                            } },
//        {"_evtSel_gThetaTransL", {  
//                                make_pair(2,criteria{"photonTheta1",2.5,11.4,11.4,fltmax}),
//                                make_pair(3,criteria{"photonTheta2",2.5,11.4,11.4,fltmax}),
//                                make_pair(4,criteria{"photonTheta3",2.5,11.4,11.4,fltmax}),
//                                make_pair(5,criteria{"photonTheta4",2.5,11.4,11.4,fltmax}),
//                            } },
//        {"_evtSel_gThetaTransT", {
//                                make_pair(2,criteria{"photonTheta1",2.5,10.0,12.5,fltmax}),
//                                make_pair(3,criteria{"photonTheta2",2.5,10.0,12.5,fltmax}),
//                                make_pair(4,criteria{"photonTheta3",2.5,10.0,12.5,fltmax}),
//                                make_pair(5,criteria{"photonTheta4",2.5,10.0,12.5,fltmax}),
//                            } },
//        {"_evtSel_gET1", {
//                                make_pair(6,criteria{"photonE1",0.11,fltmax,fltmin,fltmin}),
//                                make_pair(7,criteria{"photonE2",0.11,fltmax,fltmin,fltmin}),
//                                make_pair(8,criteria{"photonE3",0.11,fltmax,fltmin,fltmin}),
//                                make_pair(9,criteria{"photonE4",0.11,fltmax,fltmin,fltmin}),
//                            } },
//        {"_evtSel_gET2", {
//                                make_pair(6,criteria{"photonE1",0.12,fltmax,fltmin,fltmin}),
//                                make_pair(7,criteria{"photonE2",0.12,fltmax,fltmin,fltmin}),
//                                make_pair(8,criteria{"photonE3",0.12,fltmax,fltmin,fltmin}),
//                                make_pair(9,criteria{"photonE4",0.12,fltmax,fltmin,fltmin}),
//                            } },
//        {"_evtSel_pMomT1", {make_pair(10,criteria{"proton_momentum",0.35,fltmax,fltmin,fltmin})} },
//        {"_evtSel_pMomT2", {make_pair(10,criteria{"proton_momentum",0.40,fltmax,fltmin,fltmin})} },
//        {"_evtSel_pZL", {make_pair(11,criteria{"proton_z",51,79,fltmin,fltmin})} },
//        {"_evtSel_pZT", {make_pair(11,criteria{"proton_z",53,77,fltmin,fltmin})} },
//        {"_evtSel_mmsqT1", {make_pair(12,criteria{"mmsq",-0.024,0.024,fltmin,fltmin})} },
//        {"_evtSel_mmsqT2", {make_pair(12,criteria{"mmsq",-0.030,0.030,fltmin,fltmin})} },
    };

    /////////////////////
    // EVENT SELECTION SCANS - MAINLY FOR CHISQ, UE, AND PROTON MOMENTUM WHICH WAS SHOWING LARGE BARLOWS
    //    Typically we would do these scans for a particular asymmetry binning (i.e. integrated)
    /////////////////////
//    vector<string> sbTags={"_ASBS"}; 
//    string weightVar="weightASBS"; // Will not overwrite AccWeight if "AccWeight" else overwrite weightVar with AccWeight*weightBS
//    showSummary=true;
//    outputFolder="fitAsymmetryPlots_results_scan/";
//    eventSelects={
//        {"_ueLT015", {make_pair(0,criteria{"unusedEnergy",fltmin,0.15,fltmin,fltmin})} },
//        {"_ueLT025", {make_pair(0,criteria{"unusedEnergy",fltmin,0.25,fltmin,fltmin})} },
//        {"_ueLT040", {make_pair(0,criteria{"unusedEnergy",fltmin,0.40,fltmin,fltmin})} },
//        {"_ueLT060", {make_pair(0,criteria{"unusedEnergy",fltmin,0.60,fltmin,fltmin})} },
//        {"_ueLT080", {make_pair(0,criteria{"unusedEnergy",fltmin,0.80,fltmin,fltmin})} },
//        {"_ueLT100", {make_pair(0,criteria{"unusedEnergy",fltmin,1.00,fltmin,fltmin})} },
//        {"_chiLT11", {make_pair(1,criteria{"chiSq",fltmin,11,fltmin,fltmin})} }, 
//        {"_chiLT14", {make_pair(1,criteria{"chiSq",fltmin,14,fltmin,fltmin})} },
//        {"_chiLT17", {make_pair(1,criteria{"chiSq",fltmin,17,fltmin,fltmin})} },
//        {"_chiLT20", {make_pair(1,criteria{"chiSq",fltmin,20,fltmin,fltmin})} },
//        {"_chiLT25", {make_pair(1,criteria{"chiSq",fltmin,25,fltmin,fltmin})} },
//        {"_chiLT30", {make_pair(1,criteria{"chiSq",fltmin,30,fltmin,fltmin})} },
//        {"_chiLT40", {make_pair(1,criteria{"chiSq",fltmin,40,fltmin,fltmin})} },
//        {"_pMomGT035", {make_pair(10,criteria{"proton_momentum",0.35,fltmax,fltmin,fltmin})} },
//        {"_pMomGT040", {make_pair(10,criteria{"proton_momentum",0.40,fltmax,fltmin,fltmin})} },
//        {"_pMomGT045", {make_pair(10,criteria{"proton_momentum",0.45,fltmax,fltmin,fltmin})} },
//        {"_pMomGT050", {make_pair(10,criteria{"proton_momentum",0.50,fltmax,fltmin,fltmin})} },
//        {"_pMomGT055", {make_pair(10,criteria{"proton_momentum",0.55,fltmax,fltmin,fltmin})} },
//        {"_pMomGT060", {make_pair(10,criteria{"proton_momentum",0.60,fltmax,fltmin,fltmin})} },
//        {"_pMomGT065", {make_pair(10,criteria{"proton_momentum",0.65,fltmax,fltmin,fltmin})} },
//        {"_pMomGT070", {make_pair(10,criteria{"proton_momentum",0.70,fltmax,fltmin,fltmin})} },
//    };

    /////////////////////
    // SPLIT ONE SIDE OF THE ETA MASS DISTRIBUTION 
    //   NO NEED TO MODIFY WEIGHTS. SPLITTING SIGNAL AND SIDEBAND REGIONS IN HALF, CANCELS
    //   need to make the selection on Meta manually do not select on Mpi0
    //   Can check mass_plots_etaLeft_ASBS.pdf and mass_plots_etaRight_ASBS.pdf to see if the selection is working. It appears so
    /////////////////////
//    vector<string> sbTags={"_etaRight_ASBS"};  
//    string weightVar="weightASBS";

    /////////////////////
    // PHI OFFSET SYSTEMATIC 
    /////////////////////
//    vector<string> sbTags={"_freePhi_ASBS"};  
//    string weightVar="weightASBS";
//    freePhase=true; // free all phase offsets for systematics

    /////////////////////
    // FLUX NORMALIZATION SYSTEMATIC 
    //  ** CANNOT DO BOTH REGIONS WITH DIFFERENT PROCESSES IN runFitAsymPlots.py SINCE IT SELECTS ON eventSelectSingles ** 
    /////////////////////
//    vector<string> sbTags={"_fluxNormPlus5Perc_ASBS"};  
//    string weightVar="weightASBS";
//    fluxNormScaleFactor=1.05; // Scale flux normalization by +/- X% for systematics


    /////////////////////
    // SIDEBAND SYSTEMATICS 
    /////////////////////
//    sbRegions = {
//        //{2.5, 1.5, 2.0, 2.5, 2.5, 4.0}, //Narrow
//        //{3.5, 0.5, 4.0, 3.5, 1.5, 8.0} //Wider
//        {0, 0, 0, 2.5, 2.5, 4.0}, //Narrow, no longer subtracting pi0
//        {0, 0, 0, 3.5, 1.5, 8.0} //Wider, no longer subtracting pi0
//    }; 
//    vector<string> sbTags={"_sb_252540_ASBS", "_sb_351580_ASBS"}; 
//    string weightVar="weightASBS"; // Will not overwrite AccWeight if "AccWeight" else overwrite weightVar with AccWeight*weightBS

    /////////////////////
    // CHECK TO SEE IF WE SHOULD SIDEBAND SUBTRACT OR NOT
    //   Fast eta has less data so we construct only 1 bin for the test. Dont do the additional binning in u3,  s12, s23. Integrated case only
    //   Only select on eta sidebands, where basically all the background exists
    //   Need to change 'nt1bins' to reflect the binning
    /////////////////////
//    vector<string> sbTags={"_sbL_1bin"}; // {"_sbR"}, {"_sbL"}, {"_sbR_1bin"}, {"_sbL_1bin"} 
//    string weightVar="AccWeight"; // only subtract accidentals since we are testing sideband subtraction


    /////////////////////
    // END 
    /////////////////////


    map<string,vector<pair<int,criteria>>> eventSelectSingles;
    int ikey=0;
    for (auto const& it : eventSelects){
        if (ikey==single)
            eventSelectSingles[it.first] = it.second;
        ++ikey;
    }

    if (eventSelectSingles.size()==0){
        cout << "Process " << single << " did not have any systematic to run over!" << endl;
        return 1;
    }

    /////////////////////
    // EXTRACT ASYMS IN A LOOP
    /////////////////////
    for (int i=0; i<(int)sbTags.size(); ++i){ // loop mainly for sideband subtraction systematic scan but used by all scans
        vector<double> fluxRatios_90_0 = {  4.346818e+12/4.188001e+12*fluxNormScaleFactor, 0.965429*fluxNormScaleFactor, 0.918503*fluxNormScaleFactor }; 
        vector<double> fluxRatios_45_135 = {  4.076065e+12/4.095013e+12*fluxNormScaleFactor , 1.02261*fluxNormScaleFactor, 1.03254*fluxNormScaleFactor };
        constructFitArgs args = {weightVar, freePhase, fluxRatios_90_0, fluxRatios_45_135};
        percentages = extractAsymmetries(sbTags[i],sbRegions[i],args,eventSelectSingles);
    }
    if (showSummary){
        cout << endl;
        cout << "---------------------------- SUMMARY ---------------------------------" << endl;
        for (auto const& mapElement: percentages){
            //cout << mapElement.first << " has " << std::setprecision(3) << mapElement.second/percentages["_ASBS"]*100 << "%" 
            //    << " events compared to standard default selections" << endl;
            cout << mapElement.first << " has " << std::setprecision(3) << mapElement.second << " percentage of events post-to-pre selection" << endl;
        } 
        cout << "----------------------------------------------------------------------" << endl;
        cout << endl;
    }

    std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    std::cout << "Total run time = " << std::chrono::duration_cast<std::chrono::seconds>(end_time - begin_time).count() << "s" << std::endl;

    return 0;
}


void constructAndFit(
    unordered_map<string, float*> &variable_map, // holds all the variables of interest in the tree as values in a map. Used to fill array_variable_map 
    unordered_map<string, vector<float>> &array_variable_map, // contains all the data that passed the defaultSelections 
    unordered_map<string, vector<int>> &array_BeamAngles_map, 
    vector<criteria> selections,
    constructFitArgs args,
    string fitSaveLoc,
    ofstream* saveCsv
    ){
    
    TLatex *text = new TLatex();
    text->SetTextSize(0.08);
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << selections[0].minval1;
    std::string sminval1 = ss.str();
    ss.str(""); ss.clear();
    ss << std::fixed << std::setprecision(2) << selections[0].maxval1;
    std::string smaxval1 = ss.str();
    bool overlayBinInfoIntoAsymPlots=true;

    cout << "Constructing phi histograms and fitting" << endl;

    string weightVar=args.weightVar;
    bool freePhase=args.freePhase;
    vector<double> fluxRatios_90_0=args.fluxRatios_90_0;
    vector<double> fluxRatios_45_135=args.fluxRatios_45_135;

    string phis_orientation[5] =  {"phi000", "phi045", "phi090", "phi135", "phiAMO"};
    float t1HalfWidth=1.0/nt1bins/2;

    gStyle->SetOptFit(111);
    gStyle->SetStatY(1);
    gStyle->SetStatX(1);
    gStyle->SetStatW(0.18);
    gStyle->SetStatH(0.18);
    gStyle->SetFitFormat("3.2g");

    // *****************************
    // INITALIZE THE HISTOGRAMS
    // *****************************
    TH1F *phi000_eta[nDataSets+1][nt1bins]; // +1 is for the total phase 1 dataset
    TH1F *phi045_eta[nDataSets+1][nt1bins]; 
    TH1F *phi090_eta[nDataSets+1][nt1bins]; 
    TH1F *phi135_eta[nDataSets+1][nt1bins]; 
    TH1F *phiAMO_eta[nDataSets+1][nt1bins]; 
    TH1F *phi000_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi045_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi090_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi135_pi0[nDataSets+1][nt1bins]; 
    TH1F *phiAMO_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi000_090_eta[nDataSets+1][nt1bins]; 
    TH1F *phi045_135_eta[nDataSets+1][nt1bins]; 
    TH1F *phi000_090_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi045_135_pi0[nDataSets+1][nt1bins]; 
    cout << "Defined phi histograms" << endl;
    
    float minphi=-180;
    float maxphi=180;
    float nphi=30;
    int binphi=9;
    for (int iData=0; iData<nDataSets; ++iData){
        for (int it1bin=0; it1bin<nt1bins; ++it1bin){
            string yaxis="Entries / "+to_string(binphi)+" degrees";
            phi000_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi000_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi045_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi045_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi); 
            phi090_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi090_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi135_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi135_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phiAMO_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phiAMO_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi000_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi000_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi045_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi045_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi090_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi090_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi135_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi135_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phiAMO_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phiAMO_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
        }
    }
    cout << "Initialized phi histograms" << endl;
    
    ////////////////////////////////////
    // FILL THE HISTOGRAMS
    ////////////////////////////////////
    float weight;
    float selectionVariable;
    float teta;
    float tpi0;
    bool selected;
    int BeamAngle;
    float total_selected=0;
    int iteta;
    int itpi0;
    float p0; // parameter for shiftedCosine 
    int counter;
    float total_entries=0;
    for (int iData=0; iData<nDataSets; ++iData){
        cout << "Filling phi histograms for: " << dataSetTag[iData] << endl;
        Long64_t nentries=(Long64_t)array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"].size();
        total_entries+=nentries;
        for (Long64_t ientry=0; ientry<nentries; ++ientry){
            // Load the weight first so we can keep a tally of the total integral before selections
            weight=array_variable_map[dataSetTag[iData]+"_"+weightVar][ientry];
            total_entries+=weight;

            selected=true;
            // These skip bools will be used to skip filling (and fitting since we have a min entry requirement) depending
            //    on the variable we are binning in (the first vector in vector<criteria> selections)
            bool skipPi0Hists=false; 
            bool skipEtaHists=false;
            if (selections[0].variable=="Mpi0p")
                skipPi0Hists=true;
            if (selections[0].variable=="Metap")
                skipEtaHists=true;
            for (auto selection: selections){
                selectionVariable=array_variable_map[dataSetTag[iData]+"_"+selection.variable][ientry]; 
                if(selectionVariable<=0.0){
                    if (verbose) cout << "Kinematic binning variable is negative! We only expect positive variables. Exiting..." << endl;
                    exit(1);
                }
                selected *= ((selectionVariable>selection.minval1)*(selectionVariable<selection.maxval1)) || 
                            ((selectionVariable>selection.minval2)*(selectionVariable<selection.maxval2));
            } 
            if(selected){
                BeamAngle=array_BeamAngles_map[dataSetTag[iData]][ientry];
                teta=array_variable_map[dataSetTag[iData]+"_mandelstam_teta"][ientry];
                tpi0=array_variable_map[dataSetTag[iData]+"_mandelstam_tpi0"][ientry];
                //cout << BeamAngle << " " << teta << " " << tpi0 << endl;
             
                iteta=int(ceil(teta/(1.0/nt1bins))-1);
                itpi0=int(ceil(tpi0/(1.0/nt1bins))-1);   

                if(BeamAngle==0){
                    if ((iteta<nt1bins)*!skipEtaHists)
                        phi000_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if ((itpi0<nt1bins)*!skipPi0Hists)
                        phi000_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==45){
                    if ((iteta<nt1bins)*!skipEtaHists)
                        phi045_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if ((itpi0<nt1bins)*!skipPi0Hists)
                        phi045_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==90){
                    if ((iteta<nt1bins)*!skipEtaHists)
                        phi090_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if ((itpi0<nt1bins)*!skipPi0Hists)
                        phi090_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==135){
                    if ((iteta<nt1bins)*!skipEtaHists)
                        phi135_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if ((itpi0<nt1bins)*!skipPi0Hists)
                        phi135_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==-1){
                    if ((iteta<nt1bins)*!skipEtaHists)
                        phiAMO_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if ((itpi0<nt1bins)*!skipPi0Hists)
                        phiAMO_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                total_selected+=weight;
            }
        }
    }
    cout << "Total weighted entries post default selections: " << total_entries << endl;
    cout << " Total weighted entries after passing additional selections: " << total_selected << endl;
    
    TF1 *fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
    TF1 *fit_sc = new TF1("fit_flat",shiftedCos,-180,180,numDOFsig_sc); 
    fit_asym->SetLineColor(kRed);
    fit_sc->SetLineColor(kRed);
    //////////////////////////////////////
    // SCALE THE HISTOGRAMS BY FLUX RATIO
    //////////////////////////////////////
    for (int it1bin=0; it1bin<nt1bins; ++it1bin){
        phi000_eta[3][it1bin] = new TH1F(("total_phi000_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phi045_eta[3][it1bin] = new TH1F(("total_phi045_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi); 
        phi090_eta[3][it1bin] = new TH1F(("total_phi090_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phi135_eta[3][it1bin] = new TH1F(("total_phi135_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phiAMO_eta[3][it1bin] = new TH1F(("total_phiAMO_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phi000_pi0[3][it1bin] = new TH1F(("total_phi000_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phi045_pi0[3][it1bin] = new TH1F(("total_phi045_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phi090_pi0[3][it1bin] = new TH1F(("total_phi090_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phi135_pi0[3][it1bin] = new TH1F(("total_phi135_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);
        phiAMO_pi0[3][it1bin] = new TH1F(("total_phiAMO_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", nphi, minphi, maxphi);

        for (int iData=0; iData<nDataSets; ++iData){
            // First we scale all the para orientations with the flux ratios
            if (verbose) { cout << "(" << dataSetTag[iData] << " 000 eta) entries/maximum before scaling: " 
                << phi000_eta[iData][it1bin]->GetEntries() << "/" << phi000_eta[iData][it1bin]->GetMaximum() << endl;}
            phi000_eta[iData][it1bin]->Scale(fluxRatios_90_0[iData]);
            phi000_pi0[iData][it1bin]->Scale(fluxRatios_90_0[iData]);
            phi135_eta[iData][it1bin]->Scale(fluxRatios_45_135[iData]);
            phi135_pi0[iData][it1bin]->Scale(fluxRatios_45_135[iData]);
            if (verbose){ cout << "(" << dataSetTag[iData] << " 000 eta) entries/maximum after scaling: " << \
                phi000_eta[iData][it1bin]->GetEntries() << "/" << phi000_eta[iData][it1bin]->GetMaximum() << endl;}
        
            // Add up all the datasets for the total
            phi000_eta[3][it1bin]->Add(phi000_eta[iData][it1bin]);
            phi045_eta[3][it1bin]->Add(phi045_eta[iData][it1bin]); 
            phi090_eta[3][it1bin]->Add(phi090_eta[iData][it1bin]);
            phi135_eta[3][it1bin]->Add(phi135_eta[iData][it1bin]);
            phiAMO_eta[3][it1bin]->Add(phiAMO_eta[iData][it1bin]);
            phi000_pi0[3][it1bin]->Add(phi000_pi0[iData][it1bin]);
            phi045_pi0[3][it1bin]->Add(phi045_pi0[iData][it1bin]);
            phi090_pi0[3][it1bin]->Add(phi090_pi0[iData][it1bin]);
            phi135_pi0[3][it1bin]->Add(phi135_pi0[iData][it1bin]);
            phiAMO_pi0[3][it1bin]->Add(phiAMO_pi0[iData][it1bin]);
        }
        if (verbose) cout << "Phase 1 000 eta entries " << phi000_eta[3][it1bin]->GetEntries() << endl;
        
        //////////////////////////////////
        // FIT PRELIMINARIES
        //////////////////////////////////
        TCanvas *allCanvases = new TCanvas("","",1440,600);
        
        TFitResultPtr fitPointer;
        float asymmetries_000_eta;
        float asymmetries_000_eta_err;
        float asymmetries_045_eta;
        float asymmetries_045_eta_err;
        float asymmetries_000_pi0;
        float asymmetries_000_pi0_err;
        float asymmetries_045_pi0;
        float asymmetries_045_pi0_err;
        float pol_frac_000; 
        float pol_frac_090; 
        float pol_frac_045; 
        float pol_frac_135; 
        float asymmetry000_090_eta_chiPerDOF;
        float asymmetry045_135_eta_chiPerDOF;
        float asymmetry000_090_pi0_chiPerDOF;
        float asymmetry045_135_pi0_chiPerDOF;
        for(int iData=0; iData<4; ++iData){
            pol_frac_000 = polfractions[iData][0];//+-polfractions_err[iData][0]; 
            pol_frac_045 = polfractions[iData][1];//+-polfractions_err[iData][1]; 
            pol_frac_090 = polfractions[iData][2];//+-polfractions_err[iData][2]; 
            pol_frac_135 = polfractions[iData][3];//+-polfractions_err[iData][3]; 
            poloffset[4]=Phi0_0_90[iData];
            poloffset[5]=Phi0_45_135[iData];
            // *****************************
            // GETTING ASYMMETRIES
            // *****************************
            TH1* asymmetry000_090_eta = phi090_eta[iData][it1bin]->GetAsymmetry(phi000_eta[iData][it1bin]);
            TH1* asymmetry045_135_eta = phi045_eta[iData][it1bin]->GetAsymmetry(phi135_eta[iData][it1bin]);
            asymmetry000_090_eta->SetTitle("0/90 Asymmetry");
            asymmetry045_135_eta->SetTitle("45/135 Asymmetry");
            asymmetry000_090_eta->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            asymmetry045_135_eta->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            TH1* asymmetry000_090_pi0 = phi090_pi0[iData][it1bin]->GetAsymmetry(phi000_pi0[iData][it1bin]);
            TH1* asymmetry045_135_pi0 = phi045_pi0[iData][it1bin]->GetAsymmetry(phi135_pi0[iData][it1bin]);
            asymmetry000_090_pi0->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            asymmetry045_135_pi0->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            asymmetry000_090_pi0->SetTitle("0/90 Asymmetry");
            asymmetry045_135_pi0->SetTitle("45/135 Asymmetry");

            // Constructing sums of paired polarizations for systematics
            //   the perpendicular clones will be thrown away later after we sum them with 0 pols
            //   Scale everything by the cross polarization magnitude
            //   Cross weight sum the orientations
            phi000_090_eta[iData][it1bin]=(TH1F*)phi000_eta[iData][it1bin]->Clone("crossWeighted000_eta"); 
            phi045_135_eta[iData][it1bin]=(TH1F*)phi045_eta[iData][it1bin]->Clone("crossWeighted045_eta");
            phi000_090_pi0[iData][it1bin]=(TH1F*)phi000_pi0[iData][it1bin]->Clone("crossWeighted000_pi0");
            phi045_135_pi0[iData][it1bin]=(TH1F*)phi045_pi0[iData][it1bin]->Clone("crossWeighted045_pi0");
            TH1F* crossWeighted090_eta=(TH1F*)phi090_eta[iData][it1bin]->Clone("crossWeighted090_eta"); 
            TH1F* crossWeighted135_eta=(TH1F*)phi135_eta[iData][it1bin]->Clone("crossWeighted135_eta");
            TH1F* crossWeighted090_pi0=(TH1F*)phi090_pi0[iData][it1bin]->Clone("crossWeighted090_pi0");
            TH1F* crossWeighted135_pi0=(TH1F*)phi135_pi0[iData][it1bin]->Clone("crossWeighted135_pi0");
            phi000_090_eta[iData][it1bin]->Scale(pol_frac_090); 
            phi045_135_eta[iData][it1bin]->Scale(pol_frac_135);
            phi000_090_pi0[iData][it1bin]->Scale(pol_frac_090);
            phi045_135_pi0[iData][it1bin]->Scale(pol_frac_135);
            crossWeighted090_eta->Scale(pol_frac_000);
            crossWeighted135_eta->Scale(pol_frac_045);
            crossWeighted090_pi0->Scale(pol_frac_000);
            crossWeighted135_pi0->Scale(pol_frac_045);
            phi000_090_eta[iData][it1bin]->Add(crossWeighted090_eta);
            phi045_135_eta[iData][it1bin]->Add(crossWeighted135_eta);
            phi000_090_pi0[iData][it1bin]->Add(crossWeighted090_pi0);
            phi045_135_pi0[iData][it1bin]->Add(crossWeighted135_pi0);

            // Loading phi distributions into a vector
            TH1F *phis_eta[7] =  {phi000_eta[iData][it1bin], phi045_eta[iData][it1bin], phi090_eta[iData][it1bin], 
                                    phi135_eta[iData][it1bin], phi000_090_eta[iData][it1bin], phi045_135_eta[iData][it1bin], phiAMO_eta[iData][it1bin]} ;
            TH1F *phis_pi0[7] =  {phi000_pi0[iData][it1bin], phi045_pi0[iData][it1bin], phi090_pi0[iData][it1bin], 
                                    phi135_pi0[iData][it1bin], phi000_090_pi0[iData][it1bin], phi045_135_pi0[iData][it1bin], phiAMO_pi0[iData][it1bin]};
            vector<float> eta_psig;
            vector<float> eta_psig_err;
            vector<float> pi0_psig;
            vector<float> pi0_psig_err;

        
            // *****************************
            // Fitting asymmetry for eta
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(2,1);
            allCanvases->cd(1);
            fit_asym->SetParameters(pol_frac_090,pol_frac_000,0,Phi0_0_90[iData]);
            fit_asym->FixParameter(0,pol_frac_090);
            fit_asym->FixParameter(1,pol_frac_000);
            fit_asym->SetParLimits(2,-1,1);
            if (freePhase)
                fit_asym->SetParLimits(3,Phi0_0_90[iData]-10,Phi0_0_90[iData]+10);
            else
                fit_asym->FixParameter(3,Phi0_0_90[iData]);

            if ((phi090_eta[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)*(phi000_eta[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)){
                for (int i=0; i<maxiter; i++){
                    fit_asym->SetParameter(2,gRandom->Uniform()-0.5);
                    fitPointer = asymmetry000_090_eta->Fit(fit_asym,fitOption.c_str());
                    if ( ((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))//*(abs(fit_asym->GetParError(2))<1) )
                        break;
                    else{
                        if (verbose){ cout << "  Attempt " << i+1 << " failed for asym 0/90 eta... status, asym, asymerr: " << 
                            (int)fitPointer << ", " << fit_asym->GetParameter(2) << ", " << abs(fit_asym->GetParError(2)) << endl;}
                    }
                }
                if ( !(((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))){//*(abs(fit_asym->GetParError(2))<1)) ){
                    if (verbose) cout << "Asym fit 0/90 for eta was not successful after maxiter=" << maxiter << " times! exiting..." << endl;
                    exit(1);
                }
                asymmetries_000_eta = fit_asym->GetParameter(2);
                asymmetries_000_eta_err = fit_asym->GetParError(2);
                asymmetry000_090_eta_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_000_eta = 0;
                asymmetries_000_eta_err = -999;
                asymmetry000_090_eta_chiPerDOF = -999;
            }
            asymmetry000_090_eta->Draw("SAME");
            if (overlayBinInfoIntoAsymPlots)
                text->DrawTextNDC(0.12, 0.14, (sminval1+"<"+selections[0].variable+"<"+smaxval1).c_str());
            asymmetry000_090_eta->SetAxisRange(-1,1,"Y");
        
            allCanvases->cd(2);
            fit_asym->SetParameters(pol_frac_135,pol_frac_045,0.0,Phi0_45_135[iData]);
            fit_asym->FixParameter(0,pol_frac_135);
            fit_asym->FixParameter(1,pol_frac_045);
            fit_asym->SetParLimits(2,-1,1);
            if (freePhase)
                fit_asym->SetParLimits(3,Phi0_45_135[iData]-10,Phi0_45_135[iData]+10);
            else
                fit_asym->FixParameter(3,Phi0_45_135[iData]);
            if ((phi135_eta[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)*(phi045_eta[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)){
                for (int i=0; i<maxiter; i++){
                    fit_asym->SetParameter(2,gRandom->Uniform()-0.5);
                    fitPointer = asymmetry045_135_eta->Fit(fit_asym,fitOption.c_str());
                    if ( ((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))//*(abs(fit_asym->GetParError(2))<1) )
                        break;
                    else{
                        if (verbose){ cout << "  Attempt " << i+1 << " failed for asym 45/-45 eta... status, asym, asymerr: " << 
                            (int)fitPointer << ", " << fit_asym->GetParameter(2) << ", " << abs(fit_asym->GetParError(2)) << endl;}
                    }
                }
                if ( !(((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))){//*(abs(fit_asym->GetParError(2))<1)) ){
                    if (verbose) cout << "Asym fit 45/-45 for eta was not successful after maxiter=" << maxiter << " times! exiting..." << endl;
                    exit(1);
                }
                asymmetries_045_eta = fit_asym->GetParameter(2);
                asymmetries_045_eta_err = fit_asym->GetParError(2);
                asymmetry045_135_eta_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_045_eta = 0;
                asymmetries_045_eta_err = -999;
                asymmetry045_135_eta_chiPerDOF = -999;
            }
            asymmetry045_135_eta->Draw("SAME");
            if (overlayBinInfoIntoAsymPlots)
                text->DrawTextNDC(0.12, 0.14, (sminval1+"<"+selections[0].variable+"<"+smaxval1).c_str());
            asymmetry045_135_eta->SetAxisRange(-1,1,"Y");
            allCanvases->SaveAs((fitSaveLoc+"_eta"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"."+fileType).c_str()); 

            // *****************************
            // Fitting shifted cos - phi for eta
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(4,2);
            counter=0;
            for (auto phi: phis_eta){
            	allCanvases->cd(counter+1);
            	phi->SetTitle((polnames[counter]).c_str());
            	p0 = phi->Integral()/phi->GetNbinsX();
                if (phi->GetEntries()>minEntriesInPhiHists){
                    for (int i=0; i<maxiter; i++){
                        if (counter<6){
            	            fit_sc->SetParameters(p0*(gRandom->Uniform()*2),gRandom->Uniform()-0.5,poloffset[counter]); // p1 should be PSigma
                            if (freePhase)
                                fit_sc->SetParLimits(2,poloffset[counter]-10,poloffset[counter]+10);
                            else
                                fit_sc->FixParameter(2,poloffset[counter]);
                        }
                        else{ // index 6 corresponds to AMO dataset
                            fit_sc->SetParameters(p0*(gRandom->Uniform()*2),gRandom->Uniform()-0.5,0);
                            //fit_sc->SetParLimits(2,-180,180); // NOT SURE HOW TO HANDLE AMO BUT FREEING PHI LEADS TO SOME FAILED FITS
                            fit_sc->FixParameter(2,poloffset[counter]);
                        }
            	        fitPointer = phi->Fit(fit_sc,fitOption.c_str());
                        if ( fitPointer==0 ){
                            if (verbose) cout << "Shifted cos fit " << polnames[counter] << " for eta was successful!" << endl;
                            break;
                        }
                        else{
                            if (verbose){ cout << "  Attempt " << i << " shifted cos fit " << polnames[counter] << " for eta was not successful! Contained " 
                                << phi->GetEntries() << " entries with integral " << phi->Integral() << endl;}
                        }
                    }
                    if (fitPointer!=0){
                        if (verbose){ cout << "Shifted cos fit " << polnames[counter] << " for eta was not successful after maxiter=" << maxiter 
                            << " times!"; }
                        if (fitPhi){
                            cout << " exiting..." << endl;
                            exit(1);
                        }
                        else{
                            cout << " You requested not to exit, so we will continue running" << endl;}
                    }
            	    eta_psig.push_back(fit_sc->GetParameter(1));
            	    eta_psig_err.push_back(fit_sc->GetParError(1));
            	    phi->Draw("SAME");
                }
                else{
                    if (verbose) cout << "Shifted cos fit " << polnames[counter] << " for eta does not have enough events to safely fit..." << endl;
            	    eta_psig.push_back(0);
            	    eta_psig_err.push_back(-999);
                }
            	++counter;
            }
            allCanvases->SaveAs((fitSaveLoc+"_eta"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"_phis."+fileType).c_str()); 
            
            // *****************************
            // Fitting asymmetry for pi0
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(2,1);
            allCanvases->cd(1);
            fit_asym->SetParameters(pol_frac_090,pol_frac_000,0.0,Phi0_0_90[iData]);
            fit_asym->FixParameter(0,pol_frac_090);
            fit_asym->FixParameter(1,pol_frac_000);
            fit_asym->SetParLimits(2,-1,1);
            if (freePhase)
                fit_asym->SetParLimits(3,Phi0_0_90[iData]-10,Phi0_0_90[iData]+10);
            else
                fit_asym->FixParameter(3,Phi0_0_90[iData]);
            if ((phi090_pi0[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)*(phi000_pi0[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)){
                for (int i=0; i<maxiter; i++){
                    fit_asym->SetParameter(2,gRandom->Uniform()-0.5);
                    fitPointer = asymmetry000_090_pi0->Fit(fit_asym,fitOption.c_str());
                    if ( ((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))//*(abs(fit_asym->GetParError(2))<1) )
                        break;
                    else{
                        if (verbose){ cout << "  Attempt " << i+1 << " failed for asym 0/90 pi0... status, asym, asymerr: " << 
                            (int)fitPointer << ", " << fit_asym->GetParameter(2) << ", " << abs(fit_asym->GetParError(2)) << endl;}
                    }
                }
                if ( !(((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))){//*(abs(fit_asym->GetParError(2))<1)) ){
                    if (verbose) cout << "Asym fit 0/90 for pi0 was not successful after maxiter=" << maxiter << " times! exiting..." << endl;
                    exit(1);
                }
                asymmetries_000_pi0 = fit_asym->GetParameter(2);
                asymmetries_000_pi0_err = fit_asym->GetParError(2);
                asymmetry000_090_pi0_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_000_pi0 = 0;
                asymmetries_000_pi0_err = -999;
                asymmetry000_090_pi0_chiPerDOF = -999;
            }
            asymmetry000_090_pi0->Draw("SAME");
            if (overlayBinInfoIntoAsymPlots)
                text->DrawTextNDC(0.12, 0.14, (sminval1+"<"+selections[0].variable+"<"+smaxval1).c_str());
            asymmetry000_090_pi0->SetAxisRange(-1,1,"Y");
            
            allCanvases->cd(2);
            fit_asym->SetParameters(pol_frac_135,pol_frac_045,0.0,Phi0_45_135[iData]);
            fit_asym->FixParameter(0,pol_frac_135);
            fit_asym->FixParameter(1,pol_frac_045);
            fit_asym->SetParLimits(2,-1,1);
            if (freePhase)
                fit_asym->SetParLimits(3,Phi0_45_135[iData]-10,Phi0_45_135[iData]+10);
            else
                fit_asym->FixParameter(3,Phi0_45_135[iData]);
            if ((phi135_pi0[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)*(phi045_pi0[iData][it1bin]->GetEntries()>minEntriesForAsymPhiHists)){
                for (int i=0; i<maxiter; i++){
                    fit_asym->SetParameter(2,gRandom->Uniform()-0.5);
                    fitPointer = asymmetry045_135_pi0->Fit(fit_asym,fitOption.c_str());
                    if ( ((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))//*(abs(fit_asym->GetParError(2))<1) )
                        break;
                    else{
                        if (verbose){ cout << "  Attempt " << i+1 << " failed for asym -45/45 pi0... status, asym, asymerr: " << 
                            (int)fitPointer << ", " << fit_asym->GetParameter(2) << ", " << abs(fit_asym->GetParError(2)) << endl;}
                    }
                }
                if ( !(((int)fitPointer==0)*(fit_asym->GetParameter(2)!=1))){//*(abs(fit_asym->GetParError(2))<1)) ){
                    if (verbose) cout << "Asym fit 45/-45 for pi0 was not successful after maxiter=" << maxiter << " times! exiting..." << endl;
                    exit(1);
                }
                asymmetries_045_pi0 = fit_asym->GetParameter(2);
                asymmetries_045_pi0_err = fit_asym->GetParError(2);
                asymmetry045_135_pi0_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_045_pi0 = 0; 
                asymmetries_045_pi0_err = -999; 
                asymmetry045_135_pi0_chiPerDOF = -999;
            }
            asymmetry045_135_pi0->Draw("SAME");
            if (overlayBinInfoIntoAsymPlots)
                text->DrawTextNDC(0.12, 0.14, (sminval1+"<"+selections[0].variable+"<"+smaxval1).c_str());
            asymmetry045_135_pi0->SetAxisRange(-1,1,"Y");
            allCanvases->SaveAs((fitSaveLoc+"_pi0"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"."+fileType).c_str()); 

            // *****************************
            // Fitting shifted cos - phi for pi0
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(4,2);
            counter=0;
            for (auto phi: phis_pi0){
            	allCanvases->cd(counter+1);
            	phi->SetTitle((polnames[counter]).c_str());
            	p0 = phi->Integral()/phi->GetNbinsX();
                if (phi->GetEntries()>minEntriesInPhiHists){
                    for (int i=0; i<maxiter; i++){
                        if (counter<6){
            	            fit_sc->SetParameters(p0*(gRandom->Uniform()*2),gRandom->Uniform()-0.5,poloffset[counter]); // p1 is PSigma
                            if (freePhase)
                                fit_sc->SetParLimits(2,poloffset[counter]-10,poloffset[counter]+10);
                            else
                                fit_sc->FixParameter(2,poloffset[counter]);
                        }
                        else{ // index 6 corresponds to AMO dataset
                            fit_sc->SetParameters(p0*(gRandom->Uniform()*2),gRandom->Uniform()-0.5,0);
                            //fit_sc->SetParLimits(2,-180,180); // NOT SURE HOW TO HANDLE AMO BUT FREEING PHI LEADS TO SOME FAILED FITS
                            fit_sc->FixParameter(2,0);
                        }
            	        fitPointer = phi->Fit(fit_sc,fitOption.c_str());
                        if ( fitPointer==0 ){
                            if (verbose) cout << "Shifted cos fit " << polnames[counter] << " for pi0 was successful!" << endl;
                            break;
                        }
                        else{
                            if (verbose){ cout << "  Attempt " << i << " shifted cos fit  " << polnames[counter] << " for pi0 was not successful! Contained " 
                                << phi->GetEntries() << " entries with integral " << phi->Integral() << endl;}
                        }
                    }
                    if (fitPointer!=0){
                        if (verbose){ cout << "Shifted cos fit " << polnames[counter] << " for pi0 was not successful after maxiter=" << maxiter 
                            << " times!"; }
                        if (fitPhi){
                            cout << " exiting..." << endl;
                            exit(1);
                        }
                        else{
                            cout << " You requested not to exit, so we will continue running" << endl;}
                    }
            	    pi0_psig.push_back(fit_sc->GetParameter(1));
            	    pi0_psig_err.push_back(fit_sc->GetParError(1));
            	    phi->Draw("SAME");
                }
                else{
                    if (verbose) cout << "Shifted cos fit for " << polnames[counter] << " pi0 does not have enough events to safely fit..." << endl;
            	    pi0_psig.push_back(0);
            	    pi0_psig_err.push_back(-999);
                }
            	++counter;
            }
            allCanvases->SaveAs((fitSaveLoc+"_pi0"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"_phis."+fileType).c_str()); 
        
            // ********************************
            // Writing out the results
            // ********************************
            //float midt1=0.2*it1bin+0.1;
            float midt1=t1HalfWidth*2*(it1bin+0.5);
            *saveCsv << selections[0].variable << " " << dataSetTag2[iData] << " "
                     << midt1 << " " << t1HalfWidth << " "
                     << asymmetries_000_eta << " " << asymmetries_000_eta_err << " " 
                     << asymmetries_045_eta << " " << asymmetries_045_eta_err << " "
                     << asymmetries_000_pi0 << " " << asymmetries_000_pi0_err << " " 
                     << asymmetries_045_pi0 << " " << asymmetries_045_pi0_err << " "
                     << asymmetry000_090_eta_chiPerDOF << " " << asymmetry045_135_eta_chiPerDOF << " "
                     << asymmetry000_090_pi0_chiPerDOF << " " << asymmetry045_135_pi0_chiPerDOF << " "
                     << phi000_eta[iData][it1bin]->GetEntries() << " " << phi045_eta[iData][it1bin]->GetEntries() << " "
                     << phi090_eta[iData][it1bin]->GetEntries() << " " << phi135_eta[iData][it1bin]->GetEntries() << " "
                     << phi000_pi0[iData][it1bin]->GetEntries() << " " << phi045_pi0[iData][it1bin]->GetEntries() << " "
                     << phi090_pi0[iData][it1bin]->GetEntries() << " " << phi135_pi0[iData][it1bin]->GetEntries();
            for (int i=0; i<(int)eta_psig.size(); ++i){
                *saveCsv << " " << eta_psig[i] << " " << eta_psig_err[i] << " " << pi0_psig[i] << " " << pi0_psig_err[i];
            }
            for (auto selection : selections)
                *saveCsv << " " << selection.variable << " " << selection.minval1 << " " << selection.maxval1;
                // Will only save the details of the first region. The loose/nominal/tighter selections are all distint enough 
            *saveCsv << endl;

            delete crossWeighted090_eta;
            delete crossWeighted135_eta;
            delete crossWeighted090_pi0;
            delete crossWeighted135_pi0;
        } // end data loop
    } // end t loop
    
    ///////////////////////////////////
    // CLEAN UP THE POINTERS FOR NEXT ITER
    ///////////////////////////////////
    for(int iData=0; iData<4; ++iData){
        for (int it1bin=0; it1bin<nt1bins; ++it1bin){
            delete phi000_eta[iData][it1bin];
            delete phi045_eta[iData][it1bin]; 
            delete phi090_eta[iData][it1bin];
            delete phi135_eta[iData][it1bin];
            delete phiAMO_eta[iData][it1bin];
            delete phi000_pi0[iData][it1bin];
            delete phi045_pi0[iData][it1bin];
            delete phi090_pi0[iData][it1bin];
            delete phi135_pi0[iData][it1bin];
            delete phiAMO_pi0[iData][it1bin];
        }
    }
}

map<string, float> extractAsymmetries(
        string sbTag, 
        vector<float> sbRegion, 
        constructFitArgs args,
        map<string,vector<pair<int,criteria>>> eventSelects
        ){
    /// The output of this program is "percetages"
    map<string, float> percentages; // Percentage remaining Post defSels compared to the entire tree

    float pi0_sig=sbRegion[0];  
    float pi0_skip=sbRegion[1]; 
    float pi0_sb=sbRegion[2];   
    float eta_sig=sbRegion[3];  
    float eta_skip=sbRegion[4]; 
    float eta_sb=sbRegion[5];   
    ///////////////////////////////////
    // Storage variables to set branch addresses to 
    ///////////////////////////////////
    float phi_eta_lab;
    float phi_pi0_lab;
    float Mpi0;
    float Meta;
    float Mpi0p;
    float Metap;
    float mandelstam_t;
    float mandelstam_tp;
    float mandelstam_teta;
    float mandelstam_tpi0;
    float Mpi0eta;
    float AccWeight;
    float weightASBS;
    int BeamAngle;
    float unusedEnergy;
    float chiSq;
    float mmsq;
    float photonTheta1;
    float photonTheta2;
    float photonTheta3;
    float photonTheta4;
    float photonE1;
    float photonE2;
    float photonE3;
    float photonE4;
    float proton_momentum;
    float proton_z;

    unordered_map<string, float*> variable_map={
        {"phi_eta_lab",&phi_eta_lab},
        {"phi_pi0_lab",&phi_pi0_lab},
        {"Mpi0",&Mpi0},
        {"Meta",&Meta},
        {"Mpi0p",&Mpi0p},
        {"Metap",&Metap},
        {"unusedEnergy",&unusedEnergy},
        {"chiSq",&chiSq},
        {"mmsq",&mmsq},
        {"photonTheta1",&photonTheta1},
        {"photonTheta2",&photonTheta2},
        {"photonTheta3",&photonTheta3},
        {"photonTheta4",&photonTheta4},
        {"photonE1",&photonE1},
        {"photonE2",&photonE2},
        {"photonE3",&photonE3},
        {"photonE4",&photonE4},
        {"proton_momentum",&proton_momentum},
        {"proton_z",&proton_z},
        {"mandelstam_t",&mandelstam_t},
        {"mandelstam_tp",&mandelstam_tp},
        {"mandelstam_teta",&mandelstam_teta},
        {"mandelstam_tpi0",&mandelstam_tpi0},
        {"Mpi0eta",&Mpi0eta},
        {"AccWeight",&AccWeight},
        {"weightASBS",&weightASBS},
    };

    vector<string> variables={"phi_eta_lab","phi_pi0_lab","Mpi0","Meta","Mpi0p","Metap",
        "unusedEnergy","chiSq","mmsq","proton_momentum","proton_z",
        "photonTheta1","photonTheta2","photonTheta3","photonTheta4",
        "photonE1","photonE2","photonE3","photonE4",
        "mandelstam_t","mandelstam_tp","mandelstam_teta","mandelstam_tpi0","Mpi0eta","AccWeight","weightASBS"};
    unordered_map<string, vector<float>> array_variable_map_all; // Contains 'all' the values as opposed to array_variable_map, which contains selectiond
    unordered_map<string, vector<int>> array_BeamAngles_map_all;

    // **************************************************
    // Define eta and pi0 signal and sideband regions
    // **************************************************
    float pi0sigL=pi0_peak-pi0_sig*pi0_std;
    float pi0sigR=pi0_peak+pi0_sig*pi0_std;
    float pi0sbRL=pi0_peak+(pi0_sig+pi0_skip)*pi0_std;
    float pi0sbLR=pi0_peak-(pi0_sig+pi0_skip)*pi0_std;
    float pi0sbRR=pi0_peak+(pi0_sig+pi0_skip+pi0_sb)*pi0_std;
    float pi0sbLL=pi0_peak-(pi0_sig+pi0_skip+pi0_sb)*pi0_std;
    float etasigL=eta_peak-eta_sig*eta_std;
    float etasigR=eta_peak+eta_sig*eta_std;
    float etasbRL=eta_peak+(eta_sig+eta_skip)*eta_std;
    float etasbLR=eta_peak-(eta_sig+eta_skip)*eta_std;
    float etasbRR=eta_peak+(eta_sig+eta_skip+eta_sb)*eta_std;
    float etasbLL=eta_peak-(eta_sig+eta_skip+eta_sb)*eta_std;

    float eta_sbweight=-1.0*eta_sig/eta_sb;
    float pi0_sbweight=-1.0*pi0_sig/pi0_sb;

    float weightBS;
    float weightBSpi0;
    float weightBSeta;

    // **************************************************
    // LOAD ALL THE DATA INTO MEMORY TO SAVE ON FILE I/O
    // **************************************************
    TTree* tree;
    vector<Long64_t> nentries_all;
    for (int iData=0; iData<3; ++iData){
        string dataFileName=folder+"D"+dataSetTag[iData]+"_selected_acc_flat.root";
    	cout << "LOADING ROOT FILE: " << dataFileName << endl; 
    	TFile *dataFile = new TFile(dataFileName.c_str());
        dataFile->GetObject("kin",tree);
        tree->SetBranchStatus("*",0); // disable all branches and re-enable only the ones we want to look at
        //nentries_all.push_back(10000); 
        nentries_all.push_back(tree->GetEntries());
        cout << dataSetTag[iData] << " nentries: " << nentries_all[iData] << endl;
        
        for(auto variable: variables){
            array_variable_map_all[dataSetTag[iData]+"_"+variable]=vector<float>{};
            array_variable_map_all[dataSetTag[iData]+"_"+variable].reserve(nentries_all[iData]);
            tree->SetBranchStatus(variable.c_str(),1);
            tree->SetBranchAddress(variable.c_str(),variable_map[variable]);
        }
        tree->SetBranchStatus("BeamAngle",1);
        tree->SetBranchAddress("BeamAngle",&BeamAngle);
    
        auto start = std::chrono::system_clock::now();
        for (Long64_t ientry=0; ientry<nentries_all[iData]; ++ientry){
            tree->GetEntry(ientry);
            if ( (ientry%50000)==0 ){
                auto end=std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed=end-start;
                cout << ientry << "/" << nentries_all[iData] << "  elapsed time: " << elapsed.count() << "s" << endl;
            }

            //////////////////////
            /// DETERMINE WEIGHTS
            //////////////////////
            if ( Mpi0 > pi0_peak-pi0_std*pi0_sig && Mpi0 < pi0_peak+pi0_std*pi0_sig ) { weightBSpi0 = 1; } 
            else if ( Mpi0 > pi0_peak-pi0_std*(pi0_sig+pi0_skip+pi0_sb) && Mpi0 < pi0_peak-pi0_std*(pi0_sig+pi0_skip) ) { weightBSpi0 = pi0_sbweight; } 
            else if ( Mpi0 > pi0_peak+pi0_std*(pi0_sig+pi0_skip) && Mpi0 < pi0_peak+pi0_std*(pi0_sig+pi0_skip+pi0_sb) ) { weightBSpi0 = pi0_sbweight; } 
            else { weightBSpi0 = 0; }
            if ( Meta > eta_peak-eta_std*eta_sig && Meta < eta_peak+eta_std*eta_sig ) { weightBSeta = 1; } 
            else if ( Meta > eta_peak-eta_std*(eta_sig+eta_skip+eta_sb) && Meta < eta_peak-eta_std*(eta_sig+eta_skip) ) { weightBSeta = eta_sbweight; } 
            else if ( Meta > eta_peak+eta_std*(eta_sig+eta_skip) && Meta < eta_peak+eta_std*(eta_sig+eta_skip+eta_sb) ) { weightBSeta = eta_sbweight; } 
            else { weightBSeta = 0; }
            //weightBS=weightBSpi0*weightBSeta;
            weightBS=weightBSeta;
            //if ( abs(weightASBS-weightBS*AccWeight)>0.00001 ){ // CHECK IF WEIGHTING SCHEME WORKS AS EXPECTED BY COMPARING TO THE WEIGHTS IN THE TREE
            //    cout << "WEIGHTS NOT MATCHING!" << endl;
            //    exit(1);
            //}

            // Do not select on Mpi0 nor Meta if we are going to do sideband subtraction also
            for(string variable: variables){
                if ((variable==args.weightVar)*(args.weightVar!="AccWeight")) // overwrite the weights variable with the calculated
                    array_variable_map_all[dataSetTag[iData]+"_"+variable].push_back(AccWeight*weightBS); 
                else
                    array_variable_map_all[dataSetTag[iData]+"_"+variable].push_back(*variable_map[variable]); 
            }
            array_BeamAngles_map_all[dataSetTag[iData]].push_back(BeamAngle);
        }
    }

    // *************************************
    // Loading the data in from the histogram 
    // *************************************
    TH1F* dHist_Mpi0 = new TH1F("dHist_Mpi0","dHist_Mpi0;M(#gamma_{1}#gamma_{2}) GeV;Entries / 0.001",200,0.05,0.25);
    TH1F* dHist_Meta = new TH1F("dHist_Meta","dHist_Meta;M(#gamma_{3}#gamma_{4}) GeV;Entries / 0.002",300,0.25,0.85);
    TH1F* dHist_unusedEnergy = new TH1F("dHist_unusedEnergy","Unused Energy (GeV)",200,0,0.3);
    TH1F* dHist_chiSq = new TH1F("chiSq", "chiSq", 200,0,20);
    TH1F* dHist_photonTheta1 = new TH1F("photonTheta1", "Photon theta 1", 1500, 0, 15);
    TH1F* dHist_photonTheta2 = new TH1F("photonTheta2", "Photon theta 2", 1500, 0, 15);
    TH1F* dHist_photonTheta3 = new TH1F("photonTheta3", "Photon theta 3", 1500, 0, 15);
    TH1F* dHist_photonTheta4 = new TH1F("photonTheta4", "Photon theta 4", 1500, 0, 15);
    TH1F* dHist_photonE1= new TH1F("photonE1","photon E 1", 200, 0, 0.5);
    TH1F* dHist_photonE2= new TH1F("photonE2","photon E 2", 200, 0, 0.5);
    TH1F* dHist_photonE3= new TH1F("photonE3","photon E 3", 200, 0, 0.5);
    TH1F* dHist_photonE4= new TH1F("photonE4","photon E 4", 200, 0, 0.5);
    TH1F* dHist_protonMomentum = new TH1F("proton_momentum", "proton momentum", 200, 0, 1);
    TH1F* dHist_protonZ = new TH1F("proton_z", "proton z", 200, 45, 85);
    TH1F* dHist_mmsq = new TH1F("mmsq", "missing mass sq", 200, -0.1,0.1);
    unordered_map<string, TH1F*> diagnosticHists={
        {"unusedEnergy",dHist_unusedEnergy},
        {"chiSq",dHist_chiSq},
        {"photonTheta1",dHist_photonTheta1},
        {"photonTheta2",dHist_photonTheta2},
        {"photonTheta3",dHist_photonTheta3},
        {"photonTheta4",dHist_photonTheta4},
        {"photonE1",dHist_photonE1},
        {"photonE2",dHist_photonE2},
        {"photonE3",dHist_photonE3},
        {"photonE4",dHist_photonE4},
        {"proton_momentum",dHist_protonMomentum},
        {"proton_z",dHist_protonZ},
        {"mmsq",dHist_mmsq}
    };
    TBox* box = new TBox();
    TCanvas* c1 = new TCanvas("c1","c1",1440,600);
    for (auto const& ele: eventSelects){
        // //////////////////////////
        // CLEAN/RESET ALL HISTOGRAMS
        // //////////////////////////
        dHist_Mpi0->Reset();
        dHist_Meta->Reset();
        for (auto const& mapElement: diagnosticHists){
            mapElement.second->Reset();
        }
        c1->Clear();

        // //////////////////////////
        // SETUP THE SELECTIONS
        // //////////////////////////
        string evtTag=ele.first;
        defaultSelections defSels;
        if (evtTag!="default"){ // modify defaultSelections if requested
            vector<pair<int,criteria>> eventSelect=ele.second;
            for (auto eventSel: eventSelect){
                defSels.selections[eventSel.first]=eventSel.second;
            }
        }
        else{ evtTag=""; }// For the default case, we will not include a tag to the save file names

        cout << "\n\n=========================================================================================" << endl;
        cout << "=========================================================================================" << endl;
        if (evtTag!="")
            cout << " * * * * * * * * " << evtTag << " * * * * * * * * " << endl;
        else
            cout << " * * * * * * * * Default * * * * * * * * " << endl;
        cout << "=========================================================================================" << endl;
        cout << "=========================================================================================" << endl;
        cout << "Event selections:" << endl;
        for (auto eventSel: defSels.selections){
            cout << eventSel.minval1 << " < " <<eventSel.variable << " < " << eventSel.maxval1 << " OR "; 
            cout << eventSel.minval2 << " < " <<eventSel.variable << " < " << eventSel.maxval2 << endl; 
        }

        unordered_map<string, vector<float>> array_variable_map;
        unordered_map<string, vector<int>> array_BeamAngles_map;
        double signalIntegral=0;
        float total_nentries_postDefaultSelects=0;
        float total_nentries_preDefaultSelects=0;
        float weight;
        float accidental;
        for (int iData=0; iData<3; ++iData){
            cout << "Filtering Dataset: " << dataSetTag[iData] << endl;
            for(auto variable: variables){
                array_variable_map[dataSetTag[iData]+"_"+variable]=vector<float>{};
                array_variable_map[dataSetTag[iData]+"_"+variable].reserve(nentries_all[iData]);
            }
            float nentries_postDefaultSelects=0;
            float nentries_preDefaultSelects=0;
            
            bool selectPi0;
            bool selectEta;
            bool selected;
            float selectionVariable;
            auto start = std::chrono::system_clock::now();
            for (Long64_t ientry=0; ientry<nentries_all[iData]; ++ientry){
                if ( (ientry%100000)==0 ){
                    auto end=std::chrono::system_clock::now();
                    std::chrono::duration<double> elapsed=end-start;
                    cout << ientry << "/" << nentries_all[iData] << "  elapsed time: " << elapsed.count() << "s" << endl;
                }
                //////////////////////
                /// CHECK DEFAULT EVENT SELECTIONS 
                //////////////////////
                selected=1;
                for (auto eventSel: defSels.selections){
                    selectionVariable=array_variable_map_all[dataSetTag[iData]+"_"+eventSel.variable][ientry]; 
                    selected *= ((selectionVariable>eventSel.minval1)*(selectionVariable<eventSel.maxval1)) || 
                                ((selectionVariable>eventSel.minval2)*(selectionVariable<eventSel.maxval2));
//                    cout << eventSel.minval1 << " < " << eventSel.variable << " < " << eventSel.maxval1 << " OR ";
//                    cout << eventSel.minval2 << " < " << eventSel.variable << " < " << eventSel.maxval2 << endl;
//                    cout << " " << selected << " " << dataSetTag[iData]+"_"+eventSel.variable << " " << selectionVariable << endl;
                }

                weight=array_variable_map_all[dataSetTag[iData]+"_"+args.weightVar][ientry];
                Mpi0=array_variable_map_all[dataSetTag[iData]+"_Mpi0"][ientry];
                Meta=array_variable_map_all[dataSetTag[iData]+"_Meta"][ientry];
                accidental=array_variable_map_all[dataSetTag[iData]+"_AccWeight"][ientry];
                //cout << "weight, Mpi0, Meta, accidental " << weight << ", " << Mpi0 << ", " << Meta << ", " << accidental << endl;
                nentries_preDefaultSelects+=accidental; // 03/18/23 changed to accidental. For syst comparison we cant be using weightASBS

                if (!selected)
                    continue;

                //////////////////////
                // MAKE SELECTIONS
                //////////////////////
                // ******* Select Pi0 Signal // There really isn't any background. Just set to true
                //selectPi0=(Mpi0>pi0sigL)*(Mpi0<pi0sigR);
                selectPi0=true; 

                // ******* Select Eta Sidebands
                //selectEta=(Meta>etasbLL)*(Meta<etasbLR)||(Meta>etasbRL)*(Meta<etasbRR);
                // ******* Select Eta Right Sideband
                //selectEta=(Meta>etasbRL)*(Meta<etasbRR);
                // ******* Select Eta Left Sideband
                //selectEta=(Meta>etasbLL)*(Meta<etasbLR);
                // ******* Select Eta Signal
                //selectEta=(Meta>etasigL)*(Meta<etasigR);

                // ******* Select Eta Right Half Peak
                //selectEta=(Meta>eta_peak);
                // ******* Select Eta Left Half Peak
                //selectEta=(Meta<eta_peak);
                
                // Do not select on Mpi0 nor Meta if we are going to do sideband subtraction also
                selectEta=true;
                if(selectPi0*selectEta){ 
                    for(string variable: variables){
                        array_variable_map[dataSetTag[iData]+"_"+variable].push_back(array_variable_map_all[dataSetTag[iData]+"_"+variable][ientry]); 
                    }
                    array_BeamAngles_map[dataSetTag[iData]].push_back(array_BeamAngles_map_all[dataSetTag[iData]][ientry]);

                    dHist_Mpi0->Fill(Mpi0,accidental);
                    dHist_Meta->Fill(Meta,accidental);
                    // Fill the diagnostic histograms with accidental weights
                    for (auto const& ele: diagnosticHists){
                        ele.second->Fill(array_variable_map_all[dataSetTag[iData]+"_"+ele.first][ientry],accidental);
                    }
                    signalIntegral+=accidental; // 03/18/23 changed to accidental since the histogram is actually using accidental weights
                    nentries_postDefaultSelects+=accidental; // 03/18/23 changed to accidental. For syst comparison we cant be using weightASBS
                }
            }
            cout << "  remaining weighted entries percentage after default selections and possible sidebands subtraction: " 
                << std::setprecision(3) << nentries_postDefaultSelects/nentries_preDefaultSelects*100 << "%" << endl;
            total_nentries_preDefaultSelects+=nentries_preDefaultSelects;
            total_nentries_postDefaultSelects+=nentries_postDefaultSelects;
        }
        cout << "TOTAL weighted nentries: " << total_nentries_preDefaultSelects << endl;
        float percentage = (float)total_nentries_postDefaultSelects/total_nentries_preDefaultSelects*100;
        percentages[evtTag+sbTag]=percentage;
        cout << "  remaining weighted nentries percentage after default selections: " << std::setprecision(3) << percentage << "%" << endl;

        gSystem->Exec(("mkdir -p "+outputFolder+"/diagnostics").c_str());

        // Draw a diagnostic plot to denote the Meta and Mpi0 distributions
        dHist_Mpi0->SetMinimum(0);
        dHist_Meta->SetMinimum(0);
        c1->Divide(2,1);
        c1->cd(1);
        dHist_Mpi0->GetXaxis()->SetNdivisions(5);
        dHist_Mpi0->Draw("HIST");
        dHist_Mpi0->SetTitle(("Signal Integral: "+to_string((int)signalIntegral)).c_str());
        //box->SetFillColorAlpha(kGreen+2,0.3);
        //box->DrawBox(pi0sigL,0,pi0sigR,dHist_Mpi0->GetMaximum()*1.0);
        //box->SetFillColorAlpha(kRed+1,0.3);
        //box->DrawBox(pi0sbLL,0,pi0sbLR,dHist_Mpi0->GetMaximum()*1.0);
        //box->DrawBox(pi0sbRL,0,pi0sbRR,dHist_Mpi0->GetMaximum()*1.0);
        c1->cd(2);
        dHist_Meta->Draw("HIST");
        dHist_Meta->SetTitle(("Signal Integral: "+to_string((int)signalIntegral)).c_str());
        box->SetFillColorAlpha(kGreen+2,0.3);
        box->DrawBox(etasigL,0,etasigR,dHist_Meta->GetMaximum()*1.0);
        box->SetFillColorAlpha(kRed+1,0.3);
        box->DrawBox(etasbLL,0,etasbLR,dHist_Meta->GetMaximum()*1.0);
        box->DrawBox(etasbRL,0,etasbRR,dHist_Meta->GetMaximum()*1.0);
        c1->SaveAs((outputFolder+"mass_plots"+evtTag+sbTag+".pdf").c_str());

        for (auto const& ele: diagnosticHists){
            c1->Clear();
            c1->Divide(1,1);
            if (ele.first=="unusedEnergy")
                gPad->SetLogy(1);
            else
                gPad->SetLogy(0);
            ele.second->Draw("HIST");
            c1->SaveAs((outputFolder+"diagnostics/"+ele.first+evtTag+sbTag+".png").c_str());
        }

        cout << "Opening a csv to write results to..." << endl;
        ofstream saveCsv;
        saveCsv.open((outputFolder+"results"+evtTag+sbTag+".csv").c_str());
        saveCsv << "binVar data midt1 t1_err asym_000_eta asym_000_eta_err asym_045_eta asym_045_eta_err ";
        saveCsv << "asym_000_pi0 asym_000_pi0_err asym_045_pi0 asym_045_pi0_err ";
        saveCsv << "000090_eta_chiPerDOF 045135_eta_chiPerDOF 000090_pi0_chiPerDOF 045135_pi0_chiPerDOF ";
        saveCsv << "entries_000_eta entries_045_eta entries_090_eta entries_135_eta ";
        saveCsv << "entries_000_pi0 entries_045_pi0 entries_090_pi0 entries_135_pi0";
        for (int i=0; i<(int)(sizeof(polnames)/sizeof(polnames[0])); ++i){
            saveCsv << " " << polnames[i] << "_eta " << polnames[i] << "_eta_err " << polnames[i] << "_pi0 " << polnames[i] << "_pi0_err";
        }
        for (int i=1; i<3; ++i){ // the selection is the sub-binning i.e. u3,s12,s23. Second is for selecting Mpi0p region
            saveCsv << " selection" << i << " minval" << i << " maxval" << i;
        }
        saveCsv << endl;

        string fitFolder;

        ////////////////////////////////////////////////////////////
        //   BINNED IN U3
        ////////////////////////////////////////////////////////////
        fitFolder=outputFolder+"binned_u3"+evtTag+sbTag+"/";
        gSystem->Exec(("mkdir -p "+fitFolder).c_str());
        
        //const int ntbins=1;
        //float mints[ntbins]={0.000};
        //float maxts[ntbins]={100.0};
        const int ntbins=4;
        float mints[ntbins]={0.0,0.5,1.000,0.000};
        float maxts[ntbins]={0.5,1.0,100.0,100.0};
        for(int it=0; it<ntbins; ++it){
            vector<criteria> selections={
                {"mandelstam_t",mints[it],maxts[it],fltmin,fltmin},
                {"Mpi0p",1.4,100,fltmin,fltmin} 
            };
            string fitSaveLoc=fitFolder+"asym_tB"+to_string(it); 
            constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,args,fitSaveLoc,&saveCsv);
        }

        ////////////////////////////////////////////////////////////
        //   BINNED IN PI0P
        ////////////////////////////////////////////////////////////
        fitFolder=outputFolder+"binned_spi0p"+evtTag+sbTag+"/";
        gSystem->Exec(("mkdir -p "+fitFolder).c_str());
        //float minPi0P[5] = {1.15, 1.4, 1.6, 1.75, 1.95};
        //float maxPi0P[5] = {1.4, 1.6, 1.75, 1.95, 2.4};
        float minPi0P[3] = {1.15, 1.4, 2.20};
        float maxPi0P[3] = {1.40, 2.2, 10.0};
        for(int it=0; it<3; ++it){
            vector<criteria> selections={
                {"Mpi0p",minPi0P[it],maxPi0P[it],fltmin,fltmin},
                {"Mpi0p",0,100,fltmin,fltmin} 
            };
            string fitSaveLoc=fitFolder+"asym_Mpi0pB"+to_string(it); 
            constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,args,fitSaveLoc,&saveCsv);
        }

        ////////////////////////////////////////////////////////////
        //   BINNED IN METAP
        ////////////////////////////////////////////////////////////
        fitFolder=outputFolder+"binned_setap"+evtTag+sbTag+"/";
        gSystem->Exec(("mkdir -p "+fitFolder).c_str());
        //float minEtaP[5] = {1.5, 1.65, 1.9, 2.2, 2.4};
        //float maxEtaP[5] = {1.65, 1.9, 2.2, 2.4, 2.8};
        float minEtaP[3] = {0.0, 2.1, 2.60};
        float maxEtaP[3] = {2.1, 2.6, 10.0};
        for(int it=0; it<3; ++it){
            vector<criteria> selections={
                {"Metap",minEtaP[it],maxEtaP[it],fltmin,fltmin},
                {"Mpi0p",1.4,100,fltmin,fltmin} 
            };
            string fitSaveLoc=fitFolder+"asym_Mpi0pB"+to_string(it); 
            constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,args,fitSaveLoc,&saveCsv);
        }

        ////////////////////////////////////////////////////////////
        //   BINNED IN MPI0ETA
        ////////////////////////////////////////////////////////////
        fitFolder=outputFolder+"binned_s12"+evtTag+sbTag+"/";
        cout << "Making fit folder: " << fitFolder << endl;
        gSystem->Exec(("mkdir -p "+fitFolder).c_str());
        //float minMpi0eta[5] = {1.65, 1.9, 2.15, 2.4, 2.65};
        //float maxMpi0eta[5] = {1.9, 2.15, 2.4, 2.65, 2.9};
        //float minMpi0eta[7] = {1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8};
        //float maxMpi0eta[7] = {1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
        float minMpi0eta[3] = {1.6, 2.1, 2.6};
        float maxMpi0eta[3] = {2.1, 2.6, 3.1};
        for(int it=0; it<(int)(sizeof(minMpi0eta)/sizeof(minMpi0eta[0])); ++it){
            vector<criteria> selections={
                {"Mpi0eta",minMpi0eta[it],maxMpi0eta[it],fltmin,fltmin},
                {"Mpi0p",1.4,100,fltmin,fltmin} 
            };
            string fitSaveLoc=fitFolder+"asym_Mpi0etaB"+to_string(it); 
            constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,args,fitSaveLoc,&saveCsv);
        }
        saveCsv.close();
    }
    return percentages;
}














#include "DSelector_ver20.h"
bool showOutput = false;
bool showMassCalc = false;
bool onlyNamesPi0_1 = true; // true if we want to show only the histograms with _1 in their names so we can merge them with _2

int itersToRun = 0;
int protonID=0;

// Set topologyString to empty string to not select on the topology string. Oherwise we will compare the toplogy string to and return if not equal
//string topologyString="4#gammap";
string topologyString="4#gammap[2#pi^{0}]"; // for some reason the topology string in the thrown trees are 4#gammap instead of 4#gammap[2#pi^{0}]

string selectDetector="ALL";
string polarization="degALL";
string tag="_data_2017_mEllipse_8288_tLT1";

bool is_pi0eta=false;
int mcprocess=0;

void DSelector_ver20::Init(TTree *locTree)
{
        count_combos=0;
        cout << "STARTING" << endl;
	targetCenter = {0,0,65};

        if (is_pi0eta) {
                lowMass = 0.7;
                upMass = 2.0;
                etaProtonBaryonCut = 1.65;
                pi0ProtonBaryonCut = 2;
		binScale = (upMass-lowMass)/numBinsMass;
                //using the kin data
                ellipseX = 0.135784; ellipseY = 0.548036; ellipseXr = 2*0.00753584; ellipseYr = 2*0.0170809;
	        ellipseXr_loose = 5*0.00753584; ellipseYr_loose=5*0.0170809;
        }
        else {
                lowMass = 0.3;
                upMass = 1.6;
                pi0ProtonBaryonCut = 1.8;
                etaProtonBaryonCut = 1.8;
                //ellipseX = 0.134547; ellipseY = 0.134547; ellipseXr = 0.025449; ellipseYr = 0.025449;
                //using the kin data
                ellipseX = 0.135881; ellipseY = 0.135881; ellipseXr = 0.0160375; ellipseYr = 0.0160375;
                ellipseXBS1 = 0.135881; ellipseYBS1 = 0.135881; ellipseXrBS1 = 0.022; ellipseYrBS1 = 0.022;
                ellipseXBS2 = 0.135881; ellipseYBS2 = 0.135881; ellipseXrBS2 = 0.045; ellipseYrBS2 = 0.045;
		ellipseXr_loose=0.0391; ellipseYr_loose=0.0391;
		areaRatio = 0.167; //double checked these values
        }

        //USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
        dOutputFileName = polarization+tag+"_DSelector_output.root"; //"" for none
        dOutputTreeFileName = polarization+tag+"_tree_DSelector.root"; //"" for none
        dFlatTreeFileName = polarization+tag+"_treeFlat_DSelector.root"; //output flat tree (one combo per tree entry), "" for none
        dFlatTreeName = "tree_4g_flat"; //if blank, default name will be chosen
	dSaveDefaultFlatBranches = false; // False: don't save default branches, reduce disk footprint.

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	// dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	// dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/		

	// cutString will be used by alot of histograms when defining names 
	// Just defining alot of stuff before hand which will be used in the labeling of the histograms 
	std::string cutString;
	std::string cutsApplied="";
	std::string cutsBase="";
	std::vector<std::string> cutVariations;
	std::string massBin;   

        cout << "INITILIZED ACTIONS" << endl;
        countThrownEvents = new TH1F("count thrown events", "Cuts=noCut;counts", 1, 0, 1);
        //dHist_thrown_tp = new TH1F("Thrown tp", "Cuts=noCut;-t' momentum transfer of #pi^{0}+#eta;Entries / 0.06 GeV", 100, 0, 6);
        dHist_thrown_tp_selected = new TH1F("Thrown tp selected", "Cuts=noCut;-t' momentum transfer of #pi^{0}+#eta;Entries / 0.06 GeV", 100, 0, 6);

        dHist_BeamAngle = new TH1F("BeamAngle", "Beam Angle with no cuts applied;Beam Angle (GeV)", 180,-180,180);
        dHist_BeamAngle->SetYTitle("Entries / 2 Degree");
	// MPE = multiple photons per event. So instead of just filling the beam angle once per run we will fill it per combo.
        dHist_BeamAngleMPE = new TH1F("BeamAngleMPE", "Beam Angle with no cuts applied;Beam Angle (GeV)", 180,-180,180);
        dHist_BeamAngleMPE->SetYTitle("Entries / 2 Degree");
	dHist_Cuts = new TH1F("CutsPassed", "Number of times a cut has been passed", 17,0,17);
	//dHist_numCombos = new TH1I("numCombos", "" , 5,0,5);
	for (int i =0; i<3; ++i){
		if (is_pi0eta){
			dHist_checkEllipseBS[i] = new TH2F(("checkEllipseBS"+std::to_string(i)+"noCutOnlyRegionSelected").c_str(), ";#pi^{0} Mass (GeV) with Entries / 0.001 GeV;#eta Mass (GeV) with Entries / 0.0025 GeV", atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()), atof(etaBinRange[0].c_str()), atof(etaBinRange[1].c_str()), atof(etaBinRange[2].c_str()));
		}
		else {
			dHist_checkEllipseBS[i] = new TH2F(("checkEllipseBS"+std::to_string(i)+"noCutOnlyRegionSelected").c_str(), ";#pi^{0} Mass (GeV) with Entries / 0.001 GeV;#pi^{0} Mass (GeV) with Entries / 0.001 GeV", atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()), atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()));
		}
	}

        /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

        //EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
        //The type for the branch must be included in the brackets
        //1st function argument is the name of the branch
        //2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
        dTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta_thrown");
        dTreeInterface->Create_Branch_Fundamental<Int_t>("mcprocess");
        dTreeInterface->Create_Branch_FundamentalArray<Float_t>("DataWeight","NumCombos");
        /*
           dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
           dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
           dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
           dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
           */

        /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

        //EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
        //The type for the branch must be included in the brackets
        //1st function argument is the name of the branch
        //2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
        /*
           dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
           dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
           dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
           dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
           */

        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("mcprocess");
	// *** SOME THROWN VARIABLE CHECK
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tp_thrown"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_t_thrown"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Ebeam_thrown"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Ebeam"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("whichSignalRegion"); // *** 1=signal 2=sideband 3=skip
        dFlatTreeInterface->Create_Branch_Fundamental<bool>("isTruePi0Eta"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<bool>("insideEllipse"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta_thrown");
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("numSpect");
        dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("spectroscopicID","numSpect"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<bool>("beamPhotonMatchToThrown"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("rfTime"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<UChar_t>("locNumUnusedShowers");
	// **************************

        // ** Including some basic cuts  and associated variables
        // **************************
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mmsq");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_momentum");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_z");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_R");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_dEdxCDC");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_13");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_24");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_23");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_14");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("omegaCut");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pVH");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pPhotonE");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pPhotonTheta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pMagP3Proton");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pzCutmin");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pRProton");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pMissingMassSquared");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pdEdxCDCProton");

        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("DOFKinFit"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("AccWeight"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("BeamAngle"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightASBS"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightBS"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightBSpi0"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightBSeta"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("uniqueComboID"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tp"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_t"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("ptGT1");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("ptLT05");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("ptGT05LT1");
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("protonID"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("beamID"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("chiSq"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("chiSq_gpi0pi0"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("unusedEnergy"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0p");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Metap");
        if (is_pi0eta){
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0g3"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0g4"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("logMpi0g3"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("logMpi0g4"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0g3");
        	dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0g4");

        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Meta"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0_meas"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Meta_meas"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta_meas"); //fundamental = char, int, float, double, etc.
		dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_teta_meas");	
		dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_teta");	
		dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tpi0_meas");	
		dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tpi0");	
        }
        else{
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0pi0"); //fundamental = char, int, float, double, etc.
        }
        // introducing some variables to keep track of how much unique particles we have
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_eta_pi0eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0_pi0eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("baseAsymCut");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("baseAsymCut_mDelta");

	// introducing variables to show the pi0/eta was detected in
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("pi0DetectedIn"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("etaDetectedIn"); //fundamental = char, int, float, double, etc.
	
        // Introduce some angles to use in the phase space distance calculation
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("prodPlane_x"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("prodPlane_y"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("prodPlane_z"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_X_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_pi0_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_X_cm_meas"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_X_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_X_lab"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Phi"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_lab"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_pi0_lab"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_gj_meas"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_gj_meas"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_hel"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_hel"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosThetaHighestEphotonIneta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosThetaHighestEphotonInpi0_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_pi0_g3"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_pi0_g4"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_g12"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_g13"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_g14"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_g23"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_g24"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph123Rest_angle_g34"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_g12"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_g13"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_g14"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_g23"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_g24"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("ph124Rest_angle_g34"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_x");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_y");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_omega");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("pi0_energy");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta1");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta2");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta3");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta4");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE1");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE2");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE3");	
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE4");	
        //dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phiHighestEphotonIneta_etaFrame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosThetaHighestEphotonInpi0_pi0Frame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phiHighestEphotonInpi0_pi0Frame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("HighestEphotonIneta_etaFrame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("HighestEphotonIneta_pi0Frame"); //fundamental = char, int, float, double, etc.

        /************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

        //TO SAVE PROCESSING TIME
        //If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
        //By default, for each event, the data is retrieved for all branches
        //If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
        //Do this by doing something similar to the commented code below

        //dTreeInterface->Clear_GetEntryBranches(); //now get none
        //dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
	//
	if(showOutput) { cout << "Finished creating branch funamentals" <<endl; }

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);
} // end of initialization

Bool_t DSelector_ver20::Process(Long64_t locEntry)
{
    	if(showOutput){cout << "Starting next process looping" << endl;}
    	// The Process() function is called for each entry in the tree. The entry argument
    	// specifies which entry in the currently loaded tree is to be processed.
    	//
    	// This function should contain the "body" of the analysis. It can contain
    	// simple or elaborate selection criteria, run algorithms on the data
    	// of the event and typically fill histograms.
    	//
    	// The processing can be stopped by calling Abort().
    	// Use fStatus to set the return value of TTree::Process().
    	// The return value is currently not used.

    	//CALL THIS FIRST
    	DSelector::Process(locEntry); //Gets the data from the tree for the entry
    	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
    	//TLorentzVector locProductionX4 = Get_X4_Production();

    	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

    	//Only if the run number change
    	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
    	UInt_t locRunNumber = Get_RunNumber();
    	// we must have the following condition or else we run into errors, or maybe rcdb hangs due to the amount of querries or something. Anyways, it wouldn't work. 
    	if(locRunNumber != dPreviousRunNumber)
    	{
    	    dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
            // AMO RUNS HAVE VALUE OF -1 , CHECK GLUEX_ROOT_ANALYSIS
    	    hasPolarizationAngle = dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, locPolarizationAngle);
    	    if(showOutput){cout << "Getting beam polarization and filling used runs" << endl;}
    	    if (usedRuns.find(locRunNumber)==usedRuns.end()){
    	        if (hasPolarizationAngle) {
    	            dHist_BeamAngle->Fill(locPolarizationAngle);
    	        }
    	        else {
    	    	dHist_BeamAngle->Fill(-1);
    	        }

    	        usedRuns.insert(locRunNumber);
    	    }
    	    // Making a cut on the polarization angle to split the dataset up
    	    dPreviousRunNumber = locRunNumber;
    	}

        // ARE WE ONLY ASKING FOR A SPECIFIC POLARIZATION WITH THE "polarization" VARIABLE
    	keepPolarization=false;
    	if (polarization=="degALL"){ keepPolarization=true; }
    	if (locPolarizationAngle==0 && polarization=="deg000") { keepPolarization=true; }
    	if (locPolarizationAngle==45 && polarization=="deg045") { keepPolarization=true; }
    	if (locPolarizationAngle==90 && polarization=="deg090") { keepPolarization=true; }
    	if (locPolarizationAngle==135 && polarization=="deg135") { keepPolarization=true; }
    	if (!hasPolarizationAngle && polarization=="degAMO") { keepPolarization=true; }
	if (!keepPolarization) { cout << "Throwing out the event since locPolarizationAngle is " << locPolarizationAngle << " and our selection criteria requires " << polarization << endl; }
    

        // WHAT IS THE BEAM POLARIZATION IN THIS EVENT
    	keepPolarization000=false;
    	keepPolarization045=false;
    	keepPolarization090=false;
    	keepPolarization135=false;
    	keepPolarizationAMO=false;
    	if ( locPolarizationAngle==0 ) { keepPolarization000=true; }
    	if ( locPolarizationAngle==45 ) { keepPolarization045=true; }
    	if ( locPolarizationAngle==90 ) { keepPolarization090=true; }
    	if ( locPolarizationAngle==135 ) { keepPolarization135=true; }
    	if ( !hasPolarizationAngle ) { keepPolarizationAMO=true; }


	std::vector<int> parentArray;
	std::vector<int> pids;
	int locNumThrown = Get_NumThrown();

	/************************************************* PARSE THROWN TOPOLOGY ***************************************/
	TString locThrownTopology = Get_ThrownTopologyString();
        //if ((topologyString!="")*dIsMC){
        //    if (locThrownTopology != topologyString.c_str() ){
        //        cout << "incorrect locThrownTopology = " << locThrownTopology << endl;
        //        return kTRUE;
        //    }
        //    else {
        //        cout << "correct locThrownTopology = " << locThrownTopology << endl;
        //    }
        //}
	//countThrownEvents->Fill(1);

	//TLorentzVector locProtonP4_thrown;
        //Particle_t locPID;
        //Int_t locParentID;
	//for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	//{	
	//    dThrownWrapper->Set_ArrayIndex(loc_i);
	//    locPID = dThrownWrapper->Get_PID();
	//    locParentID = dThrownWrapper->Get_ParentIndex();
        //    if ( (locPID==14)*(locParentID==-1) ){ // if its a proton and it has no parent, is the primary recoiled proton
        //        mandelstam_t_thrown = -(dThrownWrapper->Get_P4()-dTargetP4).M2();
        //    }
        //}

        float locBeamE_thrown=0;
        float locT_thrown=0; 
        TLorentzVector locProtonP4_thrown;

	if(dThrownBeam != NULL)
		locBeamE_thrown = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
	    //Set branch array indices corresponding to this particle
	    dThrownWrapper->Set_ArrayIndex(loc_i);
	    Particle_t locPID = dThrownWrapper->Get_PID();
	    TLorentzVector locThrownP4_thrown = dThrownWrapper->Get_P4();
            if (locPID==14)
                locProtonP4_thrown=locThrownP4_thrown;
	}

	locMetapi0_thrown=(locPi0P4_thrown+locEtaP4_thrown).M();
	locT_thrown=-(dTargetP4-locProtonP4_thrown).M2();		
	//bool bMetapi0_thrown = (locMetapi0_thrown>1.04)*(locMetapi0_thrown<1.56);
	bool bmandelstamt_thrown=(locT_thrown<1.0)*(locT_thrown>0.1); 
        bool bBeamE_thrown = (locBeamE_thrown<8.8)*(locBeamE_thrown>8.2);
        bool bTopology = locThrownTopology==topologyString; 
        bool selection_thrown=bTopology*bBeamE_thrown*bmandelstamt_thrown//*bMpi0eta_thrown;
        if (dIsMC*!selection_thrown)
            return kTRUE;

    /********************************************* SETUP UNIQUENESS TRACKING ********************************************/
    // For uniqueness tracking we have to make a separate uniqueness tracking for each combination. This is partly because our general cuts span numerous quantities that can come from individual particles
    // from particle pairs, and from the entire combo. This makes it difficult to not apply the pi0 cut when filling the M(pi0) hist since two combos can contain the same 2 photons yet one can pass and the
    // other not (I think this could happen). If it can happen then the order in which the algorithm encouters the combo can fill the plot one or twice due to the other particles in the combo affecting the cut.
    // So we have to make filling the histogram and inserting into the uniquess set based on the same condition. 
    // There are 4Gamma + 1Proton + 1Beam in a combo. 

    //ANALYSIS ACTIONS: Reset uniqueness tracking for each action
    //For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
    Reset_Actions_NewEvent();
    //dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

    //PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
    //Sometimes, some content is the exact same between one combo and the next
    //e.g. maybe two combos have different beam particles, but the same data for the final-state
    //When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
    //So, for each quantity you histogram, keep track of what particles you used (for a givMpi0g4//Then for each combo, just compare to what you used before, and make sure it's unique


    //EXmap<Particle_t, AMPLE 2: Combo-specific info:
    //In general: Could have multiple particles with the same PID: Use a set of Int_t's
    //In general: Multiple PIDs, so multiple sets: Contain within a map
    //Multiple combos: Contain maps within a set (easier, faster to search)
    //use set<map .... > > even for reactions with the only one type of final state particles as in the case of 4 gammas
    // MAKING SURE ALL THE SETS ARE EMPTY SO WE DONT START A LOOP WITH ELEMENTS IN IT ALREADY...

    /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

    /*
       Int_t locMyInt = 7;
       dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

       TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
       dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

       for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
       dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
       */
    dTreeInterface->Fill_Fundamental<Int_t>("mcprocess", mcprocess);

    /************************************************* LOOP OVER COMBOS *************************************************/



    map<TString, vector<topology>> topologyMap;

    //Loop over combos
    //dHist_numCombos->Fill(Get_NumCombos());
    for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
    {
        if(showOutput) { cout << "\n\n\n*********************************************************\n**************************************************\n########    EventIdx, ComboIdx: " << (eventIdx) << ", " << loc_i << "    #############" << endl; }

	// Fill the polarization angle MPE = multiple photons per event or for all combinations. We will do this instead of filling per event just in case there was for whatever reason a photon in a combo can be matched with a tagged photon with a different beam polarization
	if (keepPolarization) {
    		dHist_BeamAngleMPE->Fill(locPolarizationAngle);
    	}

        //Set branch array indices for combo and all combo particles
        dComboWrapper->Set_ComboIndex(loc_i);

        // Is used to indicate when combos have been cut
        if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
            continue; // Combo has been cut previously

        /********************************************** GET PARTICLE INDICES *********************************************/

        //Used for tracking uniqueness when filling histograms, and for determining unused particles

        //Step 0
        if(showOutput){cout << "** Getting tracks p4, x4, ids" << endl;}
        Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
        if(showOutput){cout <<"Got beamID"<<endl;}
        Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
        if(showOutput){cout <<"Got ProtonID"<<endl;}

        //Step 1
        Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
        Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();
        if(showOutput){cout <<"Got photon1/2 ID"<<endl;}

        //Step 2
        Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();
        Int_t locPhoton4NeutralID = dPhoton4Wrapper->Get_NeutralID();
        if(showOutput){cout <<"Got photon3/4 ID"<<endl;}
        if(showOutput){cout << "Got Beam, proton, 4 photon IDs" << endl;}

        // ********************************************************************************************************************************
        // Get 4 Momenta!!! 
        // ********************************************************************************************************************************

        // Get P4's: //is kinfit if kinfit performed, else is measured
        //dTargetP4 is target p4
        if(showOutput){cout << "** Getting kin fitted X4, P4" << endl;}
        //Step 0
        TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4_Measured();
        TLorentzVector locBeamX4 = dComboBeamWrapper->Get_X4_Measured();


        //if(showOutput){cout << "Got kin fitted beam x4 P4" << endl;}
        TLorentzVector locProtonP4 = dProtonWrapper->Get_P4_Measured();
        TLorentzVector locProtonX4 = dProtonWrapper->Get_X4_Measured();
        protonX4[0]=locProtonX4.X();
        protonX4[1]=locProtonX4.Y();
        protonX4[2]=locProtonX4.Z();
        protonX4[3]=locProtonX4.T();

        //Step 1
        // using dDecayingPi0/EtaWrapper->Get_P4 no longer works, gives me an error now, luckily we don't use it anyways
        //TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
        TLorentzVector locPhoton1X4_Kin = dPhoton1Wrapper->Get_X4();
        TLorentzVector locPhoton2X4_Kin = dPhoton2Wrapper->Get_X4();
        if(showOutput){cout << "Got kin fitted photon 1,2 x4 and P4" << endl;}
        //Step 2
        //TLorentzVector locDecayingEtaP4 = dDecayingEtaWrapper->Get_P4();
        TLorentzVector locPhoton3X4_Kin = dPhoton3Wrapper->Get_X4();
        TLorentzVector locPhoton4X4_Kin = dPhoton4Wrapper->Get_X4();
        if(showOutput){cout << "Got kin fitted photon 3,4 x4 and P4" << endl;}


        // Get Measured P4's:
        //Step 0
        if(showOutput){cout << "** Getting measured x4 p4" << endl;}
        TLorentzVector locBeamP4_Kin = dComboBeamWrapper->Get_P4();
        TLorentzVector locBeamX4_Kin = dComboBeamWrapper->Get_X4();
        TLorentzVector locProtonP4_Kin = dProtonWrapper->Get_P4();
        TLorentzVector locProtonX4_Kin = dProtonWrapper->Get_X4();
        //Step 1
        TLorentzVector locPhoton1P4_Kin = dPhoton1Wrapper->Get_P4();
        TLorentzVector locPhoton2P4_Kin = dPhoton2Wrapper->Get_P4();
        TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4_Measured();
        TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4_Measured();
        //Step 2
        TLorentzVector locPhoton3P4_Kin = dPhoton3Wrapper->Get_P4();
        TLorentzVector locPhoton4P4_Kin = dPhoton4Wrapper->Get_P4();
        TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4_Measured();
        TLorentzVector locPhoton4P4 = dPhoton4Wrapper->Get_P4_Measured();

        // We need to use shower position since the photons do not really have a track. 
        TLorentzVector locPhoton1X4_Shower = dPhoton1Wrapper->Get_X4_Shower();
        TLorentzVector locPhoton2X4_Shower = dPhoton2Wrapper->Get_X4_Shower();
        TLorentzVector locPhoton3X4_Shower = dPhoton3Wrapper->Get_X4_Shower();
        TLorentzVector locPhoton4X4_Shower = dPhoton4Wrapper->Get_X4_Shower();
        TLorentzVector locPhoton1X4 = dPhoton1Wrapper->Get_X4_Measured();
        TLorentzVector locPhoton2X4 = dPhoton2Wrapper->Get_X4_Measured();
        TLorentzVector locPhoton3X4 = dPhoton3Wrapper->Get_X4_Measured();
        TLorentzVector locPhoton4X4 = dPhoton4Wrapper->Get_X4_Measured();


        std::vector< TLorentzVector > allPhotonP4Vectors_Kin;
        allPhotonP4Vectors_Kin.push_back(locPhoton1P4_Kin);
        allPhotonP4Vectors_Kin.push_back(locPhoton2P4_Kin);
        allPhotonP4Vectors_Kin.push_back(locPhoton3P4_Kin);
        allPhotonP4Vectors_Kin.push_back(locPhoton4P4_Kin);
        std::vector< TLorentzVector > allPhotonX4Vectors_Kin;
        allPhotonX4Vectors_Kin.push_back(locPhoton1X4_Kin);
        allPhotonX4Vectors_Kin.push_back(locPhoton2X4_Kin);
        allPhotonX4Vectors_Kin.push_back(locPhoton3X4_Kin);
        allPhotonX4Vectors_Kin.push_back(locPhoton4X4_Kin);
        std::vector< TLorentzVector > allPhotonX4Vectors;
        allPhotonX4Vectors.push_back(locPhoton1X4);
        allPhotonX4Vectors.push_back(locPhoton2X4);
        allPhotonX4Vectors.push_back(locPhoton3X4);
        allPhotonX4Vectors.push_back(locPhoton4X4);
        std::vector< TLorentzVector > allPhoton4Vectors_Shower;
        allPhoton4Vectors_Shower.push_back(locPhoton1X4_Shower);
        allPhoton4Vectors_Shower.push_back(locPhoton2X4_Shower);
        allPhoton4Vectors_Shower.push_back(locPhoton3X4_Shower);
        allPhoton4Vectors_Shower.push_back(locPhoton4X4_Shower);
        std::vector< TLorentzVector > allPhoton4Vectors_meas;
        allPhoton4Vectors_meas.push_back(locPhoton1P4);
        allPhoton4Vectors_meas.push_back(locPhoton2P4);
        allPhoton4Vectors_meas.push_back(locPhoton3P4);
        allPhoton4Vectors_meas.push_back(locPhoton4P4);
        std::vector< DNeutralParticleHypothesis* > allPhotonWrappers;
        allPhotonWrappers.push_back(dPhoton1Wrapper);
        allPhotonWrappers.push_back(dPhoton2Wrapper);
        allPhotonWrappers.push_back(dPhoton3Wrapper);
        allPhotonWrappers.push_back(dPhoton4Wrapper);
        std::vector<Int_t> photonIds;
        photonIds.push_back(locPhoton1NeutralID);
        photonIds.push_back(locPhoton2NeutralID);
        photonIds.push_back(locPhoton3NeutralID);
        photonIds.push_back(locPhoton4NeutralID);

	int nPhotons = (int) (allPhotonP4Vectors_Kin.size());
	if(showOutput){cout << "There are " << nPhotons << " photons" << endl;  }

        if(showOutput){cout << "Got measured x4 p4 with using shower X4 for photons" << endl;}


        /********************************************* COMBINE FOUR-MOMENTUM ********************************************/

        // DO YOUR STUFF HERE

        if(showOutput){cout << "** Combining 4-vectors to get pi0, eta, pi0eta, etaProton, pi0Proton" << endl;}
        // Combine 4-vectors
        TLorentzVector locMissingP4 = locBeamP4 + dTargetP4;
        locMissingP4 -= locProtonP4 + locPhoton1P4 + locPhoton2P4 + locPhoton3P4 + locPhoton4P4;

        //TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
        //TLorentzVector locEtaP4 = locPhoton3P4 + locPhoton4P4;

        TLorentzVector locPi0P4_Kin = locPhoton1P4_Kin + locPhoton2P4_Kin;
        TLorentzVector locEtaP4_Kin = locPhoton3P4_Kin + locPhoton4P4_Kin;
        TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
        TLorentzVector locEtaP4 = locPhoton3P4 + locPhoton4P4;

        TLorentzVector mismatchPair_13 = locPhoton1P4_Kin + locPhoton3P4_Kin;
        TLorentzVector mismatchPair_14 = locPhoton1P4_Kin + locPhoton4P4_Kin;
        TLorentzVector mismatchPair_23 = locPhoton2P4_Kin + locPhoton3P4_Kin;
        TLorentzVector mismatchPair_24 = locPhoton2P4_Kin + locPhoton4P4_Kin;
        TLorentzVector mismatchPair_132 = locPhoton1P4_Kin + locPhoton3P4_Kin + locPhoton2P4_Kin;
        TLorentzVector mismatchPair_134 = locPhoton1P4_Kin + locPhoton3P4_Kin + locPhoton4P4_Kin;
        TLorentzVector mismatchPair_142 = locPhoton1P4_Kin + locPhoton4P4_Kin + locPhoton2P4_Kin;
        TLorentzVector mismatchPair_143 = locPhoton1P4_Kin + locPhoton4P4_Kin + locPhoton3P4_Kin;
        TLorentzVector mismatchPair_231 = locPhoton2P4_Kin + locPhoton3P4_Kin + locPhoton1P4_Kin;
        TLorentzVector mismatchPair_234 = locPhoton2P4_Kin + locPhoton3P4_Kin + locPhoton4P4_Kin;
        TLorentzVector mismatchPair_241 = locPhoton2P4_Kin + locPhoton4P4_Kin + locPhoton1P4_Kin;
        TLorentzVector mismatchPair_243 = locPhoton2P4_Kin + locPhoton4P4_Kin + locPhoton3P4_Kin;
        TLorentzVector mismatchPair_341 = locPhoton3P4_Kin + locPhoton4P4_Kin + locPhoton1P4_Kin;
        TLorentzVector mismatchPair_342 = locPhoton3P4_Kin + locPhoton4P4_Kin + locPhoton2P4_Kin;
        TLorentzVector mismatchPair_123 = locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton3P4_Kin;
        TLorentzVector mismatchPair_124 = locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton4P4_Kin;

        TLorentzVector mixingPi0Eta_Kin = locPi0P4_Kin + locEtaP4_Kin;
        TLorentzVector mixingPi0Eta = locPi0P4 + locEtaP4;
        TLorentzVector mixingEtaProton_Kin = locEtaP4_Kin+locProtonP4_Kin;
        TLorentzVector mixingPi0Proton_Kin = locPi0P4_Kin+locProtonP4_Kin;

	massGammaPi0[0] = (locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton3P4_Kin).M();
	massGammaPi0[1] = (locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton4P4_Kin).M();
	massGammaEta[0] = (locPhoton3P4_Kin + locPhoton4P4_Kin + locPhoton1P4_Kin).M();
	massGammaEta[1] = (locPhoton3P4_Kin + locPhoton4P4_Kin + locPhoton2P4_Kin).M();


        if(showOutput){cout << "Combined 4-vectors" << endl;}

        /******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

        // Loop through the analysis actions, executing them in order for the active particle combo

        // Perform_Action gives an error  with our current code probably because it does not have any actions to perform
        //dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
        if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
            continue;


        //if you manually execute any actions, and it fails a cut, be sure to call:
        //dComboWrapper->Set_IsComboCut(true);

	TLorentzVector locPi0EtaP4 = locPhoton1P4 + locPhoton2P4 + locPhoton3P4 + locPhoton4P4;
	TLorentzVector locPi0EtaP4_Kin = locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton3P4_Kin + locPhoton4P4_Kin;

        /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

        /*
           TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
        //for arrays below: 2nd argument is value, 3rd is array index
        //NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
        //So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
        dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
        dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
        */


        // ********************************************************************************************************************************
        // Calculating some variables!!! 
        // ********************************************************************************************************************************
        if(showOutput){cout << "\nCalculating variables to fill the histograms with and to put into histVals..." << endl;}
        //Missing Mass Squared
        locMissingMassSquared = locMissingP4.M2();
        // 4 Momenta related variables  
        locBeamE = locBeamP4_Kin.E();
        if (abs(locBeamE-Ebeam_thrown)<0.0001){ beamPhotonMatchToThrown=true; } 
        else { beamPhotonMatchToThrown=false; }

        locEtaMass_Kin = locEtaP4_Kin.M();
        //double locEtaMass_KinFit = locEtaP4.M();
        locPi0Mass_Kin = locPi0P4_Kin.M();
        //double locPi0Mass_KinFit = locPi0P4.M();

        locEtaMass = locEtaP4.M();
        locPi0Mass = locPi0P4.M();

        mismatchPairMass_13 = mismatchPair_13.M();
        mismatchPairMass_14 = mismatchPair_14.M();
        mismatchPairMass_23 = mismatchPair_23.M();
        mismatchPairMass_24 = mismatchPair_24.M();
        mismatchPairMass_132 = mismatchPair_132.M();
        mismatchPairMass_134 = mismatchPair_134.M();
        mismatchPairMass_142 = mismatchPair_142.M();
        mismatchPairMass_143 = mismatchPair_143.M();
        mismatchPairMass_231 = mismatchPair_231.M();
        mismatchPairMass_234 = mismatchPair_234.M();
        mismatchPairMass_241 = mismatchPair_241.M();
        mismatchPairMass_243 = mismatchPair_243.M();
        mismatchPairMass_341 = mismatchPair_341.M();
        mismatchPairMass_342 = mismatchPair_342.M();
        mismatchPairMass_123 = mismatchPair_123.M();
        mismatchPairMass_124 = mismatchPair_124.M();
        rejectPi0PairMismatch=
            (mismatchPairMass_13<mismatchPi0Mean-2*mismatchPi0Std)*(mismatchPairMass_13>mismatchPi0Mean+2*mismatchPi0Std) || 
            (mismatchPairMass_24<mismatchPi0Mean-2*mismatchPi0Std)*(mismatchPairMass_24>mismatchPi0Mean+2*mismatchPi0Std) || 
            (mismatchPairMass_14<mismatchPi0Mean-2*mismatchPi0Std)*(mismatchPairMass_14>mismatchPi0Mean+2*mismatchPi0Std) || 
            (mismatchPairMass_23<mismatchPi0Mean-2*mismatchPi0Std)*(mismatchPairMass_23>mismatchPi0Mean+2*mismatchPi0Std);

        locEtaE_Kin = locEtaP4_Kin.E();
        locPi0E_Kin = locPi0P4_Kin.E();

        locEtaProton_Kin = mixingEtaProton_Kin.M();
        locPi0Proton_Kin = mixingPi0Proton_Kin.M();
        locPi0Eta_Kin = mixingPi0Eta_Kin.M();
	locPi0Eta_resolution = locPi0Eta_Kin-locPi0Eta_thrown;
	if(showOutput){
            cout << "locPi0Eta_Kin: " << locPi0Eta_Kin << endl;
	    cout << "locPi0Eta_thrown: " << locPi0Eta_thrown << endl;
	    cout << "locPi0Eta_resolution: " << locPi0Eta_resolution << endl;
        }
        locPi0Eta = mixingPi0Eta.M();

        // IN THE FOLLOWING SECTION WE WILL CALCUALTE BY OURSELVES THE MASS OF THE PI0 AND THE ETA WITH DIFFERENT STARTING POINTS ( USING THE PROTON X3 VS USING THE TARGET CENTER)
        //
        TLorentzVector anyVertexMomentum1;
        TLorentzVector anyVertexMomentum2;
        TVector3 vecFromVertexToShower;

        // Calculating the pi0 mass using the charged proton x3 as the initial starting point of the photon
        vecFromVertexToShower = locPhoton1X4_Shower.Vect()-locProtonX4_Kin.Vect();
        vecFromVertexToShower *= locPhoton1P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum1.SetVect(vecFromVertexToShower);
        anyVertexMomentum1.SetE(locPhoton1P4.E());
        vecFromVertexToShower = locPhoton2X4_Shower.Vect()-locProtonX4_Kin.Vect();
        vecFromVertexToShower *= locPhoton2P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum2.SetVect(vecFromVertexToShower);
        anyVertexMomentum2.SetE(locPhoton2P4.E());
        locPi0Mass_charged = (anyVertexMomentum1+anyVertexMomentum2).M();

        // Calculating the eta mass using the charged proton x3 as the initial starting point of the photon
        vecFromVertexToShower = locPhoton3X4_Shower.Vect()-locProtonX4_Kin.Vect();
        vecFromVertexToShower *= locPhoton3P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum1.SetVect(vecFromVertexToShower);
        anyVertexMomentum1.SetE(locPhoton3P4.E());
        vecFromVertexToShower = locPhoton4X4_Shower.Vect()-locProtonX4_Kin.Vect();
        vecFromVertexToShower *= locPhoton4P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum2.SetVect(vecFromVertexToShower);
        anyVertexMomentum2.SetE(locPhoton4P4.E());
        locEtaMass_charged = (anyVertexMomentum1+anyVertexMomentum2).M();

        // Calculating the pi0 mass using the target center x3 as the initial starting point of the photon
        vecFromVertexToShower = locPhoton1X4_Shower.Vect()-targetCenter;
        vecFromVertexToShower *= locPhoton1P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum1.SetVect(vecFromVertexToShower);
        anyVertexMomentum1.SetE(locPhoton1P4.E());
        vecFromVertexToShower = locPhoton2X4_Shower.Vect()-targetCenter;
        vecFromVertexToShower *= locPhoton2P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum2.SetVect(vecFromVertexToShower);
        anyVertexMomentum2.SetE(locPhoton2P4.E());
        locPi0Mass_target = (anyVertexMomentum1+anyVertexMomentum2).M();

        // Calculating the eta mass using the target center x3 as the initial starting point of the photon
        vecFromVertexToShower = locPhoton3X4_Shower.Vect()-targetCenter;
        vecFromVertexToShower *= locPhoton3P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum1.SetVect(vecFromVertexToShower);
        anyVertexMomentum1.SetE(locPhoton3P4.E());
        vecFromVertexToShower = locPhoton4X4_Shower.Vect()-targetCenter;
        vecFromVertexToShower *= locPhoton4P4.E()/vecFromVertexToShower.Mag();
        anyVertexMomentum2.SetVect(vecFromVertexToShower);
        anyVertexMomentum2.SetE(locPhoton4P4.E());
        locEtaMass_target = (anyVertexMomentum1+anyVertexMomentum2).M();

        if(showMassCalc) {
            cout << "Mass of pi0 for labeled=" << locPi0Mass << "  charged=" << locPi0Mass_charged << "  target=" << locPi0Mass_target << endl;
            cout << "Mass of eta for labeled=" << locEtaMass << "  charged=" << locEtaMass_charged << "  target=" << locEtaMass_target << endl;
            cout << endl;
        }	

        // Event Selction Variables
        //Charged Track
        locPtProton = TMath::Sqrt(TMath::Power(locProtonP4_Kin.Px(),2)+TMath::Power(locProtonP4_Kin.Py(),2));
        locPzProton = locProtonP4_Kin.Pz();
        locPolarAngleProton = locProtonP4_Kin.Theta()*radToDeg;
        locXProton = locProtonX4_Kin.X();
        locYProton = locProtonX4_Kin.Y();
        locRProton = TMath::Sqrt(TMath::Power(locXProton,2)+TMath::Power(locYProton,2));
        locdzProton = locProtonX4_Kin.Z();
        locdEdxCDCProton = dProtonWrapper->Get_dEdx_CDC();
        locdEdxFDCProton = dProtonWrapper->Get_dEdx_FDC();
        locMagP3Proton = TMath::Sqrt(TMath::Power(locProtonP4_Kin.Px(),2)+TMath::Power(locProtonP4_Kin.Py(),2)+TMath::Power(locProtonP4_Kin.Pz(),2));


	if (showOutput){ cout << "Filling photon variables" << endl; } 

        TLorentzVector newPhotonP4Vector;
        TLorentzVector newPhotonX4Vector;
        TLorentzVector newPhotonX4Vector_meas;
        TLorentzVector newPhotonVector_Shower;
        TLorentzVector newPhotonVector_meas;
        DNeutralParticleHypothesis* newPhotonWrapper;
        for (auto iPhoton=0; iPhoton<nPhotons ; ++iPhoton){
              newPhotonP4Vector = allPhotonP4Vectors_Kin[iPhoton];
              newPhotonX4Vector = allPhotonX4Vectors_Kin[iPhoton];
              newPhotonX4Vector_meas = allPhotonX4Vectors[iPhoton];
              newPhotonVector_Shower = allPhoton4Vectors_Shower[iPhoton];
              newPhotonVector_meas = allPhoton4Vectors_meas[iPhoton];
              newPhotonWrapper = allPhotonWrappers[iPhoton];
              //photonThetas[iPhoton] = newPhotonP4Vector.Theta()*radToDeg; // *radToDeg; // IF WE WANT TO BRING BACK THE OLDER CUT, double counted in pPhotonTheta
              //photonPhis[iPhoton] = newPhotonP4Vector.Phi()*radToDeg;
              photonThetas[iPhoton] = newPhotonVector_Shower.Theta()*radToDeg; // *radToDeg; // IF WE WANT TO BRING BACK THE OLDER CUT, double counted in pPhotonTheta
              photonPhis[iPhoton] = newPhotonVector_Shower.Phi()*radToDeg;
              photonEnergies[iPhoton] = newPhotonP4Vector.E();
              photonXs_Kin[iPhoton] = newPhotonX4Vector.X();
              photonYs_Kin[iPhoton] = newPhotonX4Vector.Y();
              photonZs_Kin[iPhoton] = newPhotonX4Vector.Z();
              photonTs_Kin[iPhoton] = newPhotonX4Vector.T();

              photonThetas_fromX4_meas[iPhoton] = newPhotonX4Vector_meas.Theta()*radToDeg;
              photonThetas_meas[iPhoton] = newPhotonVector_meas.Theta()*radToDeg;
              photonThetas_Shower[iPhoton] = newPhotonVector_Shower.Theta()*radToDeg;
              photonXs_Shower[iPhoton] = newPhotonVector_Shower.X();
              photonYs_Shower[iPhoton] = newPhotonVector_Shower.Y();
              photonZs_Shower[iPhoton] = newPhotonVector_Shower.Z();
              photonTs_Shower[iPhoton] = newPhotonVector_Shower.T();
              photonDeltaTs[iPhoton] = newPhotonX4Vector.T()-(dComboWrapper->Get_RFTime() + (newPhotonX4Vector.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
              photonDetectedSyss[iPhoton] = newPhotonWrapper->Get_Detector_System_Timing();
              E1E9_FCAL[iPhoton] = newPhotonWrapper->Get_E1E9_FCAL();
              E9E25_FCAL[iPhoton] = newPhotonWrapper->Get_E9E25_FCAL();
              SumU_FCAL[iPhoton] = newPhotonWrapper->Get_SumU_FCAL();
              SumV_FCAL[iPhoton] = newPhotonWrapper->Get_SumV_FCAL();
              Energy_BCALPreshower[iPhoton] = newPhotonWrapper->Get_Energy_BCALPreshower();
              Energy_BCAL[iPhoton] = newPhotonWrapper->Get_Energy_BCAL();
              SigLong_BCAL[iPhoton] = newPhotonWrapper->Get_SigLong_BCAL();
              SigTheta_BCAL[iPhoton] = newPhotonWrapper->Get_SigTheta_BCAL();
              showerQuality_FCAL[iPhoton] = newPhotonWrapper->Get_Shower_Quality();
              SigTrans_BCAL[iPhoton] = newPhotonWrapper->Get_SigTrans_BCAL();
              DeltaPhi_BCAL[iPhoton] = newPhotonWrapper->Get_TrackBCAL_DeltaPhi();
              DeltaZ_BCAL[iPhoton] = newPhotonWrapper->Get_TrackBCAL_DeltaZ();
              DOCA_FCAL[iPhoton] = newPhotonWrapper->Get_TrackFCAL_DOCA();
        }

        if(showOutput) { cout << "Getting some combo specific variables like CL" << endl; }
        //Kinematic fit variables
        locCLKinFit = dComboWrapper->Get_ConfidenceLevel_KinFit("");
        locUnusedEnergy = dComboWrapper->Get_Energy_UnusedShowers();
        locNumExtraNeutralShowers = Get_NumNeutralHypos()-4;
        locDOFKinFit = dComboWrapper->Get_NDF_KinFit("");
        locChiSqKinFit = dComboWrapper->Get_ChiSq_KinFit("");
        //cout << "locChiSqKinFit: " << locChiSqKinFit << endl;
        locChiSqKinFit_gpi0pi0 = dComboWrapper->Get_ChiSq_KinFit( "gpi0pi0" );
        //TString DOF; DOF.Form("%d", dComboWrapper->Get_NDF_KinFit(""));
        //TString title = "Num DOF: "; title+=DOF.Data();

        if(showOutput) { cout << "Getting some shower variables" << endl; }
        locE1E9_FCAL_proton = dProtonWrapper->Get_E1E9_FCAL();
        locE9E25_FCAL_proton = dProtonWrapper->Get_E9E25_FCAL();
        locSumU_FCAL_proton = dProtonWrapper->Get_SumU_FCAL();
        locSumV_FCAL_proton = dProtonWrapper->Get_SumV_FCAL();
        locEnergy_BCALPreshower_proton = dProtonWrapper->Get_Energy_BCALPreshower();
        locEnergy_BCAL_proton = dProtonWrapper->Get_Energy_BCAL();
        locSigTheta_BCAL_proton = dProtonWrapper->Get_SigTheta_BCAL();
        locSigTrans_BCAL_proton = dProtonWrapper->Get_SigTrans_BCAL();;
        locSigLong_BCAL_proton = dProtonWrapper->Get_SigLong_BCAL();;

        // Timing variables, first section is about the proton
        // Get_RFTime is quoted at the center of the target. We have to shift it outwards.
        // We shift it  by the time it takes for light to travel between the initial target and the recoiled target.  
        RFtime = dComboWrapper->Get_RFTime() + (locProtonX4_Kin.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458;
        RFtime_meas = dComboWrapper->Get_RFTime() + (locProtonX4.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458;
        RFtimeProton = locProtonX4_Kin.T() - RFtime;

        // Might not make sense to shift relative to the BeamX4.Z since we cant really track the photon in the detector. Use Protons maybe
        //double locDeltaTRF = locBeamX4.T() - (dComboWrapper->Get_RFTime() + (locBeamX4.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
        locDeltaTRF = locBeamX4_Kin.T() - RFtime;
        locDeltaTRF_meas = locBeamX4.T() - RFtime;

        // Accidental Subtracting we do a weighted histogram. Where we scale the increment to our histogram by the weight. The weight is 1 for in central beam bunch RF time and -0.5 for off central.
        weightAS = 1;
        if (18>abs(locDeltaTRF) && abs(locDeltaTRF)>2){ weightAS = -0.125;} //about 1/8 since we capture a sample of 8 outside peaks
        if (18<abs(locDeltaTRF)) { weightAS=0; }
        //if (14>abs(locDeltaTRF) && abs(locDeltaTRF) >2){ weightAS = -0.125;} //about 1/8 since we capture a sample of 8 outside peaks
        // Shifted realtive to the proton
        inBCAL = 0; inTOF = 0; inFCAL = 0; inSTART = 0; inSYS_NULL = 0;
        if (dProtonWrapper->Get_Detector_System_Timing() == SYS_BCAL){ inBCAL = 1;}
        if (dProtonWrapper->Get_Detector_System_Timing() == SYS_TOF){ inTOF = 1;}
        if (dProtonWrapper->Get_Detector_System_Timing() == SYS_FCAL){ inFCAL = 1;}
        if (dProtonWrapper->Get_Detector_System_Timing() == SYS_START){ inSTART = 1;}
        if (dProtonWrapper->Get_Detector_System_Timing() == SYS_NULL){ inSYS_NULL = 1;}

        if(showOutput) { cout << "Got timing quantites" << endl; }


        // CALCULATING PAIRWISE QUANTITES BETWEEN PHOTONS.
        // the weird thing about this is that I have to do a cut on the dij3 of all the pairwise combinations. So I have to cacluate it before I can start filling most of the other histograms (that uses that dij3Cut)
        // I also have to worry about uniquess tracking when trying to do the filling here otherwise the dij3 noncut and cut graphs would share the same tracking which is wrong. The one with cuts should only fill the 
        // unique set when the cut has been passed. SO! We have to do our uniqueness tracking out here rather than before we fill. Since we will have done the uniquenss tracking here we do not have to do it in the future
        // when we are trying to fill the histograms.
        int countBothInFCAL = 0;
        int countBothInBCAL = 0;
        int countNotInEither = 0;
        int hist_id=0;
        double dij3;// for the FCAL
        double angle_ij_fcal;
        double angle_ij;
        double deltaZ_ij;
        double deltaPhi_ij;
        double dij;// for the BCAL
        double R_BCALi;
        double z_BCALi;
        double theta_BCALi;
        double R_BCALj;
        double z_BCALj;
        double theta_BCALj;
        double delta_Rij_BCAL;
        double delta_Zij_BCAL;
        double delta_Thetaij_BCAL;
        //std::vector<double> angle_ijVec;
        //std::vector<double> deltaZ_ijVec;
        //std::vector<double> deltaPhi_ijVec;
        std::vector< double> dij3VecFCAL;
        std::vector< double> angleVecFCAL;
        std::vector< double> dijVecBCAL;
        std::vector< map<Particle_t, set<Int_t> > > usedPairIdsFCAL;// will be used to do our uniqueness trackingo over the entire event
        std::vector< map<Particle_t, set<Int_t> > > usedPairIdsBCAL;

        map<Particle_t, set<Int_t> > usingPairIdsFCAL;
        map<Particle_t, set<Int_t> > usingPairIdsBCAL; 


        if(showOutput){ cout << "\nCalculating pairwise quantities between photons, doing uniqueness tracking now instead of later since there will be a vector of data to fill\n-----------------------------------------" << endl;}
        std::vector<TLorentzVector> photonX4 = {locPhoton1X4_Shower, locPhoton2X4_Shower, locPhoton3X4_Shower, locPhoton4X4_Shower};
        if(showOutput){ cout << "Set up photonX4 and photonIds to use for pairwise quantity calculating!" << endl;}
        for (int i=0;i<4;++i){
        	for(int j=i+1; j<4;j++){
        		if(showOutput){ cout << "-- ith,jth photon: " << std::to_string(i) << ", " << std::to_string(j) << endl;}
        		// if both the photons are not in the resppective detector there is no point on continuing to check the histograms + cuts + uniquenss tracking
        		if(photonDetectedSyss[i]==SYS_FCAL && photonDetectedSyss[j]==SYS_FCAL){
        			++countBothInFCAL;
        			usingPairIdsFCAL.clear();
                                usingPairIdsFCAL[Unknown].insert(locBeamID);
        			usingPairIdsFCAL[Gamma].insert(photonIds[i]);
        			usingPairIdsFCAL[Gamma].insert(photonIds[j]);
                        	dij3 = (photonX4[i]-photonX4[j]).Vect().Mag();
        			dij3VecFCAL.push_back(dij3);
                        	angle_ij_fcal = (photonX4[i].Vect().Angle(photonX4[j].Vect()))*radToDeg;
                                angleVecFCAL.push_back(angle_ij_fcal);
                                usedPairIdsFCAL.push_back(usingPairIdsFCAL);
        		}
        		else if(photonDetectedSyss[i]==SYS_BCAL && photonDetectedSyss[j]==SYS_BCAL){
        			++countBothInBCAL;
        			usingPairIdsBCAL.clear();
                                usingPairIdsBCAL[Unknown].insert(locBeamID);
        			usingPairIdsBCAL[Gamma].insert(photonIds[i]);
        			usingPairIdsBCAL[Gamma].insert(photonIds[j]);
                                usedPairIdsBCAL.push_back(usingPairIdsBCAL);
                        	angle_ij = (photonX4[i].Vect().Angle(photonX4[j].Vect()))*radToDeg;
                        	if (angle_ij > 180) { angle_ij = 360-angle_ij; }
                        	deltaZ_ij = abs(photonX4[i].Z()-photonX4[j].Z());
                        	deltaPhi_ij = abs(photonX4[i].Phi()-photonX4[j].Phi())*radToDeg;
                        	if (deltaPhi_ij > 180) { deltaPhi_ij = 360-deltaPhi_ij; }
        			if(showOutput){ cout << "BOTH IN BCAL - Got BCAL angle, deltaZ, deltaPhi, distanceOnCylinder!" << endl;}
                        	//dij3 = (photonX4[i]-photonX4[j]).Vect().Mag();

        			// a better metric to use would be to use dij which is the distance conneting two points on the cylindrical BCAL BUT constraining to be on the surface to tell us about separation
        			// 65 < R < 90 cm is the inner and outer dimensions of the BCAL
        			R_BCALi = TMath::Sqrt(TMath::Sq(photonX4[i].X())+TMath::Sq(photonX4[i].Y()));
        			z_BCALi = photonX4[i].Z();
        			theta_BCALi = TMath::ATan(photonX4[i].Y()/photonX4[i].X());
        			R_BCALj = TMath::Sqrt(TMath::Sq(photonX4[j].X())+TMath::Sq(photonX4[j].Y()));
        			z_BCALj = photonX4[j].Z();
        			theta_BCALj = TMath::ATan(photonX4[j].Y()/photonX4[j].X());
        			delta_Rij_BCAL = R_BCALj - R_BCALi;
        			delta_Zij_BCAL = z_BCALj - z_BCALi;
        			delta_Thetaij_BCAL = theta_BCALj - theta_BCALi;
        			dij = TMath::Sqrt(TMath::Sq(delta_Rij_BCAL)+TMath::Sq((65+delta_Rij_BCAL)*delta_Thetaij_BCAL)+TMath::Sq(delta_Zij_BCAL));

                                if(showOutput)cout << "BCAL dij: " << dij << endl;
        			dijVecBCAL.push_back(dij);
        			//angle_ijVec[hist_id].push_back(angle_ij);
        			//deltaZ_ijVec[hist_id].push_back(deltaZ_ij);
        		}
        		else if((photonDetectedSyss[i]==SYS_BCAL && photonDetectedSyss[j]==SYS_FCAL) || (photonDetectedSyss[i]==SYS_FCAL && photonDetectedSyss[j]==SYS_BCAL)){
                                if(showOutput)cout << "SPLIT" << endl;
        			++countNotInEither;
        			if(showOutput){ cout << "ONE IN FCAL ONE IN BCAL" << endl;}
        		}
        		else { cout << "\n\n*************************\nERROR - Pairwise both in FCAL, both in FCAL, in either, is not complete!\n***********************************" << endl;}
        	}
        }
        if((countBothInBCAL+countNotInEither+countBothInFCAL)==6){ if(showOutput) {cout << "\n\n*********************\nPairwise photons counted correctly!\n*******************\nNum both in FCAL: " + std::to_string(countBothInFCAL) +"\nNum both in BCAL: " + std::to_string(countBothInBCAL) + "\nNum one in either: " + std::to_string(countNotInEither)<< endl; }}
        else { cout << "\n\n*********************\nERROR - Pairwise photons not counted correctly!\n*******************\n" << endl;}

        // *********************** PHOTON PAIR PLOTS ************************
        if(showOutput){ cout << endl << endl; }
        if(showOutput){ cout << "Calculating prodplanephi and cosTheta in different system" << endl;}
        // calculating phi
        //double decayPlanePhi = dAnalysisUtilities.Calc_DecayPlanePsi_Vector_3BodyDecay(locBeamE, Proton, locProtonP4_Kin, mixingPi0Eta_Kin, locPi0P4_Kin, locEtaP4_Kin, locDecayPlaneTheta);
        double prodPlanePhi_eta = locEtaP4_Kin.Phi(); //dAnalysisUtilities.Calc_ProdPlanePhi_Pseudoscalar(locBeamE, Proton, locEtaP4_Kin);
        double prodPlanePhi_pi0 = locPi0P4_Kin.Phi(); //dAnalysisUtilities.Calc_ProdPlanePhi_Pseudoscalar(locBeamE, Proton, locPi0P4_Kin);
	if(showOutput)cout << "prodPlanePhi_eta: " << prodPlanePhi_eta << endl;
        locPhi_eta = prodPlanePhi_eta*radToDeg;
        locPhi_pi0 = prodPlanePhi_pi0*radToDeg;
	if(showOutput)cout << "locPhi_eta: " << locPhi_eta << endl;

        // Calculating kinematic variables like t and cosTheta
        mandelstam_t = -(dTargetP4-locProtonP4_Kin).M2();
        mandelstam_t_pe = -(locBeamP4_Kin-mixingPi0Eta).M2();
        mandelstam_teta = -(locBeamP4_Kin-locEtaP4_Kin).M2();
        mandelstam_tpi0 = -(locBeamP4_Kin-locPi0P4_Kin).M2();
        mandelstam_teta_meas = -(locBeamP4-locEtaP4).M2();
        mandelstam_teta_meas = -(locBeamP4-locPi0P4).M2();
	idx_t_eta = (int)( (mandelstam_teta-tMin)/tStep ); 
	idx_t_pi0 = (int)( (mandelstam_tpi0-tMin)/tStep ); 
	idx_m = (int)( (locPi0Eta-mMin)/mStep );
	if ( mandelstam_teta < tMin || locPi0Eta < mMin || mandelstam_teta>tMax || locPi0Eta>mMax ) {
		teta_recCounts = -1;
	}
	else {
		teta_recCounts = idx_t_eta+num_tBins*idx_m;
	}
	if ( mandelstam_tpi0 < tMin || locPi0Eta < mMin || mandelstam_tpi0>tMax || locPi0Eta>mMax ) {
		tpi0_recCounts = -1;
	}
	else {
		tpi0_recCounts = idx_t_pi0+num_tBins*idx_m;
	}
	if(showOutput)cout << "teta_recCounts: " << teta_recCounts << endl;
	if(showOutput)cout << "tpi0_recCounts: " << tpi0_recCounts << endl;

        //////////////////////////////////////
        //////////////////////////////////////
        //////////////////////////////////////
        //////////////////////// CACULATING SOME ANGLES IN THE OMEGA REST FRAME //////////////
        TLorentzVector ph123_photon1P4 = locPhoton1P4_Kin;
        TLorentzVector ph123_photon2P4 = locPhoton2P4_Kin;
        TLorentzVector ph123_photon3P4 = locPhoton3P4_Kin;
        TLorentzVector ph123_photon4P4 = locPhoton4P4_Kin;
        TLorentzVector ph124_photon1P4 = locPhoton1P4_Kin;
        TLorentzVector ph124_photon2P4 = locPhoton2P4_Kin;
        TLorentzVector ph124_photon3P4 = locPhoton3P4_Kin;
        TLorentzVector ph124_photon4P4 = locPhoton4P4_Kin;
        TLorentzVector ph123_vec = locPhoton1P4_Kin+locPhoton2P4_Kin+locPhoton3P4_Kin; 
        TLorentzVector ph124_vec = locPhoton1P4_Kin+locPhoton2P4_Kin+locPhoton4P4_Kin; 
        ph123_photon1P4.Boost(-ph123_vec.BoostVector());
        ph123_photon2P4.Boost(-ph123_vec.BoostVector());
        ph123_photon3P4.Boost(-ph123_vec.BoostVector());
        ph123_photon4P4.Boost(-ph123_vec.BoostVector());
        ph124_photon1P4.Boost(-ph124_vec.BoostVector());
        ph124_photon2P4.Boost(-ph124_vec.BoostVector());
        ph124_photon3P4.Boost(-ph124_vec.BoostVector());
        ph124_photon4P4.Boost(-ph124_vec.BoostVector());
        ph123Rest_angle_pi0_g3=ph123_photon3P4.Vect().Angle( (ph123_photon1P4+ph123_photon2P4).Vect() )*radToDeg;
        ph124Rest_angle_pi0_g4=ph124_photon4P4.Vect().Angle( (ph124_photon1P4+ph124_photon2P4).Vect() )*radToDeg;

        ph123Rest_angle_g12=ph123_photon1P4.Vect().Angle( ph123_photon2P4.Vect() )*radToDeg;
        ph123Rest_angle_g13=ph123_photon1P4.Vect().Angle( ph123_photon3P4.Vect() )*radToDeg;
        ph123Rest_angle_g14=ph123_photon1P4.Vect().Angle( ph123_photon4P4.Vect() )*radToDeg;
        ph123Rest_angle_g23=ph123_photon2P4.Vect().Angle( ph123_photon3P4.Vect() )*radToDeg;
        ph123Rest_angle_g24=ph123_photon2P4.Vect().Angle( ph123_photon4P4.Vect() )*radToDeg;
        ph123Rest_angle_g34=ph123_photon3P4.Vect().Angle( ph123_photon4P4.Vect() )*radToDeg;

        ph124Rest_angle_g12=ph124_photon1P4.Vect().Angle( ph124_photon2P4.Vect() )*radToDeg;
        ph124Rest_angle_g13=ph124_photon1P4.Vect().Angle( ph124_photon3P4.Vect() )*radToDeg;
        ph124Rest_angle_g14=ph124_photon1P4.Vect().Angle( ph124_photon4P4.Vect() )*radToDeg;
        ph124Rest_angle_g23=ph124_photon2P4.Vect().Angle( ph124_photon3P4.Vect() )*radToDeg;
        ph124Rest_angle_g24=ph124_photon2P4.Vect().Angle( ph124_photon4P4.Vect() )*radToDeg;
        ph124Rest_angle_g34=ph124_photon3P4.Vect().Angle( ph124_photon4P4.Vect() )*radToDeg;
        //////////////////////////////////////
        //////////////////////////////////////
        //////////////////////////////////////


        // We will also calculate the cos theta in the lab frame. This is used to see if we are actually getting the right amount of events in the detector
        theta_pi0_lab = locPi0P4.Theta()*radToDeg;
        theta_eta_lab = locEtaP4.Theta()*radToDeg;
        phi_pi0eta_lab = locPi0EtaP4.Phi()*radToDeg;

        // Determine the CM vectors to be used to calculate the angles in the Hel frame
        // also calculating the cosTheta of pi0eta system, pi0, eta  in the CM frame
        TLorentzVector cm_vec = locBeamP4_Kin+dTargetP4;
        TLorentzVector cm_vec_meas = locBeamP4+dTargetP4;
        TLorentzVector mixingPi0Eta_cm = mixingPi0Eta_Kin;
        TLorentzVector mixingPi0Eta_cm_meas = mixingPi0Eta;
        TLorentzVector pi0_cm = locPi0P4_Kin;
        TLorentzVector pi0_cm_meas = locPi0P4;
        TLorentzVector eta_cm = locEtaP4_Kin;
        TLorentzVector eta_cm_meas = locEtaP4;
        TLorentzVector beam_cm = locBeamP4_Kin;
        TLorentzVector beam_cm_meas = locBeamP4;
        TLorentzVector recoil_cm = locProtonP4_Kin;
        TLorentzVector recoil_cm_meas = locProtonP4;
        TLorentzVector largestEPhotonInPi0_CM;
        if (locPhoton1P4_Kin.E() > locPhoton2P4_Kin.E()){
            largestEPhotonInPi0_CM = locPhoton1P4_Kin;
        }
        else {
            largestEPhotonInPi0_CM = locPhoton2P4_Kin;
        }
        largestEPhotonInPi0_CM.Boost(-cm_vec.BoostVector());
        cosTheta_largestEinPi0_CM = largestEPhotonInPi0_CM.CosTheta();

        mixingPi0Eta_cm.Boost(-cm_vec.BoostVector());
        mixingPi0Eta_cm_meas.Boost(-cm_vec_meas.BoostVector());
        pi0_cm.Boost(-cm_vec.BoostVector());
        eta_cm.Boost(-cm_vec.BoostVector());
        eta_cm_meas.Boost(-cm_vec_meas.BoostVector());
        beam_cm.Boost(-cm_vec.BoostVector());
        beam_cm_meas.Boost(-cm_vec_meas.BoostVector());
        recoil_cm.Boost(-cm_vec.BoostVector());
        recoil_cm_meas.Boost(-cm_vec.BoostVector());

        cosTheta_pi0eta_CM = mixingPi0Eta_cm.CosTheta();
        cosTheta_pi0eta_CM_meas = mixingPi0Eta_cm_meas.CosTheta();
        cosTheta_pi0_CM = pi0_cm.CosTheta();
        cosTheta_eta_CM = eta_cm.CosTheta();
        theta_pi0_CM = pi0_cm.Theta();
        theta_eta_CM = eta_cm.Theta();
        phi_pi0eta_CM = mixingPi0Eta_cm.Phi()*radToDeg;
        phi_pi0_CM = pi0_cm.Phi()*radToDeg;
        phi_eta_CM = eta_cm.Phi()*radToDeg;

        mandelstam_t0_oldForm = -(TMath::Power(-mixingPi0Eta_Kin.M2()/(2*(locBeamP4_Kin+dTargetP4).M()),2)-(beam_cm-mixingPi0Eta_cm).M2());
	double term1 = -mixingPi0Eta_Kin.M2()/2/cm_vec.M();
	double shift = (beam_cm.Vect().Mag()-mixingPi0Eta_cm.Vect().Mag())*(beam_cm.Vect().Mag()-mixingPi0Eta_cm.Vect().Mag());
        mandelstam_t0 = -(TMath::Power(-mixingPi0Eta_Kin.M2()/(2*(locBeamP4_Kin+dTargetP4).M()),2)-TMath::Power(beam_cm.Vect().Mag()-mixingPi0Eta_cm.Vect().Mag(),2));
        mandelstam_tp = mandelstam_t-mandelstam_t0;
        mandelstam_tp_pe = mandelstam_t_pe-mandelstam_t0;
	mandelstam_tp_oldForm = mandelstam_t-mandelstam_t0_oldForm;


	double PBeam_cm = locBeamP4_Kin.Vect().Mag()*dTargetP4.M()/cm_vec.M();
	double EX_cm = (cm_vec.M2()+mixingPi0Eta_cm.M2()-locProtonP4_Kin.M2())/2/cm_vec.M();
	double PX_cm = sqrt(EX_cm*EX_cm-mixingPi0Eta_Kin.M2());
	double shiftFromEqn = (PBeam_cm-PX_cm)*(PBeam_cm-PX_cm);
	if(showOutput)cout << "Shift, Shift from eqn = " << shift << ", " << shiftFromEqn << endl;
	if(showOutput)cout << "Mandelstam_t = " << mandelstam_t << endl;


	
	// -----------------------------------------------------------
	// -------------------- DEFINING SOME ANGLES ------------------
	// -----------------------------------------------------------

        TLorentzVector beam_res = beam_cm;
        TLorentzVector beam_res_meas = beam_cm_meas;
        TLorentzVector recoil_res = recoil_cm;
        TLorentzVector recoil_res_meas = recoil_cm_meas;
        TLorentzVector eta_res = eta_cm;
        TLorentzVector eta_res_meas = eta_cm_meas;
        TLorentzVector pi0_res = pi0_cm;
        TLorentzVector mixingPi0Eta_res = mixingPi0Eta_cm;

        // Boost using mixingPi0Eta_cm instead of res... since that is what we are changing.
        beam_res.Boost(-mixingPi0Eta_cm.BoostVector());
        recoil_res.Boost(-mixingPi0Eta_cm.BoostVector());
        eta_res.Boost(-mixingPi0Eta_cm.BoostVector());
        pi0_res.Boost(-mixingPi0Eta_cm.BoostVector());
        mixingPi0Eta_res.Boost(-mixingPi0Eta_cm.BoostVector());
        TLorentzVector largestEPhotonInEta_GJ;
        if (locPhoton3P4_Kin.E() > locPhoton4P4_Kin.E()){
            largestEPhotonInEta_GJ = locPhoton3P4_Kin;
        }
        else {
            largestEPhotonInEta_GJ = locPhoton4P4_Kin;
        }
        largestEPhotonInEta_GJ.Boost(-cm_vec.BoostVector());
        largestEPhotonInEta_GJ.Boost(-mixingPi0Eta_cm.BoostVector());
        TVector3 largestEPhotonInEta_res_unit = largestEPhotonInEta_GJ.Vect().Unit();

        TVector3 pi0_res_unit = pi0_res.Vect().Unit();	
        TVector3 eta_res_unit = eta_res.Vect().Unit();	
        TVector3 pi0eta_res_unit = mixingPi0Eta_res.Vect().Unit();	
        // Calculate cosTheta, phi in maybe the GJ axes.
        // since we already defined the x,y,z as TVector3 we don't have to do it again.
        TVector3 z = beam_res.Vect().Unit();
        // this y should be the normal of the production plane. If we do a boost in a direction in the production plane the perp direction doesn't change. We could use the beam and the recoiled proton to define the
        // production plane in this new frame. Let us define it in the CM frame. 
        // TVector3 y = mixingPi0Eta_cm.Vect().Cross(beam_cm.Vect()).Unit(); // DEFINITION I HAVE BEEN USING FOR A LONG TIME
        TVector3 y = (beam_cm.Vect().Unit().Cross(-1*recoil_cm.Vect().Unit())).Unit();
        TVector3 x = y.Cross(z).Unit();

        // recylcing our use of norm and the defined angles
        TVector3 angles_pi0;
        TVector3 angles_eta;
        TVector3 angles_pi0eta;
        angles_pi0.SetXYZ ( pi0_res_unit.Dot(x), pi0_res_unit.Dot(y), pi0_res_unit.Dot(z) );
        angles_eta.SetXYZ ( eta_res_unit.Dot(x), eta_res_unit.Dot(y), eta_res_unit.Dot(z) );
        angles_pi0eta.SetXYZ ( pi0eta_res_unit.Dot(x), pi0eta_res_unit.Dot(y), pi0eta_res_unit.Dot(z) );
        TVector3 angles_largestEinEta_GJ = largestEPhotonInEta_res_unit; 
        angles_largestEinEta_GJ.SetXYZ( angles_largestEinEta_GJ.Dot(x), angles_largestEinEta_GJ.Dot(y), angles_largestEinEta_GJ.Dot(z) );

        theta_pi0_GJ = angles_pi0.Theta()*radToDeg;
        theta_eta_GJ = angles_eta.Theta()*radToDeg;
        cosTheta_pi0_GJ = angles_pi0.CosTheta();
        cosTheta_eta_GJ = angles_eta.CosTheta();
        cosTheta_largestEinEta_GJ = angles_largestEinEta_GJ.CosTheta();
        phi_pi0_GJ = angles_pi0.Phi()*radToDeg;
        phi_eta_GJ = angles_eta.Phi()*radToDeg;	

	// some extra variables from the measured 4 vectors
        beam_res_meas.Boost(-mixingPi0Eta_cm_meas.BoostVector());
       	eta_res_meas.Boost(-mixingPi0Eta_cm_meas.BoostVector());
        z = beam_res_meas.Vect().Unit();
        // y = mixingPi0Eta_cm_meas.Vect().Cross(beam_cm_meas.Vect()).Unit(); // DEFINITION I HAVE BEEN USING FOR A LONG TIME
        y = (beam_cm_meas.Vect().Unit().Cross(-1*recoil_cm_meas.Vect().Unit())).Unit();
        x = y.Cross(z).Unit();
        TVector3 eta_res_unit_meas = eta_res_meas.Vect().Unit();	
        angles_eta.SetXYZ ( eta_res_unit_meas.Dot(x), eta_res_unit_meas.Dot(y), eta_res_unit_meas.Dot(z) );
        cosTheta_eta_GJ_meas = angles_eta.CosTheta();
        phi_eta_GJ_meas = angles_eta.Phi()*radToDeg;	

        /////////////////////// HELICITY ANGLES
        z = -1. * recoil_res.Vect().Unit();
        y = (beam_cm.Vect().Unit().Cross(-1*recoil_cm.Vect().Unit())).Unit();
        x = y.Cross(z);

        TVector3 angles( (eta_res_unit).Dot(x),
                (eta_res_unit).Dot(y),
                (eta_res_unit).Dot(z) );

        cosTheta_eta_hel = angles.CosTheta();
        theta_eta_hel = angles.Theta()*radToDeg;
        phi_eta_hel = angles.Phi()*radToDeg;
        
        TVector3 eps(TMath::Cos(locPolarizationAngle*TMath::DegToRad()), TMath::Sin(locPolarizationAngle*TMath::DegToRad()), 0.0); // beam polarization vector
        Phi = TMath::ATan2(y.Dot(eps), beam_cm.Vect().Unit().Dot(eps.Cross(y)))*radToDeg;



        //cout << "*cosTheta*" << endl;
        //cout << "Boosting first into the CM frame then boosting into the resonance frame with the resonances momentum vector" << endl;
        //cout << "The resonance's momentum vector in the boosted frame without axis redefinition:" << endl;
        //cout << mixingPi0Eta_res.Vect().X() << ","<< mixingPi0Eta_res.Vect().Y() << ","<< mixingPi0Eta_res.Vect().Z() << endl;
        //cout << "The resonance's momentum UNIT vector in the boosted frame without axis redefinition:" << endl;
        //cout << mixingPi0Eta_res.Vect().Unit().X() << ","<< mixingPi0Eta_res.Vect().Unit().Y() << ","<< mixingPi0Eta_res.Vect().Unit().Z() << endl;
        //cout << "Projecting into the the GJ angles we get the following Vector:" << endl;
        //cout << angles_pi0eta.X() << ","<< angles_pi0eta.Y() << ","<< angles_pi0eta.Z() << endl;
        //cout << angles_pi0eta.CosTheta() << endl;
        //cout << mixingPi0Eta_res.Vect().Dot(z)/mixingPi0Eta_res.Vect().Mag() << endl;
        //cout << "*cosTheta*" << endl;


        // Angle() returns the angle between two vectors where the angle is in radians.
        angleBetweenPi0Eta = pi0_res_unit.Angle(eta_res_unit)*radToDeg;
        if(showOutput) { cout << "Filling inCone bools" << endl; }
        for (int iDelta=0; iDelta<2; ++iDelta){
            // THere should be some redundancy here. since mixingPi0Eta is the sum fo the individual 4 vectors it must be that the pi0 and eta decaay back to back in the resonance rest frame.
            // Therefore largeAngle=180 always and only 1 of the mesons have to be in the cone. As a check it is good to see if the number of times it passes these 3 degenerate conditions are the same!
            pi0_inCone[iDelta] = abs(theta_pi0_GJ)>90-delta[iDelta] && abs(theta_pi0_GJ)<90+delta[iDelta] && abs(phi_pi0_GJ)>90-delta[iDelta] && abs(phi_pi0_GJ)<90+delta[iDelta]; 
            eta_inCone[iDelta] = abs(theta_eta_GJ)>90-delta[iDelta] && abs(theta_eta_GJ)<90+delta[iDelta] && abs(phi_eta_GJ)>90-delta[iDelta] && abs(phi_eta_GJ)<90+delta[iDelta]; 
            largeAngle[iDelta] = angleBetweenPi0Eta>180-2*delta[iDelta]; // that is the closest they can be
            withinCone[iDelta] = largeAngle[iDelta] && pi0_inCone[iDelta] && eta_inCone[iDelta]; 
        }
        q = sqrt(pi0_cm.Pz()*pi0_cm.Pz() + eta_cm.Pz()*eta_cm.Pz() + recoil_cm.Pz()*recoil_cm.Pz());
        pi0_cmZ = pi0_cm.Pz();
        eta_cmZ = eta_cm.Pz();
        recoil_cmZ = recoil_cm.Pz();
        omega = 0;
	// the orientation is the same as defined in https://www.dropbox.com/s/iunmarm97zd8yvu/JLab-cake-Feb19.pdf?dl=0
        if(recoil_cmZ > 0 && pi0_cmZ < 0 && eta_cmZ > 0) {
            // in quadrant I, so asin returns correct answer
            omega = asin(sqrt(3/2.)*recoil_cmZ/q);
        }
        if(recoil_cmZ > 0 && pi0_cmZ < 0 && eta_cmZ < 0) {
	    // quadrant II 
            // use pi0_cmZ as X-axis since gives simpler quadrant definition
            omega = asin(sqrt(3/2.)*pi0_cmZ/q) + 2*TMath::Pi()/3;
        }
        if(recoil_cmZ > 0 && pi0_cmZ > 0 && eta_cmZ < 0) {
            // in quadrant III, so need to subtract asin from PI to get correct angle
            omega = TMath::Pi() - asin(sqrt(3/2.)*recoil_cmZ/q);
        }
        if(recoil_cmZ < 0 && pi0_cmZ > 0 && eta_cmZ < 0) {
            // in quadrant IV, so need to subtract asin from PI to get correct angle
            omega = TMath::Pi() - asin(sqrt(3/2.)*recoil_cmZ/q);
        }
        if(recoil_cmZ < 0 && pi0_cmZ > 0 && eta_cmZ > 0) {
	    // quardrant V
            // use pi0_cmZ as X-axis since gives simpler quadrant definition
            omega = (TMath::Pi() - asin(sqrt(3/2.)*pi0_cmZ/q)) + 2*TMath::Pi()/3;
        }
        if(recoil_cmZ < 0 && pi0_cmZ < 0 && eta_cmZ > 0) {
            // in quadrant VI, so asin returns the correct answer but we dont want negative numbers. arcsin would return [-90,90]. Lets add 360 to this
            omega = asin(sqrt(3/2.)*recoil_cmZ/q)+2*TMath::Pi();
        }

        // compute vanHove_x and y before changing omega from rad to degrees
        vanHove_x = q*cos(omega);
        vanHove_y = q*sin(omega);

        omega = omega*radToDeg;
        if(showOutput) { cout << "Got Vanhove quantities" << endl; }

        // eow that we gotten the PIDs we can also do our uniqueness tracking setup now.
        if(showOutput){cout << "Setting up the sets and maps to be used in uniqueness tracking" << endl;}
 
        // ********************************************************************************************************************************
        // Recomputing some booleans since they could not have been sucessfully initialized without running though the loop first.
        // ********************************************************************************************************************************
        if(showOutput){cout << "Re-evaluating bool variables to fill into histCuts" << endl;}

        // Beam cuts
        pBeamAsymE = locBeamE > 8.0 && locBeamE < 8.7;
        pBeamE30to46 = locBeamE > 3.0 && locBeamE < 4.6;
        pBeamE46to62 = locBeamE > 4.6 && locBeamE < 6.2;
        pBeamE62to78 = locBeamE > 6.2 && locBeamE < 7.8;
        pBeamE78to94 = locBeamE > 7.8 && locBeamE < 9.4;
        pBeamE94to11 = locBeamE > 9.4 && locBeamE < 11.0;
        pBeamE8GeVPlus = locBeamE > 8.0;
        pBeamE82to88 = locBeamE > 8.2 && locBeamE< 8.8;

        dEdxCut = TMath::Power(10,-6)*(0.9+TMath::Exp(3.0-3.5*(locMagP3Proton+0.05)/.93827)); // The bigger then number multiplying MagP3 the sharper the cut.
        // pi0Eta specifc cuts
        pEtaProtonBaryonCut = locEtaProton_Kin >= etaProtonBaryonCut;
        ppi0ProtonBaryonCut = locPi0Proton_Kin >= pi0ProtonBaryonCut;
        pBeamE = locBeamE  >= beamECut;

        iLowMass = lowMass; // reinitialize
        iUpMass = lowMass+binScale; // reinitialize
        for (int bin=0; bin<numBinsMass; ++bin){
            p_phiMassBinned[bin] = locPi0Eta_Kin>iLowMass && locPi0Eta_Kin<iUpMass; 
            iLowMass+=binScale;
            iUpMass+=binScale;
        }
        iLowMass_t = lowMass_t; // reinitialize
        iUpMass_t = lowMass_t+binScale_t; // reinitialize
        for (int bin=0; bin<numBinsMass_t; ++bin){
            // we use this set of bools to bin the t distribution in bins of M(pi0eta)
            p_tMassBinned[bin] = locPi0Eta_Kin>iLowMass_t && locPi0Eta_Kin<iUpMass_t; 
            iLowMass_t+=binScale_t;
            iUpMass_t+=binScale_t;
        }
        // THIS IS FOR THE pi0Mass binned in E(pi0) selecting out the f2 peak
        for (int bin=0; bin<numRegions_E; ++bin){
            p_pi0MassPi0Eregion_1[bin] = locPi0E_Kin<iUpE[bin] && locPi0E_Kin>iLowE[bin];
            p_pi0MassPi0Eregion_2[bin] = locEtaE_Kin<iUpE[bin] && locEtaE_Kin>iLowE[bin];
        }
        pSelectf2 = locPi0Eta_Kin<1.5 && locPi0Eta_Kin>1;

        outsideEllipse_loose=0;
        if (locPi0Mass_Kin<=(ellipseX-ellipseXr_loose) || locPi0Mass_Kin>=(ellipseX+ellipseXr_loose)){ outsideEllipse_loose = 1;}
        double ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseX)/ellipseXr_loose,2))*TMath::Power(ellipseYr_loose,2));
        if (locEtaMass_Kin>=(ellipseDeltaY+ellipseY) || locEtaMass_Kin<=(-ellipseDeltaY+ellipseY)){ outsideEllipse_loose = 1; }
        pinsideEllipse_loose = !outsideEllipse_loose;

        outsideEllipse = 0;
        if (locPi0Mass_Kin<=(ellipseX-ellipseXr) || locPi0Mass_Kin>=(ellipseX+ellipseXr)){ outsideEllipse = 1;}
        ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseX)/ellipseXr,2))*TMath::Power(ellipseYr,2));
        if (locEtaMass_Kin>=(ellipseDeltaY+ellipseY) || locEtaMass_Kin<=(-ellipseDeltaY+ellipseY)){ outsideEllipse = 1; }
        pinsideEllipse = !outsideEllipse;

        outsideEllipse_13_24 = 0;
        if (mismatchPairMass_13<=(ellipseX-ellipseXr) || mismatchPairMass_13>=(ellipseX+ellipseXr)){ outsideEllipse_13_24 = 1;}
        ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((mismatchPairMass_13-ellipseX)/ellipseXr,2))*TMath::Power(ellipseXr,2));
        if (mismatchPairMass_24>=(ellipseDeltaY+ellipseX) || mismatchPairMass_24<=(-ellipseDeltaY+ellipseX)){ outsideEllipse_13_24 = 1; }

        outsideEllipse_14_23 = 0;
        if (mismatchPairMass_14<=(ellipseX-ellipseXr) || mismatchPairMass_14>=(ellipseX+ellipseXr)){ outsideEllipse_14_23 = 1;}
        ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((mismatchPairMass_14-ellipseX)/ellipseXr,2))*TMath::Power(ellipseXr,2));
        if (mismatchPairMass_23>=(ellipseDeltaY+ellipseX) || mismatchPairMass_23<=(-ellipseDeltaY+ellipseX)){ outsideEllipse_14_23 = 1; }


        omegaCut=
            (mismatchPairMass_13<0.15)*(mismatchPairMass_24<0.15) || (mismatchPairMass_14<0.15)*(mismatchPairMass_23<0.15) ||
            (mismatchPairMass_13<0.12)*(mismatchPairMass_23<0.12) || (mismatchPairMass_14<0.12)*(mismatchPairMass_24<0.12);

        // General Cuts
        pUnusedEnergy = locUnusedEnergy <= unusedEnergyCut;
	pLooseUnusedEnergy = locUnusedEnergy <= 0.5;
        pChiSq = locChiSqKinFit <= ChiSqCut;	
	pLooseChiSq = locChiSqKinFit <= 500;
        chiSq100 = locChiSqKinFit <= 100;
        pDeltaTRF = abs(locDeltaTRF) <= RFCut;
        pMissingMassSquared = locMissingMassSquared <= MMsqCut;

        // Neutral Cuts
        pdij3pass = 1;
        // we must index the dij3Vec which is a vector with size id. id_noCutDij3 would be the index of the no cut dij3 histogram
        // which we would use to calculate whether the dij3 cut passes. Probably the most fair way to calculate it rather than
        // using the one where we have already uniquenss trakced it.
        for (std::size_t i=0; i<dij3VecFCAL.size();i++){
            // since we want all the dij3FCAL to be <= 12.5 cm we can check if any of them is > 12.5 since its simpler. 
            if ( dij3VecFCAL[i] <= dijCut){
                pdij3pass = 0;
                break;
            }
        }
        pPhoton1E = photonEnergies[0] >= ECut;
        pPhoton2E = photonEnergies[1] >= ECut;
        pPhoton3E = photonEnergies[2] >= ECut;
        pPhoton4E = photonEnergies[3] >= ECut;
        pPhotonE = pPhoton1E*pPhoton2E*pPhoton3E*pPhoton4E;

        pPhoton1Theta = (photonThetas[0]>=thetaCutMin && photonThetas[0]<=thetaCutMax1) || photonThetas[0]>=thetaCutMax2;
        pPhoton2Theta = (photonThetas[1]>=thetaCutMin && photonThetas[1]<=thetaCutMax1) || photonThetas[1]>=thetaCutMax2;
        pPhoton3Theta = (photonThetas[2]>=thetaCutMin && photonThetas[2]<=thetaCutMax1) || photonThetas[2]>=thetaCutMax2;
        pPhoton4Theta = (photonThetas[3]>=thetaCutMin && photonThetas[3]<=thetaCutMax1) || photonThetas[3]>=thetaCutMax2;
        pPhotonTheta = pPhoton1Theta*pPhoton2Theta*pPhoton3Theta*pPhoton4Theta;

        pShowerQuality0 = showerQuality_FCAL[0] > 0.5;
        pShowerQuality1 = showerQuality_FCAL[1] > 0.5;
        pShowerQuality2 = showerQuality_FCAL[2] > 0.5;
        pShowerQuality3 = showerQuality_FCAL[3] > 0.5;
        // Charged Cuts
        pMagP3Proton = locMagP3Proton >= P3Cut;
        pzCutmin = zCutmin <= locdzProton && locdzProton <= zCutmax;
        pRProton = locRProton <= Rcut;
        pdEdxCDCProton = locdEdxCDCProton >= dEdxCut;
        pReg1 = locPzProton>Reg1Xmin && locPtProton>Reg1Ymin && locPzProton<Reg1Xmax && locPtProton<Reg1Ymax;
        pReg2 = locPzProton>Reg2Xmin && locPtProton>Reg2Ymin && locPzProton<Reg2Xmax && locPtProton<Reg2Ymax;
        pReg3 = locPzProton>Reg3Xmin && locPtProton>Reg3Ymin && locPzProton<Reg3Xmax && locPtProton<Reg3Ymax;
        pReg4 = locPzProton>Reg4Xmin && locPtProton>Reg4Ymin && locPzProton<Reg4Xmax && locPtProton<Reg4Ymax;

        // locWherePhoton will be set equal to photonDetectedSyss[N] 
        for (int i = 0; i<4; ++i){
            pPhotonInBCALorFCAL[i] = photonDetectedSyss[i] == SYS_BCAL || photonDetectedSyss[i] == SYS_FCAL;
            pPhotonInFCAL[i] = photonDetectedSyss[i] == SYS_FCAL;
            pPhotonInBCAL[i] = photonDetectedSyss[i] == SYS_BCAL;
        }

        pPi0InFCAL = photonDetectedSyss[0]==SYS_FCAL && photonDetectedSyss[1]==SYS_FCAL;
        pPi0InBCAL = photonDetectedSyss[0]==SYS_BCAL && photonDetectedSyss[1]==SYS_BCAL;
        pPi0InSplit = (photonDetectedSyss[0]==SYS_FCAL && photonDetectedSyss[1]==SYS_BCAL) || (photonDetectedSyss[1]==SYS_FCAL && photonDetectedSyss[0]==SYS_BCAL);
        pEtaInFCAL = photonDetectedSyss[2]==SYS_FCAL && photonDetectedSyss[3]==SYS_FCAL;
        pEtaInBCAL = photonDetectedSyss[2]==SYS_BCAL && photonDetectedSyss[3]==SYS_BCAL;
        pEtaInSplit = (photonDetectedSyss[2]==SYS_FCAL && photonDetectedSyss[3]==SYS_BCAL) || (photonDetectedSyss[3]==SYS_FCAL && photonDetectedSyss[2]==SYS_BCAL);

	if ( pPi0InFCAL ) { pi0DetectedIn = 0; } 
	if ( pPi0InBCAL ) { pi0DetectedIn = 1; } 
	if ( pPi0InSplit ) { pi0DetectedIn = 2; } 
	if ( pEtaInFCAL ) { etaDetectedIn = 0; } 
	if ( pEtaInBCAL ) { etaDetectedIn = 1; } 
	if ( pEtaInSplit ) { etaDetectedIn = 2; } 


        if (omega < 0) {
            pVanHove = omega > -120 && omega < -60; 
        }
        else{
            pVanHove = omega > 240 && omega < 300; 
        }

        ptpLT1 = mandelstam_tp<1; 
        ptpLT05 = mandelstam_tp<0.5; 
        ptGT1 = mandelstam_t>1; 
        ptLT05 = mandelstam_t<0.5; 
        ptGT05LT1 = mandelstam_t<1 && mandelstam_t>0.5; 



        // cut to accept the delta peak in M(pi0proton)
        // we will actually use it to reject the delta peak in the cuts. pi0pi0 should also have this resonance right?
        pMPi0P14 = locPi0Proton_Kin<1.4;

        // Various combinations of cuts, the majority of them will be used just for a few histograms like when showing unused energy graph we will use mUE which
        // removes the UE cut from allGeneralCutsPassed. m prefix basically stands for minus
        dzRP = pMagP3Proton*pzCutmin*pRProton;
        dzR = pzCutmin*pRProton;

        pShowerQuality=pShowerQuality0*pShowerQuality1*pShowerQuality2*pShowerQuality3;
        pShowerQuality=true;
	
	// ***********************
	// cuts applied to all
	bool cata = (locBeamE > 8.2 && locBeamE < 8.8);//*pVH;//chiSq100*pMpi0etaDoubleRegge;//(mandelstam_teta<1)*pMpi0etaDoubleRegge;
	// ***********************
	// temporary cuts
        //pinsideEllipse=pinsideEllipse_loose;
        //ptpLT1 = true;  // do not select on tp
        ptpLT1 = (mandelstam_t>0.1)*(mandelstam_t<1.0); /// Overwriting tpLT1 to use regular mandelstam t for comparisons to Malte and Colin
        pRProton=true;
        //cata=true;
        //ptpLT1=true;
	// ***********************

        mEllipse_pre = cata*ptpLT1*pUnusedEnergy*pChiSq*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;

        // Calculate weights ignoring the other dimension
        if ( locPi0Mass_Kin > pi0Mean-pi0Std*pi0Sig && locPi0Mass_Kin < pi0Mean+pi0Std*pi0Sig ) { weightBSpi0 = 1; } 
        else if ( locPi0Mass_Kin > pi0Mean-pi0Std*(pi0Sig+pi0Skip+pi0SB) && locPi0Mass_Kin < pi0Mean-pi0Std*(pi0Sig+pi0Skip) ) { weightBSpi0 = pi0SBweight; } 
        else if ( locPi0Mass_Kin > pi0Mean+pi0Std*(pi0Sig+pi0Skip) && locPi0Mass_Kin < pi0Mean+pi0Std*(pi0Sig+pi0Skip+pi0SB) ) { weightBSpi0 = pi0SBweight; } 
        else { weightBSpi0 = 0; }
        if ( locEtaMass_Kin > etaMean-etaStd*etaSig && locEtaMass_Kin < etaMean+etaStd*etaSig ) { weightBSeta = 1; } 
        else if ( locEtaMass_Kin > etaMean-etaStd*(etaSig+etaSkip+etaSB) && locEtaMass_Kin < etaMean-etaStd*(etaSig+etaSkip) ) { weightBSeta = etaSBweight; } 
        else if ( locEtaMass_Kin > etaMean+etaStd*(etaSig+etaSkip) && locEtaMass_Kin < etaMean+etaStd*(etaSig+etaSkip+etaSB) ) { weightBSeta = etaSBweight; } 
        else { weightBSeta = 0; }
        weightBS=weightBSpi0*weightBSeta;

        // now that we have defined both the weights we can multiply them together
        weightAS_BS = weightAS*weightBS;
	dTreeInterface->Fill_Fundamental<Float_t>("DataWeight", weightAS_BS, loc_i);

	if ( abs(weightAS_BS)>1){
		if(showOutput)cout << "weightBS,weightAS: " << weightBS << ", " << weightAS << endl;
	}

	++count_combos;
        if(pShowerQuality){ ++count_ShowerQuality; dHist_Cuts->Fill(cutNames[0],1);}
        if(pBeamE8GeVPlus){ ++count_BeamE8GeVPlus; dHist_Cuts->Fill(cutNames[1],1);} 
        if(pUnusedEnergy){ ++count_UnusedEnergy;dHist_Cuts->Fill(cutNames[2],1);}
        if(pChiSq){ ++count_ChiSq;dHist_Cuts->Fill(cutNames[3],1);}
        if(pDeltaTRF){ ++count_DeltaTRF;dHist_Cuts->Fill(cutNames[4],1);}
        if(pdij3pass){ ++count_dij3pass;dHist_Cuts->Fill(cutNames[5],1);}
        if(pPhotonE){ ++count_PhotonE;dHist_Cuts->Fill(cutNames[6],1);}
        if(pPhotonTheta){ ++count_PhotonTheta;dHist_Cuts->Fill(cutNames[7],1);}
        if(pMagP3Proton){ ++count_MagP3Proton;dHist_Cuts->Fill(cutNames[8],1);}
        if(pzCutmin){ ++count_zCutmin;dHist_Cuts->Fill(cutNames[9],1);}
        if(pRProton){ ++count_RProton;dHist_Cuts->Fill(cutNames[10],1);}
        if(pMissingMassSquared){ ++count_MissingMassSquared;dHist_Cuts->Fill(cutNames[11],1);}
        if(pdEdxCDCProton){ ++count_dEdxCDCProton;dHist_Cuts->Fill(cutNames[12],1);}
        if(pinsideEllipse){ ++count_insideEllipse;dHist_Cuts->Fill(cutNames[13],1);}
        if(allGeneralCutsPassed){ ++count_allGeneralCutsPassed;dHist_Cuts->Fill(cutNames[14],1);}
	if(!pMPi0P14) { ++count_MPi0P14; dHist_Cuts->Fill(cutNames[15],1);}
	if(ptpLT1){ ++count_seanResTest; dHist_Cuts->Fill(cutNames[16],1);}


        if(showOutput) {cout << "Start Filling histVals and histCuts" << endl;}

        if(inBox_noOtherCuts[4]){ whichSignalRegion=1; } // signal
        else if(inBox_noOtherCuts[10]){ whichSignalRegion=2; } // corner
        else if(inBox_noOtherCuts[11]){ whichSignalRegion=3; } // eta 
        else if(inBox_noOtherCuts[12]){ whichSignalRegion=4; } // pi0
        else { whichSignalRegion=5; } // skip

	if(showOutput){ cout << "Filling histogram's uniqueness elements" << endl; }

	if ( selectDetector == "FCAL" ) { detectorCut=pEtaInFCAL; }
	else if ( selectDetector == "BCAL" ) { detectorCut=pEtaInBCAL; }
	else if ( selectDetector == "SPLIT" ) { detectorCut=pEtaInSplit; }
	else { detectorCut=true; }

        if (!mEllipse_pre || !detectorCut){
	    if (showOutput) { cout << "Did not pass cut, moving on.... " << endl; }  
            dComboWrapper->Set_IsComboCut(true); continue; 
        }

        /****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

        /*
        //FILL ANY CUSTOM BRANCHES FIRST!!
        Int_t locMyInt_Flat = 7;
        dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

        TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

        for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
        {
        dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
        TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
        }
        */

        dFlatTreeInterface->Fill_Fundamental<Int_t>("mcprocess", mcprocess);
	// ** SOME THROWN VARIBLES **
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_t_thrown", locT_thrown); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Ebeam_thrown", locBeamE_thrown); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Ebeam", locBeamE); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<UChar_t>("locNumUnusedShowers", locNumUnusedShowers);
	// **************************

        // ** Including some basic cuts  and associated variables
        // **************************
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mmsq",locMissingMassSquared);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_momentum",locMagP3Proton);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_z",locdzProton);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_R",locRProton);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_dEdxCDC",locdEdxCDCProton);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_13",mismatchPairMass_13);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_24",mismatchPairMass_24);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_23",mismatchPairMass_23);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_14",mismatchPairMass_14);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("omegaCut",omegaCut);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pVH",pVH);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pPhotonE",pPhotonE);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pPhotonTheta",pPhotonTheta);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pMagP3Proton",pMagP3Proton);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pzCutmin",pzCutmin);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pRProton",pRProton);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pMissingMassSquared",pMissingMassSquared);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("pdEdxCDCProton",pdEdxCDCProton);

        dFlatTreeInterface->Fill_Fundamental<Float_t>("DOFKinFit", locDOFKinFit);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("AccWeight", weightAS);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("BeamAngle", locPolarizationAngle);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("weightASBS", weightAS_BS);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("weightBS", weightBS);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("weightBSpi0", weightBSpi0);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("weightBSeta", weightBSeta);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("uniqueComboID", uniqueComboID);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("protonID", locProtonTrackID);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("beamID", locBeamID);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("chiSq", locChiSqKinFit);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("chiSq_gpi0pi0", locChiSqKinFit_gpi0pi0);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("unusedEnergy", locUnusedEnergy);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0p", locPi0Proton_Kin);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Metap", locEtaProton_Kin);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("whichSignalRegion", whichSignalRegion);
        dFlatTreeInterface->Fill_Fundamental<bool>("isTruePi0Eta", isTruePi0Eta);
        dFlatTreeInterface->Fill_Fundamental<bool>("insideEllipse", pinsideEllipse);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta_thrown", Mpi0eta_thrown);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("numSpect", numSpect);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("spectroscopicID",locProtonTrackID,0); // branch name, value, array index
        dFlatTreeInterface->Fill_Fundamental<Int_t>("spectroscopicID",locPhoton1NeutralID,1); // branch name, value, array index
        dFlatTreeInterface->Fill_Fundamental<Int_t>("spectroscopicID",locPhoton2NeutralID,2); // branch name, value, array index
        dFlatTreeInterface->Fill_Fundamental<Int_t>("spectroscopicID",locPhoton3NeutralID,3); // branch name, value, array index
        dFlatTreeInterface->Fill_Fundamental<Int_t>("spectroscopicID",locPhoton4NeutralID,4); // branch name, value, array index
        dFlatTreeInterface->Fill_Fundamental<bool>("beamPhotonMatchToThrown", beamPhotonMatchToThrown);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("rfTime", locDeltaTRF);

	if(showOutput){ cout << "Filled some fundamental branches" << endl; } 

        ++uniqueComboID;
        if (is_pi0eta){
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0g3", massGammaPi0[0]);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0g4", massGammaPi0[1]);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("logMpi0g3", log(massGammaPi0[0]));
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("logMpi0g4", log(massGammaPi0[1]));
        	dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0g3",isNotRepeated_pi0g3);
        	dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0g4",isNotRepeated_pi0g4);

        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0", locPi0Mass_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Meta", locEtaMass_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta", locPi0Eta_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0_meas", locPi0Mass);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Meta_meas", locEtaMass);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta_meas", locPi0Eta);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_teta_meas", mandelstam_teta_meas);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_teta", mandelstam_teta);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tpi0_meas", mandelstam_tpi0_meas);
        	dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tpi0", mandelstam_tpi0);
        }
        else{
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0", locPi0Mass_Kin);
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0", locEtaMass_Kin);
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0pi0", locPi0Eta_Kin);
        }
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_t", mandelstam_t);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tp", mandelstam_tp);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("ptGT1", ptGT1);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("ptLT05", ptLT05);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("ptGT05LT1", ptGT05LT1);
	if(showOutput){ cout << "Filled some more fundamental branches" << endl; } 
	
	// some variables to track uniqueness
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0",isNotRepeated_pi0);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_eta",isNotRepeated_eta);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0eta",isNotRepeated_pi0eta);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_eta_pi0eta",isNotRepeated_eta_pi0eta);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0_pi0eta",isNotRepeated_pi0_pi0eta);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("baseAsymCut",baseAsymCut);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("baseAsymCut_mDelta",baseAsymCut_mDelta);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 

	// some variables to show where the pi0/eta detected
        dFlatTreeInterface->Fill_Fundamental<Float_t>("pi0DetectedIn", pi0DetectedIn); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("etaDetectedIn", etaDetectedIn); //fundamental = char, int, float, double, etc.
	
        // Introduce some angles to use in the phase space distance calculation
        dFlatTreeInterface->Fill_Fundamental<Float_t>("prodPlane_x", y.X()); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("prodPlane_y", y.Y()); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("prodPlane_z", y.Z()); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_X_cm", cosTheta_pi0eta_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_cm", cosTheta_eta_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_pi0_cm", cosTheta_pi0_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_X_cm_meas", cosTheta_pi0eta_CM_meas); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_X_cm", phi_pi0eta_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_X_lab", phi_pi0eta_lab); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Phi", Phi); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_lab", locPhi_eta); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_pi0_lab", locPhi_pi0); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_gj", cosTheta_eta_GJ); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_gj_meas", cosTheta_eta_GJ_meas); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_gj", phi_eta_GJ); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_gj_meas", phi_eta_GJ_meas); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_hel", cosTheta_eta_hel); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_hel", phi_eta_hel); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosThetaHighestEphotonIneta_gj", cosTheta_largestEinEta_GJ); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosThetaHighestEphotonInpi0_cm", cosTheta_largestEinPi0_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_pi0_g3", ph123Rest_angle_pi0_g3); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_pi0_g4", ph124Rest_angle_pi0_g4); //fundamental = char, int, float, double, etc.

        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_g12", ph123Rest_angle_g12); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_g13", ph123Rest_angle_g13); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_g14", ph123Rest_angle_g14); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_g23", ph123Rest_angle_g23); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_g24", ph123Rest_angle_g24); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph123Rest_angle_g34", ph123Rest_angle_g34); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_g12", ph124Rest_angle_g12); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_g13", ph124Rest_angle_g13); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_g14", ph124Rest_angle_g14); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_g23", ph124Rest_angle_g23); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_g24", ph124Rest_angle_g24); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Float_t>("ph124Rest_angle_g34", ph124Rest_angle_g34); //fundamental = char, int, float, double, etc.

        dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_x",vanHove_x);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_y",vanHove_y);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_omega",omega);
        // If we were to do Q-Values for the pi0pi0 system we will probably do the same thing and use one of the pi0 as the discriminator varible and check against the other
        dFlatTreeInterface->Fill_Fundamental<Float_t>("pi0_energy", locPi0E_Kin );
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta1",photonThetas[0]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta2",photonThetas[1]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta3",photonThetas[2]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta4",photonThetas[3]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE1",photonEnergies[0]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE2",photonEnergies[1]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE3",photonEnergies[2]);	
        dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE4",photonEnergies[3]);	

	if(showOutput){ cout << "Filled fundamental branches" <<endl; }
        Fill_FlatTree(); //for the active combo
    } // end of combo loop
    ++eventIdx;
    if(showOutput){cout << "\n\n **************** Finishing the combo loop ***************\n**********************************************************\n" << endl;}

    //FILL HISTOGRAMS: Num combos / events surviving actions
    Fill_NumCombosSurvivedHists();
    if(showOutput){ cout << "Fillied NumComboSurvivedHists" << endl; } 

    /******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
    /*
    //Thrown beam: just use directly
    if(dThrownBeam != NULL)
    double locEnergy = dThrownBeam->Get_P4().E();

    //Loop over throwns
    for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
    {
    //Set branch array indices corresponding to this particle
    dThrownWrapper->Set_ArrayIndex(loc_i);

    //Do stuff with the wrapper here ...
    }
    */
    /****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
    /*
    //Loop over beam particles (note, only those appearing in combos are present)
    for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
    {
    //Set branch array indices corresponding to this particle
    dBeamWrapper->Set_ArrayIndex(loc_i);

    //Do stuff with the wrapper here ...
    }

    //Loop over charged track hypotheses
    for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
    {
    //Set branch array indices corresponding to this particle
    dChargedHypoWrapper->Set_ArrayIndex(loc_i);

    //Do stuff with the wrapper here ...
    }

    //Loop over neutral particle hypotheses
    for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
    {
    //Set branch array indices corresponding to this particle
    dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

    //Do stuff with the wrapper here ...
    }
    */

    /************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
    /*
       Bool_t locIsEventCut = true;
       for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
    //Set branch array indices for combo and all combo particles
    dComboWrapper->Set_ComboIndex(loc_i);
    // Is used to indicate when combos have been cut
    if(dComboWrapper->Get_IsComboCut())
    continue;
    locIsEventCut = false; // At least one combo succeeded
    break;
    }
    if(!locIsEventCut && dOutputTreeFileName != "")
    Fill_OutputTree();
    */


    Bool_t locIsEventCut = true;
    if(showOutput) { cout << "Looping over through the comobs to check passed or not" << endl; } 
    for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
        //Set branch array indices for combo and all combo particles
        dComboWrapper->Set_ComboIndex(loc_i);
        // Is used to indicate when combos have been cut
        if(dComboWrapper->Get_IsComboCut()){
    	    if(showOutput) { cout << "Combo did not pass cuts" << endl; } 
            continue;
	}
        locIsEventCut = false; // At least one combo succeeded                                                     
        if(showOutput) { cout << "Combo passed cuts!"  << endl; } 
        break;
    }
    if(!locIsEventCut && dOutputTreeFileName != ""){ 
	    if (showOutput) {cout<<"Filling output tree" << endl; }
	    Fill_OutputTree(); 
    }
    return kTRUE; // this return should close the process loop to return false as the kTrue as the output.
}// end of process loop

void DSelector_ver20::Finalize(void)
{
    if(true){cout << "Count combos: " << std::to_string(count_combos) << endl;} 
    if(true){cout << "Percent passed ShowerQuality: " << std::to_string((double)count_ShowerQuality/count_combos)<<endl; }
    if(true){cout << "Percent passed BeamE8GeVPlus: " << std::to_string((double)count_BeamE8GeVPlus/count_combos)<<endl; }
    if(true){cout << "Percent passed UnusedEnergy: " << std::to_string((double)count_UnusedEnergy/count_combos)<<endl;}
    if(true){cout << "Percent passed CLKinFit: " << std::to_string((double)count_ChiSq/count_combos)<<endl;}
    if(true){cout << "Percent passed DeltaTRF: " << std::to_string((double)count_DeltaTRF/count_combos)<<endl;}
    if(true){cout << "Percent passed dij3pass: " << std::to_string((double)count_dij3pass/count_combos)<<endl;}
    if(true){cout << "Percent passed PhotonE: " << std::to_string((double)count_PhotonE/count_combos)<<endl;}
    if(true){cout << "Percent passed PhotonTheta: " << std::to_string((double)count_PhotonTheta/count_combos)<<endl;}
    if(true){cout << "Percent passed MagP3Proton: " << std::to_string((double)count_MagP3Proton/count_combos)<<endl;}
    if(true){cout << "Percent passed zCutmin: " << std::to_string((double)count_zCutmin/count_combos)<<endl;}
    if(true){cout << "Percent passed RProton: " << std::to_string((double)count_RProton/count_combos)<<endl;}
    if(true){cout << "Percent passed MissingMassSquared: " << std::to_string((double)count_MissingMassSquared/count_combos)<<endl;}
    if(true){cout << "Percent passed dEdxCDCProton: " << std::to_string((double)count_dEdxCDCProton/count_combos)<<endl;}
    if(true){cout << "Percent passed insideEllipse: " << std::to_string((double)count_insideEllipse/count_combos)<<endl;}
    if(true){cout << "Percent passed allGeneralCutsPassed: " << std::to_string((double)count_allGeneralCutsPassed/count_combos)<<endl;}
    if(true){cout << "Percent passed mMPi0P14: " << std::to_string((double)count_MPi0P14/count_combos)<<endl;}
    if(true){cout << "Percent pass resolutionTest: " << std::to_string((double)count_seanResTest/count_combos) << endl; }
    DSelector::Finalize(); //Saves results to the output file
};


#include "DSelector_etapi.h"

// "4#gammap[#pi^{0},#eta]"
string topologyString="4#gammap[#pi^{0},#eta]";

void DSelector_etapi::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = ""; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = "output_flat.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "kin"; //if blank, default name will be chosen
	dSaveDefaultFlatBranches = false; // False: don't save default branches, reduce disk footprint.
	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.
        
	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

        if (dFlatTreeFileName!=""){
            // Fundamental = char, int, float, double, etc.
	    // AmpTools tree output - step 2
	    // Creating new branches in the flat tree
 	    SetupAmpTools_FlatTree(); // sets most of the branches necesary for AmpTools PWA
 	    dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Target_Mass"); 
 	    dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("PID_FinalState","NumFinalState");
 	    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("BeamAngle");
            // Photon Related 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta1");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta2");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta3");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonTheta4");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE1");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE2");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE3");	
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("photonE4");	
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pPhotonE");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pPhotonTheta");
            // Proton Related
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_momentum");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_z");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_R");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("proton_dEdxCDC");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pMagP3Proton");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pzCutmin");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pRProton");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pdEdxCDCProton");
            // Exclusivity Related
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("DOFKinFit"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("chiSq"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("unusedEnergy"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mmsq");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pMissingMassSquared");
            // Kinematics Related
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_13");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_24");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_23");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mismatchPairMass_14");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("omegaCut");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0p");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Metap");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0g3"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0g4"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Meta"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Ebeam");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tp"); 
 	    dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_t"); 
	    dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_teta");	
	    dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tpi0");	
            ////// Angles related
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Phi"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_X_cm"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_X_cm"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_X_lab"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_lab"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_pi0_lab"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_gj"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_gj"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_hel"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_hel"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_omega");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pVH_pi0p");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pVH_etap");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("pVH");
            // Weighting Related
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("AccWeight"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightASBS"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightBS"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightBSpi0"); 
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("weightBSeta"); 
            // Thrown quantities
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta_thrown");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_t_thrown");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Ebeam_thrown");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isCorrectCombo");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isCorrectBeam");
            dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isCorrectSpect");
        }
	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

//	// EXAMPLE: Create deque for histogramming particle masses:
//	// // For histogramming the phi mass in phi -> K+ K-
//	// // Be sure to change this and dAnalyzeCutActions to match reaction
//	std::deque<Particle_t> MyPhi;
//	MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);
//
//	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
//	//false/true below: use measured/kinfit data
//
//	//PID
//	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
//	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
//	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));
//
//	//PIDFOM (for charged tracks)
//	dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
//	//dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
//	//dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));
//
//	//MASSES
//	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
//	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));
//
//	//KINFIT RESULTS
//	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));
//
//	//CUT MISSING MASS
//	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));
//
//	//CUT ON SHOWER QUALITY
//	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));
//
//	//BEAM ENERGY
//	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
//	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999
//
//	//KINEMATICS
//	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
//
//	// ANALYZE CUT ACTIONS
//	// // Change MyPhi to match reaction
//	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );
//
//	//INITIALIZE ACTIONS
//	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
//	Initialize_Actions();
//	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	// Best practices is to include the bin width in the axis labels
	dHist_BeamEnergy = new TH1F("BeamEnergy", ";Beam Energy (GeV);Entries / 0.1 GeV", 90, 3.0, 12.0);
	dHist_Metapi_tot = new TH1F("Metapi_tot",";M(4#gamma) (GeV);Entries / 0.02 GeV",100,0.5,2.5);
	dHist_Metapi_sig = new TH1F("Metapi_sig",";M(4#gamma) (GeV);Entries / 0.02 GeV",100,0.5,2.5);
	dHist_Metapi_bkg = new TH1F("Metapi_bkg",";M(4#gamma) (GeV);Entries / 0.02 GeV",100,0.5,2.5);
	dHist_Mpi0p = new TH1F("Mpi0p",";M(#gamma_{1}#gamma_{2}p) (GeV);Entries / 0.04 GeV",100,0,4);
	dHist_Metap = new TH1F("Metap",";M(#gamma_{3}#gamma_{4}p) (GeV);Entries / 0.04 GeV",100,0,4);
	dHist_Meta = new TH1F("Meta",";M(#gamma_{3}#gamma_{4}) (GeV);Entries / 0.05 GeV",100,0.3,0.8);
	dHist_Mpi0 = new TH1F("Mpi0",";M(#gamma_{1}#gamma_{2}) (GeV);Entries / 0.0012 GeV",100,0.08,0.2);
	dHist_t = new TH1F("mandelstam_t",";-t GeV^{2};Entries / 0.05 GeV^{2}",60,0,3);
	dHist_rf = new TH1F("rftime",";RF (ns);Entries / 0.4 ns",120,-24,24);
	dHist_mmsq = new TH1F("mmsq",";MMsq;Entries / 0.002 GeV",100,-0.1,0.1);
	dHist_chiSq = new TH1F("chiSq",";#chi^{2};Entries / 2",50,0,100);  
	dHist_photonThetaPi0 = new TH1F("photonThetaPi0",";#theta of #gamma_{1}(#gamma_{2}) GeV;Entries / 0.05 GeV^{2} ",80,0,40); 
	dHist_photonThetaEta = new TH1F("photonThetaEta",";#theta of #gamma_{3}(#gamma_{4}) GeV;Entries / 0.05 GeV^{2} ",80,0,40); 
	dHist_dEdx_momentum = new TH2F("dEdx_momentum",";Proton Momentum Entries / 0.04 GeV/c;dEdx_{CDC} Entries / 3E-7 GeV/cm",100,0,4,100,0,0.00003);
	dHist_protonZ = new TH1F("proton_z",";Proton z (cm);Entries / 1 cm", 50, 40,90); 
	dHist_cosThetaHelVsMetapi0 = new TH2F("cosThetaHelVsMetapi0",";M(#eta#pi^{0}) Entries / 0.02 GeV;cos_{hel}(#theta) #eta Entries / 0.04",100,0.5,2.5,50,-1,1);
	dHist_cosThetaGJVsMetapi0 = new TH2F("cosThetaGJVsMetapi0",";M(#eta#pi^{0}) Entries / 0.02 GeV;cos_{GJ}(#theta) #eta Entries / 0.04",100,0.5,2.5,50,-1,1);
	dHist_combosRemaining = new TH1F("combosRemaining",";# Combos / Event Passed Selections",7,0,7);

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
	// dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");

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

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_etapi::Process(Long64_t locEntry)
{
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

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	        hasPolarizationAngle = dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, locPolarizationAngle);
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
//	Reset_Actions_NewEvent();
//	dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/

	TString locThrownTopology = Get_ThrownTopologyString();
        topologies.insert(locThrownTopology);

	//Thrown beam: just use directly
        float locBeamE_thrown=0;
        float locMetapi0_thrown=0;
        float locT_thrown=0; 
        TLorentzVector locProtonP4_thrown;
        TLorentzVector locEtaP4_thrown;
        TLorentzVector locPi0P4_thrown;

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
                else if (locPID==7)
                    locPi0P4_thrown=locThrownP4_thrown;
                else if (locPID==17)
                    locEtaP4_thrown=locThrownP4_thrown;
	}

	locMetapi0_thrown=(locPi0P4_thrown+locEtaP4_thrown).M();
	locT_thrown=-(dTargetP4-locProtonP4_thrown).M2();		
	//bool bMetapi0_thrown = (locMetapi0_thrown>1.04)*(locMetapi0_thrown<1.56);
	bool bmandelstamt_thrown=(locT_thrown<1.0)*(locT_thrown>0.1); 
        bool bBeamE_thrown = (locBeamE_thrown<8.8)*(locBeamE_thrown>8.2);
        bool bMpi0eta_thrown = (locMetapi0_thrown<1.80)*(locMetapi0_thrown>0.8);
        bool bTopology = locThrownTopology==topologyString; 
        cout << locThrownTopology << endl;
        bool selection_thrown=bTopology*bBeamE_thrown;//*bmandelstamt_thrown*bMpi0eta_thrown;
        if (dIsMC*!selection_thrown)
            return kTRUE;

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	int combos_remaining=0;
	float combo_weight=0;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()){ // Is false when tree originally created
			continue;} // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		//Step 2
		Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();
		Int_t locPhoton4NeutralID = dPhoton4Wrapper->Get_NeutralID();

                // We have access to the truth information so we can construct a scheme to select the true combination
                //    in the simulated recon tree there is a branch "isTrueCombo" that exists but for some reason
                //    it is all set to zero for our trees. We can do this manually though by checking PIDs
                // Set the default values to true so that if we do not have thrown information we do not need to do any matching
                // USAGE: Filling the output tree with the information allows us to select on the truth
                //      1. Accidental subtraction statistically selects the true beam photon = isCorrectBeam
                //      2. Mass sideband subtraction statistically selects the true pi0/eta combination = isCorrectSpect
                //      The idea is to apply all selections then make a plot. Compare subtraction scheme to correct particles
                bool isCorrectCombo=true; 
                bool isCorrectBeam=true;
                bool isCorrectSpect=true;
	        vector<Int_t> thrownPIDs;
	        vector<Int_t> parentIDs;
                vector<Int_t> matchedParentPIDs;
                cout << endl;
	        if (Get_NumThrown()!=0){
	            for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	            {	
	                dThrownWrapper->Set_ArrayIndex( loc_i );
	                //cout << "thrown PID: " << dThrownWrapper->Get_PID() << " with parent thrown index: " <<  dThrownWrapper->Get_ParentIndex() << endl;
	                thrownPIDs.push_back(dThrownWrapper->Get_PID());
	                parentIDs.push_back(dThrownWrapper->Get_ParentIndex());
	            }

                    // Obtain all the thrown IDs for the photons
                    vector<Int_t> thrownID_phs = {
                        dPhoton1Wrapper->Get_ThrownIndex(), 
                        dPhoton2Wrapper->Get_ThrownIndex(), 
                        dPhoton3Wrapper->Get_ThrownIndex(), 
                        dPhoton4Wrapper->Get_ThrownIndex(),
                    };

                    for (auto thrownID: thrownID_phs){
                        if (thrownID != -1){ // if -1 then not matched to a thrown particle
                            //cout << "ph parent " << parentIDs[thrownID] << " has PID " << thrownPIDs[parentIDs[thrownID]] << endl;
                            matchedParentPIDs.push_back(thrownPIDs[parentIDs[thrownID]]);
                        }
                        else{ // if any of  the photons did not match to a thrown particle we do not have the correct combo clearly
                            isCorrectSpect=false;
                            //cout << "ph has no parent" << endl;
                            matchedParentPIDs.push_back(-1);
                        }
                    }
                    // photons 1,2 should pair to a pi0 (Geant PID=7) and photons 3,4 should pair to an eta (17)
                    //    proton should have a proton PID=14
                    if ((matchedParentPIDs[0]==7)*
                        (matchedParentPIDs[1]==7)*
                        (matchedParentPIDs[2]==17)*
                        (matchedParentPIDs[3]==17)*
                        (thrownPIDs[dProtonWrapper->Get_ThrownIndex()]==14)){
                        isCorrectSpect=true;}
                    else{
                        isCorrectSpect=false;}

                    // Checking to see if the beam photon matches the thrown by comparing the energies
                    //cout << "Thrown:Combo Beam E " << locBeamE_thrown << ":" << (dComboBeamWrapper->Get_P4()).E() << endl;
                    if ( abs(locBeamE_thrown-(dComboBeamWrapper->Get_P4()).E())<0.0001 )
                        isCorrectBeam=true;
                    else
                        isCorrectBeam=false;

                    // Finally, the true combo is when we have the true spectroscopic combination and the true beam photon
                    isCorrectCombo=isCorrectSpect*isCorrectBeam;
                    cout << "correct beam/spect/combo: " << isCorrectBeam << "/" << isCorrectSpect << "/" << isCorrectCombo <<endl;
	        }


		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		TLorentzVector locProtonX4 = dProtonWrapper->Get_X4();
		//Step 1
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		//Step 2
		TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4();
		TLorentzVector locPhoton4P4 = dPhoton4Wrapper->Get_P4();
		//Construct Intermediate Resonances
		TLorentzVector locEtaP4=locPhoton3P4+locPhoton4P4;
		TLorentzVector locPi0P4=locPhoton1P4+locPhoton2P4;

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
		//Step 2
		TLorentzVector locPhoton3P4_Measured = dPhoton3Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton4P4_Measured = dPhoton4Wrapper->Get_P4_Measured();

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		// 0 for in-time events, non-zero integer for out-of-time photons
		Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); 
		Int_t locNumOutOfTimeBunchesInTree = 4; //YOU need to specify this number
		//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times 
		//	this number for left + right bunches) 
		Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		// Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); 
		//locAccidentalScalingFactor=1; // Lawrence - testing something
		Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); 
		// Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; 
		if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
			dComboWrapper->Set_IsComboCut(true); 
			continue; 
		} 

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE
		//// CONSTRUCT SIDEBAND WEIGHTS
		// A hidden step here that required fitting M(pi0) and M(eta) to extract 
		// 	the peak and widths and associated weightings. Take these number
		// 	as given for now
		float Mpi0=locPi0P4.M();
		float Meta=locEtaP4.M();
		float pi0Mean=0.135881;
		float etaMean=0.548625;
		float pi0Std=0.0076;
		float etaStd=0.0191;
                float pi0Sig=3; // number of sigmas to the left and right of the mean
                float etaSig=3;
                float pi0Skip=1; // number of sigmas to skip right after signal region
                float etaSkip=1;
                float pi0SB=1.5; // number of sigmas that define the sideband region which is right after the skip region
                float etaSB=2;
		float pi0_sbweight;
		float eta_sbweight;
		// The signal regions are both +/- 3 sigmas around the peak the left and right sidebands 
		// 	which are some N sigmas wide with some M sigma skip region included 
		// 	between the signal and sideband regions. The weight = the ratio the lengths
		// 	spanned by the signal to that of the sideband times -1.
		if (Mpi0 > pi0Mean-pi0Sig*pi0Std && Mpi0 < pi0Mean+pi0Sig*pi0Std){ pi0_sbweight=1; }
		else if (Mpi0 > pi0Mean+(pi0Sig+pi0Skip)*pi0Std && Mpi0 < pi0Mean+(pi0Sig+pi0Skip+pi0SB)*pi0Std){ pi0_sbweight=-1*pi0Sig/pi0SB; }
		else if (Mpi0 > pi0Mean-(pi0Sig+pi0Skip+pi0SB)*pi0Std && Mpi0 < pi0Mean-(pi0Sig+pi0Skip)*pi0Std){ pi0_sbweight=-1*pi0Sig/pi0SB; }
		else { pi0_sbweight=0; }
		if (Meta > etaMean-etaSig*etaStd && Meta < etaMean+etaSig*etaStd){ eta_sbweight=1; }
		else if (Meta > etaMean+(etaSig+etaSkip)*etaStd && Meta < etaMean+(etaSig+etaSkip+etaSB)*etaStd){ eta_sbweight=-1*etaSig/etaSB; }
		else if (Meta > etaMean-(etaSig+etaSkip+etaSB)*etaStd && Meta < etaMean-(etaSig+etaSkip)*etaStd){ eta_sbweight=-1*etaSig/etaSB; }
		else { eta_sbweight=0; }
		float sbweight=pi0_sbweight*eta_sbweight;
		float weight=sbweight*locHistAccidWeightFactor;
                //weight=1;
		// Reject combinations with zero weights. Zero weights take up space and do nothing. 
		// 	Worse, it might cause the amptools unbinned likelihood fits to break
		bool bWeight=(weight==0) ? false : true; 

		// AMPTOOLS REQUIRES 4 TREES (data, bkgnd, accmc, genmc)
		//    data = selected trees of DATA where signal region has been selected, all weights = 1
		//    bkgnd = selected trees of DATA where sidebands have been selected, all weights = -weight
		//    accmc = selected trees of ACCEPTANCE MC where signal+sidebands have been selected, all weights = weight
		//    genmc = thrown trees created during simulation process
		bool bSignalRegion;
		float branchWeight;
                int choice=3;
		//---------CHOICE 1 FOR "data" RUN OVER SIGNAL/DATA-------------
                if (choice==1){
		    bSignalRegion=(pi0_sbweight==1)*(eta_sbweight==1)*(locHistAccidWeightFactor==1); // Keep combos ONLY in the signal region
		    branchWeight=1;}
		//---------CHOICE 2 FOR "bkgnd" RUN OVER SIGNAL/DATA-------------
                if (choice==2){
		    bSignalRegion=!((pi0_sbweight==1)*(eta_sbweight==1)*(locHistAccidWeightFactor==1)); // Keep combos ONLY in the sideband region
		    branchWeight=-weight;}
		//---------CHOICE 3 FOR "accmc" RUN OVER FLAT MC-------------
                if (choice==3){
		    bSignalRegion=true; // Keep combos that exist in the signal AND sideband region
		    branchWeight=weight;}
		//----------------------

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
//		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
//		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
//			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
			//dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); // Alternate version with accidental subtraction

			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();


		//////////////////////////////////////////////////////////////////
		//  DEFINING SELECTIONS AND INTERESTING VARIABLES
		//////////////////////////////////////////////////////////////////
		//// 1. NEUTRAL SHOWER RELATED (selecting good photons)
		// Low energy photons are more likely to be spurious, require a minimum E
		bool bPhotonE=(locPhoton1P4.E()>0.1)*(locPhoton2P4.E()>0.1)*(locPhoton3P4.E()>0.1)*(locPhoton4P4.E()>0.1); 
		// Working in degrees instead of radians, we remove photons near the beamline (<~2.5) and near the BCAL/FCAL transition (<~11.9, >~10.3)
		float radToDeg=180/3.14159;
		bool bPhotonTheta=
			((locPhoton1P4.Theta()*radToDeg>=2.5 && locPhoton1P4.Theta()*radToDeg<=10.3) || locPhoton1P4.Theta()*radToDeg>=11.9)*
			((locPhoton2P4.Theta()*radToDeg>=2.5 && locPhoton2P4.Theta()*radToDeg<=10.3) || locPhoton2P4.Theta()*radToDeg>=11.9)*
			((locPhoton3P4.Theta()*radToDeg>=2.5 && locPhoton3P4.Theta()*radToDeg<=10.3) || locPhoton3P4.Theta()*radToDeg>=11.9)*
			((locPhoton4P4.Theta()*radToDeg>=2.5 && locPhoton4P4.Theta()*radToDeg<=10.3) || locPhoton4P4.Theta()*radToDeg>=11.9);

		//// 2. CHARGED TRACK RELATED (selecting good protons)
		// protons need some momentum be reconstructed properly
		bool bProtonMomentum=locProtonP4.Vect().Mag()>0.3; 
		// separate proton/pi+ based on energy loss in CDC
		bool bProton_dEdx=dProtonWrapper->Get_dEdx_CDC()>=TMath::Power(10,-6)*(0.9+TMath::Exp(3.0-3.5*(locProtonP4.Vect().Mag()+0.05)/.93827)); 
		// require proton to come from ~ the target region [52,78]centimeters
		bool bProtonZ = 52 <= locProtonX4.Z() && locProtonX4.Z() <= 78;

		//// 3. EXCLUSIVITY RELATED (ensure we select exclusive gp->4gp reaction)
		// Kinematic fit for this tree only attempts to quantify how well conservation of 4-momentum is maintained
		// 	4 NDF in this fit where 13.277 corresponds to a p=0.01
		bool bChiSq=dComboWrapper->Get_ChiSq_KinFit("")<13.277;
		// Unused energy is the sum of unused neutral shower energy not used by this current combo
		// 	require no unused energy meaning the event only has 4 neutral shower hypotheses to limit final state combinatorics
		bool bUnusedEnergy=dComboWrapper->Get_Energy_UnusedShowers()<0.05;
		// No missing particles are expected = no missing mass
		bool bMMsq=abs(locMissingMassSquared)<0.05;
		
		//// 4. KINEMATICS RELATED (extra selections related to kinematics)
		float Mpi0p=(locPi0P4+locProtonP4).M();
		float Metap=(locEtaP4+locProtonP4).M();
		float Metapi0=(locPi0P4+locEtaP4).M();
		float mandelstam_t=-(dTargetP4-locProtonP4).M2();		
                float mandelstam_teta = -(locBeamP4-locEtaP4).M2();
                float mandelstam_tpi0 = -(locBeamP4-locPi0P4).M2();
		// Select on coherent peak for region of high polarization. The AMPTOOLS fit using Zlm amplitudes will use the polarization
		// 	for extra separation power (will tell us something about the production mechanism)
		bool bBeamEnergy=(locBeamP4.E()>8.2)*(locBeamP4.E()<8.8); 
		bool bMetapi0 = (Metapi0>1.04)*(Metapi0<1.56); // Select the a2(1320) mass region
		// Meson production occurs with small-t whereas baryon production occurs with large-t. This analysis cares about mesons
		bool bmandelstamt=(mandelstam_t<0.3)*(mandelstam_t>0.1); 
		// There are a few baryon resonances, the Delta+(1232) being the largest. We can reject it 
		bool bMpi0p=Mpi0p>1.4; 
		// We are interested in the eta+pi0 system where we can see multiple resonances
		// 	in the mass spectra. The spin of the resonance has a direct influence
		// 	on the decay angular distributions of the daughter particles (eta/pi)
		// 	The helicity frame and Gottfried-Jackson frame are both in the rest 
		// 	frame of the etapi system and requires a redefinition of the axes. 
		// 	For this analysis, the helicity frame will be used as the amplitudes 
		// 	that you will use are based on this frame
		// First we need to boost to center of mass frame
		TLorentzRotation cmRestBoost( -(locBeamP4+dTargetP4).BoostVector() );
		TLorentzVector pi0_cm = cmRestBoost * locPi0P4; 
		TLorentzVector eta_cm = cmRestBoost * locEtaP4; 
		TLorentzVector beam_cm = cmRestBoost * locBeamP4;
		TLorentzVector recoil_cm = cmRestBoost * locProtonP4;
   		TLorentzVector resonance = pi0_cm + eta_cm;
		// We boost again, now to the resonances rest frame
   		TLorentzRotation resRestBoost( -resonance.BoostVector() );
   		TLorentzVector beam_res   = resRestBoost * beam_cm;
   		TLorentzVector recoil_res = resRestBoost * recoil_cm;
   		TLorentzVector eta_res = resRestBoost * eta_cm;
		//// Redefinition of the axes
   		TVector3 z = -1. * recoil_res.Vect().Unit(); // CALCULATING FOR Helicity frame
   		TVector3 y = (beam_cm.Vect().Unit().Cross(-recoil_cm.Vect().Unit())).Unit();
   		TVector3 x = y.Cross(z);
   		TVector3 angles( (eta_res.Vect()).Dot(x),
   		      (eta_res.Vect()).Dot(y),
   		      (eta_res.Vect()).Dot(z) );
   		float cosTheta_hel = angles.CosTheta();
   		float phi_hel = angles.Phi();
   		z = beam_res.Vect().Unit(); // CALCULATING FOR Gottfried-Jackson frame
   		x = y.Cross(z);
   		angles.SetXYZ( (eta_res.Vect()).Dot(x),
   		      (eta_res.Vect()).Dot(y),
   		      (eta_res.Vect()).Dot(z) );
   		float cosTheta_gj = angles.CosTheta();
   		float phi_gj = angles.Phi();
                TVector3 eps(TMath::Cos(locPolarizationAngle*TMath::DegToRad()), TMath::Sin(locPolarizationAngle*TMath::DegToRad()), 0.0); // beam polarization vector
                float Phi = TMath::ATan2(y.Dot(eps), beam_cm.Vect().Unit().Dot(eps.Cross(y)))*radToDeg;
                float mandelstam_t0 = -(TMath::Power(-(locPi0P4+locEtaP4).M2()/(2*(locBeamP4+dTargetP4).M()),2)
                                        -TMath::Power(beam_cm.Vect().Mag()-(pi0_cm+eta_cm).Vect().Mag(),2));
                float mandelstam_tp = mandelstam_t-mandelstam_t0;
                std::tuple<double, double> vh = dAnalysisUtilities.Calc_vanHoveCoord(recoil_cm,pi0_cm,eta_cm);
                float q = get<0>(vh);
                float omega = get<1>(vh);
                bool bVH_pi0p = -29.61317407*atan(-0.9877663*(locPi0P4+locEtaP4).M()+2.77936736)+330.46008765 > omega;
                bool bVH_etap = 45.26878219*atan(-0.88242654*(locPi0P4+locEtaP4).M()+3.14340627)+193.59347205 < omega;
                bool bVH = bVH_pi0p*bVH_etap;

		// 5. With the above selections and sidebands stucture, the subtraction near threshold has problems
		// 	The backgrounds that populate the near threshold region are pi0pi0->4g and omega->3gamma
		// 	These backgrounds tend to populate the lower masses in the alternative (not {g1g2~pi0,g3g4~eta}) photon pairs
		float Mg1g3=(locPhoton1P4+locPhoton3P4).M();
		float Mg1g4=(locPhoton1P4+locPhoton4P4).M();
		float Mg2g3=(locPhoton2P4+locPhoton3P4).M();
		float Mg2g4=(locPhoton2P4+locPhoton4P4).M();
		bool bLowMassAltCombo=
		    !(((Mg1g3<0.15)*(Mg2g4<0.15)) || ((Mg1g4<0.15)*(Mg2g3<0.15)) ||
		    ((Mg1g3<0.12)*(Mg2g3<0.12)) || ((Mg1g4<0.12)*(Mg2g4<0.12)));

                // Turn off some selections related to M(4g) and t so that we can use another program to
                //      split the final flat trees up. This should lower our total run times
		bMetapi0=true; 
                bmandelstamt=true; 

		//// We can finally multiply all of our selections together to define our final selection criteria. 
		//	Since we have defined SELECTIONS we have to flip the boolean to get a CUT since the 
		//	FLAG used asks if the combo should be cut
		bool selection=bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*
				bmandelstamt*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion;
                //bool selection=(dComboWrapper->Get_ChiSq_KinFit("")<100)*(dComboWrapper->Get_Energy_UnusedShowers()<0.5); 
                //bool selection=true;

		// We generally do not want to apply a cut on a histogram we are trying to view. To this extent,
		// 	we will just apply all other selections that are not used in the current plot
		// 	i.e. if we want to plot MMSq we will apply all selections EXCEPT for bMMsq
		if (bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*
                    bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){ 
			dHist_Mpi0->Fill((locPhoton1P4+locPhoton2P4).M(),locHistAccidWeightFactor);
			dHist_Meta->Fill((locPhoton3P4+locPhoton4P4).M(),locHistAccidWeightFactor);
		}
//		if (!bWeight){ // We do not want to keep any combos with weight 0. Not just a waste of space, can cause problems during fitting
//			dComboWrapper->Set_IsComboCut(true);
//			continue;
//		}
		if(bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bBeamEnergy*bmandelstamt*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_mmsq->Fill(locMissingMassSquared,weight);}
		if(bPhotonE*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_photonThetaPi0->Fill(locPhoton1P4.Theta()*radToDeg,weight);
			dHist_photonThetaPi0->Fill(locPhoton2P4.Theta()*radToDeg,weight);
			dHist_photonThetaEta->Fill(locPhoton3P4.Theta()*radToDeg,weight);
			dHist_photonThetaEta->Fill(locPhoton4P4.Theta()*radToDeg,weight);}
		if(bPhotonE*bPhotonTheta*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_dEdx_momentum->Fill(locProtonP4.Vect().Mag(),dProtonWrapper->Get_dEdx_CDC(),weight);}
		if(bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_protonZ->Fill(locProtonX4.Z(),weight);}
		if(bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_t->Fill(mandelstam_t,weight);}
		if(bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_Mpi0p->Fill(Mpi0p,weight);}
		if(bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*bMpi0p*bLowMassAltCombo*bMetapi0*bSignalRegion){
			dHist_chiSq->Fill(dComboWrapper->Get_ChiSq_KinFit(""),weight);}
		if(bPhotonE*bPhotonTheta*bProtonMomentum*bProton_dEdx*bProtonZ*bChiSq*bUnusedEnergy*bMMsq*bBeamEnergy*bmandelstamt*bMpi0p*bLowMassAltCombo*bSignalRegion){
			dHist_Metapi_sig->Fill(Metapi0,weight);
			dHist_cosThetaHelVsMetapi0->Fill(Metapi0,cosTheta_hel,weight);
			dHist_cosThetaGJVsMetapi0->Fill(Metapi0,cosTheta_gj,weight);
			if ( (pi0_sbweight==1)*(eta_sbweight==1) )
				dHist_Metapi_tot->Fill(Metapi0,locHistAccidWeightFactor);
			else
				dHist_Metapi_bkg->Fill(Metapi0,locHistAccidWeightFactor*sbweight);
		}
		//E.g. Cut
		if(!selection){
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		++combos_remaining;
		combo_weight+=weight;
		dHist_rf->Fill(locDeltaT_RF);
		dHist_Metap->Fill(Metap,weight);

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		// RECOMMENDED: FILL ACCIDENTAL WEIGHT
		// dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight",locHistAccidWeightFactor);

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

                if (dFlatTreeFileName!=""){
                    // Photon Related 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta1", locPhoton1P4.Theta());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta2", locPhoton2P4.Theta());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta3", locPhoton3P4.Theta());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonTheta4", locPhoton4P4.Theta());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE1", locPhoton1P4.E());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE2", locPhoton2P4.E());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE3", locPhoton3P4.E());	
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("photonE4", locPhoton4P4.E());	
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pPhotonE", bPhotonE);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pPhotonTheta", bPhotonTheta);
                    // Proton Related
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_momentum", locProtonP4.Vect().Mag());
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_z", locProtonX4.Z());
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("proton_dEdxCDC", dProtonWrapper->Get_dEdx_CDC());
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pMagP3Proton", bProtonMomentum);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pzCutmin", bProtonZ);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pdEdxCDCProton", bProton_dEdx);
                    // Exclusivity Related
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("DOFKinFit", dComboWrapper->Get_NDF_KinFit("")); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("chiSq", dComboWrapper->Get_ChiSq_KinFit("")); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("unusedEnergy",dComboWrapper->Get_Energy_UnusedShowers()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mmsq",locMissingMassSquared);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pMissingMassSquared", bMMsq);
                    // Kinematics Related
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_13", Mg1g3);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_24", Mg2g4);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_23", Mg2g3);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mismatchPairMass_14", Mg1g4);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("omegaCut",bLowMassAltCombo);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0p",Mpi0p);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Metap",Metap);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0g3",(locPi0P4+locPhoton3P4).M()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0g4",(locPi0P4+locPhoton4P4).M()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0",Mpi0); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Meta",Meta); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta",Metapi0); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Ebeam",locBeamP4.E()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tp",mandelstam_tp); 
 	            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_t",mandelstam_t); 
	            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_teta",mandelstam_teta);	
	            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tpi0",mandelstam_tpi0);	
                    ////// Angles related
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Phi",Phi); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_X_cm",(pi0_cm+eta_cm).CosTheta()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_X_cm",(pi0_cm+eta_cm).Phi()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_X_lab",(locPi0P4+locEtaP4).Phi()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_lab",locEtaP4.Phi()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_pi0_lab",locPi0P4.Phi()); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_gj", cosTheta_gj); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_gj", phi_gj); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_hel",cosTheta_hel); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_hel",phi_hel); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_omega",omega);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pVH_pi0p", bVH_pi0p);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pVH_etap", bVH_etap);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("pVH", bVH);
                    // Weighting Related
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("AccWeight", locHistAccidWeightFactor); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("weightASBS", weight); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("weightBS", sbweight); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("weightBSpi0", pi0_sbweight); 
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("weightBSeta", eta_sbweight); 
                    // Thrown Variables
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta_thrown",locMetapi0_thrown);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_t_thrown",locT_thrown);
                    dFlatTreeInterface->Fill_Fundamental<Float_t>("Ebeam_thrown",locBeamE_thrown);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("isCorrectCombo",isCorrectCombo);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("isCorrectBeam",isCorrectBeam);
                    dFlatTreeInterface->Fill_Fundamental<Bool_t>("isCorrectSpect",isCorrectSpect);
		    // AmpTools tree output - step 3
		    // Filling the branches of the flat tree
		    vector<TLorentzVector> locFinalStateP4; // should be in the same order as PID_FinalState
		    locFinalStateP4.push_back(locProtonP4); 
		    locFinalStateP4.push_back(locPi0P4);
		    locFinalStateP4.push_back(locEtaP4);
		    dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", branchWeight);
		    dFlatTreeInterface->Fill_Fundamental<Int_t>("BeamAngle", locPolarizationAngle); // include so we can split on this branch later
		    dFlatTreeInterface->Fill_Fundamental<Float_t>("Target_Mass", 0.9382720); // Necesary for divideData.pl not for AmpTools itself (I think)
		    dFlatTreeInterface->Fill_Fundamental<Int_t>("PID_FinalState", 2212, 0); // proton
		    dFlatTreeInterface->Fill_Fundamental<Int_t>("PID_FinalState", 111, 1);  // Pi0
		    dFlatTreeInterface->Fill_Fundamental<Int_t>("PID_FinalState", 221, 2);  // Eta
		    FillAmpTools_FlatTree(locBeamP4, locFinalStateP4);
		    Fill_FlatTree(); //for the active combo
                }
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
//	Fill_NumCombosSurvivedHists();
//
	dHist_combosRemaining->Fill(combos_remaining);

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
        // If we want to output a reduced tree we can do it using this block of code.
        //      This is useful if we want to apply some pre-selections to filter a tree
        //      so that we can run over this filtered tree when making tighter selections
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()){
			continue;}
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
        // The second condition causes problems when dOutputTreeFileName is actually ""
	if(!locIsEventCut)// && dOutputTreeFileName != "") 
		Fill_OutputTree();

	return kTRUE;
}

void DSelector_etapi::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE
        cout << "Topologies in chain" << endl;
        for (auto topology: topologies)
            cout << topology << endl;

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}

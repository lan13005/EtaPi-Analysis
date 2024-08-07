#include "DSelector_thrown.h"

string topologyString="4#gammap[#pi^{0},#eta]";
//string topologyString="5#gammap[2#pi^{0}]";
float radToDeg=180/3.14159;

void DSelector_thrown::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "thrown.root"; //"" for none
	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning into separate files for AmpTools
	//dOutputTreeFileNameMap["pol000"] = ".root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["pol045"] = ".root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["pol090"] = ".root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["pol135"] = ".root"; //key is user-defined, value is output file name
        
	dFlatTreeFileName = "output_flat.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "kin"; //if blank, default name will be chosen
	dSaveDefaultFlatBranches = false; // False: don't save default branches, reduce disk footprint.

	dSkipNoTriggerEvents = false; // 04/03/24 required for now

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	dPreviousRunNumber = 0;

        if (dFlatTreeFileName!=""){
            // Fundamental = char, int, float, double, etc.
	    	// AmpTools tree output - step 2
	    	// Creating new branches in the flat tree
 	    	SetupAmpTools_FlatTree(); // sets most of the branches necesary for AmpTools PWA
 	    	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Target_Mass"); 
 	    	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("PID_FinalState","NumFinalState");
 	    	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("BeamAngle");
 	    	//dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Weight"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0p_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Metap_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_hel_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_gj_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_hel_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_gj_thrown"); 
 	    	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_t_thrown"); 
 	    	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Ebeam_thrown"); 
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_omega_thrown");
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_x_thrown");
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("vanHove_y_thrown");
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Phi_thrown");
        	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("pVH_thrown");
        }

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
	//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
	//By default, for each event, the data is retrieved for all branches
	//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
	//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_thrown::Process(Long64_t locEntry)
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

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/******************************************* LOOP OVER THROWN DATA ***************************************/

	TString locThrownTopology = Get_ThrownTopologyString();
        topologies.insert(locThrownTopology);

	//Thrown beam: just use directly
	double locBeamEnergyUsedForBinning = 0.0;
	if(dThrownBeam != NULL)
		locBeamEnergyUsedForBinning = dThrownBeam->Get_P4().E();
        TLorentzVector locBeamP4=dThrownBeam->Get_P4();

        TLorentzVector locProtonP4;
        TLorentzVector locEtaP4;
        TLorentzVector locPi0P4;

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
		Particle_t locPID = dThrownWrapper->Get_PID();
		TLorentzVector locThrownP4 = dThrownWrapper->Get_P4();

        if (locPID==14)
            locProtonP4=locThrownP4;
        else if (locPID==7)
            locPi0P4=locThrownP4;
        else if (locPID==17)
            locEtaP4=locThrownP4;
		//cout << "Thrown " << loc_i << ": " << locPID << ", " << locThrownP4.Px() << ", " << locThrownP4.Py() << ", " << locThrownP4.Pz() << ", " << locThrownP4.E() << endl;
	}
        //cout << locPi0P4.M() << ", " << locEtaP4.M() << ", " << locProtonP4.M() << endl;

	//OR Manually:
	//BEWARE: Do not expect the particles to be at the same array indices from one event to the next!!!!
	//Why? Because while your channel may be the same, the pions/kaons/etc. will decay differently each event.

	//BRANCHES: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat#TTree_Format:_Simulated_Data
	TClonesArray** locP4Array = dTreeInterface->Get_PointerToPointerTo_TClonesArray("Thrown__P4");
	TBranch* locPIDBranch = dTreeInterface->Get_Branch("Thrown__PID");
/*
	Particle_t locThrown1PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);
	TLorentzVector locThrown1P4 = *((TLorentzVector*)(*locP4Array)->At(0));
	cout << "Particle 1: " << locThrown1PID << ", " << locThrown1P4.Px() << ", " << locThrown1P4.Py() << ", " << locThrown1P4.Pz() << ", " << locThrown1P4.E() << endl;
	Particle_t locThrown2PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);
	TLorentzVector locThrown2P4 = *((TLorentzVector*)(*locP4Array)->At(1));
	cout << "Particle 2: " << locThrown2PID << ", " << locThrown2P4.Px() << ", " << locThrown2P4.Py() << ", " << locThrown2P4.Pz() << ", " << locThrown2P4.E() << endl;
*/

	float Metapi0=(locPi0P4+locEtaP4).M();
	float Metap=(locEtaP4+locProtonP4).M();
	float Mpi0p=(locPi0P4+locProtonP4).M();
	float mandelstam_t=-(dTargetP4-locProtonP4).M2();		
    float beam_e=locBeamP4.E();
	//bool bMetapi0 = (Metapi0>1.04)*(Metapi0<1.56);
	bool bmandelstamt=(mandelstam_t<1.0)*(mandelstam_t>0.1); 
    bool bMpi0eta = (Metapi0<1.80)*(Metapi0>0.8);
    bool bBeamE = (beam_e<8.8)*(beam_e>8.2);
    bool bTopology = locThrownTopology==topologyString; 
    bool selection=bTopology;//*bBeamE;//*bmandelstamt*bMpi0eta;

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
    std::tuple<double, double> vh = dAnalysisUtilities.Calc_vanHoveCoord(recoil_cm,pi0_cm,eta_cm);
    float q = get<0>(vh);
    float omega = get<1>(vh);
    float vanHove_x=q*cos(omega);
    float vanHove_y=q*sin(omega);
    float pVH=(float)filterOmega(omega*radToDeg,(locPi0P4+locEtaP4).M());

    if ((dFlatTreeFileName!="")*(selection)){
	    vector<TLorentzVector> locFinalStateP4; // should be in the same order as PID_FinalState
	    locFinalStateP4.push_back(locProtonP4); 
	    locFinalStateP4.push_back(locPi0P4);
	    locFinalStateP4.push_back(locEtaP4);
	    dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1);
	    dFlatTreeInterface->Fill_Fundamental<Int_t>("BeamAngle", locPolarizationAngle); // include so we can split on this branch later
	    dFlatTreeInterface->Fill_Fundamental<Float_t>("Target_Mass", 0.9382720); // Necesary for divideData.pl not for AmpTools itself (I think)
	    dFlatTreeInterface->Fill_Fundamental<Int_t>("PID_FinalState", 2212, 0); // proton
	    dFlatTreeInterface->Fill_Fundamental<Int_t>("PID_FinalState", 111, 1);  // Pi0
	    dFlatTreeInterface->Fill_Fundamental<Int_t>("PID_FinalState", 221, 2);  // Eta
 	    dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_t_thrown",mandelstam_t); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta_thrown",Metapi0); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Metap_thrown",Metap); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0p_thrown",Mpi0p); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_hel_thrown",cosTheta_hel); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_gj_thrown",cosTheta_gj); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_hel_thrown",phi_hel); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_gj_thrown",phi_gj); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Ebeam_thrown",beam_e); 
        dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_omega_thrown",omega*radToDeg);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_x_thrown",vanHove_x);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("vanHove_y_thrown",vanHove_y);
        dFlatTreeInterface->Fill_Fundamental<Float_t>("Phi_thrown",Phi);
	    FillAmpTools_FlatTree(locBeamP4, locFinalStateP4);
	    Fill_FlatTree(); //for the active combo
    }

	/******************************************* BIN THROWN DATA INTO SEPARATE TREES FOR AMPTOOLS ***************************************/

/*
	//THESE KEYS MUST BE DEFINED IN THE INIT SECTION (along with the output file names)
	if((locBeamEnergyUsedForBinning >= 8.0) && (locBeamEnergyUsedForBinning < 9.0))
		Fill_OutputTree("Bin1"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 9.0) && (locBeamEnergyUsedForBinning < 10.0))
		Fill_OutputTree("Bin2"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 10.0) && (locBeamEnergyUsedForBinning < 11.0))
		Fill_OutputTree("Bin3"); //your user-defined key
*/

	return kTRUE;
}

void DSelector_thrown::Finalize(void)
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

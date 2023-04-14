#include "DSelector_thrown.h"
#include "TRandom.h"

string topologyString="4#gammap"; // Dont think this is actually necessary. Its always 4g

string polarization="degALL";
string tag="_flat_2017_gen";

void DSelector_thrown::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = ""; //"" for none
	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning into separate files for AmpTools

	dFlatTreeFileName = "degALL_flat_2017_gen_flat_DSelector.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "kin"; //if blank, default name will be chosen
	dSaveDefaultFlatBranches = false; // False: don't save default branches, reduce disk footprint.

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	dPreviousRunNumber = 0;

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

        if (dFlatTreeFileName!=""){
            dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("BeamAngle"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tp"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_t"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_teta"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("mandelstam_tpi0"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0eta"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_gj"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_gj"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("cosTheta_eta_hel"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("phi_eta_hel"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Ebeam"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Metap");
            dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Mpi0p");
            dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("p_p4_kin");
            dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_p4_kin");
        }
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
	//
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//
	//

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	cout << "RunNumber = " << locRunNumber << endl;
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
        	hasPolarizationAngle = dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, locPolarizationAngle);
		dPreviousRunNumber = locRunNumber;
        	cout << "Getting beam polarization and filling used runs" << endl;
	}

    	keepPolarization=false;
    	if (polarization=="degALL"){ keepPolarization=true; }
    	else if (locPolarizationAngle==0 && polarization=="deg000") { keepPolarization=true; }
    	else if (locPolarizationAngle==45 && polarization=="deg045") { keepPolarization=true; }
    	else if (locPolarizationAngle==90 && polarization=="deg090") { keepPolarization=true; }
    	else if (locPolarizationAngle==135 && polarization=="deg135") { keepPolarization=true; }
    	else if (!hasPolarizationAngle && polarization=="degAMO") { keepPolarization=true; }
    	else {  cout << "FUDGE! THERE IS AN UNEXPECTED OUTCOME FROM THE POLARIZATION CHECK" << endl; } 

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/
	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/******************************************* LOOP OVER THROWN DATA ***************************************/

	//Thrown beam: just use directly
	bool beamEisFinite = 0;
	double locBeamEnergyUsedForBinning = 0.0;
	if(dThrownBeam != NULL)
		locBeamEnergyUsedForBinning = dThrownBeam->Get_P4().E();
	if(TMath::Finite(locBeamEnergyUsedForBinning)==1){
		++beamEfin;
		beamEisFinite = 1;	
	}
	else { ++beamEinf; }	

        // Checking thrown topology
	TString locThrownTopology = Get_ThrownTopologyString();
        if (topologyString!=""){
            if (locThrownTopology != topologyString.c_str() ){
                cout << "incorrect locThrownTopology = " << locThrownTopology << endl;
                return kTRUE;
            }
            else {
                cout << "correct locThrownTopology = " << locThrownTopology << endl;
            }
        }

        double cosTheta_eta_GJ; 
        double phi_eta_GJ; 
        double cosTheta_eta_hel;
        double phi_eta_hel;
	double locPi0EtaMass;
	double mandelstam_t;

	std::vector<TLorentzVector> allP4;
	std::vector<int> parentArray;
	std::vector<int> pids;
	TLorentzVector locEtaP4;
	TLorentzVector locPi0P4;
	TLorentzVector locProtonP4;
	TLorentzVector locTargetP4 = {0,0,0,0.938};
	TLorentzVector locBeamP4=dThrownBeam->Get_P4();
	//Loop over throwns
	// we use is_in as a way to select the events we want to show instead of dumping everything
	const bool is_in = true;//showOutput.find(eventIdx) != showOutput.end();
	if (is_in){
		cout << "########    EventIdx: " << eventIdx << "    #############" << endl;
	}

	int locNumThrown = Get_NumThrown();
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{	
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
		Particle_t locPID = dThrownWrapper->Get_PID();
		Int_t locParentID = dThrownWrapper->Get_ParentIndex();
                Int_t locParentPID;
                //////  Quickly swap over and get the PID  ///////
                if ( (int)locParentID != -1 ) {
                    dThrownWrapper->Set_ArrayIndex(locParentID);
                    locParentPID = dThrownWrapper->Get_PID();
                    dThrownWrapper->Set_ArrayIndex(loc_i);
                }
                else{
                    locParentPID=-1;
                }
                //////////////////////////////////////////////////
		allP4.push_back(dThrownWrapper->Get_P4());

		if (is_in){
			cout << "Thrown " << loc_i << " with (PID,ParentArrayLoc) = (" << locPID << " ," << locParentPID << ")"  << endl;
		}
		parentArray.push_back(locParentID);
		pids.push_back(locPID);
	}

        locProtonP4=allP4[0]; 
        locPi0P4=allP4[1]+allP4[2]; 
        locEtaP4=allP4[3]+allP4[4]; 
        cout << "Mpi0: " << locPi0P4.M() << endl;
        cout << "Meta: " << locEtaP4.M() << endl;

	// check for etas in idxInitial
	cout << "Calculated EBeam bins" << endl;
	
	TLorentzVector locPi0EtaP4 = locPi0P4+locEtaP4;

        TLorentzVector cm_vec = locBeamP4+locTargetP4;
        TLorentzVector locPi0EtaP4_cm = locPi0EtaP4;
        TLorentzVector locPi0P4_cm = locPi0P4;
        TLorentzVector locEtaP4_cm = locEtaP4;
        TLorentzVector locBeamP4_cm = locBeamP4;
        TLorentzVector locProtonP4_cm = locProtonP4;
        locPi0EtaP4_cm.Boost(-cm_vec.BoostVector());
        locPi0P4_cm.Boost(-cm_vec.BoostVector());
        locEtaP4_cm.Boost(-cm_vec.BoostVector());
        locBeamP4_cm.Boost(-cm_vec.BoostVector());
        locProtonP4_cm.Boost(-cm_vec.BoostVector());

	TLorentzVector locPi0P4_res = locPi0P4_cm;
	TLorentzVector locEtaP4_res = locEtaP4_cm;
	TLorentzVector locBeamP4_res = locBeamP4_cm;
	TLorentzVector locProtonP4_res = locProtonP4_cm;
	locPi0P4_res.Boost(-locPi0EtaP4_cm.BoostVector());
	locEtaP4_res.Boost(-locPi0EtaP4_cm.BoostVector());
	locBeamP4_res.Boost(-locPi0EtaP4_cm.BoostVector());
	locProtonP4_res.Boost(-locPi0EtaP4_cm.BoostVector());

	double radToDeg = 57.3;
        TVector3 locPi0P4_res_unit = locPi0P4_res.Vect().Unit();
        TVector3 locEtaP4_res_unit = locEtaP4_res.Vect().Unit();

        //////// 
        // GJ 
        ////////
        // Calculate cosTheta, phi in maybe the GJ axes.
        // since we already defined the x,y,z as TVector3 we don't have to do it again.
        TVector3 z = locBeamP4_res.Vect().Unit();
        // this y should be the normal of the production plane. If we do a boost in a direction in the production plane the perp direction doesn't change. We could use the beam and the recoiled proton to define the
        // production plane in this new frame. Let us define it in the CM frame. 
        //TVector3 y = locPi0EtaP4_cm.Vect().Cross(locBeamP4_cm.Vect()).Unit(); // OLD FORMULATION 
        TVector3 y = (locBeamP4_cm.Vect().Unit().Cross(-1*locProtonP4_cm.Vect().Unit())).Unit();
        TVector3 x = y.Cross(z).Unit();

	TVector3 angles_pi0;
	TVector3 angles_eta;
        angles_pi0.SetXYZ ( locPi0P4_res_unit.Dot(x), locPi0P4_res_unit.Dot(y), locPi0P4_res_unit.Dot(z) );
        angles_eta.SetXYZ ( locEtaP4_res_unit.Dot(x), locEtaP4_res_unit.Dot(y), locEtaP4_res_unit.Dot(z) );

        double cosTheta_pi0_GJ = angles_pi0.CosTheta();
        cosTheta_eta_GJ= angles_eta.CosTheta();
        double phi_pi0_GJ = angles_pi0.Phi()*radToDeg;
        phi_eta_GJ = angles_eta.Phi()*radToDeg;

        //////// 
        // HELICITY 
        ////////
        z = -1. * locProtonP4_res.Vect().Unit();
        //y = locPi0EtaP4_cm.Vect().Cross(locBeamP4_cm.Vect()).Unit(); // OLD FORMULATION
        y = (locBeamP4_cm.Vect().Unit().Cross(-1*locProtonP4_cm.Vect().Unit())).Unit();
        x = y.Cross(z).Unit();
        angles_pi0.SetXYZ ( locPi0P4_res_unit.Dot(x), locPi0P4_res_unit.Dot(y), locPi0P4_res_unit.Dot(z) );
        angles_eta.SetXYZ ( locEtaP4_res_unit.Dot(x), locEtaP4_res_unit.Dot(y), locEtaP4_res_unit.Dot(z) );

        double cosTheta_pi0_hel = angles_pi0.CosTheta();
        cosTheta_eta_hel= angles_eta.CosTheta();
        double phi_pi0_hel = angles_pi0.Phi()*radToDeg;
        phi_eta_hel = angles_eta.Phi()*radToDeg;
	cout << "Calculated the kinematic angles" << endl;

        locPi0EtaMass = locPi0EtaP4.M();
        mandelstam_t = -(locProtonP4-locTargetP4).M2();
	double mandelstam_abst = abs(mandelstam_t);
	//double mandelstam_t0 = -((locProtonP4.M2()-locPi0EtaP4.M2()-locTargetP4.M2())/(2*(locBeamP4+locTargetP4).M())-(locBeamP4_cm-locPi0EtaP4_cm).M2());
	//above formulation is wrong. the last term uses magnitue, not p4
	double mandelstam_t0 = -(TMath::Power(-locPi0EtaP4.M2()/(2*(locBeamP4+locTargetP4).M()),2)-TMath::Power(locBeamP4_cm.Vect().Mag()-locPi0EtaP4_cm.Vect().Mag(),2));
	mandelstam_tp = mandelstam_t-mandelstam_t0;
	
        double Metap = (locEtaP4+locProtonP4).M();
        double Mpi0p = (locPi0P4+locProtonP4).M();

	bool pBeamE8GeV = locBeamP4.E() > 8;
	bool pBeamE8288 = 8.2 < locBeamP4.E() &&  locBeamP4.E() < 8.8;
    
        TLorentzVector p_p4_kin = locProtonP4;
        TLorentzVector beam_p4_kin = locBeamP4;

        cout << "Passed selections! Saving" << endl;

        if (dFlatTreeFileName!=""){
            dFlatTreeInterface->Fill_Fundamental<Int_t>("BeamAngle", locPolarizationAngle); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tp", mandelstam_tp); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_t", mandelstam_t); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_teta", mandelstam_teta); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("mandelstam_tpi0", mandelstam_tpi0); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_gj",cosTheta_eta_GJ); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_gj",phi_eta_GJ); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("cosTheta_eta_hel",cosTheta_eta_hel); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("phi_eta_hel",phi_eta_hel); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Ebeam", locBeamP4.E()); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0eta", locPi0EtaMass); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Mpi0p", Mpi0p); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_Fundamental<Float_t>("Metap", Metap); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Fill_TObject<TLorentzVector>("p_p4_kin",p_p4_kin);
            dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_p4_kin",beam_p4_kin);
	    Fill_FlatTree(); //for the active combo
        }
	//Fill_OutputTree();
	++eventIdx;

	//OR Manually:
	//BEWARE: Do not expect the particles to be at the same array indices from one event to the next!!!!
	//Why? Because while your channel may be the same, the pions/kaons/etc. will decay differently each event.
	
	//cout << matchFOM << endl;

	//BRANCHES: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat#TTree_Format:_Simulated_Data
/*
	Particle_t locThrown1PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);
	TLorentzVector locThrown1P4 = *((TLorentzVector*)(*locP4Array)->At(0));
	cout << "Particle 1: " << locThrown1PID << ", " << locThrown1P4.Px() << ", " << locThrown1P4.Py() << ", " << locThrown1P4.Pz() << ", " << locThrown1P4.E() << endl;
	Particle_t locThrown2PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);
	TLorentzVector locThrown2P4 = *((TLorentzVector*)(*locP4Array)->At(1));
	cout << "Particle 2: " << locThrown2PID << ", " << locThrown2P4.Px() << ", " << locThrown2P4.Py() << ", " << locThrown2P4.Pz() << ", " << locThrown2P4.E() << endl;
*/


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

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}

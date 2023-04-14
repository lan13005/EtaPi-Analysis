// macro to process analysis TTree with DSelector
// We cannot just run this macro, the library doesnt load properly. We can run the following two lines of code
//.x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C
//.x runDSelector.C

#include <iostream> 
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"


//R__LOAD_LIBRARY(/group/halld/Software/builds/Linux_CentOS7-x86_64-gcc4.8.5/gluex_root_analysis/gluex_root_analysis-0.5/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so)
//R__LOAD_LIBRARY(/d/home/ln16/gluex_top/gluex_root_analysis/gluex_root_analysis_1.7.0/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so) 
R__LOAD_LIBRARY(libDSelector.so)
   
void runDSelector_7_17_14(bool useproof = 1, string path = "") 
{
	cout << "Loaded using R__LOAD_LIBRARY" << endl;
	// Load DSelector library
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	// change the directory that proof saves the data to
	//gEnv->SetValue("ProofLite.Sandbox", "/d/grid15/ln16/.proof");
	int proof_Nthreads = 36;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	TString nameOfTree = "pi0pi0__B3_F1_M7_Tree"; 
	TChain *chain = new TChain(nameOfTree);

	// **********************************************************************************	
	// ************************** ------ PI0PI0 BELOW ---------**************************	
	// **********************************************************************************	
        // This is probably the right one
	chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver20/tree_pi0pi0__B3_F1_M7/merged/tree_pi0pi0__B3_F1_M7_*");
        // These two data branches have an extra U1 flag. I think I dont want this one
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/tree_pi0pi0__B3_F1_U1_M7/merged/tree_pi0pi0__B3_F1_U1_M7_03*");
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/tree_pi0pi0__B3_F1_U1_M7/merged/tree_pi0pi0__B3_F1_U1_M7_03033*");
	//
	// MC with f2 resonance
	//chain->Add("/d/grid15/ln16/rootFiles/pi0pi0/Aug27Sep09/tree_pi0pi0__B3_F1_U1_M7.root");
	//chain->Add("/d/grid13/ln16/MC/pi0pi0_f2_interactive/hddm/combined/tree_pi0pi0__B3_F1_U1_M7.root");
        //
        // Flat 2017 MC
        //chain->Add("/d/grid17/ln16/rootFiles/pi0pi0_flat/trees/tree_pi0pi0__B3_F1_M7_*");
        
        
	// should change the name below from data to reco when running over MC
	string degAngle="degALL";
        string tag="_data_2017_mEllipse_8288_tLT1";

	string options = "";
	if(useproof) { // add TTree to chain and use PROOFLiteManager
		string outputHistFileName;
                string outputTreeFileName;
		outputHistFileName = degAngle+tag+"_hists_DSelector.root";
		outputTreeFileName = degAngle+tag+"_tree_DSelector.root"; 
	
		DPROOFLiteManager::Process_Chain(chain, "DSelector_ver20.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		chain->Process("DSelector_ver20.C+", options.data());
		
	}
	
	return;
}

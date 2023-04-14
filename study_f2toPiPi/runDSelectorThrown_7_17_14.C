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

R__LOAD_LIBRARY(libDSelector.so) 
   
void runDSelectorThrown_7_17_14(bool proof = 1, string path = "") 
{
	// Load DSelector library
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	int proof_Nthreads = 12;

	// open ROOT files and TTree
	TString nameOfTree = "Thrown_Tree"; // pi0eta__B4_Tree is the old one
	TChain *chain = new TChain(nameOfTree);

        chain->Add("/d/grid17/ln16/rootFiles/pi0pi0_flat/thrown/tree_thrown_gen_amp_*");

	string degAngle="degALL";
        string tag="_flat_2017_gen";

	string options = "";
	if(proof) { // add TTree to chain and use PROOFLiteManager
		string outputHistFileName = degAngle+tag+"_hists_DSelector.root";//_GEANT4.root");
                string outputTreeFileName = degAngle+tag+"_trees_DSelector.root";//_GEANT4.root");
		DPROOFLiteManager::Process_Chain(chain, "DSelector_thrown.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		chain->Process("DSelector_thrown.C+", options.data());
	}
	
	return;
}

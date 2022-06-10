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
   
void runDSelector_thrown(string inputFileLocation, string treeName, string outputFileNameTag, int proof_Nthreads) 
{
	// Load DSelector library
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	TChain *chain = new TChain(treeName.c_str());

        chain->Add(inputFileLocation.c_str());
        
	string outputHistFileName = outputFileNameTag+"_hists.root";
        string outputTreeFileName = outputFileNameTag+"_trees.root";
	DPROOFLiteManager::Process_Chain(chain, "DSelector_thrown.C++",  proof_Nthreads, outputHistFileName, outputTreeFileName, "");
	//chain->Process("DSelector_thrown.C++", "");
	
	return;
}

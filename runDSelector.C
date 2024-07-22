#include <iostream> 
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

R__LOAD_LIBRARY(libDSelector.so)
void runDSelector(string inputFileLocation, string treeName, string outputFileNameTag, int proof_Nthreads) 
{
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	TChain *chain = new TChain(treeName.c_str());

    chain->Add(inputFileLocation.c_str());

	string outputHistFileName=outputFileNameTag+"_hists.root"; // name of the output root file that contains the histograms
    string outputTreeFileName=outputFileNameTag+"_tree.root"; // name of the output root file that contains the tree

	//// Choice 1: run with proof with your desired number of threads
	DPROOFLiteManager::Process_Chain(chain, "DSelector_etapi.C++",  proof_Nthreads, outputHistFileName, outputTreeFileName, "");
	//// Choice 2: run interactively - useful for debugging sometimes when running over small (sub)samples
	// chain->Process("DSelector_etapi.C+");

	return;
}

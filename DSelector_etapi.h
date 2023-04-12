#ifndef DSelector_etapi_h
#define DSelector_etapi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TLorentzRotation.h"

bool filterOmega(float omega, float Mpi0eta){
    // omega should be in degrees and mass in GeV
    return -29.0*atan(-1.05*Mpi0eta+2.78)+328 > omega;
}

map<TString,Int_t> mapTopologyToInt={
    // Run DSelector with topologyCounter to determinethe counts of each MC topology
    //   Use mergeTopologyOutputs.py to output the counts AND the following lines that build the map
    {"4#gammap[#pi^{0},#eta]",0},
    {"4#gammap[2#pi^{0}]",1},
    {"6#gammap[3#pi^{0}]",2},
    {"5#gammap[2#pi^{0},#omega]",3},
    {"3#gammap[#pi^{0},#omega]",4},
    {"6#gammap[2#pi^{0},#eta]",5},
    {"4#gammap[#eta]",6},
    {"5#gammap[2#pi^{0}]",7},
    {"4#gammap[#pi^{0},#eta']",8},
    {"3#gammap[#eta,#phi]",9}
};

class DSelector_etapi : public DSelector
{
	public:

		DSelector_etapi(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_etapi(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO
		bool hasPolarizationAngle; // if event has a polarization angle then diamond radiator was used. Else amorphous radiator was used
		int locPolarizationAngle; // stores the polarization angle, -1=AMO

		bool dIsMC;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;
		DNeutralParticleHypothesis* dPhoton4Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1F* dHist_BeamEnergy;
		TH1F* dHist_Mpi0p;
		TH1F* dHist_Metap;
		TH1F* dHist_protonZ;
		TH1F* dHist_t;
		TH2F* dHist_cosThetaHelVsMetapi0;
		TH2F* dHist_cosThetaGJVsMetapi0;
		TH1F* dHist_Metapi_tot;
		TH1F* dHist_Metapi_sig;
		TH1F* dHist_Metapi_bkg;
		TH1F* dHist_Meta;
		TH1F* dHist_Mpi0;
		TH1F* dHist_rf;
		TH1F* dHist_mmsq;
		TH1F* dHist_chiSq;
		TH1F* dHist_photonThetaPi0;
		TH1F* dHist_photonThetaEta;
		TH2F* dHist_dEdx_momentum;
		TH1F* dHist_combosRemaining;
                
                //// Topology histograms
                //TH1F* dHistThrownTopologies;
                //map<TString, TH1I*> dHistInvariantMass_ThrownTopology;
                map<TString,int> topologyCount;

	ClassDef(DSelector_etapi, 0);
};

void DSelector_etapi::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

void print_sorted_map(map<TString, int> &my_map) { // Thank you Chat-GPT
    cout << "Printing topologies sorted by counts!" << endl;
    // Create a vector of pairs to store the key-value pairs of the map
    vector<pair<TString, int>> my_vector;

    // Copy the key-value pairs from the map to the vector
    for (auto &pair : my_map) {
        my_vector.push_back(pair);
    }

    // Sort the vector by the values in ascending order
    sort(my_vector.begin(), my_vector.end(),
         [](const pair<TString, int> &a, const pair<TString, int> &b) {
             return a.second < b.second;
         });

    // Print the sorted key-value pairs
    for (auto &pair : my_vector) {
        cout << pair.first << ": " << pair.second << endl;
    }
}

#endif // DSelector_etapi_h

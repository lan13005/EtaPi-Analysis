#ifndef DSelector_thrown_h
#define DSelector_thrown_h

#include <iostream>

#include "DSelector/DSelector.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TRandom.h"

void findDaughters( std::vector<int> parentArray, std::vector<int> &daughters, int selfLoc ) {
	cout << "Finding daughters of parent " << selfLoc << endl;
        for ( auto parent=0; parent<(int)parentArray.size(); ++parent ){
                if ( parentArray[parent] == selfLoc ){
                        daughters.push_back(parent);
			cout << "  found daughter " << parent << endl; 
                }
        }
}

void findParents( std::vector<int> parentArray, std::vector<int> &parents) {
	for (auto parent=0; parent<(int)parentArray.size(); ++parent) {
		if ( parentArray[parent] == -1 ){
			parents.push_back(parent);
		}
	}
}


class DSelector_thrown : public DSelector
{
	public:

		DSelector_thrown(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_thrown(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:
		// These variables are for the teta vs Meta plots
		int num_tBins=14;
		double tMin=0;
		double tMax=2.8;
		int num_massBins=12;
		const int numHists = num_tBins*num_massBins;
		double mMin=1.6;
		double mMax=2.8;
		double tStep=(tMax-tMin)/num_tBins;
		double mStep=(mMax-mMin)/num_massBins;
		int idx_t_eta;
		int idx_t_pi0;
		int idx_m;
		double mandelstam_teta;
		double mandelstam_tpi0;
		double teta_genCounts;
		double tpi0_genCounts;

		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO
                bool hasPolarizationAngle; // true if is polarized but false if through deg5ous radiator or no data. Under what cicumstances is the second one true.
                int locPolarizationAngle; // actual polarization angle
		bool keepPolarization;
                TH1F* dHist_BeamAngle;
                TH1F* dHist_SelectedBeamAngle;

                int equal=0;
                int notEqual=0;
                int notPass=0;
                int pass=0;
		int beamEinf=0;
		int beamEfin=0;
                vector<Int_t> vecEqual;
                vector<Int_t> vecNotEqual;
                vector<Int_t> vecNotPassCuts;
                vector<Int_t> vecHasSecondaryDaughters;
                vector<Int_t> vecNumThrownNot7;
                Int_t eventIdx=0;
		TRandom* rgen = new TRandom(1992);
		TH1F* countThrownEvents;
		TH1F* dHist_prodPlanePS_000;
		TH1F* dHist_prodPlanePS_045;
		TH1F* dHist_prodPlanePS_090;
		TH1F* dHist_prodPlanePS_135;
		TH1F* dHist_prodPlanePS_AMO;
		TH1F* dHist_prodPlanePS_000_rejSamp;
		TH1F* dHist_prodPlanePS_045_rejSamp;
		TH1F* dHist_prodPlanePS_090_rejSamp;
		TH1F* dHist_prodPlanePS_135_rejSamp;
		TH1F* dHist_prodPlanePS_AMO_rejSamp;
                TH1I* dHist_PID;
                TH1I* dHist_NumThrown;
		TH1F *mandelstam_tpAll; 
		TH1F *mandelstam_tpAll_selected;
		TH1F *mandelstam_tAll; 
		TH1F *mandelstam_tpLT1;
		TH1F *mandelstam_tpLT06;
		TH1F *mandelstam_tpGT05LT1;
		TH2F *dHist_cosThetaVsMass_tpAll;
		TH2F *dHist_cosThetaVsMass_tpLT1;
		TH2F *dHist_cosThetaVsMass_tpLT06;
		TH2F *dHist_cosThetaVsMass_tpGT05LT1;
		TH2F *dHist_phiVsMass;
		TH1F *dHist_phi;
		TH1F *dHist_cosTheta;
		TH1F *dHist_beamE;
		TH1F *dHist_beamECut;
		TH1F *dHist_numEventsOnePi0OneEta;
        	TH1F *dHist_genCounts_eta_tAll;
        	TH1F *dHist_genCounts_pi0_tAll;
        	TH1F *dHist_genCounts_eta_tGT1;
        	TH1F *dHist_genCounts_pi0_tGT1;
        	TH1F *dHist_genCounts_eta_tLT05;
        	TH1F *dHist_genCounts_pi0_tLT05;
        	TH1F *dHist_genCounts_eta_tGT05LT1;
        	TH1F *dHist_genCounts_pi0_tGT05LT1;

                TH1F* dHist_Mpi0eta;
		bool pBeamE[12];
		bool pBeamE8to9GeV;
		TH1F *dHist_pi0eta1DBeam[12];
		TH1F *dHist_pi0eta1D;
		TH1F *dHist_phi8GeVPlus;
		TH1F *dHist_cosTheta8GeVPlus;
		double mandelstam_tp;


		int maxevent;
		int ievent;
		set<Int_t> showOutput =  {6, 82,188};

	ClassDef(DSelector_thrown, 0);
};

#endif // DSelector_thrown_h

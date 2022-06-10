#ifndef DSelector_thrown_h
#define DSelector_thrown_h

#include <iostream>

#include "DSelector/DSelector.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_thrown : public DSelector
{
	public:

		DSelector_thrown(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_thrown(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO
		bool hasPolarizationAngle; // if event has a polarization angle then diamond radiator was used. Else amorphous radiator was used
		int locPolarizationAngle; // stores the polarization angle, -1=AMO

                set<TString> topologies;

	ClassDef(DSelector_thrown, 0);
};

#endif // DSelector_thrown_h

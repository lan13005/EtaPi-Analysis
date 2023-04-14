#ifndef DSelector_ver20_h
#define DSelector_ver20_h

#include <iostream>
#include <fstream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TRandom.h"

#include "DSelector_helperFuncs.h"

bool showThrownTopology=false;

// Used to look for the daughters of a specific thrown particle 
void findDaughters( std::vector<int> parentArray, std::vector<int> &daughters, int selfLoc ) {
        for ( auto parent=0; parent<(int)parentArray.size(); ++parent ){
                if ( parentArray[parent] == selfLoc ){
                        daughters.push_back(parent);
                }
        }
}
// used to look for primary particles
void findParents( std::vector<int> parentArray, std::vector<int> &parents) {
	for (auto parent=0; parent<(int)parentArray.size(); ++parent) {
		if ( parentArray[parent] == -1 ){
			parents.push_back(parent);
		}
	}
}
// Getting the thrown particle (if any) that a reconstrcuted particle, like a photon of a pi0, matched to.
Int_t getParents( Int_t* thrownID, vector<Int_t> parentIDs, vector<Int_t> thrownPIDs, TString* string_ph ){
        Int_t initialThrownID = *thrownID;
        Int_t thrownPID = thrownPIDs[*thrownID];
        *thrownID  = parentIDs[*thrownID];
        if (*thrownID != -1){
                Int_t parentPID = thrownPIDs[*thrownID];
                cout << "ThrownID: " << initialThrownID << " with PID=" << thrownPID << " has ParentID: " << *thrownID << " with PID=" << parentPID << endl;
                *string_ph += "("+to_string(parentPID)+")";
                return *thrownID;
        }
        else {
                cout << "Found primary particle" << endl;
                //*string_ph += "(-1)";
                return *thrownID;
        }
}

// saves the topology information if using simulated data. We can try to track what thrown particles make up the reconstructed particles
struct topology {
        TString locThrownTopology;
        TString composition;
        TString beamProtonID;
        TString spectroscopicID;
        double chiSq;
        double unusedEnergy;
        Int_t nUnusedShowers;
        double pi0Mass;
        double etaMass;
        double pi0etaMass;
};



//void withinBox(bool inBox[], bool inBox_noOtherCuts[],bool additionalCut, double x, double y, double xmin, double xmax, double ymin, double ymax, double xskip, double yskip){
//	// regions are:
//	//  0 1 2
//	//  3 4 5 
//	//  6 7 8
//	//  where xmin, xmax, ymin,ymax all belong to region 5
//	//  the 10th element is the interesection of the negation of all regions 
//	//  11th element is 1345 and 12th element is 0268
//	double xlength = xmax-xmin;
//	double ylength = ymax-ymin;
//	inBox[4] = x<xmax && x>xmin && y<ymax && y>ymin;
//	inBox[3] = x<(xmin-xskip) && x>(xmin-xskip-xlength/2) && y<ymax && y>ymin;
//	inBox[5] = x<(xmax+xskip+xlength/2) && x>(xmax+xskip) && y<ymax && y>ymin;
//	inBox[1] = x<xmax && x>xmin && y<(ymax+yskip+ylength) && y>(ymax+yskip);
//	inBox[7] = x<xmax && x>xmin && y>(ymin-yskip-ylength) && y<(ymin-yskip);
//	inBox[0] = x<(xmin-xskip) && x>(xmin-xskip-xlength/2) && y<(ymax+yskip+ylength) && y>(ymax+yskip);
//	inBox[2] = x<(xmax+xskip+xlength/2) && x>(xmax+xskip) && y<(ymax+yskip+ylength) && y>(ymax+yskip);
//	inBox[6] = x<(xmin-xskip) && x>(xmin-xskip-xlength/2) &&  y>(ymin-yskip-ylength) && y<(ymin-yskip);
//	inBox[8] = x<(xmax+xskip+xlength/2) && x>(xmax+xskip) && y>(ymin-yskip-ylength) && y<(ymin-yskip);
//	inBox[9] = !inBox[0] * !inBox[1] * !inBox[2] * !inBox[3] * !inBox[4] * !inBox[5] * !inBox[6] * !inBox[7] * !inBox[8];	
//	inBox[10] = inBox[0] || inBox[2] || inBox[6] || inBox[8];
//	inBox[11] = inBox[1] || inBox[7];
//       	inBox[12] = inBox[3] || inBox[5];
//
//	for (int i=0; i<13; ++i){
//		inBox_noOtherCuts[i] = inBox[i]; // save the inBox bools before modifiying them by the cut. We will use these for defining the weightBS
//		inBox[i]*=additionalCut;
//	}
//}

void withinBox(bool inBox[], bool inBox_noOtherCuts[],bool additionalCut, double pi0Mass, double etaMass, double pi0Mean, double etaMean, double pi0Std, 
                    double etaStd, double pi0Sig, double pi0Skip, double pi0SB, double etaSig, double etaSkip, double etaSB){
        // Sig/Skip/SB are all in terms of stdevs
	// regions are:
	//  0 1 2
	//  3 4 5 
	//  6 7 8
	//  where xmin, xmax, ymin,ymax all belong to region 5
	//  the 10th element is the interesection of the negation of all regions 
	//  11th element is 1345 and 12th element is 0268
        
        double xmin=pi0Mean-pi0Std*pi0Sig;
        double xmax=pi0Mean+pi0Std*pi0Sig;
        double ymin=etaMean-etaStd*etaSig;
        double ymax=etaMean+etaStd*etaSig;
        double xskip=pi0Std*pi0Skip;
        double yskip=etaStd*etaSkip;
        double xsb=pi0Std*pi0SB;
        double ysb=etaStd*etaSB;
        double x=pi0Mass;
        double y=etaMass;

	inBox[4] = x<xmax && x>xmin && y<ymax && y>ymin;
	inBox[3] = x<(xmin-xskip) && x>(xmin-xskip-xsb) && y<ymax && y>ymin;
	inBox[5] = x<(xmax+xskip+xsb) && x>(xmax+xskip) && y<ymax && y>ymin;
	inBox[1] = x<xmax && x>xmin && y<(ymax+yskip+ysb) && y>(ymax+yskip);
	inBox[7] = x<xmax && x>xmin && y>(ymin-yskip-ysb) && y<(ymin-yskip);
	inBox[0] = x<(xmin-xskip) && x>(xmin-xskip-xsb) && y<(ymax+yskip+ysb) && y>(ymax+yskip);
	inBox[2] = x<(xmax+xskip+xsb) && x>(xmax+xskip) && y<(ymax+yskip+ysb) && y>(ymax+yskip);
	inBox[6] = x<(xmin-xskip) && x>(xmin-xskip-xsb) &&  y>(ymin-yskip-ysb) && y<(ymin-yskip);
	inBox[8] = x<(xmax+xskip+xsb) && x>(xmax+xskip) && y>(ymin-yskip-ysb) && y<(ymin-yskip);
        // skip
	inBox[9] = !inBox[0] * !inBox[1] * !inBox[2] * !inBox[3] * !inBox[4] * !inBox[5] * !inBox[6] * !inBox[7] * !inBox[8];	
        // corner
	inBox[10] = inBox[0] || inBox[2] || inBox[6] || inBox[8];
        // eta sideband
	inBox[11] = inBox[1] || inBox[7];
        // pi0 sideband
       	inBox[12] = inBox[3] || inBox[5];

	for (int i=0; i<13; ++i){
		inBox_noOtherCuts[i] = inBox[i]; // save the inBox bools before modifiying them by the cut. We will use these for defining the weightBS
		inBox[i]*=additionalCut;
	}
}


struct histDef_1D{
	TH1F* hist;
	string name;
	bool* cut;
	std::vector< double* > values;
	double* weights; 

	bool saveHistValues=false;
	string baseLocation = "/d/grid15/ln16/pi0eta/092419/";

	void Write(double value, double weight){
		if (saveHistValues){
			// if it doesnt exist
			if ( gSystem->AccessPathName((baseLocation+"newGraphs_histValues").c_str() ) ){
				cout << "Making newGraphs_histValues" << endl;
				gSystem->Exec(("mkdir "+baseLocation+"newGraphs_histValues").c_str());
			}
			if ( gSystem->AccessPathName( ( (baseLocation+"newGraphs_histValues/"+name+".txt").c_str() ) ) ){
				cout << baseLocation << "newGraphs_histValues/" << name << ".txt not found, making it" << endl;
				gSystem->Exec(("touch "+baseLocation+"newGraphs_histValues/"+name+".txt").c_str());
			}
			cout << ("echo "+to_string(value)+" "+to_string(weight)+" >> "+baseLocation+"newGraphs_histValues/"+name+".txt").c_str() << endl;
			cout << "Value is actually: " << value << endl;
			gSystem->Exec(("echo "+to_string(value)+" "+to_string(weight)+" >> "+baseLocation+"newGraphs_histValues/"+name+".txt").c_str());
		}	
	}

	void clear(){
	    	values.clear();
	}
		

};

struct histDef_2D{
    TH2F* hist;
    string name;
    bool* cut;
    std::vector< double* > valuesX;
    std::vector< double* > valuesY;
    double* weights; 
    void clear(){
        valuesX.clear();
        valuesY.clear();
    }
};

class trackingGroup{
    private:
	//set< pair<Int_t,Int_t> > usedIds;
        set< pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > > usedPairMapIds;
        set< map<Particle_t, set<Int_t> > >  usedMapIds;
	//string name;

    public:
        TCanvas *anyCanvas;
        trackingGroup(){ cout << "Calling constructor" << endl; };
        ~trackingGroup(){ cout << "Calling destructor" << endl; };

        void setNameMakeCanvas(string groupName){
            anyCanvas = new TCanvas((groupName).c_str(),"",1440, 900);
        }

        std::vector<histDef_1D> allHists_1D;
        std::vector<histDef_2D> allHists_2D;
	std::vector<std::vector<double>> allValues_1D;
	std::vector<std::vector<double>> allValues_1D_weight;
	std::vector<std::vector<double>> allValues_2DX;
	std::vector<std::vector<double>> allValues_2DY;
	std::vector<std::vector<double>> allValues_2D_weight;
	std::vector<double> emptyVector;

        std::vector< set< map<Particle_t, set<Int_t> > > > allUsedMapIds_1D;
        std::vector< set< map<Particle_t, set<Int_t> > > > allUsedMapIds_2D;

        // every time we add another histogram to track we will push back another tracking set to track it.
        void insert(histDef_1D hist){
            allHists_1D.push_back(hist);
	    allValues_1D.push_back(emptyVector);
	    allValues_1D_weight.push_back(emptyVector);
            allUsedMapIds_1D.push_back(usedMapIds);
        }
        void insert_2D(histDef_2D hist){
            // ***** NOTE THAT IT IS BEST IF WE DO NOT MIX USING BOTH OF TYPES OF TRACKING TYPES ***** NOT WORKED ON YET
            allHists_2D.push_back(hist);
	    allValues_2DX.push_back(emptyVector);
	    allValues_2DY.push_back(emptyVector);
	    allValues_2D_weight.push_back(emptyVector);
            allUsedMapIds_2D.push_back(usedMapIds);
        }

        void fillHistograms_vectorMap( std::vector< map<Particle_t, set<Int_t> > > beingUsedPairIds){
            // *********************** THIS COMMENTED CODE IS TO CHECK THE OUTPUT OF THE DIJ3 FCAL FILLING TO MAKE SURE WE ARE DOING THINGS CORRECLTY ***********************************
            //  cout << "size of {tracking,value}={"<<beingUsedPairIds.size()<<","<<allHists_1D[0].values.size()<<"}"<<endl;
            // for (UInt_t iValue=0; iValue<beingUsedPairIds.size(); ++iValue){
            //    for ( auto elem : beingUsedPairIds[iValue] ){
            //        cout << elem.first << "| ";
            //        for (auto it=elem.second.begin(); it != elem.second.end(); ++it){
            //            cout << *it << " ";
            //        }
            //    }
            //bool uniqueBool= allUsedMapIds_1D[0].find(beingUsedPairIds[iValue])==allUsedMapIds_1D[0].end();
            //bool cutBool =  *(allHists_1D[0].cut);
            //cout << " -- " << *(allHists_1D[0].values[iValue]) << " Unique? " << uniqueBool << " passCut? " << cutBool << endl; 
            //}
                
             for (UInt_t iValue=0; iValue<beingUsedPairIds.size(); ++iValue){
                 for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
	              //if (allUsedMapIds_1D[iHist].find(beingUsedPairIds[iValue])==allUsedMapIds_1D[iHist].end() && *(allHists_1D[iHist].cut)  ){
	              //   allUsedMapIds_1D[iHist].insert(beingUsedPairIds[iValue]); //we get a iterator which references the element of the set so we need to dereference.              
	              if (*(allHists_1D[iHist].cut)  ){
                         allHists_1D[iHist].hist->Fill( *(allHists_1D[iHist].values[iValue]), *(allHists_1D[iHist].weights) );
                      }
                 }
             }

             for (UInt_t iValue=0; iValue<beingUsedPairIds.size(); ++iValue){
                 for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
	              //if (allUsedMapIds_2D[iHist].find(beingUsedPairIds[iValue])==allUsedMapIds_2D[iHist].end() && *(allHists_2D[iHist].cut)  ){
	              //   allUsedMapIds_2D[iHist].insert(beingUsedPairIds[iValue]); //we get a iterator which references the element of the set so we need to dereference.              
	              if ( *(allHists_2D[iHist].cut)  ){
                         allHists_2D[iHist].hist->Fill( *(allHists_2D[iHist].valuesX[iValue]),*(allHists_2D[iHist].valuesY[iValue]), *(allHists_2D[iHist].weights) );
                      }
                 }
             }
            clear_values();
        }

        // use when we are starting a new event. Only when using default uniqueness tracking
        void clear_tracking(){
            for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
                allUsedMapIds_1D[iHist].clear();
            }
            for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
                allUsedMapIds_2D[iHist].clear();
            }
        }
        

        // Similar to above but just for 1D distributions where we don't need to track a pair of maps for correlations
        void saveHistograms(){
	     //cout << "Starting to fill hists in group: " << name << endl;
             for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
	        if ( *(allHists_1D[iHist].cut)  ){
		    //    cout << "(" << allHists_1D[iHist].name << ") values in order: " << *val << endl;
                    for (UInt_t ival=0; ival<allHists_1D[iHist].values.size(); ++ival){
                        allValues_1D[iHist].push_back(*(allHists_1D[iHist].values[ival]));
		        allValues_1D_weight[iHist].push_back( *(allHists_1D[iHist].weights) );
		    //allHists_1D[iHist].Write( value, *(allHists_1D[iHist].weights) );
                    }
                }
             }
             for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
	         if ( *(allHists_2D[iHist].cut)  ){
                    for (UInt_t ival=0; ival<allHists_2D[iHist].valuesX.size(); ++ival){
                        allValues_2DX[iHist].push_back(*(allHists_2D[iHist].valuesX[ival]));
                        allValues_2DY[iHist].push_back(*(allHists_2D[iHist].valuesY[ival]));
		        allValues_2D_weight[iHist].push_back( *(allHists_2D[iHist].weights) );
                    }
                }
             }
             // Definitely dont want to clear histDef since those are the addresses of the fill variables
            //clear_values();
        }

	void fillHistograms(){
		// use back and pop_back to fill the histograms, slowly consuming the values that passed the cuts
		for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
			//cout << "Filling some histograms..." << endl;
                        //for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
	                //   if (allUsedMapIds_1D[iHist].find(beingUsedMap)==allUsedMapIds_1D[iHist].end() && *(allHists_1D[iHist].cut)  ){
	                //       allUsedMapIds_1D[iHist].insert(beingUsedMap); //we get a iterator which references the element of the set so we need to dereference 
			if (allValues_1D[iHist].size() > 0) {
				int numVals = allValues_1D[iHist].size();
				for ( int iValue=0; iValue<numVals; ++iValue){
					if (allHists_1D[iHist].name=="pi0Mass_Kin_noCut" ){
						cout << "Filling 1D hist: " << allHists_1D[iHist].name << " with " << allValues_1D[iHist].back() <<
							" with weight=" << allValues_1D_weight[iHist].back() << 
							" -- nCombos=" << numVals << endl;
					}
		    			allHists_1D[iHist].hist->Fill( allValues_1D[iHist].back(), allValues_1D_weight[iHist].back() );
		    			allValues_1D[iHist].pop_back();
		    			allValues_1D_weight[iHist].pop_back();
				}
			}
		}
		for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
			//cout << "Filling some histograms 2D..." << endl;
			if (allValues_2DX[iHist].size() > 0) {
				//cout << allValues_2DX[iHist].back() << allValues_2D_weight[iHist].back() << endl;
				int numVals = allValues_2DX[iHist].size();
				for ( int iValue=0; iValue<numVals; ++iValue){
		    			allHists_2D[iHist].hist->Fill( allValues_2DX[iHist].back(), allValues_2DY[iHist].back(), allValues_2D_weight[iHist].back() );
					//cout << "Filling 2D hist: " << allHists_2D[iHist].name << " with " << allValues_2DX[iHist].back() << ", " << allValues_2DY[iHist].back() << 
					//	" with weight=" << allValues_2D_weight[iHist].back() << endl;
		    			allValues_2DX[iHist].pop_back();
		    			allValues_2DY[iHist].pop_back();
		    			allValues_2D_weight[iHist].pop_back();
				}
			}
		}
	}

        // have to use this everytime we fillHistogram since we have a vector of values we input, i.e. photonXs
        void clear_values(){
            	for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
            	    allHists_1D[iHist].clear(); // calling histDef's clear function to remove values
            	}
            	for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
            	    allHists_2D[iHist].clear();
		}
        }
};

    
        

class DSelector_ver20 : public DSelector
{
	public:
		DSelector_ver20(TTree* locTree = NULL) : DSelector(locTree){
                }
		virtual ~DSelector_ver20(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:
                trackingGroup group_PhNB;
                trackingGroup group_pairFCAL;
                trackingGroup group_pairBCAL;
                trackingGroup groupHists;

		void Get_ComboWrappers(void);
		void Finalize(void);

                // BEAM POLARIZATION INFORMATION
                UInt_t dPreviousRunNumber;
                bool dIsPolarizedFlag; //else is AMO
                bool dIsPARAFlag; //else is PERP or AMO
                bool hasPolarizationAngle; // true if is polarized but false if through deg5ous radiator or no data. Under what cicumstances is the second one true.
                int locPolarizationAngle; // actual polarization angle
		bool keepPolarization=true;
		bool keepPolarization000=true;
		bool keepPolarization045=true;
		bool keepPolarization090=true;
		bool keepPolarization135=true;
		bool keepPolarizationAMO=true;
                set<UInt_t> usedRuns; //think we are over counting our filling of beam angles since the below condition only checks against previous run, doesn't work for recurrences.


		// These variables are for the teta vs Meta plots
		int num_tBins=14;
		double tMin=0;
		double tMax=2.8;
		static const int num_massBins=12;
		const int numHists = num_tBins*num_massBins;
		double mMin=1.6;
		double mMax=2.8;
		double tStep=(tMax-tMin)/num_tBins;
		double mStep=(mMax-mMin)/num_massBins;
		double teta_recCounts=0;
		double tpi0_recCounts=0;
		int idx_t_eta;
		int idx_t_pi0;
		int idx_m;
		bool passMassBin_tetaIntegrated[num_massBins];
		bool passMassBin_tpi0Integrated[num_massBins];


		// Variables are for binning Meta in terms of Mpi0. Hopefully to understand Q-values more
		static const int num_mpi0Bins=10;
		bool passMpi0Bin[num_mpi0Bins];


                Int_t eventIdx=0;
                bool isNotRepeated_eta=true;
                bool isNotRepeated_pi0=true;
                bool isNotRepeated_pi0g3=true;
                bool isNotRepeated_pi0g4=true;
                bool isNotRepeated_pi0eta=true;
                bool isNotRepeated_pi0_pi0eta=true;
                bool isNotRepeated_eta_pi0eta=true;
		
                string uniqueSpectroscopicPi0EtaID;
                string uniqueSpectroscopicEtaID;
                string uniqueSpectroscopicPi0ID;
		ULong64_t spectroscopicComboID=0;
		int digitsInSpectroscopicComboID;
		int digitsInEvent=0;
		int digitsInCombo=0;
		int digitsInRun=0;
		int maxDigitsInEvent=7;
		int maxDigitsInCombo=2;
		int maxDigitsInRun=5;
		string paddedCombo="";
		string paddedEvent="";
		string paddedRun="";


                bool isTruePi0Eta=false;
                Int_t numSpect=5; // 1 proton, 4 photons
                string spectroscopicID="notThrown"; // just setting this as the default for data. If MC then thrown data exists and this will be updated

		Int_t uniqueComboID=0;

                //string degAngle = "deg0";

                // ANALYZE CUT ACTIONS
                // // Automatically makes mass histograms where one cut is missing
                DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		TH1F* dHistThrownTopologies;
		map<TString, TH1I*> dHistInvariantMass_ThrownTopology;

                //CREATE REACTION-SPECIFIC PARTICLE ARRAYS

                //Step 0
                DParticleComboStep* dStep0Wrapper;
                DBeamParticle* dComboBeamWrapper;
                DChargedTrackHypothesis* dProtonWrapper;

                //Step 1
                DParticleComboStep* dStep1Wrapper;
                DKinematicData* dDecayingPi0Wrapper;
                DNeutralParticleHypothesis* dPhoton1Wrapper;
                DNeutralParticleHypothesis* dPhoton2Wrapper;

                //Step 2
                DParticleComboStep* dStep2Wrapper;
                DKinematicData* dDecayingEtaWrapper;
                DNeutralParticleHypothesis* dPhoton3Wrapper;
                DNeutralParticleHypothesis* dPhoton4Wrapper;

                // DEFINE YOUR HISTOGRAMS HERE
                // EXAMPLES:
		TVector3 targetCenter;
	
		// This is for calculating the masses with different starting points
		double locEtaMass_charged=1;
		double locPi0Mass_charged=1;
		double locEtaMass_target=1;
		double locPi0Mass_target=1;

		TH1F* dHist_mandelstam_t_thrown;
		TH1F* dHist_mandelstam_t0_thrown;
		TH1F* dHist_mandelstam_tp_thrown;
		double mandelstam_t_thrown;
		double mandelstam_t0_thrown;
		double mandelstam_tp_thrown;
		double mandelstam_t0_oldForm;
		double mandelstam_tp_oldForm;
		double Ebeam_thrown;
                double Mpi0eta_thrown;

		TH1F* countThrownEvents;
		//TH1I* dHist_numCombos;
		TH1F* dHist_thrown_tp;
		TH1F* dHist_thrown_tp_selected;
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
                TH1F* dHist_BeamAngle;
                TH1F* dHist_BeamAngleMPE;
		TH1F* dHist_Cuts;
		TH2F* dHist_checkEllipseBS[3];
		const char *cutNames[17] = {"pShowerQuality","pBeamE8GeVPlus","pUnusedEnergy","pChiSq" ,"pDeltaTRF","pdij3pass","pPhotonE","pPhotonTheta","pMagP3Proton","pzCutmin","pRProton","pMissingMassSquared","pdEdxCDCProton","pinsideEllipse", "allGeneralCutsPassed", "mMPi0P14","mMandelstamT_mBeamE8GeVPlus"};


// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//
// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//

	        // cutString will be used by alot of histograms when defining names 
	        std::string cutString;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////// ********* DEFINING VALUES TO FILL HISTS WITH ****************/////////////////////////////////////////////

		// will be used to bin the kinematic phi and t distributions
		double lowMass;
		double upMass;
		// we need to use static const since we will use this variable to make the array. This region of the code is probably the "file scope" and declares things at compile time instead of run time. 
		// The array is in this region so it must take in a static object where the preprocessor will just replace the variables. 
		static const int numBinsMass = 13;
		double binScale;

		double lowMass_t = 0.8;
		double upMass_t = 1.8;
		static const int numBinsMass_t = 10;
		double binScale_t = (upMass_t-lowMass_t)/numBinsMass_t;

		double radToDeg = 180./TMath::Pi();
		double locWherePhoton=1;

		// Accidental subtraction variables. applyAccSub will either be = to weight or noAccSub=1.
		double weightAS=1;
		double noWeight=1;
                double weightAS_BS=1;
                double weightAS_B=1;
                double appliedWeight=1;

		//************* Other variables
		double locMissingMassSquared=1;
		double locBeamE=1;
		double locCLKinFit=1;
		double locUnusedEnergy=1;
		double locNumExtraNeutralShowers=1;
                UChar_t locNumUnusedShowers;
		double locChiSqKinFit=1;
		double locChiSqKinFit_gpi0pi0=1;
		double locDOFKinFit=1;

		//************ Particle related variables  
		double locEtaE_Kin=1;
		double locPi0E_Kin=1;
		double locEtaMass_Kin=1;
		double locPi0Mass_Kin=1;
		double mismatchPairMass_13=1;
		double mismatchPairMass_14=1;
		double mismatchPairMass_23=1;
		double mismatchPairMass_24=1;
		double mismatchPairMass_34=1;
                double mismatchPairMass_132=1;
                double mismatchPairMass_134=1;
                double mismatchPairMass_142=1;
                double mismatchPairMass_143=1;
                double mismatchPairMass_231=1;
                double mismatchPairMass_234=1;
                double mismatchPairMass_241=1;
                double mismatchPairMass_243=1;
                double mismatchPairMass_341=1;
                double mismatchPairMass_342=1;
                double mismatchPairMass_123=1;
                double mismatchPairMass_124=1;
                double mismatchPi0Mean=0.135781;
                double mismatchPi0Std=0.006677;
                bool rejectPi0PairMismatch=1;
		double locEtaProton_Kin=1;
		double locPi0Proton_Kin=1;
		double locPi0Eta_Kin=1;
		double locPi0Eta_thrown=1;
		double locPi0Eta_resolution=1;
		double locPi0Eta=1;

		double locEtaMass=1;
		double locPi0Mass=1;
		//************ Shower shape variables
		std::vector<double> E1E9_FCAL={1,1,1,1};
		std::vector<double> E9E25_FCAL={1,1,1,1};
		std::vector<double> SumU_FCAL={1,1,1,1};
		std::vector<double> SumV_FCAL={1,1,1,1};
		std::vector<double> Energy_BCALPreshower={1,1,1,1};
		std::vector<double> Energy_BCAL={1,1,1,1};
		std::vector<double> SigLong_BCAL={1,1,1,1};
		std::vector<double> SigTrans_BCAL={1,1,1,1};
		std::vector<double> SigTheta_BCAL={1,1,1,1};
		std::vector<double> DeltaPhi_BCAL={1,1,1,1};
		std::vector<double> DeltaZ_BCAL={1,1,1,1};
		std::vector<double> showerQuality_FCAL={1,1,1,1};
		std::vector<double> DOCA_FCAL={1,1,1,1};
		double locSigTheta_BCAL_proton;
		double locSigTrans_BCAL_proton;
		double locSigLong_BCAL_proton;
		double locE1E9_FCAL_proton;
		double locE9E25_FCAL_proton;
		double locSumU_FCAL_proton;
		double locSumV_FCAL_proton;
		double locEnergy_BCALPreshower_proton;
		double locEnergy_BCAL_proton;

		double pi0DetectedIn;
		double etaDetectedIn;

		//************ Charged Track
		double locPtProton=1;
		double locPzProton=1;
		double locPolarAngleProton=1;
		double locXProton=1;
		double locYProton=1;
		double locRProton=1;
		double locdzProton=1;
		double locdEdxCDCProton=1;
		double locdEdxFDCProton=1;
		double locMagP3Proton=1;
                double protonX4[4]={1,1,1,1};

		std::vector<double> massGammaPi0={1,1};
		std::vector<double> massGammaEta={1,1};

		//************ Neutral Track
		std::vector<double> photonThetas={1,1,1,1};
		std::vector<double> photonEnergies={1,1,1,1};
		std::vector<double> photonPhis={1,1,1,1};
		std::vector<double> photonXs_Kin={1,1,1,1};
		std::vector<double> photonYs_Kin={1,1,1,1};
		std::vector<double> photonZs_Kin={1,1,1,1};
		std::vector<double> photonTs_Kin={1,1,1,1};
		std::vector<double> photonThetas_fromX4_meas={1,1,1,1};
		std::vector<double> photonThetas_Shower={1,1,1,1};
		std::vector<double> photonThetas_meas={1,1,1,1};
		std::vector<double> photonXs_Shower={1,1,1,1};
		std::vector<double> photonYs_Shower={1,1,1,1};
		std::vector<double> photonZs_Shower={1,1,1,1};
		std::vector<double> photonTs_Shower={1,1,1,1};
		std::vector<double> photonDeltaTs={1,1,1,1};
		std::vector<double> photonDetectedSyss={1,1,1,1};

		bool dIsMC;

		// distance, angle, Z distance, phi between photon pairs. They must be vectors since we are not sure about the size of them yet so this will allow us to push_back
		//std::vector<double> dij3Vec;
		//std::vector<double> angle_ijVec; 
		//std::vector<double> deltaZ_ijVec; 
		//std::vector<double> deltaPhi_ijVec;
		// single version of the above
		double locPhotonDijFCAL=1;
		double locPhotonDijBCAL=1;
		double locPhotonAij=1;
		double locPhotonZij=1;
		double locPhotonPij=1;
		
		//*********** Kinematic variables
		double locDecayPlaneTheta=1;
		double locPhi_eta=1;
		double locPhi_pi0=1;
		// Calculating kinematic variables like t and cosTheta
		double mandelstam_teta=1;
		double mandelstam_teta_meas=1;
		double mandelstam_tpi0=1;
		double mandelstam_tpi0_meas=1;
		double mandelstam_tp=1;
		double mandelstam_tp_pe=1;
		double mandelstam_t=1;
		double mandelstam_abst=1;
		double mandelstam_t_pe=1;
		double mandelstam_t0=1;
		// Calculate cosTheta in maybe the gottfried-jackson frame.
		double theta_pi0_lab=1;
		double theta_eta_lab=1;
		double cosTheta_decayPlane_hel=1;
		double phi_decayPlane_hel=1;
		double cosTheta_decayPlane_GJ=1;
		double cosTheta_largestEinEta_GJ=1;
		double cosTheta_largestEinPi0_CM=1;
		double phi_decayPlane_GJ=1;
		double theta_eta_hel=1;
                double Phi=1;
                double phi_pi0eta_lab=1;
		double theta_pi0_GJ=1;
		double phi_pi0_GJ=1;
		double cosTheta_eta_hel=1;
		double phi_eta_hel=1;
		double theta_eta_GJ=1;
		double phi_eta_GJ=1;
		double phi_eta_GJ_meas=1;
		double cosTheta_pi0_GJ=1;
		double cosTheta_eta_GJ=1;
		double cosTheta_eta_GJ_meas=1;
		double phi_pi0eta_GJ=1;
		double phi_X_relativeToBeamPol=1;
                // In the omega rest frame
                double ph123Rest_angle_pi0_g3=1;
                double ph124Rest_angle_pi0_g4=1;
                double ph123Rest_angle_g12=1;
                double ph123Rest_angle_g13=1;
                double ph123Rest_angle_g14=1;
                double ph123Rest_angle_g23=1;
                double ph123Rest_angle_g24=1;
                double ph123Rest_angle_g34=1;
                double ph124Rest_angle_g12=1;
                double ph124Rest_angle_g13=1;
                double ph124Rest_angle_g14=1;
                double ph124Rest_angle_g23=1;
                double ph124Rest_angle_g24=1;
                double ph124Rest_angle_g34=1;

		double vanHove_x;
		double vanHove_y;
                double q;
                double pi0_cmZ;
                double eta_cmZ;
                double recoil_cmZ;
                double omega;


		double locYDotZ_GJ=1;
		double angleBetweenPi0Eta=1;
		int delta[4] = {5,10};
		bool withinCone[4]={true,true,true,true};
		bool pi0_inCone[4]={true,true,true,true};
		bool eta_inCone[4]={true,true,true,true};
		bool largeAngle[4]={true,true,true,true};
		std::vector<double> countCone = {0,1,2,3,4};
		// calculating the cosTheta of pi0eta system, pi0, eta  in the CM framei
		double cosTheta_pi0eta_CM=1;
		double cosTheta_pi0eta_CM_meas=1;
		double cosTheta_pi0_CM=1;
		double cosTheta_eta_CM=1;
		double phi_pi0eta_CM=1;
		double phi_pi0_CM=1;
		double phi_eta_CM=1;
		double theta_pi0_CM=1;
		double theta_eta_CM=1;
		
		//*********** Timing variables, first section is about the proton
		// shifted relative to the beam
		double RFtime=1;
		double RFtime_meas=1;
		double RFtimeProton=1;
		double locDeltaTRF=1;
		double locDeltaTRF_meas=1;
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////// ********* DEFINING ADDITIONAL CUTS TO APPLY TO GRAPHS ****************/////////////////////////////////////////////
		
		
		// ******************************** DEFINING CUT THRESHOLDS ***************************************
		// Beam asymmetry stuff
		double Emin = 8.2;
		double Emax = 8.8;
		TRandom* rgen = new TRandom(1992);
		
		/////////////// Charged Track Cuts/////////
		double Rcut = 2;
                // loose window. There is an alternative that loses ~10% of events
		//double zCutmin = 42; double zCutmax = 82; // 52-78 is the new 2021 recommeneded
		double zCutmin = 52; double zCutmax = 78; // 52-78 is the new 2021 recommeneded
		double dEdxCut; // not static so we need to recompute it. TMath::Power(10,-6)*(0.9+TMath::Exp(3.0-3.3*locMagP3Proton/.93827)); // The bigger then number multiplying MagP3 the sharper the cut. 
		double P3Cut = .3; // 0.35 are the new 2021 recommeneded

                
		// zooming in on 4 specific regions in the following PzPt graphs
		double Reg1Xmin = 0.6; double Reg1Xmax = 1.1; double Reg1Ymin = 0.05; double Reg1Ymax = 0.3;
		double Reg2Xmin = 0.3; double Reg2Xmax = 0.8; double Reg2Ymin = 0.35; double Reg2Ymax = 0.65;
		double Reg3Xmin = 0.1; double Reg3Xmax = 0.5; double Reg3Ymin = 0.; double Reg3Ymax = 0.35;
		double Reg4Xmin = 0.6; double Reg4Xmax = 1.5; double Reg4Ymin = 0.2; double Reg4Ymax = 0.65;
		////////////// Photons //////////////////
		double ECut = 0.1;
                double ECuts[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
                bool allGeneralCutsPassed_pi0FCAL_ECuts[10];
		double ECutFCAL = 0.1;
		double thetaCutMin = 2.5; double thetaCutMax1 = 10.3; double thetaCutMax2 = 11.9;//11.5;
		double dijCut = 12.5;
		///////////// General ///////////////////
		double unusedEnergyCut = 0.1;
		double MMsqCut = 0.05;
		double ChiSqCut = 20.0;
		double originalChiSqCut = 13.277;
		double chiSq100 = 100;
		double RFCut = 0.5*4; // 4ns is the beam period.
		double beamECut = 6;
		double etaProtonBaryonCut;
		double pi0ProtonBaryonCut;
		// mUEChiSq
		double ellipseX; double ellipseY; double ellipseXr; double ellipseYr; 
    		double ellipseXBS1; double ellipseYBS1; double ellipseXrBS1; double ellipseYrBS1;
    		double ellipseXBS2; double ellipseYBS2; double ellipseXrBS2; double ellipseYrBS2;
		double skipX, skipY;
		double ellipseXr_loose, ellipseYr_loose;
		double weightBS=1;
		double weightBSpi0=1;
		double weightBSeta=1;
		double weightB=1;
		double areaRatio=1;
		
		// Beam cuts
		bool pBeamAsymE=true;
		bool pBeamE30to46=true;
		bool pBeamE46to62=true;
		bool pBeamE62to78=true;
		bool pBeamE78to94=true;
		bool pBeamE94to11=true;
		bool pBeamE8GeVPlus=true;
		// For the Deck beam asymmetry
		bool pBeamE82to88=true;
		bool pMpi0etaDoubleRegge=true;

		
		// pi0Eta specifc cuts
		bool pEtaProtonBaryonCut=true;
		bool ppi0ProtonBaryonCut=true;
		bool pBeamE=true;
		bool p_phiMassBinned[numBinsMass]; 
		double iLowMass;
		double iUpMass;
		bool p_tMassBinned[numBinsMass_t]; 
		double iLowMass_t;
		double iUpMass_t;

		bool p_massTBinned[10];

		// SOME VARIABLES FOR PURITY STUDY BY BINNING ON CHISQ AND MPI0ETA. WILL USE TO PLOT META AND MPI0
		double minMpi0eta = 0.6;	
		double maxMpi0eta = 2.2;
		static const int numMpi0etaBins = 8;
		double binStepMpi0eta = (maxMpi0eta-minMpi0eta)/numMpi0etaBins;
		static const int numRegions_UE=10;
		static const int numRegions_ChiSq=7;
		bool p_pi0MassEtaMassUEregion[numRegions_UE]; 
		bool p_pi0MassEtaMassChiSqregion[numRegions_ChiSq*numMpi0etaBins]; 
		double iUpUE;
		double iUpChiSq;
		
		// FOR THE pi0 MASS BINNED IN E(pi0) selecting out the f2(1270) region
		static const int numRegions_E=17;
		bool p_pi0MassPi0Eregion_1[numRegions_E]; 
		bool p_pi0MassPi0Eregion_2[numRegions_E]; 
		double iLowE[numRegions_E]  = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5};
		double iUpE[numRegions_E] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5,9.5};
		bool pSelectf2=true;
		// General Cuts
		bool pVanHove=true;
		bool mBeamE=true;
		bool pUnusedEnergy=true;
		bool pChiSq=true;
		bool pLooseChiSq=true;
		bool pLooseUnusedEnergy=true;
		bool pDeltaTRF=true ;
		bool pMissingMassSquared=true;
		
		bool pUEregion0=true;
		bool pUEregion1=true;
		bool pUEregion2=true ;
		bool pUEregion3=true;
		bool pUEregion4=true;
		bool pUEregion5=true;
		bool pUEregion6=true;
		bool pUEregion7=true;
		bool pUEregion8=true;
		bool pUEregion9=true;
		// NeuUEregiontral Cuts
		bool outsideEllipse_loose=true;
		bool pinsideEllipse_loose=true;
		bool outsideEllipse=true;
		bool pinsideEllipse=true;
		bool outsideEllipse_13_24=true;
		bool outsideEllipse_14_23=true;
		bool outsideEllipseBS1=true;
		bool pinsideEllipseBS1=true;
		bool outsideEllipseBS2=true;
		bool pinsideEllipseBS2=true;
                bool omegaCut=true;
		
		bool inBox[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};
		bool inBox_noOtherCuts[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};
		TLorentzVector pi0Momentum;
		TLorentzVector etaMomentum;
		double Mpi0eta_4;
		double Mpi0eta_10;
		double Mpi0eta_11;
		double Mpi0eta_12;

		bool pYellowBKG=true;
		bool pdij3pass=true;
		bool pPhoton1E=true;
		bool pPhoton2E=true;
		bool pPhoton3E=true;
		bool pPhoton4E=true;
		bool pPhotonE=true ;
		bool pPhoton1Theta=true;
		bool pPhoton2Theta=true;
		bool pPhoton3Theta=true;
		bool pPhoton4Theta=true;
		bool pPhotonTheta=true;
		
		// Charged Cuts
		bool pMagP3Proton=true;
		bool pzCutmin=true;
		bool pRProton=true;
		bool pdEdxCDCProton=true;
		bool pReg1=true;
		bool pReg2=true;
		bool pReg3=true;
		bool pReg4=true;
		
		//showerQuality
		bool pShowerQuality0 = true;
		bool pShowerQuality1 = true;
		bool pShowerQuality2 = true;
		bool pShowerQuality3 = true;
		bool pShowerQuality = true;


		// locWherePhoton will be set equal to photonDetectedSyss[N] 
		bool pPhotonInBCALorFCAL[4];
		bool pPhotonInFCAL[4];
		bool pPhotonInBCAL[4];
		bool pPi0InFCAL=true;
		bool pPi0InBCAL=true;
		bool pPi0InSplit=true;
		bool pEtaInFCAL=true; 
		bool pEtaInBCAL=true; 
		bool pEtaInSplit=true;
		bool pPi0InFCAL_mismatch=true;
		bool pPi0InBCAL_mismatch=true;
		bool pPi0InSplit_mismatch=true;
		bool pEtaInFCAL_mismatch=true; 
		bool pEtaInBCAL_mismatch=true; 
		bool pEtaInSplit_mismatch=true;
                bool allGeneralCutsPassed_thresh_pi0BCAL=true;
                bool allGeneralCutsPassed_thresh_pi0FCAL=true;
                bool allGeneralCutsPassed_thresh_pi0FCAL_bigThetaCut=true;
                bool allGeneralCutsPassed_thresh_pi0SPLIT=true;


		bool detectorCut=true;

		bool ptpLT1=true;
		// Times it passes a cut
		int count_events=0;
		int count_totEvents=0;
		int count_seanResTest=0;
		int count_combos=0;
		int count_MPi0P14=0;	
		int count_correctTopology=0;
		int count_ShowerQuality=0;
		int count_BeamE8GeVPlus=0;
		int count_UnusedEnergy=0;
		int count_ChiSq=0;
		int count_DeltaTRF=0;
		int count_dij3pass=0;
		int count_PhotonE=0;
		int count_PhotonTheta=0;
		int count_MagP3Proton=0;
		int count_zCutmin=0;
		int count_RProton=0;
		int count_MissingMassSquared=0;
		int count_dEdxCDCProton=0;
		int count_insideEllipse=0;
		int count_allGeneralCutsPassed=0;
		int count_allGeneralCutsPassedPlusTracked=0;
		
		// location cuts
		bool inBCAL=true;
		bool inTOF=true;
		bool inFCAL=true;
		bool inSTART=true;
		bool inSYS_NULL=true;
		
		bool medBool[5]={true,true,true,true,true};

		// baryon cuts
		bool mMandelstamT_mdelta=true;
		bool mMandelstamT_mdelta_petaProton=true;
		bool mMandelstamT_mdelta_pVanHove=true;
		bool mDelta=true;

		
		// Various combinations of cuts, the majority of them will be used just for a few histograms like when showing unused energy graph we will use mUE which
		// removes the UE cut from allGeneralCutsPassed. m prefix basically stands for minus
		bool allGeneralCutsPassed=true;
		bool allGeneralCutsPassed_tpLT05=true;
		bool allGeneralCutsSinglePolarization=true;
		
		// ------------------------------------------------------------------
		// ***************** THE ORIGINAL BA MEASUREMENT ********************
		// ------------------------------------------------------------------
		static const int numMpi0eta = 5; // there will be this many lower thesholds to calculate asymmetries vs teta/tpi0 with
		static const int numTBins = 5;
		static const int numMpi0etaRes = 9; // there will be this many bins in Mpi0eta to find the asymmetries in. Will extend down into the resonance region
		double lowerMpi0eta[numMpi0etaRes] = {0.9, 1.060, 1.24, 1.4, 1.65, 1.9, 2.15, 2.4, 2.65};
		double upperMpi0eta[numMpi0etaRes] = {1.060, 1.24, 1.4, 1.65, 1.9, 2.15, 2.4, 2.65, 2.9};
                bool ptEtaUpThreshold[4];  // Checking Mpi0eta yield asymmetry varying the upper threshold of teta/tpi.
                bool ptPi0UpThreshold[4];
                bool ptEtaBin[numMpi0etaRes]; // Plotting Mpi0eta in bins of teta/tpi. Can see contributions to the asymmetry as a function of Mpi0eta in bins of teta/tpi
                bool ptPi0Bin[numMpi0etaRes];
		bool ptEtaBeamAsym_000[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_045[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_090[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_135[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_AMO[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_backwardPi0P_000[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_backwardPi0P_045[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_backwardPi0P_090[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_backwardPi0P_135[numMpi0eta*numTBins];
		bool ptEtaBeamAsym_backwardPi0P_AMO[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_000[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_045[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_090[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_135[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_AMO[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_backwardEtaP_000[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_backwardEtaP_045[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_backwardEtaP_090[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_backwardEtaP_135[numMpi0eta*numTBins];
		bool ptPi0BeamAsym_backwardEtaP_AMO[numMpi0eta*numTBins];

		// We can also look at BA in bins of M(pi0eta) for fast pion and for fast eta
                // 6 bins in t for 5 bins between 0 and 1, and t<1 is the 6th
                // numMpi0etaRes has 9 bins which extend down into the resonance region
                static const int numt1BinRes=6;
		bool pMpi0etaBeamAsym_000_fastEtaBinned[numt1BinRes*numMpi0etaRes]; // WOULD SELECT ON FAST ETA
		bool pMpi0etaBeamAsym_045_fastEtaBinned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_090_fastEtaBinned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_135_fastEtaBinned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_AMO_fastEtaBinned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_000_fastPi0Binned[numt1BinRes*numMpi0etaRes]; // WOULD SELECT ON FAST PI0
		bool pMpi0etaBeamAsym_045_fastPi0Binned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_090_fastPi0Binned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_135_fastPi0Binned[numt1BinRes*numMpi0etaRes];
		bool pMpi0etaBeamAsym_AMO_fastPi0Binned[numt1BinRes*numMpi0etaRes];

		bool pMpi0etaBeamAsym_000_fastEta[numMpi0etaRes]; // WOULD SELECT ON FAST ETA
		bool pMpi0etaBeamAsym_045_fastEta[numMpi0etaRes];
		bool pMpi0etaBeamAsym_090_fastEta[numMpi0etaRes];
		bool pMpi0etaBeamAsym_135_fastEta[numMpi0etaRes];
		bool pMpi0etaBeamAsym_AMO_fastEta[numMpi0etaRes];
		bool pMpi0etaBeamAsym_000_fastPi0[numMpi0etaRes]; // WOULD SELECT ON FAST PI0
		bool pMpi0etaBeamAsym_045_fastPi0[numMpi0etaRes];
		bool pMpi0etaBeamAsym_090_fastPi0[numMpi0etaRes];
		bool pMpi0etaBeamAsym_135_fastPi0[numMpi0etaRes];
		bool pMpi0etaBeamAsym_AMO_fastPi0[numMpi0etaRes];

                // Will just make 4 bins of u3 from 0-2. So 0.5 GeV^2 bins
		bool ptrecoilBeamAsym_000_fastEta[4]; // WOULD SELECT ON FAST ETA
		bool ptrecoilBeamAsym_045_fastEta[4];
		bool ptrecoilBeamAsym_090_fastEta[4];
		bool ptrecoilBeamAsym_135_fastEta[4];
		bool ptrecoilBeamAsym_AMO_fastEta[4];
		bool ptrecoilBeamAsym_000_fastPi0[4]; // WOULD SELECT ON FAST PI0
		bool ptrecoilBeamAsym_045_fastPi0[4];
		bool ptrecoilBeamAsym_090_fastPi0[4];
		bool ptrecoilBeamAsym_135_fastPi0[4];
		bool ptrecoilBeamAsym_AMO_fastPi0[4];

		// We can also look at BA in bins of t_recoil for fast pion and for fast eta
		bool ptrecoilBeamAsym_000_fastEtaBinned[5*4]; // WOULD SELECT ON FAST ETA
		bool ptrecoilBeamAsym_045_fastEtaBinned[5*4];
		bool ptrecoilBeamAsym_090_fastEtaBinned[5*4];
		bool ptrecoilBeamAsym_135_fastEtaBinned[5*4];
		bool ptrecoilBeamAsym_AMO_fastEtaBinned[5*4];
		bool ptrecoilBeamAsym_000_fastPi0Binned[5*4]; // WOULD SELECT ON FAST PI0
		bool ptrecoilBeamAsym_045_fastPi0Binned[5*4];
		bool ptrecoilBeamAsym_090_fastPi0Binned[5*4];
		bool ptrecoilBeamAsym_135_fastPi0Binned[5*4];
		bool ptrecoilBeamAsym_AMO_fastPi0Binned[5*4];
		
                // Apply tighter cuts to select more strongly the teta/tpi asymmetries
		bool ptightBeamAsym_000_fastEta[5]; // WOULD SELECT ON FAST ETA
		bool ptightBeamAsym_045_fastEta[5];
		bool ptightBeamAsym_090_fastEta[5];
		bool ptightBeamAsym_135_fastEta[5];
		bool ptightBeamAsym_AMO_fastEta[5];
		bool ptightBeamAsym_000_fastPi0[5]; // WOULD SELECT ON FAST PI0
		bool ptightBeamAsym_045_fastPi0[5];
		bool ptightBeamAsym_090_fastPi0[5];
		bool ptightBeamAsym_135_fastPi0[5];
		bool ptightBeamAsym_AMO_fastPi0[5];

		// THE GOAL HERE IS TO LOOK AT THE BA CONTRIBUTIONS IN BINS OF M(PI0P) AND IN BINS OF M(ETAP). IN M(PI0P) WE SEE 3 N* RESONANCES. IF WE BIN IN THIS VARIABLE
		// WE CAN DETERMINE WHAT THE CONTRIBUTIONS ARE TO THE TOTAL ASYMMETRY WE ARE SEEING
		static const int numBaryonBins=5;
	        double minPi0P[numBaryonBins] = {1.15, 1.4, 1.6, 1.75, 1.95};
	        double maxPi0P[numBaryonBins] = {1.4, 1.6, 1.75, 1.95, 2.4};
	        double minEtaP[numBaryonBins] = {1.5, 1.65, 1.9, 2.2, 2.4};
	        double maxEtaP[numBaryonBins] = {1.65, 1.9, 2.2, 2.4, 2.8};
	        //double minPi0P[3] = {1.1, 1.5, 1.9};
	        //double maxPi0P[3] = {1.5, 1.9, 2.3};
	        //double minEtaP[3] = {1.5, 1.9, 2.3};
	        //double maxEtaP[3] = {1.9, 2.3, 2.7};
		bool pMpi0pBeamAsym_000[numBaryonBins*numTBins]; // WOULD SELECT ON FAST ETA
		bool pMpi0pBeamAsym_045[numBaryonBins*numTBins];
		bool pMpi0pBeamAsym_090[numBaryonBins*numTBins];
		bool pMpi0pBeamAsym_135[numBaryonBins*numTBins];
		bool pMpi0pBeamAsym_AMO[numBaryonBins*numTBins];
		bool pMetapBeamAsym_000[numBaryonBins*numTBins]; // WOULD SELECT ON FAST PI0
		bool pMetapBeamAsym_045[numBaryonBins*numTBins];
		bool pMetapBeamAsym_090[numBaryonBins*numTBins];
		bool pMetapBeamAsym_135[numBaryonBins*numTBins];
		bool pMetapBeamAsym_AMO[numBaryonBins*numTBins];
		bool pMpi0pFastEta[numBaryonBins];
		bool pMetapFastPi0[numBaryonBins];
                bool beamAsymVeryFastEta=true;
		

		bool mMandelstamT=true;
		bool mMandelstamT_mBeamE8GeVPlus=true;
		//bool pDiffCL=true; 
                bool mP3Proton=true;
                bool mFCALShower=true;
		bool pDiffUE=true; 
		bool mRProton=true;
		bool mRProtonZMin = true; 
		bool mdEdxCDC = true;
		bool mZMin = true;
		bool mMagP3 = true;
		bool mPhotonE = true;
		bool mPhotonTheta = true;
		bool mdij3 = true;
		bool mUE = true;
		bool mUEChiSq = true;
		bool mChiSq = true;
		bool mMMSq = true;

                // Mpi0p alternative using VH
		bool pVH = true;
		bool mMPi0P14_VH = true;
		// mEllipseRY contains both the red and yellow regions only.
		bool mEllipse = true;
		bool mEllipse_pre = true;
		bool mEllipse_pre_tGT01LT03 = true;
		bool mEllipse_pre_tGT03LT06 = true;
		bool mEllipse_pre_tGT06LT10 = true;
		bool mEllipse_pre_tAll = true;
		bool mEllipse_pre_tAll_EbeamAll = true;
                bool mEllipse_pre_tpLT05=true;
		bool ptGT1=true;
		bool ptLT05=true;
                bool ptpLT05=true;
		bool ptGT05LT1=true;
		bool mEllipse_pre_tGT1 = true;
		bool mEllipse_pre_tGT05LT1 = true;
		bool mEllipse_pre_tLT05 = true;
		bool mEllipse_pre_tAll_delta = true;
		bool mEllipseUE = true;
		bool mEllipseUE_pre = true;
		bool mEllipseUEChiSq = true;
		bool mEllipseUEChiSq_pre = true;
		bool mEllipseLooseUEChiSq_pre = true;
		bool mEllipseChiSq = true;
		bool mEllipseChiSq_pre = true;
		bool pMPi0P14=true;
		bool mMPi0P14=true;
		bool baseCuts=true;
		bool baseCuts_mChiUE=true;

		ofstream outputDeckBA;
		bool baseAsymCut_mDelta=true;
		bool baseAsymCut_fastEta=true;
		bool baseAsymCut_fastPi0=true;
		bool baseAsymCut_res=true;
		bool baseAsymCut_mDelta_fastEta=true;
		bool baseAsymCut_mDelta_fastPi0=true;
		bool baseAsymCut=true;
		bool baseAsymCut_mMpi0etaDoubleRegge=true;
		bool baseAsymCut_backwardPi0P=true;
		bool baseAsymCut_backwardPi0P_mDelta=true;
		bool baseAsymCut_backwardEtaP=true;
		bool baseAsymComparisionRegge=true;
		bool baseAsymComparisionDelta=true;
		bool looseCutsUEChiSq=true;
                bool kinematicSelected_looseCutsUEChiSq=true;
                bool combinatoricStudy=true;
		bool pBase_pT_pIE_pBE8288_pMPi0P_pDelta = true;
                bool noCut=true;

                bool dzRP=true;
                bool dzR=true;
		bool mEllipse_pre_pi0BCAL = true;
		bool mEllipse_pre_pi0FCAL = true;
		bool mEllipse_pre_pi0FCAL_thresh = true;
		bool mEllipse_pre_etaFCAL_thresh = true;
		bool mEllipse_pre_pi0SPLIT = true;
		bool mEllipse_pre_etaBCAL = true;
		bool mEllipse_pre_etaFCAL = true;
		bool mEllipse_pre_etaSPLIT = true;


                static const int nThetaBins=5;
                double thetaBinWidth=40/nThetaBins;
                //double etaBCAL_thetaBins[nThetaBins+1]={0,11,180};
                //double etaFCAL_thetaBins[nThetaBins+1]={0,3,180};
                //double etaSPLIT_thetaBins[nThetaBins+1]={0,6,180};
                //double pi0BCAL_thetaBins[nThetaBins+1]={0,17,180};
                //double pi0FCAL_thetaBins[nThetaBins+1]={0,4,180};
                //double pi0SPLIT_thetaBins[nThetaBins+1]={0,10,180};
		bool mEllipse_pre_pi0BCAL_thetaBin[nThetaBins];
		bool mEllipse_pre_pi0FCAL_thetaBin[nThetaBins];
		bool mEllipse_pre_pi0SPLIT_thetaBin[nThetaBins];
		bool mEllipse_pre_etaBCAL_thetaBin[nThetaBins];
		bool mEllipse_pre_etaFCAL_thetaBin[nThetaBins];
		bool mEllipse_pre_etaSPLIT_thetaBin[nThetaBins];




                bool allGen_barybkg=true;
                bool allGen_vanHove=true;
		bool cutsToApply = true; 

                Int_t whichSignalRegion;
                bool beamPhotonMatchToThrown;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////// ********* DEFINING ALL HISTOGRAMS ****************/////////////////////////////////////////////
		
		std::string pi0BinRange[3] = {"200", "0.05","0.25"} ;
		std::string etaBinRange[3] = {"300", "0.25","0.85"} ;
		int id;
       		int id_noCutDij3=0;
        	std::vector<std::vector<int>> vec_group_ids;
		std::string groupNames[16] = {"beam", "p1", "ph12b1", "ph1234", "ph12", "ph34", "ph34b1", "entire combo", "ph12p1", "ph34p1", "ph1234p1", "FCAL Pairs", "phN", "BCAL Pairs", "ph13", "ph24"};
        	std::vector<int> group_ids;
        	int groupVec_it; // this isnt used here but will be used later to denote the graphs we are filling.

		ofstream compositionFile;

	ClassDef(DSelector_ver20, 0);
};

void DSelector_ver20::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
        dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
        dDecayingEtaWrapper = dStep2Wrapper->Get_InitialParticle();
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_ver20_h



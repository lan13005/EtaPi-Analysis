#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"
#include "TFile.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "EtaPiPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderFilter.h"
#include "AMPTOOLS_AMPS/TwoPSHelicity.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/KMatrixPVectorAmp.h"

typedef EtaPiPlotGenerator PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( TwoPSHelicity() );
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( KMatrixPVectorAmp() );
  AmpToolsInterface::registerAmplitude( Piecewise() );
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderFilter() );
}

using namespace std;

int main( int argc, char* argv[] ){
  
  
      // ************************
      // usage
      // ************************
  
    auto help_message = []()
    {
        cout << "======\nUsage:" << endl << "======" << endl;
        cout << "etapi_plotter <results file name> "  << 
            "-o <output file name> " <<
            "-s <amp string> " << 
            "-a <doAccCorr | T/F> " << 
            "-F <L waves to gather fit fractions for> " <<
            "-var <plot all variables or not | T/F> " << 
            "-gen <plot gen histograms or not T/F>" << 
            endl << endl;
        cout << "second allowed form: <file> -g show GUI" << endl;
        exit(1);
    };
    
    if ((argc < 9)*(string(argv[2])!="-g")){
      help_message();
      return 0;
    }

    bool showGui = false;
    string outName = "etapi_plot.root";
    string resultsName(argv[1]);
    bool keepAllAmps=true; // always resort to plotting all the amps
    string ampString;
    string sdoAccCorr; // sdoAccCorr can actually be anything but will default to false if not set to "T"
    string Ls;
    string splotAllVars;
    bool plotAllVars=true;
    string splotGenData;
    bool plotGenData=true;
    bool doAccCorr=true; // setting some default that will always get overwritten - just to get rid of the compiler warnings 
    for (int i = 2; i < argc; i++){
  
      string arg(argv[i]);
  
      if (arg == "-g"){
        showGui = true;
      }
      if (arg == "-h"){
          help_message();
        exit(1);
      }
      if (arg == "-s"){
          ampString=argv[++i];
          if (ampString==""){ 
              keepAllAmps=true;
          }
          else {
              keepAllAmps=false;
          }
      }
      if (arg == "-a"){
          sdoAccCorr=argv[++i];
          if (sdoAccCorr=="true"){ 
              doAccCorr=true;
              cout << "Doing acceptance correction!" << endl;
          }
          else{ 
              doAccCorr=false;
              cout << "Not doing acceptance correction!" << endl;
          }
      }
      if (arg == "-F"){
          Ls=argv[++i];
      }
      if (arg == "-var"){
          splotAllVars=argv[++i];
          if (splotAllVars=="true"){ 
              plotAllVars=true;
              cout << "Making plots of all variables!" << endl;
          }
          else{ 
              plotAllVars=false;
              cout << "Only plotting mass plots!" << endl;
          }
      }
      if (arg == "-gen"){
          splotGenData=argv[++i];
          if (splotGenData=="true"){ 
              plotGenData=true;
              cout << "Making gen plots also!" << endl;
          }
          else{ 
              plotGenData=false;
          }
      }
  
    }
    // ************************
    // parse the command line parameters
    // ************************
    cout << "Fit results file name    = " << resultsName << endl;

    // ************************
    // load the results and display the configuration info
    // ************************
    FitResults results( resultsName );
    if( !results.valid() ){
      cout << "Invalid fit results in file:  " << resultsName << endl;
      exit( 1 );
    }

    // ************************
    // set up the plot generator
    // ************************
    atiSetup();
    cout << "Loading FitResults into PlotGen" << endl;
    PlotGenerator::Option opt;
    if (plotGenData)
        opt=PlotGenerator::kDefault;
    else
        opt=PlotGenerator::kNoGenMC;
    PlotGen plotGen( results, opt ); // second argument is defined as: enum Option { kDefault = 0, kNoGenMC }; enum automatically +1 for has no initializer

    // ************************
    // start the GUI
    // ************************
    if(showGui) {
          cout << ">> Plot generator ready, starting GUI..." << endl;
          
          int dummy_argc = 0;
          char* dummy_argv[] = {};  
          TApplication app( "app", &dummy_argc, dummy_argv );
          
          gStyle->SetFillColor(10);
          gStyle->SetCanvasColor(10);
          gStyle->SetPadColor(10);
          gStyle->SetFillStyle(1001);
          gStyle->SetPalette(1);
          gStyle->SetFrameFillColor(10);
          gStyle->SetFrameFillStyle(1001);
          
          PlotFactory factory( plotGen );	
          PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
          app.Run();
          return 0;
    }
    
    // ************************
    // Saving plots
    // ************************
    if (keepAllAmps)
        cout << "\n------------\nPlotting all amplitudes\n------------" << endl;

    string waveset; 
    string wave;
    stringstream ss_wavesets(ampString);
    while(getline(ss_wavesets, waveset, ';')){
        cout << "waveset: " << waveset << endl;
        map<string,int> selectedAmps;
        stringstream ss_waves(waveset);
        cout << "\n------------\nPlotting only with these amplitude contributions:" << endl;
        while(getline(ss_waves, wave, '_')){
            selectedAmps[wave] = -1; // Just building map. -1 is placeholder
            cout << wave << endl;
        }
        cout << "------------\n" << endl;

        outName = "etapi_plot_"+waveset+".root";
        cout << "Output file name    = " << outName << endl << endl;

        // ************************
        // set up an output ROOT file to store histograms
        // ************************

        cout << "Creating output file" << endl;
        TFile* plotfile = new TFile( outName.c_str(), "recreate");
        TH1::AddDirectory(kFALSE);
  
        vector<string> reactionList = results.reactionList();
        for (int ireact=0; ireact<(int)reactionList.size(); ++ireact){
            plotGen.enableReaction( results.reactionList()[ireact] );
        }
        vector<string> sums = plotGen.uniqueSums();
        vector<string> amps = plotGen.uniqueAmplitudes();
        cout << "\npossible sums:\n----------" << endl;
        for (auto i:sums){
          cout << i << endl;
        }
        cout << "\npossible amps:\n----------" << endl;
        for (unsigned int i=0; i<amps.size(); ++i){
          cout << amps[i];
          if (selectedAmps.find(amps[i])!=selectedAmps.end()){
              selectedAmps[amps[i]]=i;
              cout << " saving this one";
          }
          cout << endl;
        }
  
        cout << "Determining amplitues/sums to use" << endl;
        if(!keepAllAmps){
          for (unsigned int i=0; i<amps.size(); ++i){
            plotGen.disableAmp(i);
          }
          for (std::map<string,int>::iterator it=selectedAmps.begin(); it!=selectedAmps.end(); ++it){
            plotGen.enableAmp(it->second);
          }
        }
  
        // loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
        for (unsigned int isum = 0; isum <= sums.size(); isum++){
  
          // turn on all sums by default
          for (unsigned int i = 0; i < sums.size(); i++){
            plotGen.enableSum(i);
          }
  
          // for individual contributions turn off all sums but the one of interest
          if (isum < sums.size()){
            for (unsigned int i = 0; i < sums.size(); i++){
              if (i != isum) plotGen.disableSum(i);
            }
          }
  
          // loop over data, accMC, and genMC
//          cout << "Looping of datasets" << endl;
          for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
            if (isum < sums.size() && iplot == PlotGenerator::kData) continue; // only plot data once
  
            // loop over different variables
//            cout << "   looping over requested histograms" << endl;
            for (unsigned int ivar  = 0; ivar  < EtaPiPlotGenerator::kNumHists; ivar++){
                for (int ireact=0; ireact<(int)reactionList.size(); ++ireact){ // list of reactions: i.e. EtaPi0_000, EtaPi0_045...
                    string reactionName = reactionList[ireact];
                    string histname =  reactionName+"_";

                    if (!plotAllVars){
                        if (
                                (ivar == EtaPiPlotGenerator::kEtaCosTheta) ||
                                (ivar == EtaPiPlotGenerator::kPhi) || 
                                (ivar == EtaPiPlotGenerator::kphi) || 
                                (ivar == EtaPiPlotGenerator::kPsi) || 
                                (ivar == EtaPiPlotGenerator::kt) 
                           )
                           continue; 
                    }
                    if ( (!plotGenData)&&((iplot == PlotGenerator::kGenMC)) ) {
                        continue;
                    }

                    // set unique histogram name for each plot (could put in directories...)
                    if (ivar == EtaPiPlotGenerator::kEtaPiMass)  histname += "Metapi";
                    else if (ivar == EtaPiPlotGenerator::kEtaPiMass_40MeVBin)  histname += "Metapi_40MeVBin";

                    // include these above in !plotAllVars condition if we do not wish to always draw them
                    else if (ivar == EtaPiPlotGenerator::kEtaProtonMass)  histname += "Metaproton";
                    else if (ivar == EtaPiPlotGenerator::kPi0ProtonMass)  histname += "Mpi0proton";

                    else if (ivar == EtaPiPlotGenerator::kEtaCosTheta)  histname += "cosTheta";
                    else if (ivar == EtaPiPlotGenerator::kPhi)  histname += "Phi";
                    else if (ivar == EtaPiPlotGenerator::kphi)  histname += "phi";
                    else if (ivar == EtaPiPlotGenerator::kPsi)  histname += "psi";
                    else if (ivar == EtaPiPlotGenerator::kt)  histname += "t";
                    else continue;
  
                    if (iplot == PlotGenerator::kData) histname += "dat";
                    if (iplot == PlotGenerator::kBkgnd) histname += "bkg";
                    if (iplot == PlotGenerator::kAccMC) histname += "acc";
                    if (iplot == PlotGenerator::kGenMC) histname += "gen";
  
                    if (isum < sums.size()){
                      //ostringstream sdig;  sdig << (isum + 1);
                      //histname += sdig.str();
  
                    // get name of sum for naming histogram
                      string sumName = sums[isum];
                      histname += "_";
                      histname += sumName;
                    }
  
//                    cout << "      projecting histogram..." << endl;
                    Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
//                    cout << "      converting hists to ROOT format..." << endl;
                    TH1* thist = hist->toRoot();
                    thist->SetName(histname.c_str());
                    plotfile->cd();
//                    cout << "      writing..." << endl;
                    thist->Write();
                }
            }
          }
        }
        plotfile->Close();
    }
    
    ofstream outfile;
    outfile.open( "etapi_plotter_output.log" );
    outfile << "TOTAL EVENTS = " << results.intensity(doAccCorr).first << " +- " << results.intensity(doAccCorr).second << endl;
    vector<string> fullamps = plotGen.fullAmplitudes();
    map<string,vector<string>> uniqueAmps;
    for (unsigned int i = 0; i < fullamps.size(); i++){
      string amp = fullamps[i].substr(fullamps[i].find_last_of("::")+1,fullamps[i].npos);
      if (uniqueAmps.find(amp) == uniqueAmps.end())
        uniqueAmps[amp]=vector<string>{fullamps[i]};
      else
        uniqueAmps[amp].push_back(fullamps[i]);
      vector<string> useamp;  useamp.push_back(fullamps[i]);
      // THIS WOULD OUTPUT ALL THE AMPLITUDES INCLUDING NegativeRe, NegativeIm 
      //   all the pols, every intialized amplitude that is defined in the config file 
      outfile << "FIT FRACTION " << fullamps[i] << " = "
           << results.intensity(useamp, doAccCorr).first /
              results.intensity(doAccCorr).first <<  " +- "
           << results.intensity(useamp, doAccCorr).second /
              results.intensity(doAccCorr).first <<  endl;
    }
    // SUM OVER POLARIZATION AND REAL AND IMAGINARY PARTS
    for (auto ele: uniqueAmps){
      outfile << "FIT FRACTION " << ele.first << " = "
           << results.intensity(ele.second, doAccCorr).first /
              results.intensity(doAccCorr).first <<  " +- "
           << results.intensity(ele.second, doAccCorr).second /
              results.intensity(doAccCorr).first <<  endl;
    }


    // load the Ls and split along underscores to get a vector
    //   the first couple of letters is the amplitude name. i.e. D2++ and pD2+- will belong to groups D and pD
    vector<string> vLs; // vector for Ls
    stringstream ss_Ls(Ls);
    cout << "\n------------\nPlotting only with these amplitude contributions:" << endl;
    while(getline(ss_Ls, wave, '_')){
        vLs.push_back(wave);
    }
    map<string,vector<string>> merged;
    for(auto L: vLs){
        for(auto ele: uniqueAmps){
            int sizeL=(int)L.size();
            if (ele.first.substr(0,sizeL)==L){
                if (merged.find(L)!=merged.end())
                    merged[L].insert(merged[L].end(),uniqueAmps[ele.first].begin(),uniqueAmps[ele.first].end());
                else
                    merged[L]=uniqueAmps[ele.first];
            }
        }
    } 
    for (auto ele: merged){
      outfile << "FIT FRACTION " << ele.first << " = "
           << results.intensity(ele.second, doAccCorr).first /
              results.intensity(doAccCorr).first <<  " +- "
           << results.intensity(ele.second, doAccCorr).second /
              results.intensity(doAccCorr).first <<  endl;
    }

//    for (auto ele: merged["S"])
//        cout << ele << endl;
//    cout << endl;
//    for (auto ele: merged["D"])
//        cout << ele << endl;
//    cout << endl;
//    for (auto ele: merged["pD"])
//        cout << ele << endl;
//    cout << endl;
    




  
          // ************************
          // retrieve amplitudes for output
          // ************************
        /*
        // parameters to check
        vector< string > pars;
        pars.push_back("Pi+Pi-::Positive::S0+_re");
        pars.push_back("Pi+Pi-::Positive::S0+_im");
  
        // file for writing parameters (later switch to putting in ROOT file)
        ofstream outfile;
        outfile.open( "twopi_fitPars.txt" );
  
        for(unsigned int i = 0; i<pars.size(); i++) {
          double parValue = results.parValue( pars[i] );
          double parError = results.parError( pars[i] );
          outfile << parValue << "\t" << parError << "\t";
        }
  
        // Note: For twopi_amp_plotter: The following computations are nonsense for amplitudes
  
        // covariance matrix
        vector< vector< double > > covMatrix;
        covMatrix = results.errorMatrix();
  
        double SigmaN = results.parValue(pars[3]) + results.parValue(pars[6]);
        double SigmaN_err = covMatrix[5][5] + covMatrix[8][8] + 2*covMatrix[5][8];
  
        double SigmaD = 0.5*(1 - results.parValue(pars[0])) + results.parValue(pars[2]);
        double SigmaD_err = 0.5*0.5*covMatrix[2][2] + covMatrix[4][4] - 2*0.5*covMatrix[2][4];
  
        double Sigma = SigmaN/SigmaD;
        double Sigma_err = fabs(Sigma) * sqrt(SigmaN_err/SigmaN/SigmaN + SigmaD_err/SigmaD/SigmaD);
  
        double P = 2*results.parValue(pars[6]) - results.parValue(pars[4]);
        double P_err = sqrt(2*2*covMatrix[8][8] + covMatrix[6][6] - 2*2*covMatrix[6][8]);
  
        Sigma = Sigma_err = P = P_err = 0;
        outfile << Sigma << "\t" << Sigma_err << "\t";
        outfile << P << "\t" << P_err << "\t";
  
        outfile << endl;
        */
  
    return 0;

}


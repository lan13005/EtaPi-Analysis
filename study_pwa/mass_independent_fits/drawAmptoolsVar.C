#include <TVector3.h>

float PI=3.14159;
float histLineWidth=5;
float lowMass=1.04; //0.80;
float uppMass=1.56;
int nbins=13; //25; // THIS IS NOT TO BE CONFUSED WITH NM. THIS IS THE BINNING FOR THE MASS HISTOGRAM
int nbinsAng=30;
int nm=1; // NUMBER OF BINNINGS TO MAKE THE HISTOGRAMS IN
string foutTag="phase1";//"kmatrix"; 
bool doAccCorr=true; // SHOULD WE ACCEPTANCE CORRECT THE YIELDS?
//// THE FILES YOU WANT TO RUN OVER IS IN THE DRAWAMPTOOLSVAR FUNCTION. SINCE WE USE A FOR LOOP FOR 
//       INITIALIZATION WE CANNOT DO IN GLOBAL SCOPE

map<string,TH1*> getHists(string file, string tag){
    map<string,TH1*> hists;

    for (int m=0; m<nm; ++m){
        hists[("Mpi0eta_mBin"+to_string(m)).c_str()]=new TH1F( ("Mpi0eta_mBin"+to_string(m)+"-"+tag).c_str(), 
                "Invariant Mass of #eta #pi;M(#eta#pi)", nbins, lowMass, uppMass );
        hists[("cosThetaGJ_mBin"+to_string(m)).c_str()]=new TH1F(("cosTheta_mBin"+to_string(m)+"-"+tag).c_str(), "cos( #theta ) GJ;cos(#theta) GJ", 50, -1., 1. );
        hists[("cosThetaHel_mBin"+to_string(m)).c_str()]=new TH1F(("cosTheta_hel_mBin"+to_string(m)+"-"+tag).c_str(), "cos( #theta ) Hel;cos(#theta) Hel", 50, -1., 1. );
        hists[("phiGJ_mBin"+to_string(m)).c_str()]=new TH1F( ("phi_mBin"+to_string(m)+"-"+tag).c_str(), "#phi; #phi (rad)", 50, -1*PI, PI );
        hists[("phiHel_mBin"+to_string(m)).c_str()]=new TH1F( ("phi_hel_mBin"+to_string(m)+"-"+tag).c_str(), "#phi; #phi (rad)", 50, -1*PI, PI );
        hists[("t_mBin"+to_string(m)).c_str()]=new TH1F( ("t_mBin"+to_string(m)+"-"+tag).c_str(), "-t;-t (GeV^2)", 100, 0, 1 );
        hists[("ebeam_mBin"+to_string(m)).c_str()]=new TH1F( ("ebeam_mBin"+to_string(m)+"-"+tag).c_str(), "Ebeam;Ebeam (GeV^2)", 100, 4, 12 );
        hists[("Mpi0etaCosThetaGJ_mBin"+to_string(m)).c_str()]=new TH2F( ("Mpi0etaVsCosTheta_mBin"+to_string(m)+"-"+tag).c_str(), 
            ";M(#eta#pi);cos(#theta) GJ", nbins, lowMass, uppMass,  nbinsAng, -1., 1. );
        hists[("Mpi0etaCosThetaHel_mBin"+to_string(m)).c_str()]=new TH2F( ("Mpi0etaVsCosThetaHel_mBin"+to_string(m)+"-"+tag).c_str(), 
            ";M(#eta#pi);cos(#theta) Hel", nbins, lowMass, uppMass,  nbinsAng, -1., 1. );
        hists[("Mpi0etaPhiGJ_mBin"+to_string(m)).c_str()]=new TH2F( ("Mpi0etaVsPhi_mBin"+to_string(m)+"-"+tag).c_str(), 
            ";M(#eta#pi);#phi GJ", nbins, lowMass, uppMass,  nbinsAng, -1*PI, PI );
        hists[("Mpi0etaPhiHel_mBin"+to_string(m)).c_str()]=new TH2F( ("Mpi0etaVsPhiHel_mBin"+to_string(m)+"-"+tag).c_str(), 
            ";M(#eta#pi);#phi Hel", nbins, lowMass, uppMass,  nbinsAng, -1*PI, PI );
        hists[("cosThetaPhiGJ_mBin"+to_string(m)).c_str()]=new TH2F( ("cosThetaVsPhi_mBin"+to_string(m)+"-"+tag).c_str(), 
            ";cos(#theta) GJ;#phi GJ", nbinsAng, -1, 1,  nbinsAng, -1*PI, PI );
        hists[("cosThetaPhiHel_mBin"+to_string(m)).c_str()]=new TH2F( ("cosThetaVsPhiHel_mBin"+to_string(m)+"-"+tag).c_str(), 
            ";cos(#theta) Hel;#phi Hel", nbinsAng, -1, 1,  nbinsAng, -1*PI, PI );
    }

    for (auto const& x : hists){
        x.second->SetLineWidth(histLineWidth);
    }


    if (file!=""){
        TFile* data=TFile::Open(file.c_str());
        TTree* tree;
        data->GetObject("kin",tree);
        cout << "Loaded " << file << " with " << tree->GetEntries() << " entries in kin tree" << endl;

        // ******************************
        // LOAD VARS AND BRANCH ADDRESSSES
        // ******************************
        int particles;
        float particle_es[3];
        float particle_pxs[3];
        float particle_pys[3];
        float particle_pzs[3];
        float beam_e;
        float beam_px;
        float beam_py;
        float beam_pz;
        float weight;

        tree->SetBranchAddress("NumFinalState", &particles);
        tree->SetBranchAddress("E_FinalState", particle_es);
        tree->SetBranchAddress("Px_FinalState", particle_pxs);
        tree->SetBranchAddress("Py_FinalState", particle_pys);
        tree->SetBranchAddress("Pz_FinalState", particle_pzs);
        tree->SetBranchAddress("E_Beam", &beam_e);
        tree->SetBranchAddress("Px_Beam", &beam_px);
        tree->SetBranchAddress("Py_Beam", &beam_py);
        tree->SetBranchAddress("Pz_Beam", &beam_pz);
        tree->SetBranchAddress("Weight", &weight);

        TLorentzVector beam_p4;
        TLorentzVector target_p4(0,0,0,0.9382719);
        TLorentzVector recoil_p4;
        TLorentzVector eta_p4;
        TLorentzVector pi0_p4;

        // ******************************
        // READ IN DATA
        // ******************************
        Long64_t nentries=tree->GetEntries();
        cout << "There are nentries: " << nentries << endl;

        float mandelstam_t;
        float cosTheta_hel;
        float phi_hel;
        float cosTheta;
        float phi;

        for (auto i=0; i<nentries; ++i){
            tree->GetEntry(i);    
            beam_p4.SetPxPyPzE(beam_px,beam_py,beam_pz,beam_e);
            recoil_p4.SetPxPyPzE(particle_pxs[0],particle_pys[0],particle_pzs[0],particle_es[0]);
            pi0_p4.SetPxPyPzE(particle_pxs[1],particle_pys[1],particle_pzs[1],particle_es[1]);
            eta_p4.SetPxPyPzE(particle_pxs[2],particle_pys[2],particle_pzs[2],particle_es[2]);
            TLorentzVector resonance=pi0_p4+eta_p4;

            mandelstam_t = -(target_p4-recoil_p4).M2();

            TLorentzVector cm=beam_p4+target_p4;
            TLorentzRotation cmBoost( -cm.BoostVector() );
            TLorentzVector beam_cm=cmBoost*beam_p4;
            TLorentzVector recoil_cm=cmBoost*recoil_p4;
            TLorentzVector pi0_cm=cmBoost*pi0_p4;
            TLorentzVector eta_cm=cmBoost*eta_p4;
            TLorentzVector resonance_cm=pi0_cm+eta_cm;

            TLorentzRotation resonanceBoost( -resonance_cm.BoostVector() );
            TLorentzVector beam_res = resonanceBoost*beam_cm;
            TLorentzVector recoil_res = resonanceBoost*recoil_cm;
            TLorentzVector eta_res = resonanceBoost*eta_cm;
            TLorentzVector pi0_res = resonanceBoost*pi0_cm;
            
            /// ************************************************
            // This should be calculating the angles in the helictiy frame
            /// ************************************************

            // normal to the production plane
            // Actually -recoil_cm = resonance_cm. We use recoil since it is measured better
            TVector3 z_hel = -1. * recoil_res.Vect().Unit();
            TVector3 y_hel = (beam_cm.Vect().Unit().Cross(-recoil_cm.Vect().Unit())).Unit();
            TVector3 x_hel = y_hel.Cross(z_hel).Unit();
            TVector3 angles_hel((eta_res.Vect()).Dot(x_hel),
                                (eta_res.Vect()).Dot(y_hel),
                                (eta_res.Vect()).Dot(z_hel) );
            cosTheta_hel = angles_hel.CosTheta();
            phi_hel = angles_hel.Phi();

            
            /// ************************************************
            // This should be calculating the angles in the GJ frame
            /// ************************************************
            TVector3 z = beam_res.Vect().Unit();
            TVector3 y = y_hel; 
            TVector3 x = y.Cross(z).Unit();
            TVector3 angles(   (eta_res.Vect()).Dot(x),
                               (eta_res.Vect()).Dot(y),
                               (eta_res.Vect()).Dot(z) );
            cosTheta = angles.CosTheta();
            phi = angles.Phi();

            /// ************************************************
            // Calculate the rest of the variables
            /// ************************************************
            float pi0Mean=0.135881;
            float etaMean=0.548625;
            float pi0Std=0.0076;
            float etaStd=0.0191;
            float pi0Sig=3;
            float pi0Skip=1;
            float etaSig=3;
            float etaSkip=1;
            bool selectPi0=(pi0_p4.M()<(pi0Mean+pi0Std*pi0Sig))*(pi0_p4.M()>(pi0Mean-pi0Std*pi0Sig));
            bool selectEta=(eta_p4.M()<(etaMean+etaStd*etaSig))*(eta_p4.M()>(etaMean-etaStd*etaSig));
            bool select_a2=(resonance.M()>1.04)*(resonance.M()<1.56);
            bool select_t=(mandelstam_t<0.3);

            selectPi0=true;
            selectEta=true;
            select_a2=true;
            select_t=true;
            
            float width=(uppMass-lowMass)/nm;
            int m=(int)((resonance.M()-lowMass)/width);
            if ( (m<0)||(m>=nm) )
                continue;
            if (select_a2*select_t*selectEta*selectPi0){
                hists["Mpi0eta_mBin"+to_string(m)]->Fill(resonance.M(),weight);
                hists["cosThetaGJ_mBin"+to_string(m)]->Fill(cosTheta,weight);
                hists["cosThetaHel_mBin"+to_string(m)]->Fill(cosTheta_hel,weight);
                hists["phiGJ_mBin"+to_string(m)]->Fill(phi,weight);
                hists["phiHel_mBin"+to_string(m)]->Fill(phi_hel,weight);
                hists["t_mBin"+to_string(m)]->Fill(mandelstam_t,weight);
                hists["ebeam_mBin"+to_string(m)]->Fill(beam_e,weight);
                ((TH2F*)hists["Mpi0etaCosThetaGJ_mBin"+to_string(m)])->Fill(resonance.M(),cosTheta,weight);
                ((TH2F*)hists["Mpi0etaCosThetaHel_mBin"+to_string(m)])->Fill(resonance.M(),cosTheta_hel,weight);
                ((TH2F*)hists["Mpi0etaPhiGJ_mBin"+to_string(m)])->Fill(resonance.M(),phi,weight);
                ((TH2F*)hists["Mpi0etaPhiHel_mBin"+to_string(m)])->Fill(resonance.M(),phi_hel,weight);
                ((TH2F*)hists["cosThetaPhiGJ_mBin"+to_string(m)])->Fill(cosTheta,phi,weight);
                ((TH2F*)hists["cosThetaPhiHel_mBin"+to_string(m)])->Fill(cosTheta_hel,phi_hel,weight);
            }
        }
    }

    return hists; 
}


map<string, TH1*> loadAndCombine(vector<string> files, string tag){
    map<string,TH1*> hists=getHists(files[0],tag);
    for (int i=1; i<(int)files.size(); ++i){
        map<string, TH1*> tmp=getHists(files[i],tag+to_string(i));
        for (auto const& mapElement: tmp){
            hists[mapElement.first]->Add(tmp[mapElement.first]);
        }
    }
    return hists;
}



void drawAmptoolsVar(){
    gSystem->Exec("mkdir -p drawAmptoolsVar");

    vector<string> flatHistFiles;
    vector<string> flatThrownHistFiles;
    vector<string> totHistFiles;
    vector<string> sbHistFiles;
    vector<string> pols={"000","045","090","135"};
    string st="010020";
    string sm="104156";
    for (auto pol: pols){
        flatHistFiles.push_back("/d/grid17/ln16/dselector_v3/phase1_selected/t"+st+"_m"+sm+"/pol"+pol+"_t"+st+"_m"+sm+"_FTOT_selected_acc_flat.root");
        flatThrownHistFiles.push_back("/d/grid17/ln16/dselector_v3/phase1_selected/t"+st+"_m"+sm+"/pol"+pol+"_t"+st+"_m"+sm+"_FTOT_gen_data_flat.root");
        totHistFiles.push_back("/d/grid17/ln16/dselector_v3/phase1_selected/t"+st+"_m"+sm+"/pol"+pol+"_t"+st+"_m"+sm+"_DTOT_selected_data_flat.root");
        sbHistFiles.push_back("/d/grid17/ln16/dselector_v3/phase1_selected/t"+st+"_m"+sm+"/pol"+pol+"_t"+st+"_m"+sm+"_DTOT_selected_bkgnd_flat.root");

        //string base="/d//grid17/ln16/myDSelector/amptools/zPhase1_t0103061_e79828890/baseFiles_v3/010020_malte/";
        //flatHistFiles.push_back(base+"amptools_flat_phase1_t010020_e8288_sig_a2_pVHpi0p_"+pol+".root");
        //flatThrownHistFiles.push_back(base+"amptools_flat_gen_phase_1_t010020_e8288_tree_flat_a2_pol"+pol+".root");
        //totHistFiles.push_back(base+"amptools_data_phase1_t010020_e8288_tot_a2_pVHpi0p_"+pol+".root");
        //sbHistFiles.push_back(base+"amptools_data_phase1_t010020_e8288_sb_a2_pVHpi0p_"+pol+".root");

        //flatHistFiles.push_back("/d/grid17/ln16/dselector_v3/kmatrix_selected/tall_m080180/pol000_tall_m080180_F2018_8_selected_acc_flat.root");
        //flatThrownHistFiles.push_back("/d/grid17/ln16/dselector_v3/kmatrix_selected/tall_m080180/pol000_tall_m080180_F2018_8_gen_data_flat.root");
        //totHistFiles.push_back("/d/grid17/ln16/dselector_v3/kmatrix_selected/pol000_tall_m080180_kmatrix_selected_data_flat.root");
        //sbHistFiles.push_back("/d/grid17/ln16/dselector_v3/kmatrix_selected/pol000_tall_m080180_kmatrix_selected_bkgnd_flat.root");
    }
    // genHists need atleast 1 null string if you dont want it draw anything, otherwise it should crash.
    //     We still need getHists to initialize a set of empty histograms
    vector<string> genHistFiles={""};
    //vector<string> genHistFiles={"/d/grid17/ln16/dselector_v3/kmatrix_selected/pol000_tall_m080180_kmatrix_gen_halved_data_flat.root"};

    map<string, TH1*> flatHists = loadAndCombine(flatHistFiles,"flat");
    map<string, TH1*> flatThrownHists = loadAndCombine(flatThrownHistFiles,"flatThrown");
    map<string, TH1*> genHists = loadAndCombine(genHistFiles,"generator");
    map<string, TH1*> totHists = loadAndCombine(totHistFiles,"total");
    map<string, TH1*> sbHists = loadAndCombine(sbHistFiles,"sb");

    TLatex *t = new TLatex();
    t->SetTextColor(kRed+2);
    t->SetTextFont(43);
    t->SetTextSize(36);
    string baseText="Integral";
    for (auto const x: flatHists){
        string var = x.first;
        cout << "Drawing " << var << endl;

        TCanvas* c1=new TCanvas("","",2560,1600);
        c1->Divide(4,2);

        // Setting up signal/eff/effCorrected hists
        TH1* sigHist=(TH1*)totHists[var]->Clone("sigHist"); 
        sigHist->Add(sbHists[var],-1);
        TH1* efficiency=(TH1*)flatHists[var]->Clone("efficiency"); 
        efficiency->Divide(flatThrownHists[var]);
        TH1* corrected=(TH1*)sigHist->Clone("corrected");
        corrected->Divide(efficiency);

        c1->cd(1);
        totHists[var]->SetTitle("Tot Yield");
        totHists[var]->SetMinimum(0);
        totHists[var]->Draw("COLZ HIST");

        c1->cd(2);
        //sbHists[var]->SetLineColor(kRed);
        sbHists[var]->SetTitle("Sideband Yield");
        sbHists[var]->SetMinimum(0);
        sbHists[var]->Draw("COLZ HIST");

        c1->cd(3);
        sigHist->SetTitle("Signal Yield");
        sigHist->SetMinimum(0);
        sigHist->Draw("COLZ HIST");
        int sigYield=(int)sigHist->Integral();
        cout << to_string(sigYield) << endl;
        t->DrawLatexNDC(0.3,0.25,"Integral:");
        t->DrawLatexNDC(0.3,0.20,to_string(sigYield).c_str());

        c1->cd(4);
        corrected->SetTitle("Efficiency Corrected Signal Yield");
        corrected->SetMinimum(0);
        corrected->Draw("COLZ HIST");
        int correctedYield=(int)corrected->Integral();
        cout << "Corrected Hist yield: " << correctedYield << endl;
        t->DrawLatexNDC(0.3,0.25,"Integral:");
        t->DrawLatexNDC(0.3,0.20,to_string(correctedYield).c_str());

        c1->cd(5);
        genHists[var]->SetTitle("Generated");
        genHists[var]->SetMinimum(0);
        genHists[var]->Draw("COLZ HIST");
        int genYield=(int)genHists[var]->Integral();
        cout << "Gen Hist yield: " << genYield << endl;
        t->DrawLatexNDC(0.3,0.25,"Integral:");
        t->DrawLatexNDC(0.3,0.20,to_string(genYield).c_str());

        c1->cd(6);
        flatHists[var]->SetTitle("Flat Recon");
        flatHists[var]->SetMinimum(0);
        flatHists[var]->Draw("COLZ HIST");
        int flatYield=(int)flatHists[var]->Integral();
        cout << "Flat Hist yield: " << flatYield << endl;
        t->DrawLatexNDC(0.3,0.25,"Integral:");
        t->DrawLatexNDC(0.3,0.20,to_string(flatYield).c_str());

        c1->cd(7);
        flatThrownHists[var]->SetTitle("Flat Thrown");
        flatThrownHists[var]->SetMinimum(0);
        flatThrownHists[var]->Draw("COLZ HIST");
        int flatThrownYield=(int)flatThrownHists[var]->Integral();
        cout << "Flat Thrown Hist yield: " << flatThrownYield << endl;
        t->DrawLatexNDC(0.3,0.25,"Integral:");
        t->DrawLatexNDC(0.3,0.20,to_string(flatThrownYield).c_str());

        c1->cd(8);
        efficiency->SetTitle("Efficiency");
        efficiency->SetMinimum(0);
        efficiency->Draw("COLZ HIST");

        c1->SaveAs(("drawAmptoolsVar/"+foutTag+"_"+var+"-datBkngd.png").c_str());
        
        if (var=="Mpi0eta_mBin0"){
    	    ofstream logFile;
    	    logFile.open("drawAmptoolsVar/corrected_yields.txt");
            Int_t n = corrected->GetNbinsX();
            for (Int_t i=1; i<=n; i++) {
                if (doAccCorr){
                    logFile << corrected->GetBinLowEdge(i)+corrected->GetBinWidth(i)/2 << " " << 
                          corrected->GetBinContent(i) << endl; }
                else{
                    logFile << sigHist->GetBinLowEdge(i)+sigHist->GetBinWidth(i)/2 << " " << 
                          sigHist->GetBinContent(i) << endl; }
            }
        }
        delete c1;
    }
}

















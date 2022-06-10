int nBins=1;
string fitName="EtaPi0";
vector<string> pols={"000","045","090","135","allPols"};

vector<string> groups={"_S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++","_S0+-","_S0++","_D1--","_D0+-","_D1+-","_D0++","_D1++","_D2++","_D1--_pD1--","_D0+-_pD0+-","_D1+-_pD1+-","_D0++_pD0++","_D1++_pD1++","_D2++_pD2++","_pD1--","_pD0+-","_pD1+-","_pD0++","_pD1++","_pD2++","_S0+-_S0++","_D1--_D0+-_D1+-_D0++_D1++_D2++","_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++","_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++","_S0++_D0++_D1++_D2++_pD0++_pD1++_pD2++","_S0+-_D1--_D0+-_D1+-_pD1--_pD0+-_pD1+-"};
void overlaySingleBin(int iBin,int nBins, vector<string> names1D, vector<TCanvas*> allCanvases, string selectPol){
        gStyle->SetOptStat(kFALSE);
        string folder = ".";

        TH1F *any1DHist_dat;
        TH1F *any1DHist_acc;
        TH1F *any1DHist_bkg;
        TH1F *any1DHist_sig;

        TH1F *any1DHist_sig_all;
        TH1F *any1DHist_acc_all;

        cout << "Defined some variables..." << endl;

        int igroup=1;
        float tot_acc_yield;
        for (auto group:groups){
	    string outputFile = "/etapi_plot"+group+".root";
            cout << "opening: " << folder+outputFile << endl;
            TFile* infile = TFile::Open((folder+outputFile).c_str());

            TCanvas* c2=new TCanvas("","",1400,900); // plot all variables on a canvas for a given waveset

            int ncols=ceil(sqrt(groups.size()));
            int nrows=floor(sqrt(groups.size()));

            c2->Divide(ncols,nrows);
            for (int histIdx=0; histIdx<(int)names1D.size(); ++histIdx){
                vector<string> allPols;
                if (selectPol=="allPols")
                    allPols={"000","045","090","135"};
                else
                    allPols={selectPol};
                for (int ipol=0; ipol<(int)allPols.size(); ++ipol){
                    infile->GetObject((fitName+"_"+allPols[ipol]+"_"+names1D[histIdx]+"dat").c_str(),any1DHist_dat);
                    infile->GetObject((fitName+"_"+allPols[ipol]+"_"+names1D[histIdx]+"acc").c_str(),any1DHist_acc);

                    if (infile->GetListOfKeys()->Contains((fitName+"_"+allPols[ipol]+"_"+names1D[histIdx]+"bkg").c_str())){
                        cout << "Bkg file is included! Will add the distribution onto acc" << endl;
                        infile->GetObject((fitName+"_"+allPols[ipol]+"_"+names1D[histIdx]+"bkg").c_str(),any1DHist_bkg);
    	    	        any1DHist_bkg->SetFillColorAlpha( kBlue-6,0.5);
    	    	        any1DHist_bkg->SetLineColor(0);
                    }
                    any1DHist_sig=(TH1F*)any1DHist_dat->Clone("signal");
                    any1DHist_sig->Add(any1DHist_bkg,-1);
                    if (ipol==0){
                        any1DHist_sig_all=(TH1F*)any1DHist_sig->Clone("sig_total");
                        any1DHist_acc_all=(TH1F*)any1DHist_acc->Clone("acc_total");
                    }
                    else{
                        any1DHist_sig_all->Add(any1DHist_sig);
                        any1DHist_acc_all->Add(any1DHist_acc);
                    }
                }

                allCanvases[histIdx]->cd(igroup);
                any1DHist_sig_all->Draw();
                any1DHist_acc_all->Draw("HIST SAME");
                any1DHist_sig_all->SetTitle(group.c_str());
                any1DHist_sig_all->SetMinimum(0);
                allCanvases[histIdx]->Update();
    	    	any1DHist_acc_all->SetFillColorAlpha( kOrange,0.5);
    	    	any1DHist_acc_all->SetLineColor( 0);
                cout << "creating" << endl;

                auto legend = new TLegend(0.75,0.75,1.0,0.9);
                legend->AddEntry(any1DHist_sig,"signal","l");
                legend->AddEntry(any1DHist_acc,"acc","f");
                //legend->AddEntry(any1DHist_bkg,"bkg","f");
                legend->Draw();

                if (igroup==1){
                    // draw a pavetext showing the mass range for only the first pad
                }
                c2->cd(histIdx+1);
                any1DHist_sig_all->Draw();
                any1DHist_acc_all->Draw("HIST SAME");
                c2->Print(("overlayPlots/diagnostic"+group+"_"+selectPol+".pdf").c_str(),"pdf");
            }
            ++igroup;
        }
        // we could have put this into the above loop but then we would have to open the same root file a lot more times
        for (int histIdx=0; histIdx<(int)names1D.size(); ++histIdx){
            if (names1D[histIdx]=="Phi"){
                names1D[histIdx]="BigPhi";
            }
            if (iBin==(nBins-1)){
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+"_"+selectPol+".pdf)").c_str(),"pdf");
                continue;
            }
            if (iBin==0){
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+"_"+selectPol+".pdf(").c_str(),"pdf");
            }
            if (iBin==(nBins-1)){
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+"_"+selectPol+".pdf)").c_str(),"pdf");
            }
            else{
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+"_"+selectPol+".pdf").c_str(),"pdf");
            }
        }
}

void overlayBins(){
    int ngroups=(int)groups.size();
    int flooredRoot=(int)sqrt(ngroups);
    int nrows, ncols;
    if (flooredRoot*flooredRoot>ngroups){
        nrows=flooredRoot;
        ncols=flooredRoot;
    }
    else if (flooredRoot*(flooredRoot+1)>ngroups){
        nrows=flooredRoot;
        ncols=flooredRoot+1;
    }
    else {
        nrows=flooredRoot+1;
        ncols=flooredRoot+1;
    }
    cout << "Dividing pad to have ncols,nrows: " << ncols << ", " << nrows << endl;
    
    TCanvas* anyCanvas;
    vector<TCanvas*> allCanvases;
    std::vector<std::string> names1D = {"Metapi","Metapi_40MeVBin","cosTheta","Phi","phi","psi","t"};
    for (auto name: names1D){
        anyCanvas = new TCanvas(("c"+name).c_str(),"",1440,900);
        anyCanvas->Divide(ncols,nrows);
        allCanvases.push_back(anyCanvas);
    }
    cout << "Defined all the canvases" << endl;
    
    for (auto pol: pols){
        for (int iBin=0; iBin<nBins;++iBin){
            overlaySingleBin(iBin,nBins,names1D,allCanvases,pol);
        }
    }
}

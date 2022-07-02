// RooFit chisquared fitting https://root.cern/doc/v608/rf602__chi2fit_8C_source.html
// This study is for finite-bin-phi systematics 
//   look at pg 118 of will's analysis note https://halldweb.jlab.org/DocDB/0038/003850/006/Beam_Asymmetry_Analysis_Note_2gamma_responds_to_review1.pdf
int numDOFsig_sc = 3;
float pi=TMath::Pi();
Double_t shiftedCos(Double_t *x, Double_t *par){
	return par[0]*(1.0+0.4*par[1]*TMath::Cos(2*(x[0]-par[2])));
}

std::vector<float> drawAndFitSample(int nsamples, int nbins){
    float avg=(float)nsamples/nbins;

    TCanvas* c1 = new TCanvas("","",1440,900);
    gStyle->SetOptFit(111);
    gStyle->SetStatY(1);
    gStyle->SetStatX(1);
    gStyle->SetStatW(0.18);
    gStyle->SetStatH(0.18);
    //gStyle->SetFitFormat("3.2g");

    // parameters
    RooRealVar p0("p0","p0",avg,0,3*avg);
    RooRealVar p1("p1","p1",1,0,2);
    RooRealVar p2("p2","p2",0,-1*pi,pi);

    // Dependent variables
    RooRealVar phi("phi",0,-1*pi,pi);
    phi.setBins(nbins);
    RooPlot* frame = phi.frame();

    // Fit functions and generating function
    RooGenericPdf asymmetry("asymmetry", "1+0.4*cos(2*phi)", RooArgList(phi));
    //RooGenericPdf fit_fcn("fit_fcn", "p0*(1+0.4*p1*cos(2*(phi-p2)))", RooArgList(p0,p1,p2,phi));
    TF1 *fit_sc = new TF1("fit_sc",shiftedCos,-1*pi,pi,numDOFsig_sc); 

    // Sample, bin, chiSq fit 
    RooDataSet *toyMC = asymmetry.generate(RooArgSet(phi),nsamples);

    //RooDataHist* dHist_toyMC = toyMC->binnedClone();
    //fit_fcn.chi2FitTo(*dHist_toyMC,RooFit::FitOptions("E S Q"));

    // Plotting v0
    //toyMC->plotOn(frame);
    //fit_fcn.plotOn(frame, RooFit::LineColor(kRed));
    //fit_fcn.paramOn(frame, RooFit::Layout(0.725,1.0,0.825));
    //frame->Draw();
    //dHist_toyMC->Draw();
    //fit_sc->Draw("SAME");
    //c1->SaveAs("fitAsymmetryPlots_mctest.png");
    
    // Plotting v1
    TH1* dHist_toyMC = toyMC->createHistogram("dHist_toyMC",phi);
    dHist_toyMC->SetMinimum(0);
    dHist_toyMC->SetTitle(("Toy MC;#phi (radians);Events / 2#pi/"+to_string(nbins)).c_str());
    fit_sc->SetParameters(avg,1,0);
    fit_sc->SetLineColor(kRed);
    fitStatus = dHist_toyMC->Fit(fit_sc,"E S");
    //dHist_toyMC->Draw();
    //fit_sc->Draw("SAME");
    //c1->SaveAs("fitAsymmetryPlots_mctest.png");
    
    std::vector<float> pars = {(float)fit_sc->GetParameter(1), (float)fit_sc->GetParError(1)};

    delete dHist_toyMC;
    delete c1;

    return pars;
}

void fitAsymmetryPlots_binSysTest(){
    int nsamples=1000000;
    int nensembles=500;
    std::vector<float> pars;

    ofstream saveCsv;
    saveCsv.open("results_bin_systematics.csv");
    saveCsv << "nbins ith_ensemble par par_err" << endl;

    vector<int> nbinss={10,20,30,40,50,60,70,100,150};
    for (auto nbins: nbinss){
        for (int iensemble=0; iensemble<nensembles; ++iensemble){
            pars = drawAndFitSample(nsamples, nbins);
            saveCsv << nbins << " " << iensemble << " " << pars[0] << " " << pars[1] << endl; 
            cout << "nbins: " << nbins << " ith_ensemble: " << iensemble << endl;
        }
    }

    // ********
    // For multiple MC tests
    // ********
//    int nensembles=500;
//    int nbins=30;
//    TH1F* dHist_parVal = new TH1F("par_val","par_val",50,0.9,1.1);
//    TH1F* dHist_parErr = new TH1F("par_err","par_err",50,0,0.06);
//
//    for (int iensemble=0; iensemble<nensembles; ++iensemble){
//        pars = drawAndFitSample(nsamples, nbins);
//        dHist_parVal->Fill(pars[0]);
//        dHist_parErr->Fill(pars[1]);
//    }
//    
//    TCanvas* c1 = new TCanvas("","",1440,900);
//    c1->Divide(2,1);
//    c1->cd(1);
//    dHist_parVal->Draw();
//    c1->cd(2);
//    dHist_parErr->Draw();
//    c1->SaveAs("fitAsymmetryPlots_fitmcbin.png");
}

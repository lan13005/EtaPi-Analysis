#include "RooRelBreitWigner.h"
//R__ADD_LIBRARY_PATH("/d/grid17/ln16/dselector_v3/study_pwa/simple_fits/")
R__LOAD_LIBRARY(/d/grid17/ln16/dselector_v3/study_pwa/simple_fits/RooRelBreitWigner_C.so)
//#include "/d/grid17/ln16/dselector_v3/study_pwa/simple_fits/RooRelBreitWigner_C.so"

void fit(){
    vector<string> forms={"bw","voigt"};

    //gSystem->Load("/d/grid17/ln16/dselector_v3/study_pwa/simple_fits/RooRelBreitWigner_C.so");
    using namespace RooFit;
    ofstream ofile;
    ofile.open("yields.txt");
    ofile << "t\tform\tnsig\tnsig_err\tnormalization\tpeak\tpeak_err\twidth\twidth_err" << endl;
    
    TCanvas* c = new TCanvas("c1","",1200,900);
    gStyle->SetOptFit(1);

    vector<string> runs={"2017_1","2018_1","2018_8"};
    vector<string> ts={"010020","0200325","0325050","050075","075100"};
    vector<float> mints={0.1,0.2,0.325,0.5,0.75};
    vector<float> maxts={0.2,0.325,0.5,0.75,1.0};

    /// We can include extra bins to the left and right so that we can determine if extra range is needed for normalzation
    //    This is needed since nsig and nbkg only integrate in the fit region, but we care about projecting to the full range
    //    [1.04, 1.80] seems to be enough
    int extraBinsL=0; // 10
    int extraBinsR=0; // 100
    float minm=1.04-0.0076*extraBinsL;
    float maxm=1.80+0.0076*extraBinsR;
    float rangem=maxm-minm;
    int mbins=100+(extraBinsL+extraBinsR);
    float fit_min=1.08;
    float fit_max=1.52;
    int fitBins=int((fit_max-fit_min)/(maxm-minm)*mbins);
    
    float maximum=0;
    for (int it=0; it<(int)ts.size(); ++it){
        float mint=mints[it];
        float maxt=maxts[it];
        string t=ts[it];
        TH1F* hists[3];
        int i=0; // indexes runs
        for (auto run: runs){
            ////////////////////////
            //  LOAD DATA
            ////////////////////////
            string floc="/d/grid17/ln16/dselector_v3/phase1_selected_v4/D"+run+"_selected_nominal_acc_flat.root";
            cout << "loading: " << floc << endl;
            ROOT::RDataFrame d("kin", floc.c_str());
            auto df=d.Filter(Form("mandelstam_t>%f && mandelstam_t<%f",mint,maxt));
            df=df.Filter(Form("Mpi0eta>%f && Mpi0eta<%f",minm,maxm));
            df=df.Filter(Form("-29.0*atan(-1.05*Mpi0eta+2.78)+328 > vanHove_omega"));
            df=df.Filter("BeamAngle != -1"); // REJECT AMO
            auto df_hist=df.Histo1D({("M4g"+to_string(i)+"_data").c_str(),"M4g",mbins,minm,maxm},"Mpi0eta","Weight");
            if (i==0)
                hists[0]=(TH1F*)df_hist->Clone();
            else
                hists[0]->Add((TH1F*)df_hist->Clone());
            cout << "\tcurrent total entries: " << hists[0]->Integral() << endl;
            ++i;
        }

        cout << endl;
        i=0;
        for (auto run: runs){
            ////////////////////////
            //  LOAD ACC
            ////////////////////////
            string floc="/d/grid17/ln16/dselector_v3/phase1_selected_v4/F"+run+"_selected_nominal_acc_flat.root";
            cout << "loading: " << floc << endl;
            ROOT::RDataFrame dacc("kin", floc.c_str());
            auto df=dacc.Filter(Form("mandelstam_t>%f && mandelstam_t<%f",mint,maxt));
            df=df.Filter(Form("Mpi0eta>%f && Mpi0eta<%f",minm,maxm));
            df=df.Filter(Form("-29.0*atan(-1.05*Mpi0eta+2.78)+328 > vanHove_omega"));
            //df=df.Filter(Form("mandelstam_t_thrown>%f && mandelstam_t_thrown<%f",mint,maxt));
            //df=df.Filter(Form("Mpi0eta_thrown>%f && Mpi0eta_thrown<%f",minm,maxm));
            //df=df.Filter("BeamAngle != -1"); // REJECT AMO; NO NEED TO REJECT SINCE ITS NOT SIMULATED
            auto df_hist=df.Histo1D({("M4g"+to_string(i)+"_acc").c_str(),"M4g",mbins,minm,maxm},"Mpi0eta","Weight");
            if (i==0)
                hists[1]=(TH1F*)df_hist->Clone();
            else
                hists[1]->Add((TH1F*)df_hist->Clone());
            cout << "\tcurrent total entries: " << hists[1]->Integral() << endl;
            ++i;
        }

        cout << endl;
        i=0;
        for (auto run: runs){
            ////////////////////////
            //  LOAD GEN
            ////////////////////////
            string floc="/d/grid17/ln16/dselector_v3/phase1_selected_v4/F"+run+"_gen_data_flat.root";
            cout << "loading: " << floc << endl;
            ROOT::RDataFrame dgen("kin", floc.c_str(),{"Weight","mandelstam_t_thrown","Mpi0eta_thrown"});
            auto df=dgen.Filter(Form("mandelstam_t_thrown>%f && mandelstam_t_thrown<%f",mint,maxt));
            df=df.Filter(Form("Mpi0eta_thrown>%f && Mpi0eta_thrown<%f",minm,maxm));
            //df=df.Filter("BeamAngle != -1"); // REJECT AMO; NO NEED TO REJECT SINCE ITS NOT SIMULATED
            auto df_hist=df.Histo1D({("M4g"+to_string(i)+"_data").c_str(),"M4g",mbins,minm,maxm},"Mpi0eta_thrown","Weight");
            if (i==0)
                hists[2]=(TH1F*)df_hist->Clone();
            else
                hists[2]->Add((TH1F*)df_hist->Clone());
            cout << "\tcurrent total entries: " << hists[2]->Integral() << endl;
            ++i;
        }
                
        TH1F* hist=(TH1F*)hists[0]->Clone();
        hist->Multiply(hists[2]);
        hist->Divide(hists[1]);
        if (it==0)
            maximum=hist->GetMaximum()*1.05;
        hist->Draw();

        double hmin = hist->GetXaxis()->GetXmin();
        double hmax = hist->GetXaxis()->GetXmax();
        
        // Declare observable x
        RooAbsPdf* sigPDF;
        RooPlot* frame;
        RooFitResult* fitter;
        int iform=0;
        for (auto form: forms){
            // Draw all frames on a canvas
            c->Clear();
            gPad->SetLeftMargin(0.15);
                     
            RooRealVar x("x","x",hmin,hmax) ;
            RooDataHist dh("dh","dh",x,Import(*hist)) ;
            
            frame = x.frame(Title(" ")) ;
            dh.plotOn(frame,MarkerColor(1),MarkerSize(0.9),MarkerStyle(21));  //this will show histogram data points on canvas 
            //dh.statOn(frame);  //this will display hist stat on canvas
            
            //// SIGNAL 
            RooRealVar mean("mean","mean",1.315, 1.30, 1.33);
            RooRealVar width("width","width",0.11, 0, 0.2);
            RooRealVar sigma("sigma","sigma",0.005, 0, 0.05);

            RooRealVar radius("radius","radius",3.1,3.1,3.1);
            RooRealVar mass_a("mass_a","mass_a",0.548);
            RooRealVar mass_b("mass_b","mass_b",0.135);
            RooRealVar spin("spin","spin",2);
            spin.setConstant(kTRUE);
            mass_a.setConstant(kTRUE);
            mass_b.setConstant(kTRUE);
            radius.setConstant(kTRUE);
            RooRelBreitWigner relBW("relBW","relBW",x,mean,width,radius,mass_a,mass_b,spin);
            RooRealVar mg("mg","mg",0) ;
            RooRealVar sg("resol.","sg",0.0131) ;
            mg.setConstant(kTRUE);
            sg.setConstant(kTRUE);
            RooGaussian resolution("resolution","resolution",x,mg,sg) ;
            //x.setBins(10000,"cache"); // Set #bins to be used for FFT sampling to 10000
            RooFFTConvPdf voigt("voigt","RelBW (x) Gauss", x, relBW, resolution) ;
            //RooVoigtian relBW("relBW","relBW",x,mean,width,sigma);

            //// BKGND
            RooRealVar c0("bkg_c0","coefficient #0", -3,-10,1) ;
            RooExponential bkgnd("bkgnd","bkgnd",x,c0); // Chebychev assumes a constant term already! including one arg is already linear
            //RooRealVar c1("bkg_c1","coefficient #1", 0.5,-5,5) ;
            //RooChebychev bkgnd("bkgnd","bkgnd",x,RooArgList(c0,c1)); // Chebychev assumes a constant term already! including one arg is already linear
            //RooBernstein bkgnd("bkgnd","bkgnd",x,RooArgList(c0,c1));//,c2));
            //RooGenericPdf bkgnd("bkgnd",
            //        //("@1*(1-(@0-"+to_string(minm)+")/"+to_string(rangem)+")"+
            //        //"+@2*((@0-"+to_string(minm)+")/"+to_string(rangem)+")").c_str(),RooArgList(x,c1,c2));
            //        "c1*(1-x/2)+c2*(x/2)",RooArgList(x,c1,c2));

            //// TOTAL
            if (form=="bw")
                sigPDF=&relBW;
            else
                sigPDF=&voigt;
            RooRealVar nsig("nsig","#signal events",100000,0.,500000) ;
            RooRealVar nbkg("nbkg","#background events",100000,0.,1000000) ;
            RooAddPdf sum("sum","sig+bkgnd",RooArgList(*sigPDF,bkgnd),RooArgList(nsig,nbkg)) ; 
            
            x.setRange("signal",fit_min,fit_max);
            x.setRange("full",minm,maxm);
            fitter = sum.fitTo(dh,
//                    Minos(true),
                    SumW2Error(true),
                    Save(true),
                    Range("signal")
                    );
            sum.plotOn(frame,LineColor(4),LineWidth(8));//this will show fit overlay on canvas 
            sum.plotOn(frame,Components(*sigPDF),LineColor(8),LineWidth(8));
            sum.plotOn(frame,Components(bkgnd),LineColor(2),LineStyle(kDashed),LineWidth(8));
            sum.paramOn(frame,Layout(0.7,0.95,0.95)); //this will display the fit parameters on canvas
            frame->getAttText()->SetTextSize(0.03);
            //fitter->Print("v");
            
            int n_param = fitter->floatParsFinal().getSize();
            int ndf = fitBins-n_param;
            float chi_square=frame->chiSquare();
            float reduced_chi_square = chi_square/ndf;
            
            frame->GetXaxis()->SetTitle("M(#eta#pi^{0}) GeV/c^{2}");  frame->GetXaxis()->SetTitleOffset(1.2);
            float binsize = hist->GetBinWidth(1); char Bsize[50]; 
            sprintf(Bsize,"Events / %2.2f",binsize);
            frame->GetYaxis()->SetTitle(Bsize);  
            frame->GetYaxis()->SetTitleOffset(1.2);
            frame->SetMinimum(0);
            frame->SetMaximum(maximum); //frame->GetYaxis()->GetXmax());
            frame->Draw();

            TLatex *tlatex = new TLatex();
            tlatex->SetTextAlign(22);
            tlatex->SetTextColor(kRed+2);
            tlatex->SetTextFont(43);
            tlatex->SetTextSize(36);
            tlatex->DrawLatexNDC(0.7,0.55,Form("t #in [%0.4f,%0.4f]",mint,maxt));
            tlatex->DrawLatexNDC(0.7,0.50,Form("#chi^{2} = %0.3f",chi_square));
            tlatex->DrawLatexNDC(0.7,0.45,Form("NDF = %i",ndf));
            tlatex->DrawLatexNDC(0.7,0.40,Form("#chi^{2}/NDF = %0.3f",reduced_chi_square));

            RooAbsReal* integral_signal=sigPDF->createIntegral(x,NormSet(x),Range("signal"));
            RooAbsReal* integral_full=sigPDF->createIntegral(x,NormSet(x),Range("full"));
            float normalization = integral_signal->getVal()/integral_full->getVal();

            ofile << mint+(maxt-mint)/2 << "\t" << form << "\t" << nsig.getVal()*normalization  << "\t" << nsig.getError()*normalization << "\t"
                << normalization << "\t" 
                << mean.getVal() << "\t" << mean.getError() << "\t"
                << width.getVal() << "\t" << width.getError() 
                << endl;

            c->Print(("results/"+forms[iform]+"_plus_poly_t"+t+".pdf").c_str());  
            c->Print(("results/"+forms[iform]+"_plus_poly_t"+t+".png").c_str());  
            ++iform;
        }
    }
    ofile.close();
}

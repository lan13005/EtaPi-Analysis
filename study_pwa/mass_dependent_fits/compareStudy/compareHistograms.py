#!/usr/bin/python

from ROOT import TFile, TCanvas, TH1F, TLegend

rf="etapi_plot_S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++.root"
t010015=TFile.Open("./010020_cfg0_m104156/010020_cfg0_m104156_t010150/"+rf)
t010015mis=TFile.Open("./010020_cfg0_m104156/010020_cfg0_m104156_t01015mismatchThrown0/"+rf)
t010020=TFile.Open("./010020_cfg0_m104156/010020_cfg0_m104156_t010200/"+rf)
t010020mis=TFile.Open("./010020_cfg0_m104156/010020_cfg0_m104156_t01020mismatchThrown0/"+rf)

c=TCanvas("","",1920,1080)
c.Divide(2,2)

h010015=TH1F()
h010015mis=TH1F()
h010020=TH1F()
h010020mis=TH1F()

variables=["t","Metapi","cosTheta"]

for i,variable in enumerate(variables):
    c.cd(i+1)

    t010015.GetObject("EtaPi0_000_"+variable+"acc",h010015)
    t010015mis.GetObject("EtaPi0_000_"+variable+"acc",h010015mis)
    t010020.GetObject("EtaPi0_000_"+variable+"acc",h010020)
    t010020mis.GetObject("EtaPi0_000_"+variable+"acc",h010020mis)
    
    h010015.SetLineColor(1)
    h010015mis.SetLineColor(2)
    h010020.SetLineColor(8)
    h010020mis.SetLineColor(4)
    
    h010015.SetLineWidth(10)
    h010015mis.SetLineWidth(10)
    h010020.SetLineWidth(10)
    h010020mis.SetLineWidth(10)
    
    if variable=="t":
        h010020mis.GetXaxis().SetRangeUser(0.1,0.2)
    elif variable=="Metapi":
        h010020mis.GetXaxis().SetRangeUser(1.04,1.56)
    
    h010020mis.Draw()
    h010020.Draw("SAME")
    h010015.Draw("SAME")
    h010015mis.Draw("SAME")

    h010015.SetName("010015")
    h010015mis.SetName("010015mis")
    h010020.SetName("010020")
    h010020mis.SetName("010020mis")
    
    legend=TLegend(0.80,0.7,1.0,0.9);
    legend.AddEntry("010015","010015","l")
    legend.AddEntry("010015mis","010015mis","l")
    legend.AddEntry("010020","010020","l")
    legend.AddEntry("010020mis","010020mis","l")
    legend.Draw("SAME")

c.SaveAs("compare.pdf")

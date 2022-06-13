#!/usr/bin/python

import root_pandas as rp
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os
import ROOT
from decimal import Decimal
import pandas as pd
import re
import uproot3 as uproot
import scipy.stats as sts

import mplhep as hep
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])
SMALL_SIZE = 20
MEDIUM_SIZE = 22
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=17)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

baseDir="zExpectedLeakage/"
os.system("mkdir -p "+baseDir)
f=open(baseDir+"tabulate_calculations.csv","w")
f.write("Channel,Weight,signalRegion,Ecenter,Ewidth,Yield,YieldErr,XSec,XSecErr,Flux,FluxErr,Eff,EffErr,BR,FR\n")
simplify={"AccWeight":"AS", "weightASBS":"ASBS"}


def runAnalysis(channel,weightBranch,signalTag):
    getB1 = True if channel=="b1" else False
    getOmega = True if channel=="omega" else False
    getEta = True if channel=="eta" else False
    getF1 = True if channel=="f1" else False
    getF2 = True if channel=="f2" else False
    getEtaPr = True if channel=="etapr" else False
    getA2Pi = True if channel=="a2pi" else False
    
    forceGetNewFluxVals=False # should we get the flux again

    recons=[] # store all the different runs
    throwns=[] # store all the different runs

    # THIS HAS TO BE THE ORDER SINCE MOST SIMULATIONS WERE MADE WITH FA2018
    #   OTHERWISE THE ORDERING OF THE RECONS LIST NEEDS TO BE FIXED
    colors=['blue', 'red', 'black']
    labels=["Fall 2018", "Spring 2018","Spring 2017"] 
    
    #### Plotting sidebands distributions comparing data to expected b1 background leakage
    def getVal(search):
        val=-1
        with open("/d/grid17/ln16/myDSelector/DSelector_etapi.C") as selector:
            for line in selector:
                if search in line:
                    val=float(line.split("=")[1].split(";")[0])
        return val
    def getSidebands():
        with open("/d/grid17/ln16/myDSelector/DSelector_etapi.C") as selector:
            for line in selector:
                pi0Mean=getVal("pi0Mean")
                etaMean=getVal("etaMean")
                pi0Std=getVal("pi0Std")
                etaStd=getVal("etaStd")
                pi0Sig=getVal("pi0Sig")
                etaSig=getVal("etaSig")
                pi0Skip=getVal("pi0Skip")
                etaSkip=getVal("etaSkip")
                pi0SB=getVal("pi0SB")
                etaSB=getVal("etaSB")
        pi0args=[pi0Mean,pi0Std,pi0Sig,pi0Skip,pi0SB]
        etaargs=[etaMean,etaStd,etaSig,etaSkip,etaSB]
        return pi0args,etaargs
    
    fig,ax=plt.subplots(2,2,figsize=(12,9))
    ax=ax.flatten()
    target=1.22*1E-9 # its in INVERSE barns units - http://hadron.physics.fsu.edu/Theses/AErnst_FSU_thesis.pdf
    
    ## To load the data in the trees
    thrownEnergy="Ebeam_thrown"
    dataEnergy="Ebeam_thrown" # "Ebeam" if we want the reconstructed
    columns=[dataEnergy,"AccWeight","weightBSpi0","weightBSeta","weightASBS","Mpi0eta","Mpi0g3","Mpi0g4","Mpi0","Meta","pVH",
            "cosTheta_eta_gj","cosTheta_eta_hel","phi_eta_gj","phi_eta_hel"]
    


    if getB1:
        BR=1*0.0834*0.988*0.988; # assumed b1->omegapi dominant=1, omega->pi0g is 8.24, 99% decay to 2g for the 2 pi0s 
        csAxisLabel=r"$b_1$"+" Cross Section (b)"
        title=r"$b_1\rightarrow\omega\pi^0\rightarrow5\gamma$"
        outputFolder=baseDir+"zB1_efficiency/"
        binedges=np.linspace(6.6,11.4,49)
        binMin=14
        binMax=25
        sp2017=[1.414,1.430,1.398,1.378,1.372,1.367,1.336,1.341,1.341,1.268,1.304]
        sp2017_err=[0.052,0.047,0.047,0.041,0.042,0.037,0.036,0.037,0.031,0.047,0.068]
        sp2018=[1.409,1.383,1.409,1.404,1.388,1.331,1.341,1.320,1.268,1.294,1.362]
        sp2018_err=[0.084,0.083,0.073,0.068,0.073,0.062,0.063,0.058,0.047,0.073,0.131]
        fa2018=[1.524,1.529,1.466,1.482,1.461,1.440,1.393,1.409,1.383,1.367,1.325]
        fa2018_err=[0.052,0.052,0.047,0.047,0.037,0.032,0.037,0.036,0.031,0.047,0.053]
        crossSections=np.array([fa2018,sp2018,sp2017])*1000 
        crossSectionErrs=np.array([fa2018_err,sp2018_err,sp2017_err])*1000
        baseFolder="./zDSelectedBkgndSamples/b1vps_"
        for run in ["2018_8","2018_1","2017_1"]:
            recons.append(rp.read_root(baseFolder+run+"/bkgndSample_recon_acc_flat.root",columns=columns))
            throwns.append(rp.read_root(baseFolder+run+"/bkgndSample_gen_data_flat.root",columns=[thrownEnergy]))
    
    if getEta:
        BR=1*0.3628*0.988*0.988*0.988; # eta->3pi0 33% of the time and each pi0 decays to 2g 99% of the time
        csAxisLabel=r"$\eta$"+" Cross Section (nb)"
        title=r"$\eta\rightarrow 3\pi\rightarrow6\gamma$"
        outputFolder=baseDir+"zEta_efficiency/"
        binedges=np.linspace(8,9,3)
        binMin=0
        binMax=-1
        crossSections=np.array([[37.9,34.5]]) 
        crossSectionErrs=np.array([[0,0]]) 
        baseFolder="./zDSelectedBkgndSamples/eta_to_3pi/"
        recons=[rp.read_root(baseFolder+"bkgndSample_recon_acc_flat.root",columns=columns)]
        throwns=[rp.read_root(baseFolder+"bkgndSample_gen_data_flat.root",columns=[thrownEnergy])]
    
    if getF1:
        BR=0.522*0.394*0.988*0.988; # f1->etapi0pi0 52% of the time, eta->2g ~40% of the time, 2 pi0->2g ~99%
        csAxisLabel=r"$f_1(1285)$"+" Cross Section (nb)"
        title=r"$f_1(1285)\rightarrow\eta\pi\pi\rightarrow6\gamma$"
        outputFolder=baseDir+"zF1_efficiency/"
        binedges=np.linspace(8,9,3)
        binMin=0
        binMax=-1
        crossSections=[37.9,34.5]
        crossSectionErrs=np.array([[0,0]]) 
        baseFolder="./zDSelectedBkgndSamples/f1_1285_to_etapipi/"
        recons=[rp.read_root(baseFolder+"/bkgndSample_recon_acc_flat.root",columns=columns)]
        throwns=[rp.read_root(baseFolder+"/bkgndSample_gen_data_flat.root",columns=[thrownEnergy])]
    
    if getEtaPr:
        BR=0.224*0.394*0.988*0.988; # etapr->etapi0pi0 22.4% of the time, eta->2g ~40% of the time, 2 pi0->2g ~99%
        csAxisLabel=r"$\eta'$"+" Cross Section (nb)"
        title=r"$\eta'\rightarrow\eta\pi\pi\rightarrow6\gamma$"
        outputFolder=baseDir+"zEtaPr_efficiency/"
        binedges=np.linspace(8,9,3)
        binMin=0
        binMax=-1
        crossSections=np.array([[12.2,11]])/0.394 # George did not divide the cross sections by the branching ratio for eta->2g which is ~40%
        crossSectionErrs=np.array([[0.234,0.234]]) 
        baseFolder="./zDSelectedBkgndSamples/etap_to_etapipi/"
        recons=[rp.read_root(baseFolder+"/bkgndSample_recon_acc_flat.root",columns=columns)]
        throwns=[rp.read_root(baseFolder+"/bkgndSample_gen_data_flat.root",columns=[thrownEnergy])]

    if getA2Pi:
        BR=0.145*0.394*0.988*0.988; # a2->etapi 14.5% of the time, eta->2g ~40% of the time, 2 pi0->2g ~99%
        csAxisLabel=r"$a_2\pi$"+" Cross Section (nb)"
        title=r"$a_2\pi\rightarrow\eta\pi\pi\rightarrow6\gamma$"
        outputFolder=baseDir+"zA2Pi_efficiency/"
        binedges=np.linspace(8,9,3)
        binMin=0
        binMax=-1
        crossSections=np.array([[37.9,34.5]]) 
        crossSectionErrs=np.array([[0,0]]) 
        baseFolder="./zDSelectedBkgndSamples/a2pi/"
        recons=[rp.read_root(baseFolder+"/bkgndSample_recon_acc_flat.root",columns=columns)]
        throwns=[rp.read_root(baseFolder+"/bkgndSample_gen_data_flat.root",columns=[thrownEnergy])]

    if getF2:
        BR=0.842*0.988*0.988; # a2->etapi 14.5% of the time, eta->2g ~40% of the time, 2 pi0->2g ~99%
        csAxisLabel=r"$f_2$"+" Cross Section (nb)"
        title=r"$f_2\rightarrow\pi\pi\rightarrow4\gamma$"
        outputFolder=baseDir+"zF2_efficiency/"
        binedges=np.linspace(8,9,2)
        binMin=0
        binMax=-1
        crossSections=np.array([[6.0366596]]) 
        crossSectionErrs=np.array([[0.93989444]]) 
        baseFolder="./zDSelectedBkgndSamples/pi0pi0/"
        recons=[rp.read_root(baseFolder+"/bkgndSample_recon_acc_flat.root",columns=columns)]
        throwns=[rp.read_root(baseFolder+"/bkgndSample_gen_data_flat.root",columns=[thrownEnergy])]
    

    if getOmega:
        BR=0.0834*0.988; # a2->etapi 14.5% of the time, eta->2g ~40% of the time, 2 pi0->2g ~99%
        csAxisLabel=r"$\omega$"+" Cross Section (nb)"
        title=r"$\omega\rightarrow\pi\gamma\rightarrow3\gamma$"
        outputFolder=baseDir+"zOmega_efficiency/"
        binedges=np.linspace(8,9,6)
        binMin=1
        binMax=-2
        crossSections=np.array([[1.356, 1.331, 1.293, 1.275, 1.237]])*1000
        crossSectionErrs=np.array([[0.08, 0.08, 0.08, 0.08, 0.08]])*1000
        baseFolder="./zDSelectedBkgndSamples/omega_pi0g/"
        recons=[rp.read_root(baseFolder+"/bkgndSample_recon_acc_flat.root",columns=columns)]
        throwns=[rp.read_root(baseFolder+"/bkgndSample_gen_data_flat.root",columns=[thrownEnergy])]


    os.system("mkdir -p "+outputFolder)
    xerrs=(binedges[1]-binedges[0])/2
    if binMax==-1:
        includingBinMax=len(bincenters)
    else:
        includingBinMax=binMax
    emin=binedges[binMin]
    emax=binedges[binMax]
    binedges=binedges[binMin:binMax+1] # we want to include the max bin
    bincenters=(binedges[1:]+binedges[:-1])/2
    binwidth=binedges[1]-binedges[0]
    print("bin centers: {}".format(bincenters))
    crossSections=[cs[binMin:binMax+1] for cs in crossSections]
    crossSectionErrs=[err[binMin:binMax+1] for err in crossSectionErrs]
    assert((len(binedges)-1==len(crossSections[0])) and 
        (len(crossSections[0])==len(crossSectionErrs[0])))

    #############################################
    ###### END  
    #############################################



    ####################################
    # Plot Cross Sections 
    ####################################
    for i,crossSection,crossSectionErr in zip(range(len(crossSections)),crossSections,crossSectionErrs):
        ax[0].errorbar(bincenters,crossSection,yerr=crossSectionErr,xerr=xerrs,fmt="ro",ecolor=colors[i],label=labels[i])
        ax[0].set_ylabel(csAxisLabel,size=16)
        ax[0].set_xlabel("Beam Energy (GeV)",size=16)
        ax[0].set_ylim(bottom=0, top=ax[0].get_ylim()[1]*1.2)
        ax[0].axvline(emin,c="black",linestyle="--")
        ax[0].axvline(emax,c="black",linestyle="--")
        ax[0].legend()
    
    
    ####################################
    # Get/Plot reconstruction efficiency 
    ####################################
    if weightBranch=="AccWeight" and selectSignalRegion:
        for i in len(recons):
            recons[i]=recons[i][(recons[i]["weightBSpi0"]==1)&(recons[i]["weightBSeta"]==1)]


    efficiencies=[]
    efficiencyErrors=[]
    for j,recon,thrown in zip(range(len(recons)),recons,throwns):
        dat_counts, edges, _ = ax[0].hist(recon[dataEnergy],bins=binedges,weights=recon[weightBranch],label="8.2<\n$E_{beam}$\n<8.8 GeV")
        thrown_counts, edges, _ = ax[1].hist(thrown[thrownEnergy],bins=binedges,label="8.2<\n$E_{beam}$\n<8.8 GeV")
    
        # We will skip the bins that have exactly 0 entries, cant really calculate an efficiency and error from them due to div-by-zero
        skipSinceZero=[False if dat_count==0 and thrown_count==0 else True for dat_count,thrown_count in zip(dat_counts,thrown_counts)]
        dat_counts=dat_counts[skipSinceZero]
        thrown_counts=thrown_counts[skipSinceZero]
        print(dat_counts)
        print(thrown_counts)
        _bincenters=bincenters[skipSinceZero]
        efficiency=[dat_count/thrown_count for dat_count,thrown_count in zip(dat_counts,thrown_counts)]
        efficiencies.append(efficiency)
        
        # Calculate efficiencies
        dat_counts_err=np.sqrt(abs(dat_counts)) # since we oversubtract sometimes we will have negative yields
        thrown_counts_err=np.sqrt(thrown_counts)
        efficiencyErrors.append(efficiency*np.sqrt( 
            (dat_counts_err/dat_counts)*(dat_counts_err/dat_counts) + 
            (thrown_counts_err/thrown_counts)*(thrown_counts_err/thrown_counts) ))
    
        print("efficiency:")
        print(efficiencies[j])
        print("efficiencies error:")
        print(efficiencyErrors[j])
        
        ax[1].errorbar(_bincenters,efficiencies[j],yerr=efficiencyErrors[j],xerr=xerrs,fmt="ro",ecolor=colors[j],label=labels[j])
        ax[1].set_ylabel("Efficiency",size=16)
        ax[1].set_xlabel("Beam Energy (GeV)",size=16)
        ax[1].ticklabel_format(style='sci')
        ax[1].axvline(emin,c="black",linestyle="--")
        ax[1].axvline(emax,c="black",linestyle="--")

    # convert to numpy arrays
    efficiencies=np.array(efficiencies)
    efficiencyErrors=np.array(efficiencyErrors)
    
    
    ####################################
    # Obtain Tagged Flux In Specific Bins
    ####################################
    print("\n\nGetting tagged flux\n------------------------")
    
    runs=["2018_8","2018_1","2017_1"]
    runStarts=[50677,40856,30274]
    runEnds=[51768,42577,31057]
    rcdbQueries=[""," --rcdb-query='@is_2018production and @status_approved'"," --rcdb-query='@is_2018production and @status_approved and beam_on_current > 49'"]
    fluxCounts=[]
    for i in range(3):
        fluxCounts.append([])
        cmd_base="/d/grid13/gluex/gluex_top/hd_utilities/hd_utilities-1.17/psflux/plot_flux_ccdb.py --begin-run="+str(runStarts[i])+" --end-run="+str(runEnds[i])
        cmd_bins="--num-bins="+str(len(binedges)-1)
        cmd_lowE="--energy-min="+str(binedges[0])
        cmd_uppE="--energy-max="+str(binedges[-1])
        cmds=[cmd_base,cmd_bins,cmd_lowE,cmd_uppE]
        cmd=" ".join(cmds)
        cmd+=rcdbQueries[i]
        if not os.path.exists(outputFolder+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root") or forceGetNewFluxVals:
            print("Running following command:")
            print(cmd)
            os.system(cmd)
            os.system("mv flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root "+outputFolder)
        else:
            print(outputFolder+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root exists already. Lets use this one")
    
        fluxFile=ROOT.TFile.Open(outputFolder+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root")
        fluxHist=fluxFile.Get("tagged_flux")
        for j in range(fluxHist.GetNbinsX()):
            count=fluxHist.GetBinContent(j+1)
            print("Bin{0} counts: {1}".format(j,count))
            fluxCounts[i].append(count)
        fluxCounts[i]=fluxCounts[i][binMin:binMax+1]
    fluxCounts=np.array(fluxCounts)
    fluxErrors=np.sqrt(fluxCounts) # poisson errors for the flux

    ####################################
    # Obtain Tagged Flux In Fine Bins For Scaling Factor
    #   The flux binning we choose to get depends on the binning of the cross sections
    #   but we want to make predictions for ONLY the coherent peak [8.2,8.8]. We can 
    #   estimate the scale factor by getting a finely binned flux histogram with wider range
    #   and integrating between various regions
    # i.e. cross sections measured from [8,9] GeV but we want only [8.2,8.8] then the flux
    #      must be scaled down by some factor
    ####################################
    print("\n\nDetermining Flux Ratio")
    def fluxIntegrator(fluxFiles,minx,maxx):
        fluxIntegrals=[]
        for fluxFile in fluxFiles:
            h=uproot.open(fluxFile)['tagged_flux']
            edges=h.edges
            centers=edges[:-1]+(edges[1]-edges[0])/2
            values=h.values
            
            extrema=[minx,maxx]
            binRange=np.digitize(extrema,centers)
            fluxIntegral=values[binRange[0]:binRange[1]].sum()
            print("({}) Flux Integral between [{},{}] = {}".format(fluxFile,minx,maxx,fluxIntegral))
            fluxIntegrals.append(fluxIntegral)
        return np.array(fluxIntegrals)

    fluxHists=[]
    for i in range(3):
        cmd_base="/d/grid13/gluex/gluex_top/hd_utilities/hd_utilities-1.17/psflux/plot_flux_ccdb.py --begin-run="+str(runStarts[i])+" --end-run="+str(runEnds[i])
        cmd_bins="--num-bins="+str(500)
        cmd_lowE="--energy-min="+str(7)
        cmd_uppE="--energy-max="+str(9)
        cmds=[cmd_base,cmd_bins,cmd_lowE,cmd_uppE]
        cmd=" ".join(cmds)
        cmd+=rcdbQueries[i]
        if not os.path.exists(baseDir+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+"_finerBins.root"):
            print("Running following command:")
            print(cmd)
            os.system(cmd)
            os.system("mv flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root "+baseDir+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+"_fineBins.root")
        else:
            print(baseDir+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+"_finerBins.root exists already. Lets use this one")
        fluxHists.append(ROOT.TFile.Open(baseDir+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+"_finerBins.root").Get("tagged_flux"))
    
    fluxFiles=[baseDir+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+"_finerBins.root" for i in range(len(runStarts))]
    flux8288=fluxIntegrator(fluxFiles,8.2,8.8)
    fluxErange=fluxIntegrator(fluxFiles,emin,emax)
    fluxRatios=flux8288/fluxErange
    print("  Flux Ratios: {}".format(fluxRatios))
    
    ####################################
    # Calculating Yield
    ####################################
    print("\n\nGetting Expected Yields\n------------------------")
    expectedYields=[]
    expectedYieldErrs=[]
    for i in range(3): # Loop over the phase 1 datasets
        j = i if len(recons)==3 else 0
        cs=crossSections[j]
        cserr=crossSectionErrs[j]
        fc=fluxCounts[i] ## use i NOT j, flux will always have 3 runs attached but cross sections might simply reuse the same efficiency + cross section
        fcerr=fluxErrors[i]
        eff=efficiencies[j]
        efferr=efficiencyErrors[j]
        energies=bincenters[binMin:includingBinMax]
        # Currently we assume no errors on the target, branching ratio, and flux ratios
        expectedYield=cs*eff*fc*target*BR*fluxRatios
        expectedYields.append(expectedYield)
        expectedYieldErrs.append(expectedYield*np.sqrt((cserr/cs)*(cserr/cs)+(fcerr/fc)*(fcerr/fc)+(efferr/eff)*(efferr/eff)))
    expectedYields=np.array(expectedYields)
    expectedYieldErrs=np.array(expectedYieldErrs)
    totalExpectedYield=expectedYields.sum(axis=0) # total expected yield across all energy bins
    totalExpectedYieldErr=np.sqrt((expectedYieldErrs**2).sum(axis=0)) # errors add in quadruture

    ax[2].errorbar(energies,totalExpectedYield,yerr=totalExpectedYieldErr,xerr=xerrs,fmt="ro",ecolor='black')
    ax[2].set_ylabel("Phase 1 Expected Yield",size=16)
    ax[2].set_xlabel("Beam Energy (GeV)",size=16)
    #ax[5].set_ylim(bottom=0)
    
    nsigs=3
    lower3SigEstimate=sum(totalExpectedYield-nsigs*totalExpectedYieldErr)
    upper3SigEstimate=sum(totalExpectedYield+nsigs*totalExpectedYieldErr)
    s_lower3SigEstimate='%.2E' % Decimal(lower3SigEstimate)
    s_upper3SigEstimate='%.2E' % Decimal(upper3SigEstimate)
    # might swap upper and lower estimates later but maximumSeparation is always the upper. The most we will subtract or add is determined from upper
    if upper3SigEstimate<lower3SigEstimate: # if upper limit < lower limit
        s_lower3SigEstimate,s_upper3SigEstimate=s_upper3SigEstimate,s_lower3SigEstimate
        lower3SigEstimate,upper3SigEstimate=upper3SigEstimate,lower3SigEstimate
    ax[2].set_title(r"3$\sigma$ Integral: ["+s_lower3SigEstimate+", "+s_upper3SigEstimate+"]",size=16,fontweight='bold') 

    print("\n\nPlotting final plot!\n------------------------")
    plt.tight_layout()
    plt.savefig(outputFolder+"efficiency_"+weightBranch+signalTag+".png")
    
    ###################################
    # Plotting
    ###################################
    dataBaseLoc="zDSelectedBkgndSamples/gluex_"
    datas=[]
    for i,run in enumerate(["2018_8","2018_1","2017_1"]):
        tmp=rp.read_root(dataBaseFolder+run+"/bkgndSample_recon_acc_flat.root",columns=columns)
        datas.append(tmp)
        print("{} has {} entries".format(run,datas[i].AccWeight.sum()))
    data=pd.concat(datas)
    dataFull=data.copy() # make a copy of the data without any selections so we could get a standardized binning
    if weightBranch=="AccWeight" and selectSignalRegion:
        data=data[(data["weightBSpi0"]==1)&(data["weightBSeta"]==1)]
    print("phase1 has {} entries".format(data.AccWeight.sum()))
    
    for varx,label in zip(["Mpi0eta","Mpi0g3","Mpi0g4","cosTheta_eta_gj","cosTheta_eta_hel","phi_eta_gj","phi_eta_hel"],
                [r"$M(4\gamma)$ GeV",r"$M(\pi\gamma_3)$ GeV",r"$M(\pi\gamma_4)$ GeV",r"$cos\theta_{GJ}$ radians",r"$cos\theta_{hel}$ radians",
                    r"$\phi_{GJ}$i degrees",r"$\phi_{hel}$ degrees"]):
        fig,ax=plt.subplots(1,1,figsize=(8,6))
        _, edges = np.histogram(dataFull[varx],bins=100)
        centers=edges[:-1]+(edges[1]-edges[0])/2
        dcount=ax.hist(data[varx],weights=data[weightBranch],bins=edges,histtype='step',color='black',linewidth=2,label="Phase1 Data")[0]
        if len(recons)==3:
            rcounts=[]
            for recon in recons:
                rcount=np.histogram(recon[varx],weights=recon[weightBranch],bins=edges)[0]
                rcount/=rcount.sum() # make a density plot - which simply represents the shape
                assert(abs(rcount.sum()-1)<0.001) 
                rcounts.append(rcount)
            rcounts=np.array(rcounts)
            # to do proper matrix multiplication we need to transpose things a few times. The idea is to
            #   take each normalized recon distribution and construct a new density plot that is the
            #   weighted average of them (the 3 phase1 datasets)- weighted by the expectedYields
            rcounts=(rcounts.T*expectedYields).T.sum(axis=0)/expectedYields.sum()
            assert(abs(rcounts.sum()-1)<0.001) # assert things are still normalized
        else:
            rcounts=np.histogram(recons[0][varx],weights=recons[0][weightBranch],bins=edges)[0]
            rcounts/=rcounts.sum()

        ax.step(centers,rcounts*upper3SigEstimate,color='red',linewidth=2,label=r"Upper 3$\sigma$ limit"+"\n leakage")
        #ax.fill_between(centers,lower3SigEstimate*rcount,upper3SigEstimate*rcount,color='red',alpha=0.5,linewidth=1,label=r"3$\sigma$ limit"+"\n leakage")
        #ax.step(centers,dcount-rcount,color='green',linewidth=2,label="Data-MC") # plot the difference
        ax.set_xlabel(label)
        ax.set_ylabel("Entries / {0:0.3f} GeV".format(edges[1]-edges[0]))
        ax.set_title(title)
        ax.axhline(0,linewidth=1,c='black')
        ax.legend()
        plt.tight_layout()
        plt.savefig(outputFolder+"expectedYield_"+varx+"_"+weightBranch+signalTag+".png")
        ax.set_yscale('log')
        plt.savefig(outputFolder+"expectedYield_"+varx+"_"+weightBranch+signalTag+"_log.png")
    
    if not selectSignalRegion: # dont make these sideband plots if we are going to select the signal  mass region
        pi0args,etaargs=getSidebands()
        fig,axes=plt.subplots(1,2,figsize=(14,6))
        for varx,args,label,ax in zip(["Mpi0","Meta"],[pi0args,etaargs],[r"$M(\gamma_1\gamma_2)$ GeV",r"$M(\gamma_3\gamma_4)$ GeV"],axes):
            edges = np.histogram(data[varx],weights=data["AccWeight"],bins=100)[1]
            centers=edges[:-1]+(edges[1]-edges[0])/2
            dcount=ax.hist(data[varx],weights=data["AccWeight"],bins=edges,histtype='step',color='black',linewidth=2,label="Phase1 Data")[0]
            assert(abs(rcount.sum()-1)<0.001)
            if len(recons)==3:
                rcounts=[]
                for recon in recons:
                    rcount=np.histogram(recon[varx],weights=recon[weightBranch],bins=edges)[0]
                    rcount/=rcount.sum() # make a density plot - which simply represents the shape
                    assert(abs(rcount.sum()-1)<0.001) 
                    rcounts.append(rcount)
                rcounts=np.array(rcounts)
                # to do proper matrix multiplication we need to transpose things a few times. The idea is to
                #   take each normalized recon distribution and construct a new density plot that is the
                #   weighted average of them - weighted by the expectedYields
                rcounts=(rcounts.T*expectedYields).T.sum(axis=0)/expectedYields.sum()
                assert(abs(rcounts.sum()-1)<0.001) # assert things are still normalized
            else:
                rcounts=np.histogram(recons[0][varx],weights=recons[0][weightBranch],bins=edges)[0]
                rcounts/=rcounts.sum()

            ax.step(centers,rcounts*upper3SigEstimate,color='red',linewidth=2,label=r"Upper 3$\sigma$ limit"+"\n leakage")
            #ax.fill_between(centers,lower3SigEstimate*rcount,upper3SigEstimate*rcount,color='red',alpha=0.5,linewidth=1,label=r"3$\sigma$ limit"+"\n leakage")

            ax.axvspan(args[0]-args[2]*args[1],args[0]+args[2]*args[1],color='green',alpha=0.3)
            ax.axvspan(args[0]-(args[2]+args[3]+args[4])*args[1],args[0]-(args[2]+args[3])*args[1],color='red',alpha=0.3)
            ax.axvspan(args[0]+(args[2]+args[3])*args[1],args[0]+(args[2]+args[3]+args[4])*args[1],color='red',alpha=0.3)
            ax.set_xlabel(label)
            ax.set_ylabel("Entries / {0:0.3f} GeV".format(edges[1]-edges[0]))
            ax.set_title(title)
            ax.legend()
            ax.axhline(0,linewidth=1,c='black')
        plt.tight_layout()
        plt.savefig(outputFolder+"sidebands.png")
        
        # THIS SECTION NEEDS SOME WORK
#        # - Need to update these terms -
#        # fluxRatio->fluxRatios
#    for i in range(len(energies)):
#        f.write("{},{},{},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e},{:0.2e}\n".format(
#            channel,simplify[weightBranch],selectSignalRegion,
#            energies[i],binwidth,expectedYields[i],expectedYieldErrs[i],cs[i],cserr[i],fc[i],fcerr[i],eff[i],efferr[i],target,BR,fluxRatio))

def montage(filename, nrows, ncols, baseDir,outputname):
    # search subdirectories for filenames and montage them
    savedHists=[]
    dirs=os.listdir(baseDir)
    for dir in dirs:
        if os.path.isdir(baseDir+dir):
            subdir=os.listdir(baseDir+dir)
            for f in subdir:
                if filename==f:
                    fullpath=baseDir+dir+"/"+f
                    savedHists.append(fullpath)
    files=" ".join(savedHists)
    cmd="montage "+files+" -mode concatenate -tile "+str(nrows)+"x"+str(ncols)+" "+outputname
    print(cmd)
    os.system(cmd)


###    The first AccWeight does not select the mass signal region whereas the second does. weightASBS also does not select on mass region
for weightBranch, selectSignalRegion,signalTag in zip(["AccWeight", "AccWeight", "weightASBS"],[False, True, False],["","_sigRegion",""]):
    channels=["omega"]#,"omega","eta","etapr","f1","a2pi","f2"]
    for channel in channels:
        runAnalysis(channel,weightBranch,signalTag)
    for varx in ["Mpi0eta","Mpi0g3","Mpi0g4","cosTheta_eta_gj","cosTheta_eta_hel","phi_eta_gj","phi_eta_hel"]:
        montage("expectedYield_"+varx+"_"+weightBranch+signalTag+".png",3,3,baseDir,baseDir+"/montage_"+varx+"_"+weightBranch+signalTag+".png")
        montage("expectedYield_"+varx+"_"+weightBranch+signalTag+"_log.png",3,3,baseDir,baseDir+"/montage_"+varx+"_"+weightBranch+signalTag+"_log.png")

f.close()
os.system("cp "+baseDir+"tabulate_calculations.csv /d/home/ln16/notebooks/thesis_cutSelection_results/a2/")














#!/usr/bin/python3

import glob
import pickle 
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
mpl.rcParams['figure.dpi'] = 200
import numpy as np
import os
import uproot3 as uproot
import mplhep as hep
import matplotlib.patches as patches
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])
SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

pd.set_option('display.max_rows', None)

#######################################
###### Make things look prettier #######
#######################################
prettyWave={
    "S0++":"S_{0}^{+}",
    "D2-+":"D_{-2}^{+}",
    "D1-+":"D_{-1}^{+}",
    "D0++":"D_{0}^{+}",
    "D1++":"D_{1}^{+}",
    "D2++":"D_{2}^{+}",
    "S0+-":"S_{0}^{-}",
    "D2--":"D_{-2}^{-}",
    "D1--":"D_{-1}^{-}",
    "D0+-":"D_{0}^{-}",
    "D1+-":"D_{1}^{-}",
    "D2+-":"D_{2}^{-}",
    "D2-+_pD2-+":"D_{-2}^{+}",
    "D1-+_pD1-+":"D_{-1}^{+}",
    "D0++_pD0++":"D_{0}^{+}",
    "D1++_pD1++":"D_{1}^{+}",
    "D2++_pD2++":"D_{2}^{+}",
    "S0+-_pS0+-":"S_{0}^{-}",
    "D2--_pD2--":"D_{-2}^{-}",
    "D1--_pD1--":"D_{-1}^{-}",
    "D0+-_pD0+-":"D_{0}^{-}",
    "D1+-_pD1+-":"D_{1}^{-}",
    "D2+-_pD2+-":"D_{2}^{-}",
    "P0++":"P_{0}^{+}",
    "P1++":"P_{1}^{+}",
    "P0+-":"P_{0}^{-}",
    "P1+-":"P_{1}^{-}",
}

rearrangeWaveNotation={
    "S0+":"S0++",
    "S0-":"S0+-",
    "D-2+":"D2-+",
    "D-2-":"D2--",
    "D-1+":"D1-+",
    "D-1-":"D1--",
    "D0+":"D0++",
    "D0-":"D0+-",
    "D1+":"D1++",
    "D1-":"D1+-",
    "D2+":"D2++",
    "D2-":"D2+-",
}
rearrangeWaveNotationErr={k+"_err":v+"_err" for k,v in rearrangeWaveNotation.items()}
rearrangeWaveNotation={**rearrangeWaveNotation,**rearrangeWaveNotationErr}


def loadMergedPols(fileName,histBaseName,pols):
    '''
    We use amptools' plotter program to make histograms that were output into separate root files
    Use this function to load the root file and grab the histogram you want ~ integrated over some set of polarizations
    '''
    totalValues=0
    for pol in pols:
        hist=uproot.open(fileName)["EtaPi0_"+pol+"_"+histBaseName]
        edges=hist.edges
        width=edges[1]-edges[0]
        value=hist.values
        totalValues+=value
    return totalValues, edges, width
def relBW(x,particle_mass,width):
    ''' Equation for relativistic breit-wigner'''
    gamma=np.sqrt(particle_mass*particle_mass*(particle_mass*particle_mass+width*width))
    k=2*np.sqrt(2)*particle_mass*width*gamma/np.pi/np.sqrt(particle_mass*particle_mass+gamma)
    return k / ((x*x - particle_mass*particle_mass)**2 + (particle_mass*width)**2)
def loadTheory(tbin):
    ''' Load the theory curves for the different m-projections '''
    theoryCurves=pd.read_csv(baseTheoryFolder+"Waves_bin"+str(tbin)+".txt",sep='  ',
                           names=["D0++", "D1++", "D2++", "D1--", "D0+-", "D1+-"], engine='python')
    theoryCurves=theoryCurves.reset_index()
    theoryCurves.rename(columns={"index":"mass"},inplace=True)
    return theoryCurves
def constructAltTheory(tbin):
    '''
    Everytime we call this we have to load all 5 t bins since we have to keep the same proportions across t-bins
    '''
    theories=[loadTheory(i).iloc[:,1:].sum() for i in range(1,6)]
    theories=[theory/theories[0]["D2++"] for theory in theories]
    theory_masses=loadTheory(1).mass
    theory_waves=loadTheory(1).iloc[:,1:].columns

    mapWaves = {wave:relBW(loadTheory(1).mass,1.3182,0.1134)*theories[tbin-1][wave] for wave in theory_waves}
    mapWaves["mass"]=theory_masses
    mapWaves=pd.DataFrame(mapWaves)
    return mapWaves
def loadValue(floc,search="bestMinimum"):
    if not os.path.exists(floc+"etapi_result.fit"):
        return -1
    with open(floc+"etapi_result.fit") as f:
        lines=[line.rstrip().lstrip() for line in f.readlines()]
        line=[float(line.split("\t")[1]) for line in lines if line.split("\t")[0]==search][0]
    return line
def getBestFiles(floc,orderby='iteration',keepOnlyConvergedFits=True):
    ########################################################
    # OBTAIN A LIST OF ALL THE FILES WE WOULD LIKE TO PLOT
    #   Modify the list (i.e. based on best likelihood) to select what you want to show
    ########################################################
    files=np.array([])
    bestFiles=np.array([])
    ts_output=np.array([])
    best_ts_output=np.array([])
    finalStatuses=[]
    for j,t in enumerate(ts):
        fs=np.array([floc+"/"+t+"/"+f+"/" for f in os.listdir(floc+"/"+t) if t in f])
        if len(fs)==0:
            if any( [True if ".root" in f else False for f in os.listdir(floc+"/"+t)] ):
                fs=np.array([floc+"/"+t+"/"])
            else:
                print("getBestFiles cannot find the correct folder")
                print(f"  Search location: {floc}/{t}")
                print(f"  Found files: {fs}")
                exit(1)
        fs=np.array([f for f in fs if "hidden" not in f]) # selectIterations.py program can hide folders from being drawn using the hidden substr in the name
        nlls=np.array([loadValue(f) for f in fs])
        statuses=np.array([loadValue(f,search='lastMinuitCommandStatus') for f in fs])
        estatuses=np.array([loadValue(f,search='eMatrixStatus') for f in fs])

        order=np.argsort(nlls)
        if orderby=='iteration':
            order=np.argsort([int(f.split('/')[-2].split('_')[1]) for f in fs])
        fs=fs[order]
        nlls=nlls[order]
        statuses=statuses[order]
        estatuses=estatuses[order]
        finalStatuses.append([])
        for file,nll,status,estatus in zip(fs,nlls,statuses,estatuses):
            if (status==0)&(estatus==3):
                finalStatuses[j].append(0)
                x='converged'
            elif (status==0)&(estatus==1):
                finalStatuses[j].append(1)
                x='approximate'
            else:
                finalStatuses[j].append(2)
                x='failed'
            print(f"{file.split('/')[-2]} | NLL: {nll:0.0f} : DeltaBestNLL: {nll-nlls[0]:0.0f} : minuit {status:0.0f} : ematrix {estatus:0.0f} : {x}") 
        fs=list(fs)
        finalStatuses[j]=np.array(finalStatuses[j])
        convergedIdxs=np.arange(len(nlls))
        if keepOnlyConvergedFits:
            convergedIdxs=convergedIdxs[finalStatuses[j]==0]
        bestIdx=convergedIdxs[np.argmin(nlls[convergedIdxs])]
        bestFiles=np.append(bestFiles,fs[bestIdx])
        print(f'{fs[bestIdx]} was the best fit file')
        files=np.append(files,fs)
        ts_output=np.append(ts_output,np.array([t]*len(fs)))
        best_ts_output=np.append(best_ts_output,np.array([t]))

    print()
    for finalStatus in finalStatuses:
        print(f'{sum(finalStatus==0)}:{sum(finalStatus==1)}:{sum(finalStatus==2)} FITS OUT OF {len(finalStatus)} converged:approximate:failed') 
    print()

    finalStatuses=np.concatenate(finalStatuses)
    #print(f"files length: {len(files)}")
    #print(f"ts_output files length: {len(ts_output)}")
    return bestFiles, files, ts_output, best_ts_output, finalStatuses


def loadFitFractions(files):
    ''' 
    amptools' plotter program can also output yields + acceptance corrected yields
    In the current setup the acceptance corrected yields are normalized by the total corrected
    yield to obtain a fit fraction. We can just multiply by the total yield to get the corrected
    yield in a wave
    '''
    etapi_plotter_ofile="etapi_plotter_output.log"

    ### Have to prefill a list of keys since each etapi_plotter log file might not contain the same waveset
    #####  i.e. if we are scanning over different wavesets
    waves=["S0++", "S0+-",
            "D2-+", "D1-+", "D0++", "D1++", "D2++",
            "pD2-+", "pD1-+", "pD0++", "pD1++", "pD2++",
            "D2--", "D1--", "D0+-", "D1+-", "D2+-",
            "pD2--", "pD1--", "pD0+-", "pD1+-", 'pD2+-',
            "allD2--", "allD1--", "allD0+-", "allD1+-", "allD2+-",
            "allD2-+", "allD1-+", "allD0++", "allD1++", "allD2++",
#            "P0++", "P1++", "P0+-", "P1+-",
            "D","D+","D-",
            "S","S+","S-",
            "pD","pD+","pD-",
            "allD","allS",
            "allD+","allD-"
            ]
    waveInts_ts={wave:[] for wave in waves}
    waveIntErrs_ts={wave:[] for wave in waves}
    totals=[]
    waves_present_in_nominal=set()

    # Get the number of fits per t-bin so we can properly locate the nominal result
    nFitsPer_t=[]
    for it,t in enumerate(ts):
        nFitsPer_t.append(sum([1 if t in f else 0 for f in files]))
        print(f'fits in {t} bin: {nFitsPer_t[it-1]}')

    def getIdx(fname):
        ###  One way to fill is with the nominal values which is at iteration 0: bars are invisible in the syst plot and will not contribute to systematics
        which_tbin=[1 if t in fname else 0 for t in ts] 
        assert( sum(which_tbin)==1 ) # make sure fname only hints at one particular t-bin
        idx=which_tbin.index(1) # which t-bin are we interested in? 010020=0, 0200325=1 ... with the nominal scheme. Dependsn on 'ts' list
        idx_of_nominal=0 if idx==0 else sum(nFitsPer_t[:idx])
        return idx, idx_of_nominal

    for currentSize,file in zip(range(1,len(files)+1),files):
        print(file)
        fname=file if '.log' in file else file+"/"+etapi_plotter_ofile # the file might already have be a log file as in the case for bootstrapped results
        idx, idx_of_nominal = getIdx(fname)
        if os.path.exists(fname): # if the file exists load it, if it does that then we will let the fill loop include a value
            with open(fname) as fin:
                for line in fin:
                    if "TOTAL EVENTS" in line:
                        total=float(line.split("=")[1].split("+-")[0].rstrip().lstrip())
                        total_err=float(line.split("=")[1].split("+-")[1].rstrip().lstrip())
                        totals.append(total)
                    if line.startswith("FIT FRACTION") and "::" not in line:
                        wave=line.split(" ")[2]
                        waveInt=float(line.split(" ")[4].rstrip().lstrip())
                        waveInt_err=float(line.split(" ")[6].rstrip().lstrip())
                        if wave in waveInts_ts:
                            waveInts_ts[wave].append(waveInt)
                            waveIntErrs_ts[wave].append(waveInt_err)
                            if currentSize==1: # save what the nominal waveset is
                                waves_present_in_nominal.add(wave)
            ### Fill missing waves (not in the nominal) with some value. Need this loop so the next loop will work. 
            for wave in list(set(waves)-set(waves_present_in_nominal)):
                if len(waveInts_ts[wave])!=len(totals):
                    waveInts_ts[wave].append(0) # zero might give us problems since nBarlow would be 0
                    waveIntErrs_ts[wave].append(0)
            ### The following loop will not work if there are missing waves in the nominal model since the index will not have been filled yet
            for wave in waveInts_ts.keys():
                if len(waveInts_ts[wave])!=currentSize:
                    print(f'   ^-- This file is missing {wave} so use index {idx_of_nominal} in waveInts_ts which corresponds to the nominal result for {ts[idx]}')
                    #waveInts_ts[wave].append(waveInts_ts[wave][idx_of_nominal])
                    #waveIntErrs_ts[wave].append(waveIntErrs_ts[wave][idx_of_nominal])
                    waveInts_ts[wave].append(0) # when drawing the fit fraction systematic overview plot we clearly cannot just set to the nominal
                    waveIntErrs_ts[wave].append(0)
                #print(f'{wave} len waveInts: {len(waveInts_ts[wave])} len total:{len(totals)}')
                assert( len(waveInts_ts[wave])==len(totals) )
        else:
            print(f'  ^-- {"/".join(fname.split("/")[-3:])} does not exist! Loading the nominal values in its place')
            print(f'   ^-- This file suggests the use of index {idx_of_nominal} in waveInts_ts which corresponds to the nominal result for {ts[idx]}')
            for wave in waveInts_ts.keys():
                waveInts_ts[wave].append(waveInts_ts[wave][idx_of_nominal])
                waveIntErrs_ts[wave].append(waveIntErrs_ts[wave][idx_of_nominal])
            totals.append(totals[idx_of_nominal])

    ### ADDITIONAL DIAGNOSTICS, IF WE WANT TO SEE THE PAIRING IS CORRECT. ALL T-BINS MUST HAVE SAME NUMBER OF FITS OTHERWISE RESHAPE WILL NOT WORK
    #if len(set(nFitsPer_t))==1:
    #    for wave in waveInts_ts.keys():
    #        print(f'{wave}: {np.array(waveInts_ts[wave]).reshape(-1,nFitsPer_t[0])}')
    #    print(f'totals: {np.array(totals).reshape(-1,nFitsPer_t[0])}')
    #else:
    #    print("Appears that # fits / tbin is not constant. Will not print reshaped values for diagnostic")

    waveInts_ts={k:np.array(v) for k,v in waveInts_ts.items()}
    waveIntErrs_ts={k:np.array(v) for k,v in waveIntErrs_ts.items()}
    # Make sure there are the same number of elements in each wave. If you get an
    #    error here it might be due to the definition of "waves" in this function. It might
    #    be requesting for more waves than exist
    assert len(set([len(v) for k, v in waveInts_ts.items()]))==1
    assert len(set([len(v) for k, v in waveIntErrs_ts.items()]))==1
    return waveInts_ts, waveIntErrs_ts, totals


def combineBR(br1,br2,br3):
    ''' Combine 3 branching ratios by multiplication and propagate their uncertainties '''
    br=br1[0]*br2[0]*br3[0]
    brErr=br*np.sqrt((br1[1]/br1[0])**2+(br2[1]/br2[0])**2+(br3[1]/br3[0])**2)
    return [br,brErr]
def plotWaves(floc,mi_df,ofileTag,orderby,
              selectRef=[0,1],
              resonances=["","p"],
              wavesets=np.array([["S0+-","D1--","D0+-","D1+-"],["S0++","D0++","D1++","D2++"]]),
              #wavesets=np.array([["S0+-","D0+-","D2+-"],["S0++","D0++","D2++"]]),
              plotMI=False,plotTheory=False,plotData=True,plotBest=True,useBSerr_MI=False,
              histtype="fill"):
    '''
    md_floc: folder location containing the etapi_plotter root files
    mi_df: dataframe containing the fit results for mass-indep fits
    ofiletag: output file tag to include in the output file name -> positiveRef_[tag].ftype
    refs: where 0,1 represents -,+ reflectivies to plot. i.e. [0,1] will plot both [1] will just plot + ref
    resonaces=["","p"] where ""=a2 or "p"=a2prime or set resonances=["both"] to merge (for mass-indep comparison)
    plotMI: should we plot the mass-independent fit results located in mi_df?
    plotTheroy: should we overlay Vincent's predictions?
    '''
    fileName="etapi_plot_"
    wavesets=wavesets[selectRef]
    refs=np.array(["negativeRef","positiveRef"])[selectRef]
    naturalities=np.array(["Unnatural\nProduction","Natural\nProduction"])[selectRef]
    xlabel=r"$M(\eta\pi)$ $GeV^2$"
    ylabel=r"Entries / 40 MeV"
    md_curves=[]

    bestFiles, files, ts_output, best_ts_output, finalStatuses = getBestFiles(floc,orderby)
    if plotBest:
        files=bestFiles
        ts_output=best_ts_output

    ########################################################
    #### With multiple reinitializations obtaining a scale factor is a bit ambiguious.
    ####   currently we will just take the average D2++ amplitude as the scale to match
    ########################################################
    theory=constructAltTheory(1)
    fs=[f for f in files if ts[0] in f]
    assert len(fs)>0 # if there are no sub folders found you might be using the older folder scheme
    maxD2pps=[]
    for i,f in enumerate(fs):
        #### Get the maximum value of the D2++ wave in the smallest t-bin so that we can scale all theory curves to it
        histdat, edges, width = loadMergedPols(f+fileName+"D2++.root","Metapi_40MeVBindat",["000","045","090","135"])
        histacc, edgesacc, widthacc = loadMergedPols(f+fileName+"D2++.root","Metapiacc",["000","045","090","135"])
        binScaleFactor=width/widthacc
        maxD2pps.append(histacc.max()*binScaleFactor)
    xmin=edges[:-1][histdat>0][0] # edges has 1 more dimension than the counts, need to account for that
    xmax=edges[1:][histdat>0][-1]
    maxD2pp=np.mean(maxD2pps)
    if plotBest:
        assert(len(maxD2pps)==1)
    scaleFactor=maxD2pp/theory["D2++"].max()

    ########################################################
    # Begin drawing plots we have selected
    ########################################################
    for iws,waveset,ref,naturality in zip(range(len(wavesets)),wavesets,refs,naturalities):
        fig,axes=plt.subplots(4,5,figsize=(20,14),sharex=True)#,sharey=True)
        for prime in resonances:
            for iw,wave in enumerate(waveset):
                for it,t,tLabel in zip(range(len(ts)),ts,tLabels):
                    ############### Load the data to be plotted - ok to overwrite histdat and histacc now
                    if wave[0]=="D":
                        selectWave=prime+wave if prime in ["","p"] else wave+"_p"+wave
                        if prime=="":
                            histLabel=f"$a_2(1320)$ ${prettyWave[wave]}$"
                        elif prime=="p":
                            histLabel=f"$a_2(1700)$ ${prettyWave[wave]}$"
                        else:
                            histLabel=r"$a_2/a_2(1700)$"+f" ${prettyWave[wave]}$"
                    else: #All non-D waves
                        selectWave=prime+wave if prime in [""] else wave
                        histLabel=f"${prettyWave[wave]}$"

                    rootFileName=fileName+selectWave+".root"
                    fs=[floc+"/"+t+"/"+f+"/" for f in os.listdir(floc+"/"+t) if t in f]
                    maxValInT=0 if it==0 else maxValInT
                    for i,fileLoc in enumerate(fs):
                        if not os.path.exists(fileLoc) or not fileLoc in files:
                            continue
                        histdat,_,_=loadMergedPols(fileLoc+rootFileName,"Metapi_40MeVBindat",["000","045","090","135"])
                        histbkg,_,_=loadMergedPols(fileLoc+rootFileName,"Metapi_40MeVBinbkg",["000","045","090","135"])
                        values=histdat-histbkg
                        maxValInT=max(values) if max(values)>maxValInT and it==0 else maxValInT
                        histacc,_,_=loadMergedPols(fileLoc+rootFileName,"Metapiacc",["000","045","090","135"])
                        valuesacc=histacc*binScaleFactor
                        if plotData:
                            hep.histplot(values,edges,c='black',ax=axes[iw,it])

                        if (prime=="" or prime=="both")*(wave[0]=="D"):
                            color="orange"
                        elif (prime=="p")*(wave[0]=="D"):
                            color="royalblue"
                        else:
                            color='darkgray'

                        _label = histLabel if i==0 else ""
                        #_alpha = 1.0 if fileLoc in bestFiles else 0.4
                        #_linewidth = 4 if fileLoc in bestFiles else 2
                        _alpha = 1.0 if fileLoc in bestFiles else 1.0
                        _linewidth = 2 if fileLoc in bestFiles else 3
                        _linestyle = 'solid' if fileLoc in bestFiles else 'dotted'
                        hep.histplot(valuesacc,edgesacc,color=color,ax=axes[iw,it],alpha=_alpha,label=_label,linewidth=_linewidth,
                                histtype=histtype,linestyle=_linestyle)
                        if fileLoc in bestFiles:
                            md_curves.append([ref,prime,wave,t,valuesacc,edgesacc[:-1]+(edgesacc[1]-edgesacc[0])/2])

                    #### Mass independent results 
                    if plotMI:
                        mi=mi_df[mi_df.t==t]
                        if useBSerr_MI:
                            yerr=(mi[wave+"_err_bsl"],mi[wave+"_err_bsu"])
                            label="Mass Indep.\nBS uncert."
                        else:
                            yerr=mi[wave+"_err"]
                            label="Mass Indep.\nMinuit uncert."
                        if (prime=="" or prime=="both"):
                            axes[iw,it].errorbar(mi.mass,mi[wave],yerr=yerr,c='black',fmt='o',markersize=2,label=label,zorder=99)

                    #### Theory curves
                    if plotTheory:
                        theory=loadTheory(it+1)
                        theory=constructAltTheory(it+1)
                        if wave in theory.columns and prime=="":
                            axes[iw,it].plot(theory.mass,theory[wave].values*scaleFactor,c="mediumseagreen",label="Theory",linewidth=3)

                    ## Draw text with t-bin
                    if iw==0:
                        plt.text(0.15, 1.1, tLabel, fontsize=24, transform=axes[iw,it].transAxes)

                    if it==len(ts)-1:
                        axes[iw,it].legend(prop={"size":18},loc=1)#,bbox_to_anchor=(1,0.95))
                    #axes[iw,it].set_ylim(0.001,maxValInT*1.2)
                    axes[iw,it].set_ylim(0.001)
                    axes[iw,it].set_xlim(xmin,xmax)

        axes[-1,2].set_xlabel(xlabel,size=30)
        axes[1,0].set_ylabel(ylabel,size=30)
        axes[1,0].yaxis.set_label_coords(-0.25, 0.1)

        plt.text(0.35, 0.65, naturality, fontsize=24, transform=axes[0,3].transAxes, weight="bold")
        plt.tight_layout()
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.savefig(ofolder+ref+"_"+ofileTag+"."+ftype)

    md_curves=pd.DataFrame(md_curves,columns=["refl","prime","wave","t","intensity","mass"])
    return files, md_curves, ts_output



def getYieldsInLog(files):
    def getYield(floc):
        floc=f'{floc}/stdout.log'
        if os.path.exists(floc):
            option=2
            with open(floc) as f:
                lines=f.readlines()
                ## Need to incorporate ] in the following search otherwise option 1 and 2 overlap...
                if option==1:
                    lines=[line for line in lines if "] Weight Integral" in line]
                elif option==2:
                    lines=[line for line in lines if "] Accidental Subtracted Weight Integral" in line]
                else:
                    lines=[line for line in lines if "] Number of events kept" in line]
            lines=list(set(lines))
            datalines=[float(line.split(" ")[-1].rstrip()) for line in lines if "DTOT" in line and "data" in line]
            bkgndlines=[float(line.split(" ")[-1].rstrip()) for line in lines if "DTOT" in line and "bkgnd" in line]
            #accmclines=[float(line.split(" ")[-1].rstrip()) for line in lines if "FTOT" in line and "acc" in line]
            #print(f'data integrals: {datalines}')
            #print(f'bkgnd integrals: {bkgndlines}')
            #print(f'accmc integrals: {accmclines}')
            if option==1:
                integral=sum(datalines)-sum(bkgndlines) # use this if we calculate yield based on the weight integral
            else:
                integral=sum(datalines)+sum(bkgndlines) # use this if we calculate yield based on the events, or Accidental weighted
            return integral
        else:
            print(f"{floc} does not exist... yields will be set to 0")
            return 0 

    variationYields=[]
    for f in files:
        variationYields.append(getYield(f))
    return variationYields

def getCrossSection(floc,barlow_csv='',bootstrap_folder='',
        orderby='iteration',keepOnlyConvergedFits=True, variations='', dump='',returnYields=False):
    '''
    variations: should be a list containing the name sof the variations you are testing. This will only affect the dumped results
    dump: dump various results into a pickle file so it can be analyzed further later
    '''
    bestFiles, files, ts_output, best_ts_output, finalStatuses = getBestFiles(floc,orderby,keepOnlyConvergedFits)

    variationYields=getYieldsInLog(files)

    midts_output=[midts_map[t] for t in ts_output]
    tLabels_output=[tLabels_map[t] for t in ts_output]
    tBinWidths_output=[tBinWidths_map[t] for t in ts_output]
    xerrs_output=[xerrs_map[t] for t in ts_output]

    ## LOAD CORECTED YIELDS YIELDS
    waveInts_ts, waveIntErrs_ts, totals = loadFitFractions(files)

    #accCorrected_D=waveInts_ts["D"]*totals
    #accCorrectedErr_D=waveIntErrs_ts["D"]*totals
    ## NEW METHOD
    accCorrected_Dpos=waveInts_ts["D+"]*totals
    accCorrectedErr_Dpos=waveIntErrs_ts["D+"]*totals
    accCorrected_Dneg=waveInts_ts["D-"]*totals
    accCorrectedErr_Dneg=waveIntErrs_ts["D-"]*totals
    ## OLD METHOD: IS INCORRECT BUT WORKS FOR D-WAVES THAT SHARE PHASES. Seems like central values are good but errors are not
    #accCorrected_Dpos=np.array([v for k,v in waveInts_ts.items() if k[0]=="D" and k[-1]=="+" and len(k)==4]).sum(axis=0)*totals
    #accCorrectedErr_Dpos=np.sqrt(np.power(np.array([v for k,v in waveIntErrs_ts.items() if k[0]=="D" and k[-1]=="+" and len(k)==4]),2).sum(axis=0))*totals
    #accCorrected_Dneg=np.array([v for k,v in waveInts_ts.items() if k[0]=="D" and k[-1]=="-" and len(k)==4]).sum(axis=0)*totals
    #accCorrectedErr_Dneg=np.sqrt(np.power(np.array([v for k,v in waveIntErrs_ts.items() if k[0]=="D" and k[-1]=="-" and len(k)==4]),2).sum(axis=0))*totals

    bootstrap_folder=bootstrap_folder.rstrip('/')
    if bootstrap_folder!='':
        ## Load the bootstrapped uncertainties. Overwrite the above loaded uncertainties taken from Minuit
        print('overwritting minuit uncertainty (of nominal fit only!) with bootstrapped ones:')
        print(f"  ^-- FROM: {bootstrap_folder}")
        print(f"  ^-- AND FROM: {bootstrap_folder}_fitFrac")
        ##########################################
        #### 1) Update fit fraction uncertainties
        files_bs=glob.glob(f'{bootstrap_folder}_fitFrac/bootstrap_result_mean_std_t*')
        waveInts_ts_bs, waveIntErrs_ts_bs, totals_bs = loadFitFractions(files_bs)
        for value in waveInts_ts_bs.values():
            assert( len(value)==5 )
        for i in range(ntbins):
            # Overwrite all the fit fractions of all waves
            for wave,value in waveIntErrs_ts.items(): 
                nominalIdx=sum(np.array(midts_output)<midts[i]) # nominal index for a particular t-bin
                waveIntErrs_ts[wave][nominalIdx]=waveIntErrs_ts_bs[wave][i] # overwrite nominal errors
        ##########################################
        #### 2) Update Yield uncertainity (This incorporate fluctuations in the poisson constrained yields that fit fracs dont)
        files_bs=glob.glob(f'{bootstrap_folder}/bootstrap_result_mean_std_t*')
        waveInts_ts_bs, waveIntErrs_ts_bs, totals_bs = loadFitFractions(files_bs)
        for value in waveInts_ts_bs.values():
            assert( len(value)==ntbins )
        for i in range(ntbins):
            nominalIdx=sum(np.array(midts_output)<midts[i]) # nominal index for a particular t-bin
            print(f'tbin{i} {waveIntErrs_ts["D+"][i]} {accCorrectedErr_Dpos[nominalIdx]}')
            accCorrectedErr_Dpos[nominalIdx]=waveIntErrs_ts_bs["D+"][i] # Do not need to multipy by totals since files_bs already incorporates it
            accCorrectedErr_Dneg[nominalIdx]=waveIntErrs_ts_bs["D-"][i] # Do not need to multipy by totals since files_bs already incorporates it
        ##########################################

    #################
    # SETUP - LOAD FLUXES
    #################
    spring_2017_lumi=uproot.open(lumi_floc+"flux_30274_31057.root")["tagged_lumi"].values[0]
    spring_2018_1_lumi=uproot.open(lumi_floc+"flux_40856_42577.root")["tagged_lumi"].values[0]
    fall_2018_8_lumi=uproot.open(lumi_floc+"flux_50677_51768.root")["tagged_lumi"].values[0]
    spring_2017_lumi_amo=uproot.open(lumi_floc+"flux_AMO_30274_31057.root")["tagged_lumi"].values[0]
    spring_2018_1_lumi_amo=uproot.open(lumi_floc+"flux_AMO_40856_42577.root")["tagged_lumi"].values[0]
    fall_2018_8_lumi_amo=uproot.open(lumi_floc+"flux_AMO_50677_51768.root")["tagged_lumi"].values[0]
    phase1_polarized_lumi=spring_2017_lumi+spring_2018_1_lumi+fall_2018_8_lumi-spring_2017_lumi_amo-spring_2018_1_lumi_amo-fall_2018_8_lumi_amo
    print(f"2017_1 total lumi: {spring_2017_lumi:0.5e}")
    print(f"2017_1 unpolarized lumi: {spring_2017_lumi_amo:0.5e}")
    print(f"2017_1 polarized lumi: {spring_2017_lumi-spring_2017_lumi_amo:0.5e}")
    print(f"2018_1 total lumi: {spring_2018_1_lumi:0.5e}")
    print(f"2018_1 unpolarized lumi: {spring_2018_1_lumi_amo:0.5e}")
    print(f"2018_1 polarized lumi: {spring_2018_1_lumi-spring_2018_1_lumi_amo:0.5e}")
    print(f"2018_8 total lumi: {fall_2018_8_lumi:0.5e}")
    print(f"2018_8 unpolarized lumi: {fall_2018_8_lumi_amo:0.5e}")
    print(f"2018_8 polarized lumi: {fall_2018_8_lumi-fall_2018_8_lumi_amo:0.5e}")
    print(f"phase1 total lumi: {spring_2017_lumi+spring_2018_1_lumi+fall_2018_8_lumi:0.5e}")
    print(f"phase1 unpolarized lumi: {spring_2017_lumi_amo+spring_2018_1_lumi_amo+fall_2018_8_lumi_amo:0.5e}")
    print(f"phase1 polarized lumi: {phase1_polarized_lumi:0.5e}")

    ############################################################
    ######### LOAD SYSTEMATIC OUTPUT FROM ANOTHER PROGRAM ######
    ############################################################
    if barlow_csv!='': 
        barlows=pd.read_csv(barlow_csv, dtype={'tbin':str})
        barlows['absDiffPos']=barlows['diffPos'].abs()
        barlows['absDiffNeg']=barlows['diffNeg'].abs()
        barlows['threshold']=[4 if c in ['event','counting'] else 0 for c in barlows['class']]
        print_output=barlows[((barlows.nBarlowPos.abs()>barlows.threshold)|(barlows.nBarlowNeg.abs()>barlows.threshold))]
        print_output=print_output[['variation','class','nBarlowPos','nBarlowNeg','posNom','negNom','posVar','negVar']]
        print(print_output.sort_values(by=['class','variation']))
        
        crossSectionPosSystErr=[]
        crossSectionNegSystErr=[]
        for t in ts:
            #### METHOD 1: Include all significant variations, sqrt( Sum difference^2 / N-1 )
            #diffsq=barlows[(abs(barlows['nBarlowPos'])>threshold)&(barlows['tbin']==t)]['diffPos']**2
            #errPos=np.sqrt(diffsq.sum()/(len(diffsq)-1)) if len(diffsq)>1 else 0 # cannot calculate an error with less than or equal to 1 data point...
            #diffsq=barlows[(abs(barlows['nBarlowNeg'])>threshold)&(barlows['tbin']==t)]['diffNeg']**2
            #errNeg=np.sqrt(diffsq.sum()/(len(diffsq)-1)) if len(diffsq)>1 else 0
            #### METHOD 2: Include all significant variations, group by class, take max deviation, combine quadratically
            # Pos
            significantPos=barlows[(abs(barlows['nBarlowPos'])>barlows.threshold)&(barlows['tbin']==t)]
            significantPos=significantPos.groupby('class')['absDiffPos'].max().reset_index()
            corr_syst_values = [correlated_systematic*barlows.loc[barlows.tbin==t,'posNom'].iloc[0] for 
                    correlated_systematic in correlated_systematics] # [0], all values are same anyway
            errPos=np.sqrt( (np.concatenate([significantPos.absDiffPos.values,corr_syst_values])**2).sum() ) 
            # Neg
            significantNeg=barlows[(abs(barlows['nBarlowNeg'])>barlows.threshold)&(barlows['tbin']==t)]
            significantNeg=significantNeg.groupby('class')['absDiffNeg'].max().reset_index()
            corr_syst_values = [correlated_systematic*barlows.loc[barlows.tbin==t,'negNom'].iloc[0] for 
                    correlated_systematic in correlated_systematics] # [0], all values are same anyway
            errNeg=np.sqrt( (np.concatenate([significantNeg.absDiffNeg.values,corr_syst_values])**2).sum() )
            ##### FILL THE ERRORS #####
            crossSectionPosSystErr.append(errPos)
            crossSectionNegSystErr.append(errNeg)


    ############################################################
    ############## CALCULATE XSEC #####################
    ############################################################
    ## Polarized lumi basically no error bars
    crossSectionPos=accCorrected_Dpos/phase1_polarized_lumi/a2FullReact_br[0]/tBinWidths_output
    crossSectionNeg=accCorrected_Dneg/phase1_polarized_lumi/a2FullReact_br[0]/tBinWidths_output
    crossSection=crossSectionPos+crossSectionNeg
    parityAsym=(crossSectionPos-crossSectionNeg)/(crossSectionPos+crossSectionNeg)
    # pb to microbarns
    crossSection/=1e6
    crossSectionPos/=1e6
    crossSectionNeg/=1e6
    #for a,b,c,d,e,f in zip(files,accCorrected_Dpos,totals,waveInts_ts['D+'],crossSectionPos,tBinWidths_output):
    #    print(f'{a} {b} {c} {d} {e} {f}')
    # calculate errors 
    crossSectionPosErr=crossSectionPos*np.sqrt(
        (accCorrectedErr_Dpos/accCorrected_Dpos)*(accCorrectedErr_Dpos/accCorrected_Dpos)+
        (a2FullReact_br[1]/a2FullReact_br[0])*(a2FullReact_br[1]/a2FullReact_br[0])
                                        ) 
    crossSectionNegErr=crossSectionNeg*np.sqrt(
        (accCorrectedErr_Dneg/accCorrected_Dneg)*(accCorrectedErr_Dneg/accCorrected_Dneg)+
        (a2FullReact_br[1]/a2FullReact_br[0])*(a2FullReact_br[1]/a2FullReact_br[0])
                                        ) 

    ###### HERE IS OUR CHANCE TO INCORPORATE THE SYSTEMATICS ########
    crossSectionPosTotalErr=crossSectionPosErr # always initialize these but we can skip them later so we dont mess with the for loops later
    crossSectionNegTotalErr=crossSectionNegErr
    repeat_syst=int(len(crossSectionPosErr)/ntbins) # a single systematic for each t-bin, but sometimes we would like to plot additional xsec variations
    if barlow_csv!='':
        assert(len(crossSectionPosErr)%ntbins==0)
        crossSectionPosSystErr=np.array(crossSectionPosSystErr*repeat_syst)
        crossSectionNegSystErr=np.array(crossSectionNegSystErr*repeat_syst)
        crossSectionSystErr=np.sqrt(crossSectionPosSystErr**2+crossSectionNegSystErr**2)
    else: # if no systematics file is supplied, will just set to zero
        crossSectionPosSystErr=np.zeros(repeat_syst*ntbins)
        crossSectionNegSystErr=np.zeros(repeat_syst*ntbins)
        crossSectionSystErr=np.zeros(repeat_syst*ntbins)

    crossSectionPosTotalErr=np.sqrt(crossSectionPosTotalErr**2+crossSectionPosSystErr**2)
    crossSectionNegTotalErr=np.sqrt(crossSectionNegTotalErr**2+crossSectionNegSystErr**2)

    ##### Propagate erorrs 
    crossSectionErr=np.sqrt(crossSectionPosErr**2+crossSectionNegErr**2)
    crossSectionTotalErr=np.sqrt(crossSectionPosTotalErr**2+crossSectionNegTotalErr**2)
    parityAsymErr=np.sqrt(
        np.power((2*crossSectionNeg/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionPosErr),2)+
        np.power((2*crossSectionPos/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionNegErr),2)
    )
    parityAsymSystErr=np.sqrt(
        np.power((2*crossSectionNeg/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionPosSystErr),2)+
        np.power((2*crossSectionPos/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionNegSystErr),2)
    )
    parityAsymTotalErr=np.sqrt(
        np.power((2*crossSectionNeg/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionPosTotalErr),2)+
        np.power((2*crossSectionPos/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionNegTotalErr),2)
    )

    xsecpos=(crossSectionPos, crossSectionPosErr, crossSectionPosSystErr, crossSectionPosTotalErr)
    xsecneg=(crossSectionNeg, crossSectionNegErr, crossSectionNegSystErr, crossSectionNegTotalErr)
    xsec=(crossSection, crossSectionErr, crossSectionSystErr, crossSectionTotalErr)
    pa=(parityAsym, parityAsymErr, parityAsymSystErr, parityAsymTotalErr)
    additionalInfo=(tLabels_output,f'{phase1_polarized_lumi:0.3f}',f'${a2FullReact_br[0]:0.4f}\pm{a2FullReact_br[1]:0.4f}$')

    if dump!='':
        dumpfile={
                'xsecpos':xsecpos,
                'xsecneg':xsecneg,
                'pa':pa,
                'midts_output':midts_output,
                'tWidths_output':tBinWidths_output,
                'waveInts_ts':waveInts_ts,
                'waveIntErrs_ts':waveIntErrs_ts,
                'totals':totals,
                'finalStatuses':finalStatuses,
                'variations':variations}
        with open(dump, 'wb') as ofile:
            pickle.dump(dumpfile,ofile)
    if returnYields:
        return xsecpos, xsecneg, xsec, pa, midts_output, tBinWidths_output, waveInts_ts, waveIntErrs_ts, totals, finalStatuses, additionalInfo, variationYields
    else:
        return xsecpos, xsecneg, xsec, pa, midts_output, tBinWidths_output, waveInts_ts, waveIntErrs_ts, totals, finalStatuses, additionalInfo

#################
# SETUP - BASIC VARIABLES
#################
mi_waveset="D-1-_D0+_D0-_D1+_D1-_D2+_S0+_S0-"
a2_br=[0.145,0.012]
a2prime_br=[0.036,0.011]
eta_br=[0.3941,0.002]
pi0_br=[0.98823,0.00034]
target=1 # 1.22*1e-9 # set to 1 = no longer used, we use luminosity directly

ntbins=5 # select how much t-bins you will use
xerrs=np.array([0.05,0.0625,0.0875,0.125,0.125])
tBinWidths=2*xerrs
ts=["010020","0200325","0325050","050075","075100"]
#midts=[0.15,0.2625,0.4125,0.625,0.875] # simple center
midts=[0.1501,0.2560,0.4038,0.6198,0.8691] # 3/28/23 weighted center of the t-bin, dependent on weights and t-distribution
tLabels=[r"$0.1<-t<0.2~\textrm{GeV}^2$", r"$0.2<-t<0.325~\textrm{GeV}^2$", r"$0.325<-t<0.5~\textrm{GeV}^2$", r"$0.5<-t<0.75~\textrm{GeV}^2$", r"$0.75<-t<1.0~\textrm{GeV}^2$"]
tLabels_unitless=[r"$0.1<-t<0.2$", r"$0.2<-t<0.325$", r"$0.325<-t<0.5$", r"$0.5<-t<0.75$", r"$0.75<-t<1.0$"]
#ntbins=4 # select how much t-bins you will use
#xerrs=np.array([0.05625,0.075,0.10625,0.125])
#tBinWidths=2*xerrs
#ts=['01502625','0262504125','041250625','06250875']
#midts=[0.20625,0.3375,0.51875,0.75]
#tLabels=[f"0.15<t<0.2625 $GeV^2$", f"0.2625<t<0.4125 $GeV^2$", f"0.4125<t<0.625 $GeV^2$", f"0.625<t<0.875 $GeV^2$"]
#tLabels_unitless=[r"$0.15<-t<0.2625$", r"$0.2625<-t<0.4125$", r"$0.4125<-t<0.625$", r"$0.625<-t<0.875$"]
xerrs=xerrs[:ntbins]
tBinWidths=tBinWidths[:ntbins]
ts=ts[:ntbins]
midts=midts[:ntbins]
tLabels=tLabels[:ntbins]

#ntbins=2 # select how much t-bins you will use
#xerrs=np.array([0.05,0.125])
#tBinWidths=2*xerrs
#ts=["010020","075100"]
#midts=[0.15,0.875]
#tLabels=[f"0.1<t<0.2 $GeV^2$", f"0.75<t<1.0 $GeV^2$"]

mi_waveset_vec=mi_waveset.split("_")
mi_waveset_vec=[rearrangeWaveNotation[ele] for ele in mi_waveset_vec]
a2FullReact_br=combineBR(a2_br,eta_br,pi0_br)
a2primeFullReact_br=combineBR(a2prime_br,eta_br,pi0_br)
midts_map={k:v for k,v in zip(ts,midts)}
tLabels_map={k:v for k,v in zip(ts,tLabels_unitless)}
tBinWidths_map={k:v for k,v in zip(ts,tBinWidths)}
xerrs_map={k:v for k,v in zip(ts,xerrs)}

## These are generally the external systematics. I know, I know, this is a crap place to put it. Im tired though
correlated_systematics=[0.1137,0.03,0.05] # photon reconstruction efficiency, proton recon eff, flux normalization

baseTheoryFolder="theory_pred/"
lumi_floc="flux/" 



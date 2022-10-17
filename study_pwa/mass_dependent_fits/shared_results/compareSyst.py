#!/usr/bin/python3

import pickle
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200
import numpy as np
import os
import uproot3 as uproot
import mplhep as hep
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

def loadValue(floc,search="bestMinimum"):
    with open(floc+"etapi_result.fit") as f:
        lines=[line.rstrip().lstrip() for line in f.readlines()]
        line=[float(line.split("\t")[1]) for line in lines if line.split("\t")[0]==search][0]
    return line

def getBestFiles(floc,verbose=True):
    ########################################################
    # OBTAIN A LIST OF ALL THE FILES WE WOULD LIKE TO PLOT
    #   Modify the list (i.e. based on best likelihood) to select what you want to show
    ########################################################
    files=np.array([])
    bestFiles=np.array([])
    ts_output=np.array([])
    for t in ts:
        fs=np.array([floc+"/"+t+"/"+f+"/" for f in os.listdir(floc+"/"+t) if t in f])
        if len(fs)==0:
            if any( [True if ".root" in f else False for f in os.listdir(floc+"/"+t)] ):
                fs=np.array([floc+"/"+t+"/"])
            else:
                print("getBestFiles cannot find the correct folder")
                exit(1)
        fs=np.array([f for f in fs if "hidden" not in f]) # selectIterations.py program can hide folders from being drawn using the hidden substr in the name
        nlls=np.array([loadValue(f) for f in fs])
        #order=np.argsort(nlls)
        order=np.argsort([int(f.split('/')[-2].split('_')[1]) for f in fs])
        fs=fs[order]
        nlls=nlls[order]
        if verbose:
            [print(f"{file} | NLL: {nll:0.0f} : DeltaBestNLL: {nll-nlls[0]:0.0f}") for file,nll in zip(fs,nlls)]
        fs=list(fs)
        bestFiles=np.append(bestFiles,fs[0])
        files=np.append(files,fs)
        ts_output=np.append(ts_output,np.array([t]*len(fs)))
    return bestFiles,files,ts_output

def loadFitFractions(files):
    ''' 
    amptools' plotter program can also output yields + acceptance corrected yields
    In the current setup the acceptance corrected yields are normalized by the total corrected
    yield to obtain a fit fraction. We can just multiply by the total yield to get the corrected
    yield in a wave
    '''
    etapi_plotter_ofile="etapi_plotter_output.log"
    waveInts_ts={}
    waveIntErrs_ts={}
    totals=[]
    for file in files:
        fname=file+"/"+etapi_plotter_ofile
        with open(fname) as fin:
            for line in fin:
                if "TOTAL EVENTS" in line:
                    total=float(line.split("=")[1].split("+-")[0].rstrip().lstrip())
                    total_err=float(line.split("=")[1].split("+-")[1].rstrip().lstrip())
                if line.startswith("FIT FRACTION") and "::" not in line:
                    wave=line.split(" ")[2]
#                     if "+" in wave or "-" in wave:
#                         continue
                    waveInt=float(line.split(" ")[4].rstrip().lstrip())
                    waveInt_err=float(line.split(" ")[6].rstrip().lstrip())
                    
                    if wave in waveInts_ts:
                        waveInts_ts[wave].append(waveInt)
                        waveIntErrs_ts[wave].append(waveInt_err)
                    else:
                        waveInts_ts[wave]=[waveInt]
                        waveIntErrs_ts[wave]=[waveInt_err]
        totals.append(total)
    
    waveInts_ts={k:np.array(v) for k,v in waveInts_ts.items()}
    waveIntErrs_ts={k:np.array(v) for k,v in waveIntErrs_ts.items()}
    
    return waveInts_ts, waveIntErrs_ts, totals

def combineBR(br1,br2,br3):
    ''' Combine 3 branching ratios by multiplication and propagate their uncertainties '''
    br=br1[0]*br2[0]*br3[0]
    brErr=br*np.sqrt((br1[1]/br1[0])**2+(br2[1]/br2[0])**2+(br3[1]/br3[0])**2)
    return [br,brErr]


a2_br=[0.145,0.012]
a2prime_br=[0.036,0.011]
eta_br=[0.3941,0.002]
pi0_br=[0.98823,0.00034]
target=1.22*1e-9

ntbins=5 # select how much t-bins you will use
xerrs=[0.05,0.0625,0.0875,0.125,0.125]
tBinWidths=np.array([0.1,0.125,0.175,0.25,0.25])
ts=["010020","0200325","0325050","050075","075100"]
midts=[0.15,0.2625,0.4125,0.625,0.875]
lowts=["0.1","0.2",'0.325','0.5','0.75']
upts=["0.2",'0.325','0.5','0.75','1.0']
tLabels=[f"0.1<t<0.2 $GeV^2$", f"0.2<t<0.325 $GeV^2$", f"0.325<t<0.5 $GeV^2$", f"0.5<t<0.75 $GeV^2$", f"0.75<t<1.0 $GeV^2$"]
xerrs=xerrs[:ntbins]
tBinWidths=tBinWidths[:ntbins]
ts=ts[:ntbins]
midts=midts[:ntbins]
tLabels=tLabels[:ntbins]

a2FullReact_br=combineBR(a2_br,eta_br,pi0_br)
a2primeFullReact_br=combineBR(a2prime_br,eta_br,pi0_br)
midts_map={k:v for k,v in zip(ts,midts)}
tBinWidths_map={k:v for k,v in zip(ts,tBinWidths)}
xerrs_map={k:v for k,v in zip(ts,xerrs)}

flux_floc="flux/" 
pkl_output={}

#################
# SETUP - LOAD DATA
#################
parser = argparse.ArgumentParser(description='Draw diagnostic plots')
parser.add_argument('md_floc', help='folder location containing drawn results root files from PlotGen scripts - for mass dependent fits')
parser.add_argument('otag', help='output folder tag to drop results in. Folder will have format systematic_[tag]')
parser.add_argument('variations', help='variations that are tested')
args = parser.parse_args()
md_floc=args.md_floc
otag=args.otag
variations=args.variations
variations=variations.split(";")
print(variations)

if set(ts).issubset(set(os.listdir(md_floc))):
    print(f" * Folder format is good! Plotting...")
else:
    print(f" * Folder format not good! Expects {len(ts)} sub directories named {ts}")
    print(" * Function will not plot things properly! Exiting...")
#    exit(1)

ofolder=md_floc+f"/systematic_{otag}/"
os.system("mkdir -p "+ofolder)

ftype="pdf" # ["png","pdf"] to save the final images as. pdf gives an unimportant warning.

bestFiles, files, ts_output = getBestFiles(md_floc,True)
print(files)

midts_output=[midts_map[t] for t in ts_output]
tBinWidths_output=[tBinWidths_map[t] for t in ts_output]
xerrs_output=[xerrs_map[t] for t in ts_output]

## LOAD CORECTED YIELDS YIELDS
waveInts_ts, waveIntErrs_ts, totals = loadFitFractions(files)

accCorrected_D=waveInts_ts["D"]*totals
accCorrectedErr_D=waveIntErrs_ts["D"]*totals
accCorrected_Dpos=waveInts_ts["D+"]*totals
accCorrectedErr_Dpos=waveIntErrs_ts["D+"]*totals
accCorrected_Dneg=waveInts_ts["D-"]*totals
accCorrectedErr_Dneg=waveIntErrs_ts["D-"]*totals
## OLD METHOD: IS INCORRECT WBUT WORKS FOR D-WAVES THAT SHARE PHASES
accCorrected_Dpos=np.array([v for k,v in waveInts_ts.items() if k[0]=="D" and k[-1]=="+"]).sum(axis=0)*totals
accCorrectedErr_Dpos=np.sqrt(np.power(np.array([v for k,v in waveIntErrs_ts.items() if k[0]=="D" and k[-1]=="+"]),2).sum(axis=0))*totals
accCorrected_Dneg=np.array([v for k,v in waveInts_ts.items() if k[0]=="D" and k[-1]=="-"]).sum(axis=0)*totals
accCorrectedErr_Dneg=np.sqrt(np.power(np.array([v for k,v in waveIntErrs_ts.items() if k[0]=="D" and k[-1]=="-"]),2).sum(axis=0))*totals

#################
# SETUP - LOAD FLUXES
#################
spring_2017_flux=uproot.open(flux_floc+"flux_30274_31057.root")["tagged_flux"].values[0]
spring_2018_1_flux=uproot.open(flux_floc+"flux_40856_42577.root")["tagged_flux"].values[0]
fall_2018_8_flux=uproot.open(flux_floc+"flux_50677_51768.root")["tagged_flux"].values[0]
spring_2017_flux_amo=uproot.open(flux_floc+"flux_AMO_30274_31057.root")["tagged_flux"].values[0]
spring_2018_1_flux_amo=uproot.open(flux_floc+"flux_AMO_40856_42577.root")["tagged_flux"].values[0]
fall_2018_8_flux_amo=uproot.open(flux_floc+"flux_AMO_50677_51768.root")["tagged_flux"].values[0]
phase1_polarized_flux=spring_2017_flux+spring_2018_1_flux+fall_2018_8_flux-spring_2017_flux_amo-spring_2018_1_flux_amo-fall_2018_8_flux_amo
print(f"phase1 polarized flux: {phase1_polarized_flux:0.3e}")

############################################################
############## CALCULATE XSEC #####################
############################################################
## Polarized flux basically no error bars
crossSection=accCorrected_D/phase1_polarized_flux/target/a2FullReact_br[0]/tBinWidths_output
crossSectionErr=crossSection*np.sqrt(
    (accCorrectedErr_D/accCorrected_D)*(accCorrectedErr_D/accCorrected_D)+
    (a2FullReact_br[1]/a2FullReact_br[0])*(a2FullReact_br[1]/a2FullReact_br[0])
                                    )
crossSectionPos=accCorrected_Dpos/phase1_polarized_flux/target/a2FullReact_br[0]/tBinWidths_output
crossSectionPosErr=crossSectionPos*np.sqrt(
    (accCorrectedErr_Dpos/accCorrected_Dpos)*(accCorrectedErr_Dpos/accCorrected_Dpos)+
    (a2FullReact_br[1]/a2FullReact_br[0])*(a2FullReact_br[1]/a2FullReact_br[0])
                                    )
crossSectionNeg=accCorrected_Dneg/phase1_polarized_flux/target/a2FullReact_br[0]/tBinWidths_output
crossSectionNegErr=crossSectionNeg*np.sqrt(
    (accCorrectedErr_Dneg/accCorrected_Dneg)*(accCorrectedErr_Dneg/accCorrected_Dneg)+
    (a2FullReact_br[1]/a2FullReact_br[0])*(a2FullReact_br[1]/a2FullReact_br[0])
                                    )

# nb to microbarns
crossSection/=1000
crossSectionErr/=1000
crossSectionPos/=1000
crossSectionPosErr/=1000
crossSectionNeg/=1000
crossSectionNegErr/=1000

parityAsym=(crossSectionPos-crossSectionNeg)/(crossSectionPos+crossSectionNeg)
parityAsymErr=np.sqrt(
    np.power((2*crossSectionNeg/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionPosErr),2)+
    np.power((2*crossSectionPos/np.power(crossSectionPos+crossSectionNeg,2)*crossSectionNegErr),2)
)

############################################################
################### GLUEX SYS PLOT  ########################
############################################################
markersize=0 
xerrs_output=0 # set temporarily to zero to not show the x error bars

nvariations=int(len(crossSection)/len(ts))-1 # minus 1 since the first will always be the nominal!
print(f"len xsec array: {len(crossSection)}")
print(f"number of variations: {nvariations}")
assert nvariations==len(variations)
xerrs_output=0
markersize=0 

fig,axes=plt.subplots(1,5,figsize=(16,int(0.5*nvariations)+2),sharey=True)
axes=axes.flatten()

cpos='orangered'
cneg='royalblue'

for iax,ax in enumerate(axes):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
#    ax.set_xlim(-0.6,0.6)
    ax.set_xlim(-0.3,0.3)
    ax.set_ylim(-0.5,nvariations+0.5)
    ax.axvline(0,c='black',linewidth=1)

    ax.errorbar(x=0,y=0.1,yerr=0,xerr=crossSectionPosErr[iax*(nvariations+1)]/crossSectionPos[iax*(nvariations+1)],c=cpos,capsize=4)
    ax.errorbar(x=0,y=-0.1,yerr=0,xerr=crossSectionNegErr[iax*(nvariations+1)]/crossSectionNeg[iax*(nvariations+1)],c=cneg,capsize=4)

    ax.set_title(f"{lowts[iax]}< $-t$ <{upts[iax]} $GeV^2$")
    ax.minorticks_off()
    ax.tick_params(axis=u'both', which=u'both',length=0)

#    ax.set_xticks(np.arange(-0.6,0.6,0.3))
    ax.text(0.3, -0.6, r'$\sigma_{pos}$'+f'={crossSectionPos[iax*(nvariations+1)]:0.3f}nb'+'\n'+
                    r'$\sigma_{neg}$'+f'={crossSectionNeg[iax*(nvariations+1)]:0.3f}nb', horizontalalignment='center',verticalalignment='center',
                    size=12)
    ax.set_yticks(range(nvariations+1))
    ax.set_yticklabels(['Nominal Stat.']+variations)

axes[2].set_xlabel(r"$\frac{\sigma_{variation}-\sigma_{nominal}}{\sigma_{nominal}}$")

for it in range(ntbins):
    for xsec, xsec_err, j in zip([crossSectionPos,crossSectionNeg],[crossSectionPosErr,crossSectionNegErr],range(2)):
        diff=[(xsec[(nvariations+1)*it+(iv+1)]-xsec[(nvariations+1)*it])/xsec[(nvariations+1)*it] for iv in range(nvariations)]
        #diff=[(xsec[(nvariations+1)*it]-xsec[(nvariations+1)*it+(iv+1)])+xsec[(nvariations+1)*it] for iv in range(nvariations)]
        errs=[np.sqrt(abs(xsec_err[(nvariations+1)*it]**2-xsec_err[(nvariations+1)*it+(iv+1)]**2))/xsec[(nvariations+1)*it] for iv in range(nvariations)]
        #errs=[np.sqrt(abs(xsec_err[(nvariations+1)*it]**2-xsec_err[(nvariations+1)*it+(iv+1)]**2)) for iv in range(nvariations)]
        y=np.arange(1,nvariations+1)+0.05 if j==0 else np.arange(1,nvariations+1)-0.05
        c=cpos if j==0 else cneg
        axes[it].errorbar(diff,y,xerr=errs,yerr=0,c=c,linewidth=6,ls='none')

plt.tight_layout()
plt.savefig(ofolder+"systOverview."+ftype)
print("Finished drawing plot")


############################################################
############## GLUEX WAVE FRACTION PLOT ####################
############################################################

fig,axes=plt.subplots(1,5,figsize=(16,int(0.5*nvariations)+2),sharey=True)
axes=axes.flatten()

for iax,ax in enumerate(axes):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
#    ax.set_xlim(-0.3,0.3)
    ax.set_ylim(-0.5,nvariations+0.5)
    ax.axvline(0,c='black',linewidth=1)

    ax.set_title(f"{lowts[iax]}< $-t$ <{upts[iax]} $GeV^2$")
    ax.minorticks_off()
    ax.tick_params(axis=u'both', which=u'both',length=0)

    ax.set_yticks(range(nvariations+1))
    ax.set_yticklabels(['Nominal Stat.']+variations)

axes[2].set_xlabel("D-wave\nFitFractions",size=14)

map_cls={
        'D1--':[cneg,'dotted',r'$D_{-1}^{-}$'],
        'D0++':[cpos,'dotted',r'$D_{0}^{+}$'],
        'D0+-':[cneg,'dashed',r'$D_{0}^{-}$'],
        'D1++':[cpos,'dashed',r'$D_{1}^{+}$'],
        'D1+-':[cneg,'solid',r'$D_{1}^{-}$'],
        'D2++':[cpos,'solid',r'$D_{2}^{+}$'],
        #'S0+-':['black','dashed',r'$S_{0}^{-}$ x0.3'],
        #'S0++':['black','solid',r'$S_{0}^{+}$ x0.3'],
        'D':['black','dashed',r'D'],
        'S':['black','solid',r'S'],
        }

print("Finished D-wave fraction plot")
for j in range(2): # [fit fraction, acceptance corrected yields]
    hs=[] # save the plot artists so we can delete them to prepare for the new iteration
    for it in range(ntbins):
        for wave in map_cls.keys():
            dfrac=waveInts_ts[wave] if j==0 else waveInts_ts[wave]*totals
            c,l,label = map_cls[wave]
            dfrac=dfrac[(nvariations+1)*it:(nvariations+1)*(it+1)]
            if wave[0]=="S":
                dfrac=dfrac*0.3
            y=np.arange(0,nvariations+1)
            h,=axes[it].plot(dfrac,y,c=c,linestyle=l,linewidth=2,label=label,marker='o',markersize=3)
            hs.append(h)
    #axes[0].legend(bbox_to_anchor=(0.8, -0.1, 0.2, 0.08),ncol=4,prop={'size':14})
    handles, labels = axes[0].get_legend_handles_labels()
    #3 = 0.175,7 = 0.085
    yleg=0.2425-0.0225*(nvariations+1)
    fig.legend(handles, labels, bbox_to_anchor=(0.06, yleg, 0.35, 0.125),ncol=4,prop={'size':14})

    plt.tight_layout()
    #fig.subplots_adjust(bottom=0.25) # or whatever
    ftag='fitfrac' if j==0 else 'accCorrYield'
    plt.savefig(ofolder+"systDwaveFractions_"+ftag+"."+ftype)
    for h in hs:
        h.remove()

################### DUMP THE RESULTS 
################### DUMP THE RESULTS 
pkl_output["t"]=midts
pkl_output["correct_D"]=accCorrected_D
pkl_output["correct_D_err"]=accCorrectedErr_D
pkl_output["correct_Dpos"]=accCorrected_Dpos
pkl_output["correct_Dpos_err"]=accCorrectedErr_Dpos
pkl_output["correct_Dneg"]=accCorrected_Dneg
pkl_output["correct_Dneg_err"]=accCorrectedErr_Dneg
pkl_output["xsec"]=crossSection
pkl_output["xsec_err"]=crossSectionErr
pkl_output["xsec_pos"]=crossSectionPos
pkl_output["xsec_pos_err"]=crossSectionPosErr
pkl_output["xsec_neg"]=crossSectionNeg
pkl_output["xsec_neg_err"]=crossSectionNegErr
#pkl_output["parity_asym"]=parityAsym
#pkl_output["parity_asym_err"]=parityAsymErr

pickle.dump(pkl_output, open(ofolder+"source.pkl","wb"))


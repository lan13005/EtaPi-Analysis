#!/usr/bin/python3

from helper import *
import time

#################
# SETUP - LOAD DATA
#################
parser = argparse.ArgumentParser(description='Draw diagnostic plots')
parser.add_argument('md_floc', help='folder location containing drawn results root files from PlotGen scripts - for mass dependent fits')
parser.add_argument('-s', nargs='?', default='', help='csv file containing sigma barlow and variation differences')
parser.add_argument('-b', nargs='?', default='', help='folder containing bootstrapped uncertainties')
args = parser.parse_args()
md_floc=args.md_floc
barlow_csv=args.s
bootstrap_folder=args.b

if barlow_csv=='':
    print("\n***********************\nYou did not pass a barlow_csv argument. Not plotting total uncertainties with syst. err!\n***********************")
else:
    print("\n***********************\nYou requested syst. uncertainties to be included by providing barlow_csv\n***********************")
if bootstrap_folder=='':
    print("\n***********************\nYou did not pass a bootstrap_folder argument. Using default MINUIT uncertainties!\n***********************")
else:
    print("\n***********************\nYou requested bootstrapped uncertainties to be included by providing bootstrap_folder\n***********************")
print('\n')
time.sleep(3)

if set(ts).issubset(set(os.listdir(md_floc))):
    print(f" * Folder format is good! Plotting...")
else:
    print(f" * Folder format not good! Expects {len(ts)} sub directories named {ts}")
    print(" * Function will not plot things properly! Exiting...")
#    exit(1)

ofolder=md_floc+"/xsec/"
os.system("mkdir -p "+ofolder)

def plotWaves(floc,mi_df,ofileTag,orderby,
              selectRef=[0,1],
              resonances=["","p"],
              #wavesets=np.array([["S0+-","D1--","D0+-","D1+-","P0+-","P1+-"],["S0++","D0++","D1++","D2++","P0++","P1++"]]),
              wavesets=np.array([["S0+-","D1--","D0+-","D1+-"],["S0++","D0++","D1++","D2++"]]),
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
        fig,axes=plt.subplots(len(waveset),5,figsize=(20,14),sharex=True)#,sharey=True)
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
                    axes[iw,it].set_ylim(0.001,5000)
                    #axes[iw,it].set_ylim(0.001)
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


### Output the results to a csv in the output folder
odf_xsec={}
odf_pa={}

ftype="pdf" # ["png","pdf"] to save the final images as. pdf gives an unimportant warning.

# DRAW THE STAMP PLOTS WITH 5 TBINS AND TMD+S WAVESET
mis=[] # if you want to overlay mass-independent fits you can link them here
# which wave to plot? where ""=a2 or "p"=a2prime or set resonances=["both"] to merge (for mass-indep comparison)
#   plotBest will determine what files end up in files variable
#files, md_curves, ts_output = plotWaves(md_floc,mis,"nominal_a2a2p",orderby='',resonances=["","p"],histtype='step',plotData=True,plotBest=False, plotTheory=True)
#files, md_curves, ts_output = plotWaves(md_floc,mis,"nominal_Dwave",orderby='',resonances=["both"],histtype='step',plotData=True,plotBest=False, plotTheory=True)

orderby='iteration' 
xsecpos, xsecneg, xsec, pa, midts_output, tBinWidths_output, waveInts_ts, waveIntErrs_ts, totals, finalStatuses, additionalInfo = \
        getCrossSection(md_floc,barlow_csv,bootstrap_folder,orderby=orderby)
    
crossSectionPos, crossSectionPosErr, crossSectionPosSystErr, crossSectionPosTotalErr = xsecpos
crossSectionNeg, crossSectionNegErr, crossSectionNegSystErr, crossSectionNegTotalErr = xsecneg
crossSection, crossSectionErr, crossSectionSystErr, crossSectionTotalErr = xsec
parityAsym, parityAsymErr, parityAsymSystErr, parityAsymTotalErr = pa
tLabels_output, phase1_polarized_lumi, branching_frac_of_decay_chain = additionalInfo

#if barlow_csv!='':
#    nvariations=int(len(midts_output)/len(ts))
#    crossSectionPos=crossSectionPos[::nvariations]
#    crossSectionPosErr=crossSectionPosErr[::nvariations]
#    crossSectionPosTotalErr=crossSectionPosTotalErr[::nvariations]
#    crossSectionNeg=crossSectionNeg[::nvariations]
#    crossSectionNegErr=crossSectionNegErr[::nvariations]
#    crossSectionNegTotalErr=crossSectionNegTotalErr[::nvariations]
#    crossSection=crossSection[::nvariations]
#    crossSectionErr=crossSectionErr[::nvariations]
#    crossSectionTotalErr=crossSectionTotalErr[::nvariations]
#    parityAsym=parityAsym[::nvariations]
#    parityAsymErr=parityAsymErr[::nvariations]
#    parityAsymTotalErr=parityAsymTotalErr[::nvariations]
#    midts_output=midts_output[::nvariations]
#    totals=totals[::nvariations]
#    finalStatuses=finalStatuses[::nvariations]
#    for k,v in waveInts_ts.items():
#        waveInts_ts[k]=v[::nvariations]
#    for k,v in waveIntErrs_ts.items():
#        waveIntErrs_ts[k]=v[::nvariations]
#print(len(crossSectionPos))
#print(len(crossSectionPosErr))
#print(len(crossSectionPosTotalErr))
#print(len(crossSectionNeg))
#print(len(crossSectionNegErr))
#print(len(crossSectionNegTotalErr))
#print(len(crossSection))
#print(len(crossSectionErr))
#print(len(crossSectionTotalErr))
#print(len(parityAsym))
#print(len(parityAsymErr))
#print(len(parityAsymTotalErr))


############################################################
############## VINCENT' PREDICTIONS #####################
############################################################
bapa_theory=pd.read_csv(baseTheoryFolder+"Ten_A2_Fit_85GeV.txt",delimiter=";")
bapa_theory=bapa_theory[["beam asymmetry"]]
bapa_theory.columns=["ba"]

crossSection_theory=pd.read_csv(baseTheoryFolder+"Bands_TMD_A2.txt", delim_whitespace=True, header=None)
crossSection_theory=crossSection_theory.loc[:,[0,7,8,9]]
crossSection_theory.columns=["t","cs_low_err","cs","cs_up_err"]

psig_theory=pd.read_csv(baseTheoryFolder+"Psig_Bands_TMD_A2.txt", delim_whitespace=True, header=None)
psig_theory=psig_theory.loc[:,[7,8,9]]
psig_theory.columns=["pa_low_err","pa","pa_up_err"]

crossSection_theory=pd.concat([crossSection_theory,psig_theory,bapa_theory],axis=1)
crossSection_theory["cs_low_err"]=crossSection_theory["cs"]-crossSection_theory["cs_low_err"]
crossSection_theory["cs_up_err"]=crossSection_theory["cs_up_err"]-crossSection_theory["cs"]
crossSection_theory["pa_low_err"]=crossSection_theory["pa"]-crossSection_theory["pa_low_err"]
crossSection_theory["pa_up_err"]=crossSection_theory["pa_up_err"]-crossSection_theory["pa"]
crossSection_theory=crossSection_theory[(crossSection_theory.t<1)&(crossSection_theory.t>0.1)]

crossSection_theory=crossSection_theory[(crossSection_theory.t<1)&(crossSection_theory.t>0.1)]

############################################################
############## GLUEX TOTAL - POS - NEG #####################
############################################################
fig_xsec,axes_xsec=plt.subplots(1,1,figsize=(12,6))
arbitraryFactor=1
tot_alpha=1.0
tot_ewidth=2
tot_capsize=3
tot_capwidth=2
### Used to draw the final results
stat_ewidth=15 if barlow_csv!='' else 15#3
stat_alpha=0.5 if barlow_csv!='' else 0.5#1
### Used to scatter multiple errorbars to compare the resutls
#stat_ewidth=2 if barlow_csv!='' else 2
#stat_alpha=1 if barlow_csv!='' else 1
stat_capsize=0
stat_capwidth=0

def drawRectangle(xy,ax):
    rect = patches.Rectangle(xy, 0.03, 0.1, linewidth=0, edgecolor='k', facecolor='k', alpha=stat_alpha, transform=ax.transAxes)
    ax.add_patch(rect)
    ax.text(xy[0]+0.05,xy[1]+0.04, 'Statisitical Error', size=15, transform=ax.transAxes)

####### SOME CODE TO SHIFT THE POINTS IF WE HAVE MULTIPLE VARIATIONS PER T-BIN #######
if len(crossSection)>ntbins: #crossSection is flattened iterations x tbin so i0t0 i1t0 i0t1 i1t1...
    # shift by fixed amounts if the same # of iterations per t-bin. Useful when combined with selectIterations.py 
    #    to compare two sets of iterations: i.e. bestChiSq vs bestResidual
    unique,counts=np.unique(midts_output, return_counts=True)
    #map_tCounts={k:v for k,v in zip(unique,counts)} 
    i=0
    #jitters=np.linspace(-0.002*(iterations-1),0.002*(iterations-1),iterations) # 17 iterations, [-0.015,0.015] too narrow # deprecated line
    for n,iterations in enumerate(counts): 
        jitters=np.linspace(-0.0008*iterations/2,0.008*iterations/2,iterations) # 17 iterations, [-0.015,0.015] too narrow
        for j,jitter in enumerate(jitters): # add jitter to each iteartion separately
            midts_output[j+i]+=jitter
        i+=iterations
    markersize=0 
elif len(crossSection)==ntbins:
    markersize=6
else:
    print("Programs expected to find atleast 5 tbins, exiting...")
    exit()
######################################################################################

# determine parity asym
for j,poserr,negerr,toterr,ebarwidth,alpha,capsize,capwidth in zip(range(2),[crossSectionPosErr,crossSectionPosTotalErr],[crossSectionNegErr,crossSectionNegTotalErr],
        [crossSectionErr,crossSectionTotalErr],
        [stat_ewidth,tot_ewidth],[stat_alpha,tot_alpha],[stat_capsize,tot_capsize],[stat_capwidth,tot_capwidth]):
    if barlow_csv=='' and j==1:
        continue
    condition=(j==1 and barlow_csv!='') or (j==0 and barlow_csv=='')
    label="Pos. Ref" if condition else ''
    axes_xsec.errorbar(midts_output,crossSectionPos,fmt='o',xerr=0,yerr=poserr,c='tab:red',label=label,markersize=markersize,
            elinewidth=ebarwidth,alpha=alpha,capsize=capsize,markeredgewidth=capwidth)
    label="Neg. Ref" if condition else ''
    axes_xsec.errorbar(midts_output,crossSectionNeg,fmt='o',xerr=0,yerr=negerr,c='tab:blue',label=label,markersize=markersize,
            elinewidth=ebarwidth,alpha=alpha,capsize=capsize,markeredgewidth=capwidth)
    label="GlueX Phase 1\n$E_{\gamma}\in$ [8.2,8.8] GeV" if condition else ''
    axes_xsec.errorbar(midts_output,crossSection,fmt='o',xerr=0,yerr=toterr,c='black',markersize=markersize,
                     label=label,elinewidth=ebarwidth,alpha=alpha,capsize=capsize,markeredgewidth=capwidth)
axes_xsec.set_xlabel(r"-t $GeV^2$",size=30)
axes_xsec.set_ylabel(r"$\frac{d\sigma}{dt}$ [$\frac{\mu b}{GeV^2}$]",size=30)
axes_xsec.plot(crossSection_theory["t"],crossSection_theory["cs"]/arbitraryFactor,c="mediumseagreen",linewidth=3,
             label="TMD predictions\n$E_{\gamma}$=8.5 GeV",zorder=0)
axes_xsec.fill_between(crossSection_theory["t"],
                      (crossSection_theory["cs"]-crossSection_theory["cs_low_err"])/arbitraryFactor,
                      (crossSection_theory["cs"]+crossSection_theory["cs_up_err"])/arbitraryFactor,
                      color="mediumseagreen",alpha=0.6,zorder=0)
if barlow_csv!='':
   drawRectangle((0.73,0.42), axes_xsec)

axes_xsec.tick_params(axis='both',labelsize=24)
axes_xsec.set_xlim(0.1)
axes_xsec.set_ylim(0,0.4)
axes_xsec.ticklabel_format(style='plain')
l1=axes_xsec.legend(prop={"size":16},loc='upper right')
plt.savefig(ofolder+"a2_xsec."+ftype)

############################################################
############## GLUEX TO CLAS #####################
############################################################
clas_e3545={}
clas_e3545["x"]=[0.3, 0.55, 0.9, 1.3]
clas_e3545["y"]=[0.46039, 0.268257, 0.313349, 0.252332]
clas_e3545["ex"]=[0.1, 0.15, 0.2, 0.2]
clas_e3545["ey"]=[0.0645401, 0.0403087, 0.0313953, 0.0314814]

clas_e4555={}
clas_e4555["x"]=[0.3, 0.55, 0.9, 1.3, 1.75]
clas_e4555["y"]=[0.23734, 0.151709, 0.244609, 0.123577,0.0707616]
clas_e4555["ex"]=[0.1, 0.15, 0.2, 0.2, 0.25]
clas_e4555["ey"]=[0.0647159, 0.0319749, 0.0243778, 0.0177268, 0.0117884]

fig,axes=plt.subplots(1,1,figsize=(12,6))

axes.errorbar(clas_e3545['x'],clas_e3545['y'],xerr=clas_e3545['ex'],yerr=clas_e3545['ey'],fmt='o',c='orange',label="CLAS $E_{\gamma}$=[3.5,4.5] GeV")
axes.errorbar(clas_e4555['x'],clas_e4555['y'],xerr=clas_e4555['ex'],yerr=clas_e4555['ey'],fmt='o',c='magenta',label="CLAS $E_{\gamma}$=[4.5,5.5] GeV")

for j,toterr,ebarwidth,alpha,capsize,capwidth in zip(range(2),[crossSectionErr,crossSectionTotalErr],
        [stat_ewidth,tot_ewidth],[stat_alpha,tot_alpha],[stat_capsize,tot_capsize],[stat_capwidth,tot_capwidth]):
    if barlow_csv=='' and j==1:
        continue
    condition=(j==1 and barlow_csv!='') or (j==0 and barlow_csv=='')
    label="GlueX Phase 1\n$E_{\gamma}\in$ [8.2,8.8] GeV" if condition else ''
    axes.errorbar(midts_output,crossSection,fmt='o',xerr=0,yerr=toterr,c='black',markersize=markersize,
                     label=label,elinewidth=ebarwidth,alpha=alpha,capsize=capsize,markeredgewidth=capwidth)
axes.set_xlabel(r"-t $GeV^2$",size=30)
axes.set_ylabel(r"$\frac{d\sigma}{dt}$ [$\frac{\mu b}{GeV^2}$]",size=30)
axes.plot(crossSection_theory["t"],crossSection_theory["cs"]/arbitraryFactor,c="mediumseagreen",linewidth=3,
             label="TMD predictions\n$E_{\gamma}$=8.5 GeV",zorder=0)
axes.fill_between(crossSection_theory["t"],
                      (crossSection_theory["cs"]-crossSection_theory["cs_low_err"])/arbitraryFactor,
                      (crossSection_theory["cs"]+crossSection_theory["cs_up_err"])/arbitraryFactor,
                      color="mediumseagreen",alpha=0.6,zorder=0)

if barlow_csv!='':
   drawRectangle((0.73,0.42), axes)
axes.tick_params(axis='both',labelsize=24)
axes.ticklabel_format(style='plain')
axes.legend(prop={"size":16},loc='upper right')
axes.set_xlim(0.1,1.0)
axes.set_ylim(0,0.8)
plt.savefig(ofolder+"a2_xsec_clas."+ftype)

############################################################
############## PARITY ASYMMETRY #####################
############################################################
fig_pa,axes_pa=plt.subplots(1,1,figsize=(12,6))

for j,paerr,ebarwidth,alpha,capsize,capwidth in zip(range(2),[parityAsymErr,parityAsymTotalErr],
        [stat_ewidth,tot_ewidth],[stat_alpha,tot_alpha],[stat_capsize,tot_capsize],[stat_capwidth,tot_capwidth]):
    if barlow_csv=='' and j==1:
        continue
    condition=(j==1 and barlow_csv!='') or (j==0 and barlow_csv=='')
    label="GlueX Phase 1\n$E_{\gamma}\in$ [8.2,8.8] GeV" if condition else ''
    axes_pa.errorbar(midts_output,parityAsym,fmt='o',xerr=0,yerr=paerr,c='black',
            label=label,markersize=markersize,elinewidth=ebarwidth,alpha=alpha,capsize=capsize,markeredgewidth=capwidth)
axes_pa.tick_params(axis='both',labelsize=24)
axes_pa.set_xlabel(r"-t $GeV^2$",size=30)
axes_pa.set_ylabel(r"$P_{\sigma}$",size=30)
axes_pa.plot(crossSection_theory["t"],crossSection_theory["pa"],c="mediumseagreen",linewidth=4,label="TMD predictions \n$E_{\gamma}$=8.5 GeV",zorder=0)
axes_pa.fill_between(crossSection_theory["t"],
                      (crossSection_theory["pa"]-crossSection_theory["pa_low_err"])/arbitraryFactor,
                      (crossSection_theory["pa"]+crossSection_theory["pa_up_err"])/arbitraryFactor,
                      color="mediumseagreen",alpha=0.6,zorder=0)

if barlow_csv!='':
   drawRectangle((0.3,0.15), axes_pa)
axes_pa.tick_params(axis='both',labelsize=24)
axes_pa.set_xlim(0.1)
axes_pa.set_ylim(-1,1)
axes_pa.legend(prop={"size":16},loc='lower left')
axes_pa.axhline(0,c='gray',alpha=0.8,linestyle='-',zorder=-1)

plt.tight_layout()

plt.savefig(ofolder+"a2_parityAsym."+ftype)


odf_xsec['-t']=midts_output
odf_xsec[r'$-t~\mathrm{GeV}^2$']=tLabels_output
#odf['tWidth']=tBinWidths_output
odf_xsec[r'$\frac{d\sigma^+}{dt}~[\frac{\mu b}{\mathrm{GeV}^2}]$']=[f'${v:0.4f}\pm{stat:0.4f}\pm{syst:0.04f}$' for v,stat,syst in zip(crossSectionPos,crossSectionPosErr,crossSectionPosSystErr)]
#odf_xsec['xsecpos_staterr']=crossSectionPosErr
#odf_xsec['xsecpos_systerr']=crossSectionPosSystErr
#odf_xsec['xsecpos_totalerr']=crossSectionPosTotalErr
odf_xsec[r'$\frac{d\sigma^-}{dt}~[\frac{\mu b}{\mathrm{GeV}^2}]$']=[f'${v:0.4f}\pm{stat:0.4f}\pm{syst:0.04f}$' for v,stat,syst in zip(crossSectionNeg,crossSectionNegErr,crossSectionNegSystErr)]
#odf_xsec['xsecneg_staterr']=crossSectionNegErr
#odf_xsec['xsecneg_systerr']=crossSectionNegSystErr
#odf_xsec['xsecneg_totalerr']=crossSectionNegTotalErr
odf_xsec[r'$\frac{d\sigma}{dt}~[\frac{\mu b}{\mathrm{GeV}^2}]$']=[f'${v:0.4f}\pm{stat:0.4f}\pm{syst:0.04f}$' for v,stat,syst in zip(crossSection,crossSectionErr,crossSectionSystErr)]
#odf_xsec['xsec_staterr']=crossSectionErr
#odf_xsec['xsec_systerr']=crossSectionSystErr
#odf_xsec['xsec_totalerr']=crossSectionTotalErr
odf_xsec=pd.DataFrame(odf_xsec)

odf_pa[r'$-t~\mathrm{GeV}^2$']=tLabels_output
odf_pa[r'$P_{\sigma}$']=[f'${v:0.4f}\pm{stat:0.4f}\pm{syst:0.04f}$' for v,stat,syst in zip(parityAsym,parityAsymErr,parityAsymSystErr)]
#odf_pa['pa_staterr']=parityAsymErr
#odf_pa['pa_totalerr']=parityAsymTotalErr
odf_pa[r'$L~[\mathrm{pb}^{-1}]$']=[phase1_polarized_lumi]*len(midts_output)
odf_pa[r'$\Gamma~[\%]$']=[branching_frac_of_decay_chain]*len(midts_output)
odf_pa=pd.DataFrame(odf_pa)

#pd.options.display.max_rows = len(odf)
pd.set_option('display.width', 500)
pd.set_option("max_colwidth", 1000)
#pd.options.display.max_columns = len(odf.columns)
with open(f'{ofolder}/results_xsecs.txt', 'w') as f:
    print(odf_xsec.to_latex(escape=False,index=False),file=f)
with open(f'{ofolder}/results_pa.txt', 'w') as f:
    print(odf_pa.to_latex(escape=False,index=False),file=f)
#pd.DataFrame({'pos_staterr':crossSectionPosErr,'neg_staterr':crossSectionNegErr,'total_staterr':crossSectionErr}).to_csv(f'{ofolder}/results_staterrs.csv')
results_df={ # The output of this is used for integrate_xsec.py
        'xsecpos_staterr':crossSectionPosErr,
        'xsecneg_staterr':crossSectionNegErr,
        'xsecpos_systerr':crossSectionPosSystErr,
        'xsecneg_systerr':crossSectionNegSystErr,
        'xsec_staterr':crossSectionErr,
        'xsec_systerr':crossSectionSystErr,
        'xsecpos':crossSectionPos,
        'xsecneg':crossSectionNeg,
        'xsec':crossSection,
        't':midts_output,
        'tWidth':tBinWidths_output,
        }
pd.DataFrame(results_df).to_csv(f'{ofolder}/results.csv')



if orderby!='':
    print(f"Ordering the fit results by {orderby}. The results in each t-bin will be slightly shifted so that left-right is ordered accordingly")
else:
    print(f"No ordering of fit results is requested. It will be ordered according to how os.listdir() chooses")

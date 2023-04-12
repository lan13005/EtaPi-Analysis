#!/usr/bin/python3

import sys
from helper import *
from matplotlib.patches import ConnectionPatch


#################
# SETUP - LOAD DATA
#################
parser = argparse.ArgumentParser(description='Draw diagnostic plots')
parser.add_argument('md_floc', help='folder location containing drawn results root files from PlotGen scripts - for mass dependent fits')
parser.add_argument('bootstrap_folder', default='', help='folder containing bootstrap results')
parser.add_argument('otag', help='output folder tag to drop results in. Folder will have format systematic_[tag]')
parser.add_argument('variations', help='variations that are tested')
parser.add_argument('classes', help='classes, how to group the variations for maximum deviation calculation')
parser.add_argument('splits', help='# variations in each group')
parser.add_argument('keepOnlyConvergedFits', nargs='?', default='True', help='variations that are tested')
args = parser.parse_args()
md_floc=args.md_floc
bootstrap_folder=args.bootstrap_folder
otag=args.otag
variations=args.variations
variations=variations.split(";") if variations != '' else [] 
classes=args.classes
classes=classes.split(";") if classes != '' else [] 
splits=args.splits
splits=splits.split(";") if splits != '' else [] 
splits=[int(s) for s in splits]
keepOnlyConvergedFits=args.keepOnlyConvergedFits=='True'
print(f'keepOnlyConvergedFits: {keepOnlyConvergedFits}')
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

returnYields=True
result=getCrossSection(md_floc,'',bootstrap_folder,orderby='iteration', keepOnlyConvergedFits=keepOnlyConvergedFits,variations=variations,
        dump=ofolder+"input_data.pkl",returnYields=returnYields)
if returnYields:
    xsecpos, xsecneg, xsec, pa, midts_output, tBinWidths_output, waveInts_ts, waveIntErrs_ts, totals, finalStatuses, additionalInfo, variationYields = result
else:
    xsecpos, xsecneg, xsec, pa, midts_output, tBinWidths_output, waveInts_ts, waveIntErrs_ts, totals, finalStatuses, additionalInfo = result
crossSectionPos, crossSectionPosErr, crossSectionPosSystErr, crossSectionPosTotalErr = xsecpos
crossSectionNeg, crossSectionNegErr, crossSectionNegSystErr, crossSectionNegTotalErr = xsecneg
crossSection, crossSectionErr, crossSectionSystErr, crossSectionTotalErr = xsec
parityAsym, parityAsymErr, parityAsymSystErr, parityAsymTotalErr = pa
# tLabels_output, phase1_polarized_lumi, branching_frac_of_decay_chain = additionalInfo # not used for systematics plots

lowts=["0.1","0.2",'0.325','0.5','0.75']
upts=["0.2",'0.325','0.5','0.75','1.0']

############################################################
################### GLUEX SYS PLOT  ########################
###########################################################
nvariations=int(len(crossSection)/len(ts))-1 # minus 1 since the first will always be the nominal!
print(f"len xsec array: {len(crossSection)}")
print(f"number of variations: {nvariations}")
print(f'{variations}')
print(f'{len(variations)}')
assert nvariations==len(variations)
#print(crossSectionPos.reshape(-1,nvariations+1))
#print(crossSectionNeg.reshape(-1,nvariations+1))

fig,axes=plt.subplots(1,5,figsize=(20,int(0.5*nvariations)+3),sharey=True)
fig.subplots_adjust(bottom=0.07) # 0.3 might be good for small number of variations # plots might not fit properly if nvariations is small, we can make manual adjustments...
axes=axes.flatten()

fig2s=[]
ax2s=[]
ms=4
fig2_bottom=0.15 # 0.25 # for a small number of variations
for i in range(5):
    fig2,ax2=plt.subplots(1,ms,figsize=(20,int(0.5*nvariations)+3))
    fig2.subplots_adjust(bottom=fig2_bottom,left=0.15)
    fig2s.append(fig2)
    ax2s.append(ax2)

cpos='orangered'
cneg='royalblue'

#xlims=(-1.5,1.5)
xlims=(-0.5,0.5)
#xlims=(-0.25,0.25)
for it in range(ntbins):
    for k,ax in enumerate([axes[it], ax2s[it][0]]):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xlim(*xlims)
        ax.set_ylim(-0.5,nvariations+0.5)
        ax.axvline(0,c='black',linewidth=1)

        ax.errorbar(x=0,y=0.1,yerr=0,xerr=crossSectionPosErr[it*(nvariations+1)]/crossSectionPos[it*(nvariations+1)],c=cpos,capsize=4)
        ax.errorbar(x=0,y=-0.1,yerr=0,xerr=crossSectionNegErr[it*(nvariations+1)]/crossSectionNeg[it*(nvariations+1)],c=cneg,capsize=4)

        ax.set_title(f"{lowts[it]}< $-t$ <{upts[it]} $GeV^2$")
        ax.minorticks_off()
        ax.tick_params(axis=u'both', which=u'both',length=0)

        #ax.set_xticks(np.arange(-0.6,0.6,0.3))
        #ax.text(0.3, -0.6, r'$\sigma_{pos}$'+f'={crossSectionPos[it*(nvariations+1)]:0.3f}nb'+'\n'+
        #                r'$\sigma_{neg}$'+f'={crossSectionNeg[it*(nvariations+1)]:0.3f}nb', horizontalalignment='center',verticalalignment='center',
        #                size=12)
        ax.text(0.35, -0.1, r'$\sigma_{pos}$'+f'={crossSectionPos[it*(nvariations+1)]:0.3f}nb'+'\n'+
                        r'$\sigma_{neg}$'+f'={crossSectionNeg[it*(nvariations+1)]:0.3f}nb', horizontalalignment='center',verticalalignment='center',
                        size=14)
        ax.set_yticks(range(nvariations+1))
        labels=['Nominal Stat.']+variations

        #yieldSavedInLogs=sum(variationYields)!=0 
        #if yieldSavedInLogs:
        #    for iv in range(nvariations+1):
        #        for it in range(ntbins):
        #            i=it*(nvariations+1)
        #            labels[iv]+=f" | {100*float(variationYields[i+iv])/variationYields[i]:>3.0f}%"
        ax.set_yticklabels(labels)
        if k==0:
            ax.set_xlabel(r"$\frac{\sigma_{variation}-\sigma_{nominal}}{\sigma_{nominal}}$",size=26)
        else:
            ax.set_xlabel(r"$(\sigma_{variation}-\sigma_{nominal})/\sigma_{nominal}$",size=24)

#axes[2].set_xlabel(r"$\frac{\sigma_{variation}-\sigma_{nominal}}{\sigma_{nominal}}$")

mapStatuses={0:'C',1:'A',2:'F'}
tCol=[]
variationCol=[]
classCol=[]
barlowPosCol=[]
diffPosCol=[]
percDiffPosCol=[]
posNomCol=[]
posVarCol=[]
barlowNegCol=[]
diffNegCol=[]
percDiffNegCol=[]
negNomCol=[]
negVarCol=[]
statusCol=[]
yieldsCol=[]
drawXoverFailedFits=True
drawLineSplittingVariationClasses=True

def drawSplitting(_axes,_splits=splits):
    print("Drawing Splitting")
    _xlims=(_axes[0].get_xlim()[0], _axes[-1].get_xlim()[1])
    _ylims0=_axes[0].get_ylim()
    _ylims1=_axes[-1].get_ylim()
    scaleFactor=(_ylims1[1]-_ylims1[0])/(_ylims0[1]-_ylims0[0])
    for split in _splits[:-1]: # dont actually draw the last one since it is at the top of the plot already
        con = ConnectionPatch(
            xyA=(_xlims[0],split+0.5), 
            xyB=(_xlims[1],(split+0.5)*scaleFactor), 
            coordsA="data", coordsB="data",axesA=_axes[0], axesB=_axes[-1], color="black")
        _axes[-1].add_artist(con)

if drawLineSplittingVariationClasses:
    drawSplitting(axes,splits)
    #for it in range(5):
    #    drawSplitting(ax2s[it],splits) # ax2s[it] is a single plot, need to make it into a list

for it in range(ntbins):
    for k, ax in enumerate([axes[it],ax2s[it][0]]):
        for xsec, xsec_err, j in zip([crossSectionPos,crossSectionNeg],[crossSectionPosErr,crossSectionNegErr],range(2)):
            nom=np.array([xsec[(nvariations+1)*it] for iv in range(nvariations)])
            var=np.array([xsec[(nvariations+1)*it+(iv+1)] for iv in range(nvariations)])
            yields=np.array([variationYields[(nvariations+1)*it+(iv+1)] for iv in range(nvariations)])
            #perc_diff=np.array([(xsec[(nvariations+1)*it+(iv+1)]-xsec[(nvariations+1)*it])/xsec[(nvariations+1)*it] for iv in range(nvariations)])
            perc_diff=(var-nom)/nom
            diff=perc_diff*xsec[(nvariations+1)*it]
            ## Option 1: Barlow Errors sqrt(sx**2-sy**2)
            errs=np.array([np.sqrt(abs(xsec_err[(nvariations+1)*it]**2-xsec_err[(nvariations+1)*it+(iv+1)]**2))/xsec[(nvariations+1)*it] for iv in range(nvariations)])
            ## Option 2: Variation Errors
            #errs=np.array([abs(xsec_err[(nvariations+1)*it+(iv+1)])/xsec[(nvariations+1)*it] for iv in range(nvariations)])
            errs[abs(errs)<0.0001]=0.0001
            nbarlows=abs(perc_diff/errs)
            #errs=[np.sqrt(abs(xsec_err[(nvariations+1)*it]**2-xsec_err[(nvariations+1)*it+(iv+1)]**2)) for iv in range(nvariations)]
            y=np.arange(1,nvariations+1)+0.05 if j==0 else np.arange(1,nvariations+1)-0.05
            c=cpos if j==0 else cneg
            ax.errorbar(perc_diff,y,xerr=errs,yerr=0,c=c,linewidth=6,ls='none')

            finalStatus=finalStatuses[(nvariations+1)*it:(nvariations+1)*(it+1)] 
            ystatus=np.arange(len(finalStatus))[finalStatus!=0]
            if drawXoverFailedFits:
                ax.scatter([0]*len(ystatus),ystatus,marker='x',s=90,c='black',zorder=9)

            if k==0: # k is for the axis drawer. We only want to include it once.
                if j==0:
                    tCol.extend([ts[it]]*nvariations)
                    variationCol.extend(variations)
                    classCol.extend(classes)
                    barlowPosCol.extend(nbarlows)
                    diffPosCol.extend(diff)
                    percDiffPosCol.extend(perc_diff)
                    posNomCol.extend(nom)
                    posVarCol.extend(var)
                    statusCol.extend([mapStatuses[s] for s in finalStatus[1:]]) # start from 1 since we want to skip the nominal
                    yieldsCol.extend(yields)
                else:
                    barlowNegCol.extend(nbarlows)
                    diffNegCol.extend(diff)
                    percDiffNegCol.extend(perc_diff)
                    negNomCol.extend(nom)
                    negVarCol.extend(var)

def format(col):
    return [f'{x:0.05f}' for x in col]

#### Draw the plot
fig.tight_layout()
fig.savefig(ofolder+"systOverview."+ftype)

#### Nominal output scheme 
map_t={k:v for k,v in zip(ts,midts)}

barlowDF={}
barlowDF['tbin']=tCol
barlowDF['t']=format([map_t[t] for t in tCol])
barlowDF['variation']=variationCol
barlowDF['class']=classCol
barlowDF['posNom']=format(posNomCol)
barlowDF['negNom']=format(negNomCol)
barlowDF['posVar']=format(posVarCol)
barlowDF['negVar']=format(negVarCol)
barlowDF['nBarlowPos']=format(barlowPosCol)
barlowDF['diffPos']=format(diffPosCol)
barlowDF['percDiffPos']=format(percDiffPosCol)
barlowDF['nBarlowNeg']=format(barlowNegCol)
barlowDF['diffNeg']=format(diffNegCol)
barlowDF['percDiffNeg']=format(percDiffNegCol)
barlowDF['yields']=format(yieldsCol)
barlowDF['status']=statusCol
for k,v in barlowDF.items():
    print(f'{k} {len(v)}')

barlowDF=pd.DataFrame(barlowDF)
barlowDF.to_csv(ofolder+"barlows.csv",float_format='0.5f',index=False)

### Alternative output scheme that is more readable for analysis notes
map_trange={
        '010020':'0.1-0.2',
        '0200325':'0.2-0.325',
        '0325050':'0.325-0.5',
        '050075':'0.5-0.75',
        '075100':'0.75-1.0',
        }
barlowDF={}
barlowDF['trange']=[map_trange[t] for t in tCol]
#barlowDF['t']=[map_t[t] for t in tCol]
barlowDF['Variation']=variationCol
barlowDF['(+)Nom.']=format(posNomCol)
barlowDF['(+)\%Diff']=format(percDiffPosCol)
barlowDF['pNbarlow']=barlowPosCol
barlowDF['(-)Nom.']=format(negNomCol)
barlowDF['(-)\%Diff']=format(percDiffNegCol)
barlowDF['nNbarlow']=barlowNegCol
#barlowDF['dataYield']=format(yieldsCol)
#barlowDF['status']=statusCol
barlowDF=pd.DataFrame(barlowDF)
nThreshold=4
#barlowDF=barlowDF[(barlowDF.pNbarlow.abs()>nThreshold)|(barlowDF.nNbarlow.abs()>nThreshold)].reset_index(drop=True)
for bdf,tag in zip([barlowDF[(barlowDF.pNbarlow.abs()>nThreshold)].reset_index(drop=True), barlowDF[(barlowDF.nNbarlow.abs()>nThreshold)].reset_index(drop=True)],
        ['pos','neg']):
    bdf.pNbarlow=format(bdf.pNbarlow)
    bdf.nNbarlow=format(bdf.nNbarlow)
    bdf.rename(columns={'pNbarlow':'(+)N_{barlow}','nNbarlow':'(-)N_{barlow}'})
    #### OUTPUT ALTERNATIVE FORM AND LATEX TABLE
    bdf.to_csv(ofolder+f"barlows_altform_{tag}.csv",float_format='0.5f',index=False)
    stdoutOrigin=sys.stdout 
    sys.stdout = open(ofolder+f"latex_table_{tag}", "w")
    print(bdf.to_latex(float_format='0.3f',index=False,escape=False))
    sys.stdout.close()
    sys.stdout=stdoutOrigin


############################################################
############## GLUEX WAVE FRACTION PLOT ####################
############################################################

coherent_sums={
        "D+":set(['D2-+','D1-+','D0++','D1++','D2++']),
        "D-":set(['D2--','D1--','D0+-','D1+-','D2+-']),
        "pD+":set(['pD2-+','pD1-+','pD0++','pD1++','pD2++']),
        "pD-":set(['pD2--','pD1--','pD0+-','pD1+-','pD2+-']),
        "allD+":set(['D2-+','D1-+','D0++','D1++','D2++','pD2-+','pD1-+','pD0++','pD1++','pD2++']),
        "allD-":set(['D2--','D1--','D0+-','D1+-','D2+-','pD2--','pD1--','pD0+-','pD1+-','pD2+-']),
        "allD2-+":set(['D2-+','pD2-+']),
        "allD1-+":set(['D1-+','pD1-+']),
        "allD0++":set(['D0++','pD0++']),
        "allD1++":set(['D1++','pD1++']),
        "allD2++":set(['D2++','pD2++']),
        "allD2--":set(['D2--','pD2--']),
        "allD1--":set(['D1--','pD1--']),
        "allD0+-":set(['D0+-','pD0+-']),
        "allD1+-":set(['D1+-','pD1+-']),
        "allD2+-":set(['D2+-','pD2+-'])
#        "summed":set(["S0++","S0+-"]),
        }
#for wave in ["D+",'D-','pD+','pD-']:
#    coherent_sums["summed"].update(coherent_sums[wave])

incoherent_sum_totals={}

map_cls_multiple={
        'interference_D':{
        'D+':[cpos,'solid',r'$D^{+}$'],
        'D-':[cneg,'solid',r'$D^{-}$']},
        'interference_pD':{
        'pD+':[cpos,'solid',r'$pD^{+}$'],
        'pD-':[cneg,'solid',r'$pD^{-}$']},
        'interference_allD':{
        'allD+':[cpos,'solid',r'$D^{+}$'],
        'allD-':[cneg,'solid',r'$D^{-}$']},
        'interference_D2m':{
        'allD2-+':[cpos,'solid',r'$D_{-2}^{+}$'],
        'allD2--':[cneg,'solid',r'$D_{-2}^{-}$']},
        'interference_D1m':{
        'allD1-+':[cpos,'solid',r'$D_{-1}^{+}$'],
        'allD1--':[cneg,'solid',r'$D_{-1}^{-}$']},
        'interference_D0p':{
        'allD0++':[cpos,'solid',r'$D_{0}^{+}$'],
        'allD0+-':[cneg,'solid',r'$D_{0}^{-}$']},
        'interference_D1p':{
        'allD1++':[cpos,'solid',r'$D_{1}^{+}$'],
        'allD1+-':[cneg,'solid',r'$D_{1}^{-}$']},
        'interference_D2p':{
        'allD2++':[cpos,'solid',r'$D_{2}^{+}$'],
        'allD2+-':[cneg,'solid',r'$D_{2}^{-}$']},
        'interference_allD_individual':{
        'allD1--':[cneg,'dotted',r'$D_{-1}^{-}$'],
        'allD0++':[cpos,'dotted',r'$D_{0}^{+}$'],
        'allD0+-':[cneg,'dashed',r'$D_{0}^{-}$'],
        'allD1++':[cpos,'dashed',r'$D_{1}^{+}$'],
        'allD1+-':[cneg,'solid',r'$D_{1}^{-}$'],
        'allD2++':[cpos,'solid',r'$D_{2}^{+}$']},
        'invidualSx03_D':{
        'D1--':[cneg,'dotted',r'$D_{-1}^{-}$'],
        'D0++':[cpos,'dotted',r'$D_{0}^{+}$'],
        'D0+-':[cneg,'dashed',r'$D_{0}^{-}$'],
        'D1++':[cpos,'dashed',r'$D_{1}^{+}$'],
        'D1+-':[cneg,'solid',r'$D_{1}^{-}$'],
        'D2++':[cpos,'solid',r'$D_{2}^{+}$'],
        'S0+-':['black','dashed',r'$S_{0}^{-}$ x0.3'],
        'S0++':['black','solid',r'$S_{0}^{+}$ x0.3']},
        'invidualSx03_pD':{
        'pD1--':[cneg,'dotted',r'$D_{-1}^{-}$'],
        'pD0++':[cpos,'dotted',r'$D_{0}^{+}$'],
        'pD0+-':[cneg,'dashed',r'$D_{0}^{-}$'],
        'pD1++':[cpos,'dashed',r'$D_{1}^{+}$'],
        'pD1+-':[cneg,'solid',r'$D_{1}^{-}$'],
        'pD2++':[cpos,'solid',r'$D_{2}^{+}$'],
        'S0+-':['black','dashed',r'$S_{0}^{-}$ x0.3'],
        'S0++':['black','solid',r'$S_{0}^{+}$ x0.3']},
        'invidualSx03_allD':{
        'allD1--':[cneg,'dotted',r'$D_{-1}^{-}$'],
        'allD0++':[cpos,'dotted',r'$D_{0}^{+}$'],
        'allD0+-':[cneg,'dashed',r'$D_{0}^{-}$'],
        'allD1++':[cpos,'dashed',r'$D_{1}^{+}$'],
        'allD1+-':[cneg,'solid',r'$D_{1}^{-}$'],
        'allD2++':[cpos,'solid',r'$D_{2}^{+}$'],
        'S0+-':['black','dashed',r'$S_{0}^{-}$ x0.3'],
        'S0++':['black','solid',r'$S_{0}^{+}$ x0.3']},
        #'invidualSx03':{
        #'D0++':[cpos,'dotted',r'$D_{0}^{+}$'],
        #'D0+-':[cneg,'dashed',r'$D_{0}^{-}$'],
        #'D2++':[cpos,'dotted',r'$D_{2}^{+}$'],
        #'D2+-':[cneg,'dashed',r'$D_{2}^{-}$'],
        #'S0+-':['black','dashed',r'$S_{0}^{-}$ x0.3'],
        #'S0++':['black','solid',r'$S_{0}^{+}$ x0.3']}
        }

axisLocations=[-1,-1,-1,-1,-1,-1,-1,-1,3,1,-1,2] # this value relates to the axis you want to draw the values of map_cls_multiple into.
for i,nameTag2 in enumerate(map_cls_multiple.keys()):
    map_cls=map_cls_multiple[nameTag2]
    m=axisLocations[i]
    for nameTag in ["FitFraction"]:#,"CorrectedYield"]:
        fig,axes=plt.subplots(1,5,figsize=(20,int(0.5*nvariations)+2),sharey=True)
        axes=axes.flatten()
        for it in range(ntbins):
            axisList=[axes[it],ax2s[it][m]] if m!=-1 else [axes[it]]
            for iax,ax in enumerate(axisList):
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.set_ylim(-0.5,nvariations+0.5)
                ax.axvline(0,c='black',linewidth=1)
        
                ax.set_title(f"{lowts[it]}< $-t$ <{upts[it]} $GeV^2$")
                ax.minorticks_off()
                ax.tick_params(axis=u'both', which=u'both',length=0)
        
                ax.set_yticks(range(nvariations+1))
                if iax==0:
                    ax.set_yticklabels(['Nominal']+variations)
                else:
                    ax.set_yticklabels(['']*(nvariations+1))
                #xlabel=r'(Incoh. - coh.)/coh.' if nameTag2.startswith("interference") else nameTag
                xlabel=nameTag if not nameTag2.startswith('interference') else f'Interference {nameTag}'
                ax.set_xlabel(xlabel,size=18)

        y=np.arange(0,nvariations+1)
        for it in range(ntbins):
            axisList=[axes[it],ax2s[it][m]] if m!=-1 else [axes[it]]
            for ax in axisList:
                for wave in map_cls.keys():
                    scaleFactor=1 if nameTag=='FitFraction' else totals
                    c,l,label = map_cls[wave]
                    for k in coherent_sums.keys(): # reset totals for the next t-bin
                        incoherent_sum_totals[k]=np.zeros(nvariations+1)
                    if nameTag2.startswith('interference'):
                        print(f'sum {wave} with sum at: {incoherent_sum_totals[wave]}')
                        for wave2 in coherent_sums[wave]: # loop over all waves that should have belonged in a coherent sum, i.e. D+: {D0++,D1++...}
                            print(f'  adding {wave2}')
                            dfrac=waveInts_ts[wave2]*scaleFactor
                            dfrac=dfrac[(nvariations+1)*it:(nvariations+1)*(it+1)]
                            incoherent_sum_totals[wave]+=dfrac
                            print(f'  current totals: {incoherent_sum_totals[wave]}')
                            #print(f'Adding {wave2} results in total {incoherent_sum_totals}')
                        dfrac=waveInts_ts[wave]*scaleFactor
                        dfrac=dfrac[(nvariations+1)*it:(nvariations+1)*(it+1)]
                        print(f'  final totals: {incoherent_sum_totals[wave]}')
                        print(f'  coherent totals: {dfrac}')
                        #ax.axvline(dfrac[0],c=c,linewidth=1,alpha=1.0) # draw coherent sum value as a vertical line
                        diff=[incoherent_sum_totals[wave][ivary]-dfrac[ivary] for ivary in range(nvariations+1)] # Difference
                        #diff=[diff[ivary]/dfrac[ivary] for ivary in range(nvariations+1)] # Percent Difference
                        ax.plot(diff,y,c=c,linestyle=l,linewidth=2,label=label,marker='o',markersize=5) # draw difference of incoherent sum and coherent
                    else:
                        dfrac=waveInts_ts[wave]*scaleFactor
                        dfrac=dfrac[(nvariations+1)*it:(nvariations+1)*(it+1)]
                        if wave[0]=="S": # if it is of the form invidualSx03_allD then we will scale by what is after the x
                            try:
                                scaleFactor=float(label.split('x')[-1])
                            except:
                                scaleFactor=1
                            dfrac=dfrac*scaleFactor
                        ax.plot(dfrac,y,c=c,linestyle=l,linewidth=2,label=label,marker='o',markersize=5)
        ###########################################################################################################################################################
        ##### For +wave systematic scan. There exist potential waves that are not in the nominal waveset that are added for the 
        #####    waveset scan. These waves will only have ntbins number of values > 0. 
        #####    ** All this does is draw a limegreen +sign on the fit fraction **
        if not nameTag2.startswith('interference'): # summed waves do not need modification
            nWavesIncluded=0
            for wave2,values2 in waveInts_ts.items():
                # if wave2 is not explicitly requested to be drawn but there are fits that include that wave (has associated fit fraction)
                #    wave2 should look for waves with format: {D2++, pD2++, allD2++}
                wavePrefix=nameTag2.split('_')[-1]
                bWavePrefix=wave2.startswith(wavePrefix) # make sure wave2 starts with [D, pD, allD]
                bWaveFormat=len(wave2)==3+len(wavePrefix) and wave2[-1] in ['+','-'] # make sure wave2 has 4 characters where the last character is for reflectivity
                bSingleVariation=sum(values2>0)==ntbins # only 1 variation uses this wave (for +wave systematic scan)
                bWaveNotUsed=wave2 not in map_cls.keys() # make sure the wave is not part of the nominal drawn waveset
                if bWaveNotUsed and bSingleVariation and bWaveFormat and bWavePrefix:
                    #print(f"DRAWING {wave2}!")
                    yloc=np.argwhere(values2>0)[0] # values2 = [variation x tbin] where first elements are for variations performed in first tbin
                    scaleFactor=1 if nameTag=='FitFraction' else totals
                    for it in range(ntbins):
                        axisList=[axes[it],ax2s[it][m]] if m!=-1 else [axes[it]]
                        for ax in axisList:
                            dfrac=waveInts_ts[wave2]*scaleFactor
                            dfrac=dfrac[(nvariations+1)*it+yloc]
                            l='+wave' if nWavesIncluded==0 else ''
                            ax.scatter(dfrac,yloc,c='limegreen',marker='+',label=l,s=100,zorder=9) # large zorder to put on top
                            #print(f'tBin{it} for {wave2} with frac {dfrac} drawn with ylocation {yloc}')
                    nWavesIncluded+=1
        ###########################################################################################################################################################
        # only want to draw for the first t-bin on axes but for all individually separated ax2s plots
        for j,ax in enumerate([axes[0]]+[ax2s[it][m] for it in range(ntbins)]): 
            handles, labels = ax.get_legend_handles_labels()
            if j==0:
                fig.legend(handles, labels, ncol=5,prop={'size':12}, loc='lower left', frameon=True, fancybox=True, facecolor='white', 
                            framealpha=0.9, borderpad=0.1, edgecolor='black',)
            else:
                if m not in [2,3]: # only draw the legend on the first fitfraction plot. m is the index into the axis
                    continue
                loc=(1.0,-fig2_bottom+0.05) if m==2 else (0.8,-0.21) # adjust xlocation of legend
                #size=12 if m==2 else 18
                size=18 if m==2 else 24
                ax.legend(handles, labels, ncol=5,prop={'size':size}, frameon=True, fancybox=True, facecolor='white', 
                            framealpha=0.9, borderpad=0.1, edgecolor='black', bbox_to_anchor=loc)

        if drawLineSplittingVariationClasses:
            drawSplitting(axes,splits)
            for it in range(5):
                drawSplitting(ax2s[it],splits) # ax2s[it] is a single plot, need to make it into a list

        for it in range(ntbins):
            for im,mTag in zip(range(1,ms),[r'$a_2$',r"$a_2~and~a_2^{\prime}$",r"$a_2~a_2^{\prime}$"'\n'"Interference"]):
                #ax2s[it][im].set_xlabel(ax2s[it][im].get_xlabel()+mTag)
                ax2s[it][im].set_xlabel(mTag+"\nFitFraction")
    
        fig.tight_layout()
        ftag='fitfrac'
        fig.savefig(ofolder+"syst_"+nameTag2+"_"+nameTag+"_"+ftag+"."+ftype)

for it in range(5):
#    fig2s[it].tight_layout()
    fig2s[it].savefig(ofolder+"barlow_fitfrac_t"+str(it)+"."+ftype)





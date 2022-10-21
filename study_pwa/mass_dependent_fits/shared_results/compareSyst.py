#!/usr/bin/python3

from helper import *

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

xsecpos, xsecneg, xsec, pa, midts_output, waveInts_ts, waveIntErrs_ts = getCrossSection(md_floc,'',orderby='iteration')
crossSectionPos, crossSectionPosErr, crossSectionPosTotalErr = xsecpos
crossSectionNeg, crossSectionNegErr, crossSectionNegTotalErr = xsecneg
crossSection, crossSectionErr, crossSectionTotalErr = xsec
parityAsym, parityAsymErr, parityAsymTotalErr = pa

lowts=["0.1","0.2",'0.325','0.5','0.75']
upts=["0.2",'0.325','0.5','0.75','1.0']

############################################################
################### GLUEX SYS PLOT  ########################
############################################################
nvariations=int(len(crossSection)/len(ts))-1 # minus 1 since the first will always be the nominal!
print(f"len xsec array: {len(crossSection)}")
print(f"number of variations: {nvariations}")
assert nvariations==len(variations)

fig,axes=plt.subplots(1,5,figsize=(16,int(0.5*nvariations)+2),sharey=True)
axes=axes.flatten()

cpos='orangered'
cneg='royalblue'

for iax,ax in enumerate(axes):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(-0.6,0.6)
#    ax.set_xlim(-0.3,0.3)
    ax.set_ylim(-0.5,nvariations+0.5)
    ax.axvline(0,c='black',linewidth=1)

    ax.errorbar(x=0,y=0.1,yerr=0,xerr=crossSectionPosErr[iax*(nvariations+1)]/crossSectionPos[iax*(nvariations+1)],c=cpos,capsize=4)
    ax.errorbar(x=0,y=-0.1,yerr=0,xerr=crossSectionNegErr[iax*(nvariations+1)]/crossSectionNeg[iax*(nvariations+1)],c=cneg,capsize=4)

    ax.set_title(f"{lowts[iax]}< $-t$ <{upts[iax]} $GeV^2$")
    ax.minorticks_off()
    ax.tick_params(axis=u'both', which=u'both',length=0)

#    ax.set_xticks(np.arange(-0.6,0.6,0.3))
    #ax.text(0.3, -0.6, r'$\sigma_{pos}$'+f'={crossSectionPos[iax*(nvariations+1)]:0.3f}nb'+'\n'+
    #                r'$\sigma_{neg}$'+f'={crossSectionNeg[iax*(nvariations+1)]:0.3f}nb', horizontalalignment='center',verticalalignment='center',
    #                size=12)
    ax.text(0.35, -0.1, r'$\sigma_{pos}$'+f'={crossSectionPos[iax*(nvariations+1)]:0.3f}nb'+'\n'+
                    r'$\sigma_{neg}$'+f'={crossSectionNeg[iax*(nvariations+1)]:0.3f}nb', horizontalalignment='center',verticalalignment='center',
                    size=10)
    ax.set_yticks(range(nvariations+1))
    ax.set_yticklabels(['Nominal Stat.']+variations)

axes[2].set_xlabel(r"$\frac{\sigma_{variation}-\sigma_{nominal}}{\sigma_{nominal}}$")

barlowDF={}
tCol=[]
variationCol=[]
barlowPosCol=[]
diffPosCol=[]
barlowNegCol=[]
diffNegCol=[]
for it in range(ntbins):
    for xsec, xsec_err, j in zip([crossSectionPos,crossSectionNeg],[crossSectionPosErr,crossSectionNegErr],range(2)):
        perc_diff=np.array([(xsec[(nvariations+1)*it+(iv+1)]-xsec[(nvariations+1)*it])/xsec[(nvariations+1)*it] for iv in range(nvariations)])
        diff=perc_diff*xsec[(nvariations+1)*it]
        errs=np.array([np.sqrt(abs(xsec_err[(nvariations+1)*it]**2-xsec_err[(nvariations+1)*it+(iv+1)]**2))/xsec[(nvariations+1)*it] for iv in range(nvariations)])
        errs[abs(errs)<0.0001]=0.0001
        nbarlows=perc_diff/errs
        #errs=[np.sqrt(abs(xsec_err[(nvariations+1)*it]**2-xsec_err[(nvariations+1)*it+(iv+1)]**2)) for iv in range(nvariations)]
        y=np.arange(1,nvariations+1)+0.05 if j==0 else np.arange(1,nvariations+1)-0.05
        c=cpos if j==0 else cneg
        axes[it].errorbar(perc_diff,y,xerr=errs,yerr=0,c=c,linewidth=6,ls='none')

        if j==0:
            tCol.extend([ts[it]]*nvariations)
            variationCol.extend(variations)
            barlowPosCol.extend(nbarlows)
            diffPosCol.extend(diff)
        else:
            barlowNegCol.extend(nbarlows)
            diffNegCol.extend(diff)

def format(col):
    return [f'{x:0.05f}' for x in col]

barlowDF['t']=tCol
barlowDF['variation']=variationCol
barlowDF['nBarlowPos']=format(barlowPosCol)
barlowDF['diffPos']=format(diffPosCol)
barlowDF['nBarlowNeg']=format(barlowNegCol)
barlowDF['diffNeg']=format(diffNegCol)
barlowDF=pd.DataFrame(barlowDF)

plt.tight_layout()
plt.savefig(ofolder+"systOverview."+ftype)
barlowDF.to_csv(ofolder+"barlows.csv",float_format='0.5f',index=False)
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
    ax.set_ylim(-0.5,nvariations+0.5)
    ax.axvline(0,c='black',linewidth=1)

    ax.set_title(f"{lowts[iax]}< $-t$ <{upts[iax]} $GeV^2$")
    ax.minorticks_off()
    ax.tick_params(axis=u'both', which=u'both',length=0)

    ax.set_yticks(range(nvariations+1))
    ax.set_yticklabels(['Nominal Stat.']+variations)

axes[2].set_xlabel("Fit Fractions",size=14)

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

print("Finished Fit Fraction plot")
for it in range(ntbins):
    for wave in map_cls.keys():
        dfrac=waveInts_ts[wave]
        c,l,label = map_cls[wave]
        dfrac=dfrac[(nvariations+1)*it:(nvariations+1)*(it+1)]
        if wave[0]=="S":
            dfrac=dfrac*0.3
        y=np.arange(0,nvariations+1)
        axes[it].plot(dfrac,y,c=c,linestyle=l,linewidth=2,label=label,marker='o',markersize=3)
handles, labels = axes[0].get_legend_handles_labels()
yleg=(0.19+0.013*2)-0.015*nvariations
bottom=(0.35)-0.015*nvariations
fig.legend(handles, labels, bbox_to_anchor=(0.06, yleg, 0.35, 0.125),ncol=4,prop={'size':14})
#fig.legend(handles, labels, bbox_to_anchor=(0.06, 0.2, 0.35, 0.125), ncol=4,prop={'size':14})

plt.tight_layout()
#fig.subplots_adjust(bottom=0.35) # or whatever
fig.subplots_adjust(bottom=bottom) # or whatever
ftag='fitfrac'
plt.savefig(ofolder+"systDwaveFractions_"+ftag+"."+ftype)

#nvary=14
#fig.subplots_adjust(bottom=0.15)
#bbox_to_anchor=(0.06, 0, 0.35, 0.125)

#nvary=2
#0.35
#0.2









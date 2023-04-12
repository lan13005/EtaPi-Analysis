

#!/usr/bin/python3

##################
# 1. Draws some overview plots that shows the systematic variations we are performing
# 2. 
##################


###################################################################
import pandas as pd
import matplotlib.pyplot as plt
import uproot3 as uproot
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 100
import numpy as np
import seaborn as sns
import os

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
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

###################################################################

def getSelection(df,sectionStr):
    vars1=sectionStr.split(" ")[::3]
    mins1=[float(x) for x in sectionStr.split(" ")[1::3]]
    maxs1=[float(x) for x in sectionStr.split(" ")[2::3]]
    
    select=pd.Series(np.ones(len(df),dtype=bool))
    for var1,min1,max1 in zip(vars1,mins1,maxs1):
        if var1[0]=='!':
            var1=var1[1:]
            select &= ~((df[var1]>min1)&(df[var1]<max1))
        else:
            select &= (df[var1]>min1)&(df[var1]<max1)
    return select


def loadDF(fileName,treeName,cols):
    ''' Function to load a root file '''
    tree=uproot.open(fileName)[treeName]
    df=tree.arrays(cols,outputtype=pd.DataFrame).reset_index(drop=True)
    return df

cols=["unusedEnergy","chiSq",'mmsq','pVH',
      "photonTheta1","photonTheta2","photonTheta3","photonTheta4",
      "photonE1","photonE2","photonE3","photonE4",
      "proton_momentum",'proton_z',
      'mandelstam_t','Mpi0eta','Mpi0','Meta','Mpi0p','Metap',
      'Weight'
     ]

baseDir="rootFiles/"
ts=["010020","0200325","0325050","050075","075100"]
m="104180"
extraTag="_selectGenTandM"
dfs=[]
for t in ts:
    for pol in ["000","045","090","135"]:
        baseLoc=baseDir+"t"+t+"_m"+m+extraTag+"/"
        dataloc="pol"+pol+"_t"+t+"_m"+m+extraTag+"_DTOT_selected_acc_flat.root"
        df=loadDF(baseLoc+dataloc, 'kin', cols)
        dfs.append(df)


df=pd.concat(dfs).reset_index(drop=True)
df=df[getSelection(df,'pVH 0.5 999')].reset_index(drop=True)


fig,axes=plt.subplots(4,3,figsize=(15,18))
axes=axes.flatten()

def drawHist(var,label,ax,weight='Weight',xmin=None,xmax=None,nbins=30):
    tmp=df
    if xmin:
        tmp=df[df[var]>xmin]
    if xmax:
        tmp=df[df[var]<xmax]
    out=np.histogram(tmp[var],weights=tmp[weight],bins=nbins)
    binwidth=out[1][1]-out[1][0]
    hep.histplot(out,ax=ax,c='black',histtype='step')
    ax.set_ylim(0.1)
    ax.set_xlabel(label)
    ax.set_ylabel(f'Events / {binwidth:.3f}')
    if xmax:
        ax.set_xlim(ax.get_xlim()[0],xmax)
    if xmin:
        ax.set_xlim(xmin,ax.get_xlim()[1])

drawHist('Mpi0eta',r'$M(\eta\pi$ [GeV]',axes[0])
drawHist('unusedEnergy','Unused Energy [GeV]',axes[1])
drawHist('chiSq',r'$\chi^2$',axes[2])
drawHist('mmsq',r'$MM^2$ [$GeV^2$]',axes[3],xmin=-0.06,xmax=0.06)
drawHist('proton_z',r'$z_{proton}$ [cm]',axes[4],xmin=49,xmax=81)
drawHist('photonTheta1',r'$\theta_{\gamma_1}$ [degrees]',axes[5],xmax=14,nbins=30)
drawHist('photonTheta3',r'$\theta_{\gamma_3}$ [degrees]',axes[6],xmax=14,nbins=30)
drawHist('photonE1',r'$E_{\gamma_1}$ [GeV]',axes[7],nbins=30,xmin=0.05,xmax=0.5)
drawHist('photonE3',r'$E_{\gamma_3}$ [GeV]',axes[8],nbins=30,xmin=0.05,xmax=0.5)
drawHist('Mpi0p',r'$M(\pi^0p)$ [GeV]',axes[9],nbins=50,xmax=3.6)
    
axes[1].set_yscale('log')
axes[10].axis(False)
axes[11].axis(False)

cnom='red'
cvary='gray'

# Unused Energy
axes[1].axvline(0.01,c=cnom,linestyle='--')
axes[1].axvline(0.14,c=cvary,linestyle='--')
axes[1].axvline(0.17,c=cvary,linestyle='--')

# chiSq
axes[2].axvline(13.277,c=cnom,linestyle='--')
axes[2].axvline(11,c=cvary,linestyle='--')
axes[2].axvline(16,c=cvary,linestyle='--')

# mmsq
axes[3].axvline(-0.05,c=cnom,linestyle='--')
axes[3].axvline(0.05,c=cnom,linestyle='--')
axes[3].axvline(-0.013,c=cvary,linestyle='--')
axes[3].axvline(0.013,c=cvary,linestyle='--')
axes[3].axvline(-0.01,c=cvary,linestyle='--')
axes[3].axvline(0.01,c=cvary,linestyle='--')

# proton z
axes[4].axvline(52,c=cnom,linestyle='--')
axes[4].axvline(78,c=cnom,linestyle='--')
axes[4].axvline(51,c=cvary,linestyle='--')
axes[4].axvline(79,c=cvary,linestyle='--')
axes[4].axvline(53,c=cvary,linestyle='--')
axes[4].axvline(77,c=cvary,linestyle='--')

# photon theta 1/2
axes[5].axvline(2.5,c=cnom,linestyle='--')
axes[5].axvline(10.3,c=cnom,linestyle='--')
axes[5].axvline(11.9,c=cnom,linestyle='--')
axes[5].axvline(2.1,c=cvary,linestyle='--')
axes[5].axvline(2.9,c=cvary,linestyle='--')
axes[5].axvline(10.1,c=cvary,linestyle='--')
axes[5].axvline(12.1,c=cvary,linestyle='--')
axes[5].axvline(10.4,c=cvary,linestyle='--')
axes[5].axvline(11.7,c=cvary,linestyle='--')

# photon theta 3/4
axes[6].axvline(2.5,c=cnom,linestyle='--')
axes[6].axvline(10.3,c=cnom,linestyle='--')
axes[6].axvline(11.9,c=cnom,linestyle='--')
axes[6].axvline(2.1,c=cvary,linestyle='--')
axes[6].axvline(2.9,c=cvary,linestyle='--')
axes[6].axvline(10.1,c=cvary,linestyle='--')
axes[6].axvline(12.1,c=cvary,linestyle='--')
axes[6].axvline(10.4,c=cvary,linestyle='--')
axes[6].axvline(11.7,c=cvary,linestyle='--')

# photon E 1/2
axes[7].axvline(0.1,c=cnom,linestyle='--')
axes[7].axvline(0.12,c=cvary,linestyle='--')
axes[7].axvline(0.13,c=cvary,linestyle='--')

# photon E 3/4
axes[8].axvline(0.1,c=cnom,linestyle='--')
axes[8].axvline(0.12,c=cvary,linestyle='--')
axes[8].axvline(0.13,c=cvary,linestyle='--')

# Mpi0 
axes[9].axvline(1.4,c=cvary,linestyle='--')
axes[9].axvline(1.5,c=cvary,linestyle='--')

plt.tight_layout()

plt.savefig("/d/grid17/ln16/dselector_v3/study_pwa/mass_dependent_fits/systOverview.pdf")



###################################################################################
###################################################################################

## MASS SIDEBAND DRAWER. NEED TO LOAD ANOTHER DATASET SINCE THE DATA IN 
##    PHASE1_SELECTED_V4 ALL HAVE SIDEBAND SKIP REGION REMOVED

###################################################################################
###################################################################################



dfs_looseChi=[]
datasets=["2017_1","2018_1","2018_8"]
for dataset in datasets:
#     fileName="/d/grid17/ln16/myDSelector/zThesisCuts/degALL_data_"+dataset+"_stdProtonPhotonExclusivity_looseChi_for_thesis_treeFlat_DSelector.root"
    fileName="/d/grid17/ln16/dselector_v3/study_eventSelections/rootFiles/looseUE/D"+dataset+"_selected_acc_flat.root"
    treeName="kin"
    cols=["Mpi0eta","weightASBS","weightBS","cosTheta_eta_gj","cosTheta_eta_hel","phi_eta_gj","phi_eta_hel",
          "Mpi0","Meta","AccWeight","chiSq","Ebeam","mandelstam_t","mandelstam_tp","Mpi0p","Metap","omegaCut","unusedEnergy",
          "proton_z",'mmsq',
         "mandelstam_teta","mandelstam_tpi0","run","event","rfTime","phi_eta_lab","phi_pi0_lab","BeamAngle",
         "photonE1","photonE2","photonE3","photonE4","photonTheta1","photonTheta2","photonTheta3","photonTheta4",
         "pVH","pVH_etap","pVH_pi0p"]
    tree=uproot.open(fileName)[treeName]
#     df=tree.arrays(cols,entrystart=0,entrystop=100000,outputtype=pd.DataFrame).reset_index(drop=True)
    df=tree.arrays(cols,outputtype=pd.DataFrame).reset_index(drop=True)
    dfs_looseChi.append(df)
    
yields=[dfs_looseChi[i].weightASBS.sum() for i in range(len(dfs_looseChi))]
yieldsFrac=[val/sum(yields) for val in yields]
for i,val in enumerate(yieldsFrac):
    print(datasets[i]+" integral: %0.2f" % val)


phase1_looseChi=pd.concat(dfs_looseChi)
tmp=phase1_looseChi[
    (phase1_looseChi.Mpi0eta>1.04)&(phase1_looseChi.Mpi0eta<1.72)
    &(phase1_looseChi.mandelstam_t<1.0)&(phase1_looseChi.mandelstam_t>0.1)
    &(phase1_looseChi.chiSq<13.277)
    &(phase1_looseChi.omegaCut)
    &(phase1_looseChi.pVH)
].reset_index(drop=True)

pi0Mean=0.135881
etaMean=0.548625
pi0Std=0.0076
etaStd=0.0191
pi0Sig=3
pi0Skip=1
pi0SB=1.5
etaSig=3
etaSkip=1
etaSB=2

fig,axes=plt.subplots(3,2,figsize=(12,12))
axes=axes.flatten()

def drawSidebands(df,regions,axshift,text):
    for ax, var, pvar, min1, max1, mean, std, i in zip(axes[2*axshift:2*axshift+2],["Mpi0","Meta"],[r'$M(\gamma_1\gamma_2)$',r'$M(\gamma_3\gamma_4)$'],[0.06,0.3],[0.22,0.85],
                                                   [pi0Mean,etaMean],[pi0Std,etaStd],range(2)):
        out=np.histogram(df[var],weights=df['AccWeight'],bins=75)
        hep.histplot(out,ax=ax,c='black')
        ax.set_ylim(0,max(out[0])*1.15)
        ax.set_xlim(min1,max1)
        ax.axvspan(mean-regions[3*i]*std,mean+regions[0]*std,color='green',alpha=0.4)
        ax.axvspan(mean-sum(regions[3*i:3*i+3])*std,mean-sum(regions[3*i:3*i+2])*std,color='red',alpha=0.4)
        ax.axvspan(mean+sum(regions[3*i:3*i+2])*std,mean+sum(regions[3*i:3*i+3])*std,color='red',alpha=0.4)
        ax.text(0.80,0.85,text,transform=ax.transAxes,size=16)
        axes[4+i].set_xlabel(pvar)

def getYieldInRegion(df,regions,text):
    df2=df
    for var,mean,std,i in zip(['Mpi0','Meta'],[pi0Mean,etaMean],[pi0Std,etaStd],range(2)):
        df2=df2[
            ((df2[var]>mean-regions[3*i]*std)&(df2[var]<mean+regions[3*i]*std)) |
            ((df2[var]>mean-sum(regions[3*i:3*i+3]))&(df2[var]<mean-sum(regions[3*i:3*i+2])*std)) |
            ((df2[var]>mean+sum(regions[3*i:3*i+2])*std)&(df2[var]<mean+sum(regions[3*i:3*i+3])*std))
        ]
    print(f'{text} yield: {df2.AccWeight.sum():0.0f}')

regions=[3, 1, 1.5, 3, 1, 2]
drawSidebands(tmp,regions,0,'nominal')
getYieldInRegion(tmp,regions,'nominal')

regions=[2.75, 1.5, 1.25, 2.75, 1.5, 1.75]
drawSidebands(tmp,regions,1,'tighter')
getYieldInRegion(tmp,regions,'tighter')

regions=[3.25, 0.5, 1.75, 3.25, 0.5, 2.25]
drawSidebands(tmp,regions,2,'looser')
getYieldInRegion(tmp,regions,'looser')

plt.savefig("/d/grid17/ln16/dselector_v3/study_pwa/mass_dependent_fits/systOverview_sideband.pdf")















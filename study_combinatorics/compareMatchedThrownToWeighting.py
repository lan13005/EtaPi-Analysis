#!/usr/bin/python

## The purpose of this script is to generate some plots that compares the matched thrown combos
##     to that of the results from weighting. Hopefully if the weighting is good then we have
##     successfully selected good combinations

##########################
#### PRELIMINARIES
##########################
import uproot3 as uproot
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

import mplhep as hep
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])

SMALL_SIZE = 20
MEDIUM_SIZE = 24
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

##########################
##########################
def loadDF(fileName,treeName,cols):
    tree=uproot.open(fileName)[treeName]
    df=tree.arrays(cols,outputtype=pd.DataFrame).reset_index(drop=True)
    return df

selection_vars=["AccWeight","weightBS","isCorrectBeam","isCorrectSpect"]
kinematic_vars=[
        "Mpi0eta",
        #"Mpi0",
        #"Meta",
        "cosTheta_eta_hel",
        "phi_eta_hel",
        "mandelstam_t"]
pretty_vars=[
        r"$M(4\gamma)$ (GeV)",
        #r"$M(\gamma_1\gamma_2)$ (GeV)",
        #r"$M(\gamma_3\gamma_4)$ (GeV)",
        r"$cos_{hel}(\theta)$",
        r"$\phi_{hel}$ (deg)",
        r"-t (GeV^2)"]
df=loadDF("./kmatrix_selected_acc_flat.root","kin",selection_vars+kinematic_vars)

nbins=50
nkinvars=len(kinematic_vars)
fig,axes=plt.subplots(3,nkinvars,figsize=(5*nkinvars,15))
for i,kinematic_var in enumerate(kinematic_vars):
    edges=np.linspace(df[kinematic_var].min(),df[kinematic_var].max(),nbins+1)
    # HOW GOOD IS ACCIDENTAL SUBTRACTION AT SELECTING THE RIGHT BEAM PHOTON?
    hep.histplot(np.histogram(df[kinematic_var],weights=df["AccWeight"],bins=edges),c='red',ax=axes[0,i],linewidth=4,label="Acc.Sub")
    hep.histplot(np.histogram(df.loc[df.isCorrectBeam,kinematic_var],bins=edges),c='black',ax=axes[0,i],linewidth=4,label="Corr.Beam")
    # HOW GOOD IS SIDEBAND SUBTRACTION AT SELECTING THE RIGHT PHOTON PAIRING?
    hep.histplot(np.histogram(df[kinematic_var],weights=df["weightBS"],bins=edges),c='red',ax=axes[1,i],linewidth=4,label="SB.Sub")
    hep.histplot(np.histogram(df.loc[df.isCorrectSpect,kinematic_var],bins=edges),c='black',ax=axes[1,i],linewidth=4,label="Corr.Spect")
    # HOW GOOD IS SIDEBAND SUBTRACTION AT SELECTING THE RIGHT PHOTON PAIRING?
    hep.histplot(np.histogram(df[kinematic_var],weights=df["weightBS"]*df["AccWeight"],bins=edges),c='red',ax=axes[2,i],linewidth=4,label="Acc.Sub\nSB.Sub")
    hep.histplot(np.histogram(df.loc[df.isCorrectBeam&df.isCorrectSpect,kinematic_var],bins=edges),c='black',ax=axes[2,i],linewidth=4,label="Corr.Beam\nCorr.Spect")

    axes[0,i].set_xlabel(pretty_vars[i])
    axes[1,i].set_xlabel(pretty_vars[i])
    axes[2,i].set_xlabel(pretty_vars[i])


axes[0,0].legend()
axes[1,0].legend()
axes[2,0].legend()
plt.tight_layout()
plt.savefig("comparison.pdf")





























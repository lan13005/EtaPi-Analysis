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
import pickle

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

forceGetNewFluxVals=True
target=1.22*1E-9 # its in INVERSE barns units - http://hadron.physics.fsu.edu/Theses/AErnst_FSU_thesis.pdf
BR=0.8429*0.988*0.988; # assumed b1->omegapi dominant=1, omega->pi0g is 8.24, 99% decay to 2g for the 2 pi0s 
binedges=np.linspace(8.2,8.8,2)
bincenters=(binedges[1:]+binedges[:-1])/2

####################################
# Obtain Tagged Flux
####################################
print("\n\nGetting tagged flux\n------------------------")
fluxCounts=[]

runs=["2017"]#,"2018_1","2018_8"]
runStarts=[30274]#,40856,50677]
runEnds=[31057]#,42577,51768]
rcdbQueries=[""]#," --rcdb-query='@is_2018production and @status_approved'"," --rcdb-query='@is_2018production and @status_approved and beam_on_current > 49'"]
for i in range(len(runs)):
    fluxCounts.append([])
    cmd_base="/d/grid13/gluex/gluex_top/hd_utilities/hd_utilities-1.17/psflux/plot_flux_ccdb.py --begin-run="+str(runStarts[i])+" --end-run="+str(runEnds[i])
    cmd_bins="--num-bins="+str(len(binedges)-1)
    cmd_lowE="--energy-min="+str(binedges[0])
    cmd_uppE="--energy-max="+str(binedges[-1])
    cmds=[cmd_base,cmd_bins,cmd_lowE,cmd_uppE]
    cmd=" ".join(cmds)
    cmd+=rcdbQueries[i]
    if not os.path.exists("flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root") or forceGetNewFluxVals:
        print("Running following command:")
        print(cmd)
        os.system(cmd)
    else:
        print("flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root exists already. Lets use this one")

    fluxFile=ROOT.TFile.Open("flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root")
    fluxHist=fluxFile.Get("tagged_flux")
    for j in range(fluxHist.GetNbinsX()):
        count=fluxHist.GetBinContent(j+1)
        print("Bin{0} counts: {1}".format(j,count))
        fluxCounts[i].append(count)
fluxCounts=np.array(fluxCounts)
fluxCounts=fluxCounts.sum(axis=0)

fluxErrors=np.sqrt(fluxCounts)

####################################
# Calculating Cross Section
####################################
print("\n\nGetting Expected Yields\n------------------------")
yields=np.array([134851.0])
yields_err=np.array([20996.0])
fc=fluxCounts
fcerr=fluxErrors
energies=bincenters

crossSections=yields/(fc*target*BR)
crossSectionErrs=crossSections*np.sqrt((yields_err/yields)*(yields_err/yields)+(fcerr/fc)*(fcerr/fc))

print("cross sections: {}".format(crossSections))
print("cross sections errs: {}".format(crossSectionErrs))
    
    















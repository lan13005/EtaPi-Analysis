#!/usr/bin/python

import ROOT
import numpy as np
import os

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

import os
import multiprocessing
import itertools

print("\n\nGetting tagged flux\n------------------------")

pols=["AMO",""]
runs=["2017","2018_1","2018_8"]
runStarts=[30274,40856,50677]
runEnds=[31057,42577,51768]
rcdbQueries=[""," --rcdb-query='@is_2018production and @status_approved'"," --rcdb-query='@is_2018production and @status_approved and beam_on_current > 49'"]
fluxCounts=[]

def getFlux(args):
    pol,i=args
    fluxCounts.append([])
    cmd_base="/d/grid13/gluex/gluex_top/hd_utilities/hd_utilities-1.44/psflux/plot_flux_ccdb.py --begin-run="+str(runStarts[i])+" --end-run="+str(runEnds[i])
    cmd_bins="--num-bins=1"
    cmd_lowE="--energy-min=8.2"
    cmd_uppE="--energy-max=8.8"
    cmd_targetLength="--target-length=26"
    cmds=[cmd_base,cmd_bins,cmd_lowE,cmd_uppE,cmd_targetLength]
    cmd=" ".join(cmds)
    
    if pol!="":
        cmd+=" --pol="+str(pol)
    
    cmd+=rcdbQueries[i]
    print("Running following command:")
    print(cmd)
    os.system(cmd)


args=list(itertools.product(pols, range(len(runs))))
with multiprocessing.Pool(len(args)) as p:
    p.map(getFlux,args)

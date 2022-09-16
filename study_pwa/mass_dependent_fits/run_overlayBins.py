#!/usr/bin/python3

import os
import multiprocessing
import shutil
import subprocess
from checkFits import checkFits

workingDir=os.getcwd()

#ts=["010016", "016021", "021027", "027034", "034042", "042051", "051061", "061072", "072085", "085100"]
#ts=["010020", "0200325", "0325050", "050075", "075100"]
ts=["010020"]

subdir="./"

drawAllGoodFits=True
fitFileName="etapi_result.fit"
doAccCorr="true" # this should generally be true. AccCorr is chosen during the fit, we just want to extract corrected yields for cs measurements
plotAllVars="true" # should we plot all variables in etapi_plotter or just plot the mass plot
plotGenData="false" # should we plot gen mc in etapi_plotter or just plot {dat,bkgnd,accmc} trees
Ls="S_D_pD"
#Ls="S_D"

waves=[
    ############### TMD WITH D-PRIME
    "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++",
    "S0+-;S0++", # individual S waves
    "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
    "D1--_pD1--;D0+-_pD0+-;D1+-_pD1+-;D0++_pD0++;D1++_pD1++;D2++_pD2++", # Merge a2/a2prime with same M and reflectivity
    "pD1--;pD0+-;pD1+-;pD0++;pD1++;pD2++", # individual D prime
    "S0+-_S0++", #sum S waves
    "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
    "pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++", #sum D prime waves
    "D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++", # sum ALL D waves
    "S0++_D0++_D1++_D2++_pD0++_pD1++_pD2++", # all + ref
    "S0+-_D1--_D0+-_D1+-_pD1--_pD0+-_pD1+-" # all - ref
    ############### TMD NO D-PRIME
#    "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++",
#    "S0+-;S0++", # individual S waves
#    "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
#    "S0+-_S0++", #sum S waves
#    "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
#    "S0++_D0++_D1++_D2++", # all + ref
#    "S0+-_D1--_D0+-_D1+-" # all - ref
    ############## KMATRIX 
    #"S0+-_S0++_D0+-_D0++_D2+-_D2++",
    #"S0+-;S0++", # individual S waves
    #"D0+-;D0++;D2+-;D2++", # individual D waves  
    #"S0+-_S0++", #sum S waves
    #"D0+-_D0++_D2+-_D2++", #sum D waves
    #"S0++_D0++_D2++", # all + ref
    #"S0+-_D0+-_D2+-" # all - ref
]
waves=";".join(waves)

def draw(folder):
    if shutil.which("etapi_plotter") is None:
        print("HEY DUMMY! etapi_plotter program does not exist. You probably forgot to source your environment...")
        exit()

    os.chdir(folder)
    folder=workingDir+"/"+folder
    cmd="python3 "+workingDir+"/overlayBins.py 2 '"+waves+"' '"+fitFileName+"' '"+workingDir+"' '"+Ls+"' '"+doAccCorr+"' '"+plotAllVars+"' '"+plotGenData+"' '"+folder+"'"
    print(cmd)
    os.system(cmd)
    os.chdir("..")
    cmd="rsync -a "+folder+" "+"_".join(folder.split("_")[:-1])+"; rm -r "+folder
    print(cmd)
    os.system(cmd)

    return 0

### Make a bunch of folders with only one ".fit" file inside so we can run overlayBins.C in it
folders=[]
convergedFiles=[]
for t in ts:
    fs=[subdir+f for f in os.listdir(subdir) if t in f]
    #fs=[t]
    for f in fs:
        if drawAllGoodFits:
            cf=list(checkFits(f))
            dotIdx=1 if f[0]=="." else 0
            convergedIterations=[c.split(".")[dotIdx].split("result_")[1] for c in cf]
            convergedFiles+=cf
            # do not overwrite (by rerunning) drawing script if the folders already exist
            folders+=[f+"_"+i for i in convergedIterations if not os.path.exists(f+"/"+f+"_"+i)]
        else:
            if not os.path.exists(f+"/etapi_result.fit"):
                print(f'The following results file does not exist! {f+"/etapi_result.fit"}. Exiting...')
                exit(1)
            convergedFiles+=[f+"/etapi_result.fit"]
            folders+=[f+"_0"]

maxProcesses=15 # its about 5GB memory per process if we do not make the genmc plots
print(f"Running etapi_plotter over {len(folders)} files with {maxProcesses} processes...")

[os.system("mkdir -p "+folder) for folder in folders]
[os.system("cp "+convergedFile+" "+folder+"/"+fitFileName) for convergedFile,folder in zip(convergedFiles,folders)]

if len(folders)==1:
    nthreads=2
else:
    if len(folders)<maxProcesses:
        nthreads=len(folders)
    else:
        nthreads=maxProcesses
with multiprocessing.Pool(nthreads) as p:
    ### CHANGING THE WAVESET TO LOOP OVER
    groupVec=waves.split(";")
    groupsWithQuotes=['"_'+tmp+'"' if tmp!="" else '""' for tmp in groupVec]
    repStr="vector<string> groups="
    vecStr="{"
    vecStr+=",".join(groupsWithQuotes)
    vecStr+="}"
    sedArgs=["sed","-i",'s@'+repStr+'.*;@'+repStr+vecStr+';@g',workingDir+"/overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    p.map(draw,folders)

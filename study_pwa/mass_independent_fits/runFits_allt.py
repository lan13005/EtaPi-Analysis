#!/usr/bin/python

import subprocess
import os
import sys

if len(sys.argv)!=4:
    print("Requires 3 arguments")
    print("outFolder: output folder name to dump the results to")
    print("fitName: fit name")
    print("isBS: are you bootstrapping?")
    print('example: ./runFits_allt.py "phase1_m104180_noAccCorr" "EtaPi_fit" "0"')
    exit()

outFolder=sys.argv[1]
fitName=sys.argv[2]
isBS=bool(int(sys.argv[3]))

print("\n-----------------------")
print("outFolder: {}".format(outFolder))
print("fitName: {}".format(fitName))
print("isBS: {}".format(isBS))
print("-----------------------\n")

### This condition should not be used as we will not be able to bootstrap with it
#if os.path.isdir(outFolder):
#    raise ValueError("Output folder already exists! Choose a different name or delete it!")

ts=["010020","0200325","0325050","050075","075100","010020"] # supposed to cycle back to the first element to reset everything
#ts=["010020","075100","010020"] # supposed to cycle back to the first element to reset everything
for i,t in enumerate(ts):
    searchStr='$t=.*'
    replaceStr='$t="'+t+'";'
    sedArgs=["sed","-i",'s@'+searchStr+'@'+replaceStr+'@g','divideData.pl']
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    if i<len(ts)-1: 
        if isBS: 
            os.system("ln -sfn /dev/shm/"+fitName+"_"+t+" "+fitName) # need n flag for no-derference otherwise the symlink doesnt overwrite properly
        else:
            os.system("./divideData.pl")
        os.system("python runFits.py")
        os.system("mkdir -p "+outFolder+"/"+t)
        os.system("mv finalAmps* "+outFolder+"/"+t) # move any finalAmps folder (including bootstrapped) to the output destination
        if not isBS: 
            os.system("mv /dev/shm/"+fitName+" /dev/shm/"+fitName+"_"+t)





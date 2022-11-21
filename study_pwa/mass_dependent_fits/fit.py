#!/usr/bin/python3

import os
import sys
import numpy as np
import subprocess
import uproot3 as uproot
import pandas as pd
import random
import time
from checkCfgFile import checkCfgFile

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def spawnProcessChangeSetting(old,new):
    '''
    Replace the value of varName to varValue in the file called fileName. Depending on the value type we can include quotes or not
    '''
    sedArgs=["sed","-i","s@"+old+"@"+new+"@g","setup_mass_dep_fits.py"]
    print(" ".join(sedArgs))
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()


def replaceStr(search,replace,fileName):
    print("replace str: "+replace)
    sedArgs=["sed","-i",'s@'+search+'@'+replace+'@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

def main():
    ##########################################
    # Basic setup
    ##########################################
    nprocesses=9
    fitFileName="etapi_result.fit"
    niters=[1,3] # [n,m] where n is #restarts and m is #iters per restart. In total we should have n*m iterations 
    workingDir=os.getcwd()

    ts=["010020","0200325","0325050","050075","075100"]
    tmins=[0.1,0.2,0.325,0.5,0.75]
    tmaxs=[0.2,0.325,0.5,0.75,1.0]
    ms=["104172"]
    mmins=[1.04]
    mmaxs=[1.72]
    pcwsBins=[17]
    #ms=["104156","104160","104164","104168","104172","104176","104180"]
    #mmins=[1.04,1.04,1.04,1.04,1.04,1.04,1.04]
    #mmaxs=[1.56,1.60,1.64,1.68,1.72,1.76,1.80]
    #pcwsBins=[13,14,15,16,17,18,19]

    baseCfgFile="config_files/etapi_hybrid_pwave.cfg"
    
    for i in range(niters[0]):
        for j, tmin, tmax in zip(range(len(ts)),tmins,tmaxs):
            t=ts[j]
        
            os.system("rm -f "+fitFileName)
            os.system("rm -f fit*log")
        
            for ic,m,mmin,mmax,pcwsBin in zip(range(len(ms)),ms,mmins,mmaxs,pcwsBins):
                cfgFile=baseCfgFile.split("/")[-1].split(".")[0]+"-copy"
        
                setup_seed=-1
                cmd=f'./setup_mass_dep_fits.py {baseCfgFile} {setup_seed} {pcwsBin} {mmin} {mmax} {tmin} {tmax}'
                print(cmd)
                os.system(cmd) # reinitialize
        
                if not checkCfgFile(cfgFile+".cfg"): 
                    print("\n**** SOMETHING WRONG WITH CONFIG! SEE ABOVE. EXITING ****\n")
                    exit(1)

                ### THIS IS USEFUL IF YOU WANT TO WRITE CFG FILES TO THIS DIRECTORY FOR MALTE ###
#                os.system(f'mv {cfgFile+".cfg"} {"etapi_hybrid_t"+t+"_m"+m+".cfg"}')
                #################################################################################
        
                print("Starting fits")
                cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+cfgFile+".cfg -r "+str(niters[1])+" -m 1000000 -t 1.0 -x 1 -f 0.15" 
                pipeCmd=' > fit.log'
                print(cmd+pipeCmd)
                os.system(cmd+pipeCmd)
               
                # Move results to the desired folder
                ofolder=t+"_"+str(i) 
                os.system("mkdir -p "+ofolder)
                os.system("mv -f etapi_result*.fit "+ofolder+" 2>/dev/null")
                os.system("mv -f "+cfgFile+"*.cfg "+ofolder+" 2>/dev/null")
                os.system("mv -f *.log "+ofolder+" 2>/dev/null")
                os.system("mv -f *.ni "+ofolder+" 2>/dev/null")
                #os.system("ln -snfr rootFiles "+ofolder+"/rootFiles")
        
            _j=j+1 if j!=len(ts)-1 else 0
            spawnProcessChangeSetting(ts[j],ts[_j]) # prepare for new t bin
    

if __name__ == "__main__":
    start_time=time.time()
    main()
    stop_time=time.time()
    print(f"Total elapsed time: {stop_time-start_time}s or {(stop_time-start_time)/60}m")
    
